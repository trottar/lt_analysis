"""Detached, source-bound joint pion MM/delta amplitude study; no ROOT imports."""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import zipfile

import numpy as np
import scipy
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import brentq, lsq_linear

CHANNELS = ('pi_n', 'pi_delta', 'pi_sidis')
SETTINGS = ('Center_highe', 'Left_lowe', 'Center_lowe', 'Left_highe', 'Right_highe')
DELTA = np.array([-10., -7., -4., -1., 2., 5., 8., 11., 14., 17., 20.])


def require(ok, message):
    if not ok:
        raise ValueError('joint_mm_delta:' + message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def arrays(raw):
    with np.load(io.BytesIO(raw), allow_pickle=False) as f:
        return {k: f[k] for k in f.files}


def finite(a, shape=None, nonnegative=False):
    a = np.asarray(a, dtype=float)
    require(np.isfinite(a).all(), 'nonfinite_array')
    if shape is not None:
        require(a.shape == shape, 'array_shape')
    if nonnegative:
        require(np.all(a >= 0), 'negative_array')
    return a


def close(a, b, label):
    require(np.allclose(a, b, rtol=1e-10, atol=1e-12), label)


def solve(X, y, C, target, widths):
    """Existing bounded quadratic least-squares convention with retained priors."""
    X = finite(X)
    y = finite(y, (len(X),))
    C = finite(C, (len(X), len(X)))
    target = finite(target, (X.shape[1],), True)
    widths = np.asarray(widths, dtype=float)
    require(widths.shape == target.shape and np.all(widths > 0), 'prior_widths')
    require(len(X) > 0 and np.allclose(C, C.T, rtol=1e-11, atol=1e-14), 'covariance_shape_or_symmetry')
    L = np.linalg.cholesky(C)
    A = solve_triangular(L, X, lower=True)
    b = solve_triangular(L, y, lower=True)
    active_prior = np.isfinite(widths)
    R = (np.eye(X.shape[1]) / widths[:, None])[active_prior]
    t = (target / widths)[active_prior]
    result = lsq_linear(np.vstack((A, R)), np.r_[b, t], bounds=(0., np.inf), method='trf')
    require(result.success, 'solver_nonconvergence')
    beta = result.x.copy()
    beta[result.active_mask < 0] = 0.
    return dict(beta=beta, active=(result.active_mask < 0), weighted_design=A,
                weighted_target=b, prior_design=R, prior_target=t,
                data_quadratic=float(np.sum((A @ beta - b)**2)),
                prior_quadratic=float(np.sum((R @ beta - t)**2)))


def mc_terms(raw, variance, L):
    raw = finite(raw, nonnegative=True)
    V = finite(variance, raw.shape, True)
    L = finite(L, (raw.shape[1], raw.shape[1]), True)
    Q = raw @ L.T
    # All admitted current parents retain renormalize_shifted_templates=true.
    g = L[1:-1].sum(0)
    close((raw*g).sum(), 1., 'normalization_identity')
    return L, V, Q, (V*g) @ L.T, float((V*g*g).sum()), g


def mc_covariance(terms, beta, take):
    take = np.asarray(take, dtype=int)
    N = terms[0][2].shape[1]
    dd, mm = np.divmod(take, N)
    C = np.zeros((len(take), len(take)))
    for coefficient, (L, V, Q, u, alpha, g) in zip(beta, terms):
        c = np.zeros_like(C)
        for delta in np.unique(dd):
            ids = np.flatnonzero(dd == delta)
            A = L[mm[ids]]
            c[np.ix_(ids, ids)] = (A*V[delta]) @ A.T
        q = Q.ravel()[take]
        uu = u.ravel()[take]
        c += -np.outer(uu, q)-np.outer(q, uu)+alpha*np.outer(q, q)
        C += coefficient**2*c
    return (C+C.T)/2


def integrated_mc_variance(terms, beta, mask):
    variance = 0.
    for b, (L, V, Q, u, alpha, g) in zip(beta, terms):
        derivative = mask @ L-float(np.sum(mask*Q))*g
        variance += b*b*float(np.sum(V*derivative**2))
    return max(0., variance)


def sampling_covariance(X, C, data_variance, answer):
    active = answer['active']
    free = np.flatnonzero(~active)
    A, R = answer['weighted_design'], answer['prior_design']
    H = A.T@A + R.T@R
    inv = np.zeros_like(H)
    if len(free):
        require(np.linalg.matrix_rank(H[np.ix_(free, free)]) == len(free), 'unidentified_free_parameters')
        inv[np.ix_(free, free)] = np.linalg.inv(H[np.ix_(free, free)])
    K = inv @ cho_solve(cho_factor(C, lower=True), X).T
    data = (K*data_variance) @ K.T
    mc = K @ (C-np.diag(data_variance)) @ K.T
    total = (data+mc+data.T+mc.T)/2
    scale = np.sqrt(np.maximum(0, np.diag(total)))
    correlation = np.divide(total, np.outer(scale, scale), out=np.zeros_like(total), where=np.outer(scale, scale)>0)
    return dict(data=data, mc=mc, total=total, curvature=inv, correlation=correlation, response=K)


def template_parameter_derivatives(problem, X, y, C, answer, cov):
    """Differentiate the actual fixed-weight, fixed-prior normal equations.

    Includes the design-transpose times residual term, important with model
    mismatch. Boundary derivatives condition on the returned active set.
    """
    take = problem['take']
    N = problem['y'].shape[1]
    dd, mm = np.divmod(take, N)
    beta = answer['beta']
    wres = cho_solve(cho_factor(C, lower=True), y-X@beta)
    derivatives = []
    covariance = np.zeros((len(beta), len(beta)))
    for ci, (L, V, Q, u, alpha, g) in enumerate(problem['terms']):
        A = np.outer(cov['curvature'][:, ci], wres)-beta[ci]*cov['response']
        D = np.zeros((len(beta), Q.size))
        for delta in np.unique(dd):
            indices = np.flatnonzero(dd == delta)
            D[:, delta*N:(delta+1)*N] = A[:, indices]@L[mm[indices]]
        D -= np.outer(A@Q.ravel()[take], np.tile(g, Q.shape[0]))
        covariance += (D*V.ravel())@D.T
        derivatives.append(D)
    return derivatives, (covariance+covariance.T)/2


def profiles(answer):
    """Delta-penalized-objective=1 contours, not calibrated confidence limits."""
    A, R = answer['weighted_design'], answer['prior_design']
    F = np.vstack((A, R))
    z = np.r_[answer['weighted_target'], answer['prior_target']]
    beta = answer['beta']
    optimum = float(np.sum((F@beta-z)**2))
    result = []
    for j in range(len(beta)):
        other = np.array([k for k in range(len(beta)) if k != j])
        def value(x):
            fit = lsq_linear(F[:, other], z-F[:, j]*x, bounds=(0., np.inf), method='trf')
            require(fit.success, 'profile_nonconvergence')
            return float(np.sum((F[:, other]@fit.x+F[:, j]*x-z)**2))-optimum-1.
        low = 0. if value(0.) <= 0 else brentq(value, 0., beta[j], xtol=1e-12)
        high = beta[j]+max(1e-6, beta[j], 1./max(1e-12, np.linalg.norm(F[:, j])))
        for _ in range(50):
            if value(high) > 0:
                break
            high = 2*high+1e-6
        else:
            raise ValueError('joint_mm_delta:unbounded_profile')
        upper = brentq(value, beta[j], high, xtol=1e-12)
        result.append(dict(lower=float(low), upper=float(upper), lower_is_boundary=bool(low == 0.)))
    return result


def load_inputs(fixed_path, projection_path):
    if __package__:
        from . import extract_h10_delta as extraction
    else:
        import extract_h10_delta as extraction
    fixed_path, projection_path = Path(fixed_path), Path(projection_path)
    with zipfile.ZipFile(fixed_path) as fixed, zipfile.ZipFile(projection_path) as projected:
        require(fixed.testzip() is None and projected.testzip() is None, 'zip_crc')
        require(len(fixed.namelist()) == len(set(fixed.namelist())) and len(projected.namelist()) == len(set(projected.namelist())), 'duplicate_members')
        manifest = json.loads(fixed.read('members-sha256.json'))
        for n, meta in manifest.items():
            b = fixed.read(n)
            require((len(b), digest(b)) == (meta['bytes'], meta['sha256']), 'fixed_member_hash')
        verification = json.loads(fixed.read('verification.json'))
        projection_digest = digest(projection_path.read_bytes())
        require(projection_digest == verification['projection']['sha256'], 'projection_hash')
        receipt = json.loads(projected.read('receipt.json'))
        require(receipt['status'] == 'verified_projection' and receipt['parent_channel_closures'] == 45, 'projection_status')
        require(receipt['primary_reference_sha256'] == verification['historical_acquisition']['sha256'], 'reference_binding')
        require(receipt['runtime_head'] == verification['head'], 'source_binding')
        close(receipt['delta_edges'], DELTA, 'delta_edges')
        observed = arrays(fixed.read('numerical_inputs.npz'))
        extracts = json.loads(fixed.read('evidence/frozen_fit_input_extracts.json'))
        models, roles = {}, {}
        for role in receipt['roles']:
            rid = role['id']
            require(rid not in roles, 'duplicate_role')
            roles[rid] = role
            raw = projected.read(rid+'.npz')
            meta = receipt['members'][rid+'.npz']
            require((len(raw), digest(raw)) == (meta['bytes'], meta['sha256']), 'projected_member_hash')
            a = arrays(raw)
            check, parts = extraction.check_role(a, role['closure_reference'])
            require(check == role['parents'], 'projection_receipt_closure')
            for k, v in parts.items():
                close(v, a[k], 'projected_array_closure')
            models[rid] = a
        require(set(roles) == {s+'_'+c for s in SETTINGS for c in CHANNELS}, 'five_setting_roles')
    docs = {}
    for doc in extracts:
        filename = doc['archive_member'].split('/')[-1]
        setting = filename.split('_kaon_')[0]+'_'+('highe' if '_highe_' in filename else 'lowe')
        require(setting not in docs and setting in SETTINGS, 'fit_setting')
        docs[setting] = doc
    require(set(docs) == set(SETTINGS), 'missing_fit_extract')
    return observed, models, roles, docs, dict(fixed_archive_sha256=digest(fixed_path.read_bytes()), projection_sha256=projection_digest,
        projection_runtime_head=receipt['runtime_head'], historical_reference_sha256=receipt['primary_reference_sha256'],
        current_implementation_sha256=digest(Path(__file__).read_bytes()), cohort_kind=receipt['cohort_kind'])


def build_problem(setting, j, observed, models, roles, docs):
    prefix = f'{setting}/t{j+1}'
    fit = docs[setting]['parents'][j]
    dx = fit['pion_diagnostics']
    joint = dx['joint_refinement']
    B = finite(observed[prefix+'/B'], (3,), True)
    close(B, [joint['coefficients'][c] for c in CHANNELS], 'applied_B_identity')
    target = finite([dx['staged_component_amplitudes'][c] for c in CHANNELS], (3,), True)
    widths = finite([joint['regularization_widths'][c] for c in CHANNELS], (3,), True)
    require(np.all(widths > 0) and joint['regularization_enabled'], 'original_regularization')
    edges = finite(observed[prefix+'/mm_edges'])
    require(np.all(np.diff(edges) > 0), 'mm_edges')
    N = len(edges)+1
    y = finite(observed[prefix+'/nommcuts/data'], (11, N))
    variance = finite(observed[prefix+'/nommcuts/sumw2'], (11, N), True)
    count = finite(observed[prefix+'/nommcuts/count'], (11, N), True)
    require(np.all(count == np.floor(count)), 'noninteger_data_counts')
    require(np.all((count > 0) == (variance > 0)), 'empty_variance_classification')
    for field in ('data', 'sumw2', 'count'):
        close(sum(observed[f'{prefix}/nommcuts/{s}/{field}'] for s in ('prompt', 'rand', 'dummy', 'dummy_rand')),
              observed[f'{prefix}/nommcuts/{field}'], 'signed_source_closure')
    terms = []
    parent_refs = []
    for channel in CHANNELS:
        rid = setting+'_'+channel
        a = models[rid]
        ref = roles[rid]['closure_reference']['parents'][j]
        require(ref['renormalize'] is True, 'unsupported_normalization')
        for key in ('raw_sumw', 'raw_sumw2', 'count', 'alignment_operator', 'aligned'):
            close(observed[f'{prefix}/{channel}/{key}'], a[f't{j+1}_{key}'], 'fixed_projection_binding')
        terms.append(mc_terms(a[f't{j+1}_raw_sumw'], a[f't{j+1}_raw_sumw2'], a[f't{j+1}_alignment_operator']))
        parent_refs.append(ref)
    Xfull = np.stack([t[2] for t in terms], axis=-1)
    domain = np.zeros(N, bool)
    ix = np.asarray(joint['fit_bin_indices'], dtype=int)
    require(len(ix) > 3 and len(set(ix)) == len(ix) and np.all((ix >= 1)&(ix < N-1)), 'original_fit_domain')
    domain[ix] = True
    take = np.flatnonzero(((variance > 0)&domain[None, :])[:10])
    require(len(take) > 3, 'no_supported_cells')
    frozen = fit['frozen_inputs']['pion']
    close(y.sum(0), frozen['values'], 'frozen_input_content')
    close(variance.sum(0), frozen['stored_sumw2'], 'frozen_input_variance')
    # Replay the old marginal fit before using a new joint dimension.
    # Historical replay uses the historical fit's shapes, not regenerated Q.
    # New cohorts may legitimately change histogram contents as well as hashes.
    historical_X = np.column_stack([finite(fit['aligned_shapes'][c]['values'], (N,), True) for c in CHANNELS])
    for c in CHANNELS:
        close(fit['aligned_shapes'][c]['edges'], edges, 'historical_fit_edges')
    old = solve(historical_X[ix], y.sum(0)[ix], np.diag(variance.sum(0)[ix]), target, widths)
    require(np.allclose(old['beta'], B, rtol=1e-6, atol=1e-7), 'original_fit_replay')
    close(old['data_quadratic'], joint['data_chi2_contribution'], 'original_data_objective')
    return dict(setting=setting, parent=j+1, prefix=prefix, edges=edges, B=B, target=target, widths=widths,
        y=y, variance=variance, count=count, Xfull=Xfull, take=take, terms=terms, parent_refs=parent_refs,
        alignment=fit['pion_component_alignment'], original_replay_max_error=float(np.max(abs(old['beta']-B))),
        raw_parts=[models[setting+'_'+ch][f't{j+1}_raw_sumw'] for ch in CHANNELS])


def regional_summary(problem, beta, covariance, fitted, counting_covariance, mc_cov, response=None, mc_derivatives=None):
    y, variance, X = problem['y'], problem['variance'], problem['Xfull']
    prediction = X@beta
    N = y.shape[1]
    ed = problem['edges']
    selections = [('broad', np.r_[False, np.ones(N-2, bool), False]),
                  ('lambda_complete', np.r_[False, (ed[:-1]>=1.1)&(ed[1:]<=1.16), False])]
    result = []
    for name, region in selections:
        mask = np.zeros_like(y, bool)
        mask[:10] = region
        gradient = X[mask].sum(0)
        data = float(y[mask].sum())
        model = float(prediction[mask].sum())
        mc = integrated_mc_variance(problem['terms'], beta, mask)
        residual_data_variance = float(variance[mask].sum())
        residual_MC_variance = mc
        if fitted:
            take = problem['take']
            cross = response@(variance.ravel()[take]*mask.ravel()[take])
            residual_data_variance += float(gradient@counting_covariance@gradient-2*gradient@cross)
            residual_MC_variance = 0.
            for ci, ((L,V,Q,u,alpha,g), D) in enumerate(zip(problem['terms'], mc_derivatives)):
                direct = beta[ci]*(mask@L-float(np.sum(mask*Q))*g).ravel()
                propagated = direct+gradient@D
                residual_MC_variance += float(np.sum(V.ravel()*propagated**2))
        result.append(dict(region=name, data=data, model=model, residual=data-model,
            data_counting_sigma=float(np.sqrt(variance[mask].sum())), fixed_amplitude_MC_sigma=float(np.sqrt(mc)),
            parameter_data_sigma=float(np.sqrt(max(0, gradient@counting_covariance@gradient))),
            parameter_MC_sigma=float(np.sqrt(max(0, gradient@mc_cov@gradient))),
            conditional_residual_counting_sigma=float(np.sqrt(max(0, residual_data_variance))),
            conditional_residual_MC_sigma=float(np.sqrt(max(0, residual_MC_variance))),
            conditional_in_sample_residual_sigma=float(np.sqrt(max(0, residual_data_variance+residual_MC_variance))),
            parameter_uncertainty_note='Estimator derivative only; MC template and parameter terms share inputs and are not independent',
            occupied_cells=int(np.count_nonzero(mask&(variance>0))),
            empty_cell_model=float(prediction[mask&(variance==0)].sum()),
            components=(X[mask].sum(0)*beta).tolist()))
    return result


def fit_parent(problem):
    if __package__:
        from . import extract_h10_delta as extraction
    else:
        import extract_h10_delta as extraction
    take = problem['take']
    N = problem['y'].shape[1]
    X = problem['Xfull'].reshape(-1, 3)[take]
    y = problem['y'].ravel()[take]
    vd = problem['variance'].ravel()[take]
    B, target, widths = (problem[k] for k in ('B', 'target', 'widths'))
    Cmc = mc_covariance(problem['terms'], B, take)
    C = np.diag(vd)+Cmc
    answer = solve(X, y, C, target, widths)
    beta = answer['beta']
    cov = sampling_covariance(X, C, vd, answer)
    mc_derivatives, exact_mc = template_parameter_derivatives(problem, X, y, C, answer, cov)
    cov['mc'] = exact_mc
    cov['total'] = cov['data']+exact_mc
    sig = np.sqrt(np.maximum(0., np.diag(cov['total'])))
    cov['correlation'] = np.divide(cov['total'],np.outer(sig,sig),out=np.zeros_like(cov['total']),where=np.outer(sig,sig)>0)
    original_residual = y-X@B
    candidate_residual = y-X@beta
    fac = cho_factor(C, lower=True)
    # Same precision for fixed/candidate comparisons, including model variations.
    alt = []
    def alternate(label, xx, yy, cc, tar=target, wid=widths, heldout=None):
        a = solve(xx, yy, cc, tar, wid)
        entry = dict(kind=label, beta=a['beta'].tolist(), boundary=a['active'].tolist())
        if heldout is not None:
            omitted = heldout
            entry.update(omitted_data=float(y[omitted].sum()), omitted_fixed_model=float((X[omitted]@B).sum()),
                         omitted_refit_model=float((X[omitted]@a['beta']).sum()))
        alt.append(entry)
    alternate('counting_only', X, y, np.diag(vd))
    alternate('MC_covariance_at_candidate', X, y, np.diag(vd)+mc_covariance(problem['terms'], beta, take))
    alternate('without_regularization', X, y, C, wid=np.full(3, np.inf))
    ed = problem['edges']
    lam = np.r_[False, (ed[:-1]>=1.1)&(ed[1:]<=1.16), False]
    keep = ~lam[take % N]
    if np.count_nonzero(keep)>3:
        alternate('Lambda_excluded', X[keep], y[keep], C[np.ix_(keep, keep)])
    for delta in np.unique(take//N):
        keep = take//N != delta
        if np.count_nonzero(keep)>3:
            alternate('omit_delta_'+str(delta+1), X[keep], y[keep], C[np.ix_(keep, keep)], heldout=~keep)
    for ci, channel in enumerate(CHANNELS):
        info = problem['alignment']['components'][channel]
        shifts = {float(c['total_shift_gev']) for c in info.get('candidate_summaries', [])
                  if not c.get('rejection_reasons') and 'total_shift_gev' in c}
        if info.get('baseline_total_shift_gev') is not None:
            shifts.add(float(info['baseline_total_shift_gev']))
        if shifts:
            for shift in sorted({min(shifts), max(shifts)}):
                if shift == problem['parent_refs'][ci]['shift']:
                    continue
                ref = deepcopy(problem['parent_refs'][ci])
                ref['shift'] = shift
                L = extraction.transform(ref)
                raw = np.asarray(ref['raw']['values'])  # Parent identity retained by transform.
                close((L@raw)[1:-1].sum(), 1., 'alternate_alignment_normalization')
                # Q = S L.T; S is supplied explicitly, never inferred from Sumw2.
                S = problem['raw_parts'][ci]
                Q = S@L.T
                xx = X.copy()
                xx[:, ci] = Q.ravel()[take]
                alternate(f'alignment_{channel}_{shift:.15g}', xx, y, C)
    F = np.vstack((answer['weighted_design'], answer['prior_design']))
    ev = np.linalg.eigvalsh(answer['weighted_design'].T@answer['weighted_design'])
    rank = int(np.count_nonzero(ev > ev[-1]*1e-12)) if ev[-1] > 0 else 0
    delta_variants = np.array([a['beta'] for a in alt if a['kind'].startswith('omit_delta_')])
    alignment_variants = np.array([a['beta'] for a in alt if a['kind'].startswith('alignment_')])
    envelope = lambda a: dict(minimum=np.minimum(a.min(0), beta).tolist(), maximum=np.maximum(a.max(0), beta).tolist()) if len(a) else None
    fixed_regional = regional_summary(problem, B, None, False, np.zeros((3,3)), np.zeros((3,3)))
    fitted_regional = regional_summary(problem, beta, cov['total'], True, cov['data'], cov['mc'], cov['response'], mc_derivatives)
    result = dict(setting=problem['setting'], parent=problem['parent'], original_B=B.tolist(), beta=beta.tolist(),
        staged_target=target.tolist(), retained_prior_widths=widths.tolist(), boundaries=answer['active'].tolist(),
        covariance_counting=cov['data'].tolist(), covariance_MC=cov['mc'].tolist(), covariance_conditional=cov['total'].tolist(),
        penalized_curvature=cov['curvature'].tolist(), correlation=cov['correlation'].tolist(),
        profile_delta_objective_1=profiles(answer), occupied_fit_cells=len(take), data_only_rank=rank,
        data_design_condition=float(np.sqrt(ev[-1]/ev[0])) if rank==3 else None,
        original_marginal_replay_max_error=problem['original_replay_max_error'],
        fixed_data_quadratic=float(original_residual@cho_solve(fac, original_residual)),
        candidate_data_quadratic=float(candidate_residual@cho_solve(fac, candidate_residual)),
        fixed_prior_quadratic=float(np.sum(((B-target)/widths)**2)), candidate_prior_quadratic=answer['prior_quadratic'],
        fixed_regions=fixed_regional, fitted_regions=fitted_regional, sensitivities=alt,
        acceptance_influence_envelope=envelope(delta_variants), alignment_scan_envelope=envelope(alignment_variants),
        model_uncertainty_status='Empirical influence and retained scan envelopes; no calibrated probability law or quadrature combination',
        statistical_status='Conditional signed GLS with retained data-dependent regularization; not independent validation',
        MC_error_method='Derivative of bounded normal equations, fixed active set/precision/prior; includes design-transpose residual term',
        fit_data_covariance_status='No independent old-B/candidate covariance; within-candidate fixed-prior sampling derivative available')
    predictions = dict(fixed=problem['Xfull']@B, candidate=problem['Xfull']@beta,
        components=problem['Xfull']*beta, take=take,
        denominator_inverse_change=np.divide(problem['Xfull'].sum(0)@B, problem['Xfull'].sum(0)@beta,
            out=np.zeros(N), where=(problem['Xfull'].sum(0)@beta)>0),
        denominator_valid=(problem['Xfull'].sum(0)@beta)>0)
    return result, predictions


def run(fixed, projection, output, settings=SETTINGS, parents=(3,1,2)):
    output = Path(output)
    require(not output.exists(), 'output_exists')
    observed, models, roles, docs, provenance = load_inputs(fixed, projection)
    require(len(set(settings)) == len(settings) and set(settings) <= set(SETTINGS) and len(settings)>0, 'setting_selection')
    require(len(set(parents)) == len(parents) and set(parents) <= {1,2,3} and len(parents)>0, 'parent_selection')
    results, predictions = [], {}
    for setting in settings:
        for parent in parents:
            p = build_problem(setting, parent-1, observed, models, roles, docs)
            result, a = fit_parent(p)
            results.append(result)
            for k, v in a.items():
                predictions[f'{setting}/t{parent}/{k}'] = v
            print('DETACHED FIT', setting, 't'+str(parent), result['beta'], flush=True)
    payload = dict(schema='kaonlt_detached_joint_mm_delta/v1', source_provenance=provenance,
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        full_scope=(set(settings)==set(SETTINGS) and set(parents)=={1,2,3}),
        scientifically_accepted=False, production_applied=False, fit_source='signed pion control, retained prior targets/widths',
        primary_covariance='counting + correlated finite-MC, frozen at original B', parents=results,
        limitations=['Same-data development; no independent validation', 'Retained staged priors depend on fitted data; propagated conditional on them',
                     'Alignment and delta-influence envelopes are sensitivity, not physical covariance', 'Upstream cleaning/normalization/systematic covariance unavailable',
                     'Kaon-side numerator, pion weights, Method A/B and Outcome B unchanged', 'Empty-data model mass retained separately; not assigned fictitious variance'])
    encoded = json.dumps(payload, indent=2, allow_nan=False)+'\n'
    output.mkdir(parents=True, exist_ok=False)
    (output/'result.json').write_text(encoded, encoding='utf-8')
    np.savez_compressed(output/'predictions.npz', **predictions)
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixed-comparison', type=Path, required=True)
    parser.add_argument('--projection', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--settings', nargs='+', choices=SETTINGS, default=list(SETTINGS))
    parser.add_argument('--parents', nargs='+', type=int, choices=(1,2,3), default=[3,1,2])
    args = parser.parse_args(argv)
    run(args.fixed_comparison, args.projection, args.output, tuple(args.settings), tuple(args.parents))


if __name__ == '__main__':
    main()
