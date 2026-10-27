"""Detached conditional statistics and training-dependency audit; no fitting/selection."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def training_dependency_audit(parent):
    """Current packet/AB schemas contain processed full-sample evidence, not raw training inputs."""
    b = parent.get('candidate_B', {})
    observed = {
        'processed_pion_and_kaon_histograms': all(k in parent.get('frozen_inputs', {}) for k in ('pion', 'kaon')),
        'full_sample_alignment': bool(b.get('pion_component_alignment')),
        'full_sample_staged_and_joint_fit': bool(b.get('diagnostics')),
        'full_sample_protected_K_B': 'H_kaon_pion_bg_fit_total_refined' in b,
    }
    missing = [
        'raw learning populations and common physical-event grouping across PID trees',
        'training-only coordinate calibration before processed event capture',
        'training-only proton-cleaning fit and its learning population',
        'training-only setting and local alignment or independently calibrated alignment',
        'training-only staged amplitudes and derived prior widths',
        'training-only protected kaon/signal fits and K_B',
        'training-only local coefficients for every predeclared strength',
    ]
    return {'status': 'unavailable', 'independent_predictive_test_implemented': False,
        'first_missing_step': missing[0], 'missing_training_dependencies': missing,
        'observed_full_sample_fields': observed, 'disjoint_record_folds_are_sufficient': False,
        'established_fitter': {'source': 'src/cuts/pion_component_fits.py::build_particle_subtraction_component_result',
            'accepts_detached_training_histograms': True,
            'complete_training_chain_from_current_packet': False,
            'ROOT_training_refits_executed': False},
        'SIMC': 'independent selected-channel shapes may be reused with independently established calibration; no template uncertainty ensemble',
        'fitted_parameter_covariance': 'unavailable: upstream learned stages/priors and protected models not refitted',
        'repeated_training_variability': 'unavailable', 'physical_bias': 'bias_truth_unavailable',
        'MSE': 'MSE_selection_unavailable', 'selected_lambda': None,
        'truth_audit': 'selected pion SIMC is not independently truth-tagged kaon-PID leakage with verified normalization'}


def conditional_covariance(sumw2, edges, weights, domains):
    """Exact fixed-weight covariance of the same signed independent-count observations."""
    require(isinstance(edges, (list, tuple)) and len(edges) >= 2 and all(map(finite, edges)) and
            all(a < b for a, b in zip(edges, edges[1:])), 'invalid_covariance_axis')
    n = len(edges)+1
    require(isinstance(sumw2, (list, tuple)) and len(sumw2) == n and
            all(finite(v) and v >= 0 for v in sumw2), 'invalid_covariance_Sumw2')
    require(isinstance(weights, dict) and weights and 'B' in weights, 'B_reference_required')
    require(all(isinstance(w, (list, tuple)) and len(w) == n and all(map(finite, w)) and
                w[0] == w[-1] == 0 for w in weights.values()), 'invalid_covariance_weights_or_flows')
    names = list(weights)
    result = {'model_order': names, 'status': 'available', 'scope': 'conditional_fixed_weight',
        'assumption': 'same independent-count observations; source coefficients and all transfer weights fixed',
        'interpretation': 'frozen-weight plug-in count covariance, not the conditional or repeated-sample uncertainty of fitted estimators',
        'fitted_parameter_covariance': 'unavailable', 'repeated_training_variance': 'unavailable',
        'independent_model_validation': False, 'domains': {}}
    for name, window in domains.items():
        require(window is None or (len(window) == 2 and all(map(finite, window)) and window[0] <= window[1]),
                'invalid_covariance_domain')
        ix = [i+1 for i, (a, b) in enumerate(zip(edges, edges[1:]))
              if window is None or window[0] <= (a+b)/2 <= window[1]]
        matrix = [[math.fsum(sumw2[i]*weights[a][i]*weights[b][i] for i in ix)
                   for b in names] for a in names]
        require(all(finite(v) for row in matrix for v in row), 'nonfinite_covariance')
        ib = names.index('B')
        paired = [math.fsum(sumw2[i]*(weights[a][i]-weights['B'][i])*(weights[a][i]-weights['B'][i]) for i in ix) for a in names]
        naive = [matrix[j][j]+matrix[ib][ib] for j in range(len(names))]
        residual = [math.fsum((matrix[j][j], matrix[ib][ib], -2*matrix[j][ib]))-paired[j]
                    for j in range(len(names))]
        require(all(abs(v) <= 1e-11*max(naive[j], 1e-300) for j, v in enumerate(residual)),
                'paired_covariance_identity')
        correlation = [[None if matrix[i][i] == 0 or matrix[j][j] == 0 else
                        matrix[i][j]/(math.sqrt(matrix[i][i])*math.sqrt(matrix[j][j]))
                        for j in range(len(names))] for i in range(len(names))]
        result['domains'][name] = {'covariance': matrix, 'correlation': correlation,
            'paired_difference_variance_from_B': paired,
            'naive_independent_variance_from_B': naive,
            'paired_identity_roundoff': residual,
            'paired_variance_method': 'direct squared difference weights; no cancellation/clipping',
            'window_policy': 'inclusive normal-bin centers, no width weighting; flows excluded'}
    return result


def event_replay_statistics(record, domains):
    """Reuse existing signed spectra; do not reread/select events or change applications."""
    result = {}
    for label, weight_key, app_key in (('ratio', 'weights', 'application'),
            ('production_convention', 'production_convention_weights', 'production_convention_application')):
        weights = {name: model[weight_key] for name, model in record['models'].items()}
        result[label] = {}
        for cut in ('nommcuts', 'allcuts'):
            base = record['populations']['pion'][cut]
            parent = base['parent']['sumw2']
            for parts in (list(base['sources'].values()), [*base['children'].values(), base['outside_phi']]):
                require(all(len(p['sumw2']) == len(parent) for p in parts), 'covariance_partition_geometry')
                require(all(math.isclose(v, math.fsum(p['sumw2'][i] for p in parts), rel_tol=1e-11, abs_tol=1e-30)
                            for i, v in enumerate(parent)), 'covariance_source_or_phi_closure')
            stats = conditional_covariance(parent, record['edges'], weights, domains)
            for domain, values in stats['domains'].items():
                for i, name in enumerate(stats['model_order']):
                    model = record['models'][name]
                    expected = model[app_key][cut]['parent']['domains'][domain]['sumw2']
                    require(math.isclose(values['covariance'][i][i], expected, rel_tol=1e-10, abs_tol=1e-30),
                            'covariance_diagonal_application_identity')
                    if 'paired_difference_from_B' in model:
                        expected = model['paired_difference_from_B'][label][cut]['parent']['domains'][domain]['sumw2']
                        require(math.isclose(values['paired_difference_variance_from_B'][i], expected,
                                             rel_tol=1e-10, abs_tol=1e-30), 'covariance_paired_application_identity')
            stats['source_and_phi_closure'] = 'PASS: covariance linear in verified base Sumw2 partitions'
            result[label][cut] = stats
    return result


def audit_existing(documents, provenance, candidate, source_commit):
    """Use already calculated C; validate exact A/B binding and never invoke an optimizer."""
    try:
        from . import pion_local_t_transfer_variance_comparison as c
        from . import pion_local_t_event_validation_replay as replay
    except ImportError:
        import pion_local_t_transfer_variance_comparison as c
        import pion_local_t_event_validation_replay as replay
    require(candidate.get('schema_version') == c.SCHEMA and candidate.get('status') == 'available' and
            candidate.get('source_HEAD_of_input') == source_commit and candidate.get('inputs_unchanged') is True and
            candidate.get('diagnostic_only') is True and candidate.get('production_applied') is False and
            candidate.get('scientifically_accepted') is False and candidate.get('policy') == json.loads(json.dumps(c.POLICY)),
            'candidate_schema_source_policy')
    require(candidate.get('implementation_sha256') == hashlib.sha256(Path(c.__file__).read_bytes()).hexdigest(),
            'candidate_implementation_identity')
    require(set(documents) == set(c.SETTINGS), 'both_settings_required')
    for setting in c.SETTINGS:
        c.validate_document(documents[setting], setting, matched_head=source_commit)
        require(candidate['provenance'][setting]['sha256'] == provenance[setting]['sha256'], 'candidate_AB_hash_mismatch')
    parents = candidate.get('parents', [])
    require(len(parents) == 6 and {(p['setting'], p['t_bin']) for p in parents} ==
            {(s, t) for s in c.SETTINGS for t in (1, 2, 3)}, 'six_parents_required')
    output = {'schema_version': 'pion_local_t_statistical_validation/v1', 'diagnostic_only': True,
        'production_applied': False, 'scientifically_accepted': False, 'selected_lambda': None,
        'source_HEAD': source_commit, 'scope': 'histogram_only_conditional_uncertainty',
        'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_fitting_performed': False,
        'event_allcuts_validation': 'unavailable: no matched signed-event packets supplied',
        'provenance': deepcopy(provenance), 'parents': []}
    for row in parents:
        p = documents[row['setting']]['parents'][row['t_bin']-1]
        require(row['A_B_unchanged'] == p, 'candidate_saved_AB_mismatch')
        q = c.problem(p)  # Array/config validation only; no least-squares call.
        require(row['fixed_K_B'] == p['candidate_B']['H_kaon_pion_bg_fit_total_refined']['values'], 'fixed_K_B_changed')
        fits = row['C']
        require([f['lambda'] for f in fits] == list(c.LAMBDAS) and
                all(f['status'] == 'available' and f['selected'] is False and f['scientifically_accepted'] is False for f in fits),
                'unselected_seven_strengths_required')
        weights = {'A': row['A_transfer']['weights'], 'B': row['B_transfer']['weights']}
        variants = [(weights['A'], p['baseline_A']['H_pion_fit_total_refined']['values'],
                     p['baseline_A']['H_kaon_pion_bg_fit_total_refined']['values']),
                    (weights['B'], p['candidate_B']['H_pion_fit_total_refined']['values'], row['fixed_K_B'])]
        for fit in fits:
            weights['C_'+format(fit['lambda'], 'g')] = fit['transfer']['weights']
            variants.append((fit['transfer']['weights'], fit['model'], row['fixed_K_B']))
        n = len(q['edges'])+1
        for w, model, numerator in variants:
            require(len(w) == len(model) == len(numerator) == n, 'model_covariance_geometry')
            for i in range(1, n-1):
                require(all(finite(v) and v >= 0 for v in (model[i], numerator[i], w[i])) and
                        (numerator[i] == 0 or model[i] > 0), 'unsupported_covariance_transfer')
                expected = numerator[i]/model[i] if model[i] > 0 else 0.
                require(math.isclose(w[i], expected, rel_tol=1e-12, abs_tol=1e-14), 'candidate_weight_ratio_mismatch')
        h = p['frozen_inputs']['pion']
        output['parents'].append({'setting': row['setting'], 't_bin': row['t_bin'],
            'training_dependencies': training_dependency_audit(p),
            'conditional_statistics': conditional_covariance(h['stored_sumw2'], h['edges'], weights, replay.DOMAINS)})
    return output


def render_pdf(result, path):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(path) as pdf:
        for parent in result['parents']:
            stats = parent['conditional_statistics']; names = stats['model_order']
            fig, axes = plt.subplots(2, 2, figsize=(13, 10))
            fig.subplots_adjust(left=.08, right=.96, bottom=.20, top=.88, wspace=.34, hspace=.48)
            fig.suptitle('{} lowe t{} | frozen-weight covariance on genuine stored histograms\nNO independent model validation; NO selected strength'.format(parent['setting'], parent['t_bin']))
            for col, domain in enumerate(('normal_full_MM', 'signal_1.10_1.16')):
                d = stats['domains'][domain]
                a = np.array([[np.nan if v is None else v for v in row] for row in d['correlation']])
                image = axes[0, col].imshow(a, vmin=-1, vmax=1, cmap='coolwarm')
                axes[0, col].set(xticks=range(len(names)), yticks=range(len(names)), xticklabels=names,
                    yticklabels=names, title=domain+' | fixed-weight correlation')
                axes[0, col].tick_params(axis='x', labelrotation=60, labelsize=8)
                axes[0, col].tick_params(axis='y', labelsize=8)
                fig.colorbar(image, ax=axes[0, col], shrink=.85)
                axes[1, col].plot(names, np.sqrt(d['naive_independent_variance_from_B']), 'o--', color='gray',
                                 label='incorrect independence assumption (reference only)')
                axes[1, col].plot(names, np.sqrt(d['paired_difference_variance_from_B']), 'o-', color='blue',
                                 label='same-observation paired error')
                axes[1, col].set(ylabel='fixed-weight sigma of background difference from B')
                axes[1, col].tick_params(axis='x', labelrotation=60, labelsize=8)
            handles, labels = axes[1, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .06), ncol=2, fontsize=8)
            fig.text(.5, .015, 'Histogram-only; allcuts/event evidence unavailable. No fitted-parameter covariance, repeated-training variance or bias/MSE truth.',
                     ha='center', fontsize=8)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)


def main(argv=None):
    try:
        from . import pion_local_t_transfer_variance_comparison as c
        from . import pion_local_t_event_validation_replay as replay
    except ImportError:
        import pion_local_t_transfer_variance_comparison as c
        import pion_local_t_event_validation_replay as replay
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('center-json', 'left-json', 'candidate-json', 'output-json', 'output-pdf'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--candidate-sha256', required=True, help='Independently reviewed hash of the already computed C JSON; no new fit')
    args = parser.parse_args(argv)
    inputs = [args.center_json, args.left_json, args.candidate_json]
    outputs = [args.output_json, args.output_pdf]
    root = Path(__file__).resolve().parents[1]
    require(args.output_json.suffix == '.json' and args.output_pdf.suffix == '.pdf', 'output_extensions')
    require(len({p.resolve() for p in outputs}) == 2 and not {p.resolve() for p in inputs}&{p.resolve() for p in outputs}, 'input_output_collision')
    require(all(p.parent.is_dir() and not p.exists() and not p.is_symlink() and not p.resolve().is_relative_to(root) for p in outputs), 'new_external_outputs_required')
    raw = {str(p): p.read_bytes() for p in inputs}
    require(hashlib.sha256(raw[str(args.candidate_json)]).hexdigest() == args.candidate_sha256,
            'reviewed_candidate_hash_mismatch')
    docs, provenance = c.load_inputs(center_json=args.center_json, left_json=args.left_json, source_commit=args.source_commit)
    result = audit_existing(docs, provenance, replay.strict_json(raw[str(args.candidate_json)]), args.source_commit)
    result['input_sha256'] = {p: hashlib.sha256(value).hexdigest() for p, value in raw.items()}
    render_pdf(result, args.output_pdf)
    require(all(p.read_bytes() == raw[str(p)] for p in inputs), 'input_mutated_during_audit')
    result['pdf'] = {'status': 'rendered_locally', 'sha256': hashlib.sha256(args.output_pdf.read_bytes()).hexdigest(),
                     'expected_pages': 6, 'ROOT_or_farm_validated': False}
    with args.output_json.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')
    print(json.dumps({'parents': len(result['parents']), 'held_out_status': 'unavailable',
                      'conditional_covariance': 'available', 'selected_lambda': None}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
