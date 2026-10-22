"""Detached histogram-level template A/B/C feasibility; never production fitting.

CLI: python -B testing/pion_local_t_transfer_variance_comparison.py --archive ZIP
     --output-json NEW.json [--output-pdf NEW.pdf]
Alternatively supply exactly --center-json and --left-json. No ROOT imports,
event loops, farm execution, production imports, clipping or strength selection.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import platform
import re
import struct
import sys
import zipfile

import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import chi2

SOURCE = "6d13abb4681b938c56def2a1f455244b0880096c"
COMPONENTS = ("pi_n", "pi_sidis", "pi_delta")
SETTINGS = ("Center", "Left")
LAMBDAS = (0., .01, .1, 1., 10., 100., 1000.)
SCHEMA = "pion_local_t_transfer_variance_comparison/v1"
DOMAINS = {"normal_full_MM": None, "full_0.7_1.3": (.7, 1.3),
           "signal_1.10_1.16": (1.10, 1.16), "visual_1.105_1.135": (1.105, 1.135)}
POLICY = {
    "objective": "chi2_B_data + staged_prior_B + lambda * V_normal_full_MM(beta)/V_normal_full_MM(B)",
    "risk": "V_D=sum_D stored_sumw2_i*(K_B_i/P_beta_i)^2; K_B fixed protected-applied B",
    "normalization": "B risk over ALL normal bins, computed once before sweep; zero B risk uses 1",
    "support": "K>0 and S2>0 requires finite P>0; no floors, clipping, omitted troublesome bins or zero-cost unsupported bins",
    "lambda_grid": list(LAMBDAS), "selection": "none; predeclared in-sample feasibility only",
    "domains": DOMAINS, "window_rule": "inclusive bin centers, normal bins only; no fractional-bin yield integration",
    "solver": "scipy.optimize.least_squares, trf, nonnegative bounds; zero uses original staged start/default tolerances",
    "B_tolerances": {"coefficient_absolute": 1e-7, "objective_absolute": 1e-7, "objective_relative": 1e-9},
}
LIMITS = [
    "Detached local histogram feasibility, not a farm product or scientific acceptance.",
    "Template A/B/C differs from HGCer Method A/B and fitted A_n/B_n amplitudes.",
    "No independent event replay, held-out partition, RF-restored/all-cuts yield or cross section.",
    "Transfer-applied spectra below are signed no-MM-cut histogram proxies; all-cuts Lambda application unavailable.",
    "Sumw2 is a signed-control statistical proxy, not full event or fitted-coefficient covariance.",
    "Conditional Jacobian covariance/sensitivity is not a confidence interval near bounds; data and fitted parameters overlap.",
    "Historical Center A replay remains prior independent evidence, not newly verified here.",
    "Kaon model and templates are fixed; no signal/SIMC unity, purity, or production preference is inferred.",
]


class InvariantError(ValueError):
    """First failed input, provenance or numerical control invariant."""


def require(ok, reason):
    if not ok:
        raise InvariantError(reason)


def encoded(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")) + "\n").encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def basename(setting):
    return setting + "_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json"


def histogram(h, *, variance=False):
    require(isinstance(h, dict) and h.get("status") == "available", "histogram unavailable")
    e, y, err = (np.asarray(h[k], dtype=float) for k in ("edges", "values", "errors"))
    require(e.ndim == y.ndim == err.ndim == 1 and len(e) >= 2 and len(y) == len(e)+1
            and len(err) == len(y), "histogram bin layout")
    require(np.all(np.isfinite(e)) and np.all(np.diff(e)>0) and np.all(np.isfinite(y))
            and np.all(np.isfinite(err)) and np.all(err>=0), "histogram numeric validity")
    packed = [float(len(e)-1), e[0], e[-1], *e]
    packed += [v for pair in zip(y, err) for v in pair]
    fingerprint = sha(b"".join(struct.pack("!d", v) for v in packed))
    require(fingerprint == h.get("content_error_fingerprint"), "content/error fingerprint mismatch")
    if variance:
        require(isinstance(h.get("stored_sumw2"), list), "explicit Sumw2 absent")
        v = np.asarray(h["stored_sumw2"], dtype=float)
        require(v.shape == y.shape and np.all(np.isfinite(v)) and np.all(v>=0), "Sumw2 invalid")
        require(np.allclose(err*err, v, atol=1e-18, rtol=1e-10), "error/Sumw2 mismatch")
        return e, y, v
    return e, y


def load_inputs(*, archive=None, center_json=None, left_json=None):
    """Read only explicit local inputs. Never extract a ZIP or resolve farm paths."""
    require(bool(archive) != bool(center_json or left_json), "supply archive OR exact Center/Left JSON pair")
    docs, provenance = {}, {}
    if archive:
        raw = Path(archive).read_bytes()
        provenance["archive_sha256"] = sha(raw)
        with zipfile.ZipFile(archive) as z:
            require(len(z.namelist()) == len(set(z.namelist())), "duplicate archive member")
            manifest_raw = z.read("manifest.json")
            manifest = json.loads(manifest_raw)
            require(manifest.get("complete") is True and not manifest.get("errors")
                    and manifest.get("git_head") == SOURCE
                    and manifest.get("required_analysis_commit") == SOURCE, "archive provenance differs")
            entries = list(manifest.get("global_artifacts", {}).values())
            expected = {"global/"+basename(s).replace(".json", suffix)
                        for s in SETTINGS for suffix in (".json", ".pdf")}
            require({a.get("archive_path") for a in entries} == expected, "four-artifact manifest scope")
            for entry in entries:
                content = z.read(entry["archive_path"])
                require(sha(content) == entry.get("sha256") and len(content) == entry.get("byte_size"),
                        "archive member hash/size mismatch:"+entry["archive_path"])
            provenance["manifest_sha256"] = sha(manifest_raw)
            for s in SETTINGS:
                content = z.read("global/"+basename(s))
                docs[s] = json.loads(content)
                provenance[s] = {"member": "global/"+basename(s), "sha256": sha(content)}
    else:
        require(center_json and left_json, "both Center and Left JSON required")
        for s, path in zip(SETTINGS, (center_json, left_json)):
            require(Path(path).name == basename(s), "input basename/scope mismatch")
            content = Path(path).read_bytes()
            docs[s] = json.loads(content)
            provenance[s] = {"member": basename(s), "sha256": sha(content)}
        provenance["archive_verification"] = "unavailable: explicit JSON pair, no archive supplied"
    for s in SETTINGS:
        validate_document(docs[s], s)
    return docs, provenance


def validate_document(d, setting, *, matched_head=None):
    expected_head = SOURCE if matched_head is None else matched_head
    require(isinstance(expected_head, str) and re.fullmatch('[0-9a-f]{40}', expected_head) is not None and
            d.get("schema_version") == "pion_t_template_comparison/v1" and d.get("source_HEAD") == expected_head,
            "JSON source/schema provenance differs")
    require(d.get("status") == "available" and d.get("diagnostic_only") is True
            and d.get("production_applied") is False and not d.get("source_identity_rejection_reasons"),
            "JSON availability/diagnostic identity")
    require(d.get("phi_setting") == setting and d.get("setting") == {
        "EPSSET": "low", "ParticleType": "kaon", "Q2": "4p4", "W": "2p74"}, "setting scope differs")
    require(len(d.get("parents", [])) == 3 and [p.get("t_bin_index") for p in d["parents"]] == [0, 1, 2],
            "three canonical parents required per setting")
    require(len(d.get("canonical_t_edges", [])) == 4, "canonical t edges absent")
    for i, p in enumerate(d["parents"]):
        require(p.get("t_edges") == d["canonical_t_edges"][i:i+2], "canonical interval mismatch")
        validate_parent(p)


def validate_parent(p):
    require(p.get("status") == "available" and p.get("authoritative_parent_unchanged") is True
            and not p.get("rejection_reasons") and p.get("production_applied") is False,
            "parent unavailable/mutated")
    a, b = p["baseline_A"], p["candidate_B"]
    require(p.get("alignment_mode") == "fixed_parent_alignment_template_shape_test"
            and a.get("pion_component_alignment") == b.get("pion_component_alignment")
            and bool(b.get("pion_component_alignment")), "fixed alignment differs/absent")
    require(a.get("resolved_subtraction_config") == b.get("resolved_subtraction_config"), "fitting config differs")
    for side, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
        frozen = p["frozen_inputs"][side]
        histogram(frozen, variance=True)
        for fit in (a, b):
            histogram(fit[key], variance=True)
            require(fit[key] == frozen, "frozen input content/error/Sumw2 differs")
    for n in COMPONENTS:
        t = p["templates"][n]
        prov = t.get("provenance", {})
        require(t.get("status") == "genuinely_local" and t.get("normalization_count") == 1
                and all(t.get(k) is True for k in ("canonical_selection_closure", "weighted_child_to_parent_closure",
                                                  "production_global_selection_closure"))
                and t.get("t_edges") == p["t_edges"] and prov.get("status") == "available"
                and isinstance(prov.get("sha256"), str) and len(prov["sha256"]) == 64
                and prov.get("tree_name") and not prov.get("missing_branches"), "invalid local template provenance:"+n)
        histogram(t["histogram"], variance=True)
        require(t["histogram"]["values"] == t.get("normalized_values")
                and np.allclose(t["histogram"]["stored_sumw2"], t.get("normalized_variances"), rtol=1e-10, atol=1e-18)
                and abs(t.get("normalized_integral", 0)-1)<1e-10, "raw local template integrity:"+n)
        # Fitting uses these already-aligned columns, NEVER raw normalized_values.
        e, values = histogram(b["H_simc_shape_"+n])
        require(list(e) == b["H_pion_control_input"]["edges"] and np.all(values>=0), "aligned template axis/positivity")


def problem(p):
    b = p["candidate_B"]
    diag = b["diagnostics"]["pion"]
    joint = diag["joint_refinement"]
    require(joint.get("success") is True and joint.get("status") == "success", "B solver did not converge")
    names = tuple(joint["coefficients"])
    require(set(names) == set(COMPONENTS), "B joint component scope")
    e, y, s2 = histogram(b["H_pion_control_input"], variance=True)
    ix = np.asarray(joint["fit_bin_indices"], dtype=int)
    require(len(ix)>3 and len(set(ix)) == len(ix) and np.all(ix>=1) and np.all(ix<len(y)-1)
            and np.all(s2[ix]>0), "B fit-bin support invalid")
    t = np.column_stack([histogram(b["H_simc_shape_"+n])[1] for n in names])
    centers = (e[:-1]+e[1:])/2
    expected_ix = [i+1 for i, x in enumerate(centers) if diag["fit_min"] <= x <= diag["fit_max"]
                   and not any(lo<=x<=hi for lo, hi in diag.get("exclude_windows", [])) and s2[i+1]>0]
    require(ix.tolist() == expected_ix, "B fit-bin list inconsistent with recorded objective")
    initial = np.array([diag["staged_component_amplitudes"][n] for n in names])
    tau = np.array([joint["regularization_widths"][n] for n in names])
    require(joint.get("regularization_enabled") is True and np.all(np.isfinite(initial))
            and np.all(initial>=0) and np.all(np.isfinite(tau)) and np.all(tau>0), "staged/prior support")
    ke, k = histogram(b["H_kaon_pion_bg_fit_total_refined"])
    require(np.array_equal(ke, e) and np.all(k>=0), "protected-applied kaon model invalid")
    protected = sum(histogram(b["H_kaon_fit_"+n+"_scaled_refined"])[1] for n in names)
    require(np.allclose(k, protected, rtol=1e-10, atol=1e-12), "protected-applied component closure")
    return {"names": names, "edges": e, "y": y, "s2": s2, "t": t, "ix": ix,
            "initial": initial, "tau": tau, "k": k, "joint": joint}


def base_residual(q, beta):
    ix = q["ix"]
    return np.r_[(q["y"][ix]-q["t"][ix]@beta)/np.sqrt(q["s2"][ix]),
                 (beta-q["initial"])/q["tau"]]


def reconstruct_b(q):
    """Original deterministic solver settings, original initialization; no C risk."""
    r = least_squares(lambda x: base_residual(q, x), q["initial"], bounds=(np.zeros(3), np.full(3, np.inf)), method="trf")
    require(r.success, "B reconstruction nonconvergence:"+r.message)
    expected = np.array([q["joint"]["coefficients"][n] for n in q["names"]])
    delta = float(np.max(np.abs(r.x-expected)))
    require(delta<=1e-7, "B coefficient mismatch:"+str(delta))
    residual = base_residual(q, r.x)
    terms = {"data_chi2_contribution": float(residual[:-3]@residual[:-3]),
             "regularization_contribution": float(residual[-3:]@residual[-3:])}
    for key, val in terms.items():
        require(np.isclose(val, q["joint"][key], atol=1e-7, rtol=1e-9), "B objective mismatch:"+key)
    require(np.allclose(q["t"]@expected, q["recorded_model"], atol=1e-12, rtol=1e-9)
            if "recorded_model" in q else True, "B model closure mismatch")
    return r, {"status": "PASS", "max_absolute_coefficient_difference": delta, **terms}


def support_risk(k, model, s2, indices):
    """No unsafe bin is replaced by a zero weight; JSON null means unavailable."""
    ix = np.asarray(indices, dtype=int)
    active = ix[(k[ix]>0)&(s2[ix]>0)]
    bad = active[(~np.isfinite(model[active])) | (model[active]<=0)]
    if len(bad):
        return {"status": "unavailable", "unsupported_bins": bad.tolist(), "V": None}
    val = float(np.sum(s2[active]*(k[active]/model[active])**2))
    return {"status": "available" if np.isfinite(val) else "unstable", "unsupported_bins": [],
            "V": val if np.isfinite(val) else None}


def domains(q):
    c = (q["edges"][:-1]+q["edges"][1:])/2
    return {name: np.arange(1, len(c)+1) if bounds is None else np.flatnonzero((c>=bounds[0])&(c<=bounds[1]))+1
            for name, bounds in DOMAINS.items()}


def transfer_summary(q, model, k):
    normal = np.arange(1, len(model)-1)
    risk = support_risk(k, model, q["s2"], normal)
    unsupported = set(risk["unsupported_bins"])
    # A zero numerator gives exactly zero weight, including 0/0. A positive
    # numerator and zero denominator remains null even if no events occupy it.
    w = np.zeros(len(model))
    ok = np.isfinite(model)&(model>0)
    np.divide(k, model, out=w, where=ok)
    w[(k>0)&~ok] = np.nan
    applied = w*q["y"]
    stats = {}
    for name, ix in domains(q).items():
        v = support_risk(k, model, q["s2"], ix)
        supported = np.all(np.isfinite(w[ix]))
        stats[name] = {**v, "bins": ix.tolist(), "fixed_kaon_model_integral": float(np.sum(k[ix])),
                       "modeled_transfer_integral": float(np.sum(w[ix]*model[ix])) if supported else None,
                       "model_ratio_closure_max_absolute": float(np.max(np.abs(w[ix]*model[ix]-k[ix])))
                       if supported and len(ix) else None,
                       "signed_histogram_applied_integral": float(np.sum(applied[ix])) if supported else None,
                       "applied_minus_model": float(np.sum(applied[ix]-k[ix])) if supported else None,
                       "statistical_sigma_proxy": float(np.sqrt(v["V"])) if v["V"] is not None else None,
                       "max_bin_statistical_amplification_w2": float(np.max(w[ix]**2)) if supported and len(ix) else None}
    finite = w[normal][np.isfinite(w[normal])]
    return {"status": risk["status"], "unsupported_bins": sorted(unsupported),
            "unsupported_MM_intervals": [[float(q["edges"][i-1]), float(q["edges"][i])] for i in sorted(unsupported)],
            "domains": stats,
            "undefined_weight_bins": np.flatnonzero(~np.isfinite(w)).tolist(),
            "weight_bin_quantiles_unweighted_not_event_quantiles": dict(zip(("q50", "q90", "q99", "max"),
                np.quantile(finite, [.5, .9, .99, 1]).tolist())) if len(finite) else {},
            "weights": [float(v) if np.isfinite(v) else None for v in w],
            "signed_histogram_applied": [float(v) if np.isfinite(v) else None for v in applied],
            "allcuts_Lambda_application": {"status": "unavailable", "reason": "no all-cuts event replay in input JSON"},
            "flow_policy": "retained in arrays, excluded from risk/objective/integrals"}


def conditional_diagnostics(q, r):
    """Gauss-Newton curvature only; no statistical coverage assertion."""
    jac = np.asarray(r.jac)
    rank = int(np.linalg.matrix_rank(jac))
    if rank != 3:
        return {"status": "unavailable", "reason": "rank deficient curvature", "rank": rank}
    cov = np.linalg.inv(jac.T@jac)
    sig = np.sqrt(np.diag(cov))
    corr = cov/np.outer(sig, sig)
    # Differentiation of the augmented residual equation, conditional on fixed
    # templates, K, prior and risk normalization; shared events not independent.
    dy = -cov@jac[:len(q["ix"])].T/np.sqrt(q["s2"][q["ix"]])[None, :]
    return {"status": "conditional_only", "covariance": cov.tolist(), "correlation": corr.tolist(),
            "sigma": sig.tolist(), "rank": rank, "condition_number": float(np.linalg.cond(jac)),
            "d_beta_d_data_fit_bins": dy.tolist(), "bound_hits": (r.x<=1e-10).tolist(),
            "limits": "fixed-template Gauss-Newton curvature; not CI; especially invalid near bounds"}


def fit_candidate(q, lam, b_control, normalization, *, max_nfev=None):
    require(lam in LAMBDAS, "lambda outside predeclared scan")
    normal = domains(q)["normal_full_MM"]
    active = normal[(q["k"][normal]>0)&(q["s2"][normal]>0)]
    start = q["initial"] if lam == 0 else b_control.x

    def residual(beta):
        base = base_residual(q, beta)
        if lam == 0:
            return base
        p = q["t"]@beta
        require(np.all(p[active]>0) and np.all(np.isfinite(p[active])),
                "unsupported optimizer support bins:"+str(active[p[active]<=0].tolist()))
        return np.r_[base, np.sqrt(lam/normalization)*np.sqrt(q["s2"][active])*q["k"][active]/p[active]]

    try:
        r = least_squares(residual, start, bounds=(np.zeros(3), np.full(3, np.inf)), method="trf", max_nfev=max_nfev)
    except (ValueError, FloatingPointError) as exc:
        return {"lambda": lam, "status": "unavailable", "reason": str(exc), "initial_coefficients": start.tolist()}
    result = {"lambda": lam, "status": "available" if r.success else "unstable", "solver_message": r.message,
              "solver_status": int(r.status), "nfev": r.nfev, "coefficients": dict(zip(q["names"], r.x.tolist())),
              "initial_coefficients": start.tolist(), "scientifically_accepted": False, "selected": False,
              "last_residual_vector": r.fun.tolist()}
    if not r.success:
        result["reason"] = "solver nonconvergence; no subsequent lambda attempted"
        return result
    if lam == 0:
        require(np.max(np.abs(r.x-b_control.x))<=1e-7, "C zero regularization differs from reconstructed B")
    base = base_residual(q, r.x)
    data_chi = float(base[:-3]@base[:-3])
    model = q["t"]@r.x
    transfer = transfer_summary(q, model, q["k"])
    data_ndf = len(q["ix"])-3
    result.update(model=model.tolist(), components={n: (q["t"][:, i]*r.x[i]).tolist() for i, n in enumerate(q["names"])},
                  residual_sigma_fit_bins=base[:-3].tolist(), data_chi2=data_chi, ndf=data_ndf,
                  data_chi2_ndf=data_chi/data_ndf, data_chi2_p_value_descriptive=float(chi2.sf(data_chi, data_ndf)),
                  p_value_limit="data-only descriptive reference; penalized fitted estimator, not calibrated goodness-of-fit",
                  prior_contribution=float(base[-3:]@base[-3:]), risk_normalization=normalization,
                  total_objective=float(r.fun@r.fun), transfer=transfer,
                  parameter_shift_from_B=(r.x-b_control.x).tolist(), conditional=conditional_diagnostics(q, r),
                  shape_fit_worsened=data_chi>q["joint"]["data_chi2_contribution"]+1e-7)
    result["status"] = transfer["status"]
    if "k_components" in q:
        result["component_statistical_proxies"] = {
            n: {domain: support_risk(column, model, q["s2"], ix) for domain, ix in domains(q).items()}
            for n, column in q["k_components"].items()}
        result["component_statistical_proxy_limit"] = "components share control events; per-component risks cannot be summed as independent variances"
    if result["conditional"]["status"] == "conditional_only":
        cov = np.array(result["conditional"]["covariance"])
        component_uncertainty = {}
        for name, ix in domains(q).items():
            positive = ix[model[ix]>0]
            g = -np.sum((q["y"][positive]*q["k"][positive]/model[positive]**2)[:, None]*q["t"][positive], axis=0)
            component_uncertainty[name] = {"conditional_coefficient_sigma_proxy": float(np.sqrt(max(0., g@cov@g))),
                                            "not_combined_with_Sumw2": "shared-data covariance absent"}
        result["conditional_application_uncertainty"] = component_uncertainty
    return result


def matched_fit_identity(fit, authority, provenance, setting, parent, implementation_hash):
    """Cohort-specific attempt identity, never a model-quality certificate."""
    payload = {k: v for k, v in fit.items() if k not in ('fit_identity', 'model_identity')}
    return digest({'cohort': digest(authority), 'AB_sha256': provenance[setting]['sha256'],
                   'setting': setting, 't_bin_index': parent['t_bin_index'],
                   'parent_fingerprint': parent['baseline_parent_fingerprint'],
                   'policy_sha256': digest(POLICY), 'implementation_sha256': implementation_hash, 'fit': payload})


def compare(docs, provenance, *, matched_authority=None):
    """All six numerical controls must pass before ANY positive-lambda fit."""
    before = digest(docs)
    prepared = []
    for s in SETTINGS:
        validate_document(docs[s], s, matched_head=None if matched_authority is None else matched_authority['source_HEAD'])
        for p in docs[s]["parents"]:
            q = problem(p)
            q["k_components"] = {n: histogram(p["candidate_B"]["H_kaon_fit_"+n+"_scaled_refined"])[1]
                                 for n in q["names"]}
            q["recorded_model"] = histogram(p["candidate_B"]["H_pion_fit_total_refined"])[1]
            r, receipt = reconstruct_b(q)
            prepared.append((s, p, q, r, receipt))
    output = {"schema_version": SCHEMA, "status": "available", "diagnostic_only": True,
              "production_applied": False, "scientifically_accepted": False, "provenance": deepcopy(provenance),
              "source_HEAD_of_input": SOURCE, "policy": deepcopy(POLICY), "limits": LIMITS[:],
              "input_metadata": {s: {k: deepcopy(docs[s][k]) for k in (
                  "source_files", "original_analysis_input_identity", "original_input_population_accounting",
                  "canonical_t_edges", "canonical_phi_edges")} for s in SETTINGS},
              "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "implementation_sha256": sha(Path(__file__).read_bytes()), "parents": []}
    if matched_authority is not None:
        require(matched_authority.get('schema_version') == 'pion_matched_cohort_authority/v1'
                and matched_authority.get('mode') == 'matched_cohort'
                and set(matched_authority.get('inputs', {})) == set(SETTINGS)
                and matched_authority.get('private_input_policy') in (None, 'frozen_shiftprep_private_copy'), 'matched_C_authority')
        for s in SETTINGS:
            require(isinstance(provenance[s].get('sha256'), str) and
                    re.fullmatch('[0-9a-f]{64}', provenance[s]['sha256']) is not None, 'matched_C_AB_source_hash')
            require(provenance[s]['canonical_json_sha256'] == digest(docs[s]), 'matched_C_AB_identity')
            expected = matched_authority['inputs'][s]
            if matched_authority.get('private_input_policy') == 'frozen_shiftprep_private_copy':
                require(provenance[s]['original_input_identity'] == expected, 'matched_C_original_input_identity')
                expected = provenance[s]['runtime_input_identity']
            for role in ('data', 'dummy'):
                require(all(docs[s]['original_analysis_input_identity'][role][k] == expected[role][k]
                            for k in ('sha256', 'bytes')), 'matched_C_source_identity')
        output['source_HEAD_of_input'] = matched_authority['source_HEAD']
        output['cohort_authority'] = deepcopy(matched_authority)
    for s, p, q, r, receipt in prepared:
        baseline_model = histogram(p["candidate_B"]["H_pion_fit_total_refined"])[1]
        b_transfer = transfer_summary(q, baseline_model, q["k"])
        risk = b_transfer["domains"]["normal_full_MM"]["V"]
        row = {"setting": s, "t_bin": p["t_bin_index"]+1, "t_edges": p["t_edges"],
               "B_reconstruction": receipt, "A_B_unchanged": deepcopy(p), "B_transfer": b_transfer,
               "A_transfer": transfer_summary(q, histogram(p["baseline_A"]["H_pion_fit_total_refined"])[1],
                                               histogram(p["baseline_A"]["H_kaon_pion_bg_fit_total_refined"])[1]),
               "C": [], "component_order": list(q["names"]), "fixed_K_B": q["k"].tolist()}
        if risk is None:
            row.update(status="unavailable", reason="B unsupported risk normalization", unsupported_bins=b_transfer["unsupported_bins"])
            output["status"] = "unavailable"
        else:
            for lam in LAMBDAS:
                candidate = fit_candidate(q, lam, r, risk if risk>0 else 1.)
                if matched_authority is not None:
                    candidate['fit_identity'] = matched_fit_identity(candidate, matched_authority, provenance,
                                                                     s, p, output['implementation_sha256'])
                    if 'model' in candidate:
                        candidate['model_identity'] = digest([candidate['fit_identity'], candidate['model'], row['fixed_K_B']])
                row["C"].append(candidate)
                if candidate["status"] != "available":
                    output["status"] = candidate["status"]
                    output["first_failed_invariant"] = {"setting": s, "t_bin": row["t_bin"], "attempt": candidate}
                    break  # fail closed, retain attempt and exact first invariant
        output["parents"].append(row)
        if output["status"] != "available":
            break
    require(digest(docs) == before, "A/B input mutation")
    output["inputs_unchanged"] = True
    return output


def render_pdf(result, path):
    """Three pages per parent; all strengths shown, none selected."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    def visible_scale(ax, window, extra=()):
        # One scale shared by every comparator in a panel, calculated from the
        # displayed window so a large outside peak cannot erase the signal zoom.
        values = []
        for line in ax.get_lines():
            xx = np.asarray(line.get_xdata(), dtype=float)
            yy = np.asarray(line.get_ydata(), dtype=float)
            if xx.shape == yy.shape:
                values.extend(yy[(xx>=window[0])&(xx<=window[1])&np.isfinite(yy)].tolist())
        mask = (x>=window[0])&(x<=window[1])
        for array in extra:
            values.extend(np.asarray(array)[mask].tolist())
        values = [v for v in values if np.isfinite(v)]
        if values:
            lo, hi = min(0., min(values)), max(0., max(values))
            pad = max((hi-lo)*.08, 1e-15)
            ax.set_ylim(lo-pad, hi+pad)

    with PdfPages(path, metadata={"Title": "Detached A/B/C feasibility - NOT ACCEPTED", "CreationDate": None}) as pdf:
        for row in result["parents"]:
            p = row["A_B_unchanged"]
            q = problem(p)
            x = (q["edges"][:-1]+q["edges"][1:])/2
            candidates = [c for c in row["C"] if "model" in c]
            colors = plt.cm.viridis(np.linspace(.05, .95, len(candidates)))
            title = f'{row["setting"]}/lowe t{row["t_bin"]}: detached in-sample feasibility; NO selected C'
            for page in ("control", "transfer", "sweep", "components and risk"):
                fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
                fig.suptitle(title+"\n"+page, fontsize=13)
                if page == "control":
                    for col, window in enumerate(((.7, 1.3), (1.105, 1.135))):
                        ax, residual_ax = axes[0, col], axes[1, col]
                        ax.errorbar(x, q["y"][1:-1], np.sqrt(q["s2"][1:-1]), fmt="k.", ms=3, label="signed control +/- Sumw2")
                        for label, fit, color in (("A inherited", p["baseline_A"], "tab:blue"), ("B local", p["candidate_B"], "tab:red")):
                            m = np.array(fit["H_pion_fit_total_refined"]["values"])
                            ax.plot(x, m[1:-1], color=color, lw=2, label=label)
                            v = q["s2"][1:-1]>0
                            residual_ax.plot(x[v], (q["y"][1:-1][v]-m[1:-1][v])/np.sqrt(q["s2"][1:-1][v]), color=color, label=label)
                        for c, color in zip(candidates, colors):
                            ax.plot(x, c["model"][1:-1], color=color, lw=1, label=f'C lambda={c["lambda"]:g}')
                            v = q["s2"][1:-1]>0
                            residual_ax.plot(x[v], (q["y"][1:-1][v]-np.array(c["model"])[1:-1][v])/np.sqrt(q["s2"][1:-1][v]), color=color)
                        ax.set(xlim=window, ylabel="weighted bin content", title="Full MM" if col == 0 else "Lambda visual window")
                        residual_ax.set(xlim=window, xlabel="MM (GeV)", ylabel="(data-model)/sqrt(Sumw2)")
                        residual_ax.axhline(0, color="k", lw=.5)
                        visible_scale(ax, window, (q["y"][1:-1]-np.sqrt(q["s2"][1:-1]),
                                                   q["y"][1:-1]+np.sqrt(q["s2"][1:-1])))
                        visible_scale(residual_ax, window)
                    axes[0, 0].legend(fontsize=7, ncol=2)
                elif page == "transfer":
                    for col, window in enumerate(((.7, 1.3), (1.10, 1.16))):
                        ax, wax = axes[0, col], axes[1, col]
                        h = p["frozen_inputs"]["kaon"]
                        ax.errorbar(x, h["values"][1:-1], h["errors"][1:-1], fmt="k.", ms=3, label="post-proton kaon input")
                        ax.plot(x, q["k"][1:-1], "k--", label="fixed protected K_B")
                        for label, tr, color in (("A histogram proxy", row["A_transfer"], "tab:blue"), ("B histogram proxy", row["B_transfer"], "tab:red")):
                            ax.plot(x, tr["signed_histogram_applied"][1:-1], color=color, label=label)
                            wax.plot(x, tr["weights"][1:-1], color=color, label=label)
                        for c, color in zip(candidates, colors):
                            ax.plot(x, c["transfer"]["signed_histogram_applied"][1:-1], color=color, label=f'C lambda={c["lambda"]:g}')
                            wax.plot(x, c["transfer"]["weights"][1:-1], color=color)
                        ax.set(xlim=window, ylabel="weighted bin content", title="Extended no-MM-cut histogram proxy")
                        wax.set(xlim=window, xlabel="MM (GeV)", ylabel="K/P (no clipping)")
                        wax.set_yscale("symlog", linthresh=.01)
                        visible_scale(wax, window)
                        visible_scale(ax, window, (np.array(h["values"])[1:-1]-np.array(h["errors"])[1:-1],
                                                   np.array(h["values"])[1:-1]+np.array(h["errors"])[1:-1]))
                    axes[0, 0].legend(fontsize=7, ncol=2)
                elif page == "sweep":
                    ax = axes[0, 0]
                    for annotation_index, (c, color) in enumerate(zip(candidates, colors)):
                        risk = c["transfer"]["domains"]["normal_full_MM"]["V"]
                        ax.scatter(c["data_chi2"], risk, color=color)
                        ax.annotate(f'{c["lambda"]:g}', (c["data_chi2"], risk), fontsize=8,
                                    xytext=(8, -12-annotation_index*11 if annotation_index<4 else 8),
                                    textcoords="offset points", arrowprops={"arrowstyle": "-", "lw": .5})
                    ax.set(xlabel="data chi2 (original fit bins)", ylabel="V full normal MM", yscale="symlog", title="Pareto sweep: labels are lambda")
                    for c, color in zip(candidates, colors):
                        axes[0, 1].plot(range(3), [c["coefficients"][n] for n in q["names"]], "o-", color=color, label=f'{c["lambda"]:g}')
                    axes[0, 1].set(xticks=range(3), xticklabels=q["names"], ylabel="pion-control amplitude", title="Sensitivity / bounds (zero)")
                    # Components for EVERY strength; component line styles are fixed.
                    for c, color in zip(candidates, colors):
                        for n, style in zip(q["names"], ("-", "--", ":")):
                            axes[1, 0].plot(x, c["components"][n][1:-1], style, color=color, lw=1)
                    for fit, color in ((p["baseline_A"], "tab:blue"), (p["candidate_B"], "tab:red")):
                        for n, style in zip(q["names"], ("-", "--", ":")):
                            axes[1, 0].plot(x, fit["H_pion_fit_"+n+"_scaled_refined"]["values"][1:-1], style, color=color, lw=1.5)
                    axes[1, 0].set(xlim=(.7, 1.3), xlabel="MM (GeV)", ylabel="component content", title="Components: n solid / sidis dashed / delta dotted")
                    axes[1, 1].axis("off")
                    text = ["C=0 reproduces B: "+row["B_reconstruction"]["status"],
                            "V normalization: fixed full-normal-bin B risk", "All strengths unselected; shape worsening retained.",
                            "All-cuts Lambda yield/application: UNAVAILABLE", "No RF-restored yield, event replay or full covariance.",
                            "P-values are descriptive for penalized estimates.", "Conditional curvature is NOT a confidence interval."]
                    for c in row["C"]:
                        text.append(f'lambda={c["lambda"]:g}: {c["status"]}'+("; shape chi2 worsens" if c.get("shape_fit_worsened") else ""))
                    axes[1, 1].text(.02, .98, "\n".join(text), va="top", fontsize=9)
                else:
                    for col, window in enumerate(((.7, 1.3), (1.105, 1.135))):
                        for label, fit, color in (("A", p["baseline_A"], "tab:blue"), ("B", p["candidate_B"], "tab:red")):
                            for n, style in zip(q["names"], ("-", "--", ":")):
                                axes[0, col].plot(x, fit["H_pion_fit_"+n+"_scaled_refined"]["values"][1:-1], style, color=color, lw=2)
                        for c, color in zip(candidates, colors):
                            for n, style in zip(q["names"], ("-", "--", ":")):
                                axes[0, col].plot(x, c["components"][n][1:-1], style, color=color, lw=1)
                        axes[0, col].set(xlim=window, xlabel="MM (GeV)", ylabel="component content",
                                        title="n solid / sidis dashed / delta dotted; A blue, B red")
                        visible_scale(axes[0, col], window)
                    for ax in axes[1]:
                        ax.axis("off")
                    lines = ["lambda     chi2/ndf    V full      V signal     sigma signal"]
                    for c in candidates:
                        tr = c["transfer"]["domains"]
                        lines.append(f'{c["lambda"]:7g}  {c["data_chi2_ndf"]:10.4g}  '
                                     f'{tr["normal_full_MM"]["V"]:10.4g}  '
                                     f'{tr["signal_1.10_1.16"]["V"]:10.4g}  '
                                     f'{tr["signal_1.10_1.16"]["statistical_sigma_proxy"]:10.4g}')
                    axes[1, 0].text(0, .98, "\n".join(lines), va="top", fontsize=8, family="monospace")
                    notes = ["Statistical sigma: sqrt(V), fixed denominator proxy.",
                             "Bins selected by centers; NOT extracted yields.",
                             "Full MM includes 0.94-1.03 GeV.", "No missing-support clipping; all reported in JSON.",
                             "All fitted covariance/correlations and bound hits:", "retained in JSON, conditional only.",
                             "A/B staged, joint, protected-fit limitations:", "retained unchanged in JSON."]
                    axes[1, 1].text(.02, .98, "\n".join(notes), va="top", fontsize=9)
                fig.supxlabel("LOCAL HISTOGRAM DIAGNOSTIC ONLY. All-cuts signal yield unavailable; no production or farm acceptance.", fontsize=9)
                pdf.savefig(fig)
                plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--center-json", type=Path)
    parser.add_argument("--left-json", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path)
    args = parser.parse_args(argv)
    outputs = [p.resolve() for p in (args.output_json, args.output_pdf) if p]
    inputs = [p.resolve() for p in (args.archive, args.center_json, args.left_json) if p]
    require(len(set(outputs)) == len(outputs) and not set(inputs)&set(outputs), "input/output path collision")
    repo = Path(__file__).resolve().parents[1]
    require(all(not p.is_relative_to(repo) for p in outputs), "diagnostic outputs must stay outside source checkout")
    require(args.output_json.suffix == ".json" and (not args.output_pdf or args.output_pdf.suffix == ".pdf"),
            "output extensions must be JSON/PDF")
    require(all(not p.exists() and p.parent.is_dir() for p in outputs), "outputs must be new files in existing local directories")
    print(json.dumps(POLICY, sort_keys=True), flush=True)  # predeclared BEFORE loading/comparing
    try:
        docs, provenance = load_inputs(archive=args.archive, center_json=args.center_json, left_json=args.left_json)
        result = compare(docs, provenance)
    except (ValueError, TypeError, KeyError, OSError, zipfile.BadZipFile) as exc:
        result = {"schema_version": SCHEMA, "status": "unavailable", "first_failed_invariant": str(exc),
                  "policy": POLICY, "limits": LIMITS, "production_applied": False}
    if args.output_pdf and "parents" in result:
        try:
            render_pdf(result, args.output_pdf)
            result["pdf"] = {"status": "rendered_locally", "sha256": sha(args.output_pdf.read_bytes())}
        except (ImportError, ValueError, OSError, RuntimeError, AttributeError) as exc:
            result["pdf"] = {"status": "unavailable", "reason": str(exc)}
    else:
        result["pdf"] = {"status": "unavailable", "reason": "not requested or numerical gate failed"}
    args.output_json.write_bytes(encoded(result))
    print(json.dumps({"status": result["status"], "output_json_sha256": sha(args.output_json.read_bytes()),
                      "parents": len(result.get("parents", [])), "pdf": result["pdf"]}))
    return 0 if result["status"] == "available" else 2


if __name__ == "__main__":
    sys.exit(main())
