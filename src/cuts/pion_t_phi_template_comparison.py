"""Opt-in detached child fits using the production component fitter.

No ROOT or science import on import. This module never supplies a production
result. Inherited A/B curves are event-applied backgrounds, not child fits.
"""
from __future__ import annotations

from bisect import bisect_right
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import subprocess

ENABLE_ENV = "LT_ANALYSIS_PION_T_PHI_TEMPLATE_COMPARISON"
SCHEMA = "pion_t_phi_template_comparison/v1"
COMPONENTS = ("pi_n", "pi_delta", "pi_sidis")
SOURCES = ("prompt", "rand", "dummy", "dummy_rand")


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def requested(inp, setting, environ=None):
    if (os.environ if environ is None else environ).get(ENABLE_ENV) != "1":
        return False
    from pion_t_template_comparison import selected
    return ((os.environ if environ is None else environ).get(ENABLE_ENV) == "1"
            and not inp.get("bg_opt_prepass_mode")
            and str(inp.get("analysis_subtraction_context", "")).strip() != "bg_optimizer_pre_parent"
            and str(inp.get("ParticleType", "")).lower() == "kaon"
            and str(inp.get("Q2")) == "4p4" and str(inp.get("W")) == "2p74"
            and selected(inp, setting, environ))


def rows_from_caches(hist):
    """Observe scalar cache rows; preserve each population's own bin convention."""
    nd, nu, nw = (float(hist[k]) for k in ("normfac_data", "normfac_dummy", "nWindows"))
    require(all(math.isfinite(x) and x > 0 for x in (nd, nu, nw)), "source_normalization_invalid")
    coeff = dict(zip(SOURCES, (nd, -nd/nw, -nu, nu/nw)))
    cache = hist["_authoritative_pion_control_source_cache"]
    rows = {"pion": [], "kaon": []}
    for r in cache["records"]:
        require(r["source_label"] in coeff and r["coefficient"] == coeff[r["source_label"]]
                and r["coordinate_fingerprint"] == cache["coordinate_fingerprint"], "pion_row_provenance")
        rows["pion"].append(dict(source=r["source_label"], entry=int(r["entry_index"]),
            mm=float(r["adj_MM"]), ti=int(r["t_index"]), pi=r["phi_index"],
            allcuts=bool(r["allcuts"]), nommcuts=bool(r["nommcuts"]), weight=coeff[r["source_label"]]))
    fields = ("adj_MM", "source_label", "source_entry_index", "t_index", "phi_index",
              "allcuts", "nommcuts", "proton_cleaning_factor", "proton_rf_accept")
    for section in hist["_yield_data_event_cache"].values():
        n = len(section["adj_MM"])
        require(all(len(section[k]) == n for k in fields), "kaon_cache_column_length")
        for i in range(n):
            s = str(section["source_label"][i])
            s = "dummy" if s == "dummy_prompt" else s
            require(s in coeff, "kaon_source_unknown")
            factor = float(section["proton_cleaning_factor"][i])
            require(math.isfinite(factor), "nonfinite_proton_factor")
            rows["kaon"].append(dict(source=s, entry=int(section["source_entry_index"][i]),
                mm=float(section["adj_MM"][i]), ti=int(section["t_index"][i]), pi=int(section["phi_index"][i]),
                allcuts=bool(section["allcuts"][i]), nommcuts=bool(section["nommcuts"][i]),
                weight=coeff[s]*factor if bool(section["proton_rf_accept"][i]) else 0.))
    return rows, coeff


def project(rows, edges, nt, np):
    """Signed content and quadratic Sumw2, with flows/outside-phi/source closure."""
    edges = list(map(float, edges))
    require(len(edges) >= 2 and all(math.isfinite(x) for x in edges)
            and all(a < b for a, b in zip(edges, edges[1:])), "invalid_mm_axis")
    def empty():
        return {"values": [0.]*(len(edges)+1), "variances": [0.]*(len(edges)+1), "count": 0}
    parents = [empty() for _ in range(nt)]
    children = [[empty() for _ in range(np)] for _ in range(nt)]
    outside = [empty() for _ in range(nt)]
    sources = [{s: empty() for s in SOURCES} for _ in range(nt)]
    seen = set()
    for r in rows:
        identity = (r["source"], r["entry"])
        require(identity not in seen, "duplicate_source_entry")
        seen.add(identity)
        require(r["source"] in SOURCES and r["entry"] >= 0 and 0 <= r["ti"] < nt,
                "invalid_source_or_parent")
        require(math.isfinite(r["mm"]) and math.isfinite(r["weight"]), "nonfinite_event")
        pi = r["pi"]
        require(pi is None or pi == -1 or 0 <= pi < np, "invalid_phi_index")
        if not r["selected"]:
            continue
        ti = r["ti"]
        index = bisect_right(edges, r["mm"])  # ROOT flow bins, upper edge is overflow
        for cell in (parents[ti], sources[ti][r["source"]],
                     outside[ti] if pi is None or pi < 0 else children[ti][pi]):
            cell["count"] += 1
            cell["values"][index] += r["weight"]
            cell["variances"][index] += r["weight"]**2
    for ti in range(nt):
        for key in ("values", "variances"):
            for i, value in enumerate(parents[ti][key]):
                require(math.isclose(value, math.fsum(c[key][i] for c in [*children[ti], outside[ti]]),
                                     rel_tol=1e-10, abs_tol=1e-12), "child_closure_failed")
                require(math.isclose(value, math.fsum(c[key][i] for c in sources[ti].values()),
                                     rel_tol=1e-10, abs_tol=1e-12), "source_closure_failed")
    return dict(parents=parents, children=children, outside_phi=outside, sources=sources, closure="PASS", edges=edges)


def make_hist(reference, cell, name):
    h = reference.Clone(name)
    h.SetDirectory(0)
    h.Reset("ICES")
    h.Sumw2()
    for i, (value, variance) in enumerate(zip(cell["values"], cell["variances"])):
        h.SetBinContent(i, value)
        h.SetBinError(i, math.sqrt(variance))
    h.SetEntries(cell["count"])
    return h


def scope_evidence(payload, ti, pi, selected, usable):
    """Mirror the existing resolver's support decisions, not a new fallback."""
    direct = ((payload or {}).get("binned_shapes") or {}).get("t_bin{}".format(ti+1)) or {}
    if usable(direct.get("phi_bin{}".format(pi+1))):
        return "local_t_phi"
    if not usable(selected):
        return "unavailable"
    # Aggregate support is tested by the same resolver, with setting fallback
    # removed. An empty t aggregate cannot be mislabelled by value equality.
    from pion_component_fits import _resolve_single_scope_hist
    probe = dict(payload or {}, setting_shape_full=None)
    resolved = _resolve_single_scope_hist(probe, "scope_evidence", ti, pi)
    return "local_t_aggregate_fallback" if usable(resolved) else "setting_fallback"


def fit_cell(pion, kaon, shapes, inp, *, ti, pi, parent, signal, sigma0,
             sigma0_diagnostics, fitter, phi_interval=None):
    """The existing fitter owns every learned quantity and acceptance policy."""
    def clone(h, name):
        if h is None:
            return None
        result = h.Clone(name)
        result.SetDirectory(0)
        return result
    scope = "t_bin{}phi_bin{}".format(ti+1, pi+1)
    key = {"analysis_scope": scope, "t_bin": {"index": ti, "edges": list(parent["t_edges"])},
           "phi_bin": {"index": pi, "edges": list(phi_interval) if phi_interval is not None else None},
           "epsilon": inp.get("EPSSET"), "phi_setting": parent["fit_result"]["particle_subtraction_phi_setting"],
           "coordinate_fingerprint": parent.get("coordinate_fingerprint")}
    return fitter(clone(pion, scope+"_pion"), clone(kaon, scope+"_kaon"),
        {c: clone(shapes.get(c), scope+"_"+c) for c in COMPONENTS}, deepcopy(inp),
        analysis_scope=scope, kaon_signal_shape=clone(signal, scope+"_lambda"),
        kaon_sigma0_shape=clone(sigma0, scope+"_sigma0"),
        kaon_sigma0_source_diagnostics=deepcopy(sigma0_diagnostics),
        mm_offset_data=parent["fit_result"]["template_mm_offset_data"],
        phi_setting=parent["fit_result"]["particle_subtraction_phi_setting"],
        context="detached_"+scope,
        parent_alignment=deepcopy(parent["fit_result"].get("pion_component_alignment")),
        alignment_bin_key=key)


def apply_transfer(fit, cells, inp):
    """Same helper and configured policy; rejection stays unavailable."""
    from background_config import (resolve_particle_subtraction_weight_denominator_floor,
                                   resolve_particle_subtraction_weight_clip_bounds)
    from pion_component_subtraction import build_simc_shape_pion_control_weights
    from pion_background_chain_audit import snapshot
    low, high = resolve_particle_subtraction_weight_clip_bounds(inp)
    payload = build_simc_shape_pion_control_weights(fit, clip_min=low, clip_max=high,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inp))
    weights = list(map(float, payload["weights"]))
    require(len(weights) == len(cells["nommcuts"]["values"]), "weight_axis_size_mismatch")
    # Existing helper defines zero flow weights. Check rather than reinterpret.
    require(weights[0] == weights[-1] == 0 and all(math.isfinite(x) for x in weights), "invalid_transfer_weights")
    return dict(status="available", policy=snapshot(payload["diagnostics"]), weights=weights,
        signal_window=[float(inp["mm_min"]),float(inp["mm_max"])],
        uncertainty="fixed fitted-weight statistical Sumw2 only; fitted covariance and training variation unavailable",
        cuts={cut: {"values": [v*w for v,w in zip(cell["values"], weights)],
                    "variances": [v*w*w for v,w in zip(cell["variances"], weights)]}
              for cut, cell in cells.items()})


def paired_fixed_weight(applications, cells):
    """Same-event covariance conditional on the fitted weights, not fit errors."""
    available = {k:v for k,v in applications.items() if v.get("status")=="available"}
    pairs = {}
    for i,(a,aa) in enumerate(available.items()):
        for b,bb in list(available.items())[i+1:]:
            pairs[a+"__"+b] = {cut: {
                "covariance": [v*x*y for v,x,y in zip(c["variances"],aa["weights"],bb["weights"])],
                "difference_variance": [v*(x-y)**2 for v,x,y in zip(c["variances"],aa["weights"],bb["weights"])]}
                for cut,c in cells.items()}
    return pairs


def aggregate_applications(cells, method, cut):
    applications = [c["applications"].get(method,{}) for c in cells]
    if not applications or any(a.get("status")!="available" for a in applications):
        return dict(status="unavailable",reason="one_or_more_child_applications_unavailable")
    return dict(status="available",scope="populated canonical phi children only; outside-phi excluded",
        **{k:[math.fsum(a["cuts"][cut][k][i] for a in applications)
              for i in range(len(applications[0]["cuts"][cut][k]))] for k in ("values","variances")})


def observe(hist, inp, t_edges, phi_edges, outpath, signal_payload, sigma0_payload):
    if not requested(inp, hist.get("phi_setting")):
        return None
    from background_config import (BG_OPT_PREPASS_CONTEXT, resolve_particle_subtraction_mode,
                                   resolve_pion_subtraction_scope)
    if (str(inp.get("analysis_subtraction_context", "")).strip() == BG_OPT_PREPASS_CONTEXT
            or resolve_particle_subtraction_mode(inp) != "simc_shape_components"
            or resolve_pion_subtraction_scope(inp) != "t_bin"):
        return None
    import ROOT
    import apply_cuts
    import root_histogram_ownership as ownership
    from pion_background_chain_audit import snapshot, histogram, digest
    from pion_component_fits import (build_particle_subtraction_component_result,
        resolve_scope_component_shapes, resolve_scope_single_shape, _hist_has_usable_support)
    from pion_component_shapes import load_setting_pion_component_shapes, load_raw_t_pion_comparison_shapes
    from hgcer_hole import apply_HGCer_hole_cut
    from pion_t_template_comparison import compare_parent, file_identity, REQUIRED_SOURCE_FILES
    from background_config import get_pion_component_dynamic_alignment_config
    from pion_component_subtraction import evaluate_particle_subtraction_component_fit_result
    setting = hist["phi_setting"]
    epsilon = "highe" if str(inp["EPSSET"]).lower() in ("high", "highe") else "lowe"
    destination = Path(outpath) / (setting+"_kaon_Q4p4W2p74_"+epsilon+"_pion-t-phi-template-comparison.json")
    require(not os.path.lexists(destination), "child_comparison_destination_exists")
    parents = hist["_pion_t_bin_parent_results"]
    require(len(t_edges) == 4 and len(phi_edges) == 10 and len(parents) == 3,
            "child_comparison_canonical_inventory")
    before = digest(parents)
    cut_state = dict(vars(apply_cuts))
    debug_state = ownership._DEBUG_ENABLED, list(ownership._DEBUG_RECORDS)
    directory, pad = ROOT.gDirectory, ROOT.gPad
    result = dict(schema_version=SCHEMA, diagnostic_only=True, production_applied=False,
        scientifically_accepted=False, phi_setting=setting, t_edges=list(map(float,t_edges)),
        signal_window=[float(inp["mm_min"]),float(inp["mm_max"])],
        phi_edges=list(map(float,phi_edges)), status="failed", parents=[], cells=[],
        limitations=["conditional full-sample fits; no held-out evidence or physical truth",
                    "parent A/B control curves have no uniquely defined child projection",
                    "fixed-weight Sumw2 excludes fitted covariance and repeated-training uncertainty"],
        methods={"Global": "historical A: t-fit amplitudes with inherited Global templates",
                 "Local_t": "fresh B: raw-weighted t templates; fixed A alignment",
                 "Local_t_phi": "independent existing fitter and configured child alignment"})
    if epsilon == "highe":
        result["setting"] = {k: inp.get(k) for k in ("Q2", "W", "EPSSET", "ParticleType")}
    try:
        rows, result["source_coefficients"] = rows_from_caches(hist)
        edges = histogram(parents[0]["H_pion_control"])["edges"]
        projections = {pop: {cut: project([dict(r, selected=r[cut]) for r in rr], edges, 3, 9)
                            for cut in ("allcuts", "nommcuts")} for pop,rr in rows.items()}
        result["observations"] = projections
        result["selected_row_fingerprints"] = {p: digest(rr) for p,rr in rows.items()}
        result["coordinate"] = snapshot(hist["_authoritative_pion_control_source_cache"]["kaon_data_coordinate"])
        result["selection_provenance"] = snapshot({k:hist["_authoritative_pion_control_source_cache"].get(k)
            for k in ("coordinate_fingerprint","physical_pion_control_mask_fingerprint","pid_contract","source_accounting")})
        for ti, parent in enumerate(parents):
            require(parent["t_bin_index"] == ti and list(parent["t_edges"]) == list(t_edges[ti:ti+2]), "parent_identity")
            for pop,key in (("pion","H_pion_control"), ("kaon","H_proton_cleaned_final_rf")):
                saved = histogram(parent[key])
                require(saved["edges"] == edges, "parent_axis_mismatch")
                for field, original in (("values",saved["values"]), ("variances",saved["stored_sumw2"])):
                    require(original is not None and len(original)==len(edges)+1 and all(math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-12)
                        for a,b in zip(projections[pop]["nommcuts"]["parents"][ti][field], original)),
                        "parent_observation_closure:"+pop+":"+field)
        result["parent_input_closure"] = "PASS"
        # Hash the actual runtime containers, not historical filenames. These
        # may carry owner-private derived branches; no original-file equivalence
        # is inferred. Source ROOT files are read only.
        result["runtime_inputs"] = {}
        for role,key in (("data","InFile_DATA"),("dummy","InFile_DUMMY")):
            obj = hist[key]
            backing = obj.GetCurrentFile() if hasattr(obj,"GetCurrentFile") else obj
            identity = file_identity(backing.GetName())
            require(identity.get("status")=="available", "runtime_input_identity_unavailable:"+role)
            result["runtime_inputs"][role] = identity
        hole = apply_HGCer_hole_cut(inp["Q2"], inp["W"], inp["EPSSET"], setting)
        shapes = load_setting_pion_component_shapes(deepcopy(inp), setting, particle_type="kaon",
                    t_bins=t_edges, phi_bins=phi_edges, hgcer_cutg=hole, context="detached_t_phi")
        # Loader set_val has temporary module globals. Restore before raw-loader
        # selection: it explicitly uses the active, already initialized cuts.
        vars(apply_cuts).clear()
        vars(apply_cuts).update(cut_state)
        raw = load_raw_t_pion_comparison_shapes(inp, setting, t_edges, phi_edges, hgcer_cutg=hole)
        for ti,parent in enumerate(parents):
            report, local_t = compare_parent(parent, {c:raw[c]["parents"][ti] for c in COMPONENTS}, inp,
                selected_shapes=resolve_scope_component_shapes(shapes, analysis_scope=parent["analysis_scope"],t_bin_index=ti),
                global_shapes=resolve_scope_component_shapes(shapes),
                signal_shape=resolve_scope_single_shape(signal_payload,t_bin_index=ti),
                sigma0_shape=resolve_scope_single_shape(sigma0_payload,t_bin_index=ti),
                sigma0_diagnostics=(sigma0_payload or {}).get("diagnostics") or {},
                fitter=build_particle_subtraction_component_result,
                dynamic_alignment_enabled=bool(get_pion_component_dynamic_alignment_config(inp_dict=inp,phi_setting=setting).get("enabled")))
            result["parents"].append(report)
            for pi in range(9):
                scope = "t_bin{}phi_bin{}".format(ti+1,pi+1)
                cell = dict(t_index=ti,phi_index=pi,t_edges=list(map(float,t_edges[ti:ti+2])),
                    phi_edges=list(map(float,phi_edges[pi:pi+2])), status="fit_failed",fit=None,
                    fit_returned=False, fit_accepted=False, fit_acceptance=None,
                    protected_model_status={"status":"unavailable","reason":"fit_not_returned"}, applications={})
                result["cells"].append(cell)
                inputs = {pop: make_hist(parent["H_pion_control"], projections[pop]["nommcuts"]["children"][ti][pi], scope+pop)
                          for pop in rows}
                cell["inputs"] = {p:histogram(h) for p,h in inputs.items()}
                child_shapes = resolve_scope_component_shapes(shapes,analysis_scope=scope,t_bin_index=ti,phi_bin_index=pi)
                cell["template_scopes"] = {c:scope_evidence(shapes["components"].get(c),ti,pi,child_shapes[c],_hist_has_usable_support) for c in COMPONENTS}
                cell["templates"] = {c:histogram(h) for c,h in child_shapes.items()}
                signal = resolve_scope_single_shape(signal_payload,scope,ti,pi)
                sigma0 = resolve_scope_single_shape(sigma0_payload,scope,ti,pi)
                cell["signal_scopes"] = {name:scope_evidence(payload,ti,pi,chosen,_hist_has_usable_support)
                    for name,payload,chosen in (("k_lambda",signal_payload,signal),("k_sigma0",sigma0_payload,sigma0))}
                cell["protected_signal_templates"] = {"k_lambda":histogram(signal),"k_sigma0":histogram(sigma0)}
                child_fit = None
                try:
                    # A single empty population is the existing fitter's
                    # sparse/zero-solution case; do not impose a new support cut.
                    require(any(any(v>0 for v in cell["inputs"][p]["stored_sumw2"]) for p in rows), "empty_observed_cell")
                    child_fit = fit_cell(inputs["pion"],inputs["kaon"],child_shapes,inp,ti=ti,pi=pi,
                        parent=parent,signal=signal,sigma0=sigma0,
                        sigma0_diagnostics=(sigma0_payload or {}).get("diagnostics") or {},
                        fitter=build_particle_subtraction_component_result,phi_interval=phi_edges[pi:pi+2])
                    cell["fit_returned"] = True
                    cell["fit"] = snapshot(child_fit)
                    cell["protected_model_status"] = snapshot(
                        ((child_fit.get("diagnostics") or {}).get("kaon") or {}).get("pi_delta_signal_protected_fit")
                        or {"status":"unavailable","reason":"protected_diagnostics_not_retained"})
                    evaluation = evaluate_particle_subtraction_component_fit_result(child_fit, inp)
                    cell["fit_acceptance"] = snapshot(evaluation)
                    for pop,key in (("pion","H_pion_control_input"),("kaon","H_kaon_nosub_input")):
                        require(histogram(child_fit.get(key)).get("content_error_fingerprint")
                            == cell["inputs"][pop].get("content_error_fingerprint"),"child_fit_input_changed:"+pop)
                    cell["fit_accepted"] = bool(evaluation["accepted"])
                    if not cell["fit_accepted"]:
                        cell["status"] = "fit_rejected"
                        cell["reason"] = evaluation["reason"]
                    else:
                        cell["status"] = "fitted_local_templates" if all(v=="local_t_phi" for v in cell["template_scopes"].values()) else "fitted_with_template_fallback"
                except Exception as exc:
                    cell["status"] = "empty_child" if str(exc) == "empty_observed_cell" else "fit_failed"
                    cell["reason"] = str(exc)
                for name,fit in (("Global",parent["fit_result"]),("Local_t",local_t),("Local_t_phi",child_fit)):
                    try:
                        require(fit is not None,"fit_unavailable")
                        if name == "Local_t_phi":
                            require(cell["fit_accepted"], "child_fit_not_accepted:" + cell.get("reason", "fit_failed"))
                        cell["applications"][name] = apply_transfer(fit,
                            {cut:projections["pion"][cut]["children"][ti][pi] for cut in ("allcuts","nommcuts")},inp)
                    except Exception as exc:
                        cell["applications"][name] = dict(status="unavailable",reason=str(exc))
                cell["paired_fixed_weight"] = paired_fixed_weight(cell["applications"],
                    {cut:projections["pion"][cut]["children"][ti][pi] for cut in ("allcuts","nommcuts")})
        result["cell_status_counts"] = {s:sum(c["status"]==s for c in result["cells"])
            for s in sorted({c["status"] for c in result["cells"]})}
        result["applied_child_aggregates"] = [{method:{cut:aggregate_applications(
            [c for c in result["cells"] if c["t_index"]==ti],method,cut)
            for cut in ("allcuts","nommcuts")} for method in ("Global","Local_t","Local_t_phi")}
            for ti in range(3)]
        result["status"] = "captured_pending_independent_validation"
    except Exception as exc:
        result["failure_reason"] = str(exc)
        raise
    finally:
        vars(apply_cuts).clear()
        vars(apply_cuts).update(cut_state)
        ownership._DEBUG_ENABLED, ownership._DEBUG_RECORDS[:] = debug_state
        if directory:
            directory.cd()
        if pad:
            pad.cd()
        result["production_parents_unchanged"] = before == digest(parents)
        require(result["production_parents_unchanged"], "production_parents_mutated")
        root = Path(__file__).resolve().parents[2]
        result["source_HEAD"] = subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True).strip()
        files = (*REQUIRED_SOURCE_FILES,"src/cuts/pion_t_phi_template_comparison.py","src/binning/calculate_yield.py",
                 "src/cuts/pion_component_subtraction.py","src/main.py","src/binning/binning_helpers.py")
        result["source_files"] = {p:file_identity(root/p) for p in files}
        result["source_identity_complete"] = all(v.get("status")=="available" for v in result["source_files"].values())
        if not result["source_identity_complete"]:
            result["status"] = "failed"
            result["failure_reason"] = "source_identity_incomplete"
        result["template_source_diagnostics"] = snapshot(shapes.get("diagnostics")) if "shapes" in locals() else None
        with destination.open("x",encoding="utf-8") as stream:
            json.dump(snapshot(result),stream,sort_keys=True,indent=2,allow_nan=False)
            stream.write("\n")
    return result
