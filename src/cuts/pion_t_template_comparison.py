"""Opt-in fixed-alignment template experiment; never supplies production weights.

Scalar accumulation and evidence serialization need no ROOT. Runtime adapters
are lazy and use the existing loader selection and component fitter.
"""
from __future__ import annotations

import bisect
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from uuid import uuid4

from canonical_binning import find_canonical_bin
from pion_background_chain_audit import digest, histogram, snapshot
import root_histogram_ownership as ownership

ENABLE_ENV = "LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON"
COMPONENTS = ("pi_n", "pi_delta", "pi_sidis")
SCHEMA = "pion_t_template_comparison/v1"
REQUIRED_SOURCE_FILES = (
    "src/cuts/pion_t_template_comparison.py", "src/utility/pion_component_shapes.py",
    "src/cuts/rand_sub.py", "src/cuts/pion_component_fits.py", "src/cuts/pion_t_bin_parents.py",
    "src/cuts/apply_cuts.py", "src/utility/background_config.py", "src/utility/canonical_binning.py",
    "src/cuts/hgcer_hole.py", "src/cuts/data_coordinates.py",
    "src/utility/root_histogram_ownership.py", "run_Prod_Analysis.sh",
)


def enabled(environ=None):
    return (os.environ if environ is None else environ).get(ENABLE_ENV) == "1"


def file_identity(path):
    path = Path(path)
    result = {"path": str(path)}
    try:
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(block)
        result.update(status="available", bytes=path.stat().st_size, sha256=hasher.hexdigest())
    except OSError as exc:
        result.update(status="unavailable", reason=str(exc))
    return result


def _edges(values):
    values = list(map(float, values))
    if len(values) < 2 or not all(map(math.isfinite, values)) or any(
            a >= b for a, b in zip(values, values[1:])):
        raise ValueError("invalid histogram/canonical edges")
    return values


def _cell(n):
    return {"count": 0, "signed": 0., "absolute": 0., "sumw2": 0.,
            "values": [0.] * n, "variances": [0.] * n}


def _fill(cell, index, weight):
    cell["count"] += 1
    cell["signed"] += weight
    cell["absolute"] += abs(weight)
    cell["sumw2"] += weight * weight
    cell["values"][index] += weight
    cell["variances"][index] += weight * weight


def _close(a, b):
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)


def _valid_provenance(provenance):
    sha = provenance.get("sha256") or ""
    return (provenance.get("status") == "available" and len(sha) == 64
            and all(c in "0123456789abcdef" for c in sha)
            and bool(provenance.get("path")) and bool(provenance.get("tree_name"))
            and provenance.get("weight_expression") == "getattr(evt, 'iter_weight', 1.0)"
            and provenance.get("coordinate") == "missmass; -t; phipq*180/pi")


def _closure(total, parts):
    return (all(_close(total[k], math.fsum(p[k] for p in parts))
                for k in ("count", "signed", "absolute", "sumw2"))
            and all(_close(total[k][i], math.fsum(p[k][i] for p in parts))
                    for k in ("values", "variances") for i in range(len(total[k]))))


def integrate_raw_templates(records, t_edges, phi_edges, mm_edges, provenance):
    """Consume selected (raw MM, -t, phi degrees, iter_weight) once.

    Normal MM bins exclude flows. Canonical intervals include the final endpoint.
    Direct t fills and independently filled raw children must close, including
    flows and Sumw2. Children are NEVER normalized to construct a parent.
    """
    t_edges, phi_edges, mm_edges = map(_edges, (t_edges, phi_edges, mm_edges))
    n = len(mm_edges) + 1
    selected, accepted, excluded = (_cell(n) for _ in range(3))
    parents = [_cell(n) for _ in t_edges[:-1]]
    children = [[_cell(n) for _ in phi_edges[:-1]] for _ in parents]
    phi_excluded = [_cell(n) for _ in parents]
    invalid = 0
    for record in records:
        x, t, phi, weight = map(float, record)
        if not all(map(math.isfinite, (x, t, phi, weight))):
            invalid += 1
            continue
        mm_bin = bisect.bisect_right(mm_edges, x)
        _fill(selected, mm_bin, weight)
        ti, pi = find_canonical_bin(t, t_edges), find_canonical_bin(phi, phi_edges)
        if ti < 0 or pi < 0:
            _fill(excluded, mm_bin, weight)
            if ti >= 0:
                _fill(phi_excluded[ti], mm_bin, weight)
            continue
        _fill(accepted, mm_bin, weight)
        _fill(parents[ti], mm_bin, weight)
        _fill(children[ti][pi], mm_bin, weight)
    canonical_closure = _closure(selected, [accepted, excluded]) and _closure(accepted, parents)
    provenance_valid = _valid_provenance(provenance)
    results = []
    for ti, raw in enumerate(parents):
        integral = math.fsum(raw["values"][1:-1])
        normal_sumw2 = math.fsum(raw["variances"][1:-1])
        neff = integral * integral / normal_sumw2 if normal_sumw2 > 0 else 0.
        child_closure = _closure(raw, children[ti])
        reasons = []
        if not provenance_valid:
            reasons.append("invalid_source_provenance")
        if invalid:
            reasons.append("nonfinite_selected_records")
        if not canonical_closure or not child_closure:
            reasons.append("weighted_selection_or_child_closure_failed")
        if integral <= 0 or neff <= 1 or any(x < 0 for x in raw["values"][1:-1]):
            reasons.append("empty_or_insufficient_positive_template_support")
        valid = not reasons
        results.append({
            "status": "genuinely_local" if valid else "unavailable", "rejection_reasons": reasons,
            "t_bin_index": ti, "t_edges": t_edges[ti:ti+2], "mm_edges": mm_edges,
            "phi_edges": phi_edges, "raw": raw, "raw_children": children[ti],
            "phi_excluded": phi_excluded[ti], "effective_support": neff,
            "effective_support_definition": "normal-bin signed sum squared / normal-bin Sumw2; must exceed 1",
            "raw_integral": integral, "normal_sumw2": normal_sumw2,
            "underflow": {k: raw[k][0] for k in ("values", "variances")},
            "overflow": {k: raw[k][-1] for k in ("values", "variances")},
            "normalization_count": 1 if valid else 0,
            "normalized_values": [x / integral for x in raw["values"]] if valid else None,
            "normalized_variances": [x / integral**2 for x in raw["variances"]] if valid else None,
            "normalized_integral": math.fsum(x / integral for x in raw["values"][1:-1]) if valid else None,
            "weighted_child_to_parent_closure": child_closure,
            "canonical_selection_closure": canonical_closure,
            "provenance": deepcopy(provenance),
            "source_selection": selected,
        })
    return {"parents": results, "selected": selected, "accepted": accepted, "excluded": excluded,
            "invalid_records": invalid, "canonical_selection_closure": canonical_closure,
            "interval_convention": "[low, high), final canonical endpoint included; ROOT MM flows"}


def _clone(hist, name):
    if hist is None:
        return None
    result = hist.Clone(name + "_" + uuid4().hex[:10])
    result.SetDirectory(0)
    return result


def global_selection_closure(evidence, reference):
    """Verify the same selected weighted source as the production global shape."""
    selected = evidence.get("source_selection") or {}
    values, variances = selected.get("values") or [], selected.get("variances") or []
    if (reference.get("status") != "available" or evidence.get("mm_edges") != reference.get("edges")
            or len(values) != len(reference["values"]) or len(variances) != len(values)):
        return False
    integral = math.fsum(values[1:-1])
    if integral <= 0:
        return False
    return (all(_close(a / integral, b) for a, b in zip(values, reference["values"]))
            and all(_close(a / integral**2, b*b) for a, b in zip(variances, reference["errors"])))


def _fit_evidence(fit):
    # Retain every solver/configuration/validation scalar, including covariance,
    # errors, correlations, bound hits and staged histories. ROOT state is copied.
    return snapshot(fit)


def compare_parent(parent, templates, inp_dict, *, selected_shapes, global_shapes,
                   signal_shape, sigma0_shape, sigma0_diagnostics, fitter,
                   dynamic_alignment_enabled):
    """Fit clones only; baseline A stays owned by and referenced from parent."""
    before = digest(parent)
    baseline = parent["fit_result"]
    alignment = baseline.get("pion_component_alignment")
    result = {"t_bin_index": parent["t_bin_index"], "t_edges": list(parent["t_edges"]),
              "pion_parent_id": parent.get("pion_parent_id"),
              "canonical_interval_pair_id": parent.get("canonical_interval_pair_id"),
              "canonical_interval_pair_hash": parent.get("canonical_interval_pair_hash"),
              "coordinate_fingerprint": parent.get("coordinate_fingerprint"),
              "kaon_data_coordinate": snapshot(parent.get("kaon_data_coordinate")),
              "baseline_parent_fingerprint": before, "baseline_A": _fit_evidence(baseline),
              "candidate_B": None, "templates": {}, "rejection_reasons": [],
              "alignment_mode": "fixed_parent_alignment_template_shape_test",
              "protected_signal_templates": {"k_lambda": histogram(signal_shape), "k_sigma0": histogram(sigma0_shape)},
              "production_applied": False, "scientifically_accepted": False}
    reasons = result["rejection_reasons"]
    if not dynamic_alignment_enabled or not isinstance(alignment, dict) or not alignment:
        reasons.append("fixed_alignment_map_unavailable_no_rescan_permitted")
    elif not all(c in (alignment.get("components") or {}) for c in COMPONENTS):
        reasons.append("incomplete_fixed_alignment_map")
    inputs = {"pion": parent.get("H_pion_control"),
              "kaon": parent.get("H_proton_cleaned_final_rf")}
    result["frozen_inputs"] = {side: histogram(h) for side, h in inputs.items()}
    for side, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
        h = result["frozen_inputs"][side]
        retained = histogram(baseline.get(key))
        if (h.get("status") != "available" or retained.get("status") != "available"
                or h.get("content_error_fingerprint") != retained.get("content_error_fingerprint")):
            reasons.append(side + "_baseline_input_identity_mismatch")
        if not isinstance(h.get("stored_sumw2"), list):
            reasons.append(side + "_explicit_sumw2_unavailable")
    for component in COMPONENTS:
        payload = templates.get(component) or {}
        evidence = deepcopy(payload.get("evidence") or {})
        candidate_hist = histogram(payload.get("hist"))
        evidence["histogram"] = candidate_hist
        evidence["global_template"] = histogram(global_shapes.get(component))
        evidence["production_global_selection_closure"] = global_selection_closure(evidence, evidence["global_template"])
        evidence["baseline_selected_template"] = histogram(selected_shapes.get(component))
        evidence["baseline_selection_status"] = (
            "globally_inherited" if evidence["global_template"].get("status") == "available"
            and evidence["global_template"].get("content_error_fingerprint")
            == evidence["baseline_selected_template"].get("content_error_fingerprint")
            else "scope_resolved_or_unavailable")
        for label in ("global_template", "baseline_selected_template"):
            reference = evidence[label]
            evidence["difference_from_" + label] = (
                [a-b for a, b in zip(candidate_hist["values"], reference["values"])]
                if candidate_hist.get("status") == reference.get("status") == "available"
                and candidate_hist["edges"] == reference["edges"] else {"status": "unavailable", "reason": "axis/input unavailable"})
        result["templates"][component] = evidence
        if (evidence.get("status") != "genuinely_local" or payload.get("hist") is None
                or not evidence.get("canonical_selection_closure")
                or not evidence.get("weighted_child_to_parent_closure")
                or not evidence.get("production_global_selection_closure")
                or not _valid_provenance(evidence.get("provenance") or {})
                or evidence.get("t_edges") != list(parent["t_edges"])
                or candidate_hist.get("status") != "available"
                or candidate_hist.get("edges") != evidence.get("mm_edges")
                or candidate_hist.get("values") != evidence.get("normalized_values")
                or not isinstance(candidate_hist.get("stored_sumw2"), list)
                or len(candidate_hist.get("stored_sumw2", [])) != len(evidence.get("normalized_variances") or [])
                or not all(_close(a, b) for a, b in zip(candidate_hist.get("stored_sumw2", []), evidence.get("normalized_variances") or []))
                or evidence.get("normalization_count") != 1
                or not _close(evidence.get("normalized_integral", 0.), 1.)):
            reasons.append(component + ":invalid_local_template_or_provenance")
    candidate = None
    if not reasons:
        kwargs = dict(analysis_scope=parent["analysis_scope"],
                      kaon_signal_shape=_clone(signal_shape, "comparison_k_lambda"),
                      kaon_sigma0_shape=_clone(sigma0_shape, "comparison_k_sigma0"),
                      kaon_sigma0_source_diagnostics=deepcopy(sigma0_diagnostics),
                      mm_offset_data=baseline["template_mm_offset_data"],
                      phi_setting=baseline["particle_subtraction_phi_setting"],
                      context="template_comparison_t{}".format(parent["t_bin_index"] + 1),
                      pion_component_alignment=deepcopy(alignment))
        debug_enabled, debug_records = ownership._DEBUG_ENABLED, list(ownership._DEBUG_RECORDS)
        try:
            candidate = fitter(_clone(inputs["pion"], "comparison_pion_input"),
                               _clone(inputs["kaon"], "comparison_kaon_input"),
                               {c: _clone(templates[c]["hist"], "comparison_" + c) for c in COMPONENTS},
                               deepcopy(inp_dict), **kwargs)
            if candidate.get("pion_component_alignment") != alignment:
                raise RuntimeError("candidate_alignment_changed")
            for side, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
                if (histogram(candidate.get(key)).get("content_error_fingerprint")
                        != result["frozen_inputs"][side]["content_error_fingerprint"]):
                    raise RuntimeError("candidate_input_identity_changed:" + side)
            if candidate.get("resolved_subtraction_config") != baseline.get("resolved_subtraction_config"):
                raise RuntimeError("candidate_fitting_configuration_changed")
            result["candidate_B"] = _fit_evidence(candidate)
        except Exception as exc:
            reasons.append("candidate_fit_unavailable:" + str(exc))
            candidate = None
        finally:
            # The shared builder resets this optional process-global debug
            # ledger. Keep its diagnostic call separate from production logs.
            ownership._DEBUG_ENABLED = debug_enabled
            ownership._DEBUG_RECORDS[:] = debug_records
    result["authoritative_parent_unchanged"] = before == digest(parent)
    if not result["authoritative_parent_unchanged"]:
        raise RuntimeError("template_comparison_authoritative_parent_mutated")
    result["status"] = "available" if not reasons else "unavailable"
    return result, candidate


# Canvas-normalized geometry: legend never shares a scientific panel.
COMPARISON_LAYOUT = {
    "canvas": (1600, 1100),
    "models": ((.01, .51, .495, .945), (.505, .51, .99, .945)),
    "residuals": ((.01, .20, .495, .50), (.505, .20, .99, .50)),
    "legend": (.01, .01, .99, .18),
}


def render_comparison_pdf(path, rows, title_prefix=""):
    """Two detached pages per t: staged/refined models and data-model residuals."""
    import ROOT
    previous_pad, previous_directory = ROOT.gPad, ROOT.gDirectory
    canvas = ROOT.TCanvas("pion_t_template_comparison", "Fixed-alignment template comparison", *COMPARISON_LAYOUT["canvas"])
    def panel(name, rect):
        canvas.cd()
        pad = ROOT.TPad(name, name, *rect)
        pad.SetLeftMargin(.12)
        pad.SetRightMargin(.035)
        pad.SetBottomMargin(.14)
        pad.SetTopMargin(.10)
        pad.Draw()
        return pad
    canvas.Print(str(path) + "[")
    try:
        for evidence, baseline, candidate in rows:
            for side, data_key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
                canvas.Clear()
                model_pads = [panel("comparison_model_{}".format(i), rect)
                              for i, rect in enumerate(COMPARISON_LAYOUT["models"])]
                residual_pads = [panel("comparison_residual_{}".format(i), rect)
                                 for i, rect in enumerate(COMPARISON_LAYOUT["residuals"])]
                legend_pad = panel("comparison_legend", COMPARISON_LAYOUT["legend"])
                legend = ROOT.TLegend(.015, .04, .985 if candidate is not None else .48, .96)
                legend.SetNColumns(4 if candidate is not None else 2)
                legend.SetTextSize(.12)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                labels = set()
                def entry(hist, label, option):
                    if label not in labels:
                        legend.AddEntry(hist, label, option)
                        labels.add(label)
                keep = [*model_pads, *residual_pads, legend_pad, legend]
                for column, suffix in enumerate(("", "_refined")):
                    model_pads[column].cd()
                    data = _clone(baseline[data_key], "comparison_display_data")
                    data.SetTitle("t{} {} {}: fixed alignment, no subtraction, NOT ACCEPTED;MM (GeV);weighted content".format(
                        evidence["t_bin_index"] + 1, side, "staged" if not suffix else "refined"))
                    data.SetMarkerStyle(20)
                    data.SetStats(False)
                    models = []
                    for label, fit, color in (("A baseline", baseline, ROOT.kBlue), ("B local", candidate, ROOT.kRed)):
                        if fit is None:
                            continue
                        source = fit.get("H_{}_fit_total{}".format(side, suffix))
                        if source is not None:
                            model = _clone(source, "comparison_display_model")
                            model.SetLineColor(color)
                            model.SetLineWidth(3)
                            models.append((label, model, fit, color))
                    extrema = [data.GetBinContent(i) + data.GetBinError(i) for i in range(1, data.GetNbinsX()+1)]
                    for _, model, _, _ in models:
                        extrema += [model.GetBinContent(i) for i in range(1, model.GetNbinsX()+1)]
                    data.SetMaximum(max(extrema + [0.]) * 1.15 + 1e-12)
                    data.SetMinimum(min([data.GetBinContent(i)-data.GetBinError(i) for i in range(1, data.GetNbinsX()+1)] + [0.]))
                    data.Draw("E1")
                    entry(data, "same frozen data", "lep")
                    keep.append(data)
                    for label, model, fit, color in models:
                        model.Draw("HIST SAME")
                        entry(model, label + " total", "l")
                        keep.append(model)
                        for ci, component in enumerate((*COMPONENTS, "k_lambda", "k_sigma0")):
                            source = fit.get("H_{}_fit_{}_scaled{}".format(side, component, suffix))
                            if source is None:
                                continue
                            h = _clone(source, "comparison_display_component")
                            h.SetLineColor(color)
                            h.SetLineStyle(ci + 2)
                            h.Draw("HIST SAME")
                            entry(h, label + " " + component, "l")
                            keep.append(h)
                    residual_pads[column].cd()
                    residuals = []
                    for mi, (label, model, _, color) in enumerate(models):
                        residual = _clone(data, "comparison_display_residual")
                        residual.Add(model, -1.)
                        residual.SetTitle("data - {} total (display only);MM (GeV);residual".format("staged" if not suffix else "refined"))
                        residual.SetLineColor(color)
                        residual.SetMarkerColor(color)
                        residuals.append(residual)
                        keep.append(residual)
                    if residuals:
                        bounds = [h.GetBinContent(i) + sign*h.GetBinError(i) for h in residuals
                                  for i in range(1, h.GetNbinsX()+1) for sign in (-1., 1.)]
                        residuals[0].SetMinimum(min(bounds + [0.]))
                        residuals[0].SetMaximum(max(bounds + [0.]) + 1e-12)
                        for ri, residual in enumerate(residuals):
                            residual.Draw("E1" if ri == 0 else "E1 SAME")
                legend_pad.cd()
                legend.Draw()
                if candidate is None:
                    note = ROOT.TPaveText(.50, .04, .99, .96, "NDC")
                    note.SetTextSize(.09)
                    note.AddText("B unavailable; no global fallback")
                    for reason in evidence["rejection_reasons"]:
                        note.AddText(reason)
                    note.Draw()
                    keep.append(note)
                canvas.cd()
                header = ROOT.TPaveText(.05, .965, .95, .999, "NDC")
                header.SetTextSize(.025)
                header.AddText(title_prefix)
                header.Draw()
                keep.append(header)
                canvas.Print(str(path))
    finally:
        canvas.Print(str(path) + "]")
        canvas.Close()
        if previous_directory:
            previous_directory.cd()
        if previous_pad:
            previous_pad.cd()


def observe_template_comparison(hist, inp_dict, *, t_edges, phi_edges, component_payload,
                                signal_payload, sigma0_payload, outpath,
                                hgcer_cutg=None, original_inputs=None):
    """One environment opt-in; no import/load/fit/write on the default path.

    Ordinary launcher passes its environment to main.py (including plain -d).
    Limit capture to the contracted Q4p4W2p74 Center/Left low-epsilon settings.
    Each setting gets its own JSON/PDF; populations are never pooled.
    """
    if not enabled():
        return None
    if (str(inp_dict.get("Q2")) != "4p4" or str(inp_dict.get("W")) != "2p74"
            or str(inp_dict.get("EPSSET")).lower() not in ("low", "lowe")
            or hist.get("phi_setting") not in ("Center", "Left")):
        return None
    from background_config import get_pion_component_dynamic_alignment_config
    from pion_component_shapes import load_raw_t_pion_comparison_shapes
    from pion_component_fits import (build_particle_subtraction_component_result,
                                    resolve_scope_component_shapes, resolve_scope_single_shape)
    parents = hist.get("_pion_t_bin_parent_results") or ()
    if not parents:
        raise RuntimeError("template_comparison_requires_frozen_t_parents")
    if len(parents) != 3 or len(t_edges) != 4:
        raise RuntimeError("template_comparison_requires_contracted_three_canonical_t_intervals")
    raw = load_raw_t_pion_comparison_shapes(inp_dict, hist["phi_setting"], t_edges, phi_edges,
                                          hgcer_cutg=hgcer_cutg)
    rows, reports = [], []
    for parent in parents:
        ti = parent["t_bin_index"]
        scope = parent["analysis_scope"]
        report, candidate = compare_parent(
            parent, {c: raw[c]["parents"][ti] for c in COMPONENTS}, inp_dict,
            selected_shapes=resolve_scope_component_shapes(component_payload, analysis_scope=scope, t_bin_index=ti),
            global_shapes=resolve_scope_component_shapes(component_payload, analysis_scope="setting-wide"),
            signal_shape=resolve_scope_single_shape(signal_payload, analysis_scope=scope, t_bin_index=ti),
            sigma0_shape=resolve_scope_single_shape(sigma0_payload, analysis_scope=scope, t_bin_index=ti),
            sigma0_diagnostics=(sigma0_payload or {}).get("diagnostics") or {},
            fitter=build_particle_subtraction_component_result,
            dynamic_alignment_enabled=bool(get_pion_component_dynamic_alignment_config(
                inp_dict=inp_dict, phi_setting=hist["phi_setting"]).get("enabled")))
        reports.append(report)
        rows.append((report, parent["fit_result"], candidate))
    root = Path(__file__).resolve().parents[2]
    source = {p: file_identity(root / p) for p in REQUIRED_SOURCE_FILES}
    source_rejections = ["required_source_identity_unavailable:" + p for p, identity in source.items()
                         if identity.get("status") != "available" or not identity.get("sha256")]
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    result = {"schema_version": SCHEMA, "diagnostic_only": True, "production_applied": False,
              "scientifically_accepted": False, "source_HEAD": head, "source_files": source,
              "source_identity_rejection_reasons": source_rejections,
              "setting": {k: inp_dict.get(k) for k in ("Q2", "W", "EPSSET", "ParticleType")},
              "phi_setting": hist["phi_setting"], "enable_environment_variable": ENABLE_ENV,
              "original_analysis_input_identity": {role: file_identity(path) for role, path in (original_inputs or {}).items()},
              "original_input_population_accounting": snapshot({
                  "source_audit": hist.get("pion_control_source_audit"),
                  "by_t": [{key: p.get(key) for key in ("source_accounting", "pion_control_source_accounting", "input_selection", "source_target_state")}
                           for p in parents]}),
              "canonical_t_edges": list(map(float, t_edges)), "canonical_phi_edges": list(map(float, phi_edges)),
              "raw_selection": {c: {k: v for k, v in raw[c].items() if k != "parents"} for c in COMPONENTS},
              "parents": reports, "status": ("unavailable" if source_rejections else
                  "available" if all(r["status"] == "available" for r in reports) else "partial_unavailable")}
    output = Path(outpath)
    output.mkdir(parents=True, exist_ok=True)
    stem = "{}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison".format(hist["phi_setting"])
    json_path, pdf_path = output / (stem + ".json"), output / (stem + ".pdf")
    result["pdf"] = {"path": str(pdf_path), "status": "not_rendered"}
    json_path.write_text(json.dumps(snapshot(result), sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    try:
        render_comparison_pdf(pdf_path, rows, "{} Q4p4W2p74 lowe | source {} | diagnostic only".format(hist["phi_setting"], head[:12]))
    except Exception as exc:
        result["status"] = "failed"
        result["pdf"] = {"path": str(pdf_path), "status": "unavailable", "reason": str(exc)}
        json_path.write_text(json.dumps(snapshot(result), sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        raise
    result["pdf"] = file_identity(pdf_path)
    json_path.write_text(json.dumps(snapshot(result), sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result
