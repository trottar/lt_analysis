"""Observations of retained pion parents, never a fit or production decision.

All arithmetic operates on scalar copies. ROOT is imported only by the optional
companion renderer; neither capture nor rendering touches a source histogram.
The JSON is evidence for a later checker, not an acceptance receipt.
"""
from __future__ import annotations

import bisect
import hashlib
import json
import math
import numbers
import os
from pathlib import Path
import tempfile
import textwrap

from root_histogram_ownership import fingerprint_histogram_content_error

SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
COMPONENTS = ("pi_n", "pi_sidis", "pi_delta")
SCHEMA = "pion_background_chain_audit/v1"
REL_TOL, ABS_TOL = 1e-10, 1e-12  # arithmetic rounding only, no physics gate


def unavailable(reason):
    return {"status": "unavailable", "reason": str(reason)}


def finite(value):
    return isinstance(value, numbers.Real) and math.isfinite(value)


def snapshot(value):
    """Copy retained evidence, making missing/nonfinite positions explicit."""
    if value is None:
        return unavailable("not retained by producer")
    if hasattr(value, "GetNbinsX"):
        return histogram(value)
    if isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, numbers.Real):
        return float(value) if finite(value) else unavailable("nonfinite: " + str(value))
    if isinstance(value, dict):
        return {str(k): snapshot(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [snapshot(v) for v in value]
    if hasattr(value, "tolist"):
        return snapshot(value.tolist())
    return unavailable("unsupported retained type: " + type(value).__name__)


def encoded(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")) + "\n").encode()


def digest(value):
    return hashlib.sha256(encoded(snapshot(value))).hexdigest()


def histogram(hist):
    if hist is None:
        return unavailable("histogram not retained")
    try:
        n = int(hist.GetNbinsX())
        edges = [float(hist.GetXaxis().GetBinLowEdge(i)) for i in range(1, n + 2)]
        values = [float(hist.GetBinContent(i)) for i in range(n + 2)]
        errors = [float(hist.GetBinError(i)) for i in range(n + 2)]
        valid = (n > 0 and all(map(finite, edges + values + errors))
                 and all(a < b for a, b in zip(edges, edges[1:])))
        result = {"status": "available" if valid else "unavailable",
                  "edges": snapshot(edges), "values": snapshot(values), "errors": snapshot(errors),
                  "bin_layout": "underflow, normal bins 1..N, overflow; no width weighting",
                  "axis_signature": digest(edges),
                  "content_error_fingerprint": (fingerprint_histogram_content_error(hist)
                                                if valid else unavailable("nonfinite/invalid axis"))}
        if not valid:
            result["reason"] = "nonfinite histogram or invalid axis"
        try:
            result["stored_sumw2"] = ([float(hist.GetSumw2()[i]) for i in range(n + 2)]
                                      if hist.GetSumw2N() else unavailable("no explicit Sumw2 storage"))
        except (AttributeError, TypeError):
            result["stored_sumw2"] = unavailable("Sumw2 accessor not available")
        return snapshot(result)
    except Exception as exc:
        return unavailable("histogram read failed: " + str(exc))


def available(h):
    return isinstance(h, dict) and h.get("status") == "available"


def compare(left, right, indices=None):
    if left is None or right is None or len(left) != len(right):
        return unavailable("missing arrays or different bin counts")
    indices = list(range(len(left))) if indices is None else list(indices)
    if not indices:
        return unavailable("no qualifying bins")
    if not all(finite(left[i]) and finite(right[i]) for i in indices):
        return unavailable("nonfinite/missing comparison bin")
    delta = [left[i] - right[i] for i in indices]
    return {"status": "PASS" if all(math.isclose(left[i], right[i], rel_tol=REL_TOL,
                                                  abs_tol=ABS_TOL) for i in indices) else "FAIL",
            "bin_indices": indices, "differences": delta,
            "signed_difference": math.fsum(delta),
            "max_absolute_difference": max(map(abs, delta))}


def hist_compare(left, right):
    if not available(left) or not available(right):
        return unavailable("comparison histogram unavailable")
    if left["edges"] != right["edges"]:
        return unavailable("axis signature mismatch")
    return compare(left["values"], right["values"])


def source_spectra(records, reference, parent, mask="nommcuts", weights=None, final_t_endpoint=None):
    """Mirror cached scalar selection and ROOT edge/flow lookup, without Fill."""
    if records is None or not available(reference):
        return unavailable("records or original control axis unavailable")
    edges = reference["edges"]
    n = len(edges) - 1
    cells = {s: {k: [0.0] * (n + 2) for k in ("count", "signed", "absolute", "sumsq")}
             for s in SOURCES}
    invalid = []
    for ordinal, record in enumerate(records):
        try:
            if record["t_index"] != parent["t_bin_index"]:
                raise ValueError("t_index mismatch")
            if record["coordinate_fingerprint"] != parent.get("coordinate_fingerprint"):
                raise ValueError("coordinate fingerprint mismatch")
            x, coefficient, t = (float(record[k]) for k in ("adj_MM", "coefficient", "adj_t"))
            if not all(map(finite, (x, coefficient, t))):
                raise ValueError("nonfinite source coordinate/coefficient")
            in_interval = parent["t_edges"][0] <= t < parent["t_edges"][1]
            final_endpoint = final_t_endpoint is not None and t == final_t_endpoint == parent["t_edges"][1]
            if not (in_interval or final_endpoint):
                raise ValueError("record outside half-open canonical t interval")
            source = record["source_label"]
            if source not in SOURCES:
                raise ValueError("unknown source")
            selected = record[mask]
            if not selected:
                continue
            i = bisect.bisect_right(edges, x)
            weight = coefficient
            if weights is not None:
                # The producer lookup forces underflow to zero; overflow uses
                # the stored overflow entry (zero for the component proposal).
                factor = 0.0 if i == 0 else float(weights[i])
                if not finite(factor):
                    raise ValueError("nonfinite stored weight")
                weight *= factor
            cell = cells[source]
            cell["count"][i] += 1
            cell["signed"][i] += weight
            cell["absolute"][i] += abs(weight)
            cell["sumsq"][i] += weight * weight
        except (KeyError, TypeError, ValueError, IndexError) as exc:
            invalid.append({"record": ordinal, "reason": str(exc)})
    total = {k: [math.fsum(cells[s][k][i] for s in SOURCES) for i in range(n + 2)]
             for k in ("count", "signed", "absolute", "sumsq")}
    result = {"status": "available" if not invalid else "unavailable",
              "selection": mask, "edges": edges, "sources": cells, "combined": total,
              "invalid_records": invalid, "records_fingerprint": digest(records),
              "flow_semantics": "ROOT half-open MM bins; application underflow forced zero, overflow stored weight",
              "canonical_t_final_endpoint": snapshot(final_t_endpoint)}
    if invalid:
        result["reason"] = "source identity/field failure; partial sums are inadmissible"
    return result


def source_closure(sources, original):
    if not available(sources) or not available(original) or sources["edges"] != original["edges"]:
        return unavailable("source evidence unavailable or axis mismatch")
    return {"content": compare(sources["combined"]["signed"], original["values"]),
            "error_squared": compare(sources["combined"]["sumsq"], [e * e for e in original["errors"]]),
            "stored_sumw2": compare(sources["combined"]["sumsq"], original["stored_sumw2"])
            if isinstance(original.get("stored_sumw2"), list) else unavailable("explicit Sumw2 not retained")}


def fit_observation(fit):
    keys = ("resolved_subtraction_config", "diagnostics", "staged_amplitudes_raw",
            "staged_amplitudes_scaled", "refined_amplitudes", "pion_component_alignment",
            "fit_status_pion", "fit_status_kaon", "fit_validation_pion", "fit_validation_kaon",
            "fallback_used", "fallback_reason", "template_mm_offset_data")
    result = {k: snapshot(fit.get(k)) for k in keys}
    result["amplitudes"] = {k: snapshot(fit.get(k)) for k in
                            ("A_n", "A_sidis", "A_delta", "B_n", "B_sidis", "B_delta")}
    hist_keys = ["H_simc_shape_" + c for c in COMPONENTS]
    hist_keys += ["H_pion_control_input", "H_kaon_nosub_input", "H_simc_shape_k_lambda",
                  "H_simc_shape_k_sigma0", "H_kaon_pion_bg_fit_total", "H_fit_residual_pion",
                  "H_fit_residual_pion_refined", "H_kaon_full_fit_residual_refined"]
    for side in ("pion", "kaon"):
        for suffix in ("", "_refined"):
            hist_keys += ["H_{}_fit_{}_scaled{}".format(side, c, suffix) for c in COMPONENTS]
            hist_keys += ["H_{}_fit_total{}".format(side, suffix)]
    hist_keys += ["H_kaon_fit_{}_scaled{}".format(signal, suffix)
                  for signal in ("k_lambda", "k_sigma0") for suffix in ("", "_refined")]
    result["histograms"] = {k: histogram(fit.get(k)) for k in hist_keys}
    result["steps"] = {}
    for side in ("pion", "kaon"):
        steps = fit.get("H_{}_fit_step_overlays".format(side)) or []
        observed = []
        for step in steps:
            copy = snapshot(step)
            a, b, c = (histogram(step.get(k)) for k in
                       ("H_baseline_before", "H_component_scaled", "H_cumulative_after"))
            copy["cumulative_arithmetic"] = (
                compare([x + y for x, y in zip(a["values"], b["values"])], c["values"])
                if all(available(h) for h in (a, b, c)) and a["edges"] == b["edges"] == c["edges"]
                else unavailable("stored step histograms missing/mismatched"))
            residual = histogram(step.get("H_residual_input"))
            copy["output_residual_arithmetic"] = (
                {"status": "derived_arithmetic_only", "edges": residual["edges"],
                 "values": [x - y for x, y in zip(residual["values"], b["values"])]}
                if available(residual) and available(b) and residual["edges"] == b["edges"]
                else unavailable("stored residual/contribution missing/mismatched"))
            # Success/message belong to the retained pass history, not the overlay.
            history = ((fit.get("diagnostics") or {}).get(side) or {}).get("staged_pass_history") or []
            index = int(step.get("pass_index", 0)) - 1
            copy["solver_outcome"] = snapshot(history[index].get(step.get("component_name"))
                                               if 0 <= index < len(history) else None)
            observed.append(copy)
        result["steps"][side] = observed or unavailable("no stored staged overlays")
    result["signal_roles"] = "K-Lambda/K-Sigma0 are signal fit templates, excluded from pion subtraction sums"
    return result


def selected_amplitudes(fit):
    diag = fit.get("diagnostics") or {}
    chosen, provenance = {}, {}
    for side in ("pion", "kaon"):
        d = diag.get(side) or {}
        key = "refined_amplitudes" if d.get("refined_amplitudes") else "staged_amplitudes_scaled"
        source = d.get(key) or {}
        protected = d.get("pi_delta_signal_protected_fit") or {}
        if side == "kaon" and protected.get("enabled") and protected.get("protected_applied_amplitudes"):
            source = protected["protected_applied_amplitudes"]
            key = "pi_delta_signal_protected_fit.protected_applied_amplitudes"
        prefix = "B" if side == "pion" else "A"
        chosen[side] = {c: source.get(c, fit.get(prefix + "_" + c[3:])) for c in COMPONENTS}
        provenance[side] = key + " with scalar fallback"
    return chosen, provenance


def transfer_observation(payload, fit):
    hs = {k: histogram(payload.get(k)) for k in
          ("H_pion_control_model", "H_kaon_pion_model", "H_pion_weight_vs_MM")}
    p, k, w = (hs[key] for key in hs)
    result = {"histograms": hs, "weights": snapshot(payload.get("weights")),
              "producer_diagnostics": snapshot(payload.get("diagnostics")),
              "interpretation": "model algebra only; never signed-event/data agreement"}
    amps, origin = selected_amplitudes(fit)
    result.update(selected_amplitudes=snapshot(amps), amplitude_source=origin, model_variant="final proposal")
    shapes = [histogram(fit.get("H_simc_shape_" + c)) for c in COMPONENTS]
    if not all(available(h) for h in [p, k, w] + shapes):
        result["closure"] = unavailable("model/shape/weight unavailable")
        return result
    if not all(h["edges"] == p["edges"] for h in [k, w] + shapes):
        result["closure"] = unavailable("axis signature mismatch")
        return result
    if not all(finite(amps[side][c]) for side in amps for c in COMPONENTS):
        result["closure"] = unavailable("selected amplitudes not retained/finite")
        return result
    reconstructed = {side: [math.fsum(amps[side][c] * h["values"][i]
                                      for c, h in zip(COMPONENTS, shapes)) for i in range(len(p["values"]))]
                     for side in amps}
    diag = payload.get("diagnostics") or {}
    floor, lo, hi = (diag.get(key) for key in ("denom_floor", "clip_min", "clip_max"))
    bins = []
    supported = []
    for i, (den, num, weight) in enumerate(zip(p["values"], k["values"], w["values"])):
        normal = 0 < i < len(p["values"]) - 1
        known = finite(floor)
        raw = (num / den if den > floor else 0.0) if known else None
        if raw is not None and not finite(raw):
            raw = 0.0  # retained producer's finite guard, before clipping
        expected = raw
        if expected is not None and finite(lo):
            expected = max(lo, expected)
        if expected is not None and finite(hi):
            expected = min(hi, expected)
        unclipped = raw is not None and math.isclose(raw, weight, abs_tol=1e-15, rel_tol=0)
        qualifies = normal and known and den > floor and unclipped
        if qualifies:
            supported.append(i)
        bins.append({"bin": i, "flow": not normal, "denominator_supported": den > floor if known else None,
                     "unsupported_numerator": den <= floor < num if known else None,
                     "unclipped": unclipped if known else unavailable("floor absent"),
                     "configured_clipping": abs(expected - raw) > 1e-15 if known else None,
                     "stored_weight_policy_match": math.isclose(expected, weight, rel_tol=REL_TOL, abs_tol=ABS_TOL) if known else None,
                     "identity_qualified": qualifies, "P_model": den, "K_model": num, "w0": weight,
                     "raw_ratio_or_zero": snapshot(raw)})
    result.update(bins=bins, denominator_floor=snapshot(floor), clip_min=lo, clip_max=hi,
                  qualified_bins=supported,
                  closure={"P_components": compare(reconstructed["pion"], p["values"]),
                           "K_components": compare(reconstructed["kaon"], k["values"]),
                           "model_algebra": compare([a * b for a, b in zip(w["values"], p["values"])],
                                                    k["values"], supported),
                           "weight_array_vs_hist": compare(result["weights"], w["values"])
                           if isinstance(result["weights"], list) else unavailable("weights absent")})
    return result


def application_observation(payload, parent, records, reference, transfer=None, final_t_endpoint=None):
    if not payload:
        return unavailable("existing application payload absent")
    h = {k: histogram(payload.get(k)) for k in
         ("H_pion_subtraction_template_MM_nosub", "H_MM_nosub_before_pion_subtraction",
          "H_MM_nosub_after_pion_subtraction", "H_pion_subtraction_template_MM",
          "H_MM_before_pion_subtraction", "H_MM_after_pion_subtraction")}
    result = {"histograms": h, "policy": snapshot({k: v for k, v in payload.items()
               if not k.startswith("H_") and k not in ("weights", "diagnostics")}),
              "producer_diagnostics": snapshot(payload.get("diagnostics")),
              "weights": snapshot(payload.get("weights"))}
    weight_reference = histogram(payload.get("H_pion_control_model"))
    result["weight_reference"] = weight_reference
    applied = h["H_pion_subtraction_template_MM_nosub"]
    for mask, suffix in (("nommcuts", "_nosub"), ("allcuts", "")):
        weights = payload.get("weights")
        compatible = (available(weight_reference) and available(reference)
                      and weight_reference["edges"] == reference["edges"] and weights is not None
                      and len(weights) == len(reference["values"]))
        sources = (source_spectra(records, reference, parent, mask, weights, final_t_endpoint)
                   if compatible else unavailable("stored application weight reference/array absent or axis mismatch"))
        result[mask + "_weighted_sources"] = sources
        actual = h["H_pion_subtraction_template_MM" + suffix]
        result[mask + "_event_reproduction"] = source_closure(sources, actual)
        before, after = (h["H_MM" + suffix + "_" + stage + "_pion_subtraction"] for stage in ("before", "after"))
        result[mask + "_after_closure"] = (
            compare([b - a for b, a in zip(before["values"], actual["values"])], after["values"])
            if all(available(x) for x in (before, actual, after)) and before["edges"] == actual["edges"] == after["edges"]
            else unavailable("before/applied/after absent or axes mismatched"))
    if transfer and available(applied):
        p = transfer["histograms"]["H_pion_control_model"]
        k = transfer["histograms"]["H_kaon_pion_model"]
        w = transfer["histograms"]["H_pion_weight_vs_MM"]
        if all(available(x) for x in (p, k, w, reference)) and all(x["edges"] == applied["edges"] for x in (p, k, w, reference)):
            residual = [x - y for x, y in zip(reference["values"], p["values"])]
            difference = [x - y for x, y in zip(applied["values"], k["values"])]
            predicted = [x * y for x, y in zip(w["values"], residual)]
            qualified = transfer.get("qualified_bins", [])
            result.update(control_minus_model=residual, applied_minus_model=difference,
                          weighted_control_residual=predicted,
                          event_model_closure=compare(applied["values"], k["values"], qualified),
                          residual_identity=compare(difference, predicted, qualified))
            normal = list(range(1, len(applied["values"]) - 1))
            top = sorted(normal, key=lambda i: (-abs(difference[i]), i))[:8]
            edges = applied["edges"]
            result["top_discrepant_bins"] = [{"bin": i, "mm_edges": edges[i-1:i+1],
                                               "difference": difference[i], "w0": w["values"][i],
                                               "identity_qualified": i in qualified} for i in top]
            result["historical_1p03_bin"] = (bisect.bisect_right(edges, 1.03)
                                               if edges[0] <= 1.03 < edges[-1] else unavailable("1.03 outside axis"))
            result["normal_bin_integrals"] = {"applied": math.fsum(applied["values"][i] for i in normal),
                "K_model": math.fsum(k["values"][i] for i in normal),
                "absolute_bin_difference": math.fsum(abs(difference[i]) for i in normal)}
    return result


def is_target(inp, phi):
    return (inp.get("ParticleType") == "kaon" and inp.get("Q2") == "4p4"
            and inp.get("W") == "2p74" and inp.get("EPSSET") == "low" and phi == "Center"
            and inp.get("particle_subtraction_mode", "simc_shape_components") == "simc_shape_components"
            and inp.get("pion_subtraction_scope", "t_bin") == "t_bin")


def build_chain_audit(hist, inp, parents, shape_payload):
    if not is_target(inp, hist.get("phi_setting")):
        return None
    cache = hist.get("_authoritative_pion_control_source_cache") or {}
    # The cache uses find_canonical_bin: only the final overall endpoint is
    # closed. Determine that endpoint from the complete supplied collection.
    final_entry = max(cache.get("by_t", []), key=lambda p: p.get("t_index", -1), default={})
    final_parent = next((p for p in parents if p.get("t_bin_index") == final_entry.get("t_index")), {})
    final_t_endpoint = (final_entry.get("t_edges") or final_parent.get("t_edges") or [None])[-1]
    results = []
    for index in range(3):
        matches = [p for p in parents if p.get("t_bin_index") == index]
        entries = [e for e in cache.get("by_t", []) if e.get("t_index") == index]
        if len(matches) != 1 or len(entries) != 1:
            results.append({"t_bin_index": index, **unavailable("expected exactly one canonical parent/cache entry")})
            continue
        parent, entry = matches[0], entries[0]
        fit = parent.get("fit_result") or {}
        control = histogram(parent.get("H_pion_control"))
        records = entry.get("records")
        source = source_spectra(records, control, parent, final_t_endpoint=final_t_endpoint)
        cut = histogram(entry.get("H_pion_control_cut"))
        cut_source = source_spectra(records, cut, parent, "allcuts", final_t_endpoint=final_t_endpoint)
        identity = {k: snapshot(parent.get(k)) for k in
                    ("pion_parent_id", "t_bin_index", "t_edges", "analysis_scope", "coordinate_fingerprint",
                     "pion_control_input_fingerprint", "proton_producer_fingerprint", "input_selection", "source_target_state")}
        coordinate = parent.get("coordinate_fingerprint")
        identity_ok = (bool(coordinate) and coordinate == entry.get("coordinate_fingerprint") == cache.get("coordinate_fingerprint")
                       and available(control) and control["content_error_fingerprint"] == parent.get("pion_control_input_fingerprint"))
        proposal = parent.get("proposed_diagnostic_application_result") or {}
        transfer = transfer_observation(proposal, fit)
        simc = {}
        for c in COMPONENTS:
            component = ((shape_payload or {}).get("components") or {}).get(c) or {}
            shapes = (component.get("binned_shapes") or {}).get("t_bin{}".format(index + 1)) or {}
            simc[c] = {"diagnostics": snapshot(component.get("diagnostics")),
                       "phi_inputs": {str(k): histogram(v) for k, v in shapes.items()},
                       "raw_phi_integrals": unavailable("producer retains normalized phi shapes only; not reconstructed")}
        observed = {"status": "available" if identity_ok else "unavailable", "identity": identity,
                    "t_bin_index": index, "identity_closure": "PASS" if identity_ok else "FAIL",
                    "control": control, "control_cut": cut, "source_spectra": source,
                    "source_closure": source_closure(source, control),
                    "cut_source_closure": source_closure(cut_source, cut),
                    "source_accounting": snapshot(entry.get("source_accounting")),
                    "cache_control_closure": hist_compare(control, histogram(entry.get("H_pion_control"))),
                    "fit_input_closure": hist_compare(control, histogram(fit.get("H_pion_control_input"))),
                    "simc_provenance": simc, "fit": fit_observation(fit), "transfer": transfer,
                    "application_policy": snapshot(parent.get("production_application_policy")),
                    "application_status": snapshot(parent.get("diagnostic_application_status")),
                    "proposal": application_observation(proposal, parent, records, control, transfer, final_t_endpoint),
                    "final": application_observation(parent.get("final_diagnostic_application_result"), parent, records, control,
                                                     final_t_endpoint=final_t_endpoint)}
        if not identity_ok:
            observed["reason"] = "parent/cache coordinate or original histogram fingerprint mismatch"
        results.append(observed)
    audit = snapshot({"schema_version": SCHEMA, "diagnostic_only": True, "production_objects_mutated": False,
             "scope": "Q4p4W2p74/Center/lowe/canonical-t1-t2-t3", "primary_parent": "t2",
             "observed_source_HEAD": "not_provided", "parents": results,
             "cache_identity": {k: cache.get(k) for k in ("definition", "coordinate_fingerprint", "kaon_data_coordinate",
                 "physical_pion_control_mask", "physical_pion_control_mask_fingerprint", "pid_contract", "source_accounting")},
             "setting_wide_comparator": unavailable("different-scope fit not passed at parent boundary; never substituted"),
             "limits": ["Current same-run detached observations, not a historical first-bad SHA or farm acceptance.",
                        "Arithmetic tolerances are rounding checks, not scientific model/data acceptance.",
                        "Final policy/fallback is distinct from component proposal and downstream production.",
                        "Unavailable fields or failed closures require later independent checker review."]})
    audit["collection_identity"] = digest([p.get("identity", p) for p in results])
    # A retained diagnostic is never implicitly complete merely because JSON
    # can be written. Expected missing stages remain machine-readable.
    missing = []
    for p in results:
        label = "t{}".format(p["t_bin_index"] + 1)
        if "identity" not in p:
            missing.append(label + ": parent/cache unavailable")
            continue
        for key in ("control", "control_cut", "source_spectra", "proposal", "final"):
            if p[key].get("status") == "unavailable":
                missing.append(label + ": " + key)
        if p["identity_closure"] != "PASS":
            missing.append(label + ": identity mismatch")
        for side in ("pion", "kaon"):
            if not isinstance(p["fit"]["steps"][side], list):
                missing.append(label + ": " + side + " steps")
        if p["transfer"].get("closure", {}).get("status") == "unavailable":
            missing.append(label + ": transfer evidence")
    audit["capture_completeness"] = {"status": "incomplete" if missing else "captured",
                                     "missing_required_evidence": missing,
                                     "acceptance": "NOT VERIFIED; independent checker and page review required"}
    return audit


def _series(h):
    return h["values"][1:-1] if available(h) else None


def page_specs(audit):
    """Pure page plan, also usable by deterministic renderer tests."""
    pages = [{"title": "Baseline pion full-chain diagnostic (non-authoritative)", "panels": [],
              "text": [audit["scope"], "Separate t1 / t2 / t3; t2 primary", "Model closure is not event closure",
                       "capture: " + audit.get("capture_completeness", {}).get("status", "unavailable")]}]
    for parent in audit["parents"]:
        index = parent["t_bin_index"] + 1
        proposal = parent.get("proposal", {})
        model = parent.get("transfer", {}).get("closure", {}).get("model_algebra", {}).get("status", "unavailable")
        event = proposal.get("event_model_closure", {}).get("status", "unavailable")
        after = parent.get("final", {}).get("nommcuts_after_closure", {}).get("status", "unavailable")
        pages[0]["text"].append("t{}: model {} / event {} / after {}".format(index, model, event, after))
    for parent in audit["parents"]:
        label = "Center lowe t{}".format(parent["t_bin_index"] + 1)
        if "identity" not in parent:
            for stage in ("Control construction", "SIMC and staged/joint fits", "Transfer", "Application and outcome"):
                pages.append({"title": label + " - " + stage, "panels": [], "text": ["unavailable: " + parent["reason"]]})
            continue
        edges = parent["control"].get("edges", [])
        def panel(title, series, unit="signed yield"):
            return {"title": title, "edges": edges, "series": [(name, values) for name, values in series if values is not None], "unit": unit}
        sources = parent["source_spectra"].get("sources", {})
        control = [(s, sources[s]["signed"][1:-1]) for s in SOURCES if s in sources]
        control += [("original control", _series(parent["control"]))]
        def matching(h):
            return h if available(h) and h["edges"] == edges else unavailable("missing histogram or axis mismatch")
        fit = {k: matching(v) for k, v in parent["fit"]["histograms"].items()}
        transfer = {k: matching(v) for k, v in parent["transfer"]["histograms"].items()}
        proposal, final = parent["proposal"], parent["final"]
        weighted = proposal.get("nommcuts_weighted_sources", {}).get("sources", {})
        pages.append({"title": label + " - Control construction", "panels": [panel("Four signed sources and original control", control),
            panel("Uncertainty and absolute support", [("original error", parent["control"].get("errors", [])[1:-1]),
                ("absolute event support", parent["source_spectra"].get("combined", {}).get("absolute", [])[1:-1])])],
            "text": ["source content / error closure: " + str([parent["source_closure"].get(k, {}).get("status", "unavailable") for k in ("content", "error_squared")]),
                     "identity: " + parent["identity_closure"],
                     "source counts/sums/sumsq and flows retained in JSON"]})
        fit_panels = [panel("Aligned SIMC shapes", [(c, _series(fit["H_simc_shape_" + c])) for c in COMPONENTS], "unit shape"),
            panel("Staged vs joint control", [("control", _series(parent["control"]))] +
                  [(c + suffix, _series(fit["H_pion_fit_" + c + "_scaled" + suffix])) for suffix in ("", "_refined") for c in COMPONENTS]),
            panel("Protected kaon stage (signal separate)", [(c, _series(fit["H_kaon_fit_" + c + "_scaled"])) for c in COMPONENTS])]
        # Actual retained progression gets its own pads/pages; no fit is reconstructed.
        for side in ("pion", "kaon"):
            steps = parent["fit"]["steps"][side]
            if isinstance(steps, list):
                for step in steps:
                    fit_panels.append(panel(side + " step " + str(step.get("step_index")),
                        [(k, _series(matching(step.get(k, {})))) for k in ("H_residual_input", "H_component_scaled", "H_cumulative_after")]))
        for offset in range(0, len(fit_panels), 3):
            pages.append({"title": label + " - SIMC and staged/joint fits " + str(offset // 3 + 1),
                          "panels": fit_panels[offset:offset+3], "text": ["staged and joint success are distinct; see retained diagnostics",
                          "stored overlays: " + ", ".join(side + " " + ("available" if isinstance(parent["fit"]["steps"][side], list) else "unavailable") for side in ("pion", "kaon")),
                          "pion status: " + str(parent["fit"]["fit_status_pion"]),
                          "joint success: " + str(((parent["fit"]["diagnostics"].get("pion") or {}).get("joint_refinement") or {}).get("success", "unavailable")),
                          "chi2/ndf, p: " + str([((parent["fit"]["diagnostics"].get("pion") or {}).get(k, "unavailable")) for k in ("chi2_ndf", "fit_p_value")]),
                          "protected status: " + str(((parent["fit"]["diagnostics"].get("kaon") or {}).get("pi_delta_signal_protected_fit") or {}).get("status", "unavailable"))]})
        top = proposal.get("top_discrepant_bins", [])
        marker = (sum(top[0]["mm_edges"]) / 2) if top else None
        pages.append({"title": label + " - Transfer", "panels": [
            panel("Signed control and P_model", [("control", _series(parent["control"])), ("P_model", _series(transfer["H_pion_control_model"]))]),
            panel("K_model (pion components only)", [("K_model", _series(transfer["H_kaon_pion_model"]))]),
            panel("Stored w0 (separate pad)", [("w0", _series(transfer["H_pion_weight_vs_MM"]))], "weight"),
            panel("Control minus model", [("control residual", proposal.get("control_minus_model", [])[1:-1])])],
            "marker": marker, "text": ["Model algebra only; support/clip qualifications in JSON", "floor: " + str(parent["transfer"].get("denominator_floor", "unavailable"))]})
        ph = {k: matching(v) for k, v in proposal.get("histograms", {}).items()}
        fh = {k: matching(v) for k, v in final.get("histograms", {}).items()}
        weighted_source = proposal.get("nommcuts_weighted_sources", {})
        support_text = ([math.fsum(weighted_source["combined"][k]) for k in ("signed", "absolute")]
                        if available(weighted_source) else "unavailable")
        pages.append({"title": label + " - Application and outcome", "panels": [
            panel("K_model vs signed event application", [("K_model", _series(transfer["H_kaon_pion_model"])),
                ("proposal applied", _series(ph.get("H_pion_subtraction_template_MM_nosub", {})))]),
            panel("Four weighted sources", [(s, weighted[s]["signed"][1:-1]) for s in SOURCES if s in weighted]),
            panel("Application minus model", [("event-model difference", proposal.get("applied_minus_model", [])[1:-1]),
                ("weighted control residual", proposal.get("weighted_control_residual", [])[1:-1])]),
            panel("Final policy before/after (may be fallback)", [(stage, _series(fh.get("H_MM_nosub_" + stage + "_pion_subtraction", {}))) for stage in ("before", "after")])],
            "marker": marker, "text": ["event-model closure: " + str(proposal.get("event_model_closure", {}).get("status", "unavailable")),
                "final after closure: " + str(final.get("nommcuts_after_closure", {}).get("status", "unavailable")),
                "weighted signed / absolute: " + str(support_text),
                "policy: " + str(parent["application_status"].get("final_status", "unavailable"))]})
    return pages


def render_companion(path, audit, root=None):
    if root is None:
        import ROOT as root
    from array import array
    pages = page_specs(audit)
    canvas = root.TCanvas("pion_chain_" + audit["collection_identity"][:12], "Pion chain diagnostic", 1600, 1000)
    try:
        canvas.Print(str(path) + "[")
        for page_index, spec in enumerate(pages):
            canvas.Clear()
            canvas.Divide(3, 2)
            retained = []
            for pad_index, panel in enumerate(spec["panels"], 1):
                canvas.cd(pad_index)
                legend = root.TLegend(0.52, 0.67, 0.89, 0.89)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                retained.append(legend)
                series = [(name, v) for name, v in panel["series"] if len(v) == len(panel["edges"]) - 1 and all(map(finite, v))]
                if not series:
                    root.TLatex().DrawLatexNDC(.12, .5, "unavailable: retained plot inputs missing")
                    continue
                all_values = [v for _, values in series for v in values]
                low, high = min(all_values + [0.]), max(all_values + [0.])
                span = max(high - low, 1e-12)
                for ordinal, (name, values) in enumerate(series):
                    graph = root.TGraph(len(values), array("d", [(a+b)/2 for a, b in zip(panel["edges"], panel["edges"][1:])]), array("d", values))
                    graph.SetTitle(panel["title"] + ";MM [GeV];" + panel["unit"])
                    graph.SetLineColor((1, 2, 4, 6, 8, 9, 46)[ordinal % 7])
                    graph.SetLineStyle(1 + ordinal // 7)
                    graph.SetLineWidth(2)
                    if ordinal == 0:
                        graph.SetMinimum(low - .1 * span)
                        graph.SetMaximum(high + .3 * span)
                    graph.Draw("AL" if ordinal == 0 else "L SAME")
                    legend.AddEntry(graph, name, "l")
                    retained.append(graph)
                legend.Draw()
                if spec.get("marker") is not None:
                    line = root.TLine(spec["marker"], low, spec["marker"], high)
                    line.SetLineStyle(2)
                    line.Draw()
                    retained.append(line)
            canvas.cd(6)
            text = root.TPaveText(.02, .02, .98, .98, "NDC")
            text.SetFillStyle(0)
            text.SetBorderSize(0)
            text.SetTextAlign(12)
            text.SetTextSize(.045)
            for line in [spec["title"], "diagnostic only / no production mutation"] + spec["text"]:
                # Long arrays live in JSON; bounded text cannot cover plot pads.
                for wrapped in textwrap.wrap(line[:320], 52):
                    text.AddText(wrapped)
            text.Draw()
            retained.append(text)
            canvas.Print(str(path))
        canvas.Print(str(path) + "]")
    finally:
        canvas.Close()
    return len(pages)


def _publish(path, payload):
    """Atomic no-replace publication; races cannot overwrite prior artifacts."""
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".pion-chain-", delete=False) as f:
        temporary = Path(f.name)
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    try:
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.is_symlink() or path.read_bytes() != payload:
                raise FileExistsError("pion chain artifact collision: " + str(path))
    finally:
        temporary.unlink()


def write_chain_artifacts(outpath, audit, renderer=render_companion):
    directory = Path(outpath).resolve(strict=True)
    if not directory.is_dir():
        raise ValueError("alignment_outpath is not an existing directory")
    stem = "pion_background_chain_v1_Q4p4W2p74_Center_lowe_" + audit["collection_identity"]
    json_path, pdf_path = (directory / (stem + ext) for ext in (".json", ".pdf"))
    if json_path.exists() or json_path.is_symlink():
        if json_path.is_symlink():
            raise FileExistsError("pion chain JSON symlink collision: " + str(json_path))
        existing = json.loads(json_path.read_bytes())
        if (set(existing) != {"schema_version", "status", "audit", "pdf"}
                or existing.get("schema_version") != SCHEMA or existing.get("audit") != audit):
            raise FileExistsError("pion chain JSON collision: " + str(json_path))
        pdf = existing.get("pdf", {})
        if pdf.get("status") == "available":
            if (existing.get("status") != "captured_pending_independent_validation"
                    or pdf.get("basename") != pdf_path.name or pdf.get("page_count") != len(page_specs(audit))
                    or pdf_path.is_symlink() or not pdf_path.is_file() or pdf_path.stat().st_size != pdf.get("bytes")
                    or hashlib.sha256(pdf_path.read_bytes()).hexdigest() != pdf.get("sha256")):
                raise FileExistsError("pion chain PDF mismatch: " + str(pdf_path))
        elif (existing.get("status") != "unavailable" or pdf.get("status") != "unavailable"
              or pdf_path.exists() or pdf_path.is_symlink()):
            raise FileExistsError("unexpected PDF beside unavailable receipt")
        return {"status": existing["status"], "json": str(json_path), "json_sha256": hashlib.sha256(json_path.read_bytes()).hexdigest(),
                "json_bytes": json_path.stat().st_size, "pdf": pdf, "reused": True}
    if pdf_path.exists() or pdf_path.is_symlink():
        raise FileExistsError("unreceipted prior PDF collision: " + str(pdf_path))
    with tempfile.TemporaryDirectory(dir=directory, prefix=".pion-chain-render-") as tmp:
        temporary_pdf = Path(tmp) / "companion.pdf"
        try:
            count = renderer(temporary_pdf, audit)
            pdf_bytes = temporary_pdf.read_bytes()
            if not pdf_bytes.startswith(b"%PDF-") or count != len(page_specs(audit)):
                raise ValueError("renderer produced invalid PDF/header/page receipt")
            pdf = {"status": "available", "basename": pdf_path.name, "bytes": len(pdf_bytes),
                   "sha256": hashlib.sha256(pdf_bytes).hexdigest(), "page_count": count,
                   "visual_acceptance": "NOT VERIFIED"}
        except Exception as exc:
            pdf_bytes = None
            pdf = unavailable("companion rendering failed: " + str(exc))
        status = "captured_pending_independent_validation" if pdf_bytes else "unavailable"
        receipt = {"schema_version": SCHEMA, "status": status, "audit": audit, "pdf": pdf}
        if pdf_bytes is not None:
            _publish(pdf_path, pdf_bytes)
        _publish(json_path, encoded(receipt))
    return {"status": status, "json": str(json_path), "json_sha256": hashlib.sha256(json_path.read_bytes()).hexdigest(),
            "json_bytes": json_path.stat().st_size, "pdf": pdf, "reused": False}


def observe_chain_audit(hist, inp, parents, shape_payload, outpath):
    if not is_target(inp, hist.get("phi_setting")):
        return None
    try:
        audit = build_chain_audit(hist, inp, parents, shape_payload)
        receipt = write_chain_artifacts(outpath, audit)
        print("[PION BACKGROUND CHAIN AUDIT] " + encoded(receipt).decode().strip())
        return receipt
    except Exception as exc:
        # Observability never changes production callbacks, policy or exceptions.
        receipt = unavailable("capture/publication failed: " + str(exc))
        print("[PION BACKGROUND CHAIN AUDIT] " + encoded(receipt).decode().strip())
        return receipt
