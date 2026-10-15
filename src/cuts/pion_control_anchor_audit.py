"""Read-only same-run neutron-anchor observations; no ROOT import or fit calls."""

from __future__ import annotations

import bisect
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
# Allow dot-product/serial-sum rounding, never a physics acceptance threshold.
REL_TOL = 1e-10
ABS_TOL = 1e-12


def unavailable(reason):
    return {"unavailable_reason": reason}


def json_safe(value):
    """Retain unavailable/nonfinite positions explicitly, including array entries."""
    if value is None:
        return unavailable("not retained by producer")
    if isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else unavailable("nonfinite: " + str(value))
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return unavailable("unsupported scalar type: " + type(value).__name__)


def _digest(value):
    return hashlib.sha256(json.dumps(json_safe(value), sort_keys=True,
                                     allow_nan=False).encode()).hexdigest()


def is_target(inp, phi_setting):
    # main.py supplies Q2/W string tokens and EPSSET='low' (not 'lowe').
    return (inp.get("ParticleType") == "kaon" and inp.get("Q2") == "4p4"
            and inp.get("W") == "2p74" and inp.get("EPSSET") == "low"
            and phi_setting == "Center")


def _close(a, b):
    return (math.isfinite(a) and math.isfinite(b)
            and math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL))


def reconstruct_moments(diagnostics):
    try:
        bins, y, sigma, template = [diagnostics[k] for k in
            ("fit_bin_indices", "fit_y", "fit_sigma", "fit_template")]
        if not bins or not (len(bins) == len(y) == len(sigma) == len(template)):
            raise ValueError("missing or unequal arrays")
        if not all(math.isfinite(v) for row in (y, sigma, template) for v in row):
            raise ValueError("nonfinite array entry")
        if not all(s > 0 for s in sigma):
            raise ValueError("nonpositive sigma")
        numerator = math.fsum((t / s) * (v / s) for v, s, t in zip(y, sigma, template))
        denominator = math.fsum((t / s) ** 2 for s, t in zip(sigma, template))
        matches = (_close(numerator, diagnostics["weighted_numerator"])
                   and _close(denominator, diagnostics["weighted_denominator"]))
        return {"status": "PASS" if matches else "FAIL", "numerator": numerator,
                "denominator": denominator, "reason": "moments match" if matches else "moment mismatch",
                "relative_tolerance": REL_TOL, "absolute_tolerance": ABS_TOL}
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        return {"status": "FAIL", "unavailable_reason": str(exc)}


def classify(step, solve, moments):
    # History retains the original solve; overlays can carry post-fit scaling.
    amplitude = solve.get("amplitude")
    if not isinstance(amplitude, (float, int)) or not math.isfinite(amplitude):
        return "not_determined"
    if solve.get("success") is False and amplitude == 0:
        message = str(solve.get("message", "")).lower()
        if any(token in message for token in
               ("missing target histogram", "missing simc template shape", "bin-count mismatch",
                "axis-range mismatch", "non-positive integral", "template has zero support",
                "no valid fit bins", "no valid amplitude-window bins")):
            return "zero_due_to_invalid_or_unsupported_template"
    if (solve.get("success") is not True or moments.get("status") != "PASS"
            or (step.get("fit_diagnostics") or {}).get("estimator") != "least_squares"):
        return "not_determined"
    n, d = moments["numerator"], moments["denominator"]
    if d > 0 and n <= 0 and amplitude == 0:
        return "zero_due_to_nonpositive_numerator"
    if d > 0 and n > 0 and amplitude > 0 and _close(amplitude, n / d):
        return "positive_staged_amplitude"
    return "not_determined"


def decompose_anchor(histogram, records, bins):
    """Use ROOT normal-bin [low, high) semantics; final upper edge is overflow."""
    if histogram is None or records is None or not bins:
        return unavailable("histogram, cache records or exact fitted bins absent")
    count = int(histogram.GetNbinsX())
    axis = histogram.GetXaxis()
    edges = [float(axis.GetBinLowEdge(i)) for i in range(1, count + 2)]
    if any(not isinstance(i, int) or not 1 <= i <= count for i in bins):
        return unavailable("invalid fitted normal-bin indices")
    rows = {i: {s: {"count": 0, "signed_sum": 0.0, "sumw2": 0.0}
                for s in SOURCES} for i in bins}
    invalid = []
    for ordinal, record in enumerate(records):
        if not record.get("nommcuts"):
            continue
        try:
            source = record["source_label"]
            x, coefficient = float(record["adj_MM"]), float(record["coefficient"])
            if source not in SOURCES or not all(map(math.isfinite, (x, coefficient))):
                raise ValueError("unknown source or nonfinite signed input")
            index = bisect.bisect_right(edges, x)
            if index in rows:
                cell = rows[index][source]
                cell["count"] += 1
                cell["signed_sum"] += coefficient
                cell["sumw2"] += coefficient * coefficient
        except (KeyError, TypeError, ValueError) as exc:
            invalid.append({"record_index": ordinal, "unavailable_reason": str(exc)})
    result = []
    for i, sources in rows.items():
        total = math.fsum(c["signed_sum"] for c in sources.values())
        variance = math.fsum(c["sumw2"] for c in sources.values())
        content, error = float(histogram.GetBinContent(i)), float(histogram.GetBinError(i))
        passed = _close(total, content) and _close(variance, error * error)
        result.append({"bin_index": i, "edges": edges[i-1:i+1], "sources": sources,
                       "combined_sum": total, "combined_error": math.sqrt(variance),
                       "hist_content": content, "hist_error": error,
                       "closure": "PASS" if passed else "FAIL",
                       "reason": "matched content and Sumw2" if passed else "content or Sumw2 mismatch"})
    return {"axis_edges": edges, "bins": result, "invalid_records": invalid,
            "records_fingerprint": _digest(records),
            "closure": "PASS" if not invalid and all(r["closure"] == "PASS" for r in result) else "FAIL",
            "relative_tolerance": REL_TOL, "absolute_tolerance": ABS_TOL}


def _simc_provenance(payload, index):
    component = ((payload or {}).get("components") or {}).get("pi_n") or {}
    shapes = (component.get("binned_shapes") or {}).get("t_bin{}".format(index + 1)) or {}
    phi = []
    for name, h in shapes.items():
        if h is not None:
            values = [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]
            phi.append({"phi_bin": name, "normal_bin_integral": math.fsum(values),
                        "positive_support": any(v > 0 for v in values),
                        "shape_fingerprint": _digest(values)})
    return {"component": "pi_n", "populated_phi_inputs": phi,
            "component_diagnostics": component.get("diagnostics"),
            "payload_diagnostics": ((payload or {}).get("diagnostics") or {}).get("pi_n"),
            "pre_normalization_phi_integrals": unavailable("NOT AVAILABLE: producer retains normalized phi shapes only"),
            "fallback": unavailable("resolved scope histogram has no retained per-t fallback decision; no resolver rerun")}


def _parent_audit(parent, cache, payload):
    index = parent["t_bin_index"]
    fit = parent.get("fit_result") or {}
    diag = ((fit.get("diagnostics") or {}).get("pion") or {})
    history = diag.get("staged_pass_history") or []
    overlays = fit.get("H_pion_fit_step_overlays") or []
    by_t = (cache or {}).get("by_t") or []
    entry = next((e for e in by_t if e.get("t_index") == index), {})
    passes = []
    for step in overlays:
        if step.get("component_name") != "pi_n":
            continue
        number = step.get("pass_index")
        solve = (history[number - 1].get("pi_n") or {}) if isinstance(number, int) and 0 < number <= len(history) else {}
        diagnostics = step.get("fit_diagnostics") or {}
        moments = reconstruct_moments(diagnostics)
        frozen = {k: step.get(k) for k in ("pass_index", "step_index", "component_name",
            "anchor_windows", "excluded_windows", "amplitude_mode", "amplitude_windows")}
        frozen.update(amplitude=solve.get("amplitude"), sigma=solve.get("sigma"),
                      overlay_amplitude=step.get("amplitude"),
                      overlay_raw_amplitude=step.get("raw_amplitude"),
                      postfit_scale_factor=step.get("postfit_scale_factor"))
        frozen.update(solve_outcome=solve, fit_diagnostics=diagnostics,
                      moment_reconstruction=moments, decision=classify(step, solve, moments),
                      aligned_template_support=diagnostics.get("fit_template"),
                      source_decomposition=decompose_anchor(parent.get("H_pion_control"),
                          entry.get("records"), diagnostics.get("fit_bin_indices")))
        passes.append(frozen)
    first = next((p for p in passes if p["pass_index"] == 1), None)
    histogram = parent.get("H_pion_control")
    histogram_observation = unavailable("parent pion-control histogram absent")
    if histogram is not None:
        n = histogram.GetNbinsX()
        histogram_observation = {
            "axis_edges": [float(histogram.GetXaxis().GetBinLowEdge(i)) for i in range(1, n + 2)],
            "content_error_fingerprint": _digest([
                [float(histogram.GetBinContent(i)), float(histogram.GetBinError(i))]
                for i in range(n + 2)]),
            "fingerprint_definition": "SHA256 JSON content/error pairs including underflow and overflow"}
    return {"identity": {k: parent.get(k) for k in ("pion_parent_id", "t_bin_index", "t_edges",
                "analysis_scope", "coordinate_fingerprint", "kaon_data_coordinate",
                "pion_control_input_fingerprint", "proton_producer_fingerprint", "pion_input_fingerprint")},
            "fit_window": {k: diag.get(k) for k in ("fit_min", "fit_max", "include_windows", "exclude_windows")},
            "first_pass": first if first else unavailable("first pi_n overlay absent"),
            "passes": passes, "staged_pass_history": history,
            "decision": first["decision"] if first else "not_determined",
            "current_first_pass_zero": first["amplitude"] == 0 if first and first["amplitude"] is not None else None,
            "historical_comparison": "current capture only; a positive first-pass amplitude does not reproduce the reported historical zero",
            "histogram_observation": histogram_observation,
            "joint_prior": {k: diag.get(k) for k in ("prior_scales", "prior_sigmas", "joint_refinement")},
            "fit_order": diag.get("fit_order"), "configured_pass_count": diag.get("n_passes"),
            "resolved_pion_config": (fit.get("resolved_subtraction_config") or {}).get("pion_control"),
            "alignment": fit.get("pion_component_alignment"),
            "amplitudes": {k: diag.get(k) for k in ("staged_amplitudes_raw", "staged_amplitudes_scaled",
                "refined_amplitudes", "refined_amplitudes_pre_postrefine_scale")},
            "final_B_n": fit.get("B_n"), "simc": _simc_provenance(payload, index),
            "source_accounting": entry.get("source_accounting"),
            "cache_record_fingerprint": (_digest(entry["records"]) if "records" in entry else None),
            "cache_coordinate_fingerprint": entry.get("coordinate_fingerprint")}


def build_anchor_audit(hist, inp, parents, payload):
    if not is_target(inp, hist.get("phi_setting")):
        return None
    cache = hist.get("_authoritative_pion_control_source_cache") or {}
    selected = [p for p in parents if p.get("t_bin_index") in (0, 1)]
    audit = {"schema_version": "pion_control_neutron_anchor_audit/v1",
             "diagnostic_only": True, "production_objects_mutated": False,
             "scope": "Q4p4W2p74/Center/lowe/canonical-t2", "epsilon_token": inp.get("EPSSET"),
             "observed_source_HEAD": "not_provided",
             "attribution": "current executed parent; not automatically the historical August 26 parent",
             "historical_provenance_match": "NOT VERIFIED",
             "source_audit": hist.get("pion_control_source_audit"),
             "cache_identity": {k: cache.get(k) for k in ("definition", "coordinate_fingerprint",
                 "kaon_data_coordinate", "pid_contract", "physical_pion_control_mask_fingerprint")},
             "parents": [_parent_audit(p, cache, payload) for p in selected]}
    if not any(p.get("t_bin_index") == 1 for p in selected):
        audit["unavailable_reason"] = "canonical t2 parent absent"
    return json_safe(audit)


def write_anchor_audit(outpath, audit):
    """Atomically publish complete bytes without replacing any existing file."""
    directory = Path(outpath)
    encoded = (json.dumps(audit, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    identity = _digest([p["identity"] for p in audit.get("parents", [])])
    path = directory / ("neutron_anchor_Q4p4W2p74_Center_lowe_canonical-t2_" + identity + ".json")
    with tempfile.NamedTemporaryFile(dir=directory, prefix=".neutron-anchor-", delete=False) as f:
        temporary = Path(f.name)
        f.write(encoded)
        f.flush()
        os.fsync(f.fileno())
    try:
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != encoded:
                raise FileExistsError("neutron-anchor audit collision: " + str(path))
    finally:
        temporary.unlink()
    return path


def observe_anchor_audit(hist, inp, parents, payload, outpath):
    if not is_target(inp, hist.get("phi_setting")):
        return
    try:
        audit = build_anchor_audit(hist, inp, parents, payload)
    except Exception as exc:
        audit = {"schema_version": "pion_control_neutron_anchor_audit/v1",
                 "diagnostic_only": True, "production_objects_mutated": False,
                 "unavailable_reason": "capture failed: " + str(exc), "parents": []}
    try:
        path = write_anchor_audit(outpath, audit)
        print("[NEUTRON ANCHOR AUDIT] " + str(path))
    except Exception as exc:
        # Diagnostic failure cannot change the production exception/fallback path.
        print("[NEUTRON ANCHOR AUDIT] unavailable: " + str(exc))
