"""Detached Part-2B relative HGCer refinement diagnostics.

Part 2 owns the frozen zero-photoelectron response map.  This module never
changes that map: it uses its transfer only as a relative t-delta score while
preserving the existing event-level pion subtraction in every canonical t bin.
"""

from __future__ import annotations

from array import array
from copy import deepcopy
import hashlib
import json
import math

import numpy as np
from scipy.optimize import minimize

from background_config import fingerprint_hgcer_pid_contract, hgcer_mask_accepts, normalize_hgcer_mask
from pion_component_subtraction import simc_shape_pion_weight_from_value
from pion_hgcer_transfer import (
    _compound_density,
    _compound_interval_probability,
    _nbinom_logpmf,
    _poisson_fit_normalization,
    _poisson_logpmf,
    fit_zero_photoelectron_response,
    poisson_zero_probability,
)

try:  # Pure-Python contracts remain usable on hosts without PyROOT.
    import ROOT
except ImportError:  # pragma: no cover - host capability
    ROOT = None

try:
    from root_histogram_ownership import clone_root_histogram
except ImportError:  # pragma: no cover - PyROOT integration only
    clone_root_histogram = None


ROOT_SAFE_REFINEMENT_LABELS = {
    "title": "PION HGCer RELATIVE REFINEMENT - PART 2B - NON-AUTHORITATIVE",
    "unavailable": "PART 2B HGCer REFINEMENT UNAVAILABLE - EXISTING PION SUBTRACTION UNCHANGED",
    "raw_score": "Raw zero-PE leakage score L (relative input, not pion normalization)",
    "correction": "Normalized relative HGCer correction C(t,delta)",
    "source": "Part 2B usable-score source",
    "uncertainty": "Part 2B relative correction total uncertainty",
}


def pion_hgcer_refinement_display_text(kind):
    return ROOT_SAFE_REFINEMENT_LABELS[str(kind)]


def _finite(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(child) for child in value]
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _fingerprint(value):
    encoded = json.dumps(_json_ready(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _root_integral(histogram):
    if histogram is None:
        return None
    try:
        return float(histogram.Integral(0, int(histogram.GetNbinsX()) + 1))
    except Exception:
        return None


def _root_histogram_closure(expected, actual, tolerance):
    result = {
        "available": expected is not None and actual is not None,
        "tolerance": float(tolerance), "passed": False,
        "content_max_abs_difference": None, "error_max_abs_difference": None,
        "integral_expected": _root_integral(expected), "integral_actual": _root_integral(actual),
    }
    if not result["available"]:
        result["reason"] = "histogram_missing"
        return result
    try:
        if int(expected.GetNbinsX()) != int(actual.GetNbinsX()):
            result["reason"] = "bin_count_mismatch"
            return result
        content_difference = error_difference = 0.0
        for index in range(0, int(expected.GetNbinsX()) + 2):
            content_difference = max(content_difference, abs(float(expected.GetBinContent(index)) - float(actual.GetBinContent(index))))
            error_difference = max(error_difference, abs(float(expected.GetBinError(index)) - float(actual.GetBinError(index))))
        result["content_max_abs_difference"] = content_difference
        result["error_max_abs_difference"] = error_difference
        result["passed"] = bool(content_difference <= tolerance and error_difference <= tolerance)
        if not result["passed"]:
            result["reason"] = "content_or_error_mismatch"
    except Exception as exc:
        result["reason"] = "closure_exception:{}".format(exc)
    return result


def _root_integral_window(histogram, bounds):
    if histogram is None or not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
        return None
    try:
        low, high = float(bounds[0]), float(bounds[1])
        if not low < high:
            return None
        axis = histogram.GetXaxis()
        return float(histogram.Integral(int(axis.FindBin(low)), int(axis.FindBin(high))))
    except Exception:
        return None


def _part2b_component_closure_metrics(parent, histograms):
    """Display-only window metrics from existing parent alignment provenance."""
    fit = (parent or {}).get("fit_result") or {}
    alignment = fit.get("pion_component_alignment") or {}
    components = alignment.get("components") or {}
    names = {
        "pi_n": "H_kaon_fit_pi_n_scaled",
        "pi_sidis": "H_kaon_fit_pi_sidis_scaled",
        "pi_delta": "H_kaon_fit_pi_delta_scaled",
    }
    metrics, unavailable = {}, []
    for component, histogram_name in names.items():
        windows = (components.get(component) or {}).get("effective_windows") or ()
        if component == "pi_delta" and len(windows) != 2:
            unavailable.append("pi_delta_sidebands_missing")
        if component != "pi_delta" and not windows:
            unavailable.append("{}_window_missing".format(component))
        for index, bounds in enumerate(windows):
            key = "{}_{}".format(component, "low" if component == "pi_delta" and index == 0 else "high" if component == "pi_delta" else "window")
            baseline = _root_integral_window((histograms or {}).get("baseline_pion"), bounds)
            refined = _root_integral_window((histograms or {}).get("refined_pion"), bounds)
            reference = _root_integral_window(fit.get(histogram_name), bounds)
            metrics[key] = {
                "component": component, "window": list(bounds), "baseline": baseline,
                "refined": refined, "reference": reference,
                "baseline_minus_reference": baseline - reference if baseline is not None and reference is not None else None,
                "refined_minus_reference": refined - reference if refined is not None and reference is not None else None,
            }
    return {"status": "available" if metrics and not unavailable else "unavailable", "reason": ";".join(unavailable) if unavailable else None, "metrics": metrics}


def _metric_update(metrics, weight, source):
    metrics["record_count"] += 1
    metrics["signed_sum"] += float(weight)
    metrics["absolute_weight_sum"] += abs(float(weight))
    metrics["sumw2"] += float(weight) * float(weight)
    entry = metrics["source_accounting"].setdefault(str(source), {
        "record_count": 0, "signed_sum": 0.0, "absolute_weight_sum": 0.0, "sumw2": 0.0,
    })
    entry["record_count"] += 1
    entry["signed_sum"] += float(weight)
    entry["absolute_weight_sum"] += abs(float(weight))
    entry["sumw2"] += float(weight) * float(weight)


def _metric_finalize(metrics):
    result = dict(metrics)
    denominator = float(result.get("sumw2", 0.0))
    numerator = float(result.get("absolute_weight_sum", 0.0))
    result["effective_entries"] = numerator * numerator / denominator if denominator > 0.0 else 0.0
    return result


def _empty_metrics():
    return {"record_count": 0, "signed_sum": 0.0, "absolute_weight_sum": 0.0, "sumw2": 0.0, "source_accounting": {}}


def _raw_cells(raw_transfer):
    return {
        (int(cell.get("t_index", -1)), int(cell.get("delta_index", -1))): cell
        for cell in (raw_transfer or {}).get("cells") or ()
    }


def _part1_prompt_samples(part1_payload, raw_cells):
    samples = {key: [] for key in raw_cells}
    for record in ((part1_payload or {}).get("records") or {}).get("pion") or ():
        if str(record.get("source_label") or "") != "prompt" or not bool(record.get("allcuts")):
            continue
        try:
            key = (int(record.get("canonical_t_index")), int(record.get("delta_index")))
        except (TypeError, ValueError):
            continue
        value = _finite(record.get("P_hgcer_npeSum"))
        if key in samples and value is not None:
            samples[key].append(value)
    return samples


def _score_from_parameters(family, parameters, contract):
    parameters = np.asarray(parameters, dtype=float)
    if len(parameters) == 0 or not np.all(np.isfinite(parameters)):
        return None
    mu = math.exp(float(parameters[0]))
    if not math.isfinite(mu) or mu <= 0.0:
        return None
    mask = (contract or {}).get("masks") or {}
    kaon_mask, control_mask = mask.get("kaon_tree") or {}, mask.get("physical_pion_control") or {}
    p0 = poisson_zero_probability(mu)
    try:
        if family == "zero_truncated_poisson":
            threshold = int(math.floor(float(control_mask["value"])))
            denominator = 1.0 - sum(math.exp(-mu + count * math.log(mu) - math.lgamma(count + 1)) for count in range(threshold + 1))
        elif family == "zero_truncated_negative_binomial":
            shape = math.exp(float(parameters[1]))
            nb_zero = math.exp(float(_nbinom_logpmf([0], mu, shape)[0]))
            threshold = int(math.floor(float(control_mask["value"])))
            positive_tail = 1.0 - sum(math.exp(float(_nbinom_logpmf([count], mu, shape)[0])) for count in range(threshold + 1))
            denominator = (1.0 - p0) * positive_tail / max(1.0e-15, 1.0 - nb_zero)
        elif family in {"zero_truncated_compound_poisson_gamma_gain", "zero_truncated_compound_poisson_exponential_gain"}:
            shape = 1.0 if family.endswith("exponential_gain") else math.exp(float(parameters[1]))
            scale = math.exp(float(parameters[-1]))
            denominator = _compound_interval_probability(mu, shape, scale, float(control_mask["value"]), 1.0e4)
        else:
            return None
    except (ArithmeticError, IndexError, KeyError, ValueError):
        return None
    return p0 / denominator if denominator > 0.0 and math.isfinite(denominator) else None


def _response_nll(values, fit, response_mask, config):
    """Exact existing response families, exposed only for Part-2B stability audit."""
    family = str((fit or {}).get("model_family") or "")
    lower, upper = [float(value) for value in (config or {}).get("fit_range") or (0.0, 20.0)]
    values = np.asarray([value for value in (_finite(item) for item in values) if value is not None], dtype=float)
    selected = values[(values >= lower) & (values <= upper)]
    if len(selected) == 0:
        return None, ()
    if family == "zero_truncated_poisson":
        selected = np.rint(selected).astype(int)
        if np.any(np.abs(selected - values[(values >= lower) & (values <= upper)]) > float(config.get("integer_tolerance", 1.0e-6))):
            return None, ()
        def nll(parameters):
            mu = math.exp(float(parameters[0]))
            norm = _poisson_fit_normalization(mu, response_mask, (lower, upper))
            return float(-np.sum(_poisson_logpmf(selected, mu)) + len(selected) * math.log(norm)) if norm > 0.0 else float("inf")
        return nll, ((-12.0, 8.0),)
    if family == "zero_truncated_negative_binomial":
        selected = np.rint(selected).astype(int)
        def nll(parameters):
            mean, shape = math.exp(float(parameters[0])), math.exp(float(parameters[1]))
            norm = sum(math.exp(float(_nbinom_logpmf([number], mean, shape)[0])) for number in range(max(0, int(math.ceil(lower))), int(math.floor(upper)) + 1) if hgcer_mask_accepts(number, response_mask))
            return float(-np.sum(_nbinom_logpmf(selected, mean, shape)) + len(selected) * math.log(norm)) if norm > 0.0 else float("inf")
        return nll, ((-12.0, 8.0), (-8.0, 12.0))
    if family in {"zero_truncated_compound_poisson_gamma_gain", "zero_truncated_compound_poisson_exponential_gain"}:
        if np.any(selected <= 0.0):
            return None, ()
        exponential = family.endswith("exponential_gain")
        def nll(parameters):
            mu = math.exp(float(parameters[0]))
            shape = 1.0 if exponential else math.exp(float(parameters[1]))
            scale = math.exp(float(parameters[-1]))
            norm = _compound_interval_probability(mu, shape, scale, lower, upper)
            density = _compound_density(selected, mu, shape, scale)
            return float(-np.sum(np.log(density)) + len(selected) * math.log(norm)) if norm > 0.0 and np.all(density > 0.0) and np.all(np.isfinite(density)) else float("inf")
        return nll, tuple([(-10.0, 8.0)] * (2 if exponential else 3))
    return None, ()


def _audit_profile(nll, parameters, bounds, config):
    parameters = np.asarray(parameters, dtype=float)
    if len(parameters) == 0:
        return {"status": "unavailable", "reason": "fit_parameters_missing"}
    center, best_nll = float(parameters[0]), float(nll(parameters))
    count = int(config.get("profile_points", 31))
    half_range = float(config.get("profile_log_mu_half_range", 3.0))
    grid = np.linspace(center - half_range, center + half_range, count)
    center_index = int(np.argmin(np.abs(grid - center)))
    order = [center_index]
    for offset in range(1, count):
        for candidate in (center_index + offset, center_index - offset):
            if 0 <= candidate < count:
                order.append(candidate)
    values, failed, warm = {}, [], {"lower": parameters[1:], "upper": parameters[1:]}
    for index in order:
        log_mu = float(grid[index])
        if len(parameters) == 1:
            value, success = float(nll(np.asarray([log_mu]))), True
        else:
            direction = "lower" if log_mu < center else "upper"
            initial = warm[direction]
            try:
                result = minimize(lambda rest: nll(np.r_[log_mu, rest]), initial, method="L-BFGS-B", bounds=list(bounds[1:]))
                success = bool(result.success and math.isfinite(float(result.fun)))
                value = float(result.fun) if success else float("inf")
                if success:
                    warm[direction] = np.asarray(result.x, dtype=float)
            except (ArithmeticError, FloatingPointError, ValueError):
                value, success = float("inf"), False
        values[index] = value
        if not success or not math.isfinite(value):
            failed.append(index)
    ordered_values = [values.get(index, float("inf")) for index in range(count)]
    target, lower, upper = best_nll + float(config.get("seed_tie_delta_nll", 0.5)), None, None
    for index in range(count - 1):
        left, right = ordered_values[index] - target, ordered_values[index + 1] - target
        if not math.isfinite(left) or not math.isfinite(right) or left * right > 0.0:
            continue
        crossing = float(grid[index]) if right == left else float(grid[index] + (grid[index + 1] - grid[index]) * (-left) / (right - left))
        if crossing < center and lower is None:
            lower = crossing
        if crossing > center and upper is None:
            upper = crossing
    return {
        "status": "two_sided" if lower is not None and upper is not None and not failed else "one_sided_or_failed",
        "grid_log_mu": list(map(float, grid)),
        "grid_nll": [value if math.isfinite(value) else None for value in ordered_values],
        "failed_point_indices": failed,
        "log_mu_interval": [lower, upper],
    }


def audit_part2b_response_score(values, fit, contract, raw_config, refinement_config):
    """Independent deterministic multistart/profile usability audit for Part 2B."""
    if not isinstance(fit, dict) or fit.get("fit_status") != "fit_valid":
        return {"status": "ineligible", "reason": "raw_fit_not_valid"}
    parameters = np.asarray(fit.get("fit_parameters_log") or (), dtype=float)
    response_mask = normalize_hgcer_mask(fit.get("response_mask") or ((contract or {}).get("masks") or {}).get("pion_tree") or {}, name="Part-2B response mask")
    nll, bounds = _response_nll(values, fit, response_mask, raw_config)
    if nll is None or len(parameters) != len(bounds):
        return {"status": "ineligible", "reason": "response_nll_unavailable"}
    family = str(fit.get("model_family"))
    observed_mean = max(float(np.mean(values)), 0.1) if values else 0.1
    if len(parameters) == 1:
        starts = [[math.log(mu)] for mu in (0.25, 0.75, 1.5, 3.0, 6.0, 12.0)]
    elif family == "zero_truncated_negative_binomial":
        starts = [[math.log(mu), math.log(shape)] for mu in (0.25, 1.0, 4.0, 12.0) for shape in (0.5, 2.0, 8.0)]
    else:
        exponential = family.endswith("exponential_gain")
        starts = []
        for mu in (0.25, 0.75, 2.0, 6.0, 12.0):
            gain_scales = (0.15, 0.6, 2.0) if exponential else (0.2, 0.8, 3.0)
            for scale in gain_scales:
                if exponential:
                    starts.append([math.log(mu), math.log(max(observed_mean / max(mu, 0.25) * scale, 0.03))])
                else:
                    for shape in (0.4, 1.5, 6.0):
                        starts.append([math.log(mu), math.log(shape), math.log(max(observed_mean / max(mu * shape, 0.25) * scale, 0.03))])
    candidates = [{"kind": "frozen", "initial": list(map(float, parameters)), "parameters": list(map(float, parameters)), "nll": float(nll(parameters)), "success": True}]
    for start in starts:
        try:
            result = minimize(nll, np.asarray(start, dtype=float), method="L-BFGS-B", bounds=list(bounds))
            success = bool(result.success and math.isfinite(float(result.fun)))
            candidate = {"kind": "second_seed_grid", "initial": list(map(float, start)), "parameters": list(map(float, result.x)) if success else None, "nll": float(result.fun) if success else None, "success": success, "message": str(result.message)}
        except (ArithmeticError, FloatingPointError, ValueError) as exc:
            candidate = {"kind": "second_seed_grid", "initial": list(map(float, start)), "parameters": None, "nll": None, "success": False, "message": str(exc)}
        if candidate["parameters"] is not None:
            candidate["mu"] = math.exp(float(candidate["parameters"][0]))
            candidate["L"] = _score_from_parameters(family, candidate["parameters"], contract)
        candidates.append(candidate)
    finite = [candidate for candidate in candidates if candidate.get("success") and _finite(candidate.get("nll")) is not None and _finite(candidate.get("L", _score_from_parameters(family, parameters, contract))) is not None]
    best_nll = min((float(candidate["nll"]) for candidate in finite), default=None)
    tied = [candidate for candidate in finite if best_nll is not None and float(candidate["nll"]) <= best_nll + float(refinement_config["seed_tie_delta_nll"])]
    tied_scores = [float(candidate.get("L", _score_from_parameters(family, parameters, contract))) for candidate in tied]
    score_spread = (max(tied_scores) - min(tied_scores)) / max(abs(float(np.mean(tied_scores))), 1.0e-15) if len(tied_scores) >= 2 else None
    profile = _audit_profile(nll, parameters, bounds, refinement_config)
    status = "stable"
    reason = None
    if len(tied_scores) < 2:
        status, reason = "ineligible", "insufficient_converged_second_seed_candidates"
    elif score_spread is not None and score_spread > float(refinement_config["seed_ambiguity_relative_score_spread"]):
        status, reason = "ambiguous", "tied_seed_scores_materially_different"
    elif profile.get("status") != "two_sided":
        status, reason = "ineligible", "profile_not_strongly_identifiable"
    return {
        "status": status, "reason": reason, "family": family,
        "winning_seed": min(finite, key=lambda candidate: float(candidate["nll"])) if finite else None,
        "candidate_count": len(candidates), "converged_candidate_count": len(finite),
        "best_nll": best_nll, "next_best_nll": sorted(float(candidate["nll"]) for candidate in finite)[1] if len(finite) > 1 else None,
        "tied_candidate_count": len(tied), "tied_score_relative_spread": score_spread,
        "candidates": candidates, "profile": profile,
    }


def build_part2b_usable_score_map(raw_transfer, part1_payload, refinement_config):
    """Rebuild a local same-t score/fallback view without mutating Part 2."""
    raw_cells = _raw_cells(raw_transfer)
    contract = deepcopy((raw_transfer or {}).get("pid_contract") or fingerprint_hgcer_pid_contract())
    raw_config = deepcopy((raw_transfer or {}).get("config") or {})
    samples = _part1_prompt_samples(part1_payload, raw_cells)
    audits, direct, output = {}, {}, {}
    for key, raw_cell in raw_cells.items():
        fit, solution = raw_cell.get("fit") or {}, raw_cell.get("solution") or {}
        audit = audit_part2b_response_score(samples.get(key) or [], fit, contract, raw_config, refinement_config)
        audits[key] = audit
        transfer = _finite(solution.get("transfer"))
        if solution.get("solution_source") == "direct" and audit.get("status") == "stable" and transfer is not None:
            direct[key] = {"L": transfer, "source": "stable_direct", "raw_solution_source": "direct", "raw_fit": fit}
    t_indices = sorted({key[0] for key in raw_cells})
    for t_index in t_indices:
        t_keys = sorted(key for key in raw_cells if key[0] == t_index)
        stable_keys = sorted(key for key in direct if key[0] == t_index)
        pooled_values = [value for key in t_keys for value in samples.get(key) or () if (raw_cells[key].get("support_class") == "supported")]
        pooled_fit = fit_zero_photoelectron_response(pooled_values, response_mask=contract["masks"]["pion_tree"], contract=contract, config=raw_config) if pooled_values else {"fit_status": "fit_invalid", "reason": "no_same_t_prompt_samples"}
        pooled_audit = audit_part2b_response_score(pooled_values, pooled_fit, contract, raw_config, refinement_config) if pooled_fit.get("fit_status") == "fit_valid" else {"status": "ineligible", "reason": pooled_fit.get("reason")}
        for key in t_keys:
            if key in direct:
                output[key] = {**direct[key], "audit": audits[key]}
                continue
            delta_index = key[1]
            lower = [candidate for candidate in stable_keys if candidate[1] < delta_index]
            upper = [candidate for candidate in stable_keys if candidate[1] > delta_index]
            if lower and upper:
                left, right = lower[-1], upper[0]
                fraction = float(delta_index - left[1]) / float(right[1] - left[1])
                output[key] = {"L": (1.0 - fraction) * direct[left]["L"] + fraction * direct[right]["L"], "source": "same_t_delta_bracket", "bracket": [left[1], right[1]], "fraction": fraction, "audit": audits[key]}
            elif stable_keys and delta_index in (t_keys[0][1], t_keys[-1][1]):
                neighbour = stable_keys[0] if delta_index == t_keys[0][1] else stable_keys[-1]
                output[key] = {"L": direct[neighbour]["L"], "source": "same_t_delta_edge", "edge_neighbour": neighbour[1], "audit": audits[key]}
            elif pooled_fit.get("fit_status") == "fit_valid" and pooled_audit.get("status") == "stable" and _finite(pooled_fit.get("transfer")) is not None:
                output[key] = {"L": float(pooled_fit["transfer"]), "source": "same_t_pooled", "pooled_fit": pooled_fit, "pooled_audit": pooled_audit, "audit": audits[key]}
            else:
                output[key] = {"L": None, "source": "identity_no_hgcer_refinement", "reason": audits[key].get("reason") or pooled_audit.get("reason"), "audit": audits[key], "pooled_audit": pooled_audit}
    return {"raw_map_fingerprint": (raw_transfer or {}).get("map_fingerprint"), "audits": audits, "scores": output}


def calculate_part2b_relative_corrections(cell_rows, config):
    """Pure-Python canonical-t normalization with explicit identity fallback."""
    rows = [dict(row) for row in cell_rows or ()]
    usable = [row for row in rows if _finite(row.get("L")) is not None]
    if not usable:
        return {"status": "identity_fallback", "reason": "no_usable_relative_scores", "alpha": 1.0, "f_N": None, "f_D": None, "cells": [{**row, "C": 1.0, "correction_source": "identity_no_hgcer_refinement"} for row in rows]}
    numerator = sum(float(row.get("B", 0.0)) for row in usable)
    denominator = sum(float(row.get("B", 0.0)) * float(row["L"]) for row in usable)
    abs_numerator = sum(abs(float(row.get("B", 0.0))) for row in usable)
    abs_denominator = sum(abs(float(row.get("B", 0.0)) * float(row["L"])) for row in usable)
    f_n = abs(numerator) / abs_numerator if abs_numerator > 0.0 else 0.0
    f_d = abs(denominator) / abs_denominator if abs_denominator > 0.0 else 0.0
    floor = float(config["signed_support_ratio_floor"])
    stable = all(math.isfinite(value) for value in (numerator, denominator, f_n, f_d)) and denominator != 0.0 and numerator * denominator > 0.0 and f_n >= floor and f_d >= floor
    if not stable:
        return {"status": "identity_fallback", "reason": "signed_normalization_unstable", "alpha": 1.0, "N": numerator, "D": denominator, "f_N": f_n, "f_D": f_d, "cells": [{**row, "C": 1.0, "correction_source": "identity_signed_normalization_unstable"} for row in rows]}
    alpha = numerator / denominator
    result_cells = []
    for row in rows:
        if _finite(row.get("L")) is None:
            result_cells.append({**row, "C": 1.0, "correction_source": "identity_no_hgcer_refinement"})
        else:
            result_cells.append({**row, "C": alpha * float(row["L"]), "correction_source": row.get("source")})
    return {"status": "available", "reason": None, "alpha": alpha, "N": numerator, "D": denominator, "f_N": f_n, "f_D": f_d, "cells": result_cells}


def _score_uncertainties(score, raw_cell, raw_cells=None):
    solution = (raw_cell or {}).get("solution") or {}
    fit = (raw_cell or {}).get("fit") or {}
    source = score.get("source")
    if source == "same_t_pooled":
        fit = score.get("pooled_fit") or {}
        solution = fit
    if source in {"stable_direct", "same_t_pooled"}:
        return {"stat": _finite(solution.get("transfer_statistical_uncertainty")), "model": _finite(solution.get("transfer_model_uncertainty")), "fallback": 0.0}
    if source == "same_t_delta_bracket":
        try:
            t_index = int(raw_cell.get("t_index"))
            left, right = [int(value) for value in score.get("bracket") or ()]
            fraction = float(score.get("fraction"))
            left_sigma = _finite(((raw_cells or {}).get((t_index, left)) or {}).get("solution", {}).get("transfer_total_uncertainty"))
            right_sigma = _finite(((raw_cells or {}).get((t_index, right)) or {}).get("solution", {}).get("transfer_total_uncertainty"))
            if left_sigma is not None and right_sigma is not None:
                return {"stat": None, "model": None, "fallback": math.sqrt((1.0 - fraction) ** 2 * left_sigma ** 2 + fraction ** 2 * right_sigma ** 2)}
        except (TypeError, ValueError):
            pass
        return {"stat": None, "model": None, "fallback": None}
    if source == "same_t_delta_edge":
        try:
            neighbour = (raw_cells or {}).get((int(raw_cell.get("t_index")), int(score.get("edge_neighbour")))) or {}
            return {"stat": None, "model": None, "fallback": _finite((neighbour.get("solution") or {}).get("transfer_total_uncertainty"))}
        except (TypeError, ValueError):
            return {"stat": None, "model": None, "fallback": None}
    return {"stat": None, "model": None, "fallback": None}


def _propagate_correction_uncertainties(correction, raw_cells):
    cells = correction.get("cells") or []
    usable = [index for index, cell in enumerate(cells) if _finite(cell.get("L")) is not None and correction.get("status") == "available"]
    components = {"stat": [], "model": [], "fallback": []}
    for cell in cells:
        score = cell.get("score") or {}
        components_for_cell = _score_uncertainties(score, raw_cells.get((int(cell["t_index"]), int(cell["delta_index"]))) or {}, raw_cells)
        for key in components:
            components[key].append(components_for_cell[key])
    if correction.get("status") != "available":
        return {
            "cross_cell_covariance": "unavailable_identity_not_measured", "components": components,
            "covariance": None,
            "per_cell": [{
                "sigma_C_stat": None, "sigma_C_model": None, "sigma_C_fallback": None,
                "sigma_C_total": None, "uncertainty_status": "not_measured_identity",
            } for _ in cells],
        }
    B = np.asarray([float(cells[index].get("B", 0.0)) for index in usable], dtype=float)
    C = np.asarray([float(cells[index]["C"]) for index in usable], dtype=float)
    alpha, denominator = float(correction["alpha"]), float(correction["D"])
    jacobian = np.zeros((len(usable), len(usable)), dtype=float)
    for row_index in range(len(usable)):
        for column_index in range(len(usable)):
            jacobian[row_index, column_index] = alpha * (1.0 if row_index == column_index else 0.0) - C[row_index] * B[column_index] / denominator
    covariance, per_cell = {}, [dict() for _ in cells]
    for component, values in components.items():
        variances = np.asarray([0.0 if values[index] is None else float(values[index]) ** 2 for index in usable], dtype=float)
        matrix = jacobian @ np.diag(variances) @ jacobian.T
        covariance[component] = matrix.tolist()
        for output_index, cell_index in enumerate(usable):
            per_cell[cell_index]["sigma_C_{}".format(component)] = math.sqrt(max(0.0, float(matrix[output_index, output_index]))) if values[usable[output_index]] is not None else None
    for index, fields in enumerate(per_cell):
        values = [fields.get("sigma_C_{}".format(component)) for component in covariance]
        fields["sigma_C_total"] = math.sqrt(sum(value * value for value in values if value is not None)) if any(value is not None for value in values) else None
        fields["uncertainty_status"] = "available_assumed_diagonal" if fields["sigma_C_total"] is not None else "not_measured"
    return {"cross_cell_covariance": "unavailable_assumed_diagonal_raw_scores", "components": components, "covariance": covariance, "per_cell": per_cell}


def _audit_control_records(raw_transfer, cache):
    contract = (raw_transfer or {}).get("pid_contract") or fingerprint_hgcer_pid_contract()
    mask = (contract.get("masks") or {}).get("physical_pion_control") or {}
    fp, coordinate = str(contract.get("fingerprint") or ""), str((raw_transfer or {}).get("coordinate_fingerprint") or "")
    if not fp or str((cache or {}).get("physical_pion_control_mask_fingerprint") or "") != fp:
        return {"status": "unavailable", "reason": "physical_control_mask_fingerprint_mismatch"}
    if normalize_hgcer_mask((cache or {}).get("physical_pion_control_mask") or {}, name="Part-2B cache mask") != normalize_hgcer_mask(mask, name="Part-2B map mask"):
        return {"status": "unavailable", "reason": "physical_control_mask_expression_mismatch"}
    if not coordinate or str((cache or {}).get("coordinate_fingerprint") or "") != coordinate:
        return {"status": "unavailable", "reason": "coordinate_fingerprint_mismatch"}
    t_edges = tuple(float(value) for value in (raw_transfer or {}).get("t_edges") or ())
    delta_edges = tuple(float(value) for value in (raw_transfer or {}).get("delta_edges") or ())
    expected_sources = {"prompt", "rand", "dummy", "dummy_rand"}
    failures, source_checks = {}, {}
    for t_index, section in enumerate((cache or {}).get("by_t") or ()):
        observed = {}
        for record in (section or {}).get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            try:
                record_t = int(record.get("t_index"))
                delta_index = int(record.get("delta_index"))
                adj_t, delta = float(record.get("adj_t")), float(record.get("ssdelta"))
                coefficient = float(record.get("coefficient"))
            except (TypeError, ValueError):
                failures[t_index] = "physical_control_record_coordinates_or_coefficient_missing"
                break
            expected_t = t_index
            in_t = len(t_edges) == 0 or (t_edges[expected_t] <= adj_t <= t_edges[expected_t + 1])
            in_delta = len(delta_edges) == 0 or (0 <= delta_index < len(delta_edges) - 1 and delta_edges[delta_index] <= delta <= delta_edges[delta_index + 1])
            if (
                record_t != expected_t
                or not in_t
                or not in_delta
                or str(record.get("rf_state") or "") != "noRF"
                or str(record.get("coordinate_fingerprint") or "") != coordinate
                or not hgcer_mask_accepts(record.get("P_hgcer_npeSum"), mask)
                or not math.isfinite(coefficient)
            ):
                failures[t_index] = "physical_control_record_contract_failed"
                break
            source = str(record.get("source_label") or "")
            if source not in expected_sources:
                failures[t_index] = "physical_control_record_source_unknown"
                break
            entry = observed.setdefault(source, {"record_count": 0, "signed_coefficient_sum": 0.0})
            entry["record_count"] += 1
            entry["signed_coefficient_sum"] += coefficient
        declared = (section or {}).get("source_accounting") or {}
        for source, actual in observed.items():
            expected = declared.get(source) or {}
            declared_count = expected.get("allcuts_records", expected.get("record_count"))
            declared_coefficient = expected.get("coefficient")
            if declared_count is not None and int(declared_count) != int(actual["record_count"]):
                failures[t_index] = "physical_control_source_record_count_mismatch"
            if declared_coefficient is not None and not math.isclose(
                float(declared_coefficient) * actual["record_count"],
                actual["signed_coefficient_sum"], rel_tol=0.0, abs_tol=1.0e-12,
            ):
                failures[t_index] = "physical_control_source_coefficient_mismatch"
        for source, expected in declared.items():
            declared_count = expected.get("allcuts_records", expected.get("record_count", 0))
            if int(declared_count or 0) != int((observed.get(source) or {}).get("record_count", 0)):
                failures[t_index] = "physical_control_declared_source_population_mismatch"
        source_checks[t_index] = observed
    return {
        "status": "unavailable" if failures else "available",
        "reason": next(iter(failures.values()), None),
        "t_failures": failures, "source_accounting": source_checks,
    }


def _unavailable_refinement(reason, config, raw_transfer=None):
    return {
        "status": "unavailable", "reason": str(reason), "non_authoritative": True,
        "production_pion_subtraction_unchanged": True, "config": deepcopy(config or {}),
        "raw_hgcer_score_fingerprint": (raw_transfer or {}).get("map_fingerprint"),
        "Part3_eligibility": "unavailable",
    }


def pion_hgcer_refinement_fingerprint_payload(raw_transfer, parents, cache, config, cells):
    """Fingerprint only Part-2B physics/provenance, deliberately excluding phi."""
    fingerprint_config = {
        key: value for key, value in dict(config or {}).items()
        if "phi" not in str(key).lower()
    }
    return {
        "raw_hgcer_score_fingerprint": (raw_transfer or {}).get("map_fingerprint"),
        "parents": [str((parents.get(index) or {}).get("pion_parent_id") or "") for index in range(max(0, len((raw_transfer or {}).get("t_edges") or ()) - 1))],
        "t_edges": (raw_transfer or {}).get("t_edges"), "delta_edges": (raw_transfer or {}).get("delta_edges"),
        "coordinate_fingerprint": (raw_transfer or {}).get("coordinate_fingerprint"),
        "physical_control_mask_fingerprint": (cache or {}).get("physical_pion_control_mask_fingerprint"),
        "config": fingerprint_config, "cells": list(cells or ()),
    }


def build_pion_hgcer_relative_refinement(raw_transfer, part1_payload, pion_control_cache, t_bin_parent_results, proton_t_products, *, config):
    """Build detached Part-2B products; this function never changes production objects."""
    config = deepcopy(config or {})
    if not isinstance(raw_transfer, dict) or raw_transfer.get("status") != "available" or not raw_transfer.get("frozen"):
        return _unavailable_refinement("frozen_part2_raw_score_unavailable", config, raw_transfer)
    if ROOT is None or clone_root_histogram is None:
        return _unavailable_refinement("PyROOT_unavailable_for_part2b_products", config, raw_transfer)
    control_audit = _audit_control_records(raw_transfer, pion_control_cache)
    if control_audit.get("status") != "available":
        return _unavailable_refinement(control_audit.get("reason"), config, raw_transfer)
    raw_cells = _raw_cells(raw_transfer)
    score_view = build_part2b_usable_score_map(raw_transfer, part1_payload, config)
    parents = {int(parent.get("t_bin_index")): parent for parent in t_bin_parent_results or () if _finite(parent.get("t_bin_index")) is not None}
    products, cache_by_t = list(proton_t_products or ()), list((pion_control_cache or {}).get("by_t") or ())
    t_count = len((raw_transfer or {}).get("t_edges") or ()) - 1
    if len(products) != t_count or len(cache_by_t) != t_count:
        return _unavailable_refinement("canonical_t_parent_or_cache_count_mismatch", config, raw_transfer)
    output_cells, t_results, detached = [], [], {}
    for t_index in range(t_count):
        parent, section, proton_product = parents.get(t_index) or {}, cache_by_t[t_index] or {}, products[t_index] or {}
        application = parent.get("final_diagnostic_application_result") or {}
        expected = application.get("H_pion_subtraction_template_MM")
        reference, weights = application.get("H_pion_control_model"), application.get("weights")
        host = ((proton_product.get("cleaned_targets_pre_rf") or {}).get("h_mm"))
        if expected is None or reference is None or weights is None or host is None:
            t_results.append({"t_index": t_index, "status": "unavailable", "reason": "current_allcuts_event_template_or_noRF_cut_host_missing"})
            continue
        baseline = clone_root_histogram(expected, scope="pion_hgcer_part2b", role="baseline_pion_allcuts", name="H_MM_part2b_baseline_pion_t{}".format(t_index + 1), reset=True, sumw2=True)
        cell_hists, cell_metrics = {}, {}
        for key in sorted(key for key in raw_cells if key[0] == t_index):
            delta_index = key[1]
            cell_hists[delta_index] = clone_root_histogram(expected, scope="pion_hgcer_part2b", role="baseline_cell_pion", name="H_MM_part2b_baseline_pion_t{}_d{}".format(t_index + 1, delta_index + 1), reset=True, sumw2=True)
            cell_metrics[delta_index] = _empty_metrics()
        for record in section.get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            try:
                delta_index = int(record.get("delta_index"))
            except (TypeError, ValueError):
                continue
            if delta_index not in cell_metrics:
                continue
            weight = float(record["coefficient"]) * simc_shape_pion_weight_from_value(float(record["adj_MM"]), reference, weights)
            if not math.isfinite(weight):
                continue
            baseline.Fill(float(record["adj_MM"]), weight)
            cell_hists[delta_index].Fill(float(record["adj_MM"]), weight)
            _metric_update(cell_metrics[delta_index], weight, record.get("source_label"))
        baseline_closure = _root_histogram_closure(expected, baseline, config["parent_integral_closure_tolerance"])
        rows = []
        for key in sorted(key for key in raw_cells if key[0] == t_index):
            delta_index = key[1]
            score = score_view["scores"].get(key) or {"L": None, "source": "identity_no_hgcer_refinement"}
            metric = _metric_finalize(cell_metrics[delta_index])
            rows.append({"t_index": t_index, "delta_index": delta_index, "B": metric["signed_sum"], "L": score.get("L"), "source": score.get("source"), "score": score, "baseline": metric})
        correction = calculate_part2b_relative_corrections(rows, config)
        uncertainty = _propagate_correction_uncertainties(correction, raw_cells)
        if not baseline_closure.get("passed"):
            correction = {"status": "identity_fallback", "reason": "baseline_pion_event_template_closure_failed", "alpha": 1.0, "f_N": None, "f_D": None, "cells": [{**row, "C": 1.0, "correction_source": "identity_baseline_closure_failed"} for row in rows]}
            uncertainty = _propagate_correction_uncertainties(correction, raw_cells)
        refined = clone_root_histogram(expected, scope="pion_hgcer_part2b", role="refined_pion_allcuts", name="H_MM_part2b_refined_pion_t{}".format(t_index + 1), reset=True, sumw2=True)
        correction_lookup = {int(cell["delta_index"]): cell for cell in correction["cells"]}
        for record in section.get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            try:
                correction_cell = correction_lookup[int(record.get("delta_index"))]
            except (KeyError, TypeError, ValueError):
                continue
            weight = float(record["coefficient"]) * simc_shape_pion_weight_from_value(float(record["adj_MM"]), reference, weights)
            if math.isfinite(weight):
                refined.Fill(float(record["adj_MM"]), weight * float(correction_cell["C"]))
        integral_closure = _root_histogram_closure(baseline, refined, config["parent_integral_closure_tolerance"])
        if bool(config.get("require_parent_integral_preservation", True)) and not integral_closure.get("passed"):
            correction = {
                "status": "identity_fallback", "reason": "parent_refined_integral_closure_failed",
                "alpha": 1.0, "f_N": correction.get("f_N"), "f_D": correction.get("f_D"),
                "N": correction.get("N"), "D": correction.get("D"),
                "cells": [{**row, "C": 1.0, "correction_source": "identity_parent_integral_closure_failed"} for row in rows],
            }
            uncertainty = _propagate_correction_uncertainties(correction, raw_cells)
            refined = clone_root_histogram(
                baseline, scope="pion_hgcer_part2b", role="refined_pion_identity",
                name="H_MM_part2b_refined_pion_t{}".format(t_index + 1), reset=False, sumw2=True,
            )
            integral_closure = _root_histogram_closure(
                baseline, refined, config["parent_integral_closure_tolerance"],
            )
        baseline_clean = clone_root_histogram(host, scope="pion_hgcer_part2b", role="baseline_clean", name="H_MM_part2b_baseline_clean_t{}".format(t_index + 1), reset=False, sumw2=True)
        refined_clean = clone_root_histogram(host, scope="pion_hgcer_part2b", role="refined_clean", name="H_MM_part2b_refined_clean_t{}".format(t_index + 1), reset=False, sumw2=True)
        baseline_clean.Add(baseline, -1.0); refined_clean.Add(refined, -1.0)
        for cell in correction["cells"]:
            fields = (uncertainty.get("per_cell") or [{} for _ in correction["cells"]])[correction["cells"].index(cell)] if uncertainty.get("per_cell") else {}
            output_cells.append({**cell, **fields, "raw_score": score_view["scores"].get((t_index, int(cell["delta_index"]))) or {}, "baseline_histogram_role": "detached"})
        detached[t_index] = {"host": host, "baseline_pion": baseline, "refined_pion": refined, "baseline_clean": baseline_clean, "refined_clean": refined_clean, "cell_baseline": cell_hists}
        t_results.append({"t_index": t_index, "status": "available" if baseline_closure.get("passed") and integral_closure.get("passed") else "unavailable", "reason": None if baseline_closure.get("passed") and integral_closure.get("passed") else correction.get("reason"), "alpha": correction.get("alpha"), "N": correction.get("N"), "D": correction.get("D"), "f_N": correction.get("f_N"), "f_D": correction.get("f_D"), "baseline_event_template_closure": baseline_closure, "parent_integral_closure": integral_closure, "correction_status": correction.get("status"), "correction_reason": correction.get("reason"), "component_closure_metrics": _part2b_component_closure_metrics(parent, detached[t_index])})
    global_histograms = {}
    if detached:
        first = detached[sorted(detached)[0]]
        for role in ("host", "baseline_pion", "refined_pion", "baseline_clean", "refined_clean"):
            aggregate = clone_root_histogram(first[role], scope="pion_hgcer_part2b", role="global_{}_strict_sum".format(role), name="H_MM_part2b_{}_global".format(role), reset=True, sumw2=True)
            for entry in detached.values():
                aggregate.Add(entry[role])
            global_histograms[role] = aggregate
    essentials = [{"t_index": cell["t_index"], "delta_index": cell["delta_index"], "B": cell["B"], "L": cell.get("L"), "C": cell.get("C"), "source": cell.get("source"), "correction_source": cell.get("correction_source")} for cell in output_cells]
    refinement_fingerprint = _fingerprint(pion_hgcer_refinement_fingerprint_payload(
        raw_transfer, parents, pion_control_cache, config, essentials,
    ))
    return {"status": "available", "non_authoritative": True, "production_pion_subtraction_unchanged": True, "config": config, "raw_hgcer_score_fingerprint": raw_transfer.get("map_fingerprint"), "pion_parent_ids": [str((parents.get(index) or {}).get("pion_parent_id") or "") for index in range(t_count)], "coordinate_fingerprint": raw_transfer.get("coordinate_fingerprint"), "physical_control_mask_fingerprint": (pion_control_cache or {}).get("physical_pion_control_mask_fingerprint"), "t_edges": list(raw_transfer.get("t_edges") or ()), "delta_edges": list(raw_transfer.get("delta_edges") or ()), "control_audit": control_audit, "score_view": score_view, "cells": tuple(output_cells), "t_results": tuple(t_results), "histograms": detached, "global_histograms": global_histograms, "refinement_fingerprint": refinement_fingerprint, "frozen": True, "Part3_eligibility": "review_only"}


def serialize_pion_hgcer_refinement(payload):
    result = {key: value for key, value in (payload or {}).items() if key not in {"histograms", "global_histograms"}}
    return _json_ready(result)


def write_pion_hgcer_refinement_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(serialize_pion_hgcer_refinement(payload), handle, sort_keys=True, indent=2, allow_nan=False)
    return path


# This catalog is intentionally separate from the frozen Part-2 renderer.  A
# Part-2B visual is either rendered with its stated physics primitive or is an
# explicit unavailable status page; a text-only substitute is never valid for
# a map, score, correction, MM, or closure page.
PART2B_PAGE_SPECS = {
    "audit": {"kind": "text", "renderer": "audit", "required_roles": ()},
    "unavailable": {"kind": "text", "renderer": "status", "required_roles": ()},
    "raw_score_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("raw_L_map",)},
    "correction_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("correction_C_map",)},
    "source_map": {"kind": "graphical", "renderer": "categorical_map", "required_roles": ("score_source_map",)},
    "uncertainty_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("correction_uncertainty_map",)},
    "normalization_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("normalization_alpha_map",)},
    "score_series": {"kind": "graphical", "renderer": "series", "required_roles": ("raw_L_series",)},
    "correction_series": {"kind": "graphical", "renderer": "series", "required_roles": ("correction_C_series", "correction_uncertainty")},
    "baseline_support": {"kind": "graphical", "renderer": "series", "required_roles": ("baseline_B_series",)},
    "pion_mm": {"kind": "graphical", "renderer": "overlay", "required_roles": ("baseline_pion_mm", "refined_pion_mm")},
    "mm_products": {"kind": "graphical", "renderer": "overlay", "required_roles": ("host_mm", "baseline_pion_mm", "refined_pion_mm")},
    "clean_mm": {"kind": "graphical", "renderer": "overlay", "required_roles": ("baseline_clean_mm", "refined_clean_mm")},
    "resonance_closure": {"kind": "graphical", "renderer": "closure", "required_roles": ("baseline_pion_mm", "refined_pion_mm", "pi_n_reference", "pi_sidis_reference", "pi_delta_reference")},
    "signal_closure": {"kind": "graphical", "renderer": "closure", "required_roles": ("refined_clean_mm", "KLambda_reference")},
    "summary": {"kind": "text", "renderer": "summary", "required_roles": ()},
}


_PART2B_RENDER_SERIAL = 0


def _part2b_name(stem):
    global _PART2B_RENDER_SERIAL
    _PART2B_RENDER_SERIAL += 1
    return "part2b_{}_{}".format(stem, _PART2B_RENDER_SERIAL)


def _part2b_page(page_id, spec_key, *, t_index=None, status="ready", unavailable_reason=None, **detail):
    spec = PART2B_PAGE_SPECS[spec_key]
    page = {
        "page_id": str(page_id), "spec_key": str(spec_key), "t_index": t_index,
        "page_kind": spec["kind"], "renderer": spec["renderer"],
        "required_roles": list(spec["required_roles"]), "status": str(status),
        "unavailable_reason": unavailable_reason,
    }
    page.update(detail)
    return page


def _part2b_cells(payload):
    return sorted((payload or {}).get("cells") or (), key=lambda cell: (int(cell.get("t_index", -1)), int(cell.get("delta_index", -1))))


def extract_pion_hgcer_refinement_tdelta_matrix(payload, field):
    """JSON-safe t-delta matrix that preserves missing values as None."""
    t_edges, delta_edges = list((payload or {}).get("t_edges") or ()), list((payload or {}).get("delta_edges") or ())
    matrix = [[None for _ in range(max(0, len(delta_edges) - 1))] for _ in range(max(0, len(t_edges) - 1))]
    for cell in _part2b_cells(payload):
        try:
            t_index, delta_index = int(cell["t_index"]), int(cell["delta_index"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (0 <= t_index < len(matrix) and 0 <= delta_index < len(matrix[t_index])):
            continue
        if field == "alpha":
            value = next((entry.get("alpha") for entry in (payload or {}).get("t_results") or () if int(entry.get("t_index", -1)) == t_index), None)
        elif field == "source":
            value = cell.get("source")
        elif field == "uncertainty":
            value = cell.get("sigma_C_total")
        else:
            value = cell.get(field)
        matrix[t_index][delta_index] = value if isinstance(value, str) else _finite(value)
    return {"field": str(field), "t_edges": t_edges, "delta_edges": delta_edges, "values": matrix}


def extract_pion_hgcer_refinement_t_series(payload, t_index, field):
    edges = list((payload or {}).get("delta_edges") or ())
    rows = []
    for cell in _part2b_cells(payload):
        if int(cell.get("t_index", -1)) != int(t_index):
            continue
        delta_index = int(cell.get("delta_index", -1))
        if not 0 <= delta_index < len(edges) - 1:
            continue
        value = cell.get(field) if field not in {"uncertainty", "source"} else (cell.get("sigma_C_total") if field == "uncertainty" else cell.get("source"))
        rows.append({
            "delta_index": delta_index,
            "delta": 0.5 * (float(edges[delta_index]) + float(edges[delta_index + 1])),
            "value": value if isinstance(value, str) else _finite(value),
            "uncertainty": _finite(cell.get("sigma_C_total")), "source": cell.get("source"),
        })
    return {"field": str(field), "t_index": int(t_index), "points": rows}


def _part2b_histograms(payload, t_index=None):
    source = (payload or {}).get("global_histograms") if t_index is None else (payload or {}).get("histograms")
    if t_index is None:
        return dict(source or {})
    return dict((source or {}).get(t_index) or (source or {}).get(str(t_index)) or {})


def _part2b_closure_inputs(renderer_inputs, t_index):
    by_t = (renderer_inputs or {}).get("by_t") or {}
    return by_t.get(t_index) or by_t.get(str(t_index)) or {}


def _part2b_t_result(payload, t_index):
    return next((entry for entry in (payload or {}).get("t_results") or () if int(entry.get("t_index", -1)) == int(t_index)), {})


def expected_pion_hgcer_refinement_page_manifest(payload, *, renderer_inputs=None):
    """Runtime-sized append-only Part-2B page manifest."""
    if not isinstance(payload, dict) or payload.get("status") != "available":
        return [_part2b_page("part2b_unavailable", "unavailable", status="unavailable", unavailable_reason=(payload or {}).get("reason", "unknown"))]
    pages = [
        _part2b_page("part2b_audit", "audit"),
        _part2b_page("part2b_raw_L_map", "raw_score_map"),
        _part2b_page("part2b_C_map", "correction_map"),
        _part2b_page("part2b_score_source_map", "source_map"),
        _part2b_page("part2b_C_uncertainty_map", "uncertainty_map"),
        _part2b_page("part2b_normalization_alpha_map", "normalization_map"),
    ]
    t_count = max(0, len(payload.get("t_edges") or ()) - 1)
    for t_index in range(t_count):
        number, hists = t_index + 1, _part2b_histograms(payload, t_index)
        for page_id, spec, field in (
            ("part2b_raw_L_t{}".format(number), "score_series", "L"),
            ("part2b_C_t{}".format(number), "correction_series", "C"),
            ("part2b_B_t{}".format(number), "baseline_support", "B"),
        ):
            if any(point.get("value") is not None for point in extract_pion_hgcer_refinement_t_series(payload, t_index, field)["points"]):
                pages.append(_part2b_page(page_id, spec, t_index=t_index))
            else:
                pages.append(_part2b_page(page_id + "_unavailable", "unavailable", t_index=t_index, status="unavailable", unavailable_reason="no_defined_{}_values".format(field)))
        for suffix, spec in (("pion_mm", "pion_mm"), ("mm_products", "mm_products"), ("clean_mm", "clean_mm")):
            if hists:
                pages.append(_part2b_page("part2b_{}_t{}".format(suffix, number), spec, t_index=t_index))
            else:
                pages.append(_part2b_page("part2b_{}_t{}_unavailable".format(suffix, number), "unavailable", t_index=t_index, status="unavailable", unavailable_reason="detached_noRF_allcuts_products_missing"))
        closure = _part2b_closure_inputs(renderer_inputs, t_index)
        simc = closure.get("simc") or {}
        metrics = (_part2b_t_result(payload, t_index).get("component_closure_metrics") or {})
        if hists and all(simc.get(key) is not None for key in ("pi-n", "pi-SIDIS", "pi-delta")) and metrics.get("status") == "available":
            pages.append(_part2b_page("part2b_resonance_closure_t{}".format(number), "resonance_closure", t_index=t_index))
        else:
            pages.append(_part2b_page("part2b_resonance_closure_t{}_unavailable".format(number), "unavailable", t_index=t_index, status="unavailable", unavailable_reason=(metrics.get("reason") or "pi_n_pi_sidis_or_pi_delta_closure_missing")))
        signal = closure.get("signal") or {}
        lambda_hist = (signal.get("K-Lambda") or {}).get("hist") if isinstance(signal.get("K-Lambda"), dict) else signal.get("K-Lambda")
        if hists and lambda_hist is not None:
            required = list(PART2B_PAGE_SPECS["signal_closure"]["required_roles"])
            sigma_hist = (signal.get("K-Sigma0") or {}).get("hist") if isinstance(signal.get("K-Sigma0"), dict) else signal.get("K-Sigma0")
            required.append("KSigma0_reference" if sigma_hist is not None else "KSigma0_unavailable_status")
            pages.append(_part2b_page("part2b_signal_closure_t{}".format(number), "signal_closure", t_index=t_index, required_roles=required))
        else:
            pages.append(_part2b_page("part2b_signal_closure_t{}_unavailable".format(number), "unavailable", t_index=t_index, status="unavailable", unavailable_reason="KLambda_closure_missing"))
    if _part2b_histograms(payload, None):
        pages.extend([
            _part2b_page("part2b_pion_mm_global", "pion_mm", scope="global", required_roles=list(PART2B_PAGE_SPECS["pion_mm"]["required_roles"]) + ["strict_global_sum"]),
            _part2b_page("part2b_mm_products_global", "mm_products", scope="global", required_roles=list(PART2B_PAGE_SPECS["mm_products"]["required_roles"]) + ["strict_global_sum"]),
            _part2b_page("part2b_clean_mm_global", "clean_mm", scope="global", required_roles=list(PART2B_PAGE_SPECS["clean_mm"]["required_roles"]) + ["strict_global_sum"]),
        ])
    else:
        pages.append(_part2b_page("part2b_global_products_unavailable", "unavailable", status="unavailable", unavailable_reason="strict_global_noRF_products_missing"))
    pages.append(_part2b_page("part2b_final_summary", "summary"))
    return pages


def _part2b_root_clone(histogram, role, color, style=1):
    clone = histogram.Clone(_part2b_name(role))
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    clone.SetLineColor(color); clone.SetLineStyle(style); clone.SetLineWidth(2); clone.SetFillStyle(0)
    return clone


def _part2b_render_text(canvas, page, payload, title_prefix):
    text = ROOT.TPaveText(0.08, 0.12, 0.92, 0.88, "NDC")
    text.SetFillStyle(0); text.SetBorderSize(1); text.SetTextAlign(12); text.SetTextSize(0.026)
    if str(title_prefix).strip():
        text.AddText(str(title_prefix).strip())
    text.AddText(ROOT_SAFE_REFINEMENT_LABELS["unavailable"] if page["status"] == "unavailable" else ROOT_SAFE_REFINEMENT_LABELS["title"])
    text.AddText(str(page["page_id"]).replace("_", " "))
    if page["spec_key"] == "audit":
        audit = (payload or {}).get("control_audit") or {}
        text.AddText("physical-control cache audit: {}".format(audit.get("status", "unavailable")))
        text.AddText("raw score fingerprint: {}".format(str((payload or {}).get("raw_hgcer_score_fingerprint") or "unavailable")[:16]))
        text.AddText("refinement fingerprint: {}".format(str((payload or {}).get("refinement_fingerprint") or "unavailable")[:16]))
        text.AddText("C=1 when Part-2B cannot establish a trustworthy relative score.")
    elif page["spec_key"] == "summary":
        text.AddText("Existing pion event normalization and MM physics remain unchanged.")
        text.AddText("The 1e-4 signed-support ratio is numerical safety only, not physics support.")
        for entry in (payload or {}).get("t_results") or ():
            text.AddText("t{}: status={} alpha={} f_N={} f_D={}".format(int(entry.get("t_index", -1)) + 1, entry.get("status"), entry.get("alpha"), entry.get("f_N"), entry.get("f_D")))
        text.AddText("Part 3 review is required; authoritative pion subtraction is unchanged.")
    else:
        text.AddText("Reason: {}".format(page.get("unavailable_reason") or (payload or {}).get("reason") or "unknown"))
    text.Draw()
    return [("audit_text" if page["spec_key"] == "audit" else "status_text", text)]


def _part2b_make_map(page, matrix, title):
    delta_edges, t_edges = matrix["delta_edges"], matrix["t_edges"]
    if len(delta_edges) < 2 or len(t_edges) < 2:
        raise RuntimeError("Part-2B map geometry unavailable")
    hist = ROOT.TH2D(_part2b_name(page["page_id"]), title, len(delta_edges) - 1, array("d", [float(value) for value in delta_edges]), len(t_edges) - 1, array("d", [float(value) for value in t_edges]))
    hist.SetDirectory(0); hist.GetXaxis().SetTitle("SHMS delta (%)"); hist.GetYaxis().SetTitle("|t| (GeV^2)")
    return hist


def _part2b_render_numeric_map(canvas, page, payload, _inputs):
    fields = {"raw_score_map": ("L", "raw_L_map", ROOT_SAFE_REFINEMENT_LABELS["raw_score"]), "correction_map": ("C", "correction_C_map", ROOT_SAFE_REFINEMENT_LABELS["correction"]), "uncertainty_map": ("uncertainty", "correction_uncertainty_map", ROOT_SAFE_REFINEMENT_LABELS["uncertainty"]), "normalization_map": ("alpha", "normalization_alpha_map", "Part 2B alpha(t), repeated within canonical t")}
    field, role, title = fields[page["spec_key"]]
    matrix = extract_pion_hgcer_refinement_tdelta_matrix(payload, field)
    hist, masks = _part2b_make_map(page, matrix, title), []
    for t_index, row in enumerate(matrix["values"]):
        for delta_index, value in enumerate(row):
            if value is None:
                box = ROOT.TBox(matrix["delta_edges"][delta_index], matrix["t_edges"][t_index], matrix["delta_edges"][delta_index + 1], matrix["t_edges"][t_index + 1])
                box.SetFillColor(ROOT.kGray + 1); box.SetFillStyle(3004); box.SetLineColor(ROOT.kGray + 2); masks.append(box)
            else:
                hist.SetBinContent(delta_index + 1, t_index + 1, float(value))
    hist.Draw("COLZ TEXT")
    for box in masks:
        box.Draw("same")
    note = ROOT.TLatex(); note.SetNDC(True); note.SetTextSize(0.025); note.DrawLatex(0.12, 0.02, "Hatched cells are undefined, never zero-filled.")
    return [(role, hist), ("undefined_cell_note", note)] + [("undefined_cell_mask", box) for box in masks]


def _part2b_render_source_map(canvas, page, payload, _inputs):
    matrix = extract_pion_hgcer_refinement_tdelta_matrix(payload, "source")
    categories = sorted({value for row in matrix["values"] for value in row if value is not None})
    codes = {value: index + 1 for index, value in enumerate(categories)}
    hist = _part2b_make_map(page, matrix, ROOT_SAFE_REFINEMENT_LABELS["source"])
    for t_index, row in enumerate(matrix["values"]):
        for delta_index, value in enumerate(row):
            if value is not None:
                hist.SetBinContent(delta_index + 1, t_index + 1, codes[value])
    hist.Draw("COLZ TEXT")
    legend = ROOT.TPaveText(0.66, 0.60, 0.97, 0.92, "NDC")
    legend.SetFillStyle(0); legend.SetBorderSize(1); legend.SetTextSize(0.018); legend.AddText("score source codes")
    for value in categories:
        legend.AddText("{} = {}".format(codes[value], value))
    legend.Draw()
    return [("score_source_map", hist), ("source_legend", legend)]


def _part2b_render_series(canvas, page, payload, _inputs):
    field, role = {"score_series": ("L", "raw_L_series"), "correction_series": ("C", "correction_C_series"), "baseline_support": ("B", "baseline_B_series")}[page["spec_key"]]
    series = extract_pion_hgcer_refinement_t_series(payload, page["t_index"], field)
    rows = [row for row in series["points"] if row.get("value") is not None]
    if not rows:
        raise RuntimeError("Part-2B series has no defined value")
    graph = ROOT.TGraphErrors(len(rows), array("d", [float(row["delta"]) for row in rows]), array("d", [float(row["value"]) for row in rows]), array("d", [0.0] * len(rows)), array("d", [float(row.get("uncertainty") or 0.0) if field == "C" else 0.0 for row in rows]))
    graph.SetName(_part2b_name(page["page_id"])); graph.SetMarkerStyle(20); graph.SetLineWidth(2)
    graph.SetTitle("Part 2B {} - canonical t{};SHMS delta (%);{}".format(field, int(page["t_index"]) + 1, field)); graph.Draw("APZ")
    primitives = [(role, graph)]
    if field == "C":
        primitives.append(("correction_uncertainty", graph))
        unresolved = [row for row in series["points"] if row.get("value") is None]
        note = ROOT.TLatex(); note.SetNDC(True); note.SetTextSize(0.025); note.DrawLatex(0.13, 0.15, "Errors are undefined for identity/not-measured cells; they are not zero.")
        primitives.append(("uncertainty_note", note))
        if unresolved:
            primitives.append(("unresolved_cells", note))
    return primitives


def _part2b_render_overlay(canvas, page, payload, _inputs):
    hists = _part2b_histograms(payload, None if page.get("scope") == "global" else page["t_index"])
    sources = {
        "pion_mm": [("baseline_pion_mm", hists.get("baseline_pion"), "baseline pion", ROOT.kBlack), ("refined_pion_mm", hists.get("refined_pion"), "refined pion", ROOT.kRed + 1)],
        "mm_products": [("host_mm", hists.get("host"), "proton-cleaned noRF host", ROOT.kBlack), ("baseline_pion_mm", hists.get("baseline_pion"), "baseline pion", ROOT.kOrange + 7), ("refined_pion_mm", hists.get("refined_pion"), "refined pion", ROOT.kRed + 1)],
        "clean_mm": [("baseline_clean_mm", hists.get("baseline_clean"), "baseline clean", ROOT.kOrange + 7), ("refined_clean_mm", hists.get("refined_clean"), "refined clean", ROOT.kBlue + 1)],
    }[page["spec_key"]]
    clones = [(role, _part2b_root_clone(hist, role, color), label) for role, hist, label, color in sources if hist is not None]
    if len(clones) != len(sources):
        raise RuntimeError("Part-2B required detached MM primitive missing")
    first = clones[0][1]
    minimum, maximum = min(float(item[1].GetMinimum()) for item in clones), max(float(item[1].GetMaximum()) for item in clones)
    first.SetTitle("Part 2B detached allcuts noRF MM - {};MM (GeV);weighted counts".format("global strict sum" if page.get("scope") == "global" else "canonical t{}".format(int(page["t_index"]) + 1)))
    first.SetMinimum(minimum * 1.10 if minimum < 0.0 else 0.0); first.SetMaximum(maximum * 1.18 if maximum > 0.0 else 1.0); first.Draw("HIST")
    for _, clone, _ in clones[1:]: clone.Draw("HIST SAME")
    legend = ROOT.TLegend(0.56, 0.68, 0.89, 0.89); legend.SetBorderSize(0)
    for _, clone, label in clones: legend.AddEntry(clone, label, "l")
    legend.Draw()
    primitives = [(role, clone) for role, clone, _ in clones] + [("overlay_legend", legend)]
    if page.get("scope") == "global":
        strict = ROOT.TLatex(); strict.SetNDC(True); strict.SetTextSize(0.025)
        strict.DrawLatex(0.12, 0.02, "Global histogram is the strict sum of canonical t products.")
        primitives.append(("strict_global_sum", strict))
    return primitives


def _part2b_render_closure(canvas, page, payload, renderer_inputs):
    hists, closure = _part2b_histograms(payload, page["t_index"]), _part2b_closure_inputs(renderer_inputs, page["t_index"])
    if page["spec_key"] == "resonance_closure":
        simc = closure.get("simc") or {}
        sources = [("baseline_pion_mm", hists.get("baseline_pion"), "baseline pion", ROOT.kBlack), ("refined_pion_mm", hists.get("refined_pion"), "refined pion", ROOT.kRed + 1), ("pi_n_reference", simc.get("pi-n"), "pi-n", ROOT.kBlue + 1), ("pi_sidis_reference", simc.get("pi-SIDIS"), "pi-SIDIS", ROOT.kGreen + 2), ("pi_delta_reference", simc.get("pi-delta"), "pi-delta", ROOT.kMagenta + 2)]
    else:
        signal = closure.get("signal") or {}; lamb = signal.get("K-Lambda") or {}; sigma = signal.get("K-Sigma0") or {}
        lamb = lamb.get("hist") if isinstance(lamb, dict) else lamb; sigma = sigma.get("hist") if isinstance(sigma, dict) else sigma
        sources = [("refined_clean_mm", hists.get("refined_clean"), "refined clean", ROOT.kBlack), ("KLambda_reference", lamb, "K-Lambda display", ROOT.kBlue + 1)]
        if sigma is not None: sources.append(("KSigma0_reference", sigma, "K-Sigma0", ROOT.kGreen + 2))
    clones = [(role, _part2b_root_clone(hist, role, color), label) for role, hist, label, color in sources if hist is not None]
    if len(clones) < len(PART2B_PAGE_SPECS[page["spec_key"]]["required_roles"]):
        raise RuntimeError("Part-2B closure input missing")
    first = clones[0][1]; first.SetTitle("Part 2B closure - canonical t{};MM (GeV);weighted counts".format(int(page["t_index"]) + 1)); first.Draw("HIST")
    for _, clone, _ in clones[1:]: clone.Draw("HIST SAME")
    legend = ROOT.TLegend(0.56, 0.64, 0.89, 0.89); legend.SetBorderSize(0)
    for _, clone, label in clones: legend.AddEntry(clone, label, "l")
    legend.Draw(); primitives = [(role, clone) for role, clone, _ in clones] + [("closure_legend", legend)]
    if page["spec_key"] == "resonance_closure":
        summary = ROOT.TPaveText(0.12, 0.62, 0.53, 0.89, "NDC")
        summary.SetFillStyle(0); summary.SetBorderSize(1); summary.SetTextSize(0.018); summary.SetTextAlign(12)
        for label, metric in ((_part2b_t_result(payload, page["t_index"]).get("component_closure_metrics") or {}).get("metrics") or {}).items():
            summary.AddText("{}: baseline={} refined={} ref={}".format(label, metric.get("baseline"), metric.get("refined"), metric.get("reference")))
        summary.Draw(); primitives.append(("closure_window_metrics", summary))
    if page["spec_key"] == "signal_closure" and "KSigma0_unavailable_status" in page.get("required_roles", ()): primitives.append(("KSigma0_unavailable_status", legend))
    return primitives


def _part2b_render_graphical(canvas, page, payload, renderer_inputs):
    renderer = {"numeric_map": _part2b_render_numeric_map, "categorical_map": _part2b_render_source_map, "series": _part2b_render_series, "overlay": _part2b_render_overlay, "closure": _part2b_render_closure}[page["renderer"]]
    primitives = renderer(canvas, page, payload, renderer_inputs)
    roles = {role for role, _object in primitives}
    missing = [role for role in page.get("required_roles") or () if role not in roles]
    if missing:
        raise RuntimeError("Part-2B graphical page {} did not render semantic primitives {}".format(page["page_id"], ", ".join(missing)))
    return primitives


def render_pion_hgcer_refinement_pages(pdf_path, payload, *, title_prefix="", page_manifest=None, close_pdf=True, renderer_inputs=None):
    """Append Part-2B pages; frozen Part-2 payload and ROOT inputs are read-only."""
    if ROOT is None:
        return []
    manifest, pages = page_manifest if page_manifest is not None else [], expected_pion_hgcer_refinement_page_manifest(payload, renderer_inputs=renderer_inputs)
    emitted = []
    for page in pages:
        canvas = ROOT.TCanvas(_part2b_name(page["page_id"]), page["page_id"], 1200, 850)
        primitives = _part2b_render_text(canvas, page, payload, title_prefix) if page["page_kind"] == "text" else _part2b_render_graphical(canvas, page, payload, renderer_inputs)
        canvas.Print(pdf_path)
        entry = dict(page)
        entry["semantic_primitives"] = [{"role": role, "type": type(object).__name__} for role, object in primitives]
        manifest.append(_json_ready(entry)); emitted.append(_json_ready(entry))
    if close_pdf:
        ROOT.TCanvas(_part2b_name("close"), "Part2B close", 1, 1).Print(str(pdf_path) + ")")
    return emitted


def render_pion_hgcer_refinement_failure_page(pdf_path, reason, *, title_prefix="", page_manifest=None, close_pdf=True):
    payload = _unavailable_refinement(reason, {})
    return render_pion_hgcer_refinement_pages(pdf_path, payload, title_prefix=title_prefix, page_manifest=page_manifest, close_pdf=close_pdf)
