#! /usr/bin/python

from __future__ import annotations

import math
import os
import sys
from array import array
from copy import deepcopy

import numpy as np
import ROOT
from ROOT import (
    TCanvas,
    TLegend,
    TLine,
    TPaveText,
    gPad,
    gStyle,
    kBlack,
    kBlue,
    kGray,
    kGreen,
    kOrange,
    kRed,
    kViolet,
)
from scipy.optimize import least_squares, lsq_linear

sys.path.append("utility")
from background_config import (  # noqa: E402
    PROTON_CONTAMINATION_CLEANING_METHOD_CTIME_AERO_EVENT_WEIGHT,
    PROTON_CONTAMINATION_CLEANING_METHOD_DISABLED,
    get_proton_contamination_cleaning_config,
    resolve_proton_contamination_cleaning_enabled,
    resolve_proton_contamination_cleaning_method,
    resolve_proton_contamination_cleaning_rf_policy,
    resolve_proton_contamination_cleaning_tree_policy,
)
from prompt_trees import normalize_epsset  # noqa: E402


SUPPORT_UNSUPPORTED = "unsupported"
SUPPORT_MARGINAL = "marginal"
SUPPORT_SUPPORTED = "supported"
SUPPORT_CLASS_ORDER = (
    SUPPORT_UNSUPPORTED,
    SUPPORT_MARGINAL,
    SUPPORT_SUPPORTED,
)
SUPPORT_CLASS_TO_CODE = {
    SUPPORT_UNSUPPORTED: 0.0,
    SUPPORT_MARGINAL: 1.0,
    SUPPORT_SUPPORTED: 2.0,
}

DEFAULT_RF_BRANCH_CANDIDATES = (
    "RF",
    "RFTime",
    "RF_time",
    "RFTime_ROC1",
    "P_RFTime",
    "P_RF_tdcTime",
    "P_RF_adcTime",
    "RF_Dist",
    "RF_Distance",
    "P_RF_Dist",
    "P_RF_Distance",
    "P_RF_Dist_Track",
    "H_RFTime",
    "H_RF_tdcTime",
)


def _clone_hist(template_hist, name, reset=True):
    if template_hist is None:
        return None
    cloned = template_hist.Clone(name)
    if hasattr(cloned, "SetDirectory"):
        cloned.SetDirectory(0)
    if hasattr(cloned, "Sumw2"):
        cloned.Sumw2()
    if reset:
        cloned.Reset()
    return cloned


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def _hist_abs_integral(hist):
    if hist is None:
        return 0.0
    total = 0.0
    for bin_index in range(1, hist.GetNbinsX() + 1):
        total += abs(float(hist.GetBinContent(bin_index)))
    return float(total)


def _is_root_object(obj):
    try:
        return bool(obj is not None and obj.InheritsFrom("TObject"))
    except Exception:
        return False


_JSON_SKIP = object()


def _json_ready_value(value):
    if _is_root_object(value):
        return _JSON_SKIP
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        cleaned = {}
        for key, child_value in value.items():
            child = _json_ready_value(child_value)
            if child is _JSON_SKIP:
                continue
            cleaned[str(key)] = child
        return cleaned
    if isinstance(value, (list, tuple)):
        cleaned = []
        for child_value in value:
            child = _json_ready_value(child_value)
            if child is _JSON_SKIP:
                continue
            cleaned.append(child)
        return cleaned
    if isinstance(value, set):
        return sorted(_json_ready_value(list(value)))
    return deepcopy(value)


def _make_signature(evt, fields, round_digits):
    signature = []
    for field_name in fields:
        value = getattr(evt, field_name, None)
        if value is None:
            signature.append(None)
            continue
        try:
            signature.append(round(float(value), int(round_digits)))
        except Exception:
            signature.append(str(value))
    return tuple(signature)


def _tree_has_branch(tree, branch_name):
    if tree is None or not branch_name:
        return False
    try:
        return bool(tree.GetBranch(str(branch_name)))
    except Exception:
        return False


def _bundle_has_branch(source_bundle, branch_name):
    for source_spec in ((source_bundle or {}).get("sources") or {}).values():
        tree = (source_spec or {}).get("tree")
        if _tree_has_branch(tree, branch_name):
            return True
    return False


def _resolve_rf_branch_candidates(source_bundle):
    if str(os.environ.get("PROTON_CHECKER_DISABLE_RF", "")).strip() == "1":
        return []
    candidates = []
    forced = str(os.environ.get("PROTON_CHECKER_RF_BRANCH", "")).strip()
    if forced:
        candidates.append(forced)
    for candidate in DEFAULT_RF_BRANCH_CANDIDATES:
        if candidate not in candidates:
            candidates.append(candidate)
    return [candidate for candidate in candidates if _bundle_has_branch(source_bundle, candidate)]


def _preselection_passes(evt, evaluate_event, hole_contains, mm_min, mm_max):
    _, nommcuts, _ = evaluate_event(evt, mm_min, mm_max)
    hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) if hole_contains is not None else False
    return bool(nommcuts and (not hole_rejected))


def _collect_branch_values(source_bundle, evaluate_event, hole_contains, mm_min, mm_max, branch_name):
    values = []
    for source_spec in ((source_bundle or {}).get("sources") or {}).values():
        tree = (source_spec or {}).get("tree")
        if tree is None or not _tree_has_branch(tree, branch_name):
            continue
        for evt in tree:
            if not _preselection_passes(evt, evaluate_event, hole_contains, mm_min, mm_max):
                continue
            try:
                value = float(getattr(evt, branch_name))
            except Exception:
                continue
            if math.isfinite(value):
                values.append(value)
    return values


def _build_rf_membership_lookup(rf_source_bundle, signature_fields, round_digits):
    lookup = {}
    counts = {}
    for bundle_key in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
        tree = (rf_source_bundle or {}).get(bundle_key)
        accepted = set()
        if tree is not None:
            for evt in tree:
                accepted.add(_make_signature(evt, signature_fields, round_digits))
        lookup[bundle_key] = accepted
        counts[bundle_key] = len(accepted)
    return lookup, counts


def apply_low_epsilon_rf_after_proton_cleaning(cleaning_result, source_label, evt):
    if not isinstance(cleaning_result, dict):
        return True
    diagnostics = cleaning_result.get("diagnostics") or {}
    if not bool(diagnostics.get("rf_applied", False)):
        return True
    lookup = cleaning_result.get("_rf_signature_lookup") or {}
    fields = diagnostics.get("rf_signature_fields") or ()
    round_digits = int(diagnostics.get("signature_round_digits", 9) or 9)
    accepted = lookup.get(str(source_label), set())
    if not accepted:
        return False
    return _make_signature(evt, fields, round_digits) in accepted


def _find_peak_seed(histogram, x_min, x_max):
    if histogram is None:
        return 0.0, 0.5 * (float(x_min) + float(x_max))
    first_bin = max(1, histogram.GetXaxis().FindFixBin(float(x_min)))
    last_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(x_max), float(x_min))),
    )
    maximum = 0.0
    maximum_center = 0.5 * (float(x_min) + float(x_max))
    for bin_index in range(first_bin, last_bin + 1):
        content = float(histogram.GetBinContent(bin_index))
        if content > maximum:
            maximum = content
            maximum_center = float(histogram.GetXaxis().GetBinCenter(bin_index))
    return maximum, maximum_center


def _find_separated_peak_pair(histogram, x_min, x_max, minimum_separation, maximum_separation):
    result = {
        "valid": False,
        "lower_mean": 0.0,
        "upper_mean": 0.0,
        "lower_height": 0.0,
        "upper_height": 0.0,
    }
    if (
        histogram is None
        or float(x_max) <= float(x_min)
        or float(minimum_separation) <= 0.0
        or float(maximum_separation) <= float(minimum_separation)
    ):
        return result
    first_bin = max(1, histogram.GetXaxis().FindFixBin(float(x_min)))
    last_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(x_max), float(x_min))),
    )
    if (last_bin - first_bin) < 4:
        return result
    smoothed = np.zeros(histogram.GetNbinsX() + 1, dtype=float)
    smoothing_radius = 2
    for bin_index in range(first_bin, last_bin + 1):
        weighted_sum = 0.0
        weight_sum = 0.0
        for offset in range(-smoothing_radius, smoothing_radius + 1):
            neighbor = bin_index + offset
            if neighbor < first_bin or neighbor > last_bin:
                continue
            weight = float((smoothing_radius + 1) - abs(offset))
            weighted_sum += weight * max(float(histogram.GetBinContent(neighbor)), 0.0)
            weight_sum += weight
        if weight_sum > 0.0:
            smoothed[bin_index] = weighted_sum / weight_sum
    peak_bins = []
    for bin_index in range(first_bin + 1, last_bin):
        value = smoothed[bin_index]
        if value > 0.0 and value >= smoothed[bin_index - 1] and value > smoothed[bin_index + 1]:
            peak_bins.append(bin_index)
    if len(peak_bins) < 2:
        all_bins = list(range(first_bin, last_bin + 1))
        all_bins.sort(key=lambda index: smoothed[index], reverse=True)
        for bin_index in all_bins[: min(len(all_bins), 24)]:
            if bin_index not in peak_bins:
                peak_bins.append(bin_index)
    best_score = -1.0
    for first_index in range(len(peak_bins)):
        for second_index in range(first_index + 1, len(peak_bins)):
            first_bin_index = int(peak_bins[first_index])
            second_bin_index = int(peak_bins[second_index])
            first_mean = float(histogram.GetXaxis().GetBinCenter(first_bin_index))
            second_mean = float(histogram.GetXaxis().GetBinCenter(second_bin_index))
            separation = abs(second_mean - first_mean)
            if separation < float(minimum_separation) or separation > float(maximum_separation):
                continue
            first_height = float(smoothed[first_bin_index])
            second_height = float(smoothed[second_bin_index])
            if first_height <= 0.0 or second_height <= 0.0:
                continue
            score = min(first_height, second_height) * math.sqrt(first_height * second_height)
            if score <= best_score:
                continue
            best_score = score
            result["valid"] = True
            if first_mean < second_mean:
                result["lower_mean"] = first_mean
                result["upper_mean"] = second_mean
                result["lower_height"] = first_height
                result["upper_height"] = second_height
            else:
                result["lower_mean"] = second_mean
                result["upper_mean"] = first_mean
                result["lower_height"] = second_height
                result["upper_height"] = first_height
    return result


def _find_prominent_offset_peak_pair(histogram, x_min, x_max, minimum_separation, maximum_separation):
    result = {
        "valid": False,
        "lower_mean": 0.0,
        "upper_mean": 0.0,
        "lower_height": 0.0,
        "upper_height": 0.0,
    }
    if (
        histogram is None
        or float(x_max) <= float(x_min)
        or float(minimum_separation) <= 0.0
        or float(maximum_separation) <= float(minimum_separation)
    ):
        return result
    first_bin = max(1, histogram.GetXaxis().FindFixBin(float(x_min)))
    last_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(x_max), float(x_min))),
    )
    if (last_bin - first_bin) < 6:
        return result
    smoothed = np.zeros(histogram.GetNbinsX() + 1, dtype=float)
    smoothing_radius = 2
    global_maximum = 0.0
    for bin_index in range(first_bin, last_bin + 1):
        weighted_sum = 0.0
        weight_sum = 0.0
        for offset in range(-smoothing_radius, smoothing_radius + 1):
            neighbor = bin_index + offset
            if neighbor < first_bin or neighbor > last_bin:
                continue
            weight = float((smoothing_radius + 1) - abs(offset))
            weighted_sum += weight * max(float(histogram.GetBinContent(neighbor)), 0.0)
            weight_sum += weight
        if weight_sum > 0.0:
            smoothed[bin_index] = weighted_sum / weight_sum
            global_maximum = max(global_maximum, smoothed[bin_index])
    if not (global_maximum > 0.0):
        return result
    candidate_bins = []
    for bin_index in range(first_bin + 1, last_bin):
        value = smoothed[bin_index]
        if value > 0.0 and value >= smoothed[bin_index - 1] and value > smoothed[bin_index + 1]:
            candidate_bins.append(bin_index)
    ranked_bins = list(range(first_bin, last_bin + 1))
    ranked_bins.sort(key=lambda index: smoothed[index], reverse=True)
    exclusion_bins = 2
    for candidate in ranked_bins[: min(len(ranked_bins), 32)]:
        if all(abs(int(candidate) - int(existing)) > exclusion_bins for existing in candidate_bins):
            candidate_bins.append(int(candidate))
    best_score = -1.0
    for first_index in range(len(candidate_bins)):
        for second_index in range(first_index + 1, len(candidate_bins)):
            lower_bin = int(candidate_bins[first_index])
            upper_bin = int(candidate_bins[second_index])
            if lower_bin > upper_bin:
                lower_bin, upper_bin = upper_bin, lower_bin
            lower_mean = float(histogram.GetXaxis().GetBinCenter(lower_bin))
            upper_mean = float(histogram.GetXaxis().GetBinCenter(upper_bin))
            separation = upper_mean - lower_mean
            if (
                separation < float(minimum_separation)
                or separation > float(maximum_separation)
                or (upper_bin - lower_bin) < 3
            ):
                continue
            lower_height = float(smoothed[lower_bin])
            upper_height = float(smoothed[upper_bin])
            if (
                lower_height < (0.025 * global_maximum)
                or upper_height < (0.025 * global_maximum)
            ):
                continue
            valley = float("inf")
            for bin_index in range(lower_bin + 1, upper_bin):
                valley = min(valley, float(smoothed[bin_index]))
            if not math.isfinite(valley):
                continue
            lower_prominence = lower_height - valley
            upper_prominence = upper_height - valley
            minimum_prominence = max(0.012 * global_maximum, 0.025 * min(lower_height, upper_height))
            if lower_prominence < minimum_prominence or upper_prominence < minimum_prominence:
                continue
            lower_fraction = lower_prominence / max(lower_height, 1.0e-12)
            upper_fraction = upper_prominence / max(upper_height, 1.0e-12)
            score = (
                math.sqrt(lower_height * upper_height)
                * math.sqrt(lower_fraction * upper_fraction)
                * math.sqrt(max(lower_prominence * upper_prominence, 0.0))
            )
            if score <= best_score:
                continue
            best_score = score
            result["valid"] = True
            result["lower_mean"] = lower_mean
            result["upper_mean"] = upper_mean
            result["lower_height"] = lower_height
            result["upper_height"] = upper_height
    return result


def _find_central_quantile_range(histogram, requested_min, requested_max, lower_tail_fraction, upper_tail_fraction, padding_bins):
    if (
        histogram is None
        or float(requested_max) <= float(requested_min)
        or float(lower_tail_fraction) < 0.0
        or float(upper_tail_fraction) < 0.0
        or (float(lower_tail_fraction) + float(upper_tail_fraction)) >= 1.0
    ):
        return float(requested_min), float(requested_max)
    first_requested_bin = max(1, histogram.GetXaxis().FindFixBin(float(requested_min)))
    last_requested_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(requested_max), float(requested_min))),
    )
    total = 0.0
    for bin_index in range(first_requested_bin, last_requested_bin + 1):
        total += max(float(histogram.GetBinContent(bin_index)), 0.0)
    if total <= 0.0:
        return float(requested_min), float(requested_max)
    lower_target = float(lower_tail_fraction) * total
    upper_target = (1.0 - float(upper_tail_fraction)) * total
    cumulative = 0.0
    first_central_bin = first_requested_bin
    last_central_bin = last_requested_bin
    found_lower = False
    for bin_index in range(first_requested_bin, last_requested_bin + 1):
        cumulative += max(float(histogram.GetBinContent(bin_index)), 0.0)
        if (not found_lower) and cumulative >= lower_target:
            first_central_bin = bin_index
            found_lower = True
        if cumulative >= upper_target:
            last_central_bin = bin_index
            break
    first_central_bin = max(first_requested_bin, first_central_bin - max(int(padding_bins), 0))
    last_central_bin = min(last_requested_bin, last_central_bin + max(int(padding_bins), 0))
    range_min = max(float(requested_min), float(histogram.GetXaxis().GetBinLowEdge(first_central_bin)))
    range_max = min(float(requested_max), float(histogram.GetXaxis().GetBinUpEdge(last_central_bin)))
    if range_max <= range_min:
        return float(requested_min), float(requested_max)
    return range_min, range_max


def _estimate_value_central_range(values, fallback_min, fallback_max):
    finite_values = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if len(finite_values) < 25 or float(fallback_max) <= float(fallback_min):
        return float(fallback_min), float(fallback_max)
    finite_values.sort()
    lower = float(np.quantile(finite_values, 0.001))
    upper = float(np.quantile(finite_values, 0.999))
    if (not math.isfinite(lower)) or (not math.isfinite(upper)) or upper <= lower:
        return float(fallback_min), float(fallback_max)
    width = upper - lower
    padding = max(0.05 * width, 1.0e-6)
    return lower - padding, upper + padding


def _extract_weighted_hist_fit_inputs(histogram, fit_min, fit_max):
    x_values = []
    y_values = []
    sigma_values = []
    excluded_invalid_variance_bins = []
    for bin_index in range(1, histogram.GetNbinsX() + 1):
        x_value = float(histogram.GetXaxis().GetBinCenter(bin_index))
        if x_value < float(fit_min) or x_value > float(fit_max):
            continue
        y_value = float(histogram.GetBinContent(bin_index))
        variance = float(histogram.GetBinError(bin_index) ** 2)
        if (not math.isfinite(variance)) or variance <= 0.0:
            excluded_invalid_variance_bins.append(int(bin_index))
            continue
        x_values.append(x_value)
        y_values.append(y_value)
        sigma_values.append(math.sqrt(variance))
    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "sigma": np.asarray(sigma_values, dtype=float),
        "excluded_invalid_variance_bins": excluded_invalid_variance_bins,
    }


def _gaussian(x_values, amplitude, mean, sigma):
    safe_sigma = max(float(sigma), 1.0e-12)
    z_values = (np.asarray(x_values, dtype=float) - float(mean)) / safe_sigma
    return float(amplitude) * np.exp(-0.5 * np.square(z_values))


def _compute_parameter_covariance(weighted_design, parameter_names):
    covariance_matrix = {}
    correlation_matrix = {}
    uncertainties = {}
    if weighted_design is None or weighted_design.size == 0:
        return covariance_matrix, correlation_matrix, uncertainties
    try:
        normal_matrix = np.matmul(weighted_design.T, weighted_design)
        covariance = np.linalg.pinv(normal_matrix)
    except Exception:
        return covariance_matrix, correlation_matrix, uncertainties
    for i, left_name in enumerate(parameter_names):
        left_row = {}
        corr_row = {}
        left_var = float(covariance[i, i])
        uncertainties[left_name] = math.sqrt(max(left_var, 0.0))
        for j, right_name in enumerate(parameter_names):
            covariance_value = float(covariance[i, j])
            left_row[right_name] = covariance_value
            right_var = float(covariance[j, j])
            if left_var > 0.0 and right_var > 0.0:
                corr_row[right_name] = covariance_value / math.sqrt(left_var * right_var)
            else:
                corr_row[right_name] = None
        covariance_matrix[left_name] = left_row
        correlation_matrix[left_name] = corr_row
    return covariance_matrix, correlation_matrix, uncertainties


def _fit_global_timing_shape_with_bounds(
    histogram,
    function_name,
    fit_min,
    fit_max,
    kaon_mean_min,
    kaon_mean_max,
    proton_mean_min,
    proton_mean_max,
    proton_peak_is_lower,
    sigma_min,
    sigma_max,
    initial_sigma,
    minimum_separation,
    minimum_amplitude_significance,
    maximum_chi2_ndf,
    bound_fraction_tolerance,
    minimum_entries,
):
    inputs = _extract_weighted_hist_fit_inputs(histogram, fit_min, fit_max)
    x_values = inputs["x"]
    y_values = inputs["y"]
    sigma_values = inputs["sigma"]
    if x_values.size == 0 or np.sum(np.abs(y_values)) < float(minimum_entries):
        return {
            "valid": False,
            "fit_status": "insufficient_support",
            "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
            "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
            "function_name": str(function_name),
            "fit_min": float(fit_min),
            "fit_max": float(fit_max),
            "per_aero_fallback": False,
        }
    histogram_maximum = max(float(np.max(y_values)), 1.0)
    kaon_seed = _find_peak_seed(histogram, kaon_mean_min, kaon_mean_max)
    proton_seed = _find_peak_seed(histogram, proton_mean_min, proton_mean_max)
    initial = np.asarray(
        [
            max(float(kaon_seed[0]), 0.20 * histogram_maximum),
            float(np.clip(kaon_seed[1], float(kaon_mean_min), float(kaon_mean_max))),
            float(initial_sigma),
            max(float(proton_seed[0]), 0.10 * histogram_maximum),
            float(np.clip(proton_seed[1], float(proton_mean_min), float(proton_mean_max))),
            float(initial_sigma),
            0.02 * histogram_maximum,
        ],
        dtype=float,
    )
    lower_bounds = np.asarray(
        [0.0, float(kaon_mean_min), float(sigma_min), 0.0, float(proton_mean_min), float(sigma_min), 0.0],
        dtype=float,
    )
    upper_bounds = np.asarray(
        [
            100.0 * histogram_maximum,
            float(kaon_mean_max),
            float(sigma_max),
            100.0 * histogram_maximum,
            float(proton_mean_max),
            float(sigma_max),
            10.0 * histogram_maximum,
        ],
        dtype=float,
    )

    def residuals(parameters):
        kaon_model = _gaussian(x_values, parameters[0], parameters[1], parameters[2])
        proton_model = _gaussian(x_values, parameters[3], parameters[4], parameters[5])
        return (y_values - (kaon_model + proton_model + float(parameters[6]))) / sigma_values

    fit_result = least_squares(
        residuals,
        initial,
        bounds=(lower_bounds, upper_bounds),
        method="trf",
        loss="linear",
    )
    parameter_vector = np.asarray(fit_result.x, dtype=float)
    jacobian = getattr(fit_result, "jac", None)
    weighted_design = np.asarray(jacobian, dtype=float) if jacobian is not None else None
    covariance_matrix, correlation_matrix, uncertainties = _compute_parameter_covariance(
        weighted_design,
        (
            "kaon_amplitude",
            "kaon_mean",
            "kaon_sigma",
            "proton_amplitude",
            "proton_mean",
            "proton_sigma",
            "other_amplitude",
        ),
    )
    residual_vector = residuals(parameter_vector)
    active_bin_count = int(x_values.size)
    ndf = max(active_bin_count - 7, 1)
    chi2_data = float(np.sum(np.square(residual_vector)))
    chi2_ndf = float(chi2_data / float(ndf)) if ndf > 0 else float("inf")
    abs_entry_sum = float(np.sum(np.abs(y_values)))
    chi2_per_abs_entry = float(chi2_data / abs_entry_sum) if abs_entry_sum > 0.0 else float("inf")
    kaon_amp_err = float(uncertainties.get("kaon_amplitude", 0.0))
    proton_amp_err = float(uncertainties.get("proton_amplitude", 0.0))
    kaon_sigma = abs(float(parameter_vector[2]))
    proton_sigma = abs(float(parameter_vector[5]))
    separation_denominator = math.sqrt(max((kaon_sigma ** 2) + (proton_sigma ** 2), 0.0))
    separation = (
        abs(float(parameter_vector[4]) - float(parameter_vector[1])) / separation_denominator
        if separation_denominator > 0.0
        else 0.0
    )

    def near_bound(value, lower, upper):
        tolerance = float(bound_fraction_tolerance) * float(upper - lower)
        return (value <= float(lower) + tolerance) or (value >= float(upper) - tolerance)

    bound_hit = (
        near_bound(parameter_vector[1], kaon_mean_min, kaon_mean_max)
        or near_bound(parameter_vector[4], proton_mean_min, proton_mean_max)
        or near_bound(kaon_sigma, sigma_min, sigma_max)
        or near_bound(proton_sigma, sigma_min, sigma_max)
    )
    ordering_valid = (
        float(parameter_vector[4]) < float(parameter_vector[1])
        if bool(proton_peak_is_lower)
        else float(parameter_vector[1]) < float(parameter_vector[4])
    )
    kaon_significance = float(parameter_vector[0] / kaon_amp_err) if kaon_amp_err > 0.0 else 0.0
    proton_significance = float(parameter_vector[3] / proton_amp_err) if proton_amp_err > 0.0 else 0.0
    valid = (
        bool(getattr(fit_result, "success", False))
        and not bound_hit
        and ordering_valid
        and math.isfinite(chi2_ndf)
        and float(parameter_vector[0]) > 0.0
        and float(parameter_vector[3]) > 0.0
        and kaon_sigma > 0.0
        and proton_sigma > 0.0
        and separation >= float(minimum_separation)
        and kaon_significance >= float(minimum_amplitude_significance)
        and proton_significance >= float(minimum_amplitude_significance)
        and chi2_ndf <= float(maximum_chi2_ndf)
    )
    return {
        "valid": bool(valid),
        "fit_status": "success" if bool(getattr(fit_result, "success", False)) else "failure",
        "message": str(getattr(fit_result, "message", "")),
        "function_name": str(function_name),
        "kaon_amplitude": float(parameter_vector[0]),
        "kaon_amplitude_error": kaon_amp_err,
        "kaon_mean": float(parameter_vector[1]),
        "kaon_sigma": float(kaon_sigma),
        "proton_amplitude": float(parameter_vector[3]),
        "proton_amplitude_error": proton_amp_err,
        "proton_mean": float(parameter_vector[4]),
        "proton_sigma": float(proton_sigma),
        "other_amplitude": float(parameter_vector[6]),
        "separation": float(separation),
        "kaon_significance": float(kaon_significance),
        "proton_significance": float(proton_significance),
        "chi2_data": chi2_data,
        "chi2_ndf": chi2_ndf,
        "chi2_per_abs_entry": chi2_per_abs_entry,
        "active_bin_count": active_bin_count,
        "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
        "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
        "bound_hit": bool(bound_hit),
        "ordering_valid": bool(ordering_valid),
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "uncertainties": uncertainties,
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "per_aero_fallback": False,
        "peak_pair_found": False,
        "applied_mean_offset": 0.0,
        "offset_adjusted": False,
    }


def _fit_per_aero_fallback_timing_shape(
    histogram,
    function_name,
    display_min,
    display_max,
    proton_peak_is_lower,
    beam_bunch_spacing_ns,
    reference_kaon_mean,
    reference_proton_mean,
    sigma_min,
    minimum_amplitude_significance,
    minimum_entries,
):
    result = {
        "valid": False,
        "fit_status": "insufficient_support",
        "function_name": str(function_name),
        "per_aero_fallback": True,
        "fit_min": float(display_min),
        "fit_max": float(display_max),
    }
    if (
        histogram is None
        or _hist_integral(histogram) < float(minimum_entries)
        or float(display_max) <= float(display_min)
    ):
        return result
    central_min, central_max = _find_central_quantile_range(
        histogram,
        display_min,
        display_max,
        2.0e-4,
        2.0e-4,
        2,
    )
    if central_max <= central_min:
        return result

    candidate_pairs = []

    def add_candidate_pair(lower_mean, upper_mean, lower_height, upper_height):
        if (
            not math.isfinite(lower_mean)
            or not math.isfinite(upper_mean)
            or float(upper_mean) <= float(lower_mean)
            or float(lower_mean) < float(central_min)
            or float(upper_mean) > float(central_max)
        ):
            return
        separation = float(upper_mean - lower_mean)
        if separation < 0.16 or separation > 1.35:
            return
        for existing in candidate_pairs:
            if (
                abs(float(existing["lower_mean"]) - float(lower_mean)) < 0.035
                and abs(float(existing["upper_mean"]) - float(upper_mean)) < 0.035
            ):
                return
        candidate_pairs.append(
            {
                "valid": True,
                "lower_mean": float(lower_mean),
                "upper_mean": float(upper_mean),
                "lower_height": float(lower_height),
                "upper_height": float(upper_height),
            }
        )

    prominent_pair = _find_prominent_offset_peak_pair(
        histogram,
        central_min,
        central_max,
        max(0.18, 0.045 * float(beam_bunch_spacing_ns)),
        min(1.35, 0.70 * float(beam_bunch_spacing_ns)),
    )
    if prominent_pair.get("valid"):
        add_candidate_pair(
            prominent_pair["lower_mean"],
            prominent_pair["upper_mean"],
            prominent_pair["lower_height"],
            prominent_pair["upper_height"],
        )
    separated_pair = _find_separated_peak_pair(
        histogram,
        central_min,
        central_max,
        0.16,
        min(1.35, 0.70 * float(beam_bunch_spacing_ns)),
    )
    if separated_pair.get("valid"):
        add_candidate_pair(
            separated_pair["lower_mean"],
            separated_pair["upper_mean"],
            separated_pair["lower_height"],
            separated_pair["upper_height"],
        )
    dominant_height, dominant_mean = _find_peak_seed(histogram, central_min, central_max)
    for separation in (0.20, 0.28, 0.38, 0.50, 0.64, 0.78, 0.95, 1.15):
        add_candidate_pair(dominant_mean - separation, dominant_mean, 0.25 * dominant_height, dominant_height)
        add_candidate_pair(dominant_mean, dominant_mean + separation, dominant_height, 0.25 * dominant_height)
    if not candidate_pairs:
        return result

    histogram_maximum = max(float(histogram.GetMaximum()), 1.0)
    attempts = []
    best_accepted_index = -1
    best_diagnostic_index = -1
    best_accepted_score = float("inf")
    best_diagnostic_score = float("inf")

    for candidate_index, candidate in enumerate(candidate_pairs):
        seed_separation = float(candidate["upper_mean"] - candidate["lower_mean"])
        split = 0.5 * (float(candidate["lower_mean"]) + float(candidate["upper_mean"]))
        mean_half_width = float(np.clip(0.24 * seed_separation, 0.075, 0.18))
        ordering_gap = max(0.006, 0.018 * seed_separation)
        lower_mean_min = max(float(central_min), float(candidate["lower_mean"]) - mean_half_width)
        lower_mean_max = min(split - ordering_gap, float(candidate["lower_mean"]) + mean_half_width)
        upper_mean_min = max(split + ordering_gap, float(candidate["upper_mean"]) - mean_half_width)
        upper_mean_max = min(float(central_max), float(candidate["upper_mean"]) + mean_half_width)
        if lower_mean_max <= lower_mean_min or upper_mean_max <= upper_mean_min:
            continue
        local_sigma_max = float(np.clip(0.44 * seed_separation, 0.13, 0.32))
        local_sigma_initial = float(np.clip(0.21 * seed_separation, 0.075, 0.15))
        kaon_mean_min = lower_mean_min
        kaon_mean_max = lower_mean_max
        proton_mean_min = upper_mean_min
        proton_mean_max = upper_mean_max
        if bool(proton_peak_is_lower):
            proton_mean_min = lower_mean_min
            proton_mean_max = lower_mean_max
            kaon_mean_min = upper_mean_min
            kaon_mean_max = upper_mean_max
        attempt = _fit_global_timing_shape_with_bounds(
            histogram,
            "{}_try_{}".format(function_name, candidate_index),
            central_min,
            central_max,
            kaon_mean_min,
            kaon_mean_max,
            proton_mean_min,
            proton_mean_max,
            proton_peak_is_lower,
            sigma_min,
            local_sigma_max,
            local_sigma_initial,
            0.45,
            max(1.5, 0.60 * float(minimum_amplitude_significance)),
            1.0e9,
            0.001,
            minimum_entries,
        )
        attempt["fit_min"] = float(central_min)
        attempt["fit_max"] = float(central_max)
        attempt["per_aero_fallback"] = True
        attempt["peak_pair_found"] = True
        fit_status_accepted = bool(attempt.get("fit_status") == "success")
        ordering_accepted = (
            float(attempt.get("proton_mean", 0.0)) < float(attempt.get("kaon_mean", 0.0))
            if bool(proton_peak_is_lower)
            else float(attempt.get("kaon_mean", 0.0)) < float(attempt.get("proton_mean", 0.0))
        )
        smaller_amplitude_fraction = min(
            float(attempt.get("kaon_amplitude", 0.0)),
            float(attempt.get("proton_amplitude", 0.0)),
        ) / histogram_maximum
        required_significance = max(1.5, 0.60 * float(minimum_amplitude_significance))
        physically_accepted = (
            fit_status_accepted
            and ordering_accepted
            and math.isfinite(float(attempt.get("kaon_amplitude", 0.0)))
            and math.isfinite(float(attempt.get("proton_amplitude", 0.0)))
            and math.isfinite(float(attempt.get("kaon_mean", 0.0)))
            and math.isfinite(float(attempt.get("proton_mean", 0.0)))
            and math.isfinite(float(attempt.get("kaon_sigma", 0.0)))
            and math.isfinite(float(attempt.get("proton_sigma", 0.0)))
            and math.isfinite(float(attempt.get("chi2_per_abs_entry", float("inf"))))
            and float(attempt.get("kaon_amplitude", 0.0)) > 0.0
            and float(attempt.get("proton_amplitude", 0.0)) > 0.0
            and float(attempt.get("kaon_sigma", 0.0)) > 0.0
            and float(attempt.get("proton_sigma", 0.0)) > 0.0
            and float(attempt.get("separation", 0.0)) >= 0.55
            and float(attempt.get("kaon_significance", 0.0)) >= required_significance
            and float(attempt.get("proton_significance", 0.0)) >= required_significance
            and smaller_amplitude_fraction >= 0.020
            and float(attempt.get("chi2_per_abs_entry", float("inf"))) <= 0.35
        )
        population_penalty = 0.20 * (0.08 - smaller_amplitude_fraction) if smaller_amplitude_fraction < 0.08 else 0.0
        diagnostic_score = float(attempt.get("chi2_per_abs_entry", float("inf"))) + population_penalty
        attempts.append(attempt)
        stored_index = len(attempts) - 1
        if math.isfinite(diagnostic_score) and diagnostic_score < best_diagnostic_score:
            best_diagnostic_score = diagnostic_score
            best_diagnostic_index = stored_index
        if physically_accepted and diagnostic_score < best_accepted_score:
            best_accepted_score = diagnostic_score
            best_accepted_index = stored_index

    selected_index = best_accepted_index if best_accepted_index >= 0 else best_diagnostic_index
    if selected_index < 0:
        return result
    result = deepcopy(attempts[selected_index])
    result["function_name"] = str(function_name)
    reference_pair_valid = (
        math.isfinite(float(reference_kaon_mean))
        and math.isfinite(float(reference_proton_mean))
        and (
            float(reference_proton_mean) < float(reference_kaon_mean)
            if bool(proton_peak_is_lower)
            else float(reference_kaon_mean) < float(reference_proton_mean)
        )
    )
    if reference_pair_valid:
        fitted_midpoint = 0.5 * (float(result.get("kaon_mean", 0.0)) + float(result.get("proton_mean", 0.0)))
        reference_midpoint = 0.5 * (float(reference_kaon_mean) + float(reference_proton_mean))
        result["applied_mean_offset"] = float(fitted_midpoint - reference_midpoint)
        result["offset_adjusted"] = bool(abs(float(result["applied_mean_offset"])) > 1.0e-6)
    else:
        result["applied_mean_offset"] = 0.0
        result["offset_adjusted"] = False
    final_smaller_amplitude_fraction = min(
        float(result.get("kaon_amplitude", 0.0)),
        float(result.get("proton_amplitude", 0.0)),
    ) / histogram_maximum
    final_required_significance = max(1.5, 0.60 * float(minimum_amplitude_significance))
    final_ordering_accepted = (
        float(result.get("proton_mean", 0.0)) < float(result.get("kaon_mean", 0.0))
        if bool(proton_peak_is_lower)
        else float(result.get("kaon_mean", 0.0)) < float(result.get("proton_mean", 0.0))
    )
    result["valid"] = bool(
        best_accepted_index >= 0
        and result.get("fit_status") == "success"
        and final_ordering_accepted
        and float(result.get("separation", 0.0)) >= 0.55
        and float(result.get("kaon_significance", 0.0)) >= final_required_significance
        and float(result.get("proton_significance", 0.0)) >= final_required_significance
        and final_smaller_amplitude_fraction >= 0.020
        and math.isfinite(float(result.get("chi2_per_abs_entry", float("inf"))))
        and float(result.get("chi2_per_abs_entry", float("inf"))) <= 0.35
    )
    result["per_aero_fallback"] = True
    return result


def _sum_component_over_bins(histogram, values, fit_min, fit_max):
    if histogram is None:
        return 0.0
    total = 0.0
    for bin_index in range(1, histogram.GetNbinsX() + 1):
        x_value = float(histogram.GetXaxis().GetBinCenter(bin_index))
        if x_value < float(fit_min) or x_value > float(fit_max):
            continue
        total += float(values[bin_index - 1])
    return float(total)


def _fit_delta_timing_slice(
    histogram,
    global_shape,
    config,
    function_name,
):
    fit_min = float(global_shape.get("fit_min", (config.get("ctime_hist_range") or (-1.50, 1.25))[0]))
    fit_max = float(global_shape.get("fit_max", (config.get("ctime_hist_range") or (-1.50, 1.25))[1]))
    inputs = _extract_weighted_hist_fit_inputs(histogram, fit_min, fit_max)
    x_values = inputs["x"]
    y_values = inputs["y"]
    sigma_values = inputs["sigma"]
    minimum_entries = int((config.get("slice_fit") or {}).get("minimum_entries", 30))
    if (not global_shape.get("valid")) or x_values.size == 0 or np.sum(np.abs(y_values)) < float(minimum_entries):
        return {
            "valid": False,
            "fit_status": "insufficient_support",
            "function_name": str(function_name),
            "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
            "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
        }

    kaon_component = _gaussian(
        x_values,
        1.0,
        global_shape["kaon_mean"],
        global_shape["kaon_sigma"],
    )
    proton_component = _gaussian(
        x_values,
        1.0,
        global_shape["proton_mean"],
        global_shape["proton_sigma"],
    )
    weighted_design = np.column_stack(
        [
            kaon_component / sigma_values,
            proton_component / sigma_values,
            np.ones_like(x_values) / sigma_values,
        ]
    )
    weighted_target = y_values / sigma_values
    fit_result = lsq_linear(
        weighted_design,
        weighted_target,
        bounds=(0.0, np.inf),
    )
    amplitudes = np.asarray(fit_result.x, dtype=float)
    model_values = (
        amplitudes[0] * kaon_component
        + amplitudes[1] * proton_component
        + amplitudes[2]
    )
    residual_values = (y_values - model_values) / sigma_values
    chi2_data = float(np.sum(np.square(residual_values)))
    active_bin_count = int(x_values.size)
    ndf = max(active_bin_count - 3, 1)
    chi2_ndf = float(chi2_data / float(ndf)) if ndf > 0 else float("inf")
    covariance_matrix, correlation_matrix, uncertainties = _compute_parameter_covariance(
        weighted_design,
        ("kaon_amplitude", "proton_amplitude", "other_amplitude"),
    )
    full_kaon = _gaussian(
        np.asarray(
            [histogram.GetXaxis().GetBinCenter(bin_index) for bin_index in range(1, histogram.GetNbinsX() + 1)],
            dtype=float,
        ),
        amplitudes[0],
        global_shape["kaon_mean"],
        global_shape["kaon_sigma"],
    )
    full_proton = _gaussian(
        np.asarray(
            [histogram.GetXaxis().GetBinCenter(bin_index) for bin_index in range(1, histogram.GetNbinsX() + 1)],
            dtype=float,
        ),
        amplitudes[1],
        global_shape["proton_mean"],
        global_shape["proton_sigma"],
    )
    full_other = np.full(histogram.GetNbinsX(), float(amplitudes[2]), dtype=float)
    kaon_yield = _sum_component_over_bins(histogram, full_kaon, fit_min, fit_max)
    proton_yield = _sum_component_over_bins(histogram, full_proton, fit_min, fit_max)
    other_yield = _sum_component_over_bins(histogram, full_other, fit_min, fit_max)
    data_yield = float(np.sum(y_values))
    model_yield = float(kaon_yield + proton_yield + other_yield)
    model_data_ratio = float(model_yield / data_yield) if data_yield > 0.0 else None
    valid = (
        bool(getattr(fit_result, "success", False))
        and math.isfinite(chi2_ndf)
        and float(amplitudes[0]) >= 0.0
        and float(amplitudes[1]) >= 0.0
        and float(amplitudes[2]) >= 0.0
        and model_yield > 0.0
        and model_data_ratio is not None
        and model_data_ratio >= float((config.get("slice_fit") or {}).get("minimum_model_data_ratio", 0.50))
        and model_data_ratio <= float((config.get("slice_fit") or {}).get("maximum_model_data_ratio", 1.50))
        and chi2_ndf <= float((config.get("slice_fit") or {}).get("maximum_chi2_ndf", 5.0))
    )
    return {
        "valid": bool(valid),
        "fit_status": "success" if bool(getattr(fit_result, "success", False)) else "failure",
        "message": str(getattr(fit_result, "message", "")),
        "function_name": str(function_name),
        "kaon_amplitude": float(amplitudes[0]),
        "kaon_amplitude_error": float(uncertainties.get("kaon_amplitude", 0.0)),
        "proton_amplitude": float(amplitudes[1]),
        "proton_amplitude_error": float(uncertainties.get("proton_amplitude", 0.0)),
        "other_amplitude": float(amplitudes[2]),
        "other_amplitude_error": float(uncertainties.get("other_amplitude", 0.0)),
        "kaon_yield": float(kaon_yield),
        "proton_yield": float(proton_yield),
        "other_yield": float(other_yield),
        "data_yield": float(data_yield),
        "model_yield": float(model_yield),
        "model_data_ratio": float(model_data_ratio) if model_data_ratio is not None else None,
        "chi2_data": chi2_data,
        "chi2_ndf": chi2_ndf,
        "active_bin_count": active_bin_count,
        "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
        "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "uncertainties": uncertainties,
    }


def _resolve_standard_peak_bounds(
    histogram,
    display_min,
    display_max,
    proton_peak_is_lower,
    beam_bunch_spacing_ns,
):
    peak_pair = _find_prominent_offset_peak_pair(
        histogram,
        display_min,
        display_max,
        max(0.18, 0.045 * float(beam_bunch_spacing_ns)),
        min(1.35, 0.70 * float(beam_bunch_spacing_ns)),
    )
    if not peak_pair.get("valid"):
        peak_pair = _find_separated_peak_pair(
            histogram,
            display_min,
            display_max,
            0.16,
            min(1.35, 0.70 * float(beam_bunch_spacing_ns)),
        )
    if peak_pair.get("valid"):
        lower_mean = float(peak_pair["lower_mean"])
        upper_mean = float(peak_pair["upper_mean"])
    else:
        _, dominant_mean = _find_peak_seed(histogram, display_min, display_max)
        default_separation = float(np.clip(0.30 * float(beam_bunch_spacing_ns), 0.24, 0.90))
        lower_mean = float(dominant_mean - (0.5 * default_separation))
        upper_mean = float(dominant_mean + (0.5 * default_separation))
    separation = max(float(upper_mean - lower_mean), 0.18)
    split = 0.5 * (lower_mean + upper_mean)
    mean_half_width = float(np.clip(0.24 * separation, 0.075, 0.18))
    ordering_gap = max(0.006, 0.018 * separation)
    lower_mean_min = max(float(display_min), lower_mean - mean_half_width)
    lower_mean_max = min(split - ordering_gap, lower_mean + mean_half_width)
    upper_mean_min = max(split + ordering_gap, upper_mean - mean_half_width)
    upper_mean_max = min(float(display_max), upper_mean + mean_half_width)
    if lower_mean_max <= lower_mean_min:
        lower_mean_min = max(float(display_min), lower_mean - 0.10)
        lower_mean_max = min(split - 0.005, lower_mean + 0.10)
    if upper_mean_max <= upper_mean_min:
        upper_mean_min = max(split + 0.005, upper_mean - 0.10)
        upper_mean_max = min(float(display_max), upper_mean + 0.10)
    local_sigma_max = float(np.clip(0.44 * separation, 0.13, 0.32))
    local_sigma_initial = float(np.clip(0.21 * separation, 0.075, 0.15))
    kaon_mean_min = lower_mean_min
    kaon_mean_max = lower_mean_max
    proton_mean_min = upper_mean_min
    proton_mean_max = upper_mean_max
    if bool(proton_peak_is_lower):
        proton_mean_min = lower_mean_min
        proton_mean_max = lower_mean_max
        kaon_mean_min = upper_mean_min
        kaon_mean_max = upper_mean_max
    return {
        "kaon_mean_min": float(kaon_mean_min),
        "kaon_mean_max": float(kaon_mean_max),
        "proton_mean_min": float(proton_mean_min),
        "proton_mean_max": float(proton_mean_max),
        "sigma_max": float(local_sigma_max),
        "sigma_initial": float(local_sigma_initial),
        "peak_pair_found": bool(peak_pair.get("valid", False)),
    }


def _fit_global_timing_shape(
    histogram,
    config,
    function_name,
    proton_peak_is_lower=False,
    display_range=None,
    fit_mode="per_aero_multistart",
):
    global_cfg = config.get("global_fit") or {}
    display_min, display_max = [
        float(value)
        for value in (
            display_range
            or config.get("ctime_hist_range")
            or (-1.50, 1.25)
        )
    ]
    beam_bunch_spacing_ns = float(global_cfg.get("beam_bunch_spacing_ns", 2.004) or 2.004)
    sigma_min, sigma_max = [
        float(value)
        for value in (
            global_cfg.get("sigma_range")
            or (0.03, 0.45)
        )
    ]
    minimum_amplitude_significance = float(global_cfg.get("minimum_amplitude_significance", 2.0) or 2.0)
    minimum_entries = int(global_cfg.get("minimum_entries", 200) or 200)
    if str(fit_mode) == "per_aero_multistart":
        return _fit_per_aero_fallback_timing_shape(
            histogram,
            function_name,
            display_min,
            display_max,
            proton_peak_is_lower,
            beam_bunch_spacing_ns,
            float("nan"),
            float("nan"),
            sigma_min,
            minimum_amplitude_significance,
            minimum_entries,
        )
    if str(fit_mode) == "local_peak_rescue":
        return _fit_per_aero_fallback_timing_shape(
            histogram,
            function_name,
            display_min,
            display_max,
            proton_peak_is_lower,
            beam_bunch_spacing_ns,
            float("nan"),
            float("nan"),
            sigma_min,
            max(1.5, 0.75 * minimum_amplitude_significance),
            max(40, int(0.5 * minimum_entries)),
        )
    bounds = _resolve_standard_peak_bounds(
        histogram,
        display_min,
        display_max,
        proton_peak_is_lower,
        beam_bunch_spacing_ns,
    )
    shape = _fit_global_timing_shape_with_bounds(
        histogram,
        function_name,
        display_min,
        display_max,
        bounds["kaon_mean_min"],
        bounds["kaon_mean_max"],
        bounds["proton_mean_min"],
        bounds["proton_mean_max"],
        proton_peak_is_lower,
        sigma_min,
        bounds["sigma_max"],
        bounds["sigma_initial"],
        float(global_cfg.get("minimum_separation", 0.75) or 0.75),
        minimum_amplitude_significance,
        float(global_cfg.get("maximum_chi2_ndf", 5.0) or 5.0),
        float(global_cfg.get("bound_fraction_tolerance", 0.02) or 0.02),
        minimum_entries,
    )
    shape["peak_pair_found"] = bool(bounds.get("peak_pair_found", False))
    return shape


def _summarize_global_probe(
    branch_name,
    probe_kind,
    fit_mode,
    pid_payload,
    global_shapes,
    proton_peak_is_lower,
):
    valid_shapes = [shape for shape in global_shapes if bool(shape.get("valid"))]
    active_shapes = valid_shapes if valid_shapes else [shape for shape in global_shapes if shape]
    fit_min = float(pid_payload.get("time_hist_range", (0.0, 0.0))[0])
    fit_max = float(pid_payload.get("time_hist_range", (0.0, 0.0))[1])
    histogram_bins = int(pid_payload.get("time_hist_bins", 90) or 90)
    mean_separation = (
        float(np.mean([float(shape.get("separation", 0.0) or 0.0) for shape in valid_shapes]))
        if valid_shapes
        else 0.0
    )
    mean_chi2_ndf = (
        float(np.mean([float(shape.get("chi2_ndf", float("inf")) or float("inf")) for shape in valid_shapes]))
        if valid_shapes
        else float("inf")
    )
    mean_chi2_per_abs_entry = (
        float(np.mean([float(shape.get("chi2_per_abs_entry", float("inf")) or float("inf")) for shape in valid_shapes]))
        if valid_shapes
        else float("inf")
    )
    peak_pair_found = bool(any(shape.get("peak_pair_found") for shape in active_shapes))
    return {
        "branch": str(branch_name),
        "probe_kind": str(probe_kind),
        "fit_mode": str(fit_mode),
        "timing_branch": str(branch_name),
        "displayMin": float(fit_min),
        "displayMax": float(fit_max),
        "fitMin": float(fit_min),
        "fitMax": float(fit_max),
        "histogramBins": int(histogram_bins),
        "validShapes": int(len(valid_shapes)),
        "meanSeparation": float(mean_separation),
        "meanPoissonDevianceNdf": float(mean_chi2_ndf),
        "meanPoissonDeviancePerEntry": float(mean_chi2_per_abs_entry),
        "peakPairFound": bool(peak_pair_found),
        "perAeroMultistart": str(fit_mode) == "per_aero_multistart",
        "localPeakRescue": str(fit_mode) == "local_peak_rescue",
        "proton_peak_is_lower": bool(proton_peak_is_lower),
        "pid_payload": pid_payload,
        "global_shapes": global_shapes,
    }


def _compare_timing_probes(left_probe, right_probe):
    if int(left_probe.get("validShapes", 0)) != int(right_probe.get("validShapes", 0)):
        return 1 if int(left_probe.get("validShapes", 0)) > int(right_probe.get("validShapes", 0)) else -1
    if bool(left_probe.get("peakPairFound", False)) != bool(right_probe.get("peakPairFound", False)):
        return 1 if bool(left_probe.get("peakPairFound", False)) else -1
    if abs(float(left_probe.get("meanSeparation", 0.0)) - float(right_probe.get("meanSeparation", 0.0))) > 1.0e-9:
        return 1 if float(left_probe.get("meanSeparation", 0.0)) > float(right_probe.get("meanSeparation", 0.0)) else -1
    left_goodness = (
        float(left_probe.get("meanPoissonDeviancePerEntry", float("inf")))
        if bool(left_probe.get("perAeroMultistart", False)) or bool(left_probe.get("localPeakRescue", False))
        else float(left_probe.get("meanPoissonDevianceNdf", float("inf")))
    )
    right_goodness = (
        float(right_probe.get("meanPoissonDeviancePerEntry", float("inf")))
        if bool(right_probe.get("perAeroMultistart", False)) or bool(right_probe.get("localPeakRescue", False))
        else float(right_probe.get("meanPoissonDevianceNdf", float("inf")))
    )
    left_finite = math.isfinite(left_goodness)
    right_finite = math.isfinite(right_goodness)
    if left_finite != right_finite:
        return 1 if left_finite else -1
    if left_finite and abs(left_goodness - right_goodness) > 1.0e-9:
        return 1 if left_goodness < right_goodness else -1
    return 0


def _run_timing_probe(
    source_bundle,
    config,
    evaluate_event,
    shifted_mm_getter,
    hole_contains,
    mm_min,
    mm_max,
    timing_branch,
    proton_peak_is_lower,
    probe_kind,
    fit_mode,
):
    display_range = None
    if str(probe_kind) == "rf":
        branch_values = _collect_branch_values(
            source_bundle,
            evaluate_event,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch,
        )
        display_range = _estimate_value_central_range(
            branch_values,
            *(config.get("ctime_hist_range") or (-1.50, 1.25))
        )
    else:
        display_range = tuple(float(value) for value in (config.get("ctime_hist_range") or (-1.50, 1.25)))
    pid_payload = _build_signed_pid_histograms(
        source_bundle,
        config,
        evaluate_event,
        shifted_mm_getter,
        hole_contains,
        mm_min,
        mm_max,
        timing_branch=str(timing_branch),
        time_hist_range=display_range,
        time_hist_bins=90,
    )
    global_shapes = []
    for aero_index, slice_hist in enumerate(pid_payload["global_slice_hists"]):
        shape = _fit_global_timing_shape(
            slice_hist,
            config,
            "f_proton_cleaning_global_{}_aero_{}".format(
                str(timing_branch).replace(" ", "_"),
                aero_index,
            ),
            proton_peak_is_lower=proton_peak_is_lower,
            display_range=display_range,
            fit_mode=fit_mode,
        )
        global_shapes.append(shape)
    return _summarize_global_probe(
        timing_branch,
        probe_kind,
        fit_mode,
        pid_payload,
        global_shapes,
        proton_peak_is_lower,
    )


def _find_collection_bin(value, edges):
    if len(edges) < 2:
        return -1
    value = float(value)
    if value < float(edges[0]) or value > float(edges[-1]):
        return -1
    if value == float(edges[-1]):
        return len(edges) - 2
    for index in range(len(edges) - 1):
        low_edge = float(edges[index])
        high_edge = float(edges[index + 1])
        if low_edge <= value < high_edge:
            return index
    return -1


def _evaluate_event_proton_probability(
    timing,
    global_shape,
    slice_fit,
    denominator_floor,
):
    if not global_shape.get("valid") or not slice_fit.get("valid"):
        return 0.0
    proton_value = float(
        _gaussian(
            np.asarray([timing], dtype=float),
            slice_fit["proton_amplitude"],
            global_shape["proton_mean"],
            global_shape["proton_sigma"],
        )[0]
    )
    kaon_value = float(
        _gaussian(
            np.asarray([timing], dtype=float),
            slice_fit["kaon_amplitude"],
            global_shape["kaon_mean"],
            global_shape["kaon_sigma"],
        )[0]
    )
    other_value = max(float(slice_fit.get("other_amplitude", 0.0) or 0.0), 0.0)
    denominator = proton_value + kaon_value + other_value
    if (not math.isfinite(denominator)) or denominator <= float(denominator_floor):
        return 0.0
    return max(0.0, min(1.0, float(proton_value / denominator)))


def _build_signed_pid_histograms(
    source_bundle,
    config,
    evaluate_event,
    shifted_mm_getter,
    hole_contains,
    mm_min,
    mm_max,
    timing_branch,
    time_hist_range=None,
    time_hist_bins=None,
):
    aero_edges = [float(edge) for edge in (config.get("aero_slice_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0))]
    aero_min, aero_max = [float(value) for value in (config.get("aero_hist_range") or (0.0, 25.0))]
    if time_hist_range is None:
        time_hist_range = config.get("ctime_hist_range") or (-1.50, 1.25)
    time_min, time_max = [float(value) for value in time_hist_range]
    n_time_bins = int(time_hist_bins or 90)
    delta_min, delta_max = [float(value) for value in (config.get("delta_hist_range") or (-10.0, 20.0))]
    delta_bins = int(config.get("delta_bins", 10) or 10)

    h_global_pid = ROOT.TH2D(
        "H_proton_cleaning_global_pid",
        "Global {} vs P_aero NPE;P_aero NPE;{} [ns];Signed yield".format(timing_branch, timing_branch),
        75,
        aero_min,
        aero_max,
        n_time_bins,
        time_min,
        time_max,
    )
    h_global_pid.SetDirectory(0)
    h_global_pid.Sumw2()

    global_slice_hists = []
    for aero_index in range(len(aero_edges) - 1):
        hist = ROOT.TH1D(
            "H_proton_cleaning_global_time_slice_{}".format(aero_index),
            "Global timing slice {};{} [ns];Signed yield".format(aero_index + 1, timing_branch),
            n_time_bins,
            time_min,
            time_max,
        )
        hist.SetDirectory(0)
        hist.Sumw2()
        global_slice_hists.append(hist)

    delta_pid_hists = []
    delta_slice_hists = []
    for delta_index in range(delta_bins):
        pid_hist = ROOT.TH2D(
            "H_proton_cleaning_pid_delta_{}".format(delta_index),
            "PID slice {} ;P_aero NPE;{} [ns];Signed yield".format(delta_index + 1, timing_branch),
            75,
            aero_min,
            aero_max,
            n_time_bins,
            time_min,
            time_max,
        )
        pid_hist.SetDirectory(0)
        pid_hist.Sumw2()
        delta_pid_hists.append(pid_hist)

        slice_hists = []
        for aero_index in range(len(aero_edges) - 1):
            hist = ROOT.TH1D(
                "H_proton_cleaning_time_delta_{}_aero_{}".format(delta_index, aero_index),
                "Timing slice d{} a{};{} [ns];Signed yield".format(delta_index + 1, aero_index + 1, timing_branch),
                n_time_bins,
                time_min,
                time_max,
            )
            hist.SetDirectory(0)
            hist.Sumw2()
            slice_hists.append(hist)
        delta_slice_hists.append(slice_hists)

    source_stats = {}
    for source_name, source_spec in (source_bundle.get("sources") or {}).items():
        tree = source_spec.get("tree")
        coefficient = float(source_spec.get("coefficient", 0.0) or 0.0)
        source_stats[source_name] = {
            "tree_name": source_spec.get("tree_name"),
            "entries_seen": 0,
            "entries_used": 0,
        }
        if tree is None or coefficient == 0.0:
            continue
        for evt in tree:
            source_stats[source_name]["entries_seen"] += 1
            allcuts, nommcuts, _ = evaluate_event(evt, mm_min, mm_max)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) if hole_contains is not None else False
            if not (nommcuts and not hole_rejected):
                continue
            aero_value = float(getattr(evt, "P_aero_npeSum", 0.0))
            try:
                timing_value = float(getattr(evt, str(timing_branch)))
            except Exception:
                continue
            delta_value = float(getattr(evt, "ssdelta", 0.0))
            if (aero_value < aero_min) or (aero_value > aero_max) or (timing_value < time_min) or (timing_value > time_max):
                continue
            source_stats[source_name]["entries_used"] += 1
            h_global_pid.Fill(aero_value, timing_value, coefficient)
            aero_index = _find_collection_bin(aero_value, aero_edges)
            delta_width = (delta_max - delta_min) / float(delta_bins)
            delta_index = _find_collection_bin(delta_value, [delta_min + (delta_width * i) for i in range(delta_bins + 1)])
            if 0 <= aero_index < len(global_slice_hists):
                global_slice_hists[aero_index].Fill(timing_value, coefficient)
            if 0 <= delta_index < delta_bins:
                delta_pid_hists[delta_index].Fill(aero_value, timing_value, coefficient)
                if 0 <= aero_index < len(delta_slice_hists[delta_index]):
                    delta_slice_hists[delta_index][aero_index].Fill(timing_value, coefficient)

    return {
        "H_global_pid": h_global_pid,
        "global_slice_hists": global_slice_hists,
        "delta_pid_hists": delta_pid_hists,
        "delta_slice_hists": delta_slice_hists,
        "source_stats": source_stats,
        "aero_edges": aero_edges,
        "delta_edges": [delta_min + (((delta_max - delta_min) / float(delta_bins)) * i) for i in range(delta_bins + 1)],
        "timing_branch": str(timing_branch),
        "time_hist_range": (float(time_min), float(time_max)),
        "time_hist_bins": int(n_time_bins),
    }


def build_kaon_proton_cleaning_result(
    inpDict,
    phi_setting,
    source_bundle,
    evaluate_event,
    shifted_mm_getter,
    hole_contains,
    mm_min,
    mm_max,
    analysis_scope="setting-wide",
    context="",
):
    config = get_proton_contamination_cleaning_config(
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    method = resolve_proton_contamination_cleaning_method(
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    enabled = resolve_proton_contamination_cleaning_enabled(
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    tree_policy = resolve_proton_contamination_cleaning_tree_policy(
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    result = {
        "accepted": False,
        "enabled": bool(enabled),
        "method": method,
        "analysis_scope": str(analysis_scope),
        "context": str(context),
        "particle_type": "kaon",
        "phi_setting": phi_setting,
        "settings": config,
        "diagnostics": {
            "analysis_scope": str(analysis_scope),
            "context": str(context),
            "method": method,
            "resolved_setting_key": config.get("proton_contamination_setting_key"),
            "resolved_phi_setting": config.get("proton_contamination_phi_setting"),
            "override_layers": list(config.get("proton_contamination_override_layers") or []),
            "strict_mode": bool(config.get("strict_mode", False)),
            "failure_policy": str(config.get("failure_policy", "bypass")),
            "tree_policy": str(tree_policy),
            "rf_policy": str(config.get("rf_policy", "epsset_default_after_cleaning")),
            "input_tree_state": "explicit_noRF",
            "input_tree_particle_selection": "Cut_Kaon_Events_*_noRF",
            "fit_sample_definition": (
                "signed noRF prompt/random/dummy bundle with evaluate_data_event "
                "nommcuts and HGCer-hole rejection"
            ),
            "fit_sample_signed_combination": (
                "prompt - rand/nWindows - dummy + dummy_rand/nWindows"
            ),
            "fit_sample_uses_signed_random_dummy_subtraction": True,
            "fit_sample_requires_nommcuts": True,
            "fit_sample_requires_hgcer_hole_rejection": True,
        },
        "fallback_reason": "",
    }
    if method == PROTON_CONTAMINATION_CLEANING_METHOD_DISABLED or not enabled:
        result["fallback_reason"] = "proton cleaning disabled"
        return result
    if method != PROTON_CONTAMINATION_CLEANING_METHOD_CTIME_AERO_EVENT_WEIGHT:
        raise ValueError("Unsupported proton-cleaning method '{}'".format(method))

    missing_norf = []
    for source_name in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
        tree = ((source_bundle or {}).get("sources") or {}).get(source_name, {}).get("tree")
        if tree is None:
            missing_norf.append(source_name)
    if missing_norf:
        reason = "missing explicit noRF kaon trees: {}".format(", ".join(missing_norf))
        result["fallback_reason"] = reason
        if bool(config.get("strict_mode", False)) and not bool(config.get("allow_missing_norf_trees", False)):
            raise RuntimeError(reason)
        return result

    result["diagnostics"]["source_coefficients"] = {
        str(source_name): float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items()
    }
    ct_time_branch = str(config.get("ct_timing_branch", "CTime_ROC1") or "CTime_ROC1")
    rf_candidate_branches = _resolve_rf_branch_candidates(source_bundle)
    result["diagnostics"]["rf_candidate_branches"] = list(rf_candidate_branches)
    result["diagnostics"]["rf_timing_attempted"] = bool(rf_candidate_branches)
    result["diagnostics"]["ct_timing_evaluated"] = True
    result["diagnostics"]["timingFitUsedPerAeroDefault"] = True
    result["diagnostics"]["timingFitUsedStandardFallback"] = False
    result["diagnostics"]["timingFitUsedPreDiamondFallback"] = False
    result["diagnostics"]["timingFitUsedLocalPeakRescue"] = False

    def _probe_summary(probe):
        if not isinstance(probe, dict):
            return {}
        summary = {
            key: value
            for key, value in probe.items()
            if key not in ("pid_payload", "global_shapes")
        }
        return _json_ready_value(summary)

    def _best_rf_probe_or_none(probes):
        best_probe = None
        for probe in probes or []:
            if best_probe is None or _compare_timing_probes(probe, best_probe) > 0:
                best_probe = probe
        return best_probe

    current_fit_mode = "per_aero_multistart"
    rf_probes = [
        _run_timing_probe(
            source_bundle,
            config,
            evaluate_event,
            shifted_mm_getter,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch=candidate,
            proton_peak_is_lower=True,
            probe_kind="rf",
            fit_mode=current_fit_mode,
        )
        for candidate in rf_candidate_branches
    ]
    ct_probe = _run_timing_probe(
        source_bundle,
        config,
        evaluate_event,
        shifted_mm_getter,
        hole_contains,
        mm_min,
        mm_max,
        timing_branch=ct_time_branch,
        proton_peak_is_lower=False,
        probe_kind="ct",
        fit_mode=current_fit_mode,
    )
    best_rf_probe = _best_rf_probe_or_none(rf_probes)
    best_rf_valid_shapes = max((int(probe.get("validShapes", 0) or 0) for probe in rf_probes), default=0)

    if best_rf_valid_shapes == 0 and int(ct_probe.get("validShapes", 0) or 0) == 0:
        current_fit_mode = "standard"
        result["diagnostics"]["timingFitUsedPerAeroDefault"] = False
        result["diagnostics"]["timingFitUsedStandardFallback"] = True
        rf_probes = [
            _run_timing_probe(
                source_bundle,
                config,
                evaluate_event,
                shifted_mm_getter,
                hole_contains,
                mm_min,
                mm_max,
                timing_branch=candidate,
                proton_peak_is_lower=True,
                probe_kind="rf",
                fit_mode=current_fit_mode,
            )
            for candidate in rf_candidate_branches
        ]
        ct_probe = _run_timing_probe(
            source_bundle,
            config,
            evaluate_event,
            shifted_mm_getter,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch=ct_time_branch,
            proton_peak_is_lower=False,
            probe_kind="ct",
            fit_mode=current_fit_mode,
        )
        best_rf_probe = _best_rf_probe_or_none(rf_probes)
        best_rf_valid_shapes = max((int(probe.get("validShapes", 0) or 0) for probe in rf_probes), default=0)

    if best_rf_valid_shapes == 0 and int(ct_probe.get("validShapes", 0) or 0) == 0:
        current_fit_mode = "local_peak_rescue"
        result["diagnostics"]["timingFitUsedStandardFallback"] = False
        result["diagnostics"]["timingFitUsedLocalPeakRescue"] = True
        rf_probes = [
            _run_timing_probe(
                source_bundle,
                config,
                evaluate_event,
                shifted_mm_getter,
                hole_contains,
                mm_min,
                mm_max,
                timing_branch=candidate,
                proton_peak_is_lower=True,
                probe_kind="rf",
                fit_mode=current_fit_mode,
            )
            for candidate in rf_candidate_branches
        ]
        ct_probe = _run_timing_probe(
            source_bundle,
            config,
            evaluate_event,
            shifted_mm_getter,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch=ct_time_branch,
            proton_peak_is_lower=False,
            probe_kind="ct",
            fit_mode=current_fit_mode,
        )
        best_rf_probe = _best_rf_probe_or_none(rf_probes)
        best_rf_valid_shapes = max((int(probe.get("validShapes", 0) or 0) for probe in rf_probes), default=0)

    selected_probe = ct_probe
    timing_selection_reason = "ct_probe_ranked_better_than_rf"
    if best_rf_probe is None:
        timing_selection_reason = (
            "no_rf_branch_ct_selected"
            if rf_candidate_branches
            else "rf_disabled_or_unavailable_ct_selected"
        )
    else:
        if int(best_rf_probe.get("validShapes", 0) or 0) > 0 and _compare_timing_probes(best_rf_probe, ct_probe) >= 0:
            selected_probe = best_rf_probe
            timing_selection_reason = (
                "rf_probe_ranked_better_than_ct"
                if _compare_timing_probes(best_rf_probe, ct_probe) > 0
                else "rf_won_exact_probe_tie"
            )
        elif int(best_rf_probe.get("validShapes", 0) or 0) == 0 and int(ct_probe.get("validShapes", 0) or 0) == 0 and bool(best_rf_probe.get("peakPairFound", False)) and not bool(ct_probe.get("peakPairFound", False)):
            selected_probe = best_rf_probe
            timing_selection_reason = "rf_resolved_peak_pair_used_after_all_fit_validation_failed"

    if current_fit_mode == "per_aero_multistart":
        timing_selection_reason = "per_aerogel_multistart_default_{}".format(timing_selection_reason)
    elif current_fit_mode == "standard":
        timing_selection_reason = "standard_fit_fallback_{}".format(timing_selection_reason)
    elif current_fit_mode == "local_peak_rescue":
        timing_selection_reason = "local_peak_rescue_{}".format(timing_selection_reason)

    pid_payload = selected_probe.get("pid_payload") or {}
    global_shapes = selected_probe.get("global_shapes") or []
    valid_global_shape_count = int(selected_probe.get("validShapes", 0) or 0)

    result["diagnostics"]["rf_probe_summaries"] = [_probe_summary(probe) for probe in rf_probes]
    result["diagnostics"]["ct_probe_summary"] = _probe_summary(ct_probe)
    result["diagnostics"]["selected_timing_branch"] = str(selected_probe.get("timing_branch", ct_time_branch))
    result["diagnostics"]["selected_probe_kind"] = str(selected_probe.get("probe_kind", "ct"))
    result["diagnostics"]["selected_probe_fit_mode"] = str(current_fit_mode)
    result["diagnostics"]["rf_timing_selected"] = bool(selected_probe.get("probe_kind") == "rf")
    result["diagnostics"]["timingSelectionReason"] = str(timing_selection_reason)
    result["diagnostics"]["valid_global_shape_count"] = int(valid_global_shape_count)
    result["diagnostics"]["source_stats"] = pid_payload.get("source_stats") or {}
    result["selected_timing_branch"] = str(selected_probe.get("timing_branch", ct_time_branch))
    result["global_shapes"] = global_shapes
    result["H_global_pid"] = pid_payload.get("H_global_pid")
    result["H_global_timing_slices"] = pid_payload.get("global_slice_hists") or []
    result["H_delta_pid"] = pid_payload.get("delta_pid_hists") or []
    result["H_delta_timing_slices"] = pid_payload.get("delta_slice_hists") or []
    result["aero_edges"] = pid_payload.get("aero_edges") or []
    result["delta_edges"] = pid_payload.get("delta_edges") or []
    result["diagnostics"]["selected_time_hist_range"] = list(pid_payload.get("time_hist_range") or ())
    result["diagnostics"]["selected_time_hist_bins"] = int(pid_payload.get("time_hist_bins", 90) or 90)

    if valid_global_shape_count <= 0:
        result["fallback_reason"] = "no identifiable proton-kaon timing shapes"
        return result
    delta_fits = []
    support_by_delta = []
    proton_yield_by_delta = []
    kaon_yield_by_delta = []
    other_yield_by_delta = []
    data_yield_by_delta = []
    fitted_data_yield_by_delta = []
    model_yield_by_delta = []
    chi2_ndf_by_delta = []
    valid_slices_by_delta = []
    valid_coverage_by_delta = []
    for delta_index, slice_collection in enumerate(pid_payload["delta_slice_hists"]):
        slice_fits = []
        proton_total = 0.0
        kaon_total = 0.0
        other_total = 0.0
        data_total = _hist_abs_integral(pid_payload["delta_pid_hists"][delta_index])
        fitted_data_total = 0.0
        model_total = 0.0
        valid_slices = 0
        chi2_weighted_sum = 0.0
        chi2_weight = 0.0
        for aero_index, slice_hist in enumerate(slice_collection):
            global_shape = global_shapes[aero_index] if aero_index < len(global_shapes) else {"valid": False}
            slice_fit = _fit_delta_timing_slice(
                slice_hist,
                global_shape,
                config,
                "f_proton_cleaning_delta_{}_aero_{}".format(delta_index, aero_index),
            )
            slice_fits.append(slice_fit)
            if not slice_fit.get("valid"):
                continue
            valid_slices += 1
            proton_total += float(slice_fit.get("proton_yield", 0.0) or 0.0)
            kaon_total += float(slice_fit.get("kaon_yield", 0.0) or 0.0)
            other_total += float(slice_fit.get("other_yield", 0.0) or 0.0)
            fitted_data_total += abs(float(slice_fit.get("data_yield", 0.0) or 0.0))
            model_total += float(slice_fit.get("model_yield", 0.0) or 0.0)
            chi2_weighted_sum += float(slice_fit.get("chi2_ndf", 0.0) or 0.0) * abs(float(slice_fit.get("data_yield", 0.0) or 0.0))
            chi2_weight += abs(float(slice_fit.get("data_yield", 0.0) or 0.0))
        coverage = float(fitted_data_total / data_total) if data_total > 0.0 else 0.0
        support_thresholds = config.get("support_thresholds") or {}
        if (
            valid_slices >= int(support_thresholds.get("minimum_supported_slices", 2))
            and coverage >= float(support_thresholds.get("minimum_supported_coverage", 0.35))
            and (proton_total + kaon_total + other_total) >= float(support_thresholds.get("minimum_modeled_yield", 5.0))
        ):
            support_label = SUPPORT_SUPPORTED
        elif (
            valid_slices >= int(support_thresholds.get("minimum_marginal_slices", 1))
            and coverage >= float(support_thresholds.get("minimum_marginal_coverage", 0.15))
            and (proton_total + kaon_total + other_total) >= float(support_thresholds.get("minimum_modeled_yield", 5.0))
        ):
            support_label = SUPPORT_MARGINAL
        else:
            support_label = SUPPORT_UNSUPPORTED
        delta_fits.append(slice_fits)
        support_by_delta.append(support_label)
        proton_yield_by_delta.append(float(proton_total))
        kaon_yield_by_delta.append(float(kaon_total))
        other_yield_by_delta.append(float(other_total))
        data_yield_by_delta.append(float(data_total))
        fitted_data_yield_by_delta.append(float(fitted_data_total))
        model_yield_by_delta.append(float(model_total))
        chi2_ndf_by_delta.append(float(chi2_weighted_sum / chi2_weight) if chi2_weight > 0.0 else None)
        valid_slices_by_delta.append(int(valid_slices))
        valid_coverage_by_delta.append(float(coverage))
    result["delta_slice_fits"] = delta_fits
    result["support_by_delta"] = support_by_delta
    result["diagnostics"]["proton_yield_by_delta"] = proton_yield_by_delta
    result["diagnostics"]["kaon_yield_by_delta"] = kaon_yield_by_delta
    result["diagnostics"]["other_yield_by_delta"] = other_yield_by_delta
    result["diagnostics"]["data_yield_by_delta"] = data_yield_by_delta
    result["diagnostics"]["model_yield_by_delta"] = model_yield_by_delta
    result["diagnostics"]["valid_slices_by_delta"] = valid_slices_by_delta
    result["diagnostics"]["valid_coverage_by_delta"] = valid_coverage_by_delta
    result["diagnostics"]["chi2_ndf_by_delta"] = chi2_ndf_by_delta
    result["diagnostics"]["supported_delta_bins"] = int(sum(1 for label in support_by_delta if label == SUPPORT_SUPPORTED))
    result["diagnostics"]["marginal_delta_bins"] = int(sum(1 for label in support_by_delta if label == SUPPORT_MARGINAL))

    if not any(label in (SUPPORT_SUPPORTED, SUPPORT_MARGINAL) for label in support_by_delta):
        result["fallback_reason"] = "no supported or marginal delta bins for proton cleaning"
        return result

    rf_policy = resolve_proton_contamination_cleaning_rf_policy(
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    apply_rf = (
        rf_policy == "epsset_default_after_cleaning"
        and bool(config.get("apply_only_low_epsilon_rf", True))
        and normalize_epsset(inpDict.get("EPSSET")) == "low"
    )
    result["diagnostics"]["rf_policy"] = rf_policy
    result["diagnostics"]["rf_applied"] = bool(apply_rf)
    result["diagnostics"]["rf_signature_fields"] = list(config.get("rf_signature_fields") or [])
    result["diagnostics"]["signature_round_digits"] = int(config.get("signature_round_digits", 9) or 9)
    if apply_rf:
        rf_sources = (source_bundle or {}).get("rf_sources") or {}
        missing_rf = [name for name in ("prompt", "rand", "dummy_prompt", "dummy_rand") if rf_sources.get(name, {}).get("tree") is None]
        if missing_rf:
            reason = "missing explicit RF reference trees: {}".format(", ".join(missing_rf))
            result["fallback_reason"] = reason
            if bool(config.get("strict_mode", False)) and not bool(config.get("allow_missing_rf_reference", False)):
                raise RuntimeError(reason)
            return result
        lookup, lookup_counts = _build_rf_membership_lookup(
            {key: value.get("tree") for key, value in rf_sources.items()},
            config.get("rf_signature_fields") or (),
            int(config.get("signature_round_digits", 9) or 9),
        )
        result["_rf_signature_lookup"] = lookup
        result["diagnostics"]["rf_lookup_counts"] = lookup_counts
        result["diagnostics"]["rf_source"] = "existing_RF_skim_tree_membership"
    else:
        result["diagnostics"]["rf_lookup_counts"] = {}
        result["diagnostics"]["rf_source"] = "not_applied"

    result["accepted"] = True
    return result


def _fill_standard_target_hists(evt, adj_MM, adj_t, adj_hsdelta, weight, target_hists, allcuts=False, nommcuts=False, noholecuts=False):
    if not math.isfinite(weight) or weight == 0.0:
        return
    if noholecuts:
        if target_hists.get("hgcer_xy_nohole") is not None:
            target_hists["hgcer_xy_nohole"].Fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, weight * evt.P_hgcer_npeSum)
        if target_hists.get("hgcer_x_mm_nohole") is not None:
            target_hists["hgcer_x_mm_nohole"].Fill(evt.P_hgcer_xAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
        if target_hists.get("hgcer_y_mm_nohole") is not None:
            target_hists["hgcer_y_mm_nohole"].Fill(evt.P_hgcer_yAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
    if nommcuts:
        for key in ("h_mm_nosub", "h_mm_fit2sub", "h_mm_fit1sub", "h_mm_pisub", "h_mm_full"):
            if target_hists.get(key) is not None:
                target_hists[key].Fill(adj_MM, weight)
    if not allcuts:
        return
    phi_shift = float(evt.ph_q)
    fill_specs = (
        ("hgcer_xy", (evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, weight * evt.P_hgcer_npeSum)),
        ("hgcer_x_mm", (evt.P_hgcer_xAtCer, adj_MM, weight * evt.P_hgcer_npeSum)),
        ("hgcer_y_mm", (evt.P_hgcer_yAtCer, adj_MM, weight * evt.P_hgcer_npeSum)),
        ("mm_ct", (adj_MM, evt.CTime_ROC1, weight)),
        ("ct_beta", (evt.CTime_ROC1, evt.P_gtr_beta, weight)),
        ("mm_beta", (adj_MM, evt.P_gtr_beta, weight)),
        ("mm_h_cer", (adj_MM, evt.H_cer_npeSum, weight)),
        ("mm_h_cal", (adj_MM, evt.H_cal_etottracknorm, weight)),
        ("mm_p_cal", (adj_MM, evt.P_cal_etottracknorm, weight)),
        ("mm_p_hgcer", (adj_MM, evt.P_hgcer_npeSum, weight)),
        ("mm_p_aero", (adj_MM, evt.P_aero_npeSum, weight)),
        ("phiq_t", (phi_shift, adj_t, weight)),
        ("q2_w", (evt.Q2, evt.W, weight)),
        ("q2_t", (evt.Q2, adj_t, weight)),
        ("w_t", (evt.W, adj_t, weight)),
        ("eps_t", (evt.epsilon, adj_t, weight)),
        ("mm_t", (adj_MM, adj_t, weight)),
        ("h_ct", (evt.CTime_ROC1, weight)),
        ("h_ssxfp", (evt.ssxfp, weight)),
        ("h_ssyfp", (evt.ssyfp, weight)),
        ("h_ssxpfp", (evt.ssxpfp, weight)),
        ("h_ssypfp", (evt.ssypfp, weight)),
        ("h_hsxfp", (evt.hsxfp, weight)),
        ("h_hsyfp", (evt.hsyfp, weight)),
        ("h_hsxpfp", (evt.hsxpfp, weight)),
        ("h_hsypfp", (evt.hsypfp, weight)),
        ("h_ssxptar", (evt.ssxptar, weight)),
        ("h_ssyptar", (evt.ssyptar, weight)),
        ("h_hsxptar", (evt.hsxptar, weight)),
        ("h_hsyptar", (evt.hsyptar, weight)),
        ("h_ssdelta", (evt.ssdelta, weight)),
        ("h_hsdelta", (adj_hsdelta, weight)),
        ("h_ph_q", (phi_shift, weight)),
        ("h_th_q", (evt.th_q, weight)),
        ("h_ph_recoil", (evt.ph_recoil, weight)),
        ("h_th_recoil", (evt.th_recoil, weight)),
        ("h_q2", (evt.Q2, weight)),
        ("h_t", (adj_t, weight)),
        ("h_w", (evt.W, weight)),
        ("h_epsilon", (evt.epsilon, weight)),
        ("h_mm", (adj_MM, weight)),
        ("h_pmiss", (evt.pmiss, weight)),
        ("h_emiss", (evt.emiss, weight)),
        ("h_pmx", (evt.pmx, weight)),
        ("h_pmy", (evt.pmy, weight)),
        ("h_pmz", (evt.pmz, weight)),
        ("h_cal", (evt.H_cal_etottracknorm, weight)),
        ("h_cer", (evt.H_cer_npeSum, weight)),
        ("p_cal", (evt.P_cal_etottracknorm, weight)),
        ("p_hgcer", (evt.P_hgcer_npeSum, weight)),
        ("p_aero", (evt.P_aero_npeSum, weight)),
    )
    for key, args in fill_specs:
        hist = target_hists.get(key)
        if hist is None:
            continue
        hist.Fill(*args)


def _clone_target_map(target_templates, suffix):
    cloned = {}
    for key, hist in (target_templates or {}).items():
        cloned[key] = _clone_hist(hist, "{}{}".format(hist.GetName(), suffix), reset=True) if hist is not None else None
    return cloned


def apply_kaon_proton_cleaning_to_targets(
    cleaning_result,
    source_bundle,
    target_templates,
    evaluate_event,
    shifted_mm_getter,
    shifted_t_getter,
    hole_contains,
    mm_min,
    mm_max,
):
    if not isinstance(cleaning_result, dict):
        return {"accepted": False, "reason": "missing cleaning result"}
    if not bool(cleaning_result.get("accepted")):
        return {"accepted": False, "reason": cleaning_result.get("fallback_reason") or "cleaning result rejected"}

    config = cleaning_result.get("settings") or {}
    aero_edges = [float(edge) for edge in (cleaning_result.get("aero_edges") or config.get("aero_slice_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0))]
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    timing_branch = str(cleaning_result.get("selected_timing_branch") or "CTime_ROC1")
    denominator_floor = float(((config.get("weighting") or {}).get("denominator_floor", 1.0e-12)))
    raw_targets = _clone_target_map(target_templates, "_proton_clean_raw")
    proton_targets = _clone_target_map(target_templates, "_proton_clean_proton")
    cleaned_targets_pre_rf = _clone_target_map(target_templates, "_proton_clean_cleaned_prerf")
    final_targets = _clone_target_map(target_templates, "_proton_clean_final")

    h_weight_sum_delta = ROOT.TH1D(
        "H_proton_weight_sum_delta",
        "Average proton weight vs #delta;SHMS #delta [%];Average proton weight",
        len(delta_edges) - 1,
        float(delta_edges[0]),
        float(delta_edges[-1]),
    )
    h_weight_sum_delta.SetDirectory(0)
    h_weight_sum_delta.Sumw2()
    h_weight_norm_delta = _clone_hist(h_weight_sum_delta, "H_proton_weight_norm_delta", reset=True)
    h_weight_avg_delta = _clone_hist(h_weight_sum_delta, "H_proton_weight_avg_delta", reset=True)

    h_weight_sum_delta_aero = ROOT.TH2D(
        "H_proton_weight_sum_delta_aero",
        "Average proton weight vs #delta and aero;SHMS #delta [%];P_aero NPE",
        len(delta_edges) - 1,
        float(delta_edges[0]),
        float(delta_edges[-1]),
        len(aero_edges) - 1,
        array("d", aero_edges),
    )
    h_weight_sum_delta_aero.SetDirectory(0)
    h_weight_sum_delta_aero.Sumw2()
    h_weight_norm_delta_aero = _clone_hist(h_weight_sum_delta_aero, "H_proton_weight_norm_delta_aero", reset=True)
    h_weight_avg_delta_aero = _clone_hist(h_weight_sum_delta_aero, "H_proton_weight_avg_delta_aero", reset=True)

    h_mm_raw = raw_targets.get("h_mm_nosub")
    h_mm_proton = proton_targets.get("h_mm_nosub")
    h_mm_cleaned = final_targets.get("h_mm_nosub")

    rf_counts = {"accepted": 0, "rejected": 0}
    support_counts = {SUPPORT_SUPPORTED: 0, SUPPORT_MARGINAL: 0, SUPPORT_UNSUPPORTED: 0}

    for source_name, source_spec in (source_bundle.get("sources") or {}).items():
        tree = source_spec.get("tree")
        coefficient = float(source_spec.get("coefficient", 0.0) or 0.0)
        if tree is None or coefficient == 0.0:
            continue
        for evt in tree:
            allcuts, nommcuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) if hole_contains is not None else False
            base_allcuts = bool(allcuts)
            allcuts = bool(allcuts and (not hole_rejected))
            nommcuts = bool(nommcuts and (not hole_rejected))
            noholecuts = bool(base_allcuts)
            if not (allcuts or nommcuts or noholecuts):
                continue
            adj_mm = float(shifted_mm_getter(evt))
            adj_t = float(shifted_t_getter(evt))
            delta_value = float(getattr(evt, "ssdelta", 0.0))
            aero_value = float(getattr(evt, "P_aero_npeSum", 0.0))
            try:
                timing_value = float(getattr(evt, timing_branch))
            except Exception:
                timing_value = 0.0
            delta_index = _find_collection_bin(delta_value, delta_edges)
            aero_index = _find_collection_bin(aero_value, aero_edges)
            support_label = SUPPORT_UNSUPPORTED
            proton_weight = 0.0
            if 0 <= delta_index < len(cleaning_result.get("support_by_delta") or []):
                support_label = str(cleaning_result["support_by_delta"][delta_index])
                if 0 <= aero_index < len(cleaning_result.get("global_shapes") or []):
                    global_shape = cleaning_result["global_shapes"][aero_index]
                    slice_collection = (cleaning_result.get("delta_slice_fits") or [])[delta_index]
                    if 0 <= aero_index < len(slice_collection):
                        proton_weight = _evaluate_event_proton_probability(
                            timing_value,
                            global_shape,
                            slice_collection[aero_index],
                            denominator_floor,
                        )
            support_counts[support_label] = support_counts.get(support_label, 0) + 1
            raw_weight = float(coefficient)
            proton_component_weight = float(coefficient) * float(proton_weight)
            cleaned_weight = float(coefficient) * float(max(0.0, min(1.0, 1.0 - proton_weight)))
            _fill_standard_target_hists(
                evt,
                adj_mm,
                adj_t,
                adj_hsdelta,
                raw_weight,
                raw_targets,
                allcuts=allcuts,
                nommcuts=nommcuts,
                noholecuts=noholecuts,
            )
            _fill_standard_target_hists(
                evt,
                adj_mm,
                adj_t,
                adj_hsdelta,
                proton_component_weight,
                proton_targets,
                allcuts=allcuts,
                nommcuts=nommcuts,
                noholecuts=noholecuts,
            )
            _fill_standard_target_hists(
                evt,
                adj_mm,
                adj_t,
                adj_hsdelta,
                cleaned_weight,
                cleaned_targets_pre_rf,
                allcuts=allcuts,
                nommcuts=nommcuts,
                noholecuts=noholecuts,
            )
            rf_accept = apply_low_epsilon_rf_after_proton_cleaning(
                cleaning_result,
                source_name,
                evt,
            )
            if rf_accept:
                rf_counts["accepted"] += 1
                _fill_standard_target_hists(
                    evt,
                    adj_mm,
                    adj_t,
                    adj_hsdelta,
                    cleaned_weight,
                    final_targets,
                    allcuts=allcuts,
                    nommcuts=nommcuts,
                    noholecuts=noholecuts,
                )
            else:
                rf_counts["rejected"] += 1
            if 0 <= delta_index < (len(delta_edges) - 1):
                abs_norm = abs(float(coefficient))
                h_weight_sum_delta.Fill(delta_value, abs_norm * proton_weight)
                h_weight_norm_delta.Fill(delta_value, abs_norm)
                if 0 <= aero_index < (len(aero_edges) - 1):
                    h_weight_sum_delta_aero.Fill(delta_value, aero_value, abs_norm * proton_weight)
                    h_weight_norm_delta_aero.Fill(delta_value, aero_value, abs_norm)

    for bin_index in range(1, h_weight_avg_delta.GetNbinsX() + 1):
        denominator = float(h_weight_norm_delta.GetBinContent(bin_index))
        numerator = float(h_weight_sum_delta.GetBinContent(bin_index))
        if denominator > 0.0:
            h_weight_avg_delta.SetBinContent(bin_index, numerator / denominator)
    for x_bin in range(1, h_weight_avg_delta_aero.GetNbinsX() + 1):
        for y_bin in range(1, h_weight_avg_delta_aero.GetNbinsY() + 1):
            denominator = float(h_weight_norm_delta_aero.GetBinContent(x_bin, y_bin))
            numerator = float(h_weight_sum_delta_aero.GetBinContent(x_bin, y_bin))
            if denominator > 0.0:
                h_weight_avg_delta_aero.SetBinContent(x_bin, y_bin, numerator / denominator)

    h_proton_fraction_vs_mm = _clone_hist(h_mm_cleaned, "H_proton_fraction_vs_MM", reset=True) if h_mm_cleaned is not None else None
    if h_proton_fraction_vs_mm is not None and h_mm_raw is not None and h_mm_proton is not None:
        for bin_index in range(1, h_proton_fraction_vs_mm.GetNbinsX() + 1):
            raw_value = float(h_mm_raw.GetBinContent(bin_index))
            proton_value = float(h_mm_proton.GetBinContent(bin_index))
            if raw_value > 0.0:
                h_proton_fraction_vs_mm.SetBinContent(bin_index, max(0.0, min(1.0, proton_value / raw_value)))
            else:
                h_proton_fraction_vs_mm.SetBinContent(bin_index, 0.0)

    application = {
        "accepted": True,
        "raw_targets": raw_targets,
        "proton_targets": proton_targets,
        "cleaned_targets_pre_rf": cleaned_targets_pre_rf,
        "final_targets": final_targets,
        "H_MM_before_proton_cleaning": h_mm_raw,
        "H_MM_estimated_proton": h_mm_proton,
        "H_MM_after_proton_cleaning": h_mm_cleaned,
        "H_proton_fraction_vs_MM": h_proton_fraction_vs_mm,
        "H_proton_weight_vs_delta": h_weight_avg_delta,
        "H_proton_weight_vs_delta_aero": h_weight_avg_delta_aero,
        "rf_counts": rf_counts,
        "support_counts": support_counts,
        "diagnostics": {
            "rf_applied": bool((cleaning_result.get("diagnostics") or {}).get("rf_applied", False)),
            "rf_counts": rf_counts,
            "support_counts": support_counts,
            "raw_mm_integral": _hist_integral(h_mm_raw),
            "estimated_proton_mm_integral": _hist_integral(h_mm_proton),
            "cleaned_mm_integral": _hist_integral(h_mm_cleaned),
            "target_sample_definition": (
                "signed noRF kaon target bundle after proton cleaning, with RF "
                "membership optionally applied only after event-level proton weights"
            ),
            "selected_timing_branch": str(timing_branch),
        },
    }
    cleaning_result["application"] = application
    return application


def serialize_kaon_proton_cleaning_result(cleaning_result):
    if not isinstance(cleaning_result, dict):
        return {}
    payload = dict(cleaning_result)
    payload.pop("_rf_signature_lookup", None)
    return _json_ready_value(payload)


def summarize_kaon_proton_cleaning_result(cleaning_result):
    if not isinstance(cleaning_result, dict):
        return {}
    diagnostics = deepcopy(cleaning_result.get("diagnostics") or {})
    application = cleaning_result.get("application") or {}
    summary = {
        "accepted": bool(cleaning_result.get("accepted", False)),
        "enabled": bool(cleaning_result.get("enabled", False)),
        "method": cleaning_result.get("method"),
        "analysis_scope": cleaning_result.get("analysis_scope"),
        "context": cleaning_result.get("context"),
        "fallback_reason": cleaning_result.get("fallback_reason", ""),
        "diagnostics": diagnostics,
    }
    if isinstance(application, dict):
        summary["application"] = _json_ready_value(application.get("diagnostics") or {})
    return _json_ready_value(summary)


def print_kaon_proton_cleaning_terminal_summary(cleaning_result, output_pdf=None):
    if not isinstance(cleaning_result, dict):
        return
    diagnostics = cleaning_result.get("diagnostics") or {}
    application = cleaning_result.get("application") or {}
    source_stats = diagnostics.get("source_stats") or {}
    valid_global_shape_count = int(diagnostics.get("valid_global_shape_count", 0) or 0)
    n_aero_slices = max(len(cleaning_result.get("aero_edges") or []) - 1, 0)
    supported_delta_bins = int(diagnostics.get("supported_delta_bins", 0) or 0)
    marginal_delta_bins = int(diagnostics.get("marginal_delta_bins", 0) or 0)
    support_by_delta = cleaning_result.get("support_by_delta") or []
    unsupported_delta_bins = int(
        sum(1 for label in support_by_delta if str(label) == SUPPORT_UNSUPPORTED)
    )
    lines = [
        "\n===== Event-level proton-cleaning =====",
        "Scope: {}".format(cleaning_result.get("analysis_scope", "unknown")),
        "Context: {}".format(cleaning_result.get("context", "")),
        "Method: {}".format(cleaning_result.get("method", "unknown")),
        "Input tree state: {}".format(diagnostics.get("input_tree_state", "unknown")),
        "Particle selection: {}".format(diagnostics.get("input_tree_particle_selection", "unknown")),
        "Tree policy: {}".format(diagnostics.get("tree_policy", "unknown")),
        "RF policy: {}".format(diagnostics.get("rf_policy", "unknown")),
        "Global PID source usage:",
    ]
    for source_name in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
        source_payload = source_stats.get(source_name) or {}
        lines.append(
            "  {}: seen={} used={} tree={}".format(
                source_name,
                int(source_payload.get("entries_seen", 0) or 0),
                int(source_payload.get("entries_used", 0) or 0),
                source_payload.get("tree_name", "missing"),
            )
        )
    lines.extend(
        [
            "Identifiable global aerogel slices: {} / {}".format(
                valid_global_shape_count,
                n_aero_slices,
            ),
            "Selected timing branch: {}".format(diagnostics.get("selected_timing_branch", "unknown")),
            "Selected probe kind/mode: {}/{}".format(
                diagnostics.get("selected_probe_kind", "unknown"),
                diagnostics.get("selected_probe_fit_mode", "unknown"),
            ),
            "Timing selection reason: {}".format(diagnostics.get("timingSelectionReason", "unknown")),
            "Supported delta bins: {}".format(supported_delta_bins),
            "Marginal delta bins: {}".format(marginal_delta_bins),
            "Unsupported delta bins: {}".format(unsupported_delta_bins),
        ]
    )
    if bool(cleaning_result.get("accepted")):
        lines.append("Proton-cleaning fit status: accepted")
    else:
        lines.append("Proton-cleaning fit status: rejected")
        lines.append("Failure reason: {}".format(cleaning_result.get("fallback_reason", "unknown")))
        if output_pdf:
            lines.append("Writing diagnostic proton-cleaning plots to: {}".format(output_pdf))
    if isinstance(application, dict) and bool(application.get("accepted")):
        app_diag = application.get("diagnostics") or {}
        rf_counts = app_diag.get("rf_counts") or application.get("rf_counts") or {}
        support_counts = app_diag.get("support_counts") or application.get("support_counts") or {}
        lines.extend(
            [
                "Events in supported delta bins: {}".format(int(support_counts.get(SUPPORT_SUPPORTED, 0) or 0)),
                "Events in marginal delta bins: {}".format(int(support_counts.get(SUPPORT_MARGINAL, 0) or 0)),
                "Events in unsupported delta bins: {}".format(int(support_counts.get(SUPPORT_UNSUPPORTED, 0) or 0)),
                "RF accepted/rejected after cleaning: {}/{}".format(
                    int(rf_counts.get("accepted", 0) or 0),
                    int(rf_counts.get("rejected", 0) or 0),
                ),
                "Raw MM integral: {:.3e}".format(float(app_diag.get("raw_mm_integral", 0.0) or 0.0)),
                "Estimated proton MM integral: {:.3e}".format(float(app_diag.get("estimated_proton_mm_integral", 0.0) or 0.0)),
                "Proton-cleaned MM integral: {:.3e}".format(float(app_diag.get("cleaned_mm_integral", 0.0) or 0.0)),
            ]
        )
    if output_pdf and bool(cleaning_result.get("accepted")):
        lines.append("Diagnostic proton-cleaning plots: {}".format(output_pdf))
    lines.append("========================================")
    print("\n".join(lines))


def _draw_hist_overlay(canvas_name, title, histograms, legend_entries, output_pdf):
    canvas = TCanvas(canvas_name, title, 1000, 700)
    canvas.SetGrid()
    maximum = 0.0
    minimum = 0.0
    for hist in histograms:
        if hist is None:
            continue
        maximum = max(maximum, float(hist.GetMaximum()))
        minimum = min(minimum, float(hist.GetMinimum()))
    first_drawn = False
    legend = TLegend(0.62, 0.68, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for hist, style in histograms.items():
        if hist is None:
            continue
        hist.SetLineColor(style["color"])
        hist.SetLineStyle(style.get("style", 1))
        hist.SetLineWidth(style.get("width", 2))
        hist.SetMinimum(minimum * 1.10 if minimum < 0.0 else 0.0)
        hist.SetMaximum(maximum * 1.15 if maximum > 0.0 else 1.0)
        hist.Draw("hist" if not first_drawn else "hist same")
        legend.AddEntry(hist, legend_entries.get(hist.GetName(), hist.GetTitle()), "l")
        first_drawn = True
    legend.Draw()
    canvas.Print(output_pdf)


def _draw_status_pave(cleaning_result, extra_lines=None, x1=0.12, y1=0.74, x2=0.52, y2=0.90):
    pave = TPaveText(x1, y1, x2, y2, "NDC")
    pave.SetBorderSize(1)
    pave.SetFillStyle(0)
    pave.SetTextAlign(12)
    pave.SetTextSize(0.028)
    accepted = bool((cleaning_result or {}).get("accepted"))
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    fit_status = "accepted" if accepted else "rejected"
    pave.AddText("status: {}".format(fit_status))
    analysis_scope = str((cleaning_result or {}).get("analysis_scope") or diagnostics.get("analysis_scope") or "unknown")
    pave.AddText("scope: {}".format(analysis_scope))
    fallback_reason = str((cleaning_result or {}).get("fallback_reason") or "").strip()
    if fallback_reason:
        pave.AddText("reason: {}".format(fallback_reason[:120]))
    for line in extra_lines or ():
        if line:
            pave.AddText(str(line)[:140])
    pave.Draw()
    return pave


def _build_timing_shape_overlays(hist, global_shape, slice_fit):
    if hist is None or not global_shape or not slice_fit:
        return None
    if not global_shape.get("valid") or not slice_fit.get("valid"):
        return None
    x_values = np.asarray(
        [hist.GetXaxis().GetBinCenter(bin_index) for bin_index in range(1, hist.GetNbinsX() + 1)],
        dtype=float,
    )
    kaon_hist = _clone_hist(hist, "{}_k_overlay".format(hist.GetName()), reset=True)
    proton_hist = _clone_hist(hist, "{}_p_overlay".format(hist.GetName()), reset=True)
    total_hist = _clone_hist(hist, "{}_tot_overlay".format(hist.GetName()), reset=True)
    if kaon_hist is None or proton_hist is None or total_hist is None:
        return None
    for bin_index, x_value in enumerate(x_values, start=1):
        kaon_value = float(
            _gaussian(
                np.asarray([x_value]),
                slice_fit["kaon_amplitude"],
                global_shape["kaon_mean"],
                global_shape["kaon_sigma"],
            )[0]
        )
        proton_value = float(
            _gaussian(
                np.asarray([x_value]),
                slice_fit["proton_amplitude"],
                global_shape["proton_mean"],
                global_shape["proton_sigma"],
            )[0]
        )
        total_value = kaon_value + proton_value + float(slice_fit.get("other_amplitude", 0.0) or 0.0)
        kaon_hist.SetBinContent(bin_index, kaon_value)
        proton_hist.SetBinContent(bin_index, proton_value)
        total_hist.SetBinContent(bin_index, total_value)
    kaon_hist.SetLineColor(kBlue)
    proton_hist.SetLineColor(kRed)
    total_hist.SetLineColor(kGreen + 2)
    return {
        "kaon": kaon_hist,
        "proton": proton_hist,
        "total": total_hist,
    }


def print_kaon_proton_cleaning_pages(output_pdf, cleaning_result, title_prefix=""):
    if not isinstance(cleaning_result, dict):
        return
    if not bool(cleaning_result.get("enabled")):
        return
    gStyle.SetOptStat(0)
    prefix = str(title_prefix).strip() or "Kaon"
    application = cleaning_result.get("application") or {}
    diagnostics = cleaning_result.get("diagnostics") or {}
    support_by_delta = cleaning_result.get("support_by_delta") or []

    h_global_pid = cleaning_result.get("H_global_pid")
    if h_global_pid is not None:
        canvas = TCanvas("C_proton_cleaning_global_pid", "{} proton-cleaning global PID".format(prefix), 1000, 700)
        h_global_pid.Draw("colz")
        _draw_status_pave(cleaning_result)
        canvas.Print(output_pdf)

    global_slice_hists = cleaning_result.get("H_global_timing_slices") or []
    if global_slice_hists:
        canvas = TCanvas("C_proton_cleaning_global_slices", "{} proton-cleaning global timing slices".format(prefix), 1200, 800)
        canvas.Divide(3, 2)
        global_shapes = cleaning_result.get("global_shapes") or []
        for index, hist in enumerate(global_slice_hists[:6]):
            canvas.cd(index + 1)
            hist.SetLineColor(kBlack)
            hist.Draw("hist")
            if index < len(global_shapes):
                shape = global_shapes[index]
                overlays = _build_timing_shape_overlays(
                    hist,
                    shape,
                    {
                        "valid": shape.get("valid"),
                        "kaon_amplitude": shape.get("kaon_amplitude", 0.0),
                        "proton_amplitude": shape.get("proton_amplitude", 0.0),
                        "other_amplitude": shape.get("other_amplitude", 0.0),
                    },
                )
                if overlays is not None:
                    overlays["kaon"].Draw("hist same")
                    overlays["proton"].Draw("hist same")
                    overlays["total"].Draw("hist same")
                _draw_status_pave(
                    cleaning_result,
                    extra_lines=(
                        "slice {} valid={}".format(index + 1, bool(shape.get("valid"))),
                    ),
                    x1=0.14,
                    y1=0.74,
                    x2=0.58,
                    y2=0.88,
                )
        canvas.Print(output_pdf)

    h_support = ROOT.TH1D(
        "H_proton_cleaning_support_class",
        "{} delta support summary;SHMS #delta bin;support class".format(prefix),
        max(len(support_by_delta), 1),
        0.0,
        float(max(len(support_by_delta), 1)),
    )
    h_support.SetDirectory(0)
    h_support.Sumw2()
    for bin_index, label in enumerate(support_by_delta, start=1):
        h_support.SetBinContent(bin_index, SUPPORT_CLASS_TO_CODE.get(str(label), 0.0))

    h_proton = ROOT.TH1D(
        "H_proton_cleaning_proton_yield_delta",
        "{} fitted yields vs #delta;SHMS #delta bin;Yield".format(prefix),
        max(len((cleaning_result.get("diagnostics") or {}).get("proton_yield_by_delta") or []), 1),
        0.0,
        float(max(len((cleaning_result.get("diagnostics") or {}).get("proton_yield_by_delta") or []), 1)),
    )
    h_proton.SetDirectory(0)
    h_proton.Sumw2()
    h_kaon = _clone_hist(h_proton, "H_proton_cleaning_kaon_yield_delta", reset=True)
    h_other = _clone_hist(h_proton, "H_proton_cleaning_other_yield_delta", reset=True)
    for bin_index, value in enumerate((cleaning_result.get("diagnostics") or {}).get("proton_yield_by_delta") or [], start=1):
        h_proton.SetBinContent(bin_index, float(value))
    for bin_index, value in enumerate((cleaning_result.get("diagnostics") or {}).get("kaon_yield_by_delta") or [], start=1):
        h_kaon.SetBinContent(bin_index, float(value))
    for bin_index, value in enumerate((cleaning_result.get("diagnostics") or {}).get("other_yield_by_delta") or [], start=1):
        h_other.SetBinContent(bin_index, float(value))
    h_proton.SetLineColor(kRed)
    h_kaon.SetLineColor(kBlue)
    h_other.SetLineColor(kGray + 2)
    canvas = TCanvas("C_proton_cleaning_delta_summary", "{} proton-cleaning delta summary".format(prefix), 1200, 800)
    canvas.Divide(2, 2)
    canvas.cd(1)
    h_support.SetMinimum(-0.2)
    h_support.SetMaximum(2.2)
    h_support.Draw("hist")
    canvas.cd(2)
    h_proton.Draw("hist")
    h_kaon.Draw("hist same")
    h_other.Draw("hist same")
    legend = TLegend(0.60, 0.68, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_proton, "proton yield", "l")
    legend.AddEntry(h_kaon, "kaon yield", "l")
    legend.AddEntry(h_other, "other yield", "l")
    legend.Draw()
    canvas.cd(3)
    if application.get("H_proton_weight_vs_delta") is not None:
        application["H_proton_weight_vs_delta"].SetLineColor(kViolet + 1)
        application["H_proton_weight_vs_delta"].Draw("hist")
    canvas.cd(4)
    if application.get("H_proton_weight_vs_delta_aero") is not None:
        application["H_proton_weight_vs_delta_aero"].Draw("colz")
    else:
        first_delta_pid = next((hist for hist in (cleaning_result.get("H_delta_pid") or []) if hist is not None), None)
        if first_delta_pid is not None:
            first_delta_pid.Draw("colz")
    canvas.cd(1)
    _draw_status_pave(
        cleaning_result,
        extra_lines=(
            "supported={} marginal={}".format(
                diagnostics.get("supported_delta_bins", 0),
                diagnostics.get("marginal_delta_bins", 0),
            ),
        ),
        x1=0.14,
        y1=0.74,
        x2=0.54,
        y2=0.88,
    )
    canvas.Print(output_pdf)

    delta_pid_hists = cleaning_result.get("H_delta_pid") or []
    delta_slice_hists = cleaning_result.get("H_delta_timing_slices") or []
    delta_slice_fits = cleaning_result.get("delta_slice_fits") or []
    global_shapes = cleaning_result.get("global_shapes") or []
    for delta_index, pid_hist in enumerate(delta_pid_hists):
        slice_collection = delta_slice_hists[delta_index] if delta_index < len(delta_slice_hists) else []
        slice_fit_collection = delta_slice_fits[delta_index] if delta_index < len(delta_slice_fits) else []
        canvas = TCanvas(
            "C_proton_cleaning_delta_{}".format(delta_index),
            "{} proton-cleaning delta bin {}".format(prefix, delta_index + 1),
            1400,
            900,
        )
        canvas.Divide(3, 2)
        canvas.cd(1)
        if pid_hist is not None:
            pid_hist.Draw("colz")
        support_label = support_by_delta[delta_index] if delta_index < len(support_by_delta) else "unknown"
        _draw_status_pave(
            cleaning_result,
            extra_lines=(
                "delta bin {}".format(delta_index + 1),
                "support={}".format(support_label),
            ),
            x1=0.14,
            y1=0.72,
            x2=0.54,
            y2=0.88,
        )
        for slice_index, hist in enumerate(slice_collection[:5], start=2):
            canvas.cd(slice_index)
            if hist is None:
                continue
            hist.SetLineColor(kBlack)
            hist.Draw("hist")
            global_shape = global_shapes[slice_index - 2] if (slice_index - 2) < len(global_shapes) else {}
            slice_fit = slice_fit_collection[slice_index - 2] if (slice_index - 2) < len(slice_fit_collection) else {}
            overlays = _build_timing_shape_overlays(hist, global_shape, slice_fit)
            if overlays is not None:
                overlays["kaon"].Draw("hist same")
                overlays["proton"].Draw("hist same")
                overlays["total"].Draw("hist same")
            _draw_status_pave(
                cleaning_result,
                extra_lines=(
                    "aero slice {}".format(slice_index - 1),
                    "valid={}".format(bool(slice_fit.get("valid"))),
                ),
                x1=0.16,
                y1=0.74,
                x2=0.56,
                y2=0.88,
            )
        canvas.Print(output_pdf)

    h_mm_before = application.get("H_MM_before_proton_cleaning")
    h_mm_proton = application.get("H_MM_estimated_proton")
    h_mm_after = application.get("H_MM_after_proton_cleaning")
    if h_mm_before is not None and h_mm_proton is not None and h_mm_after is not None:
        canvas = TCanvas("C_proton_cleaning_mm", "{} proton-cleaning MM validation".format(prefix), 1000, 900)
        canvas.Divide(1, 2)
        canvas.cd(1)
        h_mm_before.SetLineColor(kBlack)
        h_mm_proton.SetLineColor(kRed)
        h_mm_after.SetLineColor(kGreen + 2)
        h_mm_before.Draw("hist")
        h_mm_proton.Draw("hist same")
        h_mm_after.Draw("hist same")
        legend = TLegend(0.60, 0.68, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(h_mm_before, "raw kaon noRF", "l")
        legend.AddEntry(h_mm_proton, "estimated proton", "l")
        legend.AddEntry(h_mm_after, "proton-cleaned", "l")
        legend.Draw()
        canvas.cd(2)
        if application.get("H_proton_fraction_vs_MM") is not None:
            application["H_proton_fraction_vs_MM"].SetLineColor(kViolet + 1)
            application["H_proton_fraction_vs_MM"].SetMinimum(0.0)
            application["H_proton_fraction_vs_MM"].SetMaximum(1.0)
            application["H_proton_fraction_vs_MM"].Draw("hist")
        canvas.Print(output_pdf)

        validation_windows = ((cleaning_result.get("settings") or {}).get("validation_windows") or {})
        for window_name in ("low_mm", "lambda_peak"):
            bounds = validation_windows.get(window_name)
            if not bounds:
                continue
            canvas = TCanvas(
                "C_proton_cleaning_mm_zoom_{}".format(window_name),
                "{} proton-cleaning {} validation".format(prefix, window_name),
                1000,
                700,
            )
            h_mm_before.GetXaxis().SetRangeUser(float(bounds[0]), float(bounds[1]))
            h_mm_proton.GetXaxis().SetRangeUser(float(bounds[0]), float(bounds[1]))
            h_mm_after.GetXaxis().SetRangeUser(float(bounds[0]), float(bounds[1]))
            h_mm_before.Draw("hist")
            h_mm_proton.Draw("hist same")
            h_mm_after.Draw("hist same")
            legend = TLegend(0.60, 0.68, 0.88, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(h_mm_before, "raw kaon noRF", "l")
            legend.AddEntry(h_mm_proton, "estimated proton", "l")
            legend.AddEntry(h_mm_after, "proton-cleaned", "l")
            legend.Draw()
            canvas.Print(output_pdf)
            h_mm_before.GetXaxis().UnZoom()
            h_mm_proton.GetXaxis().UnZoom()
            h_mm_after.GetXaxis().UnZoom()
