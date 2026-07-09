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


def _fit_global_timing_shape(
    histogram,
    config,
    function_name,
):
    fit_min, fit_max = config["ctime_hist_range"]
    inputs = _extract_weighted_hist_fit_inputs(histogram, fit_min, fit_max)
    x_values = inputs["x"]
    y_values = inputs["y"]
    sigma_values = inputs["sigma"]
    minimum_entries = int((config.get("global_fit") or {}).get("minimum_entries", 200))
    if x_values.size == 0 or np.sum(np.abs(y_values)) < float(minimum_entries):
        return {
            "valid": False,
            "fit_status": "insufficient_support",
            "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
            "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
        }

    kaon_mean_min, kaon_mean_max = (config.get("global_fit") or {}).get(
        "kaon_mean_range",
        (-0.45, 0.20),
    )
    proton_mean_min, proton_mean_max = (config.get("global_fit") or {}).get(
        "proton_mean_range",
        (0.20, 0.95),
    )
    sigma_min, sigma_max = (config.get("global_fit") or {}).get(
        "sigma_range",
        (0.03, 0.45),
    )
    histogram_maximum = max(float(np.max(y_values)), 1.0)
    kaon_seed = _find_peak_seed(histogram, kaon_mean_min, kaon_mean_max)
    proton_seed = _find_peak_seed(histogram, proton_mean_min, proton_mean_max)

    initial = np.asarray(
        [
            max(float(kaon_seed[0]), 0.20 * histogram_maximum),
            float(kaon_seed[1]),
            0.15,
            max(float(proton_seed[0]), 0.10 * histogram_maximum),
            float(proton_seed[1]),
            0.15,
            0.02 * histogram_maximum,
        ],
        dtype=float,
    )
    lower_bounds = np.asarray(
        [
            0.0,
            float(kaon_mean_min),
            float(sigma_min),
            0.0,
            float(proton_mean_min),
            float(sigma_min),
            0.0,
        ],
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
        model = kaon_model + proton_model + float(parameters[6])
        return (y_values - model) / sigma_values

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
    kaon_amp_err = uncertainties.get("kaon_amplitude", 0.0)
    proton_amp_err = uncertainties.get("proton_amplitude", 0.0)
    kaon_sigma = abs(float(parameter_vector[2]))
    proton_sigma = abs(float(parameter_vector[5]))
    separation_denominator = math.sqrt(max((kaon_sigma ** 2) + (proton_sigma ** 2), 0.0))
    separation = (
        abs(float(parameter_vector[4]) - float(parameter_vector[1])) / separation_denominator
        if separation_denominator > 0.0
        else 0.0
    )
    bound_fraction_tolerance = float(
        (config.get("global_fit") or {}).get("bound_fraction_tolerance", 0.02)
    )

    def near_bound(value, lower, upper):
        tolerance = bound_fraction_tolerance * float(upper - lower)
        return (value <= float(lower) + tolerance) or (value >= float(upper) - tolerance)

    bound_hit = (
        near_bound(parameter_vector[1], kaon_mean_min, kaon_mean_max)
        or near_bound(parameter_vector[4], proton_mean_min, proton_mean_max)
        or near_bound(kaon_sigma, sigma_min, sigma_max)
        or near_bound(proton_sigma, sigma_min, sigma_max)
    )
    kaon_significance = float(parameter_vector[0] / kaon_amp_err) if kaon_amp_err > 0.0 else 0.0
    proton_significance = float(parameter_vector[3] / proton_amp_err) if proton_amp_err > 0.0 else 0.0
    valid = (
        bool(getattr(fit_result, "success", False))
        and not bound_hit
        and math.isfinite(chi2_ndf)
        and float(parameter_vector[0]) > 0.0
        and float(parameter_vector[3]) > 0.0
        and kaon_sigma > 0.0
        and proton_sigma > 0.0
        and float(parameter_vector[1]) < float(parameter_vector[4])
        and separation >= float((config.get("global_fit") or {}).get("minimum_separation", 0.75))
        and kaon_significance >= float((config.get("global_fit") or {}).get("minimum_amplitude_significance", 2.0))
        and proton_significance >= float((config.get("global_fit") or {}).get("minimum_amplitude_significance", 2.0))
        and chi2_ndf <= float((config.get("global_fit") or {}).get("maximum_chi2_ndf", 5.0))
    )
    return {
        "valid": bool(valid),
        "fit_status": "success" if bool(getattr(fit_result, "success", False)) else "failure",
        "message": str(getattr(fit_result, "message", "")),
        "function_name": str(function_name),
        "kaon_amplitude": float(parameter_vector[0]),
        "kaon_amplitude_error": float(kaon_amp_err),
        "kaon_mean": float(parameter_vector[1]),
        "kaon_sigma": float(kaon_sigma),
        "proton_amplitude": float(parameter_vector[3]),
        "proton_amplitude_error": float(proton_amp_err),
        "proton_mean": float(parameter_vector[4]),
        "proton_sigma": float(proton_sigma),
        "other_amplitude": float(parameter_vector[6]),
        "separation": float(separation),
        "kaon_significance": float(kaon_significance),
        "proton_significance": float(proton_significance),
        "chi2_data": chi2_data,
        "chi2_ndf": chi2_ndf,
        "active_bin_count": active_bin_count,
        "excluded_invalid_variance_bins": inputs["excluded_invalid_variance_bins"],
        "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
        "bound_hit": bool(bound_hit),
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "uncertainties": uncertainties,
    }


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
    fit_min, fit_max = config["ctime_hist_range"]
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
):
    aero_edges = [float(edge) for edge in (config.get("aero_slice_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0))]
    aero_min, aero_max = [float(value) for value in (config.get("aero_hist_range") or (0.0, 25.0))]
    time_min, time_max = [float(value) for value in (config.get("ctime_hist_range") or (-1.50, 1.25))]
    delta_min, delta_max = [float(value) for value in (config.get("delta_hist_range") or (-10.0, 20.0))]
    delta_bins = int(config.get("delta_bins", 10) or 10)

    h_global_pid = ROOT.TH2D(
        "H_proton_cleaning_global_pid",
        "Global CTime_ROC1 vs P_aero NPE;P_aero NPE;CTime_ROC1 [ns];Signed yield",
        75,
        aero_min,
        aero_max,
        90,
        time_min,
        time_max,
    )
    h_global_pid.SetDirectory(0)
    h_global_pid.Sumw2()

    global_slice_hists = []
    for aero_index in range(len(aero_edges) - 1):
        hist = ROOT.TH1D(
            "H_proton_cleaning_global_time_slice_{}".format(aero_index),
            "Global timing slice {};CTime_ROC1 [ns];Signed yield".format(aero_index + 1),
            90,
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
            "PID slice {} ;P_aero NPE;CTime_ROC1 [ns];Signed yield".format(delta_index + 1),
            75,
            aero_min,
            aero_max,
            90,
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
                "Timing slice d{} a{};CTime_ROC1 [ns];Signed yield".format(delta_index + 1, aero_index + 1),
                90,
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
            timing_value = float(getattr(evt, "CTime_ROC1", 0.0))
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

    pid_payload = _build_signed_pid_histograms(
        source_bundle,
        config,
        evaluate_event,
        shifted_mm_getter,
        hole_contains,
        mm_min,
        mm_max,
    )
    result["diagnostics"]["source_stats"] = pid_payload["source_stats"]
    result["diagnostics"]["source_coefficients"] = {
        str(source_name): float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items()
    }
    result["H_global_pid"] = pid_payload["H_global_pid"]
    result["H_global_timing_slices"] = pid_payload["global_slice_hists"]
    result["H_delta_pid"] = pid_payload["delta_pid_hists"]
    result["H_delta_timing_slices"] = pid_payload["delta_slice_hists"]
    result["aero_edges"] = pid_payload["aero_edges"]
    result["delta_edges"] = pid_payload["delta_edges"]

    global_shapes = []
    valid_global_shape_count = 0
    for aero_index, slice_hist in enumerate(pid_payload["global_slice_hists"]):
        shape = _fit_global_timing_shape(
            slice_hist,
            config,
            "f_proton_cleaning_global_aero_{}".format(aero_index),
        )
        global_shapes.append(shape)
        if shape.get("valid"):
            valid_global_shape_count += 1
    result["global_shapes"] = global_shapes
    result["diagnostics"]["valid_global_shape_count"] = int(valid_global_shape_count)
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
            timing_value = float(getattr(evt, "CTime_ROC1", 0.0))
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
