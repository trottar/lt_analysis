#! /usr/bin/python

from __future__ import annotations

import gc
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
    kMagenta,
    kOrange,
    kRed,
    kViolet,
)
sys.path.append("utility")
from background_config import (  # noqa: E402
    PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT,
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

PROTON_CLEANING_EXACT_TIMING_BRANCH = "CTime_ROC1"
PROTON_CLEANING_EXACT_TIMING_RANGE = (-2.0, 2.0)
PROTON_CLEANING_EXACT_TIMING_BINS = 131
PROTON_CLEANING_EXACT_CT_LOW_EPSILON_RANGE = (-2.0, 2.0)
PROTON_CLEANING_EXACT_CT_HIGH_EPSILON_RANGE = (-4.0, 4.0)
PROTON_CLEANING_EXACT_CT_LOW_EPSILON_BINS = 131
PROTON_CLEANING_EXACT_CT_HIGH_EPSILON_BINS = 262
PROTON_CLEANING_EXACT_AERO_EDGES = (0.0, 3.0, 6.0, 10.0, 15.0, 25.0)
PROTON_CLEANING_EXACT_AERO_RANGE = (0.0, 25.0)
PROTON_CLEANING_EXACT_DELTA_RANGE = (-10.0, 20.0)
PROTON_CLEANING_EXACT_DELTA_BINS = 10
PROTON_CLEANING_EXACT_MM_VALIDATION_RANGE = (0.70, 1.50)
PROTON_CLEANING_EXACT_FIT_OPTIONS = "SRLIQ0"
PROTON_CLEANING_EXACT_GLOBAL_FIT = {
    "kaon_mean_range": (-0.45, 0.20),
    "proton_mean_range": (0.20, 0.95),
    "sigma_range": (0.03, 0.45),
    "initial_sigma": 0.15,
    "minimum_separation": 0.75,
    "minimum_amplitude_significance": 2.0,
    "maximum_poisson_deviance_ndf": 5.0,
    "maximum_poisson_deviance_per_entry": 0.85,
    "maximum_chi2_ndf": 5.0,
    "bound_fraction_tolerance": 0.02,
    "minimum_entries": 200,
}
PROTON_CLEANING_EXACT_SLICE_FIT = {
    "maximum_poisson_deviance_ndf": 5.0,
    "maximum_poisson_deviance_per_entry": 1.0,
    "maximum_chi2_ndf": 5.0,
    "minimum_model_data_ratio": 0.50,
    "maximum_model_data_ratio": 1.50,
    "minimum_entries": 30,
}
PROTON_CLEANING_EXACT_SUPPORT_THRESHOLDS = {
    "minimum_supported_slices": 2,
    "minimum_marginal_slices": 1,
    "minimum_supported_coverage": 0.35,
    "minimum_marginal_coverage": 0.15,
    "minimum_modeled_yield": 5.0,
}
PROTON_CLEANING_EXACT_WEIGHTING = {
    "denominator_floor": 1.0e-12,
}
PROTON_CLEANING_EXACT_VALIDATION_WINDOWS = {
    "low_mm": (0.80, 0.90),
    "lambda_peak": (1.105, 1.125),
}
SHMS_CENTRAL_PATH_CM = 1810.0
LIGHT_SPEED_CM_PER_NS = 29.9792
KAON_MASS_GEV = 0.493677
PROTON_MASS_GEV = 0.93827208
PROTON_CLEANING_TOF_OFFSET_RANGE = (-0.35, 0.35)
PROTON_CLEANING_TOF_OFFSET_VALIDATION = {
    "maximum_offset_error_ns": 0.10,
    "maximum_chi2_ndf": 5.0,
    "minimum_component_significance": 2.0,
    "minimum_smaller_component_fraction": 0.02,
    "reject_bound_hit_with_large_error": True,
    "bound_hit_large_error_fraction": 0.50,
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


def _hist_absolute_integral(hist, first_bin=None, last_bin=None):
    if hist is None:
        return 0.0
    try:
        start_bin = int(first_bin) if first_bin is not None else 1
        stop_bin = int(last_bin) if last_bin is not None else int(hist.GetNbinsX())
        total = 0.0
        for bin_index in range(max(1, start_bin), min(int(hist.GetNbinsX()), stop_bin) + 1):
            total += abs(float(hist.GetBinContent(bin_index)))
        return float(total)
    except Exception:
        return 0.0


def _build_hist_support_snapshot(histogram, support_entries, minimum_required_entries, first_fit_bin=None, last_fit_bin=None):
    signed_integral = _hist_integral(histogram)
    absolute_integral = _hist_absolute_integral(histogram)
    root_entries = 0.0
    effective_entries = 0.0
    fit_range_signed_integral = 0.0
    fit_range_absolute_integral = 0.0
    if histogram is not None:
        try:
            root_entries = float(histogram.GetEntries())
        except Exception:
            root_entries = 0.0
        try:
            effective_entries = float(histogram.GetEffectiveEntries())
        except Exception:
            effective_entries = 0.0
        if first_fit_bin is not None and last_fit_bin is not None:
            try:
                fit_range_signed_integral = float(histogram.Integral(int(first_fit_bin), int(last_fit_bin)))
            except Exception:
                fit_range_signed_integral = 0.0
            fit_range_absolute_integral = _hist_absolute_integral(
                histogram,
                first_bin=first_fit_bin,
                last_bin=last_fit_bin,
            )
    return {
        "support_entries": int(max(0, int(support_entries or 0))),
        "minimum_required_entries": int(max(0, int(minimum_required_entries or 0))),
        "signed_integral": float(signed_integral),
        "absolute_integral": float(absolute_integral),
        "root_entries": float(root_entries),
        "effective_entries": float(effective_entries),
        "fit_range_signed_integral": float(fit_range_signed_integral),
        "fit_range_absolute_integral": float(fit_range_absolute_integral),
    }


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


def _proton_debug_enabled(config):
    return bool((config or {}).get("debug_flares_enabled", False))


def _print_proton_debug(stage, **details):
    print("[DEBUG proton_cleaning] {}".format(stage), flush=True)
    for key, value in details.items():
        print("  {} = {}".format(key, value), flush=True)


def _format_debug_float(value, digits=4, scientific=False):
    if value is None:
        return "n/a"
    try:
        numeric = float(value)
    except Exception:
        return str(value)
    if not math.isfinite(numeric):
        return str(numeric)
    if scientific:
        return ("{0:." + str(int(digits)) + "e}").format(numeric)
    return ("{0:." + str(int(digits)) + "f}").format(numeric)


def _safe_event_float(evt, branch_name):
    try:
        value = float(getattr(evt, branch_name))
    except Exception:
        return None
    return value if math.isfinite(value) else None


def compute_shms_path_length_cm(p_shms_gev, ssxptar, ssyptar, ssdelta):
    del p_shms_gev, ssyptar  # Kept in signature for diagnostics and future path terms.
    path_length_cm = (
        SHMS_CENTRAL_PATH_CM
        + 0.11 * (1000.0 * float(ssxptar))
        + 0.057 * (float(ssdelta) / 100.0)
    )
    return float(path_length_cm)


def compute_proton_kaon_tof_separation_ns(p_shms_gev, path_length_cm):
    p_shms_gev = float(p_shms_gev)
    path_length_cm = float(path_length_cm)
    beta_kaon = p_shms_gev / math.sqrt((p_shms_gev ** 2) + (KAON_MASS_GEV ** 2))
    beta_proton = p_shms_gev / math.sqrt((p_shms_gev ** 2) + (PROTON_MASS_GEV ** 2))
    tof_kaon_ns = path_length_cm / (LIGHT_SPEED_CM_PER_NS * beta_kaon)
    tof_proton_ns = path_length_cm / (LIGHT_SPEED_CM_PER_NS * beta_proton)
    return {
        "beta_kaon": float(beta_kaon),
        "beta_proton": float(beta_proton),
        "tof_kaon_ns": float(tof_kaon_ns),
        "tof_proton_ns": float(tof_proton_ns),
        "delta_t_pk_ns": float(tof_proton_ns - tof_kaon_ns),
    }


def _build_shms_tof_payload(evt):
    payload = {
        "tof_valid": False,
        "P_gtr_p": None,
        "ssxptar": None,
        "ssyptar": None,
        "ssdelta": None,
        "shms_path_length_cm": None,
        "beta_kaon": None,
        "beta_proton": None,
        "tof_kaon_ns": None,
        "tof_proton_ns": None,
        "delta_t_pk_ns": None,
        "tof_invalid_reason": "",
        "tof_missing_fields": [],
    }
    missing_fields = []
    for branch_name in ("P_gtr_p", "ssxptar", "ssyptar", "ssdelta"):
        value = _safe_event_float(evt, branch_name)
        payload[branch_name] = value
        if value is None:
            missing_fields.append(branch_name)
    payload["tof_missing_fields"] = list(missing_fields)
    if missing_fields:
        payload["tof_invalid_reason"] = "missing_" + "_".join(missing_fields)
        return payload
    try:
        if float(payload["P_gtr_p"]) <= 0.0:
            payload["tof_invalid_reason"] = "nonpositive_P_gtr_p"
            return payload
        path_length_cm = compute_shms_path_length_cm(
            payload["P_gtr_p"],
            payload["ssxptar"],
            payload["ssyptar"],
            payload["ssdelta"],
        )
        if not math.isfinite(path_length_cm) or path_length_cm <= 0.0:
            payload["tof_invalid_reason"] = "invalid_shms_path_length"
            return payload
        tof_payload = compute_proton_kaon_tof_separation_ns(
            payload["P_gtr_p"],
            path_length_cm,
        )
        if (
            not math.isfinite(float(tof_payload["beta_kaon"]))
            or not math.isfinite(float(tof_payload["beta_proton"]))
            or float(tof_payload["beta_kaon"]) <= 0.0
            or float(tof_payload["beta_proton"]) <= 0.0
            or float(tof_payload["beta_kaon"]) > 1.0
            or float(tof_payload["beta_proton"]) > 1.0
            or not math.isfinite(float(tof_payload["delta_t_pk_ns"]))
            or float(tof_payload["delta_t_pk_ns"]) <= 0.0
        ):
            payload["tof_invalid_reason"] = "invalid_delta_t_pk"
            return payload
        payload.update(tof_payload)
        payload["shms_path_length_cm"] = float(path_length_cm)
        payload["tof_valid"] = True
        payload["tof_invalid_reason"] = ""
        return payload
    except Exception as exc:
        payload["tof_invalid_reason"] = "exception_{}".format(type(exc).__name__)
        return payload


def _mean_rms(values):
    clean_values = [float(value) for value in values if math.isfinite(float(value))]
    if not clean_values:
        return None, None
    array_values = np.asarray(clean_values, dtype=float)
    return float(np.mean(array_values)), float(np.std(array_values))


def _build_prompt_tof_summary_by_delta(
    source_bundle,
    delta_edges,
    selection_key,
    timing_branch=None,
    timing_range=None,
    aero_range=None,
    delta_range=None,
):
    n_delta_bins = max(len(delta_edges) - 1, 0)
    collections = [
        {
            "P_gtr_p": [],
            "ssxptar": [],
            "ssyptar": [],
            "shms_path_length_cm": [],
            "delta_t_pk_ns": [],
        }
        for _ in range(n_delta_bins)
    ]
    counters = [
        {
            "prompt_events_seen": 0,
            "prompt_events_with_valid_tof": 0,
            "prompt_events_with_selected_timing": 0,
            "prompt_events_inside_timing_range": 0,
            "prompt_events_inside_aero_range": 0,
            "prompt_events_used": 0,
        }
        for _ in range(n_delta_bins)
    ]
    timing_branch = str(timing_branch or "")
    timing_min = timing_max = None
    if timing_range and len(timing_range) >= 2:
        timing_min = float(timing_range[0])
        timing_max = float(timing_range[1])
    aero_min = aero_max = None
    if aero_range and len(aero_range) >= 2:
        aero_min = float(aero_range[0])
        aero_max = float(aero_range[1])
    delta_min = delta_max = None
    if delta_range and len(delta_range) >= 2:
        delta_min = float(delta_range[0])
        delta_max = float(delta_range[1])
    for source_name, source_spec, _, entry_payload in _iter_prepared_records(
        source_bundle,
        selection_key=selection_key,
        source_names=("prompt",),
    ):
        del source_name
        if not bool((source_spec or {}).get("is_prompt_source", False)):
            continue
        delta_value = (entry_payload or {}).get("delta_value")
        if delta_value is None:
            continue
        try:
            delta_value = float(delta_value)
        except Exception:
            continue
        if not math.isfinite(delta_value):
            continue
        if delta_min is not None and (delta_value < delta_min or delta_value > delta_max):
            continue
        delta_index = _find_collection_bin(delta_value, delta_edges)
        if not (0 <= delta_index < n_delta_bins):
            continue
        counters[delta_index]["prompt_events_seen"] += 1
        if not bool((entry_payload or {}).get("tof_valid", False)):
            continue
        counters[delta_index]["prompt_events_with_valid_tof"] += 1
        timing_values = (entry_payload or {}).get("timing_values") or {}
        selected_timing = timing_values.get(timing_branch) if timing_branch else None
        if selected_timing is None:
            continue
        try:
            selected_timing = float(selected_timing)
        except Exception:
            continue
        if not math.isfinite(selected_timing):
            continue
        counters[delta_index]["prompt_events_with_selected_timing"] += 1
        if (
            timing_min is not None
            and (selected_timing < timing_min or selected_timing > timing_max)
        ):
            continue
        counters[delta_index]["prompt_events_inside_timing_range"] += 1
        aero_value = (entry_payload or {}).get("aero_value")
        if aero_value is None:
            aero_value = (entry_payload or {}).get("P_aero_npeSum")
        if aero_value is None:
            continue
        try:
            aero_value = float(aero_value)
        except Exception:
            continue
        if not math.isfinite(aero_value):
            continue
        if aero_min is not None and (aero_value < aero_min or aero_value > aero_max):
            continue
        counters[delta_index]["prompt_events_inside_aero_range"] += 1
        counters[delta_index]["prompt_events_used"] += 1
        for key in collections[delta_index]:
            value = (entry_payload or {}).get(key)
            if value is not None and math.isfinite(float(value)):
                collections[delta_index][key].append(float(value))
    summaries = []
    for delta_index, collection in enumerate(collections):
        row = {
            "delta_index": int(delta_index),
            "delta_min": float(delta_edges[delta_index]),
            "delta_max": float(delta_edges[delta_index + 1]),
            "prompt_event_count": int(len(collection["delta_t_pk_ns"])),
            **counters[delta_index],
        }
        for key, values in collection.items():
            mean_value, rms_value = _mean_rms(values)
            row["mean_" + key] = mean_value
            row["rms_" + key] = rms_value
        row["valid"] = bool(
            row["prompt_event_count"] > 0
            and row.get("mean_delta_t_pk_ns") is not None
            and math.isfinite(float(row.get("mean_delta_t_pk_ns")))
            and float(row.get("mean_delta_t_pk_ns")) > 0.0
        )
        summaries.append(row)
    return summaries


def _wrap_rf_mean_to_selected_window(
    predicted_mean,
    reference_mean,
    display_min,
    display_max,
    bunch_spacing_ns,
):
    predicted_mean = float(predicted_mean)
    reference_mean = float(reference_mean)
    display_min = float(display_min)
    display_max = float(display_max)
    bunch_spacing_ns = float(bunch_spacing_ns)
    if (
        not math.isfinite(predicted_mean)
        or not math.isfinite(reference_mean)
        or not math.isfinite(display_min)
        or not math.isfinite(display_max)
        or not math.isfinite(bunch_spacing_ns)
        or bunch_spacing_ns <= 0.0
        or display_max <= display_min
    ):
        return {
            "valid": False,
            "raw_mean": predicted_mean,
            "wrapped_mean": predicted_mean,
            "period_shift": 0,
            "reason": "invalid_rf_wrap_inputs",
        }
    candidates = []
    for period_shift in range(-12, 13):
        candidate = predicted_mean + (float(period_shift) * bunch_spacing_ns)
        inside = display_min <= candidate <= display_max
        edge_distance = 0.0
        if candidate < display_min:
            edge_distance = display_min - candidate
        elif candidate > display_max:
            edge_distance = candidate - display_max
        candidates.append(
            (
                0 if inside else 1,
                float(edge_distance),
                abs(candidate - reference_mean),
                abs(period_shift),
                int(period_shift),
                float(candidate),
            )
        )
    candidates.sort()
    selected = candidates[0]
    return {
        "valid": True,
        "raw_mean": float(predicted_mean),
        "wrapped_mean": float(selected[5]),
        "period_shift": int(selected[4]),
        "inside_display_range": bool(selected[0] == 0),
        "reference_mean": float(reference_mean),
    }


def _format_probe_summary_for_log(label, probe):
    if not isinstance(probe, dict):
        return "{}: unavailable".format(label)
    goodness_key = (
        "meanPoissonDeviancePerEntry"
        if bool(probe.get("perAeroMultistart", False)) or bool(probe.get("localPeakRescue", False))
        else "meanPoissonDevianceNdf"
    )
    return (
        "{}: branch={} kind={} mode={} selection={} valid={} pair={} sep_sigma={} "
        "D/ndf={} D/N={} score={}"
    ).format(
        label,
        probe.get("branch", probe.get("timing_branch", "unknown")),
        probe.get("probe_kind", "unknown"),
        probe.get("fit_mode", "unknown"),
        probe.get("selection_key", "unknown"),
        int(probe.get("validShapes", 0) or 0),
        "yes" if probe.get("peakPairFound", False) else "no",
        _format_debug_float(probe.get("meanSeparation"), digits=4),
        _format_debug_float(probe.get("meanPoissonDevianceNdf"), digits=4),
        _format_debug_float(probe.get("meanPoissonDeviancePerEntry"), digits=4),
        _format_debug_float(probe.get(goodness_key), digits=4),
    )


def _join_rejection_reasons(reasons):
    cleaned = [str(reason).strip() for reason in (reasons or []) if str(reason).strip()]
    return "; ".join(cleaned)


def _build_exact_proton_cleaning_config(base_config):
    exact_config = deepcopy(base_config or {})
    exact_config["implementation"] = PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT
    exact_config.setdefault("timing_probe_policy", "rf_then_ct_best")
    exact_config.setdefault("rf_branch_candidates", tuple(DEFAULT_RF_BRANCH_CANDIDATES))
    exact_config.setdefault("allow_rf_probe", True)
    exact_config.setdefault("disable_rf_timing_env", "PROTON_CHECKER_DISABLE_RF")
    exact_config.setdefault("force_rf_branch_env", "PROTON_CHECKER_RF_BRANCH")
    exact_config["ct_timing_branch"] = PROTON_CLEANING_EXACT_TIMING_BRANCH
    exact_config["ctime_hist_range"] = tuple(PROTON_CLEANING_EXACT_TIMING_RANGE)
    exact_config["ctime_hist_bins"] = int(PROTON_CLEANING_EXACT_TIMING_BINS)
    exact_config["aero_slice_edges"] = tuple(PROTON_CLEANING_EXACT_AERO_EDGES)
    exact_config["aero_hist_range"] = tuple(PROTON_CLEANING_EXACT_AERO_RANGE)
    exact_config["delta_hist_range"] = tuple(PROTON_CLEANING_EXACT_DELTA_RANGE)
    exact_config["delta_bins"] = int(PROTON_CLEANING_EXACT_DELTA_BINS)
    exact_config["mm_validation_range"] = tuple(PROTON_CLEANING_EXACT_MM_VALIDATION_RANGE)
    exact_config["global_fit"] = deepcopy(PROTON_CLEANING_EXACT_GLOBAL_FIT)
    exact_config["slice_fit"] = deepcopy(PROTON_CLEANING_EXACT_SLICE_FIT)
    exact_config["support_thresholds"] = deepcopy(PROTON_CLEANING_EXACT_SUPPORT_THRESHOLDS)
    exact_config["weighting"] = deepcopy(PROTON_CLEANING_EXACT_WEIGHTING)
    exact_config["validation_windows"] = deepcopy(PROTON_CLEANING_EXACT_VALIDATION_WINDOWS)
    exact_config["tof_offset_validation"] = deepcopy(
        (base_config or {}).get("tof_offset_validation")
        or PROTON_CLEANING_TOF_OFFSET_VALIDATION
    )
    exact_config["fit_options"] = PROTON_CLEANING_EXACT_FIT_OPTIONS
    return exact_config


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


def _make_prepared_event_signature(source_label, entry_index):
    return "{}:{}".format(str(source_label), int(entry_index))


def get_kaon_proton_cleaning_event_payload(
    cleaning_result,
    source_label,
    entry_index,
    strict=False,
):
    if not isinstance(cleaning_result, dict) or not bool(cleaning_result.get("accepted")):
        return None
    lookup = cleaning_result.get("_prepared_event_weight_lookup") or {}
    signature = _make_prepared_event_signature(source_label, entry_index)
    payload = lookup.get(signature)
    if payload is None and strict:
        raise KeyError(
            "Missing proton-cleaning payload for source '{}' entry {}".format(
                source_label,
                int(entry_index),
            )
        )
    return payload


def get_kaon_proton_cleaning_event_scale(
    cleaning_result,
    source_label,
    entry_index,
    strict=False,
    scale_key="final_cleaned_factor",
):
    payload = get_kaon_proton_cleaning_event_payload(
        cleaning_result,
        source_label,
        entry_index,
        strict=strict,
    )
    if payload is None:
        return None
    try:
        return float(payload.get(str(scale_key), 1.0) or 0.0)
    except Exception:
        if strict:
            raise
        return None


def _tree_has_branch(tree, branch_name):
    if tree is None or not branch_name:
        return False
    try:
        return bool(tree.GetBranch(str(branch_name)))
    except Exception:
        return False


def _bundle_has_branch(source_bundle, branch_name, source_names=None):
    prepared_branches = (source_bundle or {}).get("available_timing_branches") or ()
    if branch_name in prepared_branches:
        return True
    allowed_sources = None
    if source_names is not None:
        allowed_sources = {str(source_name) for source_name in source_names}
    for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items():
        if allowed_sources is not None and str(source_name) not in allowed_sources:
            continue
        tree = (source_spec or {}).get("tree")
        if _tree_has_branch(tree, branch_name):
            return True
    return False


def _resolve_rf_branch_candidates(source_bundle, config=None):
    config = config or {}
    disable_env = str(config.get("disable_rf_timing_env") or "PROTON_CHECKER_DISABLE_RF")
    force_env = str(config.get("force_rf_branch_env") or "PROTON_CHECKER_RF_BRANCH")
    if not bool(config.get("allow_rf_probe", True)):
        return []
    if str(os.environ.get(disable_env, "")).strip() == "1":
        return []
    candidates = []
    forced = str(os.environ.get(force_env, "")).strip()
    if forced:
        candidates.append(forced)
    for candidate in tuple(config.get("rf_branch_candidates") or DEFAULT_RF_BRANCH_CANDIDATES):
        if candidate not in candidates:
            candidates.append(candidate)
    return [candidate for candidate in candidates if _bundle_has_branch(source_bundle, candidate)]


def _get_prepared_sources(source_bundle):
    return ((source_bundle or {}).get("prepared_sources") or {})


def _iter_prepared_records(source_bundle, require_nommcuts=False, selection_key=None, source_names=None):
    allowed_sources = None
    if source_names is not None:
        allowed_sources = {str(source_name) for source_name in source_names}
    active_selection_key = str(selection_key) if selection_key else ("nommcuts" if require_nommcuts else None)
    for source_name, source_spec in _get_prepared_sources(source_bundle).items():
        if allowed_sources is not None and str(source_name) not in allowed_sources:
            continue
        entry_map = (source_spec or {}).get("entries") or {}
        for entry_index, entry_payload in entry_map.items():
            if active_selection_key and not bool((entry_payload or {}).get(active_selection_key, False)):
                continue
            yield source_name, source_spec, int(entry_index), entry_payload


def _passes_exact_spectrometer_base_acceptance(evt):
    """C macro baseAcceptanceCut without the Q2/W diamond polygon."""
    try:
        return bool(
            float(getattr(evt, "ssdelta")) >= -10.0
            and float(getattr(evt, "ssdelta")) <= 20.0
            and float(getattr(evt, "ssxptar")) >= -0.06
            and float(getattr(evt, "ssxptar")) <= 0.06
            and float(getattr(evt, "ssyptar")) >= -0.04
            and float(getattr(evt, "ssyptar")) <= 0.04
            and float(getattr(evt, "hsdelta")) >= -8.0
            and float(getattr(evt, "hsdelta")) <= 8.0
            and float(getattr(evt, "hsxptar")) >= -0.08
            and float(getattr(evt, "hsxptar")) <= 0.08
            and float(getattr(evt, "hsyptar")) >= -0.045
            and float(getattr(evt, "hsyptar")) <= 0.045
        )
    except Exception:
        return False


def prepare_kaon_proton_cleaning_source_bundle(
    source_bundle,
    evaluate_event,
    shifted_mm_getter,
    shifted_t_getter,
    hole_contains,
    mm_min,
    mm_max,
    proton_cleaning_config=None,
):
    prepared_bundle = dict(source_bundle or {})
    prepared_sources = {}
    prepared_source_stats = {}
    available_timing_branches = set()
    requested_timing_branches = [PROTON_CLEANING_EXACT_TIMING_BRANCH]
    for branch_name in _resolve_rf_branch_candidates(source_bundle, proton_cleaning_config):
        if branch_name not in requested_timing_branches:
            requested_timing_branches.append(branch_name)
    prompt_source_spec = ((source_bundle or {}).get("sources") or {}).get("prompt") or {}
    prompt_physical_coefficient = float(prompt_source_spec.get("coefficient", 0.0) or 0.0)
    if (not math.isfinite(prompt_physical_coefficient)) or prompt_physical_coefficient == 0.0:
        raise RuntimeError(
            "invalid prompt physical coefficient for proton-cleaning fit scaling: {}".format(
                prompt_physical_coefficient
            )
        )

    for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items():
        tree = (source_spec or {}).get("tree")
        entry_payloads = {}
        entries_seen = 0
        entries_prepared = 0
        available_for_source = []
        physical_coefficient = float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        fit_coefficient = float(physical_coefficient / prompt_physical_coefficient)

        if tree is not None:
            available_for_source = [
                str(branch_name)
                for branch_name in requested_timing_branches
                if _tree_has_branch(tree, branch_name)
            ]
            available_timing_branches.update(available_for_source)

            for entry_index, evt in enumerate(tree):
                entries_seen += 1
                base_all_cuts, base_sub_cuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
                hole_rejected = (
                    hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
                    if hole_contains is not None
                    else False
                )
                allcuts = bool(base_all_cuts and (not hole_rejected))
                nommcuts = bool(base_sub_cuts and (not hole_rejected))
                noholecuts = bool(base_all_cuts)
                pre_diamond_nommcuts = bool(
                    _passes_exact_spectrometer_base_acceptance(evt)
                    and (not hole_rejected)
                )

                if not (allcuts or nommcuts or noholecuts or pre_diamond_nommcuts):
                    continue

                timing_values = {}
                for branch_name in available_for_source:
                    try:
                        branch_value = float(getattr(evt, branch_name))
                    except Exception:
                        continue
                    if math.isfinite(branch_value):
                        timing_values[str(branch_name)] = float(branch_value)

                tof_payload = _build_shms_tof_payload(evt)
                entry_payloads[int(entry_index)] = {
                    "allcuts": bool(allcuts),
                    "nommcuts": bool(nommcuts),
                    "noholecuts": bool(noholecuts),
                    "pre_diamond_nommcuts": bool(pre_diamond_nommcuts),
                    "adj_mm": float(shifted_mm_getter(evt)),
                    "adj_t": float(shifted_t_getter(evt)),
                    "adj_hsdelta": float(adj_hsdelta),
                    "delta_value": float(getattr(evt, "ssdelta", 0.0)),
                    "aero_value": float(getattr(evt, "P_aero_npeSum", 0.0)),
                    "timing_values": timing_values,
                    **tof_payload,
                }
                entries_prepared += 1

        prepared_sources[str(source_name)] = {
            **dict(source_spec or {}),
            "coefficient": float(physical_coefficient),
            "fit_coefficient": float(fit_coefficient),
            "is_prompt_source": bool(str(source_name) == "prompt"),
            "entries": entry_payloads,
            "available_timing_branches": list(available_for_source),
        }
        prepared_source_stats[str(source_name)] = {
            "tree_name": (source_spec or {}).get("tree_name"),
            "entries_seen": int(entries_seen),
            "entries_prepared": int(entries_prepared),
            "entries_missing_P_gtr_p": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if "P_gtr_p" in (entry_payload.get("tof_missing_fields") or [])
                )
            ),
            "entries_missing_ssxptar": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if "ssxptar" in (entry_payload.get("tof_missing_fields") or [])
                )
            ),
            "entries_missing_ssyptar": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if "ssyptar" in (entry_payload.get("tof_missing_fields") or [])
                )
            ),
            "entries_missing_ssdelta": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if "ssdelta" in (entry_payload.get("tof_missing_fields") or [])
                )
            ),
            "entries_invalid_shms_path_length": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if entry_payload.get("tof_invalid_reason") == "invalid_shms_path_length"
                )
            ),
            "entries_invalid_delta_t_pk": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if entry_payload.get("tof_invalid_reason") == "invalid_delta_t_pk"
                )
            ),
            "entries_valid_delta_t_pk": int(
                sum(
                    1
                    for entry_payload in entry_payloads.values()
                    if bool(entry_payload.get("tof_valid", False))
                )
            ),
        }

    prepared_bundle["prepared_sources"] = prepared_sources
    prepared_bundle["prepared_source_stats"] = prepared_source_stats
    prepared_bundle["available_timing_branches"] = sorted(available_timing_branches)
    return prepared_bundle


def _preselection_passes(evt, evaluate_event, hole_contains, mm_min, mm_max):
    _, nommcuts, _ = evaluate_event(evt, mm_min, mm_max)
    hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) if hole_contains is not None else False
    return bool(nommcuts and (not hole_rejected))


def _collect_branch_values(
    source_bundle,
    evaluate_event,
    hole_contains,
    mm_min,
    mm_max,
    branch_name,
    source_names=None,
    selection_key="nommcuts",
):
    if _get_prepared_sources(source_bundle):
        values = []
        for _, _, _, entry_payload in _iter_prepared_records(
            source_bundle,
            selection_key=selection_key,
            source_names=source_names,
        ):
            timing_value = ((entry_payload or {}).get("timing_values") or {}).get(str(branch_name))
            if timing_value is None:
                continue
            timing_value = float(timing_value)
            if math.isfinite(timing_value):
                values.append(timing_value)
        return values
    values = []
    allowed_sources = None
    if source_names is not None:
        allowed_sources = {str(source_name) for source_name in source_names}
    for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items():
        if allowed_sources is not None and str(source_name) not in allowed_sources:
            continue
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


def _resolve_beam_bunch_spacing_ns(source_bundle):
    epsset = normalize_epsset((source_bundle or {}).get("epsset"))
    return 4.0 if epsset == "high" else 2.0


def _resolve_ct_probe_configuration(source_bundle):
    epsset = normalize_epsset((source_bundle or {}).get("epsset"))
    if epsset == "high":
        ctime_range = tuple(PROTON_CLEANING_EXACT_CT_HIGH_EPSILON_RANGE)
        ctime_bins = int(PROTON_CLEANING_EXACT_CT_HIGH_EPSILON_BINS)
    else:
        ctime_range = tuple(PROTON_CLEANING_EXACT_CT_LOW_EPSILON_RANGE)
        ctime_bins = int(PROTON_CLEANING_EXACT_CT_LOW_EPSILON_BINS)
    global_fit = deepcopy(PROTON_CLEANING_EXACT_GLOBAL_FIT)
    global_fit["kaon_mean_range"] = (-0.45, 0.20)
    global_fit["proton_mean_range"] = (0.20, 0.95)
    global_fit["sigma_range"] = (0.03, 0.45)
    global_fit["initial_sigma"] = 0.15
    return {
        "timing_branch": PROTON_CLEANING_EXACT_TIMING_BRANCH,
        "display_range": tuple(float(value) for value in ctime_range),
        "fit_range": tuple(float(value) for value in ctime_range),
        "histogram_bins": int(ctime_bins),
        "kaonMeanMin": -0.45,
        "kaonMeanMax": 0.20,
        "protonMeanMin": 0.20,
        "protonMeanMax": 0.95,
        "sigmaMin": 0.03,
        "sigmaMax": 0.45,
        "sigmaInitial": 0.15,
        "minimumGlobalSeparation": 0.75,
        "boundFractionTolerance": 0.02,
        "useDeviancePerEntryValidation": False,
        "maximumPoissonDeviancePerEntry": 0.85,
        "global_fit": global_fit,
    }


def _build_unweighted_timing_projection(
    source_bundle,
    evaluate_event,
    hole_contains,
    mm_min,
    mm_max,
    branch_name,
    histogram_name,
    histogram_range,
    histogram_bins,
    source_names=None,
    selection_key="nommcuts",
):
    time_min, time_max = [float(value) for value in histogram_range]
    projection = ROOT.TH1D(
        str(histogram_name),
        str(histogram_name),
        int(histogram_bins),
        time_min,
        time_max,
    )
    projection.SetDirectory(0)
    projection.Sumw2()
    if _get_prepared_sources(source_bundle):
        for _, _, _, entry_payload in _iter_prepared_records(
            source_bundle,
            selection_key=selection_key,
            source_names=source_names,
        ):
            timing_value = ((entry_payload or {}).get("timing_values") or {}).get(str(branch_name))
            if timing_value is None:
                continue
            timing_value = float(timing_value)
            if (not math.isfinite(timing_value)) or timing_value < time_min or timing_value > time_max:
                continue
            projection.Fill(timing_value)
        return projection
    for source_spec in ((source_bundle or {}).get("sources") or {}).values():
        tree = (source_spec or {}).get("tree")
        if tree is None or not _tree_has_branch(tree, branch_name):
            continue
        for evt in tree:
            if not _preselection_passes(evt, evaluate_event, hole_contains, mm_min, mm_max):
                continue
            try:
                timing_value = float(getattr(evt, str(branch_name)))
            except Exception:
                continue
            if (not math.isfinite(timing_value)) or timing_value < time_min or timing_value > time_max:
                continue
            projection.Fill(timing_value)
    return projection


def _resolve_rf_probe_display_range(
    source_bundle,
    evaluate_event,
    hole_contains,
    mm_min,
    mm_max,
    timing_branch,
    selection_key="nommcuts",
    source_names=("prompt",),
    return_diagnostics=False,
):
    beam_bunch_spacing_ns = _resolve_beam_bunch_spacing_ns(source_bundle)
    fallback_range = (-float(beam_bunch_spacing_ns), float(beam_bunch_spacing_ns))
    branch_values = _collect_branch_values(
        source_bundle,
        evaluate_event,
        hole_contains,
        mm_min,
        mm_max,
        timing_branch,
        source_names=source_names,
        selection_key=selection_key,
    )
    diagnostics = {
        "rf_window_discovery_source": ",".join(str(source_name) for source_name in (source_names or ())),
        "rf_window_discovery_entry_count": int(len(branch_values)),
        "rf_window_discovery_branch": str(timing_branch),
        "rf_window_discovery_selection_key": str(selection_key),
        "rf_window_discovery_range": None,
    }
    raw_display_range = _estimate_value_central_range(
        branch_values,
        fallback_range[0],
        fallback_range[1],
    )
    raw_display_min, raw_display_max = [float(value) for value in raw_display_range]
    if raw_display_max <= raw_display_min:
        diagnostics["rf_window_discovery_range"] = [float(fallback_range[0]), float(fallback_range[1])]
        return (fallback_range, diagnostics) if return_diagnostics else fallback_range
    raw_display_width = float(raw_display_max - raw_display_min)
    search_bins = max(
        160,
        int(round(raw_display_width / 0.015)),
    )
    rf_window_search = ROOT.TH1D(
        "H_proton_cleaning_rf_window_search_{}".format(str(timing_branch).replace(" ", "_")),
        "RF window search",
        int(search_bins),
        float(raw_display_min),
        float(raw_display_max),
    )
    rf_window_search.SetDirectory(0)
    rf_window_search.Sumw2()
    for value in branch_values:
        timing_value = float(value)
        if (not math.isfinite(timing_value)) or timing_value < raw_display_min or timing_value > raw_display_max:
            continue
        rf_window_search.Fill(timing_value)
    raw_peak_pair = _find_separated_peak_pair(
        rf_window_search,
        raw_display_min,
        raw_display_max,
        max(0.18, 0.10 * beam_bunch_spacing_ns),
        0.80 * beam_bunch_spacing_ns,
    )
    if bool(raw_peak_pair.get("valid")) and raw_display_width > 1.15 * beam_bunch_spacing_ns:
        pair_center = 0.5 * (
            float(raw_peak_pair["lower_mean"]) +
            float(raw_peak_pair["upper_mean"])
        )
        selected_width = float(beam_bunch_spacing_ns)
        selected_min = float(pair_center - (0.5 * selected_width))
        selected_max = float(pair_center + (0.5 * selected_width))
        if selected_min < raw_display_min:
            selected_max += raw_display_min - selected_min
            selected_min = raw_display_min
        if selected_max > raw_display_max:
            selected_min -= selected_max - raw_display_max
            selected_max = raw_display_max
        selected_range = (
            max(float(selected_min), raw_display_min),
            min(float(selected_max), raw_display_max),
        )
        diagnostics["rf_window_discovery_range"] = [float(selected_range[0]), float(selected_range[1])]
        return (selected_range, diagnostics) if return_diagnostics else selected_range
    selected_range = (raw_display_min, raw_display_max)
    diagnostics["rf_window_discovery_range"] = [float(selected_range[0]), float(selected_range[1])]
    return (selected_range, diagnostics) if return_diagnostics else selected_range


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


def _evaluate_gaussian_scalar(x_value, amplitude, mean, sigma):
    if (
        (not math.isfinite(float(x_value)))
        or (not math.isfinite(float(amplitude)))
        or (not math.isfinite(float(mean)))
        or (not math.isfinite(float(sigma)))
        or float(amplitude) <= 0.0
        or float(sigma) <= 0.0
    ):
        return 0.0
    z_value = (float(x_value) - float(mean)) / float(sigma)
    return float(amplitude) * math.exp(-0.5 * z_value * z_value)


def _sum_gaussian_over_bins(histogram, amplitude, mean, sigma, fit_min, fit_max):
    if histogram is None or float(amplitude) <= 0.0 or float(sigma) <= 0.0:
        return 0.0
    total = 0.0
    for bin_index in range(1, histogram.GetNbinsX() + 1):
        x_value = float(histogram.GetXaxis().GetBinCenter(bin_index))
        if x_value < float(fit_min) or x_value > float(fit_max):
            continue
        total += _evaluate_gaussian_scalar(
            x_value,
            amplitude,
            mean,
            sigma,
        )
    return float(total)


def _sum_constant_over_bins(histogram, amplitude, fit_min, fit_max):
    if histogram is None or float(amplitude) <= 0.0:
        return 0.0
    total = 0.0
    for bin_index in range(1, histogram.GetNbinsX() + 1):
        x_value = float(histogram.GetXaxis().GetBinCenter(bin_index))
        if x_value < float(fit_min) or x_value > float(fit_max):
            continue
        total += float(amplitude)
    return float(total)


def _compute_poisson_goodness_of_fit(
    histogram,
    fit_function,
    fit_min,
    fit_max,
    number_of_free_parameters,
):
    result = {
        "deviance": 0.0,
        "fitted_bins": 0,
        "fitted_entries": 0.0,
        "ndf": 0,
        "deviance_ndf": None,
        "deviance_per_entry": None,
    }
    if (
        histogram is None
        or fit_function is None
        or float(fit_max) <= float(fit_min)
        or int(number_of_free_parameters) < 0
    ):
        return result
    deviance = 0.0
    fitted_entries = 0.0
    fitted_bins = 0
    for bin_index in range(1, histogram.GetNbinsX() + 1):
        bin_center = float(histogram.GetXaxis().GetBinCenter(bin_index))
        if bin_center < float(fit_min) or bin_center > float(fit_max):
            continue
        bin_low = float(histogram.GetXaxis().GetBinLowEdge(bin_index))
        bin_high = float(histogram.GetXaxis().GetBinUpEdge(bin_index))
        bin_width = float(bin_high - bin_low)
        if not (bin_width > 0.0):
            continue
        observed = max(float(histogram.GetBinContent(bin_index)), 0.0)
        expected = max(float(fit_function.Integral(bin_low, bin_high)) / bin_width, 1.0e-12)
        if observed > 0.0:
            deviance += 2.0 * (
                expected
                - observed
                + observed * math.log(observed / expected)
            )
        else:
            deviance += 2.0 * expected
        fitted_entries += observed
        fitted_bins += 1
    ndf = int(fitted_bins - int(number_of_free_parameters))
    result["deviance"] = float(deviance)
    result["fitted_bins"] = int(fitted_bins)
    result["fitted_entries"] = float(fitted_entries)
    result["ndf"] = int(ndf)
    if ndf > 0:
        result["deviance_ndf"] = float(deviance / float(ndf))
    if fitted_entries > 0.0:
        result["deviance_per_entry"] = float(deviance / float(fitted_entries))
    return result


def _is_near_bound(value, lower, upper, bound_fraction_tolerance):
    tolerance = float(bound_fraction_tolerance) * float(upper - lower)
    return (
        float(value) <= float(lower) + tolerance
        or float(value) >= float(upper) - tolerance
    )


def _extract_root_fit_matrices(fit_result, parameter_names):
    covariance_matrix = {}
    correlation_matrix = {}
    uncertainties = {}
    if fit_result is None:
        return covariance_matrix, correlation_matrix, uncertainties
    try:
        covariance = fit_result.GetCovarianceMatrix()
    except Exception:
        covariance = None
    try:
        correlation = fit_result.GetCorrelationMatrix()
    except Exception:
        correlation = None
    for index, parameter_name in enumerate(parameter_names):
        try:
            uncertainties[str(parameter_name)] = float(fit_result.ParError(index))
        except Exception:
            uncertainties[str(parameter_name)] = 0.0
        covariance_row = {}
        correlation_row = {}
        for other_index, other_name in enumerate(parameter_names):
            covariance_value = None
            correlation_value = None
            if covariance is not None:
                try:
                    covariance_value = float(covariance[index][other_index])
                except Exception:
                    covariance_value = None
            if correlation is not None:
                try:
                    correlation_value = float(correlation[index][other_index])
                except Exception:
                    correlation_value = None
            covariance_row[str(other_name)] = covariance_value
            correlation_row[str(other_name)] = correlation_value
        covariance_matrix[str(parameter_name)] = covariance_row
        correlation_matrix[str(parameter_name)] = correlation_row
    return covariance_matrix, correlation_matrix, uncertainties


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
    use_deviance_per_entry_validation=False,
    maximum_poisson_deviance_per_entry=None,
    support_entries=None,
):
    if support_entries is None:
        support_entries = int(round(float(histogram.GetEntries()))) if histogram is not None else 0
    support_snapshot = _build_hist_support_snapshot(
        histogram,
        support_entries,
        minimum_entries,
    )
    if histogram is None or int(support_entries) < int(minimum_entries):
        return {
            "valid": False,
            "fit_attempted": False,
            "fit_status": "insufficient_support",
            "fit_status_code": None,
            "excluded_invalid_variance_bins": 0,
            "invalid_bin_rule": "macro ROOT fit uses all histogram bins in the fit range",
            "function_name": str(function_name),
            "fit_min": float(fit_min),
            "fit_max": float(fit_max),
            "kaonMeanMin": float(kaon_mean_min),
            "kaonMeanMax": float(kaon_mean_max),
            "protonMeanMin": float(proton_mean_min),
            "protonMeanMax": float(proton_mean_max),
            "sigmaMin": float(sigma_min),
            "sigmaMax": float(sigma_max),
            "sigmaInitial": float(initial_sigma),
            "minimumGlobalSeparation": float(minimum_separation),
            "boundFractionTolerance": float(bound_fraction_tolerance),
            "useDeviancePerEntryValidation": bool(use_deviance_per_entry_validation),
            "maximumPoissonDeviancePerEntry": (
                float(maximum_poisson_deviance_per_entry)
                if maximum_poisson_deviance_per_entry is not None
                else None
            ),
            "per_aero_fallback": False,
            "peak_pair_found": False,
            "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
            "rejection_reasons": ["insufficient_entries"],
            "rejection_reason": "insufficient entries",
            **support_snapshot,
        }
    histogram_maximum = max(float(histogram.GetMaximum()), 1.0)
    kaon_seed = _find_peak_seed(histogram, kaon_mean_min, kaon_mean_max)
    proton_seed = _find_peak_seed(histogram, proton_mean_min, proton_mean_max)
    fit_function = ROOT.TF1(
        str(function_name),
        "[0] * exp(-0.5 * pow((x - [1]) / [2], 2))"
        " + "
        "[3] * exp(-0.5 * pow((x - [4]) / [5], 2))"
        " + [6]",
        float(fit_min),
        float(fit_max),
    )
    fit_function.SetParName(0, "K amplitude")
    fit_function.SetParName(1, "K mean")
    fit_function.SetParName(2, "K sigma")
    fit_function.SetParName(3, "p amplitude")
    fit_function.SetParName(4, "p mean")
    fit_function.SetParName(5, "p sigma")
    fit_function.SetParName(6, "other constant")
    fit_function.SetParameter(
        0,
        max(float(kaon_seed[0]), 0.20 * histogram_maximum),
    )
    fit_function.SetParLimits(
        0,
        0.0,
        100.0 * histogram_maximum,
    )
    fit_function.SetParameter(
        1,
        float(np.clip(kaon_seed[1], np.nextafter(float(kaon_mean_min), float(kaon_mean_max)), np.nextafter(float(kaon_mean_max), float(kaon_mean_min)))),
    )
    fit_function.SetParLimits(1, float(kaon_mean_min), float(kaon_mean_max))
    fit_function.SetParameter(
        2,
        float(np.clip(initial_sigma, float(sigma_min), float(sigma_max))),
    )
    fit_function.SetParLimits(2, float(sigma_min), float(sigma_max))
    fit_function.SetParameter(
        3,
        max(float(proton_seed[0]), 0.10 * histogram_maximum),
    )
    fit_function.SetParLimits(
        3,
        0.0,
        100.0 * histogram_maximum,
    )
    fit_function.SetParameter(
        4,
        float(np.clip(proton_seed[1], np.nextafter(float(proton_mean_min), float(proton_mean_max)), np.nextafter(float(proton_mean_max), float(proton_mean_min)))),
    )
    fit_function.SetParLimits(4, float(proton_mean_min), float(proton_mean_max))
    fit_function.SetParameter(
        5,
        float(np.clip(initial_sigma, float(sigma_min), float(sigma_max))),
    )
    fit_function.SetParLimits(5, float(sigma_min), float(sigma_max))
    fit_function.SetParameter(
        6,
        0.02 * histogram_maximum,
    )
    fit_function.SetParLimits(
        6,
        0.0,
        10.0 * histogram_maximum,
    )
    fit_result = histogram.Fit(fit_function, PROTON_CLEANING_EXACT_FIT_OPTIONS)
    fit_status_code = int(fit_result)
    covariance_matrix, correlation_matrix, uncertainties = _extract_root_fit_matrices(
        fit_result,
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
    kaon_amplitude = float(fit_function.GetParameter(0))
    kaon_mean = float(fit_function.GetParameter(1))
    kaon_sigma = abs(float(fit_function.GetParameter(2)))
    proton_amplitude = float(fit_function.GetParameter(3))
    proton_mean = float(fit_function.GetParameter(4))
    proton_sigma = abs(float(fit_function.GetParameter(5)))
    other_amplitude = float(fit_function.GetParameter(6))
    first_fit_bin = max(1, histogram.GetXaxis().FindFixBin(float(fit_min)))
    last_fit_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(fit_max), float(fit_min))),
    )
    active_bin_count = max(0, int(last_fit_bin - first_fit_bin + 1))
    fitted_entries = float(histogram.Integral(first_fit_bin, last_fit_bin))
    support_snapshot = _build_hist_support_snapshot(
        histogram,
        support_entries,
        minimum_entries,
        first_fit_bin=first_fit_bin,
        last_fit_bin=last_fit_bin,
    )
    chi2_data = float(fit_function.GetChisquare())
    fit_ndf = int(fit_function.GetNDF())
    chi2_ndf = (
        float(chi2_data / float(fit_ndf))
        if fit_ndf > 0 and math.isfinite(float(chi2_data))
        else None
    )
    chi2_per_abs_entry = (
        float(chi2_data / abs(fitted_entries))
        if abs(fitted_entries) > 0.0 and math.isfinite(float(chi2_data))
        else None
    )
    goodness = _compute_poisson_goodness_of_fit(
        histogram,
        fit_function,
        fit_min,
        fit_max,
        7,
    )
    poisson_deviance = float(goodness.get("deviance", 0.0) or 0.0)
    poisson_ndf = int(goodness.get("ndf", 0) or 0)
    poisson_deviance_ndf = goodness.get("deviance_ndf")
    poisson_deviance_per_entry = goodness.get("deviance_per_entry")
    kaon_amp_err = float(fit_function.GetParError(0))
    proton_amp_err = float(fit_function.GetParError(3))
    separation_denominator = math.sqrt(max((kaon_sigma ** 2) + (proton_sigma ** 2), 0.0))
    separation = (
        abs(float(proton_mean) - float(kaon_mean)) / separation_denominator
        if separation_denominator > 0.0
        else 0.0
    )
    bound_hit = (
        _is_near_bound(kaon_mean, kaon_mean_min, kaon_mean_max, bound_fraction_tolerance)
        or _is_near_bound(proton_mean, proton_mean_min, proton_mean_max, bound_fraction_tolerance)
        or _is_near_bound(kaon_sigma, sigma_min, sigma_max, bound_fraction_tolerance)
        or _is_near_bound(proton_sigma, sigma_min, sigma_max, bound_fraction_tolerance)
    )
    ordering_valid = (
        float(proton_mean) < float(kaon_mean)
        if bool(proton_peak_is_lower)
        else float(kaon_mean) < float(proton_mean)
    )
    kaon_significance = float(kaon_amplitude / kaon_amp_err) if kaon_amp_err > 0.0 else 0.0
    proton_significance = float(proton_amplitude / proton_amp_err) if proton_amp_err > 0.0 else 0.0
    poisson_ndf_valid = bool(
        poisson_deviance_ndf is not None and math.isfinite(float(poisson_deviance_ndf))
    )
    finite_parameters = (
        math.isfinite(float(kaon_amplitude))
        and math.isfinite(float(proton_amplitude))
        and math.isfinite(float(kaon_mean))
        and math.isfinite(float(proton_mean))
        and math.isfinite(float(kaon_sigma))
        and math.isfinite(float(proton_sigma))
        and math.isfinite(float(other_amplitude))
    )
    rejection_reasons = []
    if fit_status_code != 0:
        rejection_reasons.append("fit_status_{}".format(int(fit_status_code)))
    if not finite_parameters:
        rejection_reasons.append("nonfinite_parameters")
    if bound_hit:
        rejection_reasons.append("bound_hit")
    if not ordering_valid:
        rejection_reasons.append("ordering_invalid")
    if float(kaon_amplitude) <= 0.0:
        rejection_reasons.append("nonpositive_kaon_amplitude")
    if float(proton_amplitude) <= 0.0:
        rejection_reasons.append("nonpositive_proton_amplitude")
    if kaon_sigma <= 0.0:
        rejection_reasons.append("nonpositive_kaon_sigma")
    if proton_sigma <= 0.0:
        rejection_reasons.append("nonpositive_proton_sigma")
    if separation < float(minimum_separation):
        rejection_reasons.append("insufficient_separation")
    if kaon_significance < float(minimum_amplitude_significance):
        rejection_reasons.append("low_kaon_significance")
    if proton_significance < float(minimum_amplitude_significance):
        rejection_reasons.append("low_proton_significance")
    if not poisson_ndf_valid:
        rejection_reasons.append("invalid_poisson_deviance_ndf")
    else:
        if bool(use_deviance_per_entry_validation):
            if (
                poisson_deviance_per_entry is None
                or not math.isfinite(float(poisson_deviance_per_entry))
            ):
                rejection_reasons.append("invalid_poisson_deviance_per_entry")
            elif float(poisson_deviance_per_entry) > float(maximum_poisson_deviance_per_entry):
                rejection_reasons.append("poisson_deviance_per_entry_exceeds_max")
        elif float(poisson_deviance_ndf) > float(maximum_chi2_ndf):
            rejection_reasons.append("poisson_deviance_ndf_exceeds_max")
    valid = len(rejection_reasons) == 0
    return {
        "valid": bool(valid),
        "fit_attempted": True,
        "fit_status": "success" if fit_status_code == 0 else "failure",
        "fit_status_code": int(fit_status_code),
        "message": "",
        "function_name": str(function_name),
        "kaon_amplitude": float(kaon_amplitude),
        "kaon_amplitude_error": kaon_amp_err,
        "kaon_mean": float(kaon_mean),
        "kaon_sigma": float(kaon_sigma),
        "proton_amplitude": float(proton_amplitude),
        "proton_amplitude_error": proton_amp_err,
        "proton_mean": float(proton_mean),
        "proton_sigma": float(proton_sigma),
        "other_amplitude": float(other_amplitude),
        "separation": float(separation),
        "kaon_significance": float(kaon_significance),
        "proton_significance": float(proton_significance),
        "chi2_data": chi2_data,
        "chi2_ndf": float(chi2_ndf) if chi2_ndf is not None else None,
        "chi2_per_abs_entry": float(chi2_per_abs_entry) if chi2_per_abs_entry is not None else None,
        "poisson_deviance": float(poisson_deviance),
        "poisson_deviance_ndf": float(poisson_deviance_ndf) if poisson_deviance_ndf is not None else None,
        "poisson_deviance_per_entry": (
            float(poisson_deviance_per_entry)
            if poisson_deviance_per_entry is not None
            else None
        ),
        "goodness_ndf": int(poisson_ndf),
        "fitted_entries": float(fitted_entries),
        "active_bin_count": active_bin_count,
        "excluded_invalid_variance_bins": 0,
        "invalid_bin_rule": "macro ROOT fit uses all histogram bins in the fit range",
        "bound_hit": bool(bound_hit),
        "ordering_valid": bool(ordering_valid),
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "uncertainties": uncertainties,
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "kaonMeanMin": float(kaon_mean_min),
        "kaonMeanMax": float(kaon_mean_max),
        "protonMeanMin": float(proton_mean_min),
        "protonMeanMax": float(proton_mean_max),
        "sigmaMin": float(sigma_min),
        "sigmaMax": float(sigma_max),
        "sigmaInitial": float(initial_sigma),
        "minimumGlobalSeparation": float(minimum_separation),
        "boundFractionTolerance": float(bound_fraction_tolerance),
        "useDeviancePerEntryValidation": bool(use_deviance_per_entry_validation),
        "maximumPoissonDeviancePerEntry": (
            float(maximum_poisson_deviance_per_entry)
            if maximum_poisson_deviance_per_entry is not None
            else None
        ),
        "per_aero_fallback": False,
        "peak_pair_found": False,
        "applied_mean_offset": 0.0,
        "offset_adjusted": False,
        "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
        "rejection_reasons": rejection_reasons,
        "rejection_reason": _join_rejection_reasons(rejection_reasons),
        **support_snapshot,
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
            use_deviance_per_entry_validation=True,
            maximum_poisson_deviance_per_entry=1.0e9,
        )
        attempt["fit_min"] = float(central_min)
        attempt["fit_max"] = float(central_max)
        attempt["per_aero_fallback"] = True
        attempt["peak_pair_found"] = True
        fit_status_code = attempt.get("fit_status_code")
        fit_status_accepted = bool(
            fit_status_code == 0
            or (
                fit_status_code == 4
                and math.isfinite(float(attempt.get("poisson_deviance_per_entry", float("inf"))))
                and float(attempt.get("poisson_deviance_per_entry", float("inf"))) <= 0.08
            )
        )
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
            and math.isfinite(float(attempt.get("poisson_deviance_per_entry", float("inf"))))
            and float(attempt.get("kaon_amplitude", 0.0)) > 0.0
            and float(attempt.get("proton_amplitude", 0.0)) > 0.0
            and float(attempt.get("kaon_sigma", 0.0)) > 0.0
            and float(attempt.get("proton_sigma", 0.0)) > 0.0
            and float(attempt.get("separation", 0.0)) >= 0.55
            and float(attempt.get("kaon_significance", 0.0)) >= required_significance
            and float(attempt.get("proton_significance", 0.0)) >= required_significance
            and smaller_amplitude_fraction >= 0.020
            and float(attempt.get("poisson_deviance_per_entry", float("inf"))) <= 0.35
        )
        population_penalty = 0.20 * (0.08 - smaller_amplitude_fraction) if smaller_amplitude_fraction < 0.08 else 0.0
        diagnostic_score = float(attempt.get("poisson_deviance_per_entry", float("inf"))) + population_penalty
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
        and (
            result.get("fit_status_code") == 0
            or (
                result.get("fit_status_code") == 4
                and math.isfinite(float(result.get("poisson_deviance_per_entry", float("inf"))))
                and float(result.get("poisson_deviance_per_entry", float("inf"))) <= 0.08
            )
        )
        and final_ordering_accepted
        and float(result.get("separation", 0.0)) >= 0.55
        and float(result.get("kaon_significance", 0.0)) >= final_required_significance
        and float(result.get("proton_significance", 0.0)) >= final_required_significance
        and final_smaller_amplitude_fraction >= 0.020
        and math.isfinite(float(result.get("poisson_deviance_per_entry", float("inf"))))
        and float(result.get("poisson_deviance_per_entry", float("inf"))) <= 0.35
    )
    diagnostic_warnings = [
        str(warning)
        for warning in (result.get("diagnostic_warnings") or [])
        if str(warning).strip()
    ]
    if bool(result.get("bound_hit", False)) and "bound_hit" not in diagnostic_warnings:
        diagnostic_warnings.append("bound_hit")
    if result.get("fit_status_code") == 4 and bool(result.get("valid", False)):
        diagnostic_warnings.append("fit_status_4_recovered")
    result["diagnostic_warnings"] = sorted(set(diagnostic_warnings))
    if bool(result.get("valid", False)):
        result["rejection_reasons"] = []
        result["rejection_reason"] = ""
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


def _valid_global_shapes_with_weights(global_shapes, support_entries_by_aero=None):
    rows = []
    support_entries_by_aero = support_entries_by_aero or []
    for aero_index, shape in enumerate(global_shapes or []):
        if not bool((shape or {}).get("valid", False)):
            continue
        try:
            weight = float(support_entries_by_aero[aero_index])
        except Exception:
            weight = float((shape or {}).get("support_entries", 0.0) or 0.0)
        if not math.isfinite(weight) or weight <= 0.0:
            weight = 1.0
        rows.append((int(aero_index), shape, float(weight)))
    return rows


def _weighted_shape_average(valid_shape_rows, key):
    values = []
    weights = []
    for _, shape, weight in valid_shape_rows:
        value = (shape or {}).get(key)
        if value is None:
            continue
        try:
            value = float(value)
        except Exception:
            continue
        if not math.isfinite(value):
            continue
        values.append(value)
        weights.append(float(weight))
    if not values:
        return None
    return float(np.average(np.asarray(values, dtype=float), weights=np.asarray(weights, dtype=float)))


def _build_reference_shape_for_delta_offset(global_shapes, support_entries_by_aero=None):
    valid_shape_rows = _valid_global_shapes_with_weights(global_shapes, support_entries_by_aero)
    if not valid_shape_rows:
        return {"valid": False, "reason": "no_valid_global_shapes"}
    reference = {
        "valid": True,
        "source_aero_indices": [int(row[0]) for row in valid_shape_rows],
        "kaon_mean": _weighted_shape_average(valid_shape_rows, "kaon_mean"),
        "proton_mean": _weighted_shape_average(valid_shape_rows, "proton_mean"),
        "kaon_sigma": _weighted_shape_average(valid_shape_rows, "kaon_sigma"),
        "proton_sigma": _weighted_shape_average(valid_shape_rows, "proton_sigma"),
        "fit_min": _weighted_shape_average(valid_shape_rows, "fit_min"),
        "fit_max": _weighted_shape_average(valid_shape_rows, "fit_max"),
    }
    required_keys = ("kaon_mean", "proton_mean", "kaon_sigma", "proton_sigma", "fit_min", "fit_max")
    if any(reference.get(key) is None for key in required_keys):
        reference["valid"] = False
        reference["reason"] = "missing_reference_shape_values"
        return reference
    if (
        float(reference["kaon_sigma"]) <= 0.0
        or float(reference["proton_sigma"]) <= 0.0
        or float(reference["fit_max"]) <= float(reference["fit_min"])
    ):
        reference["valid"] = False
        reference["reason"] = "invalid_reference_shape_values"
        return reference
    return reference


def _fit_delta_common_timing_offset(
    histogram,
    reference_shape,
    tof_summary,
    exact_config,
    function_name,
    proton_peak_is_lower,
    probe_kind,
    display_range,
    bunch_spacing_ns,
    support_entries=0,
):
    offset_min, offset_max = [float(value) for value in PROTON_CLEANING_TOF_OFFSET_RANGE]
    validation_cfg = (
        (exact_config or {}).get("tof_offset_validation")
        or PROTON_CLEANING_TOF_OFFSET_VALIDATION
    )
    maximum_offset_error_ns = float(validation_cfg.get("maximum_offset_error_ns", 0.10))
    maximum_chi2_ndf = float(validation_cfg.get("maximum_chi2_ndf", 5.0))
    minimum_component_significance = float(
        validation_cfg.get("minimum_component_significance", 2.0)
    )
    minimum_smaller_component_fraction = float(
        validation_cfg.get("minimum_smaller_component_fraction", 0.02)
    )
    reject_bound_hit_with_large_error = bool(
        validation_cfg.get("reject_bound_hit_with_large_error", True)
    )
    bound_hit_large_error_fraction = float(
        validation_cfg.get("bound_hit_large_error_fraction", 0.50)
    )
    base_result = {
        "valid": False,
        "fit_attempted": False,
        "function_name": str(function_name),
        "delta_offset": 0.0,
        "delta_offset_error": None,
        "fit_status": "not_attempted",
        "fit_status_code": None,
        "delta_offset_bound_hit": False,
        "diagnostic_warnings": [],
        "rejection_reasons": [],
        "rejection_reason": "",
        "mean_delta_t_pk_ns": (tof_summary or {}).get("mean_delta_t_pk_ns"),
        "prompt_tof_event_count": int((tof_summary or {}).get("prompt_event_count", 0) or 0),
        "reference_kaon_mean": (reference_shape or {}).get("kaon_mean"),
        "reference_proton_mean": (reference_shape or {}).get("proton_mean"),
        "reference_kaon_sigma": (reference_shape or {}).get("kaon_sigma"),
        "reference_proton_sigma": (reference_shape or {}).get("proton_sigma"),
        "fit_min": None,
        "fit_max": None,
    }
    rejection_reasons = []
    if histogram is None:
        rejection_reasons.append("missing_delta_timing_histogram")
    if not bool((reference_shape or {}).get("valid", False)):
        rejection_reasons.append((reference_shape or {}).get("reason") or "invalid_reference_shape")
    if not bool((tof_summary or {}).get("valid", False)):
        rejection_reasons.append("invalid_prompt_tof_summary")
    mean_delta_t = (tof_summary or {}).get("mean_delta_t_pk_ns")
    if mean_delta_t is None or not math.isfinite(float(mean_delta_t)) or float(mean_delta_t) <= 0.0:
        rejection_reasons.append("invalid_mean_delta_t_pk")
    minimum_entries = int((exact_config.get("slice_fit") or {}).get("minimum_entries", 30))
    if int(support_entries or 0) < int(minimum_entries):
        rejection_reasons.append("insufficient_delta_offset_support")
    if rejection_reasons:
        base_result["rejection_reasons"] = list(rejection_reasons)
        base_result["rejection_reason"] = _join_rejection_reasons(rejection_reasons)
        return base_result

    fit_min = float(reference_shape["fit_min"])
    fit_max = float(reference_shape["fit_max"])
    base_result["fit_min"] = fit_min
    base_result["fit_max"] = fit_max
    if fit_max <= fit_min:
        base_result["rejection_reasons"] = ["invalid_fit_range"]
        base_result["rejection_reason"] = "invalid_fit_range"
        return base_result

    reference_kaon_mean = float(reference_shape["kaon_mean"])
    reference_proton_mean = float(reference_shape["proton_mean"])
    kaon_sigma = float(reference_shape["kaon_sigma"])
    proton_sigma = float(reference_shape["proton_sigma"])
    branch_sign = -1.0 if bool(proton_peak_is_lower) else 1.0
    raw_proton_reference = reference_kaon_mean + (branch_sign * float(mean_delta_t))
    wrap_info = {
        "valid": True,
        "raw_mean": float(raw_proton_reference),
        "wrapped_mean": float(raw_proton_reference),
        "period_shift": 0,
        "inside_display_range": True,
        "reference_mean": float(reference_proton_mean),
    }
    if str(probe_kind) == "rf":
        wrap_info = _wrap_rf_mean_to_selected_window(
            raw_proton_reference,
            reference_proton_mean,
            float(display_range[0]),
            float(display_range[1]),
            float(bunch_spacing_ns),
        )
    if not bool(wrap_info.get("valid", False)):
        base_result["rejection_reasons"] = [wrap_info.get("reason") or "invalid_rf_wrap"]
        base_result["rejection_reason"] = _join_rejection_reasons(base_result["rejection_reasons"])
        return base_result
    effective_proton_reference = float(wrap_info["wrapped_mean"])
    histogram_maximum = max(float(histogram.GetMaximum()), 1.0)
    fit_function = ROOT.TF1(
        str(function_name),
        "[0] * exp(-0.5 * pow((x - ([1] + [6])) / [2], 2))"
        " + "
        "[3] * exp(-0.5 * pow((x - ([4] + [6])) / [5], 2))"
        " + [7]",
        fit_min,
        fit_max,
    )
    fit_function.SetParName(0, "K amplitude")
    fit_function.SetParName(1, "K reference mean")
    fit_function.SetParName(2, "K sigma")
    fit_function.SetParName(3, "p amplitude")
    fit_function.SetParName(4, "p reference mean")
    fit_function.SetParName(5, "p sigma")
    fit_function.SetParName(6, "delta offset")
    fit_function.SetParName(7, "other constant")
    kaon_seed = _find_peak_seed(histogram, reference_kaon_mean - (2.0 * kaon_sigma), reference_kaon_mean + (2.0 * kaon_sigma))
    proton_seed = _find_peak_seed(histogram, effective_proton_reference - (2.0 * proton_sigma), effective_proton_reference + (2.0 * proton_sigma))
    fit_function.SetParameter(0, max(float(kaon_seed[0]), 0.10 * histogram_maximum))
    fit_function.SetParLimits(0, 0.0, 100.0 * histogram_maximum)
    fit_function.FixParameter(1, reference_kaon_mean)
    fit_function.FixParameter(2, kaon_sigma)
    fit_function.SetParameter(3, max(float(proton_seed[0]), 0.05 * histogram_maximum))
    fit_function.SetParLimits(3, 0.0, 100.0 * histogram_maximum)
    fit_function.FixParameter(4, effective_proton_reference)
    fit_function.FixParameter(5, proton_sigma)
    fit_function.SetParameter(6, 0.0)
    fit_function.SetParLimits(6, offset_min, offset_max)
    fit_function.SetParameter(7, 0.02 * histogram_maximum)
    fit_function.SetParLimits(7, 0.0, 10.0 * histogram_maximum)
    fit_result = histogram.Fit(fit_function, PROTON_CLEANING_EXACT_FIT_OPTIONS)
    fit_status_code = int(fit_result)
    offset = float(fit_function.GetParameter(6))
    offset_error = float(fit_function.GetParError(6))
    bound_hit = _is_near_bound(offset, offset_min, offset_max, 0.02)
    chi2_data = float(fit_function.GetChisquare())
    fit_ndf = int(fit_function.GetNDF())
    chi2_ndf = (
        float(chi2_data / float(fit_ndf))
        if fit_ndf > 0 and math.isfinite(float(chi2_data))
        else None
    )
    goodness = _compute_poisson_goodness_of_fit(
        histogram,
        fit_function,
        fit_min,
        fit_max,
        3,
    )
    poisson_deviance = float(goodness.get("deviance", 0.0) or 0.0)
    poisson_ndf = int(goodness.get("ndf", 0) or 0)
    poisson_deviance_ndf = goodness.get("deviance_ndf")
    poisson_deviance_per_entry = goodness.get("deviance_per_entry")
    kaon_amplitude = float(fit_function.GetParameter(0))
    kaon_amplitude_error = float(fit_function.GetParError(0))
    proton_amplitude = float(fit_function.GetParameter(3))
    proton_amplitude_error = float(fit_function.GetParError(3))
    other_amplitude = float(fit_function.GetParameter(7))
    other_amplitude_error = float(fit_function.GetParError(7))
    kaon_significance = (
        float(kaon_amplitude / kaon_amplitude_error)
        if kaon_amplitude_error > 0.0 and math.isfinite(kaon_amplitude_error)
        else None
    )
    proton_significance = (
        float(proton_amplitude / proton_amplitude_error)
        if proton_amplitude_error > 0.0 and math.isfinite(proton_amplitude_error)
        else None
    )
    component_denominator = kaon_amplitude + proton_amplitude + max(other_amplitude, 0.0)
    smaller_component_fraction = (
        float(min(kaon_amplitude, proton_amplitude) / component_denominator)
        if component_denominator > 0.0 and math.isfinite(component_denominator)
        else None
    )
    finite_outputs = (
        math.isfinite(offset)
        and math.isfinite(kaon_amplitude)
        and math.isfinite(proton_amplitude)
        and math.isfinite(other_amplitude)
    )
    rejection_reasons = []
    if fit_status_code != 0:
        rejection_reasons.append("fit_status_nonzero")
    if not finite_outputs:
        rejection_reasons.append("nonfinite_offset_fit_outputs")
    if not math.isfinite(offset_error) or offset_error <= 0.0:
        rejection_reasons.append("invalid_offset_error")
    elif offset_error > maximum_offset_error_ns:
        rejection_reasons.append("offset_error_exceeds_max")
    if chi2_ndf is None or not math.isfinite(float(chi2_ndf)):
        rejection_reasons.append("invalid_chi2_ndf")
    elif float(chi2_ndf) > maximum_chi2_ndf:
        rejection_reasons.append("chi2_ndf_exceeds_max")
    if not math.isfinite(kaon_amplitude) or kaon_amplitude <= 0.0:
        rejection_reasons.append("kaon_amplitude_nonpositive")
    if not math.isfinite(proton_amplitude) or proton_amplitude <= 0.0:
        rejection_reasons.append("proton_amplitude_nonpositive")
    if (
        kaon_significance is None
        or not math.isfinite(float(kaon_significance))
        or float(kaon_significance) < minimum_component_significance
    ):
        rejection_reasons.append("kaon_significance_below_min")
    if (
        proton_significance is None
        or not math.isfinite(float(proton_significance))
        or float(proton_significance) < minimum_component_significance
    ):
        rejection_reasons.append("proton_significance_below_min")
    if (
        smaller_component_fraction is None
        or not math.isfinite(float(smaller_component_fraction))
        or float(smaller_component_fraction) < minimum_smaller_component_fraction
    ):
        rejection_reasons.append("smaller_component_fraction_below_min")
    weak_bound_constraint = False
    if bound_hit:
        weak_bound_constraint = (
            (
                math.isfinite(offset_error)
                and offset_error
                > bound_hit_large_error_fraction * maximum_offset_error_ns
            )
            or (chi2_ndf is None or not math.isfinite(float(chi2_ndf)) or float(chi2_ndf) > maximum_chi2_ndf)
            or (
                kaon_significance is None
                or proton_significance is None
                or min(float(kaon_significance), float(proton_significance))
                < minimum_component_significance
            )
        )
    if bound_hit and reject_bound_hit_with_large_error and weak_bound_constraint:
        rejection_reasons.append("bound_hit_with_weak_constraint")
    valid = len(rejection_reasons) == 0
    diagnostic_warnings = []
    if bound_hit:
        diagnostic_warnings.append("bound_hit")
        if weak_bound_constraint:
            diagnostic_warnings.append("bound_hit_with_weak_constraint")
    result = dict(base_result)
    result.update(
        {
            "valid": bool(valid),
            "fit_attempted": True,
            "fit_status": "success" if fit_status_code == 0 else "failure",
            "fit_status_code": int(fit_status_code),
            "delta_offset": float(offset),
            "delta_offset_error": float(offset_error) if math.isfinite(offset_error) else None,
            "delta_offset_bound_hit": bool(bound_hit),
            "diagnostic_warnings": diagnostic_warnings,
            "rejection_reasons": rejection_reasons,
            "rejection_reason": _join_rejection_reasons(rejection_reasons),
            "chi2_data": chi2_data,
            "chi2_ndf": float(chi2_ndf) if chi2_ndf is not None else None,
            "fit_ndf": int(fit_ndf),
            "poisson_deviance": float(poisson_deviance),
            "poisson_deviance_ndf": (
                float(poisson_deviance_ndf)
                if poisson_deviance_ndf is not None
                else None
            ),
            "poisson_deviance_per_entry": (
                float(poisson_deviance_per_entry)
                if poisson_deviance_per_entry is not None
                else None
            ),
            "goodness_ndf": int(poisson_ndf),
            "reference_kaon_mean": float(reference_kaon_mean),
            "reference_proton_mean": float(reference_proton_mean),
            "raw_reference_proton_from_tof": float(raw_proton_reference),
            "wrapped_reference_proton_from_tof": float(effective_proton_reference),
            "rf_period_shift": int(wrap_info.get("period_shift", 0) or 0),
            "kaon_amplitude": float(kaon_amplitude),
            "kaon_amplitude_error": (
                float(kaon_amplitude_error)
                if math.isfinite(kaon_amplitude_error)
                else None
            ),
            "kaon_significance": (
                float(kaon_significance) if kaon_significance is not None else None
            ),
            "proton_amplitude": float(proton_amplitude),
            "proton_amplitude_error": (
                float(proton_amplitude_error)
                if math.isfinite(proton_amplitude_error)
                else None
            ),
            "proton_significance": (
                float(proton_significance) if proton_significance is not None else None
            ),
            "other_amplitude": float(other_amplitude),
            "other_amplitude_error": (
                float(other_amplitude_error)
                if math.isfinite(other_amplitude_error)
                else None
            ),
            "smaller_component_fraction": (
                float(smaller_component_fraction)
                if smaller_component_fraction is not None
                else None
            ),
            "tof_offset_validation": _json_ready_value(validation_cfg),
        }
    )
    return result


def _build_timing_constraint_for_cell(
    global_shape,
    delta_offset_fit,
    tof_summary,
    proton_peak_is_lower,
    probe_kind,
    display_range,
    bunch_spacing_ns,
):
    if not bool((global_shape or {}).get("valid", False)):
        return {"valid": False, "reason": "invalid_global_shape"}
    if not bool((delta_offset_fit or {}).get("valid", False)):
        return {"valid": False, "reason": "invalid_delta_offset_fit"}
    mean_delta_t = (tof_summary or {}).get("mean_delta_t_pk_ns")
    if mean_delta_t is None or not math.isfinite(float(mean_delta_t)) or float(mean_delta_t) <= 0.0:
        return {"valid": False, "reason": "invalid_mean_delta_t_pk"}
    delta_offset = float((delta_offset_fit or {}).get("delta_offset", 0.0) or 0.0)
    reference_kaon_mean = float((global_shape or {}).get("kaon_mean"))
    reference_proton_mean = float((global_shape or {}).get("proton_mean"))
    branch_sign = -1.0 if bool(proton_peak_is_lower) else 1.0
    predicted_kaon_mean = reference_kaon_mean + delta_offset
    raw_predicted_proton_mean = predicted_kaon_mean + (branch_sign * float(mean_delta_t))
    wrap_info = {
        "valid": True,
        "raw_mean": float(raw_predicted_proton_mean),
        "wrapped_mean": float(raw_predicted_proton_mean),
        "period_shift": 0,
        "inside_display_range": True,
        "reference_mean": float(reference_proton_mean),
    }
    if str(probe_kind) == "rf":
        wrap_info = _wrap_rf_mean_to_selected_window(
            raw_predicted_proton_mean,
            reference_proton_mean,
            float(display_range[0]),
            float(display_range[1]),
            float(bunch_spacing_ns),
        )
    if not bool(wrap_info.get("valid", False)):
        return {"valid": False, "reason": wrap_info.get("reason") or "invalid_rf_wrap"}
    return {
        "valid": True,
        "reference_global_kaon_mean": float(reference_kaon_mean),
        "reference_global_proton_mean": float(reference_proton_mean),
        "delta_timing_offset": float(delta_offset),
        "delta_timing_offset_error": (delta_offset_fit or {}).get("delta_offset_error"),
        "mean_delta_t_pk_ns": float(mean_delta_t),
        "predicted_kaon_mean": float(predicted_kaon_mean),
        "predicted_proton_mean_raw": float(raw_predicted_proton_mean),
        "predicted_proton_mean": float(wrap_info["wrapped_mean"]),
        "wrapped_predicted_proton_mean": float(wrap_info["wrapped_mean"]),
        "rf_period_shift": int(wrap_info.get("period_shift", 0) or 0),
        "kaon_sigma": float((global_shape or {}).get("kaon_sigma")),
        "proton_sigma": float((global_shape or {}).get("proton_sigma")),
    }


def _fit_delta_timing_slice(
    histogram,
    global_shape,
    config,
    function_name,
    use_deviance_per_entry_validation=False,
    maximum_poisson_deviance_per_entry=None,
    support_entries=None,
    timing_constraint=None,
):
    fit_min = float(global_shape.get("fit_min", (config.get("ctime_hist_range") or (-1.50, 1.25))[0]))
    fit_max = float(global_shape.get("fit_max", (config.get("ctime_hist_range") or (-1.50, 1.25))[1]))
    minimum_entries = int((config.get("slice_fit") or {}).get("minimum_entries", 30))
    if support_entries is None:
        support_entries = int(round(float(histogram.GetEntries()))) if histogram is not None else 0
    support_snapshot = _build_hist_support_snapshot(
        histogram,
        support_entries,
        minimum_entries,
    )
    invalid_global_shape = not global_shape.get("valid")
    invalid_timing_constraint = bool(timing_constraint is not None) and not bool(
        (timing_constraint or {}).get("valid", False)
    )
    if (
        histogram is None
        or invalid_global_shape
        or invalid_timing_constraint
        or int(support_entries) < int(minimum_entries)
    ):
        rejection_reasons = []
        if invalid_global_shape:
            rejection_reasons.append("invalid_global_shape")
        if invalid_timing_constraint:
            rejection_reasons.append(
                (timing_constraint or {}).get("reason") or "invalid_timing_constraint"
            )
        if histogram is None:
            rejection_reasons.append("missing_histogram")
        elif int(support_entries) < int(minimum_entries):
            rejection_reasons.append("insufficient_entries")
        return {
            "valid": False,
            "fit_attempted": False,
            "fit_status": "insufficient_support",
            "fit_status_code": None,
            "function_name": str(function_name),
            "excluded_invalid_variance_bins": 0,
            "invalid_bin_rule": "macro ROOT fit uses all histogram bins in the fit range",
            "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
            "rejection_reasons": rejection_reasons or ["insufficient_entries_or_invalid_global_shape"],
            "rejection_reason": _join_rejection_reasons(rejection_reasons or ["insufficient_entries_or_invalid_global_shape"]),
            **support_snapshot,
        }
    kaon_mean = float((timing_constraint or {}).get("predicted_kaon_mean", global_shape["kaon_mean"]))
    proton_mean = float(
        (timing_constraint or {}).get(
            "wrapped_predicted_proton_mean",
            (timing_constraint or {}).get("predicted_proton_mean", global_shape["proton_mean"]),
        )
    )
    kaon_sigma = float((timing_constraint or {}).get("kaon_sigma", global_shape["kaon_sigma"]))
    proton_sigma = float((timing_constraint or {}).get("proton_sigma", global_shape["proton_sigma"]))
    histogram_maximum = max(float(histogram.GetMaximum()), 1.0)
    kaon_seed = _find_peak_seed(
        histogram,
        kaon_mean - (2.0 * kaon_sigma),
        kaon_mean + (2.0 * kaon_sigma),
    )
    proton_seed = _find_peak_seed(
        histogram,
        proton_mean - (2.0 * proton_sigma),
        proton_mean + (2.0 * proton_sigma),
    )
    fit_function = ROOT.TF1(
        str(function_name),
        "[0] * exp(-0.5 * pow((x - [1]) / [2], 2))"
        " + "
        "[3] * exp(-0.5 * pow((x - [4]) / [5], 2))"
        " + [6]",
        float(fit_min),
        float(fit_max),
    )
    fit_function.SetParName(0, "K amplitude")
    fit_function.SetParName(1, "K mean")
    fit_function.SetParName(2, "K sigma")
    fit_function.SetParName(3, "p amplitude")
    fit_function.SetParName(4, "p mean")
    fit_function.SetParName(5, "p sigma")
    fit_function.SetParName(6, "other constant")
    fit_function.SetParameter(0, max(float(kaon_seed[0]), 0.10 * histogram_maximum))
    fit_function.SetParLimits(0, 0.0, 100.0 * histogram_maximum)
    fit_function.FixParameter(1, float(kaon_mean))
    fit_function.FixParameter(2, float(kaon_sigma))
    fit_function.SetParameter(3, max(float(proton_seed[0]), 0.05 * histogram_maximum))
    fit_function.SetParLimits(3, 0.0, 100.0 * histogram_maximum)
    fit_function.FixParameter(4, float(proton_mean))
    fit_function.FixParameter(5, float(proton_sigma))
    fit_function.SetParameter(6, 0.02 * histogram_maximum)
    fit_function.SetParLimits(6, 0.0, 10.0 * histogram_maximum)
    fit_result = histogram.Fit(fit_function, PROTON_CLEANING_EXACT_FIT_OPTIONS)
    fit_status_code = int(fit_result)
    covariance_matrix, correlation_matrix, uncertainties = _extract_root_fit_matrices(
        fit_result,
        ("kaon_amplitude", "proton_amplitude", "other_amplitude"),
    )
    kaon_amplitude = float(fit_function.GetParameter(0))
    proton_amplitude = float(fit_function.GetParameter(3))
    other_amplitude = float(fit_function.GetParameter(6))
    kaon_amplitude_error = float(fit_function.GetParError(0))
    proton_amplitude_error = float(fit_function.GetParError(3))
    other_amplitude_error = float(fit_function.GetParError(6))
    kaon_yield = _sum_gaussian_over_bins(
        histogram,
        kaon_amplitude,
        kaon_mean,
        kaon_sigma,
        fit_min,
        fit_max,
    )
    proton_yield = _sum_gaussian_over_bins(
        histogram,
        proton_amplitude,
        proton_mean,
        proton_sigma,
        fit_min,
        fit_max,
    )
    other_yield = _sum_constant_over_bins(
        histogram,
        other_amplitude,
        fit_min,
        fit_max,
    )
    first_fit_bin = max(1, histogram.GetXaxis().FindFixBin(float(fit_min)))
    last_fit_bin = min(
        histogram.GetNbinsX(),
        histogram.GetXaxis().FindFixBin(np.nextafter(float(fit_max), float(fit_min))),
    )
    data_yield = float(histogram.Integral(first_fit_bin, last_fit_bin))
    model_yield = float(kaon_yield + proton_yield + other_yield)
    model_data_ratio = float(model_yield / data_yield) if data_yield > 0.0 else None
    chi2_data = float(fit_function.GetChisquare())
    fit_ndf = int(fit_function.GetNDF())
    chi2_ndf = (
        float(chi2_data / float(fit_ndf))
        if fit_ndf > 0 and math.isfinite(float(chi2_data))
        else None
    )
    chi2_per_abs_entry = (
        float(chi2_data / abs(data_yield))
        if abs(data_yield) > 0.0 and math.isfinite(float(chi2_data))
        else None
    )
    goodness = _compute_poisson_goodness_of_fit(
        histogram,
        fit_function,
        fit_min,
        fit_max,
        3,
    )
    poisson_deviance = float(goodness.get("deviance", 0.0) or 0.0)
    poisson_ndf = int(goodness.get("ndf", 0) or 0)
    poisson_deviance_ndf = goodness.get("deviance_ndf")
    poisson_deviance_per_entry = goodness.get("deviance_per_entry")
    active_bin_count = max(0, int(last_fit_bin - first_fit_bin + 1))
    support_snapshot = _build_hist_support_snapshot(
        histogram,
        support_entries,
        minimum_entries,
        first_fit_bin=first_fit_bin,
        last_fit_bin=last_fit_bin,
    )
    slice_cfg = config.get("slice_fit") or {}
    finite_outputs = (
        math.isfinite(float(kaon_amplitude))
        and math.isfinite(float(proton_amplitude))
        and math.isfinite(float(other_amplitude))
        and math.isfinite(float(kaon_yield))
        and math.isfinite(float(proton_yield))
        and math.isfinite(float(other_yield))
    )
    rejection_reasons = []
    if not finite_outputs:
        rejection_reasons.append("nonfinite_fit_outputs")
    if model_data_ratio is None or not math.isfinite(float(model_data_ratio)):
        rejection_reasons.append("invalid_model_data_ratio")
    else:
        if float(model_data_ratio) < float(slice_cfg.get("minimum_model_data_ratio", 0.50)):
            rejection_reasons.append("model_data_ratio_below_min")
        if float(model_data_ratio) > float(slice_cfg.get("maximum_model_data_ratio", 1.50)):
            rejection_reasons.append("model_data_ratio_above_max")
    slice_fit_status_accepted = bool(
        fit_status_code == 0
        or (
            bool(global_shape.get("per_aero_fallback", False))
            and fit_status_code == 4
            and poisson_deviance_per_entry is not None
            and math.isfinite(float(poisson_deviance_per_entry))
            and float(poisson_deviance_per_entry) <= 0.05
        )
    )
    if not slice_fit_status_accepted:
        rejection_reasons.append("fit_status_{}".format(int(fit_status_code)))
    if (
        poisson_deviance_ndf is None
        or not math.isfinite(float(poisson_deviance_ndf))
    ):
        rejection_reasons.append("invalid_poisson_deviance_ndf")
    if (
        poisson_deviance_per_entry is None
        or not math.isfinite(float(poisson_deviance_per_entry))
    ):
        rejection_reasons.append("invalid_poisson_deviance_per_entry")
    else:
        if bool(use_deviance_per_entry_validation):
            if float(poisson_deviance_per_entry) > float(maximum_poisson_deviance_per_entry):
                rejection_reasons.append("poisson_deviance_per_entry_exceeds_max")
        elif (
            poisson_deviance_ndf is not None
            and math.isfinite(float(poisson_deviance_ndf))
            and float(poisson_deviance_ndf) > float(slice_cfg.get("maximum_poisson_deviance_ndf", 5.0))
        ):
            rejection_reasons.append("poisson_deviance_ndf_exceeds_max")
    if float(kaon_amplitude) < 0.0:
        rejection_reasons.append("negative_kaon_amplitude")
    if float(proton_amplitude) < 0.0:
        rejection_reasons.append("negative_proton_amplitude")
    if float(other_amplitude) < 0.0:
        rejection_reasons.append("negative_other_amplitude")
    if model_yield <= 0.0:
        rejection_reasons.append("nonpositive_model_yield")
    valid = len(rejection_reasons) == 0
    return {
        "valid": bool(valid),
        "fit_attempted": True,
        "fit_status": "success" if fit_status_code == 0 else "failure",
        "fit_status_code": int(fit_status_code),
        "message": "",
        "function_name": str(function_name),
        "reference_global_kaon_mean": (timing_constraint or {}).get("reference_global_kaon_mean", global_shape.get("kaon_mean")),
        "reference_global_proton_mean": (timing_constraint or {}).get("reference_global_proton_mean", global_shape.get("proton_mean")),
        "delta_timing_offset": (timing_constraint or {}).get("delta_timing_offset"),
        "delta_timing_offset_error": (timing_constraint or {}).get("delta_timing_offset_error"),
        "mean_delta_t_pk_ns": (timing_constraint or {}).get("mean_delta_t_pk_ns"),
        "predicted_kaon_mean": float(kaon_mean),
        "predicted_proton_mean": float(proton_mean),
        "wrapped_predicted_proton_mean": float(proton_mean),
        "predicted_proton_mean_raw": (timing_constraint or {}).get("predicted_proton_mean_raw", proton_mean),
        "rf_period_shift": (timing_constraint or {}).get("rf_period_shift", 0),
        "kaon_sigma": float(kaon_sigma),
        "proton_sigma": float(proton_sigma),
        "kaon_amplitude": float(kaon_amplitude),
        "kaon_amplitude_error": float(kaon_amplitude_error),
        "proton_amplitude": float(proton_amplitude),
        "proton_amplitude_error": float(proton_amplitude_error),
        "other_amplitude": float(other_amplitude),
        "other_amplitude_error": float(other_amplitude_error),
        "kaon_yield": float(kaon_yield),
        "proton_yield": float(proton_yield),
        "other_yield": float(other_yield),
        "data_yield": float(data_yield),
        "model_yield": float(model_yield),
        "model_data_ratio": float(model_data_ratio) if model_data_ratio is not None else None,
        "chi2_data": chi2_data,
        "chi2_ndf": float(chi2_ndf) if chi2_ndf is not None else None,
        "chi2_per_abs_entry": float(chi2_per_abs_entry) if chi2_per_abs_entry is not None else None,
        "poisson_deviance": float(poisson_deviance),
        "poisson_deviance_ndf": (
            float(poisson_deviance_ndf)
            if poisson_deviance_ndf is not None
            else None
        ),
        "poisson_deviance_per_entry": (
            float(poisson_deviance_per_entry)
            if poisson_deviance_per_entry is not None
            else None
        ),
        "goodness_ndf": int(poisson_ndf),
        "fitted_entries": float(data_yield),
        "active_bin_count": active_bin_count,
        "excluded_invalid_variance_bins": 0,
        "invalid_bin_rule": "macro ROOT fit uses all histogram bins in the fit range",
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "uncertainties": uncertainties,
        "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
        "rejection_reasons": rejection_reasons,
        "rejection_reason": _join_rejection_reasons(rejection_reasons),
        **support_snapshot,
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
    fixed_bounds = config.get("probe_fixed_bounds")
    if fixed_bounds:
        bounds = {
            "kaon_mean_min": float(fixed_bounds.get("kaonMeanMin")),
            "kaon_mean_max": float(fixed_bounds.get("kaonMeanMax")),
            "proton_mean_min": float(fixed_bounds.get("protonMeanMin")),
            "proton_mean_max": float(fixed_bounds.get("protonMeanMax")),
            "sigma_max": float(fixed_bounds.get("sigmaMax", sigma_max)),
            "sigma_initial": float(fixed_bounds.get("sigmaInitial", global_cfg.get("initial_sigma", 0.15) or 0.15)),
            "peak_pair_found": False,
        }
    else:
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
        use_deviance_per_entry_validation=bool(global_cfg.get("use_deviance_per_entry_validation", False)),
        maximum_poisson_deviance_per_entry=global_cfg.get("maximum_poisson_deviance_per_entry"),
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
    poisson_ndf_values = [
        float(shape.get("poisson_deviance_ndf"))
        for shape in valid_shapes
        if shape.get("poisson_deviance_ndf") is not None
        and math.isfinite(float(shape.get("poisson_deviance_ndf")))
    ]
    poisson_per_entry_values = [
        float(shape.get("poisson_deviance_per_entry"))
        for shape in valid_shapes
        if shape.get("poisson_deviance_per_entry") is not None
        and math.isfinite(float(shape.get("poisson_deviance_per_entry")))
    ]
    mean_poisson_deviance_ndf = (
        float(np.mean(poisson_ndf_values)) if poisson_ndf_values else float("inf")
    )
    mean_poisson_deviance_per_entry = (
        float(np.mean(poisson_per_entry_values)) if poisson_per_entry_values else float("inf")
    )
    peak_pair_found = bool(any(shape.get("peak_pair_found") for shape in active_shapes))
    representative_shape = active_shapes[0] if active_shapes else {}
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
        "meanPoissonDevianceNdf": float(mean_poisson_deviance_ndf),
        "meanPoissonDeviancePerEntry": float(mean_poisson_deviance_per_entry),
        "meanChi2Ndf": float(mean_chi2_ndf),
        "meanChi2PerAbsEntry": float(mean_chi2_per_abs_entry),
        "peakPairFound": bool(peak_pair_found),
        "kaonMeanMin": representative_shape.get("kaonMeanMin"),
        "kaonMeanMax": representative_shape.get("kaonMeanMax"),
        "protonMeanMin": representative_shape.get("protonMeanMin"),
        "protonMeanMax": representative_shape.get("protonMeanMax"),
        "sigmaMin": representative_shape.get("sigmaMin"),
        "sigmaMax": representative_shape.get("sigmaMax"),
        "sigmaInitial": representative_shape.get("sigmaInitial"),
        "minimumGlobalSeparation": representative_shape.get("minimumGlobalSeparation"),
        "boundFractionTolerance": representative_shape.get("boundFractionTolerance"),
        "useDeviancePerEntryValidation": representative_shape.get("useDeviancePerEntryValidation"),
        "maximumPoissonDeviancePerEntry": representative_shape.get("maximumPoissonDeviancePerEntry"),
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
    selection_key="nommcuts",
):
    beam_bunch_spacing_ns = _resolve_beam_bunch_spacing_ns(source_bundle)
    probe_config = deepcopy(config)
    global_cfg = dict(probe_config.get("global_fit") or {})
    global_cfg["beam_bunch_spacing_ns"] = float(beam_bunch_spacing_ns)
    probe_config["global_fit"] = global_cfg
    display_range = None
    rf_discovery_diagnostics = {}
    if str(probe_kind) == "rf":
        display_range, rf_discovery_diagnostics = _resolve_rf_probe_display_range(
            source_bundle,
            evaluate_event,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch,
            selection_key=selection_key,
            source_names=("prompt",),
            return_diagnostics=True,
        )
    else:
        ct_probe_config = _resolve_ct_probe_configuration(source_bundle)
        display_range = tuple(ct_probe_config["display_range"])
        probe_config["global_fit"] = {
            **dict(probe_config.get("global_fit") or {}),
            **dict(ct_probe_config.get("global_fit") or {}),
            "beam_bunch_spacing_ns": float(beam_bunch_spacing_ns),
        }
        probe_config["probe_fixed_bounds"] = {
            "kaonMeanMin": ct_probe_config["kaonMeanMin"],
            "kaonMeanMax": ct_probe_config["kaonMeanMax"],
            "protonMeanMin": ct_probe_config["protonMeanMin"],
            "protonMeanMax": ct_probe_config["protonMeanMax"],
            "sigmaMax": ct_probe_config["sigmaMax"],
            "sigmaInitial": ct_probe_config["sigmaInitial"],
        }
    if str(probe_kind) == "ct":
        time_hist_bins = int(_resolve_ct_probe_configuration(source_bundle)["histogram_bins"])
    else:
        time_hist_bins = _resolve_probe_time_histogram_bins(
            source_bundle,
            probe_kind,
            display_range,
        )
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
        time_hist_bins=time_hist_bins,
        selection_key=selection_key,
    )
    global_shapes = []
    for aero_index, slice_hist in enumerate(pid_payload["global_slice_hists"]):
        shape = _fit_global_timing_shape(
            slice_hist,
            probe_config,
            "f_proton_cleaning_global_{}_aero_{}".format(
                str(timing_branch).replace(" ", "_"),
                aero_index,
            ),
            proton_peak_is_lower=proton_peak_is_lower,
            display_range=display_range,
            fit_mode=fit_mode,
        )
        global_shapes.append(shape)
    summary = _summarize_global_probe(
        timing_branch,
        probe_kind,
        fit_mode,
        pid_payload,
        global_shapes,
        proton_peak_is_lower,
    )
    summary["selection_key"] = str(selection_key)
    if str(probe_kind) == "ct":
        ct_probe_config = _resolve_ct_probe_configuration(source_bundle)
        for key in (
            "kaonMeanMin",
            "kaonMeanMax",
            "protonMeanMin",
            "protonMeanMax",
            "sigmaMin",
            "sigmaMax",
            "sigmaInitial",
            "minimumGlobalSeparation",
            "boundFractionTolerance",
            "useDeviancePerEntryValidation",
            "maximumPoissonDeviancePerEntry",
        ):
            if summary.get(key) is None and key in ct_probe_config:
                summary[key] = ct_probe_config[key]
    for key, value in rf_discovery_diagnostics.items():
        summary[key] = value
    return summary


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
    kaon_mean = float(slice_fit.get("predicted_kaon_mean", global_shape["kaon_mean"]))
    proton_mean = float(
        slice_fit.get(
            "wrapped_predicted_proton_mean",
            slice_fit.get("predicted_proton_mean", global_shape["proton_mean"]),
        )
    )
    kaon_sigma = float(slice_fit.get("kaon_sigma", global_shape["kaon_sigma"]))
    proton_sigma = float(slice_fit.get("proton_sigma", global_shape["proton_sigma"]))
    proton_value = float(
        _gaussian(
            np.asarray([timing], dtype=float),
            slice_fit["proton_amplitude"],
            proton_mean,
            proton_sigma,
        )[0]
    )
    kaon_value = float(
        _gaussian(
            np.asarray([timing], dtype=float),
            slice_fit["kaon_amplitude"],
            kaon_mean,
            kaon_sigma,
        )[0]
    )
    other_value = max(float(slice_fit.get("other_amplitude", 0.0) or 0.0), 0.0)
    denominator = proton_value + kaon_value + other_value
    if (not math.isfinite(denominator)) or denominator <= float(denominator_floor):
        return 0.0
    return max(0.0, min(1.0, float(proton_value / denominator)))


def _build_prepared_event_weight_lookup(cleaning_result, source_bundle):
    if not isinstance(cleaning_result, dict):
        return {}
    prepared_sources = _get_prepared_sources(source_bundle)
    if not prepared_sources:
        return {}

    diagnostics = cleaning_result.get("diagnostics") or {}
    config = cleaning_result.get("settings") or {}
    timing_branch = str(cleaning_result.get("selected_timing_branch") or "CTime_ROC1")
    denominator_floor = float(
        ((config.get("weighting") or {}).get("denominator_floor", 1.0e-12))
    )
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    aero_edges = [float(edge) for edge in (cleaning_result.get("aero_edges") or [])]
    support_by_delta = list(cleaning_result.get("support_by_delta") or [])
    delta_slice_fits = list(cleaning_result.get("delta_slice_fits") or [])
    global_shapes = list(cleaning_result.get("global_shapes") or [])
    selected_selection_key = str(diagnostics.get("selected_timing_selection_key") or "nommcuts")
    lookup = {}
    closure_by_cell = {}
    closure_by_delta = {}
    for delta_index, slice_collection in enumerate(delta_slice_fits):
        for aero_index, slice_fit in enumerate(slice_collection or []):
            key = (int(delta_index), int(aero_index))
            closure_by_cell[key] = {
                "delta_index": int(delta_index),
                "aero_index": int(aero_index),
                "fitted_proton_yield": float((slice_fit or {}).get("proton_yield", 0.0) or 0.0),
                "summed_event_proton_probability": 0.0,
                "closure_ratio": None,
                "event_count": 0,
            }
        closure_by_delta[int(delta_index)] = {
            "delta_index": int(delta_index),
            "fitted_proton_yield": float(
                sum(float((fit or {}).get("proton_yield", 0.0) or 0.0) for fit in (slice_collection or []))
            ),
            "summed_event_proton_probability": 0.0,
            "closure_ratio": None,
            "event_count": 0,
        }

    active_sources = (source_bundle or {}).get("sources") or {}
    for source_name, source_spec in prepared_sources.items():
        tree = (active_sources.get(source_name) or {}).get("tree")
        entry_map = (source_spec or {}).get("entries") or {}
        if tree is None:
            for entry_index, entry_payload in entry_map.items():
                signature = _make_prepared_event_signature(source_name, entry_index)
                lookup[signature] = {
                    "source_label": str(source_name),
                    "source_entry_index": int(entry_index),
                    "delta_index": -1,
                    "aero_index": -1,
                    "support_label": SUPPORT_UNSUPPORTED,
                    "proton_weight": 0.0,
                    "cleaned_factor": 1.0,
                    "rf_accept": True,
                    "final_cleaned_factor": 1.0,
                }
            continue

        for entry_index, evt in enumerate(tree):
            entry_payload = entry_map.get(int(entry_index))
            if entry_payload is None:
                continue
            delta_value = float((entry_payload or {}).get("delta_value", 0.0) or 0.0)
            aero_value = float((entry_payload or {}).get("aero_value", 0.0) or 0.0)
            delta_index = _find_collection_bin(delta_value, delta_edges)
            aero_index = _find_collection_bin(aero_value, aero_edges)
            support_label = SUPPORT_UNSUPPORTED
            proton_weight = 0.0

            if 0 <= delta_index < len(support_by_delta):
                support_label = str(support_by_delta[delta_index])
                if (
                    support_label in (SUPPORT_SUPPORTED, SUPPORT_MARGINAL)
                    and 0 <= aero_index < len(global_shapes)
                    and 0 <= delta_index < len(delta_slice_fits)
                ):
                    slice_collection = delta_slice_fits[delta_index] or []
                    if 0 <= aero_index < len(slice_collection):
                        timing_value = ((entry_payload or {}).get("timing_values") or {}).get(
                            str(timing_branch)
                        )
                        if timing_value is not None:
                            proton_weight = _evaluate_event_proton_probability(
                                float(timing_value),
                                global_shapes[aero_index],
                                slice_collection[aero_index],
                                denominator_floor,
                            )

            proton_weight = max(0.0, min(1.0, float(proton_weight)))
            if bool((entry_payload or {}).get(selected_selection_key, False)):
                fit_coefficient = float((source_spec or {}).get("fit_coefficient", 0.0) or 0.0)
                weighted_probability = float(fit_coefficient * proton_weight)
                cell_key = (int(delta_index), int(aero_index))
                if cell_key in closure_by_cell:
                    closure_by_cell[cell_key]["summed_event_proton_probability"] += weighted_probability
                    closure_by_cell[cell_key]["event_count"] += 1
                if int(delta_index) in closure_by_delta:
                    closure_by_delta[int(delta_index)]["summed_event_proton_probability"] += weighted_probability
                    closure_by_delta[int(delta_index)]["event_count"] += 1
            cleaned_factor = max(0.0, min(1.0, 1.0 - proton_weight))
            rf_accept = apply_low_epsilon_rf_after_proton_cleaning(
                cleaning_result,
                source_name,
                evt,
            )
            final_cleaned_factor = cleaned_factor if rf_accept else 0.0
            signature = _make_prepared_event_signature(source_name, entry_index)
            lookup[signature] = {
                "source_label": str(source_name),
                "source_entry_index": int(entry_index),
                "delta_index": int(delta_index),
                "aero_index": int(aero_index),
                "support_label": str(support_label),
                "proton_weight": float(proton_weight),
                "cleaned_factor": float(cleaned_factor),
                "rf_accept": bool(rf_accept),
                "final_cleaned_factor": float(final_cleaned_factor),
            }

    diagnostics["prepared_event_lookup_count"] = int(len(lookup))
    diagnostics["event_weight_source"] = "setting_wide_immutable_prepared_lookup"
    for row in closure_by_cell.values():
        fitted = float(row.get("fitted_proton_yield", 0.0) or 0.0)
        if fitted > 0.0:
            row["closure_ratio"] = float(row["summed_event_proton_probability"] / fitted)
    for row in closure_by_delta.values():
        fitted = float(row.get("fitted_proton_yield", 0.0) or 0.0)
        if fitted > 0.0:
            row["closure_ratio"] = float(row["summed_event_proton_probability"] / fitted)
    diagnostics["event_weight_closure_by_cell"] = _json_ready_value(list(closure_by_cell.values()))
    diagnostics["event_weight_closure_by_delta"] = _json_ready_value(list(closure_by_delta.values()))
    return lookup


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
    selection_key="nommcuts",
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

    delta_edges = [delta_min + (((delta_max - delta_min) / float(delta_bins)) * i) for i in range(delta_bins + 1)]
    n_aero_slices = max(len(aero_edges) - 1, 0)
    global_prompt_support = [0 for _ in range(n_aero_slices)]
    cell_prompt_support = [
        [0 for _ in range(n_aero_slices)]
        for _ in range(delta_bins)
    ]
    source_stats = {}
    prepared_sources = _get_prepared_sources(source_bundle)
    if prepared_sources:
        prepared_source_stats = (source_bundle.get("prepared_source_stats") or {})
        for source_name, source_spec in prepared_sources.items():
            fit_coefficient = float((source_spec or {}).get("fit_coefficient", 0.0) or 0.0)
            prepared_stats = prepared_source_stats.get(source_name) or {}
            source_stats[source_name] = {
                "tree_name": (source_spec or {}).get("tree_name"),
                "entries_seen": int(prepared_stats.get("entries_seen", 0) or 0),
                "entries_prepared": int(prepared_stats.get("entries_prepared", 0) or 0),
                "entries_passing_nommcuts": 0,
                "entries_passing_selection": 0,
                "entries_missing_CTime_ROC1": 0,
                "entries_missing_timing_branch": 0,
                "entries_missing_selected_timing_branch": 0,
                "timing_branch": str(timing_branch),
                "selection_key": str(selection_key),
                "entries_outside_timing_range": 0,
                "entries_outside_aerogel_range": 0,
                "entries_outside_delta_range": 0,
                "entries_used_in_global_pid": 0,
                "entries_used_in_pid": 0,
                "entries_used": 0,
            }
            if fit_coefficient == 0.0:
                continue
            for entry_payload in (((source_spec or {}).get("entries") or {}).values()):
                if not bool((entry_payload or {}).get(selection_key, False)):
                    continue
                source_stats[source_name]["entries_passing_nommcuts"] += 1
                source_stats[source_name]["entries_passing_selection"] += 1
                aero_value = float((entry_payload or {}).get("aero_value", 0.0) or 0.0)
                timing_value = ((entry_payload or {}).get("timing_values") or {}).get(str(timing_branch))
                if timing_value is None:
                    source_stats[source_name]["entries_missing_CTime_ROC1"] += 1
                    source_stats[source_name]["entries_missing_timing_branch"] += 1
                    source_stats[source_name]["entries_missing_selected_timing_branch"] += 1
                    continue
                timing_value = float(timing_value)
                delta_value = float((entry_payload or {}).get("delta_value", 0.0) or 0.0)
                if (timing_value < time_min) or (timing_value > time_max):
                    source_stats[source_name]["entries_outside_timing_range"] += 1
                    continue
                if (aero_value < aero_min) or (aero_value > aero_max):
                    source_stats[source_name]["entries_outside_aerogel_range"] += 1
                    continue
                aero_index = _find_collection_bin(aero_value, aero_edges)
                if 0 <= aero_index < len(global_slice_hists):
                    source_stats[source_name]["entries_used_in_global_pid"] += 1
                    source_stats[source_name]["entries_used"] += 1
                    h_global_pid.Fill(aero_value, timing_value, fit_coefficient)
                    global_slice_hists[aero_index].Fill(timing_value, fit_coefficient)
                    if bool((source_spec or {}).get("is_prompt_source", False)):
                        global_prompt_support[aero_index] += 1
                delta_index = _find_collection_bin(delta_value, delta_edges)
                if not (0 <= delta_index < delta_bins):
                    source_stats[source_name]["entries_outside_delta_range"] += 1
                    continue
                source_stats[source_name]["entries_used_in_pid"] += 1
                delta_pid_hists[delta_index].Fill(aero_value, timing_value, fit_coefficient)
                if 0 <= aero_index < len(delta_slice_hists[delta_index]):
                    delta_slice_hists[delta_index][aero_index].Fill(timing_value, fit_coefficient)
                    if bool((source_spec or {}).get("is_prompt_source", False)):
                        cell_prompt_support[delta_index][aero_index] += 1
    else:
        prompt_physical_coefficient = float(
            ((((source_bundle or {}).get("sources") or {}).get("prompt") or {}).get("coefficient", 0.0) or 0.0)
        )
        if (not math.isfinite(prompt_physical_coefficient)) or prompt_physical_coefficient == 0.0:
            raise RuntimeError(
                "invalid prompt physical coefficient for proton-cleaning fit scaling: {}".format(
                    prompt_physical_coefficient
                )
            )
        for source_name, source_spec in (source_bundle.get("sources") or {}).items():
            tree = source_spec.get("tree")
            physical_coefficient = float(source_spec.get("coefficient", 0.0) or 0.0)
            fit_coefficient = float(physical_coefficient / prompt_physical_coefficient)
            source_stats[source_name] = {
                "tree_name": source_spec.get("tree_name"),
                "entries_seen": 0,
                "entries_prepared": 0,
                "entries_passing_nommcuts": 0,
                "entries_passing_selection": 0,
                "entries_missing_CTime_ROC1": 0,
                "entries_missing_timing_branch": 0,
                "entries_missing_selected_timing_branch": 0,
                "timing_branch": str(timing_branch),
                "selection_key": str(selection_key),
                "entries_outside_timing_range": 0,
                "entries_outside_aerogel_range": 0,
                "entries_outside_delta_range": 0,
                "entries_used_in_global_pid": 0,
                "entries_used_in_pid": 0,
                "entries_used": 0,
            }
            if tree is None or fit_coefficient == 0.0:
                continue
            for evt in tree:
                source_stats[source_name]["entries_seen"] += 1
                allcuts, nommcuts, _ = evaluate_event(evt, mm_min, mm_max)
                hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) if hole_contains is not None else False
                if not (nommcuts and not hole_rejected):
                    continue
                source_stats[source_name]["entries_passing_nommcuts"] += 1
                source_stats[source_name]["entries_passing_selection"] += 1
                aero_value = float(getattr(evt, "P_aero_npeSum", 0.0))
                try:
                    timing_value = float(getattr(evt, str(timing_branch)))
                except Exception:
                    source_stats[source_name]["entries_missing_CTime_ROC1"] += 1
                    source_stats[source_name]["entries_missing_timing_branch"] += 1
                    source_stats[source_name]["entries_missing_selected_timing_branch"] += 1
                    continue
                delta_value = float(getattr(evt, "ssdelta", 0.0))
                if (timing_value < time_min) or (timing_value > time_max):
                    source_stats[source_name]["entries_outside_timing_range"] += 1
                    continue
                if (aero_value < aero_min) or (aero_value > aero_max):
                    source_stats[source_name]["entries_outside_aerogel_range"] += 1
                    continue
                aero_index = _find_collection_bin(aero_value, aero_edges)
                if 0 <= aero_index < len(global_slice_hists):
                    source_stats[source_name]["entries_used_in_global_pid"] += 1
                    source_stats[source_name]["entries_used"] += 1
                    h_global_pid.Fill(aero_value, timing_value, fit_coefficient)
                    global_slice_hists[aero_index].Fill(timing_value, fit_coefficient)
                    if str(source_name) == "prompt":
                        global_prompt_support[aero_index] += 1
                delta_index = _find_collection_bin(delta_value, delta_edges)
                if not (0 <= delta_index < delta_bins):
                    source_stats[source_name]["entries_outside_delta_range"] += 1
                    continue
                source_stats[source_name]["entries_used_in_pid"] += 1
                delta_pid_hists[delta_index].Fill(aero_value, timing_value, fit_coefficient)
                if 0 <= aero_index < len(delta_slice_hists[delta_index]):
                    delta_slice_hists[delta_index][aero_index].Fill(timing_value, fit_coefficient)
                    if str(source_name) == "prompt":
                        cell_prompt_support[delta_index][aero_index] += 1

    return {
        "H_global_pid": h_global_pid,
        "global_slice_hists": global_slice_hists,
        "delta_pid_hists": delta_pid_hists,
        "delta_slice_hists": delta_slice_hists,
        "source_stats": source_stats,
        "global_prompt_support": global_prompt_support,
        "cell_prompt_support": cell_prompt_support,
        "aero_edges": aero_edges,
        "delta_edges": delta_edges,
        "timing_branch": str(timing_branch),
        "time_hist_range": (float(time_min), float(time_max)),
        "time_hist_bins": int(n_time_bins),
    }


def _resolve_probe_time_histogram_bins(source_bundle, probe_kind, display_range):
    if str(probe_kind) == "rf":
        display_min, display_max = [float(value) for value in display_range]
        candidate_width = max(float(display_max) - float(display_min), 0.0)
        return max(
            80,
            int(round(candidate_width / 0.0305)),
        )
    return int(_resolve_ct_probe_configuration(source_bundle)["histogram_bins"])


def _prepared_selection_has_entries(source_bundle, selection_key):
    for _, _, _, _ in _iter_prepared_records(source_bundle, selection_key=selection_key):
        return True
    return False


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
        "settings": {},
        "diagnostics": {
            "analysis_scope": str(analysis_scope),
            "context": str(context),
            "method": method,
            "debug_flares_enabled": bool(_proton_debug_enabled(config)),
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
                "prepared signed noRF prompt/random/dummy bundle from the upstream "
                "kaon sample immediately before pion subtraction, using prompt-"
                "relative count-scale coefficients for the proton timing-fit histograms"
            ),
            "fit_sample_signed_combination": (
                "norm_data*prompt - norm_data*rand/nWindows - norm_dummy*dummy + "
                "norm_dummy*dummy_rand/nWindows"
            ),
            "fit_sample_uses_signed_random_dummy_subtraction": True,
            "fit_sample_uses_normalized_coefficients": False,
            "fit_sample_uses_prompt_relative_fit_coefficients": True,
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
    implementation = str(
        config.get("implementation")
        or PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT
    ).strip()
    if implementation != PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT:
        raise ValueError(
            "Unsupported proton-cleaning implementation '{}'".format(implementation)
        )
    exact_config = _build_exact_proton_cleaning_config(config)
    result["settings"] = exact_config
    result["diagnostics"]["implementation"] = (
        PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT
    )
    ct_probe_base_configuration = _resolve_ct_probe_configuration(source_bundle)
    result["diagnostics"]["ct_probe_base_configuration"] = {
        "timing_branch": ct_probe_base_configuration["timing_branch"],
        "ctime_hist_range": list(ct_probe_base_configuration["display_range"]),
        "ctime_hist_bins": int(ct_probe_base_configuration["histogram_bins"]),
        "kaonMeanMin": ct_probe_base_configuration["kaonMeanMin"],
        "kaonMeanMax": ct_probe_base_configuration["kaonMeanMax"],
        "protonMeanMin": ct_probe_base_configuration["protonMeanMin"],
        "protonMeanMax": ct_probe_base_configuration["protonMeanMax"],
        "sigmaMin": ct_probe_base_configuration["sigmaMin"],
        "sigmaMax": ct_probe_base_configuration["sigmaMax"],
        "sigmaInitial": ct_probe_base_configuration["sigmaInitial"],
        "aero_slice_edges": list(PROTON_CLEANING_EXACT_AERO_EDGES),
        "delta_hist_range": list(PROTON_CLEANING_EXACT_DELTA_RANGE),
        "delta_bins": int(PROTON_CLEANING_EXACT_DELTA_BINS),
        "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
    }

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

    result["diagnostics"]["physical_source_coefficients"] = {
        str(source_name): float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items()
    }
    result["diagnostics"]["source_coefficients"] = deepcopy(
        result["diagnostics"]["physical_source_coefficients"]
    )
    result["diagnostics"]["prepared_source_stats"] = _json_ready_value(
        (source_bundle or {}).get("prepared_source_stats") or {}
    )
    result["diagnostics"]["available_timing_branches"] = list(
        (source_bundle or {}).get("available_timing_branches") or []
    )
    result["diagnostics"]["fit_source_coefficients"] = {
        str(source_name): float(
            (
                (prepared_source or {}).get("fit_coefficient", (source_spec or {}).get("coefficient", 0.0))
                or 0.0
            )
        )
        for source_name, source_spec in ((source_bundle or {}).get("sources") or {}).items()
        for prepared_source in [
            (((source_bundle or {}).get("prepared_sources") or {}).get(str(source_name)) or {})
        ]
    }
    beam_bunch_spacing_ns = _resolve_beam_bunch_spacing_ns(source_bundle)
    global_fit_cfg = exact_config.get("global_fit") or {}
    global_fit_cfg = dict(global_fit_cfg)
    global_fit_cfg["beam_bunch_spacing_ns"] = float(beam_bunch_spacing_ns)
    exact_config["global_fit"] = global_fit_cfg

    def _summary_without_root_payload(probe):
        if not probe:
            return None
        return _json_ready_value(
            {
                key: value
                for key, value in probe.items()
                if key not in ("pid_payload", "global_shapes")
            }
        )

    def _best_rf_probe(rf_probe_list):
        best = None
        for candidate_probe in rf_probe_list or []:
            if best is None or _compare_timing_probes(candidate_probe, best) > 0:
                best = candidate_probe
        return best

    def _run_probe_set(fit_mode, selection_key="nommcuts"):
        rf_probe_list = []
        for rf_branch in _resolve_rf_branch_candidates(source_bundle, exact_config):
            rf_probe_list.append(
                _run_timing_probe(
                    source_bundle,
                    exact_config,
                    evaluate_event,
                    shifted_mm_getter,
                    hole_contains,
                    mm_min,
                    mm_max,
                    timing_branch=rf_branch,
                    proton_peak_is_lower=True,
                    probe_kind="rf",
                    fit_mode=fit_mode,
                    selection_key=selection_key,
                )
            )
        ct_probe = _run_timing_probe(
            source_bundle,
            exact_config,
            evaluate_event,
            shifted_mm_getter,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch=PROTON_CLEANING_EXACT_TIMING_BRANCH,
            proton_peak_is_lower=False,
            probe_kind="ct",
            fit_mode=fit_mode,
            selection_key=selection_key,
        )
        return rf_probe_list, ct_probe, _best_rf_probe(rf_probe_list)

    selected_probe = None
    selected_stage_prefix = "per_aerogel_multistart_default"
    selected_stage_rf_probes = []
    selected_stage_ct_probe = None
    selected_stage_best_rf_probe = None
    probe_stage_summaries = []
    timing_fit_used_pre_diamond_fallback = False
    timing_fit_used_local_peak_rescue = False

    def _append_probe_stage(fit_mode, reason_prefix, selection_key):
        rf_probes, ct_probe, best_rf_probe = _run_probe_set(
            fit_mode,
            selection_key=selection_key,
        )
        acceptance_mode = (
            "pre_diamond" if str(selection_key) == "pre_diamond_nommcuts" else "inside_diamond"
        )
        stage_summary = {
            "fit_mode": str(fit_mode),
            "acceptance_mode": acceptance_mode,
            "reason_prefix": str(reason_prefix),
            "stage_trigger_reason": str(reason_prefix),
            "selection_key": str(selection_key),
            "rf_probes": [_summary_without_root_payload(probe) for probe in rf_probes],
            "ct_probe": _summary_without_root_payload(ct_probe),
            "best_rf_probe": _summary_without_root_payload(best_rf_probe),
        }
        probe_stage_summaries.append(stage_summary)
        return rf_probes, ct_probe, best_rf_probe, max(
            int((best_rf_probe or {}).get("validShapes", 0) or 0),
            int((ct_probe or {}).get("validShapes", 0) or 0),
        )

    rf_probes, ct_probe, best_rf_probe, best_valid_count = _append_probe_stage(
        "per_aero_multistart",
        "per_aerogel_multistart_default",
        "nommcuts",
    )
    selected_probe = ct_probe
    selected_stage_prefix = "per_aerogel_multistart_default"
    selected_stage_rf_probes = rf_probes
    selected_stage_ct_probe = ct_probe
    selected_stage_best_rf_probe = best_rf_probe

    if best_valid_count <= 0:
        rf_probes, ct_probe, best_rf_probe, best_valid_count = _append_probe_stage(
            "standard",
            "standard_fit_fallback",
            "nommcuts",
        )
        selected_probe = ct_probe
        selected_stage_prefix = "standard_fit_fallback"
        selected_stage_rf_probes = rf_probes
        selected_stage_ct_probe = ct_probe
        selected_stage_best_rf_probe = best_rf_probe

    if best_valid_count <= 0 and _prepared_selection_has_entries(source_bundle, "pre_diamond_nommcuts"):
        timing_fit_used_pre_diamond_fallback = True
        rf_probes, ct_probe, best_rf_probe, best_valid_count = _append_probe_stage(
            "standard",
            "pre_diamond_fit_fallback",
            "pre_diamond_nommcuts",
        )
        selected_probe = ct_probe
        selected_stage_prefix = "pre_diamond_fit_fallback"
        selected_stage_rf_probes = rf_probes
        selected_stage_ct_probe = ct_probe
        selected_stage_best_rf_probe = best_rf_probe

    if timing_fit_used_pre_diamond_fallback and best_valid_count <= 0:
        timing_fit_used_local_peak_rescue = True
        rf_probes, ct_probe, best_rf_probe, best_valid_count = _append_probe_stage(
            "local_peak_rescue",
            "local_peak_rescue",
            "pre_diamond_nommcuts",
        )
        selected_probe = ct_probe
        selected_stage_prefix = "local_peak_rescue"
        selected_stage_rf_probes = rf_probes
        selected_stage_ct_probe = ct_probe
        selected_stage_best_rf_probe = best_rf_probe

    best_rf_probe = selected_stage_best_rf_probe
    ct_probe = selected_stage_ct_probe
    rf_probe_count = len(selected_stage_rf_probes or [])
    local_rescue_active = bool(
        timing_fit_used_local_peak_rescue
        or (best_rf_probe or {}).get("localPeakRescue", False)
        or (ct_probe or {}).get("localPeakRescue", False)
    )
    rescue_rf_peak_pair_tiebreak = bool(
        best_rf_probe
        and local_rescue_active
        and int(best_rf_probe.get("validShapes", 0) or 0) == 0
        and int((ct_probe or {}).get("validShapes", 0) or 0) == 0
        and bool(best_rf_probe.get("peakPairFound", False))
        and not bool((ct_probe or {}).get("peakPairFound", False))
    )
    rf_compare_to_ct = (
        _compare_timing_probes(best_rf_probe, ct_probe)
        if best_rf_probe is not None and ct_probe is not None
        else None
    )
    select_rf = bool(
        best_rf_probe
        and (
            (
                int(best_rf_probe.get("validShapes", 0) or 0) > 0
                and rf_compare_to_ct is not None
                and int(rf_compare_to_ct) >= 0
            )
            or rescue_rf_peak_pair_tiebreak
        )
    )
    if select_rf:
        selected_probe = best_rf_probe
        if rescue_rf_peak_pair_tiebreak:
            timing_selection_reason = "rf_resolved_peak_pair_used_after_all_fit_validation_failed"
        elif int(rf_compare_to_ct or 0) == 0:
            timing_selection_reason = "rf_won_exact_probe_tie"
        else:
            timing_selection_reason = "rf_probe_ranked_better_than_ct"
    else:
        selected_probe = ct_probe
        if not bool(exact_config.get("allow_rf_probe", True)):
            timing_selection_reason = "rf_disabled_ct_selected"
        elif rf_probe_count <= 0:
            timing_selection_reason = "no_rf_branch_ct_selected"
        elif (
            int((best_rf_probe or {}).get("validShapes", 0) or 0) <= 0
            and int((ct_probe or {}).get("validShapes", 0) or 0) <= 0
        ):
            timing_selection_reason = "no_valid_rf_or_ct_shapes_ct_selected_for_diagnostics"
        else:
            timing_selection_reason = "ct_probe_ranked_better_than_rf"
    timing_selection_reason = "{}_{}".format(
        str(selected_stage_prefix),
        str(timing_selection_reason),
    )

    pid_payload = selected_probe.get("pid_payload") or {}
    global_shapes = selected_probe.get("global_shapes") or []
    selected_time_branch = str(
        selected_probe.get("timing_branch")
        or selected_probe.get("branch")
        or PROTON_CLEANING_EXACT_TIMING_BRANCH
    )
    selected_time_hist_range = tuple(pid_payload.get("time_hist_range") or (
        float(selected_probe.get("displayMin", ct_probe_base_configuration["display_range"][0])),
        float(selected_probe.get("displayMax", ct_probe_base_configuration["display_range"][1])),
    ))
    selected_time_hist_bins = int(
        pid_payload.get("time_hist_bins")
        or selected_probe.get("histogramBins")
        or ct_probe_base_configuration["histogram_bins"]
    )
    exact_config["selected_timing_branch"] = selected_time_branch
    exact_config["selected_timing_probe_kind"] = str(selected_probe.get("probe_kind", "ct"))
    exact_config["selected_timing_fit_mode"] = str(selected_probe.get("fit_mode", "unknown"))
    exact_config["selected_proton_peak_is_lower"] = bool(selected_probe.get("proton_peak_is_lower", False))
    exact_config["ctime_hist_range"] = tuple(selected_time_hist_range)
    exact_config["ctime_hist_bins"] = int(selected_time_hist_bins)
    for selected_key in (
        "kaonMeanMin",
        "kaonMeanMax",
        "protonMeanMin",
        "protonMeanMax",
        "sigmaMin",
        "sigmaMax",
        "sigmaInitial",
        "minimumGlobalSeparation",
        "boundFractionTolerance",
        "useDeviancePerEntryValidation",
        "maximumPoissonDeviancePerEntry",
    ):
        if selected_probe.get(selected_key) is not None:
            exact_config[selected_key] = selected_probe.get(selected_key)

    valid_global_shape_count = int(selected_probe.get("validShapes", 0) or 0)

    selected_probe_summary = _summary_without_root_payload(selected_probe)
    result["diagnostics"]["rf_timing_attempted"] = bool(rf_probe_count > 0)
    result["diagnostics"]["ct_timing_evaluated"] = True
    result["diagnostics"]["timingFitUsedPerAeroDefault"] = (
        str(selected_stage_prefix) == "per_aerogel_multistart_default"
    )
    result["diagnostics"]["timingFitUsedStandardFallback"] = (
        str(selected_stage_prefix) == "standard_fit_fallback"
    )
    result["diagnostics"]["timingFitUsedPreDiamondFallback"] = bool(timing_fit_used_pre_diamond_fallback)
    result["diagnostics"]["timingFitUsedLocalPeakRescue"] = bool(timing_fit_used_local_peak_rescue)
    result["diagnostics"]["timing_probe_policy"] = str(
        exact_config.get("timing_probe_policy", "rf_then_ct_best")
    )
    result["diagnostics"]["rf_branch_candidates"] = list(
        _resolve_rf_branch_candidates(source_bundle, exact_config)
    )
    result["diagnostics"]["timing_probe_stage_summaries"] = _json_ready_value(probe_stage_summaries)
    result["diagnostics"]["rf_probe_summaries"] = _json_ready_value(
        [_summary_without_root_payload(probe) for probe in selected_stage_rf_probes]
    )
    result["diagnostics"]["ct_probe_summary"] = _json_ready_value(
        _summary_without_root_payload(selected_stage_ct_probe)
    )
    result["diagnostics"]["best_rf_probe_summary"] = _json_ready_value(
        _summary_without_root_payload(selected_stage_best_rf_probe)
    )
    result["diagnostics"]["selected_probe_summary"] = _json_ready_value(selected_probe_summary)
    result["diagnostics"]["selected_timing_branch"] = selected_time_branch
    result["diagnostics"]["selected_probe_kind"] = str(selected_probe.get("probe_kind", "unknown"))
    result["diagnostics"]["selected_probe_fit_mode"] = str(selected_probe.get("fit_mode", "unknown"))
    result["diagnostics"]["selected_timing_selection_key"] = str(selected_probe.get("selection_key") or "nommcuts")
    result["diagnostics"]["implementation"] = (
        PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT
    )
    result["diagnostics"]["selected_probe_local_peak_rescue"] = bool(selected_probe.get("localPeakRescue", False))
    result["diagnostics"]["selected_probe_per_aero_multistart"] = bool(selected_probe.get("perAeroMultistart", False))
    result["diagnostics"]["selected_proton_peak_is_lower"] = bool(selected_probe.get("proton_peak_is_lower", False))
    result["diagnostics"]["rf_timing_selected"] = bool(str(selected_probe.get("probe_kind", "")) == "rf")
    result["diagnostics"]["rf_compare_to_ct"] = (
        None if rf_compare_to_ct is None else int(rf_compare_to_ct)
    )
    result["diagnostics"]["rf_rescue_peak_pair_tiebreak"] = bool(rescue_rf_peak_pair_tiebreak)
    result["diagnostics"]["timingSelectionReason"] = str(timing_selection_reason)
    result["diagnostics"]["valid_global_shape_count"] = int(valid_global_shape_count)
    result["diagnostics"]["source_stats"] = pid_payload.get("source_stats") or {}
    result["selected_timing_branch"] = selected_time_branch
    result["global_shapes"] = global_shapes
    result["H_global_pid"] = pid_payload.get("H_global_pid")
    result["H_global_timing_slices"] = pid_payload.get("global_slice_hists") or []
    result["H_delta_pid"] = pid_payload.get("delta_pid_hists") or []
    result["H_delta_timing_slices"] = pid_payload.get("delta_slice_hists") or []
    result["aero_edges"] = pid_payload.get("aero_edges") or []
    result["delta_edges"] = pid_payload.get("delta_edges") or []
    result["diagnostics"]["selected_time_hist_range"] = list(selected_time_hist_range)
    result["diagnostics"]["selected_time_hist_bins"] = int(selected_time_hist_bins)
    result["diagnostics"]["selected_kaonMeanMin"] = selected_probe.get("kaonMeanMin")
    result["diagnostics"]["selected_kaonMeanMax"] = selected_probe.get("kaonMeanMax")
    result["diagnostics"]["selected_protonMeanMin"] = selected_probe.get("protonMeanMin")
    result["diagnostics"]["selected_protonMeanMax"] = selected_probe.get("protonMeanMax")
    result["diagnostics"]["selected_sigmaMin"] = selected_probe.get("sigmaMin")
    result["diagnostics"]["selected_sigmaMax"] = selected_probe.get("sigmaMax")
    result["diagnostics"]["selected_sigmaInitial"] = selected_probe.get("sigmaInitial")
    result["diagnostics"]["selected_boundFractionTolerance"] = selected_probe.get("boundFractionTolerance")
    global_shape_debug_rows = []
    for aero_index, shape in enumerate(global_shapes):
        row = {
            "aero_index": int(aero_index),
            "valid": bool((shape or {}).get("valid", False)),
            "fit_status": (shape or {}).get("fit_status"),
            "fit_status_code": (shape or {}).get("fit_status_code"),
            "fit_attempted": bool((shape or {}).get("fit_attempted", False)),
            "support_entries": int((shape or {}).get("support_entries", 0) or 0),
            "minimum_required_entries": int((shape or {}).get("minimum_required_entries", 0) or 0),
            "fitted_entries": float((shape or {}).get("fitted_entries", 0.0) or 0.0),
            "signed_integral": float((shape or {}).get("signed_integral", 0.0) or 0.0),
            "absolute_integral": float((shape or {}).get("absolute_integral", 0.0) or 0.0),
            "root_entries": float((shape or {}).get("root_entries", 0.0) or 0.0),
            "effective_entries": float((shape or {}).get("effective_entries", 0.0) or 0.0),
            "fit_min": (shape or {}).get("fit_min"),
            "fit_max": (shape or {}).get("fit_max"),
            "kaon_mean": (shape or {}).get("kaon_mean"),
            "proton_mean": (shape or {}).get("proton_mean"),
            "kaon_sigma": (shape or {}).get("kaon_sigma"),
            "proton_sigma": (shape or {}).get("proton_sigma"),
            "separation": (shape or {}).get("separation"),
            "kaon_significance": (shape or {}).get("kaon_significance"),
            "proton_significance": (shape or {}).get("proton_significance"),
            "chi2_ndf": (shape or {}).get("chi2_ndf"),
            "chi2_per_abs_entry": (shape or {}).get("chi2_per_abs_entry"),
            "bound_hit": bool((shape or {}).get("bound_hit", False)),
            "ordering_valid": bool((shape or {}).get("ordering_valid", False)),
            "per_aero_fallback": bool((shape or {}).get("per_aero_fallback", False)),
            "peak_pair_found": bool((shape or {}).get("peak_pair_found", False)),
            "fit_options": (shape or {}).get("fit_options"),
            "rejection_reason": (shape or {}).get("rejection_reason"),
        }
        global_shape_debug_rows.append(_json_ready_value(row))
    result["diagnostics"]["global_shape_debug_rows"] = global_shape_debug_rows
    result["diagnostics"]["activeUseSliceDeviancePerEntryValidation"] = False
    result["diagnostics"]["activeMaximumSlicePoissonDeviancePerEntry"] = None

    if valid_global_shape_count <= 0:
        result["fallback_reason"] = "no identifiable proton-kaon timing shapes"
        return result
    selected_selection_key = str(selected_probe.get("selection_key") or "nommcuts")
    delta_edges = [float(edge) for edge in (pid_payload.get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    delta_tof_summaries = _build_prompt_tof_summary_by_delta(
        source_bundle,
        delta_edges,
        selected_selection_key,
        timing_branch=selected_time_branch,
        timing_range=selected_time_hist_range,
        aero_range=exact_config.get("aero_hist_range") or PROTON_CLEANING_EXACT_AERO_RANGE,
        delta_range=exact_config.get("delta_hist_range") or PROTON_CLEANING_EXACT_DELTA_RANGE,
    )
    result["diagnostics"]["tof_summary_by_delta"] = _json_ready_value(delta_tof_summaries)
    result["diagnostics"]["valid_tof_delta_bins"] = int(
        sum(1 for row in delta_tof_summaries if bool((row or {}).get("valid", False)))
    )
    delta_timing_offset_fits = []
    for delta_index, pid_hist in enumerate(pid_payload.get("delta_pid_hists") or []):
        delta_projection = None
        if pid_hist is not None:
            delta_projection = pid_hist.ProjectionY(
                "H_proton_cleaning_delta_offset_time_{}".format(delta_index),
                1,
                int(pid_hist.GetNbinsX()),
            )
            if hasattr(delta_projection, "SetDirectory"):
                delta_projection.SetDirectory(0)
            if hasattr(delta_projection, "Sumw2"):
                delta_projection.Sumw2()
        support_by_aero = (
            (pid_payload.get("cell_prompt_support") or [])[delta_index]
            if delta_index < len(pid_payload.get("cell_prompt_support") or [])
            else []
        )
        reference_shape = _build_reference_shape_for_delta_offset(
            global_shapes,
            support_by_aero,
        )
        offset_fit = _fit_delta_common_timing_offset(
            delta_projection,
            reference_shape,
            delta_tof_summaries[delta_index] if delta_index < len(delta_tof_summaries) else {},
            exact_config,
            "f_proton_cleaning_delta_offset_{}".format(delta_index),
            bool(selected_probe.get("proton_peak_is_lower", False)),
            str(selected_probe.get("probe_kind", "ct")),
            selected_time_hist_range,
            beam_bunch_spacing_ns,
            support_entries=sum(int(value or 0) for value in support_by_aero),
        )
        offset_fit["delta_index"] = int(delta_index)
        delta_timing_offset_fits.append(offset_fit)
    result["delta_timing_offset_fits"] = delta_timing_offset_fits
    result["diagnostics"]["delta_timing_offset_fits"] = _json_ready_value(delta_timing_offset_fits)
    result["diagnostics"]["valid_delta_offset_fits"] = int(
        sum(1 for row in delta_timing_offset_fits if bool((row or {}).get("valid", False)))
    )
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
    delta_support_debug_rows = []
    for delta_index, slice_collection in enumerate(pid_payload["delta_slice_hists"]):
        tof_summary = delta_tof_summaries[delta_index] if delta_index < len(delta_tof_summaries) else {}
        delta_offset_fit = (
            delta_timing_offset_fits[delta_index]
            if delta_index < len(delta_timing_offset_fits)
            else {"valid": False, "rejection_reason": "missing_delta_offset_fit"}
        )
        slice_fits = []
        slice_debug_rows = []
        proton_total = 0.0
        kaon_total = 0.0
        other_total = 0.0
        data_total = float(
            (pid_payload["delta_pid_hists"][delta_index]).Integral()
        )
        fitted_data_total = 0.0
        model_total = 0.0
        valid_slices = 0
        chi2_weighted_sum = 0.0
        chi2_weight = 0.0
        for aero_index, slice_hist in enumerate(slice_collection):
            global_shape = global_shapes[aero_index] if aero_index < len(global_shapes) else {"valid": False}
            timing_constraint = _build_timing_constraint_for_cell(
                global_shape,
                delta_offset_fit,
                tof_summary,
                bool(selected_probe.get("proton_peak_is_lower", False)),
                str(selected_probe.get("probe_kind", "ct")),
                selected_time_hist_range,
                beam_bunch_spacing_ns,
            )
            slice_fit = _fit_delta_timing_slice(
                slice_hist,
                global_shape,
                exact_config,
                "f_proton_cleaning_delta_{}_aero_{}".format(delta_index, aero_index),
                use_deviance_per_entry_validation=False,
                maximum_poisson_deviance_per_entry=None,
                support_entries=(
                    (
                        (pid_payload.get("cell_prompt_support") or [])[delta_index][aero_index]
                    )
                    if delta_index < len(pid_payload.get("cell_prompt_support") or [])
                    and aero_index < len((pid_payload.get("cell_prompt_support") or [])[delta_index])
                    else 0
                ),
                timing_constraint=timing_constraint,
            )
            slice_fits.append(slice_fit)
            slice_debug_rows.append(
                _json_ready_value(
                    {
                        "delta_index": int(delta_index),
                        "aero_index": int(aero_index),
                        "global_shape_valid": bool((global_shape or {}).get("valid", False)),
                        "valid": bool((slice_fit or {}).get("valid", False)),
                        "fit_attempted": bool((slice_fit or {}).get("fit_attempted", False)),
                        "support_entries": int((slice_fit or {}).get("support_entries", 0) or 0),
                        "minimum_required_entries": int((slice_fit or {}).get("minimum_required_entries", 0) or 0),
                        "fit_status": (slice_fit or {}).get("fit_status"),
                        "fit_status_code": (slice_fit or {}).get("fit_status_code"),
                        "data_yield": (slice_fit or {}).get("data_yield"),
                        "model_yield": (slice_fit or {}).get("model_yield"),
                        "model_data_ratio": (slice_fit or {}).get("model_data_ratio"),
                        "signed_integral": (slice_fit or {}).get("signed_integral"),
                        "absolute_integral": (slice_fit or {}).get("absolute_integral"),
                        "root_entries": (slice_fit or {}).get("root_entries"),
                        "effective_entries": (slice_fit or {}).get("effective_entries"),
                        "kaon_yield": (slice_fit or {}).get("kaon_yield"),
                        "proton_yield": (slice_fit or {}).get("proton_yield"),
                        "other_yield": (slice_fit or {}).get("other_yield"),
                        "chi2_ndf": (slice_fit or {}).get("chi2_ndf"),
                        "chi2_per_abs_entry": (slice_fit or {}).get("chi2_per_abs_entry"),
                        "active_bin_count": (slice_fit or {}).get("active_bin_count"),
                        "fit_options": (slice_fit or {}).get("fit_options"),
                        "rejection_reason": (slice_fit or {}).get("rejection_reason"),
                        "delta_timing_offset": (slice_fit or {}).get("delta_timing_offset"),
                        "mean_delta_t_pk_ns": (slice_fit or {}).get("mean_delta_t_pk_ns"),
                        "predicted_kaon_mean": (slice_fit or {}).get("predicted_kaon_mean"),
                        "predicted_proton_mean": (slice_fit or {}).get("predicted_proton_mean"),
                        "rf_period_shift": (slice_fit or {}).get("rf_period_shift"),
                    }
                )
            )
            if not slice_fit.get("valid"):
                continue
            valid_slices += 1
            proton_total += float(slice_fit.get("proton_yield", 0.0) or 0.0)
            kaon_total += float(slice_fit.get("kaon_yield", 0.0) or 0.0)
            other_total += float(slice_fit.get("other_yield", 0.0) or 0.0)
            fitted_data_total += float(slice_fit.get("data_yield", 0.0) or 0.0)
            model_total += float(slice_fit.get("model_yield", 0.0) or 0.0)
            chi2_weighted_sum += (
                float(slice_fit.get("chi2_ndf", 0.0) or 0.0)
                * float(slice_fit.get("data_yield", 0.0) or 0.0)
            )
            chi2_weight += float(slice_fit.get("data_yield", 0.0) or 0.0)
        coverage = float(fitted_data_total / data_total) if data_total > 0.0 else 0.0
        support_thresholds = exact_config.get("support_thresholds") or {}
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
        delta_support_debug_rows.append(
            _json_ready_value(
                {
                    "delta_index": int(delta_index),
                    "support_label": support_label,
                    "valid_slices": int(valid_slices),
                    "coverage": float(coverage),
                    "data_total": float(data_total),
                    "fitted_data_total": float(fitted_data_total),
                    "model_total": float(model_total),
                    "proton_total": float(proton_total),
                    "kaon_total": float(kaon_total),
                    "other_total": float(other_total),
                    "chi2_ndf_weighted": (
                        float(chi2_weighted_sum / chi2_weight) if chi2_weight > 0.0 else None
                    ),
                    "tof_summary": tof_summary,
                    "delta_offset_fit": delta_offset_fit or {},
                    "slice_rows": slice_debug_rows,
                }
            )
        )
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
    result["diagnostics"]["delta_support_debug_rows"] = delta_support_debug_rows
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

    result["_prepared_event_weight_lookup"] = _build_prepared_event_weight_lookup(
        result,
        source_bundle,
    )
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
    h_mm_cleaned = cleaned_targets_pre_rf.get("h_mm_nosub")
    h_mm_cleaned_final_rf = final_targets.get("h_mm_nosub")
    mm_fill_counters = {
        "events_sent_to_mm_filling": 0,
        "raw_mm_fills": 0,
        "proton_mm_fills": 0,
        "final_mm_fills": 0,
        "sum_physical_mm_weights": 0.0,
    }

    rf_counts = {"accepted": 0, "rejected": 0}
    support_counts = {SUPPORT_SUPPORTED: 0, SUPPORT_MARGINAL: 0, SUPPORT_UNSUPPORTED: 0}
    prepared_sources = _get_prepared_sources(source_bundle)

    for source_name, source_spec in (source_bundle.get("sources") or {}).items():
        tree = source_spec.get("tree")
        coefficient = float(source_spec.get("coefficient", 0.0) or 0.0)
        if tree is None or coefficient == 0.0:
            continue
        prepared_entry_map = (prepared_sources.get(source_name) or {}).get("entries") or {}
        for entry_index, evt in enumerate(tree):
            if prepared_sources:
                entry_payload = prepared_entry_map.get(int(entry_index))
                if entry_payload is None:
                    continue
                allcuts = bool((entry_payload or {}).get("allcuts", False))
                nommcuts = bool((entry_payload or {}).get("nommcuts", False))
                noholecuts = bool((entry_payload or {}).get("noholecuts", False))
                adj_hsdelta = float((entry_payload or {}).get("adj_hsdelta", 0.0) or 0.0)
                adj_mm = float((entry_payload or {}).get("adj_mm", 0.0) or 0.0)
                adj_t = float((entry_payload or {}).get("adj_t", 0.0) or 0.0)
                delta_value = float((entry_payload or {}).get("delta_value", 0.0) or 0.0)
                aero_value = float((entry_payload or {}).get("aero_value", 0.0) or 0.0)
                timing_value = ((entry_payload or {}).get("timing_values") or {}).get(str(timing_branch))
                timing_value = float(timing_value) if timing_value is not None else 0.0
            else:
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
            if not (allcuts or nommcuts or noholecuts):
                continue
            frozen_payload = get_kaon_proton_cleaning_event_payload(
                cleaning_result,
                source_name,
                entry_index,
                strict=True,
            )
            delta_index = int(frozen_payload.get("delta_index", -1) or -1)
            aero_index = int(frozen_payload.get("aero_index", -1) or -1)
            support_label = str(
                frozen_payload.get("support_label", SUPPORT_UNSUPPORTED)
                or SUPPORT_UNSUPPORTED
            )
            proton_weight = float(frozen_payload.get("proton_weight", 0.0) or 0.0)
            cleaned_factor = float(frozen_payload.get("cleaned_factor", 1.0) or 0.0)
            rf_accept = bool(frozen_payload.get("rf_accept", True))
            final_cleaned_factor = float(
                frozen_payload.get("final_cleaned_factor", cleaned_factor) or 0.0
            )
            support_counts[support_label] = support_counts.get(support_label, 0) + 1
            raw_weight = float(coefficient)
            proton_component_weight = float(coefficient) * float(proton_weight)
            cleaned_weight = float(coefficient) * float(cleaned_factor)
            final_cleaned_weight = float(coefficient) * float(final_cleaned_factor)
            if allcuts:
                mm_fill_counters["events_sent_to_mm_filling"] += 1
                mm_fill_counters["raw_mm_fills"] += 1
                mm_fill_counters["proton_mm_fills"] += 1
                mm_fill_counters["sum_physical_mm_weights"] += raw_weight
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
            if rf_accept:
                rf_counts["accepted"] += 1
                if allcuts:
                    mm_fill_counters["final_mm_fills"] += 1
                _fill_standard_target_hists(
                    evt,
                    adj_mm,
                    adj_t,
                    adj_hsdelta,
                    final_cleaned_weight,
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
        "H_MM_after_proton_cleaning_final_rf": h_mm_cleaned_final_rf,
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
            "cleaned_final_rf_mm_integral": _hist_integral(h_mm_cleaned_final_rf),
            "mm_validation_window_yields": _build_mm_window_yield_summary(
                h_mm_raw,
                h_mm_proton,
                h_mm_cleaned,
                (config.get("validation_windows") or PROTON_CLEANING_EXACT_VALIDATION_WINDOWS),
            ),
            "mm_fill_counters": dict(mm_fill_counters),
            "raw_mm_key_present": h_mm_raw is not None,
            "proton_mm_key_present": h_mm_proton is not None,
            "final_mm_key_present": h_mm_cleaned is not None,
            "target_template_keys": sorted(str(key) for key in (target_templates or {}).keys()),
            "raw_target_keys": sorted(str(key) for key in (raw_targets or {}).keys()),
            "proton_target_keys": sorted(str(key) for key in (proton_targets or {}).keys()),
            "final_target_keys": sorted(str(key) for key in (final_targets or {}).keys()),
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
    payload.pop("_prepared_event_weight_lookup", None)
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
    prepared_source_stats = diagnostics.get("prepared_source_stats") or {}
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
        "Implementation: {}".format(
            diagnostics.get("implementation", "unknown")
        ),
        "Input tree state: {}".format(diagnostics.get("input_tree_state", "unknown")),
        "Particle selection: {}".format(diagnostics.get("input_tree_particle_selection", "unknown")),
        "Tree policy: {}".format(diagnostics.get("tree_policy", "unknown")),
        "RF policy: {}".format(diagnostics.get("rf_policy", "unknown")),
        "Timing probe policy: {}".format(diagnostics.get("timing_probe_policy", "unknown")),
        "RF probe candidates: {}".format(
            ", ".join(diagnostics.get("rf_branch_candidates") or []) or "none"
        ),
        "Selected timing branch/range/bins: {} {} {}".format(
            diagnostics.get("selected_timing_branch", PROTON_CLEANING_EXACT_TIMING_BRANCH),
            tuple(diagnostics.get("selected_time_hist_range") or PROTON_CLEANING_EXACT_TIMING_RANGE),
            int(diagnostics.get("selected_time_hist_bins", PROTON_CLEANING_EXACT_TIMING_BINS) or PROTON_CLEANING_EXACT_TIMING_BINS),
        ),
        "CT probe branch/range/bins: {} {} {}".format(
            (diagnostics.get("ct_probe_base_configuration") or {}).get(
                "timing_branch",
                PROTON_CLEANING_EXACT_TIMING_BRANCH,
            ),
            tuple(
                (diagnostics.get("ct_probe_base_configuration") or {}).get(
                    "ctime_hist_range",
                    PROTON_CLEANING_EXACT_TIMING_RANGE,
                )
            ),
            int(
                (diagnostics.get("ct_probe_base_configuration") or {}).get(
                    "ctime_hist_bins",
                    PROTON_CLEANING_EXACT_TIMING_BINS,
                )
            ),
        ),
        "Exact aero edges: {}".format(
            tuple(
                ((cleaning_result.get("aero_edges") or ()) or PROTON_CLEANING_EXACT_AERO_EDGES)
            )
        ),
        "Exact delta range/bins: {} {}".format(
            (
                float((cleaning_result.get("delta_edges") or [PROTON_CLEANING_EXACT_DELTA_RANGE[0]])[0]),
                float((cleaning_result.get("delta_edges") or [PROTON_CLEANING_EXACT_DELTA_RANGE[1]])[-1]),
            )
            if cleaning_result.get("delta_edges")
            else PROTON_CLEANING_EXACT_DELTA_RANGE,
            max(len(cleaning_result.get("delta_edges") or []) - 1, PROTON_CLEANING_EXACT_DELTA_BINS),
        ),
        "Fit options: {}".format(PROTON_CLEANING_EXACT_FIT_OPTIONS),
        "Global PID source usage:",
    ]
    for source_name in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
        source_payload = source_stats.get(source_name) or {}
        lines.append(
            "  {}: seen={} used_global={} used_delta={} tree={}".format(
                source_name,
                int(source_payload.get("entries_seen", 0) or 0),
                int(source_payload.get("entries_used_in_global_pid", source_payload.get("entries_used", 0)) or 0),
                int(source_payload.get("entries_used_in_pid", 0) or 0),
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
            "RF timing selected: {}".format("yes" if diagnostics.get("rf_timing_selected") else "no"),
            "Best RF compare-to-CT: {}".format(diagnostics.get("rf_compare_to_ct", "n/a")),
            "Timing selection reason: {}".format(diagnostics.get("timingSelectionReason", "unknown")),
            _format_probe_summary_for_log(
                "Selected probe",
                diagnostics.get("selected_probe_summary"),
            ),
            _format_probe_summary_for_log(
                "Best RF probe",
                diagnostics.get("best_rf_probe_summary"),
            ),
            _format_probe_summary_for_log(
                "CT probe",
                diagnostics.get("ct_probe_summary"),
            ),
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
        mm_fill_counters = app_diag.get("mm_fill_counters") or {}
        lines.extend(
            [
                "MM target keys present raw/proton/final: {}/{}/{}".format(
                    "yes" if app_diag.get("raw_mm_key_present", False) else "no",
                    "yes" if app_diag.get("proton_mm_key_present", False) else "no",
                    "yes" if app_diag.get("final_mm_key_present", False) else "no",
                ),
                "MM fill events raw/proton/final: {}/{}/{}".format(
                    int(mm_fill_counters.get("raw_mm_fills", 0) or 0),
                    int(mm_fill_counters.get("proton_mm_fills", 0) or 0),
                    int(mm_fill_counters.get("final_mm_fills", 0) or 0),
                ),
                "Events sent to MM filling: {}".format(
                    int(mm_fill_counters.get("events_sent_to_mm_filling", 0) or 0)
                ),
                "Sum physical MM weights: {:.3e}".format(
                    float(mm_fill_counters.get("sum_physical_mm_weights", 0.0) or 0.0)
                ),
            ]
        )
    if output_pdf and bool(cleaning_result.get("accepted")):
        lines.append("Diagnostic proton-cleaning plots: {}".format(output_pdf))
    lines.append("========================================")
    print("\n".join(lines))

    if not bool(diagnostics.get("debug_flares_enabled", False)):
        return

    _print_proton_debug(
        "prepared sample overview",
        available_timing_branches=", ".join(diagnostics.get("available_timing_branches") or []) or "none",
        selected_time_hist_range=tuple(diagnostics.get("selected_time_hist_range") or ()),
        selected_time_hist_bins=diagnostics.get("selected_time_hist_bins"),
        physical_source_coefficients=_json_ready_value(diagnostics.get("physical_source_coefficients") or diagnostics.get("source_coefficients") or {}),
        fit_source_coefficients=_json_ready_value(diagnostics.get("fit_source_coefficients") or {}),
    )
    for source_name in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
        prep_stats = prepared_source_stats.get(source_name) or {}
        pid_stats = source_stats.get(source_name) or {}
        _print_proton_debug(
            "prepared source {}".format(source_name),
            tree_name=prep_stats.get("tree_name") or pid_stats.get("tree_name") or "missing",
            timing_branch=pid_stats.get("timing_branch", diagnostics.get("selected_timing_branch", "unknown")),
            selection_key=pid_stats.get("selection_key", "unknown"),
            entries_seen=prep_stats.get("entries_seen", pid_stats.get("entries_seen", 0)),
            entries_prepared=prep_stats.get("entries_prepared", 0),
            entries_passing_nommcuts=pid_stats.get("entries_passing_nommcuts", 0),
            entries_passing_selection=pid_stats.get("entries_passing_selection", 0),
            entries_missing_selected_timing_branch=pid_stats.get("entries_missing_selected_timing_branch", 0),
            entries_outside_timing_range=pid_stats.get("entries_outside_timing_range", 0),
            entries_outside_aerogel_range=pid_stats.get("entries_outside_aerogel_range", 0),
            entries_outside_delta_range=pid_stats.get("entries_outside_delta_range", 0),
            entries_used_in_global_pid=pid_stats.get("entries_used_in_global_pid", pid_stats.get("entries_used", 0)),
            entries_used_in_pid=pid_stats.get("entries_used_in_pid", 0),
        )

    for row in (diagnostics.get("global_shape_debug_rows") or []):
        _print_proton_debug(
            "global slice a{}".format(int(row.get("aero_index", 0)) + 1),
            valid=row.get("valid"),
            fit_status=row.get("fit_status"),
            fit_status_code=row.get("fit_status_code"),
            fit_attempted=row.get("fit_attempted"),
            support_entries=row.get("support_entries"),
            minimum_required_entries=row.get("minimum_required_entries"),
            signed_integral=_format_debug_float(row.get("signed_integral"), scientific=True),
            absolute_integral=_format_debug_float(row.get("absolute_integral"), scientific=True),
            root_entries=_format_debug_float(row.get("root_entries"), scientific=True),
            effective_entries=_format_debug_float(row.get("effective_entries"), scientific=True),
            fitted_entries=_format_debug_float(row.get("fitted_entries"), scientific=True),
            fit_window="[{}, {}]".format(
                _format_debug_float(row.get("fit_min")),
                _format_debug_float(row.get("fit_max")),
            ),
            kaon_mean=_format_debug_float(row.get("kaon_mean")),
            proton_mean=_format_debug_float(row.get("proton_mean")),
            kaon_sigma=_format_debug_float(row.get("kaon_sigma")),
            proton_sigma=_format_debug_float(row.get("proton_sigma")),
            separation=_format_debug_float(row.get("separation")),
            kaon_sig=_format_debug_float(row.get("kaon_significance")),
            proton_sig=_format_debug_float(row.get("proton_significance")),
            chi2_ndf=_format_debug_float(row.get("chi2_ndf")),
            dev_per_entry=_format_debug_float(row.get("chi2_per_abs_entry")),
            bound_hit=row.get("bound_hit"),
            ordering_valid=row.get("ordering_valid"),
            peak_pair_found=row.get("peak_pair_found"),
            rejection_reason=row.get("rejection_reason") or "none",
        )

    for row in (diagnostics.get("delta_support_debug_rows") or []):
        _print_proton_debug(
            "delta bin d{}".format(int(row.get("delta_index", 0)) + 1),
            support=row.get("support_label"),
            valid_slices=row.get("valid_slices"),
            coverage=_format_debug_float(row.get("coverage")),
            data_total=_format_debug_float(row.get("data_total"), scientific=True),
            fitted_data_total=_format_debug_float(row.get("fitted_data_total"), scientific=True),
            model_total=_format_debug_float(row.get("model_total"), scientific=True),
            proton_total=_format_debug_float(row.get("proton_total"), scientific=True),
            kaon_total=_format_debug_float(row.get("kaon_total"), scientific=True),
            other_total=_format_debug_float(row.get("other_total"), scientific=True),
            chi2_ndf_weighted=_format_debug_float(row.get("chi2_ndf_weighted")),
        )
        for slice_row in (row.get("slice_rows") or []):
            _print_proton_debug(
                "delta d{} aero a{}".format(
                    int(slice_row.get("delta_index", 0)) + 1,
                    int(slice_row.get("aero_index", 0)) + 1,
                ),
                global_shape_valid=slice_row.get("global_shape_valid"),
                valid=slice_row.get("valid"),
                fit_attempted=slice_row.get("fit_attempted"),
                support_entries=slice_row.get("support_entries"),
                minimum_required_entries=slice_row.get("minimum_required_entries"),
                fit_status=slice_row.get("fit_status"),
                fit_status_code=slice_row.get("fit_status_code"),
                signed_integral=_format_debug_float(slice_row.get("signed_integral"), scientific=True),
                absolute_integral=_format_debug_float(slice_row.get("absolute_integral"), scientific=True),
                root_entries=_format_debug_float(slice_row.get("root_entries"), scientific=True),
                effective_entries=_format_debug_float(slice_row.get("effective_entries"), scientific=True),
                data_yield=_format_debug_float(slice_row.get("data_yield"), scientific=True),
                model_yield=_format_debug_float(slice_row.get("model_yield"), scientific=True),
                model_ratio=_format_debug_float(slice_row.get("model_data_ratio")),
                kaon_yield=_format_debug_float(slice_row.get("kaon_yield"), scientific=True),
                proton_yield=_format_debug_float(slice_row.get("proton_yield"), scientific=True),
                other_yield=_format_debug_float(slice_row.get("other_yield"), scientific=True),
                chi2_ndf=_format_debug_float(slice_row.get("chi2_ndf")),
                dev_per_entry=_format_debug_float(slice_row.get("chi2_per_abs_entry")),
                active_bin_count=slice_row.get("active_bin_count"),
                rejection_reason=slice_row.get("rejection_reason") or "none",
            )


def _draw_hist_overlay(canvas_name, title, histograms, legend_entries, output_pdf):
    canvas = TCanvas(canvas_name, title, 1000, 700)
    canvas.SetGrid()
    drawn_objects = []
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
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()
    canvas.Modified()
    canvas.Update()
    gc.collect()
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
    if not bool(global_shape.get("fit_attempted", global_shape.get("valid", False))):
        return None
    if not bool(slice_fit.get("fit_attempted", slice_fit.get("valid", False))):
        return None
    kaon_mean = slice_fit.get("predicted_kaon_mean", global_shape.get("kaon_mean"))
    proton_mean = slice_fit.get(
        "wrapped_predicted_proton_mean",
        slice_fit.get("predicted_proton_mean", global_shape.get("proton_mean")),
    )
    kaon_sigma = slice_fit.get("kaon_sigma", global_shape.get("kaon_sigma"))
    proton_sigma = slice_fit.get("proton_sigma", global_shape.get("proton_sigma"))
    if not (
        math.isfinite(float(kaon_mean))
        and math.isfinite(float(proton_mean))
        and math.isfinite(float(kaon_sigma))
        and math.isfinite(float(proton_sigma))
        and math.isfinite(float(slice_fit.get("kaon_amplitude", float("nan"))))
        and math.isfinite(float(slice_fit.get("proton_amplitude", float("nan"))))
        and math.isfinite(float(slice_fit.get("other_amplitude", float("nan"))))
    ):
        return None
    x_values = np.asarray(
        [hist.GetXaxis().GetBinCenter(bin_index) for bin_index in range(1, hist.GetNbinsX() + 1)],
        dtype=float,
    )
    kaon_hist = _clone_hist(hist, "{}_k_overlay".format(hist.GetName()), reset=True)
    proton_hist = _clone_hist(hist, "{}_p_overlay".format(hist.GetName()), reset=True)
    other_hist = _clone_hist(hist, "{}_other_overlay".format(hist.GetName()), reset=True)
    total_hist = _clone_hist(hist, "{}_tot_overlay".format(hist.GetName()), reset=True)
    if kaon_hist is None or proton_hist is None or other_hist is None or total_hist is None:
        return None
    other_value = float(slice_fit.get("other_amplitude", 0.0) or 0.0)
    for bin_index, x_value in enumerate(x_values, start=1):
        kaon_value = float(
            _gaussian(
                np.asarray([x_value]),
                slice_fit["kaon_amplitude"],
                kaon_mean,
                kaon_sigma,
            )[0]
        )
        proton_value = float(
            _gaussian(
                np.asarray([x_value]),
                slice_fit["proton_amplitude"],
                proton_mean,
                proton_sigma,
            )[0]
        )
        total_value = kaon_value + proton_value + other_value
        kaon_hist.SetBinContent(bin_index, kaon_value)
        proton_hist.SetBinContent(bin_index, proton_value)
        other_hist.SetBinContent(bin_index, other_value)
        total_hist.SetBinContent(bin_index, total_value)
    line_style = 1 if (bool(global_shape.get("valid")) and bool(slice_fit.get("valid"))) else 2
    kaon_hist.SetLineColor(kBlue)
    kaon_hist.SetLineStyle(line_style)
    proton_hist.SetLineColor(kRed)
    proton_hist.SetLineStyle(line_style)
    other_hist.SetLineColor(kGray + 2)
    other_hist.SetLineStyle(line_style)
    total_hist.SetLineColor(kGreen + 2)
    total_hist.SetLineStyle(line_style)
    return {
        "kaon": kaon_hist,
        "proton": proton_hist,
        "other": other_hist,
        "total": total_hist,
    }


def _make_delta_axis_hist(name, title, delta_edges):
    edges = [float(edge) for edge in (delta_edges or [])]
    if len(edges) < 2:
        edges = [
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
        ]
    hist = ROOT.TH1D(
        str(name),
        str(title),
        len(edges) - 1,
        array("d", edges),
    )
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def _hist_integral_in_range(hist, x_min, x_max):
    if hist is None:
        return 0.0
    try:
        axis = hist.GetXaxis()
        first_bin = max(1, axis.FindFixBin(float(x_min)))
        last_bin = min(
            int(hist.GetNbinsX()),
            axis.FindFixBin(np.nextafter(float(x_max), float(x_min))),
        )
        if last_bin < first_bin:
            return 0.0
        return float(hist.Integral(first_bin, last_bin))
    except Exception:
        return 0.0


def _hist_integral_error_in_range(hist, x_min, x_max):
    if hist is None:
        return 0.0, 0.0
    try:
        axis = hist.GetXaxis()
        first_bin = max(1, axis.FindFixBin(float(x_min)))
        last_bin = min(
            int(hist.GetNbinsX()),
            axis.FindFixBin(np.nextafter(float(x_max), float(x_min))),
        )
        if last_bin < first_bin:
            return 0.0, 0.0
        total = 0.0
        variance = 0.0
        for bin_index in range(first_bin, last_bin + 1):
            total += float(hist.GetBinContent(bin_index))
            variance += float(hist.GetBinError(bin_index)) ** 2
        return float(total), math.sqrt(max(float(variance), 0.0))
    except Exception:
        return 0.0, 0.0


def _fraction_removed_with_uncertainty(raw_yield, raw_error, removed_yield, removed_error):
    raw_yield = float(raw_yield)
    raw_error = abs(float(raw_error))
    removed_yield = float(removed_yield)
    removed_error = abs(float(removed_error))
    if not math.isfinite(raw_yield) or not math.isfinite(removed_yield):
        return None, None
    if abs(raw_yield) <= max(1.0e-12, raw_error):
        return None, None
    fraction = removed_yield / raw_yield
    variance = (removed_error / raw_yield) ** 2
    variance += ((removed_yield * raw_error) / (raw_yield * raw_yield)) ** 2
    return float(fraction), math.sqrt(max(float(variance), 0.0))


def _build_mm_window_yield_summary(raw_hist, proton_hist, cleaned_hist, validation_windows):
    summaries = {}
    for window_name, bounds in (validation_windows or {}).items():
        if not bounds or len(bounds) < 2:
            continue
        x_min = float(bounds[0])
        x_max = float(bounds[1])
        raw_yield, raw_error = _hist_integral_error_in_range(raw_hist, x_min, x_max)
        proton_yield, proton_error = _hist_integral_error_in_range(proton_hist, x_min, x_max)
        cleaned_yield, cleaned_error = _hist_integral_error_in_range(cleaned_hist, x_min, x_max)
        removed_fraction, removed_fraction_error = _fraction_removed_with_uncertainty(
            raw_yield,
            raw_error,
            proton_yield,
            proton_error,
        )
        summaries[str(window_name)] = {
            "range": [float(x_min), float(x_max)],
            "raw_yield": float(raw_yield),
            "raw_yield_error": float(raw_error),
            "estimated_proton_yield": float(proton_yield),
            "estimated_proton_yield_error": float(proton_error),
            "cleaned_yield": float(cleaned_yield),
            "cleaned_yield_error": float(cleaned_error),
            "removed_fraction": (
                float(removed_fraction) if removed_fraction is not None else None
            ),
            "removed_fraction_error": (
                float(removed_fraction_error)
                if removed_fraction_error is not None
                else None
            ),
        }
    return summaries


def _set_hist_line_marker(hist, color, width=2, style=1, marker=20):
    if hist is None:
        return
    hist.SetLineColor(color)
    hist.SetLineStyle(style)
    hist.SetLineWidth(width)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)


def _probe_fit_goodness(probe_summary):
    if not isinstance(probe_summary, dict):
        return 0.0
    if bool(probe_summary.get("perAeroMultistart", False)) or bool(probe_summary.get("localPeakRescue", False)):
        value = probe_summary.get("meanPoissonDeviancePerEntry", 0.0)
    else:
        value = probe_summary.get("meanPoissonDevianceNdf", 0.0)
    try:
        value = float(value)
    except Exception:
        return 0.0
    return value if math.isfinite(value) else 0.0


def _print_timing_probe_comparison_page(output_pdf, cleaning_result, prefix):
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    rf_probes = [
        probe for probe in (diagnostics.get("rf_probe_summaries") or [])
        if isinstance(probe, dict)
    ]
    ct_probe = diagnostics.get("ct_probe_summary")
    probes = list(rf_probes)
    if isinstance(ct_probe, dict):
        probes.append(ct_probe)
    if not probes:
        return
    page_id = abs(id(cleaning_result))
    n_bins = max(len(probes), 1)
    h_valid = ROOT.TH1D(
        "H_proton_cleaning_probe_valid_{}".format(page_id),
        "{} timing probe valid slice count;probe;valid global slices".format(prefix),
        n_bins,
        0.0,
        float(n_bins),
    )
    h_sep = ROOT.TH1D(
        "H_proton_cleaning_probe_sep_{}".format(page_id),
        "{} timing probe peak separation;probe;mean standardized separation [#sigma]".format(prefix),
        n_bins,
        0.0,
        float(n_bins),
    )
    h_goodness = ROOT.TH1D(
        "H_proton_cleaning_probe_goodness_{}".format(page_id),
        "{} timing probe goodness;probe;ranking goodness".format(prefix),
        n_bins,
        0.0,
        float(n_bins),
    )
    for hist in (h_valid, h_sep, h_goodness):
        hist.SetDirectory(0)
        hist.Sumw2()
        hist.SetLineWidth(3)
        hist.SetMarkerStyle(20)
    h_valid.SetLineColor(kBlue)
    h_valid.SetMarkerColor(kBlue)
    h_sep.SetLineColor(kGreen + 2)
    h_sep.SetMarkerColor(kGreen + 2)
    h_goodness.SetLineColor(kRed)
    h_goodness.SetMarkerColor(kRed)

    selected_branch = str(diagnostics.get("selected_timing_branch", ""))
    for index, probe in enumerate(probes, start=1):
        label = "{}:{}".format(
            str(probe.get("probe_kind", "?")),
            str(probe.get("branch", probe.get("timing_branch", "?"))),
        )
        for hist in (h_valid, h_sep, h_goodness):
            hist.GetXaxis().SetBinLabel(index, label)
        h_valid.SetBinContent(index, float(probe.get("validShapes", 0) or 0))
        h_sep.SetBinContent(index, float(probe.get("meanSeparation", 0.0) or 0.0))
        h_goodness.SetBinContent(index, _probe_fit_goodness(probe))
        if str(probe.get("branch", probe.get("timing_branch", ""))) == selected_branch:
            h_valid.SetBinError(index, 0.15)
            h_sep.SetBinError(index, 0.0)
            h_goodness.SetBinError(index, 0.0)

    canvas = TCanvas(
        "C_proton_cleaning_timing_probe_comparison_{}".format(page_id),
        "{} proton-cleaning timing probe comparison".format(prefix),
        1400,
        900,
    )
    canvas.Divide(2, 2)
    drawn_objects = [h_valid, h_sep, h_goodness]
    canvas.cd(1)
    h_valid.SetMinimum(0.0)
    h_valid.Draw("hist text")
    gPad.Modified()
    gPad.Update()
    canvas.cd(2)
    h_sep.SetMinimum(0.0)
    h_sep.Draw("hist text")
    gPad.Modified()
    gPad.Update()
    canvas.cd(3)
    h_goodness.SetMinimum(0.0)
    h_goodness.Draw("hist text")
    gPad.Modified()
    gPad.Update()
    canvas.cd(4)
    info = TPaveText(0.08, 0.15, 0.92, 0.88, "NDC")
    info.SetBorderSize(1)
    info.SetFillStyle(0)
    info.SetTextAlign(12)
    info.SetTextSize(0.035)
    info.AddText("selected branch: {}".format(diagnostics.get("selected_timing_branch", "unknown")))
    info.AddText("selected kind/mode: {}/{}".format(
        diagnostics.get("selected_probe_kind", "unknown"),
        diagnostics.get("selected_probe_fit_mode", "unknown"),
    ))
    info.AddText("selection reason: {}".format(diagnostics.get("timingSelectionReason", "unknown")))
    info.AddText("RF selected: {}".format("yes" if diagnostics.get("rf_timing_selected") else "no"))
    info.AddText("RF candidates: {}".format(", ".join(diagnostics.get("rf_branch_candidates") or []) or "none"))
    info.Draw()
    drawn_objects.append(info)
    gPad.Modified()
    gPad.Update()
    canvas.Modified()
    canvas.Update()
    gc.collect()
    canvas.Print(output_pdf)


def _print_kaon_proton_cleaning_final_summary_page(output_pdf, cleaning_result, prefix):
    application = cleaning_result.get("application") or {}
    diagnostics = cleaning_result.get("diagnostics") or {}
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    aero_edges = [float(edge) for edge in (cleaning_result.get("aero_edges") or PROTON_CLEANING_EXACT_AERO_EDGES)]
    page_id = abs(id(cleaning_result))
    n_delta_bins = max(len(delta_edges) - 1, 1)
    support_by_delta = [str(label) for label in (cleaning_result.get("support_by_delta") or [])]

    proton_yields = [float(value) for value in (diagnostics.get("proton_yield_by_delta") or [])]
    kaon_yields = [float(value) for value in (diagnostics.get("kaon_yield_by_delta") or [])]
    other_yields = [float(value) for value in (diagnostics.get("other_yield_by_delta") or [])]
    coverage_values = [float(value) for value in (diagnostics.get("valid_coverage_by_delta") or [])]

    h_fit_fraction = _make_delta_axis_hist(
        "H_proton_cleaning_fit_fraction_delta_{}".format(page_id),
        "Integrated fitted proton fraction;SHMS #delta [%];w_{p}^{fit}(#delta)",
        delta_edges,
    )
    h_coverage = _make_delta_axis_hist(
        "H_proton_cleaning_fit_coverage_delta_{}".format(page_id),
        "Accepted PID-fit coverage;SHMS #delta [%];Accepted fit data / total data",
        delta_edges,
    )
    h_proton_yield = _make_delta_axis_hist(
        "H_proton_cleaning_summary_proton_yield_delta_{}".format(page_id),
        "Fitted yields versus SHMS #delta;SHMS #delta [%];Fitted yield",
        delta_edges,
    )
    h_kaon_yield = _make_delta_axis_hist(
        "H_proton_cleaning_summary_kaon_yield_delta_{}".format(page_id),
        "Fitted yields versus SHMS #delta;SHMS #delta [%];Fitted yield",
        delta_edges,
    )
    h_other_yield = _make_delta_axis_hist(
        "H_proton_cleaning_summary_other_yield_delta_{}".format(page_id),
        "Fitted yields versus SHMS #delta;SHMS #delta [%];Fitted yield",
        delta_edges,
    )

    for index in range(n_delta_bins):
        root_bin = index + 1
        proton_yield = proton_yields[index] if index < len(proton_yields) else 0.0
        kaon_yield = kaon_yields[index] if index < len(kaon_yields) else 0.0
        other_yield = other_yields[index] if index < len(other_yields) else 0.0
        total_yield = proton_yield + kaon_yield + other_yield
        h_proton_yield.SetBinContent(root_bin, proton_yield)
        h_kaon_yield.SetBinContent(root_bin, kaon_yield)
        h_other_yield.SetBinContent(root_bin, other_yield)
        if index < len(coverage_values):
            h_coverage.SetBinContent(root_bin, max(0.0, min(1.0, coverage_values[index])))
        support_label = support_by_delta[index] if index < len(support_by_delta) else SUPPORT_UNSUPPORTED
        if support_label != SUPPORT_UNSUPPORTED and total_yield > 0.0:
            proton_fraction = max(0.0, min(1.0, proton_yield / total_yield))
            h_fit_fraction.SetBinContent(root_bin, proton_fraction)
            h_fit_fraction.SetBinError(
                root_bin,
                math.sqrt(max(proton_fraction * (1.0 - proton_fraction) / total_yield, 0.0)),
            )

    h_applied_delta = _clone_hist(
        application.get("H_proton_weight_vs_delta"),
        "H_proton_cleaning_summary_applied_weight_delta_{}".format(page_id),
        reset=False,
    )
    if h_applied_delta is None:
        h_applied_delta = _make_delta_axis_hist(
            "H_proton_cleaning_summary_applied_weight_delta_empty_{}".format(page_id),
            "Mean applied event-level proton weight;SHMS #delta [%];#LTw_{p}^{event}#GT",
            delta_edges,
        )
    h_applied_delta.SetTitle("Mean applied event-level proton weight;SHMS #delta [%];#LTw_{p}^{event}#GT")

    h_applied_map = _clone_hist(
        application.get("H_proton_weight_vs_delta_aero"),
        "H_proton_cleaning_summary_applied_weight_map_{}".format(page_id),
        reset=False,
    )
    if h_applied_map is None:
        h_applied_map = ROOT.TH2D(
            "H_proton_cleaning_summary_applied_weight_map_empty_{}".format(page_id),
            "Mean applied proton probability;SHMS #delta [%];P_aero_npeSum;#LTw_{p}^{event}#GT",
            n_delta_bins,
            array("d", delta_edges),
            max(len(aero_edges) - 1, 1),
            array("d", aero_edges if len(aero_edges) >= 2 else [0.0, 25.0]),
        )
        h_applied_map.SetDirectory(0)
        h_applied_map.Sumw2()
    h_applied_map.SetTitle("Mean applied proton probability;SHMS #delta [%];P_aero_npeSum;#LTw_{p}^{event}#GT")

    h_mm_before = application.get("H_MM_before_proton_cleaning")
    h_mm_proton = application.get("H_MM_estimated_proton")
    h_mm_after = application.get("H_MM_after_proton_cleaning")
    h_mm_fraction = application.get("H_proton_fraction_vs_MM")

    _set_hist_line_marker(h_fit_fraction, kGray + 2, width=2, marker=21)
    _set_hist_line_marker(h_applied_delta, kRed, width=4, marker=20)
    _set_hist_line_marker(h_proton_yield, kRed, width=4, marker=20)
    _set_hist_line_marker(h_kaon_yield, kMagenta + 1, width=4, marker=20)
    _set_hist_line_marker(h_other_yield, kBlue, width=4, marker=20)
    _set_hist_line_marker(h_coverage, kBlue, width=3, marker=20)
    _set_hist_line_marker(h_mm_before, kBlack, width=2)
    _set_hist_line_marker(h_mm_proton, kRed, width=3)
    _set_hist_line_marker(h_mm_after, kGreen + 2, width=3)
    _set_hist_line_marker(h_mm_fraction, kRed, width=4)

    canvas = TCanvas(
        "C_proton_cleaning_final_summary_{}".format(page_id),
        "{} proton-cleaning final summary".format(prefix),
        1800,
        1100,
    )
    canvas.Divide(3, 2)
    drawn_objects = [
        h_fit_fraction,
        h_applied_delta,
        h_proton_yield,
        h_kaon_yield,
        h_other_yield,
        h_coverage,
        h_applied_map,
    ]

    canvas.cd(1)
    h_fit_fraction.SetMinimum(0.0)
    h_fit_fraction.SetMaximum(1.05)
    h_fit_fraction.Draw("E1")
    h_applied_delta.Draw("hist same")
    legend = TLegend(0.48, 0.70, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_fit_fraction, "integrated fitted fraction", "lep")
    legend.AddEntry(h_applied_delta, "mean applied event weight", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(2)
    max_yield = max(
        float(h_proton_yield.GetMaximum()),
        float(h_kaon_yield.GetMaximum()),
        float(h_other_yield.GetMaximum()),
        1.0,
    )
    h_proton_yield.SetMaximum(1.20 * max_yield)
    h_proton_yield.SetMinimum(0.0)
    h_proton_yield.Draw("hist")
    h_kaon_yield.Draw("hist same")
    h_other_yield.Draw("hist same")
    legend = TLegend(0.66, 0.68, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_proton_yield, "proton", "l")
    legend.AddEntry(h_kaon_yield, "kaon", "l")
    legend.AddEntry(h_other_yield, "other", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(3)
    gPad.SetRightMargin(0.16)
    h_applied_map.SetMinimum(0.0)
    h_applied_map.SetMaximum(1.0)
    h_applied_map.SetMarkerSize(0.8)
    h_applied_map.Draw("colz text")
    gPad.Modified()
    gPad.Update()

    canvas.cd(4)
    h_coverage.SetMinimum(0.0)
    h_coverage.SetMaximum(1.05)
    h_coverage.Draw("hist")
    x_min = float(delta_edges[0])
    x_max = float(delta_edges[-1])
    supported_threshold = float(PROTON_CLEANING_EXACT_SUPPORT_THRESHOLDS["minimum_supported_coverage"])
    marginal_threshold = float(PROTON_CLEANING_EXACT_SUPPORT_THRESHOLDS["minimum_marginal_coverage"])
    support_line = TLine(x_min, supported_threshold, x_max, supported_threshold)
    support_line.SetLineColor(kGreen + 2)
    support_line.SetLineStyle(2)
    support_line.SetLineWidth(2)
    support_line.Draw("same")
    marginal_line = TLine(x_min, marginal_threshold, x_max, marginal_threshold)
    marginal_line.SetLineColor(kOrange + 7)
    marginal_line.SetLineStyle(2)
    marginal_line.SetLineWidth(2)
    marginal_line.Draw("same")
    drawn_objects.extend([support_line, marginal_line])
    gPad.Modified()
    gPad.Update()

    canvas.cd(5)
    if h_mm_before is not None and h_mm_proton is not None and h_mm_after is not None:
        max_mm = max(
            float(h_mm_before.GetMaximum()),
            float(h_mm_proton.GetMaximum()),
            float(h_mm_after.GetMaximum()),
            1.0,
        )
        h_mm_before.SetMaximum(1.18 * max_mm)
        h_mm_before.SetMinimum(0.0)
        h_mm_before.Draw("hist")
        h_mm_proton.Draw("hist same")
        h_mm_after.Draw("hist same")
        legend = TLegend(0.50, 0.70, 0.88, 0.88)
        legend.SetBorderSize(1)
        legend.SetFillStyle(0)
        legend.AddEntry(h_mm_before, "raw kaon-selected MM", "l")
        legend.AddEntry(h_mm_proton, "estimated proton contamination", "l")
        legend.AddEntry(h_mm_after, "proton-cleaned kaon MM", "l")
        legend.Draw()
        drawn_objects.append(legend)
        validation_windows = (cleaning_result.get("settings") or {}).get("validation_windows") or {}
        low_bounds = validation_windows.get("low_mm") or PROTON_CLEANING_EXACT_VALIDATION_WINDOWS["low_mm"]
        lambda_bounds = validation_windows.get("lambda_peak") or PROTON_CLEANING_EXACT_VALIDATION_WINDOWS["lambda_peak"]
        low_raw, low_raw_err = _hist_integral_error_in_range(h_mm_before, low_bounds[0], low_bounds[1])
        low_proton, low_proton_err = _hist_integral_error_in_range(h_mm_proton, low_bounds[0], low_bounds[1])
        low_cleaned, low_cleaned_err = _hist_integral_error_in_range(h_mm_after, low_bounds[0], low_bounds[1])
        lambda_raw, lambda_raw_err = _hist_integral_error_in_range(h_mm_before, lambda_bounds[0], lambda_bounds[1])
        lambda_proton, lambda_proton_err = _hist_integral_error_in_range(h_mm_proton, lambda_bounds[0], lambda_bounds[1])
        lambda_cleaned, lambda_cleaned_err = _hist_integral_error_in_range(h_mm_after, lambda_bounds[0], lambda_bounds[1])
        low_removed, low_removed_err = _fraction_removed_with_uncertainty(
            low_raw, low_raw_err, low_proton, low_proton_err
        )
        lambda_removed, lambda_removed_err = _fraction_removed_with_uncertainty(
            lambda_raw, lambda_raw_err, lambda_proton, lambda_proton_err
        )
        pave = TPaveText(0.43, 0.31, 0.88, 0.64, "NDC")
        pave.SetBorderSize(1)
        pave.SetFillStyle(0)
        pave.SetTextAlign(12)
        pave.SetTextSize(0.026)
        pave.AddText("Low-MM raw/proton/clean = {:.3g}+/-{:.2g} / {:.3g}+/-{:.2g} / {:.3g}+/-{:.2g}".format(
            low_raw,
            low_raw_err,
            low_proton,
            low_proton_err,
            low_cleaned,
            low_cleaned_err,
        ))
        if low_removed is not None:
            pave.AddText("Low-MM removed = {:.1f}+/-{:.1f}%".format(100.0 * low_removed, 100.0 * low_removed_err))
        else:
            pave.AddText("Low-MM removed = unavailable")
        pave.AddText("#Lambda raw/proton/clean = {:.3g}+/-{:.2g} / {:.3g}+/-{:.2g} / {:.3g}+/-{:.2g}".format(
            lambda_raw,
            lambda_raw_err,
            lambda_proton,
            lambda_proton_err,
            lambda_cleaned,
            lambda_cleaned_err,
        ))
        if lambda_removed is not None:
            pave.AddText("#Lambda removed = {:.1f}+/-{:.1f}%".format(100.0 * lambda_removed, 100.0 * lambda_removed_err))
        else:
            pave.AddText("#Lambda removed = unavailable")
        pave.Draw()
        drawn_objects.append(pave)
        drawn_objects.extend([h_mm_before, h_mm_proton, h_mm_after])
    gPad.Modified()
    gPad.Update()

    canvas.cd(6)
    if h_mm_fraction is not None:
        h_mm_fraction.SetMinimum(0.0)
        h_mm_fraction.SetMaximum(max(1.0, 1.15 * float(h_mm_fraction.GetMaximum())))
        h_mm_fraction.SetTitle("Estimated proton fraction versus MM;MM [GeV];Estimated proton / raw")
        h_mm_fraction.Draw("hist")
        drawn_objects.append(h_mm_fraction)
    gPad.Modified()
    gPad.Update()

    canvas.Modified()
    canvas.Update()
    gc.collect()
    canvas.Print(output_pdf)


def _finite_float_or_none(value):
    if value is None:
        return None
    try:
        numeric = float(value)
    except Exception:
        return None
    return numeric if math.isfinite(numeric) else None


def _set_hist_bin_if_finite(hist, bin_index, value, error=None):
    numeric = _finite_float_or_none(value)
    if hist is None or numeric is None:
        return False
    hist.SetBinContent(int(bin_index), float(numeric))
    err_value = _finite_float_or_none(error)
    if err_value is not None:
        hist.SetBinError(int(bin_index), abs(float(err_value)))
    return True


def _print_proton_tof_constraint_diagnostics_page(output_pdf, cleaning_result, prefix):
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    tof_summaries = diagnostics.get("tof_summary_by_delta") or []
    offset_fits = (cleaning_result or {}).get("delta_timing_offset_fits") or diagnostics.get("delta_timing_offset_fits") or []
    if not tof_summaries and not offset_fits:
        return

    delta_edges = [float(edge) for edge in ((cleaning_result or {}).get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    page_id = abs(id(cleaning_result))
    h_p = _make_delta_axis_hist(
        "H_proton_cleaning_tof_mean_p_delta_{}".format(page_id),
        "Prompt SHMS momentum constraint;SHMS #delta [%];#LTp_{SHMS}#GT [GeV]",
        delta_edges,
    )
    h_path = _make_delta_axis_hist(
        "H_proton_cleaning_tof_path_delta_{}".format(page_id),
        "SHMS path-length constraint;SHMS #delta [%];#LTL#GT [cm]",
        delta_edges,
    )
    h_dt = _make_delta_axis_hist(
        "H_proton_cleaning_tof_deltat_delta_{}".format(page_id),
        "Predicted p-K TOF separation;SHMS #delta [%];#LT#Delta t_{pK}#GT [ns]",
        delta_edges,
    )
    h_offset = _make_delta_axis_hist(
        "H_proton_cleaning_tof_offset_delta_{}".format(page_id),
        "Fitted common timing offset;SHMS #delta [%];offset [ns]",
        delta_edges,
    )
    h_offset_invalid = _make_delta_axis_hist(
        "H_proton_cleaning_tof_offset_invalid_delta_{}".format(page_id),
        "Fitted common timing offset;SHMS #delta [%];offset [ns]",
        delta_edges,
    )
    h_offset_bound = _make_delta_axis_hist(
        "H_proton_cleaning_tof_offset_bound_delta_{}".format(page_id),
        "Fitted common timing offset;SHMS #delta [%];offset [ns]",
        delta_edges,
    )
    h_k_ref = _make_delta_axis_hist(
        "H_proton_cleaning_tof_k_ref_delta_{}".format(page_id),
        "Corrected timing centers;SHMS #delta [%];timing center [ns]",
        delta_edges,
    )
    h_p_ref = _make_delta_axis_hist(
        "H_proton_cleaning_tof_p_ref_delta_{}".format(page_id),
        "Corrected timing centers;SHMS #delta [%];timing center [ns]",
        delta_edges,
    )
    h_k_global = _make_delta_axis_hist(
        "H_proton_cleaning_tof_k_global_delta_{}".format(page_id),
        "Corrected timing centers;SHMS #delta [%];timing center [ns]",
        delta_edges,
    )
    h_p_global = _make_delta_axis_hist(
        "H_proton_cleaning_tof_p_global_delta_{}".format(page_id),
        "Corrected timing centers;SHMS #delta [%];timing center [ns]",
        delta_edges,
    )
    for row in tof_summaries:
        if not isinstance(row, dict):
            continue
        bin_index = int(row.get("delta_index", -1)) + 1
        if bin_index < 1 or bin_index > h_p.GetNbinsX():
            continue
        _set_hist_bin_if_finite(h_p, bin_index, row.get("mean_P_gtr_p"), row.get("rms_P_gtr_p"))
        _set_hist_bin_if_finite(h_path, bin_index, row.get("mean_shms_path_length_cm"), row.get("rms_shms_path_length_cm"))
        _set_hist_bin_if_finite(h_dt, bin_index, row.get("mean_delta_t_pk_ns"), row.get("rms_delta_t_pk_ns"))

    for row in offset_fits:
        if not isinstance(row, dict):
            continue
        bin_index = int(row.get("delta_index", -1)) + 1
        if bin_index < 1 or bin_index > h_offset.GetNbinsX():
            continue
        if bool(row.get("valid", False)):
            _set_hist_bin_if_finite(h_offset, bin_index, row.get("delta_offset"), row.get("delta_offset_error"))
        else:
            _set_hist_bin_if_finite(h_offset_invalid, bin_index, row.get("delta_offset"), row.get("delta_offset_error"))
        if bool(row.get("delta_offset_bound_hit", False)):
            _set_hist_bin_if_finite(h_offset_bound, bin_index, row.get("delta_offset"), row.get("delta_offset_error"))
        reference_k = _finite_float_or_none(row.get("reference_kaon_mean"))
        reference_p = _finite_float_or_none(row.get("reference_proton_mean"))
        offset = _finite_float_or_none(row.get("delta_offset"))
        wrapped_p = _finite_float_or_none(row.get("wrapped_reference_proton_from_tof"))
        if reference_k is not None:
            _set_hist_bin_if_finite(h_k_global, bin_index, reference_k)
        if reference_p is not None:
            _set_hist_bin_if_finite(h_p_global, bin_index, reference_p)
        if reference_k is not None and offset is not None:
            _set_hist_bin_if_finite(h_k_ref, bin_index, reference_k + offset, row.get("delta_offset_error"))
        if wrapped_p is not None and offset is not None:
            _set_hist_bin_if_finite(h_p_ref, bin_index, wrapped_p + offset, row.get("delta_offset_error"))

    _set_hist_line_marker(h_p, kBlue, width=2, marker=20)
    _set_hist_line_marker(h_path, kGreen + 2, width=2, marker=20)
    _set_hist_line_marker(h_dt, kMagenta + 1, width=2, marker=20)
    _set_hist_line_marker(h_offset, kRed, width=2, marker=20)
    _set_hist_line_marker(h_offset_invalid, kRed + 2, width=2, style=2, marker=24)
    _set_hist_line_marker(h_offset_bound, kOrange + 7, width=2, style=1, marker=25)
    _set_hist_line_marker(h_k_global, kBlue, width=2, style=2, marker=24)
    _set_hist_line_marker(h_p_global, kRed, width=2, style=2, marker=25)
    _set_hist_line_marker(h_k_ref, kBlue, width=3, marker=20)
    _set_hist_line_marker(h_p_ref, kRed, width=3, marker=21)

    canvas = TCanvas(
        "C_proton_cleaning_tof_constraints_{}".format(page_id),
        "{} proton-cleaning SHMS TOF constraints".format(prefix),
        1800,
        1100,
    )
    canvas.Divide(3, 2)
    drawn_objects = [
        h_p,
        h_path,
        h_dt,
        h_offset,
        h_offset_invalid,
        h_offset_bound,
        h_k_global,
        h_p_global,
        h_k_ref,
        h_p_ref,
    ]

    canvas.cd(1)
    h_p.Draw("E1")
    gPad.Modified()
    gPad.Update()

    canvas.cd(2)
    h_path.Draw("E1")
    gPad.Modified()
    gPad.Update()

    canvas.cd(3)
    h_dt.Draw("E1")
    gPad.Modified()
    gPad.Update()

    canvas.cd(4)
    offset_abs_max = max(
        abs(float(h_offset.GetMaximum())),
        abs(float(h_offset.GetMinimum())),
        abs(float(h_offset_invalid.GetMaximum())),
        abs(float(h_offset_invalid.GetMinimum())),
        abs(float(h_offset_bound.GetMaximum())),
        abs(float(h_offset_bound.GetMinimum())),
        0.05,
    )
    h_offset.SetMinimum(-1.15 * offset_abs_max)
    h_offset.SetMaximum(1.15 * offset_abs_max)
    h_offset.Draw("E1")
    h_offset_invalid.Draw("E1 same")
    h_offset_bound.Draw("E1 same")
    zero_line = TLine(float(delta_edges[0]), 0.0, float(delta_edges[-1]), 0.0)
    zero_line.SetLineColor(kGray + 2)
    zero_line.SetLineStyle(2)
    zero_line.Draw("same")
    drawn_objects.append(zero_line)
    valid_offsets = sum(1 for row in offset_fits if bool((row or {}).get("valid", False)))
    invalid_offsets = sum(1 for row in offset_fits if not bool((row or {}).get("valid", False)))
    bound_offsets = sum(1 for row in offset_fits if bool((row or {}).get("delta_offset_bound_hit", False)))
    offset_info = TPaveText(0.14, 0.73, 0.52, 0.90, "NDC")
    offset_info.SetBorderSize(1)
    offset_info.SetFillStyle(0)
    offset_info.SetTextAlign(12)
    offset_info.SetTextSize(0.032)
    offset_info.AddText("valid/invalid: {} / {}".format(valid_offsets, invalid_offsets))
    offset_info.AddText("bound-hit warnings: {}".format(bound_offsets))
    offset_info.Draw()
    drawn_objects.append(offset_info)
    gPad.Modified()
    gPad.Update()

    canvas.cd(5)
    max_center = max(
        float(h_k_global.GetMaximum()),
        float(h_p_global.GetMaximum()),
        float(h_k_ref.GetMaximum()),
        float(h_p_ref.GetMaximum()),
        1.0e-6,
    )
    min_center = min(
        float(h_k_global.GetMinimum()),
        float(h_p_global.GetMinimum()),
        float(h_k_ref.GetMinimum()),
        float(h_p_ref.GetMinimum()),
        -1.0e-6,
    )
    if max_center <= min_center:
        max_center = min_center + 1.0
    h_k_ref.SetMinimum(min_center - 0.10 * abs(max_center - min_center))
    h_k_ref.SetMaximum(max_center + 0.15 * abs(max_center - min_center))
    h_k_ref.Draw("E1")
    h_p_ref.Draw("E1 same")
    h_k_global.Draw("hist same")
    h_p_global.Draw("hist same")
    legend = TLegend(0.50, 0.62, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_k_ref, "corrected K center", "lep")
    legend.AddEntry(h_p_ref, "corrected p center", "lep")
    legend.AddEntry(h_k_global, "global K ref", "l")
    legend.AddEntry(h_p_global, "global p ref", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(6)
    closure_rows = []
    for row in diagnostics.get("event_weight_closure_by_delta") or []:
        if not isinstance(row, dict):
            continue
        ratio = _finite_float_or_none(row.get("closure_ratio"))
        delta_index = int(row.get("delta_index", -1))
        if ratio is None or not (0 <= delta_index < len(delta_edges) - 1):
            continue
        closure_rows.append((delta_index, ratio))
    graph = None
    if closure_rows:
        x_values = array(
            "d",
            [
                0.5 * (float(delta_edges[index]) + float(delta_edges[index + 1]))
                for index, _ in closure_rows
            ],
        )
        y_values = array("d", [float(value) for _, value in closure_rows])
        x_errors = array(
            "d",
            [
                0.5 * abs(float(delta_edges[index + 1]) - float(delta_edges[index]))
                for index, _ in closure_rows
            ],
        )
        y_errors = array("d", [0.0 for _ in closure_rows])
        graph = ROOT.TGraphErrors(len(closure_rows), x_values, y_values, x_errors, y_errors)
        graph.SetTitle("Event-weight closure by #delta;SHMS #delta [%];#Sigma w_{p}^{event} / fitted p yield")
        graph.SetLineColor(kViolet + 1)
        graph.SetMarkerColor(kViolet + 1)
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(3)
        graph.SetMinimum(0.0)
        graph.SetMaximum(max(1.25, 1.15 * max(float(value) for _, value in closure_rows)))
        graph.Draw("AP")
        drawn_objects.append(graph)
    else:
        empty_axis = _make_delta_axis_hist(
            "H_proton_cleaning_event_weight_closure_empty_{}".format(page_id),
            "Event-weight closure by #delta;SHMS #delta [%];#Sigma w_{p}^{event} / fitted p yield",
            delta_edges,
        )
        empty_axis.SetMinimum(0.0)
        empty_axis.SetMaximum(1.25)
        empty_axis.Draw("axis")
        drawn_objects.append(empty_axis)
    unity_line = TLine(float(delta_edges[0]), 1.0, float(delta_edges[-1]), 1.0)
    unity_line.SetLineColor(kGray + 2)
    unity_line.SetLineStyle(2)
    unity_line.Draw("same")
    drawn_objects.append(unity_line)
    closure_info = TPaveText(0.14, 0.78, 0.48, 0.90, "NDC")
    closure_info.SetBorderSize(1)
    closure_info.SetFillStyle(0)
    closure_info.SetTextAlign(12)
    closure_info.SetTextSize(0.035)
    closure_info.AddText("valid closure bins: {} / {}".format(len(closure_rows), max(len(delta_edges) - 1, 0)))
    closure_info.Draw()
    drawn_objects.append(closure_info)
    gPad.Modified()
    gPad.Update()

    canvas.Modified()
    canvas.Update()
    gc.collect()
    canvas.Print(output_pdf)


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

    _print_timing_probe_comparison_page(output_pdf, cleaning_result, prefix)
    _print_proton_tof_constraint_diagnostics_page(output_pdf, cleaning_result, prefix)

    h_global_pid = cleaning_result.get("H_global_pid")
    if h_global_pid is not None:
        canvas = TCanvas("C_proton_cleaning_global_pid", "{} proton-cleaning global PID".format(prefix), 1000, 700)
        drawn_objects = []
        h_global_pid.Draw("colz")
        drawn_objects.append(_draw_status_pave(cleaning_result))
        gPad.Modified()
        gPad.Update()
        canvas.Modified()
        canvas.Update()
        gc.collect()
        canvas.Print(output_pdf)

    global_slice_hists = cleaning_result.get("H_global_timing_slices") or []
    if global_slice_hists:
        canvas = TCanvas("C_proton_cleaning_global_slices", "{} proton-cleaning global timing slices".format(prefix), 1200, 800)
        canvas.Divide(3, 2)
        drawn_objects = []
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
                        "fit_attempted": shape.get("fit_attempted", shape.get("valid")),
                        "valid": shape.get("valid"),
                        "kaon_amplitude": shape.get("kaon_amplitude", 0.0),
                        "proton_amplitude": shape.get("proton_amplitude", 0.0),
                        "other_amplitude": shape.get("other_amplitude", 0.0),
                    },
                )
                if overlays is not None:
                    for overlay in overlays.values():
                        overlay.SetLineWidth(2)
                        overlay.Draw("hist same")
                        drawn_objects.append(overlay)
                elif shape.get("fit_attempted", shape.get("valid")):
                    _print_proton_debug(
                        "overlay construction failed",
                        slice_index=index,
                        fit_attempted=shape.get("fit_attempted"),
                        valid=shape.get("valid"),
                        kaon_mean=shape.get("kaon_mean"),
                        proton_mean=shape.get("proton_mean"),
                        kaon_sigma=shape.get("kaon_sigma"),
                        proton_sigma=shape.get("proton_sigma"),
                        kaon_amplitude=shape.get("kaon_amplitude"),
                        proton_amplitude=shape.get("proton_amplitude"),
                        other_amplitude=shape.get("other_amplitude"),
                    )
                drawn_objects.append(
                    _draw_status_pave(
                        cleaning_result,
                        extra_lines=(
                            "slice {} attempted={} valid={}".format(
                                index + 1,
                                bool(shape.get("fit_attempted", shape.get("valid"))),
                                bool(shape.get("valid")),
                            ),
	                            "status={}".format(shape.get("fit_status")),
	                            "entries={}".format(shape.get("support_entries", "n/a")),
	                            "chi2/ndf={}".format(shape.get("chi2_ndf")),
	                            "warn={}".format(
	                                ", ".join(str(w) for w in (shape.get("diagnostic_warnings") or [])) or "none"
	                            ),
	                            "reason={}".format(shape.get("rejection_reason") or "none"),
	                        ),
	                        x1=0.14,
                        y1=0.68,
                        x2=0.62,
                        y2=0.90,
                    )
                )
            gPad.Modified()
            gPad.Update()
        canvas.Modified()
        canvas.Update()
        gc.collect()
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
    drawn_objects = [h_support, h_proton, h_kaon, h_other]
    canvas.cd(1)
    h_support.SetMinimum(-0.2)
    h_support.SetMaximum(2.2)
    h_support.Draw("hist")
    gPad.Modified()
    gPad.Update()
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
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()
    canvas.cd(3)
    if application.get("H_proton_weight_vs_delta") is not None:
        application["H_proton_weight_vs_delta"].SetLineColor(kViolet + 1)
        application["H_proton_weight_vs_delta"].Draw("hist")
        drawn_objects.append(application["H_proton_weight_vs_delta"])
    gPad.Modified()
    gPad.Update()
    canvas.cd(4)
    if application.get("H_proton_weight_vs_delta_aero") is not None:
        application["H_proton_weight_vs_delta_aero"].Draw("colz")
        drawn_objects.append(application["H_proton_weight_vs_delta_aero"])
    else:
        first_delta_pid = next((hist for hist in (cleaning_result.get("H_delta_pid") or []) if hist is not None), None)
        if first_delta_pid is not None:
            first_delta_pid.Draw("colz")
            drawn_objects.append(first_delta_pid)
    gPad.Modified()
    gPad.Update()
    canvas.cd(1)
    drawn_objects.append(
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
    )
    gPad.Modified()
    gPad.Update()
    canvas.Modified()
    canvas.Update()
    gc.collect()
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
        drawn_objects = []
        canvas.cd(1)
        if pid_hist is not None:
            pid_hist.Draw("colz")
            drawn_objects.append(pid_hist)
        support_label = support_by_delta[delta_index] if delta_index < len(support_by_delta) else "unknown"
        drawn_objects.append(
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
        )
        gPad.Modified()
        gPad.Update()
        for slice_index, hist in enumerate(slice_collection[:5], start=2):
            canvas.cd(slice_index)
            if hist is None:
                continue
            hist.SetLineColor(kBlack)
            hist.Draw("hist")
            drawn_objects.append(hist)
            global_shape = global_shapes[slice_index - 2] if (slice_index - 2) < len(global_shapes) else {}
            slice_fit = slice_fit_collection[slice_index - 2] if (slice_index - 2) < len(slice_fit_collection) else {}
            overlays = _build_timing_shape_overlays(hist, global_shape, slice_fit)
            if overlays is not None:
                for overlay in overlays.values():
                    overlay.SetLineWidth(2)
                    overlay.Draw("hist same")
                    drawn_objects.append(overlay)
            elif slice_fit.get("fit_attempted", slice_fit.get("valid")):
                _print_proton_debug(
                    "overlay construction failed",
                    delta_index=delta_index,
                    aero_index=slice_index - 2,
                    fit_attempted=slice_fit.get("fit_attempted"),
                    valid=slice_fit.get("valid"),
                    kaon_mean=global_shape.get("kaon_mean"),
                    proton_mean=global_shape.get("proton_mean"),
                    kaon_sigma=global_shape.get("kaon_sigma"),
                    proton_sigma=global_shape.get("proton_sigma"),
                    kaon_amplitude=slice_fit.get("kaon_amplitude"),
                    proton_amplitude=slice_fit.get("proton_amplitude"),
                    other_amplitude=slice_fit.get("other_amplitude"),
                )
            drawn_objects.append(
                _draw_status_pave(
                    cleaning_result,
                    extra_lines=(
                        "d{} a{} attempted={} valid={}".format(
                            delta_index + 1,
                            slice_index - 1,
                            bool(slice_fit.get("fit_attempted", slice_fit.get("valid"))),
                            bool(slice_fit.get("valid")),
                        ),
                        "status={}".format(slice_fit.get("fit_status")),
	                        "entries={}".format(slice_fit.get("support_entries", "n/a")),
	                        "model/data={}".format(slice_fit.get("model_data_ratio")),
	                        "chi2/ndf={}".format(slice_fit.get("chi2_ndf")),
	                        "warn={}".format(
	                            ", ".join(str(w) for w in (slice_fit.get("diagnostic_warnings") or [])) or "none"
	                        ),
	                        "reason={}".format(slice_fit.get("rejection_reason") or "none"),
	                    ),
                    x1=0.16,
                    y1=0.66,
                    x2=0.62,
                    y2=0.90,
                )
            )
            gPad.Modified()
            gPad.Update()
        canvas.Modified()
        canvas.Update()
        gc.collect()
        canvas.Print(output_pdf)

    h_mm_before = application.get("H_MM_before_proton_cleaning")
    h_mm_proton = application.get("H_MM_estimated_proton")
    h_mm_after = application.get("H_MM_after_proton_cleaning")
    if h_mm_before is not None and h_mm_proton is not None and h_mm_after is not None:
        canvas = TCanvas("C_proton_cleaning_mm", "{} proton-cleaning MM validation".format(prefix), 1000, 900)
        canvas.Divide(1, 2)
        drawn_objects = [h_mm_before, h_mm_proton, h_mm_after]
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
        drawn_objects.append(legend)
        gPad.Modified()
        gPad.Update()
        canvas.cd(2)
        if application.get("H_proton_fraction_vs_MM") is not None:
            application["H_proton_fraction_vs_MM"].SetLineColor(kViolet + 1)
            application["H_proton_fraction_vs_MM"].SetMinimum(0.0)
            application["H_proton_fraction_vs_MM"].SetMaximum(1.0)
            application["H_proton_fraction_vs_MM"].Draw("hist")
            drawn_objects.append(application["H_proton_fraction_vs_MM"])
        gPad.Modified()
        gPad.Update()
        canvas.Modified()
        canvas.Update()
        gc.collect()
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
            drawn_objects = [h_mm_before, h_mm_proton, h_mm_after]
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
            drawn_objects.append(legend)
            gPad.Modified()
            gPad.Update()
            canvas.Modified()
            canvas.Update()
            gc.collect()
            canvas.Print(output_pdf)
            h_mm_before.GetXaxis().UnZoom()
            h_mm_proton.GetXaxis().UnZoom()
            h_mm_after.GetXaxis().UnZoom()

    _print_kaon_proton_cleaning_final_summary_page(output_pdf, cleaning_result, prefix)
