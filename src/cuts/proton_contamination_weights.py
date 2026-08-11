#! /usr/bin/python

from __future__ import annotations

import gc
import hashlib
import math
import os
import sys
from array import array
from copy import deepcopy

import numpy as np
import ROOT
from ROOT import (
    TCanvas,
    TLatex,
    TLegend,
    TLine,
    TPad,
    TPaveText,
    TH1D,
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
    PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED,
    PROTON_CONTAMINATION_CLEANING_METHOD_CTIME_AERO_EVENT_WEIGHT,
    PROTON_CONTAMINATION_CLEANING_METHOD_DISABLED,
    PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT,
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
APPLIED_ZERO_REASON_TO_CODE = {
    "applied": 0,
    "invalid_cell_fit": 1,
    "weak_proton_component": 2,
    "unsupported_delta": 3,
    "setting_gate_rejected": 4,
}
TIMING_CENTER_SOURCE_TO_CODE = {
    "delta_all_t_fit": 0,
    "delta_supported_t_fit": 1,
    "neighbor_delta_interpolation": 2,
    "neighbor_delta_nearest_fallback": 3,
    "stable_global_center_fallback": 4,
    "invalid_timing_center": 5,
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
PROTON_CLEANING_TOF_SUMMARY_VALIDATION = {
    "minimum_prompt_tof_events": 30,
    "minimum_valid_tof_fraction": 0.90,
}
PROTON_CLEANING_LOW_AERO_OFFSET_CONFIG = {
    "primary_range": (0.0, 5.0),
    "fallback_range": (0.0, 6.0),
    "full_diagnostic_range": (0.0, 25.0),
    "minimum_prompt_events": 20,
    "minimum_valid_tof_events": 10,
    "minimum_valid_tof_fraction": 0.50,
    "show_reference_npe": 5.0,
}
PROTON_CLEANING_TOF_REQUIRED_ALIASES = (
    "P_gtr_p",
    "ssxptar",
    "ssyptar",
    "ssdelta",
)
PROTON_CLEANING_TOF_FORBIDDEN_REPLAY_NAMES = (
    "P.gtr.p",
    "P.gtr.th",
    "P.gtr.ph",
    "P.gtr.dp",
    "P.gtr.xp",
    "P.gtr.yp",
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
    tof_summary_validation=None,
    reason_prefix="",
):
    validation_cfg = deepcopy(
        tof_summary_validation or PROTON_CLEANING_TOF_SUMMARY_VALIDATION
    )
    minimum_prompt_events = int(
        validation_cfg.get(
            "minimum_prompt_events",
            validation_cfg.get("minimum_prompt_tof_events", 30),
        )
        or 30
    )
    minimum_valid_tof_events = int(
        validation_cfg.get("minimum_valid_tof_events", minimum_prompt_events)
        or minimum_prompt_events
    )
    minimum_valid_tof_fraction = float(
        validation_cfg.get("minimum_valid_tof_fraction", 0.90) or 0.90
    )
    reason_prefix = str(reason_prefix or "").strip()
    def _summary_reason(kind):
        if reason_prefix == "low_aero":
            mapping = {
                "insufficient_prompt_events": "insufficient_low_aero_prompt_events",
                "insufficient_valid_tof_events": "insufficient_low_aero_valid_tof_events",
                "valid_tof_fraction_below_min": "low_aero_valid_tof_fraction_below_min",
                "invalid_mean_delta_t_pk_ns": "invalid_low_aero_mean_delta_t_pk",
                "invalid_mean_P_gtr_p": "invalid_low_aero_mean_P_gtr_p",
                "invalid_mean_shms_path_length_cm": "invalid_low_aero_mean_shms_path_length_cm",
            }
            return mapping.get(kind, "low_aero_{}".format(kind))
        return kind
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
            "prompt_events_with_selected_timing": 0,
            "prompt_events_inside_timing_range": 0,
            "prompt_events_inside_aero_range": 0,
            "prompt_events_inside_timing_and_aero_domain": 0,
            "prompt_events_with_valid_tof": 0,
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
        counters[delta_index]["prompt_events_inside_timing_and_aero_domain"] += 1
        if not bool((entry_payload or {}).get("tof_valid", False)):
            continue
        counters[delta_index]["prompt_events_with_valid_tof"] += 1
        candidate_values = {}
        for key in collections[delta_index]:
            value = (entry_payload or {}).get(key)
            try:
                value = float(value)
            except Exception:
                candidate_values = None
                break
            if not math.isfinite(float(value)):
                candidate_values = None
                break
            candidate_values[key] = float(value)
        if candidate_values is None:
            continue
        counters[delta_index]["prompt_events_used"] += 1
        for key, value in candidate_values.items():
            collections[delta_index][key].append(float(value))
    summaries = []
    for delta_index, collection in enumerate(collections):
        denominator = int(
            counters[delta_index]["prompt_events_inside_timing_and_aero_domain"] or 0
        )
        prompt_events_used = int(counters[delta_index]["prompt_events_used"] or 0)
        valid_tof_fraction = (
            float(prompt_events_used / float(denominator))
            if denominator > 0
            else 0.0
        )
        row = {
            "delta_index": int(delta_index),
            "delta_min": float(delta_edges[delta_index]),
            "delta_max": float(delta_edges[delta_index + 1]),
            "prompt_event_count": int(len(collection["delta_t_pk_ns"])),
            "usable_tof_events": int(prompt_events_used),
            "valid_tof_fraction": float(valid_tof_fraction),
            "tof_summary_validation": _json_ready_value(validation_cfg),
            **counters[delta_index],
        }
        for key, values in collection.items():
            mean_value, rms_value = _mean_rms(values)
            row["mean_" + key] = mean_value
            row["rms_" + key] = rms_value
        rejection_reasons = []
        if int(denominator) < int(minimum_prompt_events):
            rejection_reasons.append(_summary_reason("insufficient_prompt_events"))
        if int(prompt_events_used) < int(minimum_valid_tof_events):
            rejection_reasons.append(_summary_reason("insufficient_valid_tof_events"))
        if float(valid_tof_fraction) < float(minimum_valid_tof_fraction):
            rejection_reasons.append(_summary_reason("valid_tof_fraction_below_min"))
        for field_name in (
            "mean_delta_t_pk_ns",
            "mean_P_gtr_p",
            "mean_shms_path_length_cm",
        ):
            field_value = row.get(field_name)
            if (
                field_value is None
                or not math.isfinite(float(field_value))
                or float(field_value) <= 0.0
            ):
                rejection_reasons.append(_summary_reason("invalid_{}".format(field_name)))
        row["valid"] = bool(len(rejection_reasons) == 0)
        row["rejection_reasons"] = list(rejection_reasons)
        row["rejection_reason"] = _join_rejection_reasons(rejection_reasons)
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


def _count_rejection_reasons(rows):
    counts = {}
    for row in rows or []:
        if not isinstance(row, dict):
            continue
        reasons = row.get("rejection_reasons")
        if not reasons and row.get("rejection_reason"):
            reasons = [row.get("rejection_reason")]
        for reason in reasons or []:
            reason_text = str(reason).strip()
            if not reason_text:
                continue
            counts[reason_text] = int(counts.get(reason_text, 0)) + 1
    return counts


def _offset_summary_validation_from_low_aero_config(low_aero_config):
    low_aero_config = low_aero_config or {}
    return {
        "minimum_prompt_events": int(
            low_aero_config.get("minimum_prompt_events", 20) or 20
        ),
        "minimum_valid_tof_events": int(
            low_aero_config.get("minimum_valid_tof_events", 10) or 10
        ),
        "minimum_valid_tof_fraction": float(
            low_aero_config.get("minimum_valid_tof_fraction", 0.50) or 0.50
        ),
    }


def _project_delta_pid_timing_by_aero_range(pid_hist, name, aero_range, upper_inclusive=False):
    if pid_hist is None or not aero_range or len(aero_range) < 2:
        return None
    aero_min = float(aero_range[0])
    aero_max = float(aero_range[1])
    x_axis = pid_hist.GetXaxis()
    first_bin = max(1, int(x_axis.FindFixBin(aero_min)))
    upper_value = aero_max if upper_inclusive else np.nextafter(aero_max, aero_min)
    last_bin = min(int(pid_hist.GetNbinsX()), int(x_axis.FindFixBin(float(upper_value))))
    if last_bin < first_bin:
        return None
    projection = pid_hist.ProjectionY(str(name), first_bin, last_bin)
    if hasattr(projection, "SetDirectory"):
        projection.SetDirectory(0)
    if hasattr(projection, "Sumw2"):
        projection.Sumw2()
    return projection


def _clone_tof_summary_with_mode(row, mode, aero_range, low_aero_config):
    cloned = deepcopy(row or {})
    if aero_range and len(aero_range) >= 2:
        cloned["offset_fit_aero_min"] = float(aero_range[0])
        cloned["offset_fit_aero_max"] = float(aero_range[1])
    else:
        cloned["offset_fit_aero_min"] = None
        cloned["offset_fit_aero_max"] = None
    cloned["offset_fit_aero_mode"] = str(mode)
    cloned["low_aero_offset_config"] = _json_ready_value(low_aero_config or {})
    return cloned


def _select_low_aero_offset_summary(primary_row, fallback_row, low_aero_config):
    primary_range = tuple((low_aero_config or {}).get("primary_range") or (0.0, 5.0))
    fallback_range = tuple((low_aero_config or {}).get("fallback_range") or (0.0, 6.0))
    if bool((primary_row or {}).get("valid", False)):
        return _clone_tof_summary_with_mode(
            primary_row,
            "low_aero_0_5",
            primary_range,
            low_aero_config,
        )
    if bool((fallback_row or {}).get("valid", False)):
        selected = _clone_tof_summary_with_mode(
            fallback_row,
            "low_aero_0_6_fallback",
            fallback_range,
            low_aero_config,
        )
        selected["primary_low_aero_rejection_reasons"] = list(
            (primary_row or {}).get("rejection_reasons") or []
        )
        return selected
    selected = _clone_tof_summary_with_mode(
        primary_row or fallback_row or {},
        "unavailable",
        primary_range,
        low_aero_config,
    )
    rejection_reasons = []
    for row in (primary_row, fallback_row):
        for reason in (row or {}).get("rejection_reasons") or []:
            reason_text = str(reason).strip()
            if reason_text and reason_text not in rejection_reasons:
                rejection_reasons.append(reason_text)
    selected["valid"] = False
    selected["rejection_reasons"] = rejection_reasons or ["low_aero_offset_unavailable"]
    selected["rejection_reason"] = _join_rejection_reasons(selected["rejection_reasons"])
    selected["primary_low_aero_rejection_reasons"] = list(
        (primary_row or {}).get("rejection_reasons") or []
    )
    selected["fallback_low_aero_rejection_reasons"] = list(
        (fallback_row or {}).get("rejection_reasons") or []
    )
    return selected


def _make_offset_attempt_summary(offset_fit):
    offset_fit = offset_fit or {}
    return {
        "attempted": bool(offset_fit.get("fit_attempted", False)),
        "valid": bool(offset_fit.get("valid", False)),
        "rejection_reasons": list(offset_fit.get("rejection_reasons") or []),
        "rejection_reason": offset_fit.get("rejection_reason", ""),
        "delta_offset": offset_fit.get("delta_offset"),
        "delta_offset_error": offset_fit.get("delta_offset_error"),
        "fit_status": offset_fit.get("fit_status"),
        "fit_status_code": offset_fit.get("fit_status_code"),
        "chi2_ndf": offset_fit.get("chi2_ndf"),
    }


def _make_not_required_offset_fit(delta_index, mode, reason):
    return {
        "valid": False,
        "fit_attempted": False,
        "fit_status": "not_required",
        "fit_status_code": None,
        "delta_index": int(delta_index),
        "offset_fit_aero_mode": str(mode),
        "delta_offset": 0.0,
        "delta_offset_error": None,
        "rejection_reasons": [str(reason)],
        "rejection_reason": str(reason),
    }


CELL_FIT_SKIP_COUNTER_FIELDS = {
    "insufficient_support": "cell_fit_skipped_insufficient_support_count",
    "missing_histogram": "cell_fit_skipped_missing_histogram_count",
    "invalid_global_shape": "cell_fit_skipped_invalid_global_shape_count",
    "invalid_timing_center_model": "cell_fit_skipped_invalid_timing_center_model_count",
}
CELL_FIT_SKIPPED_OTHER_COUNTER_FIELD = "cell_fit_skipped_other_count"


def _cell_fit_skip_counter_key(fit_status):
    """Map an explicit pre-fit status to its diagnostic skip counter."""
    return CELL_FIT_SKIP_COUNTER_FIELDS.get(
        str(fit_status or "").strip(),
        CELL_FIT_SKIPPED_OTHER_COUNTER_FIELD,
    )


def _build_selected_tof_summary(delta_offset_fit, primary_summary, fallback_summary):
    """Attach the TOF diagnostic to the timing-center source actually selected."""
    source = str(
        (delta_offset_fit or {}).get("selected_timing_center_source")
        or (delta_offset_fit or {}).get("timing_center_source")
        or "unavailable"
    )
    if source == "low_aero_0_5_fit":
        selected = deepcopy(primary_summary or {})
        role = "offset_fit_input"
    elif source == "low_aero_0_6_fit":
        selected = deepcopy(fallback_summary or {})
        role = "offset_fit_input"
    elif source == "stable_global_center_fallback":
        selected = deepcopy(primary_summary or fallback_summary or {})
        role = "stable_fallback_diagnostic"
    else:
        selected = {}
        role = "unavailable"
    selected["selected_timing_center_source"] = source
    selected["tof_summary_role"] = role
    return selected


def _decorate_selected_timing_center_model(
    selected_fit,
    source,
    primary_fit,
    fallback_fit,
    primary_summary,
    fallback_summary,
    reference_shape,
):
    selected_fit = deepcopy(selected_fit or {})
    source = str(source)
    selected_fit["selected_timing_center_source"] = source
    selected_fit["timing_center_source"] = source
    selected_fit["selected_offset_source"] = source
    selected_fit["offset_refinement_applied"] = bool(source in ("low_aero_0_5_fit", "low_aero_0_6_fit"))
    selected_fit["offset_refinement_valid"] = bool(selected_fit.get("valid", False))
    selected_fit["timing_center_model_valid"] = bool(
        selected_fit.get("valid", False)
        or (
            source == "stable_global_center_fallback"
            and bool((reference_shape or {}).get("valid", False))
        )
    )
    selected_fit["primary_offset_attempt"] = _json_ready_value(
        _make_offset_attempt_summary(primary_fit)
    )
    selected_fit["fallback_offset_attempt"] = _json_ready_value(
        _make_offset_attempt_summary(fallback_fit)
    )
    selected_fit["primary_offset_attempted"] = bool((primary_fit or {}).get("fit_attempted", False))
    selected_fit["primary_offset_valid"] = bool((primary_fit or {}).get("valid", False))
    selected_fit["primary_offset_rejection_reasons"] = list((primary_fit or {}).get("rejection_reasons") or [])
    selected_fit["fallback_offset_attempted"] = bool((fallback_fit or {}).get("fit_attempted", False))
    selected_fit["fallback_offset_valid"] = bool((fallback_fit or {}).get("valid", False))
    selected_fit["fallback_offset_rejection_reasons"] = list((fallback_fit or {}).get("rejection_reasons") or [])
    selected_fit["primary_low_aero_tof_summary"] = _json_ready_value(primary_summary or {})
    selected_fit["fallback_low_aero_tof_summary"] = _json_ready_value(fallback_summary or {})
    selected_fit["offset_refinement_failure_reasons"] = []
    for attempt in (primary_fit, fallback_fit):
        if bool((attempt or {}).get("valid", False)):
            continue
        for reason in (attempt or {}).get("rejection_reasons") or []:
            reason_text = str(reason).strip()
            if reason_text and reason_text not in selected_fit["offset_refinement_failure_reasons"]:
                selected_fit["offset_refinement_failure_reasons"].append(reason_text)
    if source == "stable_global_center_fallback":
        selected_fit["valid"] = False
        selected_fit["offset_refinement_valid"] = False
        selected_fit["offset_refinement_applied"] = False
        selected_fit["delta_offset"] = 0.0
        selected_fit["delta_offset_error"] = None
        selected_fit["fit_status"] = "stable_global_center_fallback"
        selected_fit["fit_status_code"] = None
        selected_fit["rejection_reasons"] = list(
            selected_fit.get("offset_refinement_failure_reasons")
            or ["low_aero_offset_refinement_failed"]
        )
        selected_fit["rejection_reason"] = _join_rejection_reasons(
            selected_fit["rejection_reasons"]
        )
        selected_fit["timing_center_rejection_reasons"] = (
            []
            if selected_fit["timing_center_model_valid"]
            else [(reference_shape or {}).get("reason") or "invalid_reference_shape"]
        )
        selected_fit["reference_kaon_mean"] = (reference_shape or {}).get("kaon_mean")
        selected_fit["reference_proton_mean"] = (reference_shape or {}).get("proton_mean")
        selected_fit["reference_kaon_sigma"] = (reference_shape or {}).get("kaon_sigma")
        selected_fit["reference_proton_sigma"] = (reference_shape or {}).get("proton_sigma")
        selected_fit["mean_delta_t_pk_ns"] = (
            (primary_summary or {}).get("mean_delta_t_pk_ns")
            or (fallback_summary or {}).get("mean_delta_t_pk_ns")
        )
    return selected_fit


def _select_resolved_timing_center_model(
    delta_index,
    primary_fit,
    fallback_fit,
    primary_summary,
    fallback_summary,
    reference_shape,
):
    if bool((primary_fit or {}).get("valid", False)):
        return _decorate_selected_timing_center_model(
            primary_fit,
            "low_aero_0_5_fit",
            primary_fit,
            fallback_fit,
            primary_summary,
            fallback_summary,
            reference_shape,
        )
    if bool((fallback_fit or {}).get("valid", False)):
        return _decorate_selected_timing_center_model(
            fallback_fit,
            "low_aero_0_6_fit",
            primary_fit,
            fallback_fit,
            primary_summary,
            fallback_summary,
            reference_shape,
        )
    stable_model = {
        "delta_index": int(delta_index),
        "offset_fit_aero_mode": "stable_global_center_fallback",
    }
    return _decorate_selected_timing_center_model(
        stable_model,
        "stable_global_center_fallback",
        primary_fit,
        fallback_fit,
        primary_summary,
        fallback_summary,
        reference_shape,
    )


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
    exact_config["tof_summary_validation"] = deepcopy(
        (base_config or {}).get("tof_summary_validation")
        or PROTON_CLEANING_TOF_SUMMARY_VALIDATION
    )
    exact_config["low_aero_offset"] = deepcopy(
        (base_config or {}).get("low_aero_offset")
        or PROTON_CLEANING_LOW_AERO_OFFSET_CONFIG
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
    allowed_sources = None
    if source_names is not None:
        allowed_sources = {str(source_name) for source_name in source_names}
    if allowed_sources is None:
        prepared_branches = (source_bundle or {}).get("available_timing_branches") or ()
        if branch_name in prepared_branches:
            return True
    for source_name, source_spec in ((source_bundle or {}).get("prepared_sources") or {}).items():
        if allowed_sources is not None and str(source_name) not in allowed_sources:
            continue
        prepared_branches = (source_spec or {}).get("available_timing_branches") or ()
        if branch_name in prepared_branches:
            return True
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
    return [
        candidate
        for candidate in candidates
        if _bundle_has_branch(source_bundle, candidate, source_names=("prompt",))
    ]


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
    from apply_cuts import evaluate_pre_particle_subtraction_event

    prepared_bundle = dict(source_bundle or {})
    prepared_sources = {}
    prepared_source_stats = {}
    available_timing_branches = set()
    cross_stage_rows = []
    prepass_samples = dict((source_bundle or {}).get("canonical_t_prepass_samples") or {})
    cross_stage_sampling = dict((source_bundle or {}).get("canonical_t_prepass_sampling") or {})
    cross_stage_tolerance = float(
        ((proton_cleaning_config or {}).get("t_binning") or {}).get(
            "cross_stage_t_tolerance", 1.0e-10
        )
    )
    strict_cross_stage = bool((proton_cleaning_config or {}).get("strict_mode", False))
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
                # ``apply_cuts`` is the single owner of shifted t/mm and the
                # non-particle-subtraction selection shared with the prepass.
                selection_state = evaluate_pre_particle_subtraction_event(
                    evt, mm_min, mm_max, hole_contains=hole_contains
                )
                allcuts = bool(selection_state["allcuts"])
                nommcuts = bool(selection_state["nommcuts"])
                noholecuts = bool(selection_state["noholecuts"])
                hole_rejected = bool(selection_state["hole_rejected"])
                adj_hsdelta = float(selection_state["adj_hsdelta"])
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
                    "adj_mm": float(selection_state["adj_mm"]),
                    "adj_t": float(selection_state["adj_t"]),
                    "adj_hsdelta": float(adj_hsdelta),
                    "delta_value": float(getattr(evt, "ssdelta", 0.0)),
                    "aero_value": float(getattr(evt, "P_aero_npeSum", 0.0)),
                    "phi_value": float(getattr(evt, "ph_q", float("nan"))),
                    "timing_values": timing_values,
                    **tof_payload,
                }
                signature = _make_prepared_event_signature(source_name, entry_index)
                if signature in prepass_samples:
                    prepass_t = float(prepass_samples[signature])
                    prepared_t = float(selection_state["adj_t"])
                    difference = abs(prepass_t - prepared_t)
                    row = {
                        "event_signature": signature,
                        "prepass_t": prepass_t,
                        "prepared_proton_cleaning_adj_t": prepared_t,
                        "downstream_t": None,
                        "maximum_absolute_difference": difference,
                        "consistent": bool(difference <= cross_stage_tolerance),
                    }
                    cross_stage_rows.append(row)
                    if strict_cross_stage and not row["consistent"]:
                        raise RuntimeError(
                            "prepass/proton shifted-t mismatch for {}".format(signature)
                        )
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
    prepared_bundle["cross_stage_t_consistency"] = cross_stage_rows
    prepared_bundle["cross_stage_t_sampling"] = cross_stage_sampling
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
            "timing_model_valid": False,
            "cell_fit_valid": False,
            "proton_component_detected": False,
            "proton_component_significance": None,
            "proton_component_below_significance": True,
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
        "smaller_component_fraction_definition": "min(K_amp,p_amp)/(K_amp+p_amp)",
        "tof_offset_validation": _json_ready_value(validation_cfg),
    }
    rejection_reasons = []
    if histogram is None:
        rejection_reasons.append("missing_delta_timing_histogram")
    if not bool((reference_shape or {}).get("valid", False)):
        rejection_reasons.append((reference_shape or {}).get("reason") or "invalid_reference_shape")
    if not bool((tof_summary or {}).get("valid", False)):
        rejection_reasons.append("invalid_prompt_tof_summary")
        for reason in (tof_summary or {}).get("rejection_reasons") or []:
            reason_text = str(reason).strip()
            if reason_text and reason_text not in rejection_reasons:
                rejection_reasons.append(reason_text)
    mean_delta_t = (tof_summary or {}).get("mean_delta_t_pk_ns")
    if mean_delta_t is None or not math.isfinite(float(mean_delta_t)) or float(mean_delta_t) <= 0.0:
        rejection_reasons.append("invalid_mean_delta_t_pk")
    low_aero_config = exact_config.get("low_aero_offset") or PROTON_CLEANING_LOW_AERO_OFFSET_CONFIG
    minimum_entries = int(low_aero_config.get("minimum_prompt_events", 20) or 20)
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
        4,
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
    component_denominator = kaon_amplitude + proton_amplitude
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
            "smaller_component_fraction_definition": "min(K_amp,p_amp)/(K_amp+p_amp)",
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
    allow_stable_global_center_fallback=True,
):
    if not bool((global_shape or {}).get("valid", False)):
        return {"valid": False, "reason": "invalid_global_shape"}
    center_source = str(
        (delta_offset_fit or {}).get("timing_center_source")
        or (delta_offset_fit or {}).get("selected_timing_center_source")
        or ("low_aero_offset_fit" if bool((delta_offset_fit or {}).get("valid", False)) else "")
    )
    timing_center_model_valid = bool(
        (delta_offset_fit or {}).get(
            "timing_center_model_valid",
            bool((delta_offset_fit or {}).get("valid", False)),
        )
    )
    if not timing_center_model_valid:
        return {
            "valid": False,
            "reason": "invalid_timing_center_model",
            "timing_center_source": center_source or "unavailable",
            "offset_refinement_valid": bool((delta_offset_fit or {}).get("valid", False)),
            "offset_refinement_applied": bool((delta_offset_fit or {}).get("offset_refinement_applied", False)),
            "offset_refinement_failure_reasons": list((delta_offset_fit or {}).get("offset_refinement_failure_reasons") or []),
        }
    delta_offset = float((delta_offset_fit or {}).get("delta_offset", 0.0) or 0.0)
    reference_kaon_mean = float((global_shape or {}).get("kaon_mean"))
    reference_proton_mean = float((global_shape or {}).get("proton_mean"))
    stable_center_fallback = center_source == "stable_global_center_fallback"
    if stable_center_fallback and not bool(allow_stable_global_center_fallback):
        return {
            "valid": False,
            "reason": "local_low_aero_timing_offset_required",
            "timing_center_source": center_source,
            "offset_refinement_valid": bool((delta_offset_fit or {}).get("valid", False)),
            "offset_refinement_applied": bool((delta_offset_fit or {}).get("offset_refinement_applied", False)),
            "offset_refinement_failure_reasons": list((delta_offset_fit or {}).get("offset_refinement_failure_reasons") or []),
        }
    mean_delta_t = (tof_summary or {}).get("mean_delta_t_pk_ns")
    if stable_center_fallback:
        predicted_kaon_mean = reference_kaon_mean
        raw_predicted_proton_mean = reference_proton_mean
        mean_delta_t_out = (
            float(mean_delta_t)
            if mean_delta_t is not None and math.isfinite(float(mean_delta_t))
            else abs(reference_proton_mean - reference_kaon_mean)
        )
    else:
        if mean_delta_t is None or not math.isfinite(float(mean_delta_t)) or float(mean_delta_t) <= 0.0:
            return {
                "valid": False,
                "reason": "invalid_mean_delta_t_pk",
                "timing_center_source": center_source or "unavailable",
                "offset_refinement_valid": bool((delta_offset_fit or {}).get("valid", False)),
                "offset_refinement_applied": bool((delta_offset_fit or {}).get("offset_refinement_applied", False)),
                "offset_refinement_failure_reasons": list((delta_offset_fit or {}).get("offset_refinement_failure_reasons") or []),
            }
        branch_sign = -1.0 if bool(proton_peak_is_lower) else 1.0
        predicted_kaon_mean = reference_kaon_mean + delta_offset
        raw_predicted_proton_mean = predicted_kaon_mean + (branch_sign * float(mean_delta_t))
        mean_delta_t_out = float(mean_delta_t)
    wrap_info = {
        "valid": True,
        "raw_mean": float(raw_predicted_proton_mean),
        "wrapped_mean": float(raw_predicted_proton_mean),
        "period_shift": 0,
        "inside_display_range": True,
        "reference_mean": float(reference_proton_mean),
    }
    if str(probe_kind) == "rf" and not stable_center_fallback:
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
        "timing_center_model_valid": True,
        "timing_center_source": center_source or "low_aero_offset_fit",
        "selected_timing_center_source": center_source or "low_aero_offset_fit",
        "offset_refinement_valid": bool((delta_offset_fit or {}).get("offset_refinement_valid", (delta_offset_fit or {}).get("valid", False))),
        "offset_refinement_applied": bool((delta_offset_fit or {}).get("offset_refinement_applied", False)),
        "offset_refinement_failure_reasons": list((delta_offset_fit or {}).get("offset_refinement_failure_reasons") or []),
        "reference_global_kaon_mean": float(reference_kaon_mean),
        "reference_global_proton_mean": float(reference_proton_mean),
        "delta_timing_offset": float(delta_offset),
        "delta_timing_offset_error": (delta_offset_fit or {}).get("delta_offset_error"),
        "mean_delta_t_pk_ns": float(mean_delta_t_out),
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
        if histogram is None:
            fit_status = "missing_histogram"
        elif invalid_global_shape:
            fit_status = "invalid_global_shape"
        elif invalid_timing_constraint:
            fit_status = "invalid_timing_center_model"
        else:
            fit_status = "insufficient_support"
        if invalid_global_shape:
            rejection_reasons.append("invalid_global_shape")
        if invalid_timing_constraint:
            rejection_reasons.append("invalid_timing_center_model")
            timing_reason = str((timing_constraint or {}).get("reason") or "").strip()
            if timing_reason and timing_reason not in rejection_reasons:
                rejection_reasons.append(timing_reason)
        if histogram is None:
            rejection_reasons.append("missing_histogram")
        if int(support_entries) < int(minimum_entries):
            rejection_reasons.append("insufficient_entries")
        return {
            "valid": False,
            "timing_model_valid": bool(not invalid_global_shape and not invalid_timing_constraint),
            "timing_center_model_valid": bool((timing_constraint or {}).get("timing_center_model_valid", False)),
            "timing_center_source": (timing_constraint or {}).get("timing_center_source"),
            "offset_refinement_valid": (timing_constraint or {}).get("offset_refinement_valid"),
            "offset_refinement_applied": (timing_constraint or {}).get("offset_refinement_applied"),
            "offset_refinement_failure_reasons": (timing_constraint or {}).get("offset_refinement_failure_reasons"),
            "cell_fit_valid": False,
            "proton_component_detected": False,
            "proton_component_significance": None,
            "proton_component_below_significance": True,
            "fit_attempted": False,
            "fit_status": fit_status,
            "fit_status_code": None,
            "function_name": str(function_name),
            "excluded_invalid_variance_bins": 0,
            "invalid_bin_rule": "macro ROOT fit uses all histogram bins in the fit range",
            "fit_options": PROTON_CLEANING_EXACT_FIT_OPTIONS,
            "rejection_reasons": rejection_reasons,
            "rejection_reason": _join_rejection_reasons(rejection_reasons),
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
    proton_component_significance = (
        float(proton_amplitude / proton_amplitude_error)
        if proton_amplitude_error > 0.0 and math.isfinite(float(proton_amplitude_error))
        else None
    )
    proton_significance_threshold = float(
        (config.get("tof_offset_validation") or {}).get(
            "minimum_component_significance",
            (config.get("global_fit") or {}).get("minimum_amplitude_significance", 2.0),
        )
        or 2.0
    )
    proton_component_detected = bool(
        proton_component_significance is not None
        and math.isfinite(float(proton_component_significance))
        and float(proton_component_significance) >= float(proton_significance_threshold)
        and float(proton_amplitude) > 0.0
    )
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
    raw_proton_yield = float(proton_yield)
    raw_proton_amplitude = float(proton_amplitude)
    if not proton_component_detected:
        proton_yield = 0.0
        proton_amplitude = 0.0
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
    model_yield = float(kaon_yield + raw_proton_yield + other_yield)
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
        "timing_model_valid": True,
        "cell_fit_valid": bool(valid),
        "proton_component_detected": bool(proton_component_detected and valid),
        "proton_component_significance": (
            float(proton_component_significance)
            if proton_component_significance is not None
            else None
        ),
        "proton_component_below_significance": bool(
            not proton_component_detected
        ),
        "proton_component_significance_threshold": float(proton_significance_threshold),
        "fit_attempted": True,
        "fit_status": "success" if fit_status_code == 0 else "failure",
        "fit_status_code": int(fit_status_code),
        "message": "",
        "function_name": str(function_name),
        "reference_global_kaon_mean": (timing_constraint or {}).get("reference_global_kaon_mean", global_shape.get("kaon_mean")),
        "reference_global_proton_mean": (timing_constraint or {}).get("reference_global_proton_mean", global_shape.get("proton_mean")),
        "timing_center_model_valid": (timing_constraint or {}).get("timing_center_model_valid"),
        "timing_center_source": (timing_constraint or {}).get("timing_center_source"),
        "offset_refinement_valid": (timing_constraint or {}).get("offset_refinement_valid"),
        "offset_refinement_applied": (timing_constraint or {}).get("offset_refinement_applied"),
        "offset_refinement_failure_reasons": (timing_constraint or {}).get("offset_refinement_failure_reasons"),
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
        "raw_proton_amplitude": float(raw_proton_amplitude),
        "proton_amplitude": float(proton_amplitude),
        "proton_amplitude_error": float(proton_amplitude_error),
        "other_amplitude": float(other_amplitude),
        "other_amplitude_error": float(other_amplitude_error),
        "kaon_yield": float(kaon_yield),
        "raw_proton_yield": float(raw_proton_yield),
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


def _finite_strictly_increasing_edges(values):
    try:
        edges = [float(value) for value in (values or [])]
    except (TypeError, ValueError):
        return []
    if len(edges) < 2 or not all(math.isfinite(value) for value in edges):
        return []
    if any(right <= left for left, right in zip(edges[:-1], edges[1:])):
        return []
    return edges


def _resolve_aerogel_summary_edges(validation_cfg):
    """Use the new summary key while accepting pre-existing setting overrides."""
    cfg = dict(validation_cfg or {})
    return _finite_strictly_increasing_edges(
        cfg.get("summary_slice_edges") or cfg.get("slice_edges") or ()
    )


def _resolve_aerogel_display_edges(validation_cfg):
    cfg = dict(validation_cfg or {})
    display_range = cfg.get("display_range") or cfg.get("hist_range") or (0.0, 25.0)
    try:
        low, high = float(display_range[0]), float(display_range[1])
        bins = max(1, int(cfg.get("display_bins", 100) or 100))
    except (TypeError, ValueError, IndexError):
        return []
    if not (math.isfinite(low) and math.isfinite(high) and high > low):
        return []
    return np.linspace(low, high, bins + 1, dtype=float).tolist()


def _safe_validation_ratio(numerator, denominator):
    try:
        numerator = float(numerator)
        denominator = float(denominator)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(numerator) and math.isfinite(denominator)) or denominator == 0.0:
        return None
    return float(numerator / denominator)


# The compact canonical-t x aerogel products have one serializable source of
# truth.  ROOT objects and PDF pages are presentation layers over this payload;
# they must never independently re-accumulate a differently-defined matrix.
_T_AEROGEL_MATRIX_METRIC_SOURCES = {
    "selected_prompt_count": "raw_prompt_event_count",
    "signed_physical_yield": "signed_event_weight_sum",
    "absolute_physical_support": "absolute_event_weight_support",
    "estimated_proton_yield": "estimated_proton_missing_mass_yield",
    "cleaned_yield": "cleaned_missing_mass_yield",
    "average_proton_probability": "average_proton_probability",
    "low_mm_removed_yield": "low_mm_removed_yield",
    "low_mm_removed_fraction": "low_mm_removed_fraction",
    "lambda_removed_yield": "lambda_removed_yield",
    "lambda_removed_fraction": "lambda_removed_fraction",
}
LAMBDA_PRESERVATION_AUDIT_LIMIT = 128


def _matrix_from_t_aerogel_cells(cells, source_key):
    return [
        [
            (cell or {}).get(source_key)
            for cell in (row or [])
        ]
        for row in (cells or [])
    ]


def _sum_matrix_row(matrix, row_index):
    values = matrix[row_index] if 0 <= row_index < len(matrix or []) else []
    return float(sum(
        float(value) for value in values
        if _finite_float_or_none(value) is not None
    ))


def _build_t_aerogel_matrix_payload(cells, t_edges, aero_edges, exclusions=None):
    """Return the canonical coarse-matrix contract for all diagnostics.

    ``cells`` is already accumulated from frozen lookup rows.  Keeping the
    conversion to named matrices here prevents JSON, CSV, ROOT, and PDF code
    from accidentally selecting different source keys for the same label.
    """
    metrics = {
        name: _matrix_from_t_aerogel_cells(cells, source)
        for name, source in _T_AEROGEL_MATRIX_METRIC_SOURCES.items()
    }
    validity = {
        name: _matrix_from_t_aerogel_cells(cells, "{}_valid".format(source))
        for name, source in _T_AEROGEL_MATRIX_METRIC_SOURCES.items()
        if name in (
            "average_proton_probability",
            "low_mm_removed_fraction",
            "lambda_removed_fraction",
        )
    }
    return {
        "schema_version": 1,
        "source": "frozen_timing_t_lookup_rows",
        "t_edges": [float(edge) for edge in (t_edges or [])],
        "aero_edges": [float(edge) for edge in (aero_edges or [])],
        "metric_sources": dict(_T_AEROGEL_MATRIX_METRIC_SOURCES),
        "metrics": metrics,
        "validity_masks": validity,
        "exclusions": dict(exclusions or {}),
    }


def _matrix_metric(matrix_payload, metric_name):
    return ((matrix_payload or {}).get("metrics") or {}).get(metric_name) or []


def _matrix_has_nonzero_content(matrix, tolerance=0.0):
    tolerance = abs(float(tolerance or 0.0))
    return any(
        abs(float(value)) > tolerance
        for row in (matrix or [])
        for value in (row or [])
        if _finite_float_or_none(value) is not None
    )


def _build_timing_t_summary(cleaning_result, matrix_payload, aggregate_all):
    """Build a selected-candidate-only final timing-t diagnostic summary."""
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    selected = diagnostics.get("selected_timing_candidate") or {}
    delta_rows = list(diagnostics.get("delta_support") or [])
    applied_cells = list(diagnostics.get("applied_timing_t_cell_map") or [])
    closure_by_delta = list(diagnostics.get("event_weight_closure_by_delta") or [])
    delta_edges = [float(edge) for edge in ((cleaning_result or {}).get("delta_edges") or [])]
    closure_by_index = {
        int(row.get("delta_index")): row
        for row in closure_by_delta
        if isinstance(row, dict) and isinstance(row.get("delta_index"), int)
    }
    cells_by_delta = {}
    for cell in applied_cells:
        if not isinstance(cell, dict) or not isinstance(cell.get("delta_index"), int):
            continue
        cells_by_delta.setdefault(int(cell["delta_index"]), []).append(cell)
    per_delta = []
    for row in delta_rows:
        if not isinstance(row, dict) or not isinstance(row.get("delta_index"), int):
            continue
        index = int(row["delta_index"])
        closure = closure_by_index.get(index, {})
        cells = cells_by_delta.get(index, [])
        applied_yield = float(sum(
            float((cell or {}).get("applied_proton_yield", 0.0) or 0.0)
            for cell in cells
        ))
        data_total = float(row.get("data_total", 0.0) or 0.0)
        per_delta.append({
            "delta_index": index,
            "delta_low": delta_edges[index] if index < len(delta_edges) - 1 else None,
            "delta_high": delta_edges[index + 1] if index + 1 < len(delta_edges) else None,
            "support_label": row.get("support_label", SUPPORT_UNSUPPORTED),
            "raw_data_total": data_total,
            "fitted_data_total": float(row.get("fitted_data_total", 0.0) or 0.0),
            "model_yield": float(row.get("model_total", 0.0) or 0.0),
            "proton_yield": float(row.get("proton_total", 0.0) or 0.0),
            "kaon_yield": float(row.get("kaon_total", 0.0) or 0.0),
            "other_yield": float(row.get("other_total", 0.0) or 0.0),
            "valid_t_cell_count": int(row.get("valid_t_cells", 0) or 0),
            "coverage": float(row.get("coverage", 0.0) or 0.0),
            "applied_proton_yield": applied_yield,
            "applied_proton_fraction": (
                float(applied_yield / data_total) if data_total > 0.0 else None
            ),
            "mean_event_lookup_probability": _safe_validation_ratio(
                closure.get("summed_event_proton_probability"),
                closure.get("event_count"),
            ),
            "lookup_event_count": int(closure.get("event_count", 0) or 0),
        })
    mm_diagnostics = diagnostics.get("timing_t_mm_diagnostics") or {}
    mm_aggregate = mm_diagnostics.get("aggregate") or {}
    applied_mm_aggregate = ((mm_diagnostics.get("applied") or {}).get("aggregate") or {})
    lambda_gate = diagnostics.get("lambda_preservation_gate") or {}
    return {
        "schema_version": 1,
        "source": "selected_timing_candidate_and_frozen_lookup_only",
        "selected_candidate": {
            "timing_branch": selected.get("timing_branch") or diagnostics.get("timing_branch"),
            "selected": bool(selected.get("selected", True)),
            "candidate_selection_rank": selected.get("candidate_selection_rank"),
            "candidate_selection_tuple": diagnostics.get("candidate_selection_tuple") or [],
        },
        "delta_edges": delta_edges,
        "per_delta": per_delta,
        "setting_support": diagnostics.get("setting_support") or {},
        "matrix_source": (matrix_payload or {}).get("source"),
        "missing_mass_totals": {
            "raw": mm_aggregate.get("raw_missing_mass_yield", (aggregate_all or {}).get("raw_missing_mass_yield")),
            "estimated_proton": mm_aggregate.get(
                "estimated_proton_missing_mass_yield",
                (aggregate_all or {}).get("estimated_proton_missing_mass_yield"),
            ),
            "cleaned": mm_aggregate.get(
                "cleaned_missing_mass_yield",
                (aggregate_all or {}).get("cleaned_missing_mass_yield"),
            ),
            "proposed": {
                "raw": mm_aggregate.get("raw_missing_mass_yield"),
                "estimated_proton": mm_aggregate.get("estimated_proton_missing_mass_yield"),
                "cleaned_pre_rf": mm_aggregate.get("cleaned_missing_mass_yield"),
            },
            "applied": {
                "raw": applied_mm_aggregate.get("raw_missing_mass_yield"),
                "estimated_proton": applied_mm_aggregate.get("estimated_proton_missing_mass_yield"),
                "cleaned_pre_rf": applied_mm_aggregate.get("cleaned_missing_mass_yield"),
                "cleaned_final_rf": applied_mm_aggregate.get("final_rf_cleaned_missing_mass_yield"),
            },
        },
        "lambda_preservation_gate": lambda_gate,
        "proton_subtraction_mm": mm_diagnostics,
    }


def _resolve_timing_t_mm_display_range(config):
    """Return a valid MM display range without changing any analysis cuts."""
    mm_config = dict((config or {}).get("mm_diagnostics") or {})
    candidate = mm_config.get("display_range") or (config or {}).get("mm_validation_range")
    if not candidate or len(candidate) < 2:
        candidate = PROTON_CLEANING_EXACT_MM_VALIDATION_RANGE
    try:
        low, high = float(candidate[0]), float(candidate[1])
    except (TypeError, ValueError, IndexError):
        low, high = PROTON_CLEANING_EXACT_MM_VALIDATION_RANGE
    if not (math.isfinite(low) and math.isfinite(high) and high > low):
        low, high = PROTON_CLEANING_EXACT_MM_VALIDATION_RANGE
    return [float(low), float(high)]


def _resolve_lambda_preservation_gate_config(config):
    """Resolve the timing-t Lambda gate without introducing a second window."""
    gate = {
        "enabled": True,
        "validation_window_key": "lambda_peak",
        "maximum_lambda_removed_fraction": 0.10,
        "minimum_raw_prompt_events": 20,
        "minimum_absolute_support": None,
        "minimum_positive_signed_yield": None,
        "insufficient_support_policy": "bypass",
        "affects_event_weights": True,
        "affects_fit_acceptance": False,
    }
    gate.update(dict((config or {}).get("lambda_preservation_gate") or {}))
    policy = str(gate.get("insufficient_support_policy") or "").strip().lower()
    if policy != "bypass":
        raise ValueError(
            "Unsupported Lambda-preservation insufficient-support policy '{}'; "
            "only 'bypass' is supported".format(policy)
        )
    window_key = str(gate.get("validation_window_key") or "").strip()
    validation_windows = dict(
        (config or {}).get("validation_windows")
        or PROTON_CLEANING_EXACT_VALIDATION_WINDOWS
    )
    if window_key not in validation_windows:
        raise ValueError(
            "Lambda-preservation validation window key '{}' is unavailable in "
            "validation_windows".format(window_key)
        )
    bounds = validation_windows[window_key]
    try:
        window_low, window_high = float(bounds[0]), float(bounds[1])
    except (TypeError, ValueError, IndexError) as exc:
        raise ValueError(
            "Lambda-preservation validation window '{}' is invalid".format(window_key)
        ) from exc
    if not (
        math.isfinite(window_low)
        and math.isfinite(window_high)
        and window_high > window_low
    ):
        raise ValueError(
            "Lambda-preservation validation window '{}' is invalid".format(window_key)
        )
    maximum = _finite_float_or_none(gate.get("maximum_lambda_removed_fraction"))
    if maximum is None:
        raise ValueError("Lambda-preservation maximum removal fraction must be finite")
    prompt_minimum = int(gate.get("minimum_raw_prompt_events", 20) or 0)
    if prompt_minimum < 0:
        raise ValueError("Lambda-preservation minimum prompt count cannot be negative")
    for key in ("minimum_absolute_support", "minimum_positive_signed_yield"):
        value = gate.get(key)
        if value is None:
            continue
        finite_value = _finite_float_or_none(value)
        if finite_value is None or finite_value < 0.0:
            raise ValueError("Lambda-preservation {} must be a nonnegative finite value or None".format(key))
        gate[key] = finite_value
    gate.update(
        {
            "validation_window_key": window_key,
            "window_low": float(window_low),
            "window_high": float(window_high),
            "maximum_lambda_removed_fraction": float(maximum),
            "minimum_raw_prompt_events": int(prompt_minimum),
            "insufficient_support_policy": policy,
        }
    )
    return gate


def _lambda_gate_stage_values(row, stage):
    """Read proposed/applied frozen factors while tolerating legacy test rows."""
    prefix = "{}_".format(str(stage))
    probability = _finite_float_or_none(row.get(prefix + "proton_probability"))
    if probability is None:
        probability = _finite_float_or_none(row.get("proton_weight"))
    cleaned_factor = _finite_float_or_none(row.get(prefix + "cleaned_factor"))
    if cleaned_factor is None and probability is not None:
        cleaned_factor = 1.0 - probability
    final_cleaned_factor = _finite_float_or_none(
        row.get(prefix + "final_cleaned_factor")
    )
    if final_cleaned_factor is None:
        final_cleaned_factor = cleaned_factor
    return probability, cleaned_factor, final_cleaned_factor


def _evaluate_lambda_preservation_gate(cleaning_result, event_rows, config):
    """Make one pre-RF Lambda-preservation decision for the full setting."""
    gate = _resolve_lambda_preservation_gate_config(config)
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    result = {
        "enabled": bool(gate.get("enabled", True)),
        "validation_window_key": gate["validation_window_key"],
        "window_low": gate["window_low"],
        "window_high": gate["window_high"],
        "raw_prompt_count": 0,
        "raw_signed_yield": 0.0,
        "raw_absolute_support": 0.0,
        "raw_signed_to_absolute_support_ratio": None,
        "proposed_proton_yield": 0.0,
        "proposed_removed_fraction": None,
        "maximum_removed_fraction": gate["maximum_lambda_removed_fraction"],
        "minimum_raw_prompt_events": gate["minimum_raw_prompt_events"],
        "minimum_absolute_support": gate["minimum_absolute_support"],
        "minimum_positive_signed_yield": gate["minimum_positive_signed_yield"],
        "support_valid": False,
        "support_reasons": [],
        "observational_warnings": [],
        "status": "insufficient_support",
        "production_action": "bypass",
        "timing_fit_accepted": bool((cleaning_result or {}).get("accepted", False)),
        "setting_support_label": (diagnostics.get("setting_support") or {}).get(
            "support_label", SUPPORT_UNSUPPORTED
        ),
        "proton_cleaning_committed": False,
        "proposed_lookup_retained": True,
        "applied_lookup_zeroed": True,
        "pre_rf_evaluation": True,
        "source": "shared_nommcuts_frozen_lookup_rows",
        "lambda_rows_seen": 0,
        "lambda_rows_with_nonfinite_values": 0,
        "configuration": _json_ready_value(gate),
    }
    if not result["enabled"]:
        result.update(
            {
                "support_valid": None,
                "status": "disabled",
                "production_action": "apply",
                "proton_cleaning_committed": True,
                "applied_lookup_zeroed": False,
            }
        )
        return result

    raw_signed = raw_absolute = proposed_yield = 0.0
    prompt_count = nonfinite_rows = selected_rows = 0
    for row in event_rows or []:
        adj_mm = _finite_float_or_none((row or {}).get("adj_mm"))
        if adj_mm is None or not (gate["window_low"] <= adj_mm < gate["window_high"]):
            continue
        selected_rows += 1
        physical_weight = _finite_float_or_none((row or {}).get("physical_weight"))
        proposed_probability, _cleaned, _final = _lambda_gate_stage_values(row, "proposed")
        if physical_weight is None or proposed_probability is None:
            nonfinite_rows += 1
            continue
        raw_signed += physical_weight
        raw_absolute += abs(physical_weight)
        proposed_yield += physical_weight * proposed_probability
        prompt_count += int(bool((row or {}).get("is_prompt_source", False)))

    result.update(
        {
            "raw_prompt_count": int(prompt_count),
            "raw_signed_yield": float(raw_signed),
            "raw_absolute_support": float(raw_absolute),
            "proposed_proton_yield": float(proposed_yield),
            "lambda_rows_seen": int(selected_rows),
            "lambda_rows_with_nonfinite_values": int(nonfinite_rows),
        }
    )
    if math.isfinite(raw_signed) and math.isfinite(raw_absolute) and raw_absolute > 0.0:
        result["raw_signed_to_absolute_support_ratio"] = float(raw_signed / raw_absolute)
    if math.isfinite(raw_signed) and raw_signed > 0.0 and math.isfinite(proposed_yield):
        result["proposed_removed_fraction"] = float(proposed_yield / raw_signed)

    reasons = []
    if nonfinite_rows:
        reasons.append("nonfinite_lambda_row")
    if prompt_count < gate["minimum_raw_prompt_events"]:
        reasons.append("minimum_raw_prompt_events_not_met")
    if not (math.isfinite(raw_signed) and raw_signed > 0.0):
        reasons.append("raw_signed_yield_not_positive_or_nonfinite")
    if not (math.isfinite(raw_absolute) and raw_absolute > 0.0):
        reasons.append("raw_absolute_support_not_positive_or_nonfinite")
    if not math.isfinite(proposed_yield):
        reasons.append("proposed_proton_yield_nonfinite")
    if gate["minimum_absolute_support"] is not None and raw_absolute < gate["minimum_absolute_support"]:
        reasons.append("minimum_absolute_support_not_met")
    if gate["minimum_positive_signed_yield"] is not None and raw_signed < gate["minimum_positive_signed_yield"]:
        reasons.append("minimum_positive_signed_yield_not_met")
    fraction = result["proposed_removed_fraction"]
    if fraction is None or not math.isfinite(fraction):
        reasons.append("proposed_removed_fraction_unevaluable")
    warnings = []
    if math.isfinite(proposed_yield) and proposed_yield < 0.0:
        warnings.append("negative_proposed_lambda_yield")
    if fraction is not None and fraction < 0.0:
        warnings.append("negative_lambda_removed_fraction")
    result["support_reasons"] = sorted(set(reasons))
    result["observational_warnings"] = sorted(set(warnings))
    result["support_valid"] = not bool(reasons)
    if not result["support_valid"]:
        return result
    if fraction <= gate["maximum_lambda_removed_fraction"]:
        result.update(
            {
                "status": "pass",
                "production_action": "apply",
                "proton_cleaning_committed": True,
                "applied_lookup_zeroed": False,
            }
        )
    else:
        result.update({"status": "fail", "production_action": "bypass"})
    return result


def _finalize_lambda_gate_closure(event_rows, stage):
    raw_yield = proton_yield = cleaned_yield = final_rf_cleaned_yield = 0.0
    event_count = 0
    for row in event_rows or []:
        physical_weight = _finite_float_or_none((row or {}).get("physical_weight"))
        proton_probability, cleaned_factor, final_factor = _lambda_gate_stage_values(row, stage)
        if (
            physical_weight is None
            or proton_probability is None
            or cleaned_factor is None
            or final_factor is None
        ):
            continue
        event_count += 1
        raw_yield += physical_weight
        proton_yield += physical_weight * proton_probability
        cleaned_yield += physical_weight * cleaned_factor
        final_rf_cleaned_yield += physical_weight * final_factor
    return {
        "stage": str(stage),
        "event_count": int(event_count),
        "raw_signed_yield": float(raw_yield),
        "proton_signed_yield": float(proton_yield),
        "cleaned_pre_rf_signed_yield": float(cleaned_yield),
        "final_rf_cleaned_signed_yield": float(final_rf_cleaned_yield),
        "pre_rf_closure_difference": float(raw_yield - proton_yield - cleaned_yield),
    }


def _commit_lambda_preservation_gate(lookup, event_rows, gate_result):
    """Commit or bypass the proposed lookup without modifying its fit payload."""
    status = str((gate_result or {}).get("status") or "insufficient_support")
    commit = str((gate_result or {}).get("production_action") or "bypass") == "apply"
    by_signature = {}
    for signature, payload in (lookup or {}).items():
        proposed_probability, proposed_cleaned, proposed_final = _lambda_gate_stage_values(
            payload, "proposed"
        )
        proposed_probability = 0.0 if proposed_probability is None else proposed_probability
        proposed_cleaned = 1.0 - proposed_probability if proposed_cleaned is None else proposed_cleaned
        rf_accept = bool((payload or {}).get("rf_accept", True))
        if commit:
            applied_probability = proposed_probability
            applied_cleaned = proposed_cleaned
            applied_zero_reason = (payload or {}).get("event_lookup_zero_reason")
        else:
            applied_probability = 0.0
            applied_cleaned = 1.0
            applied_zero_reason = "lambda_preservation_gate_bypass"
        applied_final = applied_cleaned if rf_accept else 0.0
        payload.update(
            {
                "proposed_proton_probability": float(proposed_probability),
                "proposed_cleaned_factor": float(proposed_cleaned),
                "proposed_final_cleaned_factor": float(proposed_final if proposed_final is not None else (proposed_cleaned if rf_accept else 0.0)),
                "lambda_gate_status": status,
                "lambda_gate_passed": True if status == "pass" else (False if status in ("fail", "insufficient_support") else None),
                "applied_proton_probability": float(applied_probability),
                "applied_cleaned_factor": float(applied_cleaned),
                "applied_final_cleaned_factor": float(applied_final),
                "applied_zero_reason": applied_zero_reason,
                # Downstream production aliases must always be committed values.
                "proton_weight": float(applied_probability),
                "cleaned_factor": float(applied_cleaned),
                "final_cleaned_factor": float(applied_final),
            }
        )
        by_signature[str(signature)] = payload
    for row in event_rows or []:
        signature = _make_prepared_event_signature(
            row.get("source_label"), row.get("source_entry_index")
        )
        payload = by_signature.get(signature)
        if payload is None:
            continue
        row.update({
            key: payload.get(key)
            for key in (
                "proposed_proton_probability", "proposed_cleaned_factor",
                "proposed_final_cleaned_factor", "lambda_gate_status",
                "lambda_gate_passed", "applied_proton_probability",
                "applied_cleaned_factor", "applied_final_cleaned_factor",
                "applied_zero_reason", "proton_weight", "cleaned_factor",
                "final_cleaned_factor",
            )
        })
    gate_result["proposed_model_closure"] = {
        "label": "proposed timing-|t| model closure",
        "lookup_pre_rf_accounting": _finalize_lambda_gate_closure(event_rows, "proposed"),
        "cell_delta_closure_preserved": True,
    }
    gate_result["final_applied_closure"] = _finalize_lambda_gate_closure(
        event_rows, "applied"
    )
    return lookup


def _lambda_gate_audit_rows(event_rows, limit=LAMBDA_PRESERVATION_AUDIT_LIMIT):
    """Return a deterministic bounded event-level proposed/applied audit."""
    ranked = []
    for row in event_rows or []:
        source = str((row or {}).get("source_label", ""))
        try:
            entry = int((row or {}).get("source_entry_index", -1))
        except (TypeError, ValueError):
            entry = -1
        signature = _make_prepared_event_signature(source, entry)
        rank = hashlib.sha256(signature.encode("utf-8")).hexdigest()
        ranked.append((rank, signature, row))
    fields = (
        "source_label", "source_entry_index", "adj_t", "adj_mm", "selected_timing",
        "rf_accept", "proposed_proton_probability", "proposed_cleaned_factor",
        "proposed_final_cleaned_factor", "lambda_gate_status", "lambda_gate_passed",
        "applied_proton_probability", "applied_cleaned_factor",
        "applied_final_cleaned_factor", "applied_zero_reason",
    )
    rows = []
    for _rank, signature, row in sorted(ranked)[: max(0, int(limit))]:
        rows.append({
            "event_signature": signature,
            "source": row.get("source_label"),
            "entry_index": row.get("source_entry_index"),
            **{field: row.get(field) for field in fields},
        })
    return {
        "sampling_method": "sha256_ranked_source_entry_signature",
        "sample_limit": int(limit),
        "source_entry_count": int(len(event_rows or [])),
        "sampled_entry_count": int(len(rows)),
        "rows": rows,
    }


def _new_timing_t_mm_accumulator():
    return {
        "event_count": 0,
        "raw_prompt_event_count": 0,
        "signed_event_weight_sum": 0.0,
        "absolute_event_weight_support": 0.0,
        "raw_missing_mass_yield": 0.0,
        "estimated_proton_missing_mass_yield": 0.0,
        "cleaned_missing_mass_yield": 0.0,
        "final_rf_cleaned_missing_mass_yield": 0.0,
    }


def _fill_timing_t_mm_accumulator(accumulator, row, physical_weight, proton_weight, cleaned_factor, final_cleaned_factor):
    accumulator["event_count"] += 1
    accumulator["raw_prompt_event_count"] += int(bool((row or {}).get("is_prompt_source", False)))
    accumulator["signed_event_weight_sum"] += physical_weight
    accumulator["absolute_event_weight_support"] += abs(physical_weight)
    accumulator["raw_missing_mass_yield"] += physical_weight
    accumulator["estimated_proton_missing_mass_yield"] += physical_weight * proton_weight
    accumulator["cleaned_missing_mass_yield"] += physical_weight * cleaned_factor
    accumulator["final_rf_cleaned_missing_mass_yield"] += physical_weight * final_cleaned_factor


def _finalize_timing_t_mm_accumulator(accumulator, windows):
    summary = dict(accumulator)
    raw = float(summary["raw_missing_mass_yield"])
    proton = float(summary["estimated_proton_missing_mass_yield"])
    cleaned = float(summary["cleaned_missing_mass_yield"])
    summary["raw_minus_estimated_proton"] = raw - proton
    summary["pre_rf_cleaning_closure_difference"] = (raw - proton) - cleaned
    summary["windows"] = {}
    for name, values in (windows or {}).items():
        raw_yield = float(values["raw_missing_mass_yield"])
        proton_yield = float(values["estimated_proton_missing_mass_yield"])
        cleaned_yield = float(values["cleaned_missing_mass_yield"])
        removed_fraction = _safe_validation_ratio(proton_yield, raw_yield)
        summary["windows"][str(name)] = {
            "range": [float(values["range"][0]), float(values["range"][1])],
            "raw_yield": raw_yield,
            "estimated_proton_yield": proton_yield,
            "cleaned_yield": cleaned_yield,
            "final_rf_cleaned_yield": float(values["final_rf_cleaned_missing_mass_yield"]),
            "removed_yield": proton_yield,
            "removed_fraction": removed_fraction,
            "removed_fraction_valid": removed_fraction is not None,
            "raw_minus_estimated_proton": raw_yield - proton_yield,
            "pre_rf_cleaning_closure_difference": (raw_yield - proton_yield) - cleaned_yield,
        }
    return summary


def _build_timing_t_mm_root_payload(event_rows, mm_payload):
    """Create display-only MM ROOT products from frozen lookup rows."""
    payload_id = abs(id(mm_payload))
    t_edges = [float(value) for value in (mm_payload.get("t_edges") or [])]
    display_range = list(mm_payload.get("display_range") or [])
    try:
        mm_low, mm_high = float(display_range[0]), float(display_range[1])
    except (TypeError, ValueError, IndexError):
        return {}
    if len(t_edges) < 2 or not (math.isfinite(mm_low) and math.isfinite(mm_high) and mm_high > mm_low):
        return {}
    display_bins = max(1, int(mm_payload.get("display_bins", 100) or 100))
    maps = {}
    for key, title in (
        ("raw", "Raw selected MM;shifted MM [GeV];|t| [GeV^{2}]"),
        ("estimated_proton", "Estimated proton MM;shifted MM [GeV];|t| [GeV^{2}]"),
        ("cleaned_pre_rf", "Proton-cleaned MM before RF restoration;shifted MM [GeV];|t| [GeV^{2}]"),
        ("cleaned_final_rf", "Proton-cleaned MM after RF restoration;shifted MM [GeV];|t| [GeV^{2}]"),
    ):
        hist = ROOT.TH2D(
            "H_timing_t_mm_{}_{}".format(key, payload_id), title,
            display_bins, mm_low, mm_high,
            len(t_edges) - 1, array("d", t_edges),
        )
        hist.SetDirectory(0)
        hist.Sumw2()
        maps[key] = hist
    for row in event_rows or []:
        t_index = row.get("t_index")
        adj_mm = _finite_float_or_none(row.get("adj_mm"))
        physical_weight = _finite_float_or_none(row.get("physical_weight"))
        proton_weight = _finite_float_or_none(row.get("proton_weight"))
        cleaned_factor = _finite_float_or_none(row.get("cleaned_factor"))
        final_cleaned_factor = _finite_float_or_none(row.get("final_cleaned_factor"))
        if (
            not isinstance(t_index, int) or not (0 <= t_index < len(t_edges) - 1)
            or adj_mm is None or physical_weight is None or proton_weight is None
        ):
            continue
        cleaned_factor = 1.0 - proton_weight if cleaned_factor is None else cleaned_factor
        final_cleaned_factor = cleaned_factor if final_cleaned_factor is None else final_cleaned_factor
        t_value = _finite_float_or_none(row.get("adj_t"))
        if t_value is None:
            t_value = 0.5 * (t_edges[t_index] + t_edges[t_index + 1])
        maps["raw"].Fill(adj_mm, t_value, physical_weight)
        maps["estimated_proton"].Fill(adj_mm, t_value, physical_weight * proton_weight)
        maps["cleaned_pre_rf"].Fill(adj_mm, t_value, physical_weight * cleaned_factor)
        maps["cleaned_final_rf"].Fill(adj_mm, t_value, physical_weight * final_cleaned_factor)

    window_fraction_hists = {}
    for window_name in (mm_payload.get("validation_windows") or {}):
        hist = ROOT.TH1D(
            "H_timing_t_mm_{}_removed_fraction_{}".format(window_name, payload_id),
            "{} proton removal fraction;|t| [GeV^{{2}}];estimated proton / raw".format(window_name),
            len(t_edges) - 1, array("d", t_edges),
        )
        hist.SetDirectory(0)
        hist.Sumw2()
        for row in mm_payload.get("per_t_bin_summary") or []:
            t_index = row.get("t_index")
            window = (row.get("windows") or {}).get(window_name) or {}
            value = _finite_float_or_none(window.get("removed_fraction"))
            if isinstance(t_index, int) and value is not None and 0 <= t_index < hist.GetNbinsX():
                hist.SetBinContent(t_index + 1, value)
        window_fraction_hists[str(window_name)] = hist
    return {
        "maps": maps,
        "window_removed_fraction": window_fraction_hists,
        "object_count": int(len(maps) + len(window_fraction_hists)),
    }


def _build_timing_t_mm_diagnostics(cleaning_result, event_rows, config, *, include_root_payload=False):
    """Summarize proton-subtraction MM effects from frozen lookup rows only."""
    mm_config = dict((config or {}).get("mm_diagnostics") or {})
    if not bool(mm_config.get("enabled", True)):
        return {
            "enabled": False,
            "affects_event_weights": False,
            "affects_fit_acceptance": False,
            "configuration": _json_ready_value(mm_config),
        }
    t_edges = _finite_strictly_increasing_edges((cleaning_result or {}).get("t_edges") or [])
    if len(t_edges) < 2:
        return {
            "enabled": True,
            "affects_event_weights": False,
            "affects_fit_acceptance": False,
            "configuration": _json_ready_value(mm_config),
            "warnings": ["timing_t_mm_edges_unavailable"],
        }
    windows = {}
    for name, bounds in ((config or {}).get("validation_windows") or PROTON_CLEANING_EXACT_VALIDATION_WINDOWS).items():
        try:
            low, high = float(bounds[0]), float(bounds[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(low) and math.isfinite(high) and high > low:
            windows[str(name)] = [low, high]
    per_t = [_new_timing_t_mm_accumulator() for _ in range(len(t_edges) - 1)]
    per_t_windows = [
        {
            name: {**_new_timing_t_mm_accumulator(), "range": bounds}
            for name, bounds in windows.items()
        }
        for _ in per_t
    ]
    aggregate = _new_timing_t_mm_accumulator()
    aggregate_windows = {
        name: {**_new_timing_t_mm_accumulator(), "range": bounds}
        for name, bounds in windows.items()
    }
    exclusions = {
        "source_rows_seen": int(len(event_rows or [])),
        "invalid_t_index": 0,
        "missing_mm_or_weight_or_probability": 0,
        "outside_display_range": 0,
        "selected_rows": 0,
    }
    display_range = _resolve_timing_t_mm_display_range(config)
    for row in event_rows or []:
        t_index = row.get("t_index")
        adj_mm = _finite_float_or_none(row.get("adj_mm"))
        physical_weight = _finite_float_or_none(row.get("physical_weight"))
        proton_weight = _finite_float_or_none(row.get("proton_weight"))
        cleaned_factor = _finite_float_or_none(row.get("cleaned_factor"))
        final_cleaned_factor = _finite_float_or_none(row.get("final_cleaned_factor"))
        if not isinstance(t_index, int) or not (0 <= t_index < len(per_t)):
            exclusions["invalid_t_index"] += 1
            continue
        if adj_mm is None or physical_weight is None or proton_weight is None:
            exclusions["missing_mm_or_weight_or_probability"] += 1
            continue
        cleaned_factor = 1.0 - proton_weight if cleaned_factor is None else cleaned_factor
        final_cleaned_factor = cleaned_factor if final_cleaned_factor is None else final_cleaned_factor
        _fill_timing_t_mm_accumulator(
            per_t[t_index], row, physical_weight, proton_weight, cleaned_factor, final_cleaned_factor,
        )
        _fill_timing_t_mm_accumulator(
            aggregate, row, physical_weight, proton_weight, cleaned_factor, final_cleaned_factor,
        )
        for name, bounds in windows.items():
            if bounds[0] <= adj_mm < bounds[1]:
                _fill_timing_t_mm_accumulator(
                    per_t_windows[t_index][name], row, physical_weight, proton_weight,
                    cleaned_factor, final_cleaned_factor,
                )
                _fill_timing_t_mm_accumulator(
                    aggregate_windows[name], row, physical_weight, proton_weight,
                    cleaned_factor, final_cleaned_factor,
                )
        exclusions["selected_rows"] += 1
        if not (display_range[0] <= adj_mm <= display_range[1]):
            exclusions["outside_display_range"] += 1

    per_t_summary = []
    for index, accumulator in enumerate(per_t):
        per_t_summary.append({
            "t_index": int(index),
            "t_low": float(t_edges[index]),
            "t_high": float(t_edges[index + 1]),
            **_finalize_timing_t_mm_accumulator(accumulator, per_t_windows[index]),
        })
    payload = {
        "schema_version": 1,
        "enabled": True,
        "source": "frozen_timing_t_lookup_rows",
        "selection": "prepared_shared_nommcuts_without_final_missing_mass_window",
        "adj_mm_source": "shared_shifted_missing_mass",
        "affects_event_weights": False,
        "affects_fit_acceptance": False,
        "configuration": _json_ready_value(mm_config),
        "t_edges": [float(edge) for edge in t_edges],
        "display_range": display_range,
        "display_bins": max(1, int(mm_config.get("display_bins", 100) or 100)),
        "validation_windows": {name: list(bounds) for name, bounds in windows.items()},
        "aggregate": _finalize_timing_t_mm_accumulator(aggregate, aggregate_windows),
        "per_t_bin_summary": per_t_summary,
        "exclusions": exclusions,
    }
    if include_root_payload:
        payload["_root_payload"] = _build_timing_t_mm_root_payload(event_rows, payload)
    return payload


def _build_timing_t_diagnostic_integrity(
    cleaning_result, matrix_payload, coarse_by_t, validation_cfg, source_row_count,
):
    """Validate observational products without consulting production inputs."""
    tolerance = abs(float((validation_cfg or {}).get("diagnostic_integrity_tolerance", 1.0e-10) or 1.0e-10))
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    failures = []
    checks = []

    def _check(name, passed, **details):
        checks.append({"name": str(name), "passed": bool(passed), **details})
        if not passed:
            failures.append(str(name))

    t_edges = list((matrix_payload or {}).get("t_edges") or [])
    _check(
        "matrix_t_edges_match_cleaning_result",
        len(t_edges) == len((cleaning_result or {}).get("t_edges") or [])
        and all(abs(float(left) - float(right)) <= tolerance for left, right in zip(t_edges, (cleaning_result or {}).get("t_edges") or [])),
        tolerance=tolerance,
    )
    metrics = (matrix_payload or {}).get("metrics") or {}
    required = ("selected_prompt_count", "signed_physical_yield", "estimated_proton_yield", "average_proton_probability")
    for name in required:
        _check("matrix_metric_present_{}".format(name), name in metrics)
    prompt_matrix = metrics.get("selected_prompt_count") or []
    for t_index, accumulator in enumerate(coarse_by_t or []):
        observed = _sum_matrix_row(prompt_matrix, t_index)
        expected = float((accumulator or {}).get("raw_prompt_event_count", 0) or 0)
        _check(
            "prompt_total_matches_coarse_matrix_t{}".format(t_index),
            abs(observed - expected) <= tolerance,
            observed=observed, expected=expected,
        )
    for t_index, row in enumerate((matrix_payload or {}).get("metrics", {}).get("cleaned_yield") or []):
        raw_row = (metrics.get("signed_physical_yield") or [[]])[t_index] if t_index < len(metrics.get("signed_physical_yield") or []) else []
        removed_row = (metrics.get("estimated_proton_yield") or [[]])[t_index] if t_index < len(metrics.get("estimated_proton_yield") or []) else []
        for aero_index, cleaned in enumerate(row or []):
            raw = raw_row[aero_index] if aero_index < len(raw_row) else None
            removed = removed_row[aero_index] if aero_index < len(removed_row) else None
            if any(_finite_float_or_none(value) is None for value in (raw, removed, cleaned)):
                continue
            _check(
                "cleaned_equals_raw_minus_removed_t{}_a{}".format(t_index, aero_index),
                abs(float(cleaned) - (float(raw) - float(removed))) <= tolerance,
                raw=float(raw), removed=float(removed), cleaned=float(cleaned),
            )
    selected = diagnostics.get("selected_timing_candidate") or {}
    _check("selected_timing_candidate_present", bool(selected))
    if int(source_row_count or 0) > 0:
        _check(
            "nonempty_frozen_source_has_matrix_content",
            _matrix_has_nonzero_content(prompt_matrix) or _matrix_has_nonzero_content(metrics.get("signed_physical_yield") or []),
            source_row_count=int(source_row_count or 0),
        )
    return {
        "schema_version": 1,
        "source": "frozen_timing_t_lookup_diagnostic_validation",
        "strict": bool((validation_cfg or {}).get("diagnostic_strict", False)),
        "tolerance": tolerance,
        "source_row_count": int(source_row_count or 0),
        "checks": checks,
        "failures": failures,
        "valid": not failures,
        "matrix_exclusions": dict((matrix_payload or {}).get("exclusions") or {}),
    }


def _make_t_aerogel_hist(name, title, t_edges, aero_edges):
    hist = ROOT.TH2D(
        str(name), str(title), len(t_edges) - 1, array("d", t_edges),
        len(aero_edges) - 1, array("d", aero_edges),
    )
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def _signed_histogram_diagnostics(hist):
    diagnostics = {
        "positive_bin_count": 0,
        "negative_bin_count": 0,
        "positive_integral": 0.0,
        "negative_integral": 0.0,
        "absolute_negative_integral": 0.0,
        "signed_integral": 0.0,
        "absolute_integral": 0.0,
    }
    if hist is None:
        return diagnostics
    for x_bin in range(1, hist.GetNbinsX() + 1):
        for y_bin in range(1, hist.GetNbinsY() + 1):
            value = float(hist.GetBinContent(x_bin, y_bin))
            if value > 0.0:
                diagnostics["positive_bin_count"] += 1
                diagnostics["positive_integral"] += value
            elif value < 0.0:
                diagnostics["negative_bin_count"] += 1
                diagnostics["negative_integral"] += value
                diagnostics["absolute_negative_integral"] += abs(value)
            diagnostics["signed_integral"] += value
            diagnostics["absolute_integral"] += abs(value)
    return diagnostics


def _build_t_aerogel_root_payload(
    cleaning_result, event_rows, validation_cfg, t_edges, matrix_payload=None,
):
    """Create display objects from frozen lookup rows only; never refit events."""
    display_edges = _resolve_aerogel_display_edges(validation_cfg)
    if len(display_edges) < 2:
        return {"global": {}, "per_t": [], "object_count": 0}
    diagnostics = cleaning_result.get("diagnostics") or {}
    selected = diagnostics.get("selected_timing_candidate") or {}
    timing_range = selected.get("display_range") or (-2.0, 2.0)
    try:
        timing_low, timing_high = float(timing_range[0]), float(timing_range[1])
        timing_bins = max(1, int(selected.get("histogram_bins", 131) or 131))
    except (TypeError, ValueError, IndexError):
        timing_low, timing_high, timing_bins = -2.0, 2.0, 131
    if not (math.isfinite(timing_low) and math.isfinite(timing_high) and timing_high > timing_low):
        timing_low, timing_high = -2.0, 2.0

    global_hists = {}
    if bool(validation_cfg.get("write_global_aero_vs_t_pages", True)):
        global_hists = {
            "H_aero_vs_t_raw_prompt": _make_t_aerogel_hist(
                "H_aero_vs_t_raw_prompt", "Raw selected prompt events;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_signed_physical": _make_t_aerogel_hist(
                "H_aero_vs_t_signed_physical", "Signed physical source yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_signed_physical_positive": _make_t_aerogel_hist(
                "H_aero_vs_t_signed_physical_positive", "Positive physical source yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_signed_physical_negative": _make_t_aerogel_hist(
                "H_aero_vs_t_signed_physical_negative", "Absolute negative physical source yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_estimated_proton": _make_t_aerogel_hist(
                "H_aero_vs_t_estimated_proton", "Proposed estimated proton physical yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_estimated_proton_positive": _make_t_aerogel_hist(
                "H_aero_vs_t_estimated_proton_positive", "Positive estimated proton yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_estimated_proton_negative_abs": _make_t_aerogel_hist(
                "H_aero_vs_t_estimated_proton_negative_abs", "Absolute negative estimated proton yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_proton_cleaned": _make_t_aerogel_hist(
                "H_aero_vs_t_proton_cleaned", "Proposed proton-cleaned physical yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_proton_cleaned_positive": _make_t_aerogel_hist(
                "H_aero_vs_t_proton_cleaned_positive", "Positive frozen cleaned yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
            "H_aero_vs_t_proton_cleaned_negative_abs": _make_t_aerogel_hist(
                "H_aero_vs_t_proton_cleaned_negative_abs", "Absolute negative frozen cleaned yield;|t| [GeV^{2}];P_aero NPE",
                t_edges, display_edges,
            ),
        }

    per_t = []
    if bool(validation_cfg.get("write_per_t_pid_pages", True)):
        for t_index in range(len(t_edges) - 1):
            tag = "t{}".format(t_index)
            row = {
                "t_index": int(t_index),
                "raw_prompt_timing_vs_aero": ROOT.TH2D(
                    "H_aero_timing_raw_prompt_{}".format(tag),
                    "Raw prompt timing versus aerogel;P_aero NPE;selected timing [ns]",
                    len(display_edges) - 1, array("d", display_edges), timing_bins, timing_low, timing_high,
                ),
                "estimated_proton_timing_vs_aero": ROOT.TH2D(
                    "H_aero_timing_estimated_proton_{}".format(tag),
                    "Proposed estimated-proton timing versus aerogel;P_aero NPE;selected timing [ns]",
                    len(display_edges) - 1, array("d", display_edges), timing_bins, timing_low, timing_high,
                ),
                "raw_signed_projection": _make_aero_axis_hist(
                    "H_aero_raw_signed_{}".format(tag), "Signed physical yield;P_aero NPE;yield", display_edges,
                ),
                "estimated_proton_projection": _make_aero_axis_hist(
                    "H_aero_estimated_proton_{}".format(tag), "Proposed estimated proton yield;P_aero NPE;yield", display_edges,
                ),
                "cleaned_projection": _make_aero_axis_hist(
                    "H_aero_cleaned_{}".format(tag), "Proposed cleaned yield;P_aero NPE;yield", display_edges,
                ),
                "average_proton_probability": _make_aero_axis_hist(
                    "H_aero_average_proton_probability_{}".format(tag), "Average proton probability;P_aero NPE;#LTw_{p}#GT", display_edges,
                ),
                "_absolute_support": _make_aero_axis_hist(
                    "H_aero_absolute_support_{}".format(tag), "", display_edges,
                ),
                "_absolute_probability": _make_aero_axis_hist(
                    "H_aero_absolute_probability_{}".format(tag), "", display_edges,
                ),
            }
            for hist in row.values():
                if _is_root_object(hist):
                    hist.SetDirectory(0)
                    hist.Sumw2()
            if bool(validation_cfg.get("write_full_per_t_pid_pages", False)):
                row["full_raw_prompt_timing_vs_aero"] = _clone_hist(
                    row["raw_prompt_timing_vs_aero"], "H_aero_full_raw_prompt_{}".format(tag), reset=True,
                )
                row["full_signed_physical_timing_vs_aero"] = _clone_hist(
                    row["estimated_proton_timing_vs_aero"], "H_aero_full_signed_physical_{}".format(tag), reset=True,
                )
                row["full_estimated_proton_timing_vs_aero"] = _clone_hist(
                    row["estimated_proton_timing_vs_aero"], "H_aero_full_estimated_proton_{}".format(tag), reset=True,
                )
                row["full_cleaned_timing_vs_aero"] = _clone_hist(
                    row["estimated_proton_timing_vs_aero"], "H_aero_full_cleaned_{}".format(tag), reset=True,
                )
            per_t.append(row)

    for event in event_rows:
        t_index = event.get("t_index")
        if not isinstance(t_index, int) or not (0 <= t_index < len(t_edges) - 1):
            continue
        aero_value = _finite_float_or_none(event.get("aero_value"))
        physical_weight = _finite_float_or_none(event.get("physical_weight"))
        proton_weight = _finite_float_or_none(event.get("proton_weight"))
        cleaned_factor = _finite_float_or_none(event.get("cleaned_factor"))
        timing_value = _finite_float_or_none(event.get("selected_timing"))
        if aero_value is None or physical_weight is None or proton_weight is None:
            continue
        cleaned_factor = 1.0 - proton_weight if cleaned_factor is None else cleaned_factor
        estimated = physical_weight * proton_weight
        cleaned = physical_weight * cleaned_factor
        if global_hists:
            if bool(event.get("is_prompt_source", False)):
                global_hists["H_aero_vs_t_raw_prompt"].Fill(event["adj_t"], aero_value, 1.0)
            global_hists["H_aero_vs_t_signed_physical"].Fill(event["adj_t"], aero_value, physical_weight)
            if physical_weight > 0.0:
                global_hists["H_aero_vs_t_signed_physical_positive"].Fill(event["adj_t"], aero_value, physical_weight)
            elif physical_weight < 0.0:
                global_hists["H_aero_vs_t_signed_physical_negative"].Fill(event["adj_t"], aero_value, abs(physical_weight))
            global_hists["H_aero_vs_t_estimated_proton"].Fill(event["adj_t"], aero_value, estimated)
            global_hists["H_aero_vs_t_proton_cleaned"].Fill(event["adj_t"], aero_value, cleaned)
            if estimated > 0.0:
                global_hists["H_aero_vs_t_estimated_proton_positive"].Fill(event["adj_t"], aero_value, estimated)
            elif estimated < 0.0:
                global_hists["H_aero_vs_t_estimated_proton_negative_abs"].Fill(event["adj_t"], aero_value, abs(estimated))
            if cleaned > 0.0:
                global_hists["H_aero_vs_t_proton_cleaned_positive"].Fill(event["adj_t"], aero_value, cleaned)
            elif cleaned < 0.0:
                global_hists["H_aero_vs_t_proton_cleaned_negative_abs"].Fill(event["adj_t"], aero_value, abs(cleaned))
        if t_index >= len(per_t):
            continue
        row = per_t[t_index]
        row["raw_signed_projection"].Fill(aero_value, physical_weight)
        row["estimated_proton_projection"].Fill(aero_value, estimated)
        row["cleaned_projection"].Fill(aero_value, cleaned)
        row["_absolute_support"].Fill(aero_value, abs(physical_weight))
        row["_absolute_probability"].Fill(aero_value, abs(physical_weight) * proton_weight)
        if timing_value is not None:
            if bool(event.get("is_prompt_source", False)):
                row["raw_prompt_timing_vs_aero"].Fill(aero_value, timing_value, 1.0)
            row["estimated_proton_timing_vs_aero"].Fill(aero_value, timing_value, estimated)
            if "full_raw_prompt_timing_vs_aero" in row:
                if bool(event.get("is_prompt_source", False)):
                    row["full_raw_prompt_timing_vs_aero"].Fill(aero_value, timing_value, 1.0)
                row["full_signed_physical_timing_vs_aero"].Fill(aero_value, timing_value, physical_weight)
                row["full_estimated_proton_timing_vs_aero"].Fill(aero_value, timing_value, estimated)
                row["full_cleaned_timing_vs_aero"].Fill(aero_value, timing_value, cleaned)
    hidden_accumulator_objects_released = 0
    for row in per_t:
        for bin_index in range(1, row["average_proton_probability"].GetNbinsX() + 1):
            support = float(row["_absolute_support"].GetBinContent(bin_index))
            numerator = float(row["_absolute_probability"].GetBinContent(bin_index))
            if support > 0.0:
                row["average_proton_probability"].SetBinContent(bin_index, numerator / support)
        # These ROOT histograms are transient implementation accumulators;
        # the frozen average is the only retained visual product.
        row.pop("_absolute_support", None)
        row.pop("_absolute_probability", None)
        hidden_accumulator_objects_released += 2

    coarse_hists = {}
    coarse_edges = list((matrix_payload or {}).get("aero_edges") or _resolve_aerogel_summary_edges(validation_cfg))
    if len(coarse_edges) >= 2 and (matrix_payload or {}).get("metrics"):
        metric_specs = (
            ("selected_prompt_count", "H_t_aero_selected_prompt_count", False),
            ("signed_physical_yield", "H_t_aero_signed_physical_yield", False),
            ("estimated_proton_yield", "H_t_aero_estimated_proton_yield", False),
            ("average_proton_probability", "H_t_aero_average_proton_probability", True),
            ("low_mm_removed_yield", "H_t_aero_low_mm_removed_yield", False),
            ("low_mm_removed_fraction", "H_t_aero_low_mm_removed_fraction", True),
            ("lambda_removed_yield", "H_t_aero_lambda_removed_yield", False),
            ("lambda_removed_fraction", "H_t_aero_lambda_removed_fraction", True),
        )
        for metric_key, base_name, has_validity_mask in metric_specs:
            value_hist = _make_t_aerogel_hist(
                base_name, "{};|t| [GeV^{{2}}];P_aero NPE".format(metric_key),
                t_edges, coarse_edges,
            )
            valid_hist = None
            if has_validity_mask:
                valid_hist = _make_t_aerogel_hist(
                    "{}_valid".format(base_name), "{} validity;|t| [GeV^{{2}}];P_aero NPE".format(metric_key),
                    t_edges, coarse_edges,
                )
            values = _matrix_metric(matrix_payload, metric_key)
            validity = ((matrix_payload or {}).get("validity_masks") or {}).get(metric_key) or []
            for t_index, value_row in enumerate(values):
                for aero_index, value in enumerate(value_row or []):
                    if value is not None:
                        value_hist.SetBinContent(t_index + 1, aero_index + 1, float(value))
                    if valid_hist is not None:
                        is_valid = (
                            bool(validity[t_index][aero_index])
                            if t_index < len(validity) and aero_index < len(validity[t_index] or [])
                            else value is not None
                        )
                        if is_valid:
                            valid_hist.SetBinContent(t_index + 1, aero_index + 1, 1.0)
            coarse_hists[base_name] = value_hist
            if valid_hist is not None:
                coarse_hists["{}_valid".format(base_name)] = valid_hist
    root_objects = list(global_hists.values())
    root_objects.extend(
        hist for row in per_t for hist in row.values() if _is_root_object(hist)
    )
    root_objects.extend(coarse_hists.values())
    return {
        "global": global_hists,
        "per_t": per_t,
        "coarse": coarse_hists,
        "display_aero_edges": display_edges,
        "timing_range": [timing_low, timing_high],
        "timing_bins": int(timing_bins),
        "object_count": int(len(root_objects)),
        "hidden_accumulator_objects_released": int(hidden_accumulator_objects_released),
        "signed_histogram_diagnostics": {
            "signed_physical": _signed_histogram_diagnostics(global_hists.get("H_aero_vs_t_signed_physical")),
            "estimated_proton": _signed_histogram_diagnostics(global_hists.get("H_aero_vs_t_estimated_proton")),
            "proton_cleaned": _signed_histogram_diagnostics(global_hists.get("H_aero_vs_t_proton_cleaned")),
        },
    }


def _build_t_aerogel_validation(
    cleaning_result, source_bundle, event_rows, config, *, include_root_payload=False,
):
    """Summarize aerogel only after timing-t probabilities are frozen.

    This function intentionally consumes prepared lookup rows, not trees or
    fitting inputs.  Its configuration is therefore incapable of changing a
    timing candidate, cell validity, proton probability, or cleaned factor.
    """
    validation_cfg = dict(config.get("aerogel_validation") or {})
    if not bool(validation_cfg.get("enabled", True)):
        return {
            "enabled": False,
            "affects_event_weights": False,
            "affects_fit_acceptance": False,
            "configuration": _json_ready_value(validation_cfg),
        }
    t_edges = _finite_strictly_increasing_edges(cleaning_result.get("t_edges") or [])
    aero_edges = _resolve_aerogel_summary_edges(validation_cfg)
    if len(t_edges) < 2 or len(aero_edges) < 2:
        return {
            "enabled": True,
            "affects_event_weights": False,
            "affects_fit_acceptance": False,
            "warnings": ["aerogel_validation_edges_unavailable"],
            "configuration": _json_ready_value(validation_cfg),
        }
    metric_names = (
        "event_count", "raw_prompt_event_count", "signed_event_weight_sum",
        "absolute_event_weight_support", "sum_proton_probability",
        "absolute_weighted_probability_sum", "raw_missing_mass_yield",
        "estimated_proton_missing_mass_yield", "cleaned_missing_mass_yield",
        "low_mm_raw_yield", "low_mm_removed_yield", "low_mm_cleaned_yield",
        "lambda_raw_yield", "lambda_removed_yield", "lambda_cleaned_yield",
    )
    count_metrics = {"event_count", "raw_prompt_event_count"}
    cells = [
        [{name: (0 if name in count_metrics else 0.0) for name in metric_names}
         for _ in range(len(aero_edges) - 1)]
        for _ in range(len(t_edges) - 1)
    ]
    windows = dict(config.get("validation_windows") or {})
    low_mm = tuple(windows.get("low_mm") or (0.80, 0.90))
    lambda_peak = tuple(windows.get("lambda_peak") or (1.105, 1.125))

    def _new_accumulator():
        return {name: (0 if name in count_metrics else 0.0) for name in metric_names}

    def _fill_accumulator(accumulator, row, physical_weight, proton_weight, cleaned_factor, adj_mm):
        accumulator["event_count"] += 1
        if bool(row.get("is_prompt_source", False)):
            accumulator["raw_prompt_event_count"] += 1
        accumulator["signed_event_weight_sum"] += physical_weight
        accumulator["absolute_event_weight_support"] += abs(physical_weight)
        accumulator["sum_proton_probability"] += physical_weight * proton_weight
        accumulator["absolute_weighted_probability_sum"] += abs(physical_weight) * proton_weight
        accumulator["raw_missing_mass_yield"] += physical_weight
        accumulator["estimated_proton_missing_mass_yield"] += physical_weight * proton_weight
        accumulator["cleaned_missing_mass_yield"] += physical_weight * cleaned_factor
        if adj_mm is not None and low_mm[0] <= adj_mm < low_mm[1]:
            accumulator["low_mm_raw_yield"] += physical_weight
            accumulator["low_mm_removed_yield"] += physical_weight * proton_weight
            accumulator["low_mm_cleaned_yield"] += physical_weight * cleaned_factor
        if adj_mm is not None and lambda_peak[0] <= adj_mm < lambda_peak[1]:
            accumulator["lambda_raw_yield"] += physical_weight
            accumulator["lambda_removed_yield"] += physical_weight * proton_weight
            accumulator["lambda_cleaned_yield"] += physical_weight * cleaned_factor

    def _finalize_accumulator(accumulator):
        summary = dict(accumulator)
        summary["average_proton_probability"] = _safe_validation_ratio(
            summary["absolute_weighted_probability_sum"], summary["absolute_event_weight_support"],
        )
        summary["low_mm_removed_fraction"] = _safe_validation_ratio(
            summary["low_mm_removed_yield"], summary["low_mm_raw_yield"],
        )
        summary["lambda_removed_fraction"] = _safe_validation_ratio(
            summary["lambda_removed_yield"], summary["lambda_raw_yield"],
        )
        # Preserve an observed zero while making an absent denominator
        # distinguishable in both JSON and the CSV validity columns.
        for metric in (
            "average_proton_probability",
            "low_mm_removed_fraction",
            "lambda_removed_fraction",
        ):
            summary["{}_valid".format(metric)] = summary[metric] is not None
        return summary

    threshold_by_t = [
        {"low": _new_accumulator(), "high": _new_accumulator()}
        for _ in range(len(t_edges) - 1)
    ]
    threshold_setting = {"low": _new_accumulator(), "high": _new_accumulator()}
    all_by_t = [_new_accumulator() for _ in range(len(t_edges) - 1)]
    # These aggregates represent exactly the rows admitted to a coarse
    # aerogel bin.  They are the reference totals for matrix integrity, while
    # ``all_by_t`` intentionally includes rows outside display coverage.
    coarse_by_t = [_new_accumulator() for _ in range(len(t_edges) - 1)]
    all_setting = _new_accumulator()
    low_limit = float(validation_cfg.get("low_reference_max_npe", 5.0))
    high_limit = float(validation_cfg.get("high_reference_min_npe", 10.0))
    source_stats = {}
    signed_positive_events = signed_negative_events = 0
    signed_positive_integral = signed_negative_integral = 0.0
    exclusions = {
        "source_rows_seen": int(len(event_rows or [])),
        "invalid_t_index": 0,
        "missing_aerogel_or_weight_or_probability": 0,
        "outside_coarse_aerogel_edges": 0,
        "coarse_rows": 0,
    }
    for row in event_rows:
        t_index = row.get("t_index")
        if not isinstance(t_index, int) or not (0 <= t_index < len(cells)):
            exclusions["invalid_t_index"] += 1
            continue
        aero_value = _finite_float_or_none(row.get("aero_value"))
        physical_weight = _finite_float_or_none(row.get("physical_weight"))
        proton_weight = _finite_float_or_none(row.get("proton_weight"))
        cleaned_factor = _finite_float_or_none(row.get("cleaned_factor"))
        adj_mm = _finite_float_or_none(row.get("adj_mm"))
        if aero_value is None or physical_weight is None or proton_weight is None:
            exclusions["missing_aerogel_or_weight_or_probability"] += 1
            continue
        cleaned_factor = 1.0 - proton_weight if cleaned_factor is None else cleaned_factor
        # Per-t totals use every valid frozen row; display/coarse membership
        # remains an independent, optional projection below.
        _fill_accumulator(all_by_t[t_index], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
        _fill_accumulator(all_setting, row, physical_weight, proton_weight, cleaned_factor, adj_mm)
        if aero_value < low_limit:
            _fill_accumulator(threshold_by_t[t_index]["low"], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
            _fill_accumulator(threshold_setting["low"], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
        if aero_value >= high_limit:
            _fill_accumulator(threshold_by_t[t_index]["high"], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
            _fill_accumulator(threshold_setting["high"], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
        aero_index = _find_collection_bin(aero_value, aero_edges)
        if 0 <= aero_index < len(cells[t_index]):
            _fill_accumulator(cells[t_index][aero_index], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
            _fill_accumulator(coarse_by_t[t_index], row, physical_weight, proton_weight, cleaned_factor, adj_mm)
            exclusions["coarse_rows"] += 1
        else:
            exclusions["outside_coarse_aerogel_edges"] += 1
        source_name = str(row.get("source_label", "unknown"))
        source_stats.setdefault(source_name, {"event_count": 0, "raw_prompt_event_count": 0, "signed_physical_yield": 0.0})
        source_stats[source_name]["event_count"] += 1
        source_stats[source_name]["raw_prompt_event_count"] += int(bool(row.get("is_prompt_source", False)))
        source_stats[source_name]["signed_physical_yield"] += physical_weight
        if physical_weight > 0.0:
            signed_positive_events += 1
            signed_positive_integral += physical_weight
        elif physical_weight < 0.0:
            signed_negative_events += 1
            signed_negative_integral += physical_weight

    average_probability = [
        [_safe_validation_ratio(cell["absolute_weighted_probability_sum"], cell["absolute_event_weight_support"]) for cell in row]
        for row in cells
    ]
    low_removed_fraction = [
        [_safe_validation_ratio(cell["low_mm_removed_yield"], cell["low_mm_raw_yield"]) for cell in row]
        for row in cells
    ]
    lambda_removed_fraction = [
        [_safe_validation_ratio(cell["lambda_removed_yield"], cell["lambda_raw_yield"]) for cell in row]
        for row in cells
    ]
    for t_index, row in enumerate(cells):
        for aero_index, cell in enumerate(row):
            cell["average_proton_probability"] = average_probability[t_index][aero_index]
            cell["average_proton_probability_valid"] = cell["average_proton_probability"] is not None
            cell["low_mm_removed_fraction"] = low_removed_fraction[t_index][aero_index]
            cell["low_mm_removed_fraction_valid"] = cell["low_mm_removed_fraction"] is not None
            cell["lambda_removed_fraction"] = lambda_removed_fraction[t_index][aero_index]
            cell["lambda_removed_fraction_valid"] = cell["lambda_removed_fraction"] is not None
    matrix_payload = _build_t_aerogel_matrix_payload(
        cells, t_edges, aero_edges, exclusions,
    )
    matrices = {name: _matrix_from_t_aerogel_cells(cells, name) for name in metric_names}
    minimum_events = max(0, int(validation_cfg.get("minimum_events_per_t_bin", 20) or 20))
    maximum_high_weight = float(validation_cfg.get("maximum_high_aero_average_weight", 0.10))
    maximum_ratio = float(validation_cfg.get("maximum_high_to_low_weight_ratio", 1.0))
    maximum_high_lambda = float(validation_cfg.get("maximum_high_aero_lambda_removed_fraction", 0.10))
    minimum_low_probability = float(validation_cfg.get("minimum_low_aero_proton_fraction", 0.50))

    all_warnings = []
    per_t_summary = []
    warnings_by_t_bin = []
    for t_index, row in enumerate(cells):
        total = _finalize_accumulator(all_by_t[t_index])
        # Threshold accumulators are filled directly from event-level NPE.
        # They deliberately do not depend on the coarse matrix boundaries.
        low = _finalize_accumulator(threshold_by_t[t_index]["low"])
        high = _finalize_accumulator(threshold_by_t[t_index]["high"])
        labels = []
        if total["raw_prompt_event_count"] < minimum_events:
            labels.append("insufficient_aerogel_support")
        low_average = low["average_proton_probability"]
        high_average = high["average_proton_probability"]
        high_to_low = _safe_validation_ratio(high_average, low_average) if low_average not in (None, 0.0) else None
        if high_average is not None and high_average > maximum_high_weight:
            labels.extend(("high_aero_estimated_proton_fraction_large", "high_aero_average_weight_exceeds_threshold"))
        if low_average is not None and low_average < minimum_low_probability:
            labels.append("proton_probability_not_concentrated_at_low_aero")
        if high_to_low is not None and high_to_low > maximum_ratio:
            labels.append("high_to_low_aero_weight_ratio_exceeds_threshold")
        if bool(validation_cfg.get("warn_if_high_fraction_exceeds_low_fraction", True)) and (
            low_average is not None and high_average is not None and high_average > low_average
        ):
            labels.append("high_aero_average_weight_exceeds_low_aero")
        if high["lambda_removed_fraction"] is not None and high["lambda_removed_fraction"] > maximum_high_lambda:
            labels.append("high_aero_lambda_removal_exceeds_threshold")
        legacy_lambda_excess = float(validation_cfg.get("high_aero_lambda_removal_excess_threshold", 0.0))
        if (
            low["lambda_removed_fraction"] is not None and high["lambda_removed_fraction"] is not None
            and high["lambda_removed_fraction"] - low["lambda_removed_fraction"] > legacy_lambda_excess
        ):
            labels.append("lambda_region_removal_large_in_high_aero")
        if (
            low["low_mm_removed_fraction"] is not None and high["low_mm_removed_fraction"] is not None
            and high["low_mm_removed_fraction"] > low["low_mm_removed_fraction"]
        ):
            labels.append("low_mm_removal_not_correlated_with_low_aerogel")
        if labels:
            labels.append("t_bin_aerogel_validation_inconsistent")
            all_warnings.extend(labels)
        per_t_summary.append({
            "t_index": int(t_index), "t_low": float(t_edges[t_index]), "t_high": float(t_edges[t_index + 1]),
            **total,
            "low_aero": low, "high_aero": high,
            "low_aero_raw_prompt_count": low["raw_prompt_event_count"],
            "high_aero_raw_prompt_count": high["raw_prompt_event_count"],
            "low_aero_signed_yield": low["signed_event_weight_sum"],
            "high_aero_signed_yield": high["signed_event_weight_sum"],
            "low_aero_absolute_support": low["absolute_event_weight_support"],
            "high_aero_absolute_support": high["absolute_event_weight_support"],
            "low_aero_estimated_proton_yield": low["estimated_proton_missing_mass_yield"],
            "high_aero_estimated_proton_yield": high["estimated_proton_missing_mass_yield"],
            "low_aero_cleaned_yield": low["cleaned_missing_mass_yield"],
            "high_aero_cleaned_yield": high["cleaned_missing_mass_yield"],
            "low_aero_average_proton_probability": low["average_proton_probability"],
            "high_aero_average_proton_probability": high["average_proton_probability"],
            "low_aero_average_weight": low_average,
            "high_aero_average_weight": high_average,
            "high_to_low_aero_weight_ratio": high_to_low,
            "low_aero_low_mm_removed_fraction": low["low_mm_removed_fraction"],
            "high_aero_low_mm_removed_fraction": high["low_mm_removed_fraction"],
            "low_aero_lambda_removed_fraction": low["lambda_removed_fraction"],
            "high_aero_lambda_removed_fraction": high["lambda_removed_fraction"],
            "warnings": sorted(set(labels)),
        })
        warnings_by_t_bin.append({
            "t_index": int(t_index), "t_low": float(t_edges[t_index]), "t_high": float(t_edges[t_index + 1]),
            "low_aero_average_weight": low_average, "high_aero_average_weight": high_average,
            "low_aero_average_weight_valid": low_average is not None,
            "high_aero_average_weight_valid": high_average is not None,
            "low_aero_lambda_removed_fraction": low["lambda_removed_fraction"],
            "high_aero_lambda_removed_fraction": high["lambda_removed_fraction"],
            "low_aero_lambda_removed_fraction_valid": low["lambda_removed_fraction"] is not None,
            "high_aero_lambda_removed_fraction_valid": high["lambda_removed_fraction"] is not None,
            "warnings": sorted(set(labels)),
        })
    aggregate_low = _finalize_accumulator(threshold_setting["low"])
    aggregate_high = _finalize_accumulator(threshold_setting["high"])
    aggregate_all = _finalize_accumulator(all_setting)
    if (
        aggregate_low["average_proton_probability"] is not None
        and aggregate_high["average_proton_probability"] is not None
        and aggregate_high["average_proton_probability"] > aggregate_low["average_proton_probability"]
    ):
        all_warnings.append("high_aero_proton_fraction_exceeds_low_aero")

    payload = {
        "enabled": True,
        "affects_event_weights": False,
        "affects_fit_acceptance": False,
        "configuration": _json_ready_value(validation_cfg),
        "configuration_observational_only": True,
        "low_reference_max_npe": low_limit,
        "high_reference_min_npe": high_limit,
        "t_edges": t_edges,
        "aero_edges": aero_edges,
        "display_aero_edges": _resolve_aerogel_display_edges(validation_cfg),
        "source_tree_provenance": source_stats,
        "matrix_payload": matrix_payload,
        # New explicit names are the authoritative compact-matrix interface.
        "selected_prompt_count_by_t_aero": _matrix_metric(matrix_payload, "selected_prompt_count"),
        "signed_physical_yield_by_t_aero": _matrix_metric(matrix_payload, "signed_physical_yield"),
        "estimated_proton_yield_by_t_aero": _matrix_metric(matrix_payload, "estimated_proton_yield"),
        "cleaned_yield_by_t_aero": _matrix_metric(matrix_payload, "cleaned_yield"),
        "event_count_by_t_aero": matrices["event_count"],
        # Compatibility aliases intentionally reference the same named matrix
        # payload rather than a second independently-computed set of arrays.
        "raw_prompt_event_count_by_t_aero": _matrix_metric(matrix_payload, "selected_prompt_count"),
        "signed_event_weight_sum_by_t_aero": _matrix_metric(matrix_payload, "signed_physical_yield"),
        "absolute_event_weight_support_by_t_aero": _matrix_metric(matrix_payload, "absolute_physical_support"),
        "proton_probability_sum_by_t_aero": _matrix_metric(matrix_payload, "estimated_proton_yield"),
        "average_proton_probability_by_t_aero": _matrix_metric(matrix_payload, "average_proton_probability"),
        "low_mm_removed_yield_by_t_aero": _matrix_metric(matrix_payload, "low_mm_removed_yield"),
        "low_mm_removed_fraction_by_t_aero": _matrix_metric(matrix_payload, "low_mm_removed_fraction"),
        "lambda_removed_yield_by_t_aero": _matrix_metric(matrix_payload, "lambda_removed_yield"),
        "lambda_removed_fraction_by_t_aero": _matrix_metric(matrix_payload, "lambda_removed_fraction"),
        "cells": cells,
        "per_t_bin_summary": per_t_summary,
        "proton_fraction_below_5_npe": aggregate_low["average_proton_probability"],
        "proton_fraction_above_10_npe": aggregate_high["average_proton_probability"],
        "average_weight_below_5_npe": aggregate_low["average_proton_probability"],
        "average_weight_above_10_npe": aggregate_high["average_proton_probability"],
        "aggregate_low_aero": aggregate_low,
        "aggregate_high_aero": aggregate_high,
        "aggregate_all_frozen_rows": aggregate_all,
        "low_aero_raw_prompt_count": aggregate_low["raw_prompt_event_count"],
        "high_aero_raw_prompt_count": aggregate_high["raw_prompt_event_count"],
        "low_aero_signed_yield": aggregate_low["signed_event_weight_sum"],
        "high_aero_signed_yield": aggregate_high["signed_event_weight_sum"],
        "low_aero_absolute_support": aggregate_low["absolute_event_weight_support"],
        "high_aero_absolute_support": aggregate_high["absolute_event_weight_support"],
        "low_aero_estimated_proton_yield": aggregate_low["estimated_proton_missing_mass_yield"],
        "high_aero_estimated_proton_yield": aggregate_high["estimated_proton_missing_mass_yield"],
        "low_aero_cleaned_yield": aggregate_low["cleaned_missing_mass_yield"],
        "high_aero_cleaned_yield": aggregate_high["cleaned_missing_mass_yield"],
        "low_aero_average_proton_probability": aggregate_low["average_proton_probability"],
        "high_aero_average_proton_probability": aggregate_high["average_proton_probability"],
        "low_aero_low_mm_removed_fraction": aggregate_low["low_mm_removed_fraction"],
        "high_aero_low_mm_removed_fraction": aggregate_high["low_mm_removed_fraction"],
        "low_aero_lambda_removed_fraction": aggregate_low["lambda_removed_fraction"],
        "high_aero_lambda_removed_fraction": aggregate_high["lambda_removed_fraction"],
        "signed_weight_diagnostics": {
            "positive_event_count": int(signed_positive_events),
            "negative_event_count": int(signed_negative_events),
            "positive_integral": float(signed_positive_integral),
            "negative_integral": float(signed_negative_integral),
            "absolute_negative_integral": float(abs(signed_negative_integral)),
        },
        "warnings": sorted(set(all_warnings)),
        "warnings_by_t_bin": warnings_by_t_bin,
        "high_aero_weight_exceeds_low_aero_by_t_bin": [
            "high_aero_average_weight_exceeds_low_aero" in row["warnings"] for row in warnings_by_t_bin
        ],
        "high_aero_lambda_removal_exceeds_threshold_by_t_bin": [
            ("lambda_region_removal_large_in_high_aero" in row["warnings"]
             or "high_aero_lambda_removal_exceeds_threshold" in row["warnings"])
            for row in warnings_by_t_bin
        ],
    }
    coarse_summary_by_t = [
        _finalize_accumulator(accumulator) for accumulator in coarse_by_t
    ]
    payload["coarse_matrix_per_t_summary"] = coarse_summary_by_t
    integrity = _build_timing_t_diagnostic_integrity(
        cleaning_result,
        matrix_payload,
        coarse_summary_by_t,
        validation_cfg,
        exclusions["source_rows_seen"],
    )
    payload["diagnostic_integrity"] = integrity
    if bool(validation_cfg.get("diagnostic_strict", False)) and not integrity["valid"]:
        raise RuntimeError(
            "timing-t diagnostic integrity failed: {}".format(
                ", ".join(integrity["failures"])
            )
        )
    if include_root_payload:
        root_payload = _build_t_aerogel_root_payload(
            cleaning_result, event_rows, validation_cfg, t_edges, matrix_payload,
        )
        candidate_accounting = (
            (cleaning_result.get("diagnostics") or {}).get("candidate_root_object_accounting") or {}
        )
        payload["root_object_accounting"] = {
            "selected_candidate_fit_root_objects_retained": int(
                candidate_accounting.get("selected_candidate_fit_root_objects_retained", 0) or 0
            ),
            "aerogel_validation_root_objects_retained": int(root_payload.get("object_count", 0)),
            "nonselected_candidate_root_objects_retained": int(
                candidate_accounting.get("rejected_candidate_global_root_objects_retained", 0) or 0
            ),
            "hidden_accumulator_objects_released": int(
                root_payload.get("hidden_accumulator_objects_released", 0) or 0
            ),
            "nonselected_candidate_per_cell_root_objects_released": int(
                candidate_accounting.get("rejected_candidate_per_cell_root_objects_released", 0) or 0
            ),
        }
        payload["signed_weight_diagnostics"].update(
            root_payload.get("signed_histogram_diagnostics") or {}
        )
        payload["_root_payload"] = root_payload
    return payload


def _build_t_prepared_event_weight_lookup(cleaning_result, source_bundle, config):
    t_edges = [float(edge) for edge in (cleaning_result.get("t_edges") or [])]
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    cell_fits = list(cleaning_result.get("delta_t_cell_fits") or [])
    center_shapes = list(cleaning_result.get("delta_timing_center_shapes") or [])
    timing_branch = str((cleaning_result.get("diagnostics") or {}).get("timing_branch") or "")
    denominator_floor = float((config.get("weighting") or {}).get("denominator_floor", 1.0e-12))
    lookup = {}
    closure_by_cell = {}
    closure_by_delta = {}
    closure_by_t = {}
    lookup_by_t_phi = {}
    underflow = overflow = nonfinite = 0
    validation_rows = []
    prepared_sources = _get_prepared_sources(source_bundle)
    canonical_phi_edges = _finite_strictly_increasing_edges(
        ((cleaning_result.get("diagnostics") or {}).get("canonical_t_binning") or {}).get("phi_edges") or []
    )
    aerogel_summary_edges = _resolve_aerogel_summary_edges(
        (config.get("aerogel_validation") or {})
    )
    rf_accept_by_signature = {}
    if bool((cleaning_result.get("diagnostics") or {}).get("rf_applied", False)):
        for source_name, source_spec in prepared_sources.items():
            tree = (((source_bundle or {}).get("sources") or {}).get(source_name) or {}).get("tree")
            entry_map = (source_spec or {}).get("entries") or {}
            if tree is None:
                continue
            for entry_index, evt in enumerate(tree):
                if int(entry_index) not in entry_map:
                    continue
                signature = _make_prepared_event_signature(source_name, entry_index)
                rf_accept_by_signature[signature] = bool(
                    apply_low_epsilon_rf_after_proton_cleaning(cleaning_result, source_name, evt)
                )
    for delta_index, slices in enumerate(cell_fits):
        closure_by_delta[delta_index] = {"delta_index": delta_index, "applied_proton_yield": 0.0, "summed_event_proton_probability": 0.0, "event_count": 0}
        for t_index, fit in enumerate(slices or []):
            applied = float((fit or {}).get("applied_proton_yield", (fit or {}).get("proton_yield", 0.0)) or 0.0)
            closure_by_cell[(delta_index, t_index)] = {
                "delta_index": delta_index,
                "t_index": t_index,
                "t_low": float(t_edges[t_index]),
                "t_high": float(t_edges[t_index + 1]),
                "raw_proton_yield": float((fit or {}).get("raw_proton_yield", 0.0) or 0.0),
                "fitted_proton_yield": float((fit or {}).get("fitted_proton_yield", (fit or {}).get("proton_yield", 0.0)) or 0.0),
                "applied_proton_yield": applied,
                "applied_cell_enabled": bool((fit or {}).get("applied_cell_enabled", False)),
                "applied_zero_reason": (fit or {}).get("applied_zero_reason"),
                "summed_event_proton_probability": 0.0,
                "event_count": 0,
            }
            closure_by_delta[delta_index]["applied_proton_yield"] += applied
            closure_by_t.setdefault(t_index, {"t_index": t_index, "t_low": float(t_edges[t_index]), "t_high": float(t_edges[t_index + 1]), "applied_proton_yield": 0.0, "summed_event_proton_probability": 0.0, "event_count": 0})["applied_proton_yield"] += applied
    for source_name, source_spec in prepared_sources.items():
        fit_coefficient = float((source_spec or {}).get("fit_coefficient", 0.0) or 0.0)
        physical_coefficient = float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        for entry_index, entry_payload in ((source_spec or {}).get("entries") or {}).items():
            t_value = float((entry_payload or {}).get("adj_t", float("nan")))
            delta_value = float((entry_payload or {}).get("delta_value", float("nan")))
            support_label = SUPPORT_UNSUPPORTED
            proton_weight = 0.0
            t_index = -1
            lookup_zero_reason = None
            if not math.isfinite(t_value):
                nonfinite += 1
                lookup_zero_reason = "nonfinite_t"
            else:
                t_index = _find_collection_bin(t_value, t_edges)
                if t_index < 0:
                    if t_edges and t_value > float(t_edges[-1]):
                        overflow += 1
                        lookup_zero_reason = "t_overflow"
                    else:
                        underflow += 1
                        lookup_zero_reason = "t_underflow"
                elif t_index >= len(t_edges) - 1:
                    overflow += 1
                    lookup_zero_reason = "t_overflow"
            delta_index = _find_collection_bin(delta_value, delta_edges)
            selected_timing_value = ((entry_payload or {}).get("timing_values") or {}).get(
                timing_branch
            )
            if 0 <= t_index < len(t_edges) - 1 and 0 <= delta_index < len(cell_fits):
                fit = (cell_fits[delta_index] or [])[t_index] if t_index < len(cell_fits[delta_index] or []) else None
                shape = center_shapes[delta_index] if delta_index < len(center_shapes) else None
                timing_value = selected_timing_value
                support_label = str((fit or {}).get("support_label", SUPPORT_UNSUPPORTED))
                if (
                    fit and shape and timing_value is not None
                    and str((fit or {}).get("support_label", SUPPORT_UNSUPPORTED))
                    in (SUPPORT_SUPPORTED, SUPPORT_MARGINAL)
                    and bool((fit or {}).get("valid", False))
                    and bool((fit or {}).get("applied_cell_enabled", False))
                ):
                    proton_weight = _evaluate_event_proton_probability(
                        float(timing_value), shape, fit, denominator_floor
                    )
                elif timing_value is None:
                    lookup_zero_reason = "missing_timing_value"
                elif fit:
                    lookup_zero_reason = str(
                        (fit or {}).get("applied_zero_reason") or "unsupported_delta"
                    )
            proton_weight = max(0.0, min(1.0, float(proton_weight)))
            cleaned_factor = 1.0 - proton_weight
            signature = _make_prepared_event_signature(source_name, entry_index)
            rf_accept = bool(rf_accept_by_signature.get(signature, True))
            lookup[signature] = {
                "source_label": str(source_name), "source_entry_index": int(entry_index),
                "delta_index": int(delta_index), "t_index": int(t_index),
                "t_low": float(t_edges[t_index]) if 0 <= t_index < len(t_edges) - 1 else None,
                "t_high": float(t_edges[t_index + 1]) if 0 <= t_index < len(t_edges) - 1 else None,
                "support_label": support_label,
                # These values are the immutable timing-model proposal.  The
                # Lambda gate below is the only owner allowed to commit the
                # production aliases used downstream.
                "proposed_proton_probability": proton_weight,
                "proposed_cleaned_factor": cleaned_factor,
                "proposed_final_cleaned_factor": cleaned_factor if rf_accept else 0.0,
                "event_lookup_zero_reason": lookup_zero_reason,
                "rf_accept": rf_accept,
            }
            if bool((entry_payload or {}).get("nommcuts", False)):
                if (delta_index, t_index) in closure_by_cell:
                    weighted_probability = fit_coefficient * proton_weight
                    cell = closure_by_cell[(delta_index, t_index)]
                    cell["summed_event_proton_probability"] += weighted_probability
                    cell["event_count"] += 1
                    closure_by_delta[delta_index]["summed_event_proton_probability"] += weighted_probability
                    closure_by_delta[delta_index]["event_count"] += 1
                    closure_by_t[t_index]["summed_event_proton_probability"] += weighted_probability
                    closure_by_t[t_index]["event_count"] += 1
                aero_value = float((entry_payload or {}).get("aero_value", float("nan")))
                aero_index = _find_collection_bin(aero_value, aerogel_summary_edges)
                validation_rows.append({
                    "t_index": int(t_index), "aero_index": int(aero_index),
                    "source_label": str(source_name), "source_entry_index": int(entry_index),
                    "is_prompt_source": bool((source_spec or {}).get("is_prompt_source", False)),
                    "physical_weight": physical_coefficient,
                    "proposed_proton_probability": proton_weight,
                    "proposed_cleaned_factor": cleaned_factor,
                    "proposed_final_cleaned_factor": cleaned_factor if rf_accept else 0.0,
                    "rf_accept": rf_accept,
                    "adj_t": t_value,
                    "phi_value": float((entry_payload or {}).get("phi_value", float("nan"))),
                    "aero_value": aero_value,
                    "selected_timing": _finite_float_or_none(selected_timing_value),
                    "support_label": support_label,
                    "event_lookup_zero_reason": lookup_zero_reason,
                    "adj_mm": float((entry_payload or {}).get("adj_mm", float("nan"))),
                })
                phi_value = _finite_float_or_none((entry_payload or {}).get("phi_value"))
                phi_index = _find_collection_bin(phi_value, canonical_phi_edges) if phi_value is not None else -1
                if 0 <= phi_index < len(canonical_phi_edges) - 1 and 0 <= t_index < len(t_edges) - 1:
                    key = (int(t_index), int(phi_index))
                    row = lookup_by_t_phi.setdefault(key, {
                        "t_index": int(t_index), "t_low": float(t_edges[t_index]), "t_high": float(t_edges[t_index + 1]),
                        "phi_index": int(phi_index), "phi_low": float(canonical_phi_edges[phi_index]),
                        "phi_high": float(canonical_phi_edges[phi_index + 1]), "event_count": 0,
                        "signed_physical_yield": 0.0, "estimated_proton_yield": 0.0,
                        "frozen_cleaned_yield": 0.0,
                    })
                    row["event_count"] += 1
                    row["signed_physical_yield"] += physical_coefficient
                    row["estimated_proton_yield"] += physical_coefficient * proton_weight
                    row["frozen_cleaned_yield"] += physical_coefficient * cleaned_factor
    # The timing fit and all cell/delta closure terms above are deliberately
    # proposal-only.  Evaluate the physics-preservation gate only after that
    # immutable proposal exists, then commit one setting-wide applied lookup.
    lambda_gate = _evaluate_lambda_preservation_gate(
        cleaning_result, validation_rows, config
    )
    _commit_lambda_preservation_gate(lookup, validation_rows, lambda_gate)
    proposed_validation_rows = []
    for row in validation_rows:
        proposal_row = dict(row)
        proposed_probability, proposed_cleaned, proposed_final = _lambda_gate_stage_values(
            proposal_row, "proposed"
        )
        proposal_row.update(
            {
                "proton_weight": proposed_probability,
                "cleaned_factor": proposed_cleaned,
                "final_cleaned_factor": proposed_final,
            }
        )
        proposed_validation_rows.append(proposal_row)

    for collection in (closure_by_cell, closure_by_delta, closure_by_t):
        for row in collection.values():
            applied = float(row.get("applied_proton_yield", 0.0) or 0.0)
            row["closure_scope"] = "proposed_timing_t_model"
            if applied > 0.0:
                row["closure_ratio"] = float(row["summed_event_proton_probability"] / applied)
                row["closure_status"] = "evaluated_applied_yield"
            else:
                row["closure_ratio"] = None
                row["closure_status"] = "not_applicable_zero_applied_yield"
    for delta_index, slices in enumerate(cell_fits):
        for t_index, fit in enumerate(slices or []):
            closure = closure_by_cell.get((delta_index, t_index))
            if closure is not None:
                fit["closure_status"] = closure.get("closure_status")
                fit["closure_ratio"] = closure.get("closure_ratio")
    diagnostics = cleaning_result.setdefault("diagnostics", {})
    diagnostics["event_weight_source"] = "setting_wide_immutable_prepared_lookup"
    diagnostics["event_weight_closure_label"] = "proposed timing-|t| model closure"
    diagnostics["lambda_preservation_gate"] = _json_ready_value(lambda_gate)
    diagnostics["lambda_preservation_event_audit"] = _json_ready_value(
        _lambda_gate_audit_rows(validation_rows)
    )
    diagnostics["event_weight_closure_by_cell"] = _json_ready_value(list(closure_by_cell.values()))
    diagnostics["event_weight_closure_by_delta"] = _json_ready_value(list(closure_by_delta.values()))
    diagnostics["event_weight_closure_by_t"] = _json_ready_value(list(closure_by_t.values()))
    diagnostics["event_weight_closure_by_delta_t"] = _json_ready_value(list(closure_by_cell.values()))
    diagnostics["event_weight_lookup_by_t_phi"] = _json_ready_value(list(lookup_by_t_phi.values()))
    diagnostics["applied_timing_t_cell_map"] = _json_ready_value(
        [dict(fit or {}) for slices in cell_fits for fit in (slices or [])]
    )
    diagnostics["event_probability_sum_by_delta_t"] = _json_ready_value([
        [float(closure_by_cell[(delta_index, t_index)]["summed_event_proton_probability"]) for t_index in range(len(slices or []))]
        for delta_index, slices in enumerate(cell_fits)
    ])
    diagnostics["event_probability_count_by_delta_t"] = _json_ready_value([
        [int(closure_by_cell[(delta_index, t_index)]["event_count"]) for t_index in range(len(slices or []))]
        for delta_index, slices in enumerate(cell_fits)
    ])
    diagnostics["t_lookup_boundary_counts"] = {"underflow": underflow, "overflow": overflow, "nonfinite": nonfinite}
    aerogel_validation = _build_t_aerogel_validation(
        cleaning_result, source_bundle, proposed_validation_rows, config, include_root_payload=True,
    )
    cleaning_result["_aerogel_vs_t_root_payload"] = aerogel_validation.pop("_root_payload", {})
    diagnostics["aerogel_vs_t_validation"] = _json_ready_value(aerogel_validation)
    # Keep the pre-existing artifact/consumer name as a compatibility alias.
    diagnostics["aerogel_validation"] = diagnostics["aerogel_vs_t_validation"]
    diagnostics["timing_t_diagnostic_integrity"] = _json_ready_value(
        aerogel_validation.get("diagnostic_integrity") or {}
    )
    # Retain the timing-model proposal beside the committed production state.
    # The former is needed to understand a rejected model; only the latter is
    # permitted to populate the downstream kaon and pion targets.
    proposed_mm_diagnostics = _build_timing_t_mm_diagnostics(
        cleaning_result, proposed_validation_rows, config, include_root_payload=True,
    )
    applied_mm_diagnostics = _build_timing_t_mm_diagnostics(
        cleaning_result, validation_rows, config, include_root_payload=True,
    )
    proposed_root_payload = proposed_mm_diagnostics.pop("_root_payload", {})
    applied_root_payload = applied_mm_diagnostics.pop("_root_payload", {})
    proposed_public = dict(proposed_mm_diagnostics)
    applied_public = dict(applied_mm_diagnostics)
    mm_diagnostics = dict(proposed_public)
    mm_diagnostics.update(
        {
            "stage": "proposed",
            "gate_status": lambda_gate.get("status"),
            "proposed": proposed_public,
            "applied": applied_public,
        }
    )
    proposed_maps = proposed_root_payload.get("maps") or {}
    applied_maps = applied_root_payload.get("maps") or {}
    mm_root_maps = {
        "raw": proposed_maps.get("raw"),
        "proposed_estimated_proton": proposed_maps.get("estimated_proton"),
        "proposed_cleaned_pre_rf": proposed_maps.get("cleaned_pre_rf"),
        "applied_cleaned_pre_rf": applied_maps.get("cleaned_pre_rf"),
        "applied_cleaned_final_rf": applied_maps.get("cleaned_final_rf"),
    }
    cleaning_result["_timing_t_mm_root_payload"] = {
        "maps": mm_root_maps,
        "window_removed_fraction": proposed_root_payload.get("window_removed_fraction") or {},
        "object_count": int(
            sum(hist is not None for hist in mm_root_maps.values())
            + len(proposed_root_payload.get("window_removed_fraction") or {})
        ),
    }
    diagnostics["timing_t_mm_diagnostics"] = _json_ready_value(mm_diagnostics)
    diagnostics["timing_t_summary"] = _json_ready_value(
        _build_timing_t_summary(
            cleaning_result,
            aerogel_validation.get("matrix_payload") or {},
            aerogel_validation.get("aggregate_all_frozen_rows") or {},
        )
    )
    diagnostics["prepared_event_lookup_count"] = int(len(lookup))
    return lookup


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
    event_probability_sum_by_delta_aero = []
    event_probability_count_by_delta_aero = []
    for delta_index, slice_collection in enumerate(delta_slice_fits):
        probability_row = []
        count_row = []
        for aero_index, _ in enumerate(slice_collection or []):
            closure_row = closure_by_cell.get((int(delta_index), int(aero_index))) or {}
            probability_row.append(
                float(closure_row.get("summed_event_proton_probability", 0.0) or 0.0)
            )
            count_row.append(int(closure_row.get("event_count", 0) or 0))
        event_probability_sum_by_delta_aero.append(probability_row)
        event_probability_count_by_delta_aero.append(count_row)
    diagnostics["event_weight_closure_by_cell"] = _json_ready_value(list(closure_by_cell.values()))
    diagnostics["event_weight_closure_by_delta"] = _json_ready_value(list(closure_by_delta.values()))
    diagnostics["event_probability_sum_by_delta_aero"] = _json_ready_value(event_probability_sum_by_delta_aero)
    diagnostics["event_probability_count_by_delta_aero"] = _json_ready_value(event_probability_count_by_delta_aero)
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


def _build_signed_timing_t_histograms(
    source_bundle,
    config,
    t_edges,
    timing_branch,
    *,
    time_hist_range,
    time_hist_bins,
):
    """Build the production timing histograms with canonical t as the cell axis."""
    t_edges = [float(edge) for edge in (t_edges if t_edges is not None else [])]
    if len(t_edges) < 2 or not np.all(np.isfinite(t_edges)) or np.any(np.diff(t_edges) <= 0.0):
        raise ValueError("timing-t cleaning requires finite, strictly increasing canonical t edges")
    time_min, time_max = [float(value) for value in time_hist_range]
    n_time_bins = int(time_hist_bins)
    if not (math.isfinite(time_min) and math.isfinite(time_max) and time_max > time_min and n_time_bins > 0):
        raise ValueError("timing-t candidate requires a finite increasing timing range and positive bin count")
    delta_min, delta_max = [float(value) for value in (config.get("delta_hist_range") or (-10.0, 20.0))]
    delta_bins = int(config.get("delta_bins", 10) or 10)
    phi_token = str((source_bundle or {}).get("phi_setting") or "setting").lower()
    eps_token = str((source_bundle or {}).get("epsset") or "eps").lower()
    prefix = "H_proton_cleaning_t_{}_{}_{}".format(phi_token, eps_token, timing_branch)
    t_axis = array("d", t_edges)

    h_global_timing = ROOT.TH1D(
        "{}_global".format(prefix), "Global {};{} [ns];Signed yield".format(timing_branch, timing_branch),
        n_time_bins, time_min, time_max,
    )
    h_global_timing.SetDirectory(0)
    h_global_timing.Sumw2()
    h_global_timing_vs_t = ROOT.TH2D(
        "{}_vs_t".format(prefix), "Global {} versus |t|;|t| [GeV^{{2}}];{} [ns];Signed yield".format(timing_branch, timing_branch),
        len(t_edges) - 1, t_axis, n_time_bins, time_min, time_max,
    )
    h_global_timing_vs_t.SetDirectory(0)
    h_global_timing_vs_t.Sumw2()
    h_delta_timing = ROOT.TH2D(
        "{}_vs_delta".format(prefix), "{} versus SHMS delta;SHMS delta;{} [ns];Signed yield".format(timing_branch, timing_branch),
        delta_bins, delta_min, delta_max, n_time_bins, time_min, time_max,
    )
    h_delta_timing.SetDirectory(0)
    h_delta_timing.Sumw2()

    global_t_slices = []
    delta_t_cells = []
    delta_all_t_hists = []
    for t_index in range(len(t_edges) - 1):
        hist = ROOT.TH1D(
            "{}_global_t_{}".format(prefix, t_index),
            "Global |t| slice {};{} [ns];Signed yield".format(t_index + 1, timing_branch),
            n_time_bins, time_min, time_max,
        )
        hist.SetDirectory(0)
        hist.Sumw2()
        global_t_slices.append(hist)
    for delta_index in range(delta_bins):
        all_t_hist = ROOT.TH1D(
            "{}_delta_{}_all_t".format(prefix, delta_index),
            "Delta {} all t;{} [ns];Signed yield".format(delta_index + 1, timing_branch),
            n_time_bins, time_min, time_max,
        )
        all_t_hist.SetDirectory(0)
        all_t_hist.Sumw2()
        delta_all_t_hists.append(all_t_hist)
        cells = []
        for t_index in range(len(t_edges) - 1):
            hist = ROOT.TH1D(
                "{}_delta_{}_t_{}".format(prefix, delta_index, t_index),
                "Timing delta {} |t| slice {};{} [ns];Signed yield".format(
                    delta_index + 1, t_index + 1, timing_branch
                ),
                n_time_bins, time_min, time_max,
            )
            hist.SetDirectory(0)
            hist.Sumw2()
            cells.append(hist)
        delta_t_cells.append(cells)

    delta_edges = np.linspace(delta_min, delta_max, delta_bins + 1)
    global_prompt_support = [0 for _ in range(len(t_edges) - 1)]
    cell_prompt_support = [[0 for _ in range(len(t_edges) - 1)] for _ in range(delta_bins)]
    source_stats = {}
    for source_name, source_spec in _get_prepared_sources(source_bundle).items():
        fit_coefficient = float((source_spec or {}).get("fit_coefficient", 0.0) or 0.0)
        stats = {
            "tree_name": (source_spec or {}).get("tree_name"),
            "timing_branch": str(timing_branch),
            "entries_passing_selection": 0,
            "entries_missing_timing_branch": 0,
            "entries_outside_timing_range": 0,
            "entries_outside_t_range": 0,
            "entries_outside_delta_range": 0,
            "entries_used": 0,
        }
        source_stats[str(source_name)] = stats
        if fit_coefficient == 0.0:
            continue
        for entry_payload in ((source_spec or {}).get("entries") or {}).values():
            if not bool((entry_payload or {}).get("nommcuts", False)):
                continue
            stats["entries_passing_selection"] += 1
            timing_value = ((entry_payload or {}).get("timing_values") or {}).get(str(timing_branch))
            if timing_value is None:
                stats["entries_missing_timing_branch"] += 1
                continue
            timing_value = float(timing_value)
            if not (time_min <= timing_value <= time_max):
                stats["entries_outside_timing_range"] += 1
                continue
            t_index = _find_collection_bin(float((entry_payload or {}).get("adj_t", float("nan"))), t_edges)
            if not (0 <= t_index < len(global_t_slices)):
                stats["entries_outside_t_range"] += 1
                continue
            delta_index = _find_collection_bin(
                float((entry_payload or {}).get("delta_value", float("nan"))), delta_edges
            )
            if not (0 <= delta_index < delta_bins):
                stats["entries_outside_delta_range"] += 1
                continue
            stats["entries_used"] += 1
            h_global_timing.Fill(timing_value, fit_coefficient)
            h_global_timing_vs_t.Fill(float(entry_payload["adj_t"]), timing_value, fit_coefficient)
            h_delta_timing.Fill(float(entry_payload["delta_value"]), timing_value, fit_coefficient)
            global_t_slices[t_index].Fill(timing_value, fit_coefficient)
            delta_all_t_hists[delta_index].Fill(timing_value, fit_coefficient)
            delta_t_cells[delta_index][t_index].Fill(timing_value, fit_coefficient)
            if bool((source_spec or {}).get("is_prompt_source", False)):
                global_prompt_support[t_index] += 1
                cell_prompt_support[delta_index][t_index] += 1
    return {
        "H_global_timing": h_global_timing,
        "H_global_timing_vs_t": h_global_timing_vs_t,
        "H_delta_timing": h_delta_timing,
        "global_t_slices": global_t_slices,
        "delta_all_t_hists": delta_all_t_hists,
        "delta_t_cells": delta_t_cells,
        "global_prompt_support": global_prompt_support,
        "cell_prompt_support": cell_prompt_support,
        "delta_edges": [float(edge) for edge in delta_edges],
        "t_edges": t_edges,
        "source_stats": source_stats,
        "time_hist_range": (time_min, time_max),
        "time_hist_bins": n_time_bins,
        "timing_branch": str(timing_branch),
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


def resolve_timing_t_candidate_configuration(
    source_bundle,
    config,
    timing_branch,
    evaluate_event,
    hole_contains,
    mm_min,
    mm_max,
):
    """Resolve one independent RF or CT timing-t candidate.

    This is deliberately shared with the established probe helpers: CT keeps
    its fixed exact ranges and RF discovers its display range from selected
    prompt records.  The returned config is private to this candidate.
    """
    timing_branch = str(timing_branch)
    probe_kind = "ct" if timing_branch == PROTON_CLEANING_EXACT_TIMING_BRANCH else "rf"
    candidate_config = deepcopy(config)
    beam_bunch_spacing_ns = _resolve_beam_bunch_spacing_ns(source_bundle)
    global_cfg = dict(candidate_config.get("global_fit") or {})
    global_cfg["beam_bunch_spacing_ns"] = float(beam_bunch_spacing_ns)
    candidate_config["global_fit"] = global_cfg
    range_diagnostics = {}
    if probe_kind == "ct":
        ct_config = _resolve_ct_probe_configuration(source_bundle)
        display_range = tuple(float(value) for value in ct_config["display_range"])
        fit_range = tuple(float(value) for value in ct_config["fit_range"])
        candidate_config["global_fit"] = {
            **global_cfg,
            **dict(ct_config.get("global_fit") or {}),
            "beam_bunch_spacing_ns": float(beam_bunch_spacing_ns),
        }
        candidate_config["probe_fixed_bounds"] = {
            key: ct_config[key]
            for key in (
                "kaonMeanMin", "kaonMeanMax", "protonMeanMin", "protonMeanMax",
                "sigmaMax", "sigmaInitial",
            )
        }
        histogram_bins = int(ct_config["histogram_bins"])
        proton_peak_is_lower = False
        range_source = "resolve_ct_probe_configuration"
        range_diagnostics = {"ct_probe_configuration": _json_ready_value(ct_config)}
    else:
        display_range, range_diagnostics = _resolve_rf_probe_display_range(
            source_bundle,
            evaluate_event,
            hole_contains,
            mm_min,
            mm_max,
            timing_branch,
            selection_key="nommcuts",
            source_names=("prompt",),
            return_diagnostics=True,
        )
        display_range = tuple(float(value) for value in display_range)
        fit_range = display_range
        histogram_bins = int(
            _resolve_probe_time_histogram_bins(source_bundle, probe_kind, display_range)
        )
        proton_peak_is_lower = True
        range_source = "resolve_rf_probe_display_range"
    return {
        "timing_branch": timing_branch,
        "probe_kind": probe_kind,
        "candidate_config": candidate_config,
        "display_range": display_range,
        "fit_range": fit_range,
        "histogram_bins": histogram_bins,
        "proton_peak_is_lower": proton_peak_is_lower,
        "range_source": range_source,
        "range_diagnostics": _json_ready_value(range_diagnostics),
    }


def _prepared_selection_has_entries(source_bundle, selection_key):
    for _, _, _, _ in _iter_prepared_records(source_bundle, selection_key=selection_key):
        return True
    return False


def _timing_center_payload(shape, source, **provenance):
    """Create a timing-center-only object safe to reuse across delta bins."""
    shape = dict(shape or {})
    fields = ("kaon_mean", "proton_mean", "kaon_sigma", "proton_sigma")
    resolved = {}
    valid = bool(shape.get("valid", False))
    for field in fields:
        try:
            value = float(shape[field])
        except (KeyError, TypeError, ValueError):
            valid = False
            break
        if not math.isfinite(value):
            valid = False
            break
        resolved[field] = value
    payload = {
        "valid": bool(valid),
        "timing_center_model_valid": bool(valid),
        "timing_center_source": str(source),
        "fit_min": shape.get("fit_min"),
        "fit_max": shape.get("fit_max"),
        "kaon_mean": resolved.get("kaon_mean"),
        "proton_mean": resolved.get("proton_mean"),
        "kaon_sigma": resolved.get("kaon_sigma"),
        "proton_sigma": resolved.get("proton_sigma"),
        "source_parameters": dict(provenance.pop("source_parameters", resolved)),
        "resolved_parameters": dict(resolved),
        "lower_neighbor_index": provenance.pop("lower_neighbor_index", None),
        "upper_neighbor_index": provenance.pop("upper_neighbor_index", None),
        "nearest_neighbor_index": provenance.pop("nearest_neighbor_index", None),
        "target_delta_center": provenance.pop("target_delta_center", None),
        "lower_delta_center": provenance.pop("lower_delta_center", None),
        "upper_delta_center": provenance.pop("upper_delta_center", None),
        "interpolation_fraction": provenance.pop("interpolation_fraction", None),
        "rejection_reasons": list(provenance.pop("rejection_reasons", [])),
    }
    payload.update(provenance)
    if not payload["valid"] and not payload["rejection_reasons"]:
        payload["rejection_reasons"] = ["invalid_timing_center"]
    return payload


def _interpolate_delta_timing_center_shape(left_shape, right_shape, left_center, target_center, right_center):
    """Interpolate only timing locations and widths between valid delta fits."""
    if not (bool((left_shape or {}).get("valid", False)) and bool((right_shape or {}).get("valid", False))):
        return None
    left_center = float(left_center)
    right_center = float(right_center)
    target_center = float(target_center)
    if not (math.isfinite(left_center) and math.isfinite(right_center) and math.isfinite(target_center)):
        return None
    if right_center <= left_center:
        return None
    fraction = (target_center - left_center) / (right_center - left_center)
    resolved = {}
    for field in ("kaon_mean", "proton_mean", "kaon_sigma", "proton_sigma"):
        try:
            left_value = float(left_shape[field])
            right_value = float(right_shape[field])
        except (KeyError, TypeError, ValueError):
            return None
        if not (math.isfinite(left_value) and math.isfinite(right_value)):
            return None
        resolved[field] = float(left_value + fraction * (right_value - left_value))
    return _timing_center_payload(
        {
            "valid": True,
            "fit_min": left_shape.get("fit_min"),
            "fit_max": left_shape.get("fit_max"),
            **resolved,
        },
        "neighbor_delta_interpolation",
        lower_delta_center=left_center,
        target_delta_center=target_center,
        upper_delta_center=right_center,
        interpolation_fraction=float(fraction),
        source_parameters={
            "lower": {field: left_shape.get(field) for field in resolved},
            "upper": {field: right_shape.get(field) for field in resolved},
        },
    )


def _nearest_delta_timing_center_shape(neighbor_shape, neighbor_index, neighbor_center, target_center):
    return _timing_center_payload(
        neighbor_shape,
        "neighbor_delta_nearest_fallback",
        nearest_neighbor_index=int(neighbor_index),
        target_delta_center=float(target_center),
        lower_delta_center=float(neighbor_center),
    )


def _classify_timing_t_support(valid_cells, coverage, modeled_yield, thresholds):
    """Use the legacy supported/marginal hierarchy for canonical-t cells."""
    thresholds = dict(thresholds or {})
    valid_cells = int(valid_cells)
    coverage = float(coverage)
    modeled_yield = float(modeled_yield)
    if (
        valid_cells >= int(thresholds.get("minimum_supported_t_slices", thresholds.get("minimum_supported_slices", 2)))
        and coverage >= float(thresholds.get("minimum_supported_coverage", 0.35))
        and modeled_yield >= float(thresholds.get("minimum_modeled_yield", 5.0))
    ):
        return SUPPORT_SUPPORTED
    if (
        valid_cells >= int(thresholds.get("minimum_marginal_t_slices", thresholds.get("minimum_marginal_slices", 1)))
        and coverage >= float(thresholds.get("minimum_marginal_coverage", 0.15))
        and modeled_yield >= float(thresholds.get("minimum_modeled_yield", 5.0))
    ):
        return SUPPORT_MARGINAL
    return SUPPORT_UNSUPPORTED


def _classify_timing_t_setting_support(delta_rows, thresholds):
    """Aggregate delta rows without allowing one isolated cell to accept a setting."""
    thresholds = dict(thresholds or {})
    rows = list(delta_rows or [])
    valid_cells = int(sum(int(row.get("valid_t_cells", 0) or 0) for row in rows))
    data_total = float(sum(float(row.get("data_total", 0.0) or 0.0) for row in rows))
    fitted_data_total = float(sum(float(row.get("fitted_data_total", 0.0) or 0.0) for row in rows))
    modeled_yield = float(sum(float(row.get("model_total", 0.0) or 0.0) for row in rows))
    coverage = float(fitted_data_total / data_total) if data_total > 0.0 else 0.0
    min_setting_cells = int(thresholds.get("minimum_setting_valid_t_cells", 2))
    supported_deltas = int(sum(row.get("support_label") == SUPPORT_SUPPORTED for row in rows))
    marginal_deltas = int(sum(row.get("support_label") == SUPPORT_MARGINAL for row in rows))
    if valid_cells < min_setting_cells:
        label = SUPPORT_UNSUPPORTED
    elif (
        supported_deltas > 0
        and coverage >= float(thresholds.get("minimum_supported_coverage", 0.35))
        and modeled_yield >= float(thresholds.get("minimum_modeled_yield", 5.0))
    ):
        label = SUPPORT_SUPPORTED
    elif (
        (supported_deltas + marginal_deltas) > 0
        and coverage >= float(thresholds.get("minimum_marginal_coverage", 0.15))
        and modeled_yield >= float(thresholds.get("minimum_modeled_yield", 5.0))
    ):
        label = SUPPORT_MARGINAL
    else:
        label = SUPPORT_UNSUPPORTED
    return {
        "support_label": label,
        "accepted": bool(label in (SUPPORT_SUPPORTED, SUPPORT_MARGINAL)),
        "valid_t_cells": valid_cells,
        "data_total": data_total,
        "fitted_data_total": fitted_data_total,
        "model_total": modeled_yield,
        "coverage": coverage,
        "supported_delta_bins": supported_deltas,
        "marginal_delta_bins": marginal_deltas,
        "minimum_setting_valid_t_cells": min_setting_cells,
    }


def _apply_timing_t_cell_map(delta_t_cell_fits, support_by_delta, setting_support):
    """Freeze fit/support decisions into the only cell map lookup may use."""
    setting_label = str(setting_support.get("support_label", SUPPORT_UNSUPPORTED))
    for delta_index, slices in enumerate(delta_t_cell_fits or []):
        delta_label = str(
            support_by_delta[delta_index]
            if delta_index < len(support_by_delta)
            else SUPPORT_UNSUPPORTED
        )
        for fit in slices or []:
            raw_yield = float((fit or {}).get("raw_proton_yield", 0.0) or 0.0)
            fitted_yield = float((fit or {}).get("proton_yield", 0.0) or 0.0)
            fit_valid = bool((fit or {}).get("cell_fit_valid", (fit or {}).get("valid", False)))
            component_detected = bool((fit or {}).get("proton_component_detected", False))
            zero_reason = None
            if not fit_valid:
                zero_reason = "invalid_cell_fit"
            elif not component_detected:
                zero_reason = "weak_proton_component"
            elif delta_label == SUPPORT_UNSUPPORTED:
                zero_reason = "unsupported_delta"
            elif setting_label == SUPPORT_UNSUPPORTED:
                zero_reason = "setting_gate_rejected"
            enabled = zero_reason is None
            fit.update({
                "raw_proton_yield": raw_yield,
                "fitted_proton_yield": fitted_yield,
                "delta_support_label": delta_label,
                "setting_support_label": setting_label,
                "applied_cell_enabled": bool(enabled),
                "applied_proton_yield": fitted_yield if enabled else 0.0,
                "applied_zero_reason": zero_reason,
                "cell_fit_valid": fit_valid,
                "proton_component_detected": component_detected,
            })


def _timing_t_candidate_selection_score(setting_support, global_shape, branch_preference):
    """Comparable RF/CT ranking: support first, normalized fit metrics last."""
    support_rank = SUPPORT_CLASS_TO_CODE.get(
        str(setting_support.get("support_label", SUPPORT_UNSUPPORTED)), 0.0
    )
    significance = _finite_float_or_none(
        (global_shape or {}).get("proton_component_significance")
    )
    deviance = _finite_float_or_none(
        (global_shape or {}).get("poisson_deviance_per_entry")
    )
    return (
        int(bool(setting_support.get("accepted", False))),
        int(support_rank),
        int(setting_support.get("supported_delta_bins", 0) or 0),
        int(setting_support.get("marginal_delta_bins", 0) or 0),
        int(setting_support.get("valid_t_cells", 0) or 0),
        float(setting_support.get("coverage", 0.0) or 0.0),
        int(bool((global_shape or {}).get("valid", False))),
        float(significance if significance is not None else -1.0e300),
        -float(deviance if deviance is not None else 1.0e300),
        int(branch_preference),
    )


def _rank_timing_t_candidate_evaluations(candidate_evaluations):
    """Return production-ranked candidates; scores are already axis-normalized."""
    ranked = sorted(
        list(candidate_evaluations or []),
        key=lambda item: tuple(item.get("candidate_selection_score") or ()),
        reverse=True,
    )
    return ranked


def _evaluate_timing_t_candidate_production(
    candidate,
    timing_payload,
    stable_global_shape,
    t_edges,
    phi_setting,
):
    """Evaluate one RF/CT candidate through the full delta x canonical-t model."""
    selected_config = candidate["candidate_config"]
    timing_branch = str(candidate["timing_branch"])
    peak_is_lower = bool(candidate["proton_peak_is_lower"])
    fit_range = candidate["fit_range"]
    rejection_reasons = []
    if not bool((stable_global_shape or {}).get("valid", False)):
        rejection_reasons.append("global_timing_fit_invalid")
        return {
            "candidate": candidate,
            "timing_payload": timing_payload,
            "stable_global_shape": stable_global_shape,
            "global_t_shapes": [],
            "delta_timing_center_shapes": [],
            "timing_center_source_by_delta": [],
            "delta_t_cell_fits": [],
            "support_by_delta": [],
            "support_by_delta_t": [],
            "delta_support": [],
            "setting_support": {
                "support_label": SUPPORT_UNSUPPORTED,
                "accepted": False,
                "valid_t_cells": 0,
                "coverage": 0.0,
                "supported_delta_bins": 0,
                "marginal_delta_bins": 0,
                "model_total": 0.0,
            },
            "rejection_reasons": rejection_reasons,
        }

    global_t_shapes = []
    for t_index, hist in enumerate(timing_payload["global_t_slices"]):
        global_t_shapes.append(
            _fit_global_timing_shape(
                hist,
                selected_config,
                "F_proton_cleaning_t_global_{}_{}_{}".format(phi_setting, timing_branch, t_index),
                proton_peak_is_lower=peak_is_lower,
                display_range=fit_range,
                fit_mode="local_peak_rescue",
            )
        )

    delta_center_shapes = []
    center_sources = []
    for delta_index, hist in enumerate(timing_payload["delta_all_t_hists"]):
        raw_shape = _fit_global_timing_shape(
            hist,
            selected_config,
            "F_proton_cleaning_t_delta_all_{}_{}".format(phi_setting, delta_index),
            proton_peak_is_lower=peak_is_lower,
            display_range=fit_range,
            fit_mode="local_peak_rescue",
        )
        source = "delta_all_t_fit"
        if not bool(raw_shape.get("valid", False)):
            source = "invalid_timing_center"
            supported_hist = _clone_hist(hist, "{}_supported".format(hist.GetName()), reset=True)
            for t_index, cell_hist in enumerate(timing_payload["delta_t_cells"][delta_index]):
                if t_index < len(global_t_shapes) and bool(global_t_shapes[t_index].get("valid", False)):
                    supported_hist.Add(cell_hist)
            supported_shape = _fit_global_timing_shape(
                supported_hist,
                selected_config,
                "F_proton_cleaning_t_delta_supported_{}_{}".format(phi_setting, delta_index),
                proton_peak_is_lower=peak_is_lower,
                display_range=fit_range,
                fit_mode="local_peak_rescue",
            )
            if bool(supported_shape.get("valid", False)):
                raw_shape = supported_shape
                source = "delta_supported_t_fit"
        delta_center_shapes.append(_timing_center_payload(raw_shape, source))
        center_sources.append(source)

    direct_center_valid = [bool(shape.get("valid", False)) for shape in delta_center_shapes]
    delta_centers = [
        0.5 * (float(timing_payload["delta_edges"][index]) + float(timing_payload["delta_edges"][index + 1]))
        for index in range(len(timing_payload["delta_edges"]) - 1)
    ]
    stable_center = _timing_center_payload(stable_global_shape, "stable_global_center_fallback")
    for delta_index, center_shape in enumerate(delta_center_shapes):
        if bool(center_shape.get("valid", False)):
            continue
        left_index = delta_index - 1
        right_index = delta_index + 1
        left_valid = 0 <= left_index < len(delta_center_shapes) and direct_center_valid[left_index]
        right_valid = 0 <= right_index < len(delta_center_shapes) and direct_center_valid[right_index]
        if left_valid and right_valid:
            interpolated = _interpolate_delta_timing_center_shape(
                delta_center_shapes[left_index], delta_center_shapes[right_index],
                delta_centers[left_index], delta_centers[delta_index], delta_centers[right_index],
            )
            if interpolated is not None:
                interpolated["lower_neighbor_index"] = int(left_index)
                interpolated["upper_neighbor_index"] = int(right_index)
                delta_center_shapes[delta_index] = interpolated
                center_sources[delta_index] = "neighbor_delta_interpolation"
                continue
        if left_valid or right_valid:
            neighbor_index = left_index if left_valid else right_index
            delta_center_shapes[delta_index] = _nearest_delta_timing_center_shape(
                delta_center_shapes[neighbor_index], neighbor_index,
                delta_centers[neighbor_index], delta_centers[delta_index],
            )
            center_sources[delta_index] = "neighbor_delta_nearest_fallback"
        elif bool(stable_center.get("valid", False)):
            delta_center_shapes[delta_index] = _timing_center_payload(
                stable_center, "stable_global_center_fallback",
                target_delta_center=delta_centers[delta_index],
            )
            center_sources[delta_index] = "stable_global_center_fallback"
        else:
            delta_center_shapes[delta_index] = _timing_center_payload(
                {}, "invalid_timing_center",
                target_delta_center=delta_centers[delta_index],
                rejection_reasons=["no_valid_timing_center_fallback"],
            )
            center_sources[delta_index] = "invalid_timing_center"

    delta_t_cell_fits = []
    support_by_delta_t = []
    support_by_delta = []
    delta_support_rows = []
    support_thresholds = dict(selected_config.get("support_thresholds") or {})
    maximum_deviance = float(
        (selected_config.get("slice_fit") or {}).get("maximum_poisson_deviance_per_entry", float("inf"))
    )
    for delta_index, cell_hists in enumerate(timing_payload["delta_t_cells"]):
        delta_fits = []
        center_shape = delta_center_shapes[delta_index]
        data_total = float(timing_payload["delta_all_t_hists"][delta_index].Integral())
        fitted_data_total = model_total = proton_total = kaon_total = other_total = 0.0
        valid_cells = 0
        for t_index, hist in enumerate(cell_hists):
            timing_constraint = {
                "valid": bool(center_shape.get("valid", False)),
                "timing_center_model_valid": bool(center_shape.get("valid", False)),
                "timing_center_source": center_sources[delta_index],
                "selected_timing_center_source": center_sources[delta_index],
                "predicted_kaon_mean": center_shape.get("kaon_mean"),
                "predicted_proton_mean": center_shape.get("proton_mean"),
                "wrapped_predicted_proton_mean": center_shape.get("proton_mean"),
                "kaon_sigma": center_shape.get("kaon_sigma"),
                "proton_sigma": center_shape.get("proton_sigma"),
            }
            fit = _fit_delta_timing_slice(
                hist, center_shape, selected_config,
                "F_proton_cleaning_delta_{}_t_{}".format(delta_index, t_index),
                use_deviance_per_entry_validation=True,
                maximum_poisson_deviance_per_entry=maximum_deviance,
                support_entries=timing_payload["cell_prompt_support"][delta_index][t_index],
                timing_constraint=timing_constraint,
            )
            fit.update({
                "delta_index": int(delta_index),
                "t_index": int(t_index),
                "t_low": float(t_edges[t_index]),
                "t_high": float(t_edges[t_index + 1]),
                "timing_center_source": center_sources[delta_index],
                "cell_fit_valid": bool(fit.get("valid", False)),
            })
            delta_fits.append(fit)
            if not bool(fit.get("valid", False)):
                continue
            valid_cells += 1
            fitted_data_total += float(fit.get("data_yield", 0.0) or 0.0)
            model_total += float(fit.get("model_yield", 0.0) or 0.0)
            proton_total += float(fit.get("proton_yield", 0.0) or 0.0)
            kaon_total += float(fit.get("kaon_yield", 0.0) or 0.0)
            other_total += float(fit.get("other_yield", 0.0) or 0.0)
        coverage = float(fitted_data_total / data_total) if data_total > 0.0 else 0.0
        delta_label = _classify_timing_t_support(valid_cells, coverage, model_total, support_thresholds)
        for fit in delta_fits:
            fit["support_label"] = delta_label if bool(fit.get("cell_fit_valid", False)) else SUPPORT_UNSUPPORTED
        delta_t_cell_fits.append(delta_fits)
        support_by_delta_t.append([str(fit["support_label"]) for fit in delta_fits])
        support_by_delta.append(delta_label)
        delta_support_rows.append({
            "delta_index": int(delta_index), "support_label": delta_label,
            "valid_t_cells": int(valid_cells), "coverage": float(coverage),
            "data_total": float(data_total), "fitted_data_total": float(fitted_data_total),
            "model_total": float(model_total), "proton_total": float(proton_total),
            "kaon_total": float(kaon_total), "other_total": float(other_total),
            "timing_center_source": center_sources[delta_index],
        })
    setting_support = _classify_timing_t_setting_support(delta_support_rows, support_thresholds)
    _apply_timing_t_cell_map(delta_t_cell_fits, support_by_delta, setting_support)
    if not bool(setting_support.get("accepted", False)):
        rejection_reasons.append("setting_support_gate_rejected")
    return {
        "candidate": candidate,
        "timing_payload": timing_payload,
        "stable_global_shape": stable_global_shape,
        "global_t_shapes": global_t_shapes,
        "delta_timing_center_shapes": delta_center_shapes,
        "timing_center_source_by_delta": center_sources,
        "delta_t_cell_fits": delta_t_cell_fits,
        "support_by_delta": support_by_delta,
        "support_by_delta_t": support_by_delta_t,
        "delta_support": delta_support_rows,
        "setting_support": setting_support,
        "rejection_reasons": rejection_reasons,
    }

def _build_timing_t_event_weight_result(
    result,
    inp_dict,
    phi_setting,
    source_bundle,
    config,
    evaluate_event,
    hole_contains,
    mm_min,
    mm_max,
):
    """Build the non-legacy delta x canonical-t production model."""
    configured_t_edges = inp_dict.get("t_bins")
    t_edges = np.asarray(configured_t_edges if configured_t_edges is not None else [], dtype=float)
    canonical = dict(inp_dict.get("canonical_t_binning") or {})
    t_cfg = dict(config.get("t_binning") or {})
    tolerance = float(t_cfg.get("edge_tolerance", 1.0e-9))
    canonical_source = canonical.get("t_edges")
    canonical_edges = np.asarray(canonical_source if canonical_source is not None else [], dtype=float)
    if (
        t_edges.size < 2
        or canonical_edges.size != t_edges.size
        or not np.allclose(t_edges, canonical_edges, rtol=0.0, atol=tolerance)
    ):
        result["fallback_reason"] = "canonical t edges unavailable or inconsistent"
        if bool(config.get("strict_mode", False)):
            raise RuntimeError(result["fallback_reason"])
        return result

    timing_config = deepcopy(config)
    timing_config["slice_fit"] = dict(config.get("t_cell_fit") or config.get("slice_fit") or {})
    timing_config["support_thresholds"] = dict(
        config.get("t_support_thresholds") or config.get("support_thresholds") or {}
    )
    available_branches = list((source_bundle or {}).get("available_timing_branches") or [])
    candidates = [PROTON_CLEANING_EXACT_TIMING_BRANCH]
    if bool(config.get("allow_rf_probe", True)):
        candidates = [
            branch for branch in _resolve_rf_branch_candidates(source_bundle, config)
            if branch in available_branches
        ] + candidates
    candidates = list(dict.fromkeys(branch for branch in candidates if branch in available_branches))
    if not candidates:
        result["fallback_reason"] = "no configured timing branch available in prepared trees"
        return result

    candidate_evaluations = []
    preferred_branches = list(config.get("timing_t_branch_preference") or candidates)
    branch_preference = {
        str(branch): int(len(preferred_branches) - index)
        for index, branch in enumerate(preferred_branches)
    }
    for branch in candidates:
        candidate = resolve_timing_t_candidate_configuration(
            source_bundle,
            timing_config,
            branch,
            evaluate_event,
            hole_contains,
            mm_min,
            mm_max,
        )
        candidate_config = candidate["candidate_config"]
        payload = _build_signed_timing_t_histograms(
            source_bundle,
            candidate_config,
            t_edges,
            branch,
            time_hist_range=candidate["display_range"],
            time_hist_bins=candidate["histogram_bins"],
        )
        shape = _fit_global_timing_shape(
            payload["H_global_timing"],
            candidate_config,
            "F_proton_cleaning_t_global_{}_{}".format(phi_setting, branch),
            proton_peak_is_lower=bool(candidate["proton_peak_is_lower"]),
            display_range=candidate["fit_range"],
            fit_mode="local_peak_rescue",
        )
        evaluation = _evaluate_timing_t_candidate_production(
            candidate, payload, shape, t_edges, phi_setting
        )
        evaluation["candidate_selection_score"] = _timing_t_candidate_selection_score(
            evaluation["setting_support"], shape,
            branch_preference.get(str(branch), int(-len(candidates))),
        )
        candidate_evaluations.append(evaluation)
    ranked_evaluations = _rank_timing_t_candidate_evaluations(candidate_evaluations)
    selected_evaluation = ranked_evaluations[0]
    candidate_diagnostics = []
    for rank, evaluation in enumerate(ranked_evaluations, start=1):
        candidate = evaluation["candidate"]
        shape = evaluation["stable_global_shape"]
        setting_support = evaluation["setting_support"]
        summary = {
            "timing_branch": candidate["timing_branch"],
            "probe_kind": candidate["probe_kind"],
            "time_hist_range": list(candidate["display_range"]),
            "display_range": candidate["display_range"],
            "fit_range": candidate["fit_range"],
            "time_hist_bins": candidate["histogram_bins"],
            "histogram_bins": candidate["histogram_bins"],
            "range_source": candidate["range_source"],
            "range_diagnostics": candidate["range_diagnostics"],
            "global_fit_valid": bool(shape.get("valid", False)),
            "global_proton_significance": _finite_float_or_none(
                shape.get("proton_component_significance")
            ),
            "poisson_deviance_per_entry": _finite_float_or_none(
                shape.get("poisson_deviance_per_entry")
            ),
            "setting_support_label": setting_support.get("support_label"),
            "setting_accepted": bool(setting_support.get("accepted", False)),
            "supported_delta_count": int(setting_support.get("supported_delta_bins", 0) or 0),
            "marginal_delta_count": int(setting_support.get("marginal_delta_bins", 0) or 0),
            "unsupported_delta_count": int(max(0, len(evaluation["delta_support"]) - int(setting_support.get("supported_delta_bins", 0) or 0) - int(setting_support.get("marginal_delta_bins", 0) or 0))),
            "valid_delta_t_cell_count": int(setting_support.get("valid_t_cells", 0) or 0),
            "valid_coverage": float(setting_support.get("coverage", 0.0) or 0.0),
            "total_modeled_yield": float(setting_support.get("model_total", 0.0) or 0.0),
            "candidate_selection_score": list(evaluation["candidate_selection_score"]),
            "candidate_selection_rank": int(rank),
            "selected": bool(evaluation is selected_evaluation),
            "rejection_reasons": list(evaluation["rejection_reasons"]),
            "candidate_fit_metrics": {
                "global_fit_valid": bool(shape.get("valid", False)),
                "proton_significance": _finite_float_or_none(
                    shape.get("proton_component_significance")
                ),
                "poisson_deviance_per_entry": _finite_float_or_none(
                    shape.get("poisson_deviance_per_entry")
                ),
            },
        }
        if bool(evaluation is selected_evaluation):
            summary["selection_reason"] = (
                "highest production-support-first candidate_selection_score"
            )
        candidate_diagnostics.append(_json_ready_value(summary))
    # The selected candidate owns the detailed delta x t ROOT payload.  Keep
    # only detached global comparison histograms for rejected candidates, then
    # drop their per-cell payloads before this setting result escapes.
    candidate_global_comparisons = []
    rejected_cell_objects_released = 0
    rejected_global_objects_retained = 0
    for evaluation in ranked_evaluations:
        if evaluation is selected_evaluation:
            continue
        candidate = evaluation.get("candidate") or {}
        timing_payload_for_release = evaluation.get("timing_payload") or {}
        comparison = {
            "timing_branch": str(candidate.get("timing_branch", "unknown")),
            "probe_kind": str(candidate.get("probe_kind", "unknown")),
            "H_global_timing": _clone_hist(
                timing_payload_for_release.get("H_global_timing"),
                "H_timing_t_rejected_candidate_global_{}_{}".format(
                    phi_setting, candidate.get("timing_branch", "unknown"),
                ),
                reset=False,
            ),
            "H_global_timing_vs_t": _clone_hist(
                timing_payload_for_release.get("H_global_timing_vs_t"),
                "H_timing_t_rejected_candidate_global_vs_t_{}_{}".format(
                    phi_setting, candidate.get("timing_branch", "unknown"),
                ),
                reset=False,
            ),
        }
        rejected_global_objects_retained += sum(
            1 for value in comparison.values() if _is_root_object(value)
        )
        candidate_global_comparisons.append(comparison)
        rejected_cell_objects_released += len(timing_payload_for_release.get("delta_t_cells") or [])
        rejected_cell_objects_released += len(timing_payload_for_release.get("global_t_slices") or [])
        evaluation["timing_payload"] = {}
        evaluation["global_t_shapes"] = []
        evaluation["delta_timing_center_shapes"] = []
        evaluation["delta_t_cell_fits"] = []
        evaluation["stable_global_shape"] = {}
    if rejected_cell_objects_released:
        gc.collect()
    selected_candidate = selected_evaluation["candidate"]
    timing_payload = selected_evaluation["timing_payload"]
    stable_global_shape = selected_evaluation["stable_global_shape"]
    selected_config = selected_candidate["candidate_config"]
    timing_branch = str(selected_candidate["timing_branch"])
    global_t_shapes = selected_evaluation["global_t_shapes"]
    delta_center_shapes = selected_evaluation["delta_timing_center_shapes"]
    center_sources = selected_evaluation["timing_center_source_by_delta"]
    delta_t_cell_fits = selected_evaluation["delta_t_cell_fits"]
    support_by_delta = selected_evaluation["support_by_delta"]
    support_by_delta_t = selected_evaluation["support_by_delta_t"]
    delta_support_rows = selected_evaluation["delta_support"]
    setting_support = selected_evaluation["setting_support"]

    result.update(
        {
            "accepted": bool(setting_support["accepted"]),
            "implementation": PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED,
            "settings": selected_config,
            "selected_timing_branch": str(timing_branch),
            "t_edges": [float(edge) for edge in t_edges],
            "delta_edges": list(timing_payload["delta_edges"]),
            "support_by_delta": support_by_delta,
            "support_by_delta_t": support_by_delta_t,
            "delta_t_cell_fits": delta_t_cell_fits,
            "delta_timing_center_shapes": delta_center_shapes,
            "global_timing_t_shapes": global_t_shapes,
            "H_global_timing": timing_payload["H_global_timing"],
            "H_global_timing_vs_t": timing_payload["H_global_timing_vs_t"],
            "H_global_timing_t_slices": timing_payload["global_t_slices"],
            "H_delta_timing": timing_payload["H_delta_timing"],
            "H_delta_timing_t_cells": timing_payload["delta_t_cells"],
        }
    )
    result["diagnostics"].update(
        {
            "implementation": PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED,
            "timing_branch": timing_branch,
            "selected_timing_candidate": next(
                row for row in candidate_diagnostics if bool(row.get("selected", False))
            ),
            "timing_candidate_diagnostics": candidate_diagnostics,
            "candidate_selection_tuple": [
                "setting_accepted", "setting_support_rank", "supported_delta_count",
                "marginal_delta_count", "valid_delta_t_cell_count", "valid_coverage",
                "global_fit_valid", "proton_significance", "negative_poisson_deviance_per_entry",
                "configured_branch_preference",
            ],
            "canonical_t_binning": _json_ready_value(canonical),
            "source_stats": _json_ready_value(timing_payload["source_stats"]),
            "timing_center_source_by_delta": center_sources,
            "delta_support": _json_ready_value(delta_support_rows),
            "setting_support": _json_ready_value(setting_support),
            "supported_delta_bins": int(setting_support["supported_delta_bins"]),
            "marginal_delta_bins": int(setting_support["marginal_delta_bins"]),
            "supported_delta_t_cells": int(
                sum(label == SUPPORT_SUPPORTED for labels in support_by_delta_t for label in labels)
            ),
            "applied_timing_t_cell_map": _json_ready_value(
                [dict(fit or {}) for slices in delta_t_cell_fits for fit in (slices or [])]
            ),
            "candidate_root_object_accounting": {
                "selected_candidate_fit_root_objects_retained": int(
                    len(timing_payload.get("global_t_slices") or [])
                    + len(timing_payload.get("delta_all_t_hists") or [])
                    + sum(len(row or []) for row in (timing_payload.get("delta_t_cells") or []))
                    + 3
                ),
                "rejected_candidate_global_root_objects_retained": int(rejected_global_objects_retained),
                "rejected_candidate_per_cell_root_objects_released": int(rejected_cell_objects_released),
            },
        }
    )
    result["_timing_t_candidate_global_comparisons"] = candidate_global_comparisons
    if not result["accepted"]:
        result["fallback_reason"] = "timing-t setting support gate rejected the available delta/t cells"
        return result

    rf_policy = resolve_proton_contamination_cleaning_rf_policy(
        inp_dict=inp_dict, phi_setting=phi_setting
    )
    apply_rf = (
        rf_policy == "epsset_default_after_cleaning"
        and bool(config.get("apply_only_low_epsilon_rf", True))
        and normalize_epsset(inp_dict.get("EPSSET")) == "low"
    )
    result["diagnostics"]["rf_policy"] = rf_policy
    result["diagnostics"]["rf_applied"] = bool(apply_rf)
    result["diagnostics"]["rf_signature_fields"] = list(config.get("rf_signature_fields") or [])
    result["diagnostics"]["signature_round_digits"] = int(
        config.get("signature_round_digits", 9) or 9
    )
    if apply_rf:
        rf_sources = (source_bundle or {}).get("rf_sources") or {}
        missing_rf = [
            name
            for name in ("prompt", "rand", "dummy_prompt", "dummy_rand")
            if rf_sources.get(name, {}).get("tree") is None
        ]
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
    result["_prepared_event_weight_lookup"] = _build_t_prepared_event_weight_lookup(
        result, source_bundle, config
    )
    return result


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
            "kinematic_context": {
                "Q2": inpDict.get("Q2"),
                "W": inpDict.get("W"),
                "epsilon": inpDict.get("EPSVAL", inpDict.get("EPSSET")),
                "epsset": inpDict.get("EPSSET"),
                "phi_setting": phi_setting,
            },
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
            "cross_stage_t_consistency": list(
                (source_bundle or {}).get("cross_stage_t_consistency") or []
            ),
            "cross_stage_t_sampling": dict(
                (source_bundle or {}).get("cross_stage_t_sampling") or {}
            ),
        },
        "fallback_reason": "",
    }
    if method == PROTON_CONTAMINATION_CLEANING_METHOD_DISABLED or not enabled:
        result["fallback_reason"] = "proton cleaning disabled"
        return result
    if method == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT:
        implementation = str(
            config.get("implementation")
            or PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED
        ).strip()
        # The retained global default identifies the legacy method.  Selecting
        # the new method without an override intentionally resolves to its
        # own implementation rather than claiming C-macro parity.
        if implementation == PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_C_SCRIPT_EXACT:
            implementation = PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED
        if implementation != PROTON_CONTAMINATION_CLEANING_IMPLEMENTATION_TIMING_T_BINNED:
            raise ValueError(
                "Unsupported timing-t proton-cleaning implementation '{}'".format(implementation)
            )
        return _build_timing_t_event_weight_result(
            result,
            inpDict,
            phi_setting,
            source_bundle,
            config,
            evaluate_event,
            hole_contains,
            mm_min,
            mm_max,
        )
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
    result["diagnostics"]["tof_required_aliases"] = list(PROTON_CLEANING_TOF_REQUIRED_ALIASES)
    result["diagnostics"]["tof_forbidden_replay_names"] = list(PROTON_CLEANING_TOF_FORBIDDEN_REPLAY_NAMES)
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
    result["diagnostics"]["tof_required_aliases"] = list(PROTON_CLEANING_TOF_REQUIRED_ALIASES)
    result["diagnostics"]["tof_forbidden_replay_names"] = list(PROTON_CLEANING_TOF_FORBIDDEN_REPLAY_NAMES)

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
    low_aero_config = deepcopy(
        exact_config.get("low_aero_offset") or PROTON_CLEANING_LOW_AERO_OFFSET_CONFIG
    )
    primary_low_aero_range = tuple(low_aero_config.get("primary_range") or (0.0, 5.0))
    fallback_low_aero_range = tuple(low_aero_config.get("fallback_range") or (0.0, 6.0))
    full_diagnostic_aero_range = tuple(
        low_aero_config.get("full_diagnostic_range") or PROTON_CLEANING_EXACT_AERO_RANGE
    )
    low_aero_summary_validation = _offset_summary_validation_from_low_aero_config(
        low_aero_config
    )
    primary_low_aero_tof_summaries = _build_prompt_tof_summary_by_delta(
        source_bundle,
        delta_edges,
        selected_selection_key,
        timing_branch=selected_time_branch,
        timing_range=selected_time_hist_range,
        aero_range=primary_low_aero_range,
        delta_range=exact_config.get("delta_hist_range") or PROTON_CLEANING_EXACT_DELTA_RANGE,
        tof_summary_validation=low_aero_summary_validation,
        reason_prefix="low_aero",
    )
    fallback_low_aero_tof_summaries = _build_prompt_tof_summary_by_delta(
        source_bundle,
        delta_edges,
        selected_selection_key,
        timing_branch=selected_time_branch,
        timing_range=selected_time_hist_range,
        aero_range=fallback_low_aero_range,
        delta_range=exact_config.get("delta_hist_range") or PROTON_CLEANING_EXACT_DELTA_RANGE,
        tof_summary_validation=low_aero_summary_validation,
        reason_prefix="low_aero",
    )
    full_aero_tof_summaries = _build_prompt_tof_summary_by_delta(
        source_bundle,
        delta_edges,
        selected_selection_key,
        timing_branch=selected_time_branch,
        timing_range=selected_time_hist_range,
        aero_range=full_diagnostic_aero_range,
        delta_range=exact_config.get("delta_hist_range") or PROTON_CLEANING_EXACT_DELTA_RANGE,
        tof_summary_validation=exact_config.get("tof_summary_validation")
        or PROTON_CLEANING_TOF_SUMMARY_VALIDATION,
    )
    result["diagnostics"]["primary_low_aero_tof_summary_by_delta"] = _json_ready_value(primary_low_aero_tof_summaries)
    result["diagnostics"]["fallback_low_aero_tof_summary_by_delta"] = _json_ready_value(fallback_low_aero_tof_summaries)
    result["diagnostics"]["full_aero_tof_summary_by_delta"] = _json_ready_value(full_aero_tof_summaries)
    result["diagnostics"]["low_aero_offset_config"] = _json_ready_value(low_aero_config)
    result["diagnostics"]["full_aerogel_range"] = list(full_diagnostic_aero_range)
    result["diagnostics"]["aerogel_reference_npe"] = float(low_aero_config.get("show_reference_npe", 5.0) or 5.0)
    result["diagnostics"]["tof_summary_validation"] = _json_ready_value(
        exact_config.get("tof_summary_validation")
        or PROTON_CLEANING_TOF_SUMMARY_VALIDATION
    )
    delta_timing_offset_fits = []
    primary_low_aero_offset_fits = []
    fallback_low_aero_offset_fits = []
    low_aero_projection_payloads = []
    for delta_index, pid_hist in enumerate(pid_payload.get("delta_pid_hists") or []):
        primary_projection = _project_delta_pid_timing_by_aero_range(
            pid_hist,
            "H_proton_cleaning_delta_offset_time_primary_{}".format(delta_index),
            primary_low_aero_range,
            upper_inclusive=False,
        )
        fallback_projection = _project_delta_pid_timing_by_aero_range(
            pid_hist,
            "H_proton_cleaning_delta_offset_time_fallback_{}".format(delta_index),
            fallback_low_aero_range,
            upper_inclusive=False,
        )
        full_projection = _project_delta_pid_timing_by_aero_range(
            pid_hist,
            "H_proton_cleaning_delta_offset_time_full_{}".format(delta_index),
            full_diagnostic_aero_range,
            upper_inclusive=True,
        )
        primary_tof_summary = (
            primary_low_aero_tof_summaries[delta_index]
            if delta_index < len(primary_low_aero_tof_summaries)
            else {}
        )
        fallback_tof_summary = (
            fallback_low_aero_tof_summaries[delta_index]
            if delta_index < len(fallback_low_aero_tof_summaries)
            else {}
        )
        support_by_aero = (
            (pid_payload.get("cell_prompt_support") or [])[delta_index]
            if delta_index < len(pid_payload.get("cell_prompt_support") or [])
            else []
        )
        reference_shape = _build_reference_shape_for_delta_offset(
            global_shapes,
            support_by_aero,
        )
        primary_offset_fit = _fit_delta_common_timing_offset(
            primary_projection,
            reference_shape,
            primary_tof_summary,
            exact_config,
            "f_proton_cleaning_delta_offset_primary_{}".format(delta_index),
            bool(selected_probe.get("proton_peak_is_lower", False)),
            str(selected_probe.get("probe_kind", "ct")),
            selected_time_hist_range,
            beam_bunch_spacing_ns,
            support_entries=int(
                (primary_tof_summary or {}).get(
                    "prompt_events_inside_timing_and_aero_domain",
                    (primary_tof_summary or {}).get("prompt_event_count", 0),
                )
                or 0
            ),
        )
        primary_offset_fit["delta_index"] = int(delta_index)
        primary_offset_fit["offset_fit_aero_mode"] = "low_aero_0_5"
        primary_offset_fit["offset_fit_aero_min"] = float(primary_low_aero_range[0])
        primary_offset_fit["offset_fit_aero_max"] = float(primary_low_aero_range[1])
        primary_offset_fit["selected_projection_name"] = (
            primary_projection.GetName() if primary_projection is not None else ""
        )
        if bool(primary_offset_fit.get("valid", False)):
            fallback_offset_fit = _make_not_required_offset_fit(
                delta_index,
                "low_aero_0_6_fallback",
                "primary_low_aero_offset_fit_valid",
            )
        else:
            fallback_offset_fit = _fit_delta_common_timing_offset(
                fallback_projection,
                reference_shape,
                fallback_tof_summary,
                exact_config,
                "f_proton_cleaning_delta_offset_fallback_{}".format(delta_index),
                bool(selected_probe.get("proton_peak_is_lower", False)),
                str(selected_probe.get("probe_kind", "ct")),
                selected_time_hist_range,
                beam_bunch_spacing_ns,
                support_entries=int(
                    (fallback_tof_summary or {}).get(
                        "prompt_events_inside_timing_and_aero_domain",
                        (fallback_tof_summary or {}).get("prompt_event_count", 0),
                    )
                    or 0
                ),
            )
            fallback_offset_fit["delta_index"] = int(delta_index)
            fallback_offset_fit["offset_fit_aero_mode"] = "low_aero_0_6_fallback"
            fallback_offset_fit["offset_fit_aero_min"] = float(fallback_low_aero_range[0])
            fallback_offset_fit["offset_fit_aero_max"] = float(fallback_low_aero_range[1])
            fallback_offset_fit["selected_projection_name"] = (
                fallback_projection.GetName() if fallback_projection is not None else ""
            )
        offset_fit = _select_resolved_timing_center_model(
            delta_index,
            primary_offset_fit,
            fallback_offset_fit,
            primary_tof_summary,
            fallback_tof_summary,
            reference_shape,
        )
        selected_mode = str(offset_fit.get("selected_timing_center_source") or "unavailable")
        offset_fit["primary_low_aero_offset_fit"] = _json_ready_value(primary_offset_fit)
        offset_fit["fallback_low_aero_offset_fit"] = _json_ready_value(fallback_offset_fit)
        offset_fit["primary_low_aero_prompt_events"] = int(
            (primary_low_aero_tof_summaries[delta_index] if delta_index < len(primary_low_aero_tof_summaries) else {}).get(
                "prompt_events_inside_timing_and_aero_domain",
                0,
            )
            or 0
        )
        offset_fit["primary_low_aero_valid_tof_events"] = int(
            (primary_low_aero_tof_summaries[delta_index] if delta_index < len(primary_low_aero_tof_summaries) else {}).get(
                "usable_tof_events",
                0,
            )
            or 0
        )
        offset_fit["fallback_low_aero_prompt_events"] = int(
            (fallback_low_aero_tof_summaries[delta_index] if delta_index < len(fallback_low_aero_tof_summaries) else {}).get(
                "prompt_events_inside_timing_and_aero_domain",
                0,
            )
            or 0
        )
        offset_fit["fallback_low_aero_valid_tof_events"] = int(
            (fallback_low_aero_tof_summaries[delta_index] if delta_index < len(fallback_low_aero_tof_summaries) else {}).get(
                "usable_tof_events",
                0,
            )
            or 0
        )
        offset_fit["full_range_prompt_events"] = int(
            (full_aero_tof_summaries[delta_index] if delta_index < len(full_aero_tof_summaries) else {}).get(
                "prompt_events_inside_timing_and_aero_domain",
                0,
            )
            or 0
        )
        offset_fit["full_range_valid_tof_events"] = int(
            (full_aero_tof_summaries[delta_index] if delta_index < len(full_aero_tof_summaries) else {}).get(
                "usable_tof_events",
                0,
            )
            or 0
        )
        delta_timing_offset_fits.append(offset_fit)
        primary_low_aero_offset_fits.append(primary_offset_fit)
        fallback_low_aero_offset_fits.append(fallback_offset_fit)
        low_aero_projection_payloads.append(
            {
                "delta_index": int(delta_index),
                "selected_mode": str(selected_mode),
                "primary_projection": primary_projection,
                "fallback_projection": fallback_projection,
                "full_projection": full_projection,
            }
        )
    result["delta_timing_offset_fits"] = delta_timing_offset_fits
    result["primary_low_aero_offset_fits"] = primary_low_aero_offset_fits
    result["fallback_low_aero_offset_fits"] = fallback_low_aero_offset_fits
    result["H_delta_offset_low_aero_projections"] = low_aero_projection_payloads
    result["diagnostics"]["delta_timing_offset_fits"] = _json_ready_value(delta_timing_offset_fits)
    result["diagnostics"]["primary_low_aero_offset_fits_by_delta"] = _json_ready_value(primary_low_aero_offset_fits)
    result["diagnostics"]["fallback_low_aero_offset_fits_by_delta"] = _json_ready_value(fallback_low_aero_offset_fits)
    result["diagnostics"]["delta_offset_rejection_counts"] = _count_rejection_reasons(delta_timing_offset_fits)
    result["diagnostics"]["valid_delta_offset_fits"] = int(
        sum(1 for row in delta_timing_offset_fits if bool((row or {}).get("valid", False)))
    )
    selected_center_counts = {}
    for row in delta_timing_offset_fits:
        source = str((row or {}).get("selected_timing_center_source") or "unavailable")
        selected_center_counts[source] = int(selected_center_counts.get(source, 0) + 1)
    result["diagnostics"]["selected_timing_center_source_counts"] = selected_center_counts
    result["diagnostics"]["selected_timing_center_source_by_delta"] = [
        str((row or {}).get("selected_timing_center_source") or "unavailable")
        for row in delta_timing_offset_fits
    ]
    for diagnostic_key, offset_fits in (
        ("primary_low_aero_offset_fit_counts", primary_low_aero_offset_fits),
        ("fallback_low_aero_offset_fit_counts", fallback_low_aero_offset_fits),
    ):
        result["diagnostics"][diagnostic_key] = {
            "attempted": int(
                sum(1 for row in offset_fits if bool((row or {}).get("fit_attempted", False)))
            ),
            "valid": int(
                sum(1 for row in offset_fits if bool((row or {}).get("valid", False)))
            ),
        }
    result["diagnostics"]["stable_center_fallback_count"] = int(
        sum(
            1
            for row in delta_timing_offset_fits
            if str((row or {}).get("selected_timing_center_source")) == "stable_global_center_fallback"
        )
    )
    selected_tof_summary_by_delta = []
    for delta_index, delta_offset_fit in enumerate(delta_timing_offset_fits):
        primary_tof_summary = (
            primary_low_aero_tof_summaries[delta_index]
            if delta_index < len(primary_low_aero_tof_summaries)
            else {}
        )
        fallback_tof_summary = (
            fallback_low_aero_tof_summaries[delta_index]
            if delta_index < len(fallback_low_aero_tof_summaries)
            else {}
        )
        selected_tof_summary = _build_selected_tof_summary(
            delta_offset_fit,
            primary_tof_summary,
            fallback_tof_summary,
        )
        delta_offset_fit["selected_tof_summary"] = _json_ready_value(selected_tof_summary)
        selected_tof_summary_by_delta.append(selected_tof_summary)
    result["diagnostics"]["selected_tof_summary_by_delta"] = _json_ready_value(
        selected_tof_summary_by_delta
    )
    # Compatibility alias: always the final source-selected summary, never the preliminary choice.
    result["diagnostics"]["tof_summary_by_delta"] = _json_ready_value(
        selected_tof_summary_by_delta
    )
    result["diagnostics"]["tof_summary_rejection_counts"] = _count_rejection_reasons(
        selected_tof_summary_by_delta
    )
    result["diagnostics"]["valid_tof_delta_bins"] = int(
        sum(1 for row in selected_tof_summary_by_delta if bool((row or {}).get("valid", False)))
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
        delta_offset_fit = (
            delta_timing_offset_fits[delta_index]
            if delta_index < len(delta_timing_offset_fits)
            else {"valid": False, "rejection_reason": "missing_delta_offset_fit"}
        )
        tof_summary = (
            selected_tof_summary_by_delta[delta_index]
            if delta_index < len(selected_tof_summary_by_delta)
            else {}
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
        cell_fit_attempt_count = 0
        cell_fit_valid_count = 0
        cell_fit_skipped_count = 0
        cell_fit_skipped_insufficient_support_count = 0
        cell_fit_skipped_missing_histogram_count = 0
        cell_fit_skipped_invalid_global_shape_count = 0
        cell_fit_skipped_invalid_timing_center_model_count = 0
        cell_fit_skipped_other_count = 0
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
                allow_stable_global_center_fallback=bool(
                    (low_aero_config or {}).get("allow_stable_global_center_fallback", True)
                ),
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
            if bool((slice_fit or {}).get("fit_attempted", False)):
                cell_fit_attempt_count += 1
            else:
                cell_fit_skipped_count += 1
                skip_counter_key = _cell_fit_skip_counter_key(
                    (slice_fit or {}).get("fit_status")
                )
                if skip_counter_key == "cell_fit_skipped_insufficient_support_count":
                    cell_fit_skipped_insufficient_support_count += 1
                elif skip_counter_key == "cell_fit_skipped_missing_histogram_count":
                    cell_fit_skipped_missing_histogram_count += 1
                elif skip_counter_key == "cell_fit_skipped_invalid_global_shape_count":
                    cell_fit_skipped_invalid_global_shape_count += 1
                elif skip_counter_key == "cell_fit_skipped_invalid_timing_center_model_count":
                    cell_fit_skipped_invalid_timing_center_model_count += 1
                else:
                    cell_fit_skipped_other_count += 1
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
                        "rejection_reasons": (slice_fit or {}).get("rejection_reasons") or [],
                        "data_yield": (slice_fit or {}).get("data_yield"),
                        "model_yield": (slice_fit or {}).get("model_yield"),
                        "model_data_ratio": (slice_fit or {}).get("model_data_ratio"),
                        "signed_integral": (slice_fit or {}).get("signed_integral"),
                        "absolute_integral": (slice_fit or {}).get("absolute_integral"),
                        "root_entries": (slice_fit or {}).get("root_entries"),
                        "effective_entries": (slice_fit or {}).get("effective_entries"),
                        "kaon_yield": (slice_fit or {}).get("kaon_yield"),
                        "raw_proton_yield": (slice_fit or {}).get("raw_proton_yield"),
                        "proton_yield": (slice_fit or {}).get("proton_yield"),
                        "other_yield": (slice_fit or {}).get("other_yield"),
                        "timing_model_valid": (slice_fit or {}).get("timing_model_valid"),
                        "timing_center_model_valid": (slice_fit or {}).get("timing_center_model_valid"),
                        "timing_center_source": (slice_fit or {}).get("timing_center_source"),
                        "offset_refinement_valid": (slice_fit or {}).get("offset_refinement_valid"),
                        "offset_refinement_applied": (slice_fit or {}).get("offset_refinement_applied"),
                        "offset_refinement_failure_reasons": (slice_fit or {}).get("offset_refinement_failure_reasons"),
                        "cell_fit_valid": (slice_fit or {}).get("cell_fit_valid"),
                        "proton_component_detected": (slice_fit or {}).get("proton_component_detected"),
                        "proton_component_significance": (slice_fit or {}).get("proton_component_significance"),
                        "proton_component_below_significance": (slice_fit or {}).get("proton_component_below_significance"),
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
            cell_fit_valid_count += 1
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
                    "selected_timing_center_source": (delta_offset_fit or {}).get("selected_timing_center_source"),
                    "offset_refinement_applied": bool((delta_offset_fit or {}).get("offset_refinement_applied", False)),
                    "timing_center_model_valid": bool((delta_offset_fit or {}).get("timing_center_model_valid", False)),
                    "cell_fit_attempt_count": int(cell_fit_attempt_count),
                    "cell_fit_valid_count": int(cell_fit_valid_count),
                    "cell_fit_skipped_count": int(cell_fit_skipped_count),
                    "cell_fit_skipped_insufficient_support_count": int(cell_fit_skipped_insufficient_support_count),
                    "cell_fit_skipped_missing_histogram_count": int(cell_fit_skipped_missing_histogram_count),
                    "cell_fit_skipped_invalid_global_shape_count": int(cell_fit_skipped_invalid_global_shape_count),
                    "cell_fit_skipped_invalid_timing_center_model_count": int(cell_fit_skipped_invalid_timing_center_model_count),
                    "cell_fit_skipped_other_count": int(cell_fit_skipped_other_count),
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
    result["diagnostics"]["cell_fit_attempt_count"] = int(
        sum(int((row or {}).get("cell_fit_attempt_count", 0) or 0) for row in delta_support_debug_rows)
    )
    result["diagnostics"]["cell_fit_valid_count"] = int(
        sum(int((row or {}).get("cell_fit_valid_count", 0) or 0) for row in delta_support_debug_rows)
    )
    result["diagnostics"]["cell_fit_skipped_count"] = int(
        sum(int((row or {}).get("cell_fit_skipped_count", 0) or 0) for row in delta_support_debug_rows)
    )
    result["diagnostics"]["cell_fit_skipped_insufficient_support_count"] = int(
        sum(
            int((row or {}).get("cell_fit_skipped_insufficient_support_count", 0) or 0)
            for row in delta_support_debug_rows
        )
    )
    for counter_key in (
        "cell_fit_skipped_missing_histogram_count",
        "cell_fit_skipped_invalid_global_shape_count",
        "cell_fit_skipped_invalid_timing_center_model_count",
        "cell_fit_skipped_other_count",
    ):
        result["diagnostics"][counter_key] = int(
            sum(int((row or {}).get(counter_key, 0) or 0) for row in delta_support_debug_rows)
        )
    n_delta_diag = max(len(delta_fits), max(len(delta_edges) - 1, 0))
    n_aero_diag = max(len(result.get("aero_edges") or []) - 1, len(global_shapes), 0)
    proton_yield_by_delta_aero = []
    raw_proton_yield_by_delta_aero = []
    kaon_yield_by_delta_aero = []
    other_yield_by_delta_aero = []
    total_fit_yield_by_delta_aero = []
    proton_fraction_by_delta_aero = []
    proton_detected_by_delta_aero = []
    for delta_index in range(n_delta_diag):
        slice_collection = delta_fits[delta_index] if delta_index < len(delta_fits) else []
        proton_row = []
        raw_proton_row = []
        kaon_row = []
        other_row = []
        total_row = []
        fraction_row = []
        detected_row = []
        for aero_index in range(n_aero_diag):
            slice_fit = slice_collection[aero_index] if aero_index < len(slice_collection or []) else {}
            proton_value = float((slice_fit or {}).get("proton_yield", 0.0) or 0.0)
            raw_proton_value = float((slice_fit or {}).get("raw_proton_yield", proton_value) or 0.0)
            kaon_value = float((slice_fit or {}).get("kaon_yield", 0.0) or 0.0)
            other_value = float((slice_fit or {}).get("other_yield", 0.0) or 0.0)
            total_value = float(proton_value + kaon_value + other_value)
            proton_row.append(float(proton_value))
            raw_proton_row.append(float(raw_proton_value))
            kaon_row.append(float(kaon_value))
            other_row.append(float(other_value))
            total_row.append(float(total_value))
            fraction_row.append(float(proton_value / total_value) if total_value > 0.0 else 0.0)
            detected_row.append(bool((slice_fit or {}).get("proton_component_detected", False)))
        proton_yield_by_delta_aero.append(proton_row)
        raw_proton_yield_by_delta_aero.append(raw_proton_row)
        kaon_yield_by_delta_aero.append(kaon_row)
        other_yield_by_delta_aero.append(other_row)
        total_fit_yield_by_delta_aero.append(total_row)
        proton_fraction_by_delta_aero.append(fraction_row)
        proton_detected_by_delta_aero.append(detected_row)
    result["diagnostics"]["proton_yield_by_delta_aero"] = proton_yield_by_delta_aero
    result["diagnostics"]["raw_proton_yield_by_delta_aero"] = raw_proton_yield_by_delta_aero
    result["diagnostics"]["kaon_yield_by_delta_aero"] = kaon_yield_by_delta_aero
    result["diagnostics"]["other_yield_by_delta_aero"] = other_yield_by_delta_aero
    result["diagnostics"]["total_fit_yield_by_delta_aero"] = total_fit_yield_by_delta_aero
    result["diagnostics"]["proton_fraction_by_delta_aero"] = proton_fraction_by_delta_aero
    result["diagnostics"]["proton_detected_by_delta_aero"] = proton_detected_by_delta_aero
    aero_edges_for_diag = [float(edge) for edge in (result.get("aero_edges") or [])]
    low_fraction_num = 0.0
    low_fraction_den = 0.0
    high_fraction_num = 0.0
    high_fraction_den = 0.0
    for aero_index in range(n_aero_diag):
        if aero_index + 1 < len(aero_edges_for_diag):
            aero_center = 0.5 * (aero_edges_for_diag[aero_index] + aero_edges_for_diag[aero_index + 1])
        else:
            aero_center = float(aero_index)
        aero_proton = sum(row[aero_index] for row in proton_yield_by_delta_aero if aero_index < len(row))
        aero_total = sum(row[aero_index] for row in total_fit_yield_by_delta_aero if aero_index < len(row))
        if aero_total <= 0.0:
            continue
        if aero_center < 5.0:
            low_fraction_num += float(aero_proton)
            low_fraction_den += float(aero_total)
        if aero_center >= 10.0:
            high_fraction_num += float(aero_proton)
            high_fraction_den += float(aero_total)
    low_aero_proton_fraction = float(low_fraction_num / low_fraction_den) if low_fraction_den > 0.0 else None
    high_aero_proton_fraction = float(high_fraction_num / high_fraction_den) if high_fraction_den > 0.0 else None
    result["diagnostics"]["proton_fraction_below_5_npe"] = low_aero_proton_fraction
    result["diagnostics"]["proton_fraction_above_10_npe"] = high_aero_proton_fraction
    if (
        low_aero_proton_fraction is not None
        and high_aero_proton_fraction is not None
        and high_aero_proton_fraction > low_aero_proton_fraction
    ):
        warning_list = list(result["diagnostics"].get("warnings") or [])
        if "high_aero_proton_fraction_exceeds_low_aero" not in warning_list:
            warning_list.append("high_aero_proton_fraction_exceeds_low_aero")
        result["diagnostics"]["warnings"] = warning_list
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


_HGCER_DIAGNOSTIC_TARGETS = (
    "hgcer_xy",
    "hgcer_x_mm",
    "hgcer_y_mm",
    "hgcer_xy_nohole",
    "hgcer_x_mm_nohole",
    "hgcer_y_mm_nohole",
)


def _new_hgcer_fill_counter():
    return {
        "seen": 0,
        "selected": 0,
        "finite": 0,
        "in_range": 0,
        "filled": 0,
        "nonzero": 0,
        "signed_weight_sum": 0.0,
        "absolute_weight_sum": 0.0,
    }


def _record_hgcer_fill_counter(counter_map, key, hist, x_value, y_value, fill_weight, selected):
    if counter_map is None:
        return
    counter = counter_map.setdefault(str(key), _new_hgcer_fill_counter())
    counter["seen"] += 1
    if not selected:
        return
    counter["selected"] += 1
    if not all(math.isfinite(float(value)) for value in (x_value, y_value, fill_weight)):
        return
    counter["finite"] += 1
    if hist is not None:
        try:
            x_axis, y_axis = hist.GetXaxis(), hist.GetYaxis()
            if (
                float(x_axis.GetXmin()) <= float(x_value) <= float(x_axis.GetXmax())
                and float(y_axis.GetXmin()) <= float(y_value) <= float(y_axis.GetXmax())
            ):
                counter["in_range"] += 1
        except Exception:
            pass
        counter["filled"] += 1
    if float(fill_weight) != 0.0:
        counter["nonzero"] += 1
    counter["signed_weight_sum"] += float(fill_weight)
    counter["absolute_weight_sum"] += abs(float(fill_weight))


def _summarize_hgcer_display_histogram(hist):
    summary = {
        "available": bool(hist is not None), "nonzero_bin_count": 0,
        "signed_integral": 0.0, "absolute_integral": 0.0,
        "positive_bin_count": 0, "negative_bin_count": 0,
    }
    if hist is None:
        return summary
    try:
        for x_bin in range(1, int(hist.GetNbinsX()) + 1):
            for y_bin in range(1, int(hist.GetNbinsY()) + 1):
                value = float(hist.GetBinContent(x_bin, y_bin))
                summary["signed_integral"] += value
                summary["absolute_integral"] += abs(value)
                if value > 0.0:
                    summary["positive_bin_count"] += 1
                    summary["nonzero_bin_count"] += 1
                elif value < 0.0:
                    summary["negative_bin_count"] += 1
                    summary["nonzero_bin_count"] += 1
    except Exception:
        summary["read_error"] = True
    return summary


def audit_timing_t_hgcer_display_targets(cleaning_result, target_hists):
    """Attach final, display-stage HGCer audit without changing histograms."""
    if str((cleaning_result or {}).get("method") or "") != PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT:
        return {}
    application = (cleaning_result or {}).get("application") or {}
    diagnostics = (cleaning_result or {}).setdefault("diagnostics", {})
    display = {
        key: _summarize_hgcer_display_histogram((target_hists or {}).get(key))
        for key in _HGCER_DIAGNOSTIC_TARGETS
    }
    audit = {
        "source": "frozen_target_fill_counters_and_final_display_histograms",
        "adj_mm_source": "shared_shifted_missing_mass",
        "fill_counters": application.get("generic_hgcer_fill_counters") or {},
        "final_display_histograms": display,
    }
    diagnostics["generic_hgcer_diagnostic_integrity"] = audit
    return audit


def _fill_standard_target_hists(
    evt, adj_MM, adj_t, adj_hsdelta, weight, target_hists, allcuts=False,
    nommcuts=False, noholecuts=False, hgcer_counters=None,
):
    hgcer_weight = float(weight) * float(evt.P_hgcer_npeSum)
    for key, x_value, y_value, selected in (
        ("hgcer_xy", evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, allcuts),
        ("hgcer_x_mm", evt.P_hgcer_xAtCer, adj_MM, allcuts),
        ("hgcer_y_mm", evt.P_hgcer_yAtCer, adj_MM, allcuts),
        ("hgcer_xy_nohole", evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, noholecuts),
        ("hgcer_x_mm_nohole", evt.P_hgcer_xAtCer, adj_MM, noholecuts),
        ("hgcer_y_mm_nohole", evt.P_hgcer_yAtCer, adj_MM, noholecuts),
    ):
        _record_hgcer_fill_counter(
            hgcer_counters, key, (target_hists or {}).get(key), x_value, y_value,
            hgcer_weight, bool(selected),
        )
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
    method_is_t_binned = str(cleaning_result.get("method") or "") == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT
    aero_edges = [float(edge) for edge in (cleaning_result.get("aero_edges") or config.get("aero_slice_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0))]
    t_edges = [float(edge) for edge in (cleaning_result.get("t_edges") or [])]
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    secondary_edges = t_edges if method_is_t_binned else aero_edges
    secondary_label = "|t| [GeV^{2}]" if method_is_t_binned else "P_aero NPE"
    secondary_name = "t" if method_is_t_binned else "aero"
    if len(delta_edges) < 2 or len(secondary_edges) < 2:
        return {"accepted": False, "reason": "missing production bin edges"}
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

    h_weight_sum_delta_secondary = ROOT.TH2D(
        "H_proton_weight_sum_delta_{}".format(secondary_name),
        "Average proton weight vs #delta and {};SHMS #delta [%];{}".format(secondary_name, secondary_label),
        len(delta_edges) - 1,
        float(delta_edges[0]),
        float(delta_edges[-1]),
        len(secondary_edges) - 1,
        array("d", secondary_edges),
    )
    h_weight_sum_delta_secondary.SetDirectory(0)
    h_weight_sum_delta_secondary.Sumw2()
    h_weight_norm_delta_secondary = _clone_hist(
        h_weight_sum_delta_secondary,
        "H_proton_weight_norm_delta_{}".format(secondary_name),
        reset=True,
    )
    h_weight_avg_delta_secondary = _clone_hist(
        h_weight_sum_delta_secondary,
        "H_proton_weight_avg_delta_{}".format(secondary_name),
        reset=True,
    )

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
    hgcer_fill_counters = {
        "raw": {}, "estimated_proton": {}, "cleaned_pre_rf": {}, "final_cleaned": {},
    }

    rf_counts = {"accepted": 0, "rejected": 0}
    support_counts = {SUPPORT_SUPPORTED: 0, SUPPORT_MARGINAL: 0, SUPPORT_UNSUPPORTED: 0}
    prepared_sources = _get_prepared_sources(source_bundle)
    diagnostics_payload = cleaning_result.setdefault("diagnostics", {})
    cross_stage_rows = list(diagnostics_payload.get("cross_stage_t_consistency") or [])
    cross_stage_by_signature = {
        str(row.get("event_signature")): row
        for row in cross_stage_rows
        if isinstance(row, dict) and row.get("event_signature")
    }
    cross_stage_tolerance = float(
        ((config.get("t_binning") or {}).get("cross_stage_t_tolerance", 1.0e-10))
    )
    strict_cross_stage = bool(config.get("strict_mode", False))

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
            signature = _make_prepared_event_signature(source_name, entry_index)
            cross_stage_row = cross_stage_by_signature.get(signature)
            if cross_stage_row is not None:
                downstream_t = float(shifted_t_getter(evt))
                prepass_t = float(cross_stage_row.get("prepass_t"))
                prepared_t = float(cross_stage_row.get("prepared_proton_cleaning_adj_t"))
                max_difference = max(
                    abs(prepass_t - prepared_t),
                    abs(prepass_t - downstream_t),
                    abs(prepared_t - downstream_t),
                )
                cross_stage_row["downstream_t"] = downstream_t
                cross_stage_row["maximum_absolute_difference"] = max_difference
                cross_stage_row["consistent"] = bool(max_difference <= cross_stage_tolerance)
                if strict_cross_stage and not cross_stage_row["consistent"]:
                    raise RuntimeError("cross-stage shifted-t mismatch for {}".format(signature))
            frozen_payload = get_kaon_proton_cleaning_event_payload(
                cleaning_result,
                source_name,
                entry_index,
                strict=True,
            )
            delta_value = frozen_payload.get("delta_index", -1)
            secondary_value = frozen_payload.get(
                "t_index" if method_is_t_binned else "aero_index", -1
            )
            # Bin zero is a real canonical cell, not an absent value.
            delta_index = int(delta_value) if isinstance(delta_value, (int, float)) else -1
            secondary_index = int(secondary_value) if isinstance(secondary_value, (int, float)) else -1
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
                hgcer_counters=hgcer_fill_counters["raw"],
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
                hgcer_counters=hgcer_fill_counters["estimated_proton"],
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
                hgcer_counters=hgcer_fill_counters["cleaned_pre_rf"],
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
                    hgcer_counters=hgcer_fill_counters["final_cleaned"],
                )
            else:
                rf_counts["rejected"] += 1
            if 0 <= delta_index < (len(delta_edges) - 1):
                abs_norm = abs(float(coefficient))
                h_weight_sum_delta.Fill(delta_value, abs_norm * proton_weight)
                h_weight_norm_delta.Fill(delta_value, abs_norm)
                secondary_value = adj_t if method_is_t_binned else aero_value
                if 0 <= secondary_index < (len(secondary_edges) - 1):
                    h_weight_sum_delta_secondary.Fill(delta_value, secondary_value, abs_norm * proton_weight)
                    h_weight_norm_delta_secondary.Fill(delta_value, secondary_value, abs_norm)

    for bin_index in range(1, h_weight_avg_delta.GetNbinsX() + 1):
        denominator = float(h_weight_norm_delta.GetBinContent(bin_index))
        numerator = float(h_weight_sum_delta.GetBinContent(bin_index))
        if denominator > 0.0:
            h_weight_avg_delta.SetBinContent(bin_index, numerator / denominator)
    for x_bin in range(1, h_weight_avg_delta_secondary.GetNbinsX() + 1):
        for y_bin in range(1, h_weight_avg_delta_secondary.GetNbinsY() + 1):
            denominator = float(h_weight_norm_delta_secondary.GetBinContent(x_bin, y_bin))
            numerator = float(h_weight_sum_delta_secondary.GetBinContent(x_bin, y_bin))
            if denominator > 0.0:
                h_weight_avg_delta_secondary.SetBinContent(x_bin, y_bin, numerator / denominator)

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
        "generic_hgcer_fill_counters": _json_ready_value(hgcer_fill_counters),
        "H_proton_fraction_vs_MM": h_proton_fraction_vs_mm,
        "H_proton_weight_vs_delta": h_weight_avg_delta,
        "H_proton_weight_vs_delta_{}".format(secondary_name): h_weight_avg_delta_secondary,
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
            "lambda_preservation_gate": _json_ready_value(
                (cleaning_result.get("diagnostics") or {}).get("lambda_preservation_gate") or {}
            ),
        },
    }
    diagnostics_payload["cross_stage_t_consistency"] = _json_ready_value(cross_stage_rows)
    diagnostics_payload["cross_stage_t_consistency_summary"] = {
        "sampling_method": (diagnostics_payload.get("cross_stage_t_sampling") or {}).get(
            "sampling_method", "implementation_consistency_sample"
        ),
        "sample_limit": (diagnostics_payload.get("cross_stage_t_sampling") or {}).get(
            "sample_limit", len(cross_stage_rows)
        ),
        "source_entry_count": (diagnostics_payload.get("cross_stage_t_sampling") or {}).get(
            "source_entry_count", {}
        ),
        "source_selected_entry_count": (diagnostics_payload.get("cross_stage_t_sampling") or {}).get(
            "source_selected_entry_count", {}
        ),
        "sampled_entry_count": (diagnostics_payload.get("cross_stage_t_sampling") or {}).get(
            "sampled_entry_count", len(cross_stage_rows)
        ),
        "sample_count": int(len(cross_stage_rows)),
        "completed_downstream_count": int(
            sum(row.get("downstream_t") is not None for row in cross_stage_rows)
        ),
        "maximum_absolute_difference": max(
            (float(row.get("maximum_absolute_difference", 0.0) or 0.0) for row in cross_stage_rows),
            default=0.0,
        ),
        "tolerance": cross_stage_tolerance,
        "all_consistent": all(bool(row.get("consistent", False)) for row in cross_stage_rows),
    }
    cleaning_result["application"] = application
    return application


def serialize_kaon_proton_cleaning_result(cleaning_result):
    if not isinstance(cleaning_result, dict):
        return {}
    payload = dict(cleaning_result)
    payload.pop("_rf_signature_lookup", None)
    payload.pop("_prepared_event_weight_lookup", None)
    payload.pop("_aerogel_vs_t_root_payload", None)
    payload.pop("_timing_t_mm_root_payload", None)
    payload.pop("_timing_t_candidate_global_comparisons", None)
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
        "Low-aero offset primary/fallback/full: {} / {} / {}".format(
            tuple((diagnostics.get("low_aero_offset_config") or {}).get("primary_range") or (0.0, 5.0)),
            tuple((diagnostics.get("low_aero_offset_config") or {}).get("fallback_range") or (0.0, 6.0)),
            tuple((diagnostics.get("low_aero_offset_config") or {}).get("full_diagnostic_range") or (0.0, 25.0)),
        ),
        "Low-aero offset thresholds prompt/valid/fraction: {} / {} / {}".format(
            int((diagnostics.get("low_aero_offset_config") or {}).get("minimum_prompt_events", 20) or 20),
            int((diagnostics.get("low_aero_offset_config") or {}).get("minimum_valid_tof_events", 10) or 10),
            float((diagnostics.get("low_aero_offset_config") or {}).get("minimum_valid_tof_fraction", 0.50) or 0.50),
        ),
        "Aerogel reference line: {} NPE".format(
            float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0)
        ),
        "Low-aero selected centers: 0-5={} 0-6={} stable={}".format(
            int((diagnostics.get("selected_timing_center_source_counts") or {}).get("low_aero_0_5_fit", 0) or 0),
            int((diagnostics.get("selected_timing_center_source_counts") or {}).get("low_aero_0_6_fit", 0) or 0),
            int((diagnostics.get("selected_timing_center_source_counts") or {}).get("stable_global_center_fallback", 0) or 0),
        ),
        "Low-aero actual fits: 0-5 attempted/valid={}/{} 0-6 attempted/valid={}/{}".format(
            int((diagnostics.get("primary_low_aero_offset_fit_counts") or {}).get("attempted", 0) or 0),
            int((diagnostics.get("primary_low_aero_offset_fit_counts") or {}).get("valid", 0) or 0),
            int((diagnostics.get("fallback_low_aero_offset_fit_counts") or {}).get("attempted", 0) or 0),
            int((diagnostics.get("fallback_low_aero_offset_fit_counts") or {}).get("valid", 0) or 0),
        ),
        "Cell fits: attempted={} valid={} skipped={} low_support={} missing_hist={} invalid_global={} invalid_center={} other={}".format(
            int(diagnostics.get("cell_fit_attempt_count", 0) or 0),
            int(diagnostics.get("cell_fit_valid_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_insufficient_support_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_missing_histogram_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_invalid_global_shape_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_invalid_timing_center_model_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_other_count", 0) or 0),
        ),
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


def _make_aero_axis_hist(name, title, aero_edges):
    edges = [float(edge) for edge in (aero_edges or [])]
    if len(edges) < 2:
        edges = [float(PROTON_CLEANING_EXACT_AERO_RANGE[0]), float(PROTON_CLEANING_EXACT_AERO_RANGE[1])]
    hist = ROOT.TH1D(
        str(name),
        str(title),
        len(edges) - 1,
        array("d", edges),
    )
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def _draw_aerogel_reference_line(hist, reference_npe=5.0, orientation="x", label="5 NPE"):
    if hist is None:
        return []
    try:
        reference_npe = float(reference_npe)
    except Exception:
        return []
    if not math.isfinite(reference_npe):
        return []
    drawn = []
    try:
        if str(orientation).lower() == "y":
            x_axis = hist.GetXaxis()
            y_axis = hist.GetYaxis()
            x_min = float(x_axis.GetXmin())
            x_max = float(x_axis.GetXmax())
            y_min = float(y_axis.GetXmin())
            y_max = float(y_axis.GetXmax())
            if not (y_min <= reference_npe <= y_max):
                return []
            line = TLine(x_min, reference_npe, x_max, reference_npe)
            text = TLatex(x_min + 0.03 * (x_max - x_min), reference_npe + 0.03 * (y_max - y_min), str(label))
        else:
            x_axis = hist.GetXaxis()
            x_min = float(x_axis.GetXmin())
            x_max = float(x_axis.GetXmax())
            if not (x_min <= reference_npe <= x_max):
                return []
            y_min = float(hist.GetMinimum())
            y_max = float(hist.GetMaximum())
            if y_max <= y_min:
                y_min = 0.0
                y_max = 1.0
            line = TLine(reference_npe, y_min, reference_npe, y_max)
            text = TLatex(reference_npe + 0.02 * (x_max - x_min), y_max - 0.08 * (y_max - y_min), str(label))
        line.SetLineColor(kOrange + 7)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw("same")
        text.SetTextColor(kOrange + 7)
        text.SetTextSize(0.028)
        text.Draw("same")
        drawn.extend([line, text])
    except Exception:
        return []
    return drawn


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


def _make_delta_aero_hist(name, title, delta_edges, aero_edges, values=None):
    delta_edges = [float(edge) for edge in (delta_edges or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    aero_edges = [float(edge) for edge in (aero_edges or [])]
    if len(aero_edges) < 2:
        aero_edges = list(PROTON_CLEANING_EXACT_AERO_EDGES)
    hist = ROOT.TH2D(
        str(name),
        str(title),
        len(delta_edges) - 1,
        array("d", delta_edges),
        len(aero_edges) - 1,
        array("d", aero_edges),
    )
    hist.SetDirectory(0)
    hist.Sumw2()
    for delta_index, row in enumerate(values or []):
        if delta_index >= hist.GetNbinsX():
            break
        for aero_index, value in enumerate(row or []):
            if aero_index >= hist.GetNbinsY():
                break
            numeric = _finite_float_or_none(value)
            if numeric is not None:
                hist.SetBinContent(delta_index + 1, aero_index + 1, float(numeric))
    return hist


def _print_low_aero_offset_diagnostics_page(output_pdf, cleaning_result, prefix):
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    offset_fits = (cleaning_result or {}).get("delta_timing_offset_fits") or diagnostics.get("delta_timing_offset_fits") or []
    primary_rows = diagnostics.get("primary_low_aero_tof_summary_by_delta") or []
    fallback_rows = diagnostics.get("fallback_low_aero_tof_summary_by_delta") or []
    full_rows = diagnostics.get("full_aero_tof_summary_by_delta") or []
    if not offset_fits and not primary_rows and not fallback_rows and not full_rows:
        return
    delta_edges = [float(edge) for edge in ((cleaning_result or {}).get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    page_id = abs(id(cleaning_result))
    h_mode = _make_delta_axis_hist(
        "H_proton_cleaning_low_aero_offset_mode_{}".format(page_id),
        "Low-aerogel offset source;SHMS #delta [%];mode code",
        delta_edges,
    )
    h_primary_prompt = _make_delta_axis_hist(
        "H_proton_cleaning_low_aero_primary_prompt_{}".format(page_id),
        "Low-aerogel offset prompt support;SHMS #delta [%];prompt events",
        delta_edges,
    )
    h_fallback_prompt = _clone_hist(h_primary_prompt, "H_proton_cleaning_low_aero_fallback_prompt_{}".format(page_id), reset=True)
    h_full_prompt = _clone_hist(h_primary_prompt, "H_proton_cleaning_low_aero_full_prompt_{}".format(page_id), reset=True)
    h_primary_valid = _make_delta_axis_hist(
        "H_proton_cleaning_low_aero_primary_valid_{}".format(page_id),
        "Low-aerogel offset usable TOF support;SHMS #delta [%];usable TOF events",
        delta_edges,
    )
    h_fallback_valid = _clone_hist(h_primary_valid, "H_proton_cleaning_low_aero_fallback_valid_{}".format(page_id), reset=True)
    h_full_valid = _clone_hist(h_primary_valid, "H_proton_cleaning_low_aero_full_valid_{}".format(page_id), reset=True)
    h_primary_fraction = _make_delta_axis_hist(
        "H_proton_cleaning_low_aero_primary_fraction_{}".format(page_id),
        "Low-aerogel valid TOF fraction;SHMS #delta [%];usable / timing+aero",
        delta_edges,
    )
    h_fallback_fraction = _clone_hist(h_primary_fraction, "H_proton_cleaning_low_aero_fallback_fraction_{}".format(page_id), reset=True)
    h_full_fraction = _clone_hist(h_primary_fraction, "H_proton_cleaning_low_aero_full_fraction_{}".format(page_id), reset=True)
    mode_codes = {
        "unavailable": 0.0,
        "low_aero_0_5": 1.0,
        "low_aero_0_6_fallback": 2.0,
        "low_aero_0_5_fit": 1.0,
        "low_aero_0_6_fit": 2.0,
        "stable_global_center_fallback": 3.0,
    }
    for index in range(max(len(delta_edges) - 1, 0)):
        fit_row = offset_fits[index] if index < len(offset_fits) else {}
        primary = primary_rows[index] if index < len(primary_rows) else {}
        fallback = fallback_rows[index] if index < len(fallback_rows) else {}
        full = full_rows[index] if index < len(full_rows) else {}
        h_mode.SetBinContent(
            index + 1,
            float(
                mode_codes.get(
                    str(
                        fit_row.get("selected_timing_center_source")
                        or fit_row.get("offset_fit_aero_mode", "unavailable")
                    ),
                    0.0,
                )
            ),
        )
        for hist, row in (
            (h_primary_prompt, primary),
            (h_fallback_prompt, fallback),
            (h_full_prompt, full),
        ):
            hist.SetBinContent(index + 1, float((row or {}).get("prompt_events_inside_timing_and_aero_domain", 0) or 0))
        for hist, row in (
            (h_primary_valid, primary),
            (h_fallback_valid, fallback),
            (h_full_valid, full),
        ):
            hist.SetBinContent(index + 1, float((row or {}).get("usable_tof_events", 0) or 0))
        for hist, row in (
            (h_primary_fraction, primary),
            (h_fallback_fraction, fallback),
            (h_full_fraction, full),
        ):
            hist.SetBinContent(index + 1, float((row or {}).get("valid_tof_fraction", 0.0) or 0.0))

    _set_hist_line_marker(h_mode, kBlack, width=3, marker=20)
    _set_hist_line_marker(h_primary_prompt, kBlue, width=3, marker=20)
    _set_hist_line_marker(h_fallback_prompt, kOrange + 7, width=3, marker=21)
    _set_hist_line_marker(h_full_prompt, kGray + 2, width=2, style=2, marker=24)
    _set_hist_line_marker(h_primary_valid, kBlue, width=3, marker=20)
    _set_hist_line_marker(h_fallback_valid, kOrange + 7, width=3, marker=21)
    _set_hist_line_marker(h_full_valid, kGray + 2, width=2, style=2, marker=24)
    _set_hist_line_marker(h_primary_fraction, kBlue, width=3, marker=20)
    _set_hist_line_marker(h_fallback_fraction, kOrange + 7, width=3, marker=21)
    _set_hist_line_marker(h_full_fraction, kGray + 2, width=2, style=2, marker=24)

    canvas = TCanvas(
        "C_proton_cleaning_low_aero_offset_{}".format(page_id),
        "{} proton-cleaning low-aerogel offset diagnostics".format(prefix),
        1800,
        1100,
    )
    canvas.Divide(2, 2)
    drawn_objects = [
        h_mode,
        h_primary_prompt,
        h_fallback_prompt,
        h_full_prompt,
        h_primary_valid,
        h_fallback_valid,
        h_full_valid,
        h_primary_fraction,
        h_fallback_fraction,
        h_full_fraction,
    ]
    canvas.cd(1)
    h_mode.SetMinimum(-0.1)
    h_mode.SetMaximum(2.2)
    h_mode.Draw("hist text")
    info = TPaveText(0.14, 0.70, 0.60, 0.90, "NDC")
    info.SetBorderSize(1)
    info.SetFillStyle(0)
    info.SetTextAlign(12)
    info.SetTextSize(0.030)
    info.AddText("0=unavailable")
    info.AddText("1=primary 0-5 NPE")
    info.AddText("2=fallback 0-6 NPE")
    info.Draw()
    drawn_objects.append(info)
    gPad.Modified()
    gPad.Update()

    canvas.cd(2)
    max_prompt = max(float(h_primary_prompt.GetMaximum()), float(h_fallback_prompt.GetMaximum()), float(h_full_prompt.GetMaximum()), 1.0)
    h_primary_prompt.SetMaximum(1.20 * max_prompt)
    h_primary_prompt.Draw("hist")
    h_fallback_prompt.Draw("hist same")
    h_full_prompt.Draw("hist same")
    legend = TLegend(0.58, 0.68, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_primary_prompt, "0-5 prompt", "l")
    legend.AddEntry(h_fallback_prompt, "0-6 prompt", "l")
    legend.AddEntry(h_full_prompt, "0-25 prompt", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(3)
    max_valid = max(float(h_primary_valid.GetMaximum()), float(h_fallback_valid.GetMaximum()), float(h_full_valid.GetMaximum()), 1.0)
    h_primary_valid.SetMaximum(1.20 * max_valid)
    h_primary_valid.Draw("hist")
    h_fallback_valid.Draw("hist same")
    h_full_valid.Draw("hist same")
    legend = TLegend(0.58, 0.68, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_primary_valid, "0-5 usable TOF", "l")
    legend.AddEntry(h_fallback_valid, "0-6 usable TOF", "l")
    legend.AddEntry(h_full_valid, "0-25 usable TOF", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(4)
    h_primary_fraction.SetMinimum(0.0)
    h_primary_fraction.SetMaximum(1.05)
    h_primary_fraction.Draw("hist")
    h_fallback_fraction.Draw("hist same")
    h_full_fraction.Draw("hist same")
    legend = TLegend(0.58, 0.68, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_primary_fraction, "0-5 fraction", "l")
    legend.AddEntry(h_fallback_fraction, "0-6 fraction", "l")
    legend.AddEntry(h_full_fraction, "0-25 diagnostic", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()
    canvas.Modified()
    canvas.Update()
    gc.collect()
    canvas.Print(output_pdf)


def _print_proton_aerogel_diagnostics_page(output_pdf, cleaning_result, prefix):
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    aero_edges = [float(edge) for edge in ((cleaning_result or {}).get("aero_edges") or PROTON_CLEANING_EXACT_AERO_EDGES)]
    delta_edges = [float(edge) for edge in ((cleaning_result or {}).get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    proton_matrix = diagnostics.get("proton_yield_by_delta_aero") or []
    kaon_matrix = diagnostics.get("kaon_yield_by_delta_aero") or []
    other_matrix = diagnostics.get("other_yield_by_delta_aero") or []
    fraction_matrix = diagnostics.get("proton_fraction_by_delta_aero") or []
    event_sum_matrix = diagnostics.get("event_probability_sum_by_delta_aero") or []
    if not (proton_matrix or kaon_matrix or other_matrix or fraction_matrix or event_sum_matrix):
        return
    page_id = abs(id(cleaning_result))
    n_aero = max(len(aero_edges) - 1, 1)
    h_fraction_aero = _make_aero_axis_hist(
        "H_proton_cleaning_fitted_fraction_aero_{}".format(page_id),
        "Fitted proton fraction versus aerogel;P_aero_npeSum;fitted proton / total",
        aero_edges,
    )
    h_proton_aero = _make_aero_axis_hist(
        "H_proton_cleaning_fitted_proton_yield_aero_{}".format(page_id),
        "Fitted yields versus aerogel;P_aero_npeSum;Fitted yield",
        aero_edges,
    )
    h_kaon_aero = _clone_hist(h_proton_aero, "H_proton_cleaning_fitted_kaon_yield_aero_{}".format(page_id), reset=True)
    h_other_aero = _clone_hist(h_proton_aero, "H_proton_cleaning_fitted_other_yield_aero_{}".format(page_id), reset=True)
    h_event_sum_aero = _make_aero_axis_hist(
        "H_proton_cleaning_event_probability_sum_aero_{}".format(page_id),
        "Applied proton probability sum versus aerogel;P_aero_npeSum;#Sigma w_{p}^{event}",
        aero_edges,
    )
    for aero_index in range(n_aero):
        proton_total = sum(float(row[aero_index]) for row in proton_matrix if aero_index < len(row))
        kaon_total = sum(float(row[aero_index]) for row in kaon_matrix if aero_index < len(row))
        other_total = sum(float(row[aero_index]) for row in other_matrix if aero_index < len(row))
        event_total = sum(float(row[aero_index]) for row in event_sum_matrix if aero_index < len(row))
        total = proton_total + kaon_total + other_total
        h_proton_aero.SetBinContent(aero_index + 1, float(proton_total))
        h_kaon_aero.SetBinContent(aero_index + 1, float(kaon_total))
        h_other_aero.SetBinContent(aero_index + 1, float(other_total))
        h_event_sum_aero.SetBinContent(aero_index + 1, float(event_total))
        h_fraction_aero.SetBinContent(aero_index + 1, float(proton_total / total) if total > 0.0 else 0.0)
    h_fraction_map = _make_delta_aero_hist(
        "H_proton_cleaning_fitted_fraction_delta_aero_{}".format(page_id),
        "Fitted proton fraction by #delta and aerogel;SHMS #delta [%];P_aero_npeSum;fitted proton / total",
        delta_edges,
        aero_edges,
        fraction_matrix,
    )
    h_event_sum_map = _make_delta_aero_hist(
        "H_proton_cleaning_event_probability_sum_delta_aero_{}".format(page_id),
        "Applied proton probability sum by #delta and aerogel;SHMS #delta [%];P_aero_npeSum;#Sigma w_{p}^{event}",
        delta_edges,
        aero_edges,
        event_sum_matrix,
    )
    reference_npe = float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0)
    _set_hist_line_marker(h_fraction_aero, kRed, width=4, marker=20)
    _set_hist_line_marker(h_proton_aero, kRed, width=4, marker=20)
    _set_hist_line_marker(h_kaon_aero, kMagenta + 1, width=4, marker=20)
    _set_hist_line_marker(h_other_aero, kBlue, width=4, marker=20)
    _set_hist_line_marker(h_event_sum_aero, kViolet + 1, width=4, marker=20)
    canvas = TCanvas(
        "C_proton_cleaning_aerogel_diagnostics_{}".format(page_id),
        "{} proton-cleaning aerogel diagnostics".format(prefix),
        1800,
        1100,
    )
    canvas.Divide(2, 2)
    drawn_objects = [h_fraction_aero, h_proton_aero, h_kaon_aero, h_other_aero, h_event_sum_aero, h_fraction_map, h_event_sum_map]
    canvas.cd(1)
    h_fraction_aero.SetMinimum(0.0)
    h_fraction_aero.SetMaximum(1.05)
    h_fraction_aero.Draw("hist")
    drawn_objects.extend(_draw_aerogel_reference_line(h_fraction_aero, reference_npe, "x"))
    gPad.Modified()
    gPad.Update()

    canvas.cd(2)
    max_yield = max(float(h_proton_aero.GetMaximum()), float(h_kaon_aero.GetMaximum()), float(h_other_aero.GetMaximum()), 1.0)
    h_proton_aero.SetMinimum(0.0)
    h_proton_aero.SetMaximum(1.20 * max_yield)
    h_proton_aero.Draw("hist")
    h_kaon_aero.Draw("hist same")
    h_other_aero.Draw("hist same")
    drawn_objects.extend(_draw_aerogel_reference_line(h_proton_aero, reference_npe, "x"))
    legend = TLegend(0.62, 0.68, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(h_proton_aero, "proton", "l")
    legend.AddEntry(h_kaon_aero, "kaon", "l")
    legend.AddEntry(h_other_aero, "other", "l")
    legend.Draw()
    drawn_objects.append(legend)
    gPad.Modified()
    gPad.Update()

    canvas.cd(3)
    gPad.SetRightMargin(0.16)
    h_fraction_map.SetMinimum(0.0)
    h_fraction_map.SetMaximum(1.0)
    h_fraction_map.Draw("colz text")
    drawn_objects.extend(_draw_aerogel_reference_line(h_fraction_map, reference_npe, "y"))
    gPad.Modified()
    gPad.Update()

    canvas.cd(4)
    h_event_sum_aero.SetMinimum(0.0)
    h_event_sum_aero.SetMaximum(max(1.0, 1.20 * float(h_event_sum_aero.GetMaximum())))
    h_event_sum_aero.Draw("hist")
    drawn_objects.extend(_draw_aerogel_reference_line(h_event_sum_aero, reference_npe, "x"))
    info = TPaveText(0.52, 0.68, 0.88, 0.88, "NDC")
    info.SetBorderSize(1)
    info.SetFillStyle(0)
    info.SetTextAlign(12)
    info.SetTextSize(0.030)
    info.AddText("fraction <5 NPE: {}".format(_format_debug_float(diagnostics.get("proton_fraction_below_5_npe"), digits=3)))
    info.AddText("fraction >10 NPE: {}".format(_format_debug_float(diagnostics.get("proton_fraction_above_10_npe"), digits=3)))
    warnings = ", ".join(str(w) for w in (diagnostics.get("warnings") or [])) or "none"
    info.AddText("warnings: {}".format(warnings[:80]))
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
    method_is_t_binned = str(cleaning_result.get("method") or "") == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    if len(delta_edges) < 2:
        delta_edges = np.linspace(
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[0]),
            float(PROTON_CLEANING_EXACT_DELTA_RANGE[1]),
            int(PROTON_CLEANING_EXACT_DELTA_BINS) + 1,
        ).tolist()
    secondary_edges = [float(edge) for edge in (
        cleaning_result.get("t_edges") if method_is_t_binned else cleaning_result.get("aero_edges")
    ) or (PROTON_CLEANING_EXACT_AERO_EDGES if not method_is_t_binned else [0.0, 1.0])]
    secondary_label = "|t| [GeV^{2}]" if method_is_t_binned else "P_aero_npeSum"
    secondary_key = "H_proton_weight_vs_delta_t" if method_is_t_binned else "H_proton_weight_vs_delta_aero"
    closure_yield_label = "applied p yield" if method_is_t_binned else "fitted p yield"
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
        application.get(secondary_key),
        "H_proton_cleaning_summary_applied_weight_map_{}".format(page_id),
        reset=False,
    )
    if h_applied_map is None:
        h_applied_map = ROOT.TH2D(
            "H_proton_cleaning_summary_applied_weight_map_empty_{}".format(page_id),
            "Mean applied proton probability;SHMS #delta [%];{};#LTw_{{p}}^{{event}}#GT".format(secondary_label),
            n_delta_bins,
            array("d", delta_edges),
            max(len(secondary_edges) - 1, 1),
            array("d", secondary_edges if len(secondary_edges) >= 2 else [0.0, 1.0]),
        )
        h_applied_map.SetDirectory(0)
        h_applied_map.Sumw2()
    h_applied_map.SetTitle(
        "Mean applied proton probability;SHMS #delta [%];{};#LTw_{{p}}^{{event}}#GT".format(secondary_label)
    )

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
    drawn_objects.extend(
        _draw_aerogel_reference_line(
            h_applied_map,
            float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0),
            "y",
        )
    )
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
    method_is_t_binned = (
        str((cleaning_result or {}).get("method") or "")
        == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT
    )
    closure_yield_label = (
        "applied p yield" if method_is_t_binned else "fitted p yield"
    )
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
        if bool(row.get("valid", False)) and reference_k is not None and offset is not None:
            _set_hist_bin_if_finite(h_k_ref, bin_index, reference_k + offset, row.get("delta_offset_error"))
        if bool(row.get("valid", False)) and wrapped_p is not None and offset is not None:
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
    offset_counts = diagnostics.get("delta_offset_rejection_counts") or _count_rejection_reasons(offset_fits)
    tof_counts = diagnostics.get("tof_summary_rejection_counts") or _count_rejection_reasons(tof_summaries)
    weak_component = (
        int(offset_counts.get("kaon_significance_below_min", 0))
        + int(offset_counts.get("proton_significance_below_min", 0))
        + int(offset_counts.get("smaller_component_fraction_below_min", 0))
    )
    tof_support = (
        int(tof_counts.get("insufficient_prompt_tof_events", 0))
        + int(tof_counts.get("insufficient_low_aero_prompt_events", 0))
        + int(tof_counts.get("insufficient_low_aero_valid_tof_events", 0))
        + int(tof_counts.get("valid_tof_fraction_below_min", 0))
        + int(tof_counts.get("low_aero_valid_tof_fraction_below_min", 0))
        + int(tof_counts.get("invalid_mean_delta_t_pk_ns", 0))
        + int(tof_counts.get("invalid_low_aero_mean_delta_t_pk", 0))
        + int(tof_counts.get("invalid_mean_P_gtr_p", 0))
        + int(tof_counts.get("invalid_low_aero_mean_P_gtr_p", 0))
        + int(tof_counts.get("invalid_mean_shms_path_length_cm", 0))
        + int(tof_counts.get("invalid_low_aero_mean_shms_path_length_cm", 0))
    )
    offset_info = TPaveText(0.14, 0.62, 0.62, 0.90, "NDC")
    offset_info.SetBorderSize(1)
    offset_info.SetFillStyle(0)
    offset_info.SetTextAlign(12)
    offset_info.SetTextSize(0.026)
    offset_info.AddText("valid/invalid: {} / {}".format(valid_offsets, invalid_offsets))
    offset_info.AddText("bound-hit warnings: {}".format(bound_offsets))
    offset_info.AddText("large-error: {}".format(int(offset_counts.get("offset_error_exceeds_max", 0))))
    offset_info.AddText("poor-goodness: {}".format(int(offset_counts.get("chi2_ndf_exceeds_max", 0))))
    offset_info.AddText("weak-component: {}".format(int(weak_component)))
    offset_info.AddText("bound-hit weak: {}".format(int(offset_counts.get("bound_hit_with_weak_constraint", 0))))
    offset_info.AddText("TOF-support: {}".format(int(tof_support)))
    offset_info.AddText(
        "skipped low/missing/global/center: {} / {} / {} / {}".format(
            int(diagnostics.get("cell_fit_skipped_insufficient_support_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_missing_histogram_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_invalid_global_shape_count", 0) or 0),
            int(diagnostics.get("cell_fit_skipped_invalid_timing_center_model_count", 0) or 0),
        )
    )
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
        graph.SetTitle("Event-weight closure by #delta;SHMS #delta [%];#Sigma w_{{p}}^{{event}} / {}".format(closure_yield_label))
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
            "Event-weight closure by #delta;SHMS #delta [%];#Sigma w_{{p}}^{{event}} / {}".format(closure_yield_label),
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


def _timing_t_layout_for_panel_count(panel_count):
    """Return the compact, deterministic timing-t report layout."""
    count = max(1, int(panel_count or 1))
    if count == 1:
        return 1, 1
    if count == 2:
        return 2, 1
    if count <= 4:
        return 2, 2
    if count <= 6:
        return 3, 2
    if count <= 8:
        return 4, 2
    if count <= 10:
        return 5, 2
    return 4, 3


def _cross_stage_visual_state(rows, validation_cfg=None):
    threshold = abs(float((validation_cfg or {}).get("cross_stage_visual_threshold", 0.0) or 0.0))
    maximum = max(
        [abs(float((row or {}).get("maximum_absolute_difference", 0.0) or 0.0)) for row in (rows or [])] or [0.0]
    )
    return {
        "visual_threshold": threshold,
        "maximum_absolute_difference": maximum,
        "render_pass_status": bool(maximum <= threshold),
    }


def _timing_t_per_t_pid_eligible(t_summary, validation_cfg=None):
    tolerance = abs(float((validation_cfg or {}).get("per_t_absolute_support_tolerance", 0.0) or 0.0))
    prompt_count = int((t_summary or {}).get("raw_prompt_event_count", 0) or 0)
    absolute_support = abs(float((t_summary or {}).get("absolute_event_weight_support", 0.0) or 0.0))
    return {
        "eligible": bool(prompt_count > 0 and absolute_support > tolerance),
        "prompt_count": prompt_count,
        "absolute_support": absolute_support,
        "absolute_support_tolerance": tolerance,
    }


def _timing_t_signed_display_limits(minimum, maximum, padding_fraction=0.12):
    """Pad a local signed-yield range without imposing a unit-scale floor."""
    minimum = _finite_float_or_none(minimum)
    maximum = _finite_float_or_none(maximum)
    minimum = 0.0 if minimum is None else float(minimum)
    maximum = 0.0 if maximum is None else float(maximum)
    y_low = min(minimum, 0.0)
    y_high = max(maximum, 0.0)
    span = y_high - y_low
    if span <= 0.0:
        span = max(abs(y_high), 1.0e-12)
    padding = max(0.0, float(padding_fraction)) * span
    return (
        y_low - padding if y_low < 0.0 else 0.0,
        y_high + padding,
    )


def _timing_t_root_content_status(hist, *, tolerance=0.0, categorical=False):
    """Classify a ROOT payload without mistaking signed cancellation for empty."""
    if hist is None:
        return {"state": "unavailable", "reason": "missing_histogram"}
    if categorical:
        return {"state": "available", "reason": "categorical_state_map"}
    tolerance = abs(float(tolerance or 0.0))
    nonzero_bins = 0
    absolute_content = 0.0
    maximum_absolute_content = 0.0
    try:
        for x_bin in range(1, int(hist.GetNbinsX()) + 1):
            y_bins = int(hist.GetNbinsY()) if hasattr(hist, "GetNbinsY") else 0
            if y_bins > 1:
                values = (float(hist.GetBinContent(x_bin, y_bin)) for y_bin in range(1, y_bins + 1))
            else:
                values = (float(hist.GetBinContent(x_bin)),)
            for value in values:
                absolute_content += abs(value)
                maximum_absolute_content = max(maximum_absolute_content, abs(value))
                if abs(value) > tolerance:
                    nonzero_bins += 1
    except Exception:
        return {"state": "available", "reason": "histogram_content_unavailable"}
    return {
        "state": "available" if nonzero_bins else "empty",
        "reason": "nonzero_bins" if nonzero_bins else "zero_content",
        "nonzero_bins": int(nonzero_bins),
        "absolute_content": float(absolute_content),
        "maximum_absolute_content": float(maximum_absolute_content),
    }


def _draw_timing_t_status_panel(message, *, color=kRed):
    status = TPaveText(0.12, 0.35, 0.88, 0.65, "NDC")
    status.SetBorderSize(1)
    status.SetFillStyle(0)
    status.SetTextAlign(22)
    status.SetTextSize(0.040)
    status.SetTextColor(color)
    for line in str(message).split("\n"):
        status.AddText(line)
    status.Draw()
    return status


def _make_timing_t_report_canvas(
    name, title, width, height, columns=1, rows=1, header_text=None, *, panel_count=None,
):
    """Create a canvas with a dedicated, non-data header pad."""
    canvas = TCanvas(str(name), str(title), int(width), int(height))
    header = TPad("{}_header".format(name), "", 0.0, 0.91, 1.0, 1.0)
    body = TPad("{}_body".format(name), "", 0.0, 0.0, 1.0, 0.91)
    header.SetFillStyle(0)
    header.SetBorderMode(0)
    body.SetFillStyle(0)
    body.SetBorderMode(0)
    header.Draw()
    body.Draw()
    header.cd()
    label = TPaveText(0.02, 0.10, 0.98, 0.90, "NDC")
    label.SetBorderSize(0)
    label.SetFillStyle(0)
    label.SetTextAlign(12)
    label.SetTextSize(0.026)
    for line in (header_text if isinstance(header_text, (list, tuple)) else (header_text or title,)):
        label.AddText(str(line))
    label.Draw()
    body.cd()
    if panel_count is not None:
        columns, rows = _timing_t_layout_for_panel_count(panel_count)
    if int(columns) * int(rows) > 1:
        body.Divide(int(columns), int(rows))
    return canvas, body, header, label


def _draw_categorical_legend(code_map, title):
    legend = TPaveText(0.14, 0.70, 0.86, 0.90, "NDC")
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.SetTextAlign(12)
    legend.SetTextSize(0.030)
    legend.AddText(str(title))
    legend.AddText(";  ".join("{} = {}".format(code, label) for label, code in code_map.items()))
    legend.Draw()
    return legend


def _print_timing_t_closure_pages(output_pdf, cleaning_result, prefix, page_id):
    """Section F: frozen lookup closure after aerogel diagnostics."""
    diagnostics = cleaning_result.get("diagnostics") or {}
    rows = list(diagnostics.get("cross_stage_t_consistency") or [])
    summary = diagnostics.get("cross_stage_t_consistency_summary") or {}
    aero_cfg = ((diagnostics.get("aerogel_vs_t_validation") or {}).get("configuration") or {})
    visual_state = _cross_stage_visual_state(rows, aero_cfg)
    visual_threshold = visual_state["visual_threshold"]
    maximum_difference = visual_state["maximum_absolute_difference"]
    canvas, body, _header, _label = _make_timing_t_report_canvas(
        "C_timing_t_cross_stage_final_{}".format(page_id),
        "{} cross-stage shifted-t consistency".format(prefix), 1100, 750,
        header_text="Frozen lookup cross-stage shifted-t consistency", panel_count=1,
    )
    body.cd()
    if visual_state["render_pass_status"]:
        _draw_timing_t_status_panel(
            "PASS: all sampled cross-stage shifted-t differences are <= {:.6g}\n"
            "samples={} | maximum={:.6g} | strict tolerance={}".format(
                visual_threshold,
                summary.get("sample_count", len(rows)),
                maximum_difference,
                summary.get("tolerance", "n/a"),
            ),
            color=kGreen + 2,
        )
    else:
        hist = TH1D(
            "H_timing_t_cross_stage_difference_final_{}".format(page_id),
            "Cross-stage shifted-t maximum difference;sample index;max |#Delta t|",
            max(len(rows), 1), 0.0, float(max(len(rows), 1)),
        )
        hist.SetDirectory(0)
        for index, row in enumerate(rows, start=1):
            hist.SetBinContent(index, float(row.get("maximum_absolute_difference", 0.0) or 0.0))
        hist.Draw("hist")
        info = TPaveText(0.52, 0.70, 0.90, 0.90, "NDC")
        info.SetBorderSize(1)
        info.SetFillStyle(0)
        info.AddText("samples: {}".format(summary.get("sample_count", len(rows))))
        info.AddText("max difference: {}".format(maximum_difference))
        info.AddText("visual threshold: {}".format(visual_threshold))
        info.AddText("strict tolerance: {}".format(summary.get("tolerance", "n/a")))
        info.Draw()
    canvas.Print(output_pdf)

    canonical = diagnostics.get("canonical_t_binning") or {}
    t_edges = [float(edge) for edge in (cleaning_result.get("t_edges") or [])]
    phi_edges = [float(edge) for edge in (canonical.get("phi_edges") or [])]
    rows = list(diagnostics.get("event_weight_lookup_by_t_phi") or [])
    if len(t_edges) < 2 or len(phi_edges) < 2 or not rows:
        return
    count_hist = ROOT.TH2D(
        "H_timing_t_lookup_tphi_count_final_{}".format(page_id),
        "Frozen lookup event count;|t| [GeV^{2}];#phi [deg]",
        len(t_edges) - 1, array("d", t_edges), len(phi_edges) - 1, array("d", phi_edges),
    )
    count_hist.SetDirectory(0)
    proton_hist = _clone_hist(
        count_hist, "H_timing_t_lookup_tphi_proton_final_{}".format(page_id), reset=True,
    )
    for row in rows:
        t_index, phi_index = row.get("t_index"), row.get("phi_index")
        if isinstance(t_index, int) and isinstance(phi_index, int):
            if 0 <= t_index < count_hist.GetNbinsX() and 0 <= phi_index < count_hist.GetNbinsY():
                count_hist.SetBinContent(t_index + 1, phi_index + 1, float(row.get("event_count", 0) or 0))
                proton_hist.SetBinContent(t_index + 1, phi_index + 1, float(row.get("estimated_proton_yield", 0.0) or 0.0))
    canvas = TCanvas("C_timing_t_lookup_tphi_final_{}".format(page_id), "{} frozen lookup by |t| and phi".format(prefix), 1400, 700)
    canvas.Divide(2, 1)
    canvas.cd(1)
    count_hist.Draw("colz text")
    canvas.cd(2)
    proton_hist.SetTitle("Frozen lookup estimated proton yield;|t| [GeV^{2}];#phi [deg]")
    proton_hist.Draw("colz text")
    canvas.Print(output_pdf)


def _print_timing_t_mm_diagnostic_pages(output_pdf, cleaning_result, prefix, page_id):
    """Render frozen-row proton-subtraction MM diagnostics for timing-t only."""
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    mm_payload = diagnostics.get("timing_t_mm_diagnostics") or {}
    if not bool(mm_payload.get("enabled", False)):
        return
    config = dict(mm_payload.get("configuration") or {})
    root_payload = (cleaning_result or {}).get("_timing_t_mm_root_payload") or {}
    maps = root_payload.get("maps") or {}
    per_t_rows = list(mm_payload.get("per_t_bin_summary") or [])
    t_edges = [float(value) for value in (mm_payload.get("t_edges") or [])]
    if len(t_edges) < 2:
        return
    header = (
        "Frozen timing-|t| proton-subtraction MM diagnostics - proposed model versus committed applied factors"
    )
    lambda_status = str(
        (diagnostics.get("lambda_preservation_gate") or {}).get("status") or "not evaluated"
    ).upper()
    header = "{} | setting-wide Lambda gate={}".format(header, lambda_status)
    if bool(config.get("write_overview_page", True)):
        canvas, body, _header, _label = _make_timing_t_report_canvas(
            "C_timing_t_mm_overview_{}".format(page_id),
            "{} proton-subtraction MM overview".format(prefix), 1600, 1000,
            header_text=header, panel_count=4,
        )
        for index, (key, label) in enumerate((
            ("raw", "Raw selected MM"),
            ("proposed_estimated_proton", "Proposed estimated proton MM"),
            ("proposed_cleaned_pre_rf", "Proposed cleaned pre-RF MM"),
            ("applied_cleaned_pre_rf", "Final applied cleaned pre-RF MM"),
        ), start=1):
            body.cd(index)
            hist = maps.get(key)
            status = _timing_t_root_content_status(hist)
            if status["state"] in ("unavailable", "empty"):
                _draw_timing_t_status_panel(
                    "{}\n{}".format(label, status.get("reason", "unavailable")),
                )
            else:
                hist.Draw("colz")
                gPad.Modified()
                gPad.Update()
        canvas.Print(output_pdf)

    if bool(config.get("write_t_binned_pages", True)):
        raw_map = maps.get("raw")
        proton_map = maps.get("proposed_estimated_proton")
        cleaned_map = maps.get("proposed_cleaned_pre_rf")
        final_map = maps.get("applied_cleaned_pre_rf")
        panels_per_page = max(1, min(12, int(config.get("max_t_panels_per_page", 1) or 1)))
        for first in range(0, len(per_t_rows), panels_per_page):
            rows = per_t_rows[first:first + panels_per_page]
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_mm_by_t_{}_{}".format(page_id, first),
                "{} proton-subtraction MM by |t| bin".format(prefix), 1600, 1000,
                header_text=[header, "raw / proposed proton / proposed cleaned / final applied cleaned before RF"],
                panel_count=len(rows),
            )
            for panel, row in enumerate(rows, start=1):
                body.cd(panel)
                t_index = int(row.get("t_index", -1))
                support = abs(float(row.get("absolute_event_weight_support", 0.0) or 0.0))
                if (
                    support <= 0.0 or raw_map is None or proton_map is None
                    or cleaned_map is None or not (0 <= t_index < len(t_edges) - 1)
                ):
                    _draw_timing_t_status_panel(
                        "|t| bin {} [{:.6g}, {:.6g}]\nNo supported frozen MM events".format(
                            t_index,
                            float(row.get("t_low", 0.0) or 0.0),
                            float(row.get("t_high", 0.0) or 0.0),
                        ), color=kBlack,
                    )
                    continue
                projections = []
                for key, hist, color, label in (
                    ("raw", raw_map, kBlack, "raw"),
                    ("proposed_proton", proton_map, kRed, "proposed estimated proton"),
                    ("proposed_cleaned", cleaned_map, kGreen + 2, "proposed cleaned pre-RF"),
                    ("applied_cleaned", final_map, kBlue, "final applied cleaned pre-RF"),
                ):
                    if hist is None:
                        continue
                    projection = hist.ProjectionX(
                        "H_timing_t_mm_{}_projection_{}_{}".format(key, page_id, t_index),
                        t_index + 1, t_index + 1, "e",
                    )
                    projection.SetDirectory(0)
                    _set_hist_line_marker(projection, color, width=2)
                    projections.append((label, projection))
                if not projections:
                    _draw_timing_t_status_panel("No MM projection available", color=kBlack)
                    continue
                raw_status = _timing_t_root_content_status(projections[0][1])
                if raw_status["state"] == "empty":
                    display_range = mm_payload.get("display_range") or ["n/a", "n/a"]
                    _draw_timing_t_status_panel(
                        "|t| bin {} [{:.6g}, {:.6g}]\n"
                        "No populated signed-MM bins in display range [{}, {}]\n"
                        "absolute frozen support={:.5g}".format(
                            t_index,
                            float(row.get("t_low", 0.0) or 0.0),
                            float(row.get("t_high", 0.0) or 0.0),
                            display_range[0], display_range[1], support,
                        ), color=kBlack,
                    )
                    continue
                maximum = max([float(hist.GetMaximum()) for _, hist in projections] or [1.0])
                minimum = min([float(hist.GetMinimum()) for _, hist in projections] or [0.0])
                first_projection = projections[0][1]
                first_projection.SetTitle(
                    "|t| bin {} [{:.6g}, {:.6g}] GeV^{{2}};shifted MM [GeV];signed physical yield".format(
                        t_index, float(row.get("t_low", 0.0) or 0.0), float(row.get("t_high", 0.0) or 0.0),
                    )
                )
                # The signed physical yields can be much smaller than one.
                # Do not impose a unit floor: it makes a valid proton-removal
                # structure visually disappear in a small |t| interval.
                y_min, y_max = _timing_t_signed_display_limits(minimum, maximum)
                first_projection.SetMinimum(y_min)
                first_projection.SetMaximum(y_max)
                first_projection.Draw("hist")
                for _, projection in projections[1:]:
                    projection.Draw("hist same")
                legend = TLegend(0.48, 0.64, 0.88, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                for label, projection in projections:
                    legend.AddEntry(projection, label, "l")
                legend.Draw()
                gPad.Modified()
                gPad.Update()
            canvas.Print(output_pdf)

    if bool(config.get("write_window_accounting_page", True)):
        window_names = list((mm_payload.get("validation_windows") or {}).keys())
        rows_per_page = 8
        aggregate = mm_payload.get("aggregate") or {}
        for first in range(0, max(len(per_t_rows), 1), rows_per_page):
            rows = per_t_rows[first:first + rows_per_page]
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_mm_windows_{}_{}".format(page_id, first),
                "{} proton-subtraction MM window accounting".format(prefix), 1500, 900,
                header_text=header, panel_count=1,
            )
            body.cd()
            text = TPaveText(0.03, 0.04, 0.97, 0.96, "NDC")
            text.SetBorderSize(1)
            text.SetFillStyle(0)
            text.SetTextAlign(12)
            text.SetTextSize(0.024)
            if first == 0:
                text.AddText(
                    "setting total: raw={:.6g}, estimated p={:.6g}, cleaned pre-RF={:.6g}, closure(raw-p-cleaned)={:.3g}".format(
                        float(aggregate.get("raw_missing_mass_yield", 0.0) or 0.0),
                        float(aggregate.get("estimated_proton_missing_mass_yield", 0.0) or 0.0),
                        float(aggregate.get("cleaned_missing_mass_yield", 0.0) or 0.0),
                        float(aggregate.get("pre_rf_cleaning_closure_difference", 0.0) or 0.0),
                    )
                )
            for row in rows:
                text.AddText(
                    "|t|{} [{:.6g},{:.6g}]: raw={:.5g}, p={:.5g}, cleaned={:.5g}, closure={:.3g}".format(
                        row.get("t_index"), float(row.get("t_low", 0.0) or 0.0),
                        float(row.get("t_high", 0.0) or 0.0),
                        float(row.get("raw_missing_mass_yield", 0.0) or 0.0),
                        float(row.get("estimated_proton_missing_mass_yield", 0.0) or 0.0),
                        float(row.get("cleaned_missing_mass_yield", 0.0) or 0.0),
                        float(row.get("pre_rf_cleaning_closure_difference", 0.0) or 0.0),
                    )
                )
                for window_name in window_names:
                    window = (row.get("windows") or {}).get(window_name) or {}
                    fraction = window.get("removed_fraction")
                    fraction_text = "n/a" if fraction is None else "{:.4g}".format(float(fraction))
                    text.AddText(
                        "  {} [{:.6g},{:.6g}]: raw={:.5g}, p={:.5g}, cleaned={:.5g}, removed fraction={}".format(
                            window_name,
                            float((window.get("range") or [0.0, 0.0])[0]),
                            float((window.get("range") or [0.0, 0.0])[1]),
                            float(window.get("raw_yield", 0.0) or 0.0),
                            float(window.get("estimated_proton_yield", 0.0) or 0.0),
                            float(window.get("cleaned_yield", 0.0) or 0.0),
                            fraction_text,
                        )
                    )
            text.Draw()
            canvas.Print(output_pdf)


def _print_timing_t_validation_pages(output_pdf, cleaning_result, prefix):
    """Append provenance and observational t x aerogel pages for timing-t."""
    if str(cleaning_result.get("method") or "") != PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT:
        return
    diagnostics = cleaning_result.get("diagnostics") or {}
    canonical = diagnostics.get("canonical_t_binning") or {}
    page_id = abs(id(cleaning_result))
    canvas = TCanvas(
        "C_timing_t_canonical_provenance_{}".format(page_id),
        "{} |t| bin provenance".format(prefix), 1200, 800,
    )
    provenance = TPaveText(0.05, 0.06, 0.95, 0.94, "NDC")
    provenance.SetBorderSize(1)
    provenance.SetFillStyle(0)
    provenance.SetTextAlign(12)
    provenance.SetTextSize(0.026)
    provenance.AddText("|t| production-bin provenance")
    for label, value in (
        ("canonical source", canonical.get("source")),
        ("interval file", canonical.get("interval_file")),
        ("metadata file", canonical.get("metadata_file")),
        ("schema version", canonical.get("schema_version")),
        ("source epsilon", canonical.get("source_epsilon")),
        ("consumer epsilon", canonical.get("consumer_epsilon")),
        ("pre-particle-subtraction", canonical.get("source_stage")),
        ("requested / actual t bins", "{} / {}".format(canonical.get("requested_num_t_bins"), canonical.get("actual_num_t_bins"))),
        ("requested / actual phi bins", "{} / {}".format(canonical.get("requested_num_phi_bins"), canonical.get("actual_num_phi_bins"))),
        ("canonical pair ID", canonical.get("canonical_interval_pair_id")),
        ("canonical pair hash", canonical.get("canonical_interval_pair_hash")),
        ("support metric", canonical.get("support_metric")),
        ("raw counts", canonical.get("raw_event_count_by_t_bin")),
        ("signed weighted yields", canonical.get("signed_weighted_yield_by_t_bin")),
        ("absolute support", canonical.get("absolute_weighted_support_by_t_bin")),
        ("validation status", canonical.get("validation_status")),
        ("rejection reasons", canonical.get("validation_rejection_reasons")),
        ("t edges", canonical.get("t_edges")),
        ("phi edges", canonical.get("phi_edges")),
    ):
        provenance.AddText("{}: {}".format(label, value))
    provenance.Draw()
    canvas.Print(output_pdf)

    candidate_rows = list(diagnostics.get("timing_candidate_diagnostics") or [])
    canvas = TCanvas(
        "C_timing_t_candidate_summary_{}".format(page_id),
        "{} timing-t candidate summary".format(prefix), 1200, 800,
    )
    candidate_text = TPaveText(0.05, 0.06, 0.95, 0.94, "NDC")
    candidate_text.SetBorderSize(1)
    candidate_text.SetFillStyle(0)
    candidate_text.SetTextAlign(12)
    candidate_text.SetTextSize(0.025)
    candidate_text.AddText("Timing-t RF/CT candidates; selected={}".format(diagnostics.get("timing_branch")))
    candidate_text.AddText("rank tuple: {}".format(diagnostics.get("candidate_selection_tuple")))
    for row in candidate_rows:
        candidate_text.AddText(
            "{} ({}) range={} bins={} valid={} significance={} dev/entry={} support={} accepted={} rank={} selected={}".format(
                row.get("timing_branch"), row.get("probe_kind"), row.get("display_range"),
                row.get("histogram_bins"), row.get("global_fit_valid"),
                row.get("global_proton_significance"), row.get("poisson_deviance_per_entry"),
                row.get("setting_support_label"), row.get("setting_accepted"),
                row.get("candidate_selection_rank"), row.get("selected"),
            )
        )
        candidate_text.AddText("  supported/marginal/valid={} / {} / {}; reasons={}".format(
            row.get("supported_delta_count"), row.get("marginal_delta_count"),
            row.get("valid_delta_t_cell_count"), row.get("rejection_reasons"),
        ))
    setting_support = diagnostics.get("setting_support") or {}
    candidate_text.AddText("Setting support: {}".format(setting_support))
    candidate_text.Draw()
    canvas.Print(output_pdf)

    # Sections B/C: comparisons are global-only for rejected candidates;
    # detailed delta x t objects remain owned by the selected candidate.
    comparison_payloads = list(
        cleaning_result.get("_timing_t_candidate_global_comparisons") or []
    )
    selected_global = cleaning_result.get("H_global_timing")
    selected_vs_t = cleaning_result.get("H_global_timing_vs_t")
    selected_delta = cleaning_result.get("H_delta_timing")
    candidate_hists = []
    if selected_global is not None:
        candidate_hists.append(("selected {}".format(diagnostics.get("timing_branch", "?")), selected_global))
    for comparison in comparison_payloads:
        hist = comparison.get("H_global_timing")
        if hist is not None:
            candidate_hists.append(("rejected {}".format(comparison.get("timing_branch", "?")), hist))
    if candidate_hists:
        canvas = TCanvas(
            "C_timing_t_candidate_globals_{}".format(page_id),
            "{} timing-t candidate global comparisons".format(prefix), 1400, 900,
        )
        canvas.Divide(2, max(1, int(math.ceil(len(candidate_hists) / 2.0))))
        for index, (label, hist) in enumerate(candidate_hists, start=1):
            canvas.cd(index)
            hist.SetTitle("{};selected timing [ns];signed fit yield".format(label))
            hist.Draw("hist")
            gPad.Modified()
            gPad.Update()
        canvas.Print(output_pdf)
    if selected_global is not None or selected_vs_t is not None or selected_delta is not None:
        canvas = TCanvas(
            "C_timing_t_selected_model_{}".format(page_id),
            "{} selected timing-t model".format(prefix), 1500, 900,
        )
        canvas.Divide(3, 1)
        for index, (label, hist, option) in enumerate((
            ("Selected global timing", selected_global, "hist"),
            ("Selected timing versus |t|", selected_vs_t, "colz"),
            ("Selected delta timing", selected_delta, "colz"),
        ), start=1):
            canvas.cd(index)
            if hist is None:
                continue
            hist.SetTitle(label)
            hist.Draw(option)
            gPad.Modified()
            gPad.Update()
        canvas.Print(output_pdf)

    cell_rows = list(diagnostics.get("applied_timing_t_cell_map") or [])
    delta_edges = [float(edge) for edge in (cleaning_result.get("delta_edges") or [])]
    t_edges_for_cells = [float(edge) for edge in (cleaning_result.get("t_edges") or [])]
    if cell_rows and len(delta_edges) >= 2 and len(t_edges_for_cells) >= 2:
        def _delta_t_map(name, title):
            hist = ROOT.TH2D(
                "{}_{}".format(name, page_id), str(title),
                len(delta_edges) - 1, array("d", delta_edges),
                len(t_edges_for_cells) - 1, array("d", t_edges_for_cells),
            )
            hist.SetDirectory(0)
            return hist

        maps = {
            "fit_valid": _delta_t_map("H_delta_t_fit_valid", "Cell-fit validity;SHMS #delta [%];|t| [GeV^{2}]"),
            "proton_detected": _delta_t_map("H_delta_t_proton_detected", "Proton component detected;SHMS #delta [%];|t| [GeV^{2}]"),
            "delta_support": _delta_t_map("H_delta_t_delta_support", "Delta support class;SHMS #delta [%];|t| [GeV^{2}]"),
            "setting_support": _delta_t_map("H_delta_t_setting_support", "Setting support class;SHMS #delta [%];|t| [GeV^{2}]"),
            "raw_yield": _delta_t_map("H_delta_t_raw_proton_yield", "Raw fitted proton yield;SHMS #delta [%];|t| [GeV^{2}]"),
            "fitted_yield": _delta_t_map("H_delta_t_fitted_proton_yield", "Weak-component-handled fitted proton yield;SHMS #delta [%];|t| [GeV^{2}]"),
            "applied_yield": _delta_t_map("H_delta_t_applied_proton_yield", "Applied proton yield;SHMS #delta [%];|t| [GeV^{2}]"),
            "raw_minus_applied": _delta_t_map("H_delta_t_raw_minus_applied_yield", "Raw minus applied proton yield;SHMS #delta [%];|t| [GeV^{2}]"),
            "significance": _delta_t_map("H_delta_t_proton_significance", "Proton component significance;SHMS #delta [%];|t| [GeV^{2}]"),
            "applied_enabled": _delta_t_map("H_delta_t_applied_enabled", "Applied-cell enabled;SHMS #delta [%];|t| [GeV^{2}]"),
            "zero_reason": _delta_t_map("H_delta_t_applied_zero_reason", "Applied zero reason;SHMS #delta [%];|t| [GeV^{2}]"),
            "timing_center": _delta_t_map("H_delta_t_timing_center_source", "Timing-center source;SHMS #delta [%];|t| [GeV^{2}]"),
            "event_probability": _delta_t_map("H_delta_t_average_proton_probability", "Average frozen proton probability;SHMS #delta [%];|t| [GeV^{2}]"),
        }
        probability_sum = diagnostics.get("event_probability_sum_by_delta_t") or []
        probability_count = diagnostics.get("event_probability_count_by_delta_t") or []
        for row in cell_rows:
            delta_index, t_index = row.get("delta_index"), row.get("t_index")
            if not isinstance(delta_index, int) or not isinstance(t_index, int):
                continue
            if not (0 <= delta_index < maps["fit_valid"].GetNbinsX() and 0 <= t_index < maps["fit_valid"].GetNbinsY()):
                continue
            x_bin, y_bin = delta_index + 1, t_index + 1
            maps["fit_valid"].SetBinContent(x_bin, y_bin, float(bool(row.get("cell_fit_valid", row.get("valid", False)))))
            maps["proton_detected"].SetBinContent(x_bin, y_bin, float(bool(row.get("proton_component_detected", False))))
            maps["delta_support"].SetBinContent(x_bin, y_bin, SUPPORT_CLASS_TO_CODE.get(str(row.get("delta_support_label", SUPPORT_UNSUPPORTED)), 0.0))
            maps["setting_support"].SetBinContent(x_bin, y_bin, SUPPORT_CLASS_TO_CODE.get(str(row.get("setting_support_label", SUPPORT_UNSUPPORTED)), 0.0))
            raw_yield = float(row.get("raw_proton_yield", 0.0) or 0.0)
            fitted_yield = float(row.get("fitted_proton_yield", row.get("proton_yield", 0.0)) or 0.0)
            applied_yield = float(row.get("applied_proton_yield", 0.0) or 0.0)
            maps["raw_yield"].SetBinContent(x_bin, y_bin, raw_yield)
            maps["fitted_yield"].SetBinContent(x_bin, y_bin, fitted_yield)
            maps["applied_yield"].SetBinContent(x_bin, y_bin, applied_yield)
            maps["raw_minus_applied"].SetBinContent(x_bin, y_bin, raw_yield - applied_yield)
            maps["significance"].SetBinContent(x_bin, y_bin, float(row.get("proton_component_significance", 0.0) or 0.0))
            maps["applied_enabled"].SetBinContent(x_bin, y_bin, float(bool(row.get("applied_cell_enabled", False))))
            reason = "applied" if bool(row.get("applied_cell_enabled", False)) else str(row.get("applied_zero_reason") or "invalid_cell_fit")
            maps["zero_reason"].SetBinContent(x_bin, y_bin, APPLIED_ZERO_REASON_TO_CODE.get(reason, APPLIED_ZERO_REASON_TO_CODE["invalid_cell_fit"]))
            maps["timing_center"].SetBinContent(x_bin, y_bin, TIMING_CENTER_SOURCE_TO_CODE.get(str(row.get("timing_center_source") or "invalid_timing_center"), TIMING_CENTER_SOURCE_TO_CODE["invalid_timing_center"]))
            try:
                count = float(probability_count[delta_index][t_index])
                value = float(probability_sum[delta_index][t_index])
                if count > 0.0:
                    maps["event_probability"].SetBinContent(x_bin, y_bin, value / count)
            except (IndexError, TypeError, ValueError):
                pass

        canvas = TCanvas("C_timing_t_delta_t_yields_{}".format(page_id), "{} raw/fitted/applied delta-t yields".format(prefix), 1800, 700)
        canvas.Divide(4, 1)
        for index, key in enumerate(("raw_yield", "fitted_yield", "applied_yield", "raw_minus_applied"), start=1):
            canvas.cd(index)
            maps[key].Draw("colz text")
        canvas.Print(output_pdf)
        canvas = TCanvas("C_timing_t_delta_t_categories_{}".format(page_id), "{} delta-t categorical state".format(prefix), 1800, 900)
        canvas.Divide(3, 2)
        for index, key in enumerate(("fit_valid", "proton_detected", "delta_support", "setting_support", "applied_enabled", "zero_reason"), start=1):
            canvas.cd(index)
            maps[key].Draw("colz text")
            if key == "delta_support":
                _draw_categorical_legend({label: int(code) for label, code in SUPPORT_CLASS_TO_CODE.items()}, "support")
            elif key == "zero_reason":
                _draw_categorical_legend(APPLIED_ZERO_REASON_TO_CODE, "zero reason")
        canvas.Print(output_pdf)
        canvas = TCanvas("C_timing_t_delta_t_timing_{}".format(page_id), "{} delta-t timing and event state".format(prefix), 1500, 700)
        canvas.Divide(3, 1)
        for index, key in enumerate(("significance", "timing_center", "event_probability"), start=1):
            canvas.cd(index)
            maps[key].Draw("colz text")
            if key == "timing_center":
                _draw_categorical_legend(TIMING_CENTER_SOURCE_TO_CODE, "timing center")
        canvas.Print(output_pdf)

        h_applied = ROOT.TH2D(
            "H_timing_t_applied_cell_state_{}".format(page_id),
            "Applied timing-t cell state;SHMS #delta [%];|t| [GeV^{2}]",
            len(delta_edges) - 1, array("d", delta_edges),
            len(t_edges_for_cells) - 1, array("d", t_edges_for_cells),
        )
        h_closure = _clone_hist(
            h_applied, "H_timing_t_applied_cell_closure_{}".format(page_id), reset=True,
        )
        for row in cell_rows:
            delta_index, t_index = row.get("delta_index"), row.get("t_index")
            if not isinstance(delta_index, int) or not isinstance(t_index, int):
                continue
            if not (0 <= delta_index < h_applied.GetNbinsX() and 0 <= t_index < h_applied.GetNbinsY()):
                continue
            h_applied.SetBinContent(
                delta_index + 1, t_index + 1,
                1.0 if bool(row.get("applied_cell_enabled", False)) else 0.0,
            )
            closure_ratio = _finite_float_or_none(row.get("closure_ratio"))
            if closure_ratio is not None:
                h_closure.SetBinContent(delta_index + 1, t_index + 1, closure_ratio)
        canvas = TCanvas(
            "C_timing_t_applied_cell_maps_{}".format(page_id),
            "{} timing-t applied cell and closure maps".format(prefix), 1400, 700,
        )
        canvas.Divide(2, 1)
        canvas.cd(1)
        h_applied.Draw("colz text")
        canvas.cd(2)
        h_closure.SetTitle("Applied-cell closure ratio;SHMS #delta [%];|t| [GeV^{2}]")
        h_closure.Draw("colz text")
        canvas.Print(output_pdf)

    # Section F is emitted after the observational aerogel pages.  Retain the
    # existing calculations here for compatibility, but do not render them
    # before Section E.
    defer_closure_pages = True
    rows = list(diagnostics.get("cross_stage_t_consistency") or [])
    canvas = TCanvas(
        "C_timing_t_cross_stage_{}".format(page_id),
        "{} cross-stage shifted-t consistency".format(prefix), 1100, 750,
    )
    h_difference = TH1D(
        "H_timing_t_cross_stage_difference_{}".format(page_id),
        "Cross-stage shifted-t maximum difference;sample index;max |#Delta t|",
        max(len(rows), 1), 0.0, float(max(len(rows), 1)),
    )
    h_difference.SetDirectory(0)
    for index, row in enumerate(rows, start=1):
        h_difference.SetBinContent(index, float(row.get("maximum_absolute_difference", 0.0) or 0.0))
    h_difference.Draw("hist")
    consistency = diagnostics.get("cross_stage_t_consistency_summary") or {}
    text = TPaveText(0.50, 0.70, 0.90, 0.90, "NDC")
    text.SetBorderSize(1)
    text.SetFillStyle(0)
    text.AddText("prepass / proton / downstream")
    text.AddText("sampling: {}".format(consistency.get("sampling_method", "implementation_consistency_sample")))
    text.AddText("samples: {}".format(consistency.get("sample_count", len(rows))))
    text.AddText("completed downstream: {}".format(consistency.get("completed_downstream_count", 0)))
    text.AddText("max difference: {}".format(consistency.get("maximum_absolute_difference", 0.0)))
    text.AddText("tolerance: {}".format(consistency.get("tolerance", "n/a")))
    text.Draw()
    if not defer_closure_pages:
        canvas.Print(output_pdf)

    # Frozen lookup coverage in the final canonical (t, phi) coordinates.
    phi_edges_for_lookup = [float(edge) for edge in (canonical.get("phi_edges") or [])]
    t_phi_rows = list(diagnostics.get("event_weight_lookup_by_t_phi") or [])
    if len(t_edges_for_cells) >= 2 and len(phi_edges_for_lookup) >= 2 and t_phi_rows:
        h_tphi_count = ROOT.TH2D(
            "H_timing_t_lookup_tphi_count_{}".format(page_id),
            "Frozen lookup event count;|t| [GeV^{2}];#phi [deg]",
            len(t_edges_for_cells) - 1, array("d", t_edges_for_cells),
            len(phi_edges_for_lookup) - 1, array("d", phi_edges_for_lookup),
        )
        h_tphi_proton = _clone_hist(
            h_tphi_count, "H_timing_t_lookup_tphi_proton_{}".format(page_id), reset=True,
        )
        for row in t_phi_rows:
            t_index, phi_index = row.get("t_index"), row.get("phi_index")
            if isinstance(t_index, int) and isinstance(phi_index, int):
                if 0 <= t_index < h_tphi_count.GetNbinsX() and 0 <= phi_index < h_tphi_count.GetNbinsY():
                    h_tphi_count.SetBinContent(t_index + 1, phi_index + 1, float(row.get("event_count", 0) or 0))
                    h_tphi_proton.SetBinContent(t_index + 1, phi_index + 1, float(row.get("estimated_proton_yield", 0.0) or 0.0))
        canvas = TCanvas(
            "C_timing_t_lookup_tphi_{}".format(page_id),
            "{} frozen lookup by |t| and phi".format(prefix), 1400, 700,
        )
        canvas.Divide(2, 1)
        canvas.cd(1)
        h_tphi_count.Draw("colz text")
        canvas.cd(2)
        h_tphi_proton.SetTitle("Frozen lookup estimated proton yield;|t| [GeV^{2}];#phi [deg]")
        h_tphi_proton.Draw("colz text")
        if not defer_closure_pages:
            canvas.Print(output_pdf)

    aero = diagnostics.get("aerogel_vs_t_validation") or diagnostics.get("aerogel_validation") or {}
    t_edges = list(aero.get("t_edges") or [])
    aero_edges = list(aero.get("aero_edges") or [])
    if len(t_edges) < 2 or len(aero_edges) < 2:
        # Aerogel products may be disabled, but frozen-lookup closure remains
        # part of the timing-t report.
        _print_timing_t_closure_pages(output_pdf, cleaning_result, prefix, page_id)
        _print_timing_t_mm_diagnostic_pages(output_pdf, cleaning_result, prefix, page_id)
        return
    validation_cfg = dict(aero.get("configuration") or {})

    lambda_gate_status = str(
        (diagnostics.get("lambda_preservation_gate") or {}).get("status") or "not evaluated"
    ).upper()
    observational_header = (
        "Aerogel secondary PID validation only - proposed timing-|t| model; not used in production event weights"
    )

    # Section E1: exact canonical-t x fine-aerogel frozen-lookup maps.
    root_payload = cleaning_result.get("_aerogel_vs_t_root_payload") or {}
    global_hists = root_payload.get("global") or {}
    if bool(validation_cfg.get("write_global_aero_vs_t_pages", True)) and global_hists:
        canvas, body, _header, _label = _make_timing_t_report_canvas(
            "C_timing_t_aero_global_{}".format(page_id),
            "{} aerogel versus |t| validation".format(prefix), 1600, 1000,
            header_text=observational_header, panel_count=10,
        )
        for index, key in enumerate((
            "H_aero_vs_t_raw_prompt", "H_aero_vs_t_signed_physical",
            "H_aero_vs_t_signed_physical_positive", "H_aero_vs_t_signed_physical_negative",
            "H_aero_vs_t_estimated_proton", "H_aero_vs_t_estimated_proton_positive",
            "H_aero_vs_t_estimated_proton_negative_abs", "H_aero_vs_t_proton_cleaned",
            "H_aero_vs_t_proton_cleaned_positive", "H_aero_vs_t_proton_cleaned_negative_abs",
        ), start=1):
            body.cd(index)
            hist = global_hists.get(key)
            if hist is not None:
                hist.Draw("colz")
                gPad.Modified()
                gPad.Update()
        canvas.Print(output_pdf)

    # Section E2: eight compact canonical-t x summary-aerogel matrices.
    if bool(validation_cfg.get("write_t_aero_heatmaps", True)):
        matrix_payload = aero.get("matrix_payload") or {}
        coarse_hists = root_payload.get("coarse") or {}
        matrix_specs = (
            ("selected_prompt_count", "H_t_aero_selected_prompt_count", "Selected prompt count"),
            ("signed_physical_yield", "H_t_aero_signed_physical_yield", "Signed physical yield"),
            ("estimated_proton_yield", "H_t_aero_estimated_proton_yield", "Estimated proton yield"),
            ("average_proton_probability", "H_t_aero_average_proton_probability", "Average proton probability"),
            ("low_mm_removed_yield", "H_t_aero_low_mm_removed_yield", "Low-MM removed yield"),
            ("low_mm_removed_fraction", "H_t_aero_low_mm_removed_fraction", "Low-MM removal fraction"),
            ("lambda_removed_yield", "H_t_aero_lambda_removed_yield", "Lambda removed yield"),
            ("lambda_removed_fraction", "H_t_aero_lambda_removed_fraction", "Lambda removal fraction"),
        )
        for page_index, specs in enumerate((matrix_specs[:4], matrix_specs[4:]), start=1):
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_aero_matrices_{}_{}".format(page_id, page_index),
                "{} aerogel validation summary matrices".format(prefix), 1500, 1000,
                header_text=observational_header, panel_count=len(specs),
            )
            for draw_index, (metric_key, root_key, metric_label) in enumerate(specs, start=1):
                body.cd(draw_index)
                hist = coarse_hists.get(root_key)
                if hist is None:
                    hist = _make_t_aerogel_hist(
                        "H_timing_t_aero_{}_{}_{}".format(metric_key, page_id, page_index),
                        "{};|t| [GeV^{{2}}];P_aero NPE".format(metric_label),
                        t_edges, aero_edges,
                    )
                    for t_index, value_row in enumerate(_matrix_metric(matrix_payload, metric_key)):
                        for aero_index, value in enumerate(value_row or []):
                            numeric = _finite_float_or_none(value)
                            if numeric is not None and t_index < hist.GetNbinsX() and aero_index < hist.GetNbinsY():
                                hist.SetBinContent(t_index + 1, aero_index + 1, numeric)
                status = _timing_t_root_content_status(hist)
                if status["state"] == "empty":
                    _draw_timing_t_status_panel(
                        "{}\nNo populated coarse cells; see diagnostic integrity".format(metric_label),
                    )
                else:
                    hist.Draw("colz text")
                    gPad.Modified()
                    gPad.Update()
            canvas.Print(output_pdf)
        if coarse_hists:
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_aero_validity_masks_{}".format(page_id),
                "{} aerogel value and validity masks".format(prefix), 1600, 900, 3, 2,
                observational_header + " | mask: 1=evaluated, 0=undefined",
            )
            for draw_index, key in enumerate((
                "H_t_aero_average_proton_probability", "H_t_aero_average_proton_probability_valid",
                "H_t_aero_low_mm_removed_fraction", "H_t_aero_low_mm_removed_fraction_valid",
                "H_t_aero_lambda_removed_fraction", "H_t_aero_lambda_removed_fraction_valid",
            ), start=1):
                body.cd(draw_index)
                hist = coarse_hists.get(key)
                if hist is not None:
                    hist.Draw("colz text")
            canvas.Print(output_pdf)

    # Section E3: compact, per-canonical-t PID views created from the same
    # frozen rows.  Hidden numerator/support histograms never enter a page.
    if bool(validation_cfg.get("write_per_t_pid_pages", True)):
        skipped_t_summaries = []
        minimum_absolute_support = abs(float(
            validation_cfg.get("per_t_absolute_support_tolerance", 0.0) or 0.0
        ))
        for row in root_payload.get("per_t") or []:
            summary = (aero.get("per_t_bin_summary") or [])
            t_index = int(row.get("t_index", -1))
            t_summary = summary[t_index] if 0 <= t_index < len(summary) else {}
            eligibility = _timing_t_per_t_pid_eligible(t_summary, validation_cfg)
            prompt_count = eligibility["prompt_count"]
            absolute_support = eligibility["absolute_support"]
            if not eligibility["eligible"]:
                skipped_t_summaries.append({
                    "t_index": t_index,
                    "prompt_count": prompt_count,
                    "absolute_support": absolute_support,
                    "warnings": list(t_summary.get("warnings") or []),
                })
                continue
            context = diagnostics.get("kinematic_context") or {}
            local_cells = [
                cell for cell in (diagnostics.get("applied_timing_t_cell_map") or [])
                if cell.get("t_index") == t_index
            ]
            local_support = "/".join(
                "{}:{}".format(label, sum(str(cell.get("delta_support_label")) == label for cell in local_cells))
                for label in SUPPORT_CLASS_ORDER
            )
            applied_yield = sum(float(cell.get("applied_proton_yield", 0.0) or 0.0) for cell in local_cells)
            header_text = [
                "Aerogel secondary PID validation only - proposed timing-|t| model; Lambda gate={}".format(lambda_gate_status),
                "Q2={q2} GeV^2  W={w} GeV  eps={epsilon}  phi={phi} | t bin {t}: [{low:.6g}, {high:.6g}] GeV^2 | timing={branch} {timing}".format(
                q2=context.get("Q2", "n/a"), w=context.get("W", "n/a"), epsilon=context.get("epsilon", "n/a"),
                phi=context.get("phi_setting", cleaning_result.get("phi_setting", "n/a")),
                t=t_index, low=float(t_summary.get("t_low", 0.0) or 0.0), high=float(t_summary.get("t_high", 0.0) or 0.0),
                branch=diagnostics.get("timing_branch", "n/a"), timing=root_payload.get("timing_range", "n/a"),
                ),
                "physical source units: prompt={prompt} signed yield={signed:.4g} estimated proton={proton:.4g} abs support={support:.4g}".format(
                    prompt=prompt_count,
                    signed=float(t_summary.get("signed_event_weight_sum", 0.0) or 0.0),
                    proton=float(t_summary.get("estimated_proton_missing_mass_yield", 0.0) or 0.0),
                    support=absolute_support,
                ),
                "fit-model units (not compared to physical yield): applied proton yield={applied:.4g}; delta support={delta_support}; setting={setting}; warnings={warnings}".format(
                    applied=applied_yield, delta_support=local_support,
                    setting=(diagnostics.get("setting_support") or {}).get("support_label", "n/a"),
                    warnings=",".join(t_summary.get("warnings") or []) or "none",
                ),
            ]
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_aero_per_t_{}_{}".format(page_id, row.get("t_index")),
                "{} aerogel PID validation t bin {}".format(prefix, row.get("t_index")), 1500, 1000,
                header_text=header_text, panel_count=4,
            )
            body.cd(1)
            row["raw_prompt_timing_vs_aero"].Draw("colz")
            body.cd(2)
            row["estimated_proton_timing_vs_aero"].Draw("colz")
            body.cd(3)
            raw_projection = row["raw_signed_projection"]
            proton_projection = row["estimated_proton_projection"]
            cleaned_projection = row["cleaned_projection"]
            _set_hist_line_marker(raw_projection, kBlack, width=2)
            _set_hist_line_marker(proton_projection, kRed, width=3)
            _set_hist_line_marker(cleaned_projection, kGreen + 2, width=3)
            raw_projection.SetTitle("Physical aerogel projections;P_aero NPE;signed yield")
            raw_projection.Draw("hist")
            proton_projection.Draw("hist same")
            cleaned_projection.Draw("hist same")
            legend = TLegend(0.50, 0.68, 0.88, 0.88)
            legend.AddEntry(raw_projection, "raw signed", "l")
            legend.AddEntry(proton_projection, "estimated proton", "l")
            legend.AddEntry(cleaned_projection, "frozen cleaned", "l")
            legend.Draw()
            body.cd(4)
            row["average_proton_probability"].SetMinimum(0.0)
            row["average_proton_probability"].SetMaximum(1.0)
            row["average_proton_probability"].Draw("hist")
            canvas.Print(output_pdf)
            if bool(validation_cfg.get("write_full_per_t_pid_pages", False)) and "full_raw_prompt_timing_vs_aero" in row:
                full_canvas, full_body, _header, _label = _make_timing_t_report_canvas(
                    "C_timing_t_aero_per_t_full_{}_{}".format(page_id, row.get("t_index")),
                    "{} full aerogel PID detail t bin {}".format(prefix, row.get("t_index")), 1500, 1000, 2, 2,
                    header_text,
                )
                for draw_index, key in enumerate((
                    "full_raw_prompt_timing_vs_aero", "full_signed_physical_timing_vs_aero",
                    "full_estimated_proton_timing_vs_aero", "full_cleaned_timing_vs_aero",
                ), start=1):
                    full_body.cd(draw_index)
                    row[key].Draw("colz")
                full_canvas.Print(output_pdf)
        if skipped_t_summaries:
            canvas, body, _header, _label = _make_timing_t_report_canvas(
                "C_timing_t_aero_skipped_t_bins_{}".format(page_id),
                "{} skipped aerogel PID bins".format(prefix), 1300, 750,
                header_text=observational_header, panel_count=1,
            )
            body.cd()
            message = [
                "No compact PID page was rendered for the following |t| bins.",
                "Requirement: raw prompt count > 0 and absolute physical support > {:.6g}.".format(minimum_absolute_support),
            ]
            for skipped in skipped_t_summaries:
                message.append(
                    "t{t}: prompt={prompt_count}, abs support={absolute_support:.6g}, warnings={warnings}".format(
                        **{**skipped, "warnings": ",".join(skipped["warnings"]) or "none"}
                    )
                )
            _draw_timing_t_status_panel("\n".join(message), color=kBlack)
            canvas.Print(output_pdf)

    warning_rows = [row for row in (aero.get("warnings_by_t_bin") or []) if row.get("warnings")]
    integrity = diagnostics.get("timing_t_diagnostic_integrity") or aero.get("diagnostic_integrity") or {}
    per_t_details = {
        int(row.get("t_index")): row
        for row in (aero.get("per_t_bin_summary") or [])
        if isinstance(row, dict) and isinstance(row.get("t_index"), int)
    }
    canvas, body, _header, _label = _make_timing_t_report_canvas(
        "C_timing_t_aero_warnings_{}".format(page_id),
        "{} aerogel validation exceptions".format(prefix), 1400, 800,
        header_text=observational_header, panel_count=1,
    )
    warning_text = TPaveText(0.05, 0.05, 0.95, 0.95, "NDC")
    warning_text.SetBorderSize(1)
    warning_text.SetFillStyle(0)
    warning_text.SetTextAlign(12)
    warning_text.SetTextSize(0.028)
    warning_text.AddText("Warnings only; complete diagnostics are stored in JSON/CSV")
    if not bool(integrity.get("valid", True)):
        warning_text.SetTextColor(kRed)
        warning_text.AddText(
            "DIAGNOSTIC INTEGRITY WARNING: {}".format(
                ", ".join(integrity.get("failures") or []) or "unknown failure"
            )
        )
    if not warning_rows:
        warning_text.AddText("No per-t-bin aerogel validation warnings")
    for row in warning_rows:
        detail = per_t_details.get(int(row.get("t_index", -1)), {})
        warning_text.AddText(
            "t{} [{:.6g}, {:.6g}]: {} | low/high Pp={}/{} | low/high low-MM removal={}/{} | low/high Lambda removal={}/{} | abs support={:.5g}".format(
                row.get("t_index"), float(row.get("t_low", 0.0) or 0.0),
                float(row.get("t_high", 0.0) or 0.0), ", ".join(row.get("warnings") or []),
                detail.get("low_aero_average_proton_probability"),
                detail.get("high_aero_average_proton_probability"),
                detail.get("low_aero_low_mm_removed_fraction"),
                detail.get("high_aero_low_mm_removed_fraction"),
                detail.get("low_aero_lambda_removed_fraction"),
                detail.get("high_aero_lambda_removed_fraction"),
                float(detail.get("absolute_event_weight_support", 0.0) or 0.0),
            )
        )
    body.cd()
    warning_text.Draw()
    canvas.Print(output_pdf)

    # Section F: frozen lookup/cross-stage closure follows all aerogel-only
    # pages.  Section G then reports the frozen MM physics effect before the
    # selected-candidate-only final summary.
    _print_timing_t_closure_pages(output_pdf, cleaning_result, prefix, page_id)
    _print_timing_t_mm_diagnostic_pages(output_pdf, cleaning_result, prefix, page_id)


def _print_timing_t_final_summary_page(output_pdf, cleaning_result, prefix):
    """Render Section J without falling back to legacy aerogel diagnostics."""
    diagnostics = cleaning_result.get("diagnostics") or {}
    summary = diagnostics.get("timing_t_summary") or {}
    page_id = abs(id(cleaning_result))
    canvas, body, _header, _label = _make_timing_t_report_canvas(
        "C_timing_t_selected_summary_{}".format(page_id),
        "{} selected timing-t summary".format(prefix), 1450, 900,
        header_text="Timing-t final summary - selected candidate and frozen lookup only", panel_count=1,
    )
    body.cd()
    summary_text = TPaveText(0.04, 0.04, 0.96, 0.96, "NDC")
    summary_text.SetBorderSize(1)
    summary_text.SetFillStyle(0)
    summary_text.SetTextAlign(12)
    summary_text.SetTextSize(0.025)
    if not summary:
        summary_text.SetTextColor(kRed)
        summary_text.AddText("Timing-t summary unavailable; see diagnostic integrity payload")
    else:
        candidate = summary.get("selected_candidate") or {}
        setting = summary.get("setting_support") or {}
        mm = summary.get("missing_mass_totals") or {}
        gate = summary.get("lambda_preservation_gate") or {}
        per_delta = list(summary.get("per_delta") or [])
        applied = sum(float(row.get("applied_proton_yield", 0.0) or 0.0) for row in per_delta)
        data_total = sum(float(row.get("raw_data_total", 0.0) or 0.0) for row in per_delta)
        summary_text.AddText(
            "selected timing branch={} | rank={} | setting support={} accepted={}".format(
                candidate.get("timing_branch"), candidate.get("candidate_selection_rank"),
                setting.get("support_label"), setting.get("accepted"),
            )
        )
        summary_text.AddText("delta edges: {}".format(summary.get("delta_edges")))
        summary_text.AddText(
            "integrated applied proton fraction={} (applied yield={} / raw fit data={})".format(
                _safe_validation_ratio(applied, data_total), applied, data_total,
            )
        )
        proposed_mm = mm.get("proposed") or {}
        applied_mm = mm.get("applied") or {}
        summary_text.AddText(
            "proposed MM totals: raw={} estimated proton={} cleaned pre-RF={}".format(
                proposed_mm.get("raw", mm.get("raw")),
                proposed_mm.get("estimated_proton", mm.get("estimated_proton")),
                proposed_mm.get("cleaned_pre_rf", mm.get("cleaned")),
            )
        )
        summary_text.AddText(
            "final applied MM totals: raw={} estimated proton={} cleaned pre-RF={} cleaned post-RF={}".format(
                applied_mm.get("raw"), applied_mm.get("estimated_proton"),
                applied_mm.get("cleaned_pre_rf"), applied_mm.get("cleaned_final_rf"),
            )
        )
        summary_text.AddText(
            "Lambda gate={} action={} committed={}".format(
                gate.get("status"), gate.get("production_action"), gate.get("proton_cleaning_committed"),
            )
        )
        for row in per_delta:
            summary_text.AddText(
                "delta {delta_index} [{delta_low}, {delta_high}]: support={support_label} "
                "data={raw_data_total:.5g} fitted={fitted_data_total:.5g} coverage={coverage:.4g} "
                "K/p/other={kaon_yield:.5g}/{proton_yield:.5g}/{other_yield:.5g} "
                "applied={applied_proton_yield:.5g} mean lookup={mean_event_lookup_probability}".format(
                    **row
                )
            )
    summary_text.Draw()
    canvas.Print(output_pdf)


def _print_timing_t_lambda_preservation_gate_pages(output_pdf, cleaning_result, prefix):
    """Render the setting-wide gate immediately before pion diagnostics."""
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    gate = diagnostics.get("lambda_preservation_gate") or {}
    if not gate:
        return
    page_id = abs(id(cleaning_result))
    status = str(gate.get("status") or "insufficient_support").upper()
    status_color = kGreen + 2 if status == "PASS" else (kRed if status == "FAIL" else kOrange + 7)
    canvas, body, _header, _label = _make_timing_t_report_canvas(
        "C_timing_t_lambda_gate_{}".format(page_id),
        "{} Lambda-preservation gate".format(prefix), 1500, 900,
        header_text="Setting-wide Lambda-preservation gate - pre-RF proposed timing-|t| lookup",
        panel_count=1,
    )
    body.cd()
    banner = TPaveText(0.05, 0.79, 0.95, 0.94, "NDC")
    banner.SetBorderSize(2)
    banner.SetFillColor(status_color)
    banner.SetTextAlign(22)
    banner.SetTextSize(0.045)
    banner.AddText("{}: {} PROTON CLEANING".format(
        status,
        "APPLY" if str(gate.get("production_action")) == "apply" else "BYPASS",
    ))
    banner.Draw()
    text = TPaveText(0.06, 0.06, 0.94, 0.74, "NDC")
    text.SetBorderSize(1)
    text.SetFillStyle(0)
    text.SetTextAlign(12)
    text.SetTextSize(0.028)
    text.AddText(
        "timing model accepted={} | branch={} | setting support={}".format(
            gate.get("timing_fit_accepted"), diagnostics.get("timing_branch"),
            gate.get("setting_support_label"),
        )
    )
    text.AddText(
        "Lambda window {}: [{:.6g}, {:.6g}] GeV (shared validation window)".format(
            gate.get("validation_window_key"),
            float(gate.get("window_low", 0.0) or 0.0),
            float(gate.get("window_high", 0.0) or 0.0),
        )
    )
    text.AddText(
        "raw prompt count={} | raw signed yield={:.7g} | absolute support={:.7g} | signed/abs={}".format(
            gate.get("raw_prompt_count"),
            float(gate.get("raw_signed_yield", 0.0) or 0.0),
            float(gate.get("raw_absolute_support", 0.0) or 0.0),
            gate.get("raw_signed_to_absolute_support_ratio"),
        )
    )
    text.AddText(
        "proposed proton Lambda yield={:.7g} | proposed removed fraction={} | maximum allowed={:.6g}".format(
            float(gate.get("proposed_proton_yield", 0.0) or 0.0),
            gate.get("proposed_removed_fraction"),
            float(gate.get("maximum_removed_fraction", 0.10) or 0.10),
        )
    )
    text.AddText("support valid={} | reasons={}".format(
        gate.get("support_valid"), gate.get("support_reasons") or [],
    ))
    text.AddText("observational warnings={}".format(gate.get("observational_warnings") or []))
    text.AddText("proton cleaning committed={} | applied lookup zeroed={}".format(
        gate.get("proton_cleaning_committed"), gate.get("applied_lookup_zeroed"),
    ))
    final_closure = gate.get("final_applied_closure") or {}
    text.AddText(
        "final applied pre-RF closure: raw - p - cleaned = {:.5g}".format(
            float(final_closure.get("pre_rf_closure_difference", 0.0) or 0.0)
        )
    )
    text.Draw()
    canvas.Print(output_pdf)

    root_payload = (cleaning_result or {}).get("_timing_t_mm_root_payload") or {}
    maps = root_payload.get("maps") or {}
    ordered = (
        ("raw", "raw", kBlack),
        ("proposed_estimated_proton", "proposed estimated proton", kRed),
        ("proposed_cleaned_pre_rf", "proposed cleaned pre-RF", kGreen + 2),
        ("applied_cleaned_pre_rf", "final applied cleaned pre-RF", kBlue),
    )
    projections = []
    for key, label, color in ordered:
        hist = maps.get(key)
        if hist is None:
            continue
        projection = hist.ProjectionX(
            "H_timing_t_lambda_gate_{}_{}".format(key, page_id), 1, hist.GetNbinsY(), "e"
        )
        projection.SetDirectory(0)
        _set_hist_line_marker(projection, color, width=3)
        projections.append((label, projection))
    if not projections:
        return
    canvas, body, _header, _label = _make_timing_t_report_canvas(
        "C_timing_t_lambda_gate_mm_{}".format(page_id),
        "{} proposed versus applied MM".format(prefix), 1400, 900,
        header_text=(
            "Lambda gate {} - proposed model is diagnostic; final applied curve is what pion subtraction receives".format(status)
        ),
        panel_count=1,
    )
    body.cd()
    maximum = max(float(hist.GetMaximum()) for _label, hist in projections)
    minimum = min(float(hist.GetMinimum()) for _label, hist in projections)
    y_min, y_max = _timing_t_signed_display_limits(minimum, maximum)
    first = projections[0][1]
    first.SetTitle("Setting-wide shifted MM;shifted MM [GeV];signed physical yield")
    first.SetMinimum(y_min)
    first.SetMaximum(y_max)
    first.Draw("hist")
    for _label, projection in projections[1:]:
        projection.Draw("hist same")
    legend = TLegend(0.52, 0.62, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for label, projection in projections:
        legend.AddEntry(projection, label, "l")
    legend.Draw()
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

    # The opt-in timing-t route has its own ordered diagnostics report.  The
    # C-macro-reference aerogel route continues through the unchanged pages.
    if str(cleaning_result.get("method") or "") == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT:
        _print_timing_t_validation_pages(output_pdf, cleaning_result, prefix)
        _print_timing_t_final_summary_page(output_pdf, cleaning_result, prefix)
        _print_timing_t_lambda_preservation_gate_pages(output_pdf, cleaning_result, prefix)
        return

    _print_timing_probe_comparison_page(output_pdf, cleaning_result, prefix)
    _print_proton_tof_constraint_diagnostics_page(output_pdf, cleaning_result, prefix)
    _print_low_aero_offset_diagnostics_page(output_pdf, cleaning_result, prefix)
    _print_proton_aerogel_diagnostics_page(output_pdf, cleaning_result, prefix)
    _print_timing_t_validation_pages(output_pdf, cleaning_result, prefix)

    h_global_pid = cleaning_result.get("H_global_pid")
    if h_global_pid is not None:
        canvas = TCanvas("C_proton_cleaning_global_pid", "{} proton-cleaning global PID".format(prefix), 1000, 700)
        drawn_objects = []
        h_global_pid.Draw("colz")
        drawn_objects.extend(
            _draw_aerogel_reference_line(
                h_global_pid,
                float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0),
                "x",
            )
        )
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
    applied_map_key = (
        "H_proton_weight_vs_delta_t"
        if str(cleaning_result.get("method") or "") == PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT
        else "H_proton_weight_vs_delta_aero"
    )
    if application.get(applied_map_key) is not None:
        application[applied_map_key].Draw("colz")
        if applied_map_key.endswith("_aero"):
            drawn_objects.extend(
                _draw_aerogel_reference_line(
                    application[applied_map_key],
                    float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0),
                    "y",
                )
            )
        drawn_objects.append(application[applied_map_key])
    else:
        first_delta_pid = next((hist for hist in (cleaning_result.get("H_delta_pid") or []) if hist is not None), None)
        if first_delta_pid is not None:
            first_delta_pid.Draw("colz")
            drawn_objects.extend(
                _draw_aerogel_reference_line(
                    first_delta_pid,
                    float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0),
                    "x",
                )
            )
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
            drawn_objects.extend(
                _draw_aerogel_reference_line(
                    pid_hist,
                    float(diagnostics.get("aerogel_reference_npe", 5.0) or 5.0),
                    "x",
                )
            )
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
                        "support entries={}".format(slice_fit.get("support_entries", "n/a")),
                        "timing center={}".format(slice_fit.get("timing_center_source") or "unavailable"),
                        "offset applied/valid={}/{}".format(
                            bool(slice_fit.get("offset_refinement_applied", False)),
                            bool(slice_fit.get("offset_refinement_valid", False)),
                        ),
                        "model/data={}".format(slice_fit.get("model_data_ratio")),
                        "chi2/ndf={}".format(slice_fit.get("chi2_ndf")),
                        "warn={}".format(
                            ", ".join(str(w) for w in (slice_fit.get("diagnostic_warnings") or [])) or "none"
                        ),
                        "reason={}".format(slice_fit.get("rejection_reason") or "none"),
                    ),
                    x1=0.16,
                    y1=0.54,
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
