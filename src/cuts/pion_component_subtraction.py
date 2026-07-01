#! /usr/bin/python

from __future__ import annotations

import math
from copy import deepcopy

import numpy as np

from background_config import (
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    resolve_particle_subtraction_fallback_mode,
)
from mm_background_subtraction import mm_background_weight_from_value


COMPONENT_NAMES = ("pi_n", "pi_delta", "pi_sidis")


def _is_root_object(obj):
    try:
        return bool(obj is not None and obj.InheritsFrom("TObject"))
    except Exception:
        return False


def _is_root_hist(obj):
    try:
        return bool(obj is not None and obj.InheritsFrom("TH1"))
    except Exception:
        return False


_JSON_SKIP = object()


def _json_ready_particle_subtraction_value(value):
    if _is_root_object(value):
        return _JSON_SKIP
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        cleaned = {}
        for key, child_value in value.items():
            child = _json_ready_particle_subtraction_value(child_value)
            if child is _JSON_SKIP:
                continue
            cleaned[key] = child
        return cleaned
    if isinstance(value, (list, tuple)):
        cleaned = []
        for child_value in value:
            child = _json_ready_particle_subtraction_value(child_value)
            if child is _JSON_SKIP:
                continue
            cleaned.append(child)
        return cleaned
    if isinstance(value, set):
        cleaned = []
        for child_value in value:
            child = _json_ready_particle_subtraction_value(child_value)
            if child is _JSON_SKIP:
                continue
            cleaned.append(child)
        return cleaned
    return deepcopy(value)


def _clone_hist(template_hist, name, reset=False):
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


def _hist_bin_signature(hist):
    if hist is None:
        return None
    axis = hist.GetXaxis()
    return (
        int(hist.GetNbinsX()),
        float(axis.GetXmin()),
        float(axis.GetXmax()),
    )


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def _bin_index_mm_range(hist, bin_indices):
    if hist is None or not bin_indices:
        return None
    axis = hist.GetXaxis()
    sorted_bins = sorted(int(bin_index) for bin_index in bin_indices)
    return [
        float(axis.GetBinLowEdge(sorted_bins[0])),
        float(axis.GetBinUpEdge(sorted_bins[-1])),
    ]


def compute_hist_closure_metrics(reference_hist, comparison_hist):
    metrics = {
        "signature_match": False,
        "reference_integral": _hist_integral(reference_hist),
        "comparison_integral": _hist_integral(comparison_hist),
        "integral_difference": None,
        "integral_ratio": None,
        "abs_diff_integral": None,
        "max_abs_bin_diff": None,
        "max_abs_bin_center": None,
    }
    if reference_hist is None or comparison_hist is None:
        return metrics
    if _hist_bin_signature(reference_hist) != _hist_bin_signature(comparison_hist):
        return metrics

    metrics["signature_match"] = True
    ref_integral = metrics["reference_integral"]
    cmp_integral = metrics["comparison_integral"]
    metrics["integral_difference"] = float(cmp_integral - ref_integral)
    if ref_integral > 0.0:
        metrics["integral_ratio"] = float(cmp_integral / ref_integral)

    abs_diff_integral = 0.0
    max_abs_bin_diff = 0.0
    max_abs_bin_center = None
    for bin_index in range(1, reference_hist.GetNbinsX() + 1):
        ref_value = float(reference_hist.GetBinContent(bin_index))
        cmp_value = float(comparison_hist.GetBinContent(bin_index))
        abs_diff = abs(cmp_value - ref_value)
        abs_diff_integral += abs_diff
        if abs_diff > max_abs_bin_diff:
            max_abs_bin_diff = abs_diff
            max_abs_bin_center = float(reference_hist.GetBinCenter(bin_index))

    metrics["abs_diff_integral"] = float(abs_diff_integral)
    metrics["max_abs_bin_diff"] = float(max_abs_bin_diff)
    metrics["max_abs_bin_center"] = max_abs_bin_center
    return metrics


def print_particle_subtraction_weight_support_warning(weight_payload, context="", **details):
    diagnostics = (weight_payload or {}).get("diagnostics") or {}
    unsupported_count = int(diagnostics.get("pion_weight_unsupported_bin_count", 0) or 0)
    overlap_count = int(diagnostics.get("pion_weight_unsupported_overlap_bin_count", 0) or 0)
    if unsupported_count <= 0 and overlap_count <= 0:
        return

    print("WARNING: unsupported pion-weight bins encountered")
    if context:
        print("  context = {}".format(context))
    for key, value in details.items():
        print("  {} = {}".format(key, value))
    print("  unsupported_bin_count = {}".format(unsupported_count))
    unsupported_mm_range = diagnostics.get("pion_weight_unsupported_mm_range")
    if unsupported_mm_range is not None:
        print(
            "  unsupported_mm_range = [{:.4f}, {:.4f}]".format(
                float(unsupported_mm_range[0]),
                float(unsupported_mm_range[1]),
            )
        )
    if overlap_count > 0:
        print("  overlap_pion_control_bin_count = {}".format(overlap_count))
        overlap_mm_range = diagnostics.get("pion_weight_unsupported_overlap_mm_range")
        if overlap_mm_range is not None:
            print(
                "  overlap_pion_control_mm_range = [{:.4f}, {:.4f}]".format(
                    float(overlap_mm_range[0]),
                    float(overlap_mm_range[1]),
                )
            )


def _resolve_component_template_map(component_fit_result):
    return {
        "pi_n": component_fit_result.get("H_simc_shape_pi_n"),
        "pi_delta": component_fit_result.get("H_simc_shape_pi_delta"),
        "pi_sidis": component_fit_result.get("H_simc_shape_pi_sidis"),
    }


def _resolve_active_component_names(component_fit_result):
    resolved = []
    diagnostics = {}
    if isinstance(component_fit_result, dict):
        diagnostics = component_fit_result.get("diagnostics") or {}
    for scope_name in ("pion", "kaon"):
        fit_order = ((diagnostics.get(scope_name) or {}).get("fit_order")) or []
        for component_name in fit_order:
            component_name = str(component_name)
            if component_name in COMPONENT_NAMES and component_name not in resolved:
                resolved.append(component_name)

    if not resolved and isinstance(component_fit_result, dict):
        resolved_config = component_fit_result.get("resolved_subtraction_config") or {}
        for fit_target in ("pion_control", "kaon_nosub"):
            anchor_windows = ((resolved_config.get(fit_target) or {}).get("anchor_windows")) or {}
            for component_name in anchor_windows.keys():
                component_name = str(component_name)
                if component_name in COMPONENT_NAMES and component_name not in resolved:
                    resolved.append(component_name)

    if not resolved and isinstance(component_fit_result, dict):
        template_map = _resolve_component_template_map(component_fit_result)
        for component_name in COMPONENT_NAMES:
            if _is_root_hist(template_map.get(component_name)):
                resolved.append(component_name)

    return resolved or list(COMPONENT_NAMES)


def handle_particle_subtraction_fallback(payload, reason, context="particle subtraction"):
    fallback_mode = None
    if isinstance(payload, dict):
        fallback_mode = payload.get("fallback_mode")
        payload["fallback_reason"] = reason
    fallback_mode = str(fallback_mode or "").strip().lower()
    if fallback_mode == "error":
        raise RuntimeError(
            "{} failed in simc_shape_components mode: {}".format(
                context,
                reason,
            )
        )
    return payload


def evaluate_particle_subtraction_component_fit_result(component_fit_result, inpDict=None):
    fallback_mode = resolve_particle_subtraction_fallback_mode(inpDict)
    failures = []
    diagnostics = {}

    if not isinstance(component_fit_result, dict):
        failures.append("component_fit_result is missing")
        return {
            "accepted": False,
            "fallback_mode": fallback_mode,
            "reason": "; ".join(failures),
            "diagnostics": diagnostics,
        }

    if component_fit_result.get("particle_subtraction_mode") != PARTICLE_SUBTRACTION_MODE_COMPONENTS:
        failures.append("particle_subtraction_mode is not simc_shape_components")

    fit_status_pion = str(component_fit_result.get("fit_status_pion", "")).strip().lower()
    fit_status_kaon = str(component_fit_result.get("fit_status_kaon", "")).strip().lower()
    if fit_status_pion != "success":
        failures.append("fit_status_pion={}".format(component_fit_result.get("fit_status_pion")))
    if fit_status_kaon != "success":
        failures.append("fit_status_kaon={}".format(component_fit_result.get("fit_status_kaon")))

    if bool(component_fit_result.get("fallback_used")):
        failures.append(
            "component fit fallback_used ({})".format(
                component_fit_result.get("fallback_reason") or "unknown reason"
            )
        )

    nested_diagnostics = component_fit_result.get("diagnostics") or {}
    pion_validation = (nested_diagnostics.get("pion") or {}).get("validation") or {}
    kaon_validation = (nested_diagnostics.get("kaon") or {}).get("validation") or {}
    diagnostics["fit_validation_pion"] = bool(pion_validation.get("accepted"))
    diagnostics["fit_validation_kaon"] = bool(kaon_validation.get("accepted"))
    if not diagnostics["fit_validation_pion"]:
        failures.append("pion validation rejected")
    if not diagnostics["fit_validation_kaon"]:
        failures.append("kaon validation rejected")

    amplitudes = {}
    for amplitude_name in ("A_n", "A_delta", "A_sidis", "B_n", "B_delta", "B_sidis"):
        try:
            amplitude_value = float(component_fit_result.get(amplitude_name))
        except Exception:
            amplitude_value = float("nan")
        amplitudes[amplitude_name] = amplitude_value
        if not math.isfinite(amplitude_value):
            failures.append("{} is non-finite".format(amplitude_name))
        elif amplitude_value < 0.0:
            failures.append("{} is negative".format(amplitude_name))
    diagnostics["amplitudes"] = amplitudes
    active_component_names = _resolve_active_component_names(component_fit_result)
    diagnostics["active_component_names"] = list(active_component_names)

    component_templates = _resolve_component_template_map(component_fit_result)
    template_signature = None
    for component_name in active_component_names:
        template_hist = component_templates.get(component_name)
        if not _is_root_hist(template_hist):
            failures.append("missing template {}".format(component_name))
            continue
        signature = _hist_bin_signature(template_hist)
        if template_signature is None:
            template_signature = signature
        elif signature != template_signature:
            failures.append("template binning mismatch for {}".format(component_name))
        if _hist_integral(template_hist) <= 0.0:
            failures.append("non-positive template integral for {}".format(component_name))

    diagnostics["template_signature"] = template_signature
    return {
        "accepted": not failures,
        "fallback_mode": fallback_mode,
        "reason": "; ".join(failures),
        "diagnostics": diagnostics,
    }


def build_simc_shape_pion_control_weights(
    component_fit_result,
    clip_min=0.0,
    clip_max=None,
    denom_floor=1e-12,
    model_variant="final",
):
    gate_result = evaluate_particle_subtraction_component_fit_result(component_fit_result)
    if not gate_result["accepted"]:
        raise ValueError(gate_result["reason"] or "component-fit result rejected")

    template_map = _resolve_component_template_map(component_fit_result)
    active_component_names = _resolve_active_component_names(component_fit_result)
    reference_hist = None
    for component_name in active_component_names:
        reference_hist = template_map.get(component_name)
        if _is_root_hist(reference_hist):
            break
    if reference_hist is None:
        raise ValueError("No valid reference histogram is available for pion weights")

    for component_name in active_component_names:
        template_hist = template_map[component_name]
        if _hist_bin_signature(template_hist) != _hist_bin_signature(reference_hist):
            raise ValueError("Template binning mismatch for {}".format(component_name))

    h_pion_control_model = _clone_hist(
        reference_hist,
        "H_pion_control_model_{}_{}".format(
            str(model_variant or "final"),
            component_fit_result.get("analysis_scope", "scope"),
        ),
        reset=True,
    )
    h_kaon_pion_model = _clone_hist(
        reference_hist,
        "H_kaon_pion_model_{}_{}".format(
            str(model_variant or "final"),
            component_fit_result.get("analysis_scope", "scope"),
        ),
        reset=True,
    )
    h_pion_weight = _clone_hist(
        reference_hist,
        "H_pion_weight_vs_MM_{}_{}".format(
            str(model_variant or "final"),
            component_fit_result.get("analysis_scope", "scope"),
        ),
        reset=True,
    )
    h_weighted_pion_control_model = _clone_hist(
        reference_hist,
        "H_weighted_pion_control_model_{}_{}".format(
            str(model_variant or "final"),
            component_fit_result.get("analysis_scope", "scope"),
        ),
        reset=True,
    )

    mode_key = str(model_variant or "final").strip().lower()
    diagnostics = component_fit_result.get("diagnostics") or {}
    pion_diagnostics = diagnostics.get("pion") or {}
    kaon_diagnostics = diagnostics.get("kaon") or {}
    if mode_key == "staged":
        pion_amplitude_source = pion_diagnostics.get("staged_amplitudes_scaled") or {}
        kaon_amplitude_source = kaon_diagnostics.get("staged_amplitudes_scaled") or {}
    else:
        pion_amplitude_source = (
            pion_diagnostics.get("refined_amplitudes")
            or pion_diagnostics.get("staged_amplitudes_scaled")
            or {}
        )
        kaon_amplitude_source = (
            kaon_diagnostics.get("refined_amplitudes")
            or kaon_diagnostics.get("staged_amplitudes_scaled")
            or {}
        )
    pion_amplitudes = {
        "pi_n": float(pion_amplitude_source.get("pi_n", component_fit_result.get("B_n", 0.0)) or 0.0),
        "pi_delta": float(pion_amplitude_source.get("pi_delta", component_fit_result.get("B_delta", 0.0)) or 0.0),
        "pi_sidis": float(pion_amplitude_source.get("pi_sidis", component_fit_result.get("B_sidis", 0.0)) or 0.0),
    }
    kaon_amplitudes = {
        "pi_n": float(kaon_amplitude_source.get("pi_n", component_fit_result.get("A_n", 0.0)) or 0.0),
        "pi_delta": float(kaon_amplitude_source.get("pi_delta", component_fit_result.get("A_delta", 0.0)) or 0.0),
        "pi_sidis": float(kaon_amplitude_source.get("pi_sidis", component_fit_result.get("A_sidis", 0.0)) or 0.0),
    }

    for component_name in active_component_names:
        h_pion_control_model.Add(template_map[component_name], pion_amplitudes[component_name])
        h_kaon_pion_model.Add(template_map[component_name], kaon_amplitudes[component_name])

    weights = np.zeros(reference_hist.GetNbinsX() + 2, dtype=np.float64)
    unclipped_values = []
    clipped_values = []
    clipped_bin_count = 0
    unsupported_bins = []
    ratio_max_abs_error = 0.0

    for bin_index in range(1, reference_hist.GetNbinsX() + 1):
        denominator = float(h_pion_control_model.GetBinContent(bin_index))
        numerator = float(h_kaon_pion_model.GetBinContent(bin_index))

        unsupported = denominator <= float(denom_floor) and numerator > float(denom_floor)
        if unsupported:
            unsupported_bins.append(bin_index)

        if denominator <= float(denom_floor):
            weight = 0.0
        else:
            weight = numerator / denominator

        if not math.isfinite(weight):
            weight = 0.0

        unclipped_weight = float(weight)
        if clip_min is not None:
            weight = max(float(clip_min), weight)
        if clip_max is not None:
            weight = min(float(clip_max), weight)
        if abs(weight - unclipped_weight) > 1e-15:
            clipped_bin_count += 1

        weights[bin_index] = float(weight)
        h_pion_weight.SetBinContent(bin_index, float(weight))
        h_pion_weight.SetBinError(bin_index, 0.0)
        weighted_model_value = float(weight) * denominator
        h_weighted_pion_control_model.SetBinContent(bin_index, weighted_model_value)
        h_weighted_pion_control_model.SetBinError(bin_index, 0.0)
        unclipped_values.append(unclipped_weight)
        clipped_values.append(float(weight))
        ratio_max_abs_error = max(
            ratio_max_abs_error,
            abs((float(weight) * denominator) - numerator),
        )

    axis = reference_hist.GetXaxis()
    unsupported_mm_range = None
    if unsupported_bins:
        unsupported_mm_range = [
            float(axis.GetBinLowEdge(min(unsupported_bins))),
            float(axis.GetBinUpEdge(max(unsupported_bins))),
        ]

    pion_reference = component_fit_result.get("H_pion_control_input")
    unsupported_overlap_bins = []
    if pion_reference is not None and unsupported_bins:
        for bin_index in unsupported_bins:
            if float(pion_reference.GetBinContent(int(bin_index))) > 0.0:
                unsupported_overlap_bins.append(int(bin_index))

    diagnostics = {
        "model_variant": mode_key,
        "active_component_names": list(active_component_names),
        "pion_control_integral": _hist_integral(component_fit_result.get("H_pion_control_input")),
        "pion_control_model_integral": _hist_integral(h_pion_control_model),
        "kaon_pion_model_integral": _hist_integral(h_kaon_pion_model),
        "weighted_pion_control_model_integral": _hist_integral(h_weighted_pion_control_model),
        "pion_weight_min_unclipped": min(unclipped_values) if unclipped_values else 0.0,
        "pion_weight_max_unclipped": max(unclipped_values) if unclipped_values else 0.0,
        "pion_weight_min": min(clipped_values) if clipped_values else 0.0,
        "pion_weight_max": max(clipped_values) if clipped_values else 0.0,
        "pion_weight_mean": float(np.mean(clipped_values)) if clipped_values else 0.0,
        "pion_weight_rms": float(np.std(clipped_values)) if clipped_values else 0.0,
        "pion_weight_clipped_bin_count": int(clipped_bin_count),
        "pion_weight_unsupported_bin_count": int(len(unsupported_bins)),
        "pion_weight_unsupported_mm_range": unsupported_mm_range,
        "unsupported_bins": list(unsupported_bins),
        "pion_weight_unsupported_overlap_bin_count": int(len(unsupported_overlap_bins)),
        "pion_weight_unsupported_overlap_mm_range": _bin_index_mm_range(
            pion_reference,
            unsupported_overlap_bins,
        ),
        "unsupported_overlap_bins": list(unsupported_overlap_bins),
        "denom_floor": float(denom_floor),
        "ratio_max_abs_error": float(ratio_max_abs_error),
        "ratio_consistency_ok": bool(ratio_max_abs_error <= 1e-9),
        "clip_min": None if clip_min is None else float(clip_min),
        "clip_max": None if clip_max is None else float(clip_max),
    }
    diagnostics["model_closure"] = compute_hist_closure_metrics(
        h_kaon_pion_model,
        h_weighted_pion_control_model,
    )

    return {
        "weights": weights,
        "H_pion_control_model": h_pion_control_model,
        "H_kaon_pion_model": h_kaon_pion_model,
        "H_weighted_pion_control_model": h_weighted_pion_control_model,
        "H_pion_weight_vs_MM": h_pion_weight,
        "diagnostics": diagnostics,
    }


def simc_shape_pion_weight_from_value(mm_value, h_reference, pion_mm_weights):
    if h_reference is None or pion_mm_weights is None:
        return 0.0
    if not math.isfinite(float(mm_value)):
        return 0.0
    try:
        bin_index = int(h_reference.GetXaxis().FindBin(float(mm_value)))
    except Exception:
        return 0.0
    if bin_index <= 0 or bin_index >= len(pion_mm_weights):
        return 0.0
    try:
        weight = float(pion_mm_weights[bin_index])
    except Exception:
        return 0.0
    return weight if math.isfinite(weight) else 0.0


def iter_component_control_source_specs(sub_event_cache, normfac_data, normfac_dummy, nWindows, positive_template=True):
    if sub_event_cache is None:
        return []
    sign = 1.0 if positive_template else -1.0
    return [
        {
            "label": "prompt",
            "cache_section": sub_event_cache.get("prompt"),
            "coefficient": sign * float(normfac_data),
        },
        {
            "label": "rand",
            "cache_section": sub_event_cache.get("rand"),
            "coefficient": sign * (-float(normfac_data) / float(nWindows)),
        },
        {
            "label": "dummy",
            "cache_section": sub_event_cache.get("dummy"),
            "coefficient": sign * (-float(normfac_dummy)),
        },
        {
            "label": "dummy_rand",
            "cache_section": sub_event_cache.get("dummy_rand"),
            "coefficient": sign * (float(normfac_dummy) / float(nWindows)),
        },
    ]


def _resolve_selected_indices(cache_section, scope_selection, mask_key):
    if cache_section is None:
        return np.asarray([], dtype=np.int32)

    t_index = scope_selection.get("t_index")
    phi_index = scope_selection.get("phi_index")
    if phi_index is not None and t_index is not None:
        index_map_key = "{}_bin_index".format(mask_key)
        if index_map_key in cache_section:
            return cache_section[index_map_key].get((int(t_index), int(phi_index)), np.asarray([], dtype=np.int32))

    mask = np.asarray(cache_section.get(mask_key, []), dtype=bool)
    if mask.size == 0:
        return np.asarray([], dtype=np.int32)
    if t_index is not None:
        mask = np.logical_and(mask, np.asarray(cache_section.get("t_index", []), dtype=np.int32) == int(t_index))
    if phi_index is not None:
        mask = np.logical_and(mask, np.asarray(cache_section.get("phi_index", []), dtype=np.int32) == int(phi_index))
    return np.flatnonzero(mask).astype(np.int32)


def _fill_template_hist_keys(template_hists, cache_section, idx, event_weight, particle_type, pol):
    if event_weight == 0.0:
        return
    adj_t = float(cache_section["adj_t"][idx])
    adj_mm = float(cache_section["adj_MM"][idx])
    q2 = float(cache_section["Q2"][idx])
    w = float(cache_section["W"][idx])
    epsilon = float(cache_section["epsilon"][idx])
    theta_cm_deg = float(cache_section.get("theta_cm_deg", [float("nan")])[idx])
    ssxptar = float(cache_section["ssxptar"][idx])
    ssyptar = float(cache_section["ssyptar"][idx])
    hsxptar = float(cache_section["hsxptar"][idx])
    hsyptar = float(cache_section["hsyptar"][idx])

    if "t" in template_hists:
        template_hists["t"].Fill(adj_t, event_weight)
    if "Q2" in template_hists:
        template_hists["Q2"].Fill(q2, event_weight)
    if "W" in template_hists:
        template_hists["W"].Fill(w, event_weight)
    if "epsilon" in template_hists:
        template_hists["epsilon"].Fill(epsilon, event_weight)
    if "q2_w" in template_hists:
        template_hists["q2_w"].Fill(q2, w, event_weight)
    if "theta_cm" in template_hists and math.isfinite(theta_cm_deg):
        template_hists["theta_cm"].Fill(theta_cm_deg, event_weight)
    if "ssxptar" in template_hists:
        template_hists["ssxptar"].Fill(ssxptar, event_weight)
    if "ssyptar" in template_hists:
        template_hists["ssyptar"].Fill(ssyptar, event_weight)
    if "hsxptar" in template_hists:
        template_hists["hsxptar"].Fill(hsxptar, event_weight)
    if "hsyptar" in template_hists:
        template_hists["hsyptar"].Fill(hsyptar, event_weight)
    if "mm" in template_hists:
        template_hists["mm"].Fill(adj_mm, event_weight)
    if "t_vs_tmin" in template_hists:
        from theta_cm import calculate_tmin

        minus_tmin = calculate_tmin(particle_type, pol, w, q2)
        if math.isfinite(minus_tmin) and math.isfinite(adj_t):
            template_hists["t_vs_tmin"].Fill(minus_tmin, adj_t, event_weight)


def fill_simc_shape_pion_subtraction_templates(
    template_hists,
    source_specs,
    h_pion_reference,
    pion_mm_weights,
    scope_selection,
    particle_type,
    pol,
    mm_background_reference_hist=None,
    mm_background_weights=None,
    residual_weights=None,
):
    stats = {
        "n_events_allcuts": 0,
        "n_events_nommcuts": 0,
        "sum_event_weight": 0.0,
        "sum_event_weight_sq": 0.0,
    }

    for source_spec in source_specs or []:
        cache_section = source_spec.get("cache_section")
        coefficient = float(source_spec.get("coefficient", 0.0) or 0.0)
        if cache_section is None or coefficient == 0.0:
            continue

        allcut_indices = _resolve_selected_indices(cache_section, scope_selection, "allcuts")
        nommcut_indices = _resolve_selected_indices(cache_section, scope_selection, "nommcuts")

        for idx in nommcut_indices:
            adj_mm = float(cache_section["adj_MM"][idx])
            pion_weight = simc_shape_pion_weight_from_value(adj_mm, h_pion_reference, pion_mm_weights)
            event_weight = coefficient * pion_weight
            if event_weight == 0.0:
                continue
            if "mm_nosub" in template_hists:
                template_hists["mm_nosub"].Fill(adj_mm, event_weight)
            if "fit1sub" in template_hists:
                template_hists["fit1sub"].Fill(adj_mm, event_weight)
            if "pisub" in template_hists:
                template_hists["pisub"].Fill(adj_mm, event_weight)
            stats["n_events_nommcuts"] += 1

        for idx in allcut_indices:
            adj_mm = float(cache_section["adj_MM"][idx])
            pion_weight = simc_shape_pion_weight_from_value(adj_mm, h_pion_reference, pion_mm_weights)
            event_weight = coefficient * pion_weight
            if mm_background_reference_hist is not None and mm_background_weights is not None:
                event_weight *= mm_background_weight_from_value(
                    adj_mm,
                    mm_background_reference_hist,
                    mm_background_weights,
                    residual_weights=residual_weights,
                )
            if event_weight == 0.0 or not math.isfinite(event_weight):
                continue
            _fill_template_hist_keys(template_hists, cache_section, idx, event_weight, particle_type, pol)
            stats["n_events_allcuts"] += 1
            stats["sum_event_weight"] += float(event_weight)
            stats["sum_event_weight_sq"] += float(event_weight * event_weight)

    return stats


def summarize_particle_subtraction_component_payload(payload):
    if not isinstance(payload, dict):
        return {}
    summary = {}
    for key, value in payload.items():
        if key.startswith("H_") and _is_root_hist(value):
            continue
        cleaned_value = _json_ready_particle_subtraction_value(value)
        if cleaned_value is _JSON_SKIP:
            continue
        summary[key] = cleaned_value
    return summary
