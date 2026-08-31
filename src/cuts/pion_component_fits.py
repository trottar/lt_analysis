#! /usr/bin/python

from __future__ import annotations

import itertools
import csv
import hashlib
import json
import math
import os
import tempfile
from copy import deepcopy
from datetime import datetime, timezone

import numpy as np
import ROOT
from scipy.optimize import least_squares, lsq_linear
from scipy.stats import chi2 as chi2_dist

from background_config import (
    BG_OPT_MM_PLOT_MAX,
    BG_OPT_MM_PLOT_MIN,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    PION_COMPONENT_DYNAMIC_ALIGNMENT_SCHEMA_VERSION,
    get_particle_subtraction_setting_key,
    get_pion_component_dynamic_alignment_config,
    resolve_particle_subtraction_component_fit_mode,
    resolve_particle_subtraction_component_postfit_scales,
    resolve_particle_subtraction_component_postrefine_scales,
    resolve_particle_subtraction_component_prior_scales,
    resolve_particle_subtraction_component_residual_shift_settings,
    resolve_particle_subtraction_component_stage_amplitude_modes,
    resolve_particle_subtraction_component_stage_amplitude_windows,
    resolve_particle_subtraction_component_fit_excluded_windows,
    resolve_particle_subtraction_component_cleanup_validation_mm_max,
    get_particle_subtraction_component_fit_window_config,
    get_proton_contamination_cleaning_config,
    resolve_particle_subtraction_component_fit_windows,
    resolve_particle_subtraction_mode,
    resolve_particle_subtraction_root_ownership_debug,
)
from utility import normalize_hist_to_unit_area
from root_histogram_ownership import (
    clone_root_histogram,
    configure_particle_subtraction_root_ownership_debug,
    unique_root_object_name,
)


COMPONENT_NAMES = ("pi_n", "pi_delta", "pi_sidis")
KAON_SIGNAL_TEMPLATE_NAME = "k_lambda_signal"
KAON_SIGMA0_TEMPLATE_NAME = "k_sigma0_signal"
COMPONENT_PLOT_STYLE = {
    "pi_n": {"label": "pi-n", "color": ROOT.kRed + 1},
    "pi_delta": {"label": "pi-delta", "color": ROOT.kAzure + 2},
    "pi_sidis": {"label": "pi-SIDIS", "color": ROOT.kMagenta + 2},
    KAON_SIGNAL_TEMPLATE_NAME: {"label": "K-Lambda", "color": ROOT.kBlue + 1},
    KAON_SIGMA0_TEMPLATE_NAME: {"label": "K-Sigma0", "color": ROOT.kCyan + 2},
}


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


def _is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def _clone_hist(template_hist, name, title=None, reset=False):
    return clone_root_histogram(
        template_hist,
        scope="component_fit",
        role="fit_histogram",
        name=name,
        optional=True,
        reset=reset,
        sumw2=False,
        title=title,
    )


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def record_particle_subtraction_page(page_manifest, page_id, *, scope, authoritative):
    """Append one serializable, uniquely identified particle-subtraction page.

    The manifest is deliberately caller-owned: renderers only append scalar
    records after they have emitted their page, so it cannot prolong a ROOT
    object's lifetime or conceal duplicate diagnostic output.
    """
    if page_manifest is None:
        return
    if not isinstance(page_manifest, list):
        raise TypeError("particle-subtraction page manifest must be a list")
    normalized_id = str(page_id or "").strip()
    if not normalized_id:
        raise ValueError("particle-subtraction page manifest requires a page_id")
    if any(str((entry or {}).get("page_id") or "") == normalized_id for entry in page_manifest):
        raise RuntimeError("duplicate_particle_subtraction_page_id:{}".format(normalized_id))
    page_manifest.append(
        {
            "page_id": normalized_id,
            "scope": str(scope or "unknown"),
            "authoritative": bool(authoritative),
        }
    )


def _record_component_page(page_manifest, page_id_prefix, page_key, *, scope, authoritative):
    if page_id_prefix:
        record_particle_subtraction_page(
            page_manifest,
            "{}.{}".format(str(page_id_prefix).rstrip("."), page_key),
            scope=scope,
            authoritative=authoritative,
        )


def _sample_hist_value_and_variance(hist, x_value, interpolation_mode="linear"):
    if hist is None or not _is_finite_number(x_value):
        return 0.0, 0.0
    mode = str(interpolation_mode or "linear").strip().lower() or "linear"
    if mode not in ("linear", "nearest"):
        raise ValueError("Unsupported template interpolation mode '{}'".format(interpolation_mode))

    nbins = int(hist.GetNbinsX())
    if nbins <= 0:
        return 0.0, 0.0

    first_center = float(hist.GetBinCenter(1))
    last_center = float(hist.GetBinCenter(nbins))
    if x_value < first_center or x_value > last_center:
        return 0.0, 0.0

    if mode == "nearest":
        nearest_bin = 1
        nearest_distance = abs(float(x_value) - first_center)
        for bin_index in range(2, nbins + 1):
            candidate_distance = abs(float(x_value) - float(hist.GetBinCenter(bin_index)))
            if candidate_distance < nearest_distance:
                nearest_bin = int(bin_index)
                nearest_distance = candidate_distance
        error_value = float(hist.GetBinError(nearest_bin))
        return float(hist.GetBinContent(nearest_bin)), float(error_value ** 2)

    axis = hist.GetXaxis()
    anchor_bin = int(axis.FindBin(float(x_value)))
    anchor_bin = max(1, min(anchor_bin, nbins))
    anchor_center = float(hist.GetBinCenter(anchor_bin))
    if abs(float(x_value) - anchor_center) <= 1e-12:
        error_value = float(hist.GetBinError(anchor_bin))
        return float(hist.GetBinContent(anchor_bin)), float(error_value ** 2)

    if float(x_value) > anchor_center:
        left_bin = int(anchor_bin)
        right_bin = int(anchor_bin + 1)
    else:
        left_bin = int(anchor_bin - 1)
        right_bin = int(anchor_bin)
    if left_bin < 1 or right_bin > nbins:
        return 0.0, 0.0

    left_center = float(hist.GetBinCenter(left_bin))
    right_center = float(hist.GetBinCenter(right_bin))
    if abs(right_center - left_center) <= 1e-12:
        error_value = float(hist.GetBinError(anchor_bin))
        return float(hist.GetBinContent(anchor_bin)), float(error_value ** 2)

    blend = (float(x_value) - left_center) / (right_center - left_center)
    blend = min(max(blend, 0.0), 1.0)
    left_value = float(hist.GetBinContent(left_bin))
    right_value = float(hist.GetBinContent(right_bin))
    left_variance = float(hist.GetBinError(left_bin) ** 2)
    right_variance = float(hist.GetBinError(right_bin) ** 2)
    value = ((1.0 - blend) * left_value) + (blend * right_value)
    variance = (((1.0 - blend) ** 2) * left_variance) + ((blend ** 2) * right_variance)
    return float(value), float(variance)


def build_shifted_template_histogram(
    source_hist,
    delta_mm,
    sign_convention,
    output_name,
    interpolation_mode="linear",
    renormalize=True,
):
    if source_hist is None:
        return None, {}

    resolved_sign = str(sign_convention or "").strip().lower()
    if resolved_sign != "positive_moves_peak_higher_mm":
        raise ValueError(
            "Unsupported residual template shift sign convention '{}'".format(sign_convention)
        )

    shift_value = float(delta_mm) if _is_finite_number(delta_mm) else 0.0
    original_integral = _hist_integral(source_hist)
    shifted_hist = _clone_hist(source_hist, output_name, reset=True)
    if shifted_hist is None:
        return None, {}

    if abs(shift_value) <= 1e-12:
        shifted_hist.Add(source_hist)
        shifted_integral_before = _hist_integral(shifted_hist)
        renorm_factor = 1.0
        shifted_integral_after = shifted_integral_before
        if renormalize and original_integral > 0.0 and abs(shifted_integral_before - original_integral) > 1e-12:
            renorm_factor = float(original_integral / shifted_integral_before)
            shifted_hist.Scale(renorm_factor)
            shifted_integral_after = _hist_integral(shifted_hist)
        return shifted_hist, {
            "delta_mm": shift_value,
            "sign_convention": resolved_sign,
            "interpolation_mode": str(interpolation_mode or "linear"),
            "original_integral": float(original_integral),
            "shifted_integral_before_renorm": float(shifted_integral_before),
            "shifted_integral_after_renorm": float(shifted_integral_after),
            "shift_renormalization_factor": float(renorm_factor),
            "lost_integral_fraction": 0.0,
            "shift_bound_hit_flag": False,
        }

    for bin_index in range(1, shifted_hist.GetNbinsX() + 1):
        output_center = float(shifted_hist.GetBinCenter(bin_index))
        sample_x = output_center - shift_value
        sample_value, sample_variance = _sample_hist_value_and_variance(
            source_hist,
            sample_x,
            interpolation_mode=interpolation_mode,
        )
        shifted_hist.SetBinContent(bin_index, float(sample_value))
        shifted_hist.SetBinError(bin_index, math.sqrt(max(float(sample_variance), 0.0)))

    shifted_integral_before = _hist_integral(shifted_hist)
    renorm_factor = 1.0
    if renormalize and original_integral > 0.0 and shifted_integral_before > 0.0:
        renorm_factor = float(original_integral / shifted_integral_before)
        shifted_hist.Scale(renorm_factor)
    shifted_integral_after = _hist_integral(shifted_hist)
    lost_integral_fraction = 0.0
    if original_integral > 0.0 and shifted_integral_before < original_integral:
        lost_integral_fraction = float(
            max(original_integral - shifted_integral_before, 0.0) / original_integral
        )
    return shifted_hist, {
        "delta_mm": shift_value,
        "sign_convention": resolved_sign,
        "interpolation_mode": str(interpolation_mode or "linear"),
        "original_integral": float(original_integral),
        "shifted_integral_before_renorm": float(shifted_integral_before),
        "shifted_integral_after_renorm": float(shifted_integral_after),
        "shift_renormalization_factor": float(renorm_factor),
        "lost_integral_fraction": float(lost_integral_fraction),
        "shift_bound_hit_flag": False,
    }


def _build_mm_shifted_hist(template_hist, shift, name, renormalize=False):
    if template_hist is None:
        return None
    shifted_hist, _ = build_shifted_template_histogram(
        template_hist,
        shift,
        "positive_moves_peak_higher_mm",
        name,
        interpolation_mode="nearest",
        renormalize=bool(renormalize),
    )
    return shifted_hist


def _set_hist_values(hist, values):
    if hist is None:
        return
    for index, value in enumerate(values, start=1):
        hist.SetBinContent(index, float(value))
        hist.SetBinError(index, 0.0)


def _merge_template_hist_maps(*template_maps):
    merged = {}
    for template_map in template_maps:
        if not isinstance(template_map, dict):
            continue
        for template_name, template_hist in template_map.items():
            if template_hist is None:
                continue
            merged[str(template_name)] = template_hist
    return merged


def _zero_fit_result(
    target_hist,
    amplitude_prefix,
    context,
    fallback_reason,
    template_hists=None,
):
    fit_hist = _clone_hist(
        target_hist,
        "{}_fit_hist_{}".format(amplitude_prefix, context),
        reset=True,
    )
    residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_{}".format(amplitude_prefix, context),
        reset=True,
    )
    if residual_hist is not None and target_hist is not None:
        residual_hist.Add(target_hist)

    result = {
        "{}_n".format(amplitude_prefix): 0.0,
        "{}_delta".format(amplitude_prefix): 0.0,
        "{}_sidis".format(amplitude_prefix): 0.0,
        "fit_status": "fallback",
        "diagnostics": {
            "success": False,
            "status_code": None,
            "message": fallback_reason,
            "chi2": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "n_fit_bins": 0,
            "fallback_used": True,
            "fallback_reason": fallback_reason,
        },
        "fit_hist": fit_hist,
        "residual_hist": residual_hist,
        "template_hists": _merge_template_hist_maps(template_hists),
        "pi_n_scaled_hist": _clone_hist(
            target_hist,
            "{}_pi_n_scaled_{}".format(amplitude_prefix, context),
            reset=True,
        ),
        "pi_delta_scaled_hist": _clone_hist(
            target_hist,
            "{}_pi_delta_scaled_{}".format(amplitude_prefix, context),
            reset=True,
        ),
        "pi_sidis_scaled_hist": _clone_hist(
            target_hist,
            "{}_pi_sidis_scaled_{}".format(amplitude_prefix, context),
            reset=True,
        ),
        "extra_scaled_hists": {},
        "extra_component_amplitudes": {},
    }
    if amplitude_prefix == "A":
        result["pion_bg_fit_hist"] = _clone_hist(
            target_hist,
            "{}_pion_bg_fit_{}".format(amplitude_prefix, context),
            reset=True,
        )
    return result


def _is_fine_binned_analysis_scope(analysis_scope):
    scope_text = str(analysis_scope or "").strip().lower()
    if not scope_text:
        return False
    if scope_text in ("setting-wide", "setting_wide", "settingwide"):
        return False
    return ("t_bin" in scope_text) or ("phi_bin" in scope_text)


def _promote_sparse_zero_solution_if_applicable(fit_result, amplitude_prefix, analysis_scope):
    if not isinstance(fit_result, dict):
        return fit_result
    if not _is_fine_binned_analysis_scope(analysis_scope):
        return fit_result
    if str(fit_result.get("fit_status") or "").strip().lower() != "fallback":
        return fit_result

    diagnostics = deepcopy(fit_result.get("diagnostics") or {})
    message_parts = [
        diagnostics.get("fallback_reason"),
        diagnostics.get("message"),
        ((diagnostics.get("joint_refinement") or {}).get("message")),
    ]
    message_blob = " ; ".join(
        str(part).strip().lower()
        for part in message_parts
        if str(part or "").strip()
    )
    sparse_markers = (
        "no active templates available for joint refinement",
        "no valid full-range fit bins",
        "no valid template columns available for joint refinement",
        "non-positive integral for",
        "no valid fit bins",
        "insufficient fit bins",
    )
    if not any(marker in message_blob for marker in sparse_markers):
        return fit_result

    component_amplitudes = [
        float(fit_result.get("{}_n".format(amplitude_prefix), 0.0) or 0.0),
        float(fit_result.get("{}_delta".format(amplitude_prefix), 0.0) or 0.0),
        float(fit_result.get("{}_sidis".format(amplitude_prefix), 0.0) or 0.0),
    ]
    extra_component_amplitudes = [
        float(value or 0.0)
        for value in ((fit_result.get("extra_component_amplitudes") or {}).values())
    ]
    if any(abs(value) > 1e-12 for value in component_amplitudes + extra_component_amplitudes):
        return fit_result

    promoted_message = (
        "sparse fine-bin zero solution accepted: {}".format(
            diagnostics.get("fallback_reason")
            or diagnostics.get("message")
            or "no valid active templates"
        )
    )
    validation = deepcopy(diagnostics.get("validation") or {})
    validation["accepted"] = True
    validation["rejection_reasons"] = []
    if "chi2" not in validation:
        validation["chi2"] = diagnostics.get("chi2")
    if "ndf" not in validation:
        validation["ndf"] = diagnostics.get("ndf")
    if "chi2_ndf" not in validation:
        validation["chi2_ndf"] = diagnostics.get("chi2_ndf")
    if "fit_p_value" not in validation:
        validation["fit_p_value"] = diagnostics.get("fit_p_value")

    joint_refinement = deepcopy(diagnostics.get("joint_refinement") or {})
    joint_refinement["status"] = "skipped_sparse_bin"
    joint_refinement["success"] = False
    joint_refinement["message"] = promoted_message

    diagnostics["validation"] = validation
    diagnostics["success"] = True
    diagnostics["fallback_used"] = False
    diagnostics["fallback_reason"] = ""
    diagnostics["message"] = promoted_message
    diagnostics["joint_refinement"] = joint_refinement
    diagnostics["joint_refinement_status"] = "skipped_sparse_bin"
    diagnostics["sparse_bin_zero_solution_accepted"] = True
    fit_result["fit_status"] = "success"
    fit_result["diagnostics"] = diagnostics
    return fit_result


def _validate_template_hist(template_hist, target_hist, template_name):
    if target_hist is None:
        return "missing target histogram"
    if template_hist is None:
        return "missing SIMC template shape for {}".format(template_name)
    nbins = target_hist.GetNbinsX()
    xmin = target_hist.GetXaxis().GetXmin()
    xmax = target_hist.GetXaxis().GetXmax()
    if template_hist.GetNbinsX() != nbins:
        return "bin-count mismatch for {}".format(template_name)
    if (
        abs(template_hist.GetXaxis().GetXmin() - xmin) > 1e-9
        or abs(template_hist.GetXaxis().GetXmax() - xmax) > 1e-9
    ):
        return "axis-range mismatch for {}".format(template_name)
    if template_hist.Integral() <= 0.0:
        return "non-positive integral for {}".format(template_name)
    return ""


def _validate_component_shapes(component_hists, target_hist, component_names=None):
    if target_hist is None:
        return "missing target histogram"
    requested_names = [
        str(component_name)
        for component_name in (component_names or COMPONENT_NAMES)
        if str(component_name) in COMPONENT_NAMES
    ]
    if not requested_names:
        requested_names = list(COMPONENT_NAMES)
    if any(component_hists.get(name) is None for name in requested_names):
        return "missing SIMC component shape"

    for component_name in requested_names:
        message = _validate_template_hist(
            component_hists[component_name],
            target_hist,
            component_name,
        )
        if message:
            return message
    return ""


def _build_fit_inputs(
    target_hist,
    component_hists,
    fit_min,
    fit_max,
    extra_template_hists=None,
    include_windows=None,
    exclude_windows=None,
):
    if extra_template_hists is None:
        extra_template_hists = {}
    if include_windows is None:
        include_windows = []
    if exclude_windows is None:
        exclude_windows = []

    component_columns = {name: [] for name in COMPONENT_NAMES}
    extra_template_columns = {name: [] for name in extra_template_hists}
    x_values = []
    y_values = []
    sigma_values = []
    fit_bin_indices = []

    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < fit_min or x_center > fit_max:
            continue
        if include_windows and (
            not any(window_min <= x_center <= window_max for window_min, window_max in include_windows)
        ):
            continue
        if any(window_min <= x_center <= window_max for window_min, window_max in exclude_windows):
            continue

        y_value = float(target_hist.GetBinContent(bin_index))
        sigma_value = float(target_hist.GetBinError(bin_index))
        if not math.isfinite(y_value):
            continue
        if (not math.isfinite(sigma_value)) or sigma_value <= 0.0:
            sigma_value = max(math.sqrt(abs(y_value)), 1.0)

        x_values.append(x_center)
        y_values.append(y_value)
        sigma_values.append(sigma_value)
        fit_bin_indices.append(bin_index)
        for component_name in COMPONENT_NAMES:
            component_columns[component_name].append(
                float(component_hists[component_name].GetBinContent(bin_index))
            )
        for template_name, template_hist in extra_template_hists.items():
            extra_template_columns[template_name].append(
                float(template_hist.GetBinContent(bin_index))
            )

    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "sigma": np.asarray(sigma_values, dtype=float),
        "fit_bin_indices": fit_bin_indices,
        "component_columns": {
            component_name: np.asarray(values, dtype=float)
            for component_name, values in component_columns.items()
        },
        "extra_template_columns": {
            template_name: np.asarray(values, dtype=float)
            for template_name, values in extra_template_columns.items()
        },
    }


def _build_single_template_fit_inputs(
    target_hist,
    template_hist,
    fit_min,
    fit_max,
    include_windows=None,
    exclude_windows=None,
):
    if include_windows is None:
        include_windows = []
    if exclude_windows is None:
        exclude_windows = []

    x_values = []
    y_values = []
    sigma_values = []
    template_values = []
    fit_bin_indices = []

    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < fit_min or x_center > fit_max:
            continue
        if include_windows and (
            not any(window_min <= x_center <= window_max for window_min, window_max in include_windows)
        ):
            continue
        if any(window_min <= x_center <= window_max for window_min, window_max in exclude_windows):
            continue

        y_value = float(target_hist.GetBinContent(bin_index))
        sigma_value = float(target_hist.GetBinError(bin_index))
        if not math.isfinite(y_value):
            continue
        if (not math.isfinite(sigma_value)) or sigma_value <= 0.0:
            sigma_value = max(math.sqrt(abs(y_value)), 1.0)

        x_values.append(x_center)
        y_values.append(y_value)
        sigma_values.append(sigma_value)
        template_values.append(float(template_hist.GetBinContent(bin_index)))
        fit_bin_indices.append(bin_index)

    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "sigma": np.asarray(sigma_values, dtype=float),
        "template": np.asarray(template_values, dtype=float),
        "fit_bin_indices": fit_bin_indices,
    }


def _solve_nonnegative_template_amplitude(
    target_hist,
    template_hist,
    fit_min,
    fit_max,
    include_windows=None,
    exclude_windows=None,
    amplitude_windows=None,
    amplitude_mode="least_squares",
):
    template_name = getattr(template_hist, "GetName", lambda: "template")()
    validation_message = _validate_template_hist(template_hist, target_hist, template_name)
    if validation_message:
        return {
            "success": False,
            "amplitude": 0.0,
            "sigma": None,
            "chi2": None,
            "n_fit_bins": 0,
            "message": validation_message,
            "amplitude_mode": str(amplitude_mode or "least_squares"),
            "amplitude_windows": deepcopy(amplitude_windows or include_windows or []),
        }

    fit_inputs = _build_single_template_fit_inputs(
        target_hist,
        template_hist,
        fit_min,
        fit_max,
        include_windows=include_windows,
        exclude_windows=exclude_windows,
    )
    if len(fit_inputs["x"]) == 0:
        return {
            "success": False,
            "amplitude": 0.0,
            "sigma": None,
            "chi2": None,
            "n_fit_bins": 0,
            "message": "no valid fit bins",
            "amplitude_mode": str(amplitude_mode or "least_squares"),
            "amplitude_windows": deepcopy(amplitude_windows or include_windows or []),
        }

    resolved_amplitude_mode = str(amplitude_mode or "least_squares")
    resolved_amplitude_windows = deepcopy(amplitude_windows or include_windows or [])
    if resolved_amplitude_mode == "window_integral":
        amplitude_inputs = _build_single_template_fit_inputs(
            target_hist,
            template_hist,
            fit_min,
            fit_max,
            include_windows=amplitude_windows,
            exclude_windows=exclude_windows,
        )
        if len(amplitude_inputs["x"]) == 0:
            return {
                "success": False,
                "amplitude": 0.0,
                "sigma": None,
                "chi2": None,
                "n_fit_bins": 0,
                "message": "no valid amplitude-window bins",
                "amplitude_mode": resolved_amplitude_mode,
                "amplitude_windows": resolved_amplitude_windows,
                "fit_diagnostics": {},
            }
        template_sum = float(np.sum(amplitude_inputs["template"]))
        if (not math.isfinite(template_sum)) or template_sum <= 0.0:
            return {
                "success": False,
                "amplitude": 0.0,
                "sigma": None,
                "chi2": None,
                "n_fit_bins": int(len(amplitude_inputs["x"])),
                "message": "template has zero support inside amplitude window",
                "amplitude_mode": resolved_amplitude_mode,
                "amplitude_windows": resolved_amplitude_windows,
                "fit_diagnostics": {
                    "target_sum": float(np.sum(amplitude_inputs["y"])),
                    "template_sum": template_sum,
                },
            }
        target_sum = float(np.sum(amplitude_inputs["y"]))
        amplitude = max(target_sum / template_sum, 0.0)
        sigma_value = math.sqrt(float(np.sum(np.square(amplitude_inputs["sigma"])))) / template_sum
        fit_diagnostics = {
            "target_sum": target_sum,
            "template_sum": template_sum,
            "estimator": "window_integral",
            "fit_bin_indices": [int(value) for value in amplitude_inputs["fit_bin_indices"]],
            "fit_y": [float(value) for value in amplitude_inputs["y"]],
            "fit_sigma": [float(value) for value in amplitude_inputs["sigma"]],
            "fit_template": [float(value) for value in amplitude_inputs["template"]],
        }
    else:
        weighted_template = fit_inputs["template"] / fit_inputs["sigma"]
        denominator = float(np.dot(weighted_template, weighted_template))
        if (not math.isfinite(denominator)) or denominator <= 0.0:
            return {
                "success": False,
                "amplitude": 0.0,
                "sigma": None,
                "chi2": None,
                "n_fit_bins": int(len(fit_inputs["x"])),
                "message": "template has zero support inside anchor window",
                "amplitude_mode": resolved_amplitude_mode,
                "amplitude_windows": resolved_amplitude_windows,
                "fit_diagnostics": {},
            }

        weighted_target = fit_inputs["y"] / fit_inputs["sigma"]
        numerator = float(np.dot(weighted_template, weighted_target))
        amplitude = max(float(numerator / denominator), 0.0)
        sigma_value = 1.0 / math.sqrt(denominator)
        fit_diagnostics = {
            "weighted_numerator": numerator,
            "weighted_denominator": denominator,
            "estimator": "least_squares",
            "fit_bin_indices": [int(value) for value in fit_inputs["fit_bin_indices"]],
            "fit_y": [float(value) for value in fit_inputs["y"]],
            "fit_sigma": [float(value) for value in fit_inputs["sigma"]],
            "fit_template": [float(value) for value in fit_inputs["template"]],
        }

    residual = fit_inputs["y"] - amplitude * fit_inputs["template"]
    chi2_value = float(np.sum(np.square(residual / fit_inputs["sigma"])))
    return {
        "success": True,
        "amplitude": amplitude,
        "sigma": sigma_value,
        "chi2": chi2_value,
        "n_fit_bins": int(len(fit_inputs["x"])),
        "message": "",
        "amplitude_mode": resolved_amplitude_mode,
        "amplitude_windows": resolved_amplitude_windows,
        "fit_diagnostics": fit_diagnostics,
    }


def _coerce_window_map(window_map):
    coerced = {}
    for component_name, windows in (window_map or {}).items():
        if windows is None:
            continue
        if isinstance(windows, tuple) and len(windows) == 2:
            coerced[component_name] = [(float(windows[0]), float(windows[1]))]
            continue

        resolved_windows = []
        for window in windows:
            if window is None or len(window) != 2:
                continue
            resolved_windows.append((float(window[0]), float(window[1])))
        if resolved_windows:
            coerced[component_name] = resolved_windows
    return coerced


def _collect_unique_windows(window_map, ordered_names=None):
    unique_windows = []
    seen = set()
    names = ordered_names or list((window_map or {}).keys())
    for component_name in names:
        for window_min, window_max in (window_map or {}).get(component_name, []):
            key = (round(float(window_min), 8), round(float(window_max), 8))
            if key in seen:
                continue
            seen.add(key)
            unique_windows.append((float(window_min), float(window_max)))
    return unique_windows


def _compute_fit_quality(
    target_hist,
    fit_hist,
    fit_min,
    fit_max,
    include_windows=None,
    exclude_windows=None,
    n_parameters=0,
):
    fit_inputs = _build_single_template_fit_inputs(
        target_hist,
        fit_hist,
        fit_min,
        fit_max,
        include_windows=include_windows,
        exclude_windows=exclude_windows,
    )
    if len(fit_inputs["x"]) == 0:
        return {
            "chi2": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "n_fit_bins": 0,
        }

    residual = fit_inputs["y"] - fit_inputs["template"]
    chi2_value = float(np.sum(np.square(residual / fit_inputs["sigma"])))
    ndf_value = int(len(fit_inputs["x"]) - max(int(n_parameters), 0))
    chi2_ndf_value = (chi2_value / ndf_value) if ndf_value > 0 else None
    fit_p_value = float(chi2_dist.sf(chi2_value, ndf_value)) if ndf_value > 0 else None
    return {
        "chi2": chi2_value,
        "ndf": ndf_value,
        "chi2_ndf": chi2_ndf_value,
        "fit_p_value": fit_p_value,
        "n_fit_bins": int(len(fit_inputs["x"])),
    }


def _build_multi_template_fit_inputs(
    target_hist,
    template_hists,
    template_names,
    fit_min,
    fit_max,
    include_windows=None,
    exclude_windows=None,
):
    if include_windows is None:
        include_windows = []
    if exclude_windows is None:
        exclude_windows = []

    x_values = []
    y_values = []
    sigma_values = []
    fit_bin_indices = []
    template_columns = {name: [] for name in template_names}
    excluded_invalid_variance_bins = []

    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < fit_min or x_center > fit_max:
            continue
        if include_windows and (
            not any(window_min <= x_center <= window_max for window_min, window_max in include_windows)
        ):
            continue
        if any(window_min <= x_center <= window_max for window_min, window_max in exclude_windows):
            continue

        y_value = float(target_hist.GetBinContent(bin_index))
        sigma_value = float(target_hist.GetBinError(bin_index))
        if not math.isfinite(y_value):
            continue
        if (not math.isfinite(sigma_value)) or sigma_value <= 0.0:
            excluded_invalid_variance_bins.append(int(bin_index))
            continue

        x_values.append(x_center)
        y_values.append(y_value)
        sigma_values.append(sigma_value)
        fit_bin_indices.append(bin_index)
        for template_name in template_names:
            template_columns[template_name].append(
                float(template_hists[template_name].GetBinContent(bin_index))
            )

    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "sigma": np.asarray(sigma_values, dtype=float),
        "fit_bin_indices": fit_bin_indices,
        "template_columns": {
            template_name: np.asarray(values, dtype=float)
            for template_name, values in template_columns.items()
        },
        "excluded_invalid_variance_bins": excluded_invalid_variance_bins,
        "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
    }


def _build_model_hist(target_hist, template_hists, amplitude_map, hist_name):
    model_hist = _clone_hist(target_hist, hist_name, reset=True)
    if model_hist is None:
        return None
    for template_name, template_hist in (template_hists or {}).items():
        amplitude = float((amplitude_map or {}).get(template_name, 0.0) or 0.0)
        if template_hist is None or amplitude == 0.0:
            continue
        model_hist.Add(template_hist, amplitude)
    return model_hist


def _hist_axis_signature(hist):
    if not _is_root_hist(hist):
        return None
    axis = hist.GetXaxis()
    return (
        int(hist.GetNbinsX()),
        float(axis.GetXmin()),
        float(axis.GetXmax()),
    )


def _hist_shape_integrity(reference_hist, comparison_hist, tolerance=1e-12):
    """Return a shape-only comparison; never use a display normalization here."""
    metrics = {
        "available": bool(reference_hist is not None and comparison_hist is not None),
        "signature_match": False,
        "max_abs_bin_difference": None,
        "shape_identical": False,
    }
    if not metrics["available"]:
        return metrics
    if _hist_axis_signature(reference_hist) != _hist_axis_signature(comparison_hist):
        return metrics
    metrics["signature_match"] = True
    max_difference = 0.0
    for bin_index in range(1, reference_hist.GetNbinsX() + 1):
        max_difference = max(
            max_difference,
            abs(
                float(reference_hist.GetBinContent(bin_index))
                - float(comparison_hist.GetBinContent(bin_index))
            ),
        )
    metrics["max_abs_bin_difference"] = float(max_difference)
    metrics["shape_identical"] = bool(max_difference <= float(tolerance))
    return metrics


def _build_hist_from_scaled_components(target_hist, scaled_components, hist_name):
    """Assemble a model from committed scaled histograms without rescaling."""
    model_hist = _clone_hist(target_hist, hist_name, reset=True)
    if model_hist is None:
        return None
    for component_hist in (scaled_components or {}).values():
        if component_hist is not None:
            model_hist.Add(component_hist)
    return model_hist


def _hist_component_closure_metrics(model_hist, scaled_components, tolerance=1e-10):
    """Check a model against the exact scaled objects used to assemble it."""
    metrics = {
        "available": bool(model_hist is not None),
        "component_count": int(len(scaled_components or {})),
        "component_names": list((scaled_components or {}).keys()),
        "signature_match": False,
        "model_integral": _hist_integral(model_hist),
        "component_integral": None,
        "integral_difference": None,
        "max_abs_bin_difference": None,
        "tolerance": float(tolerance),
        "passed": False,
    }
    if model_hist is None:
        return metrics
    component_hists = [hist for hist in (scaled_components or {}).values() if hist is not None]
    if not component_hists:
        metrics.update({
            "component_integral": 0.0,
            "integral_difference": float(metrics["model_integral"]),
            "max_abs_bin_difference": max(abs(float(model_hist.GetBinContent(index))) for index in range(1, model_hist.GetNbinsX() + 1)) if model_hist.GetNbinsX() else 0.0,
        })
        metrics["signature_match"] = True
        metrics["passed"] = bool(
            abs(float(metrics["integral_difference"])) <= float(tolerance)
            and float(metrics["max_abs_bin_difference"]) <= float(tolerance)
        )
        return metrics
    if any(_hist_axis_signature(hist) != _hist_axis_signature(model_hist) for hist in component_hists):
        return metrics

    metrics["signature_match"] = True
    component_integral = 0.0
    max_difference = 0.0
    for bin_index in range(1, model_hist.GetNbinsX() + 1):
        expected = sum(float(hist.GetBinContent(bin_index)) for hist in component_hists)
        component_integral += expected
        max_difference = max(
            max_difference,
            abs(float(model_hist.GetBinContent(bin_index)) - expected),
        )
    metrics["component_integral"] = float(component_integral)
    metrics["integral_difference"] = float(metrics["model_integral"] - component_integral)
    metrics["max_abs_bin_difference"] = float(max_difference)
    metrics["passed"] = bool(
        abs(float(metrics["integral_difference"])) <= float(tolerance)
        and float(max_difference) <= float(tolerance)
    )
    return metrics


def _scaled_template_closure_metrics(scaled_hist, template_hist, amplitude, tolerance=1e-10):
    """Verify that a protected component received its amplitude exactly once."""
    metrics = {
        "available": bool(scaled_hist is not None and template_hist is not None),
        "signature_match": False,
        "amplitude": float(amplitude or 0.0),
        "max_abs_bin_difference": None,
        "tolerance": float(tolerance),
        "passed": False,
    }
    if not metrics["available"]:
        return metrics
    if _hist_axis_signature(scaled_hist) != _hist_axis_signature(template_hist):
        return metrics
    metrics["signature_match"] = True
    max_difference = 0.0
    for bin_index in range(1, scaled_hist.GetNbinsX() + 1):
        expected = float(amplitude) * float(template_hist.GetBinContent(bin_index))
        max_difference = max(
            max_difference,
            abs(float(scaled_hist.GetBinContent(bin_index)) - expected),
        )
    metrics["max_abs_bin_difference"] = float(max_difference)
    metrics["passed"] = bool(max_difference <= float(tolerance))
    return metrics


def _resolve_pi_delta_signal_protected_fit_config(fit_config):
    defaults = {
        "enabled": True,
        "require_k_lambda_template": True,
        "require_k_sigma0_template": True,
        "allow_lambda_only_fallback": True,
        "fit_window": (1.10, 1.30),
        "lambda_gauge_window": None,
        "lambda_gauge_constraint_mode": "gaussian",
        "lambda_gauge_min_relative_uncertainty": 0.05,
        "lambda_gauge_poor_relative_uncertainty": 0.35,
        "lambda_gauge_min_fraction": None,
        "lambda_gauge_max_fraction": None,
        "lambda_gauge_min_fit_bins": 2,
        "lambda_gauge_maximum_chi2_ndf": 10.0,
        "lambda_gauge_minimum_p_value": None,
        "lambda_gauge_min_retention_fraction": 0.90,
        "maximum_chi2_ndf": 5.0,
        "minimum_p_value": 1.0e-6,
        "nonnegative_amplitudes": True,
        "failure_policy": "zero_pi_delta",
        "template_corr_warn": 0.95,
        "affects_pi_n": False,
        "affects_pi_sidis": False,
    }
    configured = (fit_config or {}).get("pi_delta_signal_protected_fit") or {}
    if not isinstance(configured, dict):
        raise ValueError("kaon_nosub.pi_delta_signal_protected_fit must be a mapping")
    resolved = deepcopy(defaults)
    resolved.update(deepcopy(configured))
    resolved["enabled"] = bool(resolved.get("enabled", True))
    resolved["allow_lambda_only_fallback"] = bool(
        resolved.get("allow_lambda_only_fallback", True)
    )
    if not bool(resolved.get("require_k_lambda_template", True)) or not bool(
        resolved.get("require_k_sigma0_template", True)
    ):
        raise ValueError(
            "pi_delta_signal_protected_fit requires both K-Lambda and K-Sigma0 templates"
        )
    if not bool(resolved.get("nonnegative_amplitudes", True)):
        raise ValueError("pi_delta_signal_protected_fit requires nonnegative_amplitudes=True")
    policy = str(resolved.get("failure_policy") or "").strip().lower()
    if policy not in {"zero_pi_delta", "error"}:
        raise ValueError(
            "Unsupported pi_delta_signal_protected_fit failure_policy '{}'".format(policy)
        )
    resolved["failure_policy"] = policy
    try:
        resolved["template_corr_warn"] = float(resolved.get("template_corr_warn", 0.95))
    except (TypeError, ValueError):
        raise ValueError("pi_delta_signal_protected_fit.template_corr_warn must be numeric")
    if not 0.0 <= float(resolved["template_corr_warn"]) <= 1.0:
        raise ValueError(
            "pi_delta_signal_protected_fit.template_corr_warn must lie in [0, 1]"
        )
    fit_window = resolved.get("fit_window")
    if fit_window is None:
        resolved["fit_window_collection"] = []
    else:
        windows = _normalize_window_collection(fit_window)
        if not windows:
            raise ValueError("pi_delta_signal_protected_fit.fit_window must contain valid bounds")
        resolved["fit_window_collection"] = windows

    lambda_gauge_window = resolved.get("lambda_gauge_window")
    if lambda_gauge_window is None:
        resolved["lambda_gauge_window_collection"] = None
    else:
        gauge_windows = _normalize_window_collection(lambda_gauge_window)
        if len(gauge_windows) != 1:
            raise ValueError(
                "pi_delta_signal_protected_fit.lambda_gauge_window must contain one valid bounds pair"
            )
        resolved["lambda_gauge_window_collection"] = gauge_windows

    constraint_mode = str(resolved.get("lambda_gauge_constraint_mode") or "").strip().lower()
    if constraint_mode not in {"gaussian", "fixed"}:
        raise ValueError(
            "Unsupported pi_delta_signal_protected_fit.lambda_gauge_constraint_mode '{}'".format(
                constraint_mode
            )
        )
    resolved["lambda_gauge_constraint_mode"] = constraint_mode

    numeric_nonnegative = (
        "lambda_gauge_min_relative_uncertainty",
        "lambda_gauge_poor_relative_uncertainty",
        "lambda_gauge_min_retention_fraction",
    )
    for key in numeric_nonnegative:
        try:
            resolved[key] = float(resolved[key])
        except (TypeError, ValueError):
            raise ValueError("pi_delta_signal_protected_fit.{} must be numeric".format(key))
        if not math.isfinite(resolved[key]) or resolved[key] < 0.0:
            raise ValueError("pi_delta_signal_protected_fit.{} must be finite and nonnegative".format(key))
    if resolved["lambda_gauge_poor_relative_uncertainty"] <= 0.0:
        raise ValueError(
            "pi_delta_signal_protected_fit.lambda_gauge_poor_relative_uncertainty must be positive"
        )
    if resolved["lambda_gauge_min_retention_fraction"] > 1.0:
        raise ValueError(
            "pi_delta_signal_protected_fit.lambda_gauge_min_retention_fraction must lie in [0, 1]"
        )

    try:
        resolved["lambda_gauge_min_fit_bins"] = int(resolved["lambda_gauge_min_fit_bins"])
    except (TypeError, ValueError):
        raise ValueError("pi_delta_signal_protected_fit.lambda_gauge_min_fit_bins must be an integer")
    if resolved["lambda_gauge_min_fit_bins"] < 2:
        raise ValueError("pi_delta_signal_protected_fit.lambda_gauge_min_fit_bins must be at least 2")

    for key in ("maximum_chi2_ndf", "lambda_gauge_maximum_chi2_ndf"):
        value = resolved.get(key)
        if value is None:
            continue
        try:
            resolved[key] = float(value)
        except (TypeError, ValueError):
            raise ValueError("pi_delta_signal_protected_fit.{} must be numeric or None".format(key))
        if not math.isfinite(resolved[key]) or resolved[key] <= 0.0:
            raise ValueError("pi_delta_signal_protected_fit.{} must be positive or None".format(key))

    for key in ("minimum_p_value", "lambda_gauge_minimum_p_value"):
        value = resolved.get(key)
        if value is None:
            continue
        try:
            resolved[key] = float(value)
        except (TypeError, ValueError):
            raise ValueError("pi_delta_signal_protected_fit.{} must be numeric or None".format(key))
        if not 0.0 <= resolved[key] <= 1.0:
            raise ValueError("pi_delta_signal_protected_fit.{} must lie in [0, 1]".format(key))

    for key in ("lambda_gauge_min_fraction", "lambda_gauge_max_fraction"):
        value = resolved.get(key)
        if value is None:
            continue
        try:
            resolved[key] = float(value)
        except (TypeError, ValueError):
            raise ValueError("pi_delta_signal_protected_fit.{} must be numeric or None".format(key))
        if not math.isfinite(resolved[key]) or resolved[key] < 0.0:
            raise ValueError("pi_delta_signal_protected_fit.{} must be finite and nonnegative".format(key))
    lower_fraction = resolved.get("lambda_gauge_min_fraction")
    upper_fraction = resolved.get("lambda_gauge_max_fraction")
    if lower_fraction is not None and upper_fraction is not None and lower_fraction > upper_fraction:
        raise ValueError(
            "pi_delta_signal_protected_fit.lambda_gauge_min_fraction cannot exceed lambda_gauge_max_fraction"
        )
    return resolved


def _resolve_protected_signal_preservation_windows(fit_config, inp_dict, phi_setting, mm_offset_data):
    """Resolve diagnostic-only signal windows from their existing owners."""
    resolved = {
        "lambda_peak": {"available": False, "source": "validation_windows.lambda_peak", "windows": []},
        "k_sigma0_signal": {"available": False, "source": "kaon_nosub.windows.k_sigma0_signal", "windows": []},
    }
    try:
        proton_config = get_proton_contamination_cleaning_config(
            inp_dict=inp_dict,
            phi_setting=phi_setting,
        ) or {}
        lambda_window = (proton_config.get("validation_windows") or {}).get("lambda_peak")
        if isinstance(lambda_window, (list, tuple)) and len(lambda_window) == 2:
            resolved["lambda_peak"] = {
                "available": True,
                "source": "validation_windows.lambda_peak",
                "windows": [(float(lambda_window[0]), float(lambda_window[1]))],
            }
    except (TypeError, ValueError):
        pass

    sigma_bounds = ((fit_config or {}).get("windows") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    try:
        sigma_windows = _normalize_window_collection(sigma_bounds)
        if sigma_windows:
            offset = float(mm_offset_data) if bool((fit_config or {}).get("apply_mm_offset_data", False)) else 0.0
            resolved["k_sigma0_signal"] = {
                "available": True,
                "source": "kaon_nosub.windows.k_sigma0_signal",
                "windows": [(float(low) + offset, float(high) + offset) for low, high in sigma_windows],
            }
    except (TypeError, ValueError):
        pass
    return resolved


def _resolve_lambda_gauge_window(protected_config, inp_dict, phi_setting):
    """Resolve one authoritative pre-delta K-Lambda anchor window."""
    configured = (protected_config or {}).get("lambda_gauge_window_collection")
    if configured:
        return {
            "available": True,
            "window": tuple(configured[0]),
            "source": "pi_delta_signal_protected_fit.lambda_gauge_window",
        }
    try:
        proton_config = get_proton_contamination_cleaning_config(
            inp_dict=inp_dict,
            phi_setting=phi_setting,
        ) or {}
        lambda_window = (proton_config.get("validation_windows") or {}).get("lambda_peak")
        if isinstance(lambda_window, (list, tuple)) and len(lambda_window) == 2:
            low, high = float(lambda_window[0]), float(lambda_window[1])
            if math.isfinite(low) and math.isfinite(high) and high > low:
                return {
                    "available": True,
                    "window": (low, high),
                    "source": "proton_cleaning.validation_windows.lambda_peak",
                }
    except (TypeError, ValueError):
        pass
    return {
        "available": False,
        "window": None,
        "source": "proton_cleaning.validation_windows.lambda_peak",
        "reason": "lambda_gauge_window_unavailable",
    }


def _integrate_hist_windows(hist, windows):
    return float(sum(_integrate_hist_range(hist, low, high) for low, high in (windows or [])))


def _fit_lambda_gauge(h_pre_delta, lambda_template, gauge_window, protected_config):
    """Fit the immutable K-Lambda shape to R_preDelta before pi-delta fitting."""
    result = {
        "status": "missing_lambda_gauge",
        "solver_success": False,
        "quality_passed": False,
        "window": list(gauge_window) if gauge_window is not None else None,
        "normalization_method": "one_component_weighted_least_squares",
        "amplitude": None,
        "amplitude_sigma": None,
        "relative_sigma": None,
        "effective_sigma": None,
        "chi2": None,
        "ndf": None,
        "chi2_ndf": None,
        "p_value": None,
        "fit_bins": 0,
        "fit_bin_indices": [],
        "excluded_invalid_variance_bins": [],
        "template_integral_window": None,
        "data_integral_window": None,
        "gauge_predicted_yield_window": None,
        "failure_reason": None,
    }
    if h_pre_delta is None or lambda_template is None or gauge_window is None:
        result.update({"status": "missing_lambda_gauge", "failure_reason": "missing_k_lambda_template_or_pre_delta"})
        return result
    window_low, window_high = gauge_window
    result["template_integral_window"] = _integrate_hist_range(lambda_template, window_low, window_high)
    result["data_integral_window"] = _integrate_hist_range(h_pre_delta, window_low, window_high)
    if not _is_finite_number(result["template_integral_window"]) or result["template_integral_window"] <= 0.0:
        result.update({"status": "non_positive_k_lambda_gauge_support", "failure_reason": "non_positive_k_lambda_gauge_support"})
        return result
    fit_inputs = _build_multi_template_fit_inputs(
        h_pre_delta,
        {KAON_SIGNAL_TEMPLATE_NAME: lambda_template},
        [KAON_SIGNAL_TEMPLATE_NAME],
        window_low,
        window_high,
        include_windows=[(window_low, window_high)],
        exclude_windows=[],
    )
    result["fit_bins"] = int(len(fit_inputs.get("x", [])))
    result["fit_bin_indices"] = list(fit_inputs.get("fit_bin_indices") or [])
    result["excluded_invalid_variance_bins"] = list(
        fit_inputs.get("excluded_invalid_variance_bins") or []
    )
    if result["fit_bins"] < int(protected_config["lambda_gauge_min_fit_bins"]):
        result.update({"status": "insufficient_lambda_gauge_bins", "failure_reason": "insufficient_lambda_gauge_bins"})
        return result
    try:
        sigma = fit_inputs["sigma"]
        template = fit_inputs["template_columns"][KAON_SIGNAL_TEMPLATE_NAME]
        target = fit_inputs["y"]
        weighted_template = template / sigma
        weighted_target = target / sigma
        denominator = float(np.dot(weighted_template, weighted_template))
        numerator = float(np.dot(weighted_template, weighted_target))
    except (KeyError, TypeError, ValueError, FloatingPointError) as exc:
        result.update({"status": "nonfinite_lambda_gauge_design", "failure_reason": str(exc)})
        return result
    if (
        not math.isfinite(denominator)
        or denominator <= 0.0
        or not math.isfinite(numerator)
        or not np.all(np.isfinite(weighted_template))
        or not np.all(np.isfinite(weighted_target))
    ):
        result.update({"status": "nonfinite_lambda_gauge_design", "failure_reason": "nonfinite_lambda_gauge_design"})
        return result
    amplitude = max(float(numerator / denominator), 0.0)
    amplitude_sigma = float(1.0 / math.sqrt(denominator))
    if not math.isfinite(amplitude) or not math.isfinite(amplitude_sigma):
        result.update({"status": "nonfinite_lambda_gauge_amplitude", "failure_reason": "nonfinite_lambda_gauge_amplitude"})
        return result
    residual = weighted_target - amplitude * weighted_template
    chi2_value = float(np.sum(np.square(residual)))
    ndf_value = int(result["fit_bins"] - 1)
    chi2_ndf = float(chi2_value / ndf_value) if ndf_value > 0 else None
    p_value = float(chi2_dist.sf(chi2_value, ndf_value)) if ndf_value > 0 else None
    floor = float(protected_config["lambda_gauge_min_relative_uncertainty"]) * amplitude
    effective_sigma = max(amplitude_sigma, floor)
    quality_passed = bool(
        ndf_value > 0
        and math.isfinite(chi2_ndf)
        and (
            protected_config.get("lambda_gauge_maximum_chi2_ndf") is None
            or chi2_ndf <= float(protected_config["lambda_gauge_maximum_chi2_ndf"])
        )
        and (
            protected_config.get("lambda_gauge_minimum_p_value") is None
            or (p_value is not None and p_value >= float(protected_config["lambda_gauge_minimum_p_value"]))
        )
    )
    result.update(
        {
            "status": "success" if quality_passed else "poor_lambda_gauge_quality",
            "solver_success": True,
            "quality_passed": quality_passed,
            "amplitude": amplitude,
            "amplitude_sigma": amplitude_sigma,
            "relative_sigma": (float(amplitude_sigma / amplitude) if amplitude > 0.0 else None),
            "effective_sigma": effective_sigma,
            "chi2": chi2_value,
            "ndf": ndf_value,
            "chi2_ndf": chi2_ndf,
            "p_value": p_value,
            "gauge_predicted_yield_window": float(amplitude * result["template_integral_window"]),
            "failure_reason": None if quality_passed else "poor_lambda_gauge_quality",
        }
    )
    return result


def _lambda_preservation_constraint(h_pre_delta, lambda_template, pi_delta_template, lambda_gauge, retention_fraction):
    """Return the pre-solver Lambda-retention bound and its explicit status."""
    result = {
        "status": "not_evaluated",
        "lambda_gauge_predicted_yield": None,
        "lambda_pre_delta_yield": None,
        "pi_delta_template_yield": None,
        "minimum_required_retention": None,
        "a_delta_max": None,
        "bound_applied": False,
        "reason": None,
    }
    window = (lambda_gauge or {}).get("window")
    amplitude = (lambda_gauge or {}).get("amplitude")
    if window is None or not _is_finite_number(amplitude):
        result.update({"status": "unavailable", "reason": "missing_lambda_gauge"})
        return result
    low, high = window
    lambda_integral = _integrate_hist_range(lambda_template, low, high)
    pre_delta_yield = _integrate_hist_range(h_pre_delta, low, high)
    delta_integral = _integrate_hist_range(pi_delta_template, low, high)
    gauge_yield = float(amplitude) * lambda_integral
    required_yield = float(retention_fraction) * gauge_yield
    result.update(
        {
            "lambda_gauge_predicted_yield": gauge_yield,
            "lambda_pre_delta_yield": pre_delta_yield,
            "pi_delta_template_yield": delta_integral,
            "minimum_required_retention": required_yield,
        }
    )
    if not all(_is_finite_number(value) for value in (gauge_yield, pre_delta_yield, delta_integral, required_yield)):
        result.update({"status": "unavailable", "reason": "nonfinite_lambda_preservation_inputs"})
        return result
    if pre_delta_yield < required_yield - 1.0e-10:
        result.update({"status": "lambda_pre_delta_deficit", "reason": "lambda_pre_delta_deficit", "a_delta_max": 0.0})
        return result
    if delta_integral <= 0.0:
        result.update({"status": "pi_delta_no_lambda_support", "reason": "pi_delta_no_lambda_support"})
        return result
    result.update(
        {
            "status": "bounded",
            "a_delta_max": max(float((pre_delta_yield - required_yield) / delta_integral), 0.0),
            "bound_applied": True,
        }
    )
    return result


def _evaluate_lambda_preservation_gate(constraint, applied_a_delta):
    """Recheck retention after the constrained solution is proposed/applied."""
    result = deepcopy(constraint or {})
    if result.get("status") == "lambda_pre_delta_deficit":
        result.update({"gate_passed": False, "gate_reason": "lambda_pre_delta_deficit"})
        return result
    if result.get("status") == "unavailable":
        result.update({"gate_passed": False, "gate_reason": result.get("reason") or "unavailable"})
        return result
    delta_template_yield = result.get("pi_delta_template_yield")
    pre_delta_yield = result.get("lambda_pre_delta_yield")
    gauge_yield = result.get("lambda_gauge_predicted_yield")
    required_yield = result.get("minimum_required_retention")
    if not all(_is_finite_number(value) for value in (delta_template_yield, pre_delta_yield, gauge_yield, required_yield)):
        result.update({"gate_passed": False, "gate_reason": "nonfinite_lambda_preservation_inputs"})
        return result
    removed = float(applied_a_delta) * float(delta_template_yield)
    after = float(pre_delta_yield) - removed
    result.update(
        {
            "lambda_pi_delta_removed_yield": removed,
            "lambda_after_delta_yield": after,
            "lambda_after_over_gauge": (float(after / gauge_yield) if gauge_yield > 0.0 else None),
            "lambda_removed_fraction_of_pre_delta": (float(removed / pre_delta_yield) if pre_delta_yield != 0.0 else None),
            "lambda_removed_fraction_of_gauge": (float(removed / gauge_yield) if gauge_yield > 0.0 else None),
            "gate_passed": bool(after >= float(required_yield) - 1.0e-10),
        }
    )
    result["gate_reason"] = "pass" if result["gate_passed"] else "lambda_preservation_rejected"
    return result


def _signal_preservation_metrics(h_pre_delta, h_after_delta, scaled_signal_hists, windows):
    metrics = {}
    for signal_name, window_payload in (windows or {}).items():
        window_list = list(window_payload.get("windows") or [])
        if not bool(window_payload.get("available")):
            metrics[signal_name] = {
                "available": False,
                "source": window_payload.get("source"),
                "reason": "authoritative window unavailable",
            }
            continue
        before = _integrate_hist_windows(h_pre_delta, window_list)
        after = _integrate_hist_windows(h_after_delta, window_list)
        delta_removed = float(before - after)
        removed_fraction = (
            float(delta_removed / before)
            if math.isfinite(before) and math.isfinite(delta_removed) and abs(before) > 0.0
            else None
        )
        fitted_signal_hist = (scaled_signal_hists or {}).get(signal_name)
        metrics[signal_name] = {
            "available": True,
            "source": window_payload.get("source"),
            "windows": [list(window) for window in window_list],
            "pre_delta_yield": before,
            "after_delta_only_yield": after,
            "pi_delta_removed_yield": delta_removed,
            "pi_delta_removed_fraction": removed_fraction,
            "fitted_signal_template_available": bool(fitted_signal_hist is not None),
            "k_lambda_fitted_yield": (
                _integrate_hist_windows(
                    (scaled_signal_hists or {}).get(KAON_SIGNAL_TEMPLATE_NAME), window_list
                )
                if (scaled_signal_hists or {}).get(KAON_SIGNAL_TEMPLATE_NAME) is not None
                else None
            ),
            "k_sigma0_fitted_yield": (
                _integrate_hist_windows(
                    (scaled_signal_hists or {}).get(KAON_SIGMA0_TEMPLATE_NAME), window_list
                )
                if (scaled_signal_hists or {}).get(KAON_SIGMA0_TEMPLATE_NAME) is not None
                else None
            ),
            "window_diagnostic_only": bool(
                signal_name == KAON_SIGMA0_TEMPLATE_NAME and fitted_signal_hist is None
            ),
        }
    return metrics


_SIGMA0_PRE_SOLVER_SOURCE_FAILURE_REASONS = frozenset(
    {
        "no_source_configured",
        "configured_path_does_not_exist",
        "root_open_failed",
        "missing_simc_tree",
        "incompatible_tree_missing_branches",
        "zero_entry_tree",
        "no_events_passed_component_shape_cuts",
        "weighted_integral_non_positive",
        "normalization_failed",
    }
)


def _resolve_sigma0_source_availability(source_diagnostics):
    """Keep external-source provenance distinct from local template support."""
    source = deepcopy(source_diagnostics) if isinstance(source_diagnostics, dict) else {}
    source_reason = str(source.get("fallback_reason") or "").strip()
    if source_reason in _SIGMA0_PRE_SOLVER_SOURCE_FAILURE_REASONS:
        status = "unavailable"
    elif source:
        status = "available" if bool(source.get("normalized", False)) else "reported"
    else:
        status = "not_reported"
    return {
        "status": status,
        "reason": source_reason or None,
        "configured": source.get("configured"),
        "requested_environment_variable": source.get("requested_environment_variable"),
        "requested_root": source.get("requested_root"),
        "resolved_root": source.get("resolved_root"),
        "source_identity": deepcopy(source.get("source_identity") or {}),
        "loader_status": {
            "normalized": source.get("normalized"),
            "path_exists": source.get("path_exists"),
            "tree_exists": source.get("tree_exists"),
            "tree_entries": source.get("tree_entries"),
        },
    }


def _resolve_sigma0_scope_template_availability(template_hist, target_signature):
    """Describe the exact Sigma0 histogram offered to one protected scope."""
    if template_hist is None:
        return {"status": "unavailable", "reason": "k_sigma0_scope_template_missing"}
    if not _hist_has_usable_support(template_hist):
        return {
            "status": "unavailable",
            "reason": "k_sigma0_scope_template_non_positive",
        }
    if _hist_axis_signature(template_hist) != target_signature:
        return {
            "status": "unavailable",
            "reason": "k_sigma0_scope_template_binning_mismatch",
        }
    return {"status": "available", "reason": None}


def _resolve_lambda_source_availability(legacy_payload, h_kaon_signal_shape):
    """Keep the mandatory Lambda source distinct from one scope's template."""
    payload = legacy_payload if isinstance(legacy_payload, dict) else {}
    canonical_reference = (
        payload.get("H_k_lambda_simc_reference")
        or payload.get("k_lambda_simc_reference_hist")
        or payload.get("H_simc_shape_k_lambda")
    )
    input_loaded = bool(
        payload.get("k_lambda_simc_input_loaded", False)
        or h_kaon_signal_shape is not None
        or canonical_reference is not None
    )
    if canonical_reference is not None:
        return {
            "status": "available",
            "reason": None,
            "input_loaded": input_loaded,
            "reference_available": True,
            "source": "immutable_k_lambda_simc_reference",
        }
    if h_kaon_signal_shape is not None:
        return {
            "status": "available",
            "reason": None,
            "input_loaded": input_loaded,
            "reference_available": False,
            "source": "scope_k_lambda_simc_template",
        }
    return {
        "status": "unavailable",
        "reason": (
            "k_lambda_reference_build_failed"
            if input_loaded
            else "k_lambda_source_unavailable"
        ),
        "input_loaded": input_loaded,
        "reference_available": False,
        "source": None,
    }


def _resolve_lambda_scope_template_availability(template_hist, target_signature):
    """Describe the exact mandatory Lambda histogram offered to one scope."""
    if template_hist is None:
        return {"status": "unavailable", "reason": "k_lambda_scope_template_missing"}
    if not _hist_has_usable_support(template_hist):
        return {
            "status": "unavailable",
            "reason": "k_lambda_scope_template_non_positive",
        }
    if _hist_axis_signature(template_hist) != target_signature:
        return {
            "status": "unavailable",
            "reason": "k_lambda_scope_template_binning_mismatch",
        }
    return {"status": "available", "reason": None}


def _resolve_lambda_constraint_fallback_window(gauge_window_payload, protected_config):
    """Prefer the Lambda gauge window, then the protected fit-window envelope."""
    gauge_window = (gauge_window_payload or {}).get("window")
    if isinstance(gauge_window, (list, tuple)) and len(gauge_window) == 2:
        return tuple(gauge_window), (gauge_window_payload or {}).get("source") or "lambda_gauge_window"
    fit_windows = list((protected_config or {}).get("fit_window_collection") or [])
    if not fit_windows:
        return None, "protected_fit_window_unavailable"
    return (
        (
            min(float(window[0]) for window in fit_windows),
            max(float(window[1]) for window in fit_windows),
        ),
        "pi_delta_signal_protected_fit.fit_window",
    )


def _apply_signal_protected_pi_delta_fit(
    legacy_payload,
    h_kaon_nosub,
    h_pi_n_shape,
    h_pi_delta_shape,
    h_pi_sidis_shape,
    h_kaon_signal_shape,
    h_kaon_sigma0_shape,
    fit_config,
    fit_min,
    fit_max,
    excluded_windows,
    inp_dict,
    phi_setting,
    mm_offset_data,
    context="",
    kaon_sigma0_source_diagnostics=None,
):
    """Replace only the final kaon pi-delta amplitude with a protected fit."""
    protected_config = _resolve_pi_delta_signal_protected_fit_config(fit_config)
    if not bool(protected_config.get("enabled")):
        if isinstance(legacy_payload, dict):
            legacy_payload.setdefault("diagnostics", {})["pi_delta_signal_protected_fit"] = {
                "enabled": False,
                "status": "disabled",
                "fit_variant": "disabled",
                "fallback_attempted": False,
                "fallback_used": False,
                "resolved_configuration": deepcopy(protected_config),
            }
        return legacy_payload

    if not isinstance(legacy_payload, dict):
        raise ValueError("protected pi-delta fit requires a kaon staged-fit payload")
    diagnostics = legacy_payload.setdefault("diagnostics", {})
    legacy_a_delta = float(legacy_payload.get("A_delta", 0.0) or 0.0)
    early_amplitudes = {
        "pi_n": float(legacy_payload.get("A_n", 0.0) or 0.0),
        "pi_sidis": float(legacy_payload.get("A_sidis", 0.0) or 0.0),
    }
    early_templates = {"pi_n": h_pi_n_shape, "pi_sidis": h_pi_sidis_shape}
    early_scaled_hists = _build_scaled_hist_map(
        early_templates,
        early_amplitudes,
        "A_protected_early",
        context or "scope",
    )
    h_pre_delta = _clone_hist(
        h_kaon_nosub,
        "H_pi_delta_protected_fit_input_{}".format(context or "scope"),
    )
    if h_pre_delta is not None:
        for scaled_hist in early_scaled_hists.values():
            if scaled_hist is not None:
                h_pre_delta.Add(scaled_hist, -1.0)

    lambda_template = legacy_payload.get("H_k_lambda_simc_reference")
    if lambda_template is None:
        lambda_template = h_kaon_signal_shape
    sigma_template = legacy_payload.get("H_k_sigma0_simc_reference")
    if sigma_template is None:
        sigma_template = h_kaon_sigma0_shape
    template_hists = {
        KAON_SIGNAL_TEMPLATE_NAME: lambda_template,
        KAON_SIGMA0_TEMPLATE_NAME: sigma_template,
        "pi_delta": h_pi_delta_shape,
    }
    all_fit_names = [KAON_SIGNAL_TEMPLATE_NAME, KAON_SIGMA0_TEMPLATE_NAME, "pi_delta"]
    fit_names = []
    failure_status = None
    failure_reason = None
    reference_integrity = _hist_shape_integrity(h_kaon_signal_shape, lambda_template)
    template_signature = _hist_axis_signature(h_pre_delta)
    sigma0_source_availability = _resolve_sigma0_source_availability(
        kaon_sigma0_source_diagnostics
    )
    lambda_source_availability = _resolve_lambda_source_availability(
        legacy_payload,
        h_kaon_signal_shape,
    )
    lambda_scope_template_availability = {
        "status": "not_evaluated",
        "reason": None,
    }
    sigma0_scope_template_availability = {
        "status": "not_evaluated",
        "reason": None,
    }
    sigma0_availability_reason = None
    fallback_attempted = False
    fallback_used = False
    fallback_reason = None
    selected_fit_variant = None
    protected_fit_attempted = False

    if h_pre_delta is None:
        failure_status = "insufficient_support"
        failure_reason = "unable to construct pre-pi-delta residual"
    elif lambda_source_availability.get("status") != "available":
        failure_status = "missing_required_template"
        failure_reason = "mandatory K-Lambda SIMC source unavailable: {}".format(
            lambda_source_availability.get("reason") or "not reported"
        )
    else:
        lambda_scope_template_availability = _resolve_lambda_scope_template_availability(
            lambda_template,
            template_signature,
        )
    if failure_status is None and lambda_scope_template_availability.get("status") != "available":
        failure_status = "missing_required_template"
        failure_reason = "mandatory K-Lambda template unavailable in this scope: {}".format(
            lambda_scope_template_availability.get("reason") or "not reported"
        )
    elif failure_status is None and not _hist_has_usable_support(h_pi_delta_shape):
        failure_status = "missing_required_template"
        failure_reason = "missing or non-positive pi_delta template"
    elif failure_status is None and _hist_axis_signature(h_pi_delta_shape) != template_signature:
        failure_status = "missing_required_template"
        failure_reason = "template binning mismatch for pi_delta"
    elif failure_status is None:
        sigma0_scope_template_availability = _resolve_sigma0_scope_template_availability(
            sigma_template,
            template_signature,
        )
        if sigma0_scope_template_availability["status"] == "available":
            # This selection is final for the invocation.  Later design or
            # solver failures must remain conservative; they never retry with
            # Sigma0 removed.
            selected_fit_variant = "lambda_sigma0_protected"
            fit_names = list(all_fit_names)
        elif bool(protected_config.get("allow_lambda_only_fallback", True)):
            selected_fit_variant = "lambda_only_protected_fallback"
            fit_names = [KAON_SIGNAL_TEMPLATE_NAME, "pi_delta"]
            fallback_attempted = True
            fallback_reason = (
                sigma0_source_availability.get("reason")
                if sigma0_source_availability.get("status") == "unavailable"
                else sigma0_scope_template_availability.get("reason")
            )
            sigma0_availability_reason = fallback_reason
        else:
            failure_status = "missing_required_template"
            sigma0_availability_reason = sigma0_scope_template_availability.get("reason")
            failure_reason = "K-Sigma0 unavailable and Lambda-only fallback is disabled: {}".format(
                sigma0_availability_reason or "not reported"
            )

    gauge_window_payload = _resolve_lambda_gauge_window(
        protected_config,
        inp_dict,
        phi_setting,
    )
    lambda_gauge = {
        "status": "not_attempted",
        "solver_success": False,
        "quality_passed": False,
        "window": gauge_window_payload.get("window"),
        "source_histogram": "R_preDelta",
        "template_source": "immutable_aligned_k_lambda_simc",
        "failure_reason": None,
    }
    if failure_status is not None:
        # Preserve the primary template/support failure while making the
        # non-authoritative gauge state explicit for payload consumers.
        lambda_gauge.update(
            {
                "status": "missing_lambda_gauge",
                "failure_reason": failure_reason or failure_status,
            }
        )
    fit_inputs = None
    solver_result = None
    coefficients = {name: 0.0 for name in all_fit_names}
    proposed_coefficients = {name: None for name in all_fit_names}
    uncertainties = {name: None for name in all_fit_names}
    matrix_diagnostics = {}
    covariance_matrix = {}
    coefficient_correlation_matrix = {}
    constraint_metrics = {
        "mode": protected_config["lambda_gauge_constraint_mode"],
        "requested_mode": protected_config["lambda_gauge_constraint_mode"],
        "applied_mode": "not_attempted",
        "source": None,
        "anchor_amplitude": None,
        "effective_sigma": None,
        "prior_chi2": None,
        "total_chi2": None,
        "total_ndf": None,
        "total_chi2_ndf": None,
        "total_p_value": None,
        "prior_augmented_covariance_matrix": {},
        "prior_augmented_coefficient_correlation_matrix": {},
    }
    fit_metrics = {
        "chi2": None,
        "ndf": None,
        "chi2_ndf": None,
        "fit_p_value": None,
        "n_fit_bins": 0,
        "n_free_spectrum_parameters": 0,
    }
    lambda_gauge_hist = None
    lambda_constraint = {
        "status": "not_attempted",
        "source": None,
        "window": None,
        "window_source": None,
        "amplitude": None,
        "effective_sigma": None,
        "requested_mode": protected_config["lambda_gauge_constraint_mode"],
        "applied_mode": "not_attempted",
        "failure_reason": None,
    }
    preservation_constraint = {"status": "not_evaluated"}
    lambda_preservation = {"gate_passed": False, "gate_reason": "not_evaluated"}
    solver_success = False
    fit_quality_passed = False
    physics_acceptance_passed = False
    free_fit_names = []
    spectrum_design = None
    spectrum_target = None
    solver_design = None
    solver_target = None
    if failure_status is None:
        if bool(gauge_window_payload.get("available")):
            lambda_gauge = _fit_lambda_gauge(
                h_pre_delta,
                lambda_template,
                gauge_window_payload["window"],
                protected_config,
            )
            lambda_gauge.update(
                {
                    "source_histogram": "R_preDelta",
                    "template_source": "immutable_aligned_k_lambda_simc",
                    "window_source": gauge_window_payload.get("source"),
                }
            )
        else:
            lambda_gauge.update(
                {
                    "status": "missing_lambda_gauge",
                    "failure_reason": (
                        gauge_window_payload.get("reason")
                        or "lambda_gauge_window_unavailable"
                    ),
                    "window_source": gauge_window_payload.get("source"),
                }
            )

        if (
            bool(lambda_gauge.get("solver_success"))
            and _is_finite_number(lambda_gauge.get("amplitude"))
            and float(lambda_gauge.get("amplitude")) > 0.0
            and _is_finite_number(lambda_gauge.get("effective_sigma"))
            and float(lambda_gauge.get("effective_sigma")) > 0.0
        ):
            if bool(lambda_gauge.get("quality_passed")):
                lambda_constraint.update(
                    {
                        "status": "success",
                        "source": "protected_lambda_gauge",
                        "window": deepcopy(lambda_gauge.get("window")),
                        "window_source": lambda_gauge.get("window_source"),
                        "amplitude": float(lambda_gauge["amplitude"]),
                        "effective_sigma": float(lambda_gauge["effective_sigma"]),
                        "applied_mode": protected_config["lambda_gauge_constraint_mode"],
                    }
                )
            else:
                lambda_constraint.update(
                    {
                        "status": "poor_gauge_softened",
                        "source": "protected_lambda_gauge_poor_quality",
                        "window": deepcopy(lambda_gauge.get("window")),
                        "window_source": lambda_gauge.get("window_source"),
                        "amplitude": float(lambda_gauge["amplitude"]),
                        "effective_sigma": max(
                            float(lambda_gauge.get("amplitude_sigma") or 0.0),
                            float(protected_config["lambda_gauge_poor_relative_uncertainty"])
                            * float(lambda_gauge["amplitude"]),
                        ),
                        "applied_mode": "gaussian_inflated",
                        "failure_reason": "poor_lambda_gauge_quality",
                    }
                )
        else:
            fallback_window, fallback_window_source = _resolve_lambda_constraint_fallback_window(
                gauge_window_payload,
                protected_config,
            )
            fallback_reference = (
                legacy_payload.get("H_k_lambda_simc_reference")
                or legacy_payload.get("k_lambda_simc_reference_hist")
                or legacy_payload.get("H_simc_shape_k_lambda")
                or lambda_template
            )
            fallback_scale = None
            fallback_source = None
            if fallback_window is not None and fallback_reference is not None:
                low, high = fallback_window
                pre_delta_yield = _integrate_hist_range(h_pre_delta, low, high)
                reference_yield = _integrate_hist_range(fallback_reference, low, high)
                if (
                    _is_finite_number(pre_delta_yield)
                    and float(pre_delta_yield) > 0.0
                    and _is_finite_number(reference_yield)
                    and float(reference_yield) > 0.0
                ):
                    _fallback_hist, fallback_scale, fallback_source = (
                        _build_scaled_reference_hist_with_fallback(
                            h_pre_delta,
                            fallback_reference,
                            low,
                            high,
                            "H_pi_delta_lambda_constraint_fallback_{}".format(
                                context or "scope"
                            ),
                        )
                    )
            if _is_finite_number(fallback_scale) and float(fallback_scale) > 0.0:
                lambda_constraint.update(
                    {
                        "status": "canonical_pre_delta_fallback",
                        "source": "canonical_k_lambda_simc_pre_delta",
                        "window": list(fallback_window),
                        "window_source": fallback_window_source,
                        "amplitude": float(fallback_scale),
                        "effective_sigma": float(
                            protected_config["lambda_gauge_poor_relative_uncertainty"]
                        ) * float(fallback_scale),
                        "applied_mode": "gaussian_inflated",
                        "failure_reason": lambda_gauge.get("failure_reason"),
                        "normalization_source": fallback_source,
                    }
                )
            else:
                failure_status = "missing_lambda_constraint"
                failure_reason = (
                    "unable to construct canonical pre-delta K-Lambda constraint: {}".format(
                        lambda_gauge.get("failure_reason")
                        or fallback_window_source
                        or "not reported"
                    )
                )

        if failure_status is None:
            lambda_gauge_hist = _build_scaled_hist_map(
                {KAON_SIGNAL_TEMPLATE_NAME: lambda_template},
                {KAON_SIGNAL_TEMPLATE_NAME: lambda_constraint["amplitude"]},
                "A_pi_delta_lambda_gauge",
                context or "scope",
            ).get(KAON_SIGNAL_TEMPLATE_NAME)
            constraint_metrics.update(
                {
                    "mode": lambda_constraint["applied_mode"],
                    "applied_mode": lambda_constraint["applied_mode"],
                    "source": lambda_constraint["source"],
                    "anchor_amplitude": lambda_constraint["amplitude"],
                    "effective_sigma": lambda_constraint["effective_sigma"],
                }
            )
    if failure_status is None:
        preservation_constraint = _lambda_preservation_constraint(
            h_pre_delta,
            lambda_template,
            h_pi_delta_shape,
            lambda_constraint,
            protected_config["lambda_gauge_min_retention_fraction"],
        )
        if preservation_constraint.get("status") == "lambda_pre_delta_deficit":
            failure_status = "lambda_pre_delta_deficit"
            failure_reason = "lambda_pre_delta_deficit"
        elif preservation_constraint.get("status") == "unavailable":
            failure_status = "lambda_preservation_rejected"
            failure_reason = preservation_constraint.get("reason")

    if failure_status is None:
        protected_fit_attempted = True
        if lambda_constraint["applied_mode"] == "fixed":
            free_fit_names = [name for name in fit_names if name != KAON_SIGNAL_TEMPLATE_NAME]
            proposed_coefficients[KAON_SIGNAL_TEMPLATE_NAME] = float(lambda_constraint["amplitude"])
        else:
            free_fit_names = list(fit_names)
        fit_inputs = _build_multi_template_fit_inputs(
            h_pre_delta,
            template_hists,
            fit_names,
            fit_min,
            fit_max,
            include_windows=protected_config.get("fit_window_collection"),
            exclude_windows=excluded_windows,
        )
        n_fit_bins = int(len(fit_inputs.get("x", [])))
        fit_metrics["n_fit_bins"] = n_fit_bins
        fit_metrics["n_free_spectrum_parameters"] = int(len(free_fit_names))
        if n_fit_bins <= len(free_fit_names):
            failure_status = "insufficient_support"
            failure_reason = "{} fit bins for {} free protected templates".format(
                n_fit_bins,
                len(free_fit_names),
            )

    if failure_status is None:
        try:
            weighted_target = fit_inputs["y"] / fit_inputs["sigma"]
            spectrum_design = np.column_stack(
                [fit_inputs["template_columns"][name] / fit_inputs["sigma"] for name in free_fit_names]
            )
            spectrum_target = np.asarray(weighted_target, dtype=float)
            if lambda_constraint["applied_mode"] == "fixed":
                spectrum_target = spectrum_target - (
                    float(lambda_constraint["amplitude"])
                    * fit_inputs["template_columns"][KAON_SIGNAL_TEMPLATE_NAME]
                    / fit_inputs["sigma"]
                )
        except (KeyError, ValueError, TypeError, FloatingPointError) as exc:
            failure_status = "insufficient_support"
            failure_reason = "unable to build protected spectrum design: {}".format(exc)
        if failure_status is None and (
            not np.all(np.isfinite(spectrum_design)) or not np.all(np.isfinite(spectrum_target))
        ):
            failure_status = "nonfinite_solution"
            failure_reason = "non-finite protected spectrum design"

    if failure_status is None:
        # Deliberately diagnose physical spectrum columns before adding a prior.
        matrix_diagnostics = _compute_template_matrix_diagnostics(spectrum_design, free_fit_names)
        effective_rank = matrix_diagnostics.get("weighted_design_effective_rank")
        condition_number = matrix_diagnostics.get("weighted_design_condition_number")
        if effective_rank != len(free_fit_names) or not _is_finite_number(condition_number):
            failure_status = "rank_deficient"
            failure_reason = "protected physical spectrum design is rank deficient or ill-defined"

    if failure_status is None:
        solver_design = spectrum_design
        solver_target = spectrum_target
        if lambda_constraint["applied_mode"] in {"gaussian", "gaussian_inflated"}:
            lambda_index = free_fit_names.index(KAON_SIGNAL_TEMPLATE_NAME)
            effective_sigma = float(lambda_constraint["effective_sigma"])
            prior_row = np.zeros(len(free_fit_names), dtype=float)
            prior_row[lambda_index] = 1.0 / effective_sigma
            solver_design = np.vstack((spectrum_design, prior_row))
            solver_target = np.concatenate(
                (
                    spectrum_target,
                    np.asarray([float(lambda_constraint["amplitude"]) / effective_sigma]),
                )
            )
        lower_bounds = np.zeros(len(free_fit_names), dtype=float)
        upper_bounds = np.full(len(free_fit_names), np.inf, dtype=float)
        if KAON_SIGNAL_TEMPLATE_NAME in free_fit_names:
            lambda_index = free_fit_names.index(KAON_SIGNAL_TEMPLATE_NAME)
            lower_fraction = protected_config.get("lambda_gauge_min_fraction")
            upper_fraction = protected_config.get("lambda_gauge_max_fraction")
            if lower_fraction is not None:
                lower_bounds[lambda_index] = float(lower_fraction) * float(lambda_constraint["amplitude"])
            if upper_fraction is not None:
                upper_bounds[lambda_index] = float(upper_fraction) * float(lambda_constraint["amplitude"])
        if (
            preservation_constraint.get("status") == "bounded"
            and "pi_delta" in free_fit_names
        ):
            upper_bounds[free_fit_names.index("pi_delta")] = float(
                preservation_constraint["a_delta_max"]
            )
            preservation_constraint["bound_applied"] = True
        try:
            solver_result = lsq_linear(
                solver_design,
                solver_target,
                bounds=(lower_bounds, upper_bounds),
                method="trf",
            )
        except Exception as exc:
            failure_status = "solver_failure"
            failure_reason = "protected constrained solver raised: {}".format(exc)
        if failure_status is None and not bool(getattr(solver_result, "success", False)):
            failure_status = "solver_failure"
            failure_reason = str(getattr(solver_result, "message", "solver did not converge"))
        if failure_status is None and not np.all(np.isfinite(solver_result.x)):
            failure_status = "nonfinite_solution"
            failure_reason = "protected constrained solver returned non-finite amplitudes"

    if failure_status is None:
        solver_success = True
        proposed_coefficients.update(
            {
                template_name: float(solver_result.x[index])
                for index, template_name in enumerate(free_fit_names)
            }
        )
        # Physical rank, condition, and correlations must describe only the
        # measured spectrum.  Retain the prior-augmented covariance separately
        # for the constrained Gaussian solution.
        covariance_matrix, coefficient_correlation_matrix, _physical_uncertainties = _compute_parameter_covariance(
            spectrum_design,
            free_fit_names,
        )
        (
            prior_augmented_covariance,
            prior_augmented_correlation,
            fitted_uncertainties,
        ) = _compute_parameter_covariance(solver_design, free_fit_names)
        uncertainties.update(fitted_uncertainties)
        spectrum_residual = spectrum_target - np.dot(spectrum_design, solver_result.x)
        spectrum_chi2 = float(np.sum(np.square(spectrum_residual)))
        spectrum_ndf = int(len(spectrum_target) - len(free_fit_names))
        fit_metrics.update(
            {
                "chi2": spectrum_chi2,
                "ndf": spectrum_ndf,
                "chi2_ndf": (float(spectrum_chi2 / spectrum_ndf) if spectrum_ndf > 0 else None),
                "fit_p_value": (float(chi2_dist.sf(spectrum_chi2, spectrum_ndf)) if spectrum_ndf > 0 else None),
            }
        )
        total_residual = solver_target - np.dot(solver_design, solver_result.x)
        total_chi2 = float(np.sum(np.square(total_residual)))
        total_ndf = int(len(solver_target) - len(free_fit_names))
        prior_chi2 = float(total_chi2 - spectrum_chi2)
        constraint_metrics.update(
            {
                "prior_chi2": (
                    prior_chi2
                    if lambda_constraint["applied_mode"] in {"gaussian", "gaussian_inflated"}
                    else 0.0
                ),
                "total_chi2": total_chi2,
                "total_ndf": total_ndf,
                "total_chi2_ndf": (float(total_chi2 / total_ndf) if total_ndf > 0 else None),
                "total_p_value": (float(chi2_dist.sf(total_chi2, total_ndf)) if total_ndf > 0 else None),
                "prior_augmented_covariance_matrix": prior_augmented_covariance,
                "prior_augmented_coefficient_correlation_matrix": prior_augmented_correlation,
            }
        )
        fit_quality_passed = bool(
            spectrum_ndf > 0
            and _is_finite_number(fit_metrics["chi2_ndf"])
            and (
                protected_config.get("maximum_chi2_ndf") is None
                or fit_metrics["chi2_ndf"] <= float(protected_config["maximum_chi2_ndf"])
            )
            and (
                protected_config.get("minimum_p_value") is None
                or fit_metrics["fit_p_value"] >= float(protected_config["minimum_p_value"])
            )
        )
        lambda_preservation = _evaluate_lambda_preservation_gate(
            preservation_constraint,
            proposed_coefficients.get("pi_delta") or 0.0,
        )
        if not fit_quality_passed:
            failure_status = "fit_quality_rejected"
            failure_reason = "protected spectrum fit quality rejected"
        elif not bool(lambda_preservation.get("gate_passed")):
            failure_status = "lambda_preservation_rejected"
            failure_reason = lambda_preservation.get("gate_reason")
        else:
            physics_acceptance_passed = True
            coefficients.update(
                {
                    name: float(proposed_coefficients[name])
                    for name in all_fit_names
                    if proposed_coefficients.get(name) is not None
                }
            )

    selected_template_hists = {
        name: template_hists[name] for name in fit_names
    }
    scaled_signal_hists = _build_scaled_hist_map(
        selected_template_hists,
        coefficients,
        "A_pi_delta_protected",
        context or "scope",
    )
    # These are the only protected component instances.  Every downstream
    # diagnostic model reuses them directly so aligned templates and fitted
    # amplitudes cannot be applied twice.
    protected_components = {
        KAON_SIGNAL_TEMPLATE_NAME: scaled_signal_hists.get(KAON_SIGNAL_TEMPLATE_NAME),
        KAON_SIGMA0_TEMPLATE_NAME: scaled_signal_hists.get(KAON_SIGMA0_TEMPLATE_NAME),
        "pi_delta": scaled_signal_hists.get("pi_delta"),
    }
    active_protected_components = {
        name: protected_components.get(name) for name in fit_names
    }
    h_protected_total = _build_hist_from_scaled_components(
        h_pre_delta,
        active_protected_components,
        "H_pi_delta_protected_fit_total_{}".format(context or "scope"),
    )
    h_protected_residual = _clone_hist(
        h_pre_delta,
        "H_pi_delta_protected_fit_residual_{}".format(context or "scope"),
    )
    if h_protected_residual is not None and h_protected_total is not None:
        h_protected_residual.Add(h_protected_total, -1.0)
    h_after_delta = _clone_hist(
        h_pre_delta,
        "H_pi_delta_protected_after_subtraction_{}".format(context or "scope"),
    )
    if h_after_delta is not None and scaled_signal_hists.get("pi_delta") is not None:
        h_after_delta.Add(scaled_signal_hists["pi_delta"], -1.0)
    h_negative_pi_delta = _clone_hist(
        scaled_signal_hists.get("pi_delta"),
        "H_pi_delta_protected_negative_component_{}".format(context or "scope"),
    )
    if h_negative_pi_delta is not None:
        h_negative_pi_delta.Scale(-1.0)

    pion_amplitudes = {
        "pi_n": early_amplitudes["pi_n"],
        "pi_delta": coefficients["pi_delta"],
        "pi_sidis": early_amplitudes["pi_sidis"],
    }
    pion_components = {
        "pi_n": early_scaled_hists.get("pi_n"),
        "pi_sidis": early_scaled_hists.get("pi_sidis"),
        "pi_delta": scaled_signal_hists.get("pi_delta"),
    }
    h_pion_only_model = _build_hist_from_scaled_components(
        h_kaon_nosub,
        pion_components,
        "H_kaon_pion_bg_signal_protected_{}".format(context or "scope"),
    )
    full_model_components = {
        "pi_n": pion_components.get("pi_n"),
        "pi_sidis": pion_components.get("pi_sidis"),
        "pi_delta": pion_components.get("pi_delta"),
        KAON_SIGNAL_TEMPLATE_NAME: protected_components.get(KAON_SIGNAL_TEMPLATE_NAME),
        KAON_SIGMA0_TEMPLATE_NAME: protected_components.get(KAON_SIGMA0_TEMPLATE_NAME),
    }
    h_full_kaon_model = _build_hist_from_scaled_components(
        h_kaon_nosub,
        full_model_components,
        "H_kaon_full_signal_protected_model_{}".format(context or "scope"),
    )
    # The physical residual below is deliberately not data minus this model:
    # it removes pion backgrounds only.  Retain the conventional full-fit
    # residual under its own explicit diagnostic contract.
    h_full_fit_residual = _clone_hist(
        h_kaon_nosub,
        "H_kaon_full_signal_protected_fit_residual_{}".format(context or "scope"),
    )
    if h_full_fit_residual is not None and h_full_kaon_model is not None:
        h_full_fit_residual.Add(h_full_kaon_model, -1.0)
    corr_warn_pairs = []
    for left_name, row in (matrix_diagnostics.get("template_correlation_matrix") or {}).items():
        for right_name, value in (row or {}).items():
            if left_name >= right_name or not _is_finite_number(value):
                continue
            if abs(float(value)) >= float(protected_config["template_corr_warn"]):
                corr_warn_pairs.append({"left": left_name, "right": right_name, "correlation": float(value)})

    preservation_windows = _resolve_protected_signal_preservation_windows(
        fit_config,
        inp_dict,
        phi_setting,
        mm_offset_data,
    )
    fit_succeeded = bool(physics_acceptance_passed)
    fit_variant = selected_fit_variant if fit_succeeded else "zero_pi_delta_failure"
    fallback_used = bool(
        fit_succeeded and selected_fit_variant == "lambda_only_protected_fallback"
    )
    sigma0_fitted = bool(
        fit_succeeded and selected_fit_variant == "lambda_sigma0_protected"
    )
    protected_signal_closure = _hist_component_closure_metrics(
        h_protected_total,
        active_protected_components,
    )
    full_kaon_model_closure = _hist_component_closure_metrics(
        h_full_kaon_model,
        full_model_components,
    )
    protected_diagnostics = {
        "enabled": True,
        "status": "success" if fit_succeeded else failure_status,
        "fit_variant": fit_variant,
        "selected_fit_variant": selected_fit_variant,
        "protected_fit_attempted": bool(protected_fit_attempted),
        "protected_fit_succeeded": bool(fit_succeeded),
        "fallback_attempted": bool(fallback_attempted),
        "fallback_used": bool(fallback_used),
        "fallback_reason": fallback_reason,
        "sigma0_availability_reason": sigma0_availability_reason,
        "sigma0_source_availability": deepcopy(sigma0_source_availability),
        "sigma0_scope_template_availability": deepcopy(
            sigma0_scope_template_availability
        ),
        "sigma0_fitted": bool(sigma0_fitted),
        "k_lambda_source_available": bool(
            lambda_source_availability.get("status") == "available"
        ),
        "k_lambda_source_availability": deepcopy(lambda_source_availability),
        "k_lambda_scope_template_available": bool(
            lambda_scope_template_availability.get("status") == "available"
        ),
        "k_lambda_scope_template_availability": deepcopy(
            lambda_scope_template_availability
        ),
        "kaon_sigma0_source_diagnostics": deepcopy(
            kaon_sigma0_source_diagnostics or {}
        ),
        "failure_reason": failure_reason,
        "failure_policy": protected_config["failure_policy"],
        "resolved_configuration": deepcopy(protected_config),
        "fit_input_stage": "kaon_nosub_after_applied_pi_n_pi_sidis",
        "fit_window": deepcopy(protected_config.get("fit_window_collection") or []),
        "lambda_gauge_window": deepcopy(lambda_gauge.get("window")),
        "excluded_windows": deepcopy(excluded_windows or []),
        "template_sources": {
            KAON_SIGNAL_TEMPLATE_NAME: "H_k_lambda_simc_reference",
            KAON_SIGMA0_TEMPLATE_NAME: "aligned_kaon_sigma0_shape",
            "pi_delta": "selected_aligned_pi_delta_shape",
        },
        "template_availability": {
            name: bool(_hist_has_usable_support(hist)) for name, hist in template_hists.items()
        },
        "template_signature": template_signature,
        "early_amplitudes_frozen": deepcopy(early_amplitudes),
        "legacy_staged_A_delta": legacy_a_delta,
        "legacy_A_delta_role": "diagnostic_only_after_alignment_selection",
        "solver_status": (
            "success" if solver_success else (failure_status or "not_attempted")
        ),
        "solver_success": bool(solver_success),
        "fit_quality_passed": bool(fit_quality_passed),
        "lambda_preservation_passed": bool(lambda_preservation.get("gate_passed")),
        "physics_acceptance_passed": bool(physics_acceptance_passed),
        "lambda_gauge_solver_success": bool(lambda_gauge.get("solver_success")),
        "lambda_gauge_quality_passed": bool(lambda_gauge.get("quality_passed")),
        "lambda_gauge_status": lambda_gauge.get("status"),
        "lambda_gauge": deepcopy(lambda_gauge),
        "lambda_constraint": deepcopy(lambda_constraint),
        "lambda_preservation": deepcopy(lambda_preservation),
        "proposed_amplitudes": deepcopy(proposed_coefficients),
        "proposed_A_delta": proposed_coefficients.get("pi_delta"),
        "applied_A_delta": float(coefficients["pi_delta"]),
        "protected_applied_amplitudes": deepcopy(pion_amplitudes),
        "signal_amplitudes": {
            KAON_SIGNAL_TEMPLATE_NAME: (
                float(coefficients[KAON_SIGNAL_TEMPLATE_NAME])
                if fit_succeeded and KAON_SIGNAL_TEMPLATE_NAME in fit_names
                else None
            ),
            KAON_SIGMA0_TEMPLATE_NAME: (
                float(coefficients[KAON_SIGMA0_TEMPLATE_NAME])
                if sigma0_fitted
                else None
            ),
        },
        "amplitude_uncertainties": deepcopy(uncertainties),
        "fit_metrics": deepcopy(fit_metrics),
        "constraint_metrics": deepcopy(constraint_metrics),
        "fit_bin_indices": list((fit_inputs or {}).get("fit_bin_indices") or []),
        "excluded_invalid_variance_bins": list((fit_inputs or {}).get("excluded_invalid_variance_bins") or []),
        "matrix_diagnostics": deepcopy(matrix_diagnostics),
        "coefficient_covariance_matrix": deepcopy(covariance_matrix),
        "coefficient_correlation_matrix": deepcopy(coefficient_correlation_matrix),
        "bound_hit_flags": {
            name: bool(
                proposed_coefficients.get(name) is not None
                and abs(float(proposed_coefficients[name])) <= 1e-12
            )
            for name in free_fit_names
        },
        "high_template_correlation_warnings": corr_warn_pairs,
        "lambda_reference_integrity": reference_integrity,
        "early_amplitudes_frozen_integrity": {
            "before": deepcopy(early_amplitudes),
            "after": deepcopy(early_amplitudes),
            "unchanged": True,
        },
        "closure": {
            "full_five_component_model": (
                full_kaon_model_closure
                if selected_fit_variant == "lambda_sigma0_protected"
                else None
            ),
            "full_four_component_lambda_only_fallback_model": (
                full_kaon_model_closure
                if selected_fit_variant == "lambda_only_protected_fallback"
                else None
            ),
            "protected_signal_model": protected_signal_closure,
            "protected_three_component_model": (
                protected_signal_closure
                if selected_fit_variant == "lambda_sigma0_protected"
                else None
            ),
            "protected_two_component_model": (
                protected_signal_closure
                if selected_fit_variant == "lambda_only_protected_fallback"
                else None
            ),
            "pion_only_model": _hist_component_closure_metrics(
                h_pion_only_model, pion_components
            ),
            "delta_only_physics_output": _hist_component_closure_metrics(
                h_after_delta,
                {
                    "pre_delta_input": h_pre_delta,
                    "negative_pi_delta": h_negative_pi_delta,
                },
            ),
        },
        "no_double_scaling": {
            KAON_SIGNAL_TEMPLATE_NAME: _scaled_template_closure_metrics(
                protected_components.get(KAON_SIGNAL_TEMPLATE_NAME),
                lambda_template,
                coefficients[KAON_SIGNAL_TEMPLATE_NAME],
            ),
            KAON_SIGMA0_TEMPLATE_NAME: _scaled_template_closure_metrics(
                protected_components.get(KAON_SIGMA0_TEMPLATE_NAME),
                sigma_template,
                coefficients[KAON_SIGMA0_TEMPLATE_NAME],
            ),
            "pi_delta": _scaled_template_closure_metrics(
                protected_components.get("pi_delta"),
                h_pi_delta_shape,
                coefficients["pi_delta"],
            ),
        },
        "signal_templates_excluded_from_subtraction_weight": True,
        "physics_output": "H_pi_delta_protected_after_subtraction",
        "fit_residual_diagnostic_only": True,
        "signal_preservation": _signal_preservation_metrics(
            h_pre_delta,
            h_after_delta,
            scaled_signal_hists,
            preservation_windows,
        ),
    }
    if failure_status is not None and protected_config["failure_policy"] == "error":
        raise RuntimeError(
            "Signal-protected pi-delta fit failed: {}".format(failure_reason or failure_status)
        )

    # Keep early pion amplitudes untouched; only the final pi-delta source is
    # replaced.  The protected signal amplitudes remain diagnostics.
    legacy_payload["A_n"] = float(early_amplitudes["pi_n"])
    legacy_payload["A_sidis"] = float(early_amplitudes["pi_sidis"])
    legacy_payload["A_delta"] = float(coefficients["pi_delta"])
    legacy_payload["S_lambda"] = (
        float(coefficients[KAON_SIGNAL_TEMPLATE_NAME])
        if fit_succeeded and KAON_SIGNAL_TEMPLATE_NAME in fit_names
        else None
    )
    legacy_payload["S_sigma0"] = (
        float(coefficients[KAON_SIGMA0_TEMPLATE_NAME]) if sigma0_fitted else None
    )
    legacy_payload["k_lambda_fit_amplitude"] = legacy_payload["S_lambda"]
    legacy_payload["lambda_gauge_amplitude"] = lambda_gauge.get("amplitude")
    legacy_payload["lambda_gauge_amplitude_sigma"] = lambda_gauge.get("amplitude_sigma")
    legacy_payload["lambda_display_scale"] = legacy_payload.get("k_lambda_reference_scale")
    legacy_payload["pi_n_scaled_hist"] = early_scaled_hists.get("pi_n")
    legacy_payload["pi_delta_scaled_hist"] = scaled_signal_hists.get("pi_delta")
    legacy_payload["pi_sidis_scaled_hist"] = early_scaled_hists.get("pi_sidis")
    legacy_payload["k_lambda_scaled_hist"] = scaled_signal_hists.get(KAON_SIGNAL_TEMPLATE_NAME)
    legacy_payload["k_sigma0_scaled_hist"] = scaled_signal_hists.get(KAON_SIGMA0_TEMPLATE_NAME)
    legacy_payload["pion_bg_fit_hist"] = h_pion_only_model
    legacy_payload["refined_pion_bg_fit_hist"] = h_pion_only_model
    legacy_payload["refined_pi_n_scaled_hist"] = early_scaled_hists.get("pi_n")
    legacy_payload["refined_pi_delta_scaled_hist"] = scaled_signal_hists.get("pi_delta")
    legacy_payload["refined_pi_sidis_scaled_hist"] = early_scaled_hists.get("pi_sidis")
    legacy_payload["refined_k_lambda_scaled_hist"] = scaled_signal_hists.get(KAON_SIGNAL_TEMPLATE_NAME)
    legacy_payload["refined_k_sigma0_scaled_hist"] = scaled_signal_hists.get(KAON_SIGMA0_TEMPLATE_NAME)
    legacy_payload["fit_hist"] = h_full_kaon_model
    legacy_payload["refined_fit_hist"] = h_full_kaon_model
    legacy_payload["full_fit_residual_hist"] = h_full_fit_residual
    legacy_payload["refined_full_fit_residual_hist"] = h_full_fit_residual
    # This is the physical pi-only subtraction output.  The selected protected
    # signal-model residual is retained below as a diagnostic only.
    legacy_payload["residual_hist"] = h_after_delta
    legacy_payload["refined_residual_hist"] = h_after_delta
    legacy_payload["H_pi_delta_protected_fit_input"] = h_pre_delta
    legacy_payload["H_pi_delta_lambda_gauge"] = lambda_gauge_hist
    legacy_payload["H_pi_delta_protected_k_lambda"] = scaled_signal_hists.get(KAON_SIGNAL_TEMPLATE_NAME)
    legacy_payload["H_pi_delta_protected_k_sigma0"] = scaled_signal_hists.get(KAON_SIGMA0_TEMPLATE_NAME)
    legacy_payload["H_pi_delta_protected_pi_delta"] = scaled_signal_hists.get("pi_delta")
    legacy_payload["H_pi_delta_protected_fit_total"] = h_protected_total
    legacy_payload["H_pi_delta_protected_fit_residual"] = h_protected_residual
    legacy_payload["H_pi_delta_protected_after_subtraction"] = h_after_delta
    diagnostics["pi_delta_signal_protected_fit"] = protected_diagnostics
    refined_amplitudes = diagnostics.setdefault("refined_amplitudes", {})
    refined_amplitudes.update(deepcopy(pion_amplitudes))
    diagnostics.setdefault("component_amplitudes", {}).update(deepcopy(pion_amplitudes))
    diagnostics["protected_sigma0_template_available"] = bool(
        sigma0_scope_template_availability.get("status") == "available"
    )
    diagnostics["signal_templates_excluded_from_subtraction_weight"] = True
    return legacy_payload


def _integrate_hist_range(hist, x_min, x_max):
    if hist is None:
        return 0.0
    integral = 0.0
    for bin_index in range(1, hist.GetNbinsX() + 1):
        x_center = float(hist.GetBinCenter(bin_index))
        if x_min <= x_center <= x_max:
            integral += float(hist.GetBinContent(bin_index))
    return float(integral)


def _build_scaled_reference_hist(target_hist, reference_hist, x_min, x_max, hist_name):
    if target_hist is None or reference_hist is None:
        return None, None
    target_integral = _integrate_hist_range(target_hist, x_min, x_max)
    reference_integral = _integrate_hist_range(reference_hist, x_min, x_max)
    if (not math.isfinite(reference_integral)) or reference_integral <= 0.0:
        return None, None
    scale_factor = float(target_integral / reference_integral)
    scaled_hist = _clone_hist(reference_hist, hist_name)
    scaled_hist.Scale(scale_factor)
    return scaled_hist, scale_factor


def _build_scaled_reference_hist_with_fallback(target_hist, reference_hist, x_min, x_max, hist_name):
    scaled_hist, scale_factor = _build_scaled_reference_hist(
        target_hist,
        reference_hist,
        x_min,
        x_max,
        hist_name,
    )
    if scaled_hist is not None:
        return scaled_hist, scale_factor, "cut-window normalized K-Lambda SIMC"

    if target_hist is None or reference_hist is None:
        return None, None, "missing K-Lambda SIMC input"

    reference_integral = _hist_integral(reference_hist)
    if (not math.isfinite(reference_integral)) or reference_integral <= 0.0:
        return None, None, "non-positive K-Lambda SIMC integral"

    target_integral = _integrate_hist_range(target_hist, x_min, x_max)
    if (not math.isfinite(target_integral)) or target_integral <= 0.0:
        target_integral = _hist_integral(target_hist)

    if (not math.isfinite(target_integral)) or target_integral <= 0.0:
        scale_factor = 1.0
        source = "unit-normalized raw K-Lambda SIMC"
    else:
        scale_factor = float(target_integral / reference_integral)
        source = "full-range normalized K-Lambda SIMC fallback"

    scaled_hist = _clone_hist(reference_hist, hist_name)
    if scaled_hist is None:
        return None, None, "failed to clone K-Lambda SIMC"
    scaled_hist.Scale(scale_factor)
    return scaled_hist, scale_factor, source


def _clone_lambda_reference_candidate(payload, hist_key, hist_name):
    if not isinstance(payload, dict):
        return None
    reference_hist = payload.get(hist_key)
    if reference_hist is None:
        return None
    return _clone_hist(reference_hist, hist_name)


def _resolve_kaon_lambda_reference_for_plot(payload, target_hist, cut_window, scope_label, hist_name):
    if not isinstance(payload, dict):
        raise RuntimeError("K-Lambda SIMC comparison requires a retained component-fit payload")

    diagnostics = payload.get("diagnostics") or {}
    protected = diagnostics.get("pi_delta_signal_protected_fit") or (
        diagnostics.get("kaon") or {}
    ).get("pi_delta_signal_protected_fit") or {}
    gauge = protected.get("lambda_gauge") or {}
    gauge_hist = payload.get("H_pi_delta_lambda_gauge")
    if (
        bool(protected.get("enabled"))
        and gauge_hist is not None
        and bool(gauge.get("solver_success"))
        and bool(gauge.get("quality_passed"))
        and _is_finite_number(gauge.get("amplitude"))
        and float(gauge.get("amplitude")) > 0.0
    ):
        return (
            _clone_hist(gauge_hist, "{}_{}".format(hist_name, str(scope_label).replace(" ", "_"))),
            gauge.get("amplitude"),
            "canonical immutable K-Lambda SIMC",
            "protected Lambda gauge",
        )

    input_loaded = bool(payload.get("k_lambda_simc_input_loaded", False))
    canonical_reference = (
        payload.get("H_k_lambda_simc_reference")
        or payload.get("k_lambda_simc_reference_hist")
    )
    if canonical_reference is not None:
        reference_hist, scale_factor, source = _build_scaled_reference_hist_with_fallback(
            target_hist,
            canonical_reference,
            float(cut_window[0]),
            float(cut_window[1]),
            "{}_{}".format(hist_name, str(scope_label).replace(" ", "_")),
        ) if (
            target_hist is not None
            and isinstance(cut_window, (list, tuple))
            and len(cut_window) == 2
            and _is_finite_number(cut_window[0])
            and _is_finite_number(cut_window[1])
        ) else (None, None, "immutable_aligned_k_lambda_simc")
        if reference_hist is None:
            reference_hist = _clone_hist(
                canonical_reference,
                "{}_{}".format(hist_name, str(scope_label).replace(" ", "_")),
            )
            scale_factor = payload.get("k_lambda_reference_scale")
        if reference_hist is not None:
            return (
                reference_hist,
                scale_factor,
                "canonical immutable K-Lambda SIMC",
                "historical {}".format(source),
            )

    simc_reference = payload.get("H_simc_shape_k_lambda")
    if simc_reference is not None:
        reference_hist, scale_factor, source = _build_scaled_reference_hist_with_fallback(
            target_hist,
            simc_reference,
            float(cut_window[0]),
            float(cut_window[1]),
            "{}_{}".format(hist_name, str(scope_label).replace(" ", "_")),
        ) if (
            target_hist is not None
            and isinstance(cut_window, (list, tuple))
            and len(cut_window) == 2
            and _is_finite_number(cut_window[0])
            and _is_finite_number(cut_window[1])
        ) else (None, None, "H_simc_shape_k_lambda")
        if reference_hist is None:
            reference_hist = _clone_hist(
                simc_reference,
                "{}_{}".format(hist_name, str(scope_label).replace(" ", "_")),
            )
            scale_factor = payload.get("k_lambda_reference_scale")
        if reference_hist is not None:
            return (
                reference_hist,
                scale_factor,
                "canonical immutable K-Lambda SIMC (H_simc_shape_k_lambda)",
                "historical {}".format(source),
            )

    if input_loaded:
        raise RuntimeError(
            "K-Lambda SIMC input was loaded but the immutable comparison reference was lost"
        )
    raise RuntimeError(
        "K-Lambda SIMC comparison is mandatory for {} but no valid reference was retained".format(
            scope_label
        )
    )


def _mask_hist_windows_inplace(hist, windows, zero_errors=True):
    if hist is None:
        return hist
    for bin_index in range(1, hist.GetNbinsX() + 1):
        x_center = float(hist.GetBinCenter(bin_index))
        if any(window_min <= x_center <= window_max for window_min, window_max in (windows or [])):
            hist.SetBinContent(bin_index, 0.0)
            if zero_errors:
                hist.SetBinError(bin_index, 0.0)
    return hist


def _masked_hist_clone(hist, hist_name, windows, zero_errors=True):
    if hist is None:
        return None
    masked = _clone_hist(hist, hist_name)
    return _mask_hist_windows_inplace(masked, windows, zero_errors=zero_errors)


def _build_scaled_hist_map(template_hists, amplitude_map, hist_name_prefix, context):
    scaled_hists = {}
    for template_name, template_hist in (template_hists or {}).items():
        scaled_hist = _clone_hist(
            template_hist,
            "{}_{}_scaled_{}".format(hist_name_prefix, template_name, context),
        )
        if scaled_hist is not None:
            scaled_hist.Scale(float((amplitude_map or {}).get(template_name, 0.0) or 0.0))
        scaled_hists[template_name] = scaled_hist
    return scaled_hists


def _rebuild_step_overlays_with_component_scales(
    target_hist,
    step_overlays,
    template_hists,
    component_scale_map,
    context="",
):
    if target_hist is None or not step_overlays:
        return step_overlays

    resolved_scale_map = {
        str(component_name): float(scale_value or 1.0)
        for component_name, scale_value in ((component_scale_map or {}).items())
    }
    if not resolved_scale_map:
        return step_overlays
    if all(abs(scale_value - 1.0) <= 1e-12 for scale_value in resolved_scale_map.values()):
        return step_overlays

    rebuilt_overlays = []
    running_model_hist = _clone_hist(
        target_hist,
        "{}_step_postfit_running_{}".format(target_hist.GetName(), context or "scope"),
        reset=True,
    )
    for step_overlay in step_overlays:
        component_name = step_overlay.get("component_name")
        template_hist = (template_hists or {}).get(component_name)
        if template_hist is None:
            rebuilt_overlays.append(step_overlay)
            continue

        raw_amplitude = float(step_overlay.get("amplitude", 0.0) or 0.0)
        scale_factor = resolved_scale_map.get(component_name, 1.0)
        scaled_amplitude = raw_amplitude * scale_factor

        baseline_before_hist = _clone_hist(
            target_hist,
            "{}_step_postfit_baseline_{}_{}".format(
                target_hist.GetName(),
                component_name,
                step_overlay.get("step_index", 0),
            ),
            reset=True,
        )
        baseline_before_hist.Add(running_model_hist)

        residual_input_hist = _clone_hist(
            target_hist,
            "{}_step_postfit_residual_{}_{}".format(
                target_hist.GetName(),
                component_name,
                step_overlay.get("step_index", 0),
            ),
        )
        residual_input_hist.Add(baseline_before_hist, -1.0)

        component_scaled_hist = _clone_hist(
            template_hist,
            "{}_step_postfit_component_{}_{}".format(
                template_hist.GetName(),
                component_name,
                step_overlay.get("step_index", 0),
            ),
        )
        component_scaled_hist.Scale(scaled_amplitude)

        cumulative_after_hist = _clone_hist(
            baseline_before_hist,
            "{}_step_postfit_cumulative_{}_{}".format(
                baseline_before_hist.GetName(),
                component_name,
                step_overlay.get("step_index", 0),
            ),
        )
        cumulative_after_hist.Add(component_scaled_hist)

        running_model_hist = _clone_hist(
            cumulative_after_hist,
            "{}_step_postfit_running_{}_{}".format(
                cumulative_after_hist.GetName(),
                component_name,
                step_overlay.get("step_index", 0),
            ),
        )

        rebuilt_overlay = dict(step_overlay)
        rebuilt_overlay["raw_amplitude"] = raw_amplitude
        rebuilt_overlay["postfit_scale_factor"] = scale_factor
        rebuilt_overlay["amplitude"] = scaled_amplitude
        rebuilt_overlay["amplitude_mode"] = step_overlay.get("amplitude_mode", "least_squares")
        rebuilt_overlay["amplitude_windows"] = deepcopy(step_overlay.get("amplitude_windows") or [])
        rebuilt_overlay["H_baseline_before"] = baseline_before_hist
        rebuilt_overlay["H_residual_input"] = residual_input_hist
        rebuilt_overlay["H_component_scaled"] = component_scaled_hist
        rebuilt_overlay["H_cumulative_after"] = cumulative_after_hist
        rebuilt_overlays.append(rebuilt_overlay)

    return rebuilt_overlays


def _apply_component_postfit_scales(
    result,
    target_hist,
    amplitude_prefix,
    postfit_scale_map,
    fit_min,
    fit_max,
    exclude_windows=None,
    validation_options=None,
    context="",
):
    if exclude_windows is None:
        exclude_windows = []
    if validation_options is None:
        validation_options = {}
    if not isinstance(result, dict):
        return result

    resolved_scale_map = {
        str(component_name): float(scale_value or 1.0)
        for component_name, scale_value in ((postfit_scale_map or {}).items())
    }
    if not resolved_scale_map:
        return result
    if all(abs(scale_value - 1.0) <= 1e-12 for scale_value in resolved_scale_map.values()):
        return result

    raw_component_amplitudes = {
        "pi_n": float(result.get("{}_n".format(amplitude_prefix), 0.0) or 0.0),
        "pi_delta": float(result.get("{}_delta".format(amplitude_prefix), 0.0) or 0.0),
        "pi_sidis": float(result.get("{}_sidis".format(amplitude_prefix), 0.0) or 0.0),
    }
    scaled_component_amplitudes = {
        component_name: raw_component_amplitudes[component_name] * resolved_scale_map[component_name]
        for component_name in COMPONENT_NAMES
    }
    result["{}_n".format(amplitude_prefix)] = float(scaled_component_amplitudes["pi_n"])
    result["{}_delta".format(amplitude_prefix)] = float(scaled_component_amplitudes["pi_delta"])
    result["{}_sidis".format(amplitude_prefix)] = float(scaled_component_amplitudes["pi_sidis"])

    extra_component_amplitudes = {
        template_name: float(value or 0.0)
        for template_name, value in ((result.get("extra_component_amplitudes") or {}).items())
    }
    raw_extra_component_amplitudes = deepcopy(extra_component_amplitudes)
    scaled_extra_component_amplitudes = {
        template_name: (
            raw_extra_component_amplitudes[template_name]
            * float(resolved_scale_map.get(template_name, 1.0) or 1.0)
        )
        for template_name in raw_extra_component_amplitudes
    }
    result["extra_component_amplitudes"] = deepcopy(scaled_extra_component_amplitudes)
    full_amplitude_map = dict(scaled_extra_component_amplitudes)
    full_amplitude_map.update(scaled_component_amplitudes)

    template_hists = result.get("template_hists") or {}
    scaled_hist_map = _build_scaled_hist_map(
        template_hists,
        full_amplitude_map,
        amplitude_prefix,
        "{}_postfit".format(context or "scope"),
    )
    result["pi_n_scaled_hist"] = scaled_hist_map.get("pi_n")
    result["pi_delta_scaled_hist"] = scaled_hist_map.get("pi_delta")
    result["pi_sidis_scaled_hist"] = scaled_hist_map.get("pi_sidis")
    result["extra_scaled_hists"] = {
        template_name: scaled_hist_map.get(template_name)
        for template_name in scaled_extra_component_amplitudes
    }

    fit_hist = _build_model_hist(
        target_hist,
        template_hists,
        full_amplitude_map,
        "{}_fit_hist_postfit_{}".format(amplitude_prefix, context or "scope"),
    )
    residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_postfit_{}".format(amplitude_prefix, context or "scope"),
    )
    if residual_hist is not None and fit_hist is not None:
        residual_hist.Add(fit_hist, -1.0)
    result["fit_hist"] = fit_hist
    result["residual_hist"] = residual_hist

    if amplitude_prefix == "A":
        result["pion_bg_fit_hist"] = _build_model_hist(
            target_hist,
            {
                component_name: template_hists.get(component_name)
                for component_name in COMPONENT_NAMES
            },
            scaled_component_amplitudes,
            "{}_pion_bg_fit_postfit_{}".format(amplitude_prefix, context or "scope"),
        )

    diagnostics = deepcopy(result.get("diagnostics") or {})
    diagnostics["postfit_component_scales"] = deepcopy(resolved_scale_map)
    diagnostics["component_amplitudes_pre_postfit_scale"] = deepcopy(raw_component_amplitudes)
    diagnostics["component_amplitudes"] = deepcopy(scaled_component_amplitudes)
    diagnostics["extra_component_amplitudes_pre_postfit_scale"] = deepcopy(raw_extra_component_amplitudes)
    diagnostics["extra_component_amplitudes"] = deepcopy(scaled_extra_component_amplitudes)
    accepted_uncertainties = deepcopy(diagnostics.get("accepted_component_uncertainties") or {})
    for component_name, scale_value in resolved_scale_map.items():
        if accepted_uncertainties.get(component_name) is None:
            continue
        accepted_uncertainties[component_name] = (
            float(accepted_uncertainties[component_name]) * float(scale_value)
        )
    diagnostics["accepted_component_uncertainties"] = accepted_uncertainties

    validation = _evaluate_model_validation(
        target_hist,
        fit_hist,
        fit_min,
        fit_max,
        n_parameters=len(template_hists),
        exclude_windows=exclude_windows,
        oversub_sigma_tolerance=validation_options.get("oversub_sigma_tolerance", 2.0),
        max_oversub_bin_count=validation_options.get("max_oversub_bin_count"),
        max_oversub_bin_fraction=validation_options.get("max_oversub_bin_fraction"),
        max_full_range_chi2_ndf=validation_options.get("max_full_range_chi2_ndf"),
        cleanup_validation_mm_max=validation_options.get("cleanup_validation_mm_max"),
    )
    diagnostics["validation"] = deepcopy(validation)
    diagnostics["chi2"] = validation.get("chi2")
    diagnostics["ndf"] = validation.get("ndf")
    diagnostics["chi2_ndf"] = validation.get("chi2_ndf")
    diagnostics["fit_p_value"] = validation.get("fit_p_value")
    diagnostics["n_fit_bins"] = validation.get("n_fit_bins")
    result["diagnostics"] = diagnostics
    result["step_overlays"] = _rebuild_step_overlays_with_component_scales(
        target_hist,
        result.get("step_overlays") or [],
        template_hists,
        resolved_scale_map,
        context=context,
    )
    return result


def _evaluate_model_validation(
    target_hist,
    fit_hist,
    fit_min,
    fit_max,
    n_parameters=0,
    exclude_windows=None,
    oversub_sigma_tolerance=2.0,
    max_oversub_bin_count=None,
    max_oversub_bin_fraction=None,
    max_full_range_chi2_ndf=None,
    cleanup_validation_mm_max=None,
):
    if exclude_windows is None:
        exclude_windows = []

    def _evaluate_region(region_min, region_max):
        quality = _compute_fit_quality(
            target_hist,
            fit_hist,
            region_min,
            region_max,
            exclude_windows=exclude_windows,
            n_parameters=n_parameters,
        )
        oversub_bin_centers = []
        total_bins = 0
        sigma_tolerance = max(float(oversub_sigma_tolerance or 0.0), 0.0)
        for bin_index in range(1, target_hist.GetNbinsX() + 1):
            x_center = float(target_hist.GetBinCenter(bin_index))
            if x_center < region_min or x_center > region_max:
                continue
            if any(window_min <= x_center <= window_max for window_min, window_max in exclude_windows):
                continue
            total_bins += 1
            data_value = float(target_hist.GetBinContent(bin_index))
            fit_value = float(fit_hist.GetBinContent(bin_index))
            sigma_value = float(target_hist.GetBinError(bin_index))
            if (not math.isfinite(sigma_value)) or sigma_value <= 0.0:
                sigma_value = max(math.sqrt(abs(data_value)), 1.0)
            if fit_value > data_value + sigma_tolerance * sigma_value:
                oversub_bin_centers.append(x_center)

        oversub_bin_count = int(len(oversub_bin_centers))
        oversub_bin_fraction = (
            float(oversub_bin_count) / float(total_bins)
            if total_bins > 0 else 0.0
        )
        return {
            "region_min": float(region_min),
            "region_max": float(region_max),
            "oversub_bin_count": oversub_bin_count,
            "oversub_bin_fraction": oversub_bin_fraction,
            "oversub_sigma_tolerance": sigma_tolerance,
            "oversub_mm_range": (
                [float(min(oversub_bin_centers)), float(max(oversub_bin_centers))]
                if oversub_bin_centers else []
            ),
            **quality,
        }

    full_range_quality = _evaluate_region(float(fit_min), float(fit_max))
    cleanup_quality = {}
    use_cleanup_region = (
        _is_finite_number(cleanup_validation_mm_max)
        and float(cleanup_validation_mm_max) > float(fit_min)
    )
    if use_cleanup_region:
        cleanup_quality = _evaluate_region(
            float(fit_min),
            min(float(cleanup_validation_mm_max), float(fit_max)),
        )
    active_quality = cleanup_quality if cleanup_quality else full_range_quality

    accepted = True
    rejection_reasons = []
    if (
        _is_finite_number(max_oversub_bin_count)
        and int(active_quality.get("oversub_bin_count", 0) or 0) > int(max_oversub_bin_count)
    ):
        accepted = False
        rejection_reasons.append(
            "oversub_bin_count {} > {}".format(
                int(active_quality.get("oversub_bin_count", 0) or 0),
                int(max_oversub_bin_count),
            )
        )
    if (
        _is_finite_number(max_oversub_bin_fraction)
        and float(active_quality.get("oversub_bin_fraction", 0.0) or 0.0) > float(max_oversub_bin_fraction)
    ):
        accepted = False
        rejection_reasons.append(
            "oversub_bin_fraction {:.3f} > {:.3f}".format(
                float(active_quality.get("oversub_bin_fraction", 0.0) or 0.0),
                float(max_oversub_bin_fraction),
            )
        )
    if (
        _is_finite_number(max_full_range_chi2_ndf)
        and _is_finite_number(active_quality.get("chi2_ndf"))
        and float(active_quality["chi2_ndf"]) > float(max_full_range_chi2_ndf)
    ):
        accepted = False
        rejection_reasons.append(
            "chi2_ndf {:.3f} > {:.3f}".format(
                float(active_quality["chi2_ndf"]),
                float(max_full_range_chi2_ndf),
            )
        )

    return {
        "accepted": bool(accepted),
        "rejection_reasons": rejection_reasons,
        "validation_region": "cleanup" if cleanup_quality else "full_range",
        "cleanup_validation_mm_max": (
            float(cleanup_quality.get("region_max"))
            if cleanup_quality else None
        ),
        "full_range": full_range_quality,
        "cleanup_region": cleanup_quality,
        **active_quality,
    }


def _build_prior_sigma_map(fit_names, stage_uncertainties, prior_scale_map=None):
    prior_sigmas = {}
    for template_name in fit_names:
        stage_sigma = float((stage_uncertainties or {}).get(template_name, 0.0) or 0.0)
        if (not math.isfinite(stage_sigma)) or stage_sigma <= 0.0:
            stage_sigma = 1.0
        prior_scale = max(
            float(((prior_scale_map or {}).get(template_name, 1.0) or 1.0)),
            1e-6,
        )
        prior_sigmas[template_name] = max(stage_sigma * prior_scale, 1e-6)
    return prior_sigmas


def _resolve_component_fit_mode_label(mode_name):
    normalized_mode = str(mode_name or "").strip().lower()
    if not normalized_mode:
        return "staged_plus_joint"
    return normalized_mode


def _extract_component_amplitude_maps(result):
    diagnostics = deepcopy((result or {}).get("diagnostics") or {})
    raw_component_amplitudes = {
        component_name: float(
            (diagnostics.get("component_amplitudes_pre_postfit_scale") or {}).get(
                component_name,
                (diagnostics.get("component_amplitudes") or {}).get(component_name, 0.0),
            ) or 0.0
        )
        for component_name in COMPONENT_NAMES
    }
    scaled_component_amplitudes = {
        component_name: float(
            (diagnostics.get("component_amplitudes") or {}).get(component_name, 0.0) or 0.0
        )
        for component_name in COMPONENT_NAMES
    }
    raw_extra_component_amplitudes = {
        str(template_name): float(value or 0.0)
        for template_name, value in (
            diagnostics.get("extra_component_amplitudes_pre_postfit_scale")
            or result.get("extra_component_amplitudes")
            or {}
        ).items()
    }
    scaled_extra_component_amplitudes = {
        str(template_name): float(value or 0.0)
        for template_name, value in (
            diagnostics.get("extra_component_amplitudes")
            or result.get("extra_component_amplitudes")
            or {}
        ).items()
    }
    return (
        raw_component_amplitudes,
        scaled_component_amplitudes,
        raw_extra_component_amplitudes,
        scaled_extra_component_amplitudes,
    )


def _build_regularization_width_map(
    fit_names,
    stage_amplitudes,
    prior_scale_map=None,
    amplitude_floor=1e-3,
):
    floor_value = max(float(amplitude_floor or 0.0), 1e-12)
    resolved = {}
    for template_name in fit_names:
        prior_scale = max(float((prior_scale_map or {}).get(template_name, 1.0) or 1.0), 1e-12)
        reference_scale = max(
            abs(float((stage_amplitudes or {}).get(template_name, 0.0) or 0.0)),
            floor_value,
        )
        resolved[template_name] = float(prior_scale * reference_scale)
    return resolved


def _apply_component_scale_map_to_coefficients(coefficients, scale_map=None):
    resolved = {}
    for template_name, coefficient_value in (coefficients or {}).items():
        scale_value = float((scale_map or {}).get(template_name, 1.0) or 1.0)
        resolved[str(template_name)] = float(coefficient_value or 0.0) * scale_value
    return resolved


def _compute_template_matrix_diagnostics(weighted_design, fit_names):
    diagnostics = {
        "weighted_design_condition_number": None,
        "weighted_design_effective_rank": None,
        "template_correlation_matrix": {},
    }
    if weighted_design is None or len(fit_names) == 0:
        return diagnostics

    try:
        diagnostics["weighted_design_condition_number"] = float(np.linalg.cond(weighted_design))
    except Exception:
        diagnostics["weighted_design_condition_number"] = None
    try:
        diagnostics["weighted_design_effective_rank"] = int(np.linalg.matrix_rank(weighted_design))
    except Exception:
        diagnostics["weighted_design_effective_rank"] = None

    try:
        gram_matrix = np.dot(weighted_design.T, weighted_design)
        template_corr = {}
        for i, left_name in enumerate(fit_names):
            left_diag = float(gram_matrix[i, i])
            row = {}
            for j, right_name in enumerate(fit_names):
                right_diag = float(gram_matrix[j, j])
                if left_diag > 0.0 and right_diag > 0.0:
                    rho = float(gram_matrix[i, j] / math.sqrt(left_diag * right_diag))
                else:
                    rho = None
                row[right_name] = rho
            template_corr[left_name] = row
        diagnostics["template_correlation_matrix"] = template_corr
    except Exception:
        diagnostics["template_correlation_matrix"] = {}
    return diagnostics


def _compute_parameter_covariance(weighted_design, fit_names):
    covariance_matrix = {}
    correlation_matrix = {}
    uncertainties = {}
    if weighted_design is None or len(fit_names) == 0:
        return covariance_matrix, correlation_matrix, uncertainties

    try:
        normal_matrix = np.dot(weighted_design.T, weighted_design)
        covariance = np.linalg.pinv(normal_matrix)
    except Exception:
        return covariance_matrix, correlation_matrix, uncertainties

    for i, left_name in enumerate(fit_names):
        variance = float(covariance[i, i])
        uncertainties[left_name] = math.sqrt(max(variance, 0.0))
        covariance_row = {}
        correlation_row = {}
        for j, right_name in enumerate(fit_names):
            covariance_value = float(covariance[i, j])
            covariance_row[right_name] = covariance_value
            left_var = float(covariance[i, i])
            right_var = float(covariance[j, j])
            if left_var > 0.0 and right_var > 0.0:
                correlation_row[right_name] = float(
                    covariance_value / math.sqrt(left_var * right_var)
                )
            else:
                correlation_row[right_name] = None
        covariance_matrix[left_name] = covariance_row
        correlation_matrix[left_name] = correlation_row
    return covariance_matrix, correlation_matrix, uncertainties


def _run_joint_template_refinement(
    target_hist,
    template_hists,
    fit_names,
    fit_min,
    fit_max,
    initial_amplitudes,
    fit_mode,
    exclude_windows=None,
    prior_scale_map=None,
    amplitude_floor=1e-3,
):
    if exclude_windows is None:
        exclude_windows = []
    fit_names = [str(name) for name in (fit_names or []) if str(name) in (template_hists or {})]
    if not fit_names:
        return {
            "success": False,
            "status": "failure",
            "message": "no active templates available for joint refinement",
            "coefficients": {},
            "uncertainties": {},
            "active_bin_count": 0,
            "excluded_invalid_variance_bin_count": 0,
            "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
            "bound_hit_flags": {},
            "regularization_enabled": False,
            "regularization_widths": {},
            "regularization_contribution": 0.0,
            "data_chi2_contribution": None,
            "total_objective": None,
            "chi2_data": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "covariance_matrix": {},
            "correlation_matrix": {},
            "template_correlation_matrix": {},
            "weighted_design_condition_number": None,
            "weighted_design_effective_rank": None,
            "fit_bin_indices": [],
        }

    fit_inputs = _build_multi_template_fit_inputs(
        target_hist,
        template_hists,
        fit_names,
        fit_min,
        fit_max,
        exclude_windows=exclude_windows,
    )
    n_fit_bins = int(len(fit_inputs["x"]))
    if n_fit_bins == 0:
        return {
            "success": False,
            "status": "failure",
            "message": "no valid full-range fit bins",
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "active_bin_count": 0,
            "excluded_invalid_variance_bin_count": int(
                len(fit_inputs.get("excluded_invalid_variance_bins") or [])
            ),
            "invalid_bin_rule": fit_inputs.get("invalid_bin_rule"),
            "bound_hit_flags": {},
            "regularization_enabled": False,
            "regularization_widths": {},
            "regularization_contribution": 0.0,
            "data_chi2_contribution": None,
            "total_objective": None,
            "chi2_data": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "covariance_matrix": {},
            "correlation_matrix": {},
            "template_correlation_matrix": {},
            "weighted_design_condition_number": None,
            "weighted_design_effective_rank": None,
            "fit_bin_indices": [],
        }

    data_values = np.asarray(fit_inputs["y"], dtype=float)
    sigma_values = np.asarray(fit_inputs["sigma"], dtype=float)
    design_columns = [
        np.asarray(fit_inputs["template_columns"].get(name), dtype=float)
        for name in fit_names
        if fit_inputs["template_columns"].get(name) is not None
    ]
    if not design_columns:
        return {
            "success": False,
            "status": "failure",
            "message": "no valid template columns available for joint refinement",
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "active_bin_count": int(n_fit_bins),
            "excluded_invalid_variance_bin_count": int(
                len(fit_inputs.get("excluded_invalid_variance_bins") or [])
            ),
            "invalid_bin_rule": fit_inputs.get("invalid_bin_rule"),
            "bound_hit_flags": {},
            "regularization_enabled": False,
            "regularization_widths": {},
            "regularization_contribution": 0.0,
            "data_chi2_contribution": None,
            "total_objective": None,
            "chi2_data": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "covariance_matrix": {},
            "correlation_matrix": {},
            "template_correlation_matrix": {},
            "weighted_design_condition_number": None,
            "weighted_design_effective_rank": None,
            "fit_bin_indices": [int(value) for value in fit_inputs.get("fit_bin_indices") or []],
        }
    design_matrix = np.column_stack(
        design_columns
    )
    weighted_design = design_matrix / sigma_values[:, None]
    initial_vector = np.asarray(
        [
            max(float((initial_amplitudes or {}).get(template_name, 0.0) or 0.0), 0.0)
            for template_name in fit_names
        ],
        dtype=float,
    )
    regularization_enabled = (
        _resolve_component_fit_mode_label(fit_mode) == "staged_plus_regularized_joint"
    )
    regularization_widths = (
        _build_regularization_width_map(
            fit_names,
            initial_amplitudes,
            prior_scale_map=prior_scale_map,
            amplitude_floor=amplitude_floor,
        )
        if regularization_enabled else {}
    )

    def _residual_vector(parameters):
        model_values = design_matrix.dot(parameters)
        residuals = (data_values - model_values) / sigma_values
        if not regularization_enabled:
            return residuals
        reg_residuals = []
        for index, template_name in enumerate(fit_names):
            tau_value = float((regularization_widths or {}).get(template_name, 0.0) or 0.0)
            if (not math.isfinite(tau_value)) or tau_value <= 0.0:
                continue
            reg_residuals.append(
                (float(parameters[index]) - float(initial_vector[index])) / tau_value
            )
        if reg_residuals:
            residuals = np.concatenate([residuals, np.asarray(reg_residuals, dtype=float)])
        return residuals

    try:
        least_squares_result = least_squares(
            _residual_vector,
            initial_vector,
            bounds=(np.zeros(len(fit_names), dtype=float), np.full(len(fit_names), np.inf, dtype=float)),
            method="trf",
        )
    except Exception as exc:
        return {
            "success": False,
            "status": "failure",
            "message": "joint refinement exception: {}".format(exc),
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "active_bin_count": int(n_fit_bins),
            "excluded_invalid_variance_bin_count": int(
                len(fit_inputs.get("excluded_invalid_variance_bins") or [])
            ),
            "invalid_bin_rule": fit_inputs.get("invalid_bin_rule"),
            "bound_hit_flags": {},
            "regularization_enabled": bool(regularization_enabled),
            "regularization_widths": deepcopy(regularization_widths),
            "regularization_contribution": None,
            "data_chi2_contribution": None,
            "total_objective": None,
            "chi2_data": None,
            "ndf": None,
            "chi2_ndf": None,
            "fit_p_value": None,
            "covariance_matrix": {},
            "correlation_matrix": {},
            "template_correlation_matrix": {},
            "weighted_design_condition_number": None,
            "weighted_design_effective_rank": None,
            "fit_bin_indices": [int(value) for value in fit_inputs.get("fit_bin_indices") or []],
        }

    parameter_vector = np.asarray(least_squares_result.x, dtype=float)
    parameter_vector = np.clip(parameter_vector, 0.0, np.inf)
    coefficients = {
        template_name: float(parameter_vector[index])
        for index, template_name in enumerate(fit_names)
    }

    model_values = design_matrix.dot(parameter_vector)
    data_residuals = (data_values - model_values) / sigma_values
    chi2_data = float(np.sum(np.square(data_residuals)))
    regularization_contribution = 0.0
    if regularization_enabled:
        for index, template_name in enumerate(fit_names):
            tau_value = float((regularization_widths or {}).get(template_name, 0.0) or 0.0)
            if (not math.isfinite(tau_value)) or tau_value <= 0.0:
                continue
            regularization_contribution += (
                (float(parameter_vector[index]) - float(initial_vector[index])) / tau_value
            ) ** 2
    total_objective = float(chi2_data + regularization_contribution)
    n_parameters = int(len(fit_names))
    ndf = int(n_fit_bins - n_parameters)
    chi2_ndf = float(chi2_data / ndf) if ndf > 0 else None
    fit_p_value = float(chi2_dist.sf(chi2_data, ndf)) if ndf > 0 else None
    bound_hit_flags = {
        template_name: bool(abs(float(parameter_vector[index])) <= 1e-10)
        for index, template_name in enumerate(fit_names)
    }
    template_diagnostics = _compute_template_matrix_diagnostics(weighted_design, fit_names)
    covariance_matrix, correlation_matrix, uncertainties = _compute_parameter_covariance(
        weighted_design,
        fit_names,
    )

    return {
        "success": bool(getattr(least_squares_result, "success", False)),
        "status": "success" if bool(getattr(least_squares_result, "success", False)) else "failure",
        "message": str(getattr(least_squares_result, "message", "")),
        "status_code": getattr(least_squares_result, "status", None),
        "coefficients": coefficients,
        "uncertainties": uncertainties,
        "active_bin_count": int(n_fit_bins),
        "excluded_invalid_variance_bin_count": int(
            len(fit_inputs.get("excluded_invalid_variance_bins") or [])
        ),
        "invalid_bin_rule": fit_inputs.get("invalid_bin_rule"),
        "bound_hit_flags": bound_hit_flags,
        "regularization_enabled": bool(regularization_enabled),
        "regularization_widths": deepcopy(regularization_widths),
        "regularization_contribution": float(regularization_contribution),
        "data_chi2_contribution": float(chi2_data),
        "total_objective": float(total_objective),
        "chi2_data": float(chi2_data),
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
        "fit_p_value": fit_p_value,
        "covariance_matrix": covariance_matrix,
        "correlation_matrix": correlation_matrix,
        "template_correlation_matrix": template_diagnostics.get("template_correlation_matrix") or {},
        "weighted_design_condition_number": template_diagnostics.get("weighted_design_condition_number"),
        "weighted_design_effective_rank": template_diagnostics.get("weighted_design_effective_rank"),
        "fit_bin_indices": [int(value) for value in fit_inputs.get("fit_bin_indices") or []],
    }


def _apply_joint_component_refinement(
    result,
    target_hist,
    amplitude_prefix,
    fit_mode,
    fit_min,
    fit_max,
    exclude_windows=None,
    prior_scale_map=None,
    postrefine_scale_map=None,
    validation_options=None,
    context="",
    template_corr_warn=0.95,
    amplitude_floor=1e-3,
):
    if exclude_windows is None:
        exclude_windows = []
    if validation_options is None:
        validation_options = {}
    if not isinstance(result, dict):
        return result

    normalized_fit_mode = _resolve_component_fit_mode_label(fit_mode)
    diagnostics = deepcopy(result.get("diagnostics") or {})
    template_hists = result.get("template_hists") or {}
    fit_names = list(diagnostics.get("fit_order") or list(template_hists.keys()))
    fit_names = [name for name in fit_names if name in template_hists]
    stage_validation = deepcopy(diagnostics.get("validation") or {})
    diagnostics["stage_validation"] = deepcopy(diagnostics.get("stage_validation") or stage_validation)
    diagnostics["fit_mode"] = normalized_fit_mode
    diagnostics["final_solution_method"] = normalized_fit_mode
    diagnostics["stage_solution_method"] = "stage_only"

    (
        raw_component_amplitudes,
        scaled_component_amplitudes,
        raw_extra_component_amplitudes,
        scaled_extra_component_amplitudes,
    ) = _extract_component_amplitude_maps(result)
    staged_scaled_amplitudes = dict(scaled_extra_component_amplitudes)
    staged_scaled_amplitudes.update(scaled_component_amplitudes)
    staged_raw_amplitudes = dict(raw_extra_component_amplitudes)
    staged_raw_amplitudes.update(raw_component_amplitudes)

    diagnostics["staged_amplitudes_raw"] = deepcopy(staged_raw_amplitudes)
    diagnostics["staged_amplitudes_scaled"] = deepcopy(staged_scaled_amplitudes)
    diagnostics["postrefine_component_scales"] = deepcopy(postrefine_scale_map or {})

    staged_scaled_hist_map = {
        "pi_n": result.get("pi_n_scaled_hist"),
        "pi_delta": result.get("pi_delta_scaled_hist"),
        "pi_sidis": result.get("pi_sidis_scaled_hist"),
    }
    # Values in this map are ROOT histograms.  The fit result remains their
    # single owner, so a shallow map copy is the intentional handoff.
    staged_scaled_hist_map.update(dict(result.get("extra_scaled_hists") or {}))

    def _retain_staged_solution(joint_refinement_summary, message_override=None):
        retained_stage_validation = deepcopy(stage_validation)
        original_fallback_reason = str(diagnostics.get("fallback_reason") or "").strip()
        original_message = str(diagnostics.get("message") or "").strip()
        stage_failure_reason = (
            original_fallback_reason
            or original_message
            or str(message_override or "").strip()
            or str((joint_refinement_summary or {}).get("message") or "").strip()
        )
        if not retained_stage_validation and stage_failure_reason:
            retained_stage_validation = {
                "accepted": False,
                "rejection_reasons": [stage_failure_reason],
            }
        diagnostics["joint_refinement"] = deepcopy(joint_refinement_summary or {})
        diagnostics["joint_refinement_status"] = str(
            (joint_refinement_summary or {}).get("status") or "not_requested"
        )
        diagnostics["refined_amplitudes_pre_postrefine_scale"] = deepcopy(staged_scaled_amplitudes)
        diagnostics["refined_amplitudes"] = deepcopy(staged_scaled_amplitudes)
        diagnostics["amplitude_shifts"] = {
            template_name: 0.0
            for template_name in fit_names
        }
        diagnostics["amplitude_shift_fractions"] = {
            template_name: 0.0
            for template_name in fit_names
        }
        diagnostics["validation"] = deepcopy(retained_stage_validation)
        diagnostics["success"] = bool(retained_stage_validation.get("accepted"))
        diagnostics["message"] = str(
            stage_failure_reason
            if not bool(diagnostics["success"]) else (
                message_override
                or (joint_refinement_summary or {}).get("message")
                or original_message
                or ""
            )
        )
        diagnostics["chi2"] = retained_stage_validation.get("chi2")
        diagnostics["ndf"] = retained_stage_validation.get("ndf")
        diagnostics["chi2_ndf"] = retained_stage_validation.get("chi2_ndf")
        diagnostics["fit_p_value"] = retained_stage_validation.get("fit_p_value")
        diagnostics["n_fit_bins"] = (joint_refinement_summary or {}).get("active_bin_count", 0)
        diagnostics["fallback_used"] = not bool(diagnostics["success"])
        diagnostics["fallback_reason"] = (
            ""
            if bool(diagnostics["success"]) else (
                stage_failure_reason
                or "staged validation rejected: {}".format(
                    "; ".join(retained_stage_validation.get("rejection_reasons") or ["unknown"])
                )
            )
        )
        result["fit_status"] = "success" if bool(diagnostics["success"]) else "fallback"
        result["diagnostics"] = diagnostics
        result["refined_fit_hist"] = result.get("fit_hist")
        result["refined_residual_hist"] = result.get("residual_hist")
        result["refined_scaled_hist_map"] = staged_scaled_hist_map
        if amplitude_prefix == "A":
            result["refined_pion_bg_fit_hist"] = result.get("pion_bg_fit_hist")
        return result

    if normalized_fit_mode == "staged_only":
        return _retain_staged_solution({
            "requested": False,
            "status": "not_requested",
            "success": False,
            "message": "fit mode is staged_only",
            "active_bin_count": 0,
            "excluded_invalid_variance_bin_count": 0,
            "invalid_bin_rule": "exclude non-finite or non-positive Sumw2 variance bins",
            "bound_hit_flags": {},
            "regularization_enabled": False,
            "regularization_widths": {},
            "regularization_contribution": 0.0,
            "data_chi2_contribution": stage_validation.get("chi2"),
            "total_objective": stage_validation.get("chi2"),
            "chi2_data": stage_validation.get("chi2"),
            "ndf": stage_validation.get("ndf"),
            "chi2_ndf": stage_validation.get("chi2_ndf"),
            "fit_p_value": stage_validation.get("fit_p_value"),
            "covariance_matrix": {},
            "correlation_matrix": {},
            "template_correlation_matrix": {},
            "weighted_design_condition_number": None,
            "weighted_design_effective_rank": None,
        })

    refinement_result = _run_joint_template_refinement(
        target_hist,
        template_hists,
        fit_names,
        fit_min,
        fit_max,
        staged_scaled_amplitudes,
        normalized_fit_mode,
        exclude_windows=exclude_windows,
        prior_scale_map=prior_scale_map,
        amplitude_floor=amplitude_floor,
    )
    diagnostics["joint_refinement"] = deepcopy(refinement_result)
    diagnostics["joint_refinement_status"] = str(refinement_result.get("status") or "failure")

    if (
        not bool(refinement_result.get("success", False))
        and bool(stage_validation.get("accepted"))
        and (
            int(refinement_result.get("active_bin_count") or 0) <= 0
            or "no active templates" in str(refinement_result.get("message") or "").lower()
            or "no valid template columns" in str(refinement_result.get("message") or "").lower()
        )
    ):
        skipped_summary = deepcopy(refinement_result)
        skipped_summary["status"] = "skipped_no_valid_bins"
        skipped_summary["success"] = False
        skipped_summary["message"] = "joint refinement skipped: {}; staged solution retained".format(
            refinement_result.get("message") or "no valid full-range fit bins"
        )
        return _retain_staged_solution(
            skipped_summary,
            message_override=skipped_summary["message"],
        )

    refined_coefficients_raw = {
        str(template_name): float(value or 0.0)
        for template_name, value in (refinement_result.get("coefficients") or {}).items()
    }
    refined_coefficients = _apply_component_scale_map_to_coefficients(
        refined_coefficients_raw,
        scale_map=postrefine_scale_map,
    )
    refined_fit_hist_pre_postrefine = _build_model_hist(
        target_hist,
        template_hists,
        refined_coefficients_raw,
        "{}_fit_hist_refined_pre_postrefine_{}".format(amplitude_prefix, context or "scope"),
    )
    refined_fit_hist = _build_model_hist(
        target_hist,
        template_hists,
        refined_coefficients,
        "{}_fit_hist_refined_{}".format(amplitude_prefix, context or "scope"),
    )
    refined_residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_refined_{}".format(amplitude_prefix, context or "scope"),
    )
    if refined_residual_hist is not None and refined_fit_hist is not None:
        refined_residual_hist.Add(refined_fit_hist, -1.0)
    refined_scaled_hist_map = _build_scaled_hist_map(
        template_hists,
        refined_coefficients,
        amplitude_prefix,
        "{}_refined".format(context or "scope"),
    )
    refined_scaled_hist_map_pre_postrefine = _build_scaled_hist_map(
        template_hists,
        refined_coefficients_raw,
        amplitude_prefix,
        "{}_refined_pre_postrefine".format(context or "scope"),
    )
    refined_pion_bg_fit_hist_pre_postrefine = _build_model_hist(
        target_hist,
        {
            component_name: template_hists.get(component_name)
            for component_name in COMPONENT_NAMES
        },
        refined_coefficients_raw,
        "{}_pion_bg_fit_refined_pre_postrefine_{}".format(amplitude_prefix, context or "scope"),
    ) if amplitude_prefix == "A" else None
    refined_pion_bg_fit_hist = _build_model_hist(
        target_hist,
        {
            component_name: template_hists.get(component_name)
            for component_name in COMPONENT_NAMES
        },
        refined_coefficients,
        "{}_pion_bg_fit_refined_{}".format(amplitude_prefix, context or "scope"),
    ) if amplitude_prefix == "A" else None

    validation = _evaluate_model_validation(
        target_hist,
        refined_fit_hist,
        fit_min,
        fit_max,
        n_parameters=len(fit_names),
        exclude_windows=exclude_windows,
        oversub_sigma_tolerance=validation_options.get("oversub_sigma_tolerance", 2.0),
        max_oversub_bin_count=validation_options.get("max_oversub_bin_count"),
        max_oversub_bin_fraction=validation_options.get("max_oversub_bin_fraction"),
        max_full_range_chi2_ndf=validation_options.get("max_full_range_chi2_ndf"),
        cleanup_validation_mm_max=validation_options.get("cleanup_validation_mm_max"),
    )

    refined_component_amplitudes = {
        component_name: float(refined_coefficients.get(component_name, 0.0) or 0.0)
        for component_name in COMPONENT_NAMES
    }
    refined_extra_component_amplitudes = {
        template_name: float(refined_coefficients.get(template_name, 0.0) or 0.0)
        for template_name in fit_names
        if template_name not in COMPONENT_NAMES
    }
    refined_full_map = dict(refined_extra_component_amplitudes)
    refined_full_map.update(refined_component_amplitudes)
    diagnostics["refined_amplitudes_pre_postrefine_scale"] = deepcopy(refined_coefficients_raw)
    diagnostics["refined_amplitudes_postrefine_scaled"] = deepcopy(refined_full_map)
    diagnostics["staged_amplitudes_postfit_scaled"] = deepcopy(staged_scaled_amplitudes)
    diagnostics["refined_amplitudes"] = deepcopy(refined_full_map)
    diagnostics["amplitude_shifts"] = {
        template_name: float(refined_full_map.get(template_name, 0.0) - staged_scaled_amplitudes.get(template_name, 0.0))
        for template_name in fit_names
    }
    diagnostics["amplitude_shift_fractions"] = {
        template_name: (
            float(
                (refined_full_map.get(template_name, 0.0) - staged_scaled_amplitudes.get(template_name, 0.0))
                / max(abs(staged_scaled_amplitudes.get(template_name, 0.0)), 1e-12)
            )
        )
        for template_name in fit_names
    }
    diagnostics["template_corr_warn"] = float(template_corr_warn)
    diagnostics["high_template_correlations"] = [
        {
            "left": left_name,
            "right": right_name,
            "rho": float(rho_value),
        }
        for left_name, row in (refinement_result.get("template_correlation_matrix") or {}).items()
        for right_name, rho_value in (row or {}).items()
        if right_name > left_name and _is_finite_number(rho_value) and abs(float(rho_value)) > float(template_corr_warn)
    ]
    diagnostics["validation"] = deepcopy(validation)
    diagnostics["success"] = bool(
        refinement_result.get("success", False) and validation.get("accepted", False)
    )
    diagnostics["message"] = str(refinement_result.get("message") or "")
    diagnostics["chi2"] = refinement_result.get("chi2_data")
    diagnostics["ndf"] = refinement_result.get("ndf")
    diagnostics["chi2_ndf"] = refinement_result.get("chi2_ndf")
    diagnostics["fit_p_value"] = refinement_result.get("fit_p_value")
    diagnostics["n_fit_bins"] = refinement_result.get("active_bin_count")
    diagnostics["fallback_used"] = not bool(diagnostics["success"])
    diagnostics["fallback_reason"] = (
        ""
        if bool(diagnostics["success"]) else (
            "joint refinement failed validation: {}".format(
                "; ".join(validation.get("rejection_reasons") or [])
            ) if validation.get("rejection_reasons") else (
                "joint refinement failed: {}".format(refinement_result.get("message") or "unknown")
            )
        )
    )
    diagnostics["component_amplitudes"] = deepcopy(refined_component_amplitudes)
    diagnostics["extra_component_amplitudes"] = deepcopy(refined_extra_component_amplitudes)
    diagnostics["accepted_component_uncertainties"] = deepcopy(
        refinement_result.get("uncertainties") or {}
    )
    accepted_uncertainties = deepcopy(diagnostics.get("accepted_component_uncertainties") or {})
    for component_name, scale_value in (postrefine_scale_map or {}).items():
        if accepted_uncertainties.get(component_name) is None:
            continue
        accepted_uncertainties[component_name] = (
            float(accepted_uncertainties[component_name]) * float(scale_value)
        )
    diagnostics["accepted_component_uncertainties"] = accepted_uncertainties

    result["fit_status"] = "success" if bool(diagnostics["success"]) else "fallback"
    result["diagnostics"] = diagnostics
    result["{}_n".format(amplitude_prefix)] = float(refined_component_amplitudes.get("pi_n", 0.0) or 0.0)
    result["{}_delta".format(amplitude_prefix)] = float(refined_component_amplitudes.get("pi_delta", 0.0) or 0.0)
    result["{}_sidis".format(amplitude_prefix)] = float(refined_component_amplitudes.get("pi_sidis", 0.0) or 0.0)
    result["extra_component_amplitudes"] = deepcopy(refined_extra_component_amplitudes)
    result["refined_fit_hist"] = refined_fit_hist
    result["refined_fit_hist_pre_postrefine"] = refined_fit_hist_pre_postrefine
    result["refined_residual_hist"] = refined_residual_hist
    result["refined_scaled_hist_map"] = refined_scaled_hist_map
    result["refined_scaled_hist_map_pre_postrefine"] = refined_scaled_hist_map_pre_postrefine
    if amplitude_prefix == "A":
        result["refined_pion_bg_fit_hist"] = refined_pion_bg_fit_hist
        result["refined_pion_bg_fit_hist_pre_postrefine"] = refined_pion_bg_fit_hist_pre_postrefine
    return result


def _run_staged_component_pass(
    target_hist,
    template_hists,
    fit_names,
    fit_min,
    fit_max,
    anchor_window_map,
    stage_amplitude_window_map,
    stage_amplitude_mode_map,
    exclude_windows,
    amplitude_seed,
    amplitude_prefix,
    pass_index,
):
    amplitude_map = {
        template_name: max(float((amplitude_seed or {}).get(template_name, 0.0) or 0.0), 0.0)
        for template_name in template_hists
    }
    uncertainty_map = {}
    pass_summary = {}
    step_overlays = []

    for component_name in fit_names:
        baseline_before_hist = _clone_hist(
            target_hist,
            "{}_baseline_before_{}_pass{}".format(
                amplitude_prefix,
                component_name,
                pass_index + 1,
            ),
            reset=True,
        )
        for other_name in fit_names:
            if other_name == component_name:
                continue
            other_hist = template_hists.get(other_name)
            other_amplitude = float(amplitude_map.get(other_name, 0.0) or 0.0)
            if other_hist is None or other_amplitude == 0.0:
                continue
            baseline_before_hist.Add(other_hist, other_amplitude)

        residual_hist = _clone_hist(
            target_hist,
            "{}_residual_{}_pass{}".format(
                amplitude_prefix,
                component_name,
                pass_index + 1,
            ),
        )
        residual_hist.Add(baseline_before_hist, -1.0)

        stage_windows = stage_amplitude_window_map.get(component_name)
        stage_mode = str(
            (stage_amplitude_mode_map or {}).get(
                component_name,
                "window_integral" if (stage_windows or []) else "least_squares",
            )
        ).strip().lower() or "least_squares"
        include_windows = anchor_window_map.get(component_name)
        amplitude_windows = stage_windows
        if stage_mode == "least_squares" and (stage_windows or []):
            include_windows = stage_windows
            amplitude_windows = None

        solve_result = _solve_nonnegative_template_amplitude(
            residual_hist,
            template_hists[component_name],
            fit_min,
            fit_max,
            include_windows=include_windows,
            exclude_windows=exclude_windows,
            amplitude_windows=amplitude_windows,
            amplitude_mode=stage_mode,
        )
        amplitude_map[component_name] = float(solve_result.get("amplitude", 0.0) or 0.0)
        uncertainty_map[component_name] = solve_result.get("sigma")
        pass_summary[component_name] = {
            "amplitude": amplitude_map[component_name],
            "sigma": solve_result.get("sigma"),
            "chi2": solve_result.get("chi2"),
            "n_fit_bins": solve_result.get("n_fit_bins"),
            "success": bool(solve_result.get("success", False)),
            "message": solve_result.get("message", ""),
            "amplitude_mode": solve_result.get("amplitude_mode", "least_squares"),
            "amplitude_windows": deepcopy(solve_result.get("amplitude_windows") or []),
        }

        component_scaled_hist = _clone_hist(
            template_hists[component_name],
            "{}_component_scaled_{}_pass{}".format(
                amplitude_prefix,
                component_name,
                pass_index + 1,
            ),
        )
        component_scaled_hist.Scale(amplitude_map[component_name])
        cumulative_after_hist = _clone_hist(
            baseline_before_hist,
            "{}_cumulative_after_{}_pass{}".format(
                amplitude_prefix,
                component_name,
                pass_index + 1,
            ),
        )
        cumulative_after_hist.Add(component_scaled_hist)
        step_overlays.append(
            {
                "pass_index": int(pass_index + 1),
                "step_index": int(len(step_overlays) + 1),
                "component_name": component_name,
                "component_label": _component_plot_label(component_name),
                "amplitude": float(amplitude_map[component_name]),
                "sigma": solve_result.get("sigma"),
                "chi2": solve_result.get("chi2"),
                "n_fit_bins": solve_result.get("n_fit_bins"),
                "anchor_windows": deepcopy(anchor_window_map.get(component_name) or []),
                "amplitude_mode": solve_result.get("amplitude_mode", "least_squares"),
                "amplitude_windows": deepcopy(solve_result.get("amplitude_windows") or []),
                "fit_diagnostics": deepcopy(solve_result.get("fit_diagnostics") or {}),
                "excluded_windows": deepcopy(exclude_windows or []),
                "H_baseline_before": baseline_before_hist,
                "H_residual_input": residual_hist,
                "H_component_template": _clone_hist(
                    template_hists[component_name],
                    "{}_component_template_{}_pass{}".format(
                        amplitude_prefix,
                        component_name,
                        pass_index + 1,
                    ),
                ),
                "H_component_scaled": component_scaled_hist,
                "H_cumulative_after": cumulative_after_hist,
            }
        )

    return amplitude_map, uncertainty_map, pass_summary, step_overlays


def _solve_joint_template_amplitudes(
    target_hist,
    template_hists,
    fit_names,
    fit_min,
    fit_max,
    stage_amplitudes,
    prior_sigmas,
    exclude_windows=None,
):
    if exclude_windows is None:
        exclude_windows = []
    fit_names = [str(name) for name in (fit_names or []) if str(name) in (template_hists or {})]
    if not fit_names:
        return {
            "success": False,
            "coefficients": {},
            "uncertainties": {},
            "prior_sigmas": deepcopy(prior_sigmas),
            "message": "no active templates available for joint refinement",
            "n_fit_bins": 0,
        }
    fit_inputs = _build_multi_template_fit_inputs(
        target_hist,
        template_hists,
        fit_names,
        fit_min,
        fit_max,
        exclude_windows=exclude_windows,
    )
    if len(fit_inputs["x"]) == 0:
        return {
            "success": False,
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "prior_sigmas": deepcopy(prior_sigmas),
            "message": "no valid full-range fit bins",
            "n_fit_bins": 0,
        }

    design_columns = [
        fit_inputs["template_columns"].get(name) / fit_inputs["sigma"]
        for name in fit_names
        if fit_inputs["template_columns"].get(name) is not None
    ]
    if not design_columns:
        return {
            "success": False,
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "prior_sigmas": deepcopy(prior_sigmas),
            "message": "no valid template columns available for joint refinement",
            "n_fit_bins": int(len(fit_inputs["x"])),
        }

    weighted_design = np.column_stack(design_columns)
    weighted_target = fit_inputs["y"] / fit_inputs["sigma"]
    prior_rows = []
    prior_targets = []
    for index, template_name in enumerate(fit_names):
        prior_sigma = float((prior_sigmas or {}).get(template_name, 0.0) or 0.0)
        if (not math.isfinite(prior_sigma)) or prior_sigma <= 0.0:
            continue
        prior_row = np.zeros(len(fit_names), dtype=float)
        prior_row[index] = 1.0 / prior_sigma
        prior_rows.append(prior_row)
        prior_targets.append(float((stage_amplitudes or {}).get(template_name, 0.0) or 0.0) / prior_sigma)

    augmented_design = weighted_design
    augmented_target = weighted_target
    if prior_rows:
        augmented_design = np.vstack([weighted_design] + prior_rows)
        augmented_target = np.concatenate(
            [weighted_target, np.asarray(prior_targets, dtype=float)]
        )

    try:
        fit_result = lsq_linear(
            augmented_design,
            augmented_target,
            bounds=(0.0, np.inf),
            method="trf",
        )
    except Exception as exc:
        return {
            "success": False,
            "coefficients": {name: 0.0 for name in fit_names},
            "uncertainties": {},
            "prior_sigmas": deepcopy(prior_sigmas),
            "message": "joint refinement exception: {}".format(exc),
            "n_fit_bins": int(len(fit_inputs["x"])),
        }

    coefficients = {
        name: max(float(value), 0.0)
        for name, value in zip(fit_names, np.asarray(fit_result.x, dtype=float))
    }
    uncertainties = {}
    try:
        normal_matrix = np.dot(augmented_design.T, augmented_design)
        covariance = np.linalg.pinv(normal_matrix)
        for index, template_name in enumerate(fit_names):
            variance = float(covariance[index, index])
            uncertainties[template_name] = math.sqrt(max(variance, 0.0))
    except Exception:
        uncertainties = {}

    return {
        "success": bool(getattr(fit_result, "success", False)),
        "coefficients": coefficients,
        "uncertainties": uncertainties,
        "prior_sigmas": deepcopy(prior_sigmas),
        "message": getattr(fit_result, "message", ""),
        "n_fit_bins": int(len(fit_inputs["x"])),
    }


def _run_coordinate_template_updates(
    target_hist,
    template_hists,
    fit_names,
    fit_min,
    fit_max,
    initial_amplitudes,
    prior_targets=None,
    prior_sigmas=None,
    max_cycles=50,
    tolerance=1e-5,
    exclude_windows=None,
):
    if exclude_windows is None:
        exclude_windows = []
    fit_inputs = _build_multi_template_fit_inputs(
        target_hist,
        template_hists,
        fit_names,
        fit_min,
        fit_max,
        exclude_windows=exclude_windows,
    )
    if len(fit_inputs["x"]) == 0:
        return {
            "success": False,
            "converged": False,
            "coefficients": {name: 0.0 for name in fit_names},
            "message": "no valid full-range fit bins",
            "n_fit_bins": 0,
            "cycles_run": 0,
            "history": [],
        }

    weights = 1.0 / np.square(fit_inputs["sigma"])
    template_columns = fit_inputs["template_columns"]
    amplitude_map = {
        template_name: max(float((initial_amplitudes or {}).get(template_name, 0.0) or 0.0), 0.0)
        for template_name in fit_names
    }
    history = []
    cycles_allowed = max(int(max_cycles or 1), 1)
    tolerance_value = max(float(tolerance or 0.0), 0.0)
    converged = False

    for cycle_index in range(cycles_allowed):
        previous_map = deepcopy(amplitude_map)
        for template_name in fit_names:
            residual = np.asarray(fit_inputs["y"], dtype=float)
            for other_name in fit_names:
                if other_name == template_name:
                    continue
                residual = residual - float(amplitude_map.get(other_name, 0.0)) * template_columns[other_name]

            template_values = template_columns[template_name]
            numerator = float(np.dot(template_values * weights, residual))
            denominator = float(np.dot(template_values * weights, template_values))
            prior_sigma = float((prior_sigmas or {}).get(template_name, 0.0) or 0.0)
            if math.isfinite(prior_sigma) and prior_sigma > 0.0:
                numerator += float((prior_targets or {}).get(template_name, 0.0) or 0.0) / (prior_sigma ** 2)
                denominator += 1.0 / (prior_sigma ** 2)
            if (not math.isfinite(denominator)) or denominator <= 0.0:
                amplitude_map[template_name] = 0.0
                continue
            amplitude_map[template_name] = max(numerator / denominator, 0.0)

        max_rel_change = 0.0
        for template_name in fit_names:
            previous_value = float(previous_map.get(template_name, 0.0) or 0.0)
            current_value = float(amplitude_map.get(template_name, 0.0) or 0.0)
            rel_change = abs(current_value - previous_value) / max(abs(previous_value), 1e-12)
            max_rel_change = max(max_rel_change, rel_change)
        history.append(
            {
                "cycle_index": int(cycle_index + 1),
                "max_relative_change": float(max_rel_change),
                "coefficients": deepcopy(amplitude_map),
            }
        )
        if max_rel_change < tolerance_value:
            converged = True
            break

    return {
        "success": bool(converged),
        "converged": bool(converged),
        "coefficients": deepcopy(amplitude_map),
        "message": (
            "coordinate refinement converged"
            if converged else "coordinate refinement reached max cycles"
        ),
        "n_fit_bins": int(len(fit_inputs["x"])),
        "cycles_run": int(len(history)),
        "history": history,
    }


def _format_shift_point_map(shift_map):
    if not isinstance(shift_map, dict) or not shift_map:
        return "none"
    return ", ".join(
        "{}={:+.4f}".format(str(component_name), float(shift_value or 0.0))
        for component_name, shift_value in shift_map.items()
    )


def _resolve_cleanup_validation_max(fit_target, mm_offset_data=0.0, inp_dict=None, phi_setting=None):
    return resolve_particle_subtraction_component_cleanup_validation_mm_max(
        fit_target,
        mm_offset_data=mm_offset_data,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
    )


def _compute_cleanup_region_metrics(
    target_hist,
    fit_hist,
    fit_min,
    cleanup_validation_mm_max,
    n_parameters=0,
    exclude_windows=None,
):
    if target_hist is None or fit_hist is None:
        return {}
    if not _is_finite_number(cleanup_validation_mm_max):
        return {}

    cleanup_max = min(float(cleanup_validation_mm_max), float(fit_hist.GetXaxis().GetXmax()))
    if cleanup_max <= float(fit_min):
        return {}

    quality = _compute_fit_quality(
        target_hist,
        fit_hist,
        float(fit_min),
        cleanup_max,
        exclude_windows=exclude_windows,
        n_parameters=n_parameters,
    )
    residual_integral = 0.0
    pull_values = []
    fit_bin_indices = []
    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < float(fit_min) or x_center > cleanup_max:
            continue
        if any(window_min <= x_center <= window_max for window_min, window_max in (exclude_windows or [])):
            continue
        data_value = float(target_hist.GetBinContent(bin_index))
        fit_value = float(fit_hist.GetBinContent(bin_index))
        residual_integral += float(data_value - fit_value)
        sigma_value = float(target_hist.GetBinError(bin_index))
        if _is_finite_number(sigma_value) and sigma_value > 0.0:
            pull_value = float((data_value - fit_value) / sigma_value)
            pull_values.append(pull_value)
            fit_bin_indices.append(int(bin_index))

    cleanup_pull_rms = None
    cleanup_max_abs_pull = None
    if pull_values:
        cleanup_pull_rms = float(np.sqrt(np.mean(np.square(np.asarray(pull_values, dtype=float)))))
        cleanup_max_abs_pull = float(max(abs(value) for value in pull_values))
    return {
        "mm_min": float(fit_min),
        "mm_max": float(cleanup_max),
        "chi2": quality.get("chi2"),
        "ndf": quality.get("ndf"),
        "chi2_ndf": quality.get("chi2_ndf"),
        "fit_p_value": quality.get("fit_p_value"),
        "n_fit_bins": quality.get("n_fit_bins"),
        "residual_integral": float(residual_integral),
        "corrected_yield": float(residual_integral),
        "pull_rms": cleanup_pull_rms,
        "max_abs_pull": cleanup_max_abs_pull,
        "fit_bin_indices": fit_bin_indices,
    }


def _build_shift_candidate_maps(shift_settings, active_component_names):
    active_names = [str(component_name) for component_name in (active_component_names or []) if str(component_name)]
    if not active_names:
        return [{}]

    mode = str((shift_settings or {}).get("mode") or "fixed").strip().lower() or "fixed"
    configured_values = deepcopy((shift_settings or {}).get("values") or {})
    configured_grid = deepcopy((shift_settings or {}).get("scan_grid") or {})
    if mode == "fixed":
        return [
            {
                component_name: float(configured_values.get(component_name, 0.0) or 0.0)
                for component_name in active_names
            }
        ]

    grid_values = []
    for component_name in active_names:
        component_grid = list(configured_grid.get(component_name) or [])
        if not component_grid:
            component_grid = [float(configured_values.get(component_name, 0.0) or 0.0)]
        grid_values.append(component_grid)

    candidate_maps = []
    for grid_point in itertools.product(*grid_values):
        candidate_maps.append(
            {
                component_name: float(shift_value or 0.0)
                for component_name, shift_value in zip(active_names, grid_point)
            }
        )
    return candidate_maps or [{}]


def _apply_residual_shift_candidate_to_templates(
    template_hists,
    shift_map,
    bounds_map=None,
    lost_integral_warn_fraction=0.01,
    hist_name_prefix="template_shift",
):
    shifted_hists = {}
    diagnostics = {}
    warnings = []
    for template_name, template_hist in (template_hists or {}).items():
        if template_hist is None:
            diagnostics[str(template_name)] = {
                "delta_mm": 0.0,
                "component_available": False,
            }
            shifted_hists[str(template_name)] = None
            continue

        delta_mm = float((shift_map or {}).get(template_name, 0.0) or 0.0)
        shifted_hist, shift_diag = build_shifted_template_histogram(
            template_hist,
            delta_mm,
            "positive_moves_peak_higher_mm",
            "{}_{}_{}".format(hist_name_prefix, str(template_name), template_hist.GetName()),
            interpolation_mode="linear",
            renormalize=True,
        )
        component_bounds = (bounds_map or {}).get(template_name)
        if component_bounds is not None:
            bound_min, bound_max = component_bounds
            shift_diag["shift_bound_hit_flag"] = bool(
                abs(delta_mm - float(bound_min)) <= 1e-12
                or abs(delta_mm - float(bound_max)) <= 1e-12
            )
        else:
            shift_diag["shift_bound_hit_flag"] = False
        shift_diag["component_available"] = True
        shifted_hists[str(template_name)] = shifted_hist
        diagnostics[str(template_name)] = shift_diag
        if float(shift_diag.get("lost_integral_fraction", 0.0) or 0.0) > float(lost_integral_warn_fraction or 0.0):
            warnings.append(
                "{} lost_integral_fraction={:.4f}".format(
                    str(template_name),
                    float(shift_diag.get("lost_integral_fraction") or 0.0),
                )
            )
    return shifted_hists, diagnostics, warnings


def _score_residual_shift_candidate(selection_metric, fit_result, cleanup_metrics):
    metric_name = str(selection_metric or "chi2_data").strip().lower() or "chi2_data"
    diagnostics = (fit_result or {}).get("diagnostics") or {}
    validation = diagnostics.get("validation") or {}

    if metric_name == "cleanup_region_chi2":
        value = (cleanup_metrics or {}).get("chi2")
        return float(value) if _is_finite_number(value) else float("inf")
    if metric_name == "cleanup_region_yield_stability":
        value = (cleanup_metrics or {}).get("corrected_yield")
        return abs(float(value)) if _is_finite_number(value) else float("inf")

    value = diagnostics.get("chi2")
    if not _is_finite_number(value):
        value = validation.get("chi2")
    return float(value) if _is_finite_number(value) else float("inf")


def _annotate_fit_result_with_residual_shift_payload(
    fit_result,
    fit_target,
    shift_settings,
    shift_payload,
):
    if not isinstance(fit_result, dict):
        return fit_result
    if not isinstance(shift_payload, dict):
        return fit_result

    diagnostics = deepcopy(fit_result.get("diagnostics") or {})
    diagnostics["residual_component_shift"] = deepcopy(shift_payload.get("summary") or {})
    fit_result["diagnostics"] = diagnostics
    fit_result["template_shift_payload"] = shift_payload
    resolved_config_summary = deepcopy(fit_result.get("resolved_config_summary") or {})
    resolved_config_summary["residual_component_shift"] = deepcopy(shift_payload.get("summary") or {})
    fit_result["resolved_config_summary"] = resolved_config_summary
    return fit_result


def _run_component_residual_shift_selection(
    fit_target,
    target_hist,
    base_component_hists,
    base_extra_template_hists,
    shift_settings,
    fit_callback,
    fit_min,
    cleanup_validation_mm_max=None,
    exclude_windows=None,
    context="",
):
    if not bool((shift_settings or {}).get("enabled", False)):
        return None

    available_component_names = [
        component_name
        for component_name in (shift_settings.get("components") or [])
        if (base_component_hists or {}).get(component_name) is not None
        or (base_extra_template_hists or {}).get(component_name) is not None
    ]
    if not available_component_names:
        return {
            "fit_result": fit_callback(base_component_hists, base_extra_template_hists),
            "summary": {
                "enabled": False,
                "mode": str(shift_settings.get("mode") or "fixed"),
                "units": str(shift_settings.get("units") or "GeV"),
                "selection_metric": str(shift_settings.get("selection_metric") or "chi2_data"),
                "requested_components": list(shift_settings.get("components") or []),
                "active_components": [],
                "selected_shift_point": {},
                "selected_shift_reason": "no available templates for configured residual shifts",
                "candidate_count": 1,
                "candidate_summaries": [],
            },
            "original_template_hists": dict(base_component_hists or {}),
            "selected_template_hists": dict(base_component_hists or {}),
            "selected_extra_template_hists": dict(base_extra_template_hists or {}),
            "selected_component_diagnostics": {},
        }

    original_template_hists = {}
    original_template_hists.update(base_component_hists or {})
    original_template_hists.update(base_extra_template_hists or {})
    candidate_shift_maps = _build_shift_candidate_maps(shift_settings, available_component_names)
    candidate_summaries = []
    best_candidate = None
    best_key = None
    selected_reason = "best selection metric"

    for candidate_index, candidate_shift_map in enumerate(candidate_shift_maps, start=1):
        combined_template_hists = {}
        combined_template_hists.update(base_component_hists or {})
        combined_template_hists.update(base_extra_template_hists or {})
        shifted_template_hists, shifted_diags, shift_warnings = _apply_residual_shift_candidate_to_templates(
            combined_template_hists,
            candidate_shift_map,
            bounds_map=shift_settings.get("bounds") or {},
            lost_integral_warn_fraction=shift_settings.get("lost_integral_warn_fraction", 0.01),
            hist_name_prefix="{}_{}".format(fit_target, context or "scope"),
        )
        shifted_component_hists = {
            template_name: shifted_template_hists.get(template_name)
            for template_name in (base_component_hists or {}).keys()
        }
        shifted_extra_template_hists = {
            template_name: shifted_template_hists.get(template_name)
            for template_name in (base_extra_template_hists or {}).keys()
        }
        fit_result = fit_callback(shifted_component_hists, shifted_extra_template_hists)
        diagnostics = (fit_result or {}).get("diagnostics") or {}
        validation = diagnostics.get("validation") or {}
        model_hist = (
            fit_result.get("refined_fit_hist")
            or fit_result.get("fit_hist")
            or fit_result.get("refined_pion_bg_fit_hist")
            or fit_result.get("pion_bg_fit_hist")
        )
        cleanup_metrics = _compute_cleanup_region_metrics(
            target_hist,
            model_hist,
            fit_min,
            cleanup_validation_mm_max,
            n_parameters=len((fit_result.get("template_hists") or {}).keys()),
            exclude_windows=exclude_windows,
        )
        score_value = _score_residual_shift_candidate(
            shift_settings.get("selection_metric"),
            fit_result,
            cleanup_metrics,
        )
        candidate_summary = {
            "candidate_index": int(candidate_index),
            "shift_point": deepcopy(candidate_shift_map),
            "score": float(score_value) if _is_finite_number(score_value) else None,
            "fit_status": fit_result.get("fit_status"),
            "accepted": bool(validation.get("accepted")),
            "chi2": diagnostics.get("chi2"),
            "chi2_ndf": diagnostics.get("chi2_ndf"),
            "cleanup_metrics": deepcopy(cleanup_metrics),
            "warnings": list(shift_warnings),
        }
        candidate_summaries.append(candidate_summary)
        candidate_key = (
            0 if bool(validation.get("accepted")) else 1,
            float(score_value) if _is_finite_number(score_value) else float("inf"),
            int(candidate_index),
        )
        if best_candidate is None or candidate_key < best_key:
            best_candidate = {
                "fit_result": fit_result,
                "shift_point": deepcopy(candidate_shift_map),
                "shifted_template_hists": shifted_template_hists,
                "shifted_component_hists": shifted_component_hists,
                "shifted_extra_template_hists": shifted_extra_template_hists,
                "shifted_component_diagnostics": deepcopy(shifted_diags),
                "cleanup_metrics": deepcopy(cleanup_metrics),
                "warnings": list(shift_warnings),
                "accepted": bool(validation.get("accepted")),
            }
            best_key = candidate_key

    if best_candidate is None:
        return None
    if not bool(best_candidate.get("accepted")):
        selected_reason = "no accepted residual-shift candidate; retained best-scoring candidate"

    summary = {
        "enabled": True,
        "mode": str(shift_settings.get("mode") or "fixed"),
        "units": str(shift_settings.get("units") or "GeV"),
        "selection_metric": str(shift_settings.get("selection_metric") or "chi2_data"),
        "requested_components": list(shift_settings.get("components") or []),
        "active_components": list(available_component_names),
        "configured_shift_values": deepcopy(shift_settings.get("values") or {}),
        "shift_bounds": deepcopy(shift_settings.get("bounds") or {}),
        "shift_grid": deepcopy(shift_settings.get("scan_grid") or {}),
        "selected_shift_point": deepcopy(best_candidate.get("shift_point") or {}),
        "selected_shift_reason": selected_reason,
        "candidate_count": int(len(candidate_summaries)),
        "candidate_summaries": candidate_summaries,
        "per_component": deepcopy(best_candidate.get("shifted_component_diagnostics") or {}),
        "cleanup_metrics": deepcopy(best_candidate.get("cleanup_metrics") or {}),
        "warnings": list(best_candidate.get("warnings") or []),
    }
    return {
        "fit_result": best_candidate.get("fit_result"),
        "summary": summary,
        "original_template_hists": original_template_hists,
        "selected_template_hists": dict(best_candidate.get("shifted_component_hists") or {}),
        "selected_extra_template_hists": dict(best_candidate.get("shifted_extra_template_hists") or {}),
        "selected_component_diagnostics": deepcopy(best_candidate.get("shifted_component_diagnostics") or {}),
    }


def _fit_staged_anchor_templates(
    target_hist,
    component_hists,
    amplitude_prefix,
    fit_min,
    fit_max,
    anchor_windows,
    fit_order,
    stage_amplitude_windows=None,
    stage_amplitude_modes=None,
    exclude_windows=None,
    extra_positive_templates=None,
    extra_anchor_windows=None,
    n_passes=1,
    prior_scale_map=None,
    joint_refinement_enabled=True,
    max_fit_cycles=50,
    fit_tolerance=1e-5,
    validation_options=None,
    context="",
):
    if exclude_windows is None:
        exclude_windows = []
    if extra_positive_templates is None:
        extra_positive_templates = {}
    if extra_anchor_windows is None:
        extra_anchor_windows = {}
    if prior_scale_map is None:
        prior_scale_map = {}
    if validation_options is None:
        validation_options = {}

    anchor_window_map = _coerce_window_map(anchor_windows)
    stage_amplitude_window_map = _coerce_window_map(stage_amplitude_windows)
    stage_amplitude_mode_map = {
        str(component_name): str(mode_name or "").strip().lower()
        for component_name, mode_name in ((stage_amplitude_modes or {}).items())
    }
    extra_window_map = _coerce_window_map(extra_anchor_windows)
    template_hists = {
        "pi_n": component_hists["pi_n"],
        "pi_delta": component_hists["pi_delta"],
        "pi_sidis": component_hists["pi_sidis"],
    }
    template_hists.update(extra_positive_templates)

    ordered_fit_names = []
    for component_name in fit_order or COMPONENT_NAMES:
        if component_name in template_hists and component_name not in ordered_fit_names:
            ordered_fit_names.append(component_name)
    for template_name in extra_positive_templates:
        if template_name not in ordered_fit_names:
            ordered_fit_names.append(template_name)

    combined_window_map = deepcopy(anchor_window_map)
    combined_window_map.update(extra_window_map)
    fitted_names = [
        component_name for component_name in ordered_fit_names
        if (combined_window_map.get(component_name) or [])
    ]
    if not fitted_names:
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "no component anchor windows configured",
            template_hists=template_hists,
        )

    fallback_reason = _validate_component_shapes(
        component_hists,
        target_hist,
        component_names=fitted_names,
    )
    if fallback_reason:
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            fallback_reason,
            template_hists=template_hists,
        )
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(
                target_hist,
                amplitude_prefix,
                context,
                fallback_reason,
                template_hists=template_hists,
            )

    stage_amplitudes = {template_name: 0.0 for template_name in template_hists}
    stage_uncertainties = {template_name: None for template_name in template_hists}
    staged_pass_history = []
    step_overlays = []
    total_passes = max(int(n_passes or 1), 1)
    for pass_index in range(total_passes):
        pass_amplitudes, pass_uncertainties, pass_summary, pass_steps = _run_staged_component_pass(
            target_hist,
            template_hists,
            fitted_names,
            fit_min,
            fit_max,
            combined_window_map,
            stage_amplitude_window_map,
            stage_amplitude_mode_map,
            exclude_windows,
            stage_amplitudes,
            amplitude_prefix,
            pass_index,
        )
        stage_amplitudes.update(pass_amplitudes)
        stage_uncertainties.update(pass_uncertainties)
        staged_pass_history.append(pass_summary)
        step_overlays.extend(pass_steps)

    include_windows = _collect_unique_windows(
        combined_window_map,
        ordered_names=fitted_names,
    )
    validation_kwargs = {
        "oversub_sigma_tolerance": validation_options.get("oversub_sigma_tolerance", 2.0),
        "max_oversub_bin_count": validation_options.get("max_oversub_bin_count"),
        "max_oversub_bin_fraction": validation_options.get("max_oversub_bin_fraction"),
        "max_full_range_chi2_ndf": validation_options.get("max_full_range_chi2_ndf"),
        "cleanup_validation_mm_max": validation_options.get("cleanup_validation_mm_max"),
    }

    staged_template_hists = {name: template_hists[name] for name in fitted_names}
    fit_hist = _build_model_hist(
        target_hist,
        staged_template_hists,
        stage_amplitudes,
        "{}_fit_hist_{}".format(amplitude_prefix, context),
    )
    residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_{}".format(amplitude_prefix, context),
    )
    residual_hist.Add(fit_hist, -1.0)
    scaled_hist_map = _build_scaled_hist_map(
        staged_template_hists,
        stage_amplitudes,
        amplitude_prefix,
        context,
    )
    extra_scaled_hists = {
        template_name: scaled_hist_map.get(template_name)
        for template_name in extra_positive_templates
    }
    stage_validation = _evaluate_model_validation(
        target_hist,
        fit_hist,
        fit_min,
        fit_max,
        n_parameters=len(fitted_names),
        exclude_windows=exclude_windows,
        **validation_kwargs
    )
    prior_sigmas = _build_prior_sigma_map(
        fitted_names,
        stage_uncertainties,
        prior_scale_map=prior_scale_map,
    )
    pion_bg_fit_hist = _build_model_hist(
        target_hist,
        {
            component_name: staged_template_hists.get(component_name)
            for component_name in COMPONENT_NAMES
        },
        stage_amplitudes,
        "{}_pion_bg_fit_{}".format(amplitude_prefix, context),
    ) if amplitude_prefix == "A" else None

    fit_status = "success"
    fallback_used = False
    fallback_reason = ""
    if not stage_validation["accepted"]:
        fit_status = "fallback"
        fallback_used = True
        fallback_reason = (
            "ordered residual staged fit failed full-range validation: {}".format(
                "; ".join(stage_validation.get("rejection_reasons") or ["unknown"])
            )
        )

    extra_component_amplitudes = {
        template_name: float(stage_amplitudes.get(template_name, 0.0) or 0.0)
        for template_name in extra_positive_templates
    }
    diagnostics = {
        "success": (not fallback_used),
        "status_code": None,
        "message": "ordered residual staged fit",
        "chi2": stage_validation.get("chi2"),
        "ndf": stage_validation.get("ndf"),
        "chi2_ndf": stage_validation.get("chi2_ndf"),
        "fit_p_value": stage_validation.get("fit_p_value"),
        "n_fit_bins": stage_validation.get("n_fit_bins"),
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "include_windows": deepcopy(include_windows),
        "exclude_windows": deepcopy(exclude_windows),
        "fallback_used": bool(fallback_used),
        "fallback_reason": fallback_reason,
        "fit_strategy": "ordered_residual",
        "accepted_solution": "stage_only",
        "n_passes": total_passes,
        "fit_order": list(fitted_names),
        "anchor_windows": deepcopy(combined_window_map),
        "stage_amplitude_windows": deepcopy(stage_amplitude_window_map),
        "stage_amplitude_modes": deepcopy(stage_amplitude_mode_map),
        "staged_pass_history": deepcopy(staged_pass_history),
        "staged_component_amplitudes": {
            component_name: float(stage_amplitudes.get(component_name, 0.0) or 0.0)
            for component_name in fitted_names
        },
        "staged_component_uncertainties": {
            component_name: (
                None if stage_uncertainties.get(component_name) is None
                else float(stage_uncertainties.get(component_name))
            )
            for component_name in fitted_names
        },
        "accepted_component_uncertainties": {
            component_name: (
                None if stage_uncertainties.get(component_name) is None
                else float(stage_uncertainties.get(component_name))
            )
            for component_name in fitted_names
        },
        "prior_scales": deepcopy(prior_scale_map),
        "prior_sigmas": deepcopy(prior_sigmas),
        "joint_refinement": {
            "enabled": bool(joint_refinement_enabled),
            "requested": bool(joint_refinement_enabled),
            "status": "not_requested",
            "success": False,
            "message": "joint refinement not run inside staged fitter",
            "n_fit_bins": 0,
            "prior_sigmas": deepcopy(prior_sigmas),
            "coefficients": {},
            "uncertainties": {},
            "validation": {},
        },
        "coordinate_refinement": {
            "attempted": False,
            "success": False,
            "message": "coordinate refinement disabled by staged-seeded joint update",
            "n_fit_bins": 0,
            "cycles_run": 0,
            "history": [],
            "coefficients": {},
            "validation": {},
        },
        "stage_validation": deepcopy(stage_validation),
        "validation": deepcopy(stage_validation),
        "component_amplitudes": {
            component_name: float(stage_amplitudes.get(component_name, 0.0) or 0.0)
            for component_name in COMPONENT_NAMES
        },
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
    }

    result = {
        "{}_n".format(amplitude_prefix): float(stage_amplitudes.get("pi_n", 0.0) or 0.0),
        "{}_delta".format(amplitude_prefix): float(stage_amplitudes.get("pi_delta", 0.0) or 0.0),
        "{}_sidis".format(amplitude_prefix): float(stage_amplitudes.get("pi_sidis", 0.0) or 0.0),
        "fit_status": fit_status,
        "diagnostics": diagnostics,
        "fit_hist": fit_hist,
        "residual_hist": residual_hist,
        "template_hists": staged_template_hists,
        "pi_n_scaled_hist": scaled_hist_map.get("pi_n"),
        "pi_delta_scaled_hist": scaled_hist_map.get("pi_delta"),
        "pi_sidis_scaled_hist": scaled_hist_map.get("pi_sidis"),
        "extra_scaled_hists": extra_scaled_hists,
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
        "step_overlays": step_overlays,
    }
    if amplitude_prefix == "A":
        result["pion_bg_fit_hist"] = pion_bg_fit_hist
    return result


def _run_component_fit(
    target_hist,
    component_hists,
    amplitude_prefix,
    fit_min,
    fit_max,
    include_linear_background=False,
    extra_positive_templates=None,
    include_windows=None,
    exclude_windows=None,
    context="",
):
    if extra_positive_templates is None:
        extra_positive_templates = {}
    fallback_reason = _validate_component_shapes(component_hists, target_hist)
    if fallback_reason:
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            fallback_reason,
            template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
        )
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(
                target_hist,
                amplitude_prefix,
                context,
                fallback_reason,
                template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
            )

    fit_inputs = _build_fit_inputs(
        target_hist,
        component_hists,
        fit_min,
        fit_max,
        extra_template_hists=extra_positive_templates,
        include_windows=include_windows,
        exclude_windows=exclude_windows,
    )
    x_values = fit_inputs["x"]
    y_values = fit_inputs["y"]
    sigma_values = fit_inputs["sigma"]
    if len(x_values) == 0:
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "no valid fit bins",
            template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
        )

    column_names = list(COMPONENT_NAMES)
    design_columns = [fit_inputs["component_columns"][name] for name in COMPONENT_NAMES]
    lower_bounds = [0.0, 0.0, 0.0]
    upper_bounds = [np.inf, np.inf, np.inf]
    for template_name in extra_positive_templates:
        column_names.append(template_name)
        design_columns.append(fit_inputs["extra_template_columns"][template_name])
        lower_bounds.append(0.0)
        upper_bounds.append(np.inf)
    if include_linear_background:
        column_names.extend(("C0", "C1"))
        design_columns.extend((np.ones_like(x_values), x_values))
        lower_bounds.extend((-np.inf, -np.inf))
        upper_bounds.extend((np.inf, np.inf))

    design_matrix = np.column_stack(design_columns)
    weighted_design = design_matrix / sigma_values[:, None]
    weighted_target = y_values / sigma_values

    if len(x_values) <= len(column_names):
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "insufficient fit bins",
            template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
        )

    try:
        fit_result = lsq_linear(
            weighted_design,
            weighted_target,
            bounds=(np.asarray(lower_bounds, dtype=float), np.asarray(upper_bounds, dtype=float)),
            method="trf",
        )
    except Exception as exc:
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "fit solver exception: {}".format(exc),
            template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
        )

    if not getattr(fit_result, "success", False):
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "fit solver failed: {}".format(getattr(fit_result, "message", "unknown")),
            template_hists=_merge_template_hist_maps(component_hists, extra_positive_templates),
        )

    parameter_map = {
        name: float(value)
        for name, value in zip(column_names, np.asarray(fit_result.x, dtype=float))
    }
    for component_name in COMPONENT_NAMES:
        parameter_map[component_name] = max(parameter_map[component_name], 0.0)
    extra_component_amplitudes = {}
    for template_name in extra_positive_templates:
        extra_component_amplitudes[template_name] = max(parameter_map[template_name], 0.0)
        parameter_map[template_name] = extra_component_amplitudes[template_name]

    fit_values = design_matrix.dot(np.asarray(fit_result.x, dtype=float))
    residual_values = y_values - fit_values
    chi2_value = float(np.sum(np.square(residual_values / sigma_values)))
    ndf_value = int(len(x_values) - len(column_names))
    chi2_ndf_value = (chi2_value / ndf_value) if ndf_value > 0 else None
    p_value = float(chi2_dist.sf(chi2_value, ndf_value)) if ndf_value > 0 else None

    pi_n_scaled_hist = _clone_hist(
        component_hists["pi_n"],
        "{}_pi_n_scaled_{}".format(amplitude_prefix, context),
    )
    pi_delta_scaled_hist = _clone_hist(
        component_hists["pi_delta"],
        "{}_pi_delta_scaled_{}".format(amplitude_prefix, context),
    )
    pi_sidis_scaled_hist = _clone_hist(
        component_hists["pi_sidis"],
        "{}_pi_sidis_scaled_{}".format(amplitude_prefix, context),
    )
    pi_n_scaled_hist.Scale(parameter_map["pi_n"])
    pi_delta_scaled_hist.Scale(parameter_map["pi_delta"])
    pi_sidis_scaled_hist.Scale(parameter_map["pi_sidis"])
    extra_scaled_hists = {}
    for template_name, template_hist in extra_positive_templates.items():
        scaled_hist = _clone_hist(
            template_hist,
            "{}_{}_scaled_{}".format(amplitude_prefix, template_name, context),
        )
        scaled_hist.Scale(parameter_map[template_name])
        extra_scaled_hists[template_name] = scaled_hist

    fit_hist = _clone_hist(
        target_hist,
        "{}_fit_hist_{}".format(amplitude_prefix, context),
        reset=True,
    )
    pion_bg_fit_hist = _clone_hist(
        target_hist,
        "{}_pion_bg_fit_{}".format(amplitude_prefix, context),
        reset=True,
    )
    residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_{}".format(amplitude_prefix, context),
        reset=True,
    )

    fit_bin_values = []
    pion_bg_bin_values = []
    residual_bin_values = []
    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        pion_bg_value = (
            parameter_map["pi_n"] * component_hists["pi_n"].GetBinContent(bin_index)
            + parameter_map["pi_delta"] * component_hists["pi_delta"].GetBinContent(bin_index)
            + parameter_map["pi_sidis"] * component_hists["pi_sidis"].GetBinContent(bin_index)
        )
        extra_template_value = 0.0
        for template_name, template_hist in extra_positive_templates.items():
            extra_template_value += parameter_map[template_name] * template_hist.GetBinContent(bin_index)
        total_fit_value = pion_bg_value + extra_template_value
        if include_linear_background:
            total_fit_value += parameter_map["C0"] + parameter_map["C1"] * x_center
        fit_bin_values.append(total_fit_value)
        pion_bg_bin_values.append(pion_bg_value)
        residual_bin_values.append(float(target_hist.GetBinContent(bin_index)) - total_fit_value)

    _set_hist_values(fit_hist, fit_bin_values)
    _set_hist_values(pion_bg_fit_hist, pion_bg_bin_values)
    _set_hist_values(residual_hist, residual_bin_values)

    fit_status = "success"
    diagnostics = {
        "success": True,
        "status_code": getattr(fit_result, "status", None),
        "message": getattr(fit_result, "message", ""),
        "chi2": chi2_value,
        "ndf": ndf_value,
        "chi2_ndf": chi2_ndf_value,
        "fit_p_value": p_value,
        "n_fit_bins": int(len(x_values)),
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "include_windows": deepcopy(include_windows or []),
        "exclude_windows": deepcopy(exclude_windows or []),
        "fallback_used": False,
        "fallback_reason": "",
        "component_amplitudes": {
            component_name: float(parameter_map[component_name])
            for component_name in COMPONENT_NAMES
        },
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
    }
    if include_linear_background:
        diagnostics["background_terms"] = {
            "C0": float(parameter_map["C0"]),
            "C1": float(parameter_map["C1"]),
        }

    result = {
        "{}_n".format(amplitude_prefix): float(parameter_map["pi_n"]),
        "{}_delta".format(amplitude_prefix): float(parameter_map["pi_delta"]),
        "{}_sidis".format(amplitude_prefix): float(parameter_map["pi_sidis"]),
        "fit_status": fit_status,
        "diagnostics": diagnostics,
        "fit_hist": fit_hist,
        "residual_hist": residual_hist,
        "pi_n_scaled_hist": pi_n_scaled_hist,
        "pi_delta_scaled_hist": pi_delta_scaled_hist,
        "pi_sidis_scaled_hist": pi_sidis_scaled_hist,
        "extra_scaled_hists": extra_scaled_hists,
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
    }
    if amplitude_prefix == "A":
        result["pion_bg_fit_hist"] = pion_bg_fit_hist
    return result


def fit_pion_control_with_simc_shapes(
    h_pion_control,
    h_pi_n_shape,
    h_pi_delta_shape,
    h_pi_sidis_shape,
    h_kaon_sigma0_shape,
    inpDict,
    mm_offset_data=0.0,
    phi_setting=None,
    context="",
    _skip_residual_shift=False,
    alignment_windows=None,
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    fit_config = get_particle_subtraction_component_fit_window_config(
        "pion_control",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    ) or {}
    resolved_windows = resolve_particle_subtraction_component_fit_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    stage_amplitude_windows = resolve_particle_subtraction_component_stage_amplitude_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    stage_amplitude_modes = resolve_particle_subtraction_component_stage_amplitude_modes(
        "pion_control",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    if isinstance(alignment_windows, dict):
        for component_name in COMPONENT_NAMES:
            if component_name not in alignment_windows:
                continue
            resolved_windows[component_name] = _normalize_window_collection(
                alignment_windows.get(component_name)
            )
            if component_name in stage_amplitude_windows:
                stage_amplitude_windows[component_name] = _normalize_window_collection(
                    alignment_windows.get(component_name)
                )
    excluded_windows = resolve_particle_subtraction_component_fit_excluded_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    anchor_windows = {
        component_name: deepcopy(windows)
        for component_name, windows in resolved_windows.items()
    }
    prior_scale_map = resolve_particle_subtraction_component_prior_scales(
        "pion_control",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    fit_mode = resolve_particle_subtraction_component_fit_mode(
        "pion_control",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    residual_shift_settings = resolve_particle_subtraction_component_residual_shift_settings(
        "pion_control",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    postfit_scale_map = resolve_particle_subtraction_component_postfit_scales(
        "pion_control",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    postrefine_scale_map = resolve_particle_subtraction_component_postrefine_scales(
        "pion_control",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    extra_positive_templates = {}
    extra_anchor_windows = {}
    if anchor_windows.get(KAON_SIGMA0_TEMPLATE_NAME) and _hist_has_usable_support(h_kaon_sigma0_shape):
        extra_positive_templates[KAON_SIGMA0_TEMPLATE_NAME] = h_kaon_sigma0_shape
    validation_options = {
        "oversub_sigma_tolerance": fit_config.get("oversub_sigma_tolerance", 2.0),
        "max_oversub_bin_count": fit_config.get("max_oversub_bin_count"),
        "max_oversub_bin_fraction": fit_config.get("max_oversub_bin_fraction"),
        "max_full_range_chi2_ndf": fit_config.get("max_full_range_chi2_ndf"),
        "cleanup_validation_mm_max": _resolve_cleanup_validation_max(
            "pion_control",
            mm_offset_data=mm_offset_data,
            inp_dict=inpDict,
            phi_setting=phi_setting,
        ),
    }
    amplitude_floor = float(fit_config.get("joint_refinement_amplitude_floor", 1e-3) or 1e-3)
    template_corr_warn = float(fit_config.get("template_corr_warn", 0.95) or 0.95)
    if not _skip_residual_shift:
        shift_selection = _run_component_residual_shift_selection(
            "pion_control",
            h_pion_control,
            {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
            },
            {
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            residual_shift_settings,
            fit_callback=lambda shifted_components, shifted_extra_templates: fit_pion_control_with_simc_shapes(
                h_pion_control,
                shifted_components.get("pi_n"),
                shifted_components.get("pi_delta"),
                shifted_components.get("pi_sidis"),
                shifted_extra_templates.get(KAON_SIGMA0_TEMPLATE_NAME),
                inpDict,
                mm_offset_data=mm_offset_data,
                phi_setting=phi_setting,
                context=context,
                _skip_residual_shift=True,
            ),
            fit_min=fit_min,
            cleanup_validation_mm_max=validation_options.get("cleanup_validation_mm_max"),
            exclude_windows=excluded_windows,
            context="{}_pion_control".format(context or "scope"),
        )
        if shift_selection is not None:
            return _annotate_fit_result_with_residual_shift_payload(
                shift_selection.get("fit_result"),
                "pion_control",
                residual_shift_settings,
                shift_selection,
            )
    result = _fit_staged_anchor_templates(
        h_pion_control,
        {
            "pi_n": h_pi_n_shape,
            "pi_delta": h_pi_delta_shape,
            "pi_sidis": h_pi_sidis_shape,
        },
        "B",
        fit_min,
        fit_max,
        anchor_windows=anchor_windows,
        fit_order=fit_config.get("fit_order") or ("pi_n", "pi_delta", "pi_sidis"),
        stage_amplitude_windows=stage_amplitude_windows,
        stage_amplitude_modes=stage_amplitude_modes,
        exclude_windows=excluded_windows,
        extra_positive_templates=extra_positive_templates,
        extra_anchor_windows=extra_anchor_windows,
        n_passes=fit_config.get("staged_fit_passes", 1),
        prior_scale_map=prior_scale_map,
        joint_refinement_enabled=bool(fit_config.get("joint_refinement_enabled", True)),
        max_fit_cycles=fit_config.get("particle_subtraction_max_fit_cycles", 50),
        fit_tolerance=fit_config.get("particle_subtraction_fit_tolerance", 1e-5),
        validation_options=validation_options,
        context="{}_pion_control".format(context or "scope"),
    )
    result = _apply_component_postfit_scales(
        result,
        h_pion_control,
        "B",
        postfit_scale_map,
        fit_min,
        fit_max,
        exclude_windows=excluded_windows,
        validation_options=validation_options,
        context="{}_pion_control".format(context or "scope"),
    )
    result = _apply_joint_component_refinement(
        result,
        h_pion_control,
        "B",
        fit_mode,
        fit_min,
        fit_max,
        exclude_windows=excluded_windows,
        prior_scale_map=prior_scale_map,
        postrefine_scale_map=postrefine_scale_map,
        validation_options=validation_options,
        context="{}_pion_control".format(context or "scope"),
        template_corr_warn=template_corr_warn,
        amplitude_floor=amplitude_floor,
    )
    sigma0_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    sigma0_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    refined_scaled_hist_map = result.get("refined_scaled_hist_map") or {}
    return_payload = {
        "B_n": result["B_n"],
        "B_delta": result["B_delta"],
        "B_sidis": result["B_sidis"],
        "B_sigma0": None if sigma0_amplitude is None else float(sigma0_amplitude),
        "fit_mode": fit_mode,
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "template_hists": result.get("template_hists") or {},
        "fit_hist": result["fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
        "k_sigma0_scaled_hist": sigma0_scaled_hist,
        "refined_fit_hist": result.get("refined_fit_hist") or result["fit_hist"],
        "refined_residual_hist": result.get("refined_residual_hist") or result["residual_hist"],
        "refined_pi_n_scaled_hist": refined_scaled_hist_map.get("pi_n") or result["pi_n_scaled_hist"],
        "refined_pi_delta_scaled_hist": refined_scaled_hist_map.get("pi_delta") or result["pi_delta_scaled_hist"],
        "refined_pi_sidis_scaled_hist": refined_scaled_hist_map.get("pi_sidis") or result["pi_sidis_scaled_hist"],
        "refined_k_sigma0_scaled_hist": refined_scaled_hist_map.get(KAON_SIGMA0_TEMPLATE_NAME) or sigma0_scaled_hist,
        "step_overlays": result.get("step_overlays") or [],
        "resolved_config_summary": {
            "fit_target": "pion_control",
            "particle_subtraction_setting_key": fit_config.get("particle_subtraction_setting_key"),
            "particle_subtraction_phi_setting": fit_config.get("particle_subtraction_phi_setting"),
            "particle_subtraction_override_layers": deepcopy(
                fit_config.get("particle_subtraction_override_layers") or []
            ),
            "particle_subtraction_override_applied": bool(
                fit_config.get("particle_subtraction_override_applied", False)
            ),
            "fit_order": deepcopy(fit_config.get("fit_order") or []),
            "anchor_windows": deepcopy(anchor_windows),
            "excluded_windows": deepcopy(excluded_windows),
            "stage_amplitude_windows": deepcopy(stage_amplitude_windows),
            "stage_amplitude_modes": deepcopy(stage_amplitude_modes),
            "prior_scales": deepcopy(prior_scale_map),
            "postfit_component_scales": deepcopy(postfit_scale_map),
            "postrefine_component_scales": deepcopy(postrefine_scale_map),
            "fit_mode": fit_mode,
            "joint_refinement_amplitude_floor": amplitude_floor,
            "template_corr_warn": template_corr_warn,
            "cleanup_validation_mm_max": validation_options.get("cleanup_validation_mm_max"),
        },
    }
    if _skip_residual_shift:
        return return_payload

    return _annotate_fit_result_with_residual_shift_payload(
        return_payload,
        "pion_control",
        residual_shift_settings,
        {
            "summary": {
                "enabled": False,
                "mode": residual_shift_settings.get("mode"),
                "units": residual_shift_settings.get("units"),
                "selection_metric": residual_shift_settings.get("selection_metric"),
                "requested_components": list(residual_shift_settings.get("components") or []),
                "active_components": [],
                "configured_shift_values": deepcopy(residual_shift_settings.get("values") or {}),
                "shift_bounds": deepcopy(residual_shift_settings.get("bounds") or {}),
                "shift_grid": deepcopy(residual_shift_settings.get("scan_grid") or {}),
                "selected_shift_point": {},
                "selected_shift_reason": "residual shifts disabled",
                "candidate_count": 1,
                "candidate_summaries": [],
                "per_component": {},
                "cleanup_metrics": {},
                "warnings": [],
            },
            "original_template_hists": {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            "selected_template_hists": {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
            },
            "selected_extra_template_hists": {
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            "selected_component_diagnostics": {},
        },
    )


def fit_kaon_nosub_with_simc_pion_shapes(
    h_kaon_nosub,
    h_pi_n_shape,
    h_pi_delta_shape,
    h_pi_sidis_shape,
    h_kaon_signal_shape,
    h_kaon_sigma0_shape,
    inpDict,
    mm_offset_data=0.0,
    phi_setting=None,
    context="",
    _skip_residual_shift=False,
    alignment_windows=None,
    _alignment_scoring_only=False,
    kaon_sigma0_source_diagnostics=None,
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    fit_config = get_particle_subtraction_component_fit_window_config(
        "kaon_nosub",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    ) or {}
    resolved_windows = resolve_particle_subtraction_component_fit_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    stage_amplitude_windows = resolve_particle_subtraction_component_stage_amplitude_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    stage_amplitude_modes = resolve_particle_subtraction_component_stage_amplitude_modes(
        "kaon_nosub",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    if isinstance(alignment_windows, dict):
        for component_name in COMPONENT_NAMES:
            if component_name not in alignment_windows:
                continue
            resolved_windows[component_name] = _normalize_window_collection(
                alignment_windows.get(component_name)
            )
            if component_name in stage_amplitude_windows:
                stage_amplitude_windows[component_name] = _normalize_window_collection(
                    alignment_windows.get(component_name)
                )
    excluded_windows = resolve_particle_subtraction_component_fit_excluded_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    anchor_windows = {
        component_name: deepcopy(windows)
        for component_name, windows in resolved_windows.items()
    }
    mm_min = float(inpDict.get("mm_min", fit_min))
    mm_max = float(inpDict.get("mm_max", fit_max))
    extra_positive_templates = {}
    extra_anchor_windows = {}
    prior_scale_map = resolve_particle_subtraction_component_prior_scales(
        "kaon_nosub",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME, KAON_SIGNAL_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    fit_mode = resolve_particle_subtraction_component_fit_mode(
        "kaon_nosub",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    protected_fit_config = _resolve_pi_delta_signal_protected_fit_config(fit_config)
    residual_shift_settings = resolve_particle_subtraction_component_residual_shift_settings(
        "kaon_nosub",
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    postfit_scale_map = resolve_particle_subtraction_component_postfit_scales(
        "kaon_nosub",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    postrefine_scale_map = resolve_particle_subtraction_component_postrefine_scales(
        "kaon_nosub",
        component_names=("pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME, KAON_SIGNAL_TEMPLATE_NAME),
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    validation_options = {
        "oversub_sigma_tolerance": fit_config.get("oversub_sigma_tolerance", 2.0),
        "max_oversub_bin_count": fit_config.get("max_oversub_bin_count"),
        "max_oversub_bin_fraction": fit_config.get("max_oversub_bin_fraction"),
        "max_full_range_chi2_ndf": fit_config.get("max_full_range_chi2_ndf"),
        "cleanup_validation_mm_max": _resolve_cleanup_validation_max(
            "kaon_nosub",
            mm_offset_data=mm_offset_data,
            inp_dict=inpDict,
            phi_setting=phi_setting,
        ),
    }
    amplitude_floor = float(fit_config.get("joint_refinement_amplitude_floor", 1e-3) or 1e-3)
    template_corr_warn = float(fit_config.get("template_corr_warn", 0.95) or 0.95)
    include_kaon_signal_template = bool(fit_config.get("include_kaon_signal_template", False))
    if include_kaon_signal_template and _hist_has_usable_support(h_kaon_signal_shape):
        extra_positive_templates[KAON_SIGNAL_TEMPLATE_NAME] = h_kaon_signal_shape
        if mm_max > mm_min:
            tail_extension = max(float(fit_config.get("kaon_signal_tail_extension", 0.0) or 0.0), 0.0)
            signal_window_max = min(fit_max, mm_max + tail_extension)
            extra_anchor_windows[KAON_SIGNAL_TEMPLATE_NAME] = [(mm_min, signal_window_max)]
    if anchor_windows.get(KAON_SIGMA0_TEMPLATE_NAME) and _hist_has_usable_support(h_kaon_sigma0_shape):
        extra_positive_templates[KAON_SIGMA0_TEMPLATE_NAME] = h_kaon_sigma0_shape
    if not _skip_residual_shift:
        shift_selection = _run_component_residual_shift_selection(
            "kaon_nosub",
            h_kaon_nosub,
            {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
            },
            {
                KAON_SIGNAL_TEMPLATE_NAME: h_kaon_signal_shape,
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            residual_shift_settings,
            fit_callback=lambda shifted_components, shifted_extra_templates: fit_kaon_nosub_with_simc_pion_shapes(
                h_kaon_nosub,
                shifted_components.get("pi_n"),
                shifted_components.get("pi_delta"),
                shifted_components.get("pi_sidis"),
                shifted_extra_templates.get(KAON_SIGNAL_TEMPLATE_NAME),
                shifted_extra_templates.get(KAON_SIGMA0_TEMPLATE_NAME),
                inpDict,
                mm_offset_data=mm_offset_data,
                phi_setting=phi_setting,
                context=context,
                _skip_residual_shift=True,
                _alignment_scoring_only=True,
                kaon_sigma0_source_diagnostics=kaon_sigma0_source_diagnostics,
            ),
            fit_min=fit_min,
            cleanup_validation_mm_max=validation_options.get("cleanup_validation_mm_max"),
            exclude_windows=excluded_windows,
            context="{}_kaon_nosub".format(context or "scope"),
        )
        if shift_selection is not None:
            selected_result = shift_selection.get("fit_result")
            if not bool(_alignment_scoring_only):
                selected_components = shift_selection.get("selected_template_hists") or {}
                selected_extra_templates = shift_selection.get("selected_extra_template_hists") or {}
                selected_result = _apply_signal_protected_pi_delta_fit(
                    selected_result,
                    h_kaon_nosub,
                    selected_components.get("pi_n", h_pi_n_shape),
                    selected_components.get("pi_delta", h_pi_delta_shape),
                    selected_components.get("pi_sidis", h_pi_sidis_shape),
                    selected_extra_templates.get(KAON_SIGNAL_TEMPLATE_NAME, h_kaon_signal_shape),
                    selected_extra_templates.get(KAON_SIGMA0_TEMPLATE_NAME, h_kaon_sigma0_shape),
                    fit_config,
                    fit_min,
                    fit_max,
                    excluded_windows,
                    inpDict,
                    phi_setting,
                    mm_offset_data,
                    context="{}_kaon_nosub".format(context or "scope"),
                    kaon_sigma0_source_diagnostics=kaon_sigma0_source_diagnostics,
                )
            return _annotate_fit_result_with_residual_shift_payload(
                selected_result,
                "kaon_nosub",
                residual_shift_settings,
                shift_selection,
            )
    result = _fit_staged_anchor_templates(
        h_kaon_nosub,
        {
            "pi_n": h_pi_n_shape,
            "pi_delta": h_pi_delta_shape,
            "pi_sidis": h_pi_sidis_shape,
        },
        "A",
        fit_min,
        fit_max,
        anchor_windows=anchor_windows,
        exclude_windows=excluded_windows,
        fit_order=fit_config.get("fit_order") or (
            "pi_n",
            "pi_sidis",
            "pi_delta",
        ),
        stage_amplitude_windows=stage_amplitude_windows,
        stage_amplitude_modes=stage_amplitude_modes,
        extra_positive_templates=extra_positive_templates,
        extra_anchor_windows=extra_anchor_windows,
        n_passes=fit_config.get("staged_fit_passes", 1),
        prior_scale_map=prior_scale_map,
        joint_refinement_enabled=bool(fit_config.get("joint_refinement_enabled", True)),
        max_fit_cycles=fit_config.get("particle_subtraction_max_fit_cycles", 50),
        fit_tolerance=fit_config.get("particle_subtraction_fit_tolerance", 1e-5),
        validation_options=validation_options,
        context="{}_kaon_nosub".format(context or "scope"),
    )
    result = _apply_component_postfit_scales(
        result,
        h_kaon_nosub,
        "A",
        postfit_scale_map,
        fit_min,
        fit_max,
        exclude_windows=excluded_windows,
        validation_options=validation_options,
        context="{}_kaon_nosub".format(context or "scope"),
    )
    result = _apply_joint_component_refinement(
        result,
        h_kaon_nosub,
        "A",
        fit_mode,
        fit_min,
        fit_max,
        exclude_windows=excluded_windows,
        prior_scale_map=prior_scale_map,
        postrefine_scale_map=postrefine_scale_map,
        validation_options=validation_options,
        context="{}_kaon_nosub".format(context or "scope"),
        template_corr_warn=template_corr_warn,
        amplitude_floor=amplitude_floor,
    )
    signal_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    signal_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    sigma0_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    sigma0_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    refined_scaled_hist_map = result.get("refined_scaled_hist_map") or {}
    signal_reference_hist, signal_reference_scale = _build_scaled_reference_hist(
        h_kaon_nosub,
        h_kaon_signal_shape,
        mm_min,
        mm_max,
        "A_k_lambda_reference_{}_kaon_nosub".format(context or "scope"),
    )
    k_lambda_simc_reference_hist = _clone_hist(
        h_kaon_signal_shape,
        "H_k_lambda_simc_reference_{}_kaon_nosub".format(context or "scope"),
    )
    k_sigma0_simc_reference_hist = _clone_hist(
        h_kaon_sigma0_shape,
        "H_k_sigma0_simc_reference_{}_kaon_nosub".format(context or "scope"),
    )
    k_lambda_simc_input_loaded = bool(h_kaon_signal_shape is not None)
    k_lambda_simc_reference_available = bool(k_lambda_simc_reference_hist is not None)
    return_payload = {
        "A_n": result["A_n"],
        "A_delta": result["A_delta"],
        "A_sidis": result["A_sidis"],
        "S_lambda": None if signal_amplitude is None else float(signal_amplitude),
        "S_sigma0": None if sigma0_amplitude is None else float(sigma0_amplitude),
        "fit_mode": fit_mode,
        "S_lambda_reference_scale": (
            None if signal_reference_scale is None else float(signal_reference_scale)
        ),
        "k_lambda_reference_scale": (
            None if signal_reference_scale is None else float(signal_reference_scale)
        ),
        "k_lambda_fit_amplitude": None if signal_amplitude is None else float(signal_amplitude),
        "k_lambda_simc_input_loaded": bool(k_lambda_simc_input_loaded),
        "k_lambda_simc_reference_available": bool(k_lambda_simc_reference_available),
        "k_lambda_simc_reference_source": (
            "immutable_aligned_k_lambda_simc"
            if k_lambda_simc_reference_available
            else ("K-Lambda SIMC reference unavailable" if not k_lambda_simc_input_loaded else "reference_build_failed")
        ),
        "k_lambda_simc_reference_integral": (
            _hist_integral(k_lambda_simc_reference_hist)
            if k_lambda_simc_reference_hist is not None
            else None
        ),
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "template_hists": result.get("template_hists") or {},
        "fit_hist": result["fit_hist"],
        "pion_bg_fit_hist": result["pion_bg_fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
        "k_lambda_scaled_hist": signal_scaled_hist,
        "k_sigma0_scaled_hist": sigma0_scaled_hist,
        "refined_fit_hist": result.get("refined_fit_hist") or result["fit_hist"],
        "refined_fit_hist_pre_postrefine": result.get("refined_fit_hist_pre_postrefine"),
        "refined_pion_bg_fit_hist": result.get("refined_pion_bg_fit_hist") or result["pion_bg_fit_hist"],
        "refined_pion_bg_fit_hist_pre_postrefine": result.get("refined_pion_bg_fit_hist_pre_postrefine"),
        "refined_residual_hist": result.get("refined_residual_hist") or result["residual_hist"],
        "refined_pi_n_scaled_hist": refined_scaled_hist_map.get("pi_n") or result["pi_n_scaled_hist"],
        "refined_pi_delta_scaled_hist": refined_scaled_hist_map.get("pi_delta") or result["pi_delta_scaled_hist"],
        "refined_pi_sidis_scaled_hist": refined_scaled_hist_map.get("pi_sidis") or result["pi_sidis_scaled_hist"],
        "refined_k_lambda_scaled_hist": refined_scaled_hist_map.get(KAON_SIGNAL_TEMPLATE_NAME) or signal_scaled_hist,
        "refined_k_sigma0_scaled_hist": refined_scaled_hist_map.get(KAON_SIGMA0_TEMPLATE_NAME) or sigma0_scaled_hist,
        "refined_k_sigma0_scaled_hist_pre_postrefine": (
            (result.get("refined_scaled_hist_map_pre_postrefine") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
            or sigma0_scaled_hist
        ),
        "H_k_lambda_simc_reference": k_lambda_simc_reference_hist,
        "H_k_sigma0_simc_reference": k_sigma0_simc_reference_hist,
        "k_lambda_reference_hist": signal_reference_hist,
        "step_overlays": result.get("step_overlays") or [],
        "sigma0_requested": bool(anchor_windows.get(KAON_SIGMA0_TEMPLATE_NAME)),
        "sigma0_active": bool(
            sigma0_scaled_hist is not None
            or (refined_scaled_hist_map.get(KAON_SIGMA0_TEMPLATE_NAME) is not None)
        ),
        "resolved_config_summary": {
            "fit_target": "kaon_nosub",
            "particle_subtraction_setting_key": fit_config.get("particle_subtraction_setting_key"),
            "particle_subtraction_phi_setting": fit_config.get("particle_subtraction_phi_setting"),
            "particle_subtraction_override_layers": deepcopy(
                fit_config.get("particle_subtraction_override_layers") or []
            ),
            "particle_subtraction_override_applied": bool(
                fit_config.get("particle_subtraction_override_applied", False)
            ),
            "fit_order": deepcopy(fit_config.get("fit_order") or []),
            "anchor_windows": deepcopy(anchor_windows),
            "excluded_windows": deepcopy(excluded_windows),
            "stage_amplitude_windows": deepcopy(stage_amplitude_windows),
            "stage_amplitude_modes": deepcopy(stage_amplitude_modes),
            "prior_scales": deepcopy(prior_scale_map),
            "postfit_component_scales": deepcopy(postfit_scale_map),
            "postrefine_component_scales": deepcopy(postrefine_scale_map),
            "fit_mode": fit_mode,
            "sigma0_requested": bool(anchor_windows.get(KAON_SIGMA0_TEMPLATE_NAME)),
            "sigma0_active": bool(
                sigma0_scaled_hist is not None
                or (refined_scaled_hist_map.get(KAON_SIGMA0_TEMPLATE_NAME) is not None)
            ),
            "joint_refinement_amplitude_floor": amplitude_floor,
            "template_corr_warn": template_corr_warn,
            "cleanup_validation_mm_max": validation_options.get("cleanup_validation_mm_max"),
            "pi_delta_signal_protected_fit": deepcopy(protected_fit_config),
        },
    }
    if not bool(_alignment_scoring_only):
        return_payload = _apply_signal_protected_pi_delta_fit(
            return_payload,
            h_kaon_nosub,
            h_pi_n_shape,
            h_pi_delta_shape,
            h_pi_sidis_shape,
            h_kaon_signal_shape,
            h_kaon_sigma0_shape,
            fit_config,
            fit_min,
            fit_max,
            excluded_windows,
            inpDict,
            phi_setting,
            mm_offset_data,
            context="{}_kaon_nosub".format(context or "scope"),
            kaon_sigma0_source_diagnostics=kaon_sigma0_source_diagnostics,
        )
    if excluded_windows:
        for key in (
            "fit_hist",
            "full_fit_residual_hist",
            "pion_bg_fit_hist",
            "residual_hist",
            "pi_n_scaled_hist",
            "pi_delta_scaled_hist",
            "pi_sidis_scaled_hist",
            "k_lambda_scaled_hist",
            "k_sigma0_scaled_hist",
            "refined_fit_hist",
            "refined_full_fit_residual_hist",
            "refined_pion_bg_fit_hist",
            "refined_pion_bg_fit_hist_pre_postrefine",
            "refined_residual_hist",
            "refined_pi_n_scaled_hist",
            "refined_pi_delta_scaled_hist",
            "refined_pi_sidis_scaled_hist",
            "refined_k_lambda_scaled_hist",
            "refined_k_sigma0_scaled_hist",
            "refined_k_sigma0_scaled_hist_pre_postrefine",
            "refined_fit_hist_pre_postrefine",
        ):
            hist = return_payload.get(key)
            if hist is None:
                continue
            return_payload[key] = _masked_hist_clone(
                hist,
                "{}_masked".format(hist.GetName()),
                excluded_windows,
            )
    if _skip_residual_shift:
        return return_payload

    return _annotate_fit_result_with_residual_shift_payload(
        return_payload,
        "kaon_nosub",
        residual_shift_settings,
        {
            "summary": {
                "enabled": False,
                "mode": residual_shift_settings.get("mode"),
                "units": residual_shift_settings.get("units"),
                "selection_metric": residual_shift_settings.get("selection_metric"),
                "requested_components": list(residual_shift_settings.get("components") or []),
                "active_components": [],
                "configured_shift_values": deepcopy(residual_shift_settings.get("values") or {}),
                "shift_bounds": deepcopy(residual_shift_settings.get("bounds") or {}),
                "shift_grid": deepcopy(residual_shift_settings.get("scan_grid") or {}),
                "selected_shift_point": {},
                "selected_shift_reason": "residual shifts disabled",
                "candidate_count": 1,
                "candidate_summaries": [],
                "per_component": {},
                "cleanup_metrics": {},
                "warnings": [],
            },
            "original_template_hists": {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
                KAON_SIGNAL_TEMPLATE_NAME: h_kaon_signal_shape,
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            "selected_template_hists": {
                "pi_n": h_pi_n_shape,
                "pi_delta": h_pi_delta_shape,
                "pi_sidis": h_pi_sidis_shape,
            },
            "selected_extra_template_hists": {
                KAON_SIGNAL_TEMPLATE_NAME: h_kaon_signal_shape,
                KAON_SIGMA0_TEMPLATE_NAME: h_kaon_sigma0_shape,
            },
            "selected_component_diagnostics": {},
        },
    )


def _sum_hist_list_to_unit_area(hist_list, hist_name):
    template_hist = next((hist for hist in hist_list if hist is not None), None)
    if template_hist is None:
        return None

    summed_hist = _clone_hist(template_hist, hist_name, reset=True)
    for hist in hist_list:
        if hist is not None:
            summed_hist.Add(hist)
    normalize_hist_to_unit_area(summed_hist, quiet=True, context=hist_name)
    return summed_hist


def _hist_has_usable_support(hist, min_integral=1e-12):
    if hist is None:
        return False
    try:
        integral = float(hist.Integral())
    except Exception:
        return False
    return math.isfinite(integral) and integral > float(min_integral)


def _require_kaon_lambda_simc_shape(kaon_signal_shape, analysis_scope, phi_setting, context):
    if _hist_has_usable_support(kaon_signal_shape):
        return
    raise RuntimeError(
        "K-Lambda SIMC comparison is mandatory for kaon component subtraction; "
        "no usable K-Lambda template is available "
        "(phi_setting={}, analysis_scope={}, context={})".format(
            phi_setting or "unknown",
            analysis_scope or "unknown",
            context or "unknown",
        )
    )


def clone_kaon_lambda_comparison_payload(component_fit_result, context=""):
    """Return scalar K-Lambda comparison metadata for legacy callers.

    Histogram comparison state is fit-owned and is intentionally not copied
    into application payloads.  Renderers receive the fit result explicitly.
    """
    result = component_fit_result if isinstance(component_fit_result, dict) else {}
    payload = {}
    for scalar_key in (
        "S_lambda_reference_scale",
        "k_lambda_reference_scale",
        "k_lambda_fit_amplitude",
        "k_lambda_simc_input_loaded",
        "k_lambda_simc_reference_available",
        "k_lambda_simc_reference_source",
        "k_lambda_simc_reference_integral",
    ):
        payload[scalar_key] = deepcopy(result.get(scalar_key))
    return payload


def _resolve_component_scope_hist(component_entry, component_name, analysis_scope, t_bin_index=None, phi_bin_index=None):
    if not isinstance(component_entry, dict):
        return None

    setting_shape_full = component_entry.get("setting_shape_full")
    if t_bin_index is None and phi_bin_index is None:
        return setting_shape_full

    binned_shapes = component_entry.get("binned_shapes") or {}
    t_key = "t_bin{}".format(int(t_bin_index) + 1) if t_bin_index is not None else None
    if phi_bin_index is not None and t_key is not None:
        direct_hist = (binned_shapes.get(t_key) or {}).get("phi_bin{}".format(int(phi_bin_index) + 1))
        if _hist_has_usable_support(direct_hist):
            return direct_hist

        # Sparse or empty (t,phi) component bins are expected in fine slicing.
        # Fall back to the t-aggregated template before using the setting-wide shape.
        phi_shape_map = binned_shapes.get(t_key) or {}
        aggregated_hist = _sum_hist_list_to_unit_area(
            list(phi_shape_map.values()),
            "{}_{}_fallback_tbin{}".format(component_name, analysis_scope or "scope", int(t_bin_index) + 1),
        )
        if _hist_has_usable_support(aggregated_hist):
            return aggregated_hist
        return setting_shape_full

    if t_key is not None:
        phi_shape_map = binned_shapes.get(t_key) or {}
        aggregated_hist = _sum_hist_list_to_unit_area(
            list(phi_shape_map.values()),
            "{}_{}_aggregated".format(component_name, analysis_scope or "scope"),
        )
        if _hist_has_usable_support(aggregated_hist):
            return aggregated_hist
        return setting_shape_full

    return setting_shape_full


def _resolve_single_scope_hist(shape_payload, analysis_scope, t_bin_index=None, phi_bin_index=None):
    if not isinstance(shape_payload, dict):
        return None

    setting_shape_full = shape_payload.get("setting_shape_full")
    if t_bin_index is None and phi_bin_index is None:
        return setting_shape_full

    binned_shapes = shape_payload.get("binned_shapes") or {}
    t_key = "t_bin{}".format(int(t_bin_index) + 1) if t_bin_index is not None else None
    if phi_bin_index is not None and t_key is not None:
        direct_hist = (binned_shapes.get(t_key) or {}).get("phi_bin{}".format(int(phi_bin_index) + 1))
        if _hist_has_usable_support(direct_hist):
            return direct_hist

        phi_shape_map = binned_shapes.get(t_key) or {}
        aggregated_hist = _sum_hist_list_to_unit_area(
            list(phi_shape_map.values()),
            "single_shape_{}_fallback_tbin{}".format(analysis_scope or "scope", int(t_bin_index) + 1),
        )
        if _hist_has_usable_support(aggregated_hist):
            return aggregated_hist
        return setting_shape_full

    if t_key is not None:
        aggregated_hist = _sum_hist_list_to_unit_area(
            list((binned_shapes.get(t_key) or {}).values()),
            "single_shape_{}_aggregated".format(analysis_scope or "scope"),
        )
        if _hist_has_usable_support(aggregated_hist):
            return aggregated_hist
        return setting_shape_full

    return setting_shape_full


def _clone_shift_payload_hist(shift_payload, payload_keys, component_name, output_name):
    if not isinstance(payload_keys, (list, tuple)):
        payload_keys = [payload_keys]
    for payload_key in payload_keys:
        hist_map = ((shift_payload or {}).get(payload_key) or {})
        hist = hist_map.get(component_name)
        if hist is not None:
            return _clone_hist(hist, output_name)
    return None


def resolve_scope_component_shapes(
    component_payload,
    analysis_scope="setting-wide",
    t_bin_index=None,
    phi_bin_index=None,
):
    payload_components = (component_payload or {}).get("components") or {}
    resolved = {}
    for component_name in COMPONENT_NAMES:
        component_entry = payload_components.get(component_name) or {}
        resolved[component_name] = _resolve_component_scope_hist(
            component_entry,
            component_name,
            analysis_scope,
            t_bin_index=t_bin_index,
            phi_bin_index=phi_bin_index,
        )
    return resolved


def resolve_scope_single_shape(
    shape_payload,
    analysis_scope="setting-wide",
    t_bin_index=None,
    phi_bin_index=None,
):
    return _resolve_single_scope_hist(
        shape_payload,
        analysis_scope or "scope",
        t_bin_index=t_bin_index,
        phi_bin_index=phi_bin_index,
    )


ALIGNMENT_SCHEMA_VERSION = int(PION_COMPONENT_DYNAMIC_ALIGNMENT_SCHEMA_VERSION)


def _canonical_json(value):
    return json.dumps(
        _json_ready_particle_subtraction_value(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _alignment_hash(value):
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _hist_axis_specification(hist):
    if hist is None:
        return None
    axis = hist.GetXaxis()
    return {
        "nbins": int(hist.GetNbinsX()),
        "xmin": float(axis.GetXmin()),
        "xmax": float(axis.GetXmax()),
        "bin_edges": [float(axis.GetBinLowEdge(index)) for index in range(1, int(hist.GetNbinsX()) + 2)],
    }


def _hist_content_checksum(hist):
    if hist is None:
        return None
    return _alignment_hash({
        "axis": _hist_axis_specification(hist),
        "content": [float(hist.GetBinContent(index)) for index in range(0, int(hist.GetNbinsX()) + 2)],
        "error": [float(hist.GetBinError(index)) for index in range(0, int(hist.GetNbinsX()) + 2)],
    })


def _immutable_source_metadata(hist, component_name):
    hist_name = str(hist.GetName()) if hist is not None and hasattr(hist, "GetName") else None
    return {
        "immutable_source_hist_name": hist_name,
        "immutable_source_hist_identifier": {
            "component": str(component_name),
            "hist_name": hist_name,
            "axis": _hist_axis_specification(hist),
        },
        "immutable_source_hist_checksum": _hist_content_checksum(hist),
        "candidate_built_from_immutable_source": True,
        "applied_in_single_shift_operation": True,
        "shift_operation_count": 1,
    }


def _normalize_window_collection(windows):
    # A protected-fit window is commonly configured as one numeric pair,
    # while the component-window callers supply an iterable of such pairs.
    # Normalize the former into the latter before iterating so a valid
    # ``(min, max)`` configuration is not mistaken for two malformed windows.
    if (
        isinstance(windows, (list, tuple))
        and len(windows) == 2
        and all(_is_finite_number(value) for value in windows)
    ):
        windows = (windows,)
    normalized = []
    for window in windows or []:
        if not isinstance(window, (list, tuple)) or len(window) != 2:
            continue
        low_edge, high_edge = float(window[0]), float(window[1])
        if not (_is_finite_number(low_edge) and _is_finite_number(high_edge)):
            continue
        normalized.append((min(low_edge, high_edge), max(low_edge, high_edge)))
    return normalized


def transform_component_windows(canonical_windows, total_shift_gev, expansion_gev, mm_limits=None):
    """Translate pion-component windows and expand both edges symmetrically."""
    canonical = _normalize_window_collection(canonical_windows)
    shift = float(total_shift_gev) if _is_finite_number(total_shift_gev) else 0.0
    expansion = max(float(expansion_gev), 0.0) if _is_finite_number(expansion_gev) else 0.0
    limits = _normalize_window_collection([mm_limits]) if mm_limits is not None else []
    effective, clipping = [], []
    for low_edge, high_edge in canonical:
        requested = (low_edge + shift - expansion, high_edge + shift + expansion)
        resolved = requested
        if limits:
            resolved = (max(requested[0], limits[0][0]), min(requested[1], limits[0][1]))
        if resolved[1] >= resolved[0]:
            effective.append((float(resolved[0]), float(resolved[1])))
        clipping.append({"requested": requested, "effective": resolved, "clipped": resolved != requested})
    return {
        "canonical_windows": canonical,
        "effective_windows": effective,
        "total_shift_gev": shift,
        "window_expansion_gev": expansion,
        "window_clipping": clipping,
        "window_clipped": any(item["clipped"] for item in clipping),
    }


def _alignment_scan_grid(minimum, maximum, step, include_value=None):
    minimum, maximum, step = float(minimum), float(maximum), abs(float(step))
    if maximum < minimum or step <= 0.0:
        values = [minimum]
    else:
        count = int(math.floor((maximum - minimum) / step + 1e-9))
        values = [minimum + (index * step) for index in range(count + 1)]
        if values[-1] < maximum - 1e-10:
            values.append(maximum)
    if include_value is not None and _is_finite_number(include_value):
        values.append(float(include_value))
    return sorted({round(float(value), 12) for value in values})


def _alignment_windows_to_bin_indices(hist, windows):
    return [
        index for index in range(1, int(hist.GetNbinsX()) + 1)
        if any(low_edge <= float(hist.GetBinCenter(index)) <= high_edge for low_edge, high_edge in windows or [])
    ]


def _alignment_support_metrics(hist, active_bin_indices, metric_name="positive_integral"):
    """Summarize residual support without allowing positive/negative cancellation."""
    signed_integral = 0.0
    absolute_integral = 0.0
    positive_integral = 0.0
    for bin_index in active_bin_indices or []:
        value = float(hist.GetBinContent(bin_index))
        if not _is_finite_number(value):
            continue
        signed_integral += value
        absolute_integral += abs(value)
        positive_integral += max(value, 0.0)
    values = {
        "signed_integral": float(signed_integral),
        "absolute_integral": float(absolute_integral),
        "positive_integral": float(positive_integral),
    }
    resolved_metric = str(metric_name or "positive_integral")
    if resolved_metric not in values:
        resolved_metric = "positive_integral"
    return {
        "signed_support_integral": values["signed_integral"],
        "absolute_support_integral": values["absolute_integral"],
        "positive_support_integral": values["positive_integral"],
        "support_metric_used": resolved_metric,
        "support_integral_for_acceptance": values[resolved_metric],
        "support_bin_count": int(len(active_bin_indices or [])),
    }


def _build_alignment_template(source_hist, total_shift_gev, output_name, config):
    return build_shifted_template_histogram(
        source_hist,
        float(total_shift_gev),
        "positive_moves_peak_higher_mm",
        output_name,
        interpolation_mode=str((config or {}).get("interpolation_mode") or "linear"),
        renormalize=bool((config or {}).get("renormalize_shifted_templates", True)),
    )


def _score_alignment_template(target_hist, template_hist, amplitude, evaluation_indices):
    chi2_value, valid_count = 0.0, 0
    for bin_index in evaluation_indices:
        data_value = float(target_hist.GetBinContent(bin_index))
        model_value = float(template_hist.GetBinContent(bin_index)) * float(amplitude)
        sigma = float(target_hist.GetBinError(bin_index))
        if not _is_finite_number(sigma) or sigma <= 0.0:
            sigma = max(math.sqrt(abs(data_value)), 1.0)
        if not (_is_finite_number(data_value) and _is_finite_number(model_value)):
            continue
        chi2_value += ((data_value - model_value) / sigma) ** 2
        valid_count += 1
    ndf = max(valid_count - 1, 1)
    return {
        "chi2_eval": float(chi2_value),
        "ndf_eval": int(ndf),
        "chi2_eval_ndf": float(chi2_value / ndf) if valid_count else None,
        "evaluation_bin_count": int(valid_count),
    }


def _alignment_relative_improvement(baseline_score, candidate_score):
    if not (_is_finite_number(baseline_score) and _is_finite_number(candidate_score)):
        return None
    denominator = abs(float(baseline_score))
    return 0.0 if denominator <= 1e-12 else float((float(baseline_score) - float(candidate_score)) / denominator)


def _pion_control_histogram_identifier(hist):
    """Return the immutable identity used to validate a saved calibration."""
    if hist is None:
        return None
    return {
        "hist_name": str(hist.GetName()),
        "axis": _hist_axis_specification(hist),
        "checksum": _hist_content_checksum(hist),
    }


def _alignment_canonical_windows(inp_dict, phi_setting):
    resolved = resolve_particle_subtraction_component_fit_windows(
        "pion_control", mm_offset_data=0.0, inp_dict=inp_dict, phi_setting=phi_setting
    )
    return {component: _normalize_window_collection(resolved.get(component)) for component in COMPONENT_NAMES}


def _alignment_evaluation_envelope(canonical_windows, common_shift, residual_values, expansions, mm_limits, parent_windows=None):
    residual_values = [float(value) for value in residual_values if _is_finite_number(value)] or [0.0]
    expansion_values = [max(float(value), 0.0) for value in expansions if _is_finite_number(value)] or [0.0]
    maximum_expansion = max(expansion_values)
    envelope = [
        (
            low_edge + float(common_shift) + min(residual_values) - maximum_expansion,
            high_edge + float(common_shift) + max(residual_values) + maximum_expansion,
        )
        for low_edge, high_edge in _normalize_window_collection(canonical_windows)
    ]
    envelope.extend(_normalize_window_collection(parent_windows))
    return transform_component_windows(envelope, 0.0, 0.0, mm_limits)["effective_windows"]


def scan_pion_component_alignment(
    target_hist,
    residual_hist,
    source_template_hist,
    component_name,
    canonical_windows,
    common_shift_gev,
    alignment_config,
    analysis_scope,
    bin_key,
    parent_component_alignment=None,
    scan_level="parent",
):
    """Find one pion component's best offset/expansion from immutable SIMC."""
    config = alignment_config or {}
    component_config = (config.get("components") or {}).get(component_name, {}) or {}
    acceptance, ranking = config.get("acceptance") or {}, config.get("ranking") or {}
    working_hist = residual_hist or target_hist
    mm_limits = (float(working_hist.GetXaxis().GetXmin()), float(working_hist.GetXaxis().GetXmax()))
    parent = deepcopy(parent_component_alignment or {})
    parent_valid = bool(parent.get("accepted", False))
    parent_residual = float(parent.get("residual_shift_gev", 0.0) or 0.0)
    parent_expansion = float(parent.get("window_expansion_gev", 0.0) or 0.0)
    if scan_level == "fine" and not parent_valid:
        return {"accepted": False, "source": "current_common_shift_fallback", "rejection_reasons": ["parent unavailable"], **_immutable_source_metadata(source_template_hist, component_name)}

    if scan_level == "fine":
        spec = component_config.get("fine_bin_offset_scan_gev") or {}
        corrections = _alignment_scan_grid(spec.get("minimum_relative_to_parent", 0.0), spec.get("maximum_relative_to_parent", 0.0), spec.get("step", 1.0), 0.0)
        residual_values = [parent_residual + correction for correction in corrections]
        baseline = (parent_residual, parent_expansion, 0.0)
        baseline_source = "parent_setting_alignment"
    else:
        spec = component_config.get("global_offset_scan_gev") or {}
        residual_values = _alignment_scan_grid(spec.get("minimum", 0.0), spec.get("maximum", 0.0), spec.get("step", 1.0), 0.0)
        corrections = list(residual_values)
        baseline = (0.0, 0.0, 0.0)
        baseline_source = "current_common_shift"
    expansions = sorted({float(value) for value in component_config.get("window_expansion_candidates_gev") or [0.0]})
    if scan_level == "fine":
        expansions = sorted(set(expansions + [parent_expansion]))
    envelope = _alignment_evaluation_envelope(canonical_windows, common_shift_gev, residual_values + [baseline[0]], expansions + [baseline[1]], mm_limits, parent.get("effective_windows"))
    evaluation_indices = _alignment_windows_to_bin_indices(working_hist, envelope)

    candidate_points = [(baseline[0], baseline[1], baseline[2], True)]
    if scan_level == "fine":
        candidate_points.extend((parent_residual + correction, expansion, correction, False) for correction in corrections for expansion in expansions)
    else:
        candidate_points.extend((residual, expansion, residual, False) for residual in residual_values for expansion in expansions)
    candidates, seen = [], set()
    for residual, expansion, correction, is_baseline in candidate_points:
        key = (round(residual, 12), round(expansion, 12))
        if key in seen:
            continue
        seen.add(key)
        total_shift = float(common_shift_gev) + float(residual)
        template, shift_diagnostics = _build_alignment_template(
            source_template_hist,
            total_shift,
            "pion_align_{}_{}_{}_r{:+.6f}_e{:.6f}".format(analysis_scope, bin_key, component_name, residual, expansion).replace(" ", "_"),
            config,
        )
        windows = transform_component_windows(canonical_windows, total_shift, expansion, mm_limits)
        # The fine-bin parent baseline is a numerical replay of the resolved
        # parent map, including a non-zero expansion that may not be in the
        # local scan list.  It is intentionally not regenerated from a default
        # window here.
        if scan_level == "fine" and is_baseline and parent.get("effective_windows"):
            windows["effective_windows"] = _normalize_window_collection(parent.get("effective_windows"))
            windows["window_clipping"] = deepcopy(parent.get("window_clipping") or [])
        amplitude = _solve_nonnegative_template_amplitude(
            working_hist, template, mm_limits[0], mm_limits[1],
            include_windows=windows["effective_windows"], amplitude_windows=windows["effective_windows"],
        )
        active_fit_indices = _alignment_windows_to_bin_indices(
            working_hist, windows["effective_windows"]
        )
        support = _alignment_support_metrics(
            working_hist,
            active_fit_indices,
            (acceptance or {}).get("data_support_metric", "positive_integral"),
        )
        support_threshold = float(acceptance.get("minimum_data_integral", 0.0) or 0.0)
        evaluation = _score_alignment_template(working_hist, template, amplitude.get("amplitude", 0.0), evaluation_indices)
        offset_boundary = abs(residual - min(residual_values)) <= 1e-12 or abs(residual - max(residual_values)) <= 1e-12
        expansion_boundary = abs(expansion - min(expansions)) <= 1e-12 or abs(expansion - max(expansions)) <= 1e-12
        penalize_offset_boundary = bool(
            offset_boundary and bool(acceptance.get("reject_offset_boundary", False))
        )
        penalize_expansion_boundary = bool(
            expansion_boundary and bool(acceptance.get("reject_expansion_boundary", False))
        )
        offset_boundary_penalty = float(
            ranking.get("offset_boundary_penalty", ranking.get("boundary_penalty", 0.0)) or 0.0
        )
        expansion_boundary_penalty = float(
            ranking.get("expansion_boundary_penalty", ranking.get("boundary_penalty", 0.0)) or 0.0
        )
        applied_offset_boundary_penalty = 0.0
        applied_expansion_boundary_penalty = 0.0
        score = evaluation.get("chi2_eval_ndf")
        if _is_finite_number(score):
            score = float(score) + float(ranking.get("window_expansion_penalty", 0.0) or 0.0) * abs(expansion) + float(ranking.get("offset_magnitude_penalty", 0.0) or 0.0) * abs(residual)
            if not is_baseline and penalize_offset_boundary:
                applied_offset_boundary_penalty = offset_boundary_penalty
                score += applied_offset_boundary_penalty
            if not is_baseline and penalize_expansion_boundary:
                applied_expansion_boundary_penalty = expansion_boundary_penalty
                score += applied_expansion_boundary_penalty
        reasons = []
        if not amplitude.get("success"):
            reasons.append(amplitude.get("message") or "amplitude fit failed")
        if int(amplitude.get("n_fit_bins", 0) or 0) < int(acceptance.get("minimum_active_fit_bins", 0) or 0):
            reasons.append("insufficient active fit bins")
        if int(evaluation.get("evaluation_bin_count", 0) or 0) < int(acceptance.get("minimum_evaluation_bins", 0) or 0):
            reasons.append("insufficient evaluation bins")
        if support["support_integral_for_acceptance"] < support_threshold:
            reasons.append("insufficient {} support".format(support["support_metric_used"]))
        if _hist_integral(template) < float(acceptance.get("minimum_template_integral", 0.0) or 0.0):
            reasons.append("insufficient template integral")
        if float(shift_diagnostics.get("lost_integral_fraction", 0.0) or 0.0) > float(acceptance.get("maximum_lost_template_integral_fraction", 1.0) or 1.0):
            reasons.append("lost template integral exceeds threshold")
        if bool(acceptance.get("reject_offset_boundary", False)) and offset_boundary and not is_baseline:
            reasons.append("offset scan boundary")
        if bool(acceptance.get("reject_expansion_boundary", False)) and expansion_boundary and not is_baseline:
            reasons.append("expansion scan boundary")
        if not _is_finite_number(score):
            reasons.append("nonfinite score")
        warnings = []
        if expansion_boundary and not bool(acceptance.get("reject_expansion_boundary", False)):
            warnings.append("expansion scan boundary accepted without penalty")
        candidates.append({
            "residual_shift_gev": float(residual), "total_shift_gev": total_shift,
            "window_expansion_gev": float(expansion), "local_bin_correction_gev": float(correction),
            "is_baseline": bool(is_baseline), "amplitude": float(amplitude.get("amplitude", 0.0) or 0.0),
            "amplitude_uncertainty": amplitude.get("sigma"), "fit_bin_count": int(amplitude.get("n_fit_bins", 0) or 0),
            "evaluation_envelope": deepcopy(envelope), "effective_windows": deepcopy(windows["effective_windows"]),
            "window_clipping": deepcopy(windows["window_clipping"]),
            "lost_template_integral_fraction": float(shift_diagnostics.get("lost_integral_fraction", 0.0) or 0.0),
            "offset_boundary_hit": bool(offset_boundary), "expansion_boundary_hit": bool(expansion_boundary),
            "offset_boundary_penalty_applied": float(applied_offset_boundary_penalty),
            "expansion_boundary_penalty_applied": float(applied_expansion_boundary_penalty),
            "support_threshold": support_threshold, "warnings": warnings,
            **support, **evaluation, "score": score, "rejection_reasons": reasons,
        })
    baseline_candidate = next((candidate for candidate in candidates if candidate["is_baseline"]), None)
    valid = [candidate for candidate in candidates if not candidate["rejection_reasons"]]
    selected = min(valid or candidates, key=lambda candidate: float(candidate["score"]) if _is_finite_number(candidate.get("score")) else float("inf"))
    tolerance = float(acceptance.get("near_minimum_relative_score_tolerance", 0.01) or 0.0)
    near = [candidate for candidate in candidates if _is_finite_number(candidate.get("score")) and _is_finite_number(selected.get("score")) and float(candidate["score"]) <= float(selected["score"]) + max(abs(float(selected["score"])) * tolerance, 1e-12)]
    localized = (len(near) / float(len(candidates))) <= float(acceptance.get("maximum_near_minimum_candidate_fraction", 1.0) or 1.0)
    return {
        "accepted": bool(not selected.get("rejection_reasons") and (localized or not bool(acceptance.get("require_localized_minimum", True)))),
        "source": "setting_global_scan" if scan_level == "parent" else "bin_local_scan",
        "baseline_source": baseline_source, "baseline_score": baseline_candidate.get("score") if baseline_candidate else None,
        "baseline_residual_shift_gev": baseline[0], "baseline_total_shift_gev": float(common_shift_gev) + baseline[0], "baseline_window_expansion_gev": baseline[1],
        "parent_baseline_candidate_included": bool(scan_level == "fine"),
        "canonical_windows": _normalize_window_collection(canonical_windows), "evaluation_envelope": envelope,
        "selected_candidate": selected, "candidate_summaries": candidates, "candidate_count": len(candidates),
        "near_minimum_candidate_count": len(near), "near_minimum_candidate_fraction": len(near) / float(len(candidates)), "minimum_localized": bool(localized),
        **_immutable_source_metadata(source_template_hist, component_name),
    }


def _score_complete_alignment(target_hist, source_hists, components, config, evaluation_windows=None):
    residual = _clone_hist(target_hist, "{}_alignment_score".format(target_hist.GetName()))
    if residual is None:
        return {"score": None, "chi2": None, "ndf": 0}
    score_windows = _normalize_window_collection(evaluation_windows)
    use_component_envelopes = not bool(score_windows)
    for component_name in COMPONENT_NAMES:
        component = (components or {}).get(component_name) or {}
        source_hist = (source_hists or {}).get(component_name)
        if source_hist is None:
            continue
        windows = _normalize_window_collection(component.get("effective_windows"))
        if not windows:
            continue
        template, _ = _build_alignment_template(source_hist, component.get("total_shift_gev", 0.0), "{}_{}".format(residual.GetName(), component_name), config)
        amplitude = _solve_nonnegative_template_amplitude(residual, template, float(residual.GetXaxis().GetXmin()), float(residual.GetXaxis().GetXmax()), include_windows=windows, amplitude_windows=windows)
        scaled = _clone_hist(template, "{}_scaled".format(template.GetName()))
        if scaled is not None:
            scaled.Scale(float(amplitude.get("amplitude", 0.0) or 0.0))
            residual.Add(scaled, -1.0)
        if use_component_envelopes:
            score_windows.extend(_normalize_window_collection(component.get("evaluation_envelope") or windows))
    indices = _alignment_windows_to_bin_indices(target_hist, score_windows)
    chi2_value, valid_count = 0.0, 0
    for index in indices:
        sigma = float(target_hist.GetBinError(index))
        if not _is_finite_number(sigma) or sigma <= 0.0:
            sigma = max(math.sqrt(abs(float(target_hist.GetBinContent(index)))), 1.0)
        chi2_value += (float(residual.GetBinContent(index)) / sigma) ** 2
        valid_count += 1
    ndf = max(valid_count - len(COMPONENT_NAMES), 1)
    return {"score": float(chi2_value / ndf) if valid_count else None, "chi2": float(chi2_value), "ndf": int(ndf), "evaluation_bin_indices": indices}


def _alignment_parent_hash(setting_key, common_shift_gev, components, config_hash):
    return _alignment_hash({
        "setting_key": setting_key, "common_shift_gev": float(common_shift_gev), "resolved_configuration_hash": config_hash,
        "components": {name: {key: (components.get(name) or {}).get(key) for key in ("residual_shift_gev", "window_expansion_gev", "effective_windows", "immutable_source_hist_identifier", "immutable_source_hist_checksum")} for name in COMPONENT_NAMES},
    })


def build_expected_pion_alignment_metadata(
    setting_key,
    analysis_scope,
    bin_key,
    pion_control_hist,
    immutable_pion_simc_component_hists,
    parent_alignment=None,
    inp_dict=None,
    phi_setting=None,
    common_setting_shift_gev=0.0,
):
    """Build cache-validation metadata without constructing any shifted candidate."""
    if pion_control_hist is None:
        raise ValueError("pion_control_hist is required to build alignment metadata")
    config = get_pion_component_dynamic_alignment_config(
        inp_dict=inp_dict, phi_setting=phi_setting, setting_key=setting_key
    )
    sources = immutable_pion_simc_component_hists or {}
    parent = parent_alignment or {}
    return {
        "alignment_schema_version": ALIGNMENT_SCHEMA_VERSION,
        "analysis_scope": str(analysis_scope or ""),
        "bin_key": deepcopy(bin_key),
        "complete_physical_bin_identity": deepcopy(bin_key),
        "parent_setting_key": setting_key,
        "common_setting_shift_gev": (
            float(common_setting_shift_gev)
            if _is_finite_number(common_setting_shift_gev)
            else 0.0
        ),
        "parent_alignment_hash": parent.get("parent_alignment_hash"),
        "resolved_config": config,
        "resolved_configuration_hash": _alignment_hash(config),
        "pion_control_histogram_identifier": _pion_control_histogram_identifier(
            pion_control_hist
        ),
        "pion_control_histogram_provenance": "pion_control",
        "histogram_axis_specification": _hist_axis_specification(pion_control_hist),
        "immutable_source_template_identifier": {
            name: deepcopy(_immutable_source_metadata(sources.get(name), name).get("immutable_source_hist_identifier"))
            for name in COMPONENT_NAMES
        },
        "immutable_source_template_checksum": {
            name: _immutable_source_metadata(sources.get(name), name).get("immutable_source_hist_checksum")
            for name in COMPONENT_NAMES
        },
        "components": {
            name: _immutable_source_metadata(sources.get(name), name)
            for name in COMPONENT_NAMES
        },
    }


def resolve_pion_component_alignment(setting_key, analysis_scope, bin_key, pion_control_hist, immutable_pion_simc_component_hists, parent_alignment=None, inp_dict=None, phi_setting=None, common_setting_shift_gev=0.0):
    """Resolve parent/fine pion alignment without ever chaining shifted templates."""
    config = get_pion_component_dynamic_alignment_config(inp_dict=inp_dict, phi_setting=phi_setting, setting_key=setting_key)
    config_hash = _alignment_hash(config)
    scope = str(analysis_scope or "")
    is_parent = scope == "setting-wide"
    common_shift = float(common_setting_shift_gev) if _is_finite_number(common_setting_shift_gev) else 0.0
    sources = immutable_pion_simc_component_hists or {}
    canonical = _alignment_canonical_windows(inp_dict, phi_setting)
    if pion_control_hist is None:
        raise ValueError("pion_control_hist is required to resolve pion-component alignment")
    limits = (float(pion_control_hist.GetXaxis().GetXmin()), float(pion_control_hist.GetXaxis().GetXmax()))
    parent = deepcopy(parent_alignment or {})
    expected_parent_hash = _alignment_parent_hash(
        setting_key,
        common_shift,
        parent.get("components") or {},
        config_hash,
    ) if parent else None
    parent_valid = bool(
        parent.get("accepted", False)
        and parent.get("parent_alignment_hash")
        and parent.get("alignment_schema_version") == ALIGNMENT_SCHEMA_VERSION
        and parent.get("resolved_configuration_hash") == config_hash
        and parent.get("parent_setting_key") == setting_key
        and parent.get("parent_alignment_hash") == expected_parent_hash
    )
    default_components = {
        name: {
            "source": "current_common_shift_fallback", "accepted": False,
            "common_shift_gev": common_shift, "parent_residual_shift_gev": 0.0, "local_bin_correction_gev": 0.0,
            "parent_total_shift_gev": common_shift, "parent_window_expansion_gev": 0.0,
            "residual_shift_gev": 0.0, "total_shift_gev": common_shift, "window_expansion_gev": 0.0,
            **transform_component_windows(canonical.get(name), common_shift, 0.0, limits),
            "evaluation_envelope": transform_component_windows(canonical.get(name), common_shift, 0.0, limits)["effective_windows"],
            "baseline_source": "current_common_shift", "baseline_score": None,
            "baseline_residual_shift_gev": 0.0, "baseline_total_shift_gev": common_shift,
            "baseline_window_expansion_gev": 0.0, "rejection_reasons": [], "warnings": [],
            **_immutable_source_metadata(sources.get(name), name),
        }
        for name in COMPONENT_NAMES
    }
    for component in default_components.values():
        component["parent_effective_windows"] = deepcopy(component.get("effective_windows") or [])
    for component in default_components.values():
        component["calibration_source"] = "pion_control"
        component["candidate_source"] = "immutable_original_pion_simc"
        component["parent_alignment_hash"] = parent.get("parent_alignment_hash") if parent else None
        component["parent_global_score"] = parent.get("selected_score") if parent else None
        component["parent_score_in_current_bin"] = None
        component["local_candidate_score"] = None
        component["relative_improvement_over_parent"] = None
        component["scan_candidate_accepted"] = False
        component["component_applied_source"] = component["source"]
        component["component_fallback_used"] = True
        component["component_fallback_reason"] = "dynamic alignment unavailable"
    if not bool(config.get("enabled", False)) or (not is_parent and not parent_valid):
        reason = "dynamic alignment disabled" if not bool(config.get("enabled", False)) else "parent alignment unavailable, stale, incompatible, or invalid"
        return {
            "alignment_schema_version": ALIGNMENT_SCHEMA_VERSION, "analysis_scope": scope, "bin_key": deepcopy(bin_key), "complete_physical_bin_identity": deepcopy(bin_key),
            "parent_setting_key": setting_key, "source": "component_disabled" if not bool(config.get("enabled", False)) else "current_common_shift_fallback", "accepted": False,
            "common_setting_shift_gev": common_shift, "baseline_source": "current_common_shift",
            "baseline_score": None, "proposed_score": None, "proposed_relative_improvement": None,
            "applied_score": None, "applied_relative_improvement": 0.0,
            "selected_score": None, "relative_score_improvement": 0.0,
            "parent_global_score": parent.get("selected_score") if parent else None,
            "parent_score_in_current_bin": None, "parent_alignment_hash": parent.get("parent_alignment_hash") if parent else None,
            "resolved_config": config, "resolved_configuration_hash": config_hash,
            "pion_control_histogram_identifier": _pion_control_histogram_identifier(pion_control_hist),
            "pion_control_histogram_provenance": "pion_control",
            "histogram_axis_specification": _hist_axis_specification(pion_control_hist),
            "creation_timestamp": datetime.now(timezone.utc).isoformat(),
            "immutable_source_template_identifier": {
                name: deepcopy(component.get("immutable_source_hist_identifier"))
                for name, component in default_components.items()
            },
            "immutable_source_template_checksum": {
                name: component.get("immutable_source_hist_checksum")
                for name, component in default_components.items()
            },
            "proposed_components": deepcopy(default_components), "applied_components": deepcopy(default_components),
            "component_counts": {
                "locally_accepted": 0, "globally_accepted": 0,
                "parent_fallback": 0, "common_shift_fallback": len(COMPONENT_NAMES),
                "disabled": len(COMPONENT_NAMES) if not bool(config.get("enabled", False)) else 0,
            },
            "components": default_components, "persistence_status": "not_persisted", "persistence_rejection_reasons": [reason],
        }

    baseline_components = deepcopy(default_components) if is_parent else deepcopy(parent.get("components") or {})
    if not is_parent:
        # A parent contributes parameters only.  Rebind provenance to the current
        # scope's immutable raw templates without ever borrowing a shifted parent
        # histogram.
        for name in COMPONENT_NAMES:
            baseline_component = baseline_components.get(name)
            if not isinstance(baseline_component, dict):
                baseline_component = deepcopy(default_components[name])
                baseline_components[name] = baseline_component
            baseline_component.update(_immutable_source_metadata(sources.get(name), name))
            baseline_component["parent_residual_shift_gev"] = float(baseline_component.get("residual_shift_gev", 0.0) or 0.0)
            baseline_component["parent_total_shift_gev"] = float(baseline_component.get("total_shift_gev", common_shift) or common_shift)
            baseline_component["parent_window_expansion_gev"] = float(baseline_component.get("window_expansion_gev", 0.0) or 0.0)
            baseline_component["parent_effective_windows"] = deepcopy(baseline_component.get("effective_windows") or [])
            baseline_component["local_bin_correction_gev"] = 0.0
            baseline_component["baseline_source"] = "parent_setting_alignment"
            baseline_component["baseline_residual_shift_gev"] = float(baseline_component.get("residual_shift_gev", 0.0) or 0.0)
            baseline_component["baseline_total_shift_gev"] = float(baseline_component.get("total_shift_gev", common_shift) or common_shift)
            baseline_component["baseline_window_expansion_gev"] = float(baseline_component.get("window_expansion_gev", 0.0) or 0.0)
            baseline_component["parent_baseline_candidate_included"] = True

    def _candidate_summary(component):
        return {
            key: deepcopy(component.get(key))
            for key in (
                "source", "residual_shift_gev", "total_shift_gev",
                "window_expansion_gev", "local_bin_correction_gev",
                "effective_windows", "evaluation_envelope", "score",
                "offset_boundary_hit", "expansion_boundary_hit",
                "offset_boundary_penalty_applied", "expansion_boundary_penalty_applied",
                "signed_support_integral", "absolute_support_integral",
                "positive_support_integral", "support_metric_used",
                "support_integral_for_acceptance", "support_threshold",
                "rejection_reasons", "warnings",
            )
        }

    proposed_components = {}
    mixed_components = {}
    residual = _clone_hist(pion_control_hist, "{}_alignment_residual".format(pion_control_hist.GetName()))
    for name in COMPONENT_NAMES:
        component_config = (config.get("components") or {}).get(name, {}) or {}
        source_hist = sources.get(name)
        parent_component = (parent.get("components") or {}).get(name, {})
        fallback = deepcopy(baseline_components.get(name) or default_components[name])
        fallback_source = (
            "parent_setting_fallback"
            if not is_parent and parent_component
            else "current_common_shift_fallback"
        )
        fallback["source"] = fallback_source
        fallback["component_applied_source"] = fallback_source
        fallback["component_fallback_used"] = True
        fallback["component_fallback_reason"] = "component scan did not supply an accepted candidate"
        fallback["scan_candidate_accepted"] = False
        fallback.update(_immutable_source_metadata(source_hist, name))

        if source_hist is None:
            proposed = deepcopy(fallback)
            proposed["source"] = "current_common_shift_fallback"
            proposed["rejection_reasons"] = ["missing immutable source template"]
            applied = deepcopy(proposed)
            applied["component_fallback_reason"] = "missing immutable source template"
        elif not bool(component_config.get("enabled", False)):
            proposed = deepcopy(fallback)
            if not is_parent and parent_component:
                proposed["source"] = "parent_setting_fallback"
                proposed["component_fallback_reason"] = "component disabled; inherited parent alignment"
            else:
                proposed["source"] = "component_disabled"
                proposed["component_applied_source"] = "component_disabled"
                proposed["component_fallback_used"] = False
                proposed["component_fallback_reason"] = "component disabled"
            proposed["component_disabled"] = True
            applied = deepcopy(proposed)
        else:
            scan = scan_pion_component_alignment(
                pion_control_hist, residual, source_hist, name, canonical.get(name),
                common_shift, config, scope, bin_key,
                parent_component_alignment=parent_component,
                scan_level="parent" if is_parent else "fine",
            )
            candidate = deepcopy(scan.get("selected_candidate") or {})
            proposed = {
                "source": scan.get("source"), "accepted": bool(scan.get("accepted", False)),
                "common_shift_gev": common_shift,
                "parent_residual_shift_gev": float(parent_component.get("residual_shift_gev", 0.0) or 0.0),
                "local_bin_correction_gev": float(candidate.get("local_bin_correction_gev", 0.0) or 0.0),
                "parent_total_shift_gev": float(parent_component.get("total_shift_gev", common_shift) or common_shift),
                "parent_window_expansion_gev": float(parent_component.get("window_expansion_gev", 0.0) or 0.0),
                "parent_effective_windows": deepcopy(parent_component.get("effective_windows") or []),
                "residual_shift_gev": float(candidate.get("residual_shift_gev", 0.0) or 0.0),
                "total_shift_gev": float(candidate.get("total_shift_gev", common_shift) or common_shift),
                "window_expansion_gev": float(candidate.get("window_expansion_gev", 0.0) or 0.0),
                "canonical_windows": deepcopy(scan.get("canonical_windows") or []),
                "effective_windows": deepcopy(candidate.get("effective_windows") or []),
                "evaluation_envelope": deepcopy(scan.get("evaluation_envelope") or []),
                "baseline_source": scan.get("baseline_source"), "baseline_score": scan.get("baseline_score"),
                "baseline_residual_shift_gev": scan.get("baseline_residual_shift_gev"),
                "baseline_total_shift_gev": scan.get("baseline_total_shift_gev"),
                "baseline_window_expansion_gev": scan.get("baseline_window_expansion_gev"),
                "parent_baseline_candidate_included": scan.get("parent_baseline_candidate_included", False),
                "minimum_localized": scan.get("minimum_localized", False),
                "candidate_count": scan.get("candidate_count"),
                "near_minimum_candidate_count": scan.get("near_minimum_candidate_count"),
                "near_minimum_candidate_fraction": scan.get("near_minimum_candidate_fraction"),
                "candidate_summaries": deepcopy(scan.get("candidate_summaries") or []),
                "scan_candidate_accepted": bool(scan.get("accepted", False)),
                "component_applied_source": scan.get("source"),
                "component_fallback_used": False, "component_fallback_reason": None,
                "score": candidate.get("score"), "proposed_score": candidate.get("score"),
                "rejection_reasons": list(candidate.get("rejection_reasons") or []),
                "warnings": list(candidate.get("warnings") or []),
                **{
                    key: candidate.get(key)
                    for key in (
                        "offset_boundary_hit", "expansion_boundary_hit",
                        "offset_boundary_penalty_applied", "expansion_boundary_penalty_applied",
                        "signed_support_integral", "absolute_support_integral",
                        "positive_support_integral", "support_metric_used",
                        "support_integral_for_acceptance", "support_threshold",
                    )
                },
                **_immutable_source_metadata(source_hist, name),
            }
            if proposed["scan_candidate_accepted"]:
                applied = deepcopy(proposed)
            else:
                applied = deepcopy(fallback)
                applied["component_fallback_reason"] = "; ".join(
                    proposed.get("rejection_reasons") or ["candidate did not pass scan acceptance"]
                )
                applied["proposed_score"] = proposed.get("proposed_score")
                applied["score"] = scan.get("baseline_score")
                applied["candidate_summaries"] = deepcopy(proposed.get("candidate_summaries") or [])
                applied["warnings"] = list(proposed.get("warnings") or [])
                applied["scan_candidate_accepted"] = False

        proposed["proposed_candidate"] = _candidate_summary(proposed)
        proposed_components[name] = deepcopy(proposed)
        applied["proposed_candidate"] = deepcopy(proposed["proposed_candidate"])
        applied["applied_candidate"] = _candidate_summary(applied)
        applied["applied_score"] = (
            applied["applied_candidate"].get("score")
            if applied["applied_candidate"].get("score") is not None
            else applied.get("baseline_score")
        )
        mixed_components[name] = applied

        if source_hist is None:
            continue
        template, _ = _build_alignment_template(
            source_hist, applied.get("total_shift_gev", common_shift),
            "{}_applied_{}".format(residual.GetName(), name), config,
        )
        amplitude = _solve_nonnegative_template_amplitude(
            residual, template, limits[0], limits[1],
            include_windows=applied.get("effective_windows"),
            amplitude_windows=applied.get("effective_windows"),
        )
        scaled = _clone_hist(template, "{}_scaled".format(template.GetName()))
        if scaled is not None:
            scaled.Scale(float(amplitude.get("amplitude", 0.0) or 0.0))
            residual.Add(scaled, -1.0)

    # Candidate, mixed, and baseline maps share the same fixed envelope.
    fixed_evaluation_windows = []
    for component in proposed_components.values():
        fixed_evaluation_windows.extend(_normalize_window_collection(component.get("evaluation_envelope")))
    baseline_metrics = _score_complete_alignment(
        pion_control_hist, sources, baseline_components, config, fixed_evaluation_windows
    )
    proposed_metrics = _score_complete_alignment(
        pion_control_hist, sources, proposed_components, config, fixed_evaluation_windows
    )
    mixed_metrics = _score_complete_alignment(
        pion_control_hist, sources, mixed_components, config, fixed_evaluation_windows
    )
    proposed_improvement = _alignment_relative_improvement(
        baseline_metrics.get("score"), proposed_metrics.get("score")
    )
    mixed_improvement = _alignment_relative_improvement(
        baseline_metrics.get("score"), mixed_metrics.get("score")
    )
    accepted = bool(
        mixed_improvement is not None
        and mixed_improvement >= float((config.get("acceptance") or {}).get("minimum_relative_score_improvement", 0.0) or 0.0)
    )
    applied_components = deepcopy(mixed_components) if accepted else deepcopy(baseline_components)
    if not accepted:
        for name, component in applied_components.items():
            component["source"] = "current_common_shift_fallback" if is_parent else "parent_setting_fallback"
            component["component_applied_source"] = component["source"]
            component["component_fallback_used"] = True
            component["component_fallback_reason"] = "mixed map did not improve complete baseline"
            component["baseline_source"] = "current_common_shift" if is_parent else "parent_setting_alignment"
            component["baseline_score"] = baseline_metrics.get("score")
            component.update(_immutable_source_metadata(sources.get(name), name))
            component["proposed_candidate"] = _candidate_summary(proposed_components.get(name) or {})
            component["applied_candidate"] = _candidate_summary(component)
            component["proposed_score"] = component["proposed_candidate"].get("score")
            component["applied_score"] = component["applied_candidate"].get("score")

    applied_metrics = mixed_metrics if accepted else baseline_metrics
    applied_improvement = mixed_improvement if accepted else 0.0
    parent_hash = _alignment_parent_hash(setting_key, common_shift, applied_components, config_hash) if is_parent else parent.get("parent_alignment_hash")
    parent_global_score = None if is_parent else parent.get("selected_score")
    parent_score_in_current_bin = None if is_parent else baseline_metrics.get("score")
    for component in applied_components.values():
        component["calibration_source"] = "pion_control"
        component["candidate_source"] = "immutable_original_pion_simc"
        component["parent_alignment_hash"] = parent_hash
        component["parent_global_score"] = parent_global_score
        component["parent_score_in_current_bin"] = parent_score_in_current_bin
        component["local_candidate_score"] = component.get("proposed_score")
        component["relative_improvement_over_parent"] = applied_improvement
    component_counts = {
        "locally_accepted": sum(1 for component in applied_components.values() if component.get("component_applied_source") == "bin_local_scan"),
        "globally_accepted": sum(1 for component in applied_components.values() if component.get("component_applied_source") == "setting_global_scan"),
        "parent_fallback": sum(1 for component in applied_components.values() if component.get("component_applied_source") == "parent_setting_fallback"),
        "common_shift_fallback": sum(1 for component in applied_components.values() if component.get("component_applied_source") == "current_common_shift_fallback"),
        "disabled": sum(1 for component in applied_components.values() if component.get("component_applied_source") == "component_disabled"),
    }
    return {
        "alignment_schema_version": ALIGNMENT_SCHEMA_VERSION, "analysis_scope": scope, "bin_key": deepcopy(bin_key), "complete_physical_bin_identity": deepcopy(bin_key), "parent_setting_key": setting_key,
        "source": "setting_global_scan" if is_parent and accepted else ("bin_local_scan" if accepted else ("current_common_shift_fallback" if is_parent else "parent_setting_fallback")), "accepted": accepted,
        "common_setting_shift_gev": common_shift, "baseline_score": baseline_metrics.get("score"),
        "proposed_score": proposed_metrics.get("score"), "proposed_relative_improvement": proposed_improvement,
        "mixed_score": mixed_metrics.get("score"), "mixed_relative_improvement": mixed_improvement,
        "applied_score": applied_metrics.get("score"), "applied_relative_improvement": applied_improvement,
        "selected_score": applied_metrics.get("score"), "relative_score_improvement": applied_improvement,
        "parent_global_score": parent_global_score, "parent_score_in_current_bin": parent_score_in_current_bin, "parent_alignment_hash": parent_hash,
        "resolved_config": config, "resolved_configuration_hash": config_hash,
        "pion_control_histogram_identifier": _pion_control_histogram_identifier(pion_control_hist), "pion_control_histogram_provenance": "pion_control", "histogram_axis_specification": _hist_axis_specification(pion_control_hist), "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "immutable_source_template_identifier": {
            name: deepcopy(component.get("immutable_source_hist_identifier"))
            for name, component in applied_components.items()
        },
        "immutable_source_template_checksum": {
            name: component.get("immutable_source_hist_checksum")
            for name, component in applied_components.items()
        },
        "proposed_components": proposed_components, "applied_components": deepcopy(applied_components),
        "component_counts": component_counts,
        "components": applied_components, "persistence_status": "not_persisted", "persistence_rejection_reasons": [],
    }


def _alignment_safe_filename(value):
    return "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in str(value or "unknown"))


def get_pion_component_alignment_paths(outpath, setting_key, phi_setting, epsset):
    stem = "pion_component_alignment_{}_{}_{}".format(
        _alignment_safe_filename(setting_key),
        _alignment_safe_filename(phi_setting),
        _alignment_safe_filename(epsset),
    )
    return {
        "json": os.path.join(str(outpath), "{}.json".format(stem)),
        "csv": os.path.join(str(outpath), "{}.csv".format(stem)),
    }


def _alignment_record_key(payload):
    return _alignment_hash({
        "analysis_scope": payload.get("analysis_scope"),
        "identity": payload.get("complete_physical_bin_identity"),
    })


def _alignment_template_metadata_map(payload):
    return {
        name: {
            "identifier": (payload.get("components", {}).get(name) or {}).get("immutable_source_hist_identifier"),
            "checksum": (payload.get("components", {}).get(name) or {}).get("immutable_source_hist_checksum"),
        }
        for name in COMPONENT_NAMES
    }


def _alignment_compatibility_reasons(stored, expected):
    reasons = []
    for field in (
        "alignment_schema_version",
        "resolved_configuration_hash",
        "analysis_scope",
        "complete_physical_bin_identity",
        "histogram_axis_specification",
    ):
        if stored.get(field) != expected.get(field):
            reasons.append("{} mismatch".format(field))
    expected_parent_hash = expected.get("parent_alignment_hash")
    if expected_parent_hash and stored.get("parent_alignment_hash") != expected_parent_hash:
        reasons.append("parent_alignment_hash mismatch")
    if _alignment_template_metadata_map(stored) != _alignment_template_metadata_map(expected):
        reasons.append("immutable source template identifier/checksum mismatch")
    if stored.get("pion_control_histogram_identifier") != expected.get("pion_control_histogram_identifier"):
        reasons.append("pion control histogram identifier mismatch")
    return reasons


def load_pion_component_alignment(outpath, setting_key, phi_setting, epsset, expected_payload):
    """Load an exact compatible persisted map; reject stale payloads explicitly."""
    paths = get_pion_component_alignment_paths(outpath, setting_key, phi_setting, epsset)
    if not os.path.isfile(paths["json"]):
        return None, ["alignment store does not exist"], paths
    try:
        with open(paths["json"], "r") as handle:
            store = json.load(handle)
    except (OSError, ValueError) as exc:
        return None, ["alignment store unreadable: {}".format(exc)], paths
    stored = (store.get("records") or {}).get(_alignment_record_key(expected_payload))
    if not isinstance(stored, dict):
        return None, ["complete physical bin identity not found"], paths
    reasons = _alignment_compatibility_reasons(stored, expected_payload)
    if reasons:
        return None, reasons, paths
    loaded = deepcopy(stored)
    loaded["persistence_status"] = "reused"
    loaded["persistence_rejection_reasons"] = []
    return loaded, [], paths


def persist_pion_component_alignment(outpath, setting_key, phi_setting, epsset, payload):
    """Atomically persist parent/fine maps and rewrite compact CSV diagnostics."""
    if not isinstance(payload, dict) or not outpath:
        return []
    paths = get_pion_component_alignment_paths(outpath, setting_key, phi_setting, epsset)
    os.makedirs(str(outpath), exist_ok=True)
    store = {"alignment_schema_version": ALIGNMENT_SCHEMA_VERSION, "records": {}}
    if os.path.isfile(paths["json"]):
        try:
            with open(paths["json"], "r") as handle:
                existing = json.load(handle)
            if isinstance(existing, dict):
                store.update(existing)
                store.setdefault("records", {})
        except (OSError, ValueError):
            pass
    persisted = deepcopy(payload)
    requested_status = str(persisted.get("persistence_status") or "")
    persisted["persistence_status"] = (
        requested_status
        if requested_status in {"created", "rejected_stale_then_created"}
        else "created"
    )
    persisted["persistence_rejection_reasons"] = list(
        persisted.get("persistence_rejection_reasons") or []
    )
    store["alignment_schema_version"] = ALIGNMENT_SCHEMA_VERSION
    store["records"][_alignment_record_key(persisted)] = persisted
    descriptor, temporary_path = tempfile.mkstemp(prefix=".pion_alignment_", suffix=".json", dir=str(outpath))
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(store, handle, indent=2, sort_keys=True)
        os.replace(temporary_path, paths["json"])
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)

    rows = []
    for record in (store.get("records") or {}).values():
        for component_name, component in (record.get("components") or {}).items():
            rows.append({
                "alignment_schema_version": record.get("alignment_schema_version"),
                "resolved_configuration_hash": record.get("resolved_configuration_hash"),
                "setting": record.get("parent_setting_key"), "scope": record.get("analysis_scope"),
                "bin_key": _canonical_json(record.get("complete_physical_bin_identity")), "component": component_name,
                "selected_source": component.get("source"), "calibration_source": component.get("calibration_source", "pion_control"),
                "parent_alignment_hash": record.get("parent_alignment_hash"),
                "parent_global_score": component.get("parent_global_score", record.get("selected_score")),
                "parent_score_in_current_bin": component.get("parent_score_in_current_bin", record.get("parent_score_in_current_bin")),
                "baseline_source": component.get("baseline_source"), "baseline_score": component.get("baseline_score", record.get("baseline_score")),
                "proposed_score": component.get("proposed_score", record.get("proposed_score")),
                "applied_score": component.get("applied_score", record.get("applied_score")),
                "proposed_relative_improvement": record.get("proposed_relative_improvement"),
                "applied_relative_improvement": record.get("applied_relative_improvement"),
                "scan_candidate_accepted": component.get("scan_candidate_accepted"),
                "component_applied_source": component.get("component_applied_source", component.get("source")),
                "component_fallback_used": component.get("component_fallback_used"),
                "component_fallback_reason": component.get("component_fallback_reason"),
                "residual_shift_gev": component.get("residual_shift_gev"), "total_shift_gev": component.get("total_shift_gev"),
                "window_expansion_gev": component.get("window_expansion_gev"), "local_candidate_score": component.get("local_candidate_score", record.get("selected_score")),
                "relative_improvement_over_parent": component.get("relative_improvement_over_parent", record.get("relative_score_improvement")),
                "offset_boundary_hit": component.get("offset_boundary_hit"),
                "expansion_boundary_hit": component.get("expansion_boundary_hit"),
                "offset_boundary_penalty_applied": component.get("offset_boundary_penalty_applied"),
                "expansion_boundary_penalty_applied": component.get("expansion_boundary_penalty_applied"),
                "signed_support_integral": component.get("signed_support_integral"),
                "absolute_support_integral": component.get("absolute_support_integral"),
                "positive_support_integral": component.get("positive_support_integral"),
                "support_metric_used": component.get("support_metric_used"),
                "support_threshold": component.get("support_threshold"),
                "accepted": component.get("accepted", record.get("accepted")), "persistence_status": record.get("persistence_status"),
                "candidate_source": component.get("candidate_source", "immutable_original_pion_simc"),
                "immutable_source_hist_identifier": _canonical_json(component.get("immutable_source_hist_identifier")),
                "immutable_source_hist_checksum": component.get("immutable_source_hist_checksum"),
                "shift_operation_count": component.get("shift_operation_count"),
                "pion_control_histogram_identifier": _canonical_json(record.get("pion_control_histogram_identifier")),
                "pion_control_histogram_provenance": record.get("pion_control_histogram_provenance"),
                "histogram_axis_specification": _canonical_json(record.get("histogram_axis_specification")),
                "creation_timestamp": record.get("creation_timestamp"),
                "persistence_rejection_reasons": "; ".join(record.get("persistence_rejection_reasons") or []),
            })
    fieldnames = list(rows[0].keys()) if rows else ["setting", "scope", "bin_key", "component"]
    with open(paths["csv"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return [paths["json"], paths["csv"]]


def load_or_resolve_pion_component_alignment(
    outpath,
    setting_key,
    phi_setting,
    epsset,
    analysis_scope,
    bin_key,
    pion_control_hist,
    immutable_pion_simc_component_hists,
    parent_alignment=None,
    inp_dict=None,
    common_setting_shift_gev=0.0,
):
    """Load a compatible map before scanning, or resolve and persist one once."""
    expected = build_expected_pion_alignment_metadata(
        setting_key,
        analysis_scope,
        bin_key,
        pion_control_hist,
        immutable_pion_simc_component_hists,
        parent_alignment=parent_alignment,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        common_setting_shift_gev=common_setting_shift_gev,
    )
    if outpath:
        cached, rejection_reasons, artifact_paths = load_pion_component_alignment(
            outpath, setting_key, phi_setting, epsset, expected
        )
        if cached is not None:
            cached["persistence_status"] = "reused"
            cached["persistence_rejection_reasons"] = []
            return cached, "reused", [], [artifact_paths["json"], artifact_paths["csv"]]
    else:
        rejection_reasons, artifact_paths = ["persistence disabled: no OUTPATH"], {}

    resolved = resolve_pion_component_alignment(
        setting_key,
        analysis_scope,
        bin_key,
        pion_control_hist,
        immutable_pion_simc_component_hists,
        parent_alignment=parent_alignment,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        common_setting_shift_gev=common_setting_shift_gev,
    )
    if not outpath:
        resolved["persistence_status"] = "not_persisted"
        resolved["persistence_rejection_reasons"] = list(rejection_reasons)
        return resolved, "not_persisted", list(rejection_reasons), []

    stale_rejection = any(
        reason not in {"alignment store does not exist", "complete physical bin identity not found"}
        for reason in rejection_reasons
    )
    status = "rejected_stale_then_created" if stale_rejection else "created"
    resolved["persistence_status"] = status
    resolved["persistence_rejection_reasons"] = list(rejection_reasons)
    paths = persist_pion_component_alignment(
        outpath, setting_key, phi_setting, epsset, resolved
    )
    return resolved, status, list(rejection_reasons), paths


def build_particle_subtraction_component_result(
    h_pion_control,
    h_kaon_nosub,
    component_shapes,
    inpDict,
    analysis_scope,
    kaon_signal_shape=None,
    kaon_sigma0_shape=None,
    mm_offset_data=0.0,
    phi_setting=None,
    context="",
    parent_alignment=None,
    pion_component_alignment=None,
    alignment_bin_key=None,
    kaon_sigma0_source_diagnostics=None,
):
    configure_particle_subtraction_root_ownership_debug(
        resolve_particle_subtraction_root_ownership_debug(inpDict)
    )
    mode = resolve_particle_subtraction_mode(inpDict)
    if mode != PARTICLE_SUBTRACTION_MODE_COMPONENTS:
        return {
            "particle_subtraction_mode": mode,
            "analysis_scope": analysis_scope,
            "fallback_used": True,
            "fallback_reason": "particle subtraction mode is not simc_shape_components",
        }

    _require_kaon_lambda_simc_shape(
        kaon_signal_shape,
        analysis_scope,
        phi_setting,
        context,
    )
    template_mm_offset_data = float(mm_offset_data) if _is_finite_number(mm_offset_data) else 0.0
    resolved_phi_setting = phi_setting
    if resolved_phi_setting is None and isinstance(inpDict, dict):
        resolved_phi_setting = inpDict.get("phi_setting")
    setting_key = get_particle_subtraction_setting_key(inpDict)
    dynamic_config = get_pion_component_dynamic_alignment_config(
        inp_dict=inpDict,
        phi_setting=resolved_phi_setting,
        setting_key=setting_key,
    )
    dynamic_alignment_enabled = bool(dynamic_config.get("enabled", False))
    if dynamic_alignment_enabled:
        alignment_payload = pion_component_alignment or resolve_pion_component_alignment(
            setting_key,
            analysis_scope,
            alignment_bin_key or {"analysis_scope": analysis_scope, "context": context},
            h_pion_control,
            component_shapes,
            parent_alignment=parent_alignment,
            inp_dict=inpDict,
            phi_setting=resolved_phi_setting,
            common_setting_shift_gev=template_mm_offset_data,
        )
        aligned_component_shapes = {}
        alignment_windows = {}
        for component_name in COMPONENT_NAMES:
            component_alignment = (alignment_payload.get("components") or {}).get(component_name) or {}
            aligned_component_shapes[component_name], _ = _build_alignment_template(
                component_shapes.get(component_name),
                component_alignment.get("total_shift_gev", template_mm_offset_data),
                "H_MM_component_shape_{}_aligned_{}".format(component_name, context or analysis_scope),
                dynamic_config,
            )
            alignment_windows[component_name] = deepcopy(
                component_alignment.get("effective_windows") or []
            )
    else:
        alignment_payload = None
        alignment_windows = None
        aligned_component_shapes = {
            "pi_n": _build_mm_shifted_hist(
                component_shapes.get("pi_n"),
                template_mm_offset_data,
                "H_MM_component_shape_pi_n_aligned_{}".format(context or analysis_scope),
                renormalize=True,
            ),
            "pi_delta": _build_mm_shifted_hist(
                component_shapes.get("pi_delta"),
                template_mm_offset_data,
                "H_MM_component_shape_pi_delta_aligned_{}".format(context or analysis_scope),
                renormalize=True,
            ),
            "pi_sidis": _build_mm_shifted_hist(
                component_shapes.get("pi_sidis"),
                template_mm_offset_data,
                "H_MM_component_shape_pi_sidis_aligned_{}".format(context or analysis_scope),
                renormalize=True,
            ),
        }
    raw_kaon_signal_reference = _clone_hist(
        kaon_signal_shape,
        "H_k_lambda_simc_input_{}".format(context or analysis_scope),
    )
    aligned_kaon_signal_shape = _build_mm_shifted_hist(
        kaon_signal_shape,
        template_mm_offset_data,
        "H_MM_component_shape_k_lambda_aligned_{}".format(context or analysis_scope),
        renormalize=True,
    )
    aligned_kaon_sigma0_shape = _build_mm_shifted_hist(
        kaon_sigma0_shape,
        template_mm_offset_data,
        "H_MM_component_shape_k_sigma0_aligned_{}".format(context or analysis_scope),
        renormalize=True,
    )
    pion_fit = fit_pion_control_with_simc_shapes(
        h_pion_control,
        aligned_component_shapes.get("pi_n"),
        aligned_component_shapes.get("pi_delta"),
        aligned_component_shapes.get("pi_sidis"),
        aligned_kaon_sigma0_shape,
        inpDict,
        mm_offset_data=mm_offset_data,
        phi_setting=resolved_phi_setting,
        context=context,
        _skip_residual_shift=dynamic_alignment_enabled,
        alignment_windows=alignment_windows,
    )
    pion_fit = _promote_sparse_zero_solution_if_applicable(
        pion_fit,
        "B",
        analysis_scope,
    )
    kaon_fit = fit_kaon_nosub_with_simc_pion_shapes(
        h_kaon_nosub,
        aligned_component_shapes.get("pi_n"),
        aligned_component_shapes.get("pi_delta"),
        aligned_component_shapes.get("pi_sidis"),
        aligned_kaon_signal_shape,
        aligned_kaon_sigma0_shape,
        inpDict,
        mm_offset_data=mm_offset_data,
        phi_setting=resolved_phi_setting,
        context=context,
        _skip_residual_shift=dynamic_alignment_enabled,
        alignment_windows=alignment_windows,
        kaon_sigma0_source_diagnostics=kaon_sigma0_source_diagnostics,
    )
    kaon_fit = _promote_sparse_zero_solution_if_applicable(
        kaon_fit,
        "A",
        analysis_scope,
    )
    pion_shift_payload = pion_fit.get("template_shift_payload") or {}
    kaon_shift_payload = kaon_fit.get("template_shift_payload") or {}
    residual_shift_summaries = {
        "pion_control": deepcopy((pion_shift_payload.get("summary") or {})),
        "kaon_nosub": deepcopy((kaon_shift_payload.get("summary") or {})),
    }

    b_n = float(pion_fit["B_n"])
    b_delta = float(pion_fit["B_delta"])
    b_sidis = float(pion_fit["B_sidis"])
    b_sigma0 = pion_fit["B_sigma0"]
    a_n = float(kaon_fit["A_n"])
    a_delta = float(kaon_fit["A_delta"])
    a_sidis = float(kaon_fit["A_sidis"])
    s_lambda = kaon_fit["S_lambda"]
    s_sigma0 = kaon_fit["S_sigma0"]

    fallback_reasons = []
    if pion_fit["diagnostics"].get("fallback_used"):
        fallback_reasons.append("pion: {}".format(pion_fit["diagnostics"].get("fallback_reason")))
    if kaon_fit["diagnostics"].get("fallback_used"):
        fallback_reasons.append("kaon: {}".format(kaon_fit["diagnostics"].get("fallback_reason")))

    pion_diagnostics = deepcopy(pion_fit["diagnostics"])
    kaon_diagnostics = deepcopy(kaon_fit["diagnostics"])
    pion_fit_mode = str(pion_fit.get("fit_mode") or pion_diagnostics.get("fit_mode") or "staged_only")
    kaon_fit_mode = str(kaon_fit.get("fit_mode") or kaon_diagnostics.get("fit_mode") or "staged_only")
    combined_fit_mode = pion_fit_mode if pion_fit_mode == kaon_fit_mode else "mixed"

    staged_amplitudes_raw = {
        "pion_control": deepcopy(pion_diagnostics.get("staged_amplitudes_raw") or {}),
        "kaon_nosub": deepcopy(kaon_diagnostics.get("staged_amplitudes_raw") or {}),
    }
    staged_amplitudes_scaled = {
        "pion_control": deepcopy(pion_diagnostics.get("staged_amplitudes_scaled") or {}),
        "kaon_nosub": deepcopy(kaon_diagnostics.get("staged_amplitudes_scaled") or {}),
    }
    refined_amplitudes = {
        "pion_control": deepcopy(
            pion_diagnostics.get("refined_amplitudes")
            or pion_diagnostics.get("staged_amplitudes_scaled")
            or {}
        ),
        "kaon_nosub": deepcopy(
            kaon_diagnostics.get("refined_amplitudes")
            or kaon_diagnostics.get("staged_amplitudes_scaled")
            or {}
        ),
    }
    amplitude_shifts = {
        "pion_control": deepcopy(pion_diagnostics.get("amplitude_shifts") or {}),
        "kaon_nosub": deepcopy(kaon_diagnostics.get("amplitude_shifts") or {}),
    }
    amplitude_shift_fractions = {
        "pion_control": deepcopy(pion_diagnostics.get("amplitude_shift_fractions") or {}),
        "kaon_nosub": deepcopy(kaon_diagnostics.get("amplitude_shift_fractions") or {}),
    }
    k_lambda_simc_reference_hist = (
        _clone_hist(
            kaon_fit.get("H_k_lambda_simc_reference"),
            "H_k_lambda_simc_reference_{}".format(context or analysis_scope),
        )
        or _clone_hist(
            aligned_kaon_signal_shape,
            "H_k_lambda_simc_reference_{}".format(context or analysis_scope),
        )
        or _clone_hist(
            raw_kaon_signal_reference,
            "H_k_lambda_simc_reference_{}".format(context or analysis_scope),
        )
    )
    k_lambda_simc_input_loaded = bool(
        kaon_signal_shape is not None
        or aligned_kaon_signal_shape is not None
        or kaon_fit.get("k_lambda_simc_input_loaded", False)
    )
    k_lambda_simc_reference_available = bool(k_lambda_simc_reference_hist is not None)
    # The inner fitter owns a local K-Sigma0 reference.  Publish a distinct
    # shared-builder clone of that exact object, falling back only to the
    # selected aligned template when the inner reference could not be made.
    k_sigma0_fitter_reference_hist = kaon_fit.get("H_k_sigma0_simc_reference")
    k_sigma0_reference_source = (
        k_sigma0_fitter_reference_hist
        if k_sigma0_fitter_reference_hist is not None
        else aligned_kaon_sigma0_shape
    )
    k_sigma0_simc_reference_hist = _clone_hist(
        k_sigma0_reference_source,
        "H_simc_shape_k_sigma0_{}".format(context or analysis_scope),
    )

    result = {
        "particle_subtraction_mode": mode,
        "analysis_scope": analysis_scope,
        "kaon_sigma0_source_diagnostics": deepcopy(
            kaon_sigma0_source_diagnostics or {}
        ),
        "fit_mode": combined_fit_mode,
        "fit_mode_pion": pion_fit_mode,
        "fit_mode_kaon": kaon_fit_mode,
        "joint_refinement_status_pion": pion_diagnostics.get("joint_refinement_status"),
        "joint_refinement_status_kaon": kaon_diagnostics.get("joint_refinement_status"),
        "A_n": a_n,
        "A_delta": a_delta,
        "A_sidis": a_sidis,
        "S_lambda": s_lambda,
        "S_sigma0": s_sigma0,
        "S_lambda_reference_scale": kaon_fit.get("S_lambda_reference_scale"),
        "k_lambda_reference_scale": kaon_fit.get("k_lambda_reference_scale", kaon_fit.get("S_lambda_reference_scale")),
        # Compatibility-only display scale remains separate from the
        # authoritative pre-pi-delta Lambda gauge used by protected fits.
        "lambda_gauge_amplitude": kaon_fit.get("lambda_gauge_amplitude"),
        "lambda_gauge_amplitude_sigma": kaon_fit.get("lambda_gauge_amplitude_sigma"),
        "lambda_display_scale": kaon_fit.get("lambda_display_scale"),
        "k_lambda_fit_amplitude": s_lambda,
        "k_lambda_simc_input_loaded": bool(k_lambda_simc_input_loaded),
        "k_lambda_simc_reference_available": bool(k_lambda_simc_reference_available),
        "k_lambda_simc_reference_source": (
            "immutable_aligned_k_lambda_simc"
            if k_lambda_simc_reference_available
            else ("K-Lambda SIMC reference unavailable" if not k_lambda_simc_input_loaded else "reference_build_failed")
        ),
        "k_lambda_simc_reference_integral": (
            _hist_integral(k_lambda_simc_reference_hist)
            if k_lambda_simc_reference_hist is not None
            else None
        ),
        "B_n": b_n,
        "B_delta": b_delta,
        "B_sidis": b_sidis,
        "B_sigma0": b_sigma0,
        "A_over_B_n": (a_n / b_n) if b_n > 0.0 else None,
        "A_over_B_delta": (a_delta / b_delta) if b_delta > 0.0 else None,
        "A_over_B_sidis": (a_sidis / b_sidis) if b_sidis > 0.0 else None,
        "fit_status_pion": pion_fit["fit_status"],
        "fit_status_kaon": kaon_fit["fit_status"],
        "chi2_pion": pion_fit["diagnostics"].get("chi2"),
        "ndf_pion": pion_fit["diagnostics"].get("ndf"),
        "chi2_ndf_pion": pion_fit["diagnostics"].get("chi2_ndf"),
        "fit_p_value_pion": pion_fit["diagnostics"].get("fit_p_value"),
        "chi2_kaon": kaon_fit["diagnostics"].get("chi2"),
        "ndf_kaon": kaon_fit["diagnostics"].get("ndf"),
        "chi2_ndf_kaon": kaon_fit["diagnostics"].get("chi2_ndf"),
        "fit_p_value_kaon": kaon_fit["diagnostics"].get("fit_p_value"),
        "template_mm_offset_data": template_mm_offset_data,
        "template_mm_shift_applied": bool(abs(template_mm_offset_data) > 1e-12),
        "particle_subtraction_setting_key": setting_key,
        "particle_subtraction_phi_setting": resolved_phi_setting,
        "pion_component_alignment": deepcopy(alignment_payload) if alignment_payload is not None else None,
        "pion_component_alignment_source": (
            alignment_payload.get("source") if isinstance(alignment_payload, dict) else "legacy_common_shift"
        ),
        "pion_component_alignment_imported_by_kaon_fit": bool(alignment_payload is not None),
        "fallback_used": bool(fallback_reasons),
        "fallback_reason": "; ".join(fallback_reasons),
        "residual_component_shifts_enabled": bool(
            (residual_shift_summaries.get("pion_control") or {}).get("enabled")
            or (residual_shift_summaries.get("kaon_nosub") or {}).get("enabled")
        ),
        "residual_component_shift_modes": {
            "pion_control": (residual_shift_summaries.get("pion_control") or {}).get("mode"),
            "kaon_nosub": (residual_shift_summaries.get("kaon_nosub") or {}).get("mode"),
        },
        "residual_component_shift_selection_metrics": {
            "pion_control": (residual_shift_summaries.get("pion_control") or {}).get("selection_metric"),
            "kaon_nosub": (residual_shift_summaries.get("kaon_nosub") or {}).get("selection_metric"),
        },
        "residual_component_shift_units": {
            "pion_control": (residual_shift_summaries.get("pion_control") or {}).get("units"),
            "kaon_nosub": (residual_shift_summaries.get("kaon_nosub") or {}).get("units"),
        },
        "residual_component_shift_values": {
            "pion_control": deepcopy(
                (residual_shift_summaries.get("pion_control") or {}).get("selected_shift_point") or {}
            ),
            "kaon_nosub": deepcopy(
                (residual_shift_summaries.get("kaon_nosub") or {}).get("selected_shift_point") or {}
            ),
        },
        "residual_component_shift_summaries": residual_shift_summaries,
        "staged_amplitudes_raw": staged_amplitudes_raw,
        "staged_amplitudes_scaled": staged_amplitudes_scaled,
        "refined_amplitudes": refined_amplitudes,
        "amplitude_shifts": amplitude_shifts,
        "amplitude_shift_fractions": amplitude_shift_fractions,
        "resolved_subtraction_config": {
            "setting_key": setting_key,
            "phi_setting": resolved_phi_setting,
            "pion_control": deepcopy(pion_fit.get("resolved_config_summary") or {}),
            "kaon_nosub": deepcopy(kaon_fit.get("resolved_config_summary") or {}),
        },
        "diagnostics": {
            "pion": deepcopy(pion_diagnostics),
            "kaon": deepcopy(kaon_diagnostics),
        },
        "H_simc_shape_pi_n": _clone_hist(
            aligned_component_shapes.get("pi_n"),
            "H_simc_shape_pi_n_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_delta": _clone_hist(
            aligned_component_shapes.get("pi_delta"),
            "H_simc_shape_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_sidis": _clone_hist(
            aligned_component_shapes.get("pi_sidis"),
            "H_simc_shape_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_k_lambda_simc_reference": k_lambda_simc_reference_hist,
        "H_simc_shape_k_lambda": _clone_hist(
            k_lambda_simc_reference_hist,
            "H_simc_shape_k_lambda_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_k_sigma0": k_sigma0_simc_reference_hist,
        "H_pi_delta_protected_fit_input": kaon_fit.get("H_pi_delta_protected_fit_input"),
        "H_pi_delta_lambda_gauge": kaon_fit.get("H_pi_delta_lambda_gauge"),
        "H_pi_delta_protected_k_lambda": kaon_fit.get("H_pi_delta_protected_k_lambda"),
        "H_pi_delta_protected_k_sigma0": kaon_fit.get("H_pi_delta_protected_k_sigma0"),
        "H_pi_delta_protected_pi_delta": kaon_fit.get("H_pi_delta_protected_pi_delta"),
        "H_pi_delta_protected_fit_total": kaon_fit.get("H_pi_delta_protected_fit_total"),
        "H_pi_delta_protected_fit_residual": kaon_fit.get("H_pi_delta_protected_fit_residual"),
        "H_pi_delta_protected_after_subtraction": kaon_fit.get("H_pi_delta_protected_after_subtraction"),
        "H_pion_shift_original_pi_n": _clone_shift_payload_hist(
            pion_shift_payload,
            "original_template_hists",
            "pi_n",
            "H_pion_shift_original_pi_n_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_selected_pi_n": _clone_shift_payload_hist(
            pion_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_n",
            "H_pion_shift_selected_pi_n_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_original_pi_delta": _clone_shift_payload_hist(
            pion_shift_payload,
            "original_template_hists",
            "pi_delta",
            "H_pion_shift_original_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_selected_pi_delta": _clone_shift_payload_hist(
            pion_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_delta",
            "H_pion_shift_selected_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_original_pi_sidis": _clone_shift_payload_hist(
            pion_shift_payload,
            "original_template_hists",
            "pi_sidis",
            "H_pion_shift_original_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_selected_pi_sidis": _clone_shift_payload_hist(
            pion_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_sidis",
            "H_pion_shift_selected_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_original_k_sigma0": _clone_shift_payload_hist(
            pion_shift_payload,
            "original_template_hists",
            KAON_SIGMA0_TEMPLATE_NAME,
            "H_pion_shift_original_k_sigma0_{}".format(context or analysis_scope),
        ),
        "H_pion_shift_selected_k_sigma0": _clone_shift_payload_hist(
            pion_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            KAON_SIGMA0_TEMPLATE_NAME,
            "H_pion_shift_selected_k_sigma0_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_original_pi_n": _clone_shift_payload_hist(
            kaon_shift_payload,
            "original_template_hists",
            "pi_n",
            "H_kaon_shift_original_pi_n_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_selected_pi_n": _clone_shift_payload_hist(
            kaon_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_n",
            "H_kaon_shift_selected_pi_n_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_original_pi_delta": _clone_shift_payload_hist(
            kaon_shift_payload,
            "original_template_hists",
            "pi_delta",
            "H_kaon_shift_original_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_selected_pi_delta": _clone_shift_payload_hist(
            kaon_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_delta",
            "H_kaon_shift_selected_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_original_pi_sidis": _clone_shift_payload_hist(
            kaon_shift_payload,
            "original_template_hists",
            "pi_sidis",
            "H_kaon_shift_original_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_selected_pi_sidis": _clone_shift_payload_hist(
            kaon_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            "pi_sidis",
            "H_kaon_shift_selected_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_original_k_lambda": _clone_shift_payload_hist(
            kaon_shift_payload,
            "original_template_hists",
            KAON_SIGNAL_TEMPLATE_NAME,
            "H_kaon_shift_original_k_lambda_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_selected_k_lambda": _clone_shift_payload_hist(
            kaon_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            KAON_SIGNAL_TEMPLATE_NAME,
            "H_kaon_shift_selected_k_lambda_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_original_k_sigma0": _clone_shift_payload_hist(
            kaon_shift_payload,
            "original_template_hists",
            KAON_SIGMA0_TEMPLATE_NAME,
            "H_kaon_shift_original_k_sigma0_{}".format(context or analysis_scope),
        ),
        "H_kaon_shift_selected_k_sigma0": _clone_shift_payload_hist(
            kaon_shift_payload,
            ("selected_template_hists", "selected_extra_template_hists"),
            KAON_SIGMA0_TEMPLATE_NAME,
            "H_kaon_shift_selected_k_sigma0_{}".format(context or analysis_scope),
        ),
        "H_pion_fit_pi_n_scaled": pion_fit["pi_n_scaled_hist"],
        "H_pion_fit_pi_delta_scaled": pion_fit["pi_delta_scaled_hist"],
        "H_pion_fit_pi_sidis_scaled": pion_fit["pi_sidis_scaled_hist"],
        "H_pion_fit_k_sigma0_scaled": pion_fit["k_sigma0_scaled_hist"],
        "H_pion_fit_total": pion_fit["fit_hist"],
        "H_pion_fit_pi_n_scaled_refined": pion_fit.get("refined_pi_n_scaled_hist"),
        "H_pion_fit_pi_delta_scaled_refined": pion_fit.get("refined_pi_delta_scaled_hist"),
        "H_pion_fit_pi_sidis_scaled_refined": pion_fit.get("refined_pi_sidis_scaled_hist"),
        "H_pion_fit_k_sigma0_scaled_refined": pion_fit.get("refined_k_sigma0_scaled_hist"),
        "H_pion_fit_total_refined": pion_fit.get("refined_fit_hist"),
        "H_pion_fit_step_overlays": pion_fit.get("step_overlays") or [],
        "H_kaon_fit_pi_n_scaled": kaon_fit["pi_n_scaled_hist"],
        "H_kaon_fit_pi_delta_scaled": kaon_fit["pi_delta_scaled_hist"],
        "H_kaon_fit_pi_sidis_scaled": kaon_fit["pi_sidis_scaled_hist"],
        "H_kaon_fit_k_lambda_scaled": kaon_fit["k_lambda_scaled_hist"],
        "H_kaon_fit_k_sigma0_scaled": kaon_fit["k_sigma0_scaled_hist"],
        "H_kaon_fit_k_lambda_reference": kaon_fit.get("k_lambda_reference_hist"),
        "H_kaon_fit_total": kaon_fit["fit_hist"],
        "H_kaon_full_fit_residual": kaon_fit.get("full_fit_residual_hist"),
        "H_kaon_pion_bg_fit_total": kaon_fit["pion_bg_fit_hist"],
        "H_kaon_fit_pi_n_scaled_refined": kaon_fit.get("refined_pi_n_scaled_hist"),
        "H_kaon_fit_pi_delta_scaled_refined": kaon_fit.get("refined_pi_delta_scaled_hist"),
        "H_kaon_fit_pi_sidis_scaled_refined": kaon_fit.get("refined_pi_sidis_scaled_hist"),
        "H_kaon_fit_k_lambda_scaled_refined": kaon_fit.get("refined_k_lambda_scaled_hist"),
        "H_kaon_fit_k_sigma0_scaled_refined": kaon_fit.get("refined_k_sigma0_scaled_hist"),
        "H_kaon_fit_k_sigma0_scaled_refined_pre_postrefine": kaon_fit.get("refined_k_sigma0_scaled_hist_pre_postrefine"),
        "H_kaon_fit_total_refined": kaon_fit.get("refined_fit_hist"),
        "H_kaon_full_fit_residual_refined": kaon_fit.get("refined_full_fit_residual_hist"),
        "H_kaon_fit_total_refined_pre_postrefine": kaon_fit.get("refined_fit_hist_pre_postrefine"),
        "H_kaon_pion_bg_fit_total_refined": kaon_fit.get("refined_pion_bg_fit_hist"),
        "H_kaon_pion_bg_fit_total_refined_pre_postrefine": kaon_fit.get("refined_pion_bg_fit_hist_pre_postrefine"),
        "H_kaon_fit_step_overlays": kaon_fit.get("step_overlays") or [],
        "H_fit_residual_pion": pion_fit["residual_hist"],
        "H_fit_residual_kaon": kaon_fit["residual_hist"],
        "H_fit_residual_pion_refined": pion_fit.get("refined_residual_hist"),
        "H_fit_residual_kaon_refined": kaon_fit.get("refined_residual_hist"),
        "H_pion_control_input": _clone_hist(
            h_pion_control,
            "H_pion_control_input_{}".format(context or analysis_scope),
        ),
        "H_kaon_nosub_input": _clone_hist(
            h_kaon_nosub,
            "H_kaon_nosub_input_{}".format(context or analysis_scope),
        ),
    }
    result["diagnostics"]["pion"]["template_mm_offset_data"] = template_mm_offset_data
    result["diagnostics"]["pion"]["template_mm_shift_applied"] = bool(abs(template_mm_offset_data) > 1e-12)
    result["diagnostics"]["kaon"]["template_mm_offset_data"] = template_mm_offset_data
    result["diagnostics"]["kaon"]["template_mm_shift_applied"] = bool(abs(template_mm_offset_data) > 1e-12)
    result["diagnostics"]["kaon"]["sigma0_requested"] = bool(kaon_fit.get("sigma0_requested"))
    result["diagnostics"]["kaon"]["sigma0_active"] = bool(kaon_fit.get("sigma0_active"))
    result["diagnostics"]["kaon"]["k_sigma0_template_provenance"] = {
        "source_template": "kaon_sigma0_shape",
        "simc_input_available": bool(_hist_has_usable_support(kaon_sigma0_shape)),
        "aligned_template_available": bool(_hist_has_usable_support(aligned_kaon_sigma0_shape)),
        "fitter_reference_available": bool(
            _hist_has_usable_support(k_sigma0_fitter_reference_hist)
        ),
        "published_source": (
            "inner_fitter_reference"
            if k_sigma0_fitter_reference_hist is not None
            else "aligned_kaon_sigma0_shape"
        ),
        "published_hist_key": "H_simc_shape_k_sigma0",
        "published_template_available": bool(
            _hist_has_usable_support(result.get("H_simc_shape_k_sigma0"))
        ),
        "protected_fit_template_available": bool(
            ((kaon_fit.get("diagnostics") or {}).get("pi_delta_signal_protected_fit") or {})
            .get("template_availability", {})
            .get(KAON_SIGMA0_TEMPLATE_NAME, False)
        ),
    }
    result["sigma0_active_kaon"] = bool(kaon_fit.get("sigma0_active"))
    if result["fallback_used"]:
        print(
            "WARNING: pion component fit fallback used\n"
            "  analysis_scope = {}\n"
            "  context = {}\n"
            "  fit_status_pion = {}\n"
            "  fit_status_kaon = {}\n"
            "  fallback_reason = {}".format(
                analysis_scope,
                context or analysis_scope,
                result["fit_status_pion"],
                result["fit_status_kaon"],
                result["fallback_reason"] or "unknown",
            )
        )
    return result


def _style_overlay_hist(hist, color, line_style=1, line_width=2):
    if hist is None:
        return
    hist.SetLineColor(color)
    hist.SetLineStyle(line_style)
    hist.SetLineWidth(line_width)
    hist.SetFillStyle(0)
    hist.SetMarkerStyle(0)


def _format_fit_number(value):
    if not _is_finite_number(value):
        return "n/a"
    return "{:.3e}".format(float(value))


def _format_fit_metric(value):
    if not _is_finite_number(value):
        return "n/a"
    return "{:.3g}".format(float(value))


def _component_plot_label(component_name):
    return (COMPONENT_PLOT_STYLE.get(component_name) or {}).get("label", str(component_name))


def _component_plot_color(component_name):
    return (COMPONENT_PLOT_STYLE.get(component_name) or {}).get("color", ROOT.kBlack)


def _component_hist_suffix(component_name):
    if component_name == KAON_SIGNAL_TEMPLATE_NAME:
        return "k_lambda"
    if component_name == KAON_SIGMA0_TEMPLATE_NAME:
        return "k_sigma0"
    return str(component_name)


def _format_fit_strategy(diagnostics):
    if not isinstance(diagnostics, dict):
        return "n/a"
    strategy = diagnostics.get("fit_strategy")
    if not strategy:
        return "n/a"
    n_passes = diagnostics.get("n_passes")
    if _is_finite_number(n_passes):
        return "{} x{}".format(strategy, int(float(n_passes)))
    return str(strategy)


def _format_solution_method(diagnostics):
    if not isinstance(diagnostics, dict):
        return "n/a"
    method = diagnostics.get("accepted_solution")
    return str(method) if method else "n/a"


def _format_validation_status(diagnostics):
    if not isinstance(diagnostics, dict):
        return "n/a"
    validation = diagnostics.get("validation") or {}
    return "pass" if bool(validation.get("accepted")) else "fail"


def _format_window_list(windows):
    if not windows:
        return "full-range"
    formatted = []
    for window_min, window_max in windows:
        formatted.append("[{:.3f}, {:.3f}]".format(float(window_min), float(window_max)))
    return ", ".join(formatted)


def _format_excluded_window_list(windows):
    return _format_window_list(windows) if windows else "none"


def _format_component_scale_map(scale_map):
    if not isinstance(scale_map, dict) or not scale_map:
        return "n/a"
    label_map = {
        "pi_n": "n",
        "pi_delta": "delta",
        "pi_sidis": "sidis",
        KAON_SIGMA0_TEMPLATE_NAME: "sigma0",
        KAON_SIGNAL_TEMPLATE_NAME: "lambda",
    }
    ordered_names = ["pi_n", "pi_delta", "pi_sidis", KAON_SIGMA0_TEMPLATE_NAME, KAON_SIGNAL_TEMPLATE_NAME]
    ordered_names.extend(
        component_name for component_name in scale_map.keys()
        if component_name not in ordered_names
    )
    formatted = []
    for component_name in ordered_names:
        if component_name not in scale_map:
            continue
        formatted.append(
            "{}={:.2f}".format(
                label_map.get(component_name, component_name),
                float(scale_map.get(component_name, 1.0) or 1.0),
            )
        )
    return ", ".join(formatted) if formatted else "n/a"


def _component_scale_map_has_nonunity(scale_map, tolerance=1e-12):
    if not isinstance(scale_map, dict):
        return False
    for scale_value in scale_map.values():
        try:
            if abs(float(scale_value) - 1.0) > float(tolerance):
                return True
        except Exception:
            continue
    return False


def _draw_window_collection(windows, y_min, y_max, color, line_style, line_width=2):
    drawn_lines = []
    for window_min, window_max in windows or []:
        drawn_lines.extend(
            _draw_vertical_window_lines(
                window_min,
                window_max,
                y_min,
                y_max,
                color=color,
                line_style=line_style,
                line_width=line_width,
            )
        )
    return drawn_lines


def _draw_window_band_collection(
    windows,
    y_min,
    y_max,
    color,
    alpha=0.10,
    fill_style=1001,
):
    drawn_bands = []
    for window_min, window_max in windows or []:
        band = ROOT.TBox(float(window_min), float(y_min), float(window_max), float(y_max))
        band.SetLineColor(color)
        band.SetLineStyle(1)
        band.SetLineWidth(1)
        if hasattr(band, "SetFillColorAlpha"):
            band.SetFillColorAlpha(color, float(alpha))
            band.SetFillStyle(fill_style)
        else:
            band.SetFillColor(color)
            band.SetFillStyle(3002)
        band.Draw("same")
        drawn_bands.append(band)
    return drawn_bands


def _draw_vertical_window_lines(
    window_min,
    window_max,
    y_min,
    y_max,
    color=None,
    line_style=None,
    line_width=2,
):
    lines = []
    line_color = ROOT.kBlue + 1 if color is None else color
    style_value = 3 if line_style is None else line_style
    for x_value in (window_min, window_max):
        line = ROOT.TLine(float(x_value), float(y_min), float(x_value), float(y_max))
        line.SetLineColor(line_color)
        line.SetLineStyle(style_value)
        line.SetLineWidth(line_width)
        line.Draw("same")
        lines.append(line)
    return lines


def _print_component_overlay_page(
    pdf_name,
    base_hist,
    base_label,
    title,
    overlay_specs,
    stats_lines,
    cut_window=None,
):
    if base_hist is None or not overlay_specs:
        return False

    canvas = ROOT.TCanvas()
    drawn_hists = []

    base_clone = _clone_hist(base_hist, "{}_plot".format(base_hist.GetName()))
    base_clone.SetTitle(title)
    base_clone.SetLineColor(ROOT.kBlack)
    base_clone.SetLineWidth(2)
    base_clone.SetFillStyle(3001)
    base_clone.SetFillColor(ROOT.kGray + 1)
    base_clone.SetMarkerStyle(20)
    base_clone.SetMarkerSize(0.7)
    drawn_hists.append(base_clone)

    y_max = max(base_clone.GetMaximum(), 0.0)
    for hist, _, _, _ in overlay_specs:
        if hist is None:
            continue
        y_max = max(y_max, hist.GetMaximum())
    if y_max <= 0.0:
        y_max = 1.0
    base_clone.SetMaximum(1.20 * y_max)
    base_clone.SetMinimum(0.0)
    base_clone.Draw("hist")

    legend = ROOT.TLegend(0.58, 0.56, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(base_clone, base_label, "lf")

    for hist, label, color, line_style in overlay_specs:
        if hist is None:
            continue
        hist_clone = _clone_hist(hist, "{}_plot".format(hist.GetName()))
        _style_overlay_hist(hist_clone, color, line_style=line_style)
        hist_clone.Draw("hist same")
        legend.AddEntry(hist_clone, label, "l")
        drawn_hists.append(hist_clone)

    if cut_window is not None:
        window_min, window_max = cut_window
        drawn_hists.extend(
            _draw_vertical_window_lines(window_min, window_max, 0.0, 1.20 * y_max)
        )

    legend.Draw()

    if stats_lines:
        stats_box = ROOT.TPaveText(0.14, 0.64, 0.52, 0.88, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.028)
        for line in stats_lines:
            stats_box.AddText(line)
        stats_box.Draw()
        drawn_hists.append(stats_box)

    canvas.Print(pdf_name)
    canvas.Close()
    return True


def _print_protected_overlay_page(
    pdf_name,
    base_hist,
    base_label,
    title,
    overlay_specs,
    diagnostics_lines,
    window=None,
):
    """Draw protected physics curves and diagnostics on separate pads."""
    if base_hist is None:
        return
    canvas_name = unique_root_object_name(
        "protected_canvas",
        scope="protected_render",
        role="canvas",
    )
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 1000, 850)
    canvas.Divide(1, 2)
    plot_pad = canvas.cd(1)
    plot_pad.SetPad(0.0, 0.30, 1.0, 1.0)
    # ROOT pads and legends keep raw C++ pointers.  Retain every temporary
    # PyROOT proxy until after Print() so a loop-local clone cannot be
    # collected while TLegend::PaintPrimitives is dereferencing it.
    drawn_objects = [canvas, plot_pad]
    base_clone = _clone_hist(base_hist, "{}_protected_plot".format(base_hist.GetName()))
    drawn_objects.append(base_clone)
    base_clone.SetTitle(title)
    base_clone.SetLineColor(ROOT.kBlack)
    base_clone.SetLineWidth(2)
    base_clone.SetFillStyle(3001)
    base_clone.SetFillColor(ROOT.kGray + 1)
    y_max = max(float(base_clone.GetMaximum()), 0.0)
    for hist, _, _, _ in overlay_specs:
        if hist is not None:
            y_max = max(y_max, float(hist.GetMaximum()))
    base_clone.SetMaximum(1.20 * y_max if y_max > 0.0 else 1.0)
    base_clone.SetMinimum(0.0)
    base_clone.Draw("hist")
    legend = ROOT.TLegend(0.58, 0.56, 0.88, 0.88)
    drawn_objects.append(legend)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(base_clone, base_label, "lf")
    for hist, label, color, line_style in overlay_specs:
        if hist is None:
            continue
        clone = _clone_hist(hist, "{}_protected_plot".format(hist.GetName()))
        drawn_objects.append(clone)
        _style_overlay_hist(clone, color, line_style=line_style)
        clone.Draw("hist same")
        legend.AddEntry(clone, label, "l")
    if window is not None:
        drawn_objects.extend(_draw_vertical_window_lines(
            float(window[0]), float(window[1]), 0.0, float(base_clone.GetMaximum())
        ))
    legend.Draw()
    diagnostics_pad = canvas.cd(2)
    drawn_objects.append(diagnostics_pad)
    diagnostics_pad.SetPad(0.0, 0.0, 1.0, 0.30)
    diagnostics_box = ROOT.TPaveText(0.03, 0.06, 0.97, 0.94, "NDC")
    drawn_objects.append(diagnostics_box)
    diagnostics_box.SetBorderSize(0)
    diagnostics_box.SetFillStyle(0)
    diagnostics_box.SetTextAlign(12)
    diagnostics_box.SetTextSize(0.055)
    for line in diagnostics_lines:
        diagnostics_box.AddText(str(line))
    diagnostics_box.Draw()
    canvas.Print(pdf_name)
    canvas.Close()


def _format_mm_range(range_values):
    if not isinstance(range_values, (list, tuple)) or len(range_values) != 2:
        return "none"
    if not (_is_finite_number(range_values[0]) and _is_finite_number(range_values[1])):
        return "none"
    return "[{:.3f}, {:.3f}]".format(float(range_values[0]), float(range_values[1]))


def _print_single_hist_page(
    pdf_name,
    hist,
    label,
    title,
    stats_lines,
    cut_window=None,
    line_color=None,
):
    if hist is None:
        return False

    canvas = ROOT.TCanvas()
    hist_clone = _clone_hist(hist, "{}_single_plot".format(hist.GetName()))
    hist_clone.SetTitle(title)
    hist_clone.SetLineColor(ROOT.kBlack if line_color is None else line_color)
    hist_clone.SetLineWidth(2)
    hist_clone.SetFillStyle(0)
    hist_clone.SetMarkerStyle(20)
    hist_clone.SetMarkerSize(0.6)
    y_max = max(hist_clone.GetMaximum(), 0.0)
    y_min = min(hist_clone.GetMinimum(), 0.0)
    if y_max <= 0.0:
        y_max = 1.0
    hist_clone.SetMaximum(1.20 * y_max)
    hist_clone.SetMinimum(1.20 * y_min if y_min < 0.0 else 0.0)
    hist_clone.Draw("hist")

    if cut_window is not None:
        _draw_vertical_window_lines(
            cut_window[0],
            cut_window[1],
            hist_clone.GetMinimum(),
            hist_clone.GetMaximum(),
        )

    legend = ROOT.TLegend(0.62, 0.78, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_clone, label, "l")
    legend.Draw()

    if stats_lines:
        stats_box = ROOT.TPaveText(0.14, 0.58, 0.56, 0.88, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.028)
        for line in stats_lines:
            stats_box.AddText(line)
        stats_box.Draw()

    canvas.Print(pdf_name)
    canvas.Close()
    return True


def _print_component_application_status_page(
    pdf_name,
    component_payload,
    title_prefix="",
):
    if not isinstance(component_payload, dict):
        return False

    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)

    canvas = ROOT.TCanvas()
    frame = ROOT.TH1F(
        "particle_subtraction_component_application_status_frame",
        "{}Part 3 pion reweighting status".format(title_prefix),
        1,
        0.0,
        1.0,
    )
    frame.SetStats(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.0)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.GetYaxis().SetLabelSize(0.0)
    frame.Draw()

    header = ROOT.TPaveText(0.12, 0.82, 0.88, 0.92, "NDC")
    header.SetBorderSize(0)
    header.SetFillStyle(0)
    header.SetTextAlign(12)
    header.SetTextSize(0.040)
    header.AddText(
        "status: {}".format(
            "accepted" if bool(component_payload.get("accepted")) else "rejected"
        )
    )
    header.Draw()

    details = ROOT.TPaveText(0.12, 0.18, 0.88, 0.78, "NDC")
    details.SetBorderSize(0)
    details.SetFillStyle(0)
    details.SetTextAlign(12)
    details.SetTextSize(0.028)
    details.AddText(
        "analysis scope: {}".format(
            component_payload.get("analysis_scope")
            or component_payload.get("analysis_scope_label")
            or "unknown"
        )
    )
    details.AddText(
        "particle subtraction mode: {}".format(
            component_payload.get("particle_subtraction_mode") or "unknown"
        )
    )
    details.AddText(
        "fallback used: {}".format(
            "yes" if bool(component_payload.get("fallback_used")) else "no"
        )
    )
    details.AddText(
        "fallback mode: {}".format(
            component_payload.get("fallback_mode") or "unknown"
        )
    )
    details.AddText(
        "fit status pion: {}".format(
            component_payload.get("fit_status_pion") or "unknown"
        )
    )
    details.AddText(
        "fit status kaon: {}".format(
            component_payload.get("fit_status_kaon") or "unknown"
        )
    )
    details.AddText(
        "fit validation pion: {}".format(
            "pass" if bool(component_payload.get("fit_validation_pion")) else "fail"
        )
    )
    details.AddText(
        "fit validation kaon: {}".format(
            "pass" if bool(component_payload.get("fit_validation_kaon")) else "fail"
        )
    )
    if component_payload.get("pion_parent_id"):
        parent_id = str(component_payload.get("pion_parent_id"))
        details.AddText("parent id: {}".format(parent_id[:28]))
        details.AddText(
            "child application: {} ({})".format(
                component_payload.get("child_application_status") or "unknown",
                component_payload.get("child_application_mode") or "unknown",
            )
        )
        support = component_payload.get("child_control_support") or {}
        details.AddText(
            "child controls: allcuts={} nommcuts={}".format(
                support.get("n_events_allcuts", "n/a"),
                support.get("n_events_nommcuts", "n/a"),
            )
        )

    fallback_reason = str(
        component_payload.get("fallback_reason")
        or "no fallback reason recorded"
    ).strip()
    details.AddText("fallback reason:")
    for line in fallback_reason.split("; "):
        if line:
            details.AddText("  {}".format(line))
    details.Draw()

    canvas.Print(pdf_name)
    canvas.Close()
    return True


def _build_difference_hist(data_hist, model_hist, name, divide_by_sigma=False):
    if data_hist is None or model_hist is None:
        return None
    diff_hist = _clone_hist(data_hist, name, reset=True)
    if diff_hist is None:
        return None
    for bin_index in range(1, data_hist.GetNbinsX() + 1):
        data_value = float(data_hist.GetBinContent(bin_index))
        model_value = float(model_hist.GetBinContent(bin_index))
        diff_value = data_value - model_value
        if divide_by_sigma:
            sigma_value = float(data_hist.GetBinError(bin_index))
            if (not math.isfinite(sigma_value)) or sigma_value <= 0.0:
                diff_value = 0.0
            else:
                diff_value = diff_value / sigma_value
        diff_hist.SetBinContent(bin_index, float(diff_value))
        diff_hist.SetBinError(bin_index, 0.0)
    return diff_hist


def _print_joint_refinement_overlay_page(
    pdf_name,
    data_hist,
    staged_total_hist,
    refined_total_hist,
    overlay_specs,
    title,
    stats_lines,
    cut_window=None,
):
    if data_hist is None or staged_total_hist is None or refined_total_hist is None:
        return

    canvas = ROOT.TCanvas("c_joint_refine_overlay_{}".format(data_hist.GetName()), "", 1000, 800)
    canvas.SetLeftMargin(0.10)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)

    data_clone = _clone_hist(data_hist, "{}_joint_data".format(data_hist.GetName()))
    staged_clone = _clone_hist(staged_total_hist, "{}_joint_stage".format(staged_total_hist.GetName()))
    refined_clone = _clone_hist(refined_total_hist, "{}_joint_refined".format(refined_total_hist.GetName()))
    data_clone.SetTitle(title)
    data_clone.SetLineColor(ROOT.kBlack)
    data_clone.SetLineWidth(2)
    data_clone.SetFillStyle(3001)
    data_clone.SetFillColor(ROOT.kGray + 1)
    data_clone.SetMarkerStyle(20)
    data_clone.SetMarkerSize(0.7)
    _style_overlay_hist(staged_clone, ROOT.kOrange + 7, line_style=2)
    _style_overlay_hist(refined_clone, ROOT.kGreen + 2, line_style=3)
    extra_clones = []
    y_max = max(data_clone.GetMaximum(), staged_clone.GetMaximum(), refined_clone.GetMaximum(), 0.0)
    for hist, _, _, _ in overlay_specs or []:
        if hist is None:
            continue
        y_max = max(y_max, hist.GetMaximum())
    if y_max <= 0.0:
        y_max = 1.0
    data_clone.SetMaximum(1.20 * y_max)
    data_clone.SetMinimum(0.0)
    data_clone.GetXaxis().SetLabelSize(0.036)
    data_clone.GetXaxis().SetTitleSize(0.042)
    data_clone.GetYaxis().SetTitleSize(0.045)
    data_clone.GetYaxis().SetLabelSize(0.036)
    data_clone.GetYaxis().SetTitleOffset(0.95)
    data_clone.Draw("hist")
    staged_clone.Draw("hist same")
    refined_clone.Draw("hist same")
    for hist, label, color, line_style in overlay_specs or []:
        if hist is None:
            continue
        hist_clone = _clone_hist(hist, "{}_joint_overlay".format(hist.GetName()))
        _style_overlay_hist(hist_clone, color, line_style=line_style)
        hist_clone.Draw("hist same")
        extra_clones.append((hist_clone, label))
    if cut_window is not None:
        _draw_vertical_window_lines(cut_window[0], cut_window[1], 0.0, 1.20 * y_max)

    legend = ROOT.TLegend(0.50, 0.62, 0.96, 0.89)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.024)
    legend.SetNColumns(2)
    legend.AddEntry(data_clone, "data", "lf")
    legend.AddEntry(staged_clone, "staged total", "l")
    legend.AddEntry(refined_clone, "refined total", "l")
    for hist_clone, label in extra_clones:
        legend.AddEntry(hist_clone, label, "l")
    legend.Draw()

    if stats_lines:
        stats_box = ROOT.TPaveText(0.12, 0.64, 0.44, 0.90, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.022)
        for line in stats_lines:
            stats_box.AddText(line)
        stats_box.Draw()
    canvas.Print(pdf_name)
    canvas.Close()

    stage_resid = _build_difference_hist(
        data_hist,
        staged_total_hist,
        "{}_joint_stage_resid".format(data_hist.GetName()),
        divide_by_sigma=False,
    )
    refined_resid = _build_difference_hist(
        data_hist,
        refined_total_hist,
        "{}_joint_refined_resid".format(data_hist.GetName()),
        divide_by_sigma=False,
    )
    canvas = ROOT.TCanvas("c_joint_refine_resid_{}".format(data_hist.GetName()), "", 1000, 700)
    canvas.SetLeftMargin(0.10)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    stage_resid.SetTitle("Residuals: data - model")
    stage_resid.SetLineColor(ROOT.kOrange + 7)
    stage_resid.SetLineWidth(2)
    stage_resid.SetFillStyle(0)
    stage_resid.SetMarkerStyle(20)
    stage_resid.SetMarkerColor(ROOT.kOrange + 7)
    stage_resid.SetMarkerSize(0.6)
    refined_resid.SetLineColor(ROOT.kGreen + 2)
    refined_resid.SetLineWidth(2)
    refined_resid.SetMarkerStyle(24)
    refined_resid.SetMarkerColor(ROOT.kGreen + 2)
    refined_resid.SetMarkerSize(0.6)
    resid_y_min = min(stage_resid.GetMinimum(), refined_resid.GetMinimum(), 0.0)
    resid_y_max = max(stage_resid.GetMaximum(), refined_resid.GetMaximum(), 0.0)
    resid_span = max(resid_y_max - resid_y_min, 1e-3)
    stage_resid.SetMaximum(resid_y_max + 0.20 * resid_span)
    stage_resid.SetMinimum(resid_y_min - 0.20 * resid_span)
    stage_resid.GetXaxis().SetLabelSize(0.036)
    stage_resid.GetXaxis().SetTitleSize(0.042)
    stage_resid.GetYaxis().SetTitleSize(0.045)
    stage_resid.GetYaxis().SetLabelSize(0.036)
    stage_resid.GetYaxis().SetTitleOffset(1.00)
    stage_resid.Draw("hist")
    refined_resid.Draw("hist same")
    zero_line_mid = ROOT.TLine(
        float(data_hist.GetXaxis().GetXmin()),
        0.0,
        float(data_hist.GetXaxis().GetXmax()),
        0.0,
    )
    zero_line_mid.SetLineColor(ROOT.kBlack)
    zero_line_mid.SetLineStyle(3)
    zero_line_mid.Draw("same")
    if cut_window is not None:
        _draw_vertical_window_lines(
            cut_window[0],
            cut_window[1],
            stage_resid.GetMinimum(),
            stage_resid.GetMaximum(),
        )
    resid_legend = ROOT.TLegend(0.68, 0.72, 0.95, 0.88)
    resid_legend.SetBorderSize(0)
    resid_legend.SetFillStyle(0)
    resid_legend.SetTextSize(0.028)
    resid_legend.AddEntry(stage_resid, "staged residual", "l")
    resid_legend.AddEntry(refined_resid, "refined residual", "l")
    resid_legend.Draw()
    canvas.Print(pdf_name)
    canvas.Close()

    stage_pull = _build_difference_hist(
        data_hist,
        staged_total_hist,
        "{}_joint_stage_pull".format(data_hist.GetName()),
        divide_by_sigma=True,
    )
    refined_pull = _build_difference_hist(
        data_hist,
        refined_total_hist,
        "{}_joint_refined_pull".format(data_hist.GetName()),
        divide_by_sigma=True,
    )
    canvas = ROOT.TCanvas("c_joint_refine_pull_{}".format(data_hist.GetName()), "", 1000, 700)
    canvas.SetLeftMargin(0.10)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    stage_pull.SetTitle("Pulls: (data - model) / sigma")
    stage_pull.SetLineColor(ROOT.kOrange + 7)
    stage_pull.SetLineWidth(2)
    stage_pull.SetMarkerStyle(20)
    stage_pull.SetMarkerColor(ROOT.kOrange + 7)
    stage_pull.SetMarkerSize(0.6)
    refined_pull.SetLineColor(ROOT.kGreen + 2)
    refined_pull.SetLineWidth(2)
    refined_pull.SetMarkerStyle(24)
    refined_pull.SetMarkerColor(ROOT.kGreen + 2)
    refined_pull.SetMarkerSize(0.6)
    pull_y_min = min(stage_pull.GetMinimum(), refined_pull.GetMinimum(), -1.0)
    pull_y_max = max(stage_pull.GetMaximum(), refined_pull.GetMaximum(), 1.0)
    pull_span = max(pull_y_max - pull_y_min, 1.0)
    stage_pull.SetMaximum(pull_y_max + 0.20 * pull_span)
    stage_pull.SetMinimum(pull_y_min - 0.20 * pull_span)
    stage_pull.GetXaxis().SetLabelSize(0.036)
    stage_pull.GetXaxis().SetTitleSize(0.042)
    stage_pull.GetYaxis().SetTitleSize(0.045)
    stage_pull.GetYaxis().SetLabelSize(0.036)
    stage_pull.GetYaxis().SetTitleOffset(0.95)
    stage_pull.Draw("hist")
    refined_pull.Draw("hist same")
    zero_line_bot = ROOT.TLine(
        float(data_hist.GetXaxis().GetXmin()),
        0.0,
        float(data_hist.GetXaxis().GetXmax()),
        0.0,
    )
    zero_line_bot.SetLineColor(ROOT.kBlack)
    zero_line_bot.SetLineStyle(3)
    zero_line_bot.Draw("same")
    if cut_window is not None:
        _draw_vertical_window_lines(
            cut_window[0],
            cut_window[1],
            stage_pull.GetMinimum(),
            stage_pull.GetMaximum(),
        )
    pull_legend = ROOT.TLegend(0.68, 0.66, 0.95, 0.88)
    pull_legend.SetBorderSize(0)
    pull_legend.SetFillStyle(0)
    pull_legend.SetTextSize(0.028)
    pull_legend.AddEntry(stage_pull, "staged pull", "l")
    pull_legend.AddEntry(refined_pull, "refined pull", "l")
    pull_legend.Draw()
    canvas.Print(pdf_name)
    canvas.Close()


def _print_kaon_pion_bg_comparison_page(
    pdf_name,
    data_hist,
    staged_pion_bg_hist,
    refined_pion_bg_hist,
    final_pion_bg_hist,
    sigma0_hist,
    title,
    stats_lines,
    cut_window=None,
):
    if data_hist is None or staged_pion_bg_hist is None or final_pion_bg_hist is None:
        return

    canvas = ROOT.TCanvas(
        "c_kaon_pion_bg_compare_{}".format(data_hist.GetName()),
        "",
        1000,
        900,
    )
    canvas.Divide(1, 2)

    top_pad = canvas.cd(1)
    top_pad.SetBottomMargin(0.10)
    top_pad.SetLeftMargin(0.10)
    top_pad.SetRightMargin(0.04)
    top_pad.SetTopMargin(0.08)

    data_clone = _clone_hist(data_hist, "{}_kaon_compare_data".format(data_hist.GetName()))
    staged_clone = _clone_hist(staged_pion_bg_hist, "{}_kaon_compare_stage".format(staged_pion_bg_hist.GetName()))
    refined_clone = _clone_hist(
        refined_pion_bg_hist,
        "{}_kaon_compare_refined".format(refined_pion_bg_hist.GetName()),
    ) if refined_pion_bg_hist is not None else None
    final_clone = _clone_hist(final_pion_bg_hist, "{}_kaon_compare_final".format(final_pion_bg_hist.GetName()))
    sigma0_clone = _clone_hist(
        sigma0_hist,
        "{}_kaon_compare_sigma0".format(sigma0_hist.GetName()),
    ) if sigma0_hist is not None else None

    data_clone.SetTitle(title)
    data_clone.SetLineColor(ROOT.kBlack)
    data_clone.SetLineWidth(2)
    data_clone.SetFillStyle(3001)
    data_clone.SetFillColor(ROOT.kGray + 1)
    data_clone.SetMarkerStyle(20)
    data_clone.SetMarkerSize(0.7)
    _style_overlay_hist(staged_clone, ROOT.kOrange + 7, line_style=2)
    if refined_clone is not None:
        _style_overlay_hist(refined_clone, ROOT.kBlue + 1, line_style=3)
    _style_overlay_hist(final_clone, ROOT.kGreen + 2, line_style=1)
    if sigma0_clone is not None:
        _style_overlay_hist(sigma0_clone, ROOT.kCyan + 2, line_style=1)

    y_max = max(
        data_clone.GetMaximum(),
        staged_clone.GetMaximum(),
        final_clone.GetMaximum(),
        refined_clone.GetMaximum() if refined_clone is not None else 0.0,
        sigma0_clone.GetMaximum() if sigma0_clone is not None else 0.0,
        0.0,
    )
    if y_max <= 0.0:
        y_max = 1.0
    data_clone.SetMaximum(1.20 * y_max)
    data_clone.SetMinimum(0.0)
    data_clone.Draw("hist")
    staged_clone.Draw("hist same")
    if refined_clone is not None:
        refined_clone.Draw("hist same")
    final_clone.Draw("hist same")
    if sigma0_clone is not None:
        sigma0_clone.Draw("hist same")
    if cut_window is not None:
        _draw_vertical_window_lines(cut_window[0], cut_window[1], 0.0, 1.20 * y_max)

    top_legend = ROOT.TLegend(0.54, 0.60, 0.96, 0.89)
    top_legend.SetBorderSize(0)
    top_legend.SetFillStyle(0)
    top_legend.SetTextSize(0.024)
    top_legend.AddEntry(data_clone, "kaon data", "lf")
    top_legend.AddEntry(staged_clone, "staged pion-bg sum", "l")
    if refined_clone is not None:
        top_legend.AddEntry(refined_clone, "refined pion-bg sum", "l")
    top_legend.AddEntry(final_clone, "final postrefine pion-bg sum", "l")
    if sigma0_clone is not None:
        top_legend.AddEntry(sigma0_clone, "K-Sigma0 contribution", "l")
    top_legend.Draw()

    if stats_lines:
        stats_box = ROOT.TPaveText(0.12, 0.60, 0.46, 0.89, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.022)
        for line in stats_lines:
            stats_box.AddText(str(line))
        stats_box.Draw()

    bottom_pad = canvas.cd(2)
    bottom_pad.SetTopMargin(0.08)
    bottom_pad.SetBottomMargin(0.12)
    bottom_pad.SetLeftMargin(0.10)
    bottom_pad.SetRightMargin(0.04)

    stage_resid = _build_difference_hist(
        data_hist,
        staged_pion_bg_hist,
        "{}_kaon_compare_stage_resid".format(data_hist.GetName()),
        divide_by_sigma=False,
    )
    refined_resid = _build_difference_hist(
        data_hist,
        refined_pion_bg_hist,
        "{}_kaon_compare_refined_resid".format(data_hist.GetName()),
        divide_by_sigma=False,
    ) if refined_pion_bg_hist is not None else None
    final_resid = _build_difference_hist(
        data_hist,
        final_pion_bg_hist,
        "{}_kaon_compare_final_resid".format(data_hist.GetName()),
        divide_by_sigma=False,
    )

    stage_resid.SetTitle("Residuals to kaon data: data - pion-bg model")
    stage_resid.SetLineColor(ROOT.kOrange + 7)
    stage_resid.SetLineWidth(2)
    if refined_resid is not None:
        refined_resid.SetLineColor(ROOT.kBlue + 1)
        refined_resid.SetLineWidth(2)
    final_resid.SetLineColor(ROOT.kGreen + 2)
    final_resid.SetLineWidth(2)

    resid_y_min = min(
        stage_resid.GetMinimum(),
        refined_resid.GetMinimum() if refined_resid is not None else 0.0,
        final_resid.GetMinimum(),
        0.0,
    )
    resid_y_max = max(
        stage_resid.GetMaximum(),
        refined_resid.GetMaximum() if refined_resid is not None else 0.0,
        final_resid.GetMaximum(),
        0.0,
    )
    resid_span = max(resid_y_max - resid_y_min, 1e-3)
    stage_resid.SetMaximum(resid_y_max + 0.20 * resid_span)
    stage_resid.SetMinimum(resid_y_min - 0.20 * resid_span)
    stage_resid.Draw("hist")
    if refined_resid is not None:
        refined_resid.Draw("hist same")
    final_resid.Draw("hist same")
    zero_line = ROOT.TLine(
        float(data_hist.GetXaxis().GetXmin()),
        0.0,
        float(data_hist.GetXaxis().GetXmax()),
        0.0,
    )
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineStyle(3)
    zero_line.Draw("same")
    if cut_window is not None:
        _draw_vertical_window_lines(
            cut_window[0],
            cut_window[1],
            stage_resid.GetMinimum(),
            stage_resid.GetMaximum(),
        )

    bottom_legend = ROOT.TLegend(0.58, 0.70, 0.95, 0.88)
    bottom_legend.SetBorderSize(0)
    bottom_legend.SetFillStyle(0)
    bottom_legend.SetTextSize(0.024)
    bottom_legend.AddEntry(stage_resid, "data - staged pion-bg", "l")
    if refined_resid is not None:
        bottom_legend.AddEntry(refined_resid, "data - refined pion-bg", "l")
    bottom_legend.AddEntry(final_resid, "data - final postrefine pion-bg", "l")
    bottom_legend.Draw()

    canvas.Print(pdf_name)
    canvas.Close()


def _print_component_text_page(pdf_name, title, header_lines, body_lines):
    canvas = ROOT.TCanvas("c_component_text_{}".format(abs(hash(title))), "", 900, 900)
    frame = ROOT.TH1F("component_text_frame_{}".format(abs(hash(title))), title, 1, 0.0, 1.0)
    frame.SetStats(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.0)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.GetYaxis().SetLabelSize(0.0)
    frame.Draw()

    header = ROOT.TPaveText(0.10, 0.82, 0.90, 0.92, "NDC")
    header.SetBorderSize(0)
    header.SetFillStyle(0)
    header.SetTextAlign(12)
    header.SetTextSize(0.032)
    for line in header_lines or []:
        header.AddText(str(line))
    header.Draw()

    body = ROOT.TPaveText(0.10, 0.10, 0.90, 0.80, "NDC")
    body.SetBorderSize(0)
    body.SetFillStyle(0)
    body.SetTextAlign(12)
    body.SetTextSize(0.024)
    for line in body_lines or []:
        body.AddText(str(line))
    body.Draw()

    canvas.Print(pdf_name)
    canvas.Close()


def _format_matrix_text_lines(matrix_map, ordered_names):
    lines = []
    if not isinstance(matrix_map, dict) or not matrix_map:
        return ["n/a"]
    for left_name in ordered_names:
        row = matrix_map.get(left_name) or {}
        row_values = []
        for right_name in ordered_names:
            value = row.get(right_name)
            row_values.append("{}={}".format(right_name, _format_fit_metric(value)))
        lines.append("{}: {}".format(left_name, ", ".join(row_values)))
    return lines


def _print_component_step_pages(
    pdf_name,
    target_hist,
    step_overlays,
    title_prefix,
    sample_label,
):
    if target_hist is None or not step_overlays:
        return

    for step_overlay in step_overlays:
        baseline_before = step_overlay.get("H_baseline_before")
        residual_input = step_overlay.get("H_residual_input")
        component_template = step_overlay.get("H_component_template")
        component_scaled = step_overlay.get("H_component_scaled")
        cumulative_after = step_overlay.get("H_cumulative_after")
        component_name = step_overlay.get("component_name")
        component_label = step_overlay.get("component_label") or _component_plot_label(component_name)
        component_color = _component_plot_color(component_name)
        if baseline_before is None or residual_input is None or component_scaled is None or cumulative_after is None:
            continue

        canvas = ROOT.TCanvas("c_step_{}".format(step_overlay.get("step_index", 0)), "", 900, 900)
        canvas.Divide(1, 2)

        top_pad = canvas.cd(1)
        top_pad.SetBottomMargin(0.12)
        target_clone = _clone_hist(target_hist, "{}_step_target".format(target_hist.GetName()))
        baseline_clone = _clone_hist(baseline_before, "{}_step_baseline".format(baseline_before.GetName()))
        template_clone = _clone_hist(
            component_template,
            "{}_step_template".format(component_template.GetName()),
        ) if component_template is not None else None
        component_clone = _clone_hist(component_scaled, "{}_step_component".format(component_scaled.GetName()))
        cumulative_clone = _clone_hist(cumulative_after, "{}_step_cumulative".format(cumulative_after.GetName()))
        excluded_windows = step_overlay.get("excluded_windows") or []
        if excluded_windows:
            _mask_hist_windows_inplace(baseline_clone, excluded_windows)
            _mask_hist_windows_inplace(template_clone, excluded_windows)
            _mask_hist_windows_inplace(component_clone, excluded_windows)
            _mask_hist_windows_inplace(cumulative_clone, excluded_windows)
        target_clone.SetTitle(
            "{}{} step {}: {}".format(
                title_prefix,
                sample_label,
                step_overlay.get("step_index", 0),
                component_label,
            )
        )
        target_clone.SetLineColor(ROOT.kBlack)
        target_clone.SetLineWidth(2)
        target_clone.SetFillStyle(3001)
        target_clone.SetFillColor(ROOT.kGray + 1)
        target_clone.SetMarkerStyle(20)
        target_clone.SetMarkerSize(0.7)
        _style_overlay_hist(baseline_clone, ROOT.kOrange + 7, line_style=2)
        _style_overlay_hist(template_clone, ROOT.kBlue + 1, line_style=2)
        _style_overlay_hist(component_clone, component_color, line_style=1)
        _style_overlay_hist(cumulative_clone, ROOT.kGreen + 2, line_style=3)
        top_y_max = max(
            target_clone.GetMaximum(),
            baseline_clone.GetMaximum(),
            template_clone.GetMaximum() if template_clone is not None else 0.0,
            component_clone.GetMaximum(),
            cumulative_clone.GetMaximum(),
            0.0,
        )
        if top_y_max <= 0.0:
            top_y_max = 1.0
        target_clone.SetMaximum(1.20 * top_y_max)
        target_clone.SetMinimum(0.0)
        target_clone.Draw("hist")
        top_bands = _draw_window_band_collection(
            step_overlay.get("amplitude_windows") or [],
            0.0,
            1.20 * top_y_max,
            ROOT.kMagenta + 2,
            alpha=0.10,
        )
        baseline_clone.Draw("hist same")
        if template_clone is not None:
            template_clone.Draw("hist same")
        component_clone.Draw("hist same")
        cumulative_clone.Draw("hist same")
        top_anchor_lines = _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            0.0,
            1.20 * top_y_max,
            ROOT.kBlue + 1,
            3,
        )
        top_core_lines = _draw_window_collection(
            step_overlay.get("amplitude_windows") or [],
            0.0,
            1.20 * top_y_max,
            ROOT.kMagenta + 2,
            7,
        )
        top_exclude_lines = _draw_window_collection(
            step_overlay.get("excluded_windows") or [],
            0.0,
            1.20 * top_y_max,
            ROOT.kGray + 2,
            2,
        )

        top_legend = ROOT.TLegend(0.58, 0.58, 0.88, 0.88)
        top_legend.SetBorderSize(0)
        top_legend.SetFillStyle(0)
        top_legend.AddEntry(target_clone, "{} data".format(sample_label), "lf")
        top_legend.AddEntry(baseline_clone, "baseline before step", "l")
        if template_clone is not None:
            top_legend.AddEntry(template_clone, "raw SIMC template", "l")
        top_legend.AddEntry(component_clone, "{} contribution".format(component_label), "l")
        top_legend.AddEntry(cumulative_clone, "baseline after step", "l")
        top_legend.Draw()

        stats_box = ROOT.TPaveText(0.14, 0.60, 0.52, 0.88, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.028)
        fit_diagnostics = step_overlay.get("fit_diagnostics") or {}
        stats_box.AddText("pass: {}".format(step_overlay.get("pass_index", 0)))
        stats_box.AddText("component: {}".format(component_label))
        stats_box.AddText("amplitude: {}".format(_format_fit_number(step_overlay.get("amplitude"))))
        if step_overlay.get("sigma") is not None:
            stats_box.AddText("sigma: {}".format(_format_fit_number(step_overlay.get("sigma"))))
        if step_overlay.get("chi2") is not None:
            stats_box.AddText("chi2: {}".format(_format_fit_metric(step_overlay.get("chi2"))))
        stats_box.AddText("fit bins: {}".format(int(step_overlay.get("n_fit_bins", 0) or 0)))
        if str(step_overlay.get("amplitude_mode") or "least_squares") != "least_squares":
            stats_box.AddText("mode: {}".format(str(step_overlay.get("amplitude_mode"))))
        if step_overlay.get("amplitude_windows"):
            stats_box.AddText(
                "core: {}".format(_format_window_list(step_overlay.get("amplitude_windows") or []))
            )
        if fit_diagnostics.get("estimator") == "window_integral":
            stats_box.AddText(
                "core target int: {}".format(_format_fit_number(fit_diagnostics.get("target_sum")))
            )
            stats_box.AddText(
                "core template int: {}".format(_format_fit_number(fit_diagnostics.get("template_sum")))
            )
        elif fit_diagnostics.get("estimator") == "least_squares":
            stats_box.AddText(
                "LS num: {}".format(_format_fit_number(fit_diagnostics.get("weighted_numerator")))
            )
            stats_box.AddText(
                "LS den: {}".format(_format_fit_number(fit_diagnostics.get("weighted_denominator")))
            )
        if step_overlay.get("postfit_scale_factor") is not None:
            stats_box.AddText(
                "post-fit scale: {}".format(
                    _format_fit_number(step_overlay.get("postfit_scale_factor"))
                )
            )
        stats_box.AddText("anchor: {}".format(_format_window_list(step_overlay.get("anchor_windows") or [])))
        if excluded_windows:
            stats_box.AddText("exclude: {}".format(_format_window_list(excluded_windows)))
        stats_box.Draw()

        bottom_pad = canvas.cd(2)
        bottom_pad.SetTopMargin(0.08)
        bottom_pad.SetBottomMargin(0.12)
        residual_clone = _clone_hist(residual_input, "{}_step_residual".format(residual_input.GetName()))
        template_bottom_clone = _clone_hist(
            component_template,
            "{}_step_template_bottom".format(component_template.GetName()),
        ) if component_template is not None else None
        component_bottom_clone = _clone_hist(
            component_scaled,
            "{}_step_component_bottom".format(component_scaled.GetName()),
        )
        if excluded_windows:
            _mask_hist_windows_inplace(residual_clone, excluded_windows)
            _mask_hist_windows_inplace(template_bottom_clone, excluded_windows)
            _mask_hist_windows_inplace(component_bottom_clone, excluded_windows)
        residual_clone.SetTitle("Residual input to {} step".format(component_label))
        residual_clone.SetLineColor(ROOT.kBlack)
        residual_clone.SetLineWidth(2)
        residual_clone.SetFillStyle(3001)
        residual_clone.SetFillColor(ROOT.kGray + 1)
        residual_clone.SetMarkerStyle(20)
        residual_clone.SetMarkerSize(0.7)
        _style_overlay_hist(template_bottom_clone, ROOT.kBlue + 1, line_style=2)
        _style_overlay_hist(component_bottom_clone, component_color, line_style=1)
        bottom_y_max = max(
            residual_clone.GetMaximum(),
            template_bottom_clone.GetMaximum() if template_bottom_clone is not None else 0.0,
            component_bottom_clone.GetMaximum(),
            0.0,
        )
        if bottom_y_max <= 0.0:
            bottom_y_max = 1.0
        bottom_y_min = min(
            0.0,
            residual_clone.GetMinimum(),
            template_bottom_clone.GetMinimum() if template_bottom_clone is not None else 0.0,
            component_bottom_clone.GetMinimum(),
        )
        residual_clone.SetMaximum(1.20 * bottom_y_max)
        residual_clone.SetMinimum(1.20 * bottom_y_min if bottom_y_min < 0.0 else 0.0)
        residual_clone.Draw("hist")
        bottom_bands = _draw_window_band_collection(
            step_overlay.get("amplitude_windows") or [],
            residual_clone.GetMinimum(),
            residual_clone.GetMaximum(),
            ROOT.kMagenta + 2,
            alpha=0.10,
        )
        if template_bottom_clone is not None:
            template_bottom_clone.Draw("hist same")
        component_bottom_clone.Draw("hist same")
        bottom_anchor_lines = _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            residual_clone.GetMinimum(),
            residual_clone.GetMaximum(),
            ROOT.kBlue + 1,
            3,
        )
        bottom_core_lines = _draw_window_collection(
            step_overlay.get("amplitude_windows") or [],
            residual_clone.GetMinimum(),
            residual_clone.GetMaximum(),
            ROOT.kMagenta + 2,
            7,
        )
        bottom_exclude_lines = _draw_window_collection(
            step_overlay.get("excluded_windows") or [],
            residual_clone.GetMinimum(),
            residual_clone.GetMaximum(),
            ROOT.kGray + 2,
            2,
        )

        bottom_legend = ROOT.TLegend(0.58, 0.70, 0.88, 0.88)
        bottom_legend.SetBorderSize(0)
        bottom_legend.SetFillStyle(0)
        bottom_legend.AddEntry(residual_clone, "residual before step", "lf")
        if template_bottom_clone is not None:
            bottom_legend.AddEntry(template_bottom_clone, "raw SIMC template", "l")
        bottom_legend.AddEntry(component_bottom_clone, "{} fit".format(component_label), "l")
        bottom_legend.Draw()

        canvas.Print(pdf_name)
        canvas.Close()


def _format_component_template_diag_line(component_label, diagnostics):
    diagnostics = diagnostics or {}
    summary = (
        "{}: seen={} passed={} mm_passed={} norm={} fallback={}".format(
            component_label,
            int(diagnostics.get("n_events_seen", 0) or 0),
            int(diagnostics.get("n_events_passed", 0) or 0),
            int(diagnostics.get("n_events_passed_mm_window", 0) or 0),
            bool(diagnostics.get("normalized", False)),
            bool(diagnostics.get("fallback_used", False)),
        )
    )
    if component_label != "K-Sigma0":
        return summary

    identity = diagnostics.get("source_identity") or {}
    return "{} source={} Q{} W{} {} {} tree={} entries={} missing_branches={} reason={}".format(
        summary,
        diagnostics.get("resolution_source", "unknown"),
        identity.get("Q2", ""),
        identity.get("W", ""),
        identity.get("EPSSET", ""),
        identity.get("phi_setting", ""),
        diagnostics.get("tree_name", ""),
        diagnostics.get("tree_entries", ""),
        len(diagnostics.get("missing_required_branches") or []),
        diagnostics.get("fallback_reason") or "none",
    )


def _print_component_amplitude_pages(
    pdf_name,
    target_hist,
    step_overlays,
    title_prefix,
    sample_label,
):
    if target_hist is None or not step_overlays:
        return

    for step_overlay in step_overlays:
        fit_diagnostics = step_overlay.get("fit_diagnostics") or {}
        fit_bin_indices = [int(value) for value in (fit_diagnostics.get("fit_bin_indices") or [])]
        fit_y = [float(value) for value in (fit_diagnostics.get("fit_y") or [])]
        fit_sigma = [float(value) for value in (fit_diagnostics.get("fit_sigma") or [])]
        fit_template = [float(value) for value in (fit_diagnostics.get("fit_template") or [])]
        if not fit_bin_indices or not fit_y or not fit_template:
            continue
        if not (len(fit_bin_indices) == len(fit_y) == len(fit_template)):
            continue
        if len(fit_sigma) != len(fit_y):
            fit_sigma = [0.0] * len(fit_y)

        component_name = step_overlay.get("component_name")
        component_label = step_overlay.get("component_label") or _component_plot_label(component_name)
        component_color = _component_plot_color(component_name)
        amplitude = float(step_overlay.get("amplitude", 0.0) or 0.0)

        residual_fit_hist = _clone_hist(target_hist, "{}_amp_resid".format(target_hist.GetName()), reset=True)
        template_fit_hist = _clone_hist(target_hist, "{}_amp_template".format(target_hist.GetName()), reset=True)
        fitted_template_hist = _clone_hist(target_hist, "{}_amp_fitted_template".format(target_hist.GetName()), reset=True)
        ratio_raw_hist = _clone_hist(target_hist, "{}_amp_ratio_raw".format(target_hist.GetName()), reset=True)
        ratio_scaled_hist = _clone_hist(target_hist, "{}_amp_ratio_scaled".format(target_hist.GetName()), reset=True)
        ratio_raw_values = []
        ratio_scaled_values = []

        for idx, bin_index in enumerate(fit_bin_indices):
            y_value = float(fit_y[idx])
            sigma_value = float(fit_sigma[idx])
            template_value = float(fit_template[idx])
            fitted_value = amplitude * template_value

            residual_fit_hist.SetBinContent(bin_index, y_value)
            residual_fit_hist.SetBinError(bin_index, sigma_value)
            template_fit_hist.SetBinContent(bin_index, template_value)
            fitted_template_hist.SetBinContent(bin_index, fitted_value)

            if abs(template_value) > 1e-12:
                raw_ratio = y_value / template_value
                raw_ratio_err = abs(sigma_value / template_value) if sigma_value > 0.0 else 0.0
                ratio_raw_hist.SetBinContent(bin_index, raw_ratio)
                ratio_raw_hist.SetBinError(bin_index, raw_ratio_err)
                if math.isfinite(raw_ratio):
                    ratio_raw_values.append(raw_ratio)
            if abs(fitted_value) > 1e-12:
                scaled_ratio = y_value / fitted_value
                scaled_ratio_err = abs(sigma_value / fitted_value) if sigma_value > 0.0 else 0.0
                ratio_scaled_hist.SetBinContent(bin_index, scaled_ratio)
                ratio_scaled_hist.SetBinError(bin_index, scaled_ratio_err)
                if math.isfinite(scaled_ratio):
                    ratio_scaled_values.append(scaled_ratio)

        canvas = ROOT.TCanvas(
            "c_amp_step_{}".format(step_overlay.get("step_index", 0)),
            "",
            900,
            900,
        )
        canvas.Divide(1, 2)

        top_pad = canvas.cd(1)
        top_pad.SetBottomMargin(0.12)
        residual_fit_hist.SetTitle(
            "{}{} amplitude diagnostic step {}: {}".format(
                title_prefix,
                sample_label,
                step_overlay.get("step_index", 0),
                component_label,
            )
        )
        residual_fit_hist.SetLineColor(ROOT.kBlack)
        residual_fit_hist.SetLineWidth(2)
        residual_fit_hist.SetFillStyle(3001)
        residual_fit_hist.SetFillColor(ROOT.kGray + 1)
        residual_fit_hist.SetMarkerStyle(20)
        residual_fit_hist.SetMarkerSize(0.7)
        _style_overlay_hist(template_fit_hist, ROOT.kBlue + 1, line_style=2)
        _style_overlay_hist(fitted_template_hist, component_color, line_style=1)
        top_y_max = max(
            residual_fit_hist.GetMaximum(),
            template_fit_hist.GetMaximum(),
            fitted_template_hist.GetMaximum(),
            0.0,
        )
        if top_y_max <= 0.0:
            top_y_max = 1.0
        top_y_min = min(
            0.0,
            residual_fit_hist.GetMinimum(),
            template_fit_hist.GetMinimum(),
            fitted_template_hist.GetMinimum(),
        )
        residual_fit_hist.SetMaximum(1.20 * top_y_max)
        residual_fit_hist.SetMinimum(1.20 * top_y_min if top_y_min < 0.0 else 0.0)
        residual_fit_hist.Draw("E1")
        template_fit_hist.Draw("hist same")
        fitted_template_hist.Draw("hist same")
        _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            residual_fit_hist.GetMinimum(),
            residual_fit_hist.GetMaximum(),
            ROOT.kBlue + 1,
            3,
        )
        _draw_window_collection(
            step_overlay.get("amplitude_windows") or [],
            residual_fit_hist.GetMinimum(),
            residual_fit_hist.GetMaximum(),
            ROOT.kMagenta + 2,
            7,
        )
        _draw_window_collection(
            step_overlay.get("excluded_windows") or [],
            residual_fit_hist.GetMinimum(),
            residual_fit_hist.GetMaximum(),
            ROOT.kGray + 2,
            2,
        )

        top_legend = ROOT.TLegend(0.56, 0.60, 0.88, 0.88)
        top_legend.SetBorderSize(0)
        top_legend.SetFillStyle(0)
        top_legend.AddEntry(residual_fit_hist, "fit-bin residual input", "lep")
        top_legend.AddEntry(template_fit_hist, "raw template", "l")
        top_legend.AddEntry(fitted_template_hist, "fitted template", "l")
        top_legend.Draw()

        stats_box = ROOT.TPaveText(0.14, 0.56, 0.52, 0.88, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.026)
        stats_box.AddText("pass: {}".format(step_overlay.get("pass_index", 0)))
        stats_box.AddText("component: {}".format(component_label))
        stats_box.AddText("amplitude: {}".format(_format_fit_number(amplitude)))
        if step_overlay.get("sigma") is not None:
            stats_box.AddText("sigma: {}".format(_format_fit_number(step_overlay.get("sigma"))))
        if step_overlay.get("chi2") is not None:
            stats_box.AddText("chi2: {}".format(_format_fit_metric(step_overlay.get("chi2"))))
        stats_box.AddText("fit bins: {}".format(int(step_overlay.get("n_fit_bins", 0) or 0)))
        stats_box.AddText("mode: {}".format(str(step_overlay.get("amplitude_mode") or "least_squares")))
        if step_overlay.get("amplitude_windows"):
            stats_box.AddText("core: {}".format(_format_window_list(step_overlay.get("amplitude_windows") or [])))
        if fit_diagnostics.get("estimator") == "window_integral":
            stats_box.AddText("target int: {}".format(_format_fit_number(fit_diagnostics.get("target_sum"))))
            stats_box.AddText("template int: {}".format(_format_fit_number(fit_diagnostics.get("template_sum"))))
        elif fit_diagnostics.get("estimator") == "least_squares":
            stats_box.AddText("LS num: {}".format(_format_fit_number(fit_diagnostics.get("weighted_numerator"))))
            stats_box.AddText("LS den: {}".format(_format_fit_number(fit_diagnostics.get("weighted_denominator"))))
        stats_box.Draw()

        bottom_pad = canvas.cd(2)
        bottom_pad.SetTopMargin(0.08)
        bottom_pad.SetBottomMargin(0.12)
        ratio_raw_hist.SetTitle("Residual/template ratios by fit bin")
        ratio_raw_hist.SetLineColor(component_color)
        ratio_raw_hist.SetLineWidth(2)
        ratio_raw_hist.SetMarkerColor(component_color)
        ratio_raw_hist.SetMarkerStyle(20)
        ratio_raw_hist.SetMarkerSize(0.7)
        ratio_scaled_hist.SetLineColor(ROOT.kGreen + 2)
        ratio_scaled_hist.SetLineWidth(2)
        ratio_scaled_hist.SetMarkerColor(ROOT.kGreen + 2)
        ratio_scaled_hist.SetMarkerStyle(24)
        ratio_scaled_hist.SetMarkerSize(0.6)
        ratio_values = [
            float(value)
            for value in (ratio_raw_values + ratio_scaled_values + [amplitude, 1.0])
            if math.isfinite(float(value))
        ]
        ratio_y_min = min(ratio_values) if ratio_values else 0.0
        ratio_y_max = max(ratio_values) if ratio_values else 1.0
        ratio_span = max(ratio_y_max - ratio_y_min, 0.25)
        ratio_raw_hist.SetMaximum(ratio_y_max + 0.20 * ratio_span)
        ratio_raw_hist.SetMinimum(ratio_y_min - 0.20 * ratio_span)
        ratio_raw_hist.Draw("E1")
        ratio_scaled_hist.Draw("E1 same")
        _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            ratio_raw_hist.GetMinimum(),
            ratio_raw_hist.GetMaximum(),
            ROOT.kBlue + 1,
            3,
        )
        _draw_window_collection(
            step_overlay.get("amplitude_windows") or [],
            ratio_raw_hist.GetMinimum(),
            ratio_raw_hist.GetMaximum(),
            ROOT.kMagenta + 2,
            7,
        )
        _draw_window_collection(
            step_overlay.get("excluded_windows") or [],
            ratio_raw_hist.GetMinimum(),
            ratio_raw_hist.GetMaximum(),
            ROOT.kGray + 2,
            2,
        )
        amplitude_line = ROOT.TLine(
            float(target_hist.GetXaxis().GetXmin()),
            amplitude,
            float(target_hist.GetXaxis().GetXmax()),
            amplitude,
        )
        amplitude_line.SetLineColor(component_color)
        amplitude_line.SetLineStyle(2)
        amplitude_line.SetLineWidth(2)
        amplitude_line.Draw("same")
        unity_line = ROOT.TLine(
            float(target_hist.GetXaxis().GetXmin()),
            1.0,
            float(target_hist.GetXaxis().GetXmax()),
            1.0,
        )
        unity_line.SetLineColor(ROOT.kGreen + 2)
        unity_line.SetLineStyle(3)
        unity_line.SetLineWidth(2)
        unity_line.Draw("same")

        bottom_legend = ROOT.TLegend(0.50, 0.64, 0.88, 0.88)
        bottom_legend.SetBorderSize(0)
        bottom_legend.SetFillStyle(0)
        bottom_legend.AddEntry(ratio_raw_hist, "residual / raw template", "lep")
        bottom_legend.AddEntry(amplitude_line, "fitted amplitude", "l")
        bottom_legend.AddEntry(ratio_scaled_hist, "residual / fitted template", "lep")
        bottom_legend.AddEntry(unity_line, "unity closure", "l")
        bottom_legend.Draw()

        closure_mean = (
            sum(ratio_scaled_values) / len(ratio_scaled_values)
            if ratio_scaled_values else float("nan")
        )
        closure_rms = (
            math.sqrt(
                sum((value - closure_mean) ** 2 for value in ratio_scaled_values) / len(ratio_scaled_values)
            )
            if ratio_scaled_values and math.isfinite(closure_mean) else float("nan")
        )
        ratio_stats = ROOT.TPaveText(0.14, 0.74, 0.44, 0.88, "NDC")
        ratio_stats.SetBorderSize(0)
        ratio_stats.SetFillStyle(0)
        ratio_stats.SetTextAlign(12)
        ratio_stats.SetTextSize(0.024)
        ratio_stats.AddText("mean(resid/[A*T])={}".format(_format_fit_metric(closure_mean)))
        ratio_stats.AddText("rms(resid/[A*T])={}".format(_format_fit_metric(closure_rms)))
        ratio_stats.Draw()

        canvas.Print(pdf_name)
        canvas.Close()


def _residual_shift_summary_has_content(summary):
    if not isinstance(summary, dict) or not summary:
        return False
    if bool(summary.get("enabled")):
        return True
    selected_shift_point = summary.get("selected_shift_point") or {}
    if any(abs(float(delta_mm or 0.0)) > 1e-12 for delta_mm in selected_shift_point.values()):
        return True
    return int(summary.get("candidate_count", 0) or 0) > 1


def _resolve_protected_pi_delta_render_state(component_payload):
    """Resolve the PDF contract for the kaon signal-protected pi-delta stage.

    The resolved kaon configuration is the authoritative statement of whether
    protected rendering is required.  The diagnostic is deliberately only a
    fallback for older in-memory payloads that predate the resolved-config
    handoff.  Keeping this decision separate from the fit itself prevents a
    failed or incomplete protected payload from visually reviving the retired
    staged kaon pi-delta presentation.
    """
    payload = component_payload if isinstance(component_payload, dict) else {}
    resolved_config = payload.get("resolved_subtraction_config") or {}
    kaon_config = resolved_config.get("kaon_nosub") or {}
    configured = kaon_config.get("pi_delta_signal_protected_fit")

    diagnostics = payload.get("diagnostics") or {}
    kaon_diagnostics = diagnostics.get("kaon") or {}
    protected_diagnostic = kaon_diagnostics.get("pi_delta_signal_protected_fit")
    if not isinstance(protected_diagnostic, dict):
        protected_diagnostic = {}

    if not isinstance(configured, dict):
        configured = protected_diagnostic.get("resolved_configuration") or {}
    if not isinstance(configured, dict):
        configured = {}

    configured_enabled = bool(configured.get("enabled"))
    if not configured and protected_diagnostic:
        configured_enabled = bool(protected_diagnostic.get("enabled"))

    fit_variant = str(protected_diagnostic.get("fit_variant") or "").strip()
    required_histogram_keys = (
        "H_pi_delta_protected_fit_input",
        "H_pi_delta_lambda_gauge",
        "H_pi_delta_protected_k_lambda",
        "H_pi_delta_protected_pi_delta",
        "H_pi_delta_protected_fit_total",
        "H_pi_delta_protected_after_subtraction",
    )
    if fit_variant != "lambda_only_protected_fallback":
        required_histogram_keys = (
            *required_histogram_keys[:2],
            "H_pi_delta_protected_k_sigma0",
            *required_histogram_keys[2:],
        )
    missing_histograms = [
        key for key in required_histogram_keys if payload.get(key) is None
    ]
    diagnostic_available = bool(protected_diagnostic)
    diagnostic_status = str(protected_diagnostic.get("status") or "").strip().lower()
    failure_policy = (
        protected_diagnostic.get("failure_policy")
        or configured.get("failure_policy")
        or "unknown"
    )

    if not configured_enabled:
        state = "disabled"
    elif not diagnostic_available:
        state = "missing_payload"
    elif diagnostic_status != "success":
        state = "fit_failure"
    elif missing_histograms:
        state = "missing_payload"
    else:
        state = "success"

    return {
        "state": state,
        "configured_enabled": configured_enabled,
        "diagnostic_available": diagnostic_available,
        "diagnostic_status": diagnostic_status or "missing",
        "fit_variant": fit_variant or "unknown",
        "diagnostic": protected_diagnostic,
        "resolved_config": configured,
        "failure_policy": str(failure_policy),
        "missing_histograms": missing_histograms,
        "render_success_pages": state == "success",
        # This is intentionally configuration-driven: once protected mode is
        # requested, deprecated staged kaon pages must never reappear.
        "suppress_deprecated_kaon_pages": configured_enabled,
    }


def _filtered_kaon_early_component_step_overlays(step_overlays):
    """Return display-only kaon early-step overlays without changing storage."""
    allowed_components = {"pi_n", "pi_sidis"}
    return [
        dict(step_overlay)
        for step_overlay in (step_overlays or [])
        if isinstance(step_overlay, dict)
        and step_overlay.get("component_name") in allowed_components
    ]


def _print_protected_pi_delta_status_page(
    pdf_name,
    component_payload,
    render_state,
    title_prefix="",
):
    """Render an explicit unavailable page instead of a staged-kaon fallback."""
    diagnostic = render_state.get("diagnostic") or {}
    lambda_gauge = diagnostic.get("lambda_gauge") or {}
    availability = diagnostic.get("template_availability") or {}
    sigma0_source_availability = diagnostic.get("sigma0_source_availability") or {}
    sigma0_scope_availability = diagnostic.get("sigma0_scope_template_availability") or {}
    failure_policy = render_state.get("failure_policy") or "unknown"
    state = render_state.get("state") or "missing_payload"
    if state == "missing_payload":
        status_line = "PROTECTED RENDER-CONTRACT ERROR: diagnostic payload is missing or incomplete"
        detail_lines = [
            "configured protected mode has no renderable protected diagnostic payload",
            "missing histograms: {}".format(
                ", ".join(render_state.get("missing_histograms") or []) or "diagnostic payload"
            ),
        ]
    else:
        if bool(diagnostic.get("fallback_attempted")) and not bool(
            diagnostic.get("fallback_used")
        ):
            status_line = "PROTECTED FIT UNAVAILABLE: Lambda-only fallback was attempted but failed"
        else:
            status_line = "PROTECTED FIT UNAVAILABLE: staged kaon pi-delta pages remain suppressed"
        detail_lines = [
            "protected status: {}".format(render_state.get("diagnostic_status") or "unknown"),
            "failure reason: {}".format(diagnostic.get("failure_reason") or "not reported"),
            "fit variant: {} (selected={})".format(
                diagnostic.get("fit_variant") or "unknown",
                diagnostic.get("selected_fit_variant") or "none",
            ),
            "Lambda-only fallback: attempted={} used={} reason={}".format(
                bool(diagnostic.get("fallback_attempted")),
                bool(diagnostic.get("fallback_used")),
                diagnostic.get("fallback_reason") or "none",
            ),
            "K-Sigma0 source availability: {} ({})".format(
                sigma0_source_availability.get("status") or "not reported",
                sigma0_source_availability.get("reason") or "none",
            ),
            "K-Sigma0 scope-template availability: {} ({})".format(
                sigma0_scope_availability.get("status") or "not reported",
                sigma0_scope_availability.get("reason") or "none",
            ),
        ]

    sigma0_source = component_payload.get("kaon_sigma0_source_diagnostics") or {}
    if sigma0_source:
        identity = sigma0_source.get("source_identity") or {}
        detail_lines.extend(
            [
                "K-Sigma0 source: {} Q{} W{} {} {}".format(
                    sigma0_source.get("resolution_source", "unknown"),
                    identity.get("Q2", ""),
                    identity.get("W", ""),
                    identity.get("EPSSET", ""),
                    identity.get("phi_setting", ""),
                ),
                "K-Sigma0 environment: {}".format(
                    sigma0_source.get("requested_environment_variable") or "not recorded"
                ),
                "K-Sigma0 requested: {}".format(sigma0_source.get("requested_root")),
                "K-Sigma0 resolved: {}".format(sigma0_source.get("resolved_root")),
                "K-Sigma0 loader: exists={} tree_exists={} entries={} missing_branches={} reason={}".format(
                    sigma0_source.get("path_exists"),
                    sigma0_source.get("tree_exists"),
                    sigma0_source.get("tree_entries"),
                    len(sigma0_source.get("missing_required_branches") or []),
                    sigma0_source.get("fallback_reason") or "none",
                ),
            ]
        )

    if str(failure_policy).strip().lower() == "zero_pi_delta":
        detail_lines.append("No pi-delta subtraction applied for this scope")
    detail_lines.extend(
        [
            "failure policy: {}".format(failure_policy),
            "Lambda gauge: status={} solver={} quality={} reason={}".format(
                lambda_gauge.get("status") or "not attempted",
                bool(lambda_gauge.get("solver_success")),
                bool(lambda_gauge.get("quality_passed")),
                lambda_gauge.get("failure_reason") or "none",
            ),
            "applied A_delta={}".format(
                _format_fit_number(diagnostic.get("applied_A_delta"))
            ),
            "template availability: K-Lambda={} K-Sigma0={} pi-delta={}".format(
                "yes" if bool(availability.get(KAON_SIGNAL_TEMPLATE_NAME)) else "no",
                "yes" if bool(availability.get(KAON_SIGMA0_TEMPLATE_NAME)) else "no",
                "yes" if bool(availability.get("pi_delta")) else "no",
            ),
        ]
    )
    _print_component_text_page(
        pdf_name,
        "{}Signal-protected final pi-delta fit status".format(title_prefix),
        [
            "scope: {}".format(component_payload.get("analysis_scope", "unknown")),
            status_line,
        ],
        detail_lines,
    )


def _print_legacy_protected_pi_delta_pages(
    pdf_name,
    component_fit_result,
    render_state,
    title_prefix="",
    cut_window=None,
):
    """Print the two authoritative protected-kaon pages for a successful fit."""
    protected = render_state.get("diagnostic") or {}
    metrics = protected.get("fit_metrics") or {}
    amplitudes = protected.get("signal_amplitudes") or {}
    matrix = protected.get("matrix_diagnostics") or {}
    preservation = protected.get("signal_preservation") or {}
    closure = protected.get("closure") or {}
    reference_integrity = protected.get("lambda_reference_integrity") or {}
    correlations = protected.get("high_template_correlation_warnings") or []
    early_integrity = protected.get("early_amplitudes_frozen_integrity") or {}
    fit_variant = str(protected.get("fit_variant") or "").strip()

    if fit_variant == "lambda_only_protected_fallback":
        sigma0_source = protected.get("sigma0_source_availability") or {}
        sigma0_identity = sigma0_source.get("source_identity") or {}
        sigma0_source_lines = []
        if sigma0_identity or sigma0_source.get("requested_environment_variable"):
            sigma0_source_lines.extend(
                [
                    "K-Sigma0 source identity: Q{} W{} {} {}".format(
                        sigma0_identity.get("Q2", ""),
                        sigma0_identity.get("W", ""),
                        sigma0_identity.get("EPSSET", ""),
                        sigma0_identity.get("phi_setting", ""),
                    ),
                    "K-Sigma0 environment: {}".format(
                        sigma0_source.get("requested_environment_variable") or "not recorded"
                    ),
                ]
            )
        _print_component_overlay_page(
            pdf_name,
            component_fit_result.get("H_pi_delta_protected_fit_input"),
            "R_{pre#Delta}",
            "{}Signal-protected final #pi#Delta fit - K-Lambda-only fallback".format(
                title_prefix
            ),
            [
                (component_fit_result.get("H_pi_delta_protected_k_lambda"), "fitted K-Lambda", ROOT.kBlue + 1, 1),
                (component_fit_result.get("H_pi_delta_protected_pi_delta"), "protected pi-delta", ROOT.kAzure + 2, 1),
                (component_fit_result.get("H_pi_delta_protected_fit_total"), "Lambda + pi-delta total", ROOT.kGreen + 2, 2),
            ],
            [
                "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
                "K-Lambda is fitted and retained. K-Sigma0 template unavailable. ONLY pi-delta is subtracted.",
                "K-Sigma0 reason: {}".format(
                    protected.get("sigma0_availability_reason")
                    or protected.get("fallback_reason")
                    or "not reported"
                ),
                *sigma0_source_lines,
                "A_delta protected={}  legacy diagnostic={}".format(
                    _format_fit_number(protected.get("applied_A_delta")),
                    _format_fit_number(protected.get("legacy_staged_A_delta")),
                ),
                "S_lambda={}  S_sigma0=not fitted".format(
                    _format_fit_number(amplitudes.get(KAON_SIGNAL_TEMPLATE_NAME))
                ),
                "chi2/ndf={}  p={}  bins={}".format(
                    _format_fit_metric(metrics.get("chi2_ndf")),
                    _format_fit_metric(metrics.get("fit_p_value")),
                    metrics.get("n_fit_bins", 0),
                ),
                "rank={}  cond={}".format(
                    matrix.get("weighted_design_effective_rank"),
                    _format_fit_metric(matrix.get("weighted_design_condition_number")),
                ),
                "Lambda/pi-delta correlation={}".format(
                    _format_fit_metric(
                        ((matrix.get("template_correlation_matrix") or {}).get(KAON_SIGNAL_TEMPLATE_NAME) or {}).get("pi_delta")
                    )
                ),
                "K-Lambda reference integrity={}".format(
                    "pass" if bool(reference_integrity.get("shape_identical")) else "unavailable/fail"
                ),
                "early pi-n/pi-SIDIS amplitudes frozen={}".format(
                    "pass" if bool(early_integrity.get("unchanged")) else "fail"
                ),
                "two-template closure={}".format(
                    "pass" if bool((closure.get("protected_two_component_model") or {}).get("passed")) else "fail"
                ),
            ],
            cut_window=cut_window,
        )
        _print_component_overlay_page(
            pdf_name,
            component_fit_result.get("H_pi_delta_protected_fit_input"),
            "R_{pre#Delta}",
            "{}Protected #pi#Delta subtraction - K-Lambda-only fallback".format(
                title_prefix
            ),
            [
                (component_fit_result.get("H_pi_delta_protected_pi_delta"), "fitted pi-delta removed", ROOT.kAzure + 2, 1),
                (component_fit_result.get("H_pi_delta_protected_after_subtraction"), "R_{pre#Delta} - pi-delta", ROOT.kOrange + 7, 2),
            ],
            [
                "K-Lambda is fitted and retained. K-Sigma0 template unavailable. ONLY pi-delta is subtracted.",
                "physics output: R_preDelta - pi-delta only",
                "fit residual is diagnostic-only and is not a cleaned-kaon spectrum",
                "delta-only physics closure={}".format(
                    "pass" if bool((closure.get("delta_only_physics_output") or {}).get("passed")) else "fail"
                ),
            ],
            cut_window=cut_window,
        )
        return

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_pi_delta_protected_fit_input"),
        "R_{pre#Delta}",
        "{}Signal-protected final #pi#Delta fit".format(title_prefix),
        [
            (component_fit_result.get("H_pi_delta_protected_k_lambda"), "fitted K-Lambda", ROOT.kBlue + 1, 1),
            (component_fit_result.get("H_pi_delta_protected_k_sigma0"), "fitted K-Sigma0", ROOT.kCyan + 2, 1),
            (component_fit_result.get("H_pi_delta_protected_pi_delta"), "protected pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pi_delta_protected_fit_total"), "three-template total", ROOT.kGreen + 2, 2),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "K-Lambda and K-Sigma0 signal-protection templates; only pi-delta is subtracted",
            "A_delta protected={}  legacy diagnostic={}".format(
                _format_fit_number(protected.get("applied_A_delta")),
                _format_fit_number(protected.get("legacy_staged_A_delta")),
            ),
            "S_lambda={}  S_sigma0={}".format(
                _format_fit_number(amplitudes.get(KAON_SIGNAL_TEMPLATE_NAME)),
                _format_fit_number(amplitudes.get(KAON_SIGMA0_TEMPLATE_NAME)),
            ),
            "chi2/ndf={}  p={}  bins={}".format(
                _format_fit_metric(metrics.get("chi2_ndf")),
                _format_fit_metric(metrics.get("fit_p_value")),
                metrics.get("n_fit_bins", 0),
            ),
            "rank={}  cond={}".format(
                matrix.get("weighted_design_effective_rank"),
                _format_fit_metric(matrix.get("weighted_design_condition_number")),
            ),
            "template-correlation warnings={}".format(len(correlations)),
            "K-Lambda reference integrity={}".format(
                "pass" if bool(reference_integrity.get("shape_identical")) else "unavailable/fail"
            ),
            "early pi-n/pi-SIDIS amplitudes frozen={}".format(
                "pass" if bool(early_integrity.get("unchanged")) else "fail"
            ),
            "Lambda/Sigma preservation={}/{}".format(
                _format_fit_metric((preservation.get(KAON_SIGNAL_TEMPLATE_NAME) or {}).get("pi_delta_removed_fraction")),
                _format_fit_metric((preservation.get(KAON_SIGMA0_TEMPLATE_NAME) or {}).get("pi_delta_removed_fraction")),
            ),
            "three-template closure={}".format(
                "pass" if bool((closure.get("protected_three_component_model") or {}).get("passed")) else "fail"
            ),
        ],
        cut_window=cut_window,
    )

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_pi_delta_protected_fit_input"),
        "R_{pre#Delta}",
        "{}Protected #pi#Delta subtraction - only #pi#Delta removed".format(title_prefix),
        [
            (component_fit_result.get("H_pi_delta_protected_pi_delta"), "fitted pi-delta removed", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pi_delta_protected_after_subtraction"), "R_{pre#Delta} - pi-delta", ROOT.kOrange + 7, 2),
        ],
        [
            "K-Lambda and K-Sigma0 are fitted and retained; they are not subtracted",
            "physics output: R_preDelta - pi-delta only",
            "fit residual is diagnostic-only and is not a cleaned-kaon spectrum",
            "delta-only physics closure={}".format(
                "pass" if bool((closure.get("delta_only_physics_output") or {}).get("passed")) else "fail"
            ),
        ],
        cut_window=cut_window,
    )


def _print_protected_pi_delta_pages(
    pdf_name,
    component_fit_result,
    render_state,
    title_prefix="",
    cut_window=None,
):
    """Render the authoritative gauge, constrained-fit, and retention pages."""
    protected = render_state.get("diagnostic") or {}
    gauge = protected.get("lambda_gauge") or {}
    preservation = protected.get("lambda_preservation") or {}
    metrics = protected.get("fit_metrics") or {}
    constraint_metrics = protected.get("constraint_metrics") or {}
    amplitudes = protected.get("signal_amplitudes") or {}
    proposed = protected.get("proposed_amplitudes") or {}
    matrix = protected.get("matrix_diagnostics") or {}
    fit_variant = str(protected.get("fit_variant") or "").strip()
    scope = component_fit_result.get("analysis_scope", "unknown")
    title_prefix = "{} ".format(title_prefix.strip()) if title_prefix else ""
    gauge_window = gauge.get("window")
    gauge_hist = component_fit_result.get("H_pi_delta_lambda_gauge")

    _print_protected_overlay_page(
        pdf_name,
        component_fit_result.get("H_pi_delta_protected_fit_input"),
        "R_{pre#Delta}",
        "{}K-Lambda pre-pi-delta gauge".format(title_prefix),
        [(gauge_hist, "immutable K-Lambda × selected anchor", ROOT.kBlue + 1, 1)],
        [
            "scope: {}".format(scope),
            "K-Lambda template available={}".format(
                "yes" if bool(protected.get("k_lambda_scope_template_available")) else "no"
            ),
            "window: {} ({})".format(_format_window_list([gauge_window] if gauge_window else []), gauge.get("window_source") or "unknown"),
            "S_lambda,gauge={} +/- {}  status={}".format(
                _format_fit_number(gauge.get("amplitude")),
                _format_fit_number(gauge.get("amplitude_sigma")),
                gauge.get("status") or "unknown",
            ),
            "gauge chi2/ndf={}  p={}  bins={}  quality={}".format(
                _format_fit_metric(gauge.get("chi2_ndf")),
                _format_fit_metric(gauge.get("p_value")),
                gauge.get("fit_bins", 0),
                "PASS" if bool(gauge.get("quality_passed")) else "FAIL",
            ),
            "Lambda gauge quality is diagnostic only; protected anchor={}".format(
                (protected.get("lambda_constraint") or {}).get("source") or "unavailable"
            ),
            "Lambda data/gauge yield={}/{}".format(
                _format_fit_number(gauge.get("data_integral_window")),
                _format_fit_number(gauge.get("gauge_predicted_yield_window")),
            ),
        ],
        window=gauge_window,
    )

    overlays = [
        (component_fit_result.get("H_pi_delta_protected_k_lambda"), "constrained K-Lambda", ROOT.kBlue + 1, 1),
    ]
    if fit_variant == "lambda_sigma0_protected":
        overlays.append((component_fit_result.get("H_pi_delta_protected_k_sigma0"), "K-Sigma0", ROOT.kCyan + 2, 1))
    overlays.extend(
        [
            (component_fit_result.get("H_pi_delta_protected_pi_delta"), "applied pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pi_delta_protected_fit_total"), "protected fit total", ROOT.kGreen + 2, 2),
        ]
    )
    _print_protected_overlay_page(
        pdf_name,
        component_fit_result.get("H_pi_delta_protected_fit_input"),
        "R_{pre#Delta}",
        "{}Signal-protected final pi-delta fit - K-Lambda gauged{}".format(
            title_prefix,
            " + K-Sigma0" if fit_variant == "lambda_sigma0_protected" else "",
        ),
        overlays,
        [
            "scope: {}  requested/applied mode={}/{}".format(
                scope,
                constraint_metrics.get("requested_mode") or "unknown",
                constraint_metrics.get("applied_mode") or constraint_metrics.get("mode") or "unknown",
            ),
            "K-Lambda template available={}  protected fit variant={}".format(
                "yes" if bool(protected.get("k_lambda_scope_template_available")) else "no",
                fit_variant or "unknown",
            ),
            "Lambda constraint source={}  sigma={}".format(
                constraint_metrics.get("source") or "unknown",
                _format_fit_number(constraint_metrics.get("effective_sigma")),
            ),
            "S_lambda constrained/gauge={}/{}  ratio={}".format(
                _format_fit_number(amplitudes.get(KAON_SIGNAL_TEMPLATE_NAME)),
                _format_fit_number(gauge.get("amplitude")),
                _format_fit_metric(
                    float(amplitudes[KAON_SIGNAL_TEMPLATE_NAME]) / float(gauge["amplitude"])
                    if _is_finite_number(amplitudes.get(KAON_SIGNAL_TEMPLATE_NAME)) and float(gauge.get("amplitude") or 0.0) > 0.0 else None
                ),
            ),
            "A_delta proposed/applied={}/{}  legacy diagnostic={}".format(
                _format_fit_number(proposed.get("pi_delta")),
                _format_fit_number(protected.get("applied_A_delta")),
                _format_fit_number(protected.get("legacy_staged_A_delta")),
            ),
            "spectrum chi2/ndf={}  p={}  Nfree={}".format(
                _format_fit_metric(metrics.get("chi2_ndf")),
                _format_fit_metric(metrics.get("fit_p_value")),
                metrics.get("n_free_spectrum_parameters", 0),
            ),
            "rank={}  cond={}  fit quality={}".format(
                matrix.get("weighted_design_effective_rank"),
                _format_fit_metric(matrix.get("weighted_design_condition_number")),
                "PASS" if bool(protected.get("fit_quality_passed")) else "FAIL",
            ),
        ],
        window=(protected.get("fit_window") or [None])[0],
    )

    _print_protected_overlay_page(
        pdf_name,
        component_fit_result.get("H_pi_delta_protected_fit_input"),
        "R_{pre#Delta}",
        "{}Protected pi-delta subtraction - Lambda preservation".format(title_prefix),
        [
            (gauge_hist, "K-Lambda gauge", ROOT.kBlue + 1, 1),
            (component_fit_result.get("H_pi_delta_protected_pi_delta"), "pi-delta removed", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pi_delta_protected_after_subtraction"), "R_{pre#Delta} - pi-delta", ROOT.kOrange + 7, 2),
        ],
        [
            "gauge/pre/removed/after={}/{}/{}/{}".format(
                _format_fit_number(preservation.get("lambda_gauge_predicted_yield")),
                _format_fit_number(preservation.get("lambda_pre_delta_yield")),
                _format_fit_number(preservation.get("lambda_pi_delta_removed_yield")),
                _format_fit_number(preservation.get("lambda_after_delta_yield")),
            ),
            "after/gauge={}  removed/gauge={}".format(
                _format_fit_metric(preservation.get("lambda_after_over_gauge")),
                _format_fit_metric(preservation.get("lambda_removed_fraction_of_gauge")),
            ),
            "minimum retention={}  bound={}  gate={}".format(
                _format_fit_number(preservation.get("minimum_required_retention")),
                _format_fit_number(preservation.get("a_delta_max")),
                "PASS" if bool(preservation.get("gate_passed")) else "FAIL: {}".format(preservation.get("gate_reason") or "unknown"),
            ),
            "ONLY pi-delta is subtracted; K-Lambda and K-Sigma0 remain signal templates.",
        ],
        window=gauge_window,
    )


def _print_residual_shift_template_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
    protected_kaon_mode=False,
):
    shift_summaries = (component_fit_result or {}).get("residual_component_shift_summaries") or {}
    target_specs = (
        ("pion_control", "pion", "pion-control"),
        ("kaon_nosub", "kaon", "kaon no-sub"),
    )
    for summary_key, hist_prefix, sample_label in target_specs:
        summary = shift_summaries.get(summary_key) or {}
        if not _residual_shift_summary_has_content(summary):
            continue
        active_components = list(summary.get("active_components") or [])
        if not active_components:
            active_components = list((summary.get("selected_shift_point") or {}).keys())
        if not active_components:
            continue
        per_component = summary.get("per_component") or {}
        for component_name in active_components:
            suffix = _component_hist_suffix(component_name)
            original_hist = component_fit_result.get(
                "H_{}_shift_original_{}".format(hist_prefix, suffix)
            )
            shifted_hist = component_fit_result.get(
                "H_{}_shift_selected_{}".format(hist_prefix, suffix)
            )
            if original_hist is None or shifted_hist is None:
                continue

            component_label = _component_plot_label(component_name)
            component_color = _component_plot_color(component_name)
            shift_diag = per_component.get(component_name) or {}
            selected_delta = (summary.get("selected_shift_point") or {}).get(component_name)

            canvas = ROOT.TCanvas(
                "c_shift_{}_{}".format(hist_prefix, suffix),
                "",
                900,
                900,
            )
            canvas.Divide(1, 2)

            top_pad = canvas.cd(1)
            top_pad.SetBottomMargin(0.12)
            original_clone = _clone_hist(original_hist, "{}_orig".format(original_hist.GetName()))
            shifted_clone = _clone_hist(shifted_hist, "{}_shift".format(shifted_hist.GetName()))
            alignment_only_label = (
                " (alignment diagnostic only)"
                if protected_kaon_mode
                and hist_prefix == "kaon"
                and component_name == "pi_delta"
                else ""
            )
            original_clone.SetTitle(
                "{}{} residual-shift template: {}{}".format(
                    title_prefix,
                    sample_label,
                    component_label,
                    alignment_only_label,
                )
            )
            _style_overlay_hist(original_clone, ROOT.kBlue + 1, line_style=2)
            _style_overlay_hist(shifted_clone, component_color, line_style=1)
            top_y_max = max(original_clone.GetMaximum(), shifted_clone.GetMaximum(), 0.0)
            if top_y_max <= 0.0:
                top_y_max = 1.0
            original_clone.SetMaximum(1.20 * top_y_max)
            original_clone.SetMinimum(0.0)
            original_clone.Draw("hist")
            shifted_clone.Draw("hist same")

            top_legend = ROOT.TLegend(0.58, 0.70, 0.88, 0.88)
            top_legend.SetBorderSize(0)
            top_legend.SetFillStyle(0)
            top_legend.AddEntry(original_clone, "original aligned template", "l")
            top_legend.AddEntry(shifted_clone, "selected shifted template", "l")
            top_legend.Draw()

            stats_box = ROOT.TPaveText(0.14, 0.56, 0.52, 0.88, "NDC")
            stats_box.SetBorderSize(0)
            stats_box.SetFillStyle(0)
            stats_box.SetTextAlign(12)
            stats_box.SetTextSize(0.026)
            stats_box.AddText("fit target: {}".format(summary_key))
            stats_box.AddText("component: {}".format(component_label))
            if alignment_only_label:
                stats_box.AddText("legacy pi-delta alignment QA only; not applied subtraction")
            stats_box.AddText(
                "selected delta={} {}".format(
                    _format_fit_number(selected_delta),
                    str(summary.get("units") or "GeV"),
                )
            )
            stats_box.AddText("selection metric: {}".format(summary.get("selection_metric") or "n/a"))
            stats_box.AddText(
                "renorm factor={}".format(
                    _format_fit_number(shift_diag.get("shift_renormalization_factor"))
                )
            )
            stats_box.AddText(
                "lost frac={}".format(
                    _format_fit_metric(shift_diag.get("lost_integral_fraction"))
                )
            )
            stats_box.AddText(
                "bound hit={}".format(
                    "yes" if bool(shift_diag.get("shift_bound_hit_flag")) else "no"
                )
            )
            stats_box.Draw()

            bottom_pad = canvas.cd(2)
            bottom_pad.SetTopMargin(0.08)
            bottom_pad.SetBottomMargin(0.12)
            diff_hist = _clone_hist(shifted_hist, "{}_diff".format(shifted_hist.GetName()))
            diff_hist.Add(original_hist, -1.0)
            diff_hist.SetTitle("Shifted - original template")
            _style_overlay_hist(diff_hist, component_color, line_style=1)
            diff_y_max = max(abs(diff_hist.GetMaximum()), abs(diff_hist.GetMinimum()), 0.0)
            if diff_y_max <= 0.0:
                diff_y_max = 1.0
            diff_hist.SetMaximum(1.20 * diff_y_max)
            diff_hist.SetMinimum(-1.20 * diff_y_max)
            diff_hist.Draw("hist")
            zero_line = ROOT.TLine(
                float(diff_hist.GetXaxis().GetXmin()),
                0.0,
                float(diff_hist.GetXaxis().GetXmax()),
                0.0,
            )
            zero_line.SetLineColor(ROOT.kGray + 2)
            zero_line.SetLineStyle(2)
            zero_line.SetLineWidth(2)
            zero_line.Draw("same")

            bottom_legend = ROOT.TLegend(0.58, 0.78, 0.88, 0.88)
            bottom_legend.SetBorderSize(0)
            bottom_legend.SetFillStyle(0)
            bottom_legend.AddEntry(diff_hist, "shifted - original", "l")
            bottom_legend.Draw()

            canvas.Print(pdf_name)
            canvas.Close()


def _print_residual_shift_summary_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
):
    shift_summaries = (component_fit_result or {}).get("residual_component_shift_summaries") or {}
    target_specs = (
        ("pion_control", "pion-control"),
        ("kaon_nosub", "kaon no-sub"),
    )
    for summary_key, sample_label in target_specs:
        summary = shift_summaries.get(summary_key) or {}
        if not _residual_shift_summary_has_content(summary):
            continue
        candidate_summaries = list(summary.get("candidate_summaries") or [])
        body_lines = [
            "enabled={}".format(bool(summary.get("enabled"))),
            "mode={}".format(summary.get("mode") or "n/a"),
            "units={}".format(summary.get("units") or "n/a"),
            "selection_metric={}".format(summary.get("selection_metric") or "n/a"),
            "requested_components={}".format(", ".join(summary.get("requested_components") or []) or "none"),
            "active_components={}".format(", ".join(summary.get("active_components") or []) or "none"),
            "selected_shift_point={}".format(_format_shift_point_map(summary.get("selected_shift_point") or {})),
            "selected_reason={}".format(summary.get("selected_shift_reason") or "n/a"),
            "candidate_count={}".format(int(summary.get("candidate_count", 0) or 0)),
            "cleanup_chi2_ndf={}".format(
                _format_fit_metric(((summary.get("cleanup_metrics") or {}).get("chi2_ndf")))
            ),
            "cleanup_residual_integral={}".format(
                _format_fit_number(((summary.get("cleanup_metrics") or {}).get("residual_integral")))
            ),
            "cleanup_pull_rms={}".format(
                _format_fit_metric(((summary.get("cleanup_metrics") or {}).get("pull_rms")))
            ),
            "warnings={}".format(", ".join(summary.get("warnings") or []) or "none"),
        ]
        for candidate in candidate_summaries[:10]:
            cleanup_metrics = candidate.get("cleanup_metrics") or {}
            body_lines.append(
                "cand {}: shifts={} score={} acc={} cleanup chi2/ndf={} yield={}".format(
                    int(candidate.get("candidate_index", 0) or 0),
                    _format_shift_point_map(candidate.get("shift_point") or {}),
                    _format_fit_metric(candidate.get("score")),
                    "yes" if bool(candidate.get("accepted")) else "no",
                    _format_fit_metric(cleanup_metrics.get("chi2_ndf")),
                    _format_fit_number(cleanup_metrics.get("corrected_yield")),
                )
            )
        _print_component_text_page(
            pdf_name,
            "{}{} residual-shift summary".format(title_prefix, sample_label),
            [
                "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
                "fit target: {}".format(summary_key),
            ],
            body_lines,
        )


def _print_residual_shift_scan_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
    protected_kaon_mode=False,
):
    shift_summaries = (component_fit_result or {}).get("residual_component_shift_summaries") or {}
    target_specs = (
        ("pion_control", "pion-control"),
        ("kaon_nosub", "kaon no-sub"),
    )
    for summary_key, sample_label in target_specs:
        summary = shift_summaries.get(summary_key) or {}
        candidate_summaries = list(summary.get("candidate_summaries") or [])
        active_components = list(summary.get("active_components") or [])
        if len(active_components) != 1 or len(candidate_summaries) <= 1:
            continue
        component_name = active_components[0]
        graph_points = []
        yield_points = []
        for candidate in candidate_summaries:
            shift_value = (candidate.get("shift_point") or {}).get(component_name)
            score_value = candidate.get("score")
            corrected_yield = ((candidate.get("cleanup_metrics") or {}).get("corrected_yield"))
            if _is_finite_number(shift_value) and _is_finite_number(score_value):
                graph_points.append((float(shift_value), float(score_value)))
            if _is_finite_number(shift_value) and _is_finite_number(corrected_yield):
                yield_points.append((float(shift_value), float(corrected_yield)))
        if len(graph_points) <= 1:
            continue

        graph_points.sort(key=lambda item: item[0])
        yield_points.sort(key=lambda item: item[0])
        objective_graph = ROOT.TGraph(len(graph_points))
        for idx, (x_value, y_value) in enumerate(graph_points):
            objective_graph.SetPoint(idx, x_value, y_value)
        objective_graph.SetLineColor(ROOT.kBlue + 1)
        objective_graph.SetLineWidth(2)
        objective_graph.SetMarkerColor(ROOT.kBlue + 1)
        objective_graph.SetMarkerStyle(20)

        yield_graph = ROOT.TGraph(len(yield_points))
        for idx, (x_value, y_value) in enumerate(yield_points):
            yield_graph.SetPoint(idx, x_value, y_value)
        yield_graph.SetLineColor(ROOT.kMagenta + 2)
        yield_graph.SetLineWidth(2)
        yield_graph.SetMarkerColor(ROOT.kMagenta + 2)
        yield_graph.SetMarkerStyle(21)

        canvas = ROOT.TCanvas("c_shift_scan_{}".format(summary_key), "", 900, 900)
        canvas.Divide(1, 2)

        top_pad = canvas.cd(1)
        top_pad.SetBottomMargin(0.12)
        alignment_only_label = (
            " (alignment diagnostic only)"
            if protected_kaon_mode
            and summary_key == "kaon_nosub"
            and component_name == "pi_delta"
            else ""
        )
        objective_graph.SetTitle(
            "{}{} residual-shift scan: {}{}".format(
                title_prefix,
                sample_label,
                _component_plot_label(component_name),
                alignment_only_label,
            )
        )
        objective_graph.GetXaxis().SetTitle("delta MM [{}]".format(summary.get("units") or "GeV"))
        objective_graph.GetYaxis().SetTitle(str(summary.get("selection_metric") or "score"))
        objective_graph.Draw("ALP")

        bottom_pad = canvas.cd(2)
        bottom_pad.SetTopMargin(0.08)
        bottom_pad.SetBottomMargin(0.12)
        yield_graph.SetTitle("Cleanup-region corrected-yield proxy vs shift")
        yield_graph.GetXaxis().SetTitle("delta MM [{}]".format(summary.get("units") or "GeV"))
        yield_graph.GetYaxis().SetTitle("cleanup residual integral")
        yield_graph.Draw("ALP")

        canvas.Print(pdf_name)
        canvas.Close()


def print_particle_subtraction_component_template_pages(
    pdf_name,
    component_shape_payload,
    title_prefix="",
    cut_window=None,
    kaon_signal_payload=None,
    kaon_sigma0_payload=None,
    page_manifest=None,
    page_id_prefix=None,
    authoritative=False,
):
    component_shape_payload = component_shape_payload or {}
    component_map = component_shape_payload.get("components") or {}
    kaon_signal_payload = kaon_signal_payload or {}
    kaon_sigma0_payload = kaon_sigma0_payload or {}

    component_specs = [
        ("pi_n", "pi-n", ROOT.kRed + 1),
        ("pi_sidis", "pi-SIDIS", ROOT.kMagenta + 2),
        ("pi_delta", "pi-delta", ROOT.kAzure + 2),
    ]
    aux_specs = [
        ("k_lambda_signal", "K-Lambda", ROOT.kBlue + 1, kaon_signal_payload),
        ("k_sigma0_signal", "K-Sigma0", ROOT.kCyan + 2, kaon_sigma0_payload),
    ]

    def _build_page(hist_key):
        base_hist = None
        base_label = None
        overlay_specs = []
        stats_lines = [
            "mode: {}".format(component_shape_payload.get("mode", "unknown")),
            "tree: {}".format(component_shape_payload.get("tree_name", "unknown")),
        ]

        for component_name, label, color in component_specs:
            component_payload = component_map.get(component_name) or {}
            hist = component_payload.get(hist_key)
            if hist is None:
                continue
            if base_hist is None:
                base_hist = hist
                base_label = label
            else:
                overlay_specs.append((hist, label, color, 1))
            stats_lines.append(
                _format_component_template_diag_line(
                    label,
                    (component_shape_payload.get("diagnostics") or {}).get(component_name),
                )
            )

        for _, label, color, aux_payload in aux_specs:
            hist = aux_payload.get(hist_key)
            if hist is None:
                continue
            if base_hist is None:
                base_hist = hist
                base_label = label
            else:
                overlay_specs.append((hist, label, color, 2))
            stats_lines.append(
                _format_component_template_diag_line(label, aux_payload.get("diagnostics"))
            )

        return base_hist, base_label, overlay_specs, stats_lines

    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)

    full_base, full_label, full_overlays, full_stats = _build_page("setting_shape_full")
    if full_base is not None:
        if _print_component_overlay_page(
            pdf_name,
            full_base,
            full_label,
            "{}SIMC component templates (full MM range)".format(title_prefix),
            full_overlays,
            full_stats,
            cut_window=cut_window,
        ):
            _record_component_page(
                page_manifest,
                page_id_prefix,
                "component_templates_full",
                scope="setting-wide",
                authoritative=authoritative,
            )

    cut_base, cut_label, cut_overlays, cut_stats = _build_page("setting_shape")
    if cut_base is not None:
        if _print_component_overlay_page(
            pdf_name,
            cut_base,
            cut_label,
            "{}SIMC component templates (analysis MM-cut window)".format(title_prefix),
            cut_overlays,
            cut_stats,
            cut_window=cut_window,
        ):
            _record_component_page(
                page_manifest,
                page_id_prefix,
                "component_templates_cut",
                scope="setting-wide",
                authoritative=authoritative,
            )


def print_particle_subtraction_kaon_lambda_comparison_page(
    pdf_name,
    component_fit_result,
    component_payload=None,
    title_prefix="",
    cut_window=None,
    page_manifest=None,
    page_id_prefix=None,
    authoritative=False,
    proposal_payload=None,
):
    """Render the mandatory K-Lambda comparison independently of pion acceptance.

    The fit result owns the immutable SIMC reference.  The optional application
    payload supplies the actual scope's after-pion spectrum, so a per-t parent
    comparison never silently falls back to a setting-wide input.
    """
    if not isinstance(component_fit_result, dict):
        return False

    application_payload = component_payload if isinstance(component_payload, dict) else {}
    proposal_application_payload = proposal_payload if isinstance(proposal_payload, dict) else {}
    scope_label = (
        application_payload.get("analysis_scope")
        or application_payload.get("analysis_scope_label")
        or component_fit_result.get("analysis_scope")
        or "unknown"
    )
    after_hist = application_payload.get("H_MM_nosub_after_pion_subtraction")
    before_hist = application_payload.get("H_MM_nosub_before_pion_subtraction")
    proposed_after_hist = proposal_application_payload.get("H_MM_nosub_after_pion_subtraction")
    proposed_before_hist = proposal_application_payload.get("H_MM_nosub_before_pion_subtraction")
    target_hist = (
        after_hist
        or before_hist
        or proposed_after_hist
        or proposed_before_hist
        or component_fit_result.get("H_kaon_nosub_input")
    )
    if target_hist is None:
        raise RuntimeError(
            "K-Lambda SIMC comparison requires the {} scope kaon spectrum".format(
                scope_label
            )
        )

    lambda_reference_hist, lambda_reference_scale, lambda_reference_source, lambda_reference_normalization = (
        _resolve_kaon_lambda_reference_for_plot(
            component_fit_result,
            target_hist,
            cut_window,
            scope_label,
            "H_k_lambda_reference_application",
        )
    )
    protected_lambda_diagnostic = (
        (component_fit_result.get("diagnostics") or {}).get("pi_delta_signal_protected_fit")
        or ((component_fit_result.get("diagnostics") or {}).get("kaon") or {}).get(
            "pi_delta_signal_protected_fit"
        )
        or {}
    )
    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)
    has_after_spectrum = after_hist is not None
    has_proposed_after_spectrum = proposed_after_hist is not None
    final_state = application_payload.get("final_application_status") or (
        "applied_component" if bool(application_payload.get("accepted")) else "unavailable"
    )
    if has_after_spectrum:
        base_label = "final after pion subtraction ({})".format(final_state)
    elif has_proposed_after_spectrum:
        base_label = "proposed after pion subtraction (no final spectrum)"
    else:
        base_label = "kaon spectrum before pion subtraction"
    overlays = [(lambda_reference_hist, "K-Lambda SIMC comparison", ROOT.kBlue + 1, 2)]
    if proposed_after_hist is not None and proposed_after_hist is not target_hist:
        overlays.append((proposed_after_hist, "proposed component after-pion", ROOT.kRed + 1, 3))
    emitted = _print_component_overlay_page(
        pdf_name,
        target_hist,
        base_label,
        "{}Part 3 {} vs K-Lambda comparison".format(
            title_prefix,
            "after pion subtraction" if has_after_spectrum else "available kaon spectrum",
        ),
        overlays,
        [
            "scope: {}".format(scope_label),
            "application status={}".format(
                "accepted" if bool(application_payload.get("production_evaluation_accepted", application_payload.get("accepted"))) else "rejected"
            ),
            "comparison spectrum={}".format(
                base_label
            ),
            "final application state={}".format(final_state),
            "K-Lambda comparison full integral={}".format(
                _format_fit_number(_hist_integral(lambda_reference_hist))
            ),
            "K-Lambda comparison source={}".format(lambda_reference_source),
            "comparison normalization={}".format(lambda_reference_normalization),
            "K-Lambda display scale={}".format(
                _format_fit_number(lambda_reference_scale)
                if lambda_reference_scale is not None
                else "n/a"
            ),
            "Lambda gauge quality={} (diagnostic only)".format(
                "PASS"
                if bool(protected_lambda_diagnostic.get("lambda_gauge_quality_passed"))
                else "POOR/UNAVAILABLE"
            ),
            "K-Lambda fitted amplitude={}".format(
                _format_fit_number(application_payload.get("k_lambda_fit_amplitude"))
                if application_payload.get("k_lambda_fit_amplitude") is not None
                else _format_fit_number(component_fit_result.get("S_lambda"))
            ),
        ],
        cut_window=cut_window,
    )
    if emitted:
        _record_component_page(
            page_manifest,
            page_id_prefix,
            "lambda_comparison",
            scope=scope_label,
            authoritative=authoritative,
        )
    return emitted


def print_particle_subtraction_component_application_pages(
    pdf_name,
    component_payload,
    title_prefix="",
    cut_window=None,
    component_fit_result=None,
    include_lambda_page=True,
    page_manifest=None,
    page_id_prefix=None,
    authoritative=False,
):
    if not isinstance(component_payload, dict):
        return
    if not bool(component_payload.get("accepted")):
        if _print_component_application_status_page(
            pdf_name,
            component_payload,
            title_prefix=title_prefix,
        ):
            _record_component_page(
                page_manifest,
                page_id_prefix,
                "application_status",
                scope=(component_payload.get("analysis_scope") or "unknown"),
                authoritative=authoritative,
            )
        if include_lambda_page and isinstance(component_fit_result, dict):
            print_particle_subtraction_kaon_lambda_comparison_page(
                pdf_name,
                component_fit_result,
                component_payload,
                title_prefix=title_prefix,
                cut_window=cut_window,
                page_manifest=page_manifest,
                page_id_prefix=page_id_prefix,
                authoritative=authoritative,
            )
        return

    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)

    diagnostics = component_payload.get("diagnostics") or {}
    scope_label = component_payload.get("analysis_scope") or component_payload.get("analysis_scope_label") or "unknown"
    fit_mode = str(
        component_payload.get("fit_mode")
        or component_payload.get("fit_mode_kaon")
        or component_payload.get("fit_mode_pion")
        or "staged_only"
    ).strip().lower()
    joint_mode_active = fit_mode in ("staged_plus_joint", "staged_plus_regularized_joint")
    fit_render_result = (
        component_fit_result
        if isinstance(component_fit_result, dict)
        else component_payload
    )
    protected_render_state = _resolve_protected_pi_delta_render_state(fit_render_result)
    protected_kaon_mode = bool(
        protected_render_state.get("suppress_deprecated_kaon_pages")
    )
    closure_label = (
        "protected applied"
        if protected_kaon_mode
        else ("refined" if joint_mode_active else "staged")
    )
    kaon_model_label = (
        "protected applied pion-background model"
        if protected_kaon_mode
        else "kaon pion-bg model"
    )
    model_closure = diagnostics.get("model_closure") or {}
    model_closure_stage = diagnostics.get("model_closure_stage") or {}
    event_template_closure = diagnostics.get("event_template_closure") or {}
    weight_diagnostics_stage = diagnostics.get("weight_diagnostics_stage") or {}
    pion_diag = (component_payload.get("diagnostics") or {}).get("pion") or {}
    kaon_diag = (component_payload.get("diagnostics") or {}).get("kaon") or {}
    pion_postfit_scales = pion_diag.get("postfit_component_scales") or {}
    kaon_postfit_scales = kaon_diag.get("postfit_component_scales") or {}
    pion_postrefine_scales = pion_diag.get("postrefine_component_scales") or {}
    kaon_postrefine_scales = kaon_diag.get("postrefine_component_scales") or {}
    kaon_manual_scaling_active = _component_scale_map_has_nonunity(kaon_postfit_scales) or _component_scale_map_has_nonunity(kaon_postrefine_scales)
    if _print_single_hist_page(
        pdf_name,
        component_payload.get("H_pion_weight_vs_MM"),
        "w_pi(MM)",
        "{}Part 3 {} pion weight vs MM".format(title_prefix, closure_label),
        [
            "scope: {}".format(scope_label),
            "fit mode: {}".format(fit_mode or "unknown"),
            "min/max={} / {}".format(
                _format_fit_metric(diagnostics.get("pion_weight_min")),
                _format_fit_metric(diagnostics.get("pion_weight_max")),
            ),
            "mean/rms={} / {}".format(
                _format_fit_metric(diagnostics.get("pion_weight_mean")),
                _format_fit_metric(diagnostics.get("pion_weight_rms")),
            ),
            "unclipped min/max={} / {}".format(
                _format_fit_metric(diagnostics.get("pion_weight_min_unclipped")),
                _format_fit_metric(diagnostics.get("pion_weight_max_unclipped")),
            ),
            "clipped bins={}".format(int(diagnostics.get("pion_weight_clipped_bin_count", 0) or 0)),
            "unsupported bins={}".format(int(diagnostics.get("pion_weight_unsupported_bin_count", 0) or 0)),
            "unsupported MM={}".format(_format_mm_range(diagnostics.get("pion_weight_unsupported_mm_range"))),
            "denominator: pion-control model (unity-scale convention)",
            "numerator: {}{}".format(
                kaon_model_label,
                " after kaon-side scaling" if kaon_manual_scaling_active else ""
            ),
        ],
        cut_window=cut_window,
        line_color=ROOT.kViolet + 1,
    ):
        _record_component_page(page_manifest, page_id_prefix, "pion_weight", scope=scope_label, authoritative=authoritative)

    if (
        not protected_kaon_mode
        and joint_mode_active
        and component_payload.get("H_pion_weight_vs_MM_stage") is not None
    ):
        _print_component_overlay_page(
            pdf_name,
            component_payload.get("H_pion_weight_vs_MM_stage"),
            "staged w_pi(MM)",
            "{}Part 3 staged vs refined pion weight".format(title_prefix),
            [
                (component_payload.get("H_pion_weight_vs_MM"), "refined w_pi(MM)", ROOT.kViolet + 1, 1),
            ],
            [
                "scope: {}".format(scope_label),
                "staged min/max={} / {}".format(
                    _format_fit_metric(weight_diagnostics_stage.get("pion_weight_min")),
                    _format_fit_metric(weight_diagnostics_stage.get("pion_weight_max")),
                ),
                "refined min/max={} / {}".format(
                    _format_fit_metric(diagnostics.get("pion_weight_min")),
                    _format_fit_metric(diagnostics.get("pion_weight_max")),
                ),
                "staged mean/rms={} / {}".format(
                    _format_fit_metric(weight_diagnostics_stage.get("pion_weight_mean")),
                    _format_fit_metric(weight_diagnostics_stage.get("pion_weight_rms")),
                ),
                "refined mean/rms={} / {}".format(
                    _format_fit_metric(diagnostics.get("pion_weight_mean")),
                    _format_fit_metric(diagnostics.get("pion_weight_rms")),
                ),
                "kaon-side post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
                "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            ],
            cut_window=cut_window,
        )

    if _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        kaon_model_label,
        "{}Part 3 {} model closure: weighted pion model vs {}".format(
            title_prefix,
            closure_label,
            kaon_model_label,
        ),
        [
            (
                component_payload.get("H_weighted_pion_control_model"),
                "{} weighted pion-control model".format(closure_label),
                ROOT.kOrange + 7,
                2,
            ),
        ],
        [
            "scope: {}".format(scope_label),
            "signature match={}".format(
                "pass" if bool(model_closure.get("signature_match")) else "fail"
            ),
            "{} integral={}".format(
                kaon_model_label,
                _format_fit_number(model_closure.get("reference_integral"))
            ),
            "weighted pion model integral={}".format(
                _format_fit_number(model_closure.get("comparison_integral"))
            ),
            "integral ratio={}".format(
                _format_fit_metric(model_closure.get("integral_ratio"))
            ),
            "max abs bin diff={} @ MM={}".format(
                _format_fit_number(model_closure.get("max_abs_bin_diff")),
                _format_fit_metric(model_closure.get("max_abs_bin_center")),
            ),
            "{} shown{}.".format(
                kaon_model_label,
                " after kaon-side post-refine scaling" if kaon_manual_scaling_active else " without manual kaon-side scaling"
            ),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            "pion-control post-refine scales: {}".format(_format_component_scale_map(pion_postrefine_scales)),
        ],
        cut_window=cut_window,
    ):
        _record_component_page(page_manifest, page_id_prefix, "model_closure", scope=scope_label, authoritative=authoritative)

    if _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        kaon_model_label,
        "{}Part 3 {} event-template closure vs {}".format(
            title_prefix,
            closure_label,
            kaon_model_label,
        ),
        [
            (component_payload.get("H_pion_subtraction_template_MM_nosub"), "weighted pion template (full)", ROOT.kOrange + 7, 2),
        ],
        [
            "scope: {}".format(scope_label),
            "signature match={}".format(
                "pass" if bool(event_template_closure.get("signature_match")) else "fail"
            ),
            "{} integral={}".format(
                kaon_model_label,
                _format_fit_number(event_template_closure.get("reference_integral"))
            ),
            "weighted event-template integral={}".format(
                _format_fit_number(event_template_closure.get("comparison_integral"))
            ),
            "integral ratio={}".format(
                _format_fit_metric(event_template_closure.get("integral_ratio"))
            ),
            "max abs bin diff={} @ MM={}".format(
                _format_fit_number(event_template_closure.get("max_abs_bin_diff")),
                _format_fit_metric(event_template_closure.get("max_abs_bin_center")),
            ),
            "weighted event-template inherits kaon-side numerator scaling={}".format(
                "yes" if kaon_manual_scaling_active else "no"
            ),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
        ],
        cut_window=cut_window,
    ):
        _record_component_page(page_manifest, page_id_prefix, "event_template_closure", scope=scope_label, authoritative=authoritative)

    if _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_MM_nosub_before_pion_subtraction"),
        "kaon data before pion subtraction",
        "{}Part 3 kaon data vs pion-background models".format(title_prefix),
        [
            (component_payload.get("H_kaon_pion_model"), kaon_model_label, ROOT.kBlue + 1, 1),
            (component_payload.get("H_pion_subtraction_template_MM_nosub"), "weighted pion template", ROOT.kOrange + 7, 2),
        ],
        [
            "scope: {}".format(scope_label),
            "before full integral={}".format(
                _format_fit_number(component_payload.get("kaon_integral_before_pion_sub_full"))
            ),
            "{} integral={}".format(
                kaon_model_label,
                _format_fit_number(diagnostics.get("kaon_pion_model_integral"))
            ),
            "weighted template full integral={}".format(
                _format_fit_number(component_payload.get("weighted_pion_integral_full"))
            ),
            "effective scale={}".format(
                _format_fit_number(component_payload.get("particle_subtraction_effective_scale"))
            ),
            "{} shown{}.".format(
                kaon_model_label,
                " after kaon-side post-refine scaling" if kaon_manual_scaling_active else " without manual kaon-side scaling"
            ),
            "kaon-side post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
        ],
        cut_window=cut_window,
    ):
        _record_component_page(page_manifest, page_id_prefix, "data_vs_pion_model", scope=scope_label, authoritative=authoritative)

    if (
        not protected_kaon_mode
        and joint_mode_active
        and component_payload.get("H_kaon_pion_model_stage") is not None
    ):
        _print_component_overlay_page(
            pdf_name,
            component_payload.get("H_MM_nosub_before_pion_subtraction"),
            "kaon data before pion subtraction",
            "{}Part 3 staged vs refined subtraction comparison".format(title_prefix),
            [
                (component_payload.get("H_kaon_pion_model_stage"), "staged kaon pion-bg model", ROOT.kOrange + 7, 2),
                (component_payload.get("H_kaon_pion_model"), "refined kaon pion-bg model", ROOT.kBlue + 1, 1),
                (component_payload.get("H_MM_nosub_after_pion_subtraction_model_stage"), "staged model-subtracted", ROOT.kMagenta + 2, 2),
                (component_payload.get("H_MM_nosub_after_pion_subtraction_model_final"), "refined model-subtracted", ROOT.kGreen + 2, 1),
            ],
            [
                "scope: {}".format(scope_label),
                "staged model integral={}".format(
                    _format_fit_number((model_closure_stage or {}).get("reference_integral"))
                ),
                "refined model integral={}".format(
                    _format_fit_number((model_closure or {}).get("reference_integral"))
                ),
                "staged weighted model integral={}".format(
                    _format_fit_number((model_closure_stage or {}).get("comparison_integral"))
                ),
                "refined weighted model integral={}".format(
                    _format_fit_number((model_closure or {}).get("comparison_integral"))
                ),
                "kaon-side staged post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
                "kaon-side refined post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            ],
            cut_window=cut_window,
        )

    if _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_MM_nosub_before_pion_subtraction"),
        "before pion subtraction",
        "{}Part 3 full-range MM before/after pion subtraction".format(title_prefix),
        [
            (component_payload.get("H_pion_subtraction_template_MM_nosub"), "weighted pion template", ROOT.kOrange + 7, 2),
            (component_payload.get("H_MM_nosub_after_pion_subtraction"), "after pion subtraction", ROOT.kGreen + 2, 1),
        ],
        [
            "scope: {}".format(scope_label),
            "before full integral={}".format(
                _format_fit_number(component_payload.get("kaon_integral_before_pion_sub_full"))
            ),
            "template full integral={}".format(
                _format_fit_number(component_payload.get("weighted_pion_integral_full"))
            ),
            "after full integral={}".format(
                _format_fit_number(component_payload.get("kaon_integral_after_pion_sub_full"))
            ),
            "weighted template uses kaon-side scaled numerator={}".format(
                "yes" if kaon_manual_scaling_active else "no"
            ),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            "fit validation pion/kaon={}/{}".format(
                "pass" if bool(component_payload.get("fit_validation_pion")) else "fail",
                "pass" if bool(component_payload.get("fit_validation_kaon")) else "fail",
            ),
        ],
        cut_window=cut_window,
    ):
        _record_component_page(page_manifest, page_id_prefix, "before_after", scope=scope_label, authoritative=authoritative)

    if include_lambda_page:
        print_particle_subtraction_kaon_lambda_comparison_page(
            pdf_name,
            fit_render_result,
            component_payload,
            title_prefix=title_prefix,
            cut_window=cut_window,
            page_manifest=page_manifest,
            page_id_prefix=page_id_prefix,
            authoritative=authoritative,
        )


def print_particle_subtraction_component_fit_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
    cut_window=None,
    page_manifest=None,
    page_id_prefix=None,
    authoritative=False,
    refinement_detail_pdf=None,
):
    if not isinstance(component_fit_result, dict):
        return

    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)

    pion_diagnostics = ((component_fit_result.get("diagnostics") or {}).get("pion") or {})
    kaon_diagnostics = ((component_fit_result.get("diagnostics") or {}).get("kaon") or {})
    protected_render_state = _resolve_protected_pi_delta_render_state(component_fit_result)
    protected_kaon_mode = bool(
        protected_render_state.get("suppress_deprecated_kaon_pages")
    )
    pion_stage_amplitudes = pion_diagnostics.get("staged_amplitudes_scaled") or {}
    kaon_stage_amplitudes = kaon_diagnostics.get("staged_amplitudes_scaled") or {}
    pion_stage_validation = pion_diagnostics.get("stage_validation") or pion_diagnostics.get("validation") or {}
    kaon_stage_validation = kaon_diagnostics.get("stage_validation") or kaon_diagnostics.get("validation") or {}
    staged_amplitudes_raw = component_fit_result.get("staged_amplitudes_raw") or {}
    staged_amplitudes_scaled = component_fit_result.get("staged_amplitudes_scaled") or {}
    refined_amplitudes = component_fit_result.get("refined_amplitudes") or {}
    amplitude_shift_fractions = component_fit_result.get("amplitude_shift_fractions") or {}
    pion_stage_raw = (staged_amplitudes_raw.get("pion_control") or {})
    kaon_stage_raw = (staged_amplitudes_raw.get("kaon_nosub") or {})
    pion_stage_scaled = (staged_amplitudes_scaled.get("pion_control") or {})
    kaon_stage_scaled = (staged_amplitudes_scaled.get("kaon_nosub") or {})
    pion_refined_final = (refined_amplitudes.get("pion_control") or {})
    kaon_refined_final = (refined_amplitudes.get("kaon_nosub") or {})
    pion_refined_pre = (
        pion_diagnostics.get("refined_amplitudes_pre_postrefine_scale")
        or pion_refined_final
        or {}
    )
    kaon_refined_pre = (
        kaon_diagnostics.get("refined_amplitudes_pre_postrefine_scale")
        or kaon_refined_final
        or {}
    )
    pion_postfit_scales = pion_diagnostics.get("postfit_component_scales") or {}
    kaon_postfit_scales = kaon_diagnostics.get("postfit_component_scales") or {}
    pion_postrefine_scales = pion_diagnostics.get("postrefine_component_scales") or {}
    kaon_postrefine_scales = kaon_diagnostics.get("postrefine_component_scales") or {}
    pion_manual_scaling_active = _component_scale_map_has_nonunity(pion_postfit_scales) or _component_scale_map_has_nonunity(pion_postrefine_scales)
    kaon_manual_scaling_active = _component_scale_map_has_nonunity(kaon_postfit_scales) or _component_scale_map_has_nonunity(kaon_postrefine_scales)
    kaon_sigma0_active = bool(
        kaon_diagnostics.get("sigma0_active")
        or component_fit_result.get("H_kaon_fit_k_sigma0_scaled") is not None
        or component_fit_result.get("H_kaon_fit_k_sigma0_scaled_refined") is not None
    )
    alignment_payload = component_fit_result.get("pion_component_alignment") or {}
    alignment_pdf_lines = []
    if alignment_payload:
        alignment_pdf_lines = [
            "alignment baseline score={}".format(
                _format_fit_number(alignment_payload.get("baseline_score"))
            ),
            "alignment proposed score={}".format(
                _format_fit_number(alignment_payload.get("proposed_score"))
            ),
            "alignment applied score={}".format(
                _format_fit_number(alignment_payload.get("applied_score"))
            ),
            "alignment applied source={}".format(
                alignment_payload.get("source") or "unknown"
            ),
        ]
        if not bool(alignment_payload.get("accepted", False)):
            alignment_pdf_lines.append("alignment proposal rejected; baseline map applied")

    if _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_pion_control_input"),
        "pion-control data",
        "{}pion-control staged SIMC component fit".format(title_prefix),
        [
            (component_fit_result.get("H_pion_fit_pi_n_scaled"), "pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_pion_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_pion_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pion_fit_k_sigma0_scaled"), "K-Sigma0", ROOT.kCyan + 2, 1),
            (component_fit_result.get("H_pion_fit_total"), "total fit", ROOT.kGreen + 2, 2),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            *alignment_pdf_lines,
            "status: staged baseline",
            "fit mode: {}".format(component_fit_result.get("fit_mode_pion") or component_fit_result.get("fit_mode") or "unknown"),
            "strategy: {}".format(
                _format_fit_strategy(pion_diagnostics)
            ),
            "solution: stage_only",
            "validation: {}".format("pass" if bool(pion_stage_validation.get("accepted")) else "fail"),
            "template MM shift={:.6f}".format(
                float(component_fit_result.get("template_mm_offset_data") or 0.0)
            ),
            "post-fit scales: {}".format(
                _format_component_scale_map(
                    pion_diagnostics.get("postfit_component_scales")
                )
            ),
            "manual scaling active: {} (pion-control denominator convention)".format(
                "yes" if pion_manual_scaling_active else "no"
            ),
            "B_n={}  B_delta={}  B_sidis={}".format(
                _format_fit_number(pion_stage_amplitudes.get("pi_n")),
                _format_fit_number(pion_stage_amplitudes.get("pi_delta")),
                _format_fit_number(pion_stage_amplitudes.get("pi_sidis")),
            ),
            "K-Sigma0 scale={}".format(
                _format_fit_number(pion_stage_amplitudes.get(KAON_SIGMA0_TEMPLATE_NAME))
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(pion_stage_validation.get("chi2_ndf")),
                _format_fit_metric(pion_stage_validation.get("fit_p_value")),
            ),
            "anchor windows: {}".format(
                _format_window_list(
                    pion_diagnostics.get("include_windows")
                )
            ),
            "excluded windows: {}".format(
                _format_excluded_window_list(
                    pion_diagnostics.get("exclude_windows")
                )
            ),
        ],
    ):
        _record_component_page(
            page_manifest,
            page_id_prefix,
            "pion_control_fit",
            scope=component_fit_result.get("analysis_scope") or "unknown",
            authoritative=authoritative,
        )

    kaon_signal_reference_hist, kaon_signal_reference_scale, kaon_signal_reference_source, kaon_signal_reference_normalization = (
        _resolve_kaon_lambda_reference_for_plot(
            component_fit_result,
            component_fit_result.get("H_kaon_nosub_input"),
            cut_window,
            component_fit_result.get("analysis_scope", "unknown"),
            "H_k_lambda_reference_fit",
        )
    )
    has_kaon_signal_reference = kaon_signal_reference_hist is not None
    has_sigma0_component = component_fit_result.get("H_kaon_fit_k_sigma0_scaled") is not None
    kaon_title = "{}kaon no-sub SIMC pion-background fit".format(title_prefix)
    if has_kaon_signal_reference:
        kaon_title = "{}kaon no-sub SIMC pion-background fit + K-Lambda gauge".format(title_prefix)
    if has_sigma0_component:
        kaon_title = "{} + K-Sigma0".format(kaon_title)

    kaon_overlay_specs = [
        (component_fit_result.get("H_kaon_fit_pi_n_scaled"), "pi-n", ROOT.kRed + 1, 1),
        (component_fit_result.get("H_kaon_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
        (component_fit_result.get("H_kaon_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
    ]
    if has_sigma0_component:
        kaon_overlay_specs.append(
            (component_fit_result.get("H_kaon_fit_k_sigma0_scaled"), "K-Sigma0", ROOT.kCyan + 2, 1)
        )
    if has_kaon_signal_reference:
        kaon_overlay_specs.append(
            (kaon_signal_reference_hist, "K-Lambda gauge", ROOT.kBlue + 1, 2)
        )
    kaon_overlay_specs.append(
        (component_fit_result.get("H_kaon_pion_bg_fit_total"), "pion-bg sum", ROOT.kOrange + 7, 2)
    )

    if _print_component_overlay_page(
        pdf_name,
        (
            component_fit_result.get("H_kaon_nosub_input")
            if not protected_kaon_mode else None
        ),
        "kaon no-sub data",
        kaon_title.replace("SIMC pion-background fit", "staged SIMC pion-background fit"),
        kaon_overlay_specs,
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: staged baseline",
            "fit mode: {}".format(component_fit_result.get("fit_mode_kaon") or component_fit_result.get("fit_mode") or "unknown"),
            "strategy: {}".format(
                _format_fit_strategy(kaon_diagnostics)
            ),
            "solution: stage_only",
            "validation: {}".format("pass" if bool(kaon_stage_validation.get("accepted")) else "fail"),
            "template MM shift={:.6f}".format(
                float(component_fit_result.get("template_mm_offset_data") or 0.0)
            ),
            "post-fit scales: {}".format(
                _format_component_scale_map(
                    kaon_diagnostics.get("postfit_component_scales")
                )
            ),
            "manual kaon-side scaling active: {}".format(
                "yes" if kaon_manual_scaling_active else "no"
            ),
            "sigma0 active: {}".format("yes" if kaon_sigma0_active else "no"),
            "A_n={}  A_delta={}  A_sidis={}".format(
                _format_fit_number(kaon_stage_amplitudes.get("pi_n")),
                _format_fit_number(kaon_stage_amplitudes.get("pi_delta")),
                _format_fit_number(kaon_stage_amplitudes.get("pi_sidis")),
            ),
            "K-Sigma0 scale={}".format(
                _format_fit_number(kaon_stage_amplitudes.get(KAON_SIGMA0_TEMPLATE_NAME))
            ) if has_sigma0_component else "K-Sigma0 scale=n/a",
            "K-Lambda comparison scale={}".format(
                _format_fit_number(kaon_signal_reference_scale)
            ) if has_kaon_signal_reference else "K-Lambda comparison scale=n/a",
            "K-Lambda comparison source={}".format(kaon_signal_reference_source),
            "K-Lambda comparison normalization={}".format(
                kaon_signal_reference_normalization
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(kaon_stage_validation.get("chi2_ndf")),
                _format_fit_metric(kaon_stage_validation.get("fit_p_value")),
            ),
            "anchor windows: {}".format(
                _format_window_list(
                    kaon_diagnostics.get("include_windows")
                )
            ),
            "excluded windows: {}".format(
                _format_excluded_window_list(
                    kaon_diagnostics.get("exclude_windows")
                )
            ),
        ],
        cut_window=cut_window,
    ):
        _record_component_page(
            page_manifest,
            page_id_prefix,
            "kaon_fit",
            scope=component_fit_result.get("analysis_scope") or "unknown",
            authoritative=authoritative,
        )

    if protected_render_state.get("render_success_pages"):
        _print_protected_pi_delta_pages(
            pdf_name,
            component_fit_result,
            protected_render_state,
            title_prefix=title_prefix,
            cut_window=cut_window,
        )
        _record_component_page(
            page_manifest,
            page_id_prefix,
            "protected_fit",
            scope=component_fit_result.get("analysis_scope") or "unknown",
            authoritative=authoritative,
        )
    elif protected_kaon_mode:
        _print_protected_pi_delta_status_page(
            pdf_name,
            component_fit_result,
            protected_render_state,
            title_prefix=title_prefix,
        )
        _record_component_page(
            page_manifest,
            page_id_prefix,
            "protected_status",
            scope=component_fit_result.get("analysis_scope") or "unknown",
            authoritative=authoritative,
        )

    # C.4 keeps staged/protected fit evidence in the main PDF while routing
    # optional staged-versus-refined residual detail to pion-fit debug.
    refinement_pdf = refinement_detail_pdf or pdf_name
    _print_joint_refinement_overlay_page(
        refinement_pdf,
        component_fit_result.get("H_pion_control_input"),
        component_fit_result.get("H_pion_fit_total"),
        component_fit_result.get("H_pion_fit_total_refined"),
        [
            (component_fit_result.get("H_pion_fit_pi_n_scaled_refined"), "refined pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_pion_fit_pi_sidis_scaled_refined"), "refined pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_pion_fit_pi_delta_scaled_refined"), "refined pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pion_fit_k_sigma0_scaled_refined"), "refined K-Sigma0", ROOT.kCyan + 2, 1),
        ],
        "{}pion-control staged vs refined component fit".format(title_prefix),
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            *alignment_pdf_lines,
            "fit mode: {}".format(component_fit_result.get("fit_mode_pion") or component_fit_result.get("fit_mode") or "unknown"),
            "joint status: {}".format(pion_diagnostics.get("joint_refinement_status") or "unknown"),
            "refined validation: {}".format("pass" if bool((pion_diagnostics.get("validation") or {}).get("accepted")) else "fail"),
            "manual scaling active: {} (pion-control kept as denominator)".format(
                "yes" if pion_manual_scaling_active else "no"
            ),
            "post-fit scales: {}".format(_format_component_scale_map(pion_postfit_scales)),
            "post-refine scales: {}".format(_format_component_scale_map(pion_postrefine_scales)),
            "staged chi2/ndf={}  refined chi2/ndf={}".format(
                _format_fit_metric(pion_stage_validation.get("chi2_ndf")),
                _format_fit_metric(pion_diagnostics.get("chi2_ndf")),
            ),
            "B_n raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(pion_stage_raw.get("pi_n")),
                _format_fit_number(pion_refined_pre.get("pi_n")),
                _format_fit_number(pion_refined_final.get("pi_n")),
            ),
            "B_delta raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(pion_stage_raw.get("pi_delta")),
                _format_fit_number(pion_refined_pre.get("pi_delta")),
                _format_fit_number(pion_refined_final.get("pi_delta")),
            ),
            "B_sidis raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(pion_stage_raw.get("pi_sidis")),
                _format_fit_number(pion_refined_pre.get("pi_sidis")),
                _format_fit_number(pion_refined_final.get("pi_sidis")),
            ),
            "B_n stage/refined = {} / {}".format(
                _format_fit_number(pion_stage_amplitudes.get("pi_n")),
                _format_fit_number(component_fit_result.get("B_n")),
            ),
            "B_delta stage/refined = {} / {}".format(
                _format_fit_number(pion_stage_amplitudes.get("pi_delta")),
                _format_fit_number(component_fit_result.get("B_delta")),
            ),
            "B_sidis stage/refined = {} / {}".format(
                _format_fit_number(pion_stage_amplitudes.get("pi_sidis")),
                _format_fit_number(component_fit_result.get("B_sidis")),
            ),
        ],
        cut_window=cut_window,
    )

    _print_joint_refinement_overlay_page(
        refinement_pdf,
        (
            component_fit_result.get("H_kaon_nosub_input")
            if not protected_kaon_mode else None
        ),
        component_fit_result.get("H_kaon_fit_total"),
        component_fit_result.get("H_kaon_fit_total_refined"),
        [
            (component_fit_result.get("H_kaon_fit_pi_n_scaled_refined"), "refined pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_kaon_fit_pi_sidis_scaled_refined"), "refined pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_kaon_fit_pi_delta_scaled_refined"), "refined pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_kaon_fit_k_sigma0_scaled_refined"), "refined K-Sigma0", ROOT.kCyan + 2, 1),
            (component_fit_result.get("H_kaon_fit_k_lambda_scaled_refined"), "refined K-Lambda", ROOT.kBlue + 1, 1),
            (component_fit_result.get("H_kaon_pion_bg_fit_total_refined"), "refined pion-bg sum", ROOT.kOrange + 7, 2),
        ],
        "{}kaon no-sub staged vs refined component fit".format(title_prefix),
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "pion alignment: imported from pion-control calibration",
            *alignment_pdf_lines,
            "fit mode: {}".format(component_fit_result.get("fit_mode_kaon") or component_fit_result.get("fit_mode") or "unknown"),
            "joint status: {}".format(kaon_diagnostics.get("joint_refinement_status") or "unknown"),
            "refined validation: {}".format("pass" if bool((kaon_diagnostics.get("validation") or {}).get("accepted")) else "fail"),
            "manual kaon-side scaling active: {}".format(
                "yes" if kaon_manual_scaling_active else "no"
            ),
            "sigma0 active: {}".format("yes" if kaon_sigma0_active else "no"),
            "post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
            "post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            "staged chi2/ndf={}  refined chi2/ndf={}".format(
                _format_fit_metric(kaon_stage_validation.get("chi2_ndf")),
                _format_fit_metric(kaon_diagnostics.get("chi2_ndf")),
            ),
            "A_n raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(kaon_stage_raw.get("pi_n")),
                _format_fit_number(kaon_refined_pre.get("pi_n")),
                _format_fit_number(kaon_refined_final.get("pi_n")),
            ),
            "A_delta raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(kaon_stage_raw.get("pi_delta")),
                _format_fit_number(kaon_refined_pre.get("pi_delta")),
                _format_fit_number(kaon_refined_final.get("pi_delta")),
            ),
            "A_sidis raw stage / refined(pre-scale) / final = {} / {} / {}".format(
                _format_fit_number(kaon_stage_raw.get("pi_sidis")),
                _format_fit_number(kaon_refined_pre.get("pi_sidis")),
                _format_fit_number(kaon_refined_final.get("pi_sidis")),
            ),
            "A_n stage/refined = {} / {}".format(
                _format_fit_number(kaon_stage_amplitudes.get("pi_n")),
                _format_fit_number(component_fit_result.get("A_n")),
            ),
            "A_delta stage/refined = {} / {}".format(
                _format_fit_number(kaon_stage_amplitudes.get("pi_delta")),
                _format_fit_number(component_fit_result.get("A_delta")),
            ),
            "A_sidis stage/refined = {} / {}".format(
                _format_fit_number(kaon_stage_amplitudes.get("pi_sidis")),
                _format_fit_number(component_fit_result.get("A_sidis")),
            ),
        ],
        cut_window=cut_window,
    )

    _print_kaon_pion_bg_comparison_page(
        pdf_name,
        (
            component_fit_result.get("H_kaon_nosub_input")
            if not protected_kaon_mode else None
        ),
        component_fit_result.get("H_kaon_pion_bg_fit_total"),
        component_fit_result.get("H_kaon_pion_bg_fit_total_refined_pre_postrefine"),
        component_fit_result.get("H_kaon_pion_bg_fit_total_refined"),
        component_fit_result.get("H_kaon_fit_k_sigma0_scaled_refined"),
        "{}kaon no-sub pion-background comparison with sigma0 stabilizer".format(title_prefix),
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "pion alignment: imported from pion-control calibration",
            "sigma0 active: {}".format("yes" if kaon_sigma0_active else "no"),
            "staged pion-bg integral={}".format(
                _format_fit_number(_hist_integral(component_fit_result.get("H_kaon_pion_bg_fit_total")))
            ),
            "refined pion-bg integral={}".format(
                _format_fit_number(_hist_integral(component_fit_result.get("H_kaon_pion_bg_fit_total_refined_pre_postrefine")))
            ),
            "final postrefine pion-bg integral={}".format(
                _format_fit_number(_hist_integral(component_fit_result.get("H_kaon_pion_bg_fit_total_refined")))
            ),
            "K-Sigma0 integral={}".format(
                _format_fit_number(_hist_integral(component_fit_result.get("H_kaon_fit_k_sigma0_scaled_refined")))
            ),
            "post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
            "post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
        ],
        cut_window=cut_window,
    )
    _print_residual_shift_template_pages(
        pdf_name,
        component_fit_result,
        title_prefix,
        protected_kaon_mode=protected_kaon_mode,
    )
    _print_residual_shift_scan_pages(
        pdf_name,
        component_fit_result,
        title_prefix,
        protected_kaon_mode=protected_kaon_mode,
    )

    if component_fit_result.get("analysis_scope") == "setting-wide":
        kaon_display_step_overlays = (
            _filtered_kaon_early_component_step_overlays(
                component_fit_result.get("H_kaon_fit_step_overlays")
            )
            if protected_kaon_mode
            else component_fit_result.get("H_kaon_fit_step_overlays")
        )
        _print_component_step_pages(
            pdf_name,
            component_fit_result.get("H_pion_control_input"),
            component_fit_result.get("H_pion_fit_step_overlays"),
            title_prefix,
            "pion-control",
        )
        _print_component_amplitude_pages(
            pdf_name,
            component_fit_result.get("H_pion_control_input"),
            component_fit_result.get("H_pion_fit_step_overlays"),
            title_prefix,
            "pion-control",
        )
        _print_component_step_pages(
            pdf_name,
            component_fit_result.get("H_kaon_nosub_input"),
            kaon_display_step_overlays,
            title_prefix,
            "kaon no-sub",
        )
        _print_component_amplitude_pages(
            pdf_name,
            component_fit_result.get("H_kaon_nosub_input"),
            kaon_display_step_overlays,
            title_prefix,
            "kaon no-sub",
        )


def serialize_particle_subtraction_component_result(result):
    if not isinstance(result, dict):
        return {}
    serializable = {}
    for key, value in result.items():
        if _is_root_hist(value):
            continue
        if key.startswith("H_"):
            continue
        cleaned_value = _json_ready_particle_subtraction_value(value)
        if cleaned_value is _JSON_SKIP:
            continue
        serializable[key] = cleaned_value
    return serializable
