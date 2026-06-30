#! /usr/bin/python

from __future__ import annotations

import itertools
import math
from copy import deepcopy

import numpy as np
import ROOT
from scipy.optimize import least_squares, lsq_linear
from scipy.stats import chi2 as chi2_dist

from background_config import (
    BG_OPT_MM_PLOT_MAX,
    BG_OPT_MM_PLOT_MIN,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    get_particle_subtraction_setting_key,
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
    resolve_particle_subtraction_component_fit_windows,
    resolve_particle_subtraction_mode,
)
from utility import normalize_hist_to_unit_area


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


def _is_root_hist(obj):
    try:
        return bool(obj is not None and obj.InheritsFrom("TH1"))
    except Exception:
        return False


_JSON_SKIP = object()


def _json_ready_particle_subtraction_value(value):
    if _is_root_hist(value):
        return _JSON_SKIP
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
    return deepcopy(value)


def _is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def _clone_hist(template_hist, name, title=None, reset=False):
    if template_hist is None:
        return None
    cloned = template_hist.Clone(name)
    cloned.SetDirectory(0)
    if title is not None:
        cloned.SetTitle(title)
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


def _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason):
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
    design_matrix = np.column_stack(
        [np.asarray(fit_inputs["template_columns"][name], dtype=float) for name in fit_names]
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
    staged_scaled_hist_map.update(deepcopy(result.get("extra_scaled_hists") or {}))

    def _retain_staged_solution(joint_refinement_summary, message_override=None):
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
        diagnostics["validation"] = deepcopy(stage_validation)
        diagnostics["success"] = bool(stage_validation.get("accepted"))
        diagnostics["message"] = str(
            message_override
            or (joint_refinement_summary or {}).get("message")
            or ""
        )
        diagnostics["chi2"] = stage_validation.get("chi2")
        diagnostics["ndf"] = stage_validation.get("ndf")
        diagnostics["chi2_ndf"] = stage_validation.get("chi2_ndf")
        diagnostics["fit_p_value"] = stage_validation.get("fit_p_value")
        diagnostics["n_fit_bins"] = (joint_refinement_summary or {}).get("active_bin_count", 0)
        diagnostics["fallback_used"] = not bool(diagnostics["success"])
        diagnostics["fallback_reason"] = (
            ""
            if bool(diagnostics["success"]) else (
                message_override
                or "staged validation rejected: {}".format(
                    "; ".join(stage_validation.get("rejection_reasons") or ["unknown"])
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
        and int(refinement_result.get("active_bin_count") or 0) <= 0
    ):
        skipped_summary = deepcopy(refinement_result)
        skipped_summary["status"] = "skipped_no_valid_bins"
        skipped_summary["success"] = False
        skipped_summary["message"] = (
            "joint refinement skipped: no valid full-range fit bins; staged solution retained"
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
    result["refined_residual_hist"] = refined_residual_hist
    result["refined_scaled_hist_map"] = refined_scaled_hist_map
    if amplitude_prefix == "A":
        result["refined_pion_bg_fit_hist"] = refined_pion_bg_fit_hist
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

    weighted_design = np.column_stack(
        [fit_inputs["template_columns"][name] / fit_inputs["sigma"] for name in fit_names]
    )
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
        )

    fallback_reason = _validate_component_shapes(
        component_hists,
        target_hist,
        component_names=fitted_names,
    )
    if fallback_reason:
        return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)

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
        return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)

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
        )

    if not getattr(fit_result, "success", False):
        return _zero_fit_result(
            target_hist,
            amplitude_prefix,
            context,
            "fit solver failed: {}".format(getattr(fit_result, "message", "unknown")),
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
            ),
            fit_min=fit_min,
            cleanup_validation_mm_max=validation_options.get("cleanup_validation_mm_max"),
            exclude_windows=excluded_windows,
            context="{}_kaon_nosub".format(context or "scope"),
        )
        if shift_selection is not None:
            return _annotate_fit_result_with_residual_shift_payload(
                shift_selection.get("fit_result"),
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
        "refined_pion_bg_fit_hist": result.get("refined_pion_bg_fit_hist") or result["pion_bg_fit_hist"],
        "refined_residual_hist": result.get("refined_residual_hist") or result["residual_hist"],
        "refined_pi_n_scaled_hist": refined_scaled_hist_map.get("pi_n") or result["pi_n_scaled_hist"],
        "refined_pi_delta_scaled_hist": refined_scaled_hist_map.get("pi_delta") or result["pi_delta_scaled_hist"],
        "refined_pi_sidis_scaled_hist": refined_scaled_hist_map.get("pi_sidis") or result["pi_sidis_scaled_hist"],
        "refined_k_lambda_scaled_hist": refined_scaled_hist_map.get(KAON_SIGNAL_TEMPLATE_NAME) or signal_scaled_hist,
        "refined_k_sigma0_scaled_hist": refined_scaled_hist_map.get(KAON_SIGMA0_TEMPLATE_NAME) or sigma0_scaled_hist,
        "k_lambda_reference_hist": signal_reference_hist,
        "step_overlays": result.get("step_overlays") or [],
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
            "joint_refinement_amplitude_floor": amplitude_floor,
            "template_corr_warn": template_corr_warn,
            "cleanup_validation_mm_max": validation_options.get("cleanup_validation_mm_max"),
        },
    }
    if excluded_windows:
        for key in (
            "fit_hist",
            "pion_bg_fit_hist",
            "residual_hist",
            "pi_n_scaled_hist",
            "pi_delta_scaled_hist",
            "pi_sidis_scaled_hist",
            "k_lambda_scaled_hist",
            "k_sigma0_scaled_hist",
            "refined_fit_hist",
            "refined_pion_bg_fit_hist",
            "refined_residual_hist",
            "refined_pi_n_scaled_hist",
            "refined_pi_delta_scaled_hist",
            "refined_pi_sidis_scaled_hist",
            "refined_k_lambda_scaled_hist",
            "refined_k_sigma0_scaled_hist",
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
):
    mode = resolve_particle_subtraction_mode(inpDict)
    if mode != PARTICLE_SUBTRACTION_MODE_COMPONENTS:
        return {
            "particle_subtraction_mode": mode,
            "analysis_scope": analysis_scope,
            "fallback_used": True,
            "fallback_reason": "particle subtraction mode is not simc_shape_components",
        }

    template_mm_offset_data = float(mm_offset_data) if _is_finite_number(mm_offset_data) else 0.0
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
    resolved_phi_setting = phi_setting
    if resolved_phi_setting is None and isinstance(inpDict, dict):
        resolved_phi_setting = inpDict.get("phi_setting")
    setting_key = get_particle_subtraction_setting_key(inpDict)

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

    result = {
        "particle_subtraction_mode": mode,
        "analysis_scope": analysis_scope,
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
            (pion_fit.get("template_hists") or {}).get("pi_n"),
            "H_simc_shape_pi_n_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_delta": _clone_hist(
            (pion_fit.get("template_hists") or {}).get("pi_delta"),
            "H_simc_shape_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_sidis": _clone_hist(
            (pion_fit.get("template_hists") or {}).get("pi_sidis"),
            "H_simc_shape_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_k_lambda": _clone_hist(
            aligned_kaon_signal_shape,
            "H_simc_shape_k_lambda_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_k_sigma0": _clone_hist(
            aligned_kaon_sigma0_shape,
            "H_simc_shape_k_sigma0_{}".format(context or analysis_scope),
        ),
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
        "H_kaon_pion_bg_fit_total": kaon_fit["pion_bg_fit_hist"],
        "H_kaon_fit_pi_n_scaled_refined": kaon_fit.get("refined_pi_n_scaled_hist"),
        "H_kaon_fit_pi_delta_scaled_refined": kaon_fit.get("refined_pi_delta_scaled_hist"),
        "H_kaon_fit_pi_sidis_scaled_refined": kaon_fit.get("refined_pi_sidis_scaled_hist"),
        "H_kaon_fit_k_lambda_scaled_refined": kaon_fit.get("refined_k_lambda_scaled_hist"),
        "H_kaon_fit_k_sigma0_scaled_refined": kaon_fit.get("refined_k_sigma0_scaled_hist"),
        "H_kaon_fit_total_refined": kaon_fit.get("refined_fit_hist"),
        "H_kaon_pion_bg_fit_total_refined": kaon_fit.get("refined_pion_bg_fit_hist"),
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
        return

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
        return

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


def _print_component_application_status_page(
    pdf_name,
    component_payload,
    title_prefix="",
):
    if not isinstance(component_payload, dict):
        return

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
    return (
        "{}: seen={} passed={} mm_passed={} norm={} fallback={}".format(
            component_label,
            int(diagnostics.get("n_events_seen", 0) or 0),
            int(diagnostics.get("n_events_passed", 0) or 0),
            int(diagnostics.get("n_events_passed_mm_window", 0) or 0),
            bool(diagnostics.get("normalized", False)),
            bool(diagnostics.get("fallback_used", False)),
        )
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


def _print_residual_shift_template_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
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
            original_clone.SetTitle(
                "{}{} residual-shift template: {}".format(
                    title_prefix,
                    sample_label,
                    component_label,
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
        objective_graph.SetTitle(
            "{}{} residual-shift scan: {}".format(
                title_prefix,
                sample_label,
                _component_plot_label(component_name),
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
        _print_component_overlay_page(
            pdf_name,
            full_base,
            full_label,
            "{}SIMC component templates (full MM range)".format(title_prefix),
            full_overlays,
            full_stats,
            cut_window=cut_window,
        )

    cut_base, cut_label, cut_overlays, cut_stats = _build_page("setting_shape")
    if cut_base is not None:
        _print_component_overlay_page(
            pdf_name,
            cut_base,
            cut_label,
            "{}SIMC component templates (analysis MM-cut window)".format(title_prefix),
            cut_overlays,
            cut_stats,
            cut_window=cut_window,
        )


def print_particle_subtraction_component_application_pages(
    pdf_name,
    component_payload,
    title_prefix="",
    cut_window=None,
):
    if not isinstance(component_payload, dict):
        return
    if not bool(component_payload.get("accepted")):
        _print_component_application_status_page(
            pdf_name,
            component_payload,
            title_prefix=title_prefix,
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
    closure_label = "refined" if joint_mode_active else "staged"
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

    _print_single_hist_page(
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
            "numerator: kaon pion-bg model{}".format(
                " after kaon-side scaling" if kaon_manual_scaling_active else ""
            ),
        ],
        cut_window=cut_window,
        line_color=ROOT.kViolet + 1,
    )

    if joint_mode_active and component_payload.get("H_pion_weight_vs_MM_stage") is not None:
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

    _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        "kaon pion-bg model",
        "{}Part 3 {} model closure: weighted pion model vs kaon bg model".format(
            title_prefix,
            closure_label,
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
            "kaon pion model integral={}".format(
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
            "kaon bg curve shown{}.".format(
                " after kaon-side post-refine scaling" if kaon_manual_scaling_active else " without manual kaon-side scaling"
            ),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
            "pion-control post-refine scales: {}".format(_format_component_scale_map(pion_postrefine_scales)),
        ],
        cut_window=cut_window,
    )

    _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        "kaon pion-bg model",
        "{}Part 3 {} event-template closure vs kaon bg model".format(
            title_prefix,
            closure_label,
        ),
        [
            (component_payload.get("H_pion_subtraction_template_MM_nosub"), "weighted pion template (full)", ROOT.kOrange + 7, 2),
        ],
        [
            "scope: {}".format(scope_label),
            "signature match={}".format(
                "pass" if bool(event_template_closure.get("signature_match")) else "fail"
            ),
            "kaon pion model integral={}".format(
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
    )

    _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_MM_nosub_before_pion_subtraction"),
        "kaon data before pion subtraction",
        "{}Part 3 kaon data vs pion-background models".format(title_prefix),
        [
            (component_payload.get("H_kaon_pion_model"), "kaon pion-bg model", ROOT.kBlue + 1, 1),
            (component_payload.get("H_pion_subtraction_template_MM_nosub"), "weighted pion template", ROOT.kOrange + 7, 2),
        ],
        [
            "scope: {}".format(scope_label),
            "before full integral={}".format(
                _format_fit_number(component_payload.get("kaon_integral_before_pion_sub_full"))
            ),
            "kaon pion model integral={}".format(
                _format_fit_number(diagnostics.get("kaon_pion_model_integral"))
            ),
            "weighted template full integral={}".format(
                _format_fit_number(component_payload.get("weighted_pion_integral_full"))
            ),
            "effective scale={}".format(
                _format_fit_number(component_payload.get("particle_subtraction_effective_scale"))
            ),
            "kaon pion-bg model shown{}.".format(
                " after kaon-side post-refine scaling" if kaon_manual_scaling_active else " without manual kaon-side scaling"
            ),
            "kaon-side post-fit scales: {}".format(_format_component_scale_map(kaon_postfit_scales)),
            "kaon-side post-refine scales: {}".format(_format_component_scale_map(kaon_postrefine_scales)),
        ],
        cut_window=cut_window,
    )

    if joint_mode_active and component_payload.get("H_kaon_pion_model_stage") is not None:
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

    _print_component_overlay_page(
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
    )


def print_particle_subtraction_component_fit_pages(
    pdf_name,
    component_fit_result,
    title_prefix="",
    cut_window=None,
):
    if not isinstance(component_fit_result, dict):
        return

    title_prefix = (title_prefix or "").strip()
    if title_prefix:
        title_prefix = "{} ".format(title_prefix)

    pion_diagnostics = ((component_fit_result.get("diagnostics") or {}).get("pion") or {})
    kaon_diagnostics = ((component_fit_result.get("diagnostics") or {}).get("kaon") or {})
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

    _print_component_overlay_page(
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
    )

    has_kaon_signal_reference = component_fit_result.get("H_kaon_fit_k_lambda_reference") is not None
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
            (component_fit_result.get("H_kaon_fit_k_lambda_reference"), "K-Lambda gauge", ROOT.kBlue + 1, 2)
        )
    kaon_overlay_specs.append(
        (component_fit_result.get("H_kaon_pion_bg_fit_total"), "pion-bg sum", ROOT.kOrange + 7, 2)
    )

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_kaon_nosub_input"),
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
            "A_n={}  A_delta={}  A_sidis={}".format(
                _format_fit_number(kaon_stage_amplitudes.get("pi_n")),
                _format_fit_number(kaon_stage_amplitudes.get("pi_delta")),
                _format_fit_number(kaon_stage_amplitudes.get("pi_sidis")),
            ),
            "K-Sigma0 scale={}".format(
                _format_fit_number(kaon_stage_amplitudes.get(KAON_SIGMA0_TEMPLATE_NAME))
            ) if has_sigma0_component else "K-Sigma0 scale=n/a",
            "K-Lambda gauge scale={}".format(
                _format_fit_number(kaon_stage_amplitudes.get(KAON_SIGNAL_TEMPLATE_NAME))
            ) if has_kaon_signal_reference else "K-Lambda gauge scale=n/a",
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
    )

    _print_joint_refinement_overlay_page(
        pdf_name,
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
        pdf_name,
        component_fit_result.get("H_kaon_nosub_input"),
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
            "fit mode: {}".format(component_fit_result.get("fit_mode_kaon") or component_fit_result.get("fit_mode") or "unknown"),
            "joint status: {}".format(kaon_diagnostics.get("joint_refinement_status") or "unknown"),
            "refined validation: {}".format("pass" if bool((kaon_diagnostics.get("validation") or {}).get("accepted")) else "fail"),
            "manual kaon-side scaling active: {}".format(
                "yes" if kaon_manual_scaling_active else "no"
            ),
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
    _print_residual_shift_template_pages(
        pdf_name,
        component_fit_result,
        title_prefix,
    )
    _print_residual_shift_scan_pages(
        pdf_name,
        component_fit_result,
        title_prefix,
    )

    if component_fit_result.get("analysis_scope") == "setting-wide":
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
            component_fit_result.get("H_kaon_fit_step_overlays"),
            title_prefix,
            "kaon no-sub",
        )
        _print_component_amplitude_pages(
            pdf_name,
            component_fit_result.get("H_kaon_nosub_input"),
            component_fit_result.get("H_kaon_fit_step_overlays"),
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
