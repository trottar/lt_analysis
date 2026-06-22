#! /usr/bin/python

from __future__ import annotations

import math
from copy import deepcopy

import numpy as np
import ROOT
from scipy.optimize import lsq_linear
from scipy.stats import chi2 as chi2_dist

from background_config import (
    BG_OPT_MM_PLOT_MAX,
    BG_OPT_MM_PLOT_MIN,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    resolve_particle_subtraction_component_postfit_scales,
    resolve_particle_subtraction_component_stage_amplitude_windows,
    resolve_particle_subtraction_component_fit_excluded_windows,
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


def _build_mm_shifted_hist(template_hist, shift, name, renormalize=False):
    if template_hist is None:
        return None
    shift_value = float(shift) if _is_finite_number(shift) else 0.0
    shifted_hist = _clone_hist(template_hist, name, reset=True)
    if shifted_hist is None:
        return None
    if abs(shift_value) <= 1e-12:
        shifted_hist.Add(template_hist)
        return shifted_hist

    x_axis = shifted_hist.GetXaxis()
    x_min = float(x_axis.GetXmin())
    x_max = float(x_axis.GetXmax())
    for bin_index in range(1, template_hist.GetNbinsX() + 1):
        content = float(template_hist.GetBinContent(bin_index))
        error = float(template_hist.GetBinError(bin_index))
        if content == 0.0 and error == 0.0:
            continue

        shifted_x = float(template_hist.GetBinCenter(bin_index)) + shift_value
        if shifted_x < x_min or shifted_x >= x_max:
            continue

        target_bin = int(x_axis.FindBin(shifted_x))
        shifted_hist.SetBinContent(
            target_bin,
            float(shifted_hist.GetBinContent(target_bin)) + content,
        )
        shifted_hist.SetBinError(
            target_bin,
            math.sqrt(float(shifted_hist.GetBinError(target_bin)) ** 2 + error ** 2),
        )

    if renormalize:
        normalize_hist_to_unit_area(
            shifted_hist,
            quiet=True,
            context=name,
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


def _validate_component_shapes(component_hists, target_hist):
    if target_hist is None:
        return "missing target histogram"
    if any(component_hists.get(name) is None for name in COMPONENT_NAMES):
        return "missing SIMC component shape"

    for component_name in COMPONENT_NAMES:
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
            }
        target_sum = float(np.sum(amplitude_inputs["y"]))
        amplitude = max(target_sum / template_sum, 0.0)
        sigma_value = math.sqrt(float(np.sum(np.square(amplitude_inputs["sigma"])))) / template_sum
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
            }

        weighted_target = fit_inputs["y"] / fit_inputs["sigma"]
        amplitude = max(float(np.dot(weighted_template, weighted_target) / denominator), 0.0)
        sigma_value = 1.0 / math.sqrt(denominator)

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
        component_name: float((component_scale_map or {}).get(component_name, 1.0) or 1.0)
        for component_name in COMPONENT_NAMES
    }
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

    resolved_scale_map = {}
    for component_name in COMPONENT_NAMES:
        scale_value = float((postfit_scale_map or {}).get(component_name, 1.0) or 1.0)
        resolved_scale_map[component_name] = scale_value
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
    full_amplitude_map = dict(extra_component_amplitudes)
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
        for template_name in extra_component_amplitudes
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
    accepted_uncertainties = deepcopy(diagnostics.get("accepted_component_uncertainties") or {})
    for component_name in COMPONENT_NAMES:
        if accepted_uncertainties.get(component_name) is None:
            continue
        accepted_uncertainties[component_name] = (
            float(accepted_uncertainties[component_name]) * resolved_scale_map[component_name]
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
):
    if exclude_windows is None:
        exclude_windows = []
    quality = _compute_fit_quality(
        target_hist,
        fit_hist,
        fit_min,
        fit_max,
        exclude_windows=exclude_windows,
        n_parameters=n_parameters,
    )

    oversub_bin_centers = []
    total_bins = 0
    sigma_tolerance = max(float(oversub_sigma_tolerance or 0.0), 0.0)
    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < fit_min or x_center > fit_max:
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
    accepted = True
    rejection_reasons = []
    if _is_finite_number(max_oversub_bin_count) and oversub_bin_count > int(max_oversub_bin_count):
        accepted = False
        rejection_reasons.append(
            "oversub_bin_count {} > {}".format(oversub_bin_count, int(max_oversub_bin_count))
        )
    if (
        _is_finite_number(max_oversub_bin_fraction)
        and oversub_bin_fraction > float(max_oversub_bin_fraction)
    ):
        accepted = False
        rejection_reasons.append(
            "oversub_bin_fraction {:.3f} > {:.3f}".format(
                oversub_bin_fraction,
                float(max_oversub_bin_fraction),
            )
        )
    if (
        _is_finite_number(max_full_range_chi2_ndf)
        and _is_finite_number(quality.get("chi2_ndf"))
        and float(quality["chi2_ndf"]) > float(max_full_range_chi2_ndf)
    ):
        accepted = False
        rejection_reasons.append(
            "chi2_ndf {:.3f} > {:.3f}".format(
                float(quality["chi2_ndf"]),
                float(max_full_range_chi2_ndf),
            )
        )

    return {
        "accepted": bool(accepted),
        "rejection_reasons": rejection_reasons,
        "oversub_bin_count": oversub_bin_count,
        "oversub_bin_fraction": oversub_bin_fraction,
        "oversub_sigma_tolerance": sigma_tolerance,
        "oversub_mm_range": (
            [float(min(oversub_bin_centers)), float(max(oversub_bin_centers))]
            if oversub_bin_centers else []
        ),
        **quality,
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


def _run_staged_component_pass(
    target_hist,
    template_hists,
    fit_names,
    fit_min,
    fit_max,
    anchor_window_map,
    stage_amplitude_window_map,
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

        solve_result = _solve_nonnegative_template_amplitude(
            residual_hist,
            template_hists[component_name],
            fit_min,
            fit_max,
            include_windows=anchor_window_map.get(component_name),
            exclude_windows=exclude_windows,
            amplitude_windows=stage_amplitude_window_map.get(component_name),
            amplitude_mode=(
                "window_integral"
                if (stage_amplitude_window_map.get(component_name) or [])
                else "least_squares"
            ),
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
                "anchor_windows": deepcopy(anchor_window_map.get(component_name) or []),
                "amplitude_mode": solve_result.get("amplitude_mode", "least_squares"),
                "amplitude_windows": deepcopy(solve_result.get("amplitude_windows") or []),
                "excluded_windows": deepcopy(exclude_windows or []),
                "H_baseline_before": baseline_before_hist,
                "H_residual_input": residual_hist,
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


def _fit_staged_anchor_templates(
    target_hist,
    component_hists,
    amplitude_prefix,
    fit_min,
    fit_max,
    anchor_windows,
    fit_order,
    stage_amplitude_windows=None,
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

    fallback_reason = _validate_component_shapes(component_hists, target_hist)
    if fallback_reason:
        return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)

    anchor_window_map = _coerce_window_map(anchor_windows)
    stage_amplitude_window_map = _coerce_window_map(stage_amplitude_windows)
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
    }

    stage_fit_hist = _build_model_hist(
        target_hist,
        {name: template_hists[name] for name in fitted_names},
        stage_amplitudes,
        "{}_stage_fit_hist_{}".format(amplitude_prefix, context),
    )
    stage_validation = _evaluate_model_validation(
        target_hist,
        stage_fit_hist,
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

    accepted_solution = "stage_only"
    accepted_message = "ordered residual staged fit"
    accepted_amplitudes = deepcopy(stage_amplitudes)
    accepted_uncertainties = deepcopy(stage_uncertainties)
    accepted_validation = deepcopy(stage_validation)

    joint_diagnostics = {
        "enabled": bool(joint_refinement_enabled),
        "success": False,
        "message": "joint refinement disabled",
        "n_fit_bins": 0,
        "prior_sigmas": deepcopy(prior_sigmas),
        "coefficients": {},
        "uncertainties": {},
        "validation": {},
    }
    coordinate_diagnostics = {
        "attempted": False,
        "success": False,
        "message": "coordinate refinement not attempted",
        "n_fit_bins": 0,
        "cycles_run": 0,
        "history": [],
        "coefficients": {},
        "validation": {},
    }

    if joint_refinement_enabled:
        joint_result = _solve_joint_template_amplitudes(
            target_hist,
            {name: template_hists[name] for name in fitted_names},
            fitted_names,
            fit_min,
            fit_max,
            stage_amplitudes,
            prior_sigmas,
            exclude_windows=exclude_windows,
        )
        joint_diagnostics = {
            "enabled": True,
            "success": bool(joint_result.get("success", False)),
            "message": joint_result.get("message", ""),
            "n_fit_bins": joint_result.get("n_fit_bins"),
            "prior_sigmas": deepcopy(joint_result.get("prior_sigmas") or prior_sigmas),
            "coefficients": deepcopy(joint_result.get("coefficients") or {}),
            "uncertainties": deepcopy(joint_result.get("uncertainties") or {}),
            "validation": {},
        }
        if joint_result.get("coefficients"):
            joint_fit_hist = _build_model_hist(
                target_hist,
                {name: template_hists[name] for name in fitted_names},
                joint_result["coefficients"],
                "{}_joint_fit_hist_{}".format(amplitude_prefix, context),
            )
            joint_validation = _evaluate_model_validation(
                        target_hist,
                        joint_fit_hist,
                        fit_min,
                        fit_max,
                        n_parameters=len(fitted_names),
                        exclude_windows=exclude_windows,
                        **validation_kwargs
                    )
            joint_diagnostics["validation"] = deepcopy(joint_validation)
            if bool(joint_result.get("success", False)) and joint_validation["accepted"]:
                accepted_solution = "joint_prior"
                accepted_message = "ordered residual staged fit with joint prior refinement"
                accepted_amplitudes.update(joint_result["coefficients"])
                accepted_uncertainties.update(joint_result.get("uncertainties") or {})
                accepted_validation = deepcopy(joint_validation)
            else:
                coordinate_result = _run_coordinate_template_updates(
                    target_hist,
                    {name: template_hists[name] for name in fitted_names},
                    fitted_names,
                    fit_min,
                    fit_max,
                    joint_result.get("coefficients") or stage_amplitudes,
                    prior_targets=stage_amplitudes,
                    prior_sigmas=joint_diagnostics["prior_sigmas"],
                    max_cycles=max_fit_cycles,
                    tolerance=fit_tolerance,
                    exclude_windows=exclude_windows,
                )
                coordinate_diagnostics = {
                    "attempted": True,
                    "success": bool(coordinate_result.get("success", False)),
                    "message": coordinate_result.get("message", ""),
                    "n_fit_bins": coordinate_result.get("n_fit_bins"),
                    "cycles_run": coordinate_result.get("cycles_run", 0),
                    "history": deepcopy(coordinate_result.get("history") or []),
                    "coefficients": deepcopy(coordinate_result.get("coefficients") or {}),
                    "validation": {},
                }
                if coordinate_result.get("coefficients"):
                    coordinate_fit_hist = _build_model_hist(
                        target_hist,
                        {name: template_hists[name] for name in fitted_names},
                        coordinate_result["coefficients"],
                        "{}_coordinate_fit_hist_{}".format(amplitude_prefix, context),
                    )
                    coordinate_validation = _evaluate_model_validation(
                        target_hist,
                        coordinate_fit_hist,
                        fit_min,
                        fit_max,
                        n_parameters=len(fitted_names),
                        exclude_windows=exclude_windows,
                        **validation_kwargs
                    )
                    coordinate_diagnostics["validation"] = deepcopy(coordinate_validation)
                    if coordinate_validation["accepted"]:
                        accepted_solution = "coordinate_prior"
                        accepted_message = "ordered residual staged fit with coordinate refinement"
                        accepted_amplitudes.update(coordinate_result["coefficients"])
                        accepted_validation = deepcopy(coordinate_validation)

    accepted_template_hists = {name: template_hists[name] for name in fitted_names}
    fit_hist = _build_model_hist(
        target_hist,
        accepted_template_hists,
        accepted_amplitudes,
        "{}_fit_hist_{}".format(amplitude_prefix, context),
    )
    residual_hist = _clone_hist(
        target_hist,
        "{}_fit_residual_{}".format(amplitude_prefix, context),
    )
    residual_hist.Add(fit_hist, -1.0)
    scaled_hist_map = _build_scaled_hist_map(
        accepted_template_hists,
        accepted_amplitudes,
        amplitude_prefix,
        context,
    )
    extra_scaled_hists = {
        template_name: scaled_hist_map.get(template_name)
        for template_name in extra_positive_templates
    }
    pion_bg_fit_hist = _build_model_hist(
        target_hist,
        {
            component_name: accepted_template_hists.get(component_name)
            for component_name in COMPONENT_NAMES
        },
        accepted_amplitudes,
        "{}_pion_bg_fit_{}".format(amplitude_prefix, context),
    ) if amplitude_prefix == "A" else None

    fit_status = "success"
    fallback_used = False
    fallback_reason = ""
    if not accepted_validation["accepted"]:
        fit_status = "fallback"
        fallback_used = True
        fallback_reason = (
            "no ordered-fit candidate passed full-range validation: {}".format(
                "; ".join(accepted_validation.get("rejection_reasons") or ["unknown"])
            )
        )

    extra_component_amplitudes = {
        template_name: float(accepted_amplitudes.get(template_name, 0.0) or 0.0)
        for template_name in extra_positive_templates
    }
    diagnostics = {
        "success": (not fallback_used),
        "status_code": None,
        "message": accepted_message,
        "chi2": accepted_validation.get("chi2"),
        "ndf": accepted_validation.get("ndf"),
        "chi2_ndf": accepted_validation.get("chi2_ndf"),
        "fit_p_value": accepted_validation.get("fit_p_value"),
        "n_fit_bins": accepted_validation.get("n_fit_bins"),
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "include_windows": deepcopy(include_windows),
        "exclude_windows": deepcopy(exclude_windows),
        "fallback_used": bool(fallback_used),
        "fallback_reason": fallback_reason,
        "fit_strategy": "ordered_residual",
        "accepted_solution": accepted_solution,
        "n_passes": total_passes,
        "fit_order": list(fitted_names),
        "anchor_windows": deepcopy(combined_window_map),
        "stage_amplitude_windows": deepcopy(stage_amplitude_window_map),
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
                None if accepted_uncertainties.get(component_name) is None
                else float(accepted_uncertainties.get(component_name))
            )
            for component_name in fitted_names
        },
        "prior_sigmas": deepcopy(prior_sigmas),
        "joint_refinement": deepcopy(joint_diagnostics),
        "coordinate_refinement": deepcopy(coordinate_diagnostics),
        "validation": deepcopy(accepted_validation),
        "component_amplitudes": {
            component_name: float(accepted_amplitudes.get(component_name, 0.0) or 0.0)
            for component_name in COMPONENT_NAMES
        },
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
    }

    result = {
        "{}_n".format(amplitude_prefix): float(accepted_amplitudes.get("pi_n", 0.0) or 0.0),
        "{}_delta".format(amplitude_prefix): float(accepted_amplitudes.get("pi_delta", 0.0) or 0.0),
        "{}_sidis".format(amplitude_prefix): float(accepted_amplitudes.get("pi_sidis", 0.0) or 0.0),
        "fit_status": fit_status,
        "diagnostics": diagnostics,
        "fit_hist": fit_hist,
        "residual_hist": residual_hist,
        "template_hists": accepted_template_hists,
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
    context="",
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    fit_config = get_particle_subtraction_component_fit_window_config("pion_control") or {}
    resolved_windows = resolve_particle_subtraction_component_fit_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
    )
    stage_amplitude_windows = resolve_particle_subtraction_component_stage_amplitude_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
    )
    excluded_windows = resolve_particle_subtraction_component_fit_excluded_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
    )
    anchor_windows = {
        component_name: [window]
        for component_name, window in resolved_windows.items()
    }
    prior_scale_map = {
        "pi_n": float(fit_config.get("particle_subtraction_prior_scale_pi_n", 1.0) or 1.0),
        "pi_delta": float(fit_config.get("particle_subtraction_prior_scale_pi_delta", 1.5) or 1.5),
        "pi_sidis": float(fit_config.get("particle_subtraction_prior_scale_pi_sidis", 2.0) or 2.0),
    }
    postfit_scale_map = resolve_particle_subtraction_component_postfit_scales("pion_control")
    extra_positive_templates = {}
    extra_anchor_windows = {}
    include_sigma0_signal_template = bool(fit_config.get("include_sigma0_signal_template", False))
    if include_sigma0_signal_template and _hist_has_usable_support(h_kaon_sigma0_shape):
        extra_positive_templates[KAON_SIGMA0_TEMPLATE_NAME] = h_kaon_sigma0_shape
        prior_scale_map[KAON_SIGMA0_TEMPLATE_NAME] = float(
            fit_config.get("particle_subtraction_prior_scale_k_sigma0_signal", 1.0) or 1.0
        )
        sigma0_anchor_window = fit_config.get("sigma0_signal_anchor_window")
        if sigma0_anchor_window is not None and len(sigma0_anchor_window) == 2:
            sigma0_window_min = float(sigma0_anchor_window[0])
            sigma0_window_max = float(sigma0_anchor_window[1])
            if bool(fit_config.get("apply_mm_offset_data", False)):
                sigma0_window_min += float(mm_offset_data)
                sigma0_window_max += float(mm_offset_data)
            extra_anchor_windows[KAON_SIGMA0_TEMPLATE_NAME] = [
                (sigma0_window_min, sigma0_window_max)
            ]
    validation_options = {
        "oversub_sigma_tolerance": fit_config.get("oversub_sigma_tolerance", 2.0),
        "max_oversub_bin_count": fit_config.get("max_oversub_bin_count"),
        "max_oversub_bin_fraction": fit_config.get("max_oversub_bin_fraction"),
        "max_full_range_chi2_ndf": fit_config.get("max_full_range_chi2_ndf"),
    }
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
    sigma0_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    sigma0_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    return {
        "B_n": result["B_n"],
        "B_delta": result["B_delta"],
        "B_sidis": result["B_sidis"],
        "B_sigma0": None if sigma0_amplitude is None else float(sigma0_amplitude),
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "template_hists": result.get("template_hists") or {},
        "fit_hist": result["fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
        "k_sigma0_scaled_hist": sigma0_scaled_hist,
        "step_overlays": result.get("step_overlays") or [],
    }


def fit_kaon_nosub_with_simc_pion_shapes(
    h_kaon_nosub,
    h_pi_n_shape,
    h_pi_delta_shape,
    h_pi_sidis_shape,
    h_kaon_signal_shape,
    h_kaon_sigma0_shape,
    inpDict,
    mm_offset_data=0.0,
    context="",
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    fit_config = get_particle_subtraction_component_fit_window_config("kaon_nosub") or {}
    resolved_windows = resolve_particle_subtraction_component_fit_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
    )
    stage_amplitude_windows = resolve_particle_subtraction_component_stage_amplitude_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
    )
    excluded_windows = resolve_particle_subtraction_component_fit_excluded_windows(
        "kaon_nosub",
        mm_offset_data=mm_offset_data,
    )
    anchor_windows = {
        component_name: [window]
        for component_name, window in resolved_windows.items()
    }
    mm_min = float(inpDict.get("mm_min", fit_min))
    mm_max = float(inpDict.get("mm_max", fit_max))
    extra_positive_templates = {}
    extra_anchor_windows = {}
    prior_scale_map = {
        "pi_n": float(fit_config.get("particle_subtraction_prior_scale_pi_n", 1.0) or 1.0),
        "pi_delta": float(fit_config.get("particle_subtraction_prior_scale_pi_delta", 1.5) or 1.5),
        "pi_sidis": float(fit_config.get("particle_subtraction_prior_scale_pi_sidis", 2.0) or 2.0),
    }
    postfit_scale_map = resolve_particle_subtraction_component_postfit_scales("kaon_nosub")
    validation_options = {
        "oversub_sigma_tolerance": fit_config.get("oversub_sigma_tolerance", 2.0),
        "max_oversub_bin_count": fit_config.get("max_oversub_bin_count"),
        "max_oversub_bin_fraction": fit_config.get("max_oversub_bin_fraction"),
        "max_full_range_chi2_ndf": fit_config.get("max_full_range_chi2_ndf"),
    }
    include_kaon_signal_template = bool(fit_config.get("include_kaon_signal_template", False))
    if include_kaon_signal_template and _hist_has_usable_support(h_kaon_signal_shape):
        extra_positive_templates[KAON_SIGNAL_TEMPLATE_NAME] = h_kaon_signal_shape
        prior_scale_map[KAON_SIGNAL_TEMPLATE_NAME] = float(
            fit_config.get("particle_subtraction_prior_scale_k_lambda_signal", 1.0) or 1.0
        )
        if mm_max > mm_min:
            tail_extension = max(float(fit_config.get("kaon_signal_tail_extension", 0.0) or 0.0), 0.0)
            signal_window_max = min(fit_max, mm_max + tail_extension)
            extra_anchor_windows[KAON_SIGNAL_TEMPLATE_NAME] = [(mm_min, signal_window_max)]
    include_sigma0_signal_template = bool(fit_config.get("include_sigma0_signal_template", False))
    if include_sigma0_signal_template and _hist_has_usable_support(h_kaon_sigma0_shape):
        extra_positive_templates[KAON_SIGMA0_TEMPLATE_NAME] = h_kaon_sigma0_shape
        prior_scale_map[KAON_SIGMA0_TEMPLATE_NAME] = float(
            fit_config.get("particle_subtraction_prior_scale_k_sigma0_signal", 1.0) or 1.0
        )
        sigma0_anchor_window = fit_config.get("sigma0_signal_anchor_window")
        if sigma0_anchor_window is not None and len(sigma0_anchor_window) == 2:
            sigma0_window_min = float(sigma0_anchor_window[0])
            sigma0_window_max = float(sigma0_anchor_window[1])
            if bool(fit_config.get("apply_mm_offset_data", False)):
                sigma0_window_min += float(mm_offset_data)
                sigma0_window_max += float(mm_offset_data)
            extra_anchor_windows[KAON_SIGMA0_TEMPLATE_NAME] = [
                (sigma0_window_min, sigma0_window_max)
            ]
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
    signal_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    signal_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    sigma0_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
    sigma0_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGMA0_TEMPLATE_NAME)
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
        "k_lambda_reference_hist": signal_reference_hist,
        "step_overlays": result.get("step_overlays") or [],
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
        ):
            hist = return_payload.get(key)
            if hist is None:
                continue
            return_payload[key] = _masked_hist_clone(
                hist,
                "{}_masked".format(hist.GetName()),
                excluded_windows,
            )
    return return_payload


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

    pion_fit = fit_pion_control_with_simc_shapes(
        h_pion_control,
        aligned_component_shapes.get("pi_n"),
        aligned_component_shapes.get("pi_delta"),
        aligned_component_shapes.get("pi_sidis"),
        aligned_kaon_sigma0_shape,
        inpDict,
        mm_offset_data=mm_offset_data,
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
        context=context,
    )

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

    result = {
        "particle_subtraction_mode": mode,
        "analysis_scope": analysis_scope,
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
        "fallback_used": bool(fallback_reasons),
        "fallback_reason": "; ".join(fallback_reasons),
        "diagnostics": {
            "pion": deepcopy(pion_fit["diagnostics"]),
            "kaon": deepcopy(kaon_fit["diagnostics"]),
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
        "H_pion_fit_pi_n_scaled": pion_fit["pi_n_scaled_hist"],
        "H_pion_fit_pi_delta_scaled": pion_fit["pi_delta_scaled_hist"],
        "H_pion_fit_pi_sidis_scaled": pion_fit["pi_sidis_scaled_hist"],
        "H_pion_fit_k_sigma0_scaled": pion_fit["k_sigma0_scaled_hist"],
        "H_pion_fit_total": pion_fit["fit_hist"],
        "H_pion_fit_step_overlays": pion_fit.get("step_overlays") or [],
        "H_kaon_fit_pi_n_scaled": kaon_fit["pi_n_scaled_hist"],
        "H_kaon_fit_pi_delta_scaled": kaon_fit["pi_delta_scaled_hist"],
        "H_kaon_fit_pi_sidis_scaled": kaon_fit["pi_sidis_scaled_hist"],
        "H_kaon_fit_k_lambda_scaled": kaon_fit["k_lambda_scaled_hist"],
        "H_kaon_fit_k_sigma0_scaled": kaon_fit["k_sigma0_scaled_hist"],
        "H_kaon_fit_k_lambda_reference": kaon_fit.get("k_lambda_reference_hist"),
        "H_kaon_fit_total": kaon_fit["fit_hist"],
        "H_kaon_pion_bg_fit_total": kaon_fit["pion_bg_fit_hist"],
        "H_kaon_fit_step_overlays": kaon_fit.get("step_overlays") or [],
        "H_fit_residual_pion": pion_fit["residual_hist"],
        "H_fit_residual_kaon": kaon_fit["residual_hist"],
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
    return "n={:.2f}, delta={:.2f}, sidis={:.2f}".format(
        float(scale_map.get("pi_n", 1.0) or 1.0),
        float(scale_map.get("pi_delta", 1.0) or 1.0),
        float(scale_map.get("pi_sidis", 1.0) or 1.0),
    )


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
        component_clone = _clone_hist(component_scaled, "{}_step_component".format(component_scaled.GetName()))
        cumulative_clone = _clone_hist(cumulative_after, "{}_step_cumulative".format(cumulative_after.GetName()))
        excluded_windows = step_overlay.get("excluded_windows") or []
        if excluded_windows:
            _mask_hist_windows_inplace(baseline_clone, excluded_windows)
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
        _style_overlay_hist(component_clone, component_color, line_style=1)
        _style_overlay_hist(cumulative_clone, ROOT.kGreen + 2, line_style=3)
        top_y_max = max(
            target_clone.GetMaximum(),
            baseline_clone.GetMaximum(),
            component_clone.GetMaximum(),
            cumulative_clone.GetMaximum(),
            0.0,
        )
        if top_y_max <= 0.0:
            top_y_max = 1.0
        target_clone.SetMaximum(1.20 * top_y_max)
        target_clone.SetMinimum(0.0)
        target_clone.Draw("hist")
        baseline_clone.Draw("hist same")
        component_clone.Draw("hist same")
        cumulative_clone.Draw("hist same")
        _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            0.0,
            1.20 * top_y_max,
            ROOT.kBlue + 1,
            3,
        )
        _draw_window_collection(
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
        top_legend.AddEntry(component_clone, "{} contribution".format(component_label), "l")
        top_legend.AddEntry(cumulative_clone, "baseline after step", "l")
        top_legend.Draw()

        stats_box = ROOT.TPaveText(0.14, 0.60, 0.52, 0.88, "NDC")
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.028)
        stats_box.AddText("pass: {}".format(step_overlay.get("pass_index", 0)))
        stats_box.AddText("component: {}".format(component_label))
        stats_box.AddText("amplitude: {}".format(_format_fit_number(step_overlay.get("amplitude"))))
        if str(step_overlay.get("amplitude_mode") or "least_squares") != "least_squares":
            stats_box.AddText("mode: {}".format(str(step_overlay.get("amplitude_mode"))))
        if step_overlay.get("amplitude_windows"):
            stats_box.AddText(
                "core: {}".format(_format_window_list(step_overlay.get("amplitude_windows") or []))
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
        component_bottom_clone = _clone_hist(
            component_scaled,
            "{}_step_component_bottom".format(component_scaled.GetName()),
        )
        if excluded_windows:
            _mask_hist_windows_inplace(residual_clone, excluded_windows)
            _mask_hist_windows_inplace(component_bottom_clone, excluded_windows)
        residual_clone.SetTitle("Residual input to {} step".format(component_label))
        residual_clone.SetLineColor(ROOT.kBlack)
        residual_clone.SetLineWidth(2)
        residual_clone.SetFillStyle(3001)
        residual_clone.SetFillColor(ROOT.kGray + 1)
        residual_clone.SetMarkerStyle(20)
        residual_clone.SetMarkerSize(0.7)
        _style_overlay_hist(component_bottom_clone, component_color, line_style=1)
        bottom_y_max = max(
            residual_clone.GetMaximum(),
            component_bottom_clone.GetMaximum(),
            0.0,
        )
        if bottom_y_max <= 0.0:
            bottom_y_max = 1.0
        bottom_y_min = min(
            0.0,
            residual_clone.GetMinimum(),
            component_bottom_clone.GetMinimum(),
        )
        residual_clone.SetMaximum(1.20 * bottom_y_max)
        residual_clone.SetMinimum(1.20 * bottom_y_min if bottom_y_min < 0.0 else 0.0)
        residual_clone.Draw("hist")
        component_bottom_clone.Draw("hist same")
        _draw_window_collection(
            step_overlay.get("anchor_windows") or [],
            residual_clone.GetMinimum(),
            residual_clone.GetMaximum(),
            ROOT.kBlue + 1,
            3,
        )
        _draw_window_collection(
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
        bottom_legend.AddEntry(component_bottom_clone, "{} fit".format(component_label), "l")
        bottom_legend.Draw()

        canvas.Print(pdf_name)
        canvas.Close()


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
    model_closure = diagnostics.get("model_closure") or {}
    event_template_closure = diagnostics.get("event_template_closure") or {}

    _print_single_hist_page(
        pdf_name,
        component_payload.get("H_pion_weight_vs_MM"),
        "w_pi(MM)",
        "{}Part 3 pion weight vs MM".format(title_prefix),
        [
            "scope: {}".format(scope_label),
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
        ],
        cut_window=cut_window,
        line_color=ROOT.kViolet + 1,
    )

    _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        "kaon pion-bg model",
        "{}Part 3 model closure: weighted pion model vs kaon bg model".format(title_prefix),
        [
            (component_payload.get("H_weighted_pion_control_model"), "weighted pion-control model", ROOT.kOrange + 7, 2),
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
        ],
        cut_window=cut_window,
    )

    _print_component_overlay_page(
        pdf_name,
        component_payload.get("H_kaon_pion_model"),
        "kaon pion-bg model",
        "{}Part 3 event-template closure vs kaon bg model".format(title_prefix),
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

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_pion_control_input"),
        "pion-control data",
        "{}pion-control SIMC component fit".format(title_prefix),
        [
            (component_fit_result.get("H_pion_fit_pi_n_scaled"), "pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_pion_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_pion_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pion_fit_k_sigma0_scaled"), "K-Sigma0", ROOT.kCyan + 2, 1),
            (component_fit_result.get("H_pion_fit_total"), "total fit", ROOT.kGreen + 2, 2),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: {}".format(component_fit_result.get("fit_status_pion", "unknown")),
            "strategy: {}".format(
                _format_fit_strategy(((component_fit_result.get("diagnostics") or {}).get("pion") or {}))
            ),
            "solution: {}".format(
                _format_solution_method(((component_fit_result.get("diagnostics") or {}).get("pion") or {}))
            ),
            "validation: {}".format(
                _format_validation_status(((component_fit_result.get("diagnostics") or {}).get("pion") or {}))
            ),
            "template MM shift={:.6f}".format(
                float(component_fit_result.get("template_mm_offset_data") or 0.0)
            ),
            "post-fit scales: {}".format(
                _format_component_scale_map(
                    ((component_fit_result.get("diagnostics") or {}).get("pion") or {}).get("postfit_component_scales")
                )
            ),
            "B_n={}  B_delta={}  B_sidis={}".format(
                _format_fit_number(component_fit_result.get("B_n")),
                _format_fit_number(component_fit_result.get("B_delta")),
                _format_fit_number(component_fit_result.get("B_sidis")),
            ),
            "K-Sigma0 scale={}".format(
                _format_fit_number(component_fit_result.get("B_sigma0"))
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(component_fit_result.get("chi2_ndf_pion")),
                _format_fit_metric(component_fit_result.get("fit_p_value_pion")),
            ),
            "anchor windows: {}".format(
                _format_window_list(
                    ((component_fit_result.get("diagnostics") or {}).get("pion") or {}).get("include_windows")
                )
            ),
            "excluded windows: {}".format(
                _format_excluded_window_list(
                    ((component_fit_result.get("diagnostics") or {}).get("pion") or {}).get("exclude_windows")
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
        kaon_title,
        kaon_overlay_specs,
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: {}".format(component_fit_result.get("fit_status_kaon", "unknown")),
            "strategy: {}".format(
                _format_fit_strategy(((component_fit_result.get("diagnostics") or {}).get("kaon") or {}))
            ),
            "solution: {}".format(
                _format_solution_method(((component_fit_result.get("diagnostics") or {}).get("kaon") or {}))
            ),
            "validation: {}".format(
                _format_validation_status(((component_fit_result.get("diagnostics") or {}).get("kaon") or {}))
            ),
            "template MM shift={:.6f}".format(
                float(component_fit_result.get("template_mm_offset_data") or 0.0)
            ),
            "post-fit scales: {}".format(
                _format_component_scale_map(
                    ((component_fit_result.get("diagnostics") or {}).get("kaon") or {}).get("postfit_component_scales")
                )
            ),
            "A_n={}  A_delta={}  A_sidis={}".format(
                _format_fit_number(component_fit_result.get("A_n")),
                _format_fit_number(component_fit_result.get("A_delta")),
                _format_fit_number(component_fit_result.get("A_sidis")),
            ),
            "K-Sigma0 scale={}".format(
                _format_fit_number(component_fit_result.get("S_sigma0"))
            ) if has_sigma0_component else "K-Sigma0 scale=n/a",
            "K-Lambda gauge scale={}".format(
                _format_fit_number(component_fit_result.get("S_lambda_reference_scale"))
            ) if has_kaon_signal_reference else "K-Lambda gauge scale=n/a",
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(component_fit_result.get("chi2_ndf_kaon")),
                _format_fit_metric(component_fit_result.get("fit_p_value_kaon")),
            ),
            "anchor windows: {}".format(
                _format_window_list(
                    ((component_fit_result.get("diagnostics") or {}).get("kaon") or {}).get("include_windows")
                )
            ),
            "excluded windows: {}".format(
                _format_excluded_window_list(
                    ((component_fit_result.get("diagnostics") or {}).get("kaon") or {}).get("exclude_windows")
                )
            ),
        ],
        cut_window=cut_window,
    )

    if component_fit_result.get("analysis_scope") == "setting-wide":
        _print_component_step_pages(
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


def serialize_particle_subtraction_component_result(result):
    if not isinstance(result, dict):
        return {}
    serializable = {}
    for key, value in result.items():
        if _is_root_hist(value):
            continue
        if key.startswith("H_"):
            continue
        serializable[key] = deepcopy(value)
    return serializable
