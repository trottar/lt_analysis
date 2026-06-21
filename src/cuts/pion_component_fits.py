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
    get_particle_subtraction_component_fit_window_config,
    resolve_particle_subtraction_component_fit_windows,
    resolve_particle_subtraction_mode,
)
from utility import normalize_hist_to_unit_area


COMPONENT_NAMES = ("pi_n", "pi_delta", "pi_sidis")
KAON_SIGNAL_TEMPLATE_NAME = "k_lambda_signal"


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
):
    template_name = getattr(template_hist, "GetName", lambda: "template")()
    validation_message = _validate_template_hist(template_hist, target_hist, template_name)
    if validation_message:
        return {
            "success": False,
            "amplitude": 0.0,
            "chi2": None,
            "n_fit_bins": 0,
            "message": validation_message,
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
            "chi2": None,
            "n_fit_bins": 0,
            "message": "no valid fit bins",
        }

    weighted_template = fit_inputs["template"] / fit_inputs["sigma"]
    denominator = float(np.dot(weighted_template, weighted_template))
    if (not math.isfinite(denominator)) or denominator <= 0.0:
        return {
            "success": False,
            "amplitude": 0.0,
            "chi2": None,
            "n_fit_bins": int(len(fit_inputs["x"])),
            "message": "template has zero support inside anchor window",
        }

    weighted_target = fit_inputs["y"] / fit_inputs["sigma"]
    amplitude = max(float(np.dot(weighted_template, weighted_target) / denominator), 0.0)
    residual = fit_inputs["y"] - amplitude * fit_inputs["template"]
    chi2_value = float(np.sum(np.square(residual / fit_inputs["sigma"])))
    return {
        "success": True,
        "amplitude": amplitude,
        "chi2": chi2_value,
        "n_fit_bins": int(len(fit_inputs["x"])),
        "message": "",
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


def _fit_staged_anchor_templates(
    target_hist,
    component_hists,
    amplitude_prefix,
    fit_min,
    fit_max,
    anchor_windows,
    fit_order,
    extra_positive_templates=None,
    extra_anchor_windows=None,
    n_passes=3,
    context="",
):
    if extra_positive_templates is None:
        extra_positive_templates = {}
    if extra_anchor_windows is None:
        extra_anchor_windows = {}

    fallback_reason = _validate_component_shapes(component_hists, target_hist)
    if fallback_reason:
        return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)
    for template_name, template_hist in extra_positive_templates.items():
        fallback_reason = _validate_template_hist(template_hist, target_hist, template_name)
        if fallback_reason:
            return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)

    anchor_window_map = _coerce_window_map(anchor_windows)
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

    amplitude_map = {component_name: 0.0 for component_name in template_hists}
    staged_pass_history = []
    total_passes = max(int(n_passes or 1), 1)
    for pass_index in range(total_passes):
        pass_summary = {}
        for component_name in fitted_names:
            residual_hist = _clone_hist(
                target_hist,
                "{}_residual_{}_pass{}".format(
                    amplitude_prefix,
                    component_name,
                    pass_index + 1,
                ),
            )
            for other_name, other_hist in template_hists.items():
                if other_name == component_name:
                    continue
                other_amplitude = float(amplitude_map.get(other_name, 0.0) or 0.0)
                if other_hist is None or other_amplitude == 0.0:
                    continue
                residual_hist.Add(other_hist, -other_amplitude)

            solve_result = _solve_nonnegative_template_amplitude(
                residual_hist,
                template_hists[component_name],
                fit_min,
                fit_max,
                include_windows=combined_window_map.get(component_name),
            )
            amplitude_map[component_name] = float(solve_result.get("amplitude", 0.0) or 0.0)
            pass_summary[component_name] = {
                "amplitude": amplitude_map[component_name],
                "chi2": solve_result.get("chi2"),
                "n_fit_bins": solve_result.get("n_fit_bins"),
                "success": bool(solve_result.get("success", False)),
                "message": solve_result.get("message", ""),
            }
        staged_pass_history.append(pass_summary)

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
    pi_n_scaled_hist.Scale(float(amplitude_map.get("pi_n", 0.0)))
    pi_delta_scaled_hist.Scale(float(amplitude_map.get("pi_delta", 0.0)))
    pi_sidis_scaled_hist.Scale(float(amplitude_map.get("pi_sidis", 0.0)))

    extra_scaled_hists = {}
    for template_name, template_hist in extra_positive_templates.items():
        scaled_hist = _clone_hist(
            template_hist,
            "{}_{}_scaled_{}".format(amplitude_prefix, template_name, context),
        )
        scaled_hist.Scale(float(amplitude_map.get(template_name, 0.0)))
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
    )

    for component_name, scaled_hist in (
        ("pi_n", pi_n_scaled_hist),
        ("pi_delta", pi_delta_scaled_hist),
        ("pi_sidis", pi_sidis_scaled_hist),
    ):
        if scaled_hist is None:
            continue
        fit_hist.Add(scaled_hist)
        if amplitude_prefix == "A":
            pion_bg_fit_hist.Add(scaled_hist)
        residual_hist.Add(scaled_hist, -1.0)

    for scaled_hist in extra_scaled_hists.values():
        if scaled_hist is None:
            continue
        fit_hist.Add(scaled_hist)
        residual_hist.Add(scaled_hist, -1.0)

    include_windows = _collect_unique_windows(
        combined_window_map,
        ordered_names=fitted_names,
    )
    quality = _compute_fit_quality(
        target_hist,
        fit_hist,
        fit_min,
        fit_max,
        include_windows=include_windows,
        n_parameters=len(fitted_names),
    )
    extra_component_amplitudes = {
        template_name: float(amplitude_map.get(template_name, 0.0) or 0.0)
        for template_name in extra_positive_templates
    }
    diagnostics = {
        "success": True,
        "status_code": None,
        "message": "staged anchored component fit",
        "chi2": quality["chi2"],
        "ndf": quality["ndf"],
        "chi2_ndf": quality["chi2_ndf"],
        "fit_p_value": quality["fit_p_value"],
        "n_fit_bins": quality["n_fit_bins"],
        "fit_min": float(fit_min),
        "fit_max": float(fit_max),
        "include_windows": deepcopy(include_windows),
        "exclude_windows": [],
        "fallback_used": False,
        "fallback_reason": "",
        "fit_strategy": "staged_anchor",
        "n_passes": total_passes,
        "fit_order": list(fitted_names),
        "anchor_windows": deepcopy(combined_window_map),
        "staged_pass_history": deepcopy(staged_pass_history),
        "component_amplitudes": {
            component_name: float(amplitude_map.get(component_name, 0.0) or 0.0)
            for component_name in COMPONENT_NAMES
        },
        "extra_component_amplitudes": deepcopy(extra_component_amplitudes),
    }

    result = {
        "{}_n".format(amplitude_prefix): float(amplitude_map.get("pi_n", 0.0) or 0.0),
        "{}_delta".format(amplitude_prefix): float(amplitude_map.get("pi_delta", 0.0) or 0.0),
        "{}_sidis".format(amplitude_prefix): float(amplitude_map.get("pi_sidis", 0.0) or 0.0),
        "fit_status": "success",
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
    anchor_windows = {
        component_name: [window]
        for component_name, window in resolved_windows.items()
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
        fit_order=fit_config.get("fit_order") or ("pi_n", "pi_sidis", "pi_delta"),
        n_passes=fit_config.get("staged_fit_passes", 3),
        context="{}_pion_control".format(context or "scope"),
    )
    return {
        "B_n": result["B_n"],
        "B_delta": result["B_delta"],
        "B_sidis": result["B_sidis"],
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "fit_hist": result["fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
    }


def fit_kaon_nosub_with_simc_pion_shapes(
    h_kaon_nosub,
    h_pi_n_shape,
    h_pi_delta_shape,
    h_pi_sidis_shape,
    h_kaon_signal_shape,
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
    anchor_windows = {
        component_name: [window]
        for component_name, window in resolved_windows.items()
    }
    mm_min = float(inpDict.get("mm_min", fit_min))
    mm_max = float(inpDict.get("mm_max", fit_max))
    extra_positive_templates = {}
    extra_anchor_windows = {}
    if h_kaon_signal_shape is not None:
        extra_positive_templates[KAON_SIGNAL_TEMPLATE_NAME] = h_kaon_signal_shape
        if mm_max > mm_min:
            tail_extension = max(float(fit_config.get("kaon_signal_tail_extension", 0.0) or 0.0), 0.0)
            signal_window_max = min(fit_max, mm_max + tail_extension)
            extra_anchor_windows[KAON_SIGNAL_TEMPLATE_NAME] = [(mm_min, signal_window_max)]
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
        fit_order=fit_config.get("fit_order") or ("pi_n", "pi_sidis", "pi_delta"),
        extra_positive_templates=extra_positive_templates,
        extra_anchor_windows=extra_anchor_windows,
        n_passes=fit_config.get("staged_fit_passes", 3),
        context="{}_kaon_nosub".format(context or "scope"),
    )
    signal_scaled_hist = (result.get("extra_scaled_hists") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    signal_amplitude = (result.get("extra_component_amplitudes") or {}).get(KAON_SIGNAL_TEMPLATE_NAME)
    return {
        "A_n": result["A_n"],
        "A_delta": result["A_delta"],
        "A_sidis": result["A_sidis"],
        "S_lambda": None if signal_amplitude is None else float(signal_amplitude),
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "fit_hist": result["fit_hist"],
        "pion_bg_fit_hist": result["pion_bg_fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
        "k_lambda_scaled_hist": signal_scaled_hist,
    }


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
        if t_bin_index is None and phi_bin_index is None:
            resolved[component_name] = component_entry.get("setting_shape_full")
            continue

        binned_shapes = component_entry.get("binned_shapes") or {}
        t_key = "t_bin{}".format(int(t_bin_index) + 1) if t_bin_index is not None else None
        if phi_bin_index is not None and t_key is not None:
            resolved[component_name] = (
                (binned_shapes.get(t_key) or {}).get("phi_bin{}".format(int(phi_bin_index) + 1))
            )
            continue

        if t_key is not None:
            phi_shape_map = binned_shapes.get(t_key) or {}
            resolved[component_name] = _sum_hist_list_to_unit_area(
                list(phi_shape_map.values()),
                "{}_{}_aggregated".format(component_name, analysis_scope),
            )
            continue

        resolved[component_name] = component_entry.get("setting_shape_full")
    return resolved


def resolve_scope_single_shape(
    shape_payload,
    analysis_scope="setting-wide",
    t_bin_index=None,
    phi_bin_index=None,
):
    if not isinstance(shape_payload, dict):
        return None
    if t_bin_index is None and phi_bin_index is None:
        return shape_payload.get("setting_shape_full")

    binned_shapes = shape_payload.get("binned_shapes") or {}
    t_key = "t_bin{}".format(int(t_bin_index) + 1) if t_bin_index is not None else None
    if phi_bin_index is not None and t_key is not None:
        return (binned_shapes.get(t_key) or {}).get("phi_bin{}".format(int(phi_bin_index) + 1))
    if t_key is not None:
        return _sum_hist_list_to_unit_area(
            list((binned_shapes.get(t_key) or {}).values()),
            "single_shape_{}_aggregated".format(analysis_scope or "scope"),
        )
    return shape_payload.get("setting_shape_full")


def build_particle_subtraction_component_result(
    h_pion_control,
    h_kaon_nosub,
    component_shapes,
    inpDict,
    analysis_scope,
    kaon_signal_shape=None,
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

    pion_fit = fit_pion_control_with_simc_shapes(
        h_pion_control,
        component_shapes.get("pi_n"),
        component_shapes.get("pi_delta"),
        component_shapes.get("pi_sidis"),
        inpDict,
        mm_offset_data=mm_offset_data,
        context=context,
    )
    kaon_fit = fit_kaon_nosub_with_simc_pion_shapes(
        h_kaon_nosub,
        component_shapes.get("pi_n"),
        component_shapes.get("pi_delta"),
        component_shapes.get("pi_sidis"),
        kaon_signal_shape,
        inpDict,
        mm_offset_data=mm_offset_data,
        context=context,
    )

    b_n = float(pion_fit["B_n"])
    b_delta = float(pion_fit["B_delta"])
    b_sidis = float(pion_fit["B_sidis"])
    a_n = float(kaon_fit["A_n"])
    a_delta = float(kaon_fit["A_delta"])
    a_sidis = float(kaon_fit["A_sidis"])
    s_lambda = kaon_fit["S_lambda"]

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
        "B_n": b_n,
        "B_delta": b_delta,
        "B_sidis": b_sidis,
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
        "fallback_used": bool(fallback_reasons),
        "fallback_reason": "; ".join(fallback_reasons),
        "diagnostics": {
            "pion": deepcopy(pion_fit["diagnostics"]),
            "kaon": deepcopy(kaon_fit["diagnostics"]),
        },
        "H_simc_shape_pi_n": _clone_hist(
            component_shapes.get("pi_n"),
            "H_simc_shape_pi_n_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_delta": _clone_hist(
            component_shapes.get("pi_delta"),
            "H_simc_shape_pi_delta_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_pi_sidis": _clone_hist(
            component_shapes.get("pi_sidis"),
            "H_simc_shape_pi_sidis_{}".format(context or analysis_scope),
        ),
        "H_simc_shape_k_lambda": _clone_hist(
            kaon_signal_shape,
            "H_simc_shape_k_lambda_{}".format(context or analysis_scope),
        ),
        "H_pion_fit_pi_n_scaled": pion_fit["pi_n_scaled_hist"],
        "H_pion_fit_pi_delta_scaled": pion_fit["pi_delta_scaled_hist"],
        "H_pion_fit_pi_sidis_scaled": pion_fit["pi_sidis_scaled_hist"],
        "H_pion_fit_total": pion_fit["fit_hist"],
        "H_kaon_fit_pi_n_scaled": kaon_fit["pi_n_scaled_hist"],
        "H_kaon_fit_pi_delta_scaled": kaon_fit["pi_delta_scaled_hist"],
        "H_kaon_fit_pi_sidis_scaled": kaon_fit["pi_sidis_scaled_hist"],
        "H_kaon_fit_k_lambda_scaled": kaon_fit["k_lambda_scaled_hist"],
        "H_kaon_fit_total": kaon_fit["fit_hist"],
        "H_kaon_pion_bg_fit_total": kaon_fit["pion_bg_fit_hist"],
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


def _format_window_list(windows):
    if not windows:
        return "full-range"
    formatted = []
    for window_min, window_max in windows:
        formatted.append("[{:.3f}, {:.3f}]".format(float(window_min), float(window_max)))
    return ", ".join(formatted)


def _draw_vertical_window_lines(window_min, window_max, y_min, y_max):
    lines = []
    for x_value in (window_min, window_max):
        line = ROOT.TLine(float(x_value), float(y_min), float(x_value), float(y_max))
        line.SetLineColor(ROOT.kBlue + 1)
        line.SetLineStyle(3)
        line.SetLineWidth(2)
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
            (component_fit_result.get("H_pion_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_pion_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_pion_fit_total"), "total fit", ROOT.kGreen + 2, 2),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: {}".format(component_fit_result.get("fit_status_pion", "unknown")),
            "strategy: {}".format(
                _format_fit_strategy(((component_fit_result.get("diagnostics") or {}).get("pion") or {}))
            ),
            "B_n={}  B_delta={}  B_sidis={}".format(
                _format_fit_number(component_fit_result.get("B_n")),
                _format_fit_number(component_fit_result.get("B_delta")),
                _format_fit_number(component_fit_result.get("B_sidis")),
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
        ],
    )

    has_kaon_signal = component_fit_result.get("H_kaon_fit_k_lambda_scaled") is not None
    kaon_title = "{}kaon no-sub SIMC pion-background fit".format(title_prefix)
    if has_kaon_signal:
        kaon_title = "{}kaon no-sub SIMC decomposition fit".format(title_prefix)

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_kaon_nosub_input"),
        "kaon no-sub data",
        kaon_title,
        [
            (component_fit_result.get("H_kaon_fit_pi_n_scaled"), "pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_kaon_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_kaon_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_kaon_fit_k_lambda_scaled"), "K-Lambda", ROOT.kBlue + 1, 1),
            (component_fit_result.get("H_kaon_pion_bg_fit_total"), "pion-bg sum", ROOT.kOrange + 7, 2),
            (component_fit_result.get("H_kaon_fit_total"), "total fit", ROOT.kGreen + 2, 3),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: {}".format(component_fit_result.get("fit_status_kaon", "unknown")),
            "strategy: {}".format(
                _format_fit_strategy(((component_fit_result.get("diagnostics") or {}).get("kaon") or {}))
            ),
            "A_n={}  A_delta={}  A_sidis={}  S_lambda={}".format(
                _format_fit_number(component_fit_result.get("A_n")),
                _format_fit_number(component_fit_result.get("A_delta")),
                _format_fit_number(component_fit_result.get("A_sidis")),
                _format_fit_number(component_fit_result.get("S_lambda")),
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(component_fit_result.get("chi2_ndf_kaon")),
                _format_fit_metric(component_fit_result.get("fit_p_value_kaon")),
            ),
            "anchor windows: {}".format(
                _format_window_list(
                    ((component_fit_result.get("diagnostics") or {}).get("kaon") or {}).get("include_windows")
                )
            ),
        ],
        cut_window=cut_window,
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
