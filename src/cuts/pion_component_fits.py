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
    resolve_particle_subtraction_mode,
)
from utility import normalize_hist_to_unit_area


COMPONENT_NAMES = ("pi_n", "pi_delta", "pi_sidis")


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
    }
    if amplitude_prefix == "A":
        result["pion_bg_fit_hist"] = _clone_hist(
            target_hist,
            "{}_pion_bg_fit_{}".format(amplitude_prefix, context),
            reset=True,
        )
    return result


def _validate_component_shapes(component_hists, target_hist):
    if target_hist is None:
        return "missing target histogram"
    if any(component_hists.get(name) is None for name in COMPONENT_NAMES):
        return "missing SIMC component shape"

    nbins = target_hist.GetNbinsX()
    xmin = target_hist.GetXaxis().GetXmin()
    xmax = target_hist.GetXaxis().GetXmax()
    for component_name in COMPONENT_NAMES:
        hist = component_hists[component_name]
        if hist.GetNbinsX() != nbins:
            return "bin-count mismatch for {}".format(component_name)
        if abs(hist.GetXaxis().GetXmin() - xmin) > 1e-9 or abs(hist.GetXaxis().GetXmax() - xmax) > 1e-9:
            return "axis-range mismatch for {}".format(component_name)
        if hist.Integral() <= 0.0:
            return "non-positive integral for {}".format(component_name)
    return ""


def _build_fit_inputs(
    target_hist,
    component_hists,
    fit_min,
    fit_max,
    exclude_windows=None,
):
    if exclude_windows is None:
        exclude_windows = []

    component_columns = {name: [] for name in COMPONENT_NAMES}
    x_values = []
    y_values = []
    sigma_values = []
    fit_bin_indices = []

    for bin_index in range(1, target_hist.GetNbinsX() + 1):
        x_center = float(target_hist.GetBinCenter(bin_index))
        if x_center < fit_min or x_center > fit_max:
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

    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "sigma": np.asarray(sigma_values, dtype=float),
        "fit_bin_indices": fit_bin_indices,
        "component_columns": {
            component_name: np.asarray(values, dtype=float)
            for component_name, values in component_columns.items()
        },
    }


def _run_component_fit(
    target_hist,
    component_hists,
    amplitude_prefix,
    fit_min,
    fit_max,
    include_linear_background=False,
    exclude_windows=None,
    context="",
):
    fallback_reason = _validate_component_shapes(component_hists, target_hist)
    if fallback_reason:
        return _zero_fit_result(target_hist, amplitude_prefix, context, fallback_reason)

    fit_inputs = _build_fit_inputs(
        target_hist,
        component_hists,
        fit_min,
        fit_max,
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
        total_fit_value = pion_bg_value
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
        "exclude_windows": deepcopy(exclude_windows or []),
        "fallback_used": False,
        "fallback_reason": "",
        "component_amplitudes": {
            component_name: float(parameter_map[component_name])
            for component_name in COMPONENT_NAMES
        },
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
    context="",
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    result = _run_component_fit(
        h_pion_control,
        {
            "pi_n": h_pi_n_shape,
            "pi_delta": h_pi_delta_shape,
            "pi_sidis": h_pi_sidis_shape,
        },
        "B",
        fit_min,
        fit_max,
        include_linear_background=False,
        exclude_windows=None,
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
    inpDict,
    context="",
):
    fit_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    fit_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    exclude_windows = []
    mm_min = float(inpDict.get("mm_min", fit_min))
    mm_max = float(inpDict.get("mm_max", fit_max))
    if mm_max > mm_min:
        exclude_windows.append((mm_min, mm_max))
    result = _run_component_fit(
        h_kaon_nosub,
        {
            "pi_n": h_pi_n_shape,
            "pi_delta": h_pi_delta_shape,
            "pi_sidis": h_pi_sidis_shape,
        },
        "A",
        fit_min,
        fit_max,
        include_linear_background=True,
        exclude_windows=exclude_windows,
        context="{}_kaon_nosub".format(context or "scope"),
    )
    return {
        "A_n": result["A_n"],
        "A_delta": result["A_delta"],
        "A_sidis": result["A_sidis"],
        "fit_status": result["fit_status"],
        "diagnostics": result["diagnostics"],
        "fit_hist": result["fit_hist"],
        "pion_bg_fit_hist": result["pion_bg_fit_hist"],
        "residual_hist": result["residual_hist"],
        "pi_n_scaled_hist": result["pi_n_scaled_hist"],
        "pi_delta_scaled_hist": result["pi_delta_scaled_hist"],
        "pi_sidis_scaled_hist": result["pi_sidis_scaled_hist"],
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


def build_particle_subtraction_component_result(
    h_pion_control,
    h_kaon_nosub,
    component_shapes,
    inpDict,
    analysis_scope,
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
        context=context,
    )
    kaon_fit = fit_kaon_nosub_with_simc_pion_shapes(
        h_kaon_nosub,
        component_shapes.get("pi_n"),
        component_shapes.get("pi_delta"),
        component_shapes.get("pi_sidis"),
        inpDict,
        context=context,
    )

    b_n = float(pion_fit["B_n"])
    b_delta = float(pion_fit["B_delta"])
    b_sidis = float(pion_fit["B_sidis"])
    a_n = float(kaon_fit["A_n"])
    a_delta = float(kaon_fit["A_delta"])
    a_sidis = float(kaon_fit["A_sidis"])

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
        "H_pion_fit_pi_n_scaled": pion_fit["pi_n_scaled_hist"],
        "H_pion_fit_pi_delta_scaled": pion_fit["pi_delta_scaled_hist"],
        "H_pion_fit_pi_sidis_scaled": pion_fit["pi_sidis_scaled_hist"],
        "H_pion_fit_total": pion_fit["fit_hist"],
        "H_kaon_fit_pi_n_scaled": kaon_fit["pi_n_scaled_hist"],
        "H_kaon_fit_pi_delta_scaled": kaon_fit["pi_delta_scaled_hist"],
        "H_kaon_fit_pi_sidis_scaled": kaon_fit["pi_sidis_scaled_hist"],
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
            "B_n={}  B_delta={}  B_sidis={}".format(
                _format_fit_number(component_fit_result.get("B_n")),
                _format_fit_number(component_fit_result.get("B_delta")),
                _format_fit_number(component_fit_result.get("B_sidis")),
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(component_fit_result.get("chi2_ndf_pion")),
                _format_fit_metric(component_fit_result.get("fit_p_value_pion")),
            ),
        ],
    )

    _print_component_overlay_page(
        pdf_name,
        component_fit_result.get("H_kaon_nosub_input"),
        "kaon no-sub data",
        "{}kaon no-sub SIMC pion-background fit".format(title_prefix),
        [
            (component_fit_result.get("H_kaon_fit_pi_n_scaled"), "pi-n", ROOT.kRed + 1, 1),
            (component_fit_result.get("H_kaon_fit_pi_delta_scaled"), "pi-delta", ROOT.kAzure + 2, 1),
            (component_fit_result.get("H_kaon_fit_pi_sidis_scaled"), "pi-SIDIS", ROOT.kMagenta + 2, 1),
            (component_fit_result.get("H_kaon_pion_bg_fit_total"), "pion-bg sum", ROOT.kOrange + 7, 2),
            (component_fit_result.get("H_kaon_fit_total"), "total fit", ROOT.kGreen + 2, 3),
        ],
        [
            "scope: {}".format(component_fit_result.get("analysis_scope", "unknown")),
            "status: {}".format(component_fit_result.get("fit_status_kaon", "unknown")),
            "A_n={}  A_delta={}  A_sidis={}".format(
                _format_fit_number(component_fit_result.get("A_n")),
                _format_fit_number(component_fit_result.get("A_delta")),
                _format_fit_number(component_fit_result.get("A_sidis")),
            ),
            "chi2/ndf={}  p={}".format(
                _format_fit_metric(component_fit_result.get("chi2_ndf_kaon")),
                _format_fit_metric(component_fit_result.get("fit_p_value_kaon")),
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
