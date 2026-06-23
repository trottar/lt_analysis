#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 15:19:50 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
from collections import defaultdict
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

from binning_helpers import find_bin_index

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
sys.path.append("cuts")
from utility import (
    remove_bad_bins,
    get_centroid,
    prune_hist,
    compute_staged_particle_subtraction_scales,
)
from prompt_trees import get_prompt_tree_name, get_rand_tree_name
from background_config import (
    resolve_particle_subtraction_mode,
    resolve_particle_subtraction_windows,
    resolve_bg_stat_scale1,
    resolve_bg_stat_scale2,
    resolve_particle_subtraction_weight_clip_bounds,
    resolve_particle_subtraction_weight_denominator_floor,
    resolve_particle_subtraction_weight_warn_max,
)
from pion_component_shapes import (
    load_kaon_simc_signal_shape,
    load_kaon_simc_sigma0_shape,
    load_setting_pion_component_shapes,
)
from pion_component_fits import (
    build_particle_subtraction_component_result,
    print_particle_subtraction_component_application_pages,
    print_particle_subtraction_component_fit_pages,
    resolve_scope_component_shapes,
    resolve_scope_single_shape,
    serialize_particle_subtraction_component_result,
)
from mm_background_subtraction import (
    build_mm_background_weights,
    build_mm_residual_weights,
    clone_reset_hist,
    mm_background_weight_from_value,
)
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    compute_hist_closure_metrics,
    evaluate_particle_subtraction_component_fit_result,
    fill_simc_shape_pion_subtraction_templates,
    handle_particle_subtraction_fallback,
    iter_component_control_source_specs,
    print_particle_subtraction_weight_support_warning,
    simc_shape_pion_weight_from_value,
    summarize_particle_subtraction_component_payload,
)

##################################################################################################################################################


def _init_ave_hist_group(names, n_t):
    return {name: [None for _ in range(n_t)] for name in names}


def _fill_cached_ave_section(cache_section, hist_group, t_bins, update_mm_offset=False):
    mm_offset_data = None

    for idx, adj_t in enumerate(cache_section["adj_t"]):
        t_index = find_bin_index(adj_t, t_bins)
        if t_index is None:
            continue

        adj_mm = cache_section["adj_MM"][idx]

        if cache_section["nommcuts"][idx]:
            hist_group["fit1sub"][t_index].Fill(adj_mm)
            hist_group["pisub"][t_index].Fill(adj_mm)
            hist_group["nosub"][t_index].Fill(adj_mm)

        if cache_section["allcuts"][idx]:
            hist_group["t"][t_index].Fill(adj_t)
            hist_group["Q2"][t_index].Fill(cache_section["Q2"][idx])
            hist_group["W"][t_index].Fill(cache_section["W"][idx])
            hist_group["epsilon"][t_index].Fill(cache_section["epsilon"][idx])
            hist_group["mm"][t_index].Fill(adj_mm)
            if update_mm_offset:
                mm_offset_data = cache_section["mm_offset"][idx]

    return mm_offset_data


def _resolve_hist_background_sample_root(hist, inpDict, background_name):
    phi_setting = hist.get("phi_setting")
    background_samples = hist.get("background_samples") or {}
    sample_entry = background_samples.get(background_name) or {}
    if not sample_entry:
        background_samples = (((inpDict.get("background_samples") or {}).get("by_phi") or {}).get(phi_setting, {}) or {})
        sample_entry = background_samples.get(background_name) or {}
    return sample_entry.get("root")


def _get_cached_kaon_signal_shape_payload(hist, inpDict, t_bins, phi_bins, context):
    if resolve_particle_subtraction_mode(inpDict) != "simc_shape_components":
        return None
    if str(inpDict.get("ParticleType", "")).strip().lower() != "kaon":
        return None

    cache_key = "_simc_kaon_signal_shape_payload"
    cached_payload = hist.get(cache_key)
    if cached_payload is not None:
        return cached_payload

    root_filename = hist.get("rootFileSimc")
    if not root_filename:
        simc_name = hist.get("InSIMCFilename")
        if simc_name:
            root_filename = os.path.join(OUTPATH, simc_name)

    cached_payload = load_kaon_simc_signal_shape(
        root_filename,
        inpDict,
        hist["phi_setting"],
        t_bins=t_bins,
        phi_bins=phi_bins,
        context=context,
    )
    hist[cache_key] = cached_payload
    return cached_payload


def _get_cached_kaon_sigma0_shape_payload(hist, inpDict, t_bins, phi_bins, context):
    if resolve_particle_subtraction_mode(inpDict) != "simc_shape_components":
        return None
    if str(inpDict.get("ParticleType", "")).strip().lower() != "kaon":
        return None

    cache_key = "_simc_kaon_sigma0_shape_payload"
    cached_payload = hist.get(cache_key)
    if cached_payload is not None:
        return cached_payload

    root_filename = _resolve_hist_background_sample_root(hist, inpDict, "sigma0")
    cached_payload = load_kaon_simc_sigma0_shape(
        root_filename,
        inpDict,
        hist["phi_setting"],
        t_bins=t_bins,
        phi_bins=phi_bins,
        context=context,
    )
    hist[cache_key] = cached_payload
    return cached_payload


def _fill_cached_ave_simc(cache_section, hist_group, t_bins):
    for idx, minus_t in enumerate(cache_section["minus_t"]):
        t_index = find_bin_index(minus_t, t_bins)
        if t_index is None:
            continue

        iter_weight = cache_section["iter_weight"][idx]
        hist_group["t"][t_index].Fill(minus_t, iter_weight)
        hist_group["Q2"][t_index].Fill(cache_section["Q2"][idx], iter_weight)
        hist_group["W"][t_index].Fill(cache_section["W"][idx], iter_weight)
        hist_group["epsilon"][t_index].Fill(cache_section["epsilon"][idx], iter_weight)


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def _iter_ave_source_specs(
    data_event_cache,
    sub_event_cache,
    normfac_data,
    normfac_dummy,
    nWindows,
    scale_factor,
    pion_component_payload=None,
):
    yield {
        "cache_section": data_event_cache["prompt"],
        "coefficient": float(normfac_data),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": data_event_cache["rand"],
        "coefficient": -float(normfac_data) / float(nWindows),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": data_event_cache["dummy"],
        "coefficient": -float(normfac_dummy),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": data_event_cache["dummy_rand"],
        "coefficient": float(normfac_dummy) / float(nWindows),
        "apply_pion_weight": False,
    }

    if sub_event_cache is None:
        return

    if pion_component_payload is not None:
        for source_spec in iter_component_control_source_specs(
            sub_event_cache,
            normfac_data,
            normfac_dummy,
            nWindows,
            positive_template=False,
        ):
            source_spec = dict(source_spec)
            source_spec["apply_pion_weight"] = True
            yield source_spec
        return

    if scale_factor == 0.0:
        return

    yield {
        "cache_section": sub_event_cache["prompt"],
        "coefficient": -float(scale_factor) * float(normfac_data),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": sub_event_cache["rand"],
        "coefficient": float(scale_factor) * float(normfac_data) / float(nWindows),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": sub_event_cache["dummy"],
        "coefficient": float(scale_factor) * float(normfac_dummy),
        "apply_pion_weight": False,
    }
    yield {
        "cache_section": sub_event_cache["dummy_rand"],
        "coefficient": -float(scale_factor) * float(normfac_dummy) / float(nWindows),
        "apply_pion_weight": False,
    }


def _fill_ave_background_templates_for_tbin(
    template_hists,
    source_specs,
    mm_reference_hist,
    mm_background_weights,
    t_bins,
    t_index,
    pion_component_payload=None,
    residual_weights=None,
):
    for source_spec in source_specs:
        cache_section = source_spec.get("cache_section")
        coeff = float(source_spec.get("coefficient", 0.0) or 0.0)
        apply_pion_weight = bool(source_spec.get("apply_pion_weight"))
        if cache_section is None or coeff == 0.0:
            continue

        bin_indices = cache_section.get("allcut_bin_index", {}).get((t_index, None))
        if bin_indices is None:
            bin_indices = range(len(cache_section["adj_MM"]))

        for idx in bin_indices:
            if not cache_section["allcuts"][idx]:
                continue

            if "t_index" in cache_section:
                event_t_index = int(cache_section["t_index"][idx])
            else:
                event_t_index = find_bin_index(cache_section["adj_t"][idx], t_bins)

            if event_t_index != t_index:
                continue

            adj_mm = cache_section["adj_MM"][idx]
            event_weight = coeff * mm_background_weight_from_value(
                adj_mm,
                mm_reference_hist,
                mm_background_weights,
                residual_weights=residual_weights,
            )
            if apply_pion_weight:
                event_weight *= simc_shape_pion_weight_from_value(
                    adj_mm,
                    pion_component_payload["H_pion_control_model"],
                    pion_component_payload["weights"],
                )
            if event_weight == 0.0:
                continue

            template_hists["Q2"].Fill(cache_section["Q2"][idx], event_weight)
            template_hists["W"].Fill(cache_section["W"][idx], event_weight)
            template_hists["t"].Fill(cache_section["adj_t"][idx], event_weight)
            template_hists["epsilon"].Fill(cache_section["epsilon"][idx], event_weight)


def _subtract_ave_mm_background_for_tbin(
    hist_bin_dict,
    j,
    mm_reference_hist,
    background_hist,
    data_event_cache,
    sub_event_cache,
    normfac_data,
    normfac_dummy,
    nWindows,
    scale_factor,
    t_bins,
    pion_component_payload=None,
    residual_weights=None,
):
    mm_background_weights = build_mm_background_weights(mm_reference_hist, background_hist)
    template_hists = {
        "Q2": clone_reset_hist(hist_bin_dict["H_Q2_DATA_{}".format(j)], "_bg_template"),
        "W": clone_reset_hist(hist_bin_dict["H_W_DATA_{}".format(j)], "_bg_template"),
        "t": clone_reset_hist(hist_bin_dict["H_t_DATA_{}".format(j)], "_bg_template"),
        "epsilon": clone_reset_hist(hist_bin_dict["H_epsilon_DATA_{}".format(j)], "_bg_template"),
    }

    _fill_ave_background_templates_for_tbin(
        template_hists,
        _iter_ave_source_specs(
            data_event_cache,
            sub_event_cache,
            normfac_data,
            normfac_dummy,
            nWindows,
            scale_factor,
            pion_component_payload=pion_component_payload,
        ),
        mm_reference_hist,
        mm_background_weights,
        t_bins,
        j,
        pion_component_payload=pion_component_payload,
        residual_weights=residual_weights,
    )

    hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(template_hists["Q2"], -1)
    hist_bin_dict["H_W_DATA_{}".format(j)].Add(template_hists["W"], -1)
    hist_bin_dict["H_t_DATA_{}".format(j)].Add(template_hists["t"], -1)
    hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(template_hists["epsilon"], -1)

    return mm_background_weights


def _build_ave_component_template_hists(hist_bin_dict, j):
    return {
        "Q2": clone_reset_hist(hist_bin_dict["H_Q2_DATA_{}".format(j)], "_pi_component_template"),
        "W": clone_reset_hist(hist_bin_dict["H_W_DATA_{}".format(j)], "_pi_component_template"),
        "t": clone_reset_hist(hist_bin_dict["H_t_DATA_{}".format(j)], "_pi_component_template"),
        "epsilon": clone_reset_hist(hist_bin_dict["H_epsilon_DATA_{}".format(j)], "_pi_component_template"),
        "mm": clone_reset_hist(hist_bin_dict["H_MM_DATA_{}".format(j)], "_pi_component_template"),
        "mm_nosub": clone_reset_hist(hist_bin_dict["H_MM_nosub_DATA_{}".format(j)], "_pi_component_template"),
    }


def _apply_component_pion_subtraction_for_tbin(
    hist_bin_dict,
    j,
    component_fit_result,
    sub_event_cache,
    normfac_data,
    normfac_dummy,
    nWindows,
    particle_type,
    pol,
    inpDict,
):
    gate_result = evaluate_particle_subtraction_component_fit_result(component_fit_result, inpDict)
    payload = {
        "accepted": False,
        "fallback_used": True,
        "fallback_mode": gate_result.get("fallback_mode"),
        "fallback_reason": gate_result.get("reason") or "component-fit result rejected",
        "analysis_scope": component_fit_result.get("analysis_scope") if isinstance(component_fit_result, dict) else None,
        "particle_subtraction_mode": component_fit_result.get("particle_subtraction_mode") if isinstance(component_fit_result, dict) else None,
        "fit_status_pion": component_fit_result.get("fit_status_pion") if isinstance(component_fit_result, dict) else None,
        "fit_status_kaon": component_fit_result.get("fit_status_kaon") if isinstance(component_fit_result, dict) else None,
        "fit_validation_pion": bool((gate_result.get("diagnostics") or {}).get("fit_validation_pion")),
        "fit_validation_kaon": bool((gate_result.get("diagnostics") or {}).get("fit_validation_kaon")),
    }
    if not gate_result["accepted"]:
        return handle_particle_subtraction_fallback(
            payload,
            payload["fallback_reason"],
            context="ave_per_bin component pion subtraction ({}, t{})".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
            ),
        )
    if sub_event_cache is None:
        return handle_particle_subtraction_fallback(
            payload,
            "missing sub_event_cache for component-weight subtraction",
            context="ave_per_bin component pion subtraction ({}, t{})".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
            ),
        )

    clip_min, clip_max = resolve_particle_subtraction_weight_clip_bounds(inpDict)
    weight_payload = build_simc_shape_pion_control_weights(
        component_fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inpDict),
    )
    print_particle_subtraction_weight_support_warning(
        weight_payload,
        context="ave_per_bin component pion subtraction",
        phi_setting=inpDict.get("phi_setting", ""),
        t_bin=int(j) + 1,
    )

    if weight_payload["diagnostics"]["pion_weight_max"] > resolve_particle_subtraction_weight_warn_max(inpDict):
        print(
            "WARNING: pion component weight exceeded threshold\n"
            "  phi_setting = {}\n"
            "  t_bin = {}\n"
            "  max_weight = {:.4f}".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
                float(weight_payload["diagnostics"]["pion_weight_max"]),
            )
        )

    template_hists = _build_ave_component_template_hists(hist_bin_dict, j)
    fill_stats = fill_simc_shape_pion_subtraction_templates(
        template_hists,
        iter_component_control_source_specs(
            sub_event_cache,
            normfac_data,
            normfac_dummy,
            nWindows,
            positive_template=True,
        ),
        weight_payload["H_pion_control_model"],
        weight_payload["weights"],
        {"t_index": j},
        particle_type,
        pol,
    )

    h_mm_before = clone_reset_hist(hist_bin_dict["H_MM_DATA_{}".format(j)], "_before_pion_sub")
    h_mm_before.Add(hist_bin_dict["H_MM_DATA_{}".format(j)])
    h_mm_nosub_before = clone_reset_hist(hist_bin_dict["H_MM_nosub_DATA_{}".format(j)], "_before_pion_sub_full")
    h_mm_nosub_before.Add(hist_bin_dict["H_MM_nosub_DATA_{}".format(j)])

    hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(template_hists["Q2"], -1.0)
    hist_bin_dict["H_W_DATA_{}".format(j)].Add(template_hists["W"], -1.0)
    hist_bin_dict["H_t_DATA_{}".format(j)].Add(template_hists["t"], -1.0)
    hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(template_hists["epsilon"], -1.0)
    hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Add(template_hists["mm_nosub"], -1.0)
    hist_bin_dict["H_MM_pisub_DATA_{}".format(j)].Add(template_hists["mm_nosub"], -1.0)
    hist_bin_dict["H_MM_DATA_{}".format(j)].Add(template_hists["mm"], -1.0)

    h_mm_after = clone_reset_hist(hist_bin_dict["H_MM_DATA_{}".format(j)], "_after_pion_sub")
    h_mm_after.Add(hist_bin_dict["H_MM_DATA_{}".format(j)])
    h_mm_nosub_after = clone_reset_hist(h_mm_nosub_before, "_after_pion_sub_full")
    h_mm_nosub_after.Add(h_mm_nosub_before)
    h_mm_nosub_after.Add(template_hists["mm_nosub"], -1.0)

    pion_control_integral = _hist_integral(component_fit_result.get("H_pion_control_input"))
    weighted_pion_integral_cut = _hist_integral(template_hists["mm"])
    weighted_pion_integral_full = _hist_integral(template_hists["mm_nosub"])
    effective_scale = weighted_pion_integral_full / pion_control_integral if pion_control_integral > 0.0 else 0.0
    event_template_closure = compute_hist_closure_metrics(
        weight_payload.get("H_kaon_pion_model"),
        template_hists["mm_nosub"],
    )

    payload.update(
        {
            "accepted": True,
            "fallback_used": False,
            "fallback_reason": "",
            "particle_subtraction_effective_scale": float(effective_scale),
            "weighted_pion_integral": float(weighted_pion_integral_full),
            "weighted_pion_integral_cut": float(weighted_pion_integral_cut),
            "weighted_pion_integral_full": float(weighted_pion_integral_full),
            "kaon_integral_before_pion_sub": _hist_integral(h_mm_before),
            "kaon_integral_after_pion_sub": _hist_integral(h_mm_after),
            "kaon_integral_before_pion_sub_full": _hist_integral(h_mm_nosub_before),
            "kaon_integral_after_pion_sub_full": _hist_integral(h_mm_nosub_after),
            "H_pion_control_model": weight_payload["H_pion_control_model"],
            "H_kaon_pion_model": weight_payload["H_kaon_pion_model"],
            "H_weighted_pion_control_model": weight_payload.get("H_weighted_pion_control_model"),
            "H_pion_weight_vs_MM": weight_payload["H_pion_weight_vs_MM"],
            "weights": weight_payload["weights"],
            "H_pion_control_input": component_fit_result.get("H_pion_control_input").Clone(
                "{}_clone".format(component_fit_result.get("H_pion_control_input").GetName())
            ) if component_fit_result.get("H_pion_control_input") is not None else None,
            "H_kaon_nosub_input": component_fit_result.get("H_kaon_nosub_input").Clone(
                "{}_clone".format(component_fit_result.get("H_kaon_nosub_input").GetName())
            ) if component_fit_result.get("H_kaon_nosub_input") is not None else None,
            "H_pion_control_unscaled": component_fit_result.get("H_pion_control_input").Clone(
                "{}_clone".format(component_fit_result.get("H_pion_control_input").GetName())
            ) if component_fit_result.get("H_pion_control_input") is not None else None,
            "H_pion_subtraction_template_MM": template_hists["mm"],
            "H_pion_subtraction_template_MM_nosub": template_hists["mm_nosub"],
            "H_MM_before_pion_subtraction": h_mm_before,
            "H_MM_after_pion_subtraction": h_mm_after,
            "H_MM_nosub_before_pion_subtraction": h_mm_nosub_before,
            "H_MM_nosub_after_pion_subtraction": h_mm_nosub_after,
            "H_pion_fit_step_overlays": deepcopy(component_fit_result.get("H_pion_fit_step_overlays") or []),
            "H_kaon_fit_step_overlays": deepcopy(component_fit_result.get("H_kaon_fit_step_overlays") or []),
            "diagnostics": {
                **dict(weight_payload["diagnostics"]),
                **dict(fill_stats),
                "event_template_closure": event_template_closure,
            },
        }
    )
    return payload


def _process_ave_data_tree(
    tree,
    hist_group,
    shifted_mm_getter,
    shifted_t_getter,
    t_bins,
    particle_type,
    hole_contains,
    evaluate_event,
    mm_min,
    mm_max,
    progress_bar,
    update_mm_offset=False,
):
    mm_offset_data = None
    total_entries = tree.GetEntries()

    for i, evt in enumerate(tree):
        progress_bar(i, total_entries, bar_length=25)

        if particle_type == "kaon":
            base_allcuts, base_nommcuts, _ = evaluate_event(evt, mm_min, mm_max)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            allcuts = base_allcuts and not hole_rejected
            nommcuts = base_nommcuts and not hole_rejected
        else:
            allcuts, nommcuts, _ = evaluate_event(evt, mm_min, mm_max)

        if not (allcuts or nommcuts):
            continue

        adj_MM = shifted_mm_getter(evt)
        adj_t = shifted_t_getter(evt)
        t_index = find_bin_index(adj_t, t_bins)
        if t_index is None:
            continue

        if nommcuts:
            hist_group["fit1sub"][t_index].Fill(adj_MM)
            hist_group["pisub"][t_index].Fill(adj_MM)
            hist_group["nosub"][t_index].Fill(adj_MM)

        if allcuts:
            hist_group["t"][t_index].Fill(adj_t)
            hist_group["Q2"][t_index].Fill(evt.Q2)
            hist_group["W"][t_index].Fill(evt.W)
            hist_group["epsilon"][t_index].Fill(evt.epsilon)
            hist_group["mm"][t_index].Fill(adj_MM)
            if update_mm_offset:
                mm_offset_data = adj_MM - evt.MM

    return mm_offset_data


def _process_ave_simc_tree(
    tree,
    hist_group,
    particle_type,
    hole_contains,
    apply_simc_cuts,
    mm_min,
    mm_max,
    t_bins,
    progress_bar,
):
    total_entries = tree.GetEntries()

    for i, evt in enumerate(tree):
        progress_bar(i, total_entries, bar_length=25)

        if particle_type == "kaon":
            allcuts = apply_simc_cuts(evt, mm_min, mm_max) and not hole_contains(evt.phgcer_x_det, evt.phgcer_y_det)
        else:
            allcuts = apply_simc_cuts(evt, mm_min, mm_max)

        if not allcuts:
            continue

        minus_t = -evt.t
        t_index = find_bin_index(minus_t, t_bins)
        if t_index is None:
            continue

        hist_group["t"][t_index].Fill(minus_t, evt.iter_weight)
        hist_group["Q2"][t_index].Fill(evt.Q2, evt.iter_weight)
        hist_group["W"][t_index].Fill(evt.W, evt.iter_weight)
        hist_group["epsilon"][t_index].Fill(evt.epsilon, evt.iter_weight)


def process_hist_data(
    tree_data,
    tree_dummy,
    t_bins,
    phi_bins,
    nWindows,
    phi_setting,
    inpDict,
    event_cache=None,
    sub_event_cache=None,
    particle_subtraction_scale_factor=None,
    kaon_signal_shape_payload=None,
    kaon_sigma0_shape_payload=None,
):

    processed_dict = {}

    OutFilename = inpDict["OutFilename"]    
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"]
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]

    norm_factor_data = inpDict["normfac_data"]
    norm_factor_dummy = inpDict["normfac_dummy"]    

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import evaluate_data_event, get_shifted_mm, get_shifted_t, set_shift_context, set_val
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=inpDict.get("shift_mode", "raw"))

    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, phi_setting)
    
    ################################################################################################################################################
        
    prompt_tree_name = get_prompt_tree_name(ParticleType, EPSSET)
    rand_tree_name = get_rand_tree_name(ParticleType, EPSSET)

    TBRANCH_DATA  = tree_data.Get(prompt_tree_name)
    TBRANCH_RAND  = tree_data.Get(rand_tree_name)
    
    TBRANCH_DUMMY  = tree_dummy.Get(prompt_tree_name)
    TBRANCH_DUMMY_RAND  = tree_dummy.Get(rand_tree_name)

    hist_bin_dict = {}
    n_t = len(t_bins) - 1
    particle_subtraction_mode = resolve_particle_subtraction_mode(inpDict)
    component_shape_payload = None
    component_fit_results = [None for _ in range(n_t)]
    component_subtraction_payloads = [None for _ in range(n_t)]
    data_hists = _init_ave_hist_group(("Q2", "W", "t", "epsilon", "mm", "fit1sub", "pisub", "nosub"), n_t)
    rand_hists = _init_ave_hist_group(("Q2", "W", "t", "epsilon", "mm", "fit1sub", "pisub", "nosub"), n_t)
    dummy_hists = _init_ave_hist_group(("Q2", "W", "t", "epsilon", "mm", "fit1sub", "pisub", "nosub"), n_t)
    dummy_rand_hists = _init_ave_hist_group(("Q2", "W", "t", "epsilon", "mm", "fit1sub", "pisub", "nosub"), n_t)

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        from particle_subtraction import particle_subtraction_ave
        SubtractedParticle = "pion"
        subDict = {}
        fitDict = {}
        if particle_subtraction_mode == "simc_shape_components":
            component_shape_payload = load_setting_pion_component_shapes(
                inpDict,
                phi_setting,
                particle_type=ParticleType,
                t_bins=t_bins,
                phi_bins=phi_bins,
                hgcer_cutg=hgcer_cutg,
                context="ave_per_bin_scope_fits",
            )

    # Fit background and subtract
    from background_fit import bg_fit
        
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):

        hist_bin_dict["H_Q2_DATA_{}".format(j)]       = TH1D("H_Q2_DATA_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DATA_{}".format(j)]  = TH1D("H_W_DATA_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DATA_{}".format(j)]       = TH1D("H_t_DATA_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DATA_{}".format(j)]  = TH1D("H_epsilon_DATA_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DATA_{}".format(j)]       = TH1D("H_MM_DATA_{}".format(j),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)]       = TH1D("H_MM_fit1sub_DATA_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_pisub_DATA_{}".format(j)]       = TH1D("H_MM_pisub_DATA_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)]       = TH1D("H_MM_nosub_DATA_{}".format(j),"MM", 100, 0.7, 1.5)

        hist_bin_dict["H_Q2_RAND_{}".format(j)]       = TH1D("H_Q2_RAND_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_RAND_{}".format(j)]  = TH1D("H_W_RAND_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_RAND_{}".format(j)]       = TH1D("H_t_RAND_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_RAND_{}".format(j)]  = TH1D("H_epsilon_RAND_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_RAND_{}".format(j)]       = TH1D("H_MM_RAND_{}".format(j),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_fit1sub_RAND_{}".format(j)]       = TH1D("H_MM_fit1sub_RAND_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_pisub_RAND_{}".format(j)]       = TH1D("H_MM_pisub_RAND_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_nosub_RAND_{}".format(j)]       = TH1D("H_MM_nosub_RAND_{}".format(j),"MM", 100, 0.7, 1.5)

        hist_bin_dict["H_Q2_DUMMY_{}".format(j)]       = TH1D("H_Q2_DUMMY_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DUMMY_{}".format(j)]  = TH1D("H_W_DUMMY_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DUMMY_{}".format(j)]       = TH1D("H_t_DUMMY_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DUMMY_{}".format(j)]  = TH1D("H_epsilon_DUMMY_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DUMMY_{}".format(j)]       = TH1D("H_MM_DUMMY_{}".format(j),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_fit1sub_DUMMY_{}".format(j)]       = TH1D("H_MM_fit1sub_DUMMY_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_pisub_DUMMY_{}".format(j)]       = TH1D("H_MM_pisub_DUMMY_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)]       = TH1D("H_MM_nosub_DUMMY_{}".format(j),"MM", 100, 0.7, 1.5)

        hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)]       = TH1D("H_Q2_DUMMY_RAND_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)]  = TH1D("H_W_DUMMY_RAND_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)]       = TH1D("H_t_DUMMY_RAND_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)]  = TH1D("H_epsilon_DUMMY_RAND_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_DUMMY_RAND_{}".format(j),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_fit1sub_DUMMY_RAND_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_pisub_DUMMY_RAND_{}".format(j),"MM", 100, 0.7, 1.5)
        hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_nosub_DUMMY_RAND_{}".format(j),"MM", 100, 0.7, 1.5)

        data_hists["Q2"][j] = hist_bin_dict["H_Q2_DATA_{}".format(j)]
        data_hists["W"][j] = hist_bin_dict["H_W_DATA_{}".format(j)]
        data_hists["t"][j] = hist_bin_dict["H_t_DATA_{}".format(j)]
        data_hists["epsilon"][j] = hist_bin_dict["H_epsilon_DATA_{}".format(j)]
        data_hists["mm"][j] = hist_bin_dict["H_MM_DATA_{}".format(j)]
        data_hists["fit1sub"][j] = hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)]
        data_hists["pisub"][j] = hist_bin_dict["H_MM_pisub_DATA_{}".format(j)]
        data_hists["nosub"][j] = hist_bin_dict["H_MM_nosub_DATA_{}".format(j)]

        rand_hists["Q2"][j] = hist_bin_dict["H_Q2_RAND_{}".format(j)]
        rand_hists["W"][j] = hist_bin_dict["H_W_RAND_{}".format(j)]
        rand_hists["t"][j] = hist_bin_dict["H_t_RAND_{}".format(j)]
        rand_hists["epsilon"][j] = hist_bin_dict["H_epsilon_RAND_{}".format(j)]
        rand_hists["mm"][j] = hist_bin_dict["H_MM_RAND_{}".format(j)]
        rand_hists["fit1sub"][j] = hist_bin_dict["H_MM_fit1sub_RAND_{}".format(j)]
        rand_hists["pisub"][j] = hist_bin_dict["H_MM_pisub_RAND_{}".format(j)]
        rand_hists["nosub"][j] = hist_bin_dict["H_MM_nosub_RAND_{}".format(j)]

        dummy_hists["Q2"][j] = hist_bin_dict["H_Q2_DUMMY_{}".format(j)]
        dummy_hists["W"][j] = hist_bin_dict["H_W_DUMMY_{}".format(j)]
        dummy_hists["t"][j] = hist_bin_dict["H_t_DUMMY_{}".format(j)]
        dummy_hists["epsilon"][j] = hist_bin_dict["H_epsilon_DUMMY_{}".format(j)]
        dummy_hists["mm"][j] = hist_bin_dict["H_MM_DUMMY_{}".format(j)]
        dummy_hists["fit1sub"][j] = hist_bin_dict["H_MM_fit1sub_DUMMY_{}".format(j)]
        dummy_hists["pisub"][j] = hist_bin_dict["H_MM_pisub_DUMMY_{}".format(j)]
        dummy_hists["nosub"][j] = hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)]

        dummy_rand_hists["Q2"][j] = hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["W"][j] = hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["t"][j] = hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["epsilon"][j] = hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["mm"][j] = hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["fit1sub"][j] = hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["pisub"][j] = hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}".format(j)]
        dummy_rand_hists["nosub"][j] = hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)]

        # Pion subtraction by scaling simc to peak size
        if ParticleType == "kaon":
            
            subDict["H_Q2_SUB_DATA_{}".format(j)]       = TH1D("H_Q2_SUB_DATA_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DATA_{}".format(j)]  = TH1D("H_W_SUB_DATA_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DATA_{}".format(j)]       = TH1D("H_t_SUB_DATA_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DATA_{}".format(j)]  = TH1D("H_epsilon_SUB_DATA_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DATA_{}".format(j)]  = TH1D("H_MM_SUB_DATA_{}".format(j),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DATA_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DATA_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)

            subDict["H_Q2_SUB_RAND_{}".format(j)]       = TH1D("H_Q2_SUB_RAND_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_RAND_{}".format(j)]  = TH1D("H_W_SUB_RAND_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_RAND_{}".format(j)]       = TH1D("H_t_SUB_RAND_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_RAND_{}".format(j)]  = TH1D("H_epsilon_SUB_RAND_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_RAND_{}".format(j)]  = TH1D("H_MM_SUB_RAND_{}".format(j),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_RAND_{}".format(j)]  = TH1D("H_MM_nosub_SUB_RAND_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)

            subDict["H_Q2_SUB_DUMMY_{}".format(j)]       = TH1D("H_Q2_SUB_DUMMY_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DUMMY_{}".format(j)]  = TH1D("H_W_SUB_DUMMY_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DUMMY_{}".format(j)]       = TH1D("H_t_SUB_DUMMY_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DUMMY_{}".format(j)]  = TH1D("H_epsilon_SUB_DUMMY_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DUMMY_{}".format(j)]  = TH1D("H_MM_SUB_DUMMY_{}".format(j),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DUMMY_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DUMMY_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)

            subDict["H_Q2_SUB_DUMMY_RAND_{}".format(j)]       = TH1D("H_Q2_SUB_DUMMY_RAND_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_W_SUB_DUMMY_RAND_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DUMMY_RAND_{}".format(j)]       = TH1D("H_t_SUB_DUMMY_RAND_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_epsilon_SUB_DUMMY_RAND_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_MM_SUB_DUMMY_RAND_{}".format(j),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DUMMY_RAND_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)
            
    hole_contains = hgcer_cutg.IsInside if ParticleType == "kaon" else None

    print("\nBinning data...")
    if event_cache is not None:
        mm_offset_from_cache = _fill_cached_ave_section(event_cache["prompt"], data_hists, t_bins, update_mm_offset=True)
        if mm_offset_from_cache is not None:
            MM_offset_DATA = mm_offset_from_cache
    else:
        mm_offset_from_tree = _process_ave_data_tree(
            TBRANCH_DATA,
            data_hists,
            get_shifted_mm,
            get_shifted_t,
            t_bins,
            ParticleType,
            hole_contains,
            evaluate_data_event,
            mm_min,
            mm_max,
            Misc.progressBar,
            update_mm_offset=True,
        )
        if mm_offset_from_tree is not None:
            MM_offset_DATA = mm_offset_from_tree

    print("\nBinning dummy...")
    if event_cache is not None:
        _fill_cached_ave_section(event_cache["dummy"], dummy_hists, t_bins)
    else:
        _process_ave_data_tree(
            TBRANCH_DUMMY,
            dummy_hists,
            get_shifted_mm,
            get_shifted_t,
            t_bins,
            ParticleType,
            hole_contains,
            evaluate_data_event,
            mm_min,
            mm_max,
            Misc.progressBar,
        )

    print("\nBinning rand...")
    if event_cache is not None:
        _fill_cached_ave_section(event_cache["rand"], rand_hists, t_bins)
    else:
        _process_ave_data_tree(
            TBRANCH_RAND,
            rand_hists,
            get_shifted_mm,
            get_shifted_t,
            t_bins,
            ParticleType,
            hole_contains,
            evaluate_data_event,
            mm_min,
            mm_max,
            Misc.progressBar,
        )

    print("\nBinning dummy_rand...")
    if event_cache is not None:
        _fill_cached_ave_section(event_cache["dummy_rand"], dummy_rand_hists, t_bins)
    else:
        _process_ave_data_tree(
            TBRANCH_DUMMY_RAND,
            dummy_rand_hists,
            get_shifted_mm,
            get_shifted_t,
            t_bins,
            ParticleType,
            hole_contains,
            evaluate_data_event,
            mm_min,
            mm_max,
            Misc.progressBar,
        )

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        subDict["MM_offset_DATA"] = MM_offset_DATA
        particle_subtraction_ave(t_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg)
                    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
                    
        hist_bin_dict["H_Q2_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_W_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_bin_dict["H_t_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_epsilon_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_fit1sub_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_pisub_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_nosub_RAND_{}".format(j)].Scale(1/nWindows)        

        hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(hist_bin_dict["H_Q2_RAND_{}".format(j)],-1)
        hist_bin_dict["H_W_DATA_{}".format(j)].Add(hist_bin_dict["H_W_RAND_{}".format(j)],-1)
        hist_bin_dict["H_t_DATA_{}".format(j)].Add(hist_bin_dict["H_t_RAND_{}".format(j)],-1)
        hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(hist_bin_dict["H_epsilon_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_fit1sub_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_pisub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_pisub_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_nosub_RAND_{}".format(j)],-1)        

        hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows) 
        hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)        

        hist_bin_dict["H_Q2_DUMMY_{}".format(j)].Add(hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_W_DUMMY_{}".format(j)].Add(hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_t_DUMMY_{}".format(j)].Add(hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_epsilon_DUMMY_{}".format(j)].Add(hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_fit1sub_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_pisub_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)],-1)   
        
        # Data Normalization
        hist_bin_dict["H_Q2_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_W_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_t_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_epsilon_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_MM_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_MM_pisub_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].Scale(norm_factor_data)

        # Dummy Normalization
        hist_bin_dict["H_Q2_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_W_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_t_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_epsilon_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_MM_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_MM_fit1sub_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_MM_pisub_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)].Scale(norm_factor_dummy)   

        # Dummy subtraction
        hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(hist_bin_dict["H_Q2_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_W_DATA_{}".format(j)].Add(hist_bin_dict["H_W_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_t_DATA_{}".format(j)].Add(hist_bin_dict["H_t_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(hist_bin_dict["H_epsilon_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_MM_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_fit1sub_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_MM_pisub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_pisub_DUMMY_{}".format(j)], -1)
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)], -1)  

        # Remove histograms with less than event_threshold entries and negative integrals
        event_threshold = 1
        prune_hist(
            hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)],
            event_threshold
        )     
        prune_hist(
            hist_bin_dict["H_MM_pisub_DATA_{}".format(j)],
            event_threshold
        )  
        prune_hist(
            hist_bin_dict["H_MM_nosub_DATA_{}".format(j)],
            event_threshold
        )                     
        prune_hist(
            hist_bin_dict["H_MM_DATA_{}".format(j)],
            event_threshold
        )
        prune_hist(
            hist_bin_dict["H_t_DATA_{}".format(j)],
            event_threshold
        )

        # Pion subtraction by scaling simc to peak size
        if ParticleType == "kaon":
            scale_factor = 0.0
            component_payload = None
            use_legacy_scalar_subtraction = True
            if component_shape_payload is not None:
                component_fit_results[j] = build_particle_subtraction_component_result(
                    subDict["H_MM_nosub_SUB_DATA_{}".format(j)],
                    hist_bin_dict["H_MM_nosub_DATA_{}".format(j)],
                    resolve_scope_component_shapes(
                        component_shape_payload,
                        analysis_scope="t_bin{}".format(j + 1),
                        t_bin_index=j,
                    ),
                        inpDict,
                        analysis_scope="t_bin{}".format(j + 1),
                        kaon_signal_shape=resolve_scope_single_shape(
                            kaon_signal_shape_payload,
                            analysis_scope="t_bin{}".format(j + 1),
                            t_bin_index=j,
                        ),
                        kaon_sigma0_shape=resolve_scope_single_shape(
                            kaon_sigma0_shape_payload,
                            analysis_scope="t_bin{}".format(j + 1),
                            t_bin_index=j,
                        ),
                        mm_offset_data=MM_offset_DATA,
                        context="ave_{}_t{}".format(phi_setting, j + 1),
                    )
                component_payload = _apply_component_pion_subtraction_for_tbin(
                    hist_bin_dict,
                    j,
                    component_fit_results[j],
                    sub_event_cache,
                    norm_factor_data,
                    norm_factor_dummy,
                    nWindows,
                    ParticleType,
                    inpDict["POL"],
                    {**inpDict, "phi_setting": phi_setting},
                )
                component_subtraction_payloads[j] = component_payload
                if component_payload.get("accepted"):
                    use_legacy_scalar_subtraction = False
                    scale_factor = float(component_payload.get("particle_subtraction_effective_scale", 0.0) or 0.0)
                else:
                    fallback_mode = component_payload.get("fallback_mode") or "single_scale"
                    if fallback_mode == "error":
                        raise RuntimeError(
                            "ave_per_bin component pion subtraction ({}, t{}) rejected: {}".format(
                                phi_setting,
                                int(j) + 1,
                                component_payload.get("fallback_reason") or "unknown reason",
                            )
                        )
                    if fallback_mode == "single_scale":
                        use_legacy_scalar_subtraction = True
                    elif fallback_mode in ("zero", "skip_bin"):
                        use_legacy_scalar_subtraction = False
                        scale_factor = 0.0
                    else:
                        use_legacy_scalar_subtraction = True

            if use_legacy_scalar_subtraction:
                if particle_subtraction_scale_factor is not None:
                    scale_factor = float(particle_subtraction_scale_factor)
                else:
                    try:
                        subtraction_windows = resolve_particle_subtraction_windows(
                            ParticleType,
                            SubtractedParticle,
                            MM_offset_DATA,
                        )
                        scale_components = compute_staged_particle_subtraction_scales(
                            hist_bin_dict[f"H_MM_nosub_DATA_{j}"],
                            subDict[f"H_MM_nosub_SUB_DATA_{j}"],
                            subtraction_windows,
                            context="pion subtraction (t-bin {})".format(j),
                        )
                        scale_factor = scale_components["total_scale_factor"]
                    except ZeroDivisionError:
                        scale_factor = 0.0
                '''
                if scale_factor > 10.0:
                    print("\n\nWARNING: Pion scaling factor too large, likely no pion peak. Setting to zero....")
                    scale_factor = 0.0
                '''

                if phi_setting == "Center":
                    phi_scale = 0.95
                elif phi_setting == "Left":
                    phi_scale = 0.65
                elif phi_setting == "Right":
                    phi_scale = 0.65
                else:
                    raise ValueError("Invalid phi_setting: {}".format(phi_setting))

                scale_factor = scale_factor #* phi_scale

                # Scale pion to subtraction proper peak
                subDict["H_Q2_SUB_DATA_{}".format(j)].Scale(scale_factor)
                subDict["H_W_SUB_DATA_{}".format(j)].Scale(scale_factor)
                subDict["H_t_SUB_DATA_{}".format(j)].Scale(scale_factor)
                subDict["H_epsilon_SUB_DATA_{}".format(j)].Scale(scale_factor)
                subDict["H_MM_SUB_DATA_{}".format(j)].Scale(scale_factor)
                subDict["H_MM_nosub_SUB_DATA_{}".format(j)].Scale(scale_factor)

                # Apply pion subtraction
                hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(subDict["H_Q2_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_W_DATA_{}".format(j)].Add(subDict["H_W_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_t_DATA_{}".format(j)].Add(subDict["H_t_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(subDict["H_epsilon_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Add(subDict["H_MM_nosub_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_MM_pisub_DATA_{}".format(j)].Add(subDict["H_MM_nosub_SUB_DATA_{}".format(j)],-1)
                hist_bin_dict["H_MM_DATA_{}".format(j)].Add(subDict["H_MM_SUB_DATA_{}".format(j)],-1)

        # Fit background and subtract
        # ---- Statistic‑scale for this (t,φ) bin ----------------
        inpDict["bg_stat_scale1"] = resolve_bg_stat_scale1(inpDict, phi_setting)
        residual_bg_weights1 = None
        
        if inpDict["bg_stat_scale1"] > 0.0:
            fitDict["background_data_fit1_{}".format(j)] = bg_fit(
                phi_setting,
                inpDict,
                hist_bin_dict[f"H_MM_pisub_DATA_{j}"],   # wide / no–MM-cut
                hist_bin_dict[f"H_MM_DATA_{j}"],          # cut-window axis
                scaling=inpDict["bg_stat_scale1"],
                model_key=f"fixquad_{phi_setting}_{EPSSET}e",
                fit_name="Fit 1"
            )
            # ----------------------------------------------------------------

            can_subtract_with_templates = event_cache is not None and (ParticleType != "kaon" or sub_event_cache is not None)
            if can_subtract_with_templates:
                mm_stage1_input = hist_bin_dict["H_MM_DATA_{}".format(j)].Clone("{}_stage1_input".format(hist_bin_dict["H_MM_DATA_{}".format(j)].GetName()))
                if hasattr(mm_stage1_input, "SetDirectory"):
                    mm_stage1_input.SetDirectory(0)
                active_component_payload = component_subtraction_payloads[j]
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                bg_weights1 = _subtract_ave_mm_background_for_tbin(
                    hist_bin_dict,
                    j,
                    mm_stage1_input,
                    fitDict["background_data_fit1_{}".format(j)][0],
                    event_cache,
                    sub_event_cache,
                    norm_factor_data,
                    norm_factor_dummy,
                    nWindows,
                    scale_factor if ParticleType == "kaon" else 0.0,
                    t_bins,
                    pion_component_payload=active_component_payload,
                )
                residual_bg_weights1 = build_mm_residual_weights(bg_weights1)
            else:
                raise RuntimeError(
                    "Missing cached yield-event templates for ave_per_bin stage-1 MM background subtraction. "
                    "Legacy histogram rescaling fallback is disabled to preserve physics normalization."
                )
            hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)].Add(fitDict["background_data_fit1_{}".format(j)][1], -1)
            hist_bin_dict["H_MM_DATA_{}".format(j)].Add(fitDict["background_data_fit1_{}".format(j)][0], -1)  

            # Remove histograms with less than event_threshold entries and negative integrals  
            prune_hist(
                hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)],
                event_threshold
            ) 
            prune_hist(
                hist_bin_dict["H_MM_pisub_DATA_{}".format(j)],
                event_threshold
            )  
            prune_hist(
                hist_bin_dict["H_MM_nosub_DATA_{}".format(j)],
                event_threshold
            )                  
            prune_hist(
                hist_bin_dict["H_MM_DATA_{}".format(j)],
                event_threshold
            )
            prune_hist(
                hist_bin_dict["H_t_DATA_{}".format(j)],
                event_threshold
            )

        # Fit background and subtract
        # ---- Statistic‑scale for this (t,φ) bin ----------------
        inpDict["bg_stat_scale2"] = resolve_bg_stat_scale2(inpDict, phi_setting)
        
        if inpDict["bg_stat_scale2"] > 0.0:  
            fitDict["background_data_fit2_{}".format(j)] = bg_fit(
                phi_setting,
                inpDict,
                hist_bin_dict[f"H_MM_fit1sub_DATA_{j}"],   # wide / no–MM-cut
                hist_bin_dict[f"H_MM_DATA_{j}"],          # cut-window axis
                scaling=inpDict["bg_stat_scale2"],
                model_key=f"cheb2_{phi_setting}_{EPSSET}e",
                fit_name="Fit 2"
            )
            # ----------------------------------------------------------------

            can_subtract_with_templates = event_cache is not None and (ParticleType != "kaon" or sub_event_cache is not None)
            if can_subtract_with_templates:
                mm_stage2_input = hist_bin_dict["H_MM_DATA_{}".format(j)].Clone("{}_stage2_input".format(hist_bin_dict["H_MM_DATA_{}".format(j)].GetName()))
                if hasattr(mm_stage2_input, "SetDirectory"):
                    mm_stage2_input.SetDirectory(0)
                active_component_payload = component_subtraction_payloads[j]
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                _subtract_ave_mm_background_for_tbin(
                    hist_bin_dict,
                    j,
                    mm_stage2_input,
                    fitDict["background_data_fit2_{}".format(j)][0],
                    event_cache,
                    sub_event_cache,
                    norm_factor_data,
                    norm_factor_dummy,
                    nWindows,
                    scale_factor if ParticleType == "kaon" else 0.0,
                    t_bins,
                    pion_component_payload=active_component_payload,
                    residual_weights=residual_bg_weights1,
                )
            else:
                raise RuntimeError(
                    "Missing cached yield-event templates for ave_per_bin stage-2 MM background subtraction. "
                    "Legacy histogram rescaling fallback is disabled to preserve physics normalization."
                )
            hist_bin_dict["H_MM_DATA_{}".format(j)].Add(fitDict["background_data_fit2_{}".format(j)][0], -1)  

            # Remove histograms with less than event_threshold entries and negative integrals
            prune_hist(
                hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)],
                event_threshold
            )          
            prune_hist(
                hist_bin_dict["H_MM_pisub_DATA_{}".format(j)],
                event_threshold
            )          
            prune_hist(
                hist_bin_dict["H_MM_nosub_DATA_{}".format(j)],
                event_threshold
            )          
            prune_hist(
                hist_bin_dict["H_MM_DATA_{}".format(j)],
                event_threshold
            )
            prune_hist(
                hist_bin_dict["H_t_DATA_{}".format(j)],
                event_threshold
            )        

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):

        processed_dict["t_bin{}".format(j+1)] = {
            "H_Q2_DATA" : hist_bin_dict["H_Q2_DATA_{}".format(j)],
            "H_W_DATA" : hist_bin_dict["H_W_DATA_{}".format(j)],
            "H_t_DATA" : hist_bin_dict["H_t_DATA_{}".format(j)],
            "H_epsilon_DATA" : hist_bin_dict["H_epsilon_DATA_{}".format(j)],
            "H_MM_DATA" : hist_bin_dict["H_MM_DATA_{}".format(j)],
            "H_Q2_DUMMY" : hist_bin_dict["H_Q2_DUMMY_{}".format(j)],
            "H_W_DUMMY" : hist_bin_dict["H_W_DUMMY_{}".format(j)],
            "H_t_DUMMY" : hist_bin_dict["H_t_DUMMY_{}".format(j)],
            "H_epsilon_DUMMY" : hist_bin_dict["H_epsilon_DUMMY_{}".format(j)],
            "H_MM_DUMMY" : hist_bin_dict["H_MM_DUMMY_{}".format(j)],
            "H_MM_pisub_DATA" : hist_bin_dict["H_MM_pisub_DATA_{}".format(j)],
            "H_MM_fit1sub_DATA" : hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)],
            "H_MM_fit1sub_DATA" : hist_bin_dict["H_MM_fit1sub_DATA_{}".format(j)],
            "particle_subtraction_component_fit" : component_fit_results[j],
            "particle_subtraction_component_fit_summary" : serialize_particle_subtraction_component_result(
                component_fit_results[j]
            ) if component_fit_results[j] is not None else None,
            "particle_subtraction_component_payload" : component_subtraction_payloads[j],
            "particle_subtraction_component_payload_summary" : summarize_particle_subtraction_component_payload(
                component_subtraction_payloads[j]
            ) if component_subtraction_payloads[j] is not None else None,
            "particle_subtraction_component_applied" : bool(
                isinstance(component_subtraction_payloads[j], dict)
                and component_subtraction_payloads[j].get("accepted")
            ),
        }

        # Sort dictionary keys alphabetically
        processed_dict["t_bin{}".format(j+1)] = {key : processed_dict["t_bin{}".format(j+1)][key] \
                                                 for key in sorted(processed_dict["t_bin{}".format(j+1)].keys())}

        # Include Stat box
        ROOT.gStyle.SetOptStat(1)
        for i, (key,val) in enumerate(processed_dict["t_bin{}".format(j+1)].items()):
            if not hasattr(val, "Draw"):
                continue
            canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
            centroid = get_centroid(val, val.GetXaxis().GetXmin(), val.GetXaxis().GetXmax()) # [centroid, centroid_error]
            val.Draw()
            val.SetTitle("{} {}".format(phi_setting, val.GetName()))
            # Create a TLatex object to add text to the plot
            text = ROOT.TLatex()
            text.SetNDC();
            text.SetTextSize(0.04);
            text.SetTextAlign(22); # Centered alignment
            text.SetTextColor(ROOT.kRed); # Set text color to red            
            # Add the centroid value to the plot
            text.DrawLatex(0.7, 0.65, "Centroid: {:.2f}+/-{:.3f}".format(centroid[0], centroid[1]))
            if i==0 and j==0:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType))+'(')
            elif i==len(processed_dict["t_bin{}".format(j+1)].items())-1 and j==len(t_bins)-2:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType))+')')
            else:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType)))

            if key == "H_MM_DATA" and component_fit_results[j] is not None:
                print_particle_subtraction_component_fit_pages(
                    outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType)),
                    component_fit_results[j],
                    title_prefix="{} t{}".format(phi_setting, j + 1),
                    cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
                )
                if isinstance(component_subtraction_payloads[j], dict):
                    print_particle_subtraction_component_application_pages(
                        outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType)),
                        component_subtraction_payloads[j],
                        title_prefix="{} t{}".format(phi_setting, j + 1),
                        cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
                    )
            del canvas
            
    return processed_dict

def bin_data(
    kinematic_types,
    tree_data,
    tree_dummy,
    t_bins,
    phi_bins,
    nWindows,
    phi_setting,
    inpDict,
    event_cache=None,
    sub_event_cache=None,
    particle_subtraction_scale_factor=None,
    kaon_signal_shape_payload=None,
    kaon_sigma0_shape_payload=None,
):

    processed_dict = process_hist_data(
        tree_data,
        tree_dummy,
        t_bins,
        phi_bins,
        nWindows,
        phi_setting,
        inpDict,
        event_cache=event_cache,
        sub_event_cache=sub_event_cache,
        particle_subtraction_scale_factor=particle_subtraction_scale_factor,
        kaon_signal_shape_payload=kaon_signal_shape_payload,
        kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
    )
    
    binned_dict = {}

    for kin_type in kinematic_types:

        # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
        binned_t_data = []
        binned_hist_data = []
        binned_hist_dummy = []
        kin_hist_data = []
        kin_hist_dummy = []
        
        # Loop through bins in t_data and identify events in specified bins
        for j in range(len(t_bins)-1):

            H_Q2_DATA = processed_dict["t_bin{}".format(j+1)]["H_Q2_DATA"]
            H_W_DATA = processed_dict["t_bin{}".format(j+1)]["H_W_DATA"]
            H_t_DATA = processed_dict["t_bin{}".format(j+1)]["H_t_DATA"]
            H_epsilon_DATA = processed_dict["t_bin{}".format(j+1)]["H_epsilon_DATA"]

            H_Q2_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_Q2_DUMMY"]
            H_W_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_W_DUMMY"]
            H_t_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_t_DUMMY"]
            H_epsilon_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_epsilon_DUMMY"]

            # Initialize lists for tmp_binned_t_data, tmp_binned_hist_data, and tmp_binned_hist_dummy
            tmp_binned_t_data = []
            tmp_binned_hist_data = []
            tmp_binned_hist_dummy = []

            tmp_hist_data = [[],[]]
            for i in range(1, H_t_DATA.GetNbinsX() + 1):
                tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                
            tmp_binned_t_data.append(tmp_hist_data)

            if kin_type == "t":
                kin_hist_data.append(H_t_DATA.Clone())
                tmp_hist_data = [[],[]]                
                for i in range(1, H_t_DATA.GetNbinsX() + 1):        
                    tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)
            if kin_type == "Q2":
                kin_hist_data.append(H_Q2_DATA.Clone())
                tmp_hist_data = [[],[]]                
                for i in range(1, H_Q2_DATA.GetNbinsX() + 1):        
                    tmp_hist_data[0].append(H_Q2_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_Q2_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)
            if kin_type == "W":
                kin_hist_data.append(H_W_DATA.Clone())            
                tmp_hist_data = [[],[]]                
                for i in range(1, H_W_DATA.GetNbinsX() + 1):
                    tmp_hist_data[0].append(H_W_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_W_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)        
            if kin_type == "epsilon":
                kin_hist_data.append(H_epsilon_DATA.Clone())
                tmp_hist_data = [[],[]]                
                for i in range(1, H_epsilon_DATA.GetNbinsX() + 1):
                    tmp_hist_data[0].append(H_epsilon_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_epsilon_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)

            if kin_type == "t":
                kin_hist_dummy.append(H_t_DUMMY.Clone())
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_t_DUMMY.GetNbinsX() + 1):        
                    tmp_hist_dummy[0].append(H_t_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_t_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)
            if kin_type == "Q2":
                kin_hist_dummy.append(H_Q2_DUMMY.Clone())
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_Q2_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_Q2_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_Q2_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)
            if kin_type == "W":
                kin_hist_dummy.append(H_W_DUMMY.Clone())
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_W_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_W_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_W_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)        
            if kin_type == "epsilon":
                kin_hist_dummy.append(H_epsilon_DUMMY.Clone())
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_epsilon_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_epsilon_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_epsilon_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)

            binned_t_data.append(tmp_binned_t_data[0]) # Save a list of hists where each one is a t-bin
            binned_hist_data.append(tmp_binned_hist_data[0])
            binned_hist_dummy.append(tmp_binned_hist_dummy[0])

        binned_dict[kin_type] = {
            "binned_t_data" : binned_t_data,
            "binned_hist_data" : binned_hist_data,
            "binned_hist_dummy" : binned_hist_dummy,
            "kin_hist_data" : kin_hist_data,
            "kin_hist_dummy" : kin_hist_dummy
        }
        
    return binned_dict
    
def calculate_ave_data(kinematic_types, hist, t_bins, phi_bins, inpDict):

    tree_data, tree_dummy = hist["InFile_DATA"], hist["InFile_DUMMY"]
    nWindows, phi_setting = hist["nWindows"], hist["phi_setting"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]    
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    event_cache = hist.get("_yield_data_event_cache")
    sub_event_cache = hist.get("_yield_sub_event_cache")
    particle_subtraction_scale_factor = hist.get("particle_subtraction_scale_factor")
    if resolve_particle_subtraction_mode(inpDict) == "simc_shape_components":
        particle_subtraction_scale_factor = None
    kaon_signal_shape_payload = _get_cached_kaon_signal_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "ave_per_bin_signal",
    )
    kaon_sigma0_shape_payload = _get_cached_kaon_sigma0_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "ave_per_bin_sigma0",
    )
    binned_dict = bin_data(
        kinematic_types,
        tree_data,
        tree_dummy,
        t_bins,
        phi_bins,
        nWindows,
        phi_setting,
        inpDict,
        event_cache=event_cache,
        sub_event_cache=sub_event_cache,
        particle_subtraction_scale_factor=particle_subtraction_scale_factor,
        kaon_signal_shape_payload=kaon_signal_shape_payload,
        kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
    )

    group_dict = {}
    bad_bins = set()  # <-- collect (tbin_index, phibin_index) that have sem==1000.0 anywhere
    
    for kin_type in kinematic_types:
        binned_t_data = binned_dict[kin_type]["binned_t_data"]
        binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
        binned_hist_dummy = binned_dict[kin_type]["binned_hist_dummy"]
        kin_hist_data = binned_dict[kin_type]["kin_hist_data"]
        kin_hist_dummy = binned_dict[kin_type]["kin_hist_dummy"]
        
        ave_hist = []
        ave_err_hist = []
        binned_sub_data = [[],[]]
        i=0 # iter
        print("-"*25)
        # Subtract binned_hist_dummy from binned_hist_data element-wise
        for data, dummy in zip(binned_hist_data, binned_hist_dummy):
            bin_val_data, hist_val_data = data
            bin_val_dummy, hist_val_dummy = dummy
            arr_data = np.array(hist_val_data)
            try:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(arr_data * bin_val_data)
                total_count = np.sum(arr_data)
                average = weighted_sum / total_count
                if math.isnan(average) or math.isinf(average) or average <= 0.0:
                    print("Empty binning for data {} (t-bin={})... ".format(kin_type, i+1))
                    average = 0.0
                ave_hist.append(average)

                x = np.array(bin_val_data)

                # Weighted variance using the same weights (arr_data)
                # var = sum(w * (x - mu)^2) / sum(w)
                if total_count > 0:
                    var = np.sum(arr_data * (x - average)**2) / total_count
                    std_dev = np.sqrt(var)

                    # Use total_count (number of entries) as N, not number of bins
                    n = total_count
                    sem = std_dev / np.sqrt(n)
                else:
                    std_dev = 0.0
                    sem = -1000.0

                if average == 0.0:
                    sem = -1000.0 # Assign a large error if average is zero

                ave_err_hist.append(sem)
                binned_sub_data[0].append(bin_val_data)
                binned_sub_data[1].append(arr_data)
            except ZeroDivisionError:
                print("Empty binning for data {} (t-bin={})... ".format(kin_type, i+1))
                sys.exit(2)
            i+=1

        # Print statements to check sizes
        #print("Size of binned_t_data:", len(binned_t_data))
        #print("Size of binned_hist_data:", len(binned_hist_data))
        #print("Size of binned_hist_dummy:", len(binned_hist_dummy))
        #print("Size of binned_sub_data:", len(binned_sub_data[1]))
        #print("Size of ave_hist:", len(ave_hist))
        #print("Size of t_bins:", len(t_bins))
        #print("Size of phi_bins:", len(phi_bins), "\n")

        dict_lst = []
        for j in range(len(t_bins) - 1):
            tbin_index = j
            for k in range(len(phi_bins) - 1):
                phibin_index = k
                hist_val = [binned_sub_data[0][j], binned_sub_data[1][j]]
                ave_val = ave_hist[j]
                ave_err_val = ave_err_hist[j]
                print("Data average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, j+1, k+1, ave_val, ave_err_val))
                dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))

                # mark bad bins for global override later
                if ave_err_val == -1000.0:
                    bad_bins.add((tbin_index, phibin_index))                
                
        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_ave".format(kin_type) : tup[2],
                "{}_ave_err".format(kin_type) : tup[3],
            }

        group_dict[kin_type] = groups 
    
    return group_dict

def ave_per_bin_data(histlist, inpDict):

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data averages for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        binned_dict = calculate_ave_data(kinematic_types, hist, t_bins, phi_bins, inpDict)
        hist.pop("_yield_data_event_cache", None)
        hist.pop("_yield_sub_event_cache", None)
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
                
    return {"binned_DATA" : aveDict}

##################################################################################################################################################

def process_hist_simc(tree_simc, t_bins, phi_setting, inpDict, iteration, event_cache=None):

    processed_dict = {}

    OutFilename = inpDict["OutFilename"]    
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]    
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]    

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_simc_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, phi_setting)
    
    ################################################################################################################################################
        
    TBRANCH_SIMC  = tree_simc.Get("h10")
    
    hist_bin_dict = {}
    n_t = len(t_bins) - 1
    simc_hists = _init_ave_hist_group(("Q2", "W", "t", "epsilon"), n_t)
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):

        hist_bin_dict["H_Q2_SIMC_{}".format(j)]       = TH1D("H_Q2_SIMC_{}".format(j),"Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_SIMC_{}".format(j)]  = TH1D("H_W_SIMC_{}".format(j),"W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_SIMC_{}".format(j)]       = TH1D("H_t_SIMC_{}".format(j),"-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_SIMC_{}".format(j)]  = TH1D("H_epsilon_SIMC_{}".format(j),"epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        simc_hists["Q2"][j] = hist_bin_dict["H_Q2_SIMC_{}".format(j)]
        simc_hists["W"][j] = hist_bin_dict["H_W_SIMC_{}".format(j)]
        simc_hists["t"][j] = hist_bin_dict["H_t_SIMC_{}".format(j)]
        simc_hists["epsilon"][j] = hist_bin_dict["H_epsilon_SIMC_{}".format(j)]

    print("\nBinning simc...")
    hole_contains = hgcer_cutg.IsInside if ParticleType == "kaon" else None
    if event_cache is not None:
        _fill_cached_ave_simc(event_cache, simc_hists, t_bins)
    else:
        _process_ave_simc_tree(
            TBRANCH_SIMC,
            simc_hists,
            ParticleType,
            hole_contains,
            apply_simc_cuts,
            mm_min,
            mm_max,
            t_bins,
            Misc.progressBar,
        )

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
                        
        processed_dict["t_bin{}".format(j+1)] = {
            "H_Q2_SIMC" : hist_bin_dict["H_Q2_SIMC_{}".format(j)],
            "H_W_SIMC" : hist_bin_dict["H_W_SIMC_{}".format(j)],
            "H_t_SIMC" : hist_bin_dict["H_t_SIMC_{}".format(j)],
            "H_epsilon_SIMC" : hist_bin_dict["H_epsilon_SIMC_{}".format(j)],
        }
        
    return processed_dict                    
        
def bin_simc(kinematic_types, tree_simc, t_bins, phi_setting, inpDict, iteration, event_cache=None):

    processed_dict = process_hist_simc(tree_simc, t_bins, phi_setting, inpDict, iteration, event_cache=event_cache)
    
    binned_dict = {}
    
    for kin_type in kinematic_types:
        
        # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
        binned_t_simc = []
        binned_hist_simc = []
        kin_hist_data = []
        
        # Loop through bins in t_simc and identify events in specified bins
        for j in range(len(t_bins)-1):

            H_Q2_SIMC = processed_dict["t_bin{}".format(j+1)]["H_Q2_SIMC"]
            H_W_SIMC = processed_dict["t_bin{}".format(j+1)]["H_W_SIMC"]
            H_t_SIMC = processed_dict["t_bin{}".format(j+1)]["H_t_SIMC"]
            H_epsilon_SIMC = processed_dict["t_bin{}".format(j+1)]["H_epsilon_SIMC"]
            
            # Initialize lists for tmp_binned_t_simc, tmp_binned_hist_simc, and tmp_binned_hist_dummy
            tmp_binned_t_simc = []
            tmp_binned_hist_simc = []

            tmp_hist_simc = [[],[]]
            for i in range(1, H_t_SIMC.GetNbinsX() + 1):
                tmp_hist_simc[0].append(H_t_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_t_SIMC.GetBinContent(i))                
            tmp_binned_t_simc.append(tmp_hist_simc)

            if kin_type == "t":
                kin_hist_simc.append(H_t_SIMC.Clone())
                tmp_binned_hist_simc.append(tmp_hist_simc)
            if kin_type == "Q2":
                kin_hist_simc.append(H_Q2_SIMC.Clone())
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_Q2_SIMC.GetNbinsX() + 1):        
                    tmp_hist_simc[0].append(H_Q2_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_Q2_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)
            if kin_type == "W":
                kin_hist_simc.append(H_W_SIMC.Clone())                
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_W_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc[0].append(H_W_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_W_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)        
            if kin_type == "epsilon":
                kin_hist_simc.append(H_epsilon_SIMC.Clone())                
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_epsilon_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc[0].append(H_epsilon_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_epsilon_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)

            binned_t_simc.append(tmp_binned_t_simc[0]) # Save a list of hists where each one is a t-bin
            binned_hist_simc.append(tmp_binned_hist_simc[0])

            if j+1 == len(t_bins)-1:
                binned_dict[kin_type] = {
                    "binned_t_simc" : binned_t_simc,
                    "binned_hist_simc" : binned_hist_simc,
                    "kin_hist_simc" : kin_hist_simc                    
                }
        
    return binned_dict                

def calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict, iteration):

    tree_simc = hist["InFile_SIMC"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Initialize lists for binned_t_data, binned_hist_data
    event_cache = hist.get("_yield_simc_event_cache")
    binned_dict = bin_simc(kinematic_types, tree_simc, t_bins, hist["phi_setting"], inpDict, iteration, event_cache=event_cache)
    
    kin_hist_simc = binned_dict[kin_type]["kin_hist_simc"]
    
    group_dict = {}
    
    for kin_type in kinematic_types:
        binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
        binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]
    
        ave_hist = []
        ave_err_hist = []
        binned_sub_simc = [[],[]]
        i=0 # iter
        print("-"*25)
        for simc in binned_hist_simc:
            bin_val_simc, hist_val_simc = simc
            sub_val = np.array(hist_val_simc) # No dummy subtraction for simc
            try:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(sub_val * bin_val_simc)
                total_count = np.sum(sub_val)
                average = weighted_sum / total_count
                if math.isnan(average) or math.isinf(average):
                    print("Empty binning for simc {} (t-bin={})... ".format(kin_type, i+1))
                    #sys.exit(2)
                    average = 0.0
                ave_hist.append(average)
                # Calculate the standard deviation of the simc points within the bin
                std_dev = np.std(bin_val_simc)
                # Determine the number of simc points within the bin
                n = len(bin_val_simc)
                # Calculate the standard error of the mean (SEM) for the bin
                sem = std_dev / np.sqrt(n)
                # Append the uncertainty (SEM) to the list
                ave_err_hist.append(sem)                
                #print("Weighted Sum:",weighted_sum)
                #print("Total Count:",total_count)
                #print("Average for t-bin {}:".format(i+1),average)
                binned_sub_simc[0].append(bin_val_simc)
                binned_sub_simc[1].append(sub_val)
            except ZeroDivisionError:
                print("Empty binning for simc {} (t-bin={})... ".format(kin_type, i+1))
                sys.exit(2)
            i+=1

        # Print statements to check sizes
        #print("Size of binned_t_simc:", len(binned_t_simc))
        #print("Size of binned_hist_simc:", len(binned_hist_simc))
        #print("Size of binned_sub_simc:", len(binned_sub_simc[1]))
        #print("Size of ave_hist:", len(ave_hist))
        #print("Size of t_bins:", len(t_bins))
        #print("Size of phi_bins:", len(phi_bins), "\n")

        dict_lst = []
        for j in range(len(t_bins) - 1):
            tbin_index = j
            for k in range(len(phi_bins) - 1):
                phibin_index = k
                hist_val = [binned_sub_simc[0][j], binned_sub_simc[1][j]]
                ave_val = ave_hist[j]
                ave_err_val = ave_err_hist[j]
                print("Simc average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, j+1, k+1, ave_val, ave_err_val))
                dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))

        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_ave".format(kin_type) : tup[2],
                "{}_ave_err".format(kin_type) : tup[3],
            }                    
    
        group_dict[kin_type] = groups
            
    return group_dict

def ave_per_bin_simc(histlist, inpDict, iteration=False):

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding simc averages for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        binned_dict = calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict, iteration)
        hist.pop("_yield_simc_event_cache", None)
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
        
    return {"binned_SIMC" : aveDict}

##################################################################################################################################################

def grab_ave_data(histlist, inpDict):

    OutFilename = inpDict["OutFilename"]    
    
    NumtBins = inpDict["NumtBins"]
    W = inpDict["W"]
    Q2 = inpDict["Q2"]
    EPSVAL = float(inpDict["EPSVAL"] )    
    ParticleType = inpDict["ParticleType"]

    POL = float(inpDict["POL"])
        
    if POL > 0:
        polID = 'pl'
    else:
        polID = 'mn'
    
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]

    pThetaValRight = inpDict["pThetaValRight"]
    EbeamValRight = inpDict["EbeamValRight"]
    pThetaValLeft = inpDict["pThetaValLeft"]
    EbeamValLeft = inpDict["EbeamValLeft"]
    pThetaValCenter = inpDict["pThetaValCenter"]
    EbeamValCenter = inpDict["EbeamValCenter"]
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }

    # Define thpq vector relative to middle setting
    if hist["phi_setting"] == "Right":
        runNums = np.array([int(x) for x in runNumRight.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = f"{LTANAPATH}/log/{hist['phi_setting']}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                ebeam_right = float(EbeamValRight[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_-{}.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))                
                break
            else:
                continue

    if hist["phi_setting"] == "Left":
        runNums = np.array([int(x) for x in runNumLeft.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = f"{LTANAPATH}/log/{hist['phi_setting']}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                ebeam_left = float(EbeamValLeft[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))                
                break
            else:
                continue

    if hist["phi_setting"] == "Center":
        runNums = np.array([int(x) for x in runNumCenter.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = f"{LTANAPATH}/log/{hist['phi_setting']}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                ebeam_center = float(EbeamValCenter[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_center*1000))                
                break
            else:
                continue
    
    #f_avek = '{}/averages/avek.Q{}W{}.dat'.format(ParticleType, Q2.replace("p",""), W.replace("p",""))
    
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data averages for {}...".format(hist["phi_setting"]))
        print("\nIteration, therefore grabbing data from {}...".format(f_kindata))        
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        group_dict = {}
        for kin_type in kinematic_types:
            with open(f_kindata, 'r') as f:
                lines = f.readlines()
            dict_lst = []
            for j, line in enumerate(lines): # j is t-bin
                for k in range(len(phi_bins) - 1):
                    line_lst = line.split(" ") # aveQ2, errQ2, aveW, errW, avett, errtt
                    if kin_type == "Q2":
                        ave_val = float(line_lst[0])
                        ave_err_val = float(line_lst[1])                        
                    if kin_type == "W":
                        ave_val = float(line_lst[2])
                        ave_err_val = float(line_lst[3])
                    if kin_type == "t":
                        ave_val = float(line_lst[4])
                        ave_err_val = float(line_lst[5])
                    tbin_index = j
                    phibin_index = k
                    print("Data average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, tbin_index+1, phibin_index+1, ave_val, ave_err_val))
                    dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))

            # Group the tuples by the first two elements using defaultdict
            groups = defaultdict(list)
            for tup in dict_lst:
                key = (tup[0], tup[1])
                groups[key] = {
                    "{}_ave".format(kin_type) : tup[2],
                    "{}_ave_err".format(kin_type) : tup[3],                    
                }
                
            group_dict[kin_type] = groups
        
        binned_dict = group_dict
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
                
    return {"binned_DATA" : aveDict}
