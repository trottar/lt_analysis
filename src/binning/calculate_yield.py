#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 15:20:09 trottar"
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
from copy import deepcopy
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

from binning_helpers import find_bin_index, find_2d_bin_indices
from theta_cm import calculate_theta_cm_deg, calculate_tmin

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
    is_hist,
    remove_bad_bins,
    integrate_hist_range,
    prune_hist,
    compute_staged_particle_subtraction_scales,
)
from prompt_trees import get_prompt_tree_name, get_rand_tree_name
from background_config import (
    BG_OVERSUB_WARN_FRACTION,
    BG_OVERSUB_WARN_MAX_RATIO,
    get_bg_scale_setting_key,
    resolve_particle_subtraction_fallback_mode,
    resolve_particle_subtraction_mode,
    resolve_particle_subtraction_weight_clip_bounds,
    resolve_particle_subtraction_weight_denominator_floor,
    resolve_particle_subtraction_weight_warn_max,
    resolve_particle_subtraction_windows,
    resolve_bg_stat_scale1,
    resolve_bg_stat_scale2,
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
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    compute_hist_closure_metrics,
    evaluate_particle_subtraction_component_fit_result,
    fill_simc_shape_pion_subtraction_templates,
    handle_particle_subtraction_fallback,
    iter_component_control_source_specs,
    simc_shape_pion_weight_from_value,
    summarize_particle_subtraction_component_payload,
)
from mm_background_subtraction import (
    build_mm_background_weights,
    build_mm_background_weights_with_diagnostics,
    build_mm_residual_weights,
    clone_reset_hist,
    mm_background_weight_from_value,
)

##################################################################################################################################################

def integral_with_stat_error(hist):
    total = 0.0
    variance = 0.0
    for bin_index in range(1, hist.GetNbinsX() + 1):
        total += hist.GetBinContent(bin_index)
        bin_error = hist.GetBinError(bin_index)
        variance += bin_error * bin_error
    return total, math.sqrt(max(variance, 0.0))


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def _clone_hist_for_plot(hist):
    cloned_hist = hist.Clone()
    if hasattr(cloned_hist, "SetDirectory"):
        cloned_hist.SetDirectory(0)
    return cloned_hist


def _clone_processed_entry(entry):
    cloned = {}
    for key, val in entry.items():
        if is_hist(val):
            cloned[key] = _clone_hist_for_plot(val)
        else:
            cloned[key] = val
    return cloned


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


def _clone_processed_dict(processed_dict):
    return {
        key: _clone_processed_entry(entry)
        for key, entry in processed_dict.items()
    }


def _make_bg_opt_bin_cache_key(t_bins, phi_bins, shift_mode):
    return (
        tuple(round(float(val), 8) for val in np.asarray(t_bins, dtype=float)),
        tuple(round(float(val), 8) for val in np.asarray(phi_bins, dtype=float)),
        str(shift_mode),
    )


def _warn_if_oversub_diagnostics(inpDict, diagnostics, phi_setting, t_bin_index, phi_bin_index, fit_stage):
    if not diagnostics or bool(inpDict.get("suppress_bg_opt_warnings", False)):
        return
    fraction = float(diagnostics.get("affected_lambda_fraction", 0.0) or 0.0)
    max_ratio = float(diagnostics.get("max_unclamped_ratio", 0.0) or 0.0)
    if fraction > float(BG_OVERSUB_WARN_FRACTION) or max_ratio > float(BG_OVERSUB_WARN_MAX_RATIO):
        affected_mm_range = diagnostics.get("affected_mm_range")
        print(
            "WARNING: empirical MM background over-subtraction diagnostic exceeded threshold\n"
            "  epsset = {}\n"
            "  phi_setting = {}\n"
            "  t_bin = {}\n"
            "  phi_bin = {}\n"
            "  fit_stage = {}\n"
            "  affected_lambda_fraction = {:.4f}\n"
            "  max_unclamped_ratio = {:.4f}\n"
            "  oversub_bin_count = {}\n"
            "  affected_mm_range = {}".format(
                inpDict.get("EPSSET", ""),
                phi_setting,
                int(t_bin_index) + 1,
                int(phi_bin_index) + 1,
                fit_stage,
                fraction,
                max_ratio,
                int(diagnostics.get("oversub_bin_count", 0) or 0),
                affected_mm_range,
            )
        )


def _extract_stage_window_yields(entry, mm_min, mm_max):
    stage_map = {
        "raw_prompt": "H_MM_RAW_PROMPT_DATA",
        "after_random_subtraction": "H_MM_AFTER_RANDOM_DATA",
        "after_dummy_subtraction": "H_MM_nosub_DATA",
        "after_pion_subtraction": "H_MM_pisub_DATA",
        "after_empirical_fit1": "H_MM_fit1sub_DATA",
        "after_empirical_fit2": "H_MM_DATA",
        "final_lambda_window": "H_MM_DATA",
    }
    yields = {}
    for stage_name, hist_key in stage_map.items():
        hist = entry.get(hist_key)
        if hist is None:
            yields[stage_name] = None
            continue
        try:
            yields[stage_name] = float(integrate_hist_range(hist, mm_min, mm_max))
        except Exception:
            yields[stage_name] = None
    return yields


def _fill_t_vs_tmin_hist(hist, particle_type, pol, w, q2, minus_t, weight=None):
    minus_tmin = calculate_tmin(particle_type, pol, w, q2)
    if not (math.isfinite(minus_tmin) and math.isfinite(minus_t)):
        return

    if weight is None:
        hist.Fill(minus_tmin, minus_t)
    else:
        hist.Fill(minus_tmin, minus_t, weight)


def _get_simc_true_theta_cm_deg(evt, particle_type, pol):
    """Prefer the native SIMC truth theta_cm branch; fall back to kinematic reconstruction."""
    try:
        theta_cm_true = float(evt.thetacm)
    except AttributeError:
        theta_cm_true = float("nan")

    if math.isfinite(theta_cm_true):
        return math.degrees(abs(theta_cm_true))

    try:
        minus_t_true = -evt.ti
        return calculate_theta_cm_deg(particle_type, pol, evt.Wi, evt.Q2i, minus_t_true)
    except AttributeError:
        return float("nan")


def _init_ave_event_cache():
    cache_template = {
        "adj_t": [],
        "adj_MM": [],
        "theta_cm_deg": [],
        "Q2": [],
        "W": [],
        "epsilon": [],
        "ssxptar": [],
        "ssyptar": [],
        "hsxptar": [],
        "hsyptar": [],
        "allcuts": [],
        "nommcuts": [],
        "mm_offset": [],
        "t_index": [],
        "phi_index": [],
        "phi_shift_deg": [],
    }
    return {
        key: {name: values.copy() for name, values in cache_template.items()}
        for key in ("prompt", "dummy", "rand", "dummy_rand")
    }


def _append_ave_event(
    cache_section,
    adj_t,
    adj_MM,
    theta_cm_deg,
    q2,
    w,
    epsilon,
    ssxptar,
    ssyptar,
    hsxptar,
    hsyptar,
    allcuts,
    nommcuts,
    mm_offset=0.0,
    t_index=-1,
    phi_index=-1,
    phi_shift_deg=float("nan"),
):
    if not (allcuts or nommcuts):
        return

    cache_section["adj_t"].append(adj_t)
    cache_section["adj_MM"].append(adj_MM)
    cache_section["theta_cm_deg"].append(theta_cm_deg)
    cache_section["Q2"].append(q2)
    cache_section["W"].append(w)
    cache_section["epsilon"].append(epsilon)
    cache_section["ssxptar"].append(ssxptar)
    cache_section["ssyptar"].append(ssyptar)
    cache_section["hsxptar"].append(hsxptar)
    cache_section["hsyptar"].append(hsyptar)
    cache_section["allcuts"].append(bool(allcuts))
    cache_section["nommcuts"].append(bool(nommcuts))
    cache_section["mm_offset"].append(mm_offset)
    cache_section["t_index"].append(int(t_index))
    cache_section["phi_index"].append(int(phi_index))
    cache_section["phi_shift_deg"].append(phi_shift_deg)


def _freeze_ave_event_cache(event_cache):
    frozen_cache = {}
    for cache_key, cache_section in event_cache.items():
        frozen_section = {
            "adj_t": np.asarray(cache_section["adj_t"], dtype=np.float64),
            "adj_MM": np.asarray(cache_section["adj_MM"], dtype=np.float64),
            "theta_cm_deg": np.asarray(cache_section["theta_cm_deg"], dtype=np.float64),
            "Q2": np.asarray(cache_section["Q2"], dtype=np.float64),
            "W": np.asarray(cache_section["W"], dtype=np.float64),
            "epsilon": np.asarray(cache_section["epsilon"], dtype=np.float64),
            "ssxptar": np.asarray(cache_section["ssxptar"], dtype=np.float64),
            "ssyptar": np.asarray(cache_section["ssyptar"], dtype=np.float64),
            "hsxptar": np.asarray(cache_section["hsxptar"], dtype=np.float64),
            "hsyptar": np.asarray(cache_section["hsyptar"], dtype=np.float64),
            "allcuts": np.asarray(cache_section["allcuts"], dtype=bool),
            "nommcuts": np.asarray(cache_section["nommcuts"], dtype=bool),
            "mm_offset": np.asarray(cache_section["mm_offset"], dtype=np.float64),
            "t_index": np.asarray(cache_section["t_index"], dtype=np.int32),
            "phi_index": np.asarray(cache_section["phi_index"], dtype=np.int32),
            "phi_shift_deg": np.asarray(cache_section["phi_shift_deg"], dtype=np.float64),
        }
        frozen_section["allcut_bin_index"] = _build_allcut_bin_index(frozen_section)
        frozen_section["nommcut_bin_index"] = _build_nommcut_bin_index(frozen_section)
        frozen_cache[cache_key] = frozen_section
    return frozen_cache


def _build_allcut_bin_index(cache_section):
    bin_index = {}
    for idx in np.flatnonzero(cache_section["allcuts"]):
        t_bin = int(cache_section["t_index"][idx])
        phi_bin = int(cache_section["phi_index"][idx])
        if t_bin < 0 or phi_bin < 0:
            continue
        bin_index.setdefault((t_bin, phi_bin), []).append(int(idx))
    return {
        key: np.asarray(indices, dtype=np.int32)
        for key, indices in bin_index.items()
    }


def _build_nommcut_bin_index(cache_section):
    bin_index = {}
    for idx in np.flatnonzero(cache_section["nommcuts"]):
        t_bin = int(cache_section["t_index"][idx])
        phi_bin = int(cache_section["phi_index"][idx])
        if t_bin < 0 or phi_bin < 0:
            continue
        bin_index.setdefault((t_bin, phi_bin), []).append(int(idx))
    return {
        key: np.asarray(indices, dtype=np.int32)
        for key, indices in bin_index.items()
    }


def _noop_progress_bar(*args, **kwargs):
    return None


def _init_ave_simc_event_cache():
    return {
        "minus_t": [],
        "Q2": [],
        "W": [],
        "epsilon": [],
        "ssxptar": [],
        "ssyptar": [],
        "hsxptar": [],
        "hsyptar": [],
        "iter_weight": [],
    }


def _append_ave_simc_event(cache_section, minus_t, q2, w, epsilon, ssxptar, ssyptar, hsxptar, hsyptar, iter_weight):
    cache_section["minus_t"].append(minus_t)
    cache_section["Q2"].append(q2)
    cache_section["W"].append(w)
    cache_section["epsilon"].append(epsilon)
    cache_section["ssxptar"].append(ssxptar)
    cache_section["ssyptar"].append(ssyptar)
    cache_section["hsxptar"].append(hsxptar)
    cache_section["hsyptar"].append(hsyptar)
    cache_section["iter_weight"].append(iter_weight)


def _freeze_ave_simc_event_cache(event_cache):
    return {
        "minus_t": np.asarray(event_cache["minus_t"], dtype=np.float64),
        "Q2": np.asarray(event_cache["Q2"], dtype=np.float64),
        "W": np.asarray(event_cache["W"], dtype=np.float64),
        "epsilon": np.asarray(event_cache["epsilon"], dtype=np.float64),
        "ssxptar": np.asarray(event_cache["ssxptar"], dtype=np.float64),
        "ssyptar": np.asarray(event_cache["ssyptar"], dtype=np.float64),
        "hsxptar": np.asarray(event_cache["hsxptar"], dtype=np.float64),
        "hsyptar": np.asarray(event_cache["hsyptar"], dtype=np.float64),
        "iter_weight": np.asarray(event_cache["iter_weight"], dtype=np.float64),
    }


def _iter_yield_source_specs(
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


def _fill_yield_background_templates_for_bin(
    template_hists,
    source_specs,
    mm_reference_hist,
    mm_background_weights,
    t_index,
    phi_index,
    particle_type,
    pol,
    pion_component_payload=None,
    residual_weights=None,
):
    for source_spec in source_specs:
        cache_section = source_spec.get("cache_section")
        coeff = float(source_spec.get("coefficient", 0.0) or 0.0)
        apply_pion_weight = bool(source_spec.get("apply_pion_weight"))
        if cache_section is None or coeff == 0.0:
            continue

        bin_indices = cache_section.get("allcut_bin_index", {}).get((t_index, phi_index))
        if bin_indices is None:
            bin_indices = range(len(cache_section["adj_MM"]))

        for idx in bin_indices:
            if not cache_section["allcuts"][idx]:
                continue
            if int(cache_section["t_index"][idx]) != t_index or int(cache_section["phi_index"][idx]) != phi_index:
                continue

            adj_mm = cache_section["adj_MM"][idx]
            event_weight = coeff
            if apply_pion_weight:
                event_weight *= simc_shape_pion_weight_from_value(
                    adj_mm,
                    pion_component_payload["H_pion_control_model"],
                    pion_component_payload["weights"],
                )
            event_weight *= mm_background_weight_from_value(
                adj_mm,
                mm_reference_hist,
                mm_background_weights,
                residual_weights=residual_weights,
            )
            if event_weight == 0.0:
                continue

            adj_t = cache_section["adj_t"][idx]
            q2 = cache_section["Q2"][idx]
            w = cache_section["W"][idx]
            ssxptar = cache_section["ssxptar"][idx]
            ssyptar = cache_section["ssyptar"][idx]
            hsxptar = cache_section["hsxptar"][idx]
            hsyptar = cache_section["hsyptar"][idx]
            theta_cm_vals = cache_section.get("theta_cm_deg")
            if theta_cm_vals is None:
                theta_cm_deg = calculate_theta_cm_deg(particle_type, pol, w, q2, adj_t)
            else:
                theta_cm_deg = theta_cm_vals[idx]

            template_hists["t"].Fill(adj_t, event_weight)
            template_hists["Q2"].Fill(q2, event_weight)
            template_hists["W"].Fill(w, event_weight)
            template_hists["q2_w"].Fill(q2, w, event_weight)
            _fill_t_vs_tmin_hist(
                template_hists["t_vs_tmin"],
                particle_type,
                pol,
                w,
                q2,
                adj_t,
                weight=event_weight,
            )

            if math.isfinite(theta_cm_deg):
                template_hists["theta_cm"].Fill(theta_cm_deg, event_weight)
            template_hists["ssxptar"].Fill(ssxptar, event_weight)
            template_hists["ssyptar"].Fill(ssyptar, event_weight)
            template_hists["hsxptar"].Fill(hsxptar, event_weight)
            template_hists["hsyptar"].Fill(hsyptar, event_weight)


def _subtract_yield_mm_background_for_bin(
    hist_bin_dict,
    j,
    k,
    mm_reference_hist,
    background_hist,
    data_event_cache,
    sub_event_cache,
    normfac_data,
    normfac_dummy,
    nWindows,
    scale_factor,
    particle_type,
    pol,
    pion_component_payload=None,
    residual_weights=None,
):
    mm_background_weights, diagnostics = build_mm_background_weights_with_diagnostics(
        mm_reference_hist,
        background_hist,
    )
    template_hists = {
        "t": clone_reset_hist(hist_bin_dict["H_t_DATA_{}_{}".format(j, k)], "_bg_template"),
        "Q2": clone_reset_hist(hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)], "_bg_template"),
        "W": clone_reset_hist(hist_bin_dict["H_W_DATA_{}_{}".format(j, k)], "_bg_template"),
        "q2_w": clone_reset_hist(hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)], "_bg_template"),
        "theta_cm": clone_reset_hist(hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)], "_bg_template"),
        "ssxptar": clone_reset_hist(hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)], "_bg_template"),
        "ssyptar": clone_reset_hist(hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)], "_bg_template"),
        "hsxptar": clone_reset_hist(hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)], "_bg_template"),
        "hsyptar": clone_reset_hist(hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)], "_bg_template"),
        "t_vs_tmin": clone_reset_hist(hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)], "_bg_template"),
    }

    _fill_yield_background_templates_for_bin(
        template_hists,
        _iter_yield_source_specs(
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
        j,
        k,
        particle_type,
        pol,
        pion_component_payload=pion_component_payload,
        residual_weights=residual_weights,
    )

    hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(template_hists["t"], -1)
    hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(template_hists["Q2"], -1)
    hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Add(template_hists["W"], -1)
    hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(template_hists["q2_w"], -1)
    hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(template_hists["theta_cm"], -1)
    hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(template_hists["ssxptar"], -1)
    hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(template_hists["ssyptar"], -1)
    hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(template_hists["hsxptar"], -1)
    hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(template_hists["hsyptar"], -1)
    hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(template_hists["t_vs_tmin"], -1)

    return mm_background_weights, diagnostics


def _build_yield_component_template_hists(hist_bin_dict, j, k):
    return {
        "Q2": clone_reset_hist(hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "W": clone_reset_hist(hist_bin_dict["H_W_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "q2_w": clone_reset_hist(hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "theta_cm": clone_reset_hist(hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "t": clone_reset_hist(hist_bin_dict["H_t_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "hsxptar": clone_reset_hist(hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "hsyptar": clone_reset_hist(hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "ssxptar": clone_reset_hist(hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "ssyptar": clone_reset_hist(hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "t_vs_tmin": clone_reset_hist(hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "mm": clone_reset_hist(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)], "_pi_component_template"),
        "mm_nosub": clone_reset_hist(hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)], "_pi_component_template"),
    }


def _apply_component_pion_subtraction_for_bin(
    hist_bin_dict,
    j,
    k,
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
            context="calculate_yield component pion subtraction ({}, t{}, phi{})".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
                int(k) + 1,
            ),
        )
    if sub_event_cache is None:
        return handle_particle_subtraction_fallback(
            payload,
            "missing sub_event_cache for component-weight subtraction",
            context="calculate_yield component pion subtraction ({}, t{}, phi{})".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
                int(k) + 1,
            ),
        )

    clip_min, clip_max = resolve_particle_subtraction_weight_clip_bounds(inpDict)
    weight_payload = build_simc_shape_pion_control_weights(
        component_fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inpDict),
    )
    unsupported_bins = set(weight_payload["diagnostics"].get("unsupported_bins") or [])
    pion_reference = component_fit_result.get("H_pion_control_input")
    unsupported_overlap = 0
    if pion_reference is not None and unsupported_bins:
        for bin_index in unsupported_bins:
            if float(pion_reference.GetBinContent(int(bin_index))) > 0.0:
                unsupported_overlap += 1
    if unsupported_overlap > 0:
        return handle_particle_subtraction_fallback(
            payload,
            "unsupported pion-weight bins overlap pion-control content",
            context="calculate_yield component pion subtraction ({}, t{}, phi{})".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
                int(k) + 1,
            ),
        )

    if weight_payload["diagnostics"]["pion_weight_max"] > resolve_particle_subtraction_weight_warn_max(inpDict):
        print(
            "WARNING: pion component weight exceeded threshold\n"
            "  phi_setting = {}\n"
            "  t_bin = {}\n"
            "  phi_bin = {}\n"
            "  max_weight = {:.4f}".format(
                inpDict.get("phi_setting", ""),
                int(j) + 1,
                int(k) + 1,
                float(weight_payload["diagnostics"]["pion_weight_max"]),
            )
        )

    template_hists = _build_yield_component_template_hists(hist_bin_dict, j, k)
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
        {"t_index": j, "phi_index": k},
        particle_type,
        pol,
    )

    h_mm_before = _clone_hist_for_plot(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)])
    h_mm_nosub_before = _clone_hist_for_plot(hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)])
    h_mm_after = None

    hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(template_hists["Q2"], -1.0)
    hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Add(template_hists["W"], -1.0)
    hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(template_hists["q2_w"], -1.0)
    hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(template_hists["theta_cm"], -1.0)
    hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(template_hists["t"], -1.0)
    hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(template_hists["hsxptar"], -1.0)
    hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(template_hists["hsyptar"], -1.0)
    hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(template_hists["ssxptar"], -1.0)
    hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(template_hists["ssyptar"], -1.0)
    hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(template_hists["t_vs_tmin"], -1.0)
    hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Add(template_hists["mm_nosub"], -1.0)
    hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)].Add(template_hists["mm_nosub"], -1.0)
    hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(template_hists["mm"], -1.0)
    h_mm_after = _clone_hist_for_plot(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)])
    h_mm_nosub_after = _clone_hist_for_plot(h_mm_nosub_before)
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
            "H_pion_control_unscaled": _clone_hist_for_plot(component_fit_result.get("H_pion_control_input")),
            "H_pion_subtraction_template_MM": template_hists["mm"],
            "H_pion_subtraction_template_MM_nosub": template_hists["mm_nosub"],
            "H_pion_subtraction_template_t": template_hists["t"],
            "H_MM_before_pion_subtraction": h_mm_before,
            "H_MM_after_pion_subtraction": h_mm_after,
            "H_MM_nosub_before_pion_subtraction": h_mm_nosub_before,
            "H_MM_nosub_after_pion_subtraction": h_mm_nosub_after,
            "diagnostics": {
                **deepcopy(weight_payload["diagnostics"]),
                **deepcopy(fill_stats),
                "event_template_closure": event_template_closure,
            },
        }
    )
    return payload


def _init_hist_group_matrices(names, n_t, n_phi):
    return {
        name: [[None for _ in range(n_phi)] for _ in range(n_t)]
        for name in names
    }


def _process_yield_data_tree(
    tree,
    cache_section,
    hist_group,
    shifted_mm_getter,
    shifted_t_getter,
    t_bins,
    phi_bins,
    particle_type,
    pol,
    hole_contains,
    evaluate_event,
    mm_min,
    mm_max,
    progress_bar,
    update_mm_offset=False,
):
    mm_offset_data = None
    total_entries = tree.GetEntries()
    phi_scale = 180.0 / math.pi

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
        phi_shift = evt.ph_q * phi_scale
        theta_cm_deg = calculate_theta_cm_deg(particle_type, pol, evt.W, evt.Q2, adj_t)
        t_index, phi_index = find_2d_bin_indices(adj_t, phi_shift, t_bins, phi_bins)

        if t_index is not None:
            mm_offset = 0.0
            if allcuts:
                mm_offset = adj_MM - evt.MM
            _append_ave_event(
                cache_section,
                adj_t,
                adj_MM,
                theta_cm_deg,
                evt.Q2,
                evt.W,
                evt.epsilon,
                evt.ssxptar,
                evt.ssyptar,
                evt.hsxptar,
                evt.hsyptar,
                allcuts,
                nommcuts,
                mm_offset,
                t_index=t_index if t_index is not None else -1,
                phi_index=phi_index if phi_index is not None else -1,
                phi_shift_deg=phi_shift,
            )

        if t_index is None or phi_index is None:
            continue

        if nommcuts:
            hist_group["fit1sub"][t_index][phi_index].Fill(adj_MM)
            hist_group["pisub"][t_index][phi_index].Fill(adj_MM)
            hist_group["nosub"][t_index][phi_index].Fill(adj_MM)

        if allcuts:
            hist_group["t"][t_index][phi_index].Fill(adj_t)
            hist_group["Q2"][t_index][phi_index].Fill(evt.Q2)
            hist_group["W"][t_index][phi_index].Fill(evt.W)
            hist_group["q2_w"][t_index][phi_index].Fill(evt.Q2, evt.W)
            if math.isfinite(theta_cm_deg):
                hist_group["theta_cm"][t_index][phi_index].Fill(theta_cm_deg)
            hist_group["ssxptar"][t_index][phi_index].Fill(evt.ssxptar)
            hist_group["ssyptar"][t_index][phi_index].Fill(evt.ssyptar)
            hist_group["hsxptar"][t_index][phi_index].Fill(evt.hsxptar)
            hist_group["hsyptar"][t_index][phi_index].Fill(evt.hsyptar)
            _fill_t_vs_tmin_hist(
                hist_group["t_vs_tmin"][t_index][phi_index],
                particle_type,
                pol,
                evt.W,
                evt.Q2,
                adj_t,
            )
            hist_group["mm"][t_index][phi_index].Fill(adj_MM)
            if update_mm_offset:
                mm_offset_data = adj_MM - evt.MM

    return mm_offset_data


def _process_yield_simc_tree(
    tree,
    hist_group,
    ave_simc_event_cache,
    particle_type,
    pol,
    hole_contains,
    apply_simc_cuts,
    mm_min,
    mm_max,
    t_bins,
    phi_bins,
    progress_bar,
):
    total_entries = tree.GetEntries()
    phi_scale = 180.0 / math.pi

    for i, evt in enumerate(tree):
        progress_bar(i, total_entries, bar_length=25)

        if particle_type == "kaon":
            allcuts = apply_simc_cuts(evt, mm_min, mm_max) and not hole_contains(evt.phgcer_x_det, evt.phgcer_y_det)
        else:
            allcuts = apply_simc_cuts(evt, mm_min, mm_max)

        if not allcuts:
            continue

        minus_t = -evt.t
        phi_shift = evt.phipq * phi_scale
        theta_cm_deg = calculate_theta_cm_deg(particle_type, pol, evt.W, evt.Q2, minus_t)
        theta_cm_true_deg = _get_simc_true_theta_cm_deg(evt, particle_type, pol)
        t_index, phi_index = find_2d_bin_indices(minus_t, phi_shift, t_bins, phi_bins)

        if t_index is not None:
            _append_ave_simc_event(
                ave_simc_event_cache,
                minus_t,
                evt.Q2,
                evt.W,
                evt.epsilon,
                evt.ssxptar,
                evt.ssyptar,
                evt.hsxptar,
                evt.hsyptar,
                evt.iter_weight,
            )

        if t_index is None or phi_index is None:
            continue

        hist_group["t"][t_index][phi_index].Fill(minus_t, evt.iter_weight)
        hist_group["Q2"][t_index][phi_index].Fill(evt.Q2, evt.iter_weight)
        hist_group["W"][t_index][phi_index].Fill(evt.W, evt.iter_weight)
        hist_group["q2_w"][t_index][phi_index].Fill(evt.Q2, evt.W, evt.iter_weight)
        if math.isfinite(theta_cm_deg):
            hist_group["theta_cm"][t_index][phi_index].Fill(theta_cm_deg, evt.iter_weight)
        if math.isfinite(theta_cm_true_deg):
            hist_group["theta_cm_true"][t_index][phi_index].Fill(theta_cm_true_deg, evt.iter_weight)
        hist_group["ssxptar"][t_index][phi_index].Fill(evt.ssxptar, evt.iter_weight)
        hist_group["ssyptar"][t_index][phi_index].Fill(evt.ssyptar, evt.iter_weight)
        hist_group["hsxptar"][t_index][phi_index].Fill(evt.hsxptar, evt.iter_weight)
        hist_group["hsyptar"][t_index][phi_index].Fill(evt.hsyptar, evt.iter_weight)
        _fill_t_vs_tmin_hist(
            hist_group["t_vs_tmin"][t_index][phi_index],
            particle_type,
            pol,
            evt.W,
            evt.Q2,
            minus_t,
            weight=evt.iter_weight,
        )
        hist_group["mm"][t_index][phi_index].Fill(evt.missmass, evt.iter_weight)
        hist_group["mm_unweighted"][t_index][phi_index].Fill(evt.missmass)


def process_hist_data(
    tree_data,
    tree_dummy,
    normfac_data,
    normfac_dummy,
    t_bins,
    phi_bins,
    nWindows,
    phi_setting,
    inpDict,
    particle_subtraction_scale_factor=None,
    kaon_signal_shape_payload=None,
    kaon_sigma0_shape_payload=None,
):
    emit_plots = inpDict.get("yield_emit_plots", True)
    suppress_scale_warnings = bool(inpDict.get("suppress_bg_opt_warnings", False))
    progress_bar = Misc.progressBar if bool(inpDict.get("yield_show_progress", True)) else _noop_progress_bar

    processed_dict = {}
    support_hist_dict = _init_hist_group_matrices(
        ("Q2", "W", "q2_w", "theta_cm", "theta_cm_true", "mm", "t", "t_vs_tmin", "xptar", "yptar", "ssxptar", "ssyptar", "hsxptar", "hsyptar"),
        len(t_bins) - 1,
        len(phi_bins) - 1,
    )
    ave_event_cache = _init_ave_event_cache()
    sub_event_cache = None
    
    OutFilename = inpDict["OutFilename"] 

    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"]
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]    
    POL = inpDict["POL"]
    theta_cm_min = 0.0
    theta_cm_max = 180.0
    tmin_plot_min = 0.0
    tmin_plot_max = max(inpDict["tmax"], 0.0)
    
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
    n_phi = len(phi_bins) - 1
    particle_subtraction_mode = resolve_particle_subtraction_mode(inpDict)
    component_shape_payload = None
    component_fit_results = [[None for _ in range(n_phi)] for _ in range(n_t)]
    component_subtraction_payloads = [[None for _ in range(n_phi)] for _ in range(n_t)]
    dummy_norm_hist_dict = {}
    yield_hist_names = ("Q2", "W", "q2_w", "theta_cm", "mm", "fit1sub", "pisub", "nosub", "t", "t_vs_tmin", "ssxptar", "ssyptar", "hsxptar", "hsyptar")
    data_hists = _init_hist_group_matrices(yield_hist_names, n_t, n_phi)
    dummy_hists = _init_hist_group_matrices(yield_hist_names, n_t, n_phi)
    rand_hists = _init_hist_group_matrices(yield_hist_names, n_t, n_phi)
    dummy_rand_hists = _init_hist_group_matrices(yield_hist_names, n_t, n_phi)

    # Pion subtraction by scaling pion background to peak size
    if ParticleType == "kaon":
        from particle_subtraction import particle_subtraction_yield
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
                context="calculate_yield_scope_fits",
            )

    # Fit background and subtract
    from background_fit import bg_fit
        
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)] = TH1D("H_Q2_DATA_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            hist_bin_dict["H_W_DATA_{}_{}".format(j, k)] = TH1D("H_W_DATA_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_DATA_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)] = TH1D("H_theta_cm_DATA_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_DATA_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_fit1sub_DATA_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_pisub_DATA_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DATA_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)]       = TH1D("H_t_DATA_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)] = TH1D("H_hsxptar_DATA_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)] = TH1D("H_hsyptar_DATA_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)] = TH1D("H_ssxptar_DATA_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)] = TH1D("H_ssyptar_DATA_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
            hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_DATA_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_Q2_RAND_{}_{}".format(j, k)] = TH1D("H_Q2_RAND_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            hist_bin_dict["H_W_RAND_{}_{}".format(j, k)] = TH1D("H_W_RAND_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_RAND_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_theta_cm_RAND_{}_{}".format(j, k)] = TH1D("H_theta_cm_RAND_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_RAND_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_fit1sub_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_fit1sub_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_pisub_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_pisub_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)]       = TH1D("H_t_RAND_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_hsxptar_RAND_{}_{}".format(j, k)] = TH1D("H_hsxptar_RAND_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_hsyptar_RAND_{}_{}".format(j, k)] = TH1D("H_hsyptar_RAND_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssxptar_RAND_{}_{}".format(j, k)] = TH1D("H_ssxptar_RAND_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssyptar_RAND_{}_{}".format(j, k)] = TH1D("H_ssyptar_RAND_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
            hist_bin_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_RAND_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_Q2_DUMMY_{}_{}".format(j, k)] = TH1D("H_Q2_DUMMY_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            hist_bin_dict["H_W_DUMMY_{}_{}".format(j, k)] = TH1D("H_W_DUMMY_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_DUMMY_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)] = TH1D("H_theta_cm_DUMMY_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_fit1sub_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_fit1sub_DUMMY_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_pisub_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_pisub_DUMMY_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DUMMY_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)] = TH1D("H_hsxptar_DUMMY_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)] = TH1D("H_hsyptar_DUMMY_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)] = TH1D("H_ssxptar_DUMMY_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)] = TH1D("H_ssyptar_DUMMY_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
            hist_bin_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_DUMMY_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_Q2_DUMMY_RAND_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            hist_bin_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_W_DUMMY_RAND_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_theta_cm_DUMMY_RAND_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_RAND_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_fit1sub_DUMMY_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_pisub_DUMMY_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k),"MM", 100, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_RAND_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_hsxptar_DUMMY_RAND_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_hsyptar_DUMMY_RAND_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_ssxptar_DUMMY_RAND_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_ssyptar_DUMMY_RAND_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
            hist_bin_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])

            data_hists["Q2"][j][k] = hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)]
            data_hists["W"][j][k] = hist_bin_dict["H_W_DATA_{}_{}".format(j, k)]
            data_hists["q2_w"][j][k] = hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)]
            data_hists["theta_cm"][j][k] = hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)]
            data_hists["mm"][j][k] = hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]
            data_hists["fit1sub"][j][k] = hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)]
            data_hists["pisub"][j][k] = hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)]
            data_hists["nosub"][j][k] = hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)]
            data_hists["t"][j][k] = hist_bin_dict["H_t_DATA_{}_{}".format(j, k)]
            data_hists["hsxptar"][j][k] = hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)]
            data_hists["hsyptar"][j][k] = hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)]
            data_hists["ssxptar"][j][k] = hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)]
            data_hists["ssyptar"][j][k] = hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)]
            data_hists["t_vs_tmin"][j][k] = hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)]

            rand_hists["Q2"][j][k] = hist_bin_dict["H_Q2_RAND_{}_{}".format(j, k)]
            rand_hists["W"][j][k] = hist_bin_dict["H_W_RAND_{}_{}".format(j, k)]
            rand_hists["q2_w"][j][k] = hist_bin_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)]
            rand_hists["theta_cm"][j][k] = hist_bin_dict["H_theta_cm_RAND_{}_{}".format(j, k)]
            rand_hists["mm"][j][k] = hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)]
            rand_hists["fit1sub"][j][k] = hist_bin_dict["H_MM_fit1sub_RAND_{}_{}".format(j, k)]
            rand_hists["pisub"][j][k] = hist_bin_dict["H_MM_pisub_RAND_{}_{}".format(j, k)]
            rand_hists["nosub"][j][k] = hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)]
            rand_hists["t"][j][k] = hist_bin_dict["H_t_RAND_{}_{}".format(j, k)]
            rand_hists["hsxptar"][j][k] = hist_bin_dict["H_hsxptar_RAND_{}_{}".format(j, k)]
            rand_hists["hsyptar"][j][k] = hist_bin_dict["H_hsyptar_RAND_{}_{}".format(j, k)]
            rand_hists["ssxptar"][j][k] = hist_bin_dict["H_ssxptar_RAND_{}_{}".format(j, k)]
            rand_hists["ssyptar"][j][k] = hist_bin_dict["H_ssyptar_RAND_{}_{}".format(j, k)]
            rand_hists["t_vs_tmin"][j][k] = hist_bin_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)]

            dummy_hists["Q2"][j][k] = hist_bin_dict["H_Q2_DUMMY_{}_{}".format(j, k)]
            dummy_hists["W"][j][k] = hist_bin_dict["H_W_DUMMY_{}_{}".format(j, k)]
            dummy_hists["q2_w"][j][k] = hist_bin_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)]
            dummy_hists["theta_cm"][j][k] = hist_bin_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)]
            dummy_hists["mm"][j][k] = hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)]
            dummy_hists["fit1sub"][j][k] = hist_bin_dict["H_MM_fit1sub_DUMMY_{}_{}".format(j, k)]
            dummy_hists["pisub"][j][k] = hist_bin_dict["H_MM_pisub_DUMMY_{}_{}".format(j, k)]
            dummy_hists["nosub"][j][k] = hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)]
            dummy_hists["t"][j][k] = hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)]
            dummy_hists["hsxptar"][j][k] = hist_bin_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)]
            dummy_hists["hsyptar"][j][k] = hist_bin_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)]
            dummy_hists["ssxptar"][j][k] = hist_bin_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)]
            dummy_hists["ssyptar"][j][k] = hist_bin_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)]
            dummy_hists["t_vs_tmin"][j][k] = hist_bin_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)]

            dummy_rand_hists["Q2"][j][k] = hist_bin_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["W"][j][k] = hist_bin_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["q2_w"][j][k] = hist_bin_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["theta_cm"][j][k] = hist_bin_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["mm"][j][k] = hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["fit1sub"][j][k] = hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["pisub"][j][k] = hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["nosub"][j][k] = hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["t"][j][k] = hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["hsxptar"][j][k] = hist_bin_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["hsyptar"][j][k] = hist_bin_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["ssxptar"][j][k] = hist_bin_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["ssyptar"][j][k] = hist_bin_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)]
            dummy_rand_hists["t_vs_tmin"][j][k] = hist_bin_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)]

            # Pion subtraction by scaling simc to peak size
            if ParticleType == "kaon":

                subDict["H_Q2_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_Q2_SUB_DATA_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
                subDict["H_W_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_W_SUB_DATA_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_Q2_vs_W_SUB_DATA_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_SUB_DATA_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_theta_cm_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_theta_cm_SUB_DATA_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
                subDict["H_t_SUB_DATA_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DATA_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_hsxptar_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_hsxptar_SUB_DATA_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
                subDict["H_hsyptar_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_hsyptar_SUB_DATA_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
                subDict["H_ssxptar_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_ssxptar_SUB_DATA_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
                subDict["H_ssyptar_SUB_DATA_{}_{}".format(j, k)] = TH1D("H_ssyptar_SUB_DATA_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
                subDict["H_t_vs_tmin_SUB_DATA_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_SUB_DATA_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DATA_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DATA_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DATA_{}_{}".format(j, k),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)

                subDict["H_Q2_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_Q2_SUB_RAND_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
                subDict["H_W_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_W_SUB_RAND_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_Q2_vs_W_SUB_RAND_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_SUB_RAND_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_theta_cm_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_theta_cm_SUB_RAND_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
                subDict["H_t_SUB_RAND_{}_{}".format(j, k)]       = TH1D("H_t_SUB_RAND_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_hsxptar_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_hsxptar_SUB_RAND_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
                subDict["H_hsyptar_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_hsyptar_SUB_RAND_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
                subDict["H_ssxptar_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_ssxptar_SUB_RAND_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
                subDict["H_ssyptar_SUB_RAND_{}_{}".format(j, k)] = TH1D("H_ssyptar_SUB_RAND_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
                subDict["H_t_vs_tmin_SUB_RAND_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_SUB_RAND_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_RAND_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_RAND_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_RAND_{}_{}".format(j, k),"MM_nosub_{}".format(SubtractedParticle), 100, 0.7, 1.5)

                subDict["H_Q2_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_Q2_SUB_DUMMY_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
                subDict["H_W_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_W_SUB_DUMMY_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_Q2_vs_W_SUB_DUMMY_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_SUB_DUMMY_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_theta_cm_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_theta_cm_SUB_DUMMY_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
                subDict["H_t_SUB_DUMMY_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DUMMY_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_hsxptar_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_hsxptar_SUB_DUMMY_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
                subDict["H_hsyptar_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_hsyptar_SUB_DUMMY_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
                subDict["H_ssxptar_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_ssxptar_SUB_DUMMY_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
                subDict["H_ssyptar_SUB_DUMMY_{}_{}".format(j, k)] = TH1D("H_ssyptar_SUB_DUMMY_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
                subDict["H_t_vs_tmin_SUB_DUMMY_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_SUB_DUMMY_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DUMMY_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DUMMY_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)

                subDict["H_Q2_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_Q2_SUB_DUMMY_RAND_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
                subDict["H_W_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_W_SUB_DUMMY_RAND_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_Q2_vs_W_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_SUB_DUMMY_RAND_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
                subDict["H_theta_cm_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_theta_cm_SUB_DUMMY_RAND_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
                subDict["H_t_SUB_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DUMMY_RAND_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_hsxptar_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_hsxptar_SUB_DUMMY_RAND_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
                subDict["H_hsyptar_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_hsyptar_SUB_DUMMY_RAND_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
                subDict["H_ssxptar_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_ssxptar_SUB_DUMMY_RAND_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
                subDict["H_ssyptar_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH1D("H_ssyptar_SUB_DUMMY_RAND_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
                subDict["H_t_vs_tmin_SUB_DUMMY_RAND_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_SUB_DUMMY_RAND_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
                
    hole_contains = hgcer_cutg.IsInside if ParticleType == "kaon" else None

    print("\nBinning data...")
    MM_offset_DATA = _process_yield_data_tree(
        TBRANCH_DATA,
        ave_event_cache["prompt"],
        data_hists,
        get_shifted_mm,
        get_shifted_t,
        t_bins,
        phi_bins,
        ParticleType,
        POL,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        progress_bar,
        update_mm_offset=True,
    )

    print("\nBinning dummy...")
    _process_yield_data_tree(
        TBRANCH_DUMMY,
        ave_event_cache["dummy"],
        dummy_hists,
        get_shifted_mm,
        get_shifted_t,
        t_bins,
        phi_bins,
        ParticleType,
        POL,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        progress_bar,
    )

    print("\nBinning rand...")
    _process_yield_data_tree(
        TBRANCH_RAND,
        ave_event_cache["rand"],
        rand_hists,
        get_shifted_mm,
        get_shifted_t,
        t_bins,
        phi_bins,
        ParticleType,
        POL,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        progress_bar,
    )

    print("\nBinning dummy_rand...")
    _process_yield_data_tree(
        TBRANCH_DUMMY_RAND,
        ave_event_cache["dummy_rand"],
        dummy_rand_hists,
        get_shifted_mm,
        get_shifted_t,
        t_bins,
        phi_bins,
        ParticleType,
        POL,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        progress_bar,
    )

    # Pion subtraction by scaling pion background to peak size
    if ParticleType == "kaon":
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        subDict["MM_offset_DATA"] = MM_offset_DATA
        particle_subtraction_yield(t_bins, phi_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg)        
        sub_event_cache = subDict.get("_sub_event_cache")
        
    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_phi = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]
    stage_snapshot_dict = {}

    # Per-(t,phi) fractional uncertainty from the background fits (background_fit1/2)
    bg_fit1_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]    
    bg_fit2_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]  
    bg_oversub_diagnostics = [
        [{"fit1": {}, "fit2": {}} for _ in range(n_phi)]
        for _ in range(n_t)
    ]

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            stage_snapshot_dict["H_MM_RAW_PROMPT_DATA_{}_{}".format(j, k)] = _clone_hist_for_plot(
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]
            )
                            
            hist_bin_dict["H_Q2_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_W_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_theta_cm_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_fit1sub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_pisub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_hsxptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_hsyptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_ssxptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_ssyptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_W_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_theta_cm_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_fit1sub_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_pisub_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)],-1)            
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsxptar_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsyptar_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssxptar_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssyptar_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)],-1)

            hist_bin_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_Q2_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_W_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_fit1sub_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_fit1sub_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_pisub_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_pisub_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)],-1)            
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)],-1)        
            hist_bin_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)],-1)

            # Data Normalization
            hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            stage_snapshot_dict["H_MM_RAW_PROMPT_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Scale(normfac_data)

            # Dummy Normalization
            hist_bin_dict["H_Q2_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_W_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_MM_fit1sub_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_MM_pisub_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)            
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)        
            hist_bin_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)

            dummy_norm_hist_dict["H_MM_DUMMY_NORM_{}_{}".format(j, k)] = hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Clone(
                "H_MM_DUMMY_NORM_{}_{}".format(j, k)
            )
            stage_snapshot_dict["H_MM_AFTER_RANDOM_DATA_{}_{}".format(j, k)] = _clone_hist_for_plot(
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]
            )
            
            # Dummy subtraction            
            hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_W_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_fit1sub_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_pisub_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)], -1) 
            hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)], -1)

            # Remove histograms with less than event_threshold entries and negative integrals
            event_threshold = 1
            prune_hist(
                hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)],
                event_threshold
            )  
            prune_hist(
                hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)],
                event_threshold
            )   
            prune_hist(
                hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)],
                event_threshold
            )                               
            prune_hist(
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                event_threshold
            )
            prune_hist(
                hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                event_threshold
            )

            # Pion subtraction by scaling pion background to peak size
            if ParticleType == "kaon":
                scale_factor = 0.0
                component_payload = None
                use_legacy_scalar_subtraction = True

                if component_shape_payload is not None:
                    scope_result = build_particle_subtraction_component_result(
                        subDict[f"H_MM_nosub_SUB_DATA_{j}_{k}"],
                        hist_bin_dict[f"H_MM_nosub_DATA_{j}_{k}"],
                        resolve_scope_component_shapes(
                            component_shape_payload,
                            analysis_scope="t_bin{}phi_bin{}".format(j + 1, k + 1),
                            t_bin_index=j,
                            phi_bin_index=k,
                        ),
                        inpDict,
                        analysis_scope="t_bin{}phi_bin{}".format(j + 1, k + 1),
                        kaon_signal_shape=resolve_scope_single_shape(
                            kaon_signal_shape_payload,
                            analysis_scope="t_bin{}phi_bin{}".format(j + 1, k + 1),
                            t_bin_index=j,
                            phi_bin_index=k,
                        ),
                        kaon_sigma0_shape=resolve_scope_single_shape(
                            kaon_sigma0_shape_payload,
                            analysis_scope="t_bin{}phi_bin{}".format(j + 1, k + 1),
                            t_bin_index=j,
                            phi_bin_index=k,
                        ),
                        mm_offset_data=MM_offset_DATA,
                        context="yield_{}_t{}_phi{}".format(phi_setting, j + 1, k + 1),
                    )
                    component_fit_results[j][k] = scope_result
                    component_payload = _apply_component_pion_subtraction_for_bin(
                        hist_bin_dict,
                        j,
                        k,
                        scope_result,
                        sub_event_cache,
                        normfac_data,
                        normfac_dummy,
                        nWindows,
                        ParticleType,
                        POL,
                        {**inpDict, "phi_setting": phi_setting},
                    )
                    component_subtraction_payloads[j][k] = component_payload
                    if component_payload.get("accepted"):
                        use_legacy_scalar_subtraction = False
                        scale_factor = float(component_payload.get("particle_subtraction_effective_scale", 0.0) or 0.0)
                    else:
                        fallback_mode = component_payload.get("fallback_mode") or "single_scale"
                        if fallback_mode == "error":
                            raise RuntimeError(
                                "calculate_yield component pion subtraction ({}, t{}, phi{}) rejected: {}".format(
                                    phi_setting,
                                    int(j) + 1,
                                    int(k) + 1,
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
                                hist_bin_dict[f"H_MM_nosub_DATA_{j}_{k}"],
                                subDict[f"H_MM_nosub_SUB_DATA_{j}_{k}"],
                                subtraction_windows,
                                context="pion subtraction (t-bin {}, phi-bin {})".format(j, k),
                                quiet=suppress_scale_warnings,
                            )
                            scale_factor = scale_components["total_scale_factor"]

                            # Check that pion background is not over subtracting within kaon MM range
                            kaon_range_check = integrate_hist_range(
                                hist_bin_dict[f"H_MM_nosub_DATA_{j}_{k}"],
                                mm_min, mm_max
                            )

                            pion_range_check = integrate_hist_range(
                                subDict[f"H_MM_nosub_SUB_DATA_{j}_{k}"],
                                mm_min, mm_max
                            )

                            if pion_range_check > kaon_range_check:
                                if not suppress_scale_warnings:
                                    print("\n\nWARNING: Pion background larger than kaon peak in t-bin {}, phi-bin {}. Setting scaling factor to zero....".format(j, k))
                                scale_factor = 0.0
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
                    subDict["H_Q2_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_W_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_Q2_vs_W_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_theta_cm_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_t_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_hsxptar_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_hsyptar_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_ssxptar_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_ssyptar_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_t_vs_tmin_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_MM_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                    subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)

                    # Apply pion subtraction
                    hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(subDict["H_Q2_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_W_DATA_{}_{}".format(j, k)].Add(subDict["H_W_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(subDict["H_Q2_vs_W_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(subDict["H_theta_cm_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(subDict["H_t_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(subDict["H_hsxptar_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(subDict["H_hsyptar_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(subDict["H_ssxptar_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(subDict["H_ssyptar_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(subDict["H_t_vs_tmin_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Add(subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)].Add(subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)],-1)
                    hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(subDict["H_MM_SUB_DATA_{}_{}".format(j, k)],-1)

                arr_scale_factor[j][k] = scale_factor

            # Fit background and subtract
            # ---- Statistic‑scale for this (t,phi) bin ----------------
            inpDict["bg_stat_scale1"] = resolve_bg_stat_scale1(inpDict, phi_setting)
            # ----------------------------------------------------------------
            residual_bg_weights1 = None

            if inpDict["bg_stat_scale1"] > 0.0:
                # Fit background and subtract
                fitDict["background_fit1_{}_{}".format(j, k)] = bg_fit(
                    phi_setting,
                    inpDict,
                    hist_bin_dict[f"H_MM_pisub_DATA_{j}_{k}"],   # wide / no-MM-cut
                    hist_bin_dict[f"H_MM_DATA_{j}_{k}"],          # cut-window axis 
                    scaling=inpDict["bg_stat_scale1"],
                    model_key=f"fixquad_{phi_setting}_{EPSSET}e",
                    fit_name="Fit 1"
                )

                mm_stage1_input = _clone_hist_for_plot(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)])
                active_component_payload = component_subtraction_payloads[j][k]
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                active_component_payload = entry.get("particle_subtraction_component_payload")
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                bg_weights1, fit1_diagnostics = _subtract_yield_mm_background_for_bin(
                    hist_bin_dict,
                    j,
                    k,
                    mm_stage1_input,
                    fitDict["background_fit1_{}_{}".format(j, k)][0],
                    ave_event_cache,
                    sub_event_cache,
                    normfac_data,
                    normfac_dummy,
                    nWindows,
                    arr_scale_factor[j][k],
                    ParticleType,
                    POL,
                    pion_component_payload=active_component_payload,
                )
                bg_oversub_diagnostics[j][k]["fit1"] = fit1_diagnostics
                _warn_if_oversub_diagnostics(
                    inpDict,
                    fit1_diagnostics,
                    phi_setting,
                    j,
                    k,
                    "Fit 1",
                )
                residual_bg_weights1 = build_mm_residual_weights(bg_weights1)
                hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)].Add(fitDict["background_fit1_{}_{}".format(j, k)][1], -1)
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(fitDict["background_fit1_{}_{}".format(j, k)][0], -1)                       

                # Estimate fractional yield uncertainty from background_fit1 (covariance-propagated)
                try:
                    # This is dN_bg_norm returned by bg_fit(): (fit_hist_inrange, fit_vis, bg_par, f_sig, N_bg_norm_err)
                    dN_bg_norm = float(fitDict[f"background_fit1_{j}_{k}"][4])

                    mm_min = float(inpDict["mm_min"])
                    mm_max = float(inpDict["mm_max"])

                    hmm = hist_bin_dict[f"H_MM_DATA_{j}_{k}"]
                    ax  = hmm.GetXaxis()
                    ib_lo = max(1, ax.FindBin(mm_min))
                    ib_hi = min(ax.GetNbins(), ax.FindBin(mm_max))

                    # normalize BG uncertainty to signal-region counts only (same window as bg_fit() uses for N_bg_norm_err)
                    N_sig_norm = float(hmm.Integral(ib_lo, ib_hi))
                    if N_sig_norm < 0.0:
                        N_sig_norm = 0.0

                    if N_sig_norm > 0.0:
                        if (not math.isfinite(dN_bg_norm)) or (dN_bg_norm < 0.0):
                            # missing/invalid covariance => flag bin with large (100%) fractional BG uncertainty
                            bg_fit1_frac_err[j][k] = 1.0
                        else:
                            bg_fit1_frac_err[j][k] = abs(dN_bg_norm) / N_sig_norm
                    else:
                        bg_fit1_frac_err[j][k] = 0.0

                except KeyError:
                    bg_fit1_frac_err[j][k] = 0.0
                except Exception:
                    # Any unexpected failure: make it obvious in the error budget (but don’t change yields here)
                    bg_fit1_frac_err[j][k] = 1.0        

                # Remove histograms with less than event_threshold entries and negative integrals
                prune_hist(
                    hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)],
                    event_threshold
                )          
                prune_hist(
                    hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)],
                    event_threshold
                )   
                prune_hist(
                    hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)],
                    event_threshold
                )                                
                prune_hist(
                    hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                    event_threshold
                )
                prune_hist(
                    hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                    event_threshold
                )

            # Fit background and subtract
            # ---- Statistic‑scale for this (t,phi) bin ----------------
            inpDict["bg_stat_scale2"] = resolve_bg_stat_scale2(inpDict, phi_setting)
            # ----------------------------------------------------------------

            if inpDict["bg_stat_scale2"] > 0.0:
                # Fit background and subtract
                fitDict["background_fit2_{}_{}".format(j, k)] = bg_fit(
                    phi_setting,
                    inpDict,
                    hist_bin_dict[f"H_MM_fit1sub_DATA_{j}_{k}"],   # wide / no-MM-cut
                    hist_bin_dict[f"H_MM_DATA_{j}_{k}"],          # cut-window axis 
                    scaling=inpDict["bg_stat_scale2"],
                    model_key=f"cheb2_{phi_setting}_{EPSSET}e",
                    fit_name="Fit 2"
                )

                mm_stage2_input = _clone_hist_for_plot(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)])
                active_component_payload = component_subtraction_payloads[j][k]
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                _, fit2_diagnostics = _subtract_yield_mm_background_for_bin(
                    hist_bin_dict,
                    j,
                    k,
                    mm_stage2_input,
                    fitDict["background_fit2_{}_{}".format(j, k)][0],
                    ave_event_cache,
                    sub_event_cache,
                    normfac_data,
                    normfac_dummy,
                    nWindows,
                    arr_scale_factor[j][k],
                    ParticleType,
                    POL,
                    pion_component_payload=active_component_payload,
                    residual_weights=residual_bg_weights1,
                )
                bg_oversub_diagnostics[j][k]["fit2"] = fit2_diagnostics
                _warn_if_oversub_diagnostics(
                    inpDict,
                    fit2_diagnostics,
                    phi_setting,
                    j,
                    k,
                    "Fit 2",
                )
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(fitDict["background_fit2_{}_{}".format(j, k)][0], -1)            

                # Estimate fractional yield uncertainty from background_fit2 (covariance-propagated)
                try:
                    dN_bg_norm = float(fitDict[f"background_fit2_{j}_{k}"][4])

                    mm_min = float(inpDict["mm_min"])
                    mm_max = float(inpDict["mm_max"])

                    hmm = hist_bin_dict[f"H_MM_DATA_{j}_{k}"]
                    ax  = hmm.GetXaxis()
                    ib_lo = max(1, ax.FindBin(mm_min))
                    ib_hi = min(ax.GetNbins(), ax.FindBin(mm_max))

                    N_sig_norm = float(hmm.Integral(ib_lo, ib_hi))
                    if N_sig_norm < 0.0:
                        N_sig_norm = 0.0

                    if N_sig_norm > 0.0:
                        if (not math.isfinite(dN_bg_norm)) or (dN_bg_norm < 0.0):
                            bg_fit2_frac_err[j][k] = 1.0
                        else:
                            bg_fit2_frac_err[j][k] = abs(dN_bg_norm) / N_sig_norm
                    else:
                        bg_fit2_frac_err[j][k] = 0.0

                except KeyError:
                    bg_fit2_frac_err[j][k] = 0.0
                except Exception:
                    bg_fit2_frac_err[j][k] = 1.0

                # Remove histograms with less than event_threshold entries and negative integrals
                prune_hist(
                    hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)],
                    event_threshold
                ) 
                prune_hist(
                    hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)],
                    event_threshold
                ) 
                prune_hist(
                    hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)],
                    event_threshold
                )                                     
                prune_hist(
                    hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                    event_threshold
                )
                prune_hist(
                    hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                    event_threshold
                )            

            support_hist_dict["Q2"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)])
            support_hist_dict["W"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_W_DATA_{}_{}".format(j, k)])
            support_hist_dict["q2_w"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)])
            support_hist_dict["theta_cm"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)])
            support_hist_dict["mm"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)])
            support_hist_dict["t"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_t_DATA_{}_{}".format(j, k)])
            support_hist_dict["xptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["yptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["ssxptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["ssyptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["hsxptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["hsyptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)])
            support_hist_dict["t_vs_tmin"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)])

    # Checks for first plots and calls +'(' to Print
    canvas_iter = 0
    total_plots = (len(t_bins)-1) * (len(phi_bins)-1) * len(list(["H_MM_DATA_{}_{}".format(j, k), "H_t_DATA_{}_{}".format(j, k), "H_MM_DUMMY_{}_{}".format(j, k), "H_t_DUMMY_{}_{}".format(j, k)]))-1 # '-1' to remove t-phi bin edges
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1): 
            component_payload = component_subtraction_payloads[j][k]
            display_mm_sub_hist = subDict["H_MM_SUB_DATA_{}_{}".format(j, k)]
            display_t_sub_hist = subDict["H_t_SUB_DATA_{}_{}".format(j, k)]
            if isinstance(component_payload, dict):
                if component_payload.get("accepted"):
                    display_mm_sub_hist = component_payload.get("H_pion_subtraction_template_MM") or display_mm_sub_hist
                    display_t_sub_hist = component_payload.get("H_pion_subtraction_template_t") or display_t_sub_hist
                elif component_payload.get("fallback_mode") in ("zero", "skip_bin"):
                    display_mm_sub_hist = clone_reset_hist(subDict["H_MM_SUB_DATA_{}_{}".format(j, k)], "_component_empty_display")
                    display_t_sub_hist = clone_reset_hist(subDict["H_t_SUB_DATA_{}_{}".format(j, k)], "_component_empty_display")
            
            processed_entry = {
                "H_MM_DATA" : hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                "H_MM_DUMMY_NORM" : dummy_norm_hist_dict["H_MM_DUMMY_NORM_{}_{}".format(j, k)],
                "H_MM_RAW_PROMPT_DATA" : stage_snapshot_dict["H_MM_RAW_PROMPT_DATA_{}_{}".format(j, k)],
                "H_MM_AFTER_RANDOM_DATA" : stage_snapshot_dict["H_MM_AFTER_RANDOM_DATA_{}_{}".format(j, k)],
                "H_t_DATA" : hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                "H_MM_SUB_DATA" : display_mm_sub_hist,
                "H_t_SUB_DATA" : display_t_sub_hist,
                "scale_factor" : arr_scale_factor[j][k],
                "particle_subtraction_component_fit" : component_fit_results[j][k],
                "particle_subtraction_component_fit_summary" : serialize_particle_subtraction_component_result(
                    component_fit_results[j][k]
                ) if component_fit_results[j][k] is not None else None,
                "particle_subtraction_component_payload" : component_subtraction_payloads[j][k],
                "particle_subtraction_component_payload_summary" : summarize_particle_subtraction_component_payload(
                    component_subtraction_payloads[j][k]
                ) if component_subtraction_payloads[j][k] is not None else None,
                "particle_subtraction_component_applied" : bool(
                    isinstance(component_subtraction_payloads[j][k], dict)
                    and component_subtraction_payloads[j][k].get("accepted")
                ),
                # Fractional background-fit error for this bin
                "bg_fit1_frac_err" : bg_fit1_frac_err[j][k],        
                "bg_fit2_frac_err" : bg_fit2_frac_err[j][k],
                "oversub_diagnostics" : bg_oversub_diagnostics[j][k],
            }
            processed_entry.update({
                "H_Q2_DATA" : hist_bin_dict["H_Q2_DATA_{}_{}".format(j, k)],
                "H_W_DATA" : hist_bin_dict["H_W_DATA_{}_{}".format(j, k)],
                "H_Q2_vs_W_DATA" : hist_bin_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)],
                "H_theta_cm_DATA" : hist_bin_dict["H_theta_cm_DATA_{}_{}".format(j, k)],
                "H_hsxptar_DATA" : hist_bin_dict["H_hsxptar_DATA_{}_{}".format(j, k)],
                "H_hsyptar_DATA" : hist_bin_dict["H_hsyptar_DATA_{}_{}".format(j, k)],
                "H_ssxptar_DATA" : hist_bin_dict["H_ssxptar_DATA_{}_{}".format(j, k)],
                "H_ssyptar_DATA" : hist_bin_dict["H_ssyptar_DATA_{}_{}".format(j, k)],
                "H_MM_fit1sub_DATA" : hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)],
                "H_MM_pisub_DATA" : hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)],
                "H_MM_nosub_DATA" : hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)],
                "H_t_vs_tmin_DATA" : hist_bin_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)],
            })
            processed_entry["stage_window_yields"] = _extract_stage_window_yields(
                processed_entry,
                float(inpDict["mm_min"]),
                float(inpDict["mm_max"]),
            )
            processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)] = processed_entry

            # Sort dictionary keys alphabetically
            processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)] = {key : processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)][key] \
                                                                  for key in sorted(processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)].keys())}

            if not emit_plots:
                continue
            
            # Include Stat box
            ROOT.gStyle.SetOptStat(1)
            plot_entry = processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)]
            for i, (key,val) in enumerate(plot_entry.items()):

                # Track the absolute first and last plots across all iterations
                is_absolute_first = (canvas_iter == 0)
                is_absolute_last = (canvas_iter == total_plots)

                #print("Processing plot: {}, Canvas iter: {}".format(key, canvas_iter))

                if is_hist(val):
                    
                    if "MM_DATA" in key:
                        # Create a new canvas for each plot
                        canvas2 = ROOT.TCanvas("canvas2_{}".format(canvas_iter), "Canvas", 800, 600)

                        mm_nosub_plot = _clone_hist_for_plot(hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)])
                        mm_nosub_plot.SetLineColor(1)
                        mm_nosub_plot.SetFillStyle(3001)  # Set fill style to dots
                        mm_nosub_plot.SetFillColor(kBlack)  # Set fill color to black
                        mm_nosub_plot.Draw()
                        if ParticleType == "kaon":
                            if isinstance(component_payload, dict) and component_payload.get("accepted"):
                                mm_sub_plot = _clone_hist_for_plot(component_payload.get("H_pion_subtraction_template_MM_nosub"))
                            elif isinstance(component_payload, dict) and component_payload.get("fallback_mode") in ("zero", "skip_bin"):
                                mm_sub_plot = clone_reset_hist(subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)], "_component_empty_plot")
                            else:
                                mm_sub_plot = _clone_hist_for_plot(subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)])
                            mm_sub_plot.SetLineColor(2)
                            mm_sub_plot.Draw("same, E1")
                        mm_nosub_plot.SetTitle(mm_nosub_plot.GetName())
                        
                        # Ensure correct PDF opening and closing
                        pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                        if is_absolute_first:
                            print("(")
                            canvas2.Print(pdf_name + '(')
                        elif is_absolute_last:
                            print(")")
                            canvas2.Print(pdf_name + ')')
                        else:
                            canvas2.Print(pdf_name)
                            
                        # Close the canvas2 to free up memory
                        canvas2.Close()

                        # Create a new canvas for each plot
                        canvas3 = ROOT.TCanvas("canvas3_{}".format(canvas_iter), "Canvas", 800, 600)

                        mm_pisub_plot = _clone_hist_for_plot(hist_bin_dict["H_MM_pisub_DATA_{}_{}".format(j, k)])
                        mm_pisub_plot.SetLineColor(1)
                        mm_pisub_plot.SetFillStyle(3001)  # Set fill style to dots
                        mm_pisub_plot.SetFillColor(kBlack)  # Set fill color to black
                        mm_pisub_plot.Draw("hist same")
                        if inpDict["bg_stat_scale1"] > 0.0:
                            fit_plot = _clone_hist_for_plot(fitDict["background_fit1_{}_{}".format(j, k)][1])
                            fit_plot.SetLineColor(3)
                            fit_plot.Draw("same")
                        mm_pisub_plot.SetTitle(mm_pisub_plot.GetName())
                        
                        # Ensure correct PDF opening and closing
                        pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                        if is_absolute_first:
                            print("(")
                            canvas3.Print(pdf_name + '(')
                        elif is_absolute_last:
                            print(")")
                            canvas3.Print(pdf_name + ')')
                        else:
                            canvas3.Print(pdf_name)
                            
                        # Close the canvas2 to free up memory
                        canvas3.Close()      

                        # Create a new canvas for each plot
                        canvas4 = ROOT.TCanvas("canvas4_{}".format(canvas_iter), "Canvas", 800, 600)

                        mm_fit1sub_plot = _clone_hist_for_plot(hist_bin_dict["H_MM_fit1sub_DATA_{}_{}".format(j, k)])
                        mm_fit1sub_plot.SetLineColor(1)
                        mm_fit1sub_plot.SetFillStyle(3001)  # Set fill style to dots
                        mm_fit1sub_plot.SetFillColor(kBlack)  # Set fill color to black
                        mm_fit1sub_plot.Draw("hist same")
                        if inpDict["bg_stat_scale2"] > 0.0:
                            fit_plot = _clone_hist_for_plot(fitDict["background_fit2_{}_{}".format(j, k)][1])
                            fit_plot.SetLineColor(3)
                            fit_plot.Draw("same")
                        mm_fit1sub_plot.SetTitle(mm_fit1sub_plot.GetName())
                        
                        # Ensure correct PDF opening and closing
                        pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                        if is_absolute_first:
                            print("(")
                            canvas4.Print(pdf_name + '(')
                        elif is_absolute_last:
                            print(")")
                            canvas4.Print(pdf_name + ')')
                        else:
                            canvas4.Print(pdf_name)
                            
                        # Close the canvas2 to free up memory
                        canvas4.Close()                                          

                        if component_fit_results[j][k] is not None:
                            print_particle_subtraction_component_fit_pages(
                                pdf_name,
                                component_fit_results[j][k],
                                title_prefix="{} t{} phi{}".format(phi_setting, j + 1, k + 1),
                                cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
                            )
                        if isinstance(component_payload, dict):
                            print_particle_subtraction_component_application_pages(
                                pdf_name,
                                component_payload,
                                title_prefix="{} t{} phi{}".format(phi_setting, j + 1, k + 1),
                                cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
                            )

                    # Create a new canvas for each plot
                    canvas = ROOT.TCanvas("canvas_{}".format(canvas_iter), "Canvas", 800, 600)
                    plot_hist = _clone_hist_for_plot(val)
                    plot_hist.Draw()
                    plot_hist.SetTitle(plot_hist.GetName())
                    
                    if "MM_DATA" in key:
                        
                        # Create a TLatex object to add text to the plot
                        text = ROOT.TLatex()
                        text.SetNDC()
                        text.SetTextSize(0.04)
                        text.SetTextAlign(22) # Centered alignment
                        text.SetTextColor(ROOT.kBlack)

                        # Add the number of mesons to the plot
                        text.DrawLatex(0.7, 0.65, "{} Yield: {:.3e}".format(ParticleType.capitalize(), plot_hist.Integral()))
                    
                    # Ensure correct PDF opening and closing
                    pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                    if is_absolute_first:
                        print("(")
                        canvas.Print(pdf_name + '(')
                    elif is_absolute_last:
                        print(")")
                        canvas.Print(pdf_name + ')')
                    else:
                        canvas.Print(pdf_name)

                    # Close the canvas to free up memory
                    canvas.Close()
                        
                    # Increment canvas iterator AFTER printing
                    canvas_iter += 1
            
    return processed_dict, support_hist_dict, _freeze_ave_event_cache(ave_event_cache), sub_event_cache


def prepare_bg_opt_data_base_cache(hist, inpDict, t_bins, phi_bins):
    cache_key = _make_bg_opt_bin_cache_key(t_bins, phi_bins, inpDict.get("shift_mode", "raw"))
    cache_store = hist.setdefault("_bg_opt_data_base_cache_store", {})
    if cache_key in cache_store:
        return cache_store[cache_key]

    base_inp = dict(inpDict)
    base_inp["yield_emit_plots"] = False
    base_inp["yield_show_progress"] = False
    base_inp["suppress_bg_opt_warnings"] = True
    base_inp["bg_stat_scale1"] = 0.0
    base_inp["bg_stat_scale1_by_setting"] = {
        get_bg_scale_setting_key(base_inp["EPSSET"], hist["phi_setting"]): 0.0
    }
    base_inp["bg_stat_scale2"] = 0.0
    base_inp["bg_stat_scale2_by_setting"] = {
        get_bg_scale_setting_key(base_inp["EPSSET"], hist["phi_setting"]): 0.0
    }

    kaon_signal_shape_payload = _get_cached_kaon_signal_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "calculate_yield_bg_opt_signal",
    )
    kaon_sigma0_shape_payload = _get_cached_kaon_sigma0_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "calculate_yield_bg_opt_sigma0",
    )
    processed_dict, _, ave_event_cache, sub_event_cache = process_hist_data(
        hist["InFile_DATA"],
        hist["InFile_DUMMY"],
        hist["normfac_data"],
        hist["normfac_dummy"],
        t_bins,
        phi_bins,
        hist["nWindows"],
        hist["phi_setting"],
        base_inp,
        kaon_signal_shape_payload=kaon_signal_shape_payload,
        kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
    )

    base_cache = {
        "processed_dict": _clone_processed_dict(processed_dict),
        "ave_event_cache": ave_event_cache,
        "sub_event_cache": sub_event_cache,
    }
    cache_store[cache_key] = base_cache
    return base_cache


def _process_hist_data_from_base_cache(data_base_cache, t_bins, phi_bins, phi_setting, normfac_data, normfac_dummy, nWindows, inpDict):
    processed_dict = _clone_processed_dict(data_base_cache["processed_dict"])
    ave_event_cache = data_base_cache["ave_event_cache"]
    sub_event_cache = data_base_cache["sub_event_cache"]
    support_hist_dict = _init_hist_group_matrices(
        ("Q2", "W", "q2_w", "theta_cm", "theta_cm_true", "mm", "t", "t_vs_tmin", "xptar", "yptar", "ssxptar", "ssyptar", "hsxptar", "hsyptar"),
        len(t_bins) - 1,
        len(phi_bins) - 1,
    )

    ParticleType = inpDict["ParticleType"]
    EPSSET = inpDict["EPSSET"]
    POL = inpDict["POL"]
    mm_min = float(inpDict["mm_min"])
    mm_max = float(inpDict["mm_max"])
    event_threshold = 1

    from background_fit import bg_fit

    bg_scale1 = resolve_bg_stat_scale1(inpDict, phi_setting)
    bg_scale2 = resolve_bg_stat_scale2(inpDict, phi_setting)
    for j in range(len(t_bins) - 1):
        for k in range(len(phi_bins) - 1):
            entry = processed_dict["t_bin{}phi_bin{}".format(j + 1, k + 1)]
            residual_bg_weights1 = None
            entry.setdefault("oversub_diagnostics", {})
            if bg_scale1 > 0.0 or bg_scale2 > 0.0:
                hist_bin_dict = {
                    "H_Q2_DATA_{}_{}".format(j, k): entry["H_Q2_DATA"],
                    "H_W_DATA_{}_{}".format(j, k): entry["H_W_DATA"],
                    "H_Q2_vs_W_DATA_{}_{}".format(j, k): entry["H_Q2_vs_W_DATA"],
                    "H_theta_cm_DATA_{}_{}".format(j, k): entry["H_theta_cm_DATA"],
                    "H_MM_DATA_{}_{}".format(j, k): entry["H_MM_DATA"],
                    "H_t_DATA_{}_{}".format(j, k): entry["H_t_DATA"],
                    "H_hsxptar_DATA_{}_{}".format(j, k): entry["H_hsxptar_DATA"],
                    "H_hsyptar_DATA_{}_{}".format(j, k): entry["H_hsyptar_DATA"],
                    "H_ssxptar_DATA_{}_{}".format(j, k): entry["H_ssxptar_DATA"],
                    "H_ssyptar_DATA_{}_{}".format(j, k): entry["H_ssyptar_DATA"],
                    "H_t_vs_tmin_DATA_{}_{}".format(j, k): entry["H_t_vs_tmin_DATA"],
                }

            if bg_scale1 > 0.0:
                fit_result1 = bg_fit(
                    phi_setting,
                    inpDict,
                    entry["H_MM_pisub_DATA"],
                    entry["H_MM_DATA"],
                    scaling=bg_scale1,
                    model_key=f"fixquad_{phi_setting}_{EPSSET}e",
                    fit_name="Fit 1",
                )

                mm_stage1_input = _clone_hist_for_plot(entry["H_MM_DATA"])
                bg_weights1, fit1_diagnostics = _subtract_yield_mm_background_for_bin(
                    hist_bin_dict,
                    j,
                    k,
                    mm_stage1_input,
                    fit_result1[0],
                    ave_event_cache,
                    sub_event_cache,
                    normfac_data,
                    normfac_dummy,
                    nWindows,
                    entry["scale_factor"],
                    ParticleType,
                    POL,
                    pion_component_payload=active_component_payload,
                )
                entry["oversub_diagnostics"]["fit1"] = fit1_diagnostics
                _warn_if_oversub_diagnostics(
                    inpDict,
                    fit1_diagnostics,
                    phi_setting,
                    j,
                    k,
                    "Fit 1",
                )
                residual_bg_weights1 = build_mm_residual_weights(bg_weights1)
                entry["H_MM_fit1sub_DATA"].Add(fit_result1[1], -1)
                entry["H_MM_DATA"].Add(fit_result1[0], -1)

                try:
                    dN_bg_norm = float(fit_result1[4])
                    hmm = entry["H_MM_DATA"]
                    ax = hmm.GetXaxis()
                    ib_lo = max(1, ax.FindBin(mm_min))
                    ib_hi = min(ax.GetNbins(), ax.FindBin(mm_max))
                    N_sig_norm = float(hmm.Integral(ib_lo, ib_hi))
                    if N_sig_norm < 0.0:
                        N_sig_norm = 0.0
                    if N_sig_norm > 0.0:
                        if (not math.isfinite(dN_bg_norm)) or (dN_bg_norm < 0.0):
                            entry["bg_fit1_frac_err"] = 1.0
                        else:
                            entry["bg_fit1_frac_err"] = abs(dN_bg_norm) / N_sig_norm
                    else:
                        entry["bg_fit1_frac_err"] = 0.0
                except Exception:
                    entry["bg_fit1_frac_err"] = 1.0
            else:
                entry["bg_fit1_frac_err"] = 0.0

            if bg_scale2 > 0.0:
                fit_result2 = bg_fit(
                    phi_setting,
                    inpDict,
                    entry["H_MM_fit1sub_DATA"],
                    entry["H_MM_DATA"],
                    scaling=bg_scale2,
                    model_key=f"cheb2_{phi_setting}_{EPSSET}e",
                    fit_name="Fit 2",
                )

                mm_stage2_input = _clone_hist_for_plot(entry["H_MM_DATA"])
                active_component_payload = entry.get("particle_subtraction_component_payload")
                if not (isinstance(active_component_payload, dict) and active_component_payload.get("accepted")):
                    active_component_payload = None
                _, fit2_diagnostics = _subtract_yield_mm_background_for_bin(
                    hist_bin_dict,
                    j,
                    k,
                    mm_stage2_input,
                    fit_result2[0],
                    ave_event_cache,
                    sub_event_cache,
                    normfac_data,
                    normfac_dummy,
                    nWindows,
                    entry["scale_factor"],
                    ParticleType,
                    POL,
                    pion_component_payload=active_component_payload,
                    residual_weights=residual_bg_weights1,
                )
                entry["oversub_diagnostics"]["fit2"] = fit2_diagnostics
                _warn_if_oversub_diagnostics(
                    inpDict,
                    fit2_diagnostics,
                    phi_setting,
                    j,
                    k,
                    "Fit 2",
                )
                entry["H_MM_DATA"].Add(fit_result2[0], -1)

                try:
                    dN_bg_norm = float(fit_result2[4])
                    hmm = entry["H_MM_DATA"]
                    ax = hmm.GetXaxis()
                    ib_lo = max(1, ax.FindBin(mm_min))
                    ib_hi = min(ax.GetNbins(), ax.FindBin(mm_max))
                    N_sig_norm = float(hmm.Integral(ib_lo, ib_hi))
                    if N_sig_norm < 0.0:
                        N_sig_norm = 0.0
                    if N_sig_norm > 0.0:
                        if (not math.isfinite(dN_bg_norm)) or (dN_bg_norm < 0.0):
                            entry["bg_fit2_frac_err"] = 1.0
                        else:
                            entry["bg_fit2_frac_err"] = abs(dN_bg_norm) / N_sig_norm
                    else:
                        entry["bg_fit2_frac_err"] = 0.0
                except Exception:
                    entry["bg_fit2_frac_err"] = 1.0
            else:
                entry["bg_fit2_frac_err"] = 0.0

            if bg_scale1 > 0.0 or bg_scale2 > 0.0:
                prune_hist(entry["H_MM_fit1sub_DATA"], event_threshold)
                prune_hist(entry["H_MM_pisub_DATA"], event_threshold)
                prune_hist(entry["H_MM_nosub_DATA"], event_threshold)
                prune_hist(entry["H_MM_DATA"], event_threshold)
                prune_hist(entry["H_t_DATA"], event_threshold)

            entry["stage_window_yields"] = _extract_stage_window_yields(entry, mm_min, mm_max)

            support_hist_dict["Q2"][j][k] = _clone_hist_for_plot(entry["H_Q2_DATA"])
            support_hist_dict["W"][j][k] = _clone_hist_for_plot(entry["H_W_DATA"])
            support_hist_dict["q2_w"][j][k] = _clone_hist_for_plot(entry["H_Q2_vs_W_DATA"])
            support_hist_dict["theta_cm"][j][k] = _clone_hist_for_plot(entry["H_theta_cm_DATA"])
            support_hist_dict["mm"][j][k] = _clone_hist_for_plot(entry["H_MM_DATA"])
            support_hist_dict["t"][j][k] = _clone_hist_for_plot(entry["H_t_DATA"])
            support_hist_dict["xptar"][j][k] = _clone_hist_for_plot(entry["H_ssxptar_DATA"])
            support_hist_dict["yptar"][j][k] = _clone_hist_for_plot(entry["H_ssyptar_DATA"])
            support_hist_dict["ssxptar"][j][k] = _clone_hist_for_plot(entry["H_ssxptar_DATA"])
            support_hist_dict["ssyptar"][j][k] = _clone_hist_for_plot(entry["H_ssyptar_DATA"])
            support_hist_dict["hsxptar"][j][k] = _clone_hist_for_plot(entry["H_hsxptar_DATA"])
            support_hist_dict["hsyptar"][j][k] = _clone_hist_for_plot(entry["H_hsyptar_DATA"])
            support_hist_dict["t_vs_tmin"][j][k] = _clone_hist_for_plot(entry["H_t_vs_tmin_DATA"])

    return processed_dict, support_hist_dict, ave_event_cache, sub_event_cache

def bin_data(
    kin_type,
    tree_data,
    tree_dummy,
    normfac_data,
    normfac_dummy,
    t_bins,
    phi_bins,
    nWindows,
    phi_setting,
    inpDict,
    data_base_cache=None,
    particle_subtraction_scale_factor=None,
    kaon_signal_shape_payload=None,
    kaon_sigma0_shape_payload=None,
):

    if data_base_cache is None:
        processed_dict, support_hist_dict, ave_event_cache, sub_event_cache = process_hist_data(
            tree_data,
            tree_dummy,
            normfac_data,
            normfac_dummy,
            t_bins,
            phi_bins,
            nWindows,
            phi_setting,
            inpDict,
            particle_subtraction_scale_factor=particle_subtraction_scale_factor,
            kaon_signal_shape_payload=kaon_signal_shape_payload,
            kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
        )
    else:
        processed_dict, support_hist_dict, ave_event_cache, sub_event_cache = _process_hist_data_from_base_cache(
            data_base_cache,
            t_bins,
            phi_bins,
            phi_setting,
            normfac_data,
            normfac_dummy,
            nWindows,
            inpDict,
        )
    
    binned_dict = {}

    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_sub
    binned_t_data = []
    binned_hist_data = []
    binned_hist_sub = []
    mm_hist_data = []
    mm_hist_sub = []
    mm_hist_dummy_norm = []

    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_phi = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]

    # Background-fit fractional errors
    arr_bg_fit1_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)] 
    arr_bg_fit2_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]     

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_DATA"]
            H_MM_DUMMY_NORM = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_DUMMY_NORM"]
            H_t_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_DATA"]

            H_MM_SUB_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_SUB_DATA"]
            H_t_SUB_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_SUB_DATA"]

            arr_scale_factor[j][k] = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["scale_factor"]

            # Background-fit fractional error for this (t,phi) bin
            if "bg_fit1_frac_err" in processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]:
                arr_bg_fit1_frac_err[j][k] = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["bg_fit1_frac_err"]
            else:
                arr_bg_fit1_frac_err[j][k] = 0.0        
            if "bg_fit2_frac_err" in processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]:
                arr_bg_fit2_frac_err[j][k] = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["bg_fit2_frac_err"]
            else:
                arr_bg_fit2_frac_err[j][k] = 0.0

            mm_hist_data.append(H_MM_DATA.Clone())
            mm_hist_dummy_norm.append(H_MM_DUMMY_NORM.Clone())
            mm_hist_sub.append(H_MM_SUB_DATA.Clone())

            # Initialize lists for tmp_binned_t_data, tmp_binned_hist_data, and tmp_binned_hist_sub
            tmp_binned_t_data = []
            tmp_binned_hist_data = []
            tmp_binned_hist_sub = []

            tmp_hist_data = [[],[]]
            for i in range(1, H_t_DATA.GetNbinsX() + 1):
                tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                
            tmp_binned_t_data.append(tmp_hist_data)

            tmp_hist_data = [[],[]]                
            for i in range(1, H_MM_DATA.GetNbinsX() + 1):        
                tmp_hist_data[0].append(H_MM_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_MM_DATA.GetBinContent(i))
            tmp_binned_hist_data.append(tmp_hist_data)

            tmp_hist_sub = [[],[]]                
            for i in range(1, H_MM_SUB_DATA.GetNbinsX() + 1):
                tmp_hist_sub[0].append(H_MM_SUB_DATA.GetBinCenter(i))
                tmp_hist_sub[1].append(H_MM_SUB_DATA.GetBinContent(i))                    
            tmp_binned_hist_sub.append(tmp_hist_sub)

            binned_t_data.append(tmp_binned_t_data[0]) # Save a list of hists where each one is a t-bin
            binned_hist_data.append(tmp_binned_hist_data[0])
            binned_hist_sub.append(tmp_binned_hist_sub[0])

    binned_dict[kin_type] = {
        "processed_dict": processed_dict,
        "binned_t_data" : binned_t_data,
        "binned_hist_data" : binned_hist_data,
        "binned_hist_sub" : binned_hist_sub,
        "mm_hist_data" : mm_hist_data,
        "mm_hist_dummy_norm" : mm_hist_dummy_norm,
        "mm_hist_sub" : mm_hist_sub,
        "scale_factor" : arr_scale_factor,                                  
    }
    if inpDict["bg_stat_scale1"] > 0.0:                 
        binned_dict[kin_type]["bg_fit1_frac_err"] = arr_bg_fit1_frac_err                  
    if inpDict["bg_stat_scale2"] > 0.0:
        binned_dict[kin_type]["bg_fit2_frac_err"] = arr_bg_fit2_frac_err 
    binned_dict[kin_type]["support_hist_dict"] = support_hist_dict

    return binned_dict, ave_event_cache, sub_event_cache

def calculate_yield_data(kin_type, hist, t_bins, phi_bins, inpDict):

    tree_data, tree_dummy = hist["InFile_DATA"], hist["InFile_DUMMY"]
    normfac_data, normfac_dummy = hist["normfac_data"], hist["normfac_dummy"]
    nWindows, phi_setting = hist["nWindows"], hist["phi_setting"]
    particle_subtraction_scale_factor = hist.get("particle_subtraction_scale_factor")
    if resolve_particle_subtraction_mode(inpDict) == "simc_shape_components":
        particle_subtraction_scale_factor = None

    # Grab the setting by setting normalized error
    data_charge_err = inpDict["data_charge_err_{}".format(hist["phi_setting"].lower())]
    dummy_charge_err = inpDict["dummy_charge_err_{}".format(hist["phi_setting"].lower())]
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]

    data_base_cache = None
    if inpDict.get("bg_opt_use_data_cache"):
        data_base_cache = hist.get("_bg_opt_data_base_cache")
    kaon_signal_shape_payload = _get_cached_kaon_signal_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "calculate_yield_signal",
    )
    kaon_sigma0_shape_payload = _get_cached_kaon_sigma0_shape_payload(
        hist,
        inpDict,
        t_bins,
        phi_bins,
        "calculate_yield_sigma0",
    )
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_dict, ave_event_cache, sub_event_cache = bin_data(
        kin_type,
        tree_data,
        tree_dummy,
        normfac_data,
        normfac_dummy,
        t_bins,
        phi_bins,
        nWindows,
        phi_setting,
        inpDict,
        data_base_cache=data_base_cache,
        particle_subtraction_scale_factor=particle_subtraction_scale_factor,
        kaon_signal_shape_payload=kaon_signal_shape_payload,
        kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
    )
    hist["_yield_data_event_cache"] = ave_event_cache
    hist["_yield_sub_event_cache"] = sub_event_cache
    hist["_yield_data_processed_dict"] = binned_dict[kin_type]["processed_dict"]
    hist["_xsect_support_data"] = binned_dict[kin_type]["support_hist_dict"]

    binned_t_data = binned_dict[kin_type]["binned_t_data"]
    binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
    binned_hist_sub = binned_dict[kin_type]["binned_hist_sub"]
    mm_hist_data = binned_dict[kin_type]["mm_hist_data"]
    mm_hist_dummy_norm = binned_dict[kin_type]["mm_hist_dummy_norm"]
    mm_hist_sub = binned_dict[kin_type]["mm_hist_sub"]

    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_phi = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]

    # Background-fit fractional error per (t,phi) bin
    arr_bg_fit1_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]  
    arr_bg_fit2_frac_err = [[0.0 for _ in range(n_phi)] for _ in range(n_t)]  

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            arr_scale_factor[j][k] = binned_dict[kin_type]["scale_factor"][j][k]
            # Protect older files in case bg_fit1_frac_err is not present
            if "bg_fit1_frac_err" in binned_dict[kin_type]:
                arr_bg_fit1_frac_err[j][k] = binned_dict[kin_type]["bg_fit1_frac_err"][j][k]
            else:
                arr_bg_fit1_frac_err[j][k] = 0.0     
            if "bg_fit2_frac_err" in binned_dict[kin_type]:
                arr_bg_fit2_frac_err[j][k] = binned_dict[kin_type]["bg_fit2_frac_err"][j][k]
            else:
                arr_bg_fit2_frac_err[j][k] = 0.0                          

    nphi = len(phi_bins) - 1
    yield_hist = []
    yield_err_hist = []
    binned_sub_data = [[],[]]
    i=0 # iter
    print("-"*25)
    for data, final_hist, dummy_hist in zip(binned_hist_data, mm_hist_data, mm_hist_dummy_norm):
        j = i // nphi
        k = i %  nphi
        # Data is dummy, background (if used) subtracted and normalized
        bin_val_data, hist_val_data = data
        arr_data = np.array(hist_val_data)
        try:
            yld, yld_stat_err = integral_with_stat_error(final_hist)
            dummy_yld, _ = integral_with_stat_error(dummy_hist)

            bg_fit1_err = arr_bg_fit1_frac_err[j][k]        
            bg_fit2_err = arr_bg_fit2_frac_err[j][k]     

            # Normalization uncertainty from the data effective charge applies to the
            # final extracted yield, while the dummy normalization uncertainty applies
            # to the normalized dummy component that was subtracted.
            yld_data_norm_err = abs(yld) * data_charge_err
            yld_dummy_norm_err = abs(dummy_yld) * dummy_charge_err

            # Convert background-fit fractional terms to absolute uncertainties.
            yld_err = np.sqrt(
                yld_stat_err**2 +
                yld_data_norm_err**2 +
                yld_dummy_norm_err**2 +
                (abs(yld) * bg_fit1_err)**2 + 
                (abs(yld) * bg_fit2_err)**2
            )
        except ZeroDivisionError:
            yld = 0.0
            yld_err = -1000.0
        if not math.isfinite(yld):
            yld = 0.0
            yld_err = -1000.0
        if not math.isfinite(yld_err):
            yld_err = -1000.0
        yield_hist.append(yld)
        yield_err_hist.append(yld_err)
        binned_sub_data[0].append(bin_val_data)
        binned_sub_data[1].append(arr_data)
        i+=1

    # Print statements to check sizes
    #print("Size of binned_t_data:", len(binned_t_data))
    #print("Size of binned_phi_data:", len(binned_phi_data))
    #print("Size of binned_hist_data:", len(binned_hist_data))
    #print("Size of binned_hist_dummy:", len(binned_hist_dummy))
    #print("Size of binned_sub_data:", len(binned_sub_data[1]))
    #print("Size of yield_hist:", len(yield_hist))
    #print("Size of t_bins:", len(t_bins)-1)
    #print("Size of phi_bins:", len(phi_bins)-1, "\n")

    i = 0
    dict_lst = []
    for j in range(len(t_bins) - 1):
        tbin_index = j
        for k in range(len(phi_bins) - 1):
            phibin_index = k
            yield_val = yield_hist[i]
            yield_err_val = yield_err_hist[i]
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
            i+=1

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}".format(kin_type) : tup[2],
            "{}_err".format(kin_type) : tup[3],
        }            
            
    return groups

def find_yield_data(histlist, inpDict):
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_data("yield", hist, t_bins, phi_bins, inpDict)

    return {"binned_DATA" : yieldDict}

##################################################################################################################################################

def process_hist_simc(tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration):
    emit_plots = inpDict.get("yield_emit_plots", True)
    progress_bar = Misc.progressBar if bool(inpDict.get("yield_show_progress", True)) else _noop_progress_bar

    processed_dict = {}
    support_hist_dict = _init_hist_group_matrices(
        ("Q2", "W", "q2_w", "theta_cm", "theta_cm_true", "mm", "t", "t_vs_tmin", "xptar", "yptar", "ssxptar", "ssyptar", "hsxptar", "hsyptar"),
        len(t_bins) - 1,
        len(phi_bins) - 1,
    )
    ave_simc_event_cache = _init_ave_simc_event_cache()
    
    OutFilename = inpDict["OutFilename"] 
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]    
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]
    POL = inpDict["POL"]
    theta_cm_min = 0.0
    theta_cm_max = 180.0
    tmin_plot_min = 0.0
    tmin_plot_max = max(inpDict["tmax"], 0.0)
    
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
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, phi_setting)
        
    ################################################################################################################################################
        
    TBRANCH_SIMC  = tree_simc.Get("h10")
    
    hist_bin_dict = {}
    n_t = len(t_bins) - 1
    n_phi = len(phi_bins) - 1
    simc_hists = _init_hist_group_matrices(("Q2", "W", "q2_w", "theta_cm", "theta_cm_true", "mm", "t", "mm_unweighted", "t_vs_tmin", "ssxptar", "ssyptar", "hsxptar", "hsyptar"), n_t, n_phi)
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_Q2_SIMC_{}_{}".format(j, k)] = TH1D("H_Q2_SIMC_{}_{}".format(j, k), "Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
            hist_bin_dict["H_W_SIMC_{}_{}".format(j, k)] = TH1D("H_W_SIMC_{}_{}".format(j, k), "W", 100, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_Q2_vs_W_SIMC_{}_{}".format(j, k)] = TH2D("H_Q2_vs_W_SIMC_{}_{}".format(j, k), "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
            hist_bin_dict["H_theta_cm_SIMC_{}_{}".format(j, k)] = TH1D("H_theta_cm_SIMC_{}_{}".format(j, k), "theta_cm", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_theta_cm_true_SIMC_{}_{}".format(j, k)] = TH1D("H_theta_cm_true_SIMC_{}_{}".format(j, k), "theta_cm true", 100, theta_cm_min, theta_cm_max)
            hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)]       = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)]       = TH1D("H_t_SIMC_{}_{}".format(j, k),"-t", 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_hsxptar_SIMC_{}_{}".format(j, k)] = TH1D("H_hsxptar_SIMC_{}_{}".format(j, k), "HMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_hsyptar_SIMC_{}_{}".format(j, k)] = TH1D("H_hsyptar_SIMC_{}_{}".format(j, k), "HMS yptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssxptar_SIMC_{}_{}".format(j, k)] = TH1D("H_ssxptar_SIMC_{}_{}".format(j, k), "SHMS xptar", 100, -0.1, 0.1)
            hist_bin_dict["H_ssyptar_SIMC_{}_{}".format(j, k)] = TH1D("H_ssyptar_SIMC_{}_{}".format(j, k), "SHMS yptar", 100, -0.04, 0.04)
            hist_bin_dict["H_t_vs_tmin_SIMC_{}_{}".format(j, k)] = TH2D("H_t_vs_tmin_SIMC_{}_{}".format(j, k), "; -t_{min}; -t", 100, tmin_plot_min, tmin_plot_max, 100, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)] = TH1D("H_MM_SIMC_unweighted_{}_{}".format(j, k),"MM", 100, inpDict["mm_min"], inpDict["mm_max"])            
            #hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)]       = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 1000, 0.0, 2.0)
            #hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)]       = TH1D("H_t_SIMC_{}_{}".format(j, k),"-t", 100, 0.0, 2.0)
            #hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)] = TH1D("H_MM_SIMC_unweighted_{}_{}".format(j, k),"MM", 100, 0.0, 2.0)
            simc_hists["Q2"][j][k] = hist_bin_dict["H_Q2_SIMC_{}_{}".format(j, k)]
            simc_hists["W"][j][k] = hist_bin_dict["H_W_SIMC_{}_{}".format(j, k)]
            simc_hists["q2_w"][j][k] = hist_bin_dict["H_Q2_vs_W_SIMC_{}_{}".format(j, k)]
            simc_hists["theta_cm"][j][k] = hist_bin_dict["H_theta_cm_SIMC_{}_{}".format(j, k)]
            simc_hists["theta_cm_true"][j][k] = hist_bin_dict["H_theta_cm_true_SIMC_{}_{}".format(j, k)]
            simc_hists["mm"][j][k] = hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)]
            simc_hists["t"][j][k] = hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)]
            simc_hists["hsxptar"][j][k] = hist_bin_dict["H_hsxptar_SIMC_{}_{}".format(j, k)]
            simc_hists["hsyptar"][j][k] = hist_bin_dict["H_hsyptar_SIMC_{}_{}".format(j, k)]
            simc_hists["ssxptar"][j][k] = hist_bin_dict["H_ssxptar_SIMC_{}_{}".format(j, k)]
            simc_hists["ssyptar"][j][k] = hist_bin_dict["H_ssyptar_SIMC_{}_{}".format(j, k)]
            simc_hists["t_vs_tmin"][j][k] = hist_bin_dict["H_t_vs_tmin_SIMC_{}_{}".format(j, k)]
            simc_hists["mm_unweighted"][j][k] = hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)]

    for hist in hist_bin_dict.values():
        hist.Sumw2()

    print("\nBinning simc...")
    hole_contains = hgcer_cutg.IsInside if ParticleType == "kaon" else None
    _process_yield_simc_tree(
        TBRANCH_SIMC,
        simc_hists,
        ave_simc_event_cache,
        ParticleType,
        POL,
        hole_contains,
        apply_simc_cuts,
        mm_min,
        mm_max,
        t_bins,
        phi_bins,
        progress_bar,
    )

    # Checks for first plots and calls +'(' to Print
    canvas_iter = 0
    total_plots = (len(t_bins)-1) * (len(phi_bins)-1) * len(list(["H_MM_SIMC", "H_t_SIMC"]))-1 # '-1' to remove t-phi bin edges and NumEvts_bin_MM_SIMC_unweighted

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            # Normalize for yields
            hist_bin_dict["H_Q2_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_W_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_Q2_vs_W_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_theta_cm_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_theta_cm_true_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_hsxptar_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_hsyptar_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_ssxptar_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_ssyptar_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_t_vs_tmin_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)

            #print("SIMC yield for t-bin {} phi-bin {}: {:.3e}".format(j+1, k+1, hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Integral()))
            
            processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)] = {
                "H_MM_SIMC" : hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)],
                "H_t_SIMC" : hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)],
                "NumEvts_bin_MM_SIMC_unweighted" : hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)].Integral(),
            }

            support_hist_dict["Q2"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_Q2_SIMC_{}_{}".format(j, k)])
            support_hist_dict["W"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_W_SIMC_{}_{}".format(j, k)])
            support_hist_dict["q2_w"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_Q2_vs_W_SIMC_{}_{}".format(j, k)])
            support_hist_dict["theta_cm"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_theta_cm_SIMC_{}_{}".format(j, k)])
            support_hist_dict["theta_cm_true"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_theta_cm_true_SIMC_{}_{}".format(j, k)])
            support_hist_dict["mm"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)])
            support_hist_dict["t"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)])
            support_hist_dict["xptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssxptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["yptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssyptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["ssxptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssxptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["ssyptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_ssyptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["hsxptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_hsxptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["hsyptar"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_hsyptar_SIMC_{}_{}".format(j, k)])
            support_hist_dict["t_vs_tmin"][j][k] = _clone_hist_for_plot(hist_bin_dict["H_t_vs_tmin_SIMC_{}_{}".format(j, k)])

            # Sort dictionary keys alphabetically
            processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)] = {key : processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)][key] \
                                                                      for key in sorted(processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)].keys())}

            if not emit_plots:
                continue

            # Include Stat box
            ROOT.gStyle.SetOptStat(1)
            plot_entry = processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)]
            for i, (key,val) in enumerate(plot_entry.items()):
                
                # Create a new canvas for each plot
                canvas = ROOT.TCanvas("canvas_{}".format(canvas_iter), "Canvas", 800, 600)

                # Track the absolute first and last plots across all iterations
                is_absolute_first = (canvas_iter == 0)
                is_absolute_last = (canvas_iter == total_plots)

                #print("Processing plot: {}, Canvas iter: {}".format(key, canvas_iter))

                if is_hist(val):
                    
                    if "MM_SIMC" in key:
                        simc_plot = _clone_hist_for_plot(hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)])
                        simc_plot.SetLineColor(1)
                        simc_plot.Draw()
                        simc_plot.SetTitle(simc_plot.GetName())

                        # Create a TLatex object to add text to the plot
                        text = ROOT.TLatex()
                        text.SetNDC()
                        text.SetTextSize(0.04)
                        text.SetTextAlign(22) # Centered alignment
                        text.SetTextColor(ROOT.kBlack)

                        # Add the number of mesons to the plot
                        text.DrawLatex(0.7, 0.65, "{} Yield: {:.3e}".format(ParticleType.capitalize(), simc_plot.Integral()))

                    else:
                        plot_hist = _clone_hist_for_plot(val)
                        plot_hist.Draw()
                        plot_hist.SetTitle(plot_hist.GetName())

                    # Ensure correct PDF opening and closing
                    pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_simc_".format(phi_setting, ParticleType))

                    if is_absolute_first:
                        print("(")
                        canvas.Print(pdf_name + '(')
                    elif is_absolute_last:
                        print(")")
                        canvas.Print(pdf_name + ')')
                    else:
                        canvas.Print(pdf_name)

                    # Increment canvas iterator AFTER printing
                    canvas_iter += 1

                    # Close the canvas to free up memory
                    canvas.Close()
                
    return processed_dict, support_hist_dict, _freeze_ave_simc_event_cache(ave_simc_event_cache)

def bin_simc(kin_type, tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration):

    processed_dict, support_hist_dict, ave_simc_event_cache = process_hist_simc(tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration)
    
    binned_dict = {}

    # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
    binned_t_simc = []
    binned_hist_simc = []
    mm_hist_simc = []
    
    binned_unweighted_NumEvts_simc = []

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_SIMC"]
            H_t_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_SIMC"]
            NumEvts_bin_MM_SIMC_unweighted = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["NumEvts_bin_MM_SIMC_unweighted"]

            mm_hist_simc.append(H_MM_SIMC.Clone())            
            
            # Initialize lists for tmp_binned_t_simc, tmp_binned_hist_simc, and tmp_binned_hist_dummy
            tmp_binned_t_simc = []
            tmp_binned_hist_simc = []

            tmp_hist_simc = [[],[]]
            for i in range(1, H_t_SIMC.GetNbinsX() + 1):
                tmp_hist_simc[0].append(H_t_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_t_SIMC.GetBinContent(i))                
            tmp_binned_t_simc.append(tmp_hist_simc)

            tmp_hist_simc = [[],[]]                
            for i in range(1, H_MM_SIMC.GetNbinsX() + 1):        
                tmp_hist_simc[0].append(H_MM_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_MM_SIMC.GetBinContent(i))                    
            tmp_binned_hist_simc.append(tmp_hist_simc)

            binned_t_simc.append(tmp_binned_t_simc[0]) # Save a list of hists where each one is a t-bin
            binned_hist_simc.append(tmp_binned_hist_simc[0])
            
            binned_unweighted_NumEvts_simc.append(NumEvts_bin_MM_SIMC_unweighted)

    binned_dict[kin_type] = {
        "binned_t_simc" : binned_t_simc,
        "binned_hist_simc" : binned_hist_simc,
        "mm_hist_simc" : mm_hist_simc,
        "binned_unweighted_NumEvts_simc" : binned_unweighted_NumEvts_simc,
        "support_hist_dict" : support_hist_dict,
    }
        
    return binned_dict, ave_simc_event_cache

def calculate_yield_simc(kin_type, hist, t_bins, phi_bins, inpDict, iteration):

    tree_simc, normfac_simc = hist["InFile_SIMC"], hist["normfac_simc"]
    phi_setting = hist["phi_setting"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Initialize lists for binned_t_data, binned_hist_data
    binned_dict, ave_simc_event_cache = bin_simc(kin_type, tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration)
    hist["_yield_simc_event_cache"] = ave_simc_event_cache
    hist["_xsect_support_simc"] = binned_dict[kin_type]["support_hist_dict"]

    binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
    binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]
    mm_hist_simc = binned_dict[kin_type]["mm_hist_simc"]
    
    yield_hist = []
    yield_err_hist = []
    binned_sub_simc = [[],[]]
    i=0 # iter
    print("-"*25)
    for simc, simc_hist in zip(binned_hist_simc, mm_hist_simc):
        bin_val_simc, hist_val_simc = simc
        arr_simc = np.array(hist_val_simc)
        try:
            yld, yld_err = integral_with_stat_error(simc_hist)
        except ZeroDivisionError:
            yld = 0.0
            yld_err = 0.0
        if yld < 0.0:
            yld = 0.0
            yld_err = 0.0
        if math.isnan(yld) or math.isnan(yld_err):
            yld = 0.0
            yld_err = 0.0
        if math.isinf(yld) or math.isinf(yld_err):
            yld = 0.0
            yld_err = 0.0            
        yield_hist.append(yld)
        yield_err_hist.append(yld_err)
        binned_sub_simc[0].append(bin_val_simc)
        binned_sub_simc[1].append(arr_simc)
        i+=1

    # Print statements to check sizes
    #print("Size of binned_t_simc:", len(binned_t_simc))
    #print("Size of binned_phi_simc:", len(binned_phi_simc))
    #print("Size of binned_hist_simc:", len(binned_hist_simc))
    #print("Size of binned_sub_simc:", len(binned_sub_simc[1]))
    #print("Size of yield_hist:", len(yield_hist))
    #print("Size of t_bins:", len(t_bins)-1)
    #print("Size of phi_bins:", len(phi_bins)-1, "\n")

    i = 0
    dict_lst = []
    for j in range(len(t_bins) - 1):
        tbin_index = j
        for k in range(len(phi_bins) - 1):
            phibin_index = k
            yield_val = yield_hist[i]
            yield_err_val = yield_err_hist[i]
            print("Simc yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
            i+=1

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}".format(kin_type) : tup[2],
            "{}_err".format(kin_type) : tup[3],
        }            
            
    return groups

def find_yield_simc(histlist, inpDict, iteration=False):
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding simc yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_simc("yield", hist, t_bins, phi_bins, inpDict, iteration)
            
    return {"binned_SIMC" : yieldDict}

##################################################################################################################################################

def grab_yield_data(histlist, phisetlist, inpDict):

    OutFilename = inpDict["OutFilename"]
    
    Ws = inpDict["W"]
    Qs = inpDict["Q2"]
    Q2 = float(Qs.replace("p","."))
    W = float(Ws.replace("p","."))
    EPSVAL = float(inpDict["EPSVAL"] )
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]    
    ParticleType = inpDict["ParticleType"]
    pThetaValRight = inpDict["pThetaValRight"]
    EbeamValRight = inpDict["EbeamValRight"]
    pThetaValLeft = inpDict["pThetaValLeft"]
    EbeamValLeft = inpDict["EbeamValLeft"]
    pThetaValCenter = inpDict["pThetaValCenter"]
    EbeamValCenter = inpDict["EbeamValCenter"]    
    POL = float(inpDict["POL"])
    
    if POL > 0:
        polID = 'pl'
    else:
        polID = 'mn'
            
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:

        phiset = hist["phi_setting"]
        
        if phiset == "Right":
            runNums = np.array([int(x) for x in runNumRight.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                    ebeam_right = float(EbeamValRight[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                    Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))

        if phiset == "Left":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                    ebeam_left = float(EbeamValLeft[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                    Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))

        if phiset == "Center":
            runNums = np.array([int(x) for x in runNumCenter.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                    ebeam_center = float(EbeamValCenter[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                      Ws.replace("p",""), float(EPSVAL)*100)
            
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("\nIteration, therefore grabbing data from {}...".format(f_yield))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        with open(f_yield, 'r') as f:
            lines = f.readlines()
        dict_lst = []
        for line in lines:
            line_lst = line.split(" ") # yield, yield_err, phibin, tbin
            yield_val = float(line_lst[0])
            yield_err_val = float(line_lst[1])
            phibin_index = int(line_lst[2])-1
            tbin_index = int(line_lst[3])-1
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(tbin_index+1, phibin_index+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
        
        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}".format("yield") : tup[2],
                "{}_err".format("yield") : tup[3],
            }

        yieldDict[hist["phi_setting"]]["yield"] = groups
        
    return {"binned_DATA" : yieldDict}
