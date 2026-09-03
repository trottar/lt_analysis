#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 14:32:45 trottar"
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
import matplotlib.pyplot as plt
from copy import deepcopy
import sys, math, os, subprocess, csv, json
import traceback
import array
from time import perf_counter
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import kBlue, kBlack, kCyan, kRed, kGreen, kMagenta, kGray, kOrange, kAzure, kViolet
from functools import reduce

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
from utility import (
    open_root_file,
    remove_bad_bins,
    create_polar_plot,
    compute_staged_particle_subtraction_scales,
)
from lambda_gate_regression import write_csv_atomically, upsert_sorted_csv
from prompt_trees import get_prompt_tree_name, get_rand_tree_name
from background_config import (
    BG_OPT_MM_PLOT_MAX,
    BG_OPT_MM_PLOT_MIN,
    BG_OPT_MM_PLOT_NBINS,
    BG_OVERSUB_WARN_FRACTION,
    BG_OVERSUB_WARN_MAX_RATIO,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    resolve_particle_subtraction_mode,
    resolve_pion_subtraction_scope,
    resolve_particle_subtraction_windows,
    resolve_bg_stat_scale1,
    resolve_bg_stat_scale2,
    resolve_particle_subtraction_weight_clip_bounds,
    resolve_particle_subtraction_weight_denominator_floor,
    resolve_particle_subtraction_weight_warn_max,
    get_proton_contamination_cleaning_config,
    get_pion_hgcer_diagnostic_config,
    get_pion_hgcer_transfer_config,
    fingerprint_hgcer_pid_contract,
    hgcer_mask_accepts,
    get_particle_subtraction_setting_key,
    get_pion_component_dynamic_alignment_config,
)
from pion_component_fits import (
    _resolve_kaon_lambda_reference_for_plot,
    build_particle_subtraction_component_result,
    print_particle_subtraction_component_application_pages,
    print_particle_subtraction_kaon_lambda_comparison_page,
    print_particle_subtraction_component_template_pages,
    print_particle_subtraction_component_fit_pages,
    resolve_scope_component_shapes,
    resolve_scope_single_shape,
    serialize_particle_subtraction_component_result,
    load_or_resolve_pion_component_alignment,
)
from pion_component_shapes import (
    load_kaon_simc_signal_shape,
    load_setting_pion_component_shapes,
)
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    assert_component_subtraction_payload_ownership,
    compute_hist_closure_metrics,
    evaluate_component_pion_application_proposal,
    evaluate_particle_subtraction_component_fit_result,
    resolve_frozen_parent_application_policy,
    fingerprint_histogram_content_error,
    handle_particle_subtraction_fallback,
    print_particle_subtraction_weight_support_warning,
    simc_shape_pion_weight_from_value,
    summarize_particle_subtraction_component_payload,
)
from pion_t_bin_parents import (
    build_setting_t_bin_pion_parents,
    render_setting_t_bin_pion_parent_pages,
)
from root_histogram_ownership import clone_root_histogram
from canonical_binning import find_canonical_bin
from proton_contamination_weights import (
    audit_timing_t_hgcer_display_targets,
    apply_kaon_proton_cleaning_to_targets,
    build_kaon_proton_cleaning_result,
    prepare_kaon_proton_cleaning_source_bundle,
    print_kaon_proton_cleaning_terminal_summary,
    print_kaon_proton_cleaning_pages,
    serialize_kaon_proton_cleaning_result,
    summarize_kaon_proton_cleaning_result,
)
from mm_background_subtraction import (
    build_mm_background_weights,
    build_mm_background_weights_with_diagnostics,
    build_mm_residual_weights,
    clone_reset_hist,
    mm_background_weight_from_value,
)
from apply_cuts import (
    get_kaon_data_coordinate,
    get_shift_mode,
    get_shifted_mm,
    get_shifted_t,
    set_shift_context,
)
from data_coordinates import raw_event_coordinates, validate_kaon_data_coordinate_contract
from pion_hgcer_diagnostics import (
    build_pion_hgcer_tdelta_diagnostic,
    pion_hgcer_display_text,
    render_pion_hgcer_tdelta_pages,
    resolve_pion_hgcer_delta_edges,
    serialize_pion_hgcer_tdelta_diagnostic,
    write_pion_hgcer_tdelta_json,
)
from pion_hgcer_event_contract import (
    _build_identity_no_proton_cleaning_application,
    build_pion_hgcer_event_contract,
    finalize_committed_host_application,
    summarize_pion_hgcer_event_contract,
)
from pion_hgcer_refinement_method_a import (
    build_pion_hgcer_method_a,
    summarize_pion_hgcer_method_a,
)
from pion_hgcer_refinement_method_b import (
    build_pion_hgcer_method_b,
    resolve_pion_hgcer_method_b_config,
    summarize_pion_hgcer_method_b,
)
from pion_hgcer_refinement_checkpoint import (
    build_pion_hgcer_refinement_checkpoint,
    pion_hgcer_refinement_checkpoint_filename,
    write_pion_hgcer_refinement_checkpoint_json,
)
from pion_hgcer_refinement_comparison import (
    build_pion_hgcer_ab_comparison,
    build_pion_hgcer_comparison_input_contract,
    build_pion_hgcer_method_a_comparison,
    build_pion_hgcer_method_b_comparison,
)
from pion_hgcer_phase_d_checkpoint import (
    build_pion_hgcer_phase_d_checkpoint,
    pion_hgcer_phase_d_checkpoint_filename,
    write_pion_hgcer_phase_d_checkpoint_json,
)
from full_background_subtraction_plots import (
    build_full_background_subtraction_d6_payload,
    build_full_background_subtraction_d7_payload,
    build_full_background_subtraction_d8_payload,
    build_full_background_subtraction_d9_payload,
    build_full_background_subtraction_d10_payload,
    close_full_background_subtraction_pdf,
    full_background_subtraction_pdf_path,
    open_full_background_subtraction_pdf,
    render_full_background_subtraction_procedure_pages,
)
from pion_hgcer_refinement_plots import (
    build_pdf_destinations,
    build_pdf_route_manifest,
    close_diagnostic_pdf,
    method_b_display_payload,
    method_b_display_source_parity,
    open_diagnostic_pdf,
    render_pion_hgcer_refinement_pages,
    render_pion_hgcer_ab_comparison_pages,
    render_proton_main_summary_pages,
    render_setting_warning_page,
)
from pion_hgcer_transfer import (
    audit_pion_hgcer_control_population,
    apply_frozen_pion_hgcer_transfer_map,
    build_pion_hgcer_zerope_transfer_map,
    render_pion_hgcer_zerope_transfer_failure_page,
    render_pion_hgcer_zerope_transfer_pages,
    serialize_pion_hgcer_zerope_transfer,
    write_pion_hgcer_zerope_transfer_json,
)

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
# Disable statistics box by default
#ROOT.gStyle.SetOptStat(0)
################################################################################################################################################

def _format_elapsed(seconds):
    if seconds < 60.0:
        return "{:.2f} s".format(seconds)
    minutes, remainder = divmod(seconds, 60.0)
    if minutes < 60.0:
        return "{:.0f} m {:.2f} s".format(minutes, remainder)
    hours, minutes = divmod(minutes, 60.0)
    return "{:.0f} h {:.0f} m {:.2f} s".format(hours, minutes, remainder)


def _print_rand_timer(label, elapsed, total_events=None):
    if total_events and total_events > 0:
        per_event_ms = (elapsed / total_events) * 1000.0
        print("[TIMER] {}: {} ({:.3f} ms/event)".format(label, _format_elapsed(elapsed), per_event_ms))
    else:
        print("[TIMER] {}: {}".format(label, _format_elapsed(elapsed)))


def _print_rand_debug(stage, **details):
    print("[DEBUG rand_sub] {}".format(stage), flush=True)
    for key, value in details.items():
        print("  {} = {}".format(key, value), flush=True)


def _safe_tree_entries(tree):
    if tree is None:
        return None
    try:
        return int(tree.GetEntries())
    except Exception:
        return "unavailable"


def _debug_tree_status(label, tree):
    tree_name = None
    if tree is not None:
        try:
            tree_name = tree.GetName()
        except Exception:
            tree_name = "unavailable"
    _print_rand_debug(
        "tree status",
        label=label,
        exists=bool(tree is not None),
        tree_name=tree_name,
        entries=_safe_tree_entries(tree),
    )


def _warn_if_oversub_diagnostics(inpDict, diagnostics, phi_setting, fit_stage):
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
            "  fit_stage = {}\n"
            "  affected_lambda_fraction = {:.4f}\n"
            "  max_unclamped_ratio = {:.4f}\n"
            "  oversub_bin_count = {}\n"
            "  affected_mm_range = {}".format(
                inpDict.get("EPSSET", ""),
                phi_setting,
                fit_stage,
                fraction,
                max_ratio,
                int(diagnostics.get("oversub_bin_count", 0) or 0),
                affected_mm_range,
            )
        )


def _draw_hgcer_signed_display(hist, title, status_payload=None):
    """Draw a display clone without suppressing negative signed content."""
    status_payload = dict(status_payload or {})
    if hist is None or not bool(status_payload.get("available", hist is not None)):
        text = TPaveText(0.14, 0.38, 0.86, 0.62, "NDC")
        text.SetBorderSize(1)
        text.SetFillStyle(0)
        text.SetTextAlign(22)
        text.SetTextSize(0.040)
        text.SetTextColor(kRed)
        text.AddText("HGCer diagnostic unavailable")
        text.Draw()
        return text
    if int(status_payload.get("nonzero_bin_count", 0) or 0) == 0:
        text = TPaveText(0.14, 0.38, 0.86, 0.62, "NDC")
        text.SetBorderSize(1)
        text.SetFillStyle(0)
        text.SetTextAlign(22)
        text.SetTextSize(0.035)
        text.SetTextColor(kRed)
        text.AddText("No nonzero final-display bins")
        text.AddText("See timing-t HGCer fill/audit diagnostics")
        text.Draw()
        return text
    display = clone_root_histogram(
        hist,
        scope="hgcer_display",
        role="signed_display",
        title=title,
        sumw2=False,
    )
    maximum = max(abs(float(display.GetMinimum())), abs(float(display.GetMaximum())), 1.0e-12)
    if int(status_payload.get("negative_bin_count", 0) or 0) > 0:
        # Symmetric limits retain random/dummy-subtracted negative structure.
        display.SetMinimum(-maximum)
        display.SetMaximum(maximum)
    else:
        display.SetMinimum(0.0)
    display.Draw("colz")
    return display


def _fill_rand_sub_allcuts(evt, adj_MM, adj_t, adj_hsdelta, fills):
    fills["hgcer_xy"](evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
    fills["hgcer_x_mm"](evt.P_hgcer_xAtCer, adj_MM, evt.P_hgcer_npeSum)
    fills["hgcer_y_mm"](evt.P_hgcer_yAtCer, adj_MM, evt.P_hgcer_npeSum)

    phi_shift = evt.ph_q

    fills["mm_ct"](adj_MM, evt.CTime_ROC1)
    fills["ct_beta"](evt.CTime_ROC1, evt.P_gtr_beta)
    fills["mm_beta"](adj_MM, evt.P_gtr_beta)
    fills["mm_h_cer"](adj_MM, evt.H_cer_npeSum)
    fills["mm_h_cal"](adj_MM, evt.H_cal_etottracknorm)
    fills["mm_p_cal"](adj_MM, evt.P_cal_etottracknorm)
    fills["mm_p_hgcer"](adj_MM, evt.P_hgcer_npeSum)
    fills["mm_p_aero"](adj_MM, evt.P_aero_npeSum)
    fills["phiq_t"](phi_shift, adj_t)
    fills["q2_w"](evt.Q2, evt.W)
    fills["q2_t"](evt.Q2, adj_t)
    fills["w_t"](evt.W, adj_t)
    fills["eps_t"](evt.epsilon, adj_t)
    fills["mm_t"](adj_MM, adj_t)

    polar_graph = fills["polar_graph"]
    if polar_graph is not None:
        polar_graph.SetPoint(polar_graph.GetN(), (phi_shift) * (180 / math.pi), adj_t)

    fills["h_ct"](evt.CTime_ROC1)
    fills["h_ssxfp"](evt.ssxfp)
    fills["h_ssyfp"](evt.ssyfp)
    fills["h_ssxpfp"](evt.ssxpfp)
    fills["h_ssypfp"](evt.ssypfp)
    fills["h_ssdelta"](evt.ssdelta)
    fills["h_ssxptar"](evt.ssxptar)
    fills["h_ssyptar"](evt.ssyptar)
    fills["h_hsxfp"](evt.hsxfp)
    fills["h_hsyfp"](evt.hsyfp)
    fills["h_hsxpfp"](evt.hsxpfp)
    fills["h_hsypfp"](evt.hsypfp)
    fills["h_hsdelta"](adj_hsdelta)
    fills["h_hsxptar"](evt.hsxptar)
    fills["h_hsyptar"](evt.hsyptar)
    fills["h_ph_q"](phi_shift)
    fills["h_th_q"](evt.th_q)
    fills["h_ph_recoil"](evt.ph_recoil)
    fills["h_th_recoil"](evt.th_recoil)
    fills["h_pmiss"](evt.pmiss)
    fills["h_emiss"](evt.emiss)
    fills["h_pmx"](evt.pmx)
    fills["h_pmy"](evt.pmy)
    fills["h_pmz"](evt.pmz)
    fills["h_q2"](evt.Q2)
    fills["h_t"](adj_t)
    fills["h_w"](evt.W)
    fills["h_epsilon"](evt.epsilon)
    fills["h_mm"](adj_MM)
    if fills["h_cal"] is not None:
        fills["h_cal"](evt.H_cal_etottracknorm)
    if fills["h_cer"] is not None:
        fills["h_cer"](evt.H_cer_npeSum)
    if fills["p_cal"] is not None:
        fills["p_cal"](evt.P_cal_etottracknorm)
    if fills["p_hgcer"] is not None:
        fills["p_hgcer"](evt.P_hgcer_npeSum)
    if fills["p_aero"] is not None:
        fills["p_aero"](evt.P_aero_npeSum)


def _fill_rand_sub_allcuts_weighted(evt, adj_MM, adj_t, adj_hsdelta, weight, hists):
    hists["hgcer_xy"].Fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, weight * evt.P_hgcer_npeSum)
    if "hgcer_xy_nohole" in hists:
        hists["hgcer_xy_nohole"].Fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, weight * evt.P_hgcer_npeSum)
    hists["hgcer_x_mm"].Fill(evt.P_hgcer_xAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
    if "hgcer_x_mm_nohole" in hists:
        hists["hgcer_x_mm_nohole"].Fill(evt.P_hgcer_xAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
    hists["hgcer_y_mm"].Fill(evt.P_hgcer_yAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
    if "hgcer_y_mm_nohole" in hists:
        hists["hgcer_y_mm_nohole"].Fill(evt.P_hgcer_yAtCer, adj_MM, weight * evt.P_hgcer_npeSum)

    phi_shift = evt.ph_q

    hists["mm_ct"].Fill(adj_MM, evt.CTime_ROC1, weight)
    hists["ct_beta"].Fill(evt.CTime_ROC1, evt.P_gtr_beta, weight)
    hists["mm_beta"].Fill(adj_MM, evt.P_gtr_beta, weight)
    hists["mm_h_cer"].Fill(adj_MM, evt.H_cer_npeSum, weight)
    hists["mm_h_cal"].Fill(adj_MM, evt.H_cal_etottracknorm, weight)
    hists["mm_p_cal"].Fill(adj_MM, evt.P_cal_etottracknorm, weight)
    hists["mm_p_hgcer"].Fill(adj_MM, evt.P_hgcer_npeSum, weight)
    hists["mm_p_aero"].Fill(adj_MM, evt.P_aero_npeSum, weight)
    hists["phiq_t"].Fill(phi_shift, adj_t, weight)
    hists["q2_w"].Fill(evt.Q2, evt.W, weight)
    hists["q2_t"].Fill(evt.Q2, adj_t, weight)
    hists["w_t"].Fill(evt.W, adj_t, weight)
    hists["eps_t"].Fill(evt.epsilon, adj_t, weight)
    hists["mm_t"].Fill(adj_MM, adj_t, weight)

    hists["h_ct"].Fill(evt.CTime_ROC1, weight)
    hists["h_ssxfp"].Fill(evt.ssxfp, weight)
    hists["h_ssyfp"].Fill(evt.ssyfp, weight)
    hists["h_ssxpfp"].Fill(evt.ssxpfp, weight)
    hists["h_ssypfp"].Fill(evt.ssypfp, weight)
    hists["h_ssdelta"].Fill(evt.ssdelta, weight)
    hists["h_ssxptar"].Fill(evt.ssxptar, weight)
    hists["h_ssyptar"].Fill(evt.ssyptar, weight)
    hists["h_hsxfp"].Fill(evt.hsxfp, weight)
    hists["h_hsyfp"].Fill(evt.hsyfp, weight)
    hists["h_hsxpfp"].Fill(evt.hsxpfp, weight)
    hists["h_hsypfp"].Fill(evt.hsypfp, weight)
    hists["h_hsdelta"].Fill(adj_hsdelta, weight)
    hists["h_hsxptar"].Fill(evt.hsxptar, weight)
    hists["h_hsyptar"].Fill(evt.hsyptar, weight)
    hists["h_ph_q"].Fill(phi_shift, weight)
    hists["h_th_q"].Fill(evt.th_q, weight)
    hists["h_ph_recoil"].Fill(evt.ph_recoil, weight)
    hists["h_th_recoil"].Fill(evt.th_recoil, weight)
    hists["h_pmiss"].Fill(evt.pmiss, weight)
    hists["h_emiss"].Fill(evt.emiss, weight)
    hists["h_pmx"].Fill(evt.pmx, weight)
    hists["h_pmy"].Fill(evt.pmy, weight)
    hists["h_pmz"].Fill(evt.pmz, weight)
    hists["h_q2"].Fill(evt.Q2, weight)
    hists["h_t"].Fill(adj_t, weight)
    hists["h_w"].Fill(evt.W, weight)
    hists["h_epsilon"].Fill(evt.epsilon, weight)
    if "h_cal" in hists:
        hists["h_cal"].Fill(evt.H_cal_etottracknorm, weight)
    if "h_cer" in hists:
        hists["h_cer"].Fill(evt.H_cer_npeSum, weight)
    if "p_cal" in hists:
        hists["p_cal"].Fill(evt.P_cal_etottracknorm, weight)
    if "p_hgcer" in hists:
        hists["p_hgcer"].Fill(evt.P_hgcer_npeSum, weight)
    if "p_aero" in hists:
        hists["p_aero"].Fill(evt.P_aero_npeSum, weight)


def _create_rand_sub_bg_templates(target_hists):
    return {
        key: clone_reset_hist(hist, "_bg_template")
        for key, hist in target_hists.items()
        if hist is not None
    }


def _hist_integral(hist):
    if hist is None:
        return 0.0
    try:
        return float(hist.Integral())
    except Exception:
        return 0.0


def _clone_hist_detached(hist, name=None):
    return clone_root_histogram(
        hist,
        scope="rand_sub",
        role="working_histogram",
        name=name,
        optional=True,
        sumw2=False,
    )


def _clone_component_targets_for_setting_wide_diagnostic(data_targets):
    """Clone the committed production target set for read-only diagnostics.

    ``active_component_targets`` has already passed no-RF proton cleaning and
    any configured post-proton RF restoration when this helper is called.  The
    common application builder may therefore perform its usual in-place
    subtraction on these detached clones without changing production histograms.
    """
    if not isinstance(data_targets, dict):
        raise RuntimeError("setting-wide diagnostic requires active component targets")
    diagnostic_targets = {}
    for key, histogram in data_targets.items():
        diagnostic_targets[key] = clone_root_histogram(
            histogram,
            scope="setting_wide_diagnostic",
            role="application_target_{}".format(key),
            name="{}_setting_wide_diagnostic".format(key),
            optional=True,
            sumw2=False,
        )
    return diagnostic_targets


def _post_proton_cleaning_input_metadata(proton_cleaning_application):
    """Describe the committed noRF host consumed by pion subtraction."""
    if str((proton_cleaning_application or {}).get("host_state") or "") == (
        "identity_no_proton_cleaning"
    ):
        return {
            "input_selection": "no_rf_identity_no_proton_cleaning",
            "source_target_state": "post_proton_noRF",
        }
    return {
        "input_selection": "no_rf_proton_cleaning_without_rf_restoration",
        "source_target_state": "post_proton_noRF",
    }


def _open_subtracted_particle_tree_bundle(outpath, phi_setting, subtracted_particle, data_filename, dummy_filename, epsset):
    sub_data_path = f"{outpath}/{phi_setting}_{subtracted_particle}_{data_filename}.root"
    sub_dummy_path = f"{outpath}/{phi_setting}_{subtracted_particle}_{dummy_filename}.root"
    _print_rand_debug(
        "opening subtraction ROOT files",
        phi_setting=phi_setting,
        epsset=epsset,
        subtracted_particle=subtracted_particle,
        data_path=sub_data_path,
        data_exists=os.path.exists(sub_data_path),
        dummy_path=sub_dummy_path,
        dummy_exists=os.path.exists(sub_dummy_path),
    )
    sub_root_data = open_root_file(sub_data_path)
    sub_root_dummy = open_root_file(sub_dummy_path)
    sub_prompt_tree_name = get_prompt_tree_name(
        subtracted_particle, epsset, rf_state="noRF"
    )
    sub_rand_tree_name = get_rand_tree_name(
        subtracted_particle, epsset, rf_state="noRF"
    )
    bundle = {
        "sub_root_data": sub_root_data,
        "sub_root_dummy": sub_root_dummy,
        "prompt_tree_name": sub_prompt_tree_name,
        "rand_tree_name": sub_rand_tree_name,
        "prompt_tree": sub_root_data.Get(sub_prompt_tree_name),
        "rand_tree": sub_root_data.Get(sub_rand_tree_name),
        "dummy_prompt_tree": sub_root_dummy.Get(sub_prompt_tree_name),
        "dummy_rand_tree": sub_root_dummy.Get(sub_rand_tree_name),
    }
    _print_rand_debug(
        "resolved subtraction tree names",
        prompt_tree_name=sub_prompt_tree_name,
        rand_tree_name=sub_rand_tree_name,
    )
    _debug_tree_status("sub_data_prompt", bundle["prompt_tree"])
    _debug_tree_status("sub_data_rand", bundle["rand_tree"])
    _debug_tree_status("sub_dummy_prompt", bundle["dummy_prompt_tree"])
    _debug_tree_status("sub_dummy_rand", bundle["dummy_rand_tree"])
    return bundle


def _open_kaon_proton_cleaning_tree_bundle(
    infile_data,
    infile_dummy,
    particle_type,
    pol,
    epsset,
    phi_setting,
    norm_factor_data,
    norm_factor_dummy,
    n_windows,
):
    no_rf_prompt_name = get_prompt_tree_name(particle_type, epsset, rf_state="noRF")
    no_rf_rand_name = get_rand_tree_name(particle_type, epsset, rf_state="noRF")
    rf_prompt_name = get_prompt_tree_name(particle_type, epsset, rf_state="RF")
    rf_rand_name = get_rand_tree_name(particle_type, epsset, rf_state="RF")
    bundle = {
        "sources": {
            "prompt": {
                "tree_name": no_rf_prompt_name,
                "tree": infile_data.Get(no_rf_prompt_name) if infile_data is not None else None,
                "coefficient": float(norm_factor_data),
            },
            "rand": {
                "tree_name": no_rf_rand_name,
                "tree": infile_data.Get(no_rf_rand_name) if infile_data is not None else None,
                "coefficient": -float(norm_factor_data) / float(n_windows),
            },
            "dummy_prompt": {
                "tree_name": no_rf_prompt_name,
                "tree": infile_dummy.Get(no_rf_prompt_name) if infile_dummy is not None else None,
                "coefficient": -float(norm_factor_dummy),
            },
            "dummy_rand": {
                "tree_name": no_rf_rand_name,
                "tree": infile_dummy.Get(no_rf_rand_name) if infile_dummy is not None else None,
                "coefficient": float(norm_factor_dummy) / float(n_windows),
            },
        },
        "rf_sources": {
            "prompt": {
                "tree_name": rf_prompt_name,
                "tree": infile_data.Get(rf_prompt_name) if infile_data is not None else None,
            },
            "rand": {
                "tree_name": rf_rand_name,
                "tree": infile_data.Get(rf_rand_name) if infile_data is not None else None,
            },
            "dummy_prompt": {
                "tree_name": rf_prompt_name,
                "tree": infile_dummy.Get(rf_prompt_name) if infile_dummy is not None else None,
            },
            "dummy_rand": {
                "tree_name": rf_rand_name,
                "tree": infile_dummy.Get(rf_rand_name) if infile_dummy is not None else None,
            },
        },
        "phi_setting": phi_setting,
        "epsset": epsset,
        "particle_type": particle_type,
    }
    _print_rand_debug(
        "resolved proton-cleaning tree names",
        particle_type=particle_type,
        phi_setting=phi_setting,
        epsset=epsset,
        prompt_noRF=no_rf_prompt_name,
        rand_noRF=no_rf_rand_name,
        prompt_RF=rf_prompt_name,
        rand_RF=rf_rand_name,
    )
    _debug_tree_status("kaon_clean_prompt_noRF", bundle["sources"]["prompt"]["tree"])
    _debug_tree_status("kaon_clean_rand_noRF", bundle["sources"]["rand"]["tree"])
    _debug_tree_status("kaon_clean_dummy_prompt_noRF", bundle["sources"]["dummy_prompt"]["tree"])
    _debug_tree_status("kaon_clean_dummy_rand_noRF", bundle["sources"]["dummy_rand"]["tree"])
    _debug_tree_status("kaon_clean_prompt_RF", bundle["rf_sources"]["prompt"]["tree"])
    _debug_tree_status("kaon_clean_rand_RF", bundle["rf_sources"]["rand"]["tree"])
    _debug_tree_status("kaon_clean_dummy_prompt_RF", bundle["rf_sources"]["dummy_prompt"]["tree"])
    _debug_tree_status("kaon_clean_dummy_rand_RF", bundle["rf_sources"]["dummy_rand"]["tree"])
    return bundle


def _create_rand_sub_component_templates(target_hists):
    return {
        key: clone_reset_hist(hist, "_pi_component_template")
        for key, hist in target_hists.items()
        if hist is not None
    }


def _process_component_weighted_subtracted_particle_tree(
    tree,
    mm_offset_data,
    template_hists,
    particle_type,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    source_coeff,
    pion_reference_hist,
    pion_mm_weights,
    stats=None,
    tree_label=None,
    t_edges=None,
    mm_only=False,
):
    if tree is None:
        raise RuntimeError("Subtracted-particle tree '{}' is None".format(tree_label or "unnamed"))

    mm_offset_correction = 0.0 if get_shift_mode() == "shifted" else mm_offset_data

    for evt in tree:
        base_all_cuts, base_nomm_cuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max, mm_offset=mm_offset_correction)

        hole_rejected = (
            hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            if hole_contains is not None
            else False
        )
        pid_pass = evt.P_hgcer_npeSum > 2.0 if particle_type == "kaon" else True
        allcuts = base_all_cuts and not hole_rejected and pid_pass
        nommcuts = base_nomm_cuts and not hole_rejected and pid_pass

        if not (allcuts or nommcuts):
            continue

        adj_MM = get_shifted_mm(evt, mm_offset=mm_offset_correction)
        adj_t = shifted_t_getter(evt)
        if t_edges is not None:
            t_low, t_high = float(t_edges[0]), float(t_edges[1])
            if adj_t < t_low or adj_t >= t_high:
                continue
        pion_weight = source_coeff * simc_shape_pion_weight_from_value(
            adj_MM,
            pion_reference_hist,
            pion_mm_weights,
        )
        if pion_weight == 0.0 or not math.isfinite(pion_weight):
            continue

        if nommcuts:
            if "h_mm_nosub" in template_hists:
                template_hists["h_mm_nosub"].Fill(adj_MM, pion_weight)
            if "h_mm_fit2sub" in template_hists:
                template_hists["h_mm_fit2sub"].Fill(adj_MM, pion_weight)
            if "h_mm_fit1sub" in template_hists:
                template_hists["h_mm_fit1sub"].Fill(adj_MM, pion_weight)
            if "h_mm_pisub" in template_hists:
                template_hists["h_mm_pisub"].Fill(adj_MM, pion_weight)
            if "h_mm_full" in template_hists:
                template_hists["h_mm_full"].Fill(adj_MM, pion_weight)
            if stats is not None:
                stats["n_events_nommcuts"] += 1

        if not allcuts:
            continue

        if mm_only:
            # Parent diagnostics intentionally own only MM targets.  Do not
            # route this sparse map through the full HGCer/detector filler.
            if template_hists.get("h_mm") is not None:
                template_hists["h_mm"].Fill(adj_MM, pion_weight)
        else:
            _fill_rand_sub_allcuts_weighted(
                evt, adj_MM, adj_t, adj_hsdelta, pion_weight, template_hists
            )
        if stats is not None:
            stats["n_events_allcuts"] += 1
            stats["sum_event_weight"] += float(pion_weight)
            stats["sum_event_weight_sq"] += float(pion_weight * pion_weight)


def _canonical_t_index(value, t_bins):
    index = find_canonical_bin(value, t_bins)
    return int(index) if index >= 0 else None


def _build_authoritative_pion_control_source_cache(
    sub_tree_bundle,
    *,
    proton_t_products,
    t_bins,
    phi_bins,
    particle_type,
    pol,
    mm_offset_data,
    coordinate_contract,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    norm_factor_data,
    norm_factor_dummy,
    n_windows,
    delta_edges=None,
):
    """Build the one signed pion-control population used by t-bin parents.

    The records deliberately contain only immutable scalar source information.
    A parent fit consumes the cached t projection; later child consumers may
    project the same record definition locally without applying any proton
    factor to pion-control events.
    """
    coordinate_contract = validate_kaon_data_coordinate_contract(coordinate_contract)
    coordinate_fingerprint = coordinate_contract["coordinate_fingerprint"]
    coordinate_name_suffix = coordinate_fingerprint[-12:]
    products = list(proton_t_products or ())
    if len(products) != max(0, len(t_bins) - 1):
        raise RuntimeError("pion_control_cache_requires_one_proton_product_per_canonical_t")
    if not isinstance(sub_tree_bundle, dict):
        raise RuntimeError("pion_control_cache_requires_subtracted_particle_tree_bundle")
    pid_contract = fingerprint_hgcer_pid_contract()
    physical_pion_control_mask = dict(pid_contract["masks"]["physical_pion_control"])
    resolved_delta_edges = tuple(float(edge) for edge in (delta_edges or ()))
    if resolved_delta_edges and (
        len(resolved_delta_edges) < 2
        or any(resolved_delta_edges[index] >= resolved_delta_edges[index + 1]
               for index in range(len(resolved_delta_edges) - 1))
    ):
        raise RuntimeError("pion_control_cache_requires_increasing_delta_edges")

    child_cache_fields = (
        "source_label", "entry_index", "coefficient", "coordinate_fingerprint",
        "raw_t", "raw_MM", "adj_t", "adj_MM", "theta_cm_deg", "Q2", "W", "epsilon",
        "ssxptar", "ssyptar", "hsxptar", "hsyptar", "allcuts",
        "nommcuts", "t_index", "phi_index",
    )
    child_cache = {
        source: {field: [] for field in child_cache_fields}
        for source in ("prompt", "rand", "dummy", "dummy_rand")
    }
    cache = {
        "by_t": [], "source_accounting": {}, "records": [],
        "child_event_cache": child_cache,
        "kaon_data_coordinate": dict(coordinate_contract),
        "coordinate_fingerprint": coordinate_fingerprint,
        # This is the one source of truth for the already-authoritative
        # downstream control selection and for the Part-2 denominator.
        "physical_pion_control_mask": physical_pion_control_mask,
        "physical_pion_control_mask_fingerprint": pid_contract["fingerprint"],
        "pid_contract": pid_contract,
        "delta_edges": list(resolved_delta_edges),
    }
    for t_index, product in enumerate(products):
        final_targets = product.get("final_targets") or {}
        full_source = final_targets.get("h_mm_nosub")
        cut_source = final_targets.get("h_mm")
        if full_source is None or cut_source is None:
            raise RuntimeError("pion_control_cache_missing_proton_final_mm_target:t{}".format(t_index + 1))
        cache["by_t"].append({
            "t_index": int(t_index),
            "t_edges": [float(t_bins[t_index]), float(t_bins[t_index + 1])],
            "H_pion_control": clone_root_histogram(
                full_source,
                scope="pion_control_t{}".format(t_index + 1),
                role="authoritative_parent_full",
                name="H_MM_pion_control_parent_t{}_full".format(t_index + 1),
                reset=True,
                sumw2=True,
            ),
            "H_pion_control_cut": clone_root_histogram(
                cut_source,
                scope="pion_control_t{}".format(t_index + 1),
                role="authoritative_parent_cut",
                name="H_MM_pion_control_parent_t{}_cut".format(t_index + 1),
                reset=True,
                sumw2=True,
            ),
            "records": [],
            "source_accounting": {},
            "kaon_data_coordinate": dict(coordinate_contract),
            "coordinate_fingerprint": coordinate_fingerprint,
            "physical_pion_control_mask": dict(physical_pion_control_mask),
            "physical_pion_control_mask_fingerprint": pid_contract["fingerprint"],
        })

    # Diagnostic-only source overlays are detached from every production fit.
    # They make a raw-versus-analysis-frame mismatch visible without adding a
    # second control loop or changing any source normalization.
    source_coordinate_hists = {}
    if cache["by_t"]:
        for source_label in ("prompt", "rand", "dummy", "dummy_rand"):
            source_coordinate_hists[source_label] = {
                "H_MM_raw": clone_root_histogram(
                    cache["by_t"][0]["H_pion_control"],
                    scope="pion_coordinate_{}".format(source_label),
                    role="raw_mm_display",
                    name="H_MM_pion_coordinate_{}_raw_{}".format(
                        source_label, coordinate_name_suffix
                    ),
                    reset=True,
                    sumw2=True,
                ),
                "H_MM_analysis": clone_root_histogram(
                    cache["by_t"][0]["H_pion_control"],
                    scope="pion_coordinate_{}".format(source_label),
                    role="analysis_mm_display",
                    name="H_MM_pion_coordinate_{}_analysis_{}".format(
                        source_label, coordinate_name_suffix
                    ),
                    reset=True,
                    sumw2=True,
                ),
                "H_t_raw": TH1D(
                    "H_t_pion_coordinate_{}_raw_{}".format(
                        source_label, coordinate_name_suffix
                    ),
                    "Pion control raw |t|;|t| [GeV^2];signed yield",
                    240, float(t_bins[0]), float(t_bins[-1]),
                ),
                "H_t_analysis": TH1D(
                    "H_t_pion_coordinate_{}_analysis_{}".format(
                        source_label, coordinate_name_suffix
                    ),
                    "Pion control analysis |t|;|t| [GeV^2];signed yield",
                    240, float(t_bins[0]), float(t_bins[-1]),
                ),
            }
            for histogram in source_coordinate_hists[source_label].values():
                histogram.SetDirectory(0)
                histogram.Sumw2()
    cache["coordinate_diagnostics"] = source_coordinate_hists

    source_specs = (
        ("prompt", sub_tree_bundle.get("prompt_tree"), float(norm_factor_data), sub_tree_bundle.get("prompt_tree_name")),
        ("rand", sub_tree_bundle.get("rand_tree"), -float(norm_factor_data) / float(n_windows), sub_tree_bundle.get("rand_tree_name")),
        ("dummy", sub_tree_bundle.get("dummy_prompt_tree"), -float(norm_factor_dummy), sub_tree_bundle.get("prompt_tree_name")),
        ("dummy_rand", sub_tree_bundle.get("dummy_rand_tree"), float(norm_factor_dummy) / float(n_windows), sub_tree_bundle.get("rand_tree_name")),
    )
    mm_offset_correction = 0.0 if get_shift_mode() == "shifted" else float(mm_offset_data)
    for source_label, tree, coefficient, source_tree_name in source_specs:
        source_audit = cache["source_accounting"].setdefault(source_label, {
            "tree_entries_seen": 0,
            "selected_records": 0,
            "records_inside_canonical_t": 0,
            "allcuts_records": 0,
            "nommcuts_records": 0,
            "signed_weight_sum": 0.0,
            "absolute_weight_support": 0.0,
            "coefficient": float(coefficient),
            "tree_name": source_tree_name,
            "rf_state": "noRF" if str(source_tree_name).endswith("_noRF") else "RF_or_unknown",
            "pid_role": "pion_pid",
            "proton_factor_scope": "none",
            "coordinate_fingerprint": coordinate_fingerprint,
            "coordinate_closure": {
                "checked_records": 0,
                "mm_residual_sum": 0.0,
                "t_residual_sum": 0.0,
                "maximum_abs_mm_residual": 0.0,
                "maximum_abs_t_residual": 0.0,
            },
            "t_bin_migration": {},
            "canonical_t_counts": {
                str(t_index): 0 for t_index in range(max(0, len(t_bins) - 1))
            },
        })
        if tree is None:
            continue
        try:
            source_audit["tree_entries_seen"] = int(tree.GetEntries())
        except Exception:
            source_audit["tree_entries_seen"] = 0
        for entry_index, evt in enumerate(tree):
            base_allcuts, base_nommcuts, _ = evaluate_event(
                evt, mm_min, mm_max, mm_offset=mm_offset_correction
            )
            hole_rejected = (
                hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
                if hole_contains is not None else False
            )
            pid_pass = (
                hgcer_mask_accepts(getattr(evt, "P_hgcer_npeSum", None), physical_pion_control_mask)
                if str(particle_type).strip().lower() == "kaon" else True
            )
            allcuts = bool(base_allcuts and not hole_rejected and pid_pass)
            nommcuts = bool(base_nommcuts and not hole_rejected and pid_pass)
            if not (allcuts or nommcuts):
                continue
            source_audit["selected_records"] += 1
            raw_mm, raw_t = raw_event_coordinates(evt)
            adj_mm = float(get_shifted_mm(evt, mm_offset=mm_offset_correction))
            adj_t = float(shifted_t_getter(evt))
            expected_mm = raw_mm + float(coordinate_contract["mm_shift"])
            expected_t = raw_t + float(coordinate_contract["t_shift"])
            closure = source_audit["coordinate_closure"]
            closure["checked_records"] += 1
            closure["mm_residual_sum"] += adj_mm - expected_mm
            closure["t_residual_sum"] += adj_t - expected_t
            closure["maximum_abs_mm_residual"] = max(
                closure["maximum_abs_mm_residual"], abs(adj_mm - expected_mm)
            )
            closure["maximum_abs_t_residual"] = max(
                closure["maximum_abs_t_residual"], abs(adj_t - expected_t)
            )
            raw_t_index = _canonical_t_index(raw_t, t_bins)
            t_index = _canonical_t_index(adj_t, t_bins)
            delta_index = _canonical_t_index(
                float(getattr(evt, "ssdelta", float("nan"))), resolved_delta_edges
            ) if resolved_delta_edges else None
            migration_key = "{}->{}".format(
                "out" if raw_t_index is None else raw_t_index,
                "out" if t_index is None else t_index,
            )
            source_audit["t_bin_migration"][migration_key] = (
                source_audit["t_bin_migration"].get(migration_key, 0) + 1
            )
            if nommcuts:
                diagnostic_hists = source_coordinate_hists.get(source_label) or {}
                diagnostic_hists.get("H_MM_raw").Fill(raw_mm, coefficient)
                diagnostic_hists.get("H_MM_analysis").Fill(adj_mm, coefficient)
                diagnostic_hists.get("H_t_raw").Fill(raw_t, coefficient)
                diagnostic_hists.get("H_t_analysis").Fill(adj_t, coefficient)
            if t_index is None:
                continue
            source_audit["records_inside_canonical_t"] += 1
            source_audit["canonical_t_counts"][str(t_index)] += 1
            source_audit["allcuts_records"] += int(allcuts)
            source_audit["nommcuts_records"] += int(nommcuts)
            source_audit["signed_weight_sum"] += float(coefficient)
            source_audit["absolute_weight_support"] += abs(float(coefficient))
            t_source_audit = cache["by_t"][t_index]["source_accounting"].setdefault(
                source_label,
                {
                    "selected_records": 0,
                    "allcuts_records": 0,
                    "nommcuts_records": 0,
                    "signed_weight_sum": 0.0,
                    "absolute_weight_support": 0.0,
                    "coefficient": float(coefficient),
                },
            )
            t_source_audit["selected_records"] += 1
            t_source_audit["allcuts_records"] += int(allcuts)
            t_source_audit["nommcuts_records"] += int(nommcuts)
            t_source_audit["signed_weight_sum"] += float(coefficient)
            t_source_audit["absolute_weight_support"] += abs(float(coefficient))
            record = {
                "source_label": source_label,
                "entry_index": int(entry_index),
                "coefficient": float(coefficient),
                "raw_MM": raw_mm,
                "raw_t": raw_t,
                "adj_MM": adj_mm,
                "adj_t": adj_t,
                "t_index": int(t_index),
                "coordinate_fingerprint": coordinate_fingerprint,
                "source_tree_name": source_tree_name,
                "rf_state": "noRF" if str(source_tree_name).endswith("_noRF") else "RF_or_unknown",
                "pid_role": "pion_pid",
                "phi": float(evt.ph_q),
                "ssdelta": float(getattr(evt, "ssdelta", float("nan"))),
                "delta_index": int(delta_index) if delta_index is not None else None,
                "P_hgcer_npeSum": float(getattr(evt, "P_hgcer_npeSum", float("nan"))),
                "P_hgcer_xAtCer": float(getattr(evt, "P_hgcer_xAtCer", float("nan"))),
                "P_hgcer_yAtCer": float(getattr(evt, "P_hgcer_yAtCer", float("nan"))),
                "allcuts": bool(allcuts),
                "nommcuts": bool(nommcuts),
                "proton_cleaning_factor": None,
            }
            cache["records"].append(record)
            cache["by_t"][t_index]["records"].append(record)
            if nommcuts:
                cache["by_t"][t_index]["H_pion_control"].Fill(adj_mm, coefficient)
            if allcuts:
                cache["by_t"][t_index]["H_pion_control_cut"].Fill(adj_mm, coefficient)
            phi_degrees = float(evt.ph_q) * (180.0 / math.pi)
            phi_index = _canonical_t_index(phi_degrees, phi_bins)
            if phi_index is not None:
                try:
                    from theta_cm import calculate_theta_cm_deg
                    theta_cm_deg = float(
                        calculate_theta_cm_deg(
                            particle_type, pol, evt.W, evt.Q2, adj_t
                        )
                    )
                except Exception:
                    theta_cm_deg = float("nan")
                child_record = {
                    "source_label": source_label,
                    "entry_index": int(entry_index),
                    "coefficient": float(coefficient),
                    "raw_MM": raw_mm,
                    "raw_t": raw_t,
                    "adj_t": adj_t,
                    "adj_MM": adj_mm,
                    "coordinate_fingerprint": coordinate_fingerprint,
                    "theta_cm_deg": theta_cm_deg,
                    "Q2": float(evt.Q2),
                    "W": float(evt.W),
                    "epsilon": float(evt.epsilon),
                    "ssxptar": float(evt.ssxptar),
                    "ssyptar": float(evt.ssyptar),
                    "hsxptar": float(evt.hsxptar),
                    "hsyptar": float(evt.hsyptar),
                    "allcuts": bool(allcuts),
                    "nommcuts": bool(nommcuts),
                    "t_index": int(t_index),
                    "phi_index": int(phi_index),
                }
                for field, value in child_record.items():
                    child_cache[source_label][field].append(value)

    cache["records"] = tuple(cache["records"])
    for entry in cache["by_t"]:
        entry["records"] = tuple(entry["records"])
    cache["by_t"] = tuple(cache["by_t"])
    frozen_child_cache = {}
    for source_label, section in child_cache.items():
        frozen = {
            field: np.asarray(
                values,
                dtype=(str if field in ("source_label", "coordinate_fingerprint") else bool if field in ("allcuts", "nommcuts") else np.int32
                       if field in ("t_index", "phi_index") else np.float64),
            )
            for field, values in section.items()
        }
        allcut_index, nommcut_index = {}, {}
        for index, (t_index, phi_index) in enumerate(
            zip(frozen["t_index"], frozen["phi_index"])
        ):
            key = (int(t_index), int(phi_index))
            if bool(frozen["allcuts"][index]):
                allcut_index.setdefault(key, []).append(int(index))
            if bool(frozen["nommcuts"][index]):
                nommcut_index.setdefault(key, []).append(int(index))
        frozen["allcut_bin_index"] = {
            key: np.asarray(indices, dtype=np.int32)
            for key, indices in allcut_index.items()
        }
        frozen["nommcut_bin_index"] = {
            key: np.asarray(indices, dtype=np.int32)
            for key, indices in nommcut_index.items()
        }
        frozen_child_cache[source_label] = frozen
    cache["child_event_cache"] = frozen_child_cache
    for source_audit in cache["source_accounting"].values():
        closure = source_audit.get("coordinate_closure") or {}
        closure["tolerance"] = 1.0e-12
        closure["passed"] = bool(
            closure.get("maximum_abs_mm_residual", float("inf")) <= closure["tolerance"]
            and closure.get("maximum_abs_t_residual", float("inf")) <= closure["tolerance"]
        )
        closure["mm_shift"] = float(coordinate_contract["mm_shift"])
        closure["t_shift"] = float(coordinate_contract["t_shift"])
    if cache["by_t"]:
        global_full = clone_root_histogram(
            cache["by_t"][0]["H_pion_control"],
            scope="pion_control_global",
            role="authoritative_setting_wide_full",
            name="H_MM_pion_control_authoritative_canonical_t_global_full",
            reset=True,
            sumw2=True,
        )
        global_cut = clone_root_histogram(
            cache["by_t"][0]["H_pion_control_cut"],
            scope="pion_control_global",
            role="authoritative_setting_wide_cut",
            name="H_MM_pion_control_authoritative_canonical_t_global_cut",
            reset=True,
            sumw2=True,
        )
        for entry in cache["by_t"]:
            global_full.Add(entry["H_pion_control"])
            global_cut.Add(entry["H_pion_control_cut"])
        cache["H_pion_control_global"] = global_full
        cache["H_pion_control_global_cut"] = global_cut
    cache["definition"] = "one_signed_pion_control_source_cache_no_proton_weight_physical_hgcer_mask_resolved"
    cache["immutable_record_contract"] = True
    return cache


def _fill_mm_only_authoritative_pion_control_templates(
    records,
    template_hists,
    pion_reference_hist,
    pion_mm_weights,
):
    """Apply a parent weight to the frozen signed pion-control records.

    This is deliberately MM-only.  Parent diagnostics receive detached MM
    products from the proton stage, so invoking the generic detector filler
    would incorrectly require HGCer and unrelated target maps.
    """
    stats = {
        "n_events_allcuts": 0,
        "n_events_nommcuts": 0,
        "sum_event_weight": 0.0,
        "sum_event_weight_sq": 0.0,
        "source_weight_sums": {},
    }
    for record in records or ():
        adj_mm = float(record["adj_MM"])
        event_weight = float(record["coefficient"]) * simc_shape_pion_weight_from_value(
            adj_mm,
            pion_reference_hist,
            pion_mm_weights,
        )
        if event_weight == 0.0 or not math.isfinite(event_weight):
            continue
        source_label = str(record.get("source_label") or "unknown")
        source_stats = stats["source_weight_sums"].setdefault(
            source_label,
            {"signed_sum": 0.0, "absolute_support": 0.0, "allcuts_records": 0, "nommcuts_records": 0},
        )
        source_stats["signed_sum"] += event_weight
        source_stats["absolute_support"] += abs(event_weight)
        if bool(record.get("nommcuts")):
            if template_hists.get("h_mm_full") is not None:
                template_hists["h_mm_full"].Fill(adj_mm, event_weight)
            stats["n_events_nommcuts"] += 1
            source_stats["nommcuts_records"] += 1
        if bool(record.get("allcuts")):
            if template_hists.get("h_mm") is not None:
                template_hists["h_mm"].Fill(adj_mm, event_weight)
            stats["n_events_allcuts"] += 1
            stats["sum_event_weight"] += event_weight
            stats["sum_event_weight_sq"] += event_weight * event_weight
            source_stats["allcuts_records"] += 1
    return stats


def _build_authoritative_parent_mm_diagnostic_proposal(
    fit_result,
    parent_input,
    *,
    inp_dict,
    scope,
):
    """Build a proposal from the cuts-layer cache without re-looping trees."""
    proposal_check = evaluate_component_pion_application_proposal(fit_result)
    if not proposal_check.get("available"):
        raise ValueError(proposal_check.get("reason") or "proposal model unavailable")
    cut_before = parent_input.get("H_proton_cleaned_final_rf_cut")
    full_before = parent_input.get("H_proton_cleaned_final_rf")
    records = parent_input.get("pion_control_records")
    if cut_before is None or full_before is None or records is None:
        raise RuntimeError("authoritative_parent_mm_diagnostic_missing_input")

    clip_min, clip_max = resolve_particle_subtraction_weight_clip_bounds(inp_dict)
    weight_payload = build_simc_shape_pion_control_weights(
        fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inp_dict),
        proposal_mode=True,
    )
    stage_weight_payload = build_simc_shape_pion_control_weights(
        fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inp_dict),
        model_variant="staged",
        proposal_mode=True,
    )
    templates = {
        "h_mm": clone_root_histogram(
            cut_before, scope=scope, role="authoritative_parent_template_cut", reset=True, sumw2=True
        ),
        "h_mm_full": clone_root_histogram(
            full_before, scope=scope, role="authoritative_parent_template_full", reset=True, sumw2=True
        ),
    }
    stats = _fill_mm_only_authoritative_pion_control_templates(
        records,
        templates,
        weight_payload["H_pion_control_model"],
        weight_payload["weights"],
    )
    before = clone_root_histogram(cut_before, scope=scope, role="authoritative_parent_before_cut")
    full_before_copy = clone_root_histogram(
        full_before, scope=scope, role="authoritative_parent_before_full"
    )
    after = clone_root_histogram(before, scope=scope, role="authoritative_parent_after_cut")
    full_after = clone_root_histogram(
        full_before_copy, scope=scope, role="authoritative_parent_after_full"
    )
    after.Add(templates["h_mm"], -1.0)
    full_after.Add(templates["h_mm_full"], -1.0)
    model_after_stage = clone_root_histogram(
        full_before_copy, scope=scope, role="authoritative_parent_after_stage_model"
    )
    model_after_final = clone_root_histogram(
        full_before_copy, scope=scope, role="authoritative_parent_after_final_model"
    )
    if stage_weight_payload.get("H_kaon_pion_model") is not None:
        model_after_stage.Add(stage_weight_payload["H_kaon_pion_model"], -1.0)
    if weight_payload.get("H_kaon_pion_model") is not None:
        model_after_final.Add(weight_payload["H_kaon_pion_model"], -1.0)
    weighted_full = _hist_integral(templates["h_mm_full"])
    control_integral = _hist_integral(parent_input.get("H_pion_control"))
    payload = {
        "accepted": True,
        "fallback_used": False,
        "fallback_reason": "",
        "proposal_available": True,
        "diagnostic_only": True,
        "application_authoritative": False,
        "production_applied": False,
        "diagnostic_role": "proposal",
        "analysis_scope": fit_result.get("analysis_scope"),
        "input_selection": str(
            parent_input.get("input_selection")
            or "no_rf_proton_cleaning_then_rf_restored"
        ),
        "source_target_state": str(
            parent_input.get("source_target_state") or "post_proton_post_rf"
        ),
        "H_pion_control_model": weight_payload["H_pion_control_model"],
        "H_kaon_pion_model": weight_payload["H_kaon_pion_model"],
        "H_weighted_pion_control_model": weight_payload.get("H_weighted_pion_control_model"),
        "H_pion_weight_vs_MM": weight_payload["H_pion_weight_vs_MM"],
        "H_pion_control_model_stage": stage_weight_payload["H_pion_control_model"],
        "H_kaon_pion_model_stage": stage_weight_payload["H_kaon_pion_model"],
        "H_weighted_pion_control_model_stage": stage_weight_payload.get("H_weighted_pion_control_model"),
        "H_pion_weight_vs_MM_stage": stage_weight_payload["H_pion_weight_vs_MM"],
        "weights": weight_payload["weights"],
        "H_pion_subtraction_template_MM": templates["h_mm"],
        "H_pion_subtraction_template_MM_nosub": templates["h_mm_full"],
        "H_MM_before_pion_subtraction": before,
        "H_MM_after_pion_subtraction": after,
        "H_MM_nosub_before_pion_subtraction": full_before_copy,
        "H_MM_nosub_after_pion_subtraction": full_after,
        "H_MM_nosub_after_pion_subtraction_model_stage": model_after_stage,
        "H_MM_nosub_after_pion_subtraction_model_final": model_after_final,
        "weighted_pion_integral": weighted_full,
        "weighted_pion_integral_cut": _hist_integral(templates["h_mm"]),
        "weighted_pion_integral_full": weighted_full,
        "particle_subtraction_effective_scale": (
            weighted_full / control_integral if control_integral > 0.0 else 0.0
        ),
        "kaon_integral_before_pion_sub": _hist_integral(before),
        "kaon_integral_after_pion_sub": _hist_integral(after),
        "kaon_integral_before_pion_sub_full": _hist_integral(full_before_copy),
        "kaon_integral_after_pion_sub_full": _hist_integral(full_after),
        "diagnostics": {
            **dict(weight_payload.get("diagnostics") or {}),
            "weight_diagnostics_stage": dict(stage_weight_payload.get("diagnostics") or {}),
            **stats,
            "source_definition": "authoritative_pion_control_cache",
            "event_template_closure": compute_hist_closure_metrics(
                weight_payload.get("H_kaon_pion_model"), templates["h_mm_full"]
            ),
            "final_closure": compute_hist_closure_metrics(
                full_before_copy, full_after
            ),
        },
    }
    assert_component_subtraction_payload_ownership(payload)
    return payload


def _build_authoritative_parent_single_scale_final(
    proposal_payload,
    production_evaluation,
    parent_input,
    *,
    inp_dict,
    phi_setting,
    mm_offset_data,
    scope,
):
    """Apply the established single-scale fallback to the same frozen records."""
    before = proposal_payload.get("H_MM_before_pion_subtraction")
    full_before = proposal_payload.get("H_MM_nosub_before_pion_subtraction")
    reference = proposal_payload.get("H_pion_control_model")
    if before is None or full_before is None or reference is None:
        raise RuntimeError("missing_source_histogram")
    templates = {
        "h_mm": clone_root_histogram(before, scope=scope, role="single_scale_cached_template_cut", reset=True, sumw2=True),
        "h_mm_full": clone_root_histogram(full_before, scope=scope, role="single_scale_cached_template_full", reset=True, sumw2=True),
    }
    unit_weights = np.ones(int(reference.GetNbinsX()) + 2, dtype=np.float64)
    stats = _fill_mm_only_authoritative_pion_control_templates(
        parent_input.get("pion_control_records"), templates, reference, unit_weights
    )
    windows = resolve_particle_subtraction_windows(
        "kaon", "pion", mm_offset_data, inp_dict=inp_dict, phi_setting=phi_setting
    )
    try:
        scale_components = compute_staged_particle_subtraction_scales(
            full_before,
            templates["h_mm_full"],
            windows,
            context="parent diagnostic single-scale ({})".format(phi_setting),
        )
        scale_factor = float(scale_components["total_scale_factor"])
    except ZeroDivisionError:
        scale_components = None
        scale_factor = 0.0
    for histogram in templates.values():
        histogram.Scale(scale_factor)
    after = clone_root_histogram(before, scope=scope, role="single_scale_cached_after_cut")
    full_after = clone_root_histogram(full_before, scope=scope, role="single_scale_cached_after_full")
    after.Add(templates["h_mm"], -1.0)
    full_after.Add(templates["h_mm_full"], -1.0)
    final_weight = clone_root_histogram(reference, scope=scope, role="single_scale_cached_weight", reset=True)
    for bin_index in range(1, final_weight.GetNbinsX() + 1):
        final_weight.SetBinContent(bin_index, scale_factor)
        final_weight.SetBinError(bin_index, 0.0)
    expected_after = clone_root_histogram(
        full_before, scope=scope, role="single_scale_cached_expected_after"
    )
    expected_after.Add(templates["h_mm_full"], -1.0)
    final = {
        "accepted": True,
        "proposal_available": False,
        "diagnostic_role": "final",
        "final_application_status": "applied_fallback",
        "production_evaluation_accepted": False,
        "production_rejection_reasons": [
            reason.strip() for reason in str(production_evaluation.get("reason") or "").split(";") if reason.strip()
        ],
        "fallback_used": True,
        "fallback_mode": "single_scale",
        "fallback_reason": production_evaluation.get("reason") or "component-fit result rejected",
        "diagnostic_only": True,
        "application_authoritative": False,
        "production_applied": False,
        "analysis_scope": proposal_payload.get("analysis_scope"),
        "input_selection": proposal_payload.get("input_selection"),
        "source_target_state": proposal_payload.get("source_target_state"),
        "H_pion_control_model": clone_root_histogram(reference, scope=scope, role="single_scale_cached_control"),
        "H_kaon_pion_model": clone_root_histogram(templates["h_mm_full"], scope=scope, role="single_scale_cached_model"),
        "H_weighted_pion_control_model": clone_root_histogram(templates["h_mm_full"], scope=scope, role="single_scale_cached_weighted_model"),
        "H_pion_weight_vs_MM": final_weight,
        "weights": np.full(int(reference.GetNbinsX()) + 2, scale_factor, dtype=np.float64),
        "H_pion_subtraction_template_MM": templates["h_mm"],
        "H_pion_subtraction_template_MM_nosub": templates["h_mm_full"],
        "H_MM_before_pion_subtraction": clone_root_histogram(before, scope=scope, role="single_scale_cached_before"),
        "H_MM_after_pion_subtraction": after,
        "H_MM_nosub_before_pion_subtraction": clone_root_histogram(full_before, scope=scope, role="single_scale_cached_full_before"),
        "H_MM_nosub_after_pion_subtraction": full_after,
        "H_MM_nosub_after_pion_subtraction_model_stage": clone_root_histogram(full_after, scope=scope, role="single_scale_cached_stage_model"),
        "H_MM_nosub_after_pion_subtraction_model_final": clone_root_histogram(full_after, scope=scope, role="single_scale_cached_final_model"),
        "particle_subtraction_effective_scale": scale_factor,
        "weighted_pion_integral": _hist_integral(templates["h_mm_full"]),
        "weighted_pion_integral_cut": _hist_integral(templates["h_mm"]),
        "weighted_pion_integral_full": _hist_integral(templates["h_mm_full"]),
        "kaon_integral_before_pion_sub": _hist_integral(before),
        "kaon_integral_after_pion_sub": _hist_integral(after),
        "kaon_integral_before_pion_sub_full": _hist_integral(full_before),
        "kaon_integral_after_pion_sub_full": _hist_integral(full_after),
        "diagnostics": {
            **stats,
            "source_definition": "authoritative_pion_control_cache",
            "fallback_mode": "single_scale",
            "scale_components": scale_components,
            "event_template_closure": compute_hist_closure_metrics(templates["h_mm_full"], templates["h_mm_full"]),
            "final_closure": compute_hist_closure_metrics(expected_after, full_after),
        },
    }
    assert_component_subtraction_payload_ownership(final)
    return final


def _process_rand_sub_background_tree(
    tree,
    tmin,
    tmax,
    template_hists,
    particle_type,
    hole_contains,
    evaluate_event,
    mm_min,
    mm_max,
    mm_reference_hist,
    mm_background_weights,
    source_coeff,
    residual_weights=None,
    tree_label=None,
):
    if tree is None:
        raise RuntimeError("Background tree '{}' is None".format(tree_label or "unnamed"))

    for evt in tree:
        base_all_cuts, _, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
        hole_rejected = (
            hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            if hole_contains is not None
            else False
        )
        allcuts = base_all_cuts and not hole_rejected

        if not allcuts:
            continue

        adj_MM = get_shifted_mm(evt)
        adj_t = get_shifted_t(evt)

        event_weight = source_coeff * mm_background_weight_from_value(
            adj_MM,
            mm_reference_hist,
            mm_background_weights,
            residual_weights=residual_weights,
        )
        if event_weight == 0.0:
            continue

        _fill_rand_sub_allcuts_weighted(evt, adj_MM, adj_t, adj_hsdelta, event_weight, template_hists)


def _process_subtracted_particle_background_tree(
    tree,
    mm_offset_data,
    template_hists,
    particle_type,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    mm_reference_hist,
    mm_background_weights,
    source_coeff,
    pion_reference_hist=None,
    pion_mm_weights=None,
    residual_weights=None,
    tree_label=None,
):
    if tree is None:
        raise RuntimeError("Subtracted-particle background tree '{}' is None".format(tree_label or "unnamed"))

    mm_offset_correction = 0.0 if get_shift_mode() == "shifted" else mm_offset_data

    for evt in tree:
        base_all_cuts, _, adj_hsdelta = evaluate_event(evt, mm_min, mm_max, mm_offset=mm_offset_correction)

        hole_rejected = (
            hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            if hole_contains is not None
            else False
        )
        pid_pass = evt.P_hgcer_npeSum > 2.0 if particle_type == "kaon" else True
        allcuts = base_all_cuts and not hole_rejected and pid_pass

        if not allcuts:
            continue

        adj_MM = get_shifted_mm(evt, mm_offset=mm_offset_correction)
        adj_t = shifted_t_getter(evt)

        event_weight = source_coeff * mm_background_weight_from_value(
            adj_MM,
            mm_reference_hist,
            mm_background_weights,
            residual_weights=residual_weights,
        )
        if pion_reference_hist is not None and pion_mm_weights is not None:
            event_weight *= simc_shape_pion_weight_from_value(
                adj_MM,
                pion_reference_hist,
                pion_mm_weights,
            )
        if event_weight == 0.0:
            continue

        _fill_rand_sub_allcuts_weighted(evt, adj_MM, adj_t, adj_hsdelta, event_weight, template_hists)


def _apply_component_pion_subtraction_setting(
    component_fit_result,
    sub_tree_bundle,
    phi_setting,
    inpDict,
    particle_type,
    mm_offset_data,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    norm_factor_data,
    norm_factor_dummy,
    nWindows,
    data_targets,
    diagnostic_only=False,
    input_selection="no_rf_proton_cleaning_then_rf_restored",
    source_target_state="post_proton_post_rf",
    t_edges=None,
    proposal_only=False,
    mm_only=False,
):
    gate_result = evaluate_particle_subtraction_component_fit_result(component_fit_result, inpDict)
    payload = {
        "accepted": False,
        "fallback_used": True,
        "fallback_mode": gate_result.get("fallback_mode"),
        "fallback_reason": gate_result.get("reason") or "component-fit result rejected",
        "analysis_scope": component_fit_result.get("analysis_scope") if isinstance(component_fit_result, dict) else None,
        "particle_subtraction_mode": component_fit_result.get("particle_subtraction_mode") if isinstance(component_fit_result, dict) else None,
        "particle_subtraction_setting_key": component_fit_result.get("particle_subtraction_setting_key") if isinstance(component_fit_result, dict) else None,
        "particle_subtraction_phi_setting": component_fit_result.get("particle_subtraction_phi_setting") if isinstance(component_fit_result, dict) else None,
        "resolved_subtraction_config": deepcopy(
            component_fit_result.get("resolved_subtraction_config") or {}
        ) if isinstance(component_fit_result, dict) else {},
        "fit_mode": component_fit_result.get("fit_mode") if isinstance(component_fit_result, dict) else None,
        "fit_mode_pion": component_fit_result.get("fit_mode_pion") if isinstance(component_fit_result, dict) else None,
        "fit_mode_kaon": component_fit_result.get("fit_mode_kaon") if isinstance(component_fit_result, dict) else None,
        "fit_status_pion": component_fit_result.get("fit_status_pion") if isinstance(component_fit_result, dict) else None,
        "fit_status_kaon": component_fit_result.get("fit_status_kaon") if isinstance(component_fit_result, dict) else None,
        "fit_validation_pion": bool((gate_result.get("diagnostics") or {}).get("fit_validation_pion")),
        "fit_validation_kaon": bool((gate_result.get("diagnostics") or {}).get("fit_validation_kaon")),
        "diagnostic_only": bool(diagnostic_only),
        "application_authoritative": not bool(diagnostic_only),
        "production_applied": False,
        "proposal_available": False,
        "diagnostic_role": "proposal" if proposal_only else "production",
        "production_evaluation_accepted": bool(gate_result.get("accepted")),
        "production_rejection_reasons": [
            reason.strip()
            for reason in str(gate_result.get("reason") or "").split(";")
            if reason.strip()
        ],
        "input_selection": str(input_selection),
        "source_target_state": str(source_target_state),
    }
    if not gate_result["accepted"] and not proposal_only:
        if diagnostic_only:
            assert_component_subtraction_payload_ownership(payload)
            return payload
        return handle_particle_subtraction_fallback(
            payload,
            payload["fallback_reason"],
            context="rand_sub component pion subtraction ({})".format(phi_setting),
        )
    if not isinstance(sub_tree_bundle, dict):
        if diagnostic_only:
            payload["fallback_reason"] = "missing subtraction-tree bundle for component-weight subtraction"
            assert_component_subtraction_payload_ownership(payload)
            return payload
        return handle_particle_subtraction_fallback(
            payload,
            "missing subtraction-tree bundle for component-weight subtraction",
            context="rand_sub component pion subtraction ({})".format(phi_setting),
        )

    clip_min, clip_max = resolve_particle_subtraction_weight_clip_bounds(inpDict)
    weight_payload = build_simc_shape_pion_control_weights(
        component_fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inpDict),
        proposal_mode=bool(proposal_only),
    )
    stage_weight_payload = build_simc_shape_pion_control_weights(
        component_fit_result,
        clip_min=clip_min,
        clip_max=clip_max,
        denom_floor=resolve_particle_subtraction_weight_denominator_floor(inpDict),
        model_variant="staged",
        proposal_mode=bool(proposal_only),
    )
    print_particle_subtraction_weight_support_warning(
        weight_payload,
        context="rand_sub component pion subtraction",
        epsset=inpDict.get("EPSSET", ""),
        phi_setting=phi_setting,
    )

    if weight_payload["diagnostics"]["pion_weight_max"] > resolve_particle_subtraction_weight_warn_max(inpDict):
        print(
            "WARNING: pion component weight exceeded threshold\n"
            "  epsset = {}\n"
            "  phi_setting = {}\n"
            "  max_weight = {:.4f}".format(
                inpDict.get("EPSSET", ""),
                phi_setting,
                float(weight_payload["diagnostics"]["pion_weight_max"]),
            )
        )

    template_hists = _create_rand_sub_component_templates(data_targets)
    stats = {
        "n_events_allcuts": 0,
        "n_events_nommcuts": 0,
        "sum_event_weight": 0.0,
        "sum_event_weight_sq": 0.0,
    }
    source_specs = [
        ("prompt", sub_tree_bundle.get("prompt_tree"), float(norm_factor_data)),
        ("rand", sub_tree_bundle.get("rand_tree"), -float(norm_factor_data) / float(nWindows)),
        ("dummy", sub_tree_bundle.get("dummy_prompt_tree"), -float(norm_factor_dummy)),
        ("dummy_rand", sub_tree_bundle.get("dummy_rand_tree"), float(norm_factor_dummy) / float(nWindows)),
    ]

    for label, tree, coeff in source_specs:
        _process_component_weighted_subtracted_particle_tree(
            tree,
            mm_offset_data,
            template_hists,
            particle_type,
            hole_contains,
            evaluate_event,
            shifted_t_getter,
            mm_min,
            mm_max,
            coeff,
            weight_payload["H_pion_control_model"],
            weight_payload["weights"],
            stats=stats,
            tree_label="component {}".format(label),
            t_edges=t_edges,
            mm_only=bool(mm_only),
        )

    h_mm_before = clone_reset_hist(data_targets["h_mm"], "_before_pion_sub")
    h_mm_before.Add(data_targets["h_mm"])
    h_mm_full_before = clone_reset_hist(data_targets["h_mm_full"], "_before_pion_sub_full")
    h_mm_full_before.Add(data_targets["h_mm_full"])
    h_mm_full_after_stage_model = clone_reset_hist(h_mm_full_before, "_after_pion_sub_stage_model_full")
    h_mm_full_after_stage_model.Add(h_mm_full_before)
    if stage_weight_payload.get("H_kaon_pion_model") is not None:
        h_mm_full_after_stage_model.Add(stage_weight_payload["H_kaon_pion_model"], -1.0)
    h_mm_full_after_final_model = clone_reset_hist(h_mm_full_before, "_after_pion_sub_final_model_full")
    h_mm_full_after_final_model.Add(h_mm_full_before)
    if weight_payload.get("H_kaon_pion_model") is not None:
        h_mm_full_after_final_model.Add(weight_payload["H_kaon_pion_model"], -1.0)

    for key, target_hist in data_targets.items():
        template_hist = template_hists.get(key)
        if target_hist is None or template_hist is None:
            continue
        target_hist.Add(template_hist, -1.0)

    h_mm_after = clone_reset_hist(data_targets["h_mm"], "_after_pion_sub")
    h_mm_after.Add(data_targets["h_mm"])
    h_mm_full_after = clone_reset_hist(data_targets["h_mm_full"], "_after_pion_sub_full")
    h_mm_full_after.Add(data_targets["h_mm_full"])

    pion_control_integral = _hist_integral(component_fit_result.get("H_pion_control_input"))
    weighted_pion_integral_cut = _hist_integral(template_hists.get("h_mm"))
    weighted_pion_integral_full = _hist_integral(template_hists.get("h_mm_full"))
    effective_scale = weighted_pion_integral_full / pion_control_integral if pion_control_integral > 0.0 else 0.0
    event_template_closure = compute_hist_closure_metrics(
        weight_payload.get("H_kaon_pion_model"),
        template_hists.get("h_mm_full"),
    )

    payload.update(
        {
            "accepted": True,
            "fallback_used": False,
            "fallback_reason": "",
            "proposal_available": True,
            "production_applied": not bool(diagnostic_only) and not bool(proposal_only),
            "particle_subtraction_effective_scale": float(effective_scale),
            "weighted_pion_integral": float(weighted_pion_integral_full),
            "weighted_pion_integral_cut": float(weighted_pion_integral_cut),
            "weighted_pion_integral_full": float(weighted_pion_integral_full),
            "kaon_integral_before_pion_sub": _hist_integral(h_mm_before),
            "kaon_integral_after_pion_sub": _hist_integral(h_mm_after),
            "kaon_integral_before_pion_sub_full": _hist_integral(h_mm_full_before),
            "kaon_integral_after_pion_sub_full": _hist_integral(h_mm_full_after),
            "H_pion_control_model": weight_payload["H_pion_control_model"],
            "H_kaon_pion_model": weight_payload["H_kaon_pion_model"],
            "H_weighted_pion_control_model": weight_payload.get("H_weighted_pion_control_model"),
            "H_pion_weight_vs_MM": weight_payload["H_pion_weight_vs_MM"],
            "H_pion_control_model_stage": stage_weight_payload["H_pion_control_model"],
            "H_kaon_pion_model_stage": stage_weight_payload["H_kaon_pion_model"],
            "H_weighted_pion_control_model_stage": stage_weight_payload.get("H_weighted_pion_control_model"),
            "H_pion_weight_vs_MM_stage": stage_weight_payload["H_pion_weight_vs_MM"],
            "weights": weight_payload["weights"],
            "H_pion_subtraction_template_MM": template_hists.get("h_mm"),
            "H_pion_subtraction_template_MM_nosub": template_hists.get("h_mm_full"),
            "H_MM_before_pion_subtraction": h_mm_before,
            "H_MM_after_pion_subtraction": h_mm_after,
            "H_MM_nosub_before_pion_subtraction": h_mm_full_before,
            "H_MM_nosub_after_pion_subtraction": h_mm_full_after,
            "H_MM_nosub_after_pion_subtraction_model_stage": h_mm_full_after_stage_model,
            "H_MM_nosub_after_pion_subtraction_model_final": h_mm_full_after_final_model,
            "diagnostics": {
                **dict(weight_payload["diagnostics"]),
                "weight_diagnostics_stage": dict(stage_weight_payload.get("diagnostics") or {}),
                "model_closure_stage": dict((stage_weight_payload.get("diagnostics") or {}).get("model_closure") or {}),
                **dict(stats),
                "event_template_closure": event_template_closure,
            },
        }
    )
    assert_component_subtraction_payload_ownership(payload)
    return payload


def build_component_pion_application_proposal(
    component_fit_result,
    sub_tree_bundle,
    phi_setting,
    inpDict,
    particle_type,
    mm_offset_data,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    norm_factor_data,
    norm_factor_dummy,
    nWindows,
    detached_targets,
    *,
    input_selection="no_rf_proton_cleaning_then_rf_restored",
    source_target_state="post_proton_post_rf",
    t_edges=None,
    mm_only=False,
):
    """Build an evaluable, detached parent diagnostic proposal.

    Production acceptance is deliberately not consulted here.  The proposal
    still rejects malformed/non-finite component models, but a fit-quality
    rejection remains visible as a proposed correction rather than becoming an
    empty diagnostic payload.
    """
    proposal_check = evaluate_component_pion_application_proposal(component_fit_result)
    if not proposal_check.get("available"):
        raise ValueError(proposal_check.get("reason") or "proposal model unavailable")
    return _apply_component_pion_subtraction_setting(
        component_fit_result,
        sub_tree_bundle,
        phi_setting,
        inpDict,
        particle_type,
        mm_offset_data,
        hole_contains,
        evaluate_event,
        shifted_t_getter,
        mm_min,
        mm_max,
        norm_factor_data,
        norm_factor_dummy,
        nWindows,
        detached_targets,
        diagnostic_only=True,
        input_selection=input_selection,
        source_target_state=source_target_state,
        t_edges=t_edges,
        proposal_only=True,
        mm_only=bool(mm_only),
    )


def _clone_parent_diagnostic_histogram(histogram, scope, role, *, reset=False):
    if histogram is None:
        return None
    return clone_root_histogram(
        histogram,
        scope=scope,
        role=role,
        reset=reset,
    )


def _build_zero_parent_diagnostic_final(proposal_payload, production_evaluation, scope):
    """Represent the configured zero fallback on detached before spectra."""
    before = proposal_payload.get("H_MM_before_pion_subtraction")
    full_before = proposal_payload.get("H_MM_nosub_before_pion_subtraction")
    if before is None or full_before is None:
        raise RuntimeError("missing_source_histogram")
    zero_cut = _clone_parent_diagnostic_histogram(
        before, scope, "zero_fallback_template_cut", reset=True
    )
    zero_full = _clone_parent_diagnostic_histogram(
        full_before, scope, "zero_fallback_template_full", reset=True
    )
    final_before = _clone_parent_diagnostic_histogram(before, scope, "zero_fallback_before")
    final_full_before = _clone_parent_diagnostic_histogram(
        full_before, scope, "zero_fallback_full_before")
    final_after = _clone_parent_diagnostic_histogram(before, scope, "zero_fallback_after")
    final_full_after = _clone_parent_diagnostic_histogram(full_before, scope, "zero_fallback_full_after")
    final_weight = _clone_parent_diagnostic_histogram(
        proposal_payload.get("H_pion_weight_vs_MM") or before,
        scope,
        "zero_fallback_weight",
        reset=True,
    )
    final_control_model = _clone_parent_diagnostic_histogram(
        proposal_payload.get("H_pion_control_model") or before,
        scope,
        "zero_fallback_control_model",
        reset=True,
    )
    final_kaon_model = _clone_parent_diagnostic_histogram(
        proposal_payload.get("H_kaon_pion_model") or before,
        scope,
        "zero_fallback_kaon_model",
        reset=True,
    )
    final_weighted_model = _clone_parent_diagnostic_histogram(
        proposal_payload.get("H_kaon_pion_model") or before,
        scope,
        "zero_fallback_weighted_model",
        reset=True,
    )
    final = {
        "accepted": True,
        "proposal_available": False,
        "diagnostic_role": "final",
        "final_application_status": "zero",
        "production_evaluation_accepted": False,
        "production_rejection_reasons": [
            reason.strip()
            for reason in str(production_evaluation.get("reason") or "").split(";")
            if reason.strip()
        ],
        "fallback_used": True,
        "fallback_mode": "zero",
        "fallback_reason": production_evaluation.get("reason") or "component-fit result rejected",
        "diagnostic_only": True,
        "application_authoritative": False,
        "production_applied": False,
        "analysis_scope": proposal_payload.get("analysis_scope"),
        "input_selection": proposal_payload.get("input_selection"),
        "source_target_state": proposal_payload.get("source_target_state"),
        "H_pion_control_model": final_control_model,
        "H_kaon_pion_model": final_kaon_model,
        "H_weighted_pion_control_model": final_weighted_model,
        "H_pion_weight_vs_MM": final_weight,
        "weights": np.zeros(int(final_weight.GetNbinsX()) + 2, dtype=np.float64),
        "H_pion_subtraction_template_MM": zero_cut,
        "H_pion_subtraction_template_MM_nosub": zero_full,
        "H_MM_before_pion_subtraction": final_before,
        "H_MM_after_pion_subtraction": final_after,
        "H_MM_nosub_before_pion_subtraction": final_full_before,
        "H_MM_nosub_after_pion_subtraction": final_full_after,
        "H_MM_nosub_after_pion_subtraction_model_stage": _clone_parent_diagnostic_histogram(
            full_before, scope, "zero_fallback_stage_model"
        ),
        "H_MM_nosub_after_pion_subtraction_model_final": _clone_parent_diagnostic_histogram(
            full_before, scope, "zero_fallback_final_model"
        ),
        "weighted_pion_integral": 0.0,
        "weighted_pion_integral_cut": 0.0,
        "weighted_pion_integral_full": 0.0,
        "kaon_integral_before_pion_sub": _hist_integral(final_before),
        "kaon_integral_after_pion_sub": _hist_integral(final_after),
        "kaon_integral_before_pion_sub_full": _hist_integral(final_full_before),
        "kaon_integral_after_pion_sub_full": _hist_integral(final_full_after),
        "diagnostics": {
            "final_closure": compute_hist_closure_metrics(final_full_before, final_full_after),
            "event_template_closure": compute_hist_closure_metrics(zero_full, zero_full),
            "fallback_mode": "zero",
        },
    }
    assert_component_subtraction_payload_ownership(final)
    return final


def _build_single_scale_parent_diagnostic_final(
    proposal_payload,
    production_evaluation,
    sub_tree_bundle,
    phi_setting,
    inpDict,
    particle_type,
    mm_offset_data,
    hole_contains,
    evaluate_event,
    shifted_t_getter,
    mm_min,
    mm_max,
    norm_factor_data,
    norm_factor_dummy,
    nWindows,
    t_edges,
    scope,
):
    """Build the real legacy scalar fallback for one detached parent scope."""
    before = proposal_payload.get("H_MM_before_pion_subtraction")
    full_before = proposal_payload.get("H_MM_nosub_before_pion_subtraction")
    reference = proposal_payload.get("H_pion_control_model")
    if before is None or full_before is None or reference is None:
        raise RuntimeError("missing_source_histogram")
    targets = {
        "h_mm": _clone_parent_diagnostic_histogram(before, scope, "single_scale_target_cut"),
        "h_mm_full": _clone_parent_diagnostic_histogram(full_before, scope, "single_scale_target_full"),
    }
    templates = _create_rand_sub_component_templates(targets)
    unit_weights = np.ones(int(reference.GetNbinsX()) + 2, dtype=np.float64)
    stats = {"n_events_allcuts": 0, "n_events_nommcuts": 0, "sum_event_weight": 0.0, "sum_event_weight_sq": 0.0}
    source_specs = (
        ("prompt", sub_tree_bundle.get("prompt_tree"), float(norm_factor_data)),
        ("rand", sub_tree_bundle.get("rand_tree"), -float(norm_factor_data) / float(nWindows)),
        ("dummy", sub_tree_bundle.get("dummy_prompt_tree"), -float(norm_factor_dummy)),
        ("dummy_rand", sub_tree_bundle.get("dummy_rand_tree"), float(norm_factor_dummy) / float(nWindows)),
    )
    for label, tree, coefficient in source_specs:
        _process_component_weighted_subtracted_particle_tree(
            tree,
            mm_offset_data,
            templates,
            particle_type,
            hole_contains,
            evaluate_event,
            shifted_t_getter,
            mm_min,
            mm_max,
            coefficient,
            reference,
            unit_weights,
            stats=stats,
            tree_label="single-scale {}".format(label),
            t_edges=t_edges,
            mm_only=True,
        )
    windows = resolve_particle_subtraction_windows(
        particle_type,
        "pion",
        mm_offset_data,
        inp_dict=inpDict,
        phi_setting=phi_setting,
    )
    try:
        scale_components = compute_staged_particle_subtraction_scales(
            full_before,
            templates["h_mm_full"],
            windows,
            context="parent diagnostic single-scale ({})".format(phi_setting),
        )
        scale_factor = float(scale_components["total_scale_factor"])
    except ZeroDivisionError:
        scale_components = None
        scale_factor = 0.0
    for histogram in templates.values():
        histogram.Scale(scale_factor)
    for key, target in targets.items():
        target.Add(templates[key], -1.0)
    final_weight = _clone_parent_diagnostic_histogram(reference, scope, "single_scale_weight", reset=True)
    for bin_index in range(1, final_weight.GetNbinsX() + 1):
        final_weight.SetBinContent(bin_index, scale_factor)
        final_weight.SetBinError(bin_index, 0.0)
    final_control = _clone_parent_diagnostic_histogram(reference, scope, "single_scale_control")
    final_model = _clone_parent_diagnostic_histogram(templates["h_mm_full"], scope, "single_scale_model")
    final_weighted_model = _clone_parent_diagnostic_histogram(
        templates["h_mm_full"], scope, "single_scale_weighted_model"
    )
    expected_after = _clone_parent_diagnostic_histogram(
        full_before, scope, "single_scale_expected_after"
    )
    expected_after.Add(templates["h_mm_full"], -1.0)
    final = {
        "accepted": True,
        "proposal_available": False,
        "diagnostic_role": "final",
        "final_application_status": "applied_fallback",
        "production_evaluation_accepted": False,
        "production_rejection_reasons": [
            reason.strip()
            for reason in str(production_evaluation.get("reason") or "").split(";")
            if reason.strip()
        ],
        "fallback_used": True,
        "fallback_mode": "single_scale",
        "fallback_reason": production_evaluation.get("reason") or "component-fit result rejected",
        "diagnostic_only": True,
        "application_authoritative": False,
        "production_applied": False,
        "analysis_scope": proposal_payload.get("analysis_scope"),
        "input_selection": proposal_payload.get("input_selection"),
        "source_target_state": proposal_payload.get("source_target_state"),
        "H_pion_control_model": final_control,
        "H_kaon_pion_model": final_model,
        "H_weighted_pion_control_model": final_weighted_model,
        "H_pion_weight_vs_MM": final_weight,
        "weights": np.full(int(reference.GetNbinsX()) + 2, scale_factor, dtype=np.float64),
        "H_pion_subtraction_template_MM": templates["h_mm"],
        "H_pion_subtraction_template_MM_nosub": templates["h_mm_full"],
        "H_MM_before_pion_subtraction": _clone_parent_diagnostic_histogram(before, scope, "single_scale_before"),
        "H_MM_after_pion_subtraction": targets["h_mm"],
        "H_MM_nosub_before_pion_subtraction": _clone_parent_diagnostic_histogram(full_before, scope, "single_scale_full_before"),
        "H_MM_nosub_after_pion_subtraction": targets["h_mm_full"],
        "H_MM_nosub_after_pion_subtraction_model_stage": _clone_parent_diagnostic_histogram(targets["h_mm_full"], scope, "single_scale_stage_model"),
        "H_MM_nosub_after_pion_subtraction_model_final": _clone_parent_diagnostic_histogram(targets["h_mm_full"], scope, "single_scale_final_model"),
        "particle_subtraction_effective_scale": scale_factor,
        "weighted_pion_integral": _hist_integral(templates["h_mm_full"]),
        "weighted_pion_integral_cut": _hist_integral(templates["h_mm"]),
        "weighted_pion_integral_full": _hist_integral(templates["h_mm_full"]),
        "kaon_integral_before_pion_sub": _hist_integral(before),
        "kaon_integral_after_pion_sub": _hist_integral(targets["h_mm"]),
        "kaon_integral_before_pion_sub_full": _hist_integral(full_before),
        "kaon_integral_after_pion_sub_full": _hist_integral(targets["h_mm_full"]),
        "diagnostics": {
            **stats,
            "fallback_mode": "single_scale",
            "scale_components": scale_components,
            "event_template_closure": compute_hist_closure_metrics(final_model, templates["h_mm_full"]),
            "final_closure": compute_hist_closure_metrics(
                expected_after, targets["h_mm_full"]
            ),
        },
    }
    assert_component_subtraction_payload_ownership(final)
    return final


def resolve_parent_diagnostic_final_application(
    proposal_payload,
    production_evaluation,
    *,
    fallback_context=None,
):
    """Return the final diagnostic state without changing production policy.

    The caller supplies a fallback factory only for the one policy that needs
    event-level reconstruction (``single_scale``).  ``zero`` is represented
    here from detached proposal-before histograms; skip/error intentionally
    retain no fabricated final spectrum.
    """
    policy = {
        "fit_accepted": bool(production_evaluation.get("accepted")),
        "fallback_mode": str(production_evaluation.get("fallback_mode") or "error").strip().lower(),
        "reason": production_evaluation.get("reason") or None,
    }
    policy["action"] = (
        "component_weight" if policy["fit_accepted"] else
        "single_scale" if policy["fallback_mode"] == "single_scale" else
        "zero" if policy["fallback_mode"] == "zero" else
        "skip_bin" if policy["fallback_mode"] == "skip_bin" else "error"
    )
    policy["child_valid"] = policy["action"] not in ("skip_bin", "error")
    if not isinstance(proposal_payload, dict):
        return None, {
            "status": "unavailable",
            "final_status": "unavailable",
            "final_reason": "proposal_unavailable",
            "application_policy": policy,
        }
    if policy["action"] == "component_weight":
        final = dict(proposal_payload)
        final.update(
            {
                "diagnostic_role": "final",
                "final_application_status": "applied_component",
                "production_evaluation_accepted": True,
                "production_rejection_reasons": [],
                "production_applied": False,
            }
        )
        return final, {
            "status": "available",
            "final_status": "applied_component",
            "final_reason": None,
            "application_policy": policy,
        }

    fallback_mode = policy["fallback_mode"]
    scope = "pion_parent_{}".format(proposal_payload.get("analysis_scope") or "unknown")
    if fallback_mode == "zero":
        return _build_zero_parent_diagnostic_final(proposal_payload, production_evaluation, scope), {
            "status": "available",
            "final_status": "zero",
            "final_reason": production_evaluation.get("reason"),
            "application_policy": policy,
        }
    if fallback_mode == "single_scale" and callable(fallback_context):
        final = fallback_context()
        return final, {
            "status": "available",
            "final_status": "applied_fallback",
            "final_reason": production_evaluation.get("reason"),
            "application_policy": policy,
        }
    if fallback_mode == "single_scale":
        return None, {
            "status": "unavailable",
            "final_status": "unavailable",
            "final_reason": "single_scale_fallback_builder_missing",
            "application_policy": policy,
        }
    return None, {
        "status": "partial",
        "final_status": fallback_mode or "error",
        "final_reason": production_evaluation.get("reason") or "component-fit result rejected",
        "application_policy": policy,
    }


def _process_rand_sub_tree(
    tree,
    print_label,
    timer_label,
    tmin,
    tmax,
    nomm_fills,
    fills,
    particle_type,
    hole_contains,
    evaluate_event,
    mm_min,
    mm_max,
    progress_bar,
):
    # Keep downstream filling on the same no-particle-subtraction selection and
    # shifted-variable implementation used by the canonical prepass and the
    # proton-cleaning source preparation.
    from apply_cuts import evaluate_pre_particle_subtraction_event

    print(print_label)
    entries = tree.GetEntries()
    progress_time = 0.0
    loop_start = perf_counter()
    nohole_xy_fill = fills["nohole_xy"]
    nohole_x_mm_fill = fills["nohole_x_mm"]
    nohole_y_mm_fill = fills["nohole_y_mm"]

    for i, evt in enumerate(tree):
        progress_start = perf_counter()
        progress_bar(i, entries, bar_length=25)
        progress_time += perf_counter() - progress_start

        selection_state = evaluate_pre_particle_subtraction_event(
            evt, mm_min, mm_max, hole_contains=hole_contains
        )
        allcuts = bool(selection_state["allcuts"])
        nommcuts = bool(selection_state["nommcuts"])
        noholecuts = bool(selection_state["noholecuts"]) if particle_type == "kaon" else False
        adj_hsdelta = float(selection_state["adj_hsdelta"])

        if not (noholecuts or nommcuts or allcuts):
            continue

        adj_MM = float(selection_state["adj_mm"])
        adj_t = float(selection_state["adj_t"])

        if noholecuts and nohole_xy_fill is not None:
            nohole_xy_fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
            nohole_x_mm_fill(evt.P_hgcer_xAtCer, adj_MM, evt.P_hgcer_npeSum)
            nohole_y_mm_fill(evt.P_hgcer_yAtCer, adj_MM, evt.P_hgcer_npeSum)

        if nommcuts:
            for fill in nomm_fills:
                fill(adj_MM)

        if allcuts:
            _fill_rand_sub_allcuts(evt, adj_MM, adj_t, adj_hsdelta, fills)

    loop_elapsed = perf_counter() - loop_start
    _print_rand_timer(timer_label, loop_elapsed, entries)
    _print_rand_timer("{} progressBar".format(timer_label), progress_time, entries)
    _print_rand_timer("{} other".format(timer_label), max(loop_elapsed - progress_time, 0.0), entries)


def _resolve_prepass_random_window_count(inp_dict, phi_setting):
    """Read the same timing-table random-window count used by ``rand_sub``."""
    run_key = {
        "Right": "runNumRight",
        "Left": "runNumLeft",
        "Center": "runNumCenter",
    }.get(phi_setting)
    run_tokens = str(inp_dict.get(run_key, "")).split()
    if not run_tokens:
        raise RuntimeError("No run number configured for {} prepass".format(phi_setting))
    try:
        run_number = int(run_tokens[-1])
    except ValueError as exc:
        raise RuntimeError("Invalid run number for {} prepass".format(phi_setting)) from exc
    timing_path = os.path.join(UTILPATH, "DB", "PARAM", "Timing_Parameters.csv")
    matches = []
    with open(timing_path, "r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            fields = line.partition("#")[0].strip().split(",")
            if len(fields) < 6:
                continue
            try:
                if int(fields[0]) <= run_number <= int(fields[1]):
                    matches.append(float(fields[5]))
            except ValueError:
                continue
    if not matches:
        raise RuntimeError("No random-window timing entry for run {}".format(run_number))
    return float(matches[-1])


def build_pre_particle_subtraction_binning_payload(phi_setting, inp_dict, *, source_bundle=None):
    """Collect only canonical-binning records before particle subtraction.

    The payload is intentionally record based: its raw support is the number
    of selected records, while signed support uses the physical source
    coefficients.  It does not allocate final histograms or touch pion/proton
    fitting code.
    """
    from apply_cuts import (
        evaluate_pre_particle_subtraction_event,
        set_shift_context,
        set_val,
    )

    particle_type = inp_dict["ParticleType"]
    epsset = inp_dict["EPSSET"]
    set_val(inp_dict)
    set_shift_context(phi_setting=phi_setting, shift_mode=inp_dict.get("shift_mode", "raw"))

    owned_files = []
    if source_bundle is None:
        sys.path.append("normalize")
        from get_eff_charge import get_eff_charge
        from hgcer_hole import apply_HGCer_hole_cut

        # The charge/normalization resolver is the existing owner of the
        # prompt/random/dummy physical coefficients.
        get_eff_charge({"phi_setting": phi_setting}, inp_dict, all_data=False)
        data_path = "{}/{}_{}_{}.root".format(
            OUTPATH, phi_setting, particle_type, inp_dict["InDATAFilename"]
        )
        dummy_path = "{}/{}_{}_{}.root".format(
            OUTPATH, phi_setting, particle_type, inp_dict["InDUMMYFilename"]
        )
        if not (os.path.isfile(data_path) and os.path.isfile(dummy_path)):
            raise FileNotFoundError("prepass source tree file missing for {}".format(phi_setting))
        data_file = open_root_file(data_path)
        dummy_file = open_root_file(dummy_path)
        owned_files.extend((data_file, dummy_file))
        source_bundle = _open_kaon_proton_cleaning_tree_bundle(
            data_file,
            dummy_file,
            particle_type,
            inp_dict.get("POL"),
            epsset,
            phi_setting,
            float(inp_dict["normfac_data"]),
            float(inp_dict["normfac_dummy"]),
            _resolve_prepass_random_window_count(inp_dict, phi_setting),
        )
        hole_cut = (
            apply_HGCer_hole_cut(inp_dict["Q2"], inp_dict["W"], epsset, phi_setting)
            if particle_type in ("kaon", "pion")
            else None
        )
        hole_contains = hole_cut.IsInside if hole_cut is not None else None
    else:
        hole_contains = source_bundle.get("hole_contains")

    records = []
    source_stats = {}
    for source_label, source_spec in ((source_bundle or {}).get("sources") or {}).items():
        tree = (source_spec or {}).get("tree")
        coefficient = float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        selected = 0
        seen = 0
        if tree is not None:
            for entry_index, event in enumerate(tree):
                seen += 1
                state = evaluate_pre_particle_subtraction_event(
                    event,
                    float(inp_dict["mm_min"]),
                    float(inp_dict["mm_max"]),
                    hole_contains=hole_contains,
                )
                if not state["nommcuts"]:
                    continue
                try:
                    phi_value = float(event.ph_q)
                except (AttributeError, TypeError, ValueError):
                    continue
                if not (math.isfinite(state["adj_t"]) and math.isfinite(phi_value)):
                    continue
                records.append(
                    {
                        "source_label": str(source_label),
                        "entry_index": int(entry_index),
                        "adj_t": float(state["adj_t"]),
                        "adj_mm": float(state["adj_mm"]),
                        "phi_value": phi_value,
                        "physical_coefficient": coefficient,
                    }
                )
                selected += 1
        source_stats[str(source_label)] = {
            "tree_name": (source_spec or {}).get("tree_name"),
            "entries_seen": int(seen),
            "entries_selected": int(selected),
            "physical_coefficient": coefficient,
        }

    # ROOT owns the tree objects; do not close files supplied by a caller.
    for owned_file in owned_files:
        if owned_file is not None and hasattr(owned_file, "Close"):
            owned_file.Close()
    return {
        "phi_setting": phi_setting,
        "records": records,
        "t_values": [record["adj_t"] for record in records],
        "phi_values": [record["phi_value"] for record in records],
        "signed_weights": [record["physical_coefficient"] for record in records],
        "raw_event_support": int(len(records)),
        "source_stats": source_stats,
        "selection_provenance": {
            "stage": "pre_particle_subtraction",
            "selection": "shared_nommcuts_with_hgcer_hole_rejection",
            "shifted_t_getter": "apply_cuts.get_shifted_t",
            "source_bundle": "noRF_prompt_random_dummy",
        },
    }


LAMBDA_GATE_SUMMARY_FIELDS = (
    "kinematic", "epsilon", "phi_setting", "selected_timing_branch",
    "timing_fit_accepted", "setting_support", "lambda_validation_window_key",
    "lambda_window_low", "lambda_window_high", "lambda_raw_prompt_count",
    "lambda_raw_signed_yield", "lambda_raw_absolute_support",
    "lambda_signed_to_absolute_support_ratio", "lambda_proposed_proton_yield",
    "lambda_removed_fraction", "lambda_removed_fraction_limit",
    "lambda_support_valid", "lambda_support_reasons",
    "lambda_observational_warnings", "lambda_gate_status", "production_action",
    "proton_cleaning_committed", "closure_tolerance",
    "proposed_pre_rf_closure_difference", "proposed_pre_rf_closure_passed",
    "final_applied_pre_rf_closure_difference", "final_applied_closure_passed",
    "lookup_rows_checked", "lookup_commit_mismatch_count",
    "lookup_commit_integrity_passed", "lookup_commit_mismatch_categories",
)


def _build_lambda_gate_summary_row(
    lambda_gate, diagnostics, *, outfilename, epsset, phi_setting,
):
    """Flatten one setting-wide Lambda-gate result for CSV artifacts."""
    candidate = diagnostics.get("selected_timing_candidate") or {}
    return {
        "kinematic": outfilename,
        "epsilon": epsset,
        "phi_setting": phi_setting,
        "selected_timing_branch": (
            candidate.get("timing_branch") or diagnostics.get("timing_branch")
        ),
        "timing_fit_accepted": lambda_gate.get("timing_fit_accepted"),
        "setting_support": lambda_gate.get("setting_support_label"),
        "lambda_validation_window_key": lambda_gate.get("validation_window_key"),
        "lambda_window_low": lambda_gate.get("window_low"),
        "lambda_window_high": lambda_gate.get("window_high"),
        "lambda_raw_prompt_count": lambda_gate.get("raw_prompt_count"),
        "lambda_raw_signed_yield": lambda_gate.get("raw_signed_yield"),
        "lambda_raw_absolute_support": lambda_gate.get("raw_absolute_support"),
        "lambda_signed_to_absolute_support_ratio": lambda_gate.get(
            "raw_signed_to_absolute_support_ratio"
        ),
        "lambda_proposed_proton_yield": lambda_gate.get("proposed_proton_yield"),
        "lambda_removed_fraction": lambda_gate.get("proposed_removed_fraction"),
        "lambda_removed_fraction_limit": lambda_gate.get("maximum_removed_fraction"),
        "lambda_support_valid": lambda_gate.get("support_valid"),
        "lambda_support_reasons": json.dumps(
            lambda_gate.get("support_reasons") or [], separators=(",", ":")
        ),
        "lambda_observational_warnings": json.dumps(
            lambda_gate.get("observational_warnings") or [], separators=(",", ":")
        ),
        "lambda_gate_status": lambda_gate.get("status"),
        "production_action": lambda_gate.get("production_action"),
        "proton_cleaning_committed": lambda_gate.get("proton_cleaning_committed"),
        "closure_tolerance": lambda_gate.get("closure_tolerance"),
        "proposed_pre_rf_closure_difference": lambda_gate.get(
            "proposed_pre_rf_closure_difference"
        ),
        "proposed_pre_rf_closure_passed": lambda_gate.get(
            "proposed_pre_rf_closure_passed"
        ),
        "final_applied_pre_rf_closure_difference": lambda_gate.get(
            "final_applied_pre_rf_closure_difference"
        ),
        "final_applied_closure_passed": lambda_gate.get("final_applied_closure_passed"),
        "lookup_rows_checked": lambda_gate.get("lookup_rows_checked"),
        "lookup_commit_mismatch_count": lambda_gate.get("lookup_commit_mismatch_count"),
        "lookup_commit_integrity_passed": lambda_gate.get(
            "lookup_commit_integrity_passed"
        ),
        "lookup_commit_mismatch_categories": json.dumps(
            lambda_gate.get("lookup_commit_mismatch_categories") or {},
            sort_keys=True,
            separators=(",", ":"),
        ),
    }


def _write_csv_atomically(path, fieldnames, rows):
    """Write a compact CSV without exposing a partially-written artifact."""
    return write_csv_atomically(path, fieldnames, rows)


def _upsert_lambda_gate_regression_summary(outpath, row):
    """Maintain one deterministic cross-setting Lambda-gate review table."""
    path = os.path.join(outpath, "proton_cleaning_lambda_gate_regression_summary.csv")
    key_fields = ("kinematic", "epsilon", "phi_setting")
    return upsert_sorted_csv(path, row, LAMBDA_GATE_SUMMARY_FIELDS, key_fields)


def _write_timing_t_validation_artifacts(
    cleaning_result, *, outpath, particle_type, outfilename, epsset, phi_setting
):
    """Persist serializable, validation-only timing-t diagnostics."""
    if str((cleaning_result or {}).get("method") or "") != "timing_t_event_weight":
        return []
    diagnostics = (cleaning_result or {}).get("diagnostics") or {}
    aero_validation = (
        diagnostics.get("aerogel_vs_t_validation")
        or diagnostics.get("aerogel_validation")
        or {}
    )
    cross_stage = diagnostics.get("cross_stage_t_consistency") or []
    base = "{}_{}_{}_{}_timing_t".format(
        particle_type, outfilename, phi_setting, epsset
    )
    artifacts = []
    cross_json = os.path.join(outpath, "{}_cross_stage_t_consistency.json".format(base))
    with open(cross_json, "w", encoding="utf-8") as handle:
        json.dump(cross_stage, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
    artifacts.append(cross_json)
    cross_summary_json = os.path.join(outpath, "{}_cross_stage_t_consistency_summary.json".format(base))
    with open(cross_summary_json, "w", encoding="utf-8") as handle:
        json.dump(
            diagnostics.get("cross_stage_t_consistency_summary") or {},
            handle,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    artifacts.append(cross_summary_json)
    cross_csv = os.path.join(outpath, "{}_cross_stage_t_consistency.csv".format(base))
    with open(cross_csv, "w", newline="", encoding="utf-8") as handle:
        fields = (
            "event_signature", "prepass_t", "prepared_proton_cleaning_adj_t",
            "downstream_t", "maximum_absolute_difference", "consistent",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in cross_stage:
            writer.writerow({field: row.get(field) for field in fields})
    artifacts.append(cross_csv)
    timing_diagnostics = {
        "canonical_interval_pair": diagnostics.get("canonical_t_binning") or {},
        "candidate_production_support": diagnostics.get("timing_candidate_diagnostics") or [],
        "selected_candidate": diagnostics.get("selected_timing_candidate") or {},
        "setting_support": diagnostics.get("setting_support") or {},
        "delta_support": diagnostics.get("delta_support") or [],
        "applied_cell_state": diagnostics.get("applied_timing_t_cell_map") or [],
        "event_lookup_state": {
            "boundary_counts": diagnostics.get("t_lookup_boundary_counts") or {},
            "lookup_count": diagnostics.get("prepared_event_lookup_count", 0),
        },
        "frozen_lookup_stage_summary": diagnostics.get("frozen_lookup_stage_summary") or {},
        "closure_state": {
            "by_cell": diagnostics.get("event_weight_closure_by_cell") or [],
            "by_delta": diagnostics.get("event_weight_closure_by_delta") or [],
            "by_t": diagnostics.get("event_weight_closure_by_t") or [],
            "lookup_by_t_phi": diagnostics.get("event_weight_lookup_by_t_phi") or [],
        },
        "aerogel_vs_t_validation": aero_validation,
        "timing_t_mm_diagnostics": diagnostics.get("timing_t_mm_diagnostics") or {},
        "timing_t_summary": diagnostics.get("timing_t_summary") or {},
        "lambda_preservation_gate": diagnostics.get("lambda_preservation_gate") or {},
        "lambda_preservation_event_audit": diagnostics.get("lambda_preservation_event_audit") or {},
        "timing_t_diagnostic_integrity": diagnostics.get("timing_t_diagnostic_integrity") or {},
        "generic_hgcer_diagnostic_integrity": diagnostics.get("generic_hgcer_diagnostic_integrity") or {},
    }
    timing_diagnostics_json = os.path.join(outpath, "{}_candidate_cell_diagnostics.json".format(base))
    with open(timing_diagnostics_json, "w", encoding="utf-8") as handle:
        json.dump(timing_diagnostics, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
    artifacts.append(timing_diagnostics_json)
    for suffix, payload in (
        ("timing_t_summary", timing_diagnostics["timing_t_summary"]),
        ("timing_t_mm_diagnostics", timing_diagnostics["timing_t_mm_diagnostics"]),
        ("lambda_preservation_gate", timing_diagnostics["lambda_preservation_gate"]),
        ("lambda_preservation_event_audit", timing_diagnostics["lambda_preservation_event_audit"]),
        ("timing_t_diagnostic_integrity", timing_diagnostics["timing_t_diagnostic_integrity"]),
        ("generic_hgcer_diagnostic_integrity", timing_diagnostics["generic_hgcer_diagnostic_integrity"]),
    ):
        path = os.path.join(outpath, "{}_{}.json".format(base, suffix))
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        artifacts.append(path)
    lambda_gate = dict(timing_diagnostics["lambda_preservation_gate"] or {})
    if lambda_gate:
        lambda_csv = os.path.join(
            outpath, "{}_proton_cleaning_lambda_gate_summary.csv".format(base)
        )
        lambda_row = _build_lambda_gate_summary_row(
            lambda_gate,
            diagnostics,
            outfilename=outfilename,
            epsset=epsset,
            phi_setting=phi_setting,
        )
        _write_csv_atomically(lambda_csv, LAMBDA_GATE_SUMMARY_FIELDS, [lambda_row])
        artifacts.append(lambda_csv)
        artifacts.append(_upsert_lambda_gate_regression_summary(outpath, lambda_row))
    audit_rows = list((timing_diagnostics["lambda_preservation_event_audit"] or {}).get("rows") or [])
    if audit_rows:
        audit_csv = os.path.join(
            outpath, "{}_proton_cleaning_lambda_gate_event_audit.csv".format(base)
        )
        audit_fields = sorted({key for row in audit_rows for key in (row or {})})
        with open(audit_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=audit_fields)
            writer.writeheader()
            writer.writerows(audit_rows)
        artifacts.append(audit_csv)
    for suffix, rows in (
        ("candidate_support", timing_diagnostics["candidate_production_support"]),
        ("applied_cells", timing_diagnostics["applied_cell_state"]),
    ):
        csv_path = os.path.join(outpath, "{}_{}.csv".format(base, suffix))
        fields = sorted({key for row in rows for key in (row or {})})
        with open(csv_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    field: (
                        json.dumps(value, sort_keys=True, separators=(",", ":"))
                        if isinstance(value, (dict, list, tuple)) else value
                    )
                    for field, value in (row or {}).items()
                })
        artifacts.append(csv_path)
    mm_rows = list((timing_diagnostics["timing_t_mm_diagnostics"] or {}).get("per_t_bin_summary") or [])
    if mm_rows:
        mm_csv = os.path.join(outpath, "{}_proton_cleaning_tbin_mm_summary.csv".format(base))
        with open(mm_csv, "w", newline="", encoding="utf-8") as handle:
            fields = [
                "setting", "epsilon", "phi_setting", "t_index", "t_low", "t_high",
                "event_count", "raw_prompt_event_count", "signed_event_weight_sum",
                "absolute_event_weight_support", "raw_missing_mass_yield",
                "estimated_proton_missing_mass_yield", "cleaned_missing_mass_yield",
                "final_rf_cleaned_missing_mass_yield", "raw_minus_estimated_proton",
                "pre_rf_cleaning_closure_difference", "windows",
            ]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in mm_rows:
                writer.writerow({
                    "setting": outfilename,
                    "epsilon": epsset,
                    "phi_setting": phi_setting,
                    **{
                        field: (
                            json.dumps(value, sort_keys=True, separators=(",", ":"))
                            if isinstance(value, (dict, list, tuple)) else value
                        )
                        for field, value in (row or {}).items()
                        if field in fields
                    },
                })
        artifacts.append(mm_csv)
    if aero_validation:
        aero_json = os.path.join(outpath, "{}_t_aerogel_validation.json".format(base))
        with open(aero_json, "w", encoding="utf-8") as handle:
            json.dump(aero_validation, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        artifacts.append(aero_json)
        t_edges = list(aero_validation.get("t_edges") or [])
        aero_edges = list(aero_validation.get("aero_edges") or [])
        matrix_payload = dict(aero_validation.get("matrix_payload") or {})
        matrix_metrics = dict(matrix_payload.get("metrics") or {})
        matrix_validity = dict(matrix_payload.get("validity_masks") or {})

        def _matrix_cell(metric, t_index, aero_index):
            rows = matrix_metrics.get(metric) or []
            if t_index >= len(rows) or aero_index >= len(rows[t_index] or []):
                return None
            return rows[t_index][aero_index]

        def _matrix_valid(metric, t_index, aero_index):
            rows = matrix_validity.get(metric) or []
            if t_index >= len(rows) or aero_index >= len(rows[t_index] or []):
                return False
            return bool(rows[t_index][aero_index])
        def _csv_value(value):
            return (
                json.dumps(value, sort_keys=True, separators=(",", ":"))
                if isinstance(value, (dict, list, tuple)) else value
            )

        cells_csv = os.path.join(
            outpath, "{}_proton_cleaning_t_aerogel_cells.csv".format(base)
        )
        with open(cells_csv, "w", newline="", encoding="utf-8") as handle:
            fields = [
                "setting", "epsilon", "phi_setting", "t_index", "t_low", "t_high",
                "aero_index", "aero_low", "aero_high", "raw_prompt_count", "signed_yield",
                "absolute_support", "estimated_proton_yield", "cleaned_yield",
                "average_proton_probability", "average_proton_probability_valid",
                "low_mm_removed_yield", "low_mm_removed_fraction", "low_mm_removed_fraction_valid",
                "lambda_removed_yield", "lambda_removed_fraction", "lambda_removed_fraction_valid",
            ]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for t_index in range(max(len(t_edges) - 1, 0)):
                for aero_index in range(max(len(aero_edges) - 1, 0)):
                    writer.writerow({
                        "setting": outfilename,
                        "epsilon": epsset,
                        "phi_setting": phi_setting,
                        "t_index": t_index,
                        "t_low": t_edges[t_index] if t_index < len(t_edges) else None,
                        "t_high": t_edges[t_index + 1] if t_index + 1 < len(t_edges) else None,
                        "aero_index": aero_index,
                        "aero_low": aero_edges[aero_index] if aero_index < len(aero_edges) else None,
                        "aero_high": aero_edges[aero_index + 1] if aero_index + 1 < len(aero_edges) else None,
                        "raw_prompt_count": _matrix_cell("selected_prompt_count", t_index, aero_index),
                        "signed_yield": _matrix_cell("signed_physical_yield", t_index, aero_index),
                        "absolute_support": _matrix_cell("absolute_physical_support", t_index, aero_index),
                        "estimated_proton_yield": _matrix_cell("estimated_proton_yield", t_index, aero_index),
                        "cleaned_yield": _matrix_cell("cleaned_yield", t_index, aero_index),
                        "average_proton_probability": _csv_value(_matrix_cell("average_proton_probability", t_index, aero_index)),
                        "average_proton_probability_valid": _matrix_valid("average_proton_probability", t_index, aero_index),
                        "low_mm_removed_yield": _matrix_cell("low_mm_removed_yield", t_index, aero_index),
                        "low_mm_removed_fraction": _csv_value(_matrix_cell("low_mm_removed_fraction", t_index, aero_index)),
                        "low_mm_removed_fraction_valid": _matrix_valid("low_mm_removed_fraction", t_index, aero_index),
                        "lambda_removed_yield": _matrix_cell("lambda_removed_yield", t_index, aero_index),
                        "lambda_removed_fraction": _csv_value(_matrix_cell("lambda_removed_fraction", t_index, aero_index)),
                        "lambda_removed_fraction_valid": _matrix_valid("lambda_removed_fraction", t_index, aero_index),
                    })
        artifacts.append(cells_csv)

        tbin_rows = list(aero_validation.get("per_t_bin_summary") or [])
        tbin_csv = os.path.join(
            outpath, "{}_proton_cleaning_tbin_pid_summary.csv".format(base)
        )
        with open(tbin_csv, "w", newline="", encoding="utf-8") as handle:
            fields = ["setting", "epsilon", "phi_setting"] + sorted(
                {key for row in tbin_rows for key in (row or {})}
            )
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in tbin_rows:
                writer.writerow({
                    "setting": outfilename, "epsilon": epsset, "phi_setting": phi_setting,
                    **{key: _csv_value(value) for key, value in (row or {}).items()},
                })
        artifacts.append(tbin_csv)

        warning_rows = []
        for row in list(aero_validation.get("warnings_by_t_bin") or []):
            labels = list((row or {}).get("warnings") or [])
            if not labels:
                continue
            for label in labels:
                warning_rows.append({
                    "t_index": row.get("t_index"), "t_low": row.get("t_low"),
                    "t_high": row.get("t_high"), "warning": str(label),
                    "setting": outfilename, "epsilon": epsset, "phi_setting": phi_setting,
                    "low_aero_average_weight": row.get("low_aero_average_weight"),
                    "high_aero_average_weight": row.get("high_aero_average_weight"),
                    "low_aero_average_weight_valid": row.get("low_aero_average_weight_valid"),
                    "high_aero_average_weight_valid": row.get("high_aero_average_weight_valid"),
                    "low_aero_lambda_removed_fraction": row.get("low_aero_lambda_removed_fraction"),
                    "high_aero_lambda_removed_fraction": row.get("high_aero_lambda_removed_fraction"),
                    "low_aero_lambda_removed_fraction_valid": row.get("low_aero_lambda_removed_fraction_valid"),
                    "high_aero_lambda_removed_fraction_valid": row.get("high_aero_lambda_removed_fraction_valid"),
                })
        warning_csv = os.path.join(
            outpath, "{}_proton_cleaning_tbin_pid_warnings.csv".format(base)
        )
        with open(warning_csv, "w", newline="", encoding="utf-8") as handle:
            fields = [
                "setting", "epsilon", "phi_setting", "t_index", "t_low", "t_high", "warning", "low_aero_average_weight",
                "high_aero_average_weight", "low_aero_average_weight_valid", "high_aero_average_weight_valid",
                "low_aero_lambda_removed_fraction", "high_aero_lambda_removed_fraction",
                "low_aero_lambda_removed_fraction_valid", "high_aero_lambda_removed_fraction_valid",
            ]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(warning_rows)
        artifacts.append(warning_csv)
    return list(dict.fromkeys(artifacts))


def rand_sub(
    phi_setting,
    inpDict,
    shift_mode="raw",
    emit_plots=True,
    component_payload=None,
    kaon_signal_shape_payload=None,
    kaon_sigma0_shape_payload=None,
):
    total_start = perf_counter()
    setup_start = perf_counter()

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = float(inpDict["tmin"])
    tmax = float(inpDict["tmax"])
    mm_min = float(inpDict["mm_min"])
    mm_max = float(inpDict["mm_max"])
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]
    data_charge_right = inpDict["data_charge_right"] 
    data_charge_left = inpDict["data_charge_left"] 
    data_charge_center = inpDict["data_charge_center"] 
    dummy_charge_right = inpDict["dummy_charge_right"] 
    dummy_charge_left = inpDict["dummy_charge_left"] 
    dummy_charge_center = inpDict["dummy_charge_center"]
    data_charge_err_right = inpDict["data_charge_err_right"] 
    data_charge_err_left = inpDict["data_charge_err_left"] 
    data_charge_err_center = inpDict["data_charge_err_center"] 
    dummy_charge_err_right = inpDict["dummy_charge_err_right"] 
    dummy_charge_err_left = inpDict["dummy_charge_err_left"] 
    dummy_charge_err_center = inpDict["dummy_charge_err_center"]     
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"]
    InData_error_efficiency_right = inpDict["InData_error_efficiency_right"] 
    InData_error_efficiency_left = inpDict["InData_error_efficiency_left"] 
    InData_error_efficiency_center = inpDict["InData_error_efficiency_center"]    
    efficiency_table = inpDict["efficiency_table"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]
    mm_plot_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    mm_plot_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    mm_plot_nbins = int(inpDict.get("bg_opt_mm_plot_nbins", BG_OPT_MM_PLOT_NBINS))

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    histDict["phi_setting"] = phi_setting  

    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import (
        apply_data_cuts,
        apply_data_sub_cuts,
        evaluate_data_cut_bools,
        evaluate_data_event,
        get_shift_mode,
        get_kaon_data_coordinate,
        get_shifted_mm,
        get_shifted_t,
        set_shift_context,
        set_val,
    )
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=shift_mode)
    coordinate_contract = get_kaon_data_coordinate(
        required=(
            ParticleType == "kaon"
            and resolve_particle_subtraction_mode(inpDict) == "simc_shape_components"
            and resolve_pion_subtraction_scope(inpDict) == "t_bin"
        )
    )
    if coordinate_contract is not None:
        histDict["kaon_data_coordinate"] = dict(coordinate_contract)
        histDict["kaon_data_coordinate_fingerprint"] = coordinate_contract[
            "coordinate_fingerprint"
        ]
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        histDict.update({ "phi_setting" : phi_setting})
        return histDict

    InFile_DATA = open_root_file(rootFileData)

    prompt_tree_name = get_prompt_tree_name(
        ParticleType, EPSSET, rf_state="noRF"
    )
    rand_tree_name = get_rand_tree_name(
        ParticleType, EPSSET, rf_state="noRF"
    )

    TBRANCH_DATA  = InFile_DATA.Get(prompt_tree_name)

    TBRANCH_RAND  = InFile_DATA.Get(rand_tree_name)

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        return histDict

    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    hgcer_cutg = None
    if ParticleType in ("kaon", "pion"):
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, phi_setting)

    ################################################################################################################################################

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get(prompt_tree_name)
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get(rand_tree_name)

    ##############
    # HARD CODED #
    ##############
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1.0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

    if ParticleType == "kaon":
        for c0, p in zip(c0_list, h_momentum_list):
            if p == 0.889:
                c0_dict["Q2p1W2p95_lowe"] = c0 # Proper value 0.888
            elif p == 0.968:
                c0_dict["Q0p5W2p40_lowe"] = c0
                c0_dict["Q3p0W3p14_lowe"] = c0 # Proper value 1.821
                c0_dict["Q5p5W3p02_lowe"] = c0 # Proper value 0.962
            elif p == 2.185:
                c0_dict["Q0p5W2p40_highe"] = c0 # Proper value 2.066
                c0_dict["Q3p0W2p32_lowe"] = c0
            elif p == 2.328:
                c0_dict["Q4p4W2p74_lowe"] = c0
            elif p == 3.266:
                c0_dict["Q5p5W3p02_highe"] = c0            
            elif p == 4.2:
                c0_dict["Q3p0W3p14_highe"] = c0 # Proper value 4.204
            elif p == 4.712:
                c0_dict["Q4p4W2p74_highe"] = c0            
            elif p == 5.292:
                c0_dict["Q2p1W2p95_highe"] = c0
            elif p == 6.59:
                c0_dict["Q3p0W2p32_highe"] = c0
    else:
        c0_dict["Q0p4W2p20_lowe"] = 0.0
        c0_dict["Q0p4W2p20_highe"] = 0.0
        
    ##############
    ##############        
    ##############
    
    ################################################################################################################################################
    # Grabs PID cut string

    if phi_setting == "Right":
        runNums= runNumRight
        for run in runNumRight.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue

        InData_efficiency = InData_efficiency_right
    if phi_setting == "Left":
        runNums= runNumLeft
        for run in runNumLeft.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = InData_efficiency_left
    if phi_setting == "Center":
        runNums= runNumCenter
        for run in runNumCenter.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = InData_efficiency_center

    if 'pid_text' in locals():
        print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
    else:
        print("ERROR: Invalid {} log file {}!".format(phi_setting.lower(),pid_log))
        pid_text = "\nNo {} cuts file found in logs...".format(phi_setting.lower())

    ###############################################################################################################################################
    # Grab windows for random subtraction

    # Section for grabing Prompt/Random selection parameters from PARAM file
    PARAMPATH = "%s/DB/PARAM" % UTILPATH
    print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, LTANAPATH))
    TimingCutFile = "%s/Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!
    TimingCutf = open(TimingCutFile)
    try:
        TimingCutFile
    except NameError:
        print("!!!!! ERRROR !!!!!\n One (or more) of the cut files not found!\n!!!!! ERRORR !!!!!")
        sys.exit(2)
    print("Reading timing cuts from %s" % TimingCutFile)
    PromptWindow = [0, 0]
    RandomWindows = [0, 0, 0, 0]
    linenum = 0 # Count line number we're on
    TempPar = -1 # To check later
    for line in TimingCutf: # Read all lines in the cut file
        linenum += 1 # Add one to line number at start of loop
        if(linenum > 1): # Skip first line
            line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
            line = line.rstrip()
            array = line.split(",") # Convert line into an array, anything after a comma is a new entry 
            if(int(run) in range (int(array[0]), int(array[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
                TempPar += 2 # If run number is in range, set to non -1 value
                BunchSpacing = float(array[2])
                CoinOffset = float(array[3]) # Coin offset value
                nSkip = float(array[4]) # Number of random windows skipped 
                nWindows = float(array[5]) # Total number of random windows
                PromptPeak = float(array[6]) # Pion CT prompt peak positon 
    TimingCutf.close() # After scanning all lines in file, close file

    if(TempPar == -1): # If value is still -1, run number provided din't match any ranges specified so exit 
        print("!!!!! ERROR !!!!!\n Run number specified does not fall within a set of runs for which cuts are defined in %s\n!!!!! ERROR !!!!!" % TimingCutFile)
        sys.exit(3)
    elif(TempPar > 1):
        print("!!! WARNING!!! Run number was found within the range of two (or more) line entries of %s !!! WARNING !!!" % TimingCutFile)
        print("The last matching entry will be treated as the input, you should ensure this is what you want")

    # From our values from the file, reconstruct our windows 
    PromptWindow[0] = PromptPeak - (BunchSpacing/2) - CoinOffset
    PromptWindow[1] = PromptPeak + (BunchSpacing/2) + CoinOffset
    RandomWindows[0] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing) - ((nWindows/2)*BunchSpacing)
    RandomWindows[1] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing)
    RandomWindows[2] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing)
    RandomWindows[3] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) + ((nWindows/2)*BunchSpacing)

    sys.path.append("normalize")
    from get_eff_charge import get_eff_charge

    # Upate hist dictionary with effective charge
    get_eff_charge(histDict, inpDict, all_data=False)

    norm_factor_data = inpDict["normfac_data"]
    norm_factor_dummy = inpDict["normfac_dummy"]  

    ################################################################################################################################################
    # Plot definitions

    H_hsdelta_DATA  = TH1D("H_hsdelta_DATA","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DATA  = TH1D("H_hsxptar_DATA","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DATA  = TH1D("H_hsyptar_DATA","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DATA    = TH1D("H_ssxfp_DATA","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DATA    = TH1D("H_ssyfp_DATA","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DATA   = TH1D("H_ssxpfp_DATA","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DATA   = TH1D("H_ssypfp_DATA","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DATA    = TH1D("H_hsxfp_DATA","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DATA    = TH1D("H_hsyfp_DATA","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DATA   = TH1D("H_hsxpfp_DATA","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DATA   = TH1D("H_hsypfp_DATA","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DATA  = TH1D("H_ssdelta_DATA","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DATA  = TH1D("H_ssxptar_DATA","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DATA  = TH1D("H_ssyptar_DATA","SHMS yptar", 100, -0.04, 0.04)
    H_q_DATA        = TH1D("H_q_DATA","q", 100, 0.0, 10.0)
    H_Q2_DATA       = TH1D("H_Q2_DATA","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DATA  = TH1D("H_W_DATA","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DATA       = TH1D("H_t_DATA","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DATA  = TH1D("H_epsilon_DATA","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DATA  = TH1D("H_MM_DATA",f"MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_rand_dummy_DATA  = TH1D("H_MM_rand_dummy_DATA",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_dummy_DATA  = TH1D("H_MM_dummy_DATA",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_full_DATA  = TH1D("H_MM_full_DATA",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit2sub_DATA  = TH1D("H_MM_fit2sub_DATA",f"MM_fit2sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit1sub_DATA  = TH1D("H_MM_fit1sub_DATA",f"MM_fit1sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_pisub_DATA  = TH1D("H_MM_pisub_DATA",f"MM_pisub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_nosub_DATA  = TH1D("H_MM_nosub_DATA",f"MM_nosub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_th_DATA  = TH1D("H_th_DATA","X' tar", 100, -0.1, 0.1)
    H_ph_DATA  = TH1D("H_ph_DATA","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DATA  = TH1D("H_ph_q_DATA","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DATA  = TH1D("H_th_q_DATA","Theta Detected (th_xq)", 100, -math.pi, math.pi)
    H_ph_recoil_DATA  = TH1D("H_ph_recoil_DATA","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
    H_th_recoil_DATA  = TH1D("H_th_recoil_DATA","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
    H_pmiss_DATA  = TH1D("H_pmiss_DATA","pmiss", 100, 0.0, 2.0)
    H_emiss_DATA  = TH1D("H_emiss_DATA","emiss", 100, 0.0, 2.0)
    H_pmx_DATA  = TH1D("H_pmx_DATA","pmx", 100, -10.0, 10.0)
    H_pmy_DATA  = TH1D("H_pmy_DATA","pmy ", 100, -10.0, 10.0)
    H_pmz_DATA  = TH1D("H_pmz_DATA","pmz", 100, -10.0, 10.0)
    H_ct_DATA = TH1D("H_ct_DATA", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)
    H_cal_etottracknorm_DATA = TH1D("H_cal_etottracknorm_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
    H_cer_npeSum_DATA = TH1D("H_cer_npeSum_DATA", "HMS Cer Npe Sum", 100, 0, 30)
    P_cal_etottracknorm_DATA = TH1D("P_cal_etottracknorm_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
    P_hgcer_npeSum_DATA = TH1D("P_hgcer_npeSum_DATA", "SHMS HGCer Npe Sum", 100, 0, 10)
    P_aero_npeSum_DATA = TH1D("P_aero_npeSum_DATA", "SHMS Aero Npe Sum", 100, 0, 30)

    H_hsdelta_DUMMY  = TH1D("H_hsdelta_DUMMY","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DUMMY  = TH1D("H_hsxptar_DUMMY","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DUMMY  = TH1D("H_hsyptar_DUMMY","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DUMMY    = TH1D("H_ssxfp_DUMMY","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DUMMY    = TH1D("H_ssyfp_DUMMY","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DUMMY   = TH1D("H_ssxpfp_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DUMMY   = TH1D("H_ssypfp_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DUMMY    = TH1D("H_hsxfp_DUMMY","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DUMMY    = TH1D("H_hsyfp_DUMMY","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DUMMY   = TH1D("H_hsxpfp_DUMMY","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DUMMY   = TH1D("H_hsypfp_DUMMY","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DUMMY  = TH1D("H_ssdelta_DUMMY","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DUMMY  = TH1D("H_ssxptar_DUMMY","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DUMMY  = TH1D("H_ssyptar_DUMMY","SHMS yptar", 100, -0.04, 0.04)
    H_q_DUMMY        = TH1D("H_q_DUMMY","q", 100, 0.0, 10.0)
    H_Q2_DUMMY       = TH1D("H_Q2_DUMMY","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY  = TH1D("H_W_DUMMY","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY       = TH1D("H_t_DUMMY","-t", 100, inpDict["tmin"], inpDict["tmax"])  
    H_epsilon_DUMMY  = TH1D("H_epsilon_DUMMY","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY  = TH1D("H_MM_DUMMY",f"MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_full_DUMMY  = TH1D("H_MM_full_DUMMY",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit2sub_DUMMY  = TH1D("H_MM_fit2sub_DUMMY",f"MM_fit2sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit1sub_DUMMY  = TH1D("H_MM_fit1sub_DUMMY",f"MM_fit1sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_pisub_DUMMY  = TH1D("H_MM_pisub_DUMMY",f"MM_pisub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_nosub_DUMMY  = TH1D("H_MM_nosub_DUMMY",f"MM_nosub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_th_DUMMY  = TH1D("H_th_DUMMY","X' tar", 100, -0.1, 0.1)
    H_ph_DUMMY  = TH1D("H_ph_DUMMY","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DUMMY  = TH1D("H_ph_q_DUMMY","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DUMMY  = TH1D("H_th_q_DUMMY","Theta Detected (th_xq)", 100, -math.pi, math.pi)
    H_ph_recoil_DUMMY  = TH1D("H_ph_recoil_DUMMY","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
    H_th_recoil_DUMMY  = TH1D("H_th_recoil_DUMMY","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
    H_pmiss_DUMMY  = TH1D("H_pmiss_DUMMY","pmiss", 100, 0.0, 2.0)
    H_emiss_DUMMY  = TH1D("H_emiss_DUMMY","emiss", 100, 0.0, 2.0)
    H_pmx_DUMMY  = TH1D("H_pmx_DUMMY","pmx", 100, -10.0, 10.0)
    H_pmy_DUMMY  = TH1D("H_pmy_DUMMY","pmy ", 100, -10.0, 10.0)
    H_pmz_DUMMY  = TH1D("H_pmz_DUMMY","pmz", 100, -10.0, 10.0)
    H_ct_DUMMY = TH1D("H_ct_DUMMY", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)

    H_hsdelta_RAND  = TH1D("H_hsdelta_RAND","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_RAND  = TH1D("H_hsxptar_RAND","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_RAND  = TH1D("H_hsyptar_RAND","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_RAND    = TH1D("H_ssxfp_RAND","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_RAND    = TH1D("H_ssyfp_RAND","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_RAND   = TH1D("H_ssxpfp_RAND","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_RAND   = TH1D("H_ssypfp_RAND","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_RAND    = TH1D("H_hsxfp_RAND","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_RAND    = TH1D("H_hsyfp_RAND","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_RAND   = TH1D("H_hsxpfp_RAND","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_RAND   = TH1D("H_hsypfp_RAND","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_RAND  = TH1D("H_ssdelta_RAND","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_RAND  = TH1D("H_ssxptar_RAND","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_RAND  = TH1D("H_ssyptar_RAND","SHMS yptar", 100, -0.04, 0.04)
    H_q_RAND        = TH1D("H_q_RAND","q", 100, 0.0, 10.0)
    H_Q2_RAND       = TH1D("H_Q2_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_RAND  = TH1D("H_W_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_RAND       = TH1D("H_t_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_RAND  = TH1D("H_epsilon_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_RAND  = TH1D("H_MM_RAND",f"MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_rand_dummy_RAND  = TH1D("H_MM_rand_dummy_RAND",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_dummy_RAND  = TH1D("H_MM_dummy_RAND",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_full_RAND  = TH1D("H_MM_full_RAND",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit2sub_RAND  = TH1D("H_MM_fit2sub_RAND",f"MM_fit2sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit1sub_RAND  = TH1D("H_MM_fit1sub_RAND",f"MM_fit1sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_pisub_RAND  = TH1D("H_MM_pisub_RAND",f"MM_pisub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_nosub_RAND  = TH1D("H_MM_nosub_RAND",f"MM_nosub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_th_RAND  = TH1D("H_th_RAND","X' tar", 100, -0.1, 0.1)
    H_ph_RAND  = TH1D("H_ph_RAND","Y' tar", 100, -0.1, 0.1)
    H_ph_q_RAND  = TH1D("H_ph_q_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_RAND  = TH1D("H_th_q_RAND","Theta Detected (th_xq)", 100, -math.pi, math.pi)
    H_ph_recoil_RAND  = TH1D("H_ph_recoil_RAND","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
    H_th_recoil_RAND  = TH1D("H_th_recoil_RAND","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
    H_pmiss_RAND  = TH1D("H_pmiss_RAND","pmiss", 100, 0.0, 2.0)
    H_emiss_RAND  = TH1D("H_emiss_RAND","emiss", 100, 0.0, 2.0)
    H_pmx_RAND  = TH1D("H_pmx_RAND","pmx", 100, -10.0, 10.0)
    H_pmy_RAND  = TH1D("H_pmy_RAND","pmy ", 100, -10.0, 10.0)
    H_pmz_RAND  = TH1D("H_pmz_RAND","pmz", 100, -10.0, 10.0)
    H_ct_RAND = TH1D("H_ct_RAND", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)

    H_hsdelta_DUMMY_RAND  = TH1D("H_hsdelta_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DUMMY_RAND  = TH1D("H_hsxptar_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DUMMY_RAND  = TH1D("H_hsyptar_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DUMMY_RAND    = TH1D("H_ssxfp_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DUMMY_RAND    = TH1D("H_ssyfp_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DUMMY_RAND   = TH1D("H_ssxpfp_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DUMMY_RAND   = TH1D("H_ssypfp_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DUMMY_RAND    = TH1D("H_hsxfp_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DUMMY_RAND    = TH1D("H_hsyfp_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DUMMY_RAND   = TH1D("H_hsxpfp_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DUMMY_RAND   = TH1D("H_hsypfp_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DUMMY_RAND  = TH1D("H_ssdelta_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DUMMY_RAND  = TH1D("H_ssxptar_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DUMMY_RAND  = TH1D("H_ssyptar_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
    H_q_DUMMY_RAND        = TH1D("H_q_DUMMY_RAND","q", 100, 0.0, 10.0)
    H_Q2_DUMMY_RAND       = TH1D("H_Q2_DUMMY_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY_RAND  = TH1D("H_W_DUMMY_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY_RAND       = TH1D("H_t_DUMMY_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DUMMY_RAND  = TH1D("H_epsilon_DUMMY_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY_RAND  = TH1D("H_MM_DUMMY_RAND",f"MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_full_DUMMY_RAND  = TH1D("H_MM_full_DUMMY_RAND",f"MM_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit2sub_DUMMY_RAND  = TH1D("H_MM_fit2sub_DUMMY_RAND",f"MM_fit2sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_fit1sub_DUMMY_RAND  = TH1D("H_MM_fit1sub_DUMMY_RAND",f"MM_fit1sub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_pisub_DUMMY_RAND  = TH1D("H_MM_pisub_DUMMY_RAND",f"MM_pisub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_MM_nosub_DUMMY_RAND  = TH1D("H_MM_nosub_DUMMY_RAND",f"MM_nosub_{ParticleType[0].upper()}", mm_plot_nbins, mm_plot_min, mm_plot_max)
    H_th_DUMMY_RAND  = TH1D("H_th_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
    H_ph_DUMMY_RAND  = TH1D("H_ph_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DUMMY_RAND  = TH1D("H_ph_q_DUMMY_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DUMMY_RAND  = TH1D("H_th_q_DUMMY_RAND","Theta Detected (th_xq)", 100, -math.pi, math.pi)
    H_ph_recoil_DUMMY_RAND  = TH1D("H_ph_recoil_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
    H_th_recoil_DUMMY_RAND  = TH1D("H_th_recoil_DUMMY_RAND","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
    H_pmiss_DUMMY_RAND  = TH1D("H_pmiss_DUMMY_RAND","pmiss", 100, 0.0, 2.0)
    H_emiss_DUMMY_RAND  = TH1D("H_emiss_DUMMY_RAND","emiss", 100, 0.0, 2.0)
    H_pmx_DUMMY_RAND  = TH1D("H_pmx_DUMMY_RAND","pmx", 100, -10.0, 10.0)
    H_pmy_DUMMY_RAND  = TH1D("H_pmy_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
    H_pmz_DUMMY_RAND  = TH1D("H_pmz_DUMMY_RAND","pmz", 100, -10.0, 10.0)
    H_ct_DUMMY_RAND = TH1D("H_ct_DUMMY_RAND", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)

    ################################################################################################################################################
    # 2D histograms

    MM_vs_CoinTime_DATA = TH2D("MM_vs_CoinTime_DATA","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DATA = TH2D("CoinTime_vs_beta_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DATA = TH2D("MM_vs_beta_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DATA = TH2D("MM_vs_H_cer_DATA", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_DATA = TH2D("MM_vs_H_cal_DATA", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_DATA = TH2D("MM_vs_P_cal_DATA", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
    MM_vs_P_hgcer_DATA = TH2D("MM_vs_P_hgcer_DATA", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DATA = TH2D("MM_vs_P_aero_DATA", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    phiq_vs_t_DATA = TH2D("phiq_vs_t_DATA","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DATA = TGraphPolar()
    polar_phiq_vs_t_DATA.SetName("polar_phiq_vs_t_DATA")
    Q2_vs_W_DATA = TH2D("Q2_vs_W_DATA", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DATA = TH2D("Q2_vs_t_DATA", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DATA = TH2D("W_vs_t_DATA", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DATA = TH2D("EPS_vs_t_DATA", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DATA = TH2D("MM_vs_t_DATA", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DATA = TH2D("P_hgcer_xAtCer_vs_yAtCer_DATA", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DATA", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DATA = TH2D("P_hgcer_xAtCer_vs_MM_DATA", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DATA", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DATA = TH2D("P_hgcer_yAtCer_vs_MM_DATA", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DATA = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DATA", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_DUMMY = TH2D("MM_vs_CoinTime_DUMMY","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DUMMY = TH2D("CoinTime_vs_beta_DUMMY", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DUMMY = TH2D("MM_vs_beta_DUMMY", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DUMMY = TH2D("MM_vs_H_cer_DUMMY", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_P_cal_DUMMY = TH2D("MM_vs_P_cal_DUMMY", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)    
    MM_vs_H_cal_DUMMY = TH2D("MM_vs_H_cal_DUMMY", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
    MM_vs_P_hgcer_DUMMY = TH2D("MM_vs_P_hgcer_DUMMY", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DUMMY = TH2D("MM_vs_P_aero_DUMMY", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_DUMMY = TH2D("phiq_vs_t_DUMMY","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DUMMY = TGraphPolar()
    polar_phiq_vs_t_DUMMY.SetName("polar_phiq_vs_t_DUMMY")
    Q2_vs_W_DUMMY = TH2D("Q2_vs_W_DUMMY", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DUMMY = TH2D("Q2_vs_t_DUMMY", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DUMMY = TH2D("W_vs_t_DUMMY", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DUMMY = TH2D("EPS_vs_t_DUMMY", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DUMMY = TH2D("MM_vs_t_DUMMY", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DUMMY = TH2D("P_hgcer_xAtCer_vs_yAtCer_DUMMY", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DUMMY = TH2D("P_hgcer_xAtCer_vs_MM_DUMMY", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DUMMY", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DUMMY = TH2D("P_hgcer_yAtCer_vs_MM_DUMMY", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DUMMY", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_RAND = TH2D("MM_vs_CoinTime_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_RAND = TH2D("CoinTime_vs_beta_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_RAND = TH2D("MM_vs_beta_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_RAND = TH2D("MM_vs_H_cer_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_RAND = TH2D("MM_vs_H_cal_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_RAND = TH2D("MM_vs_P_cal_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)    
    MM_vs_P_hgcer_RAND = TH2D("MM_vs_P_hgcer_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_RAND = TH2D("MM_vs_P_aero_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_RAND = TH2D("phiq_vs_t_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_RAND = TH2D("Q2_vs_W_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_RAND = TH2D("Q2_vs_t_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_RAND = TH2D("W_vs_t_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_RAND = TH2D("EPS_vs_t_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_RAND = TH2D("MM_vs_t_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_RAND = TH2D("P_hgcer_xAtCer_vs_yAtCer_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_RAND = TH2D("P_hgcer_xAtCer_vs_MM_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_MM_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_RAND = TH2D("P_hgcer_yAtCer_vs_MM_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_RAND = TH2D("P_hgcer_nohole_yAtCer_vs_MM_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_DUMMY_RAND = TH2D("MM_vs_CoinTime_DUMMY_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DUMMY_RAND = TH2D("CoinTime_vs_beta_DUMMY_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DUMMY_RAND = TH2D("MM_vs_beta_DUMMY_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DUMMY_RAND = TH2D("MM_vs_H_cer_DUMMY_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_DUMMY_RAND = TH2D("MM_vs_H_cal_DUMMY_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_DUMMY_RAND = TH2D("MM_vs_P_cal_DUMMY_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)    
    MM_vs_P_hgcer_DUMMY_RAND = TH2D("MM_vs_P_hgcer_DUMMY_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DUMMY_RAND = TH2D("MM_vs_P_aero_DUMMY_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_DUMMY_RAND = TH2D("phiq_vs_t_DUMMY_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_DUMMY_RAND = TH2D("Q2_vs_W_DUMMY_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DUMMY_RAND = TH2D("Q2_vs_t_DUMMY_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DUMMY_RAND = TH2D("W_vs_t_DUMMY_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DUMMY_RAND = TH2D("EPS_vs_t_DUMMY_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DUMMY_RAND = TH2D("MM_vs_t_DUMMY_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND = TH2D("P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_xAtCer_vs_MM_DUMMY_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_yAtCer_vs_MM_DUMMY_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)        

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":        
        
        from particle_subtraction import particle_subtraction_cuts
        SubtractedParticle = "pion"
        subDict = {}
        
        subDict["H_hsdelta_SUB_DATA"]  = TH1D("H_hsdelta_SUB_DATA","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DATA"]  = TH1D("H_hsxptar_SUB_DATA","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DATA"]  = TH1D("H_hsyptar_SUB_DATA","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DATA"]    = TH1D("H_ssxfp_SUB_DATA","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DATA"]    = TH1D("H_ssyfp_SUB_DATA","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DATA"]   = TH1D("H_ssxpfp_SUB_DATA","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DATA"]   = TH1D("H_ssypfp_SUB_DATA","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DATA"]    = TH1D("H_hsxfp_SUB_DATA","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DATA"]    = TH1D("H_hsyfp_SUB_DATA","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DATA"]   = TH1D("H_hsxpfp_SUB_DATA","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DATA"]   = TH1D("H_hsypfp_SUB_DATA","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DATA"]  = TH1D("H_ssdelta_SUB_DATA","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DATA"]  = TH1D("H_ssxptar_SUB_DATA","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DATA"]  = TH1D("H_ssyptar_SUB_DATA","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DATA"]        = TH1D("H_q_SUB_DATA","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DATA"]       = TH1D("H_Q2_SUB_DATA","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DATA"]  = TH1D("H_W_SUB_DATA","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DATA"]       = TH1D("H_t_SUB_DATA","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DATA"]  = TH1D("H_epsilon_SUB_DATA","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DATA"]  = TH1D("H_MM_SUB_DATA",f"MM_{SubtractedParticle}", 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_full_SUB_DATA"]  = TH1D("H_MM_full_SUB_DATA",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_MM_nosub_SUB_DATA"]  = TH1D("H_MM_nosub_SUB_DATA",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_th_SUB_DATA"]  = TH1D("H_th_SUB_DATA","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DATA"]  = TH1D("H_ph_SUB_DATA","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DATA"]  = TH1D("H_ph_q_SUB_DATA","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DATA"]  = TH1D("H_th_q_SUB_DATA","Theta Detected (th_xq)", 100, -math.pi, math.pi)
        subDict["H_ph_recoil_SUB_DATA"]  = TH1D("H_ph_recoil_SUB_DATA","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
        subDict["H_th_recoil_SUB_DATA"]  = TH1D("H_th_recoil_SUB_DATA","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
        subDict["H_pmiss_SUB_DATA"]  = TH1D("H_pmiss_SUB_DATA","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DATA"]  = TH1D("H_emiss_SUB_DATA","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DATA"]  = TH1D("H_pmx_SUB_DATA","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DATA"]  = TH1D("H_pmy_SUB_DATA","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DATA"]  = TH1D("H_pmz_SUB_DATA","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DATA"] = TH1D("H_ct_SUB_DATA", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DATA"] = TH1D("H_cal_etottracknorm_SUB_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DATA"] = TH1D("H_cer_npeSum_SUB_DATA", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DATA"] = TH1D("P_cal_etottracknorm_SUB_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DATA"] = TH1D("P_hgcer_npeSum_SUB_DATA", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DATA"] = TH1D("P_aero_npeSum_SUB_DATA", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["H_hsdelta_SUB_RAND"]  = TH1D("H_hsdelta_SUB_RAND","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_RAND"]  = TH1D("H_hsxptar_SUB_RAND","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_RAND"]  = TH1D("H_hsyptar_SUB_RAND","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_RAND"]    = TH1D("H_ssxfp_SUB_RAND","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_RAND"]    = TH1D("H_ssyfp_SUB_RAND","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_RAND"]   = TH1D("H_ssxpfp_SUB_RAND","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_RAND"]   = TH1D("H_ssypfp_SUB_RAND","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_RAND"]    = TH1D("H_hsxfp_SUB_RAND","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_RAND"]    = TH1D("H_hsyfp_SUB_RAND","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_RAND"]   = TH1D("H_hsxpfp_SUB_RAND","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_RAND"]   = TH1D("H_hsypfp_SUB_RAND","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_RAND"]  = TH1D("H_ssdelta_SUB_RAND","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_RAND"]  = TH1D("H_ssxptar_SUB_RAND","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_RAND"]  = TH1D("H_ssyptar_SUB_RAND","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_RAND"]        = TH1D("H_q_SUB_RAND","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_RAND"]       = TH1D("H_Q2_SUB_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_RAND"]  = TH1D("H_W_SUB_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_RAND"]       = TH1D("H_t_SUB_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_RAND"]  = TH1D("H_epsilon_SUB_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_RAND"]  = TH1D("H_MM_SUB_RAND",f"MM_{SubtractedParticle}", 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_full_SUB_RAND"]  = TH1D("H_MM_full_SUB_RAND",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_MM_nosub_SUB_RAND"]  = TH1D("H_MM_nosub_SUB_RAND",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_th_SUB_RAND"]  = TH1D("H_th_SUB_RAND","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_RAND"]  = TH1D("H_ph_SUB_RAND","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_RAND"]  = TH1D("H_ph_q_SUB_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_RAND"]  = TH1D("H_th_q_SUB_RAND","Theta Detected (th_xq)", 100, -math.pi, math.pi)
        subDict["H_ph_recoil_SUB_RAND"]  = TH1D("H_ph_recoil_SUB_RAND","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
        subDict["H_th_recoil_SUB_RAND"]  = TH1D("H_th_recoil_SUB_RAND","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
        subDict["H_pmiss_SUB_RAND"]  = TH1D("H_pmiss_SUB_RAND","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_RAND"]  = TH1D("H_emiss_SUB_RAND","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_RAND"]  = TH1D("H_pmx_SUB_RAND","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_RAND"]  = TH1D("H_pmy_SUB_RAND","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_RAND"]  = TH1D("H_pmz_SUB_RAND","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_RAND"] = TH1D("H_ct_SUB_RAND", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_RAND"] = TH1D("H_cal_etottracknorm_SUB_RAND", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_RAND"] = TH1D("H_cer_npeSum_SUB_RAND", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_RAND"] = TH1D("P_cal_etottracknorm_SUB_RAND", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_RAND"] = TH1D("P_hgcer_npeSum_SUB_RAND", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_RAND"] = TH1D("P_aero_npeSum_SUB_RAND", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["H_hsdelta_SUB_DUMMY"]  = TH1D("H_hsdelta_SUB_DUMMY","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DUMMY"]  = TH1D("H_hsxptar_SUB_DUMMY","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DUMMY"]  = TH1D("H_hsyptar_SUB_DUMMY","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DUMMY"]    = TH1D("H_ssxfp_SUB_DUMMY","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DUMMY"]    = TH1D("H_ssyfp_SUB_DUMMY","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DUMMY"]   = TH1D("H_ssxpfp_SUB_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DUMMY"]   = TH1D("H_ssypfp_SUB_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DUMMY"]    = TH1D("H_hsxfp_SUB_DUMMY","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DUMMY"]    = TH1D("H_hsyfp_SUB_DUMMY","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DUMMY"]   = TH1D("H_hsxpfp_SUB_DUMMY","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DUMMY"]   = TH1D("H_hsypfp_SUB_DUMMY","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DUMMY"]  = TH1D("H_ssdelta_SUB_DUMMY","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DUMMY"]  = TH1D("H_ssxptar_SUB_DUMMY","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DUMMY"]  = TH1D("H_ssyptar_SUB_DUMMY","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DUMMY"]        = TH1D("H_q_SUB_DUMMY","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DUMMY"]       = TH1D("H_Q2_SUB_DUMMY","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DUMMY"]  = TH1D("H_W_SUB_DUMMY","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DUMMY"]       = TH1D("H_t_SUB_DUMMY","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DUMMY"]  = TH1D("H_epsilon_SUB_DUMMY","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DUMMY"]  = TH1D("H_MM_SUB_DUMMY",f"MM_{SubtractedParticle}", 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_full_SUB_DUMMY"]  = TH1D("H_MM_full_SUB_DUMMY",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_MM_nosub_SUB_DUMMY"]  = TH1D("H_MM_nosub_SUB_DUMMY",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_th_SUB_DUMMY"]  = TH1D("H_th_SUB_DUMMY","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DUMMY"]  = TH1D("H_ph_SUB_DUMMY","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DUMMY"]  = TH1D("H_ph_q_SUB_DUMMY","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DUMMY"]  = TH1D("H_th_q_SUB_DUMMY","Theta Detected (th_xq)", 100, -math.pi, math.pi)
        subDict["H_ph_recoil_SUB_DUMMY"]  = TH1D("H_ph_recoil_SUB_DUMMY","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
        subDict["H_th_recoil_SUB_DUMMY"]  = TH1D("H_th_recoil_SUB_DUMMY","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
        subDict["H_pmiss_SUB_DUMMY"]  = TH1D("H_pmiss_SUB_DUMMY","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DUMMY"]  = TH1D("H_emiss_SUB_DUMMY","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DUMMY"]  = TH1D("H_pmx_SUB_DUMMY","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DUMMY"]  = TH1D("H_pmy_SUB_DUMMY","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DUMMY"]  = TH1D("H_pmz_SUB_DUMMY","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DUMMY"] = TH1D("H_ct_SUB_DUMMY", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DUMMY"] = TH1D("H_cal_etottracknorm_SUB_DUMMY", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DUMMY"] = TH1D("H_cer_npeSum_SUB_DUMMY", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DUMMY"] = TH1D("P_cal_etottracknorm_SUB_DUMMY", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DUMMY"] = TH1D("P_hgcer_npeSum_SUB_DUMMY", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DUMMY"] = TH1D("P_aero_npeSum_SUB_DUMMY", "SHMS Aero Npe Sum", 100, 0, 30)        

        subDict["H_hsdelta_SUB_DUMMY_RAND"]  = TH1D("H_hsdelta_SUB_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DUMMY_RAND"]  = TH1D("H_hsxptar_SUB_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DUMMY_RAND"]  = TH1D("H_hsyptar_SUB_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DUMMY_RAND"]    = TH1D("H_ssxfp_SUB_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DUMMY_RAND"]    = TH1D("H_ssyfp_SUB_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DUMMY_RAND"]   = TH1D("H_ssxpfp_SUB_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DUMMY_RAND"]   = TH1D("H_ssypfp_SUB_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DUMMY_RAND"]    = TH1D("H_hsxfp_SUB_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DUMMY_RAND"]    = TH1D("H_hsyfp_SUB_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DUMMY_RAND"]   = TH1D("H_hsxpfp_SUB_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DUMMY_RAND"]   = TH1D("H_hsypfp_SUB_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DUMMY_RAND"]  = TH1D("H_ssdelta_SUB_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DUMMY_RAND"]  = TH1D("H_ssxptar_SUB_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DUMMY_RAND"]  = TH1D("H_ssyptar_SUB_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DUMMY_RAND"]        = TH1D("H_q_SUB_DUMMY_RAND","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DUMMY_RAND"]       = TH1D("H_Q2_SUB_DUMMY_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DUMMY_RAND"]  = TH1D("H_W_SUB_DUMMY_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DUMMY_RAND"]       = TH1D("H_t_SUB_DUMMY_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DUMMY_RAND"]  = TH1D("H_epsilon_SUB_DUMMY_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DUMMY_RAND"]  = TH1D("H_MM_SUB_DUMMY_RAND",f"MM_{SubtractedParticle}", 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_full_SUB_DUMMY_RAND"]  = TH1D("H_MM_full_SUB_DUMMY_RAND",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_MM_nosub_SUB_DUMMY_RAND"]  = TH1D("H_MM_nosub_SUB_DUMMY_RAND",f"MM_{SubtractedParticle}", mm_plot_nbins, mm_plot_min, mm_plot_max)
        subDict["H_th_SUB_DUMMY_RAND"]  = TH1D("H_th_SUB_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DUMMY_RAND"]  = TH1D("H_ph_SUB_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DUMMY_RAND"]  = TH1D("H_ph_q_SUB_DUMMY_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DUMMY_RAND"]  = TH1D("H_th_q_SUB_DUMMY_RAND","Theta Detected (th_xq)", 100, -math.pi, math.pi)
        subDict["H_ph_recoil_SUB_DUMMY_RAND"]  = TH1D("H_ph_recoil_SUB_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -math.pi, math.pi)
        subDict["H_th_recoil_SUB_DUMMY_RAND"]  = TH1D("H_th_recoil_SUB_DUMMY_RAND","Theta Recoil (th_bq)", 100, -math.pi, math.pi)
        subDict["H_pmiss_SUB_DUMMY_RAND"]  = TH1D("H_pmiss_SUB_DUMMY_RAND","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DUMMY_RAND"]  = TH1D("H_emiss_SUB_DUMMY_RAND","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DUMMY_RAND"]  = TH1D("H_pmx_SUB_DUMMY_RAND","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DUMMY_RAND"]  = TH1D("H_pmy_SUB_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DUMMY_RAND"]  = TH1D("H_pmz_SUB_DUMMY_RAND","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DUMMY_RAND"] = TH1D("H_ct_SUB_DUMMY_RAND", f"Electron-{ParticleType.capitalize()} CTime", 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DUMMY_RAND"] = TH1D("H_cal_etottracknorm_SUB_DUMMY_RAND", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DUMMY_RAND"] = TH1D("H_cer_npeSum_SUB_DUMMY_RAND", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DUMMY_RAND"] = TH1D("P_cal_etottracknorm_SUB_DUMMY_RAND", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DUMMY_RAND"] = TH1D("P_hgcer_npeSum_SUB_DUMMY_RAND", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DUMMY_RAND"] = TH1D("P_aero_npeSum_SUB_DUMMY_RAND", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["MM_vs_CoinTime_SUB_DATA"] = TH2D("MM_vs_CoinTime_SUB_DATA","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DATA"] = TH2D("CoinTime_vs_beta_SUB_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DATA"] = TH2D("MM_vs_beta_SUB_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DATA"] = TH2D("MM_vs_H_cer_SUB_DATA", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DATA"] = TH2D("MM_vs_H_cal_SUB_DATA", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DATA"] = TH2D("MM_vs_P_cal_SUB_DATA", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DATA"] = TH2D("MM_vs_P_hgcer_SUB_DATA", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DATA"] = TH2D("MM_vs_P_aero_SUB_DATA", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DATA"] = TH2D("phiq_vs_t_SUB_DATA","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DATA"] = TH2D("Q2_vs_W_SUB_DATA", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DATA"] = TH2D("Q2_vs_t_SUB_DATA", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DATA"] = TH2D("W_vs_t_SUB_DATA", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DATA"] = TH2D("EPS_vs_t_SUB_DATA", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DATA"] = TH2D("MM_vs_t_SUB_DATA", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DATA", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DATA", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DATA", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_DUMMY"] = TH2D("MM_vs_CoinTime_SUB_DUMMY","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DUMMY"] = TH2D("CoinTime_vs_beta_SUB_DUMMY", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DUMMY"] = TH2D("MM_vs_beta_SUB_DUMMY", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DUMMY"] = TH2D("MM_vs_H_cer_SUB_DUMMY", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DUMMY"] = TH2D("MM_vs_H_cal_SUB_DUMMY", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DUMMY"] = TH2D("MM_vs_P_cal_SUB_DUMMY", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DUMMY"] = TH2D("MM_vs_P_hgcer_SUB_DUMMY", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DUMMY"] = TH2D("MM_vs_P_aero_SUB_DUMMY", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DUMMY"] = TH2D("phiq_vs_t_SUB_DUMMY","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DUMMY"] = TH2D("Q2_vs_W_SUB_DUMMY", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DUMMY"] = TH2D("Q2_vs_t_SUB_DUMMY", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DUMMY"] = TH2D("W_vs_t_SUB_DUMMY", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DUMMY"] = TH2D("EPS_vs_t_SUB_DUMMY", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DUMMY"] = TH2D("MM_vs_t_SUB_DUMMY", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DUMMY", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DUMMY", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_RAND"] = TH2D("MM_vs_CoinTime_SUB_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_RAND"] = TH2D("CoinTime_vs_beta_SUB_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_RAND"] = TH2D("MM_vs_beta_SUB_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_RAND"] = TH2D("MM_vs_H_cer_SUB_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_RAND"] = TH2D("MM_vs_H_cal_SUB_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_RAND"] = TH2D("MM_vs_P_cal_SUB_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_RAND"] = TH2D("MM_vs_P_hgcer_SUB_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_RAND"] = TH2D("MM_vs_P_aero_SUB_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_RAND"] = TH2D("phiq_vs_t_SUB_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_RAND"] = TH2D("Q2_vs_W_SUB_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_RAND"] = TH2D("Q2_vs_t_SUB_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_RAND"] = TH2D("W_vs_t_SUB_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_RAND"] = TH2D("EPS_vs_t_SUB_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_RAND"] = TH2D("MM_vs_t_SUB_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_RAND"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_DUMMY_RAND"] = TH2D("MM_vs_CoinTime_SUB_DUMMY_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DUMMY_RAND"] = TH2D("CoinTime_vs_beta_SUB_DUMMY_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DUMMY_RAND"] = TH2D("MM_vs_beta_SUB_DUMMY_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DUMMY_RAND"] = TH2D("MM_vs_H_cer_SUB_DUMMY_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DUMMY_RAND"] = TH2D("MM_vs_H_cal_SUB_DUMMY_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_cal_SUB_DUMMY_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_hgcer_SUB_DUMMY_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_aero_SUB_DUMMY_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DUMMY_RAND"] = TH2D("phiq_vs_t_SUB_DUMMY_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DUMMY_RAND"] = TH2D("Q2_vs_W_SUB_DUMMY_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DUMMY_RAND"] = TH2D("Q2_vs_t_SUB_DUMMY_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DUMMY_RAND"] = TH2D("W_vs_t_SUB_DUMMY_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DUMMY_RAND"] = TH2D("EPS_vs_t_SUB_DUMMY_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DUMMY_RAND"] = TH2D("MM_vs_t_SUB_DUMMY_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

    # Fit background and subtract
    from background_fit import bg_fit
    _print_rand_timer("rand_sub setup {}".format(phi_setting), perf_counter() - setup_start)
        
    ################################################################################################################################################
    # Fill histograms for various trees called above

    hole_contains = hgcer_cutg.IsInside if hgcer_cutg is not None else None

    if ParticleType == "kaon":
        data_nohole_xy_fill = P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Fill
        data_nohole_x_mm_fill = P_hgcer_nohole_xAtCer_vs_MM_DATA.Fill
        data_nohole_y_mm_fill = P_hgcer_nohole_yAtCer_vs_MM_DATA.Fill
    else:
        data_nohole_xy_fill = None
        data_nohole_x_mm_fill = None
        data_nohole_y_mm_fill = None

    data_fills = {
        "nohole_xy": data_nohole_xy_fill,
        "nohole_x_mm": data_nohole_x_mm_fill,
        "nohole_y_mm": data_nohole_y_mm_fill,
        "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DATA.Fill,
        "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DATA.Fill,
        "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DATA.Fill,
        "mm_ct": MM_vs_CoinTime_DATA.Fill,
        "ct_beta": CoinTime_vs_beta_DATA.Fill,
        "mm_beta": MM_vs_beta_DATA.Fill,
        "mm_h_cer": MM_vs_H_cer_DATA.Fill,
        "mm_h_cal": MM_vs_H_cal_DATA.Fill,
        "mm_p_cal": MM_vs_P_cal_DATA.Fill,
        "mm_p_hgcer": MM_vs_P_hgcer_DATA.Fill,
        "mm_p_aero": MM_vs_P_aero_DATA.Fill,
        "phiq_t": phiq_vs_t_DATA.Fill,
        "q2_w": Q2_vs_W_DATA.Fill,
        "q2_t": Q2_vs_t_DATA.Fill,
        "w_t": W_vs_t_DATA.Fill,
        "eps_t": EPS_vs_t_DATA.Fill,
        "mm_t": MM_vs_t_DATA.Fill,
        "polar_graph": polar_phiq_vs_t_DATA,
        "h_ct": H_ct_DATA.Fill,
        "h_ssxfp": H_ssxfp_DATA.Fill,
        "h_ssyfp": H_ssyfp_DATA.Fill,
        "h_ssxpfp": H_ssxpfp_DATA.Fill,
        "h_ssypfp": H_ssypfp_DATA.Fill,
        "h_ssdelta": H_ssdelta_DATA.Fill,
        "h_ssxptar": H_ssxptar_DATA.Fill,
        "h_ssyptar": H_ssyptar_DATA.Fill,
        "h_hsxfp": H_hsxfp_DATA.Fill,
        "h_hsyfp": H_hsyfp_DATA.Fill,
        "h_hsxpfp": H_hsxpfp_DATA.Fill,
        "h_hsypfp": H_hsypfp_DATA.Fill,
        "h_hsdelta": H_hsdelta_DATA.Fill,
        "h_hsxptar": H_hsxptar_DATA.Fill,
        "h_hsyptar": H_hsyptar_DATA.Fill,
        "h_ph_q": H_ph_q_DATA.Fill,
        "h_th_q": H_th_q_DATA.Fill,
        "h_ph_recoil": H_ph_recoil_DATA.Fill,
        "h_th_recoil": H_th_recoil_DATA.Fill,
        "h_pmiss": H_pmiss_DATA.Fill,
        "h_emiss": H_emiss_DATA.Fill,
        "h_pmx": H_pmx_DATA.Fill,
        "h_pmy": H_pmy_DATA.Fill,
        "h_pmz": H_pmz_DATA.Fill,
        "h_q2": H_Q2_DATA.Fill,
        "h_t": H_t_DATA.Fill,
        "h_w": H_W_DATA.Fill,
        "h_epsilon": H_epsilon_DATA.Fill,
        "h_mm": H_MM_DATA.Fill,
        "h_cal": H_cal_etottracknorm_DATA.Fill,
        "h_cer": H_cer_npeSum_DATA.Fill,
        "p_cal": P_cal_etottracknorm_DATA.Fill,
        "p_hgcer": P_hgcer_npeSum_DATA.Fill,
        "p_aero": P_aero_npeSum_DATA.Fill,
    }
    data_nomm_fills = (
        H_MM_rand_dummy_DATA.Fill,
        H_MM_dummy_DATA.Fill,
        H_MM_full_DATA.Fill,
        H_MM_fit2sub_DATA.Fill,
        H_MM_fit1sub_DATA.Fill,
        H_MM_pisub_DATA.Fill,
        H_MM_nosub_DATA.Fill,
    )

    _process_rand_sub_tree(
        TBRANCH_DATA,
        "\nGrabbing {} {} data...".format(phi_setting, ParticleType),
        "rand_sub data loop {}".format(phi_setting),
        tmin,
        tmax,
        data_nomm_fills,
        data_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        Misc.progressBar,
    )
    # All experimental values now arrive in the kaon analysis coordinate.
    # Keep the historical value for audit only; passing it again would shift
    # prompt trees twice and leave pion control trees inconsistent.
    event_derived_mm_offset = None
    MM_offset_DATA = 0.0

    ################################################################################################################################################
    # Fill dummy histograms for various trees called above

    if ParticleType == "kaon":
        dummy_nohole_xy_fill = P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY.Fill
        dummy_nohole_x_mm_fill = P_hgcer_nohole_xAtCer_vs_MM_DUMMY.Fill
        dummy_nohole_y_mm_fill = P_hgcer_nohole_yAtCer_vs_MM_DUMMY.Fill
    else:
        dummy_nohole_xy_fill = None
        dummy_nohole_x_mm_fill = None
        dummy_nohole_y_mm_fill = None

    dummy_fills = {
        "nohole_xy": dummy_nohole_xy_fill,
        "nohole_x_mm": dummy_nohole_x_mm_fill,
        "nohole_y_mm": dummy_nohole_y_mm_fill,
        "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DUMMY.Fill,
        "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DUMMY.Fill,
        "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DUMMY.Fill,
        "mm_ct": MM_vs_CoinTime_DUMMY.Fill,
        "ct_beta": CoinTime_vs_beta_DUMMY.Fill,
        "mm_beta": MM_vs_beta_DUMMY.Fill,
        "mm_h_cer": MM_vs_H_cer_DUMMY.Fill,
        "mm_h_cal": MM_vs_H_cal_DUMMY.Fill,
        "mm_p_cal": MM_vs_P_cal_DUMMY.Fill,
        "mm_p_hgcer": MM_vs_P_hgcer_DUMMY.Fill,
        "mm_p_aero": MM_vs_P_aero_DUMMY.Fill,
        "phiq_t": phiq_vs_t_DUMMY.Fill,
        "q2_w": Q2_vs_W_DUMMY.Fill,
        "q2_t": Q2_vs_t_DUMMY.Fill,
        "w_t": W_vs_t_DUMMY.Fill,
        "eps_t": EPS_vs_t_DUMMY.Fill,
        "mm_t": MM_vs_t_DUMMY.Fill,
        "polar_graph": polar_phiq_vs_t_DUMMY,
        "h_ct": H_ct_DUMMY.Fill,
        "h_ssxfp": H_ssxfp_DUMMY.Fill,
        "h_ssyfp": H_ssyfp_DUMMY.Fill,
        "h_ssxpfp": H_ssxpfp_DUMMY.Fill,
        "h_ssypfp": H_ssypfp_DUMMY.Fill,
        "h_ssdelta": H_ssdelta_DUMMY.Fill,
        "h_ssxptar": H_ssxptar_DUMMY.Fill,
        "h_ssyptar": H_ssyptar_DUMMY.Fill,
        "h_hsxfp": H_hsxfp_DUMMY.Fill,
        "h_hsyfp": H_hsyfp_DUMMY.Fill,
        "h_hsxpfp": H_hsxpfp_DUMMY.Fill,
        "h_hsypfp": H_hsypfp_DUMMY.Fill,
        "h_hsdelta": H_hsdelta_DUMMY.Fill,
        "h_hsxptar": H_hsxptar_DUMMY.Fill,
        "h_hsyptar": H_hsyptar_DUMMY.Fill,
        "h_ph_q": H_ph_q_DUMMY.Fill,
        "h_th_q": H_th_q_DUMMY.Fill,
        "h_ph_recoil": H_ph_recoil_DUMMY.Fill,
        "h_th_recoil": H_th_recoil_DUMMY.Fill,
        "h_pmiss": H_pmiss_DUMMY.Fill,
        "h_emiss": H_emiss_DUMMY.Fill,
        "h_pmx": H_pmx_DUMMY.Fill,
        "h_pmy": H_pmy_DUMMY.Fill,
        "h_pmz": H_pmz_DUMMY.Fill,
        "h_q2": H_Q2_DUMMY.Fill,
        "h_t": H_t_DUMMY.Fill,
        "h_w": H_W_DUMMY.Fill,
        "h_epsilon": H_epsilon_DUMMY.Fill,
        "h_mm": H_MM_DUMMY.Fill,
        "h_cal": None,
        "h_cer": None,
        "p_cal": None,
        "p_hgcer": None,
        "p_aero": None,
    }
    dummy_nomm_fills = (
        H_MM_full_DUMMY.Fill,
        H_MM_fit2sub_DUMMY.Fill,
        H_MM_fit1sub_DUMMY.Fill,
        H_MM_pisub_DUMMY.Fill,
        H_MM_nosub_DUMMY.Fill,
    )

    _process_rand_sub_tree(
        TBRANCH_DUMMY,
        "\nGrabbing {} {} dummy...".format(phi_setting, ParticleType),
        "rand_sub dummy loop {}".format(phi_setting),
        tmin,
        tmax,
        dummy_nomm_fills,
        dummy_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        Misc.progressBar,
    )

    ###################################################################################################################################################    
    # Fill random histograms for various trees called above

    if ParticleType == "kaon":
        rand_nohole_xy_fill = P_hgcer_nohole_xAtCer_vs_yAtCer_RAND.Fill
        rand_nohole_x_mm_fill = P_hgcer_nohole_xAtCer_vs_MM_RAND.Fill
        rand_nohole_y_mm_fill = P_hgcer_nohole_yAtCer_vs_MM_RAND.Fill
    else:
        rand_nohole_xy_fill = None
        rand_nohole_x_mm_fill = None
        rand_nohole_y_mm_fill = None

    rand_fills = {
        "nohole_xy": rand_nohole_xy_fill,
        "nohole_x_mm": rand_nohole_x_mm_fill,
        "nohole_y_mm": rand_nohole_y_mm_fill,
        "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_RAND.Fill,
        "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_RAND.Fill,
        "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_RAND.Fill,
        "mm_ct": MM_vs_CoinTime_RAND.Fill,
        "ct_beta": CoinTime_vs_beta_RAND.Fill,
        "mm_beta": MM_vs_beta_RAND.Fill,
        "mm_h_cer": MM_vs_H_cer_RAND.Fill,
        "mm_h_cal": MM_vs_H_cal_RAND.Fill,
        "mm_p_cal": MM_vs_P_cal_RAND.Fill,
        "mm_p_hgcer": MM_vs_P_hgcer_RAND.Fill,
        "mm_p_aero": MM_vs_P_aero_RAND.Fill,
        "phiq_t": phiq_vs_t_RAND.Fill,
        "q2_w": Q2_vs_W_RAND.Fill,
        "q2_t": Q2_vs_t_RAND.Fill,
        "w_t": W_vs_t_RAND.Fill,
        "eps_t": EPS_vs_t_RAND.Fill,
        "mm_t": MM_vs_t_RAND.Fill,
        "polar_graph": None,
        "h_ct": H_ct_RAND.Fill,
        "h_ssxfp": H_ssxfp_RAND.Fill,
        "h_ssyfp": H_ssyfp_RAND.Fill,
        "h_ssxpfp": H_ssxpfp_RAND.Fill,
        "h_ssypfp": H_ssypfp_RAND.Fill,
        "h_ssdelta": H_ssdelta_RAND.Fill,
        "h_ssxptar": H_ssxptar_RAND.Fill,
        "h_ssyptar": H_ssyptar_RAND.Fill,
        "h_hsxfp": H_hsxfp_RAND.Fill,
        "h_hsyfp": H_hsyfp_RAND.Fill,
        "h_hsxpfp": H_hsxpfp_RAND.Fill,
        "h_hsypfp": H_hsypfp_RAND.Fill,
        "h_hsdelta": H_hsdelta_RAND.Fill,
        "h_hsxptar": H_hsxptar_RAND.Fill,
        "h_hsyptar": H_hsyptar_RAND.Fill,
        "h_ph_q": H_ph_q_RAND.Fill,
        "h_th_q": H_th_q_RAND.Fill,
        "h_ph_recoil": H_ph_recoil_RAND.Fill,
        "h_th_recoil": H_th_recoil_RAND.Fill,
        "h_pmiss": H_pmiss_RAND.Fill,
        "h_emiss": H_emiss_RAND.Fill,
        "h_pmx": H_pmx_RAND.Fill,
        "h_pmy": H_pmy_RAND.Fill,
        "h_pmz": H_pmz_RAND.Fill,
        "h_q2": H_Q2_RAND.Fill,
        "h_t": H_t_RAND.Fill,
        "h_w": H_W_RAND.Fill,
        "h_epsilon": H_epsilon_RAND.Fill,
        "h_mm": H_MM_RAND.Fill,
        "h_cal": None,
        "h_cer": None,
        "p_cal": None,
        "p_hgcer": None,
        "p_aero": None,
    }
    rand_nomm_fills = (
        H_MM_rand_dummy_RAND.Fill,
        H_MM_dummy_RAND.Fill,
        H_MM_full_RAND.Fill,
        H_MM_fit2sub_RAND.Fill,
        H_MM_fit1sub_RAND.Fill,
        H_MM_pisub_RAND.Fill,
        H_MM_nosub_RAND.Fill,
    )

    _process_rand_sub_tree(
        TBRANCH_RAND,
        "\nGrabbing {} {} random data...".format(phi_setting, ParticleType),
        "rand_sub random loop {}".format(phi_setting),
        tmin,
        tmax,
        rand_nomm_fills,
        rand_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        Misc.progressBar,
    )

    ###################################################################################################################################################    
    # Fill dummy random histograms for various trees called above

    if ParticleType == "kaon":
        dummy_rand_nohole_xy_fill = P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND.Fill
        dummy_rand_nohole_x_mm_fill = P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND.Fill
        dummy_rand_nohole_y_mm_fill = P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND.Fill
    else:
        dummy_rand_nohole_xy_fill = None
        dummy_rand_nohole_x_mm_fill = None
        dummy_rand_nohole_y_mm_fill = None

    dummy_rand_fills = {
        "nohole_xy": dummy_rand_nohole_xy_fill,
        "nohole_x_mm": dummy_rand_nohole_x_mm_fill,
        "nohole_y_mm": dummy_rand_nohole_y_mm_fill,
        "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND.Fill,
        "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DUMMY_RAND.Fill,
        "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DUMMY_RAND.Fill,
        "mm_ct": MM_vs_CoinTime_DUMMY_RAND.Fill,
        "ct_beta": CoinTime_vs_beta_DUMMY_RAND.Fill,
        "mm_beta": MM_vs_beta_DUMMY_RAND.Fill,
        "mm_h_cer": MM_vs_H_cer_DUMMY_RAND.Fill,
        "mm_h_cal": MM_vs_H_cal_DUMMY_RAND.Fill,
        "mm_p_cal": MM_vs_P_cal_DUMMY_RAND.Fill,
        "mm_p_hgcer": MM_vs_P_hgcer_DUMMY_RAND.Fill,
        "mm_p_aero": MM_vs_P_aero_DUMMY_RAND.Fill,
        "phiq_t": phiq_vs_t_DUMMY_RAND.Fill,
        "q2_w": Q2_vs_W_DUMMY_RAND.Fill,
        "q2_t": Q2_vs_t_DUMMY_RAND.Fill,
        "w_t": W_vs_t_DUMMY_RAND.Fill,
        "eps_t": EPS_vs_t_DUMMY_RAND.Fill,
        "mm_t": MM_vs_t_DUMMY_RAND.Fill,
        "polar_graph": None,
        "h_ct": H_ct_DUMMY_RAND.Fill,
        "h_ssxfp": H_ssxfp_DUMMY_RAND.Fill,
        "h_ssyfp": H_ssyfp_DUMMY_RAND.Fill,
        "h_ssxpfp": H_ssxpfp_DUMMY_RAND.Fill,
        "h_ssypfp": H_ssypfp_DUMMY_RAND.Fill,
        "h_ssdelta": H_ssdelta_DUMMY_RAND.Fill,
        "h_ssxptar": H_ssxptar_DUMMY_RAND.Fill,
        "h_ssyptar": H_ssyptar_DUMMY_RAND.Fill,
        "h_hsxfp": H_hsxfp_DUMMY_RAND.Fill,
        "h_hsyfp": H_hsyfp_DUMMY_RAND.Fill,
        "h_hsxpfp": H_hsxpfp_DUMMY_RAND.Fill,
        "h_hsypfp": H_hsypfp_DUMMY_RAND.Fill,
        "h_hsdelta": H_hsdelta_DUMMY_RAND.Fill,
        "h_hsxptar": H_hsxptar_DUMMY_RAND.Fill,
        "h_hsyptar": H_hsyptar_DUMMY_RAND.Fill,
        "h_ph_q": H_ph_q_DUMMY_RAND.Fill,
        "h_th_q": H_th_q_DUMMY_RAND.Fill,
        "h_ph_recoil": H_ph_recoil_DUMMY_RAND.Fill,
        "h_th_recoil": H_th_recoil_DUMMY_RAND.Fill,
        "h_pmiss": H_pmiss_DUMMY_RAND.Fill,
        "h_emiss": H_emiss_DUMMY_RAND.Fill,
        "h_pmx": H_pmx_DUMMY_RAND.Fill,
        "h_pmy": H_pmy_DUMMY_RAND.Fill,
        "h_pmz": H_pmz_DUMMY_RAND.Fill,
        "h_q2": H_Q2_DUMMY_RAND.Fill,
        "h_t": H_t_DUMMY_RAND.Fill,
        "h_w": H_W_DUMMY_RAND.Fill,
        "h_epsilon": H_epsilon_DUMMY_RAND.Fill,
        "h_mm": H_MM_DUMMY_RAND.Fill,
        "h_cal": None,
        "h_cer": None,
        "p_cal": None,
        "p_hgcer": None,
        "p_aero": None,
    }
    dummy_rand_nomm_fills = (
        H_MM_full_DUMMY_RAND.Fill,
        H_MM_fit2sub_DUMMY_RAND.Fill,
        H_MM_fit1sub_DUMMY_RAND.Fill,
        H_MM_pisub_DUMMY_RAND.Fill,
        H_MM_nosub_DUMMY_RAND.Fill,
    )

    _process_rand_sub_tree(
        TBRANCH_DUMMY_RAND,
        "\nGrabbing {} {} dummy random data...".format(phi_setting, ParticleType),
        "rand_sub dummy random loop {}".format(phi_setting),
        tmin,
        tmax,
        dummy_rand_nomm_fills,
        dummy_rand_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        mm_min,
        mm_max,
        Misc.progressBar,
    )
          
    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge    

    # Data Random subtraction window
    stage_start = perf_counter()
    P_hgcer_xAtCer_vs_yAtCer_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_RAND.Scale(1/nWindows)
    P_hgcer_xAtCer_vs_MM_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_RAND.Scale(1/nWindows)
    P_hgcer_yAtCer_vs_MM_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_RAND.Scale(1/nWindows)        
    MM_vs_CoinTime_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_RAND.Scale(1/nWindows)
    MM_vs_beta_RAND.Scale(1/nWindows)
    MM_vs_H_cer_RAND.Scale(1/nWindows)
    MM_vs_H_cal_RAND.Scale(1/nWindows)
    MM_vs_P_cal_RAND.Scale(1/nWindows)
    MM_vs_P_hgcer_RAND.Scale(1/nWindows)
    MM_vs_P_aero_RAND.Scale(1/nWindows)
    phiq_vs_t_RAND.Scale(1/nWindows)
    Q2_vs_W_RAND.Scale(1/nWindows)
    Q2_vs_t_RAND.Scale(1/nWindows)
    W_vs_t_RAND.Scale(1/nWindows)
    EPS_vs_t_RAND.Scale(1/nWindows)
    MM_vs_t_RAND.Scale(1/nWindows)    
    H_ct_RAND.Scale(1/nWindows)
    H_ssxfp_RAND.Scale(1/nWindows)
    H_ssyfp_RAND.Scale(1/nWindows)
    H_ssxpfp_RAND.Scale(1/nWindows)
    H_ssypfp_RAND.Scale(1/nWindows)
    H_hsxfp_RAND.Scale(1/nWindows)
    H_hsyfp_RAND.Scale(1/nWindows)
    H_hsxpfp_RAND.Scale(1/nWindows)
    H_hsypfp_RAND.Scale(1/nWindows)
    H_ssxptar_RAND.Scale(1/nWindows)
    H_ssyptar_RAND.Scale(1/nWindows)
    H_hsxptar_RAND.Scale(1/nWindows)
    H_hsyptar_RAND.Scale(1/nWindows)
    H_ssdelta_RAND.Scale(1/nWindows)
    H_hsdelta_RAND.Scale(1/nWindows)
    H_ph_q_RAND.Scale(1/nWindows)
    H_th_q_RAND.Scale(1/nWindows)
    H_ph_recoil_RAND.Scale(1/nWindows)
    H_th_recoil_RAND.Scale(1/nWindows)
    H_Q2_RAND.Scale(1/nWindows)
    H_W_RAND.Scale(1/nWindows)    
    H_t_RAND.Scale(1/nWindows)
    H_epsilon_RAND.Scale(1/nWindows)
    H_MM_RAND.Scale(1/nWindows)
    H_MM_full_RAND.Scale(1/nWindows)
    H_MM_rand_dummy_RAND.Scale(1/nWindows)
    H_MM_dummy_RAND.Scale(1/nWindows)
    H_MM_fit2sub_RAND.Scale(1/nWindows)
    H_MM_fit1sub_RAND.Scale(1/nWindows)
    H_MM_pisub_RAND.Scale(1/nWindows)
    H_MM_nosub_RAND.Scale(1/nWindows)
    H_pmiss_RAND.Scale(1/nWindows)
    H_emiss_RAND.Scale(1/nWindows)
    H_pmx_RAND.Scale(1/nWindows)
    H_pmy_RAND.Scale(1/nWindows)
    H_pmz_RAND.Scale(1/nWindows)

    # Data Dummy_Random subtraction window
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND.Scale(1/nWindows)
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)            
    MM_vs_CoinTime_DUMMY_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_H_cer_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_H_cal_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_P_cal_DUMMY_RAND.Scale(1/nWindows)    
    MM_vs_P_hgcer_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_P_aero_DUMMY_RAND.Scale(1/nWindows)    
    phiq_vs_t_DUMMY_RAND.Scale(1/nWindows)
    Q2_vs_W_DUMMY_RAND.Scale(1/nWindows)
    Q2_vs_t_DUMMY_RAND.Scale(1/nWindows)
    W_vs_t_DUMMY_RAND.Scale(1/nWindows)
    EPS_vs_t_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_t_DUMMY_RAND.Scale(1/nWindows)
    H_ssxfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssyfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssxpfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssypfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsxfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsyfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsxpfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsypfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssxptar_DUMMY_RAND.Scale(1/nWindows)
    H_ssyptar_DUMMY_RAND.Scale(1/nWindows)
    H_hsxptar_DUMMY_RAND.Scale(1/nWindows)
    H_hsyptar_DUMMY_RAND.Scale(1/nWindows)
    H_ssdelta_DUMMY_RAND.Scale(1/nWindows)
    H_hsdelta_DUMMY_RAND.Scale(1/nWindows)
    H_ph_q_DUMMY_RAND.Scale(1/nWindows)
    H_th_q_DUMMY_RAND.Scale(1/nWindows)
    H_ph_recoil_DUMMY_RAND.Scale(1/nWindows)
    H_th_recoil_DUMMY_RAND.Scale(1/nWindows)    
    H_Q2_DUMMY_RAND.Scale(1/nWindows)
    H_W_DUMMY_RAND.Scale(1/nWindows)
    H_t_DUMMY_RAND.Scale(1/nWindows)
    H_epsilon_DUMMY_RAND.Scale(1/nWindows)
    H_MM_DUMMY_RAND.Scale(1/nWindows)
    H_MM_full_DUMMY_RAND.Scale(1/nWindows)    
    H_MM_fit2sub_DUMMY_RAND.Scale(1/nWindows)
    H_MM_fit1sub_DUMMY_RAND.Scale(1/nWindows)
    H_MM_pisub_DUMMY_RAND.Scale(1/nWindows)
    H_MM_nosub_DUMMY_RAND.Scale(1/nWindows)
    H_pmiss_DUMMY_RAND.Scale(1/nWindows)
    H_emiss_DUMMY_RAND.Scale(1/nWindows)
    H_pmx_DUMMY_RAND.Scale(1/nWindows)
    H_pmy_DUMMY_RAND.Scale(1/nWindows)
    H_pmz_DUMMY_RAND.Scale(1/nWindows)
    #H_ct_DUMMY_RAND.Scale(1/nWindows)

    print("\n\n{} data total number of events (no subtraction): {:.3e}".format(phi_setting, H_MM_DATA.Integral()))
    print("{} dummy total number of events (no subtraction): {:.3e}".format(phi_setting, H_MM_DUMMY.Integral()))     

    ###
    # Data Random subtraction
    P_hgcer_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_xAtCer_vs_yAtCer_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_nohole_xAtCer_vs_yAtCer_RAND,-1)
    P_hgcer_xAtCer_vs_MM_DATA.Add(P_hgcer_xAtCer_vs_MM_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(P_hgcer_nohole_xAtCer_vs_MM_RAND,-1)
    P_hgcer_yAtCer_vs_MM_DATA.Add(P_hgcer_yAtCer_vs_MM_RAND,-1)
    if ParticleType == "kaon":    
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(P_hgcer_nohole_yAtCer_vs_MM_RAND,-1)        
    MM_vs_CoinTime_DATA.Add(MM_vs_CoinTime_RAND,-1)
    CoinTime_vs_beta_DATA.Add(CoinTime_vs_beta_RAND,-1)
    MM_vs_beta_DATA.Add(MM_vs_beta_RAND,-1)
    MM_vs_H_cer_DATA.Add(MM_vs_H_cer_RAND,-1)
    MM_vs_H_cal_DATA.Add(MM_vs_H_cal_RAND,-1)
    MM_vs_P_cal_DATA.Add(MM_vs_P_cal_RAND,-1)    
    MM_vs_P_hgcer_DATA.Add(MM_vs_P_hgcer_RAND,-1)
    MM_vs_P_aero_DATA.Add(MM_vs_P_aero_RAND,-1)
    phiq_vs_t_DATA.Add(phiq_vs_t_RAND,-1)
    Q2_vs_W_DATA.Add(Q2_vs_W_RAND,-1)
    Q2_vs_t_DATA.Add(Q2_vs_t_RAND,-1)
    W_vs_t_DATA.Add(W_vs_t_RAND,-1)
    EPS_vs_t_DATA.Add(EPS_vs_t_RAND,-1)
    MM_vs_t_DATA.Add(MM_vs_t_RAND,-1)    
    H_ssxfp_DATA.Add(H_ssxfp_RAND,-1)
    H_ssyfp_DATA.Add(H_ssyfp_RAND,-1)
    H_ssxpfp_DATA.Add(H_ssxpfp_RAND,-1)
    H_ssypfp_DATA.Add(H_ssypfp_RAND,-1)
    H_hsxfp_DATA.Add(H_hsxfp_RAND,-1)
    H_hsyfp_DATA.Add(H_hsyfp_RAND,-1)
    H_hsxpfp_DATA.Add(H_hsxpfp_RAND,-1)
    H_hsypfp_DATA.Add(H_hsypfp_RAND,-1)
    H_ssxptar_DATA.Add(H_ssxptar_RAND,-1)
    H_ssyptar_DATA.Add(H_ssyptar_RAND,-1)
    H_hsxptar_DATA.Add(H_hsxptar_RAND,-1)
    H_hsyptar_DATA.Add(H_hsyptar_RAND,-1)
    H_ssdelta_DATA.Add(H_ssdelta_RAND,-1)
    H_hsdelta_DATA.Add(H_hsdelta_RAND,-1)
    H_ph_q_DATA.Add(H_ph_q_RAND,-1)
    H_th_q_DATA.Add(H_th_q_RAND,-1)
    H_ph_recoil_DATA.Add(H_ph_recoil_RAND,-1)
    H_th_recoil_DATA.Add(H_th_recoil_RAND,-1)
    H_Q2_DATA.Add(H_Q2_RAND,-1)
    H_W_DATA.Add(H_W_RAND,-1)
    H_t_DATA.Add(H_t_RAND,-1)
    H_epsilon_DATA.Add(H_epsilon_RAND,-1)
    H_MM_DATA.Add(H_MM_RAND,-1)
    H_MM_full_DATA.Add(H_MM_full_RAND,-1)
    H_MM_dummy_DATA.Add(H_MM_dummy_RAND,-1)
    H_MM_fit2sub_DATA.Add(H_MM_fit2sub_RAND,-1)
    H_MM_fit1sub_DATA.Add(H_MM_fit1sub_RAND,-1)
    H_MM_pisub_DATA.Add(H_MM_pisub_RAND,-1)
    H_MM_nosub_DATA.Add(H_MM_nosub_RAND,-1)
    H_pmiss_DATA.Add(H_pmiss_RAND,-1)
    H_emiss_DATA.Add(H_emiss_RAND,-1)
    H_pmx_DATA.Add(H_pmx_RAND,-1)
    H_pmy_DATA.Add(H_pmy_RAND,-1)
    H_pmz_DATA.Add(H_pmz_RAND,-1)
    H_ct_DATA.Add(H_ct_RAND,-1)

    ###
    # Dummy Random subtraction
    P_hgcer_xAtCer_vs_yAtCer_DUMMY.Add(P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY.Add(P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND,-1)
    P_hgcer_xAtCer_vs_MM_DUMMY.Add(P_hgcer_xAtCer_vs_MM_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY.Add(P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND,-1)
    P_hgcer_yAtCer_vs_MM_DUMMY.Add(P_hgcer_yAtCer_vs_MM_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY.Add(P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND,-1)                
    MM_vs_CoinTime_DUMMY.Add(MM_vs_CoinTime_DUMMY_RAND,-1)
    CoinTime_vs_beta_DUMMY.Add(CoinTime_vs_beta_DUMMY_RAND,-1)
    MM_vs_beta_DUMMY.Add(MM_vs_beta_DUMMY_RAND,-1)
    MM_vs_H_cer_DUMMY.Add(MM_vs_H_cer_DUMMY_RAND,-1)
    MM_vs_H_cal_DUMMY.Add(MM_vs_H_cal_DUMMY_RAND,-1)
    MM_vs_P_cal_DUMMY.Add(MM_vs_P_cal_DUMMY_RAND,-1)    
    MM_vs_P_hgcer_DUMMY.Add(MM_vs_P_hgcer_DUMMY_RAND,-1)
    MM_vs_P_aero_DUMMY.Add(MM_vs_P_aero_DUMMY_RAND,-1)    
    phiq_vs_t_DUMMY.Add(phiq_vs_t_DUMMY_RAND,-1)
    Q2_vs_W_DUMMY.Add(Q2_vs_W_DUMMY_RAND,-1)
    Q2_vs_t_DUMMY.Add(Q2_vs_t_DUMMY_RAND,-1)
    W_vs_t_DUMMY.Add(W_vs_t_DUMMY_RAND,-1)
    EPS_vs_t_DUMMY.Add(EPS_vs_t_DUMMY_RAND,-1)
    MM_vs_t_DUMMY.Add(MM_vs_t_DUMMY_RAND,-1)
    H_ssxfp_DUMMY.Add(H_ssxfp_DUMMY_RAND,-1)
    H_ssyfp_DUMMY.Add(H_ssyfp_DUMMY_RAND,-1)
    H_ssxpfp_DUMMY.Add(H_ssxpfp_DUMMY_RAND,-1)
    H_ssypfp_DUMMY.Add(H_ssypfp_DUMMY_RAND,-1)
    H_hsxfp_DUMMY.Add(H_hsxfp_DUMMY_RAND,-1)
    H_hsyfp_DUMMY.Add(H_hsyfp_DUMMY_RAND,-1)
    H_hsxpfp_DUMMY.Add(H_hsxpfp_DUMMY_RAND,-1)
    H_hsypfp_DUMMY.Add(H_hsypfp_DUMMY_RAND,-1)
    H_ssxptar_DUMMY.Add(H_ssxptar_DUMMY_RAND,-1)
    H_ssyptar_DUMMY.Add(H_ssyptar_DUMMY_RAND,-1)
    H_hsxptar_DUMMY.Add(H_hsxptar_DUMMY_RAND,-1)
    H_hsyptar_DUMMY.Add(H_hsyptar_DUMMY_RAND,-1)
    H_ssdelta_DUMMY.Add(H_ssdelta_DUMMY_RAND,-1)
    H_hsdelta_DUMMY.Add(H_hsdelta_DUMMY_RAND,-1)
    H_ph_q_DUMMY.Add(H_ph_q_DUMMY_RAND,-1)
    H_th_q_DUMMY.Add(H_th_q_DUMMY_RAND,-1)
    H_ph_recoil_DUMMY.Add(H_ph_recoil_DUMMY_RAND,-1)
    H_th_recoil_DUMMY.Add(H_th_recoil_DUMMY_RAND,-1)    
    H_Q2_DUMMY.Add(H_Q2_DUMMY_RAND,-1)
    H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
    H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
    H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)
    H_MM_DUMMY.Add(H_MM_DUMMY_RAND,-1)
    H_MM_full_DUMMY.Add(H_MM_full_DUMMY_RAND,-1)
    H_MM_fit2sub_DUMMY.Add(H_MM_fit2sub_DUMMY_RAND,-1)
    H_MM_fit1sub_DUMMY.Add(H_MM_fit1sub_DUMMY_RAND,-1)
    H_MM_pisub_DUMMY.Add(H_MM_pisub_DUMMY_RAND,-1)
    H_MM_nosub_DUMMY.Add(H_MM_nosub_DUMMY_RAND,-1)
    H_pmiss_DUMMY.Add(H_pmiss_DUMMY_RAND,-1)
    H_emiss_DUMMY.Add(H_emiss_DUMMY_RAND,-1)
    H_pmx_DUMMY.Add(H_pmx_DUMMY_RAND,-1)
    H_pmy_DUMMY.Add(H_pmy_DUMMY_RAND,-1)
    H_pmz_DUMMY.Add(H_pmz_DUMMY_RAND,-1)
    H_ct_DUMMY.Add(H_ct_DUMMY_RAND,-1)
    _print_rand_timer("rand_sub random-window subtraction {}".format(phi_setting), perf_counter() - stage_start)

    print("\n\n{} data total number of events (random subtraction only!): {:.3e}".format(phi_setting, H_MM_DATA.Integral()))
    print("{} dummy total number of events (random subtraction only!): {:.3e}".format(phi_setting, H_MM_DUMMY.Integral()))  

    ###################################################################################################################################
    # Apply the setting-level normalization/background treatment used by the
    # Step 5 data-vs-SIMC overlay plots.
    #
    # Important: this is not the same path used for the final ratios in Step 6.
    # The ratio code performs a separate per-(t,phi) yield extraction/background
    # correction in binning/calculate_yield.py.
    ###################################################################################################################################
    ###
    # Data Normalization
    stage_start = perf_counter()
    P_hgcer_xAtCer_vs_yAtCer_DATA.Scale(norm_factor_data)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Scale(norm_factor_data)
    P_hgcer_xAtCer_vs_MM_DATA.Scale(norm_factor_data)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Scale(norm_factor_data)
    P_hgcer_yAtCer_vs_MM_DATA.Scale(norm_factor_data)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Scale(norm_factor_data)
    MM_vs_CoinTime_DATA.Scale(norm_factor_data)
    CoinTime_vs_beta_DATA.Scale(norm_factor_data)
    MM_vs_beta_DATA.Scale(norm_factor_data)
    MM_vs_H_cer_DATA.Scale(norm_factor_data)
    MM_vs_H_cal_DATA.Scale(norm_factor_data)
    MM_vs_P_cal_DATA.Scale(norm_factor_data)
    MM_vs_P_hgcer_DATA.Scale(norm_factor_data)
    MM_vs_P_aero_DATA.Scale(norm_factor_data)
    phiq_vs_t_DATA.Scale(norm_factor_data)
    Q2_vs_W_DATA.Scale(norm_factor_data)
    Q2_vs_t_DATA.Scale(norm_factor_data)
    W_vs_t_DATA.Scale(norm_factor_data)
    EPS_vs_t_DATA.Scale(norm_factor_data)
    MM_vs_t_DATA.Scale(norm_factor_data)
    H_ssxfp_DATA.Scale(norm_factor_data)
    H_ssyfp_DATA.Scale(norm_factor_data)
    H_ssxpfp_DATA.Scale(norm_factor_data)
    H_ssypfp_DATA.Scale(norm_factor_data)
    H_hsxfp_DATA.Scale(norm_factor_data)
    H_hsyfp_DATA.Scale(norm_factor_data)
    H_hsxpfp_DATA.Scale(norm_factor_data)
    H_hsypfp_DATA.Scale(norm_factor_data)
    H_ssxptar_DATA.Scale(norm_factor_data)
    H_ssyptar_DATA.Scale(norm_factor_data)
    H_hsxptar_DATA.Scale(norm_factor_data)
    H_hsyptar_DATA.Scale(norm_factor_data)
    H_ssdelta_DATA.Scale(norm_factor_data)
    H_hsdelta_DATA.Scale(norm_factor_data)
    H_ph_q_DATA.Scale(norm_factor_data)
    H_th_q_DATA.Scale(norm_factor_data)
    H_ph_recoil_DATA.Scale(norm_factor_data)
    H_th_recoil_DATA.Scale(norm_factor_data)
    H_Q2_DATA.Scale(norm_factor_data)
    H_W_DATA.Scale(norm_factor_data)
    H_t_DATA.Scale(norm_factor_data)
    H_epsilon_DATA.Scale(norm_factor_data)
    H_MM_DATA.Scale(norm_factor_data)
    H_MM_full_DATA.Scale(norm_factor_data)
    H_MM_rand_dummy_DATA.Scale(norm_factor_data)
    H_MM_dummy_DATA.Scale(norm_factor_data)
    H_MM_fit2sub_DATA.Scale(norm_factor_data)
    H_MM_fit1sub_DATA.Scale(norm_factor_data)
    H_MM_pisub_DATA.Scale(norm_factor_data)
    H_MM_nosub_DATA.Scale(norm_factor_data)
    H_pmiss_DATA.Scale(norm_factor_data)
    H_emiss_DATA.Scale(norm_factor_data)
    H_pmx_DATA.Scale(norm_factor_data)
    H_pmy_DATA.Scale(norm_factor_data)
    H_pmz_DATA.Scale(norm_factor_data)
    H_ct_DATA.Scale(norm_factor_data)

    ###
    # Dummy Normalization
    P_hgcer_xAtCer_vs_yAtCer_DUMMY.Scale(norm_factor_dummy)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY.Scale(norm_factor_dummy)
    P_hgcer_xAtCer_vs_MM_DUMMY.Scale(norm_factor_dummy)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY.Scale(norm_factor_dummy)
    P_hgcer_yAtCer_vs_MM_DUMMY.Scale(norm_factor_dummy)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY.Scale(norm_factor_dummy)
    MM_vs_CoinTime_DUMMY.Scale(norm_factor_dummy)
    CoinTime_vs_beta_DUMMY.Scale(norm_factor_dummy)
    MM_vs_beta_DUMMY.Scale(norm_factor_dummy)
    MM_vs_H_cer_DUMMY.Scale(norm_factor_dummy)
    MM_vs_H_cal_DUMMY.Scale(norm_factor_dummy)
    MM_vs_P_cal_DUMMY.Scale(norm_factor_dummy)
    MM_vs_P_hgcer_DUMMY.Scale(norm_factor_dummy)
    MM_vs_P_aero_DUMMY.Scale(norm_factor_dummy)
    phiq_vs_t_DUMMY.Scale(norm_factor_dummy)
    Q2_vs_W_DUMMY.Scale(norm_factor_dummy)
    Q2_vs_t_DUMMY.Scale(norm_factor_dummy)
    W_vs_t_DUMMY.Scale(norm_factor_dummy)
    EPS_vs_t_DUMMY.Scale(norm_factor_dummy)
    MM_vs_t_DUMMY.Scale(norm_factor_dummy)
    H_ssxfp_DUMMY.Scale(norm_factor_dummy)
    H_ssyfp_DUMMY.Scale(norm_factor_dummy)
    H_ssxpfp_DUMMY.Scale(norm_factor_dummy)
    H_ssypfp_DUMMY.Scale(norm_factor_dummy)
    H_hsxfp_DUMMY.Scale(norm_factor_dummy)
    H_hsyfp_DUMMY.Scale(norm_factor_dummy)
    H_hsxpfp_DUMMY.Scale(norm_factor_dummy)
    H_hsypfp_DUMMY.Scale(norm_factor_dummy)
    H_ssxptar_DUMMY.Scale(norm_factor_dummy)
    H_ssyptar_DUMMY.Scale(norm_factor_dummy)
    H_hsxptar_DUMMY.Scale(norm_factor_dummy)
    H_hsyptar_DUMMY.Scale(norm_factor_dummy)
    H_ssdelta_DUMMY.Scale(norm_factor_dummy)
    H_hsdelta_DUMMY.Scale(norm_factor_dummy)
    H_ph_q_DUMMY.Scale(norm_factor_dummy)
    H_th_q_DUMMY.Scale(norm_factor_dummy)
    H_ph_recoil_DUMMY.Scale(norm_factor_dummy)
    H_th_recoil_DUMMY.Scale(norm_factor_dummy)
    H_Q2_DUMMY.Scale(norm_factor_dummy)
    H_W_DUMMY.Scale(norm_factor_dummy)
    H_t_DUMMY.Scale(norm_factor_dummy)
    H_epsilon_DUMMY.Scale(norm_factor_dummy)
    H_MM_DUMMY.Scale(norm_factor_dummy)
    H_MM_full_DUMMY.Scale(norm_factor_dummy)    
    H_MM_fit2sub_DUMMY.Scale(norm_factor_dummy)
    H_MM_fit1sub_DUMMY.Scale(norm_factor_dummy)
    H_MM_pisub_DUMMY.Scale(norm_factor_dummy)
    H_MM_nosub_DUMMY.Scale(norm_factor_dummy)
    H_pmiss_DUMMY.Scale(norm_factor_dummy)
    H_emiss_DUMMY.Scale(norm_factor_dummy)
    H_pmx_DUMMY.Scale(norm_factor_dummy)
    H_pmy_DUMMY.Scale(norm_factor_dummy)
    H_pmz_DUMMY.Scale(norm_factor_dummy)
    H_ct_DUMMY.Scale(norm_factor_dummy)

    ###
    # Dummy subtraction
    P_hgcer_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_xAtCer_vs_yAtCer_DUMMY,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY,-1)
    P_hgcer_xAtCer_vs_MM_DATA.Add(P_hgcer_xAtCer_vs_MM_DUMMY,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(P_hgcer_nohole_xAtCer_vs_MM_DUMMY,-1)
    P_hgcer_yAtCer_vs_MM_DATA.Add(P_hgcer_yAtCer_vs_MM_DUMMY,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(P_hgcer_nohole_yAtCer_vs_MM_DUMMY,-1)                
    MM_vs_CoinTime_DATA.Add(MM_vs_CoinTime_DUMMY,-1)
    CoinTime_vs_beta_DATA.Add(CoinTime_vs_beta_DUMMY,-1)
    MM_vs_beta_DATA.Add(MM_vs_beta_DUMMY,-1)
    MM_vs_H_cer_DATA.Add(MM_vs_H_cer_DUMMY,-1)
    MM_vs_H_cal_DATA.Add(MM_vs_H_cal_DUMMY,-1)
    MM_vs_P_cal_DATA.Add(MM_vs_P_cal_DUMMY,-1)    
    MM_vs_P_hgcer_DATA.Add(MM_vs_P_hgcer_DUMMY,-1)
    MM_vs_P_aero_DATA.Add(MM_vs_P_aero_DUMMY,-1)    
    phiq_vs_t_DATA.Add(phiq_vs_t_DUMMY,-1)
    Q2_vs_W_DATA.Add(Q2_vs_W_DUMMY,-1)
    Q2_vs_t_DATA.Add(Q2_vs_t_DUMMY,-1)
    W_vs_t_DATA.Add(W_vs_t_DUMMY,-1)
    EPS_vs_t_DATA.Add(EPS_vs_t_DUMMY,-1)
    MM_vs_t_DATA.Add(MM_vs_t_DUMMY,-1)
    H_ssxfp_DATA.Add(H_ssxfp_DUMMY,-1)
    H_ssyfp_DATA.Add(H_ssyfp_DUMMY,-1)
    H_ssxpfp_DATA.Add(H_ssxpfp_DUMMY,-1)
    H_ssypfp_DATA.Add(H_ssypfp_DUMMY,-1)
    H_hsxfp_DATA.Add(H_hsxfp_DUMMY,-1)
    H_hsyfp_DATA.Add(H_hsyfp_DUMMY,-1)
    H_hsxpfp_DATA.Add(H_hsxpfp_DUMMY,-1)
    H_hsypfp_DATA.Add(H_hsypfp_DUMMY,-1)
    H_ssxptar_DATA.Add(H_ssxptar_DUMMY,-1)
    H_ssyptar_DATA.Add(H_ssyptar_DUMMY,-1)
    H_hsxptar_DATA.Add(H_hsxptar_DUMMY,-1)
    H_hsyptar_DATA.Add(H_hsyptar_DUMMY,-1)
    H_ssdelta_DATA.Add(H_ssdelta_DUMMY,-1)
    H_hsdelta_DATA.Add(H_hsdelta_DUMMY,-1)
    H_ph_q_DATA.Add(H_ph_q_DUMMY,-1)
    H_th_q_DATA.Add(H_th_q_DUMMY,-1)
    H_ph_recoil_DATA.Add(H_ph_recoil_DUMMY,-1)
    H_th_recoil_DATA.Add(H_th_recoil_DUMMY,-1)    
    H_Q2_DATA.Add(H_Q2_DUMMY,-1)
    H_W_DATA.Add(H_W_DUMMY,-1)
    H_t_DATA.Add(H_t_DUMMY,-1)
    H_epsilon_DATA.Add(H_epsilon_DUMMY,-1)
    H_MM_fit2sub_DATA.Add(H_MM_fit2sub_DUMMY,-1)
    H_MM_fit1sub_DATA.Add(H_MM_fit1sub_DUMMY,-1)
    H_MM_pisub_DATA.Add(H_MM_pisub_DUMMY,-1)
    H_MM_nosub_DATA.Add(H_MM_nosub_DUMMY,-1)
    H_MM_DATA.Add(H_MM_DUMMY,-1)
    H_MM_full_DATA.Add(H_MM_full_DUMMY,-1)        
    H_pmiss_DATA.Add(H_pmiss_DUMMY,-1)
    H_emiss_DATA.Add(H_emiss_DUMMY,-1)
    H_pmx_DATA.Add(H_pmx_DUMMY,-1)
    H_pmy_DATA.Add(H_pmy_DUMMY,-1)
    H_pmz_DATA.Add(H_pmz_DUMMY,-1)
    H_ct_DATA.Add(H_ct_DUMMY,-1)      
    _print_rand_timer("rand_sub norm/dummy subtraction {}".format(phi_setting), perf_counter() - stage_start)

    print("\n\n{} data total number of events (dummy & random subtraction): {:.3e}".format(phi_setting, H_MM_DATA.Integral()))
    print("{} dummy total number of events (dummy & random subtraction): {:.3e}".format(phi_setting, H_MM_DUMMY.Integral()))      

    binning_t_hist = None
    binning_phi_hist = None
    if ParticleType == "kaon":
        # Freeze the bin-counting inputs before slow-proton and pion subtraction.
        binning_t_hist = clone_root_histogram(
            H_t_DATA,
            scope="{}_setting".format(phi_setting),
            role="pre_particle_subtraction_t",
            name="H_t_DATA_pre_particle_subtraction_{}".format(phi_setting),
            sumw2=False,
        )
        binning_phi_hist = clone_root_histogram(
            H_ph_q_DATA,
            scope="{}_setting".format(phi_setting),
            role="pre_particle_subtraction_phi",
            name="H_ph_q_DATA_pre_particle_subtraction_{}".format(phi_setting),
            sumw2=False,
        )

    component_fit_result = None
    component_subtraction_payload = None
    component_diagnostic_payload = None
    sub_tree_bundle = None
    proton_cleaning_result = None
    proton_cleaning_application = None
    pion_hgcer_tdelta_diagnostic = None
    pion_hgcer_tdelta_json = None
    pion_hgcer_zerope_transfer = None
    pion_hgcer_zerope_transfer_json = None
    pion_hgcer_method_a = None
    pion_hgcer_method_b = None
    pion_hgcer_refinement_checkpoint = None
    pion_hgcer_refinement_checkpoint_json = None
    pion_hgcer_comparison_input = None
    pion_hgcer_method_a_comparison = None
    pion_hgcer_method_b_comparison = None
    pion_hgcer_ab_comparison = None
    pion_hgcer_phase_d_checkpoint = None
    pion_hgcer_phase_d_checkpoint_json = None

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        stage_start = perf_counter()
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        subDict["MM_offset_DATA"] = MM_offset_DATA
        subDict["legacy_MM_offset_DATA"] = {
            "effective_value": float(MM_offset_DATA or 0.0),
            "status": "disabled_by_kaon_data_coordinate"
            if coordinate_contract is not None else "legacy_raw_coordinate",
            "event_derived_value": (
                None if event_derived_mm_offset is None else float(event_derived_mm_offset)
            ),
            "coordinate_fingerprint": (
                coordinate_contract.get("coordinate_fingerprint")
                if coordinate_contract is not None else None
            ),
        }
        authoritative_component_t_bin = (
            resolve_particle_subtraction_mode(inpDict) == "simc_shape_components"
            and resolve_pion_subtraction_scope(inpDict) == "t_bin"
        )
        if not authoritative_component_t_bin:
            particle_subtraction_cuts(
                histDict, subDict, inpDict, SubtractedParticle, hgcer_cutg
            )

        if resolve_particle_subtraction_mode(inpDict) == "simc_shape_components":
            proton_cleaning_output_pdf = outputpdf.replace(
                "{}_FullAnalysis_".format(ParticleType),
                "{}_{}_rand_sub_".format(phi_setting, ParticleType),
            )
            print(
                "\nRunning integrated proton-cleaning diagnostics for {} {}...".format(
                    phi_setting,
                    ParticleType,
                )
            )
            proton_cleaning_tree_bundle = _open_kaon_proton_cleaning_tree_bundle(
                InFile_DATA,
                InFile_DUMMY,
                ParticleType,
                inpDict.get("POL"),
                EPSSET,
                phi_setting,
                norm_factor_data,
                norm_factor_dummy,
                nWindows,
            )
            proton_cleaning_tree_bundle["canonical_t_prepass_samples"] = dict(
                (inpDict.get("canonical_t_prepass_samples") or {}).get(phi_setting, {})
            )
            proton_cleaning_tree_bundle["canonical_t_prepass_sampling"] = dict(
                (inpDict.get("canonical_t_prepass_sampling") or {}).get(phi_setting, {})
            )
            proton_cleaning_tree_bundle = prepare_kaon_proton_cleaning_source_bundle(
                proton_cleaning_tree_bundle,
                evaluate_data_event,
                get_shifted_mm,
                get_shifted_t,
                hole_contains,
                mm_min,
                mm_max,
                proton_cleaning_config=get_proton_contamination_cleaning_config(
                    inp_dict=inpDict,
                    phi_setting=phi_setting,
                ),
            )
            proton_cleaning_result = build_kaon_proton_cleaning_result(
                inpDict,
                phi_setting,
                proton_cleaning_tree_bundle,
                evaluate_data_event,
                get_shifted_mm,
                hole_contains,
                mm_min,
                mm_max,
                analysis_scope="setting-wide",
                context="{}_{}_setting".format(phi_setting, EPSSET),
            )
            scope_payload = component_payload
            if scope_payload is None:
                scope_payload = load_setting_pion_component_shapes(
                    inpDict,
                    phi_setting,
                    particle_type=ParticleType,
                    context="rand_sub_setting_fit",
                )
            if kaon_signal_shape_payload is None:
                k_lambda_simc_root = os.path.join(
                    OUTPATH,
                    "{}_kaon_Simc_Q{}W{}_{}e.root".format(
                        phi_setting,
                        Q2,
                        W,
                        EPSSET,
                    ),
                )
                print(
                    "[SIMC K-Lambda] loading mandatory fallback payload for {} from {}".format(
                        phi_setting,
                        k_lambda_simc_root,
                    )
                )
                kaon_signal_shape_payload = load_kaon_simc_signal_shape(
                    k_lambda_simc_root,
                    inpDict,
                    phi_setting,
                    context="rand_sub_mandatory_k_lambda",
                )
            scope_shapes = resolve_scope_component_shapes(
                scope_payload,
                analysis_scope="setting-wide",
            )
            sub_tree_bundle = _open_subtracted_particle_tree_bundle(
                OUTPATH,
                phi_setting,
                SubtractedParticle,
                InDATAFilename,
                InDUMMYFilename,
                EPSSET,
            )
            component_targets = {
                "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DATA,
                "hgcer_xy_nohole": P_hgcer_nohole_xAtCer_vs_yAtCer_DATA if ParticleType == "kaon" else None,
                "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DATA,
                "hgcer_x_mm_nohole": P_hgcer_nohole_xAtCer_vs_MM_DATA if ParticleType == "kaon" else None,
                "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DATA,
                "hgcer_y_mm_nohole": P_hgcer_nohole_yAtCer_vs_MM_DATA if ParticleType == "kaon" else None,
                "mm_ct": MM_vs_CoinTime_DATA,
                "ct_beta": CoinTime_vs_beta_DATA,
                "mm_beta": MM_vs_beta_DATA,
                "mm_h_cer": MM_vs_H_cer_DATA,
                "mm_h_cal": MM_vs_H_cal_DATA,
                "mm_p_cal": MM_vs_P_cal_DATA,
                "mm_p_hgcer": MM_vs_P_hgcer_DATA,
                "mm_p_aero": MM_vs_P_aero_DATA,
                "phiq_t": phiq_vs_t_DATA,
                "q2_w": Q2_vs_W_DATA,
                "q2_t": Q2_vs_t_DATA,
                "w_t": W_vs_t_DATA,
                "eps_t": EPS_vs_t_DATA,
                "mm_t": MM_vs_t_DATA,
                "h_ct": H_ct_DATA,
                "h_ssxfp": H_ssxfp_DATA,
                "h_ssyfp": H_ssyfp_DATA,
                "h_ssxpfp": H_ssxpfp_DATA,
                "h_ssypfp": H_ssypfp_DATA,
                "h_hsxfp": H_hsxfp_DATA,
                "h_hsyfp": H_hsyfp_DATA,
                "h_hsxpfp": H_hsxpfp_DATA,
                "h_hsypfp": H_hsypfp_DATA,
                "h_ssxptar": H_ssxptar_DATA,
                "h_ssyptar": H_ssyptar_DATA,
                "h_hsxptar": H_hsxptar_DATA,
                "h_hsyptar": H_hsyptar_DATA,
                "h_ssdelta": H_ssdelta_DATA,
                "h_hsdelta": H_hsdelta_DATA,
                "h_ph_q": H_ph_q_DATA,
                "h_th_q": H_th_q_DATA,
                "h_ph_recoil": H_ph_recoil_DATA,
                "h_th_recoil": H_th_recoil_DATA,
                "h_q2": H_Q2_DATA,
                "h_t": H_t_DATA,
                "h_w": H_W_DATA,
                "h_epsilon": H_epsilon_DATA,
                "h_mm": H_MM_DATA,
                "h_mm_nosub": H_MM_nosub_DATA,
                "h_mm_fit2sub": H_MM_fit2sub_DATA,
                "h_mm_fit1sub": H_MM_fit1sub_DATA,
                "h_mm_pisub": H_MM_pisub_DATA,
                "h_mm_full": H_MM_full_DATA,
                "h_pmiss": H_pmiss_DATA,
                "h_emiss": H_emiss_DATA,
                "h_pmx": H_pmx_DATA,
                "h_pmy": H_pmy_DATA,
                "h_pmz": H_pmz_DATA,
            }
            active_component_targets = component_targets
            component_fit_kaon_input = H_MM_nosub_DATA
            component_input_metadata = {
                "input_selection": "no_rf_identity_no_proton_cleaning",
                "source_target_state": "post_proton_noRF",
            }
            pion_control_cache = None
            frozen_t_bins = None
            frozen_phi_bins = None

            if isinstance(proton_cleaning_result, dict):
                if bool(proton_cleaning_result.get("accepted")):
                    proton_cleaning_application = apply_kaon_proton_cleaning_to_targets(
                        proton_cleaning_result,
                        proton_cleaning_tree_bundle,
                        component_targets,
                        evaluate_data_event,
                        get_shifted_mm,
                        get_shifted_t,
                        hole_contains,
                        mm_min,
                        mm_max,
                    )
                    if bool((proton_cleaning_application or {}).get("accepted")):
                        if coordinate_contract is not None:
                            proton_cleaning_application["kaon_data_coordinate"] = dict(
                                coordinate_contract
                            )
                            proton_cleaning_application["coordinate_fingerprint"] = (
                                coordinate_contract["coordinate_fingerprint"]
                            )
                            for product in (
                                proton_cleaning_application.get("canonical_t_products") or ()
                            ):
                                product["kaon_data_coordinate"] = dict(coordinate_contract)
                                product["coordinate_fingerprint"] = coordinate_contract[
                                    "coordinate_fingerprint"
                                ]
                        # Retain the producer's detached application record for
                        # proton closure diagnostics and its ordered PDF pages.
                        proton_cleaning_result["application"] = proton_cleaning_application
                        production_map_key = (
                            "H_proton_weight_vs_delta_t"
                            if str(proton_cleaning_result.get("method") or "") == "timing_t_event_weight"
                            else "H_proton_weight_vs_delta_aero"
                        )
                        production_map_name = (
                            "H_proton_weight_vs_delta_t_DATA"
                            if production_map_key.endswith("_t")
                            else "H_proton_weight_vs_delta_aero_DATA"
                        )
                        for key, clone_name in (
                            ("H_MM_before_proton_cleaning", "H_MM_before_proton_cleaning_DATA"),
                            ("H_MM_estimated_proton", "H_MM_estimated_proton_DATA"),
                            ("H_MM_after_proton_cleaning", "H_MM_nosub_proton_cleaned_DATA"),
                            ("H_MM_after_proton_cleaning_final_rf", "H_MM_nosub_proton_cleaned_final_RF_DATA"),
                            ("H_proton_fraction_vs_MM", "H_proton_fraction_vs_MM_DATA"),
                            ("H_proton_weight_vs_delta", "H_proton_weight_vs_delta_DATA"),
                            (production_map_key, production_map_name),
                        ):
                            proton_cleaning_application[key] = _clone_hist_detached(
                                proton_cleaning_application.get(key),
                                clone_name,
                            )
                            histDict[clone_name] = proton_cleaning_application.get(key)

                        active_component_targets = proton_cleaning_application.get("final_targets") or component_targets
                        component_fit_kaon_input = active_component_targets.get("h_mm_nosub") or H_MM_nosub_DATA
                        component_input_metadata = _post_proton_cleaning_input_metadata(
                            proton_cleaning_application
                        )

                        P_hgcer_xAtCer_vs_yAtCer_DATA = active_component_targets.get("hgcer_xy")
                        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA = active_component_targets.get("hgcer_xy_nohole")
                        P_hgcer_xAtCer_vs_MM_DATA = active_component_targets.get("hgcer_x_mm")
                        P_hgcer_nohole_xAtCer_vs_MM_DATA = active_component_targets.get("hgcer_x_mm_nohole")
                        P_hgcer_yAtCer_vs_MM_DATA = active_component_targets.get("hgcer_y_mm")
                        P_hgcer_nohole_yAtCer_vs_MM_DATA = active_component_targets.get("hgcer_y_mm_nohole")
                        MM_vs_CoinTime_DATA = active_component_targets.get("mm_ct")
                        CoinTime_vs_beta_DATA = active_component_targets.get("ct_beta")
                        MM_vs_beta_DATA = active_component_targets.get("mm_beta")
                        MM_vs_H_cer_DATA = active_component_targets.get("mm_h_cer")
                        MM_vs_H_cal_DATA = active_component_targets.get("mm_h_cal")
                        MM_vs_P_cal_DATA = active_component_targets.get("mm_p_cal")
                        MM_vs_P_hgcer_DATA = active_component_targets.get("mm_p_hgcer")
                        MM_vs_P_aero_DATA = active_component_targets.get("mm_p_aero")
                        phiq_vs_t_DATA = active_component_targets.get("phiq_t")
                        Q2_vs_W_DATA = active_component_targets.get("q2_w")
                        Q2_vs_t_DATA = active_component_targets.get("q2_t")
                        W_vs_t_DATA = active_component_targets.get("w_t")
                        EPS_vs_t_DATA = active_component_targets.get("eps_t")
                        MM_vs_t_DATA = active_component_targets.get("mm_t")
                        H_ct_DATA = active_component_targets.get("h_ct")
                        H_ssxfp_DATA = active_component_targets.get("h_ssxfp")
                        H_ssyfp_DATA = active_component_targets.get("h_ssyfp")
                        H_ssxpfp_DATA = active_component_targets.get("h_ssxpfp")
                        H_ssypfp_DATA = active_component_targets.get("h_ssypfp")
                        H_hsxfp_DATA = active_component_targets.get("h_hsxfp")
                        H_hsyfp_DATA = active_component_targets.get("h_hsyfp")
                        H_hsxpfp_DATA = active_component_targets.get("h_hsxpfp")
                        H_hsypfp_DATA = active_component_targets.get("h_hsypfp")
                        H_ssxptar_DATA = active_component_targets.get("h_ssxptar")
                        H_ssyptar_DATA = active_component_targets.get("h_ssyptar")
                        H_hsxptar_DATA = active_component_targets.get("h_hsxptar")
                        H_hsyptar_DATA = active_component_targets.get("h_hsyptar")
                        H_ssdelta_DATA = active_component_targets.get("h_ssdelta")
                        H_hsdelta_DATA = active_component_targets.get("h_hsdelta")
                        H_ph_q_DATA = active_component_targets.get("h_ph_q")
                        H_th_q_DATA = active_component_targets.get("h_th_q")
                        H_ph_recoil_DATA = active_component_targets.get("h_ph_recoil")
                        H_th_recoil_DATA = active_component_targets.get("h_th_recoil")
                        H_Q2_DATA = active_component_targets.get("h_q2")
                        H_t_DATA = active_component_targets.get("h_t")
                        H_W_DATA = active_component_targets.get("h_w")
                        H_epsilon_DATA = active_component_targets.get("h_epsilon")
                        H_MM_DATA = active_component_targets.get("h_mm") or H_MM_DATA
                        H_MM_nosub_DATA = active_component_targets.get("h_mm_nosub") or H_MM_nosub_DATA
                        H_MM_fit2sub_DATA = active_component_targets.get("h_mm_fit2sub") or H_MM_fit2sub_DATA
                        H_MM_fit1sub_DATA = active_component_targets.get("h_mm_fit1sub") or H_MM_fit1sub_DATA
                        H_MM_pisub_DATA = active_component_targets.get("h_mm_pisub") or H_MM_pisub_DATA
                        H_MM_full_DATA = active_component_targets.get("h_mm_full") or H_MM_full_DATA
                        H_pmiss_DATA = active_component_targets.get("h_pmiss")
                        H_emiss_DATA = active_component_targets.get("h_emiss")
                        H_pmx_DATA = active_component_targets.get("h_pmx")
                        H_pmy_DATA = active_component_targets.get("h_pmy")
                        H_pmz_DATA = active_component_targets.get("h_pmz")

                if (
                    resolve_pion_subtraction_scope(inpDict) == "t_bin"
                    and str(ParticleType).strip().lower() == "kaon"
                ):
                    canonical_binning = inpDict.get("canonical_t_binning") or {}
                    frozen_t_bins = inpDict.get("t_bins")
                    frozen_phi_bins = inpDict.get("phi_bins")
                    if frozen_t_bins is None:
                        frozen_t_bins = canonical_binning.get("t_edges")
                    if frozen_phi_bins is None:
                        frozen_phi_bins = canonical_binning.get("phi_edges")
                    if frozen_t_bins is None or frozen_phi_bins is None:
                        raise RuntimeError(
                            "authoritative_t_bin_pion_parents_require_frozen_canonical_bins"
                        )
                    pion_hgcer_diagnostic_config = get_pion_hgcer_diagnostic_config(
                        inp_dict=inpDict,
                        phi_setting=phi_setting,
                    )
                    pion_hgcer_delta_edges, pion_hgcer_delta_edge_source = resolve_pion_hgcer_delta_edges(
                        pion_hgcer_diagnostic_config,
                        proton_cleaning_result,
                    )
                    if not bool((proton_cleaning_application or {}).get("accepted")):
                        if bool(proton_cleaning_result.get("accepted")):
                            raise RuntimeError(
                                "accepted_proton_result_missing_application"
                            )
                        proton_cleaning_application = (
                            _build_identity_no_proton_cleaning_application(
                                proton_source_bundle=proton_cleaning_tree_bundle,
                                target_templates=component_targets,
                                t_edges=frozen_t_bins,
                                delta_edges=pion_hgcer_delta_edges,
                                coordinate_fingerprint=(
                                    (coordinate_contract or {}).get(
                                        "coordinate_fingerprint"
                                    )
                                ),
                                proton_cleaning_result=proton_cleaning_result,
                            )
                        )
                    committed_host = finalize_committed_host_application(
                        proton_cleaning_result,
                        proton_cleaning_application,
                        component_targets,
                    )
                    if coordinate_contract is not None:
                        proton_cleaning_application["kaon_data_coordinate"] = dict(
                            coordinate_contract
                        )
                        proton_cleaning_application["coordinate_fingerprint"] = (
                            coordinate_contract["coordinate_fingerprint"]
                        )
                        for product in (
                            proton_cleaning_application.get("canonical_t_products")
                            or ()
                        ):
                            product["kaon_data_coordinate"] = dict(
                                coordinate_contract
                            )
                            product["coordinate_fingerprint"] = coordinate_contract[
                                "coordinate_fingerprint"
                            ]
                    proton_cleaning_result["application"] = (
                        proton_cleaning_application
                    )
                    component_input_metadata = _post_proton_cleaning_input_metadata(
                        proton_cleaning_application
                    )
                    proton_t_products = tuple(
                        proton_cleaning_application.get("canonical_t_products") or ()
                    )
                    if len(proton_t_products) != len(frozen_t_bins) - 1:
                        raise RuntimeError(
                            "authoritative_t_bin_pion_parents_require_direct_proton_products"
                        )
                    pion_control_cache = _build_authoritative_pion_control_source_cache(
                        sub_tree_bundle,
                        proton_t_products=proton_t_products,
                        t_bins=frozen_t_bins,
                        phi_bins=frozen_phi_bins,
                        particle_type=ParticleType,
                        pol=inpDict.get("POL"),
                        mm_offset_data=MM_offset_DATA,
                        coordinate_contract=coordinate_contract,
                        hole_contains=hole_contains,
                        evaluate_event=evaluate_data_event,
                        shifted_t_getter=get_shifted_t,
                        mm_min=mm_min,
                        mm_max=mm_max,
                        norm_factor_data=norm_factor_data,
                        norm_factor_dummy=norm_factor_dummy,
                        n_windows=nWindows,
                        delta_edges=pion_hgcer_delta_edges,
                    )
                    histDict["_authoritative_pion_control_source_cache"] = pion_control_cache
                    histDict["pion_control_source_audit"] = {
                        "definition": pion_control_cache.get("definition"),
                        "source_accounting": pion_control_cache.get("source_accounting"),
                        "coordinate_fingerprint": pion_control_cache.get("coordinate_fingerprint"),
                        "kaon_data_coordinate": pion_control_cache.get("kaon_data_coordinate"),
                        "physical_pion_control_mask": pion_control_cache.get("physical_pion_control_mask"),
                        "physical_pion_control_mask_fingerprint": pion_control_cache.get("physical_pion_control_mask_fingerprint"),
                    }
                    if bool(pion_hgcer_diagnostic_config.get("enabled", False)):
                        try:
                            pion_hgcer_tdelta_diagnostic = (
                                build_pion_hgcer_tdelta_diagnostic(
                                    kaon_source_bundle=proton_cleaning_tree_bundle,
                                    pion_tree_bundle=sub_tree_bundle,
                                    proton_cleaning_result=proton_cleaning_result,
                                    proton_coordinate_fingerprint=(
                                        (proton_cleaning_application or {}).get(
                                            "coordinate_fingerprint"
                                        )
                                    ),
                                    pion_control_cache=pion_control_cache,
                                    coordinate_contract=coordinate_contract,
                                    t_edges=frozen_t_bins,
                                    config=pion_hgcer_diagnostic_config,
                                    hole_contains=hole_contains,
                                    evaluate_pion_event=evaluate_data_event,
                                    mm_min=mm_min,
                                    mm_max=mm_max,
                                )
                            )
                        except Exception as exc:
                            # Part 1 is intentionally observational.  A
                            # missing branch or an unreadable diagnostic tree
                            # is reported, never promoted into a production
                            # pion/proton decision.
                            original_exception = getattr(
                                exc, "original_exception", exc
                            )
                            pion_hgcer_tdelta_diagnostic = {
                                "status": "unavailable",
                                "diagnostic_label": pion_hgcer_display_text("part1"),
                                "non_authoritative": True,
                                "production_side_effect_free": True,
                                "production_hgcer_pid_unchanged": True,
                                "reason": str(original_exception),
                                "exception_type": type(
                                    original_exception
                                ).__name__,
                                "exception_message": str(original_exception),
                                "diagnostic_stage": getattr(
                                    exc, "diagnostic_stage", "diagnostic_build"
                                ),
                                "config": pion_hgcer_diagnostic_config,
                                "source_provenance": getattr(
                                    exc, "source_provenance", None
                                ) or {},
                                "coordinate_fingerprint": (
                                    (coordinate_contract or {}).get("coordinate_fingerprint")
                                ),
                                "t_edges": (
                                    list(frozen_t_bins)
                                    if frozen_t_bins is not None
                                    else []
                                ),
                            }
                        pion_hgcer_tdelta_json = os.path.join(
                            OUTPATH,
                            "{}_{}_pion_hgcer_tdelta_diagnostic.json".format(
                                phi_setting, OutFilename
                            ),
                        )
                        write_pion_hgcer_tdelta_json(
                            pion_hgcer_tdelta_json,
                            pion_hgcer_tdelta_diagnostic,
                        )
                        histDict["pion_hgcer_tdelta_diagnostic"] = (
                            serialize_pion_hgcer_tdelta_diagnostic(
                                pion_hgcer_tdelta_diagnostic,
                                include_records=False,
                            )
                        )
                        histDict["pion_hgcer_tdelta_diagnostic_artifacts"] = [
                            pion_hgcer_tdelta_json
                        ]
                        # Part 2 consumes only the frozen Part-1 records and
                        # the already-authoritative cache.  Its every failure
                        # stays diagnostic-only and cannot alter the component
                        # fit, production pion weight, or proton products.
                        try:
                            pion_hgcer_transfer_config = get_pion_hgcer_transfer_config(
                                inp_dict=inpDict,
                                phi_setting=phi_setting,
                            )
                            if bool(pion_hgcer_transfer_config.get("enabled", False)):
                                pion_hgcer_zerope_transfer = build_pion_hgcer_zerope_transfer_map(
                                    pion_hgcer_tdelta_diagnostic,
                                    pion_control_cache,
                                    config=pion_hgcer_transfer_config,
                                )
                                pion_hgcer_zerope_transfer["application"] = (
                                    apply_frozen_pion_hgcer_transfer_map(
                                        pion_hgcer_zerope_transfer,
                                        pion_control_cache,
                                        proton_t_products,
                                    )
                                )
                        except Exception as exc:
                            pion_hgcer_zerope_transfer = {
                                "status": "unavailable",
                                "reason": str(exc),
                                "exception_type": type(exc).__name__,
                                "diagnostic_stage": "part2_build",
                                "non_authoritative": True,
                                "production_side_effect_free": True,
                                "production_pion_subtraction_unchanged": True,
                                "noRF_host_terminology": True,
                                "rf_restoration_applied": False,
                            }
                        if isinstance(pion_hgcer_zerope_transfer, dict):
                            pion_hgcer_zerope_transfer_json = os.path.join(
                                OUTPATH,
                                "{}_{}_pion_hgcer_zerope_transfer.json".format(
                                    phi_setting, OutFilename
                                ),
                            )
                            write_pion_hgcer_zerope_transfer_json(
                                pion_hgcer_zerope_transfer_json,
                                pion_hgcer_zerope_transfer,
                            )
                            histDict["pion_hgcer_zerope_transfer"] = (
                                serialize_pion_hgcer_zerope_transfer(
                                    pion_hgcer_zerope_transfer
                                )
                            )
                            histDict["pion_hgcer_zerope_transfer_artifacts"] = [
                                pion_hgcer_zerope_transfer_json
                            ]

                histDict["proton_contamination_cleaning_result_setting"] = (
                    serialize_kaon_proton_cleaning_result(proton_cleaning_result)
                )
                histDict["_proton_contamination_cleaning_result_setting"] = (
                    proton_cleaning_result
                )
                histDict["proton_contamination_cleaning_setting"] = (
                    summarize_kaon_proton_cleaning_result(proton_cleaning_result)
                )
                histDict["proton_contamination_cleaning_artifacts"] = (
                    _write_timing_t_validation_artifacts(
                        proton_cleaning_result,
                        outpath=OUTPATH,
                        particle_type=ParticleType,
                        outfilename=OutFilename,
                        epsset=EPSSET,
                        phi_setting=phi_setting,
                    )
                )
                print_kaon_proton_cleaning_terminal_summary(
                    proton_cleaning_result,
                    output_pdf=proton_cleaning_output_pdf,
                )

            setting_kaon_signal_shape = resolve_scope_single_shape(
                kaon_signal_shape_payload,
                analysis_scope="setting-wide",
            )
            component_pion_control_input = subDict["H_MM_nosub_SUB_DATA"]
            if resolve_pion_subtraction_scope(inpDict) == "t_bin":
                component_pion_control_input = (pion_control_cache or {}).get(
                    "H_pion_control_global"
                )
                if component_pion_control_input is None:
                    raise RuntimeError(
                        "authoritative_t_bin_setting_wide_diagnostic_requires_pion_control_cache"
                    )
            signal_shape_diagnostics = (
                dict(kaon_signal_shape_payload.get("diagnostics") or {})
                if isinstance(kaon_signal_shape_payload, dict)
                else {}
            )
            print(
                "[SIMC K-Lambda] setting-wide handoff phi={} root={} "
                "hist_present={} integral={:.6e} normalized={} fallback_reason={}".format(
                    phi_setting,
                    signal_shape_diagnostics.get("root_filename"),
                    setting_kaon_signal_shape is not None,
                    _hist_integral(setting_kaon_signal_shape),
                    signal_shape_diagnostics.get("normalized"),
                    signal_shape_diagnostics.get("fallback_reason") or "none",
                )
            )
            alignment_bin_key = {
                "kinematic_setting": get_particle_subtraction_setting_key(inpDict),
                "epsilon": EPSSET,
                "phi_setting": phi_setting,
                "analysis_scope": "setting-wide",
                "kaon_data_coordinate_fingerprint": str(
                    (coordinate_contract or {}).get("coordinate_fingerprint") or ""
                ),
                "t_bin": None,
                "phi_bin": None,
                "active_dimensions": {
                    "particle_type": ParticleType,
                    "polarization": inpDict.get("POL"),
                    "target": inpDict.get("target") or inpDict.get("Target"),
                },
            }
            setting_alignment = None
            alignment_status = "not_persisted"
            alignment_reasons = []
            alignment_paths = []
            dynamic_alignment_config = get_pion_component_dynamic_alignment_config(
                inp_dict=inpDict,
                phi_setting=phi_setting,
                setting_key=get_particle_subtraction_setting_key(inpDict),
            )
            if (
                bool(dynamic_alignment_config.get("enabled", False))
                and resolve_particle_subtraction_mode(inpDict) == PARTICLE_SUBTRACTION_MODE_COMPONENTS
            ):
                setting_alignment, alignment_status, alignment_reasons, alignment_paths = load_or_resolve_pion_component_alignment(
                    OUTPATH,
                    get_particle_subtraction_setting_key(inpDict),
                    phi_setting,
                    EPSSET,
                    "setting-wide",
                    alignment_bin_key,
                    component_pion_control_input,
                    scope_shapes,
                    inp_dict=inpDict,
                    # Model-only SIMC alignment retains the explicit kaon
                    # correction; experimental histograms are already in the
                    # analysis coordinate and receive no legacy offset.
                    common_setting_shift_gev=float(
                        (coordinate_contract or {}).get("mm_shift", 0.0)
                    ),
                )
            component_fit_result = build_particle_subtraction_component_result(
                component_pion_control_input,
                component_fit_kaon_input,
                scope_shapes,
                inpDict,
                analysis_scope="setting-wide",
                kaon_signal_shape=setting_kaon_signal_shape,
                kaon_sigma0_shape=resolve_scope_single_shape(
                    kaon_sigma0_shape_payload,
                    analysis_scope="setting-wide",
                ),
                kaon_sigma0_source_diagnostics=(
                    (kaon_sigma0_shape_payload or {}).get("diagnostics") or {}
                ),
                mm_offset_data=MM_offset_DATA,
                phi_setting=phi_setting,
                context="{}_{}_setting".format(phi_setting, EPSSET),
                pion_component_alignment=setting_alignment,
                alignment_bin_key=alignment_bin_key,
            )
            setting_wide_diagnostic_only = (
                resolve_pion_subtraction_scope(inpDict) == "t_bin"
            )
            component_fit_result["diagnostic_only"] = setting_wide_diagnostic_only
            component_fit_result["application_authoritative"] = not setting_wide_diagnostic_only
            component_fit_result["input_selection"] = component_input_metadata["input_selection"]
            component_fit_result["source_target_state"] = component_input_metadata["source_target_state"]
            alignment_payload = component_fit_result.get("pion_component_alignment")
            if isinstance(alignment_payload, dict):
                alignment_payload["persistence_status"] = alignment_status
                alignment_payload["persistence_rejection_reasons"] = list(alignment_reasons)
                alignment_payload["artifact_paths"] = list(alignment_paths)
                for path in alignment_paths:
                    if path not in inpDict.setdefault("pion_component_alignment_artifacts", []):
                        inpDict["pion_component_alignment_artifacts"].append(path)
            if setting_wide_diagnostic_only:
                # In t-bin production this fit remains a comparison diagnostic.
                # It must never alter a later child pion weight or spectrum.
                component_subtraction_payload = None
                if bool(inpDict.get("emit_setting_wide_pion_diagnostic", True)):
                    setting_parent_input = {
                        "H_proton_cleaned_final_rf_cut": active_component_targets.get("h_mm"),
                        "H_proton_cleaned_final_rf": active_component_targets.get("h_mm_nosub"),
                        "H_pion_control": (pion_control_cache or {}).get("H_pion_control_global"),
                        "pion_control_records": (pion_control_cache or {}).get("records"),
                        **component_input_metadata,
                    }
                    setting_scope = "pion_setting_wide_diagnostic"
                    proposal = _build_authoritative_parent_mm_diagnostic_proposal(
                        component_fit_result,
                        setting_parent_input,
                        inp_dict=inpDict,
                        scope=setting_scope,
                    )
                    setting_evaluation = evaluate_particle_subtraction_component_fit_result(
                        component_fit_result, inpDict
                    )
                    component_diagnostic_payload, setting_final_status = (
                        resolve_parent_diagnostic_final_application(
                            proposal,
                            setting_evaluation,
                            fallback_context=lambda: _build_authoritative_parent_single_scale_final(
                                proposal,
                                setting_evaluation,
                                setting_parent_input,
                                inp_dict=inpDict,
                                phi_setting=phi_setting,
                                mm_offset_data=MM_offset_DATA,
                                scope=setting_scope,
                            ),
                        )
                    )
                    if component_diagnostic_payload is None:
                        component_diagnostic_payload = proposal
                    component_diagnostic_payload["setting_wide_final_status"] = setting_final_status
                    component_diagnostic_payload["setting_wide_label"] = (
                        "SETTING-WIDE DIAGNOSTIC FIT - NON-AUTHORITATIVE"
                    )
            else:
                component_subtraction_payload = _apply_component_pion_subtraction_setting(
                    component_fit_result,
                    sub_tree_bundle,
                    phi_setting,
                    inpDict,
                    ParticleType,
                    MM_offset_DATA,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    norm_factor_data,
                    norm_factor_dummy,
                    nWindows,
                    active_component_targets,
                )
            histDict["_particle_subtraction_component_fit_setting"] = component_fit_result
            histDict["particle_subtraction_component_fit_setting"] = (
                serialize_particle_subtraction_component_result(component_fit_result)
            )
            histDict["_particle_subtraction_component_payload_setting"] = component_subtraction_payload
            histDict["particle_subtraction_component_payload_setting"] = summarize_particle_subtraction_component_payload(
                component_subtraction_payload
            )
            histDict["_particle_subtraction_component_diagnostic_payload_setting"] = component_diagnostic_payload
            histDict["particle_subtraction_component_diagnostic_payload_setting"] = (
                summarize_particle_subtraction_component_payload(component_diagnostic_payload)
            )
            if (
                resolve_pion_subtraction_scope(inpDict) == "t_bin"
                and str(ParticleType).strip().lower() == "kaon"
            ):
                canonical_binning = inpDict.get("canonical_t_binning") or {}
                proton_t_products = tuple(
                    (proton_cleaning_application or {}).get("canonical_t_products") or ()
                )
                if len(proton_t_products) != len(frozen_t_bins) - 1:
                    raise RuntimeError(
                        "authoritative_t_bin_pion_parents_require_direct_proton_products"
                    )
                if pion_control_cache is None:
                    raise RuntimeError(
                        "authoritative_t_bin_pion_parent_cache_not_built_at_cuts_stage"
                    )
                parent_inputs = []
                for t_index, proton_product in enumerate(proton_t_products):
                    final_t_targets = proton_product.get("final_targets") or {}
                    raw_t_targets = proton_product.get("raw_targets") or {}
                    proton_t_targets = proton_product.get("proton_targets") or {}
                    cleaned_t_targets = proton_product.get("cleaned_targets_pre_rf") or {}
                    control_t = pion_control_cache["by_t"][t_index]
                    final_proton_output = final_t_targets.get("h_mm_nosub")
                    if final_proton_output is None:
                        raise RuntimeError(
                            "authoritative_t_bin_parent_missing_final_proton_output:t{}".format(
                                t_index + 1
                            )
                        )
                    parent_inputs.append({
                        "t_index": int(t_index),
                        "t_edges": list(proton_product.get("t_edges") or ()),
                        "H_random_dummy_subtracted_kaon": raw_t_targets.get("h_mm_nosub"),
                        "H_proton_estimate": proton_t_targets.get("h_mm_nosub"),
                        "H_proton_cleaned": cleaned_t_targets.get("h_mm_nosub"),
                        "H_proton_cleaned_final_rf": final_proton_output,
                        "H_proton_cleaned_final_rf_cut": final_t_targets.get("h_mm"),
                        "H_pion_control": control_t.get("H_pion_control"),
                        "H_pion_control_cut": control_t.get("H_pion_control_cut"),
                        "pion_control_records": control_t.get("records"),
                        "source_accounting": proton_product.get("source_accounting"),
                        "pion_control_source_accounting": control_t.get("source_accounting"),
                        "kaon_data_coordinate": dict(
                            proton_product.get("kaon_data_coordinate") or {}
                        ),
                        "coordinate_fingerprint": str(
                            proton_product.get("coordinate_fingerprint") or ""
                        ),
                        "pion_control_coordinate_fingerprint": str(
                            control_t.get("coordinate_fingerprint")
                            or pion_control_cache.get("coordinate_fingerprint")
                            or ""
                        ),
                        # Propagate the proton-stage producer value unchanged;
                        # the parent boundary recomputes it and fails on drift.
                        "proton_final_output_fingerprint": str(
                            proton_product.get("final_output_fingerprint") or ""
                        ),
                        "pion_control_input_fingerprint": fingerprint_histogram_content_error(
                            control_t.get("H_pion_control")
                        ),
                        **component_input_metadata,
                        "source_epsilon": str(inpDict.get("EPSSET", "")).strip().lower(),
                        "consumer_epsilon": str(inpDict.get("EPSSET", "")).strip().lower(),
                        "canonical_interval_pair_id": canonical_binning.get("canonical_interval_pair_id"),
                        "canonical_interval_pair_hash": canonical_binning.get("canonical_interval_pair_hash"),
                    })

                def _build_parent_diagnostic_application(
                    *, fit_result, parent_input, t_index, t_edges
                ):
                    """Build proposal and final policy state on detached t sources."""
                    cut_source = parent_input.get("H_proton_cleaned_final_rf_cut")
                    full_source = parent_input.get("H_proton_cleaned_final_rf")
                    if cut_source is None or full_source is None:
                        raise RuntimeError("missing_t_integrated_parent_diagnostic_source")
                    diagnostic_policy = resolve_frozen_parent_application_policy(
                        {"fit_result": fit_result}, inpDict
                    )
                    production_evaluation = dict(diagnostic_policy.get("evaluation") or {})
                    parent_scope = "pion_parent_t{}".format(int(t_index) + 1)
                    proposal = _build_authoritative_parent_mm_diagnostic_proposal(
                        fit_result,
                        parent_input,
                        inp_dict=inpDict,
                        scope=parent_scope,
                    )
                    try:
                        final_payload, final_status = resolve_parent_diagnostic_final_application(
                            proposal,
                            production_evaluation,
                            fallback_context=lambda: _build_authoritative_parent_single_scale_final(
                                proposal,
                                production_evaluation,
                                parent_input,
                                inp_dict=inpDict,
                                phi_setting=phi_setting,
                                mm_offset_data=MM_offset_DATA,
                                scope=parent_scope,
                            ),
                        )
                    except Exception as exc:
                        unexpected = not any(
                            token in str(exc).lower()
                            for token in ("hist", "template", "tree", "event", "clone")
                        )
                        if bool(inpDict.get("pion_parent_diagnostic_strict", False)) and unexpected:
                            raise
                        final_payload = None
                        final_status = {
                            "status": "unavailable",
                            "final_status": "unavailable",
                            "final_reason": "fallback_build_failed",
                            "detail": str(exc),
                        }
                    final_status.setdefault("application_policy", diagnostic_policy)
                    return {
                        "proposed_diagnostic_application_result": proposal,
                        "proposed_diagnostic_application_status": {
                            "status": "available",
                            "reason": None,
                            "detail": None,
                        },
                        "final_diagnostic_application_result": final_payload,
                        "final_diagnostic_application_status": final_status,
                        "production_evaluation": production_evaluation,
                    }

                build_setting_t_bin_pion_parents(
                    histDict,
                    inpDict,
                    t_bins=frozen_t_bins,
                    parent_inputs=parent_inputs,
                    pion_component_shape_payload=scope_payload,
                    kaon_signal_shape_payload=kaon_signal_shape_payload,
                    kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
                    parent_pion_alignment=component_fit_result.get("pion_component_alignment"),
                    alignment_outpath=OUTPATH,
                    mm_offset_data=MM_offset_DATA,
                    coordinate_contract=coordinate_contract,
                    diagnostic_application_builder=_build_parent_diagnostic_application,
                )
                try:
                    pion_hgcer_event_contract = build_pion_hgcer_event_contract(
                        pion_control_cache=pion_control_cache,
                        pion_parents=histDict.get("_pion_t_bin_parent_results") or (),
                        canonical_t_global=histDict.get(
                            "_pion_authoritative_canonical_t_global"
                        ),
                        proton_source_bundle=proton_cleaning_tree_bundle,
                        proton_cleaning_result=proton_cleaning_result,
                        proton_cleaning_application=proton_cleaning_application,
                        inp_dict=inpDict,
                        canonical_binning=canonical_binning,
                        delta_edge_source=pion_hgcer_delta_edge_source,
                    )
                except Exception as exc:
                    # Phase A is a detached read-only proof.  An unexpected
                    # contract failure is retained as diagnostics and cannot
                    # disable or modify the already-established subtraction.
                    pion_hgcer_event_contract = {
                        "schema_version": "pion_hgcer_event_contract/v1",
                        "status": "unavailable",
                        "available": False,
                        "reason": str(exc),
                        "diagnostic_stage": "contract_build_exception",
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc),
                        "contract_fingerprint": None,
                        "pion_records": [],
                        "kaon_host_records": [],
                        "pion_closure": {"passed": False},
                        "host_closure": {"passed": False},
                        "production_objects_mutated": False,
                        "refinement_applied": False,
                    }
                histDict["_pion_hgcer_event_contract"] = pion_hgcer_event_contract
                histDict["pion_hgcer_event_contract_summary"] = (
                    summarize_pion_hgcer_event_contract(
                        pion_hgcer_event_contract
                    )
                )
                try:
                    pion_hgcer_method_a = build_pion_hgcer_method_a(
                        pion_hgcer_tdelta_diagnostic,
                        pion_hgcer_event_contract,
                    )
                except Exception as exc:
                    # Phase B is detached and diagnostic-only.  Its failure
                    # cannot modify or disable the frozen Phase-A subtraction.
                    pion_hgcer_method_a = {
                        "schema_version": "pion_hgcer_method_a/v1",
                        "method": "observed_positive_hgcer_response",
                        "status": "unavailable",
                        "available": False,
                        "reason": str(exc),
                        "diagnostic_stage": "runtime_build_exception",
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc),
                        "non_authoritative": True,
                        "production_objects_mutated": False,
                        "refinement_applied": False,
                        "rf_ct_required": False,
                        "zerope_model_used": False,
                        "t_edges": [],
                        "delta_edges": [],
                        "configuration": {},
                        "fingerprint": None,
                        "cells": [],
                        "summary": {},
                    }
                histDict["_pion_hgcer_method_a"] = pion_hgcer_method_a
                histDict["pion_hgcer_method_a_summary"] = (
                    summarize_pion_hgcer_method_a(pion_hgcer_method_a)
                )
                try:
                    pion_hgcer_method_b_config = (
                        resolve_pion_hgcer_method_b_config(
                            inp_dict=inpDict,
                            phi_setting=phi_setting,
                            mm_offset_data=MM_offset_DATA,
                        )
                    )
                    pion_hgcer_method_b = build_pion_hgcer_method_b(
                        pion_hgcer_event_contract,
                        config=pion_hgcer_method_b_config,
                    )
                except Exception as exc:
                    # Phase C is detached and diagnostic-only.  Preserve a
                    # structured unavailable result while production continues.
                    pion_hgcer_method_b = build_pion_hgcer_method_b(
                        pion_hgcer_event_contract,
                        config={},
                    )
                    pion_hgcer_method_b.update({
                        "reason": "runtime_method_b_configuration_exception",
                        "diagnostic_stage": "runtime_build_exception",
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc),
                    })
                histDict["_pion_hgcer_method_b"] = pion_hgcer_method_b
                histDict["pion_hgcer_method_b_summary"] = (
                    summarize_pion_hgcer_method_b(pion_hgcer_method_b)
                )
                try:
                    checkpoint_kinematic_token = (
                        get_particle_subtraction_setting_key(inpDict)
                    )
                    checkpoint_particle_type = str(ParticleType).strip().lower()
                    checkpoint_epsilon_setting = str(EPSSET).strip().lower()
                    checkpoint_epsilon_filename_token = "{}e".format(
                        checkpoint_epsilon_setting
                    )
                    checkpoint_setting = {
                        "kinematic_token": checkpoint_kinematic_token,
                        "Q2": Q2,
                        "W": W,
                        "epsilon_setting": checkpoint_epsilon_setting,
                        "epsilon_filename_token": checkpoint_epsilon_filename_token,
                        "phi_setting": phi_setting,
                        "particle_type": checkpoint_particle_type,
                    }
                    pion_hgcer_refinement_checkpoint = (
                        build_pion_hgcer_refinement_checkpoint(
                            setting=checkpoint_setting,
                            phase_a=pion_hgcer_event_contract,
                            phase_a_summary=histDict[
                                "pion_hgcer_event_contract_summary"
                            ],
                            method_a=pion_hgcer_method_a,
                            method_a_summary=histDict[
                                "pion_hgcer_method_a_summary"
                            ],
                            method_b=pion_hgcer_method_b,
                            method_b_summary=histDict[
                                "pion_hgcer_method_b_summary"
                            ],
                        )
                    )
                    pion_hgcer_refinement_checkpoint_json = os.path.join(
                        OUTPATH,
                        pion_hgcer_refinement_checkpoint_filename(
                            pion_hgcer_refinement_checkpoint["setting"]["phi_setting"],
                            pion_hgcer_refinement_checkpoint["setting"]["particle_type"],
                            pion_hgcer_refinement_checkpoint["setting"]["kinematic_token"],
                            pion_hgcer_refinement_checkpoint["setting"]["epsilon_filename_token"],
                        ),
                    )
                    write_pion_hgcer_refinement_checkpoint_json(
                        pion_hgcer_refinement_checkpoint_json,
                        pion_hgcer_refinement_checkpoint,
                    )
                    histDict["pion_hgcer_refinement_checkpoint"] = (
                        pion_hgcer_refinement_checkpoint
                    )
                    histDict["pion_hgcer_refinement_checkpoint_artifacts"] = [
                        pion_hgcer_refinement_checkpoint_json
                    ]
                    # Phase D.4 is a terminal detached diagnostic chain.  It
                    # consumes the in-memory Phase-C checkpoint only after the
                    # frozen Phase-C artifact has been written.
                    try:
                        pion_hgcer_comparison_input = (
                            build_pion_hgcer_comparison_input_contract(
                                pion_hgcer_refinement_checkpoint
                            )
                        )
                        pion_hgcer_method_a_comparison = (
                            build_pion_hgcer_method_a_comparison(
                                pion_hgcer_comparison_input
                            )
                        )
                        pion_hgcer_method_b_comparison = (
                            build_pion_hgcer_method_b_comparison(
                                pion_hgcer_comparison_input
                            )
                        )
                        pion_hgcer_ab_comparison = build_pion_hgcer_ab_comparison(
                            pion_hgcer_method_a_comparison,
                            pion_hgcer_method_b_comparison,
                        )
                        pion_hgcer_phase_d_checkpoint = (
                            build_pion_hgcer_phase_d_checkpoint(
                                setting=pion_hgcer_refinement_checkpoint["setting"],
                                method_a_comparison=pion_hgcer_method_a_comparison,
                                method_b_comparison=pion_hgcer_method_b_comparison,
                                ab_comparison=pion_hgcer_ab_comparison,
                            )
                        )
                        pion_hgcer_phase_d_checkpoint_json = os.path.join(
                            OUTPATH,
                            pion_hgcer_phase_d_checkpoint_filename(
                                pion_hgcer_phase_d_checkpoint["setting"]["phi_setting"],
                                pion_hgcer_phase_d_checkpoint["setting"]["particle_type"],
                                pion_hgcer_phase_d_checkpoint["setting"]["kinematic_token"],
                                pion_hgcer_phase_d_checkpoint["setting"]["epsilon_filename_token"],
                            ),
                        )
                        write_pion_hgcer_phase_d_checkpoint_json(
                            pion_hgcer_phase_d_checkpoint_json,
                            pion_hgcer_phase_d_checkpoint,
                        )
                        histDict["pion_hgcer_comparison_input"] = pion_hgcer_comparison_input
                        histDict["pion_hgcer_method_a_comparison"] = pion_hgcer_method_a_comparison
                        histDict["pion_hgcer_method_b_comparison"] = pion_hgcer_method_b_comparison
                        histDict["pion_hgcer_ab_comparison"] = pion_hgcer_ab_comparison
                        histDict["pion_hgcer_phase_d_checkpoint"] = pion_hgcer_phase_d_checkpoint
                        histDict["pion_hgcer_phase_d_checkpoint_artifacts"] = [
                            pion_hgcer_phase_d_checkpoint_json
                        ]
                    except Exception as exc:
                        histDict["pion_hgcer_phase_d_checkpoint_status"] = {
                            "status": "unavailable",
                            "available": False,
                            "reason": "phase_d_checkpoint_write_exception",
                            "diagnostic_stage": "runtime_phase_d_exception",
                            "exception_type": type(exc).__name__,
                            "exception_message": str(exc),
                            "non_authoritative": True,
                            "production_objects_mutated": False,
                            "refinement_applied": False,
                        }
                except Exception as exc:
                    # A checkpoint is a persistent diagnostic artifact only;
                    # failure to serialize it cannot alter production.
                    histDict["pion_hgcer_refinement_checkpoint_status"] = {
                        "status": "unavailable",
                        "available": False,
                        "reason": (
                            str(exc)
                            if str(exc).startswith("checkpoint_")
                            else "checkpoint_write_exception"
                        ),
                        "diagnostic_stage": "runtime_checkpoint_exception",
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc),
                        "non_authoritative": True,
                        "production_objects_mutated": False,
                        "refinement_applied": False,
                    }
            histDict["H_simc_shape_pi_n_SIMC"] = component_fit_result.get("H_simc_shape_pi_n")
            histDict["H_simc_shape_pi_delta_SIMC"] = component_fit_result.get("H_simc_shape_pi_delta")
            histDict["H_simc_shape_pi_sidis_SIMC"] = component_fit_result.get("H_simc_shape_pi_sidis")
            histDict["H_k_lambda_simc_reference_SIMC"] = component_fit_result.get("H_k_lambda_simc_reference")
            histDict["H_simc_shape_k_lambda_SIMC"] = component_fit_result.get("H_simc_shape_k_lambda")
            histDict["H_simc_shape_k_sigma0_SIMC"] = component_fit_result.get("H_simc_shape_k_sigma0")
            histDict["H_pion_fit_pi_n_scaled_DATA"] = component_fit_result.get("H_pion_fit_pi_n_scaled")
            histDict["H_pion_fit_pi_delta_scaled_DATA"] = component_fit_result.get("H_pion_fit_pi_delta_scaled")
            histDict["H_pion_fit_pi_sidis_scaled_DATA"] = component_fit_result.get("H_pion_fit_pi_sidis_scaled")
            histDict["H_pion_fit_total_DATA"] = component_fit_result.get("H_pion_fit_total")
            histDict["H_kaon_fit_pi_n_scaled_DATA"] = component_fit_result.get("H_kaon_fit_pi_n_scaled")
            histDict["H_kaon_fit_pi_delta_scaled_DATA"] = component_fit_result.get("H_kaon_fit_pi_delta_scaled")
            histDict["H_kaon_fit_pi_sidis_scaled_DATA"] = component_fit_result.get("H_kaon_fit_pi_sidis_scaled")
            histDict["H_kaon_fit_k_lambda_scaled_DATA"] = component_fit_result.get("H_kaon_fit_k_lambda_scaled")
            histDict["H_kaon_fit_k_sigma0_scaled_DATA"] = component_fit_result.get("H_kaon_fit_k_sigma0_scaled")
            histDict["H_kaon_fit_total_DATA"] = component_fit_result.get("H_kaon_fit_total")
            histDict["H_kaon_full_fit_residual_DATA"] = component_fit_result.get("H_kaon_full_fit_residual")
            histDict["H_kaon_pion_bg_fit_total_DATA"] = component_fit_result.get("H_kaon_pion_bg_fit_total")
            histDict["H_fit_residual_pion_DATA"] = component_fit_result.get("H_fit_residual_pion")
            histDict["H_fit_residual_kaon_DATA"] = component_fit_result.get("H_fit_residual_kaon")
            if isinstance(component_subtraction_payload, dict):
                histDict["H_pion_control_model_DATA"] = component_subtraction_payload.get("H_pion_control_model")
                histDict["H_kaon_pion_model_DATA"] = component_subtraction_payload.get("H_kaon_pion_model")
                histDict["H_pion_weight_vs_MM_DATA"] = component_subtraction_payload.get("H_pion_weight_vs_MM")
                histDict["H_pion_subtraction_template_MM_DATA"] = component_subtraction_payload.get("H_pion_subtraction_template_MM")
                histDict["H_pion_subtraction_template_MM_nosub_DATA"] = component_subtraction_payload.get("H_pion_subtraction_template_MM_nosub")
                histDict["H_MM_before_pion_subtraction_DATA"] = component_subtraction_payload.get("H_MM_before_pion_subtraction")
                histDict["H_MM_after_pion_subtraction_DATA"] = component_subtraction_payload.get("H_MM_after_pion_subtraction")

        subtraction_windows = None
        scale_components = None
        scale_factor = 0.0
        use_legacy_scalar_subtraction = True
        if component_fit_result is not None and bool(component_fit_result.get("diagnostic_only")):
            use_legacy_scalar_subtraction = False
        if isinstance(component_subtraction_payload, dict):
            if component_subtraction_payload.get("accepted"):
                use_legacy_scalar_subtraction = False
                scale_factor = float(component_subtraction_payload.get("particle_subtraction_effective_scale", 0.0) or 0.0)
            else:
                fallback_mode = component_subtraction_payload.get("fallback_mode") or "single_scale"
                if fallback_mode == "error":
                    raise RuntimeError(
                        "rand_sub component pion subtraction ({}) rejected: {}".format(
                            phi_setting,
                            component_subtraction_payload.get("fallback_reason") or "unknown reason",
                        )
                    )
                if fallback_mode == "single_scale":
                    use_legacy_scalar_subtraction = True
                elif fallback_mode in ("zero", "skip_bin"):
                    use_legacy_scalar_subtraction = False
                else:
                    use_legacy_scalar_subtraction = True

        if use_legacy_scalar_subtraction:
            try:
                subtraction_windows = resolve_particle_subtraction_windows(
                    ParticleType,
                    SubtractedParticle,
                    MM_offset_DATA,
                    inp_dict=inpDict,
                    phi_setting=phi_setting,
                )
                scale_components = compute_staged_particle_subtraction_scales(
                    H_MM_nosub_DATA,
                    subDict["H_MM_nosub_SUB_DATA"],
                    subtraction_windows,
                    context="pion subtraction ({})".format(phi_setting),
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
        histDict["particle_subtraction_scale_factor"] = scale_factor
        histDict["particle_subtraction_scale_components"] = scale_components
        histDict["particle_subtraction_windows"] = subtraction_windows
        if scale_components is not None:
            _print_rand_debug(
                "particle subtraction normalization",
                phi_setting=phi_setting,
                epsset=EPSSET,
                pi_n_window=scale_components["pi_n"]["window"],
                pi_n_data_amp=scale_components["pi_n"]["data_amp"],
                pi_n_background_amp=scale_components["pi_n"]["background_amp"],
                pi_n_scale_factor=scale_components["pi_n"]["scale_factor"],
                pi_delta_window=scale_components["pi_delta"]["window"],
                pi_delta_data_amp=scale_components["pi_delta"]["data_amp"],
                pi_delta_background_amp=scale_components["pi_delta"]["background_amp"],
                pi_delta_residual_amp=scale_components["pi_delta"]["residual_amp"],
                pi_delta_scale_factor=scale_components["pi_delta"]["scale_factor"],
                total_scale_factor=scale_factor,
            )

        if use_legacy_scalar_subtraction:
            # Apply scale factor
            subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"].Scale(scale_factor)
            subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"].Scale(scale_factor)
            subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
            subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
            subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
            subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_CoinTime_SUB_DATA"].Scale(scale_factor)
            subDict["CoinTime_vs_beta_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_beta_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_H_cer_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_H_cal_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_P_cal_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_P_hgcer_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_P_aero_SUB_DATA"].Scale(scale_factor)
            subDict["phiq_vs_t_SUB_DATA"].Scale(scale_factor)
            subDict["Q2_vs_W_SUB_DATA"].Scale(scale_factor)
            subDict["Q2_vs_t_SUB_DATA"].Scale(scale_factor)
            subDict["W_vs_t_SUB_DATA"].Scale(scale_factor)
            subDict["EPS_vs_t_SUB_DATA"].Scale(scale_factor)
            subDict["MM_vs_t_SUB_DATA"].Scale(scale_factor)
            subDict["H_ct_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssxfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssyfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssxpfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssypfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsxfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsyfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsxpfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsypfp_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssxptar_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssyptar_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsxptar_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsyptar_SUB_DATA"].Scale(scale_factor)
            subDict["H_ssdelta_SUB_DATA"].Scale(scale_factor)
            subDict["H_hsdelta_SUB_DATA"].Scale(scale_factor)
            subDict["H_ph_q_SUB_DATA"].Scale(scale_factor)
            subDict["H_th_q_SUB_DATA"].Scale(scale_factor)
            subDict["H_ph_recoil_SUB_DATA"].Scale(scale_factor)
            subDict["H_th_recoil_SUB_DATA"].Scale(scale_factor)
            subDict["H_Q2_SUB_DATA"].Scale(scale_factor)
            subDict["H_W_SUB_DATA"].Scale(scale_factor)
            subDict["H_t_SUB_DATA"].Scale(scale_factor)
            subDict["H_epsilon_SUB_DATA"].Scale(scale_factor)
            subDict["H_MM_SUB_DATA"].Scale(scale_factor)
            subDict["H_MM_nosub_SUB_DATA"].Scale(scale_factor)
            subDict["H_pmiss_SUB_DATA"].Scale(scale_factor)
            subDict["H_emiss_SUB_DATA"].Scale(scale_factor)
            subDict["H_pmx_SUB_DATA"].Scale(scale_factor)
            subDict["H_pmy_SUB_DATA"].Scale(scale_factor)
            subDict["H_pmz_SUB_DATA"].Scale(scale_factor)
            histDict["H_MM_SUB_DATA"] = subDict["H_MM_SUB_DATA"]
            histDict["H_MM_nosub_SUB_DATA"] = subDict["H_MM_nosub_SUB_DATA"]

            # Subtract pion
            P_hgcer_xAtCer_vs_yAtCer_DATA.Add(subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"],-1)
            P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"],-1)
            P_hgcer_xAtCer_vs_MM_DATA.Add(subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"],-1)
            P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"],-1)
            P_hgcer_yAtCer_vs_MM_DATA.Add(subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"],-1)
            P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"],-1)
            MM_vs_CoinTime_DATA.Add(subDict["MM_vs_CoinTime_SUB_DATA"],-1)
            CoinTime_vs_beta_DATA.Add(subDict["CoinTime_vs_beta_SUB_DATA"],-1)
            MM_vs_beta_DATA.Add(subDict["MM_vs_beta_SUB_DATA"],-1)
            MM_vs_H_cer_DATA.Add(subDict["MM_vs_H_cer_SUB_DATA"],-1)
            MM_vs_H_cal_DATA.Add(subDict["MM_vs_H_cal_SUB_DATA"],-1)
            MM_vs_P_cal_DATA.Add(subDict["MM_vs_P_cal_SUB_DATA"],-1)
            MM_vs_P_hgcer_DATA.Add(subDict["MM_vs_P_hgcer_SUB_DATA"],-1)
            MM_vs_P_aero_DATA.Add(subDict["MM_vs_P_aero_SUB_DATA"],-1)
            phiq_vs_t_DATA.Add(subDict["phiq_vs_t_SUB_DATA"],-1)
            Q2_vs_W_DATA.Add(subDict["Q2_vs_W_SUB_DATA"],-1)
            Q2_vs_t_DATA.Add(subDict["Q2_vs_t_SUB_DATA"],-1)
            W_vs_t_DATA.Add(subDict["W_vs_t_SUB_DATA"],-1)
            EPS_vs_t_DATA.Add(subDict["EPS_vs_t_SUB_DATA"],-1)
            MM_vs_t_DATA.Add(subDict["MM_vs_t_SUB_DATA"],-1)
            H_ssxfp_DATA.Add(subDict["H_ssxfp_SUB_DATA"],-1)
            H_ssyfp_DATA.Add(subDict["H_ssyfp_SUB_DATA"],-1)
            H_ssxpfp_DATA.Add(subDict["H_ssxpfp_SUB_DATA"],-1)
            H_ssypfp_DATA.Add(subDict["H_ssypfp_SUB_DATA"],-1)
            H_hsxfp_DATA.Add(subDict["H_hsxfp_SUB_DATA"],-1)
            H_hsyfp_DATA.Add(subDict["H_hsyfp_SUB_DATA"],-1)
            H_hsxpfp_DATA.Add(subDict["H_hsxpfp_SUB_DATA"],-1)
            H_hsypfp_DATA.Add(subDict["H_hsypfp_SUB_DATA"],-1)
            H_ssxptar_DATA.Add(subDict["H_ssxptar_SUB_DATA"],-1)
            H_ssyptar_DATA.Add(subDict["H_ssyptar_SUB_DATA"],-1)
            H_hsxptar_DATA.Add(subDict["H_hsxptar_SUB_DATA"],-1)
            H_hsyptar_DATA.Add(subDict["H_hsyptar_SUB_DATA"],-1)
            H_ssdelta_DATA.Add(subDict["H_ssdelta_SUB_DATA"],-1)
            H_hsdelta_DATA.Add(subDict["H_hsdelta_SUB_DATA"],-1)
            H_ph_q_DATA.Add(subDict["H_ph_q_SUB_DATA"],-1)
            H_th_q_DATA.Add(subDict["H_th_q_SUB_DATA"],-1)
            H_ph_recoil_DATA.Add(subDict["H_ph_recoil_SUB_DATA"],-1)
            H_th_recoil_DATA.Add(subDict["H_th_recoil_SUB_DATA"],-1)
            H_Q2_DATA.Add(subDict["H_Q2_SUB_DATA"],-1)
            H_W_DATA.Add(subDict["H_W_SUB_DATA"],-1)
            H_t_DATA.Add(subDict["H_t_SUB_DATA"],-1)
            H_epsilon_DATA.Add(subDict["H_epsilon_SUB_DATA"],-1)
            H_MM_fit2sub_DATA.Add(subDict["H_MM_nosub_SUB_DATA"],-1)
            H_MM_fit1sub_DATA.Add(subDict["H_MM_nosub_SUB_DATA"],-1)
            H_MM_pisub_DATA.Add(subDict["H_MM_nosub_SUB_DATA"],-1)
            H_MM_DATA.Add(subDict["H_MM_SUB_DATA"],-1)
            H_MM_full_DATA.Add(subDict["H_MM_nosub_SUB_DATA"],-1)
            H_pmiss_DATA.Add(subDict["H_pmiss_SUB_DATA"],-1)
            H_emiss_DATA.Add(subDict["H_emiss_SUB_DATA"],-1)
            H_pmx_DATA.Add(subDict["H_pmx_SUB_DATA"],-1)
            H_pmy_DATA.Add(subDict["H_pmy_SUB_DATA"],-1)
            H_pmz_DATA.Add(subDict["H_pmz_SUB_DATA"],-1)
            H_ct_DATA.Add(subDict["H_ct_SUB_DATA"],-1)
        _print_rand_timer("rand_sub pion subtraction {}".format(phi_setting), perf_counter() - stage_start)

    data_bg_targets = {
        "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DATA,
        "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DATA,
        "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DATA,
        "mm_ct": MM_vs_CoinTime_DATA,
        "ct_beta": CoinTime_vs_beta_DATA,
        "mm_beta": MM_vs_beta_DATA,
        "mm_h_cer": MM_vs_H_cer_DATA,
        "mm_h_cal": MM_vs_H_cal_DATA,
        "mm_p_cal": MM_vs_P_cal_DATA,
        "mm_p_hgcer": MM_vs_P_hgcer_DATA,
        "mm_p_aero": MM_vs_P_aero_DATA,
        "phiq_t": phiq_vs_t_DATA,
        "q2_w": Q2_vs_W_DATA,
        "q2_t": Q2_vs_t_DATA,
        "w_t": W_vs_t_DATA,
        "eps_t": EPS_vs_t_DATA,
        "mm_t": MM_vs_t_DATA,
        "h_ct": H_ct_DATA,
        "h_ssxfp": H_ssxfp_DATA,
        "h_ssyfp": H_ssyfp_DATA,
        "h_ssxpfp": H_ssxpfp_DATA,
        "h_ssypfp": H_ssypfp_DATA,
        "h_ssdelta": H_ssdelta_DATA,
        "h_ssxptar": H_ssxptar_DATA,
        "h_ssyptar": H_ssyptar_DATA,
        "h_hsxfp": H_hsxfp_DATA,
        "h_hsyfp": H_hsyfp_DATA,
        "h_hsxpfp": H_hsxpfp_DATA,
        "h_hsypfp": H_hsypfp_DATA,
        "h_hsdelta": H_hsdelta_DATA,
        "h_hsxptar": H_hsxptar_DATA,
        "h_hsyptar": H_hsyptar_DATA,
        "h_ph_q": H_ph_q_DATA,
        "h_th_q": H_th_q_DATA,
        "h_ph_recoil": H_ph_recoil_DATA,
        "h_th_recoil": H_th_recoil_DATA,
        "h_pmiss": H_pmiss_DATA,
        "h_emiss": H_emiss_DATA,
        "h_pmx": H_pmx_DATA,
        "h_pmy": H_pmy_DATA,
        "h_pmz": H_pmz_DATA,
        "h_q2": H_Q2_DATA,
        "h_t": H_t_DATA,
        "h_w": H_W_DATA,
        "h_epsilon": H_epsilon_DATA,
        "h_cal": H_cal_etottracknorm_DATA,
        "h_cer": H_cer_npeSum_DATA,
        "p_cal": P_cal_etottracknorm_DATA,
        "p_hgcer": P_hgcer_npeSum_DATA,
        "p_aero": P_aero_npeSum_DATA,
    }

    rand_debug_stage = "post-pion-subtraction"
    bg_diag1 = None
    bg_diag2 = None
    TBRANCH_SUB_DATA = sub_tree_bundle.get("prompt_tree") if isinstance(sub_tree_bundle, dict) else None
    TBRANCH_SUB_RAND = sub_tree_bundle.get("rand_tree") if isinstance(sub_tree_bundle, dict) else None
    TBRANCH_SUB_DUMMY = sub_tree_bundle.get("dummy_prompt_tree") if isinstance(sub_tree_bundle, dict) else None
    TBRANCH_SUB_DUMMY_RAND = sub_tree_bundle.get("dummy_rand_tree") if isinstance(sub_tree_bundle, dict) else None

    try:
        if ParticleType == "kaon" and TBRANCH_SUB_DATA is None:
            rand_debug_stage = "open subtracted-particle ROOT files"
            sub_tree_bundle = _open_subtracted_particle_tree_bundle(
                OUTPATH,
                phi_setting,
                SubtractedParticle,
                InDATAFilename,
                InDUMMYFilename,
                EPSSET,
            )
            TBRANCH_SUB_DATA = sub_tree_bundle.get("prompt_tree")
            TBRANCH_SUB_RAND = sub_tree_bundle.get("rand_tree")
            TBRANCH_SUB_DUMMY = sub_tree_bundle.get("dummy_prompt_tree")
            TBRANCH_SUB_DUMMY_RAND = sub_tree_bundle.get("dummy_rand_tree")

        # Fit background and subtract
        # --------------------------------------------------------------
        # Stat-scale: events that survive ALL subtractions & MM-cuts
        # --------------------------------------------------------------
        rand_debug_stage = "resolve fit 1 scale"
        inpDict["bg_stat_scale1"] = resolve_bg_stat_scale1(inpDict, phi_setting)
        residual_bg_weights1 = None
        background_fit1 = None
        background_fit2 = None
        _print_rand_debug(
            "fit 1 scale resolved",
            phi_setting=phi_setting,
            epsset=EPSSET,
            bg_stat_scale1=inpDict["bg_stat_scale1"],
            scale_factor=scale_factor if "scale_factor" in locals() else None,
        )
        active_component_payload = component_subtraction_payload if (
            isinstance(component_subtraction_payload, dict) and component_subtraction_payload.get("accepted")
        ) else None

        if inpDict["bg_stat_scale1"] > 0.0:
            rand_debug_stage = "fit 1 function build"
            background_fit1 = bg_fit(phi_setting,
                                    inpDict,
                                    H_MM_pisub_DATA,   # wide / no-cut
                                    H_MM_DATA,         # cut-window axis
                                    scaling=inpDict["bg_stat_scale1"],
                                    model_key=f"fixquad_{phi_setting}_{EPSSET}e",
                                    fit_name="Fit 1")
            mm_stage1_input = clone_reset_hist(H_MM_DATA, "_stage1_input")
            mm_stage1_input.Add(H_MM_DATA)
            bg_templates1 = _create_rand_sub_bg_templates(data_bg_targets)
            bg_weights1, bg_diag1 = build_mm_background_weights_with_diagnostics(
                mm_stage1_input,
                background_fit1[0],
            )
            _warn_if_oversub_diagnostics(
                inpDict,
                bg_diag1,
                phi_setting,
                "Fit 1",
            )

            rand_debug_stage = "fit 1 prompt data background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DATA,
                tmin,
                tmax,
                bg_templates1,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage1_input,
                bg_weights1,
                norm_factor_data,
                tree_label="prompt data fit1",
            )
            rand_debug_stage = "fit 1 random data background pass"
            _process_rand_sub_background_tree(
                TBRANCH_RAND,
                tmin,
                tmax,
                bg_templates1,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage1_input,
                bg_weights1,
                -norm_factor_data / nWindows,
                tree_label="random data fit1",
            )
            rand_debug_stage = "fit 1 prompt dummy background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DUMMY,
                tmin,
                tmax,
                bg_templates1,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage1_input,
                bg_weights1,
                -norm_factor_dummy,
                tree_label="prompt dummy fit1",
            )
            rand_debug_stage = "fit 1 random dummy background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DUMMY_RAND,
                tmin,
                tmax,
                bg_templates1,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage1_input,
                bg_weights1,
                norm_factor_dummy / nWindows,
                tree_label="random dummy fit1",
            )

            if ParticleType == "kaon" and active_component_payload is not None:
                rand_debug_stage = "fit 1 subtracted prompt data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DATA,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    -norm_factor_data,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    tree_label="subtracted prompt data fit1",
                )
                rand_debug_stage = "fit 1 subtracted random data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_RAND,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    norm_factor_data / nWindows,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    tree_label="subtracted random data fit1",
                )
                rand_debug_stage = "fit 1 subtracted prompt dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    norm_factor_dummy,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    tree_label="subtracted prompt dummy fit1",
                )
                rand_debug_stage = "fit 1 subtracted random dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY_RAND,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    -norm_factor_dummy / nWindows,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    tree_label="subtracted random dummy fit1",
                )
            elif ParticleType == "kaon" and scale_factor != 0.0:
                rand_debug_stage = "fit 1 subtracted prompt data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DATA,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    -scale_factor * norm_factor_data,
                    tree_label="subtracted prompt data fit1",
                )
                rand_debug_stage = "fit 1 subtracted random data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_RAND,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    scale_factor * norm_factor_data / nWindows,
                    tree_label="subtracted random data fit1",
                )
                rand_debug_stage = "fit 1 subtracted prompt dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    scale_factor * norm_factor_dummy,
                    tree_label="subtracted prompt dummy fit1",
                )
                rand_debug_stage = "fit 1 subtracted random dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY_RAND,
                    MM_offset_DATA,
                    bg_templates1,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage1_input,
                    bg_weights1,
                    -scale_factor * norm_factor_dummy / nWindows,
                    tree_label="subtracted random dummy fit1",
                )

            rand_debug_stage = "fit 1 histogram subtraction"
            for key, hist in data_bg_targets.items():
                hist.Add(bg_templates1[key], -1)

            residual_bg_weights1 = build_mm_residual_weights(bg_weights1)
            H_MM_fit2sub_DATA.Add(background_fit1[1], -1)
            H_MM_fit1sub_DATA.Add(background_fit1[1], -1)
            H_MM_DATA.Add(background_fit1[0], -1)
            H_MM_full_DATA.Add(background_fit1[1], -1)

        # Fit background and subtract
        # --------------------------------------------------------------
        # Stat-scale: events that survive ALL subtractions & MM-cuts
        # --------------------------------------------------------------
        rand_debug_stage = "resolve fit 2 scale"
        inpDict["bg_stat_scale2"] = resolve_bg_stat_scale2(inpDict, phi_setting)
        _print_rand_debug(
            "fit 2 scale resolved",
            phi_setting=phi_setting,
            epsset=EPSSET,
            bg_stat_scale2=inpDict["bg_stat_scale2"],
        )

        if inpDict["bg_stat_scale2"] > 0.0:
            rand_debug_stage = "fit 2 function build"
            background_fit2 = bg_fit(phi_setting,
                                    inpDict,
                                    H_MM_fit1sub_DATA,   # wide / no-cut
                                    H_MM_DATA,         # cut-window axis
                                    scaling=inpDict["bg_stat_scale2"],
                                    model_key=f"cheb2_{phi_setting}_{EPSSET}e",
                                    fit_name="Fit 2")
            mm_stage2_input = clone_reset_hist(H_MM_DATA, "_stage2_input")
            mm_stage2_input.Add(H_MM_DATA)
            bg_templates2 = _create_rand_sub_bg_templates(data_bg_targets)
            bg_weights2, bg_diag2 = build_mm_background_weights_with_diagnostics(
                mm_stage2_input,
                background_fit2[0],
            )
            _warn_if_oversub_diagnostics(
                inpDict,
                bg_diag2,
                phi_setting,
                "Fit 2",
            )

            rand_debug_stage = "fit 2 prompt data background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DATA,
                tmin,
                tmax,
                bg_templates2,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage2_input,
                bg_weights2,
                norm_factor_data,
                residual_weights=residual_bg_weights1,
                tree_label="prompt data fit2",
            )
            rand_debug_stage = "fit 2 random data background pass"
            _process_rand_sub_background_tree(
                TBRANCH_RAND,
                tmin,
                tmax,
                bg_templates2,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage2_input,
                bg_weights2,
                -norm_factor_data / nWindows,
                residual_weights=residual_bg_weights1,
                tree_label="random data fit2",
            )
            rand_debug_stage = "fit 2 prompt dummy background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DUMMY,
                tmin,
                tmax,
                bg_templates2,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage2_input,
                bg_weights2,
                -norm_factor_dummy,
                residual_weights=residual_bg_weights1,
                tree_label="prompt dummy fit2",
            )
            rand_debug_stage = "fit 2 random dummy background pass"
            _process_rand_sub_background_tree(
                TBRANCH_DUMMY_RAND,
                tmin,
                tmax,
                bg_templates2,
                ParticleType,
                hole_contains,
                evaluate_data_event,
                mm_min,
                mm_max,
                mm_stage2_input,
                bg_weights2,
                norm_factor_dummy / nWindows,
                residual_weights=residual_bg_weights1,
                tree_label="random dummy fit2",
            )

            if ParticleType == "kaon" and active_component_payload is not None:
                rand_debug_stage = "fit 2 subtracted prompt data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DATA,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    -norm_factor_data,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted prompt data fit2",
                )
                rand_debug_stage = "fit 2 subtracted random data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_RAND,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    norm_factor_data / nWindows,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted random data fit2",
                )
                rand_debug_stage = "fit 2 subtracted prompt dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    norm_factor_dummy,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted prompt dummy fit2",
                )
                rand_debug_stage = "fit 2 subtracted random dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY_RAND,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    -norm_factor_dummy / nWindows,
                    pion_reference_hist=active_component_payload["H_pion_control_model"],
                    pion_mm_weights=active_component_payload["weights"],
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted random dummy fit2",
                )
            elif ParticleType == "kaon" and scale_factor != 0.0:
                rand_debug_stage = "fit 2 subtracted prompt data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DATA,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    -scale_factor * norm_factor_data,
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted prompt data fit2",
                )
                rand_debug_stage = "fit 2 subtracted random data background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_RAND,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    scale_factor * norm_factor_data / nWindows,
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted random data fit2",
                )
                rand_debug_stage = "fit 2 subtracted prompt dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    scale_factor * norm_factor_dummy,
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted prompt dummy fit2",
                )
                rand_debug_stage = "fit 2 subtracted random dummy background pass"
                _process_subtracted_particle_background_tree(
                    TBRANCH_SUB_DUMMY_RAND,
                    MM_offset_DATA,
                    bg_templates2,
                    ParticleType,
                    hole_contains,
                    evaluate_data_event,
                    get_shifted_t,
                    mm_min,
                    mm_max,
                    mm_stage2_input,
                    bg_weights2,
                    -scale_factor * norm_factor_dummy / nWindows,
                    residual_weights=residual_bg_weights1,
                    tree_label="subtracted random dummy fit2",
                )

            rand_debug_stage = "fit 2 histogram subtraction"
            for key, hist in data_bg_targets.items():
                hist.Add(bg_templates2[key], -1)

            H_MM_fit2sub_DATA.Add(background_fit2[1], -1)
            H_MM_DATA.Add(background_fit2[0], -1)
            H_MM_full_DATA.Add(background_fit2[1], -1)
    except Exception:
        _print_rand_debug(
            "failure after pion subtraction",
            stage=rand_debug_stage,
            phi_setting=phi_setting,
            epsset=EPSSET,
            particle_type=ParticleType,
            bg_stat_scale1=inpDict.get("bg_stat_scale1"),
            bg_stat_scale2=inpDict.get("bg_stat_scale2"),
            scale_factor=scale_factor if "scale_factor" in locals() else None,
        )
        traceback.print_exc()
        raise

    stage_start = perf_counter()
    histDict["InFile_DATA"] = InFile_DATA
    histDict["InFile_DUMMY"] = InFile_DUMMY
    histDict["phi_setting"] = phi_setting
    histDict["pid_text"] = pid_text
    histDict["runNums"] = runNums.split(' ')
    histDict["nWindows"] = nWindows
    histDict["H_hsdelta_DUMMY"] =     H_hsdelta_DUMMY
    histDict["H_hsxptar_DUMMY"] =     H_hsxptar_DUMMY
    histDict["H_hsyptar_DUMMY"] =     H_hsyptar_DUMMY
    histDict["H_ssxfp_DUMMY"] =     H_ssxfp_DUMMY
    histDict["H_ssyfp_DUMMY"] =     H_ssyfp_DUMMY
    histDict["H_ssxpfp_DUMMY"] =     H_ssxpfp_DUMMY
    histDict["H_ssypfp_DUMMY"] =     H_ssypfp_DUMMY
    histDict["H_hsxfp_DUMMY"] =     H_hsxfp_DUMMY
    histDict["H_hsyfp_DUMMY"] =     H_hsyfp_DUMMY
    histDict["H_hsxpfp_DUMMY"] =     H_hsxpfp_DUMMY
    histDict["H_hsypfp_DUMMY"] =     H_hsypfp_DUMMY
    histDict["H_ssdelta_DUMMY"] =     H_ssdelta_DUMMY
    histDict["H_ssxptar_DUMMY"] =     H_ssxptar_DUMMY
    histDict["H_ssyptar_DUMMY"] =     H_ssyptar_DUMMY
    histDict["H_q_DUMMY"] =     H_q_DUMMY
    histDict["H_Q2_DUMMY"] =     H_Q2_DUMMY
    histDict["H_t_DUMMY"] =     H_t_DUMMY
    histDict["H_epsilon_DUMMY"] =     H_epsilon_DUMMY
    histDict["H_MM_DUMMY"] =     H_MM_DUMMY
    histDict["H_MM_full_DUMMY"] =     H_MM_full_DUMMY
    histDict["H_MM_fit2sub_DUMMY"] =     H_MM_fit2sub_DUMMY
    histDict["H_MM_fit1sub_DUMMY"] =     H_MM_fit1sub_DUMMY
    histDict["H_MM_pisub_DUMMY"] =     H_MM_pisub_DUMMY
    histDict["H_MM_nosub_DUMMY"] =     H_MM_nosub_DUMMY
    histDict["H_th_DUMMY"] =     H_th_DUMMY
    histDict["H_ph_DUMMY"] =     H_ph_DUMMY
    histDict["H_ph_q_DUMMY"] =     H_ph_q_DUMMY
    histDict["H_th_q_DUMMY"] =     H_th_q_DUMMY
    histDict["H_ph_recoil_DUMMY"] =     H_ph_recoil_DUMMY
    histDict["H_th_recoil_DUMMY"] =     H_th_recoil_DUMMY
    histDict["H_pmiss_DUMMY"] =     H_pmiss_DUMMY
    histDict["H_emiss_DUMMY"] =     H_emiss_DUMMY
    histDict["H_pmx_DUMMY"] =     H_pmx_DUMMY
    histDict["H_pmy_DUMMY"] =     H_pmy_DUMMY
    histDict["H_pmz_DUMMY"] =     H_pmz_DUMMY
    histDict["H_W_DUMMY"] =     H_W_DUMMY
    histDict["H_ct_DUMMY"] =     H_ct_DUMMY
    histDict["MM_vs_CoinTime_DUMMY"] = MM_vs_CoinTime_DUMMY
    histDict["CoinTime_vs_beta_DUMMY"] = CoinTime_vs_beta_DUMMY
    histDict["MM_vs_beta_DUMMY"] = MM_vs_beta_DUMMY
    histDict["phiq_vs_t_DUMMY"] = phiq_vs_t_DUMMY
    histDict["polar_phiq_vs_t_DUMMY"] = polar_phiq_vs_t_DUMMY
    histDict["H_hsdelta_DATA"] =     H_hsdelta_DATA
    histDict["H_hsxptar_DATA"] =     H_hsxptar_DATA
    histDict["H_hsyptar_DATA"] =     H_hsyptar_DATA
    histDict["H_ssxfp_DATA"] =     H_ssxfp_DATA
    histDict["H_ssyfp_DATA"] =     H_ssyfp_DATA
    histDict["H_ssxpfp_DATA"] =     H_ssxpfp_DATA
    histDict["H_ssypfp_DATA"] =     H_ssypfp_DATA
    histDict["H_hsxfp_DATA"] =     H_hsxfp_DATA
    histDict["H_hsyfp_DATA"] =     H_hsyfp_DATA
    histDict["H_hsxpfp_DATA"] =     H_hsxpfp_DATA
    histDict["H_hsypfp_DATA"] =     H_hsypfp_DATA
    histDict["H_ssdelta_DATA"] =     H_ssdelta_DATA
    histDict["H_ssxptar_DATA"] =     H_ssxptar_DATA
    histDict["H_ssyptar_DATA"] =     H_ssyptar_DATA
    histDict["H_q_DATA"] =     H_q_DATA
    histDict["H_Q2_DATA"] =     H_Q2_DATA
    histDict["H_t_DATA"] =     H_t_DATA
    histDict["_binning_H_t_DATA_pre_particle_subtraction"] = binning_t_hist
    histDict["H_epsilon_DATA"] =     H_epsilon_DATA
    histDict["H_MM_DATA"] =     H_MM_DATA
    histDict["H_MM_rand_dummy_DATA"] =     H_MM_rand_dummy_DATA
    histDict["H_MM_dummy_DATA"] =     H_MM_dummy_DATA
    histDict["H_MM_full_DATA"] =     H_MM_full_DATA        
    histDict["H_MM_fit2sub_DATA"] =     H_MM_fit2sub_DATA
    histDict["H_MM_fit1sub_DATA"] =     H_MM_fit1sub_DATA
    histDict["H_MM_pisub_DATA"] =     H_MM_pisub_DATA
    histDict["H_MM_nosub_DATA"] =     H_MM_nosub_DATA
    histDict["BG_FIT1_VIS_DATA"] = background_fit1[1] if background_fit1 is not None else None
    histDict["BG_FIT2_VIS_DATA"] = background_fit2[1] if background_fit2 is not None else None
    histDict["bg_oversub_diagnostics"] = {
        "fit1": bg_diag1 if "bg_diag1" in locals() else None,
        "fit2": bg_diag2 if "bg_diag2" in locals() else None,
    }
    if "H_MM_SUB_DATA" not in histDict:
        if isinstance(component_subtraction_payload, dict) and component_subtraction_payload.get("H_pion_subtraction_template_MM") is not None:
            histDict["H_MM_SUB_DATA"] = component_subtraction_payload.get("H_pion_subtraction_template_MM")
        else:
            histDict["H_MM_SUB_DATA"] = clone_reset_hist(subDict["H_MM_SUB_DATA"], "_component_empty")
    if "H_MM_nosub_SUB_DATA" not in histDict:
        if isinstance(component_subtraction_payload, dict) and component_subtraction_payload.get("H_pion_subtraction_template_MM_nosub") is not None:
            histDict["H_MM_nosub_SUB_DATA"] = component_subtraction_payload.get("H_pion_subtraction_template_MM_nosub")
        else:
            histDict["H_MM_nosub_SUB_DATA"] = clone_reset_hist(subDict["H_MM_nosub_SUB_DATA"], "_component_empty")
    if "particle_subtraction_scale_factor" not in histDict:
        histDict["particle_subtraction_scale_factor"] = 0.0
    if "particle_subtraction_scale_components" not in histDict:
        histDict["particle_subtraction_scale_components"] = None
    histDict["H_th_DATA"] =     H_th_DATA
    histDict["H_ph_DATA"] =     H_ph_DATA
    histDict["H_ph_q_DATA"] =     H_ph_q_DATA
    histDict["_binning_H_ph_q_DATA_pre_particle_subtraction"] = binning_phi_hist
    histDict["_binning_input_stage"] = "post_random_dummy_pre_particle_subtraction"
    histDict["H_th_q_DATA"] =     H_th_q_DATA
    histDict["H_ph_recoil_DATA"] =     H_ph_recoil_DATA
    histDict["H_th_recoil_DATA"] =     H_th_recoil_DATA
    histDict["H_pmiss_DATA"] =     H_pmiss_DATA
    histDict["H_emiss_DATA"] =     H_emiss_DATA
    histDict["H_pmx_DATA"] =     H_pmx_DATA
    histDict["H_pmy_DATA"] =     H_pmy_DATA
    histDict["H_pmz_DATA"] =     H_pmz_DATA
    histDict["H_W_DATA"] =     H_W_DATA
    histDict["H_ct_DATA"] =     H_ct_DATA
    histDict["H_cal_etottracknorm_DATA"] =     H_cal_etottracknorm_DATA
    histDict["H_cer_npeSum_DATA"] =     H_cer_npeSum_DATA
    histDict["P_cal_etottracknorm_DATA"] =     P_cal_etottracknorm_DATA
    histDict["P_hgcer_npeSum_DATA"] =     P_hgcer_npeSum_DATA
    histDict["P_aero_npeSum_DATA"] =     P_aero_npeSum_DATA
    histDict["MM_vs_CoinTime_DATA"] = MM_vs_CoinTime_DATA
    histDict["CoinTime_vs_beta_DATA"] = CoinTime_vs_beta_DATA
    histDict["MM_vs_beta_DATA"] = MM_vs_beta_DATA
    histDict["phiq_vs_t_DATA"] = phiq_vs_t_DATA
    histDict["polar_phiq_vs_t_DATA"] = polar_phiq_vs_t_DATA
    histDict["Q2_vs_W_DATA"] = Q2_vs_W_DATA
    histDict["Q2_vs_t_DATA"] = Q2_vs_t_DATA
    histDict["W_vs_t_DATA"] = W_vs_t_DATA
    histDict["EPS_vs_t_DATA"] = EPS_vs_t_DATA
    histDict["MM_vs_t_DATA"] = MM_vs_t_DATA
    histDict["MM_vs_H_cer_DATA"] = MM_vs_H_cer_DATA
    histDict["MM_vs_H_cal_DATA"] = MM_vs_H_cal_DATA
    histDict["MM_vs_P_cal_DATA"] = MM_vs_P_cal_DATA
    histDict["MM_vs_P_hgcer_DATA"] = MM_vs_P_hgcer_DATA
    histDict["MM_vs_P_aero_DATA"] = MM_vs_P_aero_DATA
    histDict["NumEvts_MM_DUMMY"] = H_MM_DUMMY.Integral()
    histDict["NumEvts_MM_DATA"] = H_MM_DATA.Integral()
    _print_rand_timer("rand_sub histDict pack {}".format(phi_setting), perf_counter() - stage_start)

    if not emit_plots:
        _print_rand_timer("rand_sub total {}".format(phi_setting), perf_counter() - total_start)
        return histDict

    # C.4 only reroutes presentation of the existing post_proton_pre_rf
    # products; their production stage, coordinates, and event factors remain
    # frozen before this terminal PDF block.
    
    ###
    # CT plots
    stage_start = perf_counter()
    ct = TCanvas()
    l_ct = TLegend(0.115,0.65,0.33,0.95)
    l_ct.SetTextSize(0.0235)
    H_ct_DATA.SetLineColor(1)
    H_ct_DATA.Draw("same, HIST")
    l_ct.AddEntry(H_ct_DATA,"{}".format(ParticleType.capitalize()))
    l_ct.Draw()

    main_pdf = outputpdf.replace(
        "{}_FullAnalysis_".format(ParticleType),
        "{}_{}_rand_sub_".format(phi_setting, ParticleType),
    )
    ct.Print(main_pdf + '(')
    pdf_destinations = build_pdf_destinations(main_pdf)
    pdf_route_manifest = build_pdf_route_manifest(main_pdf)
    supplement_manifests = {}
    # Detached main-PDF render failures are retained for terminal QA; they
    # must never interrupt the event-analysis or PDF-close paths.
    setting_renderer_failures = []
    checkpoint_for_plots = (
        pion_hgcer_refinement_checkpoint
        if isinstance(pion_hgcer_refinement_checkpoint, dict)
        else histDict.get("pion_hgcer_refinement_checkpoint")
    )
    phase_d_checkpoint_for_plots = (
        pion_hgcer_phase_d_checkpoint
        if isinstance(pion_hgcer_phase_d_checkpoint, dict)
        else histDict.get("pion_hgcer_phase_d_checkpoint")
    )
    # C.4.1 presentation context is copied from already-established runtime
    # records.  It deliberately performs no gate, closure, or host-state work.
    proton_display_diagnostics = (
        proton_cleaning_result.get("diagnostics")
        if isinstance(proton_cleaning_result, dict)
        else {}
    ) or {}
    proton_display_gate = proton_display_diagnostics.get("lambda_preservation_gate") or {}
    phase_a_display_context = {
        "lambda_gate_status": proton_display_gate.get(
            "status",
            (pion_hgcer_event_contract or {}).get("lambda_gate_status")
            if isinstance(pion_hgcer_event_contract, dict) else "not recorded",
        ),
        "production_action": proton_display_gate.get(
            "production_action",
            (pion_hgcer_event_contract or {}).get("lambda_gate_production_action")
            if isinstance(pion_hgcer_event_contract, dict) else "not recorded",
        ),
        "proton_cleaning_committed": proton_display_gate.get(
            "proton_cleaning_committed",
            (proton_cleaning_application or {}).get("accepted", "not recorded")
            if isinstance(proton_cleaning_application, dict) else "not recorded",
        ),
        "host_state": (
            (proton_cleaning_application or {}).get("host_state")
            if isinstance(proton_cleaning_application, dict) else None
        ) or (
            (pion_hgcer_event_contract or {}).get("host_state")
            if isinstance(pion_hgcer_event_contract, dict) else "not recorded"
        ),
    }
    if ParticleType == "kaon":
        for supplement_key, role in (
            ("proton_debug", "proton-debug"),
            ("pion_fit_debug", "pion-fit-debug"),
            ("hgcer_debug", "hgcer-debug"),
        ):
            supplement_manifests[supplement_key] = open_diagnostic_pdf(
                pdf_destinations[supplement_key],
                checkpoint_for_plots,
                role=role,
                main_pdf=main_pdf,
            )
        if isinstance(phase_d_checkpoint_for_plots, dict):
            supplement_manifests["phase_d_ab"] = open_diagnostic_pdf(
                pdf_destinations["phase_d_ab"],
                phase_d_checkpoint_for_plots,
                role="hgcer-ab-comparison",
                main_pdf=main_pdf,
            )
    histDict["pdf_cleanup_route_manifest"] = pdf_route_manifest

    ###
    # Q2 plots    
    CQ2 = TCanvas()

    histDict["H_Q2_DATA"].SetLineColor(1)
    histDict["H_Q2_DATA"].Draw("same, E1")

    CQ2.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # W plots    
    CW = TCanvas()

    histDict["H_W_DATA"].SetLineColor(1)
    histDict["H_W_DATA"].Draw("same, E1")

    CW.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))    
    
    ###
    # MM plots    
    CMM = TCanvas()

    histDict["H_MM_DATA"].SetLineColor(1)
    histDict["H_MM_DATA"].Draw("same, E1")

    CMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # MM full plots    

    gStyle.SetOptStat(0)

    def style_hist(h, line_col, fill_col=None, alpha=0.25, lstyle=1, lwidth=2, fstyle=1001):
        h.SetLineColor(line_col)
        h.SetLineStyle(lstyle)
        h.SetLineWidth(lwidth)

        if fill_col is None:
            h.SetFillStyle(0)  # no fill
        else:
            if hasattr(h, "SetFillColorAlpha"):
                h.SetFillColorAlpha(fill_col, alpha)
            else:
                h.SetFillColor(fill_col)  # fallback (no true alpha)
            h.SetFillStyle(fstyle)

    # --- Canvas
    CMMfull = TCanvas("CMMfull", "MM full (subtraction steps)", 1000, 700)
    CMMfull.SetGrid()

    # --- Ordered steps (earliest -> latest), increasing visual priority
    steps = [
        ("H_MM_rand_dummy_DATA", "rand_dummy", kGray+2,   kGray+2,   0.12, 1, 2),
        ("H_MM_dummy_DATA",      "dummy",      kOrange+7, kOrange+7, 0.14, 1, 2),
        ("H_MM_pisub_DATA",      "pi-sub",     kGreen+2,  kGreen+2,  0.16, 1, 2),
        ("H_MM_fit1sub_DATA",    "fit1-sub",   kAzure+2,  kAzure+2,  0.18, 2, 2),
        ("H_MM_fit2sub_DATA",    "fit2-sub",   kViolet+1, kViolet+1, 0.20, 3, 2),
        ("H_MM_full_DATA",       "full",       kBlack,    None,      0.00, 1, 3),  # final: bold line, no fill
    ]

    hlist = [histDict[k] for k, *_ in steps]
    ymax = max(h.GetMaximum() for h in hlist)
    hlist[0].SetMaximum(1.15 * ymax)
    hlist[0].SetMinimum(0)

    # --- Draw + legend
    leg = TLegend(0.62, 0.62, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for i, (key, label, lcol, fcol, a, ls, lw) in enumerate(steps):
        h = histDict[key]
        style_hist(h, lcol, fcol, alpha=a, lstyle=ls, lwidth=lw, fstyle=1001)

        opt = "hist" if i == 0 else "hist same"
        h.Draw(opt)
        leg.AddEntry(h, label, "lf" if fcol is not None else "l")

    leg.Draw()

    # --- Cut lines (blue, dotted)
    gPad.Update()
    ymin, ymax = gPad.GetUymin(), gPad.GetUymax()
    x1, x2 = float(inpDict["mm_min"]), float(inpDict["mm_max"])

    line_min = TLine(x1, ymin, x1, ymax); line_min.SetLineColor(kBlue); line_min.SetLineStyle(3); line_min.SetLineWidth(2); line_min.Draw("same")
    line_max = TLine(x2, ymin, x2, ymax); line_max.SetLineColor(kBlue); line_max.SetLineStyle(3); line_max.SetLineWidth(2); line_max.Draw("same")

    gPad.Modified(); gPad.Update()

    CMMfull.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))    
    
    ###
    # MM sub plots    
    CMMfit2sub = TCanvas()

    histDict["H_MM_fit2sub_DATA"].SetLineColor(1)
    histDict["H_MM_fit2sub_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_fit2sub_DATA"].SetFillColor(kBlack)  # Set fill color to black
    #histDict["H_MM_fit2sub_DATA"].Draw("same, E1")
    histDict["H_MM_fit2sub_DATA"].Draw("hist same")

    # C.4: individual intermediate subtraction snapshots remain computed but
    # are intentionally not emitted to the reader-facing main PDF.

    ###
    # MM sub plots    
    CMMfit1sub = TCanvas()

    histDict["H_MM_fit1sub_DATA"].SetLineColor(1)
    histDict["H_MM_fit1sub_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_fit1sub_DATA"].SetFillColor(kBlack)  # Set fill color to black
    #histDict["H_MM_fit1sub_DATA"].Draw("same, E1")
    histDict["H_MM_fit1sub_DATA"].Draw("hist same")
    if inpDict["bg_stat_scale2"] > 0.0:
        background_fit2[1].SetLineColor(3)
        background_fit2[1].Draw("same")

    # C.4: intermediate fit-1 snapshot suppressed from PDF output.

    ###
    # MM sub plots    
    CMMpisub = TCanvas()

    histDict["H_MM_pisub_DATA"].SetLineColor(1)
    histDict["H_MM_pisub_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_pisub_DATA"].SetFillColor(kBlack)  # Set fill color to black
    #histDict["H_MM_pisub_DATA"].Draw("same, E1")
    histDict["H_MM_pisub_DATA"].Draw("hist same")    
    if inpDict["bg_stat_scale1"] > 0.0:
        background_fit1[1].SetLineColor(3)
        background_fit1[1].Draw("same")

    # C.4: intermediate pion-subtraction snapshot suppressed from PDF output.

    ###
    # MM sub plots    
    CMMsub = TCanvas()

    histDict["H_MM_nosub_DATA"].SetLineColor(1)
    histDict["H_MM_nosub_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_nosub_DATA"].SetFillColor(kBlack)  # Set fill color to black
    #histDict["H_MM_nosub_DATA"].Draw("same, E1")
    histDict["H_MM_nosub_DATA"].Draw("hist same")
    if ParticleType == "kaon":        
        histDict["H_MM_nosub_SUB_DATA"].SetLineColor(2)
        histDict["H_MM_nosub_SUB_DATA"].Draw("same, E1")

    # C.4: intermediate no-subtraction snapshot suppressed from PDF output.

    if isinstance(proton_cleaning_result, dict):
        if str(proton_cleaning_result.get("method") or "") == "timing_t_event_weight":
            display_audit = audit_timing_t_hgcer_display_targets(
                proton_cleaning_result,
                {
                    "hgcer_xy": P_hgcer_xAtCer_vs_yAtCer_DATA,
                    "hgcer_xy_nohole": P_hgcer_nohole_xAtCer_vs_yAtCer_DATA,
                    "hgcer_x_mm": P_hgcer_xAtCer_vs_MM_DATA,
                    "hgcer_x_mm_nohole": P_hgcer_nohole_xAtCer_vs_MM_DATA,
                    "hgcer_y_mm": P_hgcer_yAtCer_vs_MM_DATA,
                    "hgcer_y_mm_nohole": P_hgcer_nohole_yAtCer_vs_MM_DATA,
                },
            )
            histDict["proton_contamination_cleaning_result_setting"] = (
                serialize_kaon_proton_cleaning_result(proton_cleaning_result)
            )
            histDict["proton_contamination_cleaning_setting"] = (
                summarize_kaon_proton_cleaning_result(proton_cleaning_result)
            )
            refreshed_artifacts = _write_timing_t_validation_artifacts(
                proton_cleaning_result,
                outpath=OUTPATH,
                particle_type=ParticleType,
                outfilename=OutFilename,
                epsset=EPSSET,
                phi_setting=phi_setting,
            )
            histDict["proton_contamination_cleaning_artifacts"] = list(dict.fromkeys(
                list(histDict.get("proton_contamination_cleaning_artifacts") or [])
                + list(refreshed_artifacts or [])
            ))
            _print_rand_debug(
                "timing-t HGCer final display audit",
                populated=sum(
                    int((entry or {}).get("nonzero_bin_count", 0) or 0) > 0
                    for entry in (display_audit.get("final_display_histograms") or {}).values()
                ),
            )
        try:
            render_proton_main_summary_pages(
                main_pdf,
                checkpoint_for_plots,
                proton_cleaning_result,
                proton_cleaning_application,
                phase_a_display_context=phase_a_display_context,
                page_manifest=histDict.setdefault("proton_cleaning_main_page_manifest", []),
            )
        except Exception as exc:
            setting_renderer_failures.append(
                "proton main summary: {}: {}".format(type(exc).__name__, exc)
            )
            _print_rand_debug(
                "detached proton main PDF pages unavailable",
                renderer="pion_hgcer_refinement_plots",
                exception_type=type(exc).__name__,
                exception=str(exc),
            )
        print_kaon_proton_cleaning_pages(
            pdf_destinations["proton_debug"],
            proton_cleaning_result,
            title_prefix="{} {}".format(phi_setting, ParticleType),
        )

    component_page_manifest = histDict.setdefault("pion_component_page_manifest", [])
    setting_wide_pages_enabled = bool(component_fit_result is not None) and (
        not bool(component_fit_result.get("diagnostic_only"))
        or bool(inpDict.get("emit_setting_wide_pion_diagnostic", True))
    )
    setting_wide_render_payload = (
        component_subtraction_payload
        if isinstance(component_subtraction_payload, dict)
        else component_diagnostic_payload
    )
    setting_wide_title = "{} {} SETTING-WIDE DIAGNOSTIC FIT - NON-AUTHORITATIVE".format(
        phi_setting, ParticleType
    )

    if setting_wide_pages_enabled and component_payload is not None:
        print_particle_subtraction_component_template_pages(
            pdf_destinations["pion_fit_debug"],
            component_payload,
            title_prefix=setting_wide_title,
            cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
            kaon_signal_payload=kaon_signal_shape_payload,
            kaon_sigma0_payload=kaon_sigma0_shape_payload,
            page_manifest=component_page_manifest,
            page_id_prefix="pion.setting_wide",
            authoritative=not bool(component_fit_result.get("diagnostic_only")),
        )

    if setting_wide_pages_enabled:
        print_particle_subtraction_component_fit_pages(
            pdf_destinations["pion_fit_debug"],
            component_fit_result,
            title_prefix=setting_wide_title,
            cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
            page_manifest=component_page_manifest,
            page_id_prefix="pion.setting_wide",
            authoritative=not bool(component_fit_result.get("diagnostic_only")),
        )
    if setting_wide_pages_enabled and isinstance(setting_wide_render_payload, dict):
        print_particle_subtraction_component_application_pages(
            pdf_destinations["pion_fit_debug"],
            setting_wide_render_payload,
            title_prefix=setting_wide_title,
            cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
            component_fit_result=component_fit_result,
            include_lambda_page=False,
            page_manifest=component_page_manifest,
            page_id_prefix="pion.setting_wide",
            authoritative=not bool(component_fit_result.get("diagnostic_only")),
        )
    if setting_wide_pages_enabled:
        print_particle_subtraction_kaon_lambda_comparison_page(
            pdf_destinations["pion_fit_debug"],
            component_fit_result,
            setting_wide_render_payload,
            title_prefix=setting_wide_title,
            cut_window=(float(inpDict["mm_min"]), float(inpDict["mm_max"])),
            page_manifest=component_page_manifest,
            page_id_prefix="pion.setting_wide",
            authoritative=not bool(component_fit_result.get("diagnostic_only")),
        )
    t_bin_parent_results = histDict.get("_pion_t_bin_parent_results") or []
    pion_parent_plot_contract = {}
    if t_bin_parent_results:
        pion_parent_plot_contract = render_setting_t_bin_pion_parent_pages(
            main_pdf,
            t_bin_parent_results,
            inpDict,
            title_prefix="{} {}".format(phi_setting, ParticleType),
            page_manifest=component_page_manifest,
            setting_wide_summary=(histDict.get("pion_t_parent_diagnostics") or {}).get("setting_wide"),
            canonical_t_global=histDict.get("_pion_authoritative_canonical_t_global"),
            setting_wide_enabled=setting_wide_pages_enabled,
            coordinate_audit=histDict.get("pion_control_source_audit"),
            coordinate_diagnostics=(
                (histDict.get("_authoritative_pion_control_source_cache") or {}).get(
                    "coordinate_diagnostics"
                )
            ),
            coordinate_debug_pdf=pdf_destinations["pion_fit_debug"],
        )
        histDict["pion_component_plot_contract"] = pion_parent_plot_contract
    # C.4.3 consumes only the scalar result of the actual canonical-parent
    # K-Lambda renderer.  Page-manifest presence and pre-render fit provenance
    # cannot distinguish a comparison from its unavailable-status fallback.
    canonical_parent_k_lambda_qa = list(
        (pion_parent_plot_contract or {}).get(
            "canonical_parent_k_lambda_render", []
        )
    )
    histDict["canonical_parent_k_lambda_qa"] = canonical_parent_k_lambda_qa

    ###
    # MM dummy plots    
    CMMdummy = TCanvas()

    histDict["H_MM_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_dummy_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_dummy_DATA"].SetFillColor(kBlack)  # Set fill color to black
    histDict["H_MM_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_dummy_DATA"].Draw("hist same")

    CMMdummy.Print(pdf_destinations["pion_fit_debug"])

    ###
    # MM rand dummy plots    
    CMMranddummy = TCanvas()

    histDict["H_MM_rand_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_rand_dummy_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_rand_dummy_DATA"].SetFillColor(kBlack)  # Set fill color to black    
    histDict["H_MM_rand_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_rand_dummy_DATA"].Draw("hist same")

    CMMranddummy.Print(pdf_destinations["pion_fit_debug"])
    
    ###
    # t-Phi plots        
    Cpht_data = TCanvas()

    # Create the polar plot using the function
    polar_plot = create_polar_plot(histDict["polar_phiq_vs_t_DATA"])
    # Draw the plot
    polar_plot.Draw("AP")
    
    Cpht_data.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # t plots            
    Ct = TCanvas()
    l_t = TLegend(0.115,0.45,0.33,0.95)
    l_t.SetTextSize(0.0135)

    histDict["H_t_DATA"].SetLineColor(1)
    l_t.AddEntry(histDict["H_t_DATA"],histDict["phi_setting"])
    histDict["H_t_DATA"].Draw("same, E1")

    Ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # phi plots            
    Cphi = TCanvas()
    l_phi = TLegend(0.115,0.45,0.33,0.95)
    l_phi.SetTextSize(0.0135)
    histDict["H_ph_q_DATA"].SetLineColor(1)
    l_phi.AddEntry(histDict["H_ph_q_DATA"],histDict["phi_setting"])
    histDict["H_ph_q_DATA"].Draw("same, E1")    

    Cphi.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    
    ###
    # PID Plots
    c_pid = TCanvas()

    c_pid.Divide(2,3)

    c_pid.cd(1)
    gPad.SetLogy()

    H_cal_etottracknorm_DATA.SetLineColor(1)
    H_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(2)
    gPad.SetLogy()

    H_cer_npeSum_DATA.SetLineColor(1)
    H_cer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(3)
    gPad.SetLogy()

    P_cal_etottracknorm_DATA.SetLineColor(1)
    P_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(4)
    gPad.SetLogy()

    P_hgcer_npeSum_DATA.SetLineColor(1)
    P_hgcer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(5)
    gPad.SetLogy()

    P_aero_npeSum_DATA.SetLineColor(1)
    P_aero_npeSum_DATA.Draw("same, HIST")

    c_pid.Draw()

    if ParticleType == "kaon":
        c_pid.Print(pdf_destinations["hgcer_debug"])
    else:
        c_pid.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')

    if ParticleType == "kaon":
        
        ##
        # HGCer Hole Plots
        c_hgcervsMM = TCanvas()

        c_hgcervsMM.Divide(2,2)

        hgc_display_audit = (
            ((proton_cleaning_result or {}).get("diagnostics") or {})
            .get("generic_hgcer_diagnostic_integrity", {})
            .get("final_display_histograms", {})
            if isinstance(proton_cleaning_result, dict)
            and str(proton_cleaning_result.get("method") or "") == "timing_t_event_weight"
            else {}
        )
        if hgc_display_audit:
            for pad_index, key, hist, title in (
                (1, "hgcer_x_mm", P_hgcer_xAtCer_vs_MM_DATA, "HGCer X versus shifted missing mass;HGCer xAtCer;adj_MM [GeV]"),
                (2, "hgcer_x_mm_nohole", P_hgcer_nohole_xAtCer_vs_MM_DATA, "HGCer X versus shifted missing mass (no hole);HGCer xAtCer;adj_MM [GeV]"),
                (3, "hgcer_y_mm", P_hgcer_yAtCer_vs_MM_DATA, "HGCer Y versus shifted missing mass;HGCer yAtCer;adj_MM [GeV]"),
                (4, "hgcer_y_mm_nohole", P_hgcer_nohole_yAtCer_vs_MM_DATA, "HGCer Y versus shifted missing mass (no hole);HGCer yAtCer;adj_MM [GeV]"),
            ):
                c_hgcervsMM.cd(pad_index)
                _draw_hgcer_signed_display(hist, title, hgc_display_audit.get(key))
        else:
            # Legacy C-macro-reference diagnostic rendering is intentionally
            # unchanged; adaptive signed display applies only to timing-t.
            c_hgcervsMM.cd(1)
            P_hgcer_xAtCer_vs_MM_DATA.SetMinimum(1e-6)
            P_hgcer_xAtCer_vs_MM_DATA.Draw("colz")
            c_hgcervsMM.cd(2)
            P_hgcer_nohole_xAtCer_vs_MM_DATA.SetMinimum(1e-6)
            P_hgcer_nohole_xAtCer_vs_MM_DATA.Draw("colz")
            c_hgcervsMM.cd(3)
            P_hgcer_yAtCer_vs_MM_DATA.SetMinimum(1e-6)
            P_hgcer_yAtCer_vs_MM_DATA.Draw("colz")
            c_hgcervsMM.cd(4)
            P_hgcer_nohole_yAtCer_vs_MM_DATA.SetMinimum(1e-6)
            P_hgcer_nohole_yAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.Draw()

        c_hgcervsMM.Print(pdf_destinations["hgcer_debug"])

        ##
        # HGCer Hole Plots
        c_hgcer_hole = TCanvas()

        c_hgcer_hole.Divide(2,2)

        c_hgcer_hole.cd(1)
        P_hgcer_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_yAtCer_DATA.Draw("colz")

        c_hgcer_hole.cd(2)
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Draw("colz")

        c_hgcer_hole.cd(3)
        P_hgcer_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_yAtCer_DATA.Draw("colz")
        hgcer_cutg.SetLineColor(7)
        hgcer_cutg.Draw("same")

        c_hgcer_hole.cd(4)
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Draw("colz")
        hgcer_cutg.SetLineColor(7)
        hgcer_cutg.Draw("same")    

        c_hgcer_hole.Draw()

        diagnostic_pdf = pdf_destinations["hgcer_debug"]
        if isinstance(pion_hgcer_tdelta_diagnostic, dict):
            # C.4 routes detector diagnostics to their automatic supplement.
            # The supplement is closed after every terminal renderer returns.
            c_hgcer_hole.Print(diagnostic_pdf)
            emitted_hgcer_pages = render_pion_hgcer_tdelta_pages(
                diagnostic_pdf,
                pion_hgcer_tdelta_diagnostic,
                title_prefix="{} {}".format(phi_setting, ParticleType),
                page_manifest=histDict.setdefault(
                    "pion_hgcer_tdelta_diagnostic_page_manifest", []
                ),
                close_pdf=False,
            )
            histDict["pion_hgcer_tdelta_diagnostic"] = (
                serialize_pion_hgcer_tdelta_diagnostic(
                    pion_hgcer_tdelta_diagnostic,
                    include_records=False,
                )
            )
            if pion_hgcer_tdelta_json is not None:
                write_pion_hgcer_tdelta_json(
                    pion_hgcer_tdelta_json,
                    pion_hgcer_tdelta_diagnostic,
                )
            _print_rand_debug(
                "pion HGCer t-delta diagnostic pages",
                emitted=len(emitted_hgcer_pages),
                status=pion_hgcer_tdelta_diagnostic.get("status"),
            )
            if isinstance(pion_hgcer_zerope_transfer, dict):
                # These ROOT objects already exist at this terminal display
                # stage.  They are passed ephemerally to Part 2 only for
                # closure overlays; the renderer never reloads, refits,
                # normalizes, or serializes them.
                def _part2_display_target(fit_result, application_payload, proposal_payload=None):
                    application_payload = application_payload if isinstance(application_payload, dict) else {}
                    proposal_payload = proposal_payload if isinstance(proposal_payload, dict) else {}
                    return (
                        application_payload.get("H_MM_nosub_after_pion_subtraction")
                        or application_payload.get("H_MM_nosub_before_pion_subtraction")
                        or proposal_payload.get("H_MM_nosub_after_pion_subtraction")
                        or proposal_payload.get("H_MM_nosub_before_pion_subtraction")
                        or (fit_result or {}).get("H_kaon_nosub_input")
                    )

                def _part2_lambda_display(fit_result, application_payload, scope_label, proposal_payload=None):
                    target = _part2_display_target(fit_result, application_payload, proposal_payload)
                    diagnostics = (fit_result or {}).get("diagnostics") or {}
                    protected = diagnostics.get("pi_delta_signal_protected_fit") or (diagnostics.get("kaon") or {}).get("pi_delta_signal_protected_fit") or {}
                    reference_hists = [
                        (fit_result or {}).get("H_k_lambda_simc_reference"),
                        (fit_result or {}).get("H_simc_shape_k_lambda"),
                        (fit_result or {}).get("H_pi_delta_lambda_gauge") if protected else None,
                    ]
                    snapshots = {
                        id(reference): [
                            (float(reference.GetBinContent(index)), float(reference.GetBinError(index)))
                            for index in range(0, int(reference.GetNbinsX()) + 2)
                        ]
                        for reference in reference_hists if reference is not None
                    }
                    try:
                        display_hist, scale, source, normalization = _resolve_kaon_lambda_reference_for_plot(
                            fit_result, target,
                            (float(inpDict["mm_min"]), float(inpDict["mm_max"])),
                            scope_label, "H_part2_k_lambda_display",
                        )
                        for reference in reference_hists:
                            if reference is not None and snapshots.get(id(reference)) != [
                                (float(reference.GetBinContent(index)), float(reference.GetBinError(index)))
                                for index in range(0, int(reference.GetNbinsX()) + 2)
                            ]:
                                raise RuntimeError("canonical K-Lambda template changed during display normalization")
                        return {"hist": display_hist, "scale": scale, "source": source, "normalization": normalization, "status": "available"}
                    except Exception as exc:
                        return {"hist": None, "status": "unavailable", "reason": str(exc)}

                # The control audit is intentionally performed only for the
                # detached presentation comparison.  It cannot alter the
                # already-frozen map or its proposed application products.
                try:
                    part2_control_audit = audit_pion_hgcer_control_population(
                        pion_hgcer_zerope_transfer, pion_control_cache
                    )
                except Exception as exc:
                    part2_control_audit = {
                        "status": "unavailable",
                        "reason": "control_population_audit_exception: {}".format(exc),
                        "by_t": {}, "global": None,
                    }
                # Retain only scalar audit evidence in the Part-2 sidecar.
                # The detached ROOT projections remain terminal renderer
                # inputs, never serialized analysis products.
                pion_hgcer_zerope_transfer["control_population_audit"] = {
                    "status": part2_control_audit.get("status"),
                    "reason": part2_control_audit.get("reason"),
                    "by_t": {
                        str(t_index): {
                            key: value for key, value in (entry or {}).items()
                            if key != "before"
                        }
                        for t_index, entry in (part2_control_audit.get("by_t") or {}).items()
                    },
                    "global": {
                        key: value for key, value in (part2_control_audit.get("global") or {}).items()
                        if key != "before"
                    } if isinstance(part2_control_audit.get("global"), dict) else None,
                }
                part2_renderer_inputs = {"by_t": {}, "global": {}, "control": part2_control_audit}
                for parent in t_bin_parent_results:
                    parent_fit = (parent or {}).get("fit_result") or {}
                    try:
                        parent_t_index = int((parent or {}).get("t_bin_index"))
                    except (TypeError, ValueError):
                        continue
                    part2_renderer_inputs["by_t"][parent_t_index] = {
                        "simc": {
                            "pi-n": parent_fit.get("H_kaon_fit_pi_n_scaled"),
                            "pi-delta": parent_fit.get("H_kaon_fit_pi_delta_scaled"),
                            "pi-SIDIS": parent_fit.get("H_kaon_fit_pi_sidis_scaled"),
                            "current total pion model": parent_fit.get("H_kaon_pion_bg_fit_total"),
                        },
                        "signal": {
                            "K-Lambda": _part2_lambda_display(
                                parent_fit,
                                (parent or {}).get("final_diagnostic_application_result"),
                                (parent or {}).get("analysis_scope") or "Part2_t{}".format(parent_t_index + 1),
                                (parent or {}).get("proposed_diagnostic_application_result"),
                            ),
                            "K-Sigma0": {"hist": parent_fit.get("H_kaon_fit_k_sigma0_scaled"), "status": "available" if parent_fit.get("H_kaon_fit_k_sigma0_scaled") is not None else "unavailable"},
                        },
                    }
                if isinstance(component_fit_result, dict):
                    part2_renderer_inputs["global"] = {
                        "simc": {
                            "pi-n": component_fit_result.get("H_kaon_fit_pi_n_scaled"),
                            "pi-delta": component_fit_result.get("H_kaon_fit_pi_delta_scaled"),
                            "pi-SIDIS": component_fit_result.get("H_kaon_fit_pi_sidis_scaled"),
                            "current total pion model": component_fit_result.get("H_kaon_pion_bg_fit_total"),
                        },
                        "signal": {
                            "K-Lambda": _part2_lambda_display(
                                component_fit_result, setting_wide_render_payload,
                                component_fit_result.get("analysis_scope") or "Part2_setting_wide",
                            ),
                            "K-Sigma0": {"hist": component_fit_result.get("H_kaon_fit_k_sigma0_scaled"), "status": "available" if component_fit_result.get("H_kaon_fit_k_sigma0_scaled") is not None else "unavailable"},
                        },
                    }
                part2_manifest = histDict.setdefault(
                    "pion_hgcer_zerope_transfer_page_manifest", []
                )
                try:
                    emitted_part2_pages = render_pion_hgcer_zerope_transfer_pages(
                        diagnostic_pdf,
                        pion_hgcer_zerope_transfer,
                        title_prefix="{} {}".format(phi_setting, ParticleType),
                        page_manifest=part2_manifest,
                        close_pdf=False,
                        renderer_inputs=part2_renderer_inputs,
                    )
                    pion_hgcer_zerope_transfer["rendering_status"] = "available"
                except Exception as exc:
                    # Part 2 is terminal diagnostic presentation only.  A
                    # failed page contract must still leave a closed PDF and
                    # cannot propagate into Part 1/1.5 or production state.
                    reason = "{}: {}".format(type(exc).__name__, exc)
                    emitted_part2_pages = render_pion_hgcer_zerope_transfer_failure_page(
                        diagnostic_pdf, reason,
                        title_prefix="{} {}".format(phi_setting, ParticleType),
                        page_manifest=part2_manifest,
                        close_pdf=False,
                    )
                    pion_hgcer_zerope_transfer["rendering_status"] = "unavailable"
                    pion_hgcer_zerope_transfer["rendering_failure_reason"] = reason
                pion_hgcer_zerope_transfer["render_manifest"] = emitted_part2_pages
                histDict["pion_hgcer_zerope_transfer"] = (
                    serialize_pion_hgcer_zerope_transfer(
                        pion_hgcer_zerope_transfer
                    )
                )
                if pion_hgcer_zerope_transfer_json is not None:
                    write_pion_hgcer_zerope_transfer_json(
                        pion_hgcer_zerope_transfer_json,
                        pion_hgcer_zerope_transfer,
                    )
                _print_rand_debug(
                    "pion HGCer zero-PE transfer pages",
                    emitted=len(emitted_part2_pages),
                    status=pion_hgcer_zerope_transfer.get("status"),
                )
        else:
            c_hgcer_hole.Print(diagnostic_pdf)

    if ParticleType == "kaon":
        hgcer_refinement_manifest = histDict.setdefault(
            "pion_hgcer_refinement_page_manifest", []
        )
        method_b_display = method_b_display_payload(
            pion_hgcer_method_b, checkpoint_for_plots
        )
        if (
            method_b_display.get("source") == "checkpoint_method_b"
            and not method_b_display.get("source_complete")
        ):
            setting_renderer_failures.append(
                "Method-B checkpoint display payload incomplete"
            )
        method_b_source_parity = method_b_display_source_parity(
            pion_hgcer_method_b, checkpoint_for_plots
        )
        if (
            method_b_source_parity.get("checked")
            and not method_b_source_parity.get("passed")
        ):
            setting_renderer_failures.append(
                "Method-B display-source parity mismatch: {}".format(
                    ", ".join(method_b_source_parity.get("differences") or ())
                )
            )
        method_b_coverage_parity = method_b_display.get("coverage_parity") or {}
        if (
            method_b_display.get("source") == "checkpoint_method_b"
            and method_b_coverage_parity.get("checked")
            and not method_b_coverage_parity.get("passed")
        ):
            setting_renderer_failures.append(
                "Method-B checkpoint cell/summary coverage mismatch: {}".format(
                    ", ".join(method_b_coverage_parity.get("differences") or ())
                )
            )
        aerogel_validation_for_qa = proton_display_diagnostics.get(
            "aerogel_vs_t_validation"
        ) or {}
        setting_qa_context = {
            "aerogel_warnings": aerogel_validation_for_qa.get(
                "warnings_by_t_bin", aerogel_validation_for_qa.get("warnings", "not available")
            ),
            "proton_warnings": proton_display_gate.get(
                "observational_warnings", "not available"
            ),
            "canonical_parent_k_lambda": histDict.get(
                "canonical_parent_k_lambda_qa", []
            ),
            "k_lambda_comparison": "canonical-parent summary",
            "k_sigma0_protected_region": (
                "active" if t_bin_parent_results else "not recorded"
            ),
            "k_sigma0_availability": (
                "available"
                if isinstance(component_fit_result, dict)
                and component_fit_result.get("H_kaon_fit_k_sigma0_scaled") is not None
                else "not available"
            ),
            "hgcer_diagnostic_availability": (
                pion_hgcer_tdelta_diagnostic.get("status", "not available")
                if isinstance(pion_hgcer_tdelta_diagnostic, dict) else "not available"
            ),
            "renderer_failures": setting_renderer_failures,
            "method_b_display_payload": method_b_display,
        }
        try:
            render_pion_hgcer_refinement_pages(
                main_pdf,
                checkpoint_for_plots,
                phase_a=pion_hgcer_event_contract,
                method_a=pion_hgcer_method_a,
                method_b=pion_hgcer_method_b,
                method_b_display=method_b_display,
                phase_a_display_context=phase_a_display_context,
                page_manifest=hgcer_refinement_manifest,
            )
        except Exception as exc:
            setting_renderer_failures.append(
                "refinement pages: {}: {}".format(type(exc).__name__, exc)
            )
            _print_rand_debug(
                "detached HGCer refinement PDF pages unavailable",
                renderer="pion_hgcer_refinement_plots",
                exception_type=type(exc).__name__,
                exception=str(exc),
            )
        try:
            render_setting_warning_page(
                main_pdf,
                checkpoint_for_plots,
                phase_a=pion_hgcer_event_contract,
                method_a=pion_hgcer_method_a,
                method_b=pion_hgcer_method_b,
                method_b_display=method_b_display,
                part2=pion_hgcer_zerope_transfer,
                phase_a_display_context=phase_a_display_context,
                runtime_qa_context=setting_qa_context,
                page_manifest=hgcer_refinement_manifest,
                close_pdf=True,
            )
        except Exception as exc:
            _print_rand_debug(
                "detached HGCer warning page unavailable",
                renderer="pion_hgcer_refinement_plots",
                page_id="qa.setting_warnings",
                exception_type=type(exc).__name__,
                exception=str(exc),
            )
            # Preserve a closed primary PDF even if a detached renderer fails.
            ct.Print(main_pdf + ')')
        if isinstance(phase_d_checkpoint_for_plots, dict):
            try:
                render_pion_hgcer_ab_comparison_pages(
                    pdf_destinations["phase_d_ab"],
                    phase_d_checkpoint_for_plots,
                    page_manifest=supplement_manifests.get("phase_d_ab", []),
                )
            except Exception as exc:
                setting_renderer_failures.append(
                    "Phase-D A/B pages: {}: {}".format(type(exc).__name__, exc)
                )
                _print_rand_debug(
                    "detached Phase-D A/B PDF pages unavailable",
                    renderer="pion_hgcer_refinement_plots",
                    exception_type=type(exc).__name__,
                    exception=str(exc),
                )
        # Phases D.6 through D.10 are terminal presentation only.  They receive the
        # already constructed proton-cleaning products and retain no
        # ROOT-bearing object in histDict or any downstream physics input.
        full_background_subtraction_manifest = []
        full_background_subtraction_failures = []
        full_background_subtraction_pdf = full_background_subtraction_pdf_path(
            main_pdf
        )
        full_background_subtraction_open = False
        try:
            full_background_subtraction_d6_payload = (
                build_full_background_subtraction_d6_payload(
                    proton_cleaning_result,
                    proton_cleaning_application,
                )
            )
            full_background_subtraction_d7_payload = (
                build_full_background_subtraction_d7_payload(
                    proton_cleaning_result,
                    proton_cleaning_application,
                    proton_cleaning_tree_bundle,
                )
            )
            full_background_subtraction_d8_payload = (
                build_full_background_subtraction_d8_payload(
                    histDict.get("_pion_t_bin_parent_results"),
                    histDict.get("_authoritative_pion_control_source_cache"),
                )
            )
            full_background_subtraction_d9_payload = (
                build_full_background_subtraction_d9_payload(
                    pion_hgcer_tdelta_diagnostic,
                    pion_hgcer_method_a,
                    pion_hgcer_method_a_comparison,
                )
            )
            full_background_subtraction_d10_payload = (
                build_full_background_subtraction_d10_payload(
                    pion_hgcer_event_contract,
                    pion_hgcer_method_b,
                    pion_hgcer_method_b_comparison,
                )
            )
            if (
                not full_background_subtraction_d6_payload.get("available")
                and not full_background_subtraction_d7_payload.get("available")
                and not full_background_subtraction_d8_payload.get("available")
                and not full_background_subtraction_d9_payload.get("available")
                and not full_background_subtraction_d10_payload.get("available")
            ):
                full_background_subtraction_failures.extend((
                    "D.6 procedure input unavailable: {}".format(
                        full_background_subtraction_d6_payload.get("reason")
                    ),
                    "D.7 procedure input unavailable: {}".format(
                        full_background_subtraction_d7_payload.get("reason")
                    ),
                    "D.8 procedure input unavailable: {}".format(
                        full_background_subtraction_d8_payload.get("reason")
                    ),
                    "D.9 procedure input unavailable: {}".format(
                        full_background_subtraction_d9_payload.get("reason")
                    ),
                    "D.10 procedure input unavailable: {}".format(
                        full_background_subtraction_d10_payload.get("reason")
                    ),
                ))
            else:
                full_background_subtraction_open = (
                    open_full_background_subtraction_pdf(
                        full_background_subtraction_pdf
                    )
                )
                if not full_background_subtraction_open:
                    full_background_subtraction_failures.append(
                        "full background-subtraction rendering unavailable: PyROOT not available"
                    )
                else:
                    rendered_procedure = (
                        render_full_background_subtraction_procedure_pages(
                            full_background_subtraction_pdf,
                            full_background_subtraction_d6_payload,
                            full_background_subtraction_d7_payload,
                            full_background_subtraction_d8_payload,
                            full_background_subtraction_d9_payload,
                            full_background_subtraction_d10_payload,
                            page_manifest=full_background_subtraction_manifest,
                        )
                    )
                    full_background_subtraction_failures.extend(
                        rendered_procedure.get("failures") or ()
                    )
        except Exception as exc:
            full_background_subtraction_failures.append(
                "full background-subtraction procedure pages: {}: {}".format(
                    type(exc).__name__, exc
                )
            )
            _print_rand_debug(
                "detached full background-subtraction D.6 through D.10 pages unavailable",
                renderer="full_background_subtraction_plots",
                exception_type=type(exc).__name__,
                exception=str(exc),
            )
        finally:
            if full_background_subtraction_open:
                try:
                    close_full_background_subtraction_pdf(
                        full_background_subtraction_pdf
                    )
                except Exception as exc:
                    full_background_subtraction_failures.append(
                        "full background-subtraction procedure PDF close: {}: {}".format(
                            type(exc).__name__, exc
                        )
                    )
        histDict["full_background_subtraction_page_manifest"] = [
            dict(page) for page in full_background_subtraction_manifest
        ]
        histDict["full_background_subtraction_renderer_failures"] = list(
            full_background_subtraction_failures
        )
        setting_renderer_failures.extend(full_background_subtraction_failures)
        for supplement_key, role in (
            ("proton_debug", "proton-debug"),
            ("pion_fit_debug", "pion-fit-debug"),
            ("hgcer_debug", "hgcer-debug"),
        ):
            close_diagnostic_pdf(
                pdf_destinations[supplement_key],
                checkpoint_for_plots,
                role=role,
                manifest=supplement_manifests.get(supplement_key, []),
            )
        if "phase_d_ab" in supplement_manifests:
            close_diagnostic_pdf(
                pdf_destinations["phase_d_ab"],
                phase_d_checkpoint_for_plots,
                role="hgcer-ab-comparison",
                manifest=supplement_manifests["phase_d_ab"],
            )

    _print_rand_timer("rand_sub plotting {}".format(phi_setting), perf_counter() - stage_start)
    _print_rand_timer("rand_sub total {}".format(phi_setting), perf_counter() - total_start)
    return histDict
