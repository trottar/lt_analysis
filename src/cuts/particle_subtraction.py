#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-09 14:04:49 trottar"
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
from ROOT import TH1D, TH2D, TCutG, TFile
import sys, os, math
from time import perf_counter

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
from utility import open_root_file
from prompt_trees import get_prompt_tree_name, get_rand_tree_name

################################################################################################################################################

def _format_elapsed(seconds):
    if seconds < 60.0:
        return "{:.2f} s".format(seconds)
    minutes, remainder = divmod(seconds, 60.0)
    if minutes < 60.0:
        return "{:.0f} m {:.2f} s".format(minutes, remainder)
    hours, minutes = divmod(minutes, 60.0)
    return "{:.0f} h {:.0f} m {:.2f} s".format(hours, minutes, remainder)


def _print_sub_timer(label, elapsed, total_events=None):
    if total_events and total_events > 0:
        per_event_ms = (elapsed / total_events) * 1000.0
        print("[TIMER] {}: {} ({:.3f} ms/event)".format(label, _format_elapsed(elapsed), per_event_ms))
    else:
        print("[TIMER] {}: {}".format(label, _format_elapsed(elapsed)))


def _init_sub_event_cache():
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
        "t_index": [],
        "phi_index": [],
    }
    return {
        key: {name: values.copy() for name, values in cache_template.items()}
        for key in ("prompt", "dummy", "rand", "dummy_rand")
    }


def _build_sub_allcut_bin_index(cache_section):
    if len(cache_section["allcuts"]) == 0:
        return {}

    index_map = {}
    allcut_indices = np.flatnonzero(cache_section["allcuts"])
    for idx in allcut_indices:
        key = (int(cache_section["t_index"][idx]), int(cache_section["phi_index"][idx]))
        if key not in index_map:
            index_map[key] = []
        index_map[key].append(int(idx))
    return {
        key: np.asarray(indices, dtype=np.int32)
        for key, indices in index_map.items()
    }


def _append_sub_event(cache_section, adj_t, adj_MM, theta_cm_deg, q2, w, epsilon, ssxptar, ssyptar, hsxptar, hsyptar, allcuts, nommcuts, t_index=-1, phi_index=-1):
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
    cache_section["t_index"].append(int(t_index))
    cache_section["phi_index"].append(int(phi_index))


def _freeze_sub_event_cache(event_cache):
    frozen_cache = {}
    for cache_key, cache_section in event_cache.items():
        frozen_cache[cache_key] = {
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
            "t_index": np.asarray(cache_section["t_index"], dtype=np.int32),
            "phi_index": np.asarray(cache_section["phi_index"], dtype=np.int32),
        }
        frozen_cache[cache_key]["allcut_bin_index"] = _build_sub_allcut_bin_index(frozen_cache[cache_key])
    return frozen_cache


def _fill_particle_subtraction_allcuts(evt, adj_MM, adj_t, adj_hsdelta, fills):
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
    fills["h_cal"](evt.H_cal_etottracknorm)
    fills["h_cer"](evt.H_cer_npeSum)
    fills["p_cal"](evt.P_cal_etottracknorm)
    fills["p_hgcer"](evt.P_hgcer_npeSum)
    fills["p_aero"](evt.P_aero_npeSum)


def _process_particle_subtraction_tree(
    tree,
    print_label,
    timer_label,
    mm_offset_data,
    fills,
    particle_type,
    hole_contains,
    evaluate_event,
    shifted_mm_getter,
    shifted_t_getter,
    effective_mm_offset_getter,
    mm_min,
    mm_max,
    progress_bar,
    kaon_nomm_enabled=True,
):
    print(print_label)
    entries = tree.GetEntries()
    progress_time = 0.0
    loop_start = perf_counter()
    nohole_xy_fill = fills["nohole_xy"]
    nohole_x_mm_fill = fills["nohole_x_mm"]
    nohole_y_mm_fill = fills["nohole_y_mm"]
    nomm_fill = fills["nomm"]
    mm_offset_correction = effective_mm_offset_getter(mm_offset_data)

    for i, evt in enumerate(tree):
        progress_start = perf_counter()
        progress_bar(i, entries, bar_length=25)
        progress_time += perf_counter() - progress_start

        if particle_type == "kaon":
            base_all_cuts, base_sub_cuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max, mm_offset=mm_offset_correction)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            pid_pass = evt.P_hgcer_npeSum > 2.0
            allcuts = base_all_cuts and not hole_rejected and pid_pass
            noholecuts = base_sub_cuts
            nommcuts = base_sub_cuts and not hole_rejected and pid_pass
        else:
            allcuts, nommcuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max, mm_offset=mm_offset_correction)
            noholecuts = False

        if not (noholecuts or nommcuts or allcuts):
            continue

        adj_MM = shifted_mm_getter(evt, mm_offset=mm_offset_data)

        if noholecuts and nohole_xy_fill is not None:
            nohole_xy_fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
            nohole_x_mm_fill(evt.P_hgcer_xAtCer, adj_MM, evt.P_hgcer_npeSum)
            nohole_y_mm_fill(evt.P_hgcer_yAtCer, adj_MM, evt.P_hgcer_npeSum)

        if nommcuts and (particle_type != "kaon" or kaon_nomm_enabled):
            nomm_fill(adj_MM)

        if allcuts:
            adj_t = shifted_t_getter(evt)
            _fill_particle_subtraction_allcuts(evt, adj_MM, adj_t, adj_hsdelta, fills)

    loop_elapsed = perf_counter() - loop_start
    _print_sub_timer(timer_label, loop_elapsed, entries)
    _print_sub_timer("{} progressBar".format(timer_label), progress_time, entries)
    _print_sub_timer("{} other".format(timer_label), max(loop_elapsed - progress_time, 0.0), entries)


def particle_subtraction_cuts(histDict, subDict, inpDict, SubtractedParticle, hgcer_cutg=None):
    total_start = perf_counter()
    setup_start = perf_counter()

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    sys.path.append("normalize")
    from get_eff_charge import get_eff_charge

    # Upate hist dictionary with effective charge
    get_eff_charge(histDict, inpDict, all_data=False)

    norm_factor_data = inpDict["normfac_data"]
    norm_factor_dummy = inpDict["normfac_dummy"]

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, evaluate_data_cut_bools, evaluate_data_event, get_effective_mm_offset, get_shifted_mm, get_shifted_t, set_shift_context, set_val
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=inpDict.get("shift_mode", "raw"))
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    prompt_tree_name = get_prompt_tree_name(SubtractedParticle, EPSSET)
    rand_tree_name = get_rand_tree_name(SubtractedParticle, EPSSET)

    TBRANCH_DATA  = InFile_DATA.Get(prompt_tree_name)

    TBRANCH_RAND  = InFile_DATA.Get(rand_tree_name)

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get(prompt_tree_name)
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get(rand_tree_name)

    ################################################################################################################################################
    
    H_hsdelta_DATA = subDict["H_hsdelta_SUB_DATA"]
    H_hsxptar_DATA = subDict["H_hsxptar_SUB_DATA"]
    H_hsyptar_DATA = subDict["H_hsyptar_SUB_DATA"]
    H_ssxfp_DATA = subDict["H_ssxfp_SUB_DATA"]
    H_ssyfp_DATA = subDict["H_ssyfp_SUB_DATA"]
    H_ssxpfp_DATA = subDict["H_ssxpfp_SUB_DATA"]
    H_ssypfp_DATA = subDict["H_ssypfp_SUB_DATA"]
    H_hsxfp_DATA = subDict["H_hsxfp_SUB_DATA"]
    H_hsyfp_DATA = subDict["H_hsyfp_SUB_DATA"]
    H_hsxpfp_DATA = subDict["H_hsxpfp_SUB_DATA"]
    H_hsypfp_DATA = subDict["H_hsypfp_SUB_DATA"]
    H_ssdelta_DATA = subDict["H_ssdelta_SUB_DATA"]
    H_ssxptar_DATA = subDict["H_ssxptar_SUB_DATA"]
    H_ssyptar_DATA = subDict["H_ssyptar_SUB_DATA"]
    H_q_DATA = subDict["H_q_SUB_DATA"]
    H_Q2_DATA = subDict["H_Q2_SUB_DATA"]
    H_W_DATA = subDict["H_W_SUB_DATA"]
    H_t_DATA = subDict["H_t_SUB_DATA"]
    H_epsilon_DATA = subDict["H_epsilon_SUB_DATA"]
    H_MM_DATA = subDict["H_MM_SUB_DATA"]
    H_MM_nosub_DATA = subDict["H_MM_nosub_SUB_DATA"]
    H_th_DATA = subDict["H_th_SUB_DATA"]
    H_ph_DATA = subDict["H_ph_SUB_DATA"]
    H_ph_q_DATA = subDict["H_ph_q_SUB_DATA"]
    H_th_q_DATA = subDict["H_th_q_SUB_DATA"]
    H_ph_recoil_DATA = subDict["H_ph_recoil_SUB_DATA"]
    H_th_recoil_DATA = subDict["H_th_recoil_SUB_DATA"]
    H_pmiss_DATA = subDict["H_pmiss_SUB_DATA"]
    H_emiss_DATA = subDict["H_emiss_SUB_DATA"]
    H_pmx_DATA = subDict["H_pmx_SUB_DATA"]
    H_pmy_DATA = subDict["H_pmy_SUB_DATA"]
    H_pmz_DATA = subDict["H_pmz_SUB_DATA"]
    H_ct_DATA = subDict["H_ct_SUB_DATA"]
    H_cal_etottracknorm_DATA = subDict["H_cal_etottracknorm_SUB_DATA"]
    H_cer_npeSum_DATA = subDict["H_cer_npeSum_SUB_DATA"]
    P_cal_etottracknorm_DATA = subDict["P_cal_etottracknorm_SUB_DATA"]
    P_hgcer_npeSum_DATA = subDict["P_hgcer_npeSum_SUB_DATA"]
    P_aero_npeSum_DATA = subDict["P_aero_npeSum_SUB_DATA"]

    H_hsdelta_DUMMY = subDict["H_hsdelta_SUB_DUMMY"]
    H_hsxptar_DUMMY = subDict["H_hsxptar_SUB_DUMMY"]
    H_hsyptar_DUMMY = subDict["H_hsyptar_SUB_DUMMY"]
    H_ssxfp_DUMMY = subDict["H_ssxfp_SUB_DUMMY"]
    H_ssyfp_DUMMY = subDict["H_ssyfp_SUB_DUMMY"]
    H_ssxpfp_DUMMY = subDict["H_ssxpfp_SUB_DUMMY"]
    H_ssypfp_DUMMY = subDict["H_ssypfp_SUB_DUMMY"]
    H_hsxfp_DUMMY = subDict["H_hsxfp_SUB_DUMMY"]
    H_hsyfp_DUMMY = subDict["H_hsyfp_SUB_DUMMY"]
    H_hsxpfp_DUMMY = subDict["H_hsxpfp_SUB_DUMMY"]
    H_hsypfp_DUMMY = subDict["H_hsypfp_SUB_DUMMY"]
    H_ssdelta_DUMMY = subDict["H_ssdelta_SUB_DUMMY"]
    H_ssxptar_DUMMY = subDict["H_ssxptar_SUB_DUMMY"]
    H_ssyptar_DUMMY = subDict["H_ssyptar_SUB_DUMMY"]
    H_q_DUMMY = subDict["H_q_SUB_DUMMY"]
    H_Q2_DUMMY = subDict["H_Q2_SUB_DUMMY"]
    H_W_DUMMY = subDict["H_W_SUB_DUMMY"]
    H_t_DUMMY = subDict["H_t_SUB_DUMMY"]
    H_epsilon_DUMMY = subDict["H_epsilon_SUB_DUMMY"]
    H_MM_DUMMY = subDict["H_MM_SUB_DUMMY"]
    H_MM_nosub_DUMMY = subDict["H_MM_nosub_SUB_DUMMY"]    
    H_th_DUMMY = subDict["H_th_SUB_DUMMY"]
    H_ph_DUMMY = subDict["H_ph_SUB_DUMMY"]
    H_ph_q_DUMMY = subDict["H_ph_q_SUB_DUMMY"]
    H_th_q_DUMMY = subDict["H_th_q_SUB_DUMMY"]
    H_ph_recoil_DUMMY = subDict["H_ph_recoil_SUB_DUMMY"]
    H_th_recoil_DUMMY = subDict["H_th_recoil_SUB_DUMMY"]
    H_pmiss_DUMMY = subDict["H_pmiss_SUB_DUMMY"]
    H_emiss_DUMMY = subDict["H_emiss_SUB_DUMMY"]
    H_pmx_DUMMY = subDict["H_pmx_SUB_DUMMY"]
    H_pmy_DUMMY = subDict["H_pmy_SUB_DUMMY"]
    H_pmz_DUMMY = subDict["H_pmz_SUB_DUMMY"]
    H_ct_DUMMY = subDict["H_ct_SUB_DUMMY"]
    H_cal_etottracknorm_DUMMY = subDict["H_cal_etottracknorm_SUB_DUMMY"]
    H_cer_npeSum_DUMMY = subDict["H_cer_npeSum_SUB_DUMMY"]
    P_cal_etottracknorm_DUMMY = subDict["P_cal_etottracknorm_SUB_DUMMY"]
    P_hgcer_npeSum_DUMMY = subDict["P_hgcer_npeSum_SUB_DUMMY"]
    P_aero_npeSum_DUMMY = subDict["P_aero_npeSum_SUB_DUMMY"]

    H_hsdelta_RAND = subDict["H_hsdelta_SUB_RAND"]
    H_hsxptar_RAND = subDict["H_hsxptar_SUB_RAND"]
    H_hsyptar_RAND = subDict["H_hsyptar_SUB_RAND"]
    H_ssxfp_RAND = subDict["H_ssxfp_SUB_RAND"]
    H_ssyfp_RAND = subDict["H_ssyfp_SUB_RAND"]
    H_ssxpfp_RAND = subDict["H_ssxpfp_SUB_RAND"]
    H_ssypfp_RAND = subDict["H_ssypfp_SUB_RAND"]
    H_hsxfp_RAND = subDict["H_hsxfp_SUB_RAND"]
    H_hsyfp_RAND = subDict["H_hsyfp_SUB_RAND"]
    H_hsxpfp_RAND = subDict["H_hsxpfp_SUB_RAND"]
    H_hsypfp_RAND = subDict["H_hsypfp_SUB_RAND"]
    H_ssdelta_RAND = subDict["H_ssdelta_SUB_RAND"]
    H_ssxptar_RAND = subDict["H_ssxptar_SUB_RAND"]
    H_ssyptar_RAND = subDict["H_ssyptar_SUB_RAND"]
    H_q_RAND = subDict["H_q_SUB_RAND"]
    H_Q2_RAND = subDict["H_Q2_SUB_RAND"]
    H_W_RAND = subDict["H_W_SUB_RAND"]
    H_t_RAND = subDict["H_t_SUB_RAND"]
    H_epsilon_RAND = subDict["H_epsilon_SUB_RAND"]
    H_MM_RAND = subDict["H_MM_SUB_RAND"]
    H_MM_nosub_RAND = subDict["H_MM_nosub_SUB_RAND"]    
    H_th_RAND = subDict["H_th_SUB_RAND"]
    H_ph_RAND = subDict["H_ph_SUB_RAND"]
    H_ph_q_RAND = subDict["H_ph_q_SUB_RAND"]
    H_th_q_RAND = subDict["H_th_q_SUB_RAND"]
    H_ph_recoil_RAND = subDict["H_ph_recoil_SUB_RAND"]
    H_th_recoil_RAND = subDict["H_th_recoil_SUB_RAND"]
    H_pmiss_RAND = subDict["H_pmiss_SUB_RAND"]
    H_emiss_RAND = subDict["H_emiss_SUB_RAND"]
    H_pmx_RAND = subDict["H_pmx_SUB_RAND"]
    H_pmy_RAND = subDict["H_pmy_SUB_RAND"]
    H_pmz_RAND = subDict["H_pmz_SUB_RAND"]
    H_ct_RAND = subDict["H_ct_SUB_RAND"]
    H_cal_etottracknorm_RAND = subDict["H_cal_etottracknorm_SUB_RAND"]
    H_cer_npeSum_RAND = subDict["H_cer_npeSum_SUB_RAND"]
    P_cal_etottracknorm_RAND = subDict["P_cal_etottracknorm_SUB_RAND"]
    P_hgcer_npeSum_RAND = subDict["P_hgcer_npeSum_SUB_RAND"]
    P_aero_npeSum_RAND = subDict["P_aero_npeSum_SUB_RAND"]

    H_hsdelta_DUMMY_RAND = subDict["H_hsdelta_SUB_DUMMY_RAND"]
    H_hsxptar_DUMMY_RAND = subDict["H_hsxptar_SUB_DUMMY_RAND"]
    H_hsyptar_DUMMY_RAND = subDict["H_hsyptar_SUB_DUMMY_RAND"]
    H_ssxfp_DUMMY_RAND = subDict["H_ssxfp_SUB_DUMMY_RAND"]
    H_ssyfp_DUMMY_RAND = subDict["H_ssyfp_SUB_DUMMY_RAND"]
    H_ssxpfp_DUMMY_RAND = subDict["H_ssxpfp_SUB_DUMMY_RAND"]
    H_ssypfp_DUMMY_RAND = subDict["H_ssypfp_SUB_DUMMY_RAND"]
    H_hsxfp_DUMMY_RAND = subDict["H_hsxfp_SUB_DUMMY_RAND"]
    H_hsyfp_DUMMY_RAND = subDict["H_hsyfp_SUB_DUMMY_RAND"]
    H_hsxpfp_DUMMY_RAND = subDict["H_hsxpfp_SUB_DUMMY_RAND"]
    H_hsypfp_DUMMY_RAND = subDict["H_hsypfp_SUB_DUMMY_RAND"]
    H_ssdelta_DUMMY_RAND = subDict["H_ssdelta_SUB_DUMMY_RAND"]
    H_ssxptar_DUMMY_RAND = subDict["H_ssxptar_SUB_DUMMY_RAND"]
    H_ssyptar_DUMMY_RAND = subDict["H_ssyptar_SUB_DUMMY_RAND"]
    H_q_DUMMY_RAND = subDict["H_q_SUB_DUMMY_RAND"]
    H_Q2_DUMMY_RAND = subDict["H_Q2_SUB_DUMMY_RAND"]
    H_W_DUMMY_RAND = subDict["H_W_SUB_DUMMY_RAND"]
    H_t_DUMMY_RAND = subDict["H_t_SUB_DUMMY_RAND"]
    H_epsilon_DUMMY_RAND = subDict["H_epsilon_SUB_DUMMY_RAND"]
    H_MM_DUMMY_RAND = subDict["H_MM_SUB_DUMMY_RAND"]
    H_MM_nosub_DUMMY_RAND = subDict["H_MM_nosub_SUB_DUMMY_RAND"]    
    H_th_DUMMY_RAND = subDict["H_th_SUB_DUMMY_RAND"]
    H_ph_DUMMY_RAND = subDict["H_ph_SUB_DUMMY_RAND"]
    H_ph_q_DUMMY_RAND = subDict["H_ph_q_SUB_DUMMY_RAND"]
    H_th_q_DUMMY_RAND = subDict["H_th_q_SUB_DUMMY_RAND"]
    H_ph_recoil_DUMMY_RAND = subDict["H_ph_recoil_SUB_DUMMY_RAND"]
    H_th_recoil_DUMMY_RAND = subDict["H_th_recoil_SUB_DUMMY_RAND"]
    H_pmiss_DUMMY_RAND = subDict["H_pmiss_SUB_DUMMY_RAND"]
    H_emiss_DUMMY_RAND = subDict["H_emiss_SUB_DUMMY_RAND"]
    H_pmx_DUMMY_RAND = subDict["H_pmx_SUB_DUMMY_RAND"]
    H_pmy_DUMMY_RAND = subDict["H_pmy_SUB_DUMMY_RAND"]
    H_pmz_DUMMY_RAND = subDict["H_pmz_SUB_DUMMY_RAND"]
    H_ct_DUMMY_RAND = subDict["H_ct_SUB_DUMMY_RAND"]
    H_cal_etottracknorm_DUMMY_RAND = subDict["H_cal_etottracknorm_SUB_DUMMY_RAND"]
    H_cer_npeSum_DUMMY_RAND = subDict["H_cer_npeSum_SUB_DUMMY_RAND"]
    P_cal_etottracknorm_DUMMY_RAND = subDict["P_cal_etottracknorm_SUB_DUMMY_RAND"]
    P_hgcer_npeSum_DUMMY_RAND = subDict["P_hgcer_npeSum_SUB_DUMMY_RAND"]
    P_aero_npeSum_DUMMY_RAND = subDict["P_aero_npeSum_SUB_DUMMY_RAND"]

    MM_vs_CoinTime_DATA = subDict["MM_vs_CoinTime_SUB_DATA"]
    CoinTime_vs_beta_DATA = subDict["CoinTime_vs_beta_SUB_DATA"]
    MM_vs_beta_DATA = subDict["MM_vs_beta_SUB_DATA"]
    MM_vs_H_cer_DATA = subDict["MM_vs_H_cer_SUB_DATA"]
    MM_vs_H_cal_DATA = subDict["MM_vs_H_cal_SUB_DATA"]
    MM_vs_P_cal_DATA = subDict["MM_vs_P_cal_SUB_DATA"]
    MM_vs_P_hgcer_DATA = subDict["MM_vs_P_hgcer_SUB_DATA"]
    MM_vs_P_aero_DATA = subDict["MM_vs_P_aero_SUB_DATA"]
    phiq_vs_t_DATA = subDict["phiq_vs_t_SUB_DATA"]
    Q2_vs_W_DATA = subDict["Q2_vs_W_SUB_DATA"]
    Q2_vs_t_DATA = subDict["Q2_vs_t_SUB_DATA"]
    W_vs_t_DATA = subDict["W_vs_t_SUB_DATA"]
    EPS_vs_t_DATA = subDict["EPS_vs_t_SUB_DATA"]
    MM_vs_t_DATA = subDict["MM_vs_t_SUB_DATA"]
    P_hgcer_xAtCer_vs_yAtCer_DATA = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DATA = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"]
    P_hgcer_xAtCer_vs_MM_DATA = subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"]
    P_hgcer_nohole_xAtCer_vs_MM_DATA = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"]
    P_hgcer_yAtCer_vs_MM_DATA = subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"]
    P_hgcer_nohole_yAtCer_vs_MM_DATA = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"]

    MM_vs_CoinTime_DUMMY = subDict["MM_vs_CoinTime_SUB_DUMMY"]
    CoinTime_vs_beta_DUMMY = subDict["CoinTime_vs_beta_SUB_DUMMY"]
    MM_vs_beta_DUMMY = subDict["MM_vs_beta_SUB_DUMMY"]
    MM_vs_H_cer_DUMMY = subDict["MM_vs_H_cer_SUB_DUMMY"]
    MM_vs_H_cal_DUMMY = subDict["MM_vs_H_cal_SUB_DUMMY"]
    MM_vs_P_cal_DUMMY = subDict["MM_vs_P_cal_SUB_DUMMY"]
    MM_vs_P_hgcer_DUMMY = subDict["MM_vs_P_hgcer_SUB_DUMMY"]
    MM_vs_P_aero_DUMMY = subDict["MM_vs_P_aero_SUB_DUMMY"]
    phiq_vs_t_DUMMY = subDict["phiq_vs_t_SUB_DUMMY"]
    Q2_vs_W_DUMMY = subDict["Q2_vs_W_SUB_DUMMY"]
    Q2_vs_t_DUMMY = subDict["Q2_vs_t_SUB_DUMMY"]
    W_vs_t_DUMMY = subDict["W_vs_t_SUB_DUMMY"]
    EPS_vs_t_DUMMY = subDict["EPS_vs_t_SUB_DUMMY"]
    MM_vs_t_DUMMY = subDict["MM_vs_t_SUB_DUMMY"]
    P_hgcer_xAtCer_vs_yAtCer_DUMMY = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY"]
    P_hgcer_xAtCer_vs_MM_DUMMY = subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_nohole_xAtCer_vs_MM_DUMMY = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_yAtCer_vs_MM_DUMMY = subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_nohole_yAtCer_vs_MM_DUMMY = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY"]

    MM_vs_CoinTime_RAND = subDict["MM_vs_CoinTime_SUB_RAND"]
    CoinTime_vs_beta_RAND = subDict["CoinTime_vs_beta_SUB_RAND"]
    MM_vs_beta_RAND = subDict["MM_vs_beta_SUB_RAND"]
    MM_vs_H_cer_RAND = subDict["MM_vs_H_cer_SUB_RAND"]
    MM_vs_H_cal_RAND = subDict["MM_vs_H_cal_SUB_RAND"]
    MM_vs_P_cal_RAND = subDict["MM_vs_P_cal_SUB_RAND"]
    MM_vs_P_hgcer_RAND = subDict["MM_vs_P_hgcer_SUB_RAND"]
    MM_vs_P_aero_RAND = subDict["MM_vs_P_aero_SUB_RAND"]
    phiq_vs_t_RAND = subDict["phiq_vs_t_SUB_RAND"]
    Q2_vs_W_RAND = subDict["Q2_vs_W_SUB_RAND"]
    Q2_vs_t_RAND = subDict["Q2_vs_t_SUB_RAND"]
    W_vs_t_RAND = subDict["W_vs_t_SUB_RAND"]
    EPS_vs_t_RAND = subDict["EPS_vs_t_SUB_RAND"]
    MM_vs_t_RAND = subDict["MM_vs_t_SUB_RAND"]
    P_hgcer_xAtCer_vs_yAtCer_RAND = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_RAND"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_RAND = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND"]
    P_hgcer_xAtCer_vs_MM_RAND = subDict["P_hgcer_xAtCer_vs_MM_SUB_RAND"]
    P_hgcer_nohole_xAtCer_vs_MM_RAND = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND"]
    P_hgcer_yAtCer_vs_MM_RAND = subDict["P_hgcer_yAtCer_vs_MM_SUB_RAND"]
    P_hgcer_nohole_yAtCer_vs_MM_RAND = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND"]

    MM_vs_CoinTime_DUMMY_RAND = subDict["MM_vs_CoinTime_SUB_DUMMY_RAND"]
    CoinTime_vs_beta_DUMMY_RAND = subDict["CoinTime_vs_beta_SUB_DUMMY_RAND"]
    MM_vs_beta_DUMMY_RAND = subDict["MM_vs_beta_SUB_DUMMY_RAND"]
    MM_vs_H_cer_DUMMY_RAND = subDict["MM_vs_H_cer_SUB_DUMMY_RAND"]
    MM_vs_H_cal_DUMMY_RAND = subDict["MM_vs_H_cal_SUB_DUMMY_RAND"]
    MM_vs_P_cal_DUMMY_RAND = subDict["MM_vs_P_cal_SUB_DUMMY_RAND"]
    MM_vs_P_hgcer_DUMMY_RAND = subDict["MM_vs_P_hgcer_SUB_DUMMY_RAND"]
    MM_vs_P_aero_DUMMY_RAND = subDict["MM_vs_P_aero_SUB_DUMMY_RAND"]
    phiq_vs_t_DUMMY_RAND = subDict["phiq_vs_t_SUB_DUMMY_RAND"]
    Q2_vs_W_DUMMY_RAND = subDict["Q2_vs_W_SUB_DUMMY_RAND"]
    Q2_vs_t_DUMMY_RAND = subDict["Q2_vs_t_SUB_DUMMY_RAND"]
    W_vs_t_DUMMY_RAND = subDict["W_vs_t_SUB_DUMMY_RAND"]
    EPS_vs_t_DUMMY_RAND = subDict["EPS_vs_t_SUB_DUMMY_RAND"]
    MM_vs_t_DUMMY_RAND = subDict["MM_vs_t_SUB_DUMMY_RAND"]
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"]
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND"]    
    
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
    # Fill histograms for various trees called above

    _print_sub_timer("particle_subtraction setup {}".format(phi_setting), perf_counter() - setup_start)

    hole_contains = hgcer_cutg.IsInside if ParticleType == "kaon" else None

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
        "nomm": H_MM_nosub_DATA.Fill,
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

    _process_particle_subtraction_tree(
        TBRANCH_DATA,
        "\nGrabbing {} {} subtraction data...".format(phi_setting, SubtractedParticle),
        "particle_subtraction data loop {}".format(phi_setting),
        MM_offset_DATA,
        data_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        get_shifted_mm,
        get_shifted_t,
        get_effective_mm_offset,
        mm_min,
        mm_max,
        Misc.progressBar,
    )

    ################################################################################################################################################
    # Fill histograms for various trees called above

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
        "nomm": H_MM_nosub_DUMMY.Fill,
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
        "h_cal": H_cal_etottracknorm_DUMMY.Fill,
        "h_cer": H_cer_npeSum_DUMMY.Fill,
        "p_cal": P_cal_etottracknorm_DUMMY.Fill,
        "p_hgcer": P_hgcer_npeSum_DUMMY.Fill,
        "p_aero": P_aero_npeSum_DUMMY.Fill,
    }

    _process_particle_subtraction_tree(
        TBRANCH_DUMMY,
        "\nGrabbing {} {} subtraction dummy...".format(phi_setting, SubtractedParticle),
        "particle_subtraction dummy loop {}".format(phi_setting),
        MM_offset_DATA,
        dummy_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        get_shifted_mm,
        get_shifted_t,
        get_effective_mm_offset,
        mm_min,
        mm_max,
        Misc.progressBar,
    )

    ################################################################################################################################################
    # Fill histograms for various trees called above

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
        "nomm": H_MM_nosub_RAND.Fill,
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
        "h_cal": H_cal_etottracknorm_RAND.Fill,
        "h_cer": H_cer_npeSum_RAND.Fill,
        "p_cal": P_cal_etottracknorm_RAND.Fill,
        "p_hgcer": P_hgcer_npeSum_RAND.Fill,
        "p_aero": P_aero_npeSum_RAND.Fill,
    }

    _process_particle_subtraction_tree(
        TBRANCH_RAND,
        "\nGrabbing {} {} subtraction random...".format(phi_setting, SubtractedParticle),
        "particle_subtraction random loop {}".format(phi_setting),
        MM_offset_DATA,
        rand_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        get_shifted_mm,
        get_shifted_t,
        get_effective_mm_offset,
        mm_min,
        mm_max,
        Misc.progressBar,
    )
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

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
        "nomm": H_MM_nosub_DUMMY_RAND.Fill,
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
        "h_cal": H_cal_etottracknorm_DUMMY_RAND.Fill,
        "h_cer": H_cer_npeSum_DUMMY_RAND.Fill,
        "p_cal": P_cal_etottracknorm_DUMMY_RAND.Fill,
        "p_hgcer": P_hgcer_npeSum_DUMMY_RAND.Fill,
        "p_aero": P_aero_npeSum_DUMMY_RAND.Fill,
    }

    _process_particle_subtraction_tree(
        TBRANCH_DUMMY_RAND,
        "\nGrabbing {} {} subtraction dummy random...".format(phi_setting, SubtractedParticle),
        "particle_subtraction dummy random loop {}".format(phi_setting),
        MM_offset_DATA,
        dummy_rand_fills,
        ParticleType,
        hole_contains,
        evaluate_data_event,
        get_shifted_mm,
        get_shifted_t,
        get_effective_mm_offset,
        mm_min,
        mm_max,
        Misc.progressBar,
        kaon_nomm_enabled=False,
    )

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
    H_t_DUMMY_RAND.Scale(1/nWindows)
    H_epsilon_DUMMY_RAND.Scale(1/nWindows)
    H_MM_DUMMY_RAND.Scale(1/nWindows)
    H_MM_nosub_DUMMY_RAND.Scale(1/nWindows)    
    H_pmiss_DUMMY_RAND.Scale(1/nWindows)
    H_emiss_DUMMY_RAND.Scale(1/nWindows)
    H_pmx_DUMMY_RAND.Scale(1/nWindows)
    H_pmy_DUMMY_RAND.Scale(1/nWindows)
    H_pmz_DUMMY_RAND.Scale(1/nWindows)
    H_W_DUMMY_RAND.Scale(1/nWindows)
    #H_ct_DUMMY_RAND.Scale(1/nWindows)
    
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
    H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
    H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)
    H_MM_DUMMY.Add(H_MM_DUMMY_RAND,-1)
    H_MM_nosub_DUMMY.Add(H_MM_nosub_DUMMY_RAND,-1)    
    H_pmiss_DUMMY.Add(H_pmiss_DUMMY_RAND,-1)
    H_emiss_DUMMY.Add(H_emiss_DUMMY_RAND,-1)
    H_pmx_DUMMY.Add(H_pmx_DUMMY_RAND,-1)
    H_pmy_DUMMY.Add(H_pmy_DUMMY_RAND,-1)
    H_pmz_DUMMY.Add(H_pmz_DUMMY_RAND,-1)
    H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
    H_ct_DUMMY.Add(H_ct_DUMMY_RAND,-1)
    _print_sub_timer("particle_subtraction random-window subtraction {}".format(phi_setting), perf_counter() - stage_start)

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
    H_MM_nosub_DUMMY.Scale(norm_factor_dummy)
    H_pmiss_DUMMY.Scale(norm_factor_dummy)
    H_emiss_DUMMY.Scale(norm_factor_dummy)
    H_pmx_DUMMY.Scale(norm_factor_dummy)
    H_pmy_DUMMY.Scale(norm_factor_dummy)
    H_pmz_DUMMY.Scale(norm_factor_dummy)
    H_ct_DUMMY.Scale(norm_factor_dummy)

    ###
    # Data dummy subtraction
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
    H_MM_DATA.Add(H_MM_DUMMY,-1)
    H_MM_nosub_DATA.Add(H_MM_nosub_DUMMY,-1)    
    H_pmiss_DATA.Add(H_pmiss_DUMMY,-1)
    H_emiss_DATA.Add(H_emiss_DUMMY,-1)
    H_pmx_DATA.Add(H_pmx_DUMMY,-1)
    H_pmy_DATA.Add(H_pmy_DUMMY,-1)
    H_pmz_DATA.Add(H_pmz_DUMMY,-1)
    H_ct_DATA.Add(H_ct_DUMMY,-1)
    _print_sub_timer("particle_subtraction norm/dummy subtraction {}".format(phi_setting), perf_counter() - stage_start)
    _print_sub_timer("particle_subtraction total {}".format(phi_setting), perf_counter() - total_start)

################################################################################################################################################

def particle_subtraction_ave(t_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg=None):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    norm_factor_data = inpDict["normfac_data"]
    norm_factor_dummy = inpDict["normfac_dummy"]   

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, get_effective_mm_offset, get_shifted_mm, get_shifted_t, set_shift_context, set_val
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=inpDict.get("shift_mode", "raw"))
    mm_offset_for_cuts = get_effective_mm_offset(MM_offset_DATA)
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    prompt_tree_name = get_prompt_tree_name(SubtractedParticle, EPSSET)
    rand_tree_name = get_rand_tree_name(SubtractedParticle, EPSSET)

    TBRANCH_DATA  = InFile_DATA.Get(prompt_tree_name)

    TBRANCH_RAND  = InFile_DATA.Get(rand_tree_name)

    ################################################################################################################################################
    # Defin of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get(prompt_tree_name)
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get(rand_tree_name)

    ################################################################################################################################################

    hist_dict = {}
    
    for j in range(len(t_bins)-1):
        
        hist_dict["H_Q2_DATA_{}".format(j)] = subDict["H_Q2_SUB_DATA_{}".format(j)]
        hist_dict["H_W_DATA_{}".format(j)] = subDict["H_W_SUB_DATA_{}".format(j)]
        hist_dict["H_t_DATA_{}".format(j)] = subDict["H_t_SUB_DATA_{}".format(j)]
        hist_dict["H_epsilon_DATA_{}".format(j)] = subDict["H_epsilon_SUB_DATA_{}".format(j)]
        hist_dict["H_MM_DATA_{}".format(j)] = subDict["H_MM_SUB_DATA_{}".format(j)]
        hist_dict["H_MM_nosub_DATA_{}".format(j)] = subDict["H_MM_nosub_SUB_DATA_{}".format(j)]        

        hist_dict["H_Q2_DUMMY_{}".format(j)] = subDict["H_Q2_SUB_DUMMY_{}".format(j)]
        hist_dict["H_W_DUMMY_{}".format(j)] = subDict["H_W_SUB_DUMMY_{}".format(j)]
        hist_dict["H_t_DUMMY_{}".format(j)] = subDict["H_t_SUB_DUMMY_{}".format(j)]
        hist_dict["H_epsilon_DUMMY_{}".format(j)] = subDict["H_epsilon_SUB_DUMMY_{}".format(j)]
        hist_dict["H_MM_DUMMY_{}".format(j)] = subDict["H_MM_SUB_DUMMY_{}".format(j)]
        hist_dict["H_MM_nosub_DUMMY_{}".format(j)] = subDict["H_MM_nosub_SUB_DUMMY_{}".format(j)]        

        hist_dict["H_Q2_RAND_{}".format(j)] = subDict["H_Q2_SUB_RAND_{}".format(j)]
        hist_dict["H_W_RAND_{}".format(j)] = subDict["H_W_SUB_RAND_{}".format(j)]
        hist_dict["H_t_RAND_{}".format(j)] = subDict["H_t_SUB_RAND_{}".format(j)]
        hist_dict["H_epsilon_RAND_{}".format(j)] = subDict["H_epsilon_SUB_RAND_{}".format(j)]
        hist_dict["H_MM_RAND_{}".format(j)] = subDict["H_MM_SUB_RAND_{}".format(j)]
        hist_dict["H_MM_nosub_RAND_{}".format(j)] = subDict["H_MM_nosub_SUB_RAND_{}".format(j)]        

        hist_dict["H_Q2_DUMMY_RAND_{}".format(j)] = subDict["H_Q2_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_W_DUMMY_RAND_{}".format(j)] = subDict["H_W_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_t_DUMMY_RAND_{}".format(j)] = subDict["H_t_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)] = subDict["H_epsilon_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_MM_DUMMY_RAND_{}".format(j)] = subDict["H_MM_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)] = subDict["H_MM_nosub_SUB_DUMMY_RAND_{}".format(j)]        
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1.0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

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
            
    ##############
    ##############        
    ##############

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + MM_offset_DATA

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_MM_nosub_DATA_{}".format(j)].Fill(adj_MM)            
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_Q2_DATA_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DATA_{}".format(j)].Fill(adj_t)
                    hist_dict["H_W_DATA_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DATA_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DATA_{}".format(j)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_MM_nosub_DUMMY_{}".format(j)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_Q2_DUMMY_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DUMMY_{}".format(j)].Fill(adj_t)
                    hist_dict["H_W_DUMMY_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DUMMY_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DUMMY_{}".format(j)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_MM_nosub_RAND_{}".format(j)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_Q2_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_RAND_{}".format(j)].Fill(adj_t)
                    hist_dict["H_W_RAND_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_RAND_{}".format(j)].Fill(adj_MM)
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:                
                    hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Fill(adj_MM)                                
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= adj_t < t_bins[j+1]:
                    hist_dict["H_Q2_DUMMY_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DUMMY_RAND_{}".format(j)].Fill(adj_t)
                    hist_dict["H_W_DUMMY_RAND_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DUMMY_RAND_{}".format(j)].Fill(adj_MM)
                  
    for j in range(len(t_bins)-1):

        # Data Random subtraction window
        hist_dict["H_Q2_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_W_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_dict["H_t_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_epsilon_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_nosub_RAND_{}".format(j)].Scale(1/nWindows)        

        # Data Dummy_Random subtraction window
        hist_dict["H_Q2_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_W_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_t_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)        

        ###
        # Data Random subtraction
        hist_dict["H_Q2_DATA_{}".format(j)].Add(hist_dict["H_Q2_RAND_{}".format(j)],-1)
        hist_dict["H_W_DATA_{}".format(j)].Add(hist_dict["H_W_RAND_{}".format(j)],-1)
        hist_dict["H_t_DATA_{}".format(j)].Add(hist_dict["H_t_RAND_{}".format(j)],-1)
        hist_dict["H_epsilon_DATA_{}".format(j)].Add(hist_dict["H_epsilon_RAND_{}".format(j)],-1)
        hist_dict["H_MM_DATA_{}".format(j)].Add(hist_dict["H_MM_RAND_{}".format(j)],-1)
        hist_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_dict["H_MM_nosub_RAND_{}".format(j)],-1)        

        ###
        # Dummy Random subtraction
        hist_dict["H_Q2_DUMMY_{}".format(j)].Add(hist_dict["H_Q2_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_W_DUMMY_{}".format(j)].Add(hist_dict["H_W_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_t_DUMMY_{}".format(j)].Add(hist_dict["H_t_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_epsilon_DUMMY_{}".format(j)].Add(hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_MM_DUMMY_{}".format(j)].Add(hist_dict["H_MM_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_MM_nosub_DUMMY_{}".format(j)].Add(hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)],-1)    

        ###
        # Data Normalization
        hist_dict["H_Q2_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_dict["H_W_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_dict["H_t_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_dict["H_epsilon_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_dict["H_MM_DATA_{}".format(j)].Scale(norm_factor_data)
        hist_dict["H_MM_nosub_DATA_{}".format(j)].Scale(norm_factor_data)

        ###
        # Dummy Normalization
        hist_dict["H_Q2_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_dict["H_W_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_dict["H_t_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_dict["H_epsilon_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_dict["H_MM_DUMMY_{}".format(j)].Scale(norm_factor_dummy)
        hist_dict["H_MM_nosub_DUMMY_{}".format(j)].Scale(norm_factor_dummy)

        ###
        # Data dummy subtraction
        hist_dict["H_Q2_DATA_{}".format(j)].Add(hist_dict["H_Q2_DUMMY_{}".format(j)],-1)
        hist_dict["H_W_DATA_{}".format(j)].Add(hist_dict["H_W_DUMMY_{}".format(j)],-1)
        hist_dict["H_t_DATA_{}".format(j)].Add(hist_dict["H_t_DUMMY_{}".format(j)],-1)
        hist_dict["H_epsilon_DATA_{}".format(j)].Add(hist_dict["H_epsilon_DUMMY_{}".format(j)],-1)
        hist_dict["H_MM_DATA_{}".format(j)].Add(hist_dict["H_MM_DUMMY_{}".format(j)],-1)
        hist_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_dict["H_MM_nosub_DUMMY_{}".format(j)],-1)

################################################################################################################################################

def particle_subtraction_yield(t_bins, phi_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg=None):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]
    POL = inpDict["POL"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    norm_factor_data = inpDict["normfac_data"]
    norm_factor_dummy = inpDict["normfac_dummy"]

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, get_effective_mm_offset, get_shifted_mm, get_shifted_t, set_shift_context, set_val
    sys.path.append("binning")
    from theta_cm import calculate_theta_cm_deg, calculate_tmin
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=inpDict.get("shift_mode", "raw"))
    mm_offset_for_cuts = get_effective_mm_offset(MM_offset_DATA)
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    prompt_tree_name = get_prompt_tree_name(SubtractedParticle, EPSSET)
    rand_tree_name = get_rand_tree_name(SubtractedParticle, EPSSET)

    TBRANCH_DATA  = InFile_DATA.Get(prompt_tree_name)

    TBRANCH_RAND  = InFile_DATA.Get(rand_tree_name)

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get(prompt_tree_name)
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get(rand_tree_name)

    ################################################################################################################################################

    hist_dict = {}
    sub_event_cache = _init_sub_event_cache()
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
    
            hist_dict["H_Q2_DATA_{}_{}".format(j, k)] = subDict["H_Q2_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_W_DATA_{}_{}".format(j, k)] = subDict["H_W_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)] = subDict["H_Q2_vs_W_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_theta_cm_DATA_{}_{}".format(j, k)] = subDict["H_theta_cm_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_t_DATA_{}_{}".format(j, k)] = subDict["H_t_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_hsxptar_DATA_{}_{}".format(j, k)] = subDict["H_hsxptar_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_hsyptar_DATA_{}_{}".format(j, k)] = subDict["H_hsyptar_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_ssxptar_DATA_{}_{}".format(j, k)] = subDict["H_ssxptar_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_ssyptar_DATA_{}_{}".format(j, k)] = subDict["H_ssyptar_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)] = subDict["H_t_vs_tmin_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_MM_DATA_{}_{}".format(j, k)] = subDict["H_MM_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)]            

            hist_dict["H_Q2_DUMMY_{}_{}".format(j, k)] = subDict["H_Q2_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_W_DUMMY_{}_{}".format(j, k)] = subDict["H_W_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)] = subDict["H_Q2_vs_W_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)] = subDict["H_theta_cm_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_t_DUMMY_{}_{}".format(j, k)] = subDict["H_t_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)] = subDict["H_hsxptar_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)] = subDict["H_hsyptar_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)] = subDict["H_ssxptar_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)] = subDict["H_ssyptar_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)] = subDict["H_t_vs_tmin_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)] = subDict["H_MM_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k)]            

            hist_dict["H_Q2_RAND_{}_{}".format(j, k)] = subDict["H_Q2_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_W_RAND_{}_{}".format(j, k)] = subDict["H_W_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)] = subDict["H_Q2_vs_W_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_theta_cm_RAND_{}_{}".format(j, k)] = subDict["H_theta_cm_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_t_RAND_{}_{}".format(j, k)] = subDict["H_t_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_hsxptar_RAND_{}_{}".format(j, k)] = subDict["H_hsxptar_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_hsyptar_RAND_{}_{}".format(j, k)] = subDict["H_hsyptar_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_ssxptar_RAND_{}_{}".format(j, k)] = subDict["H_ssxptar_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_ssyptar_RAND_{}_{}".format(j, k)] = subDict["H_ssyptar_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)] = subDict["H_t_vs_tmin_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_RAND_{}_{}".format(j, k)] = subDict["H_MM_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_RAND_{}_{}".format(j, k)]            

            hist_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_Q2_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_W_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_Q2_vs_W_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_theta_cm_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_t_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_hsxptar_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_hsyptar_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_ssxptar_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_ssyptar_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_t_vs_tmin_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k)]            
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1.0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

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
            
    ##############
    ##############        
    ##############

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + MM_offset_DATA

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        theta_cm_deg = calculate_theta_cm_deg(ParticleType, POL, evt.W, evt.Q2, adj_t)
        minus_tmin = calculate_tmin(ParticleType, POL, evt.W, evt.Q2)
        
        ##############
        ##############        
        ##############
        
        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):            
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            _append_sub_event(
                                sub_event_cache["prompt"],
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
                                True,
                                NOMMCUTS,
                                t_index=j,
                                phi_index=k,
                            )
                            hist_dict["H_Q2_DATA_{}_{}".format(j, k)].Fill(evt.Q2)
                            hist_dict["H_W_DATA_{}_{}".format(j, k)].Fill(evt.W)
                            hist_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Fill(evt.Q2, evt.W)
                            if math.isfinite(theta_cm_deg):
                                hist_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Fill(theta_cm_deg)
                            hist_dict["H_t_DATA_{}_{}".format(j, k)].Fill(adj_t)
                            hist_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Fill(evt.hsxptar)
                            hist_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Fill(evt.hsyptar)
                            hist_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Fill(evt.ssxptar)
                            hist_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Fill(evt.ssyptar)
                            if math.isfinite(minus_tmin):
                                hist_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Fill(minus_tmin, adj_t)
                            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        theta_cm_deg = calculate_theta_cm_deg(ParticleType, POL, evt.W, evt.Q2, adj_t)
        minus_tmin = calculate_tmin(ParticleType, POL, evt.W, evt.Q2)

        ##############
        ##############        
        ##############

        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:                
                            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            _append_sub_event(
                                sub_event_cache["dummy"],
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
                                True,
                                NOMMCUTS,
                                t_index=j,
                                phi_index=k,
                            )
                            hist_dict["H_Q2_DUMMY_{}_{}".format(j, k)].Fill(evt.Q2)
                            hist_dict["H_W_DUMMY_{}_{}".format(j, k)].Fill(evt.W)
                            hist_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)].Fill(evt.Q2, evt.W)
                            if math.isfinite(theta_cm_deg):
                                hist_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)].Fill(theta_cm_deg)
                            hist_dict["H_t_DUMMY_{}_{}".format(j, k)].Fill(adj_t)
                            hist_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)].Fill(evt.hsxptar)
                            hist_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)].Fill(evt.hsyptar)
                            hist_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)].Fill(evt.ssxptar)
                            hist_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)].Fill(evt.ssyptar)
                            if math.isfinite(minus_tmin):
                                hist_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)].Fill(minus_tmin, adj_t)
                            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        theta_cm_deg = calculate_theta_cm_deg(ParticleType, POL, evt.W, evt.Q2, adj_t)
        minus_tmin = calculate_tmin(ParticleType, POL, evt.W, evt.Q2)

        ##############
        ##############        
        ##############

        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:                
                            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            _append_sub_event(
                                sub_event_cache["rand"],
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
                                True,
                                NOMMCUTS,
                                t_index=j,
                                phi_index=k,
                            )
                            hist_dict["H_Q2_RAND_{}_{}".format(j, k)].Fill(evt.Q2)
                            hist_dict["H_W_RAND_{}_{}".format(j, k)].Fill(evt.W)
                            hist_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)].Fill(evt.Q2, evt.W)
                            if math.isfinite(theta_cm_deg):
                                hist_dict["H_theta_cm_RAND_{}_{}".format(j, k)].Fill(theta_cm_deg)
                            hist_dict["H_t_RAND_{}_{}".format(j, k)].Fill(adj_t)
                            hist_dict["H_hsxptar_RAND_{}_{}".format(j, k)].Fill(evt.hsxptar)
                            hist_dict["H_hsyptar_RAND_{}_{}".format(j, k)].Fill(evt.hsyptar)
                            hist_dict["H_ssxptar_RAND_{}_{}".format(j, k)].Fill(evt.ssxptar)
                            hist_dict["H_ssyptar_RAND_{}_{}".format(j, k)].Fill(evt.ssyptar)
                            if math.isfinite(minus_tmin):
                                hist_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)].Fill(minus_tmin, adj_t)
                            hist_dict["H_MM_RAND_{}_{}".format(j, k)].Fill(adj_MM)
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = get_shifted_mm(evt, mm_offset=MM_offset_DATA)
        adj_t = get_shifted_t(evt)
        theta_cm_deg = calculate_theta_cm_deg(ParticleType, POL, evt.W, evt.Q2, adj_t)
        minus_tmin = calculate_tmin(ParticleType, POL, evt.W, evt.Q2)

        ##############
        ##############        
        ##############


        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)
            NOMMCUTS = apply_data_sub_cuts(evt, mm_min, mm_max, mm_offset=mm_offset_for_cuts)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= adj_t < t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) < phi_bins[k+1]:
                            _append_sub_event(
                                sub_event_cache["dummy_rand"],
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
                                True,
                                NOMMCUTS,
                                t_index=j,
                                phi_index=k,
                            )
                            hist_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.Q2)
                            hist_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.W)
                            hist_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.Q2, evt.W)
                            if math.isfinite(theta_cm_deg):
                                hist_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)].Fill(theta_cm_deg)
                            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_t)
                            hist_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.hsxptar)
                            hist_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.hsyptar)
                            hist_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.ssxptar)
                            hist_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.ssyptar)
                            if math.isfinite(minus_tmin):
                                hist_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)].Fill(minus_tmin, adj_t)
                            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)

    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            
            # Data Random subtraction window    
            hist_dict["H_Q2_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_W_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_theta_cm_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_t_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_hsxptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_hsyptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_ssxptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_ssyptar_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            

            # Data Dummy_Random subtraction window
            hist_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            

            ###
            # Data Random subtraction
            hist_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(hist_dict["H_Q2_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_W_DATA_{}_{}".format(j, k)].Add(hist_dict["H_W_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(hist_dict["H_Q2_vs_W_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(hist_dict["H_theta_cm_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_dict["H_t_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_hsxptar_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_hsyptar_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_ssxptar_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_ssyptar_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(hist_dict["H_t_vs_tmin_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)],-1)            

            ###
            # Dummy Random subtraction
            hist_dict["H_Q2_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_Q2_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_W_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_W_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_Q2_vs_W_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_theta_cm_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_t_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_hsxptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_hsyptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_ssxptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_ssyptar_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_t_vs_tmin_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)],-1)    

            ###
            # Data Normalization
            hist_dict["H_Q2_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_W_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_t_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Scale(norm_factor_data)

            ###
            # Dummy Normalization
            hist_dict["H_Q2_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_W_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_t_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Scale(norm_factor_dummy)
                                   
            ###
            # Data Dummy subtraction
            hist_dict["H_Q2_DATA_{}_{}".format(j, k)].Add(hist_dict["H_Q2_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_W_DATA_{}_{}".format(j, k)].Add(hist_dict["H_W_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_Q2_vs_W_DATA_{}_{}".format(j, k)].Add(hist_dict["H_Q2_vs_W_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_theta_cm_DATA_{}_{}".format(j, k)].Add(hist_dict["H_theta_cm_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_dict["H_t_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_hsxptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_hsxptar_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_hsyptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_hsyptar_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_ssxptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_ssxptar_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_ssyptar_DATA_{}_{}".format(j, k)].Add(hist_dict["H_ssyptar_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_t_vs_tmin_DATA_{}_{}".format(j, k)].Add(hist_dict["H_t_vs_tmin_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_DUMMY_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)],-1)

    subDict["_sub_event_cache"] = _freeze_sub_event_cache(sub_event_cache)
