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
import sys, math, os, subprocess
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
from utility import open_root_file, remove_bad_bins, create_polar_plot, integrate_hist_range
from mm_background_subtraction import (
    build_mm_background_weights,
    build_mm_residual_weights,
    clone_reset_hist,
    mm_background_weight_from_value,
)
from apply_cuts import get_shift_mode, get_shifted_mm, get_shifted_t, set_shift_context

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
    hists["hgcer_x_mm"].Fill(evt.P_hgcer_xAtCer, adj_MM, weight * evt.P_hgcer_npeSum)
    hists["hgcer_y_mm"].Fill(evt.P_hgcer_yAtCer, adj_MM, weight * evt.P_hgcer_npeSum)

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
    hists["h_cal"].Fill(evt.H_cal_etottracknorm, weight)
    hists["h_cer"].Fill(evt.H_cer_npeSum, weight)
    hists["p_cal"].Fill(evt.P_cal_etottracknorm, weight)
    hists["p_hgcer"].Fill(evt.P_hgcer_npeSum, weight)
    hists["p_aero"].Fill(evt.P_aero_npeSum, weight)


def _create_rand_sub_bg_templates(target_hists):
    return {
        key: clone_reset_hist(hist, "_bg_template")
        for key, hist in target_hists.items()
        if hist is not None
    }


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
):
    for evt in tree:
        if particle_type == "kaon":
            base_all_cuts, _, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            allcuts = base_all_cuts and not hole_rejected
        else:
            base_all_cuts, _, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
            allcuts = base_all_cuts

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
    residual_weights=None,
):
    mm_offset_correction = 0.0 if get_shift_mode() == "shifted" else mm_offset_data

    for evt in tree:
        base_all_cuts, _, adj_hsdelta = evaluate_event(evt, mm_min, mm_max, mm_offset=mm_offset_correction)

        if particle_type == "kaon":
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            pid_pass = evt.P_hgcer_npeSum > 2.0
            allcuts = base_all_cuts and not hole_rejected and pid_pass
        else:
            allcuts = base_all_cuts

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
        if event_weight == 0.0:
            continue

        _fill_rand_sub_allcuts_weighted(evt, adj_MM, adj_t, adj_hsdelta, event_weight, template_hists)

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
    update_mm_offset=False,
):
    print(print_label)
    entries = tree.GetEntries()
    progress_time = 0.0
    loop_start = perf_counter()
    nohole_xy_fill = fills["nohole_xy"]
    nohole_x_mm_fill = fills["nohole_x_mm"]
    nohole_y_mm_fill = fills["nohole_y_mm"]
    mm_offset_value = None

    for i, evt in enumerate(tree):
        progress_start = perf_counter()
        progress_bar(i, entries, bar_length=25)
        progress_time += perf_counter() - progress_start

        if particle_type == "kaon":
            base_all_cuts, base_sub_cuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
            hole_rejected = hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
            allcuts = base_all_cuts and not hole_rejected
            noholecuts = base_all_cuts
            nommcuts = base_sub_cuts and not hole_rejected
        else:
            base_all_cuts, base_sub_cuts, adj_hsdelta = evaluate_event(evt, mm_min, mm_max)
            allcuts = base_all_cuts
            nommcuts = base_sub_cuts
            noholecuts = False

        if not (noholecuts or nommcuts or allcuts):
            continue

        adj_MM = get_shifted_mm(evt)
        adj_t = get_shifted_t(evt)

        if noholecuts and nohole_xy_fill is not None:
            nohole_xy_fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
            nohole_x_mm_fill(evt.P_hgcer_xAtCer, adj_MM, evt.P_hgcer_npeSum)
            nohole_y_mm_fill(evt.P_hgcer_yAtCer, adj_MM, evt.P_hgcer_npeSum)

        if nommcuts:
            for fill in nomm_fills:
                fill(adj_MM)

        if allcuts:
            _fill_rand_sub_allcuts(evt, adj_MM, adj_t, adj_hsdelta, fills)
            if update_mm_offset:
                mm_offset_value = adj_MM - evt.MM

    loop_elapsed = perf_counter() - loop_start
    _print_rand_timer(timer_label, loop_elapsed, entries)
    _print_rand_timer("{} progressBar".format(timer_label), progress_time, entries)
    _print_rand_timer("{} other".format(timer_label), max(loop_elapsed - progress_time, 0.0), entries)
    return mm_offset_value

def rand_sub(phi_setting, inpDict, shift_mode="raw", emit_plots=True):
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
        get_shifted_mm,
        get_shifted_t,
        set_shift_context,
        set_val,
    )
    set_val(inpDict) # Set global variables for optimization
    set_shift_context(phi_setting=phi_setting, shift_mode=shift_mode)
    
    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        histDict.update({ "phi_setting" : phi_setting})
        return histDict

    InFile_DATA = open_root_file(rootFileData)

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        return histDict

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))

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
    H_MM_rand_dummy_DATA  = TH1D("H_MM_rand_dummy_DATA",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_dummy_DATA  = TH1D("H_MM_dummy_DATA",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_full_DATA  = TH1D("H_MM_full_DATA",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit2sub_DATA  = TH1D("H_MM_fit2sub_DATA",f"MM_fit2sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit1sub_DATA  = TH1D("H_MM_fit1sub_DATA",f"MM_fit1sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_pisub_DATA  = TH1D("H_MM_pisub_DATA",f"MM_pisub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_nosub_DATA  = TH1D("H_MM_nosub_DATA",f"MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
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
    H_MM_full_DUMMY  = TH1D("H_MM_full_DUMMY",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit2sub_DUMMY  = TH1D("H_MM_fit2sub_DUMMY",f"MM_fit2sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit1sub_DUMMY  = TH1D("H_MM_fit1sub_DUMMY",f"MM_fit1sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_pisub_DUMMY  = TH1D("H_MM_pisub_DUMMY",f"MM_pisub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_nosub_DUMMY  = TH1D("H_MM_nosub_DUMMY",f"MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
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
    H_MM_rand_dummy_RAND  = TH1D("H_MM_rand_dummy_RAND",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_dummy_RAND  = TH1D("H_MM_dummy_RAND",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_full_RAND  = TH1D("H_MM_full_RAND",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)    
    H_MM_fit2sub_RAND  = TH1D("H_MM_fit2sub_RAND",f"MM_fit2sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit1sub_RAND  = TH1D("H_MM_fit1sub_RAND",f"MM_fit1sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_pisub_RAND  = TH1D("H_MM_pisub_RAND",f"MM_pisub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_nosub_RAND  = TH1D("H_MM_nosub_RAND",f"MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
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
    H_MM_full_DUMMY_RAND  = TH1D("H_MM_full_DUMMY_RAND",f"MM_{ParticleType[0].upper()}", 100, 0.7, 1.5)    
    H_MM_fit2sub_DUMMY_RAND  = TH1D("H_MM_fit2sub_DUMMY_RAND",f"MM_fit2sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_fit1sub_DUMMY_RAND  = TH1D("H_MM_fit1sub_DUMMY_RAND",f"MM_fit1sub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_pisub_DUMMY_RAND  = TH1D("H_MM_pisub_DUMMY_RAND",f"MM_pisub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_MM_nosub_DUMMY_RAND  = TH1D("H_MM_nosub_DUMMY_RAND",f"MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
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
        subDict["H_MM_full_SUB_DATA"]  = TH1D("H_MM_full_SUB_DATA",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)        
        subDict["H_MM_nosub_SUB_DATA"]  = TH1D("H_MM_nosub_SUB_DATA",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)
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
        subDict["H_MM_full_SUB_RAND"]  = TH1D("H_MM_full_SUB_RAND",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)        
        subDict["H_MM_nosub_SUB_RAND"]  = TH1D("H_MM_nosub_SUB_RAND",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)
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
        subDict["H_MM_full_SUB_DUMMY"]  = TH1D("H_MM_full_SUB_DUMMY",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)
        subDict["H_MM_nosub_SUB_DUMMY"]  = TH1D("H_MM_nosub_SUB_DUMMY",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)
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
        subDict["H_MM_full_SUB_DUMMY_RAND"]  = TH1D("H_MM_full_SUB_DUMMY_RAND",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)        
        subDict["H_MM_nosub_SUB_DUMMY_RAND"]  = TH1D("H_MM_nosub_SUB_DUMMY_RAND",f"MM_{SubtractedParticle}", 100, 0.7, 1.5)
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

    MM_offset_DATA = _process_rand_sub_tree(
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
        update_mm_offset=True,
    )

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
    H_MM_dummy_DATA.Add(H_MM_RAND,-1)
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

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        stage_start = perf_counter()
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        subDict["MM_offset_DATA"] = MM_offset_DATA
        particle_subtraction_cuts(histDict, subDict, inpDict, SubtractedParticle, hgcer_cutg)
        
        try:
            ##############
            # HARD CODED #
            ##############
            pi_mm_min = 0.88 + MM_offset_DATA
            pi_mm_max = 0.94 + MM_offset_DATA 
            ###pi_mm_min = 0.91 + MM_offset_DATA
            ###pi_mm_max = 0.98 + MM_offset_DATA                       
            # Scale pion to kaon data
            kaon_amp = integrate_hist_range(
                H_MM_nosub_DATA,
                pi_mm_min, pi_mm_max
            )

            pion_background_amp = integrate_hist_range(
                subDict["H_MM_nosub_SUB_DATA"],
                pi_mm_min, pi_mm_max
            )
            scale_factor = (kaon_amp / pion_background_amp)
            ##############
            ##############
            ##############
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

    if ParticleType == "kaon":
        sub_root_data = open_root_file(f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root")
        sub_root_dummy = open_root_file(f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root")
        TBRANCH_SUB_DATA = sub_root_data.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))
        TBRANCH_SUB_RAND = sub_root_data.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))
        TBRANCH_SUB_DUMMY = sub_root_dummy.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))
        TBRANCH_SUB_DUMMY_RAND = sub_root_dummy.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    # Fit background and subtract
    # --------------------------------------------------------------
    # Stat‑scale: events that survive ALL subtractions & MM‑cuts
    # --------------------------------------------------------------
    inpDict["bg_stat_scale1"] = 0.50
    residual_bg_weights1 = None

    if  inpDict["bg_stat_scale1"] > 0.0:
        background_fit1 = bg_fit(phi_setting,
                                inpDict,
                                H_MM_pisub_DATA,   # wide / no-cut
                                H_MM_DATA,         # cut-window axis
                                scaling=inpDict["bg_stat_scale1"],
                                model_key=f"fixquad_{phi_setting}_{EPSSET}e",
                                fit_name="Fit 1")
        # background_fit1[0] : scaled function   (use for subtraction)
        # background_fit1[1] : original function (use for drawing only)
        mm_stage1_input = clone_reset_hist(H_MM_DATA, "_stage1_input")
        mm_stage1_input.Add(H_MM_DATA)
        bg_templates1 = _create_rand_sub_bg_templates(data_bg_targets)
        bg_weights1 = build_mm_background_weights(mm_stage1_input, background_fit1[0])

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
        )
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
        )
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
        )
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
        )

        if ParticleType == "kaon" and scale_factor != 0.0:
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
            )
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
            )
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
            )
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
            )

        for key, hist in data_bg_targets.items():
            hist.Add(bg_templates1[key], -1)

        residual_bg_weights1 = build_mm_residual_weights(bg_weights1)
        H_MM_fit2sub_DATA.Add(background_fit1[1], -1)
        H_MM_fit1sub_DATA.Add(background_fit1[1], -1)
        H_MM_DATA.Add(background_fit1[0], -1)
        H_MM_full_DATA.Add(background_fit1[1], -1)        

    # Fit background and subtract
    # --------------------------------------------------------------
    # Stat‑scale: events that survive ALL subtractions & MM‑cuts
    # --------------------------------------------------------------
    inpDict["bg_stat_scale2"] = 0.50

    if inpDict["bg_stat_scale2"] > 0.0:
        background_fit2 = bg_fit(phi_setting,
                                inpDict,
                                H_MM_fit1sub_DATA,   # wide / no-cut
                                H_MM_DATA,         # cut-window axis
                                scaling=inpDict["bg_stat_scale2"],
                                model_key=f"cheb2_{phi_setting}_{EPSSET}e",
                                fit_name="Fit 2")
        # background_fit2[0] : scaled function   (use for subtraction)
        # background_fit2[1] : original function (use for drawing only)
        mm_stage2_input = clone_reset_hist(H_MM_DATA, "_stage2_input")
        mm_stage2_input.Add(H_MM_DATA)
        bg_templates2 = _create_rand_sub_bg_templates(data_bg_targets)
        bg_weights2 = build_mm_background_weights(mm_stage2_input, background_fit2[0])

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
        )
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
        )
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
        )
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
        )

        if ParticleType == "kaon" and scale_factor != 0.0:
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
            )
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
            )
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
            )
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
            )

        for key, hist in data_bg_targets.items():
            hist.Add(bg_templates2[key], -1)

        H_MM_fit2sub_DATA.Add(background_fit2[1], -1)
        H_MM_DATA.Add(background_fit2[0], -1)
        H_MM_full_DATA.Add(background_fit2[1], -1)
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
    histDict["H_epsilon_DATA"] =     H_epsilon_DATA
    histDict["H_MM_DATA"] =     H_MM_DATA
    histDict["H_MM_rand_dummy_DATA"] =     H_MM_rand_dummy_DATA
    histDict["H_MM_dummy_DATA"] =     H_MM_dummy_DATA
    histDict["H_MM_full_DATA"] =     H_MM_full_DATA        
    histDict["H_MM_fit2sub_DATA"] =     H_MM_fit2sub_DATA
    histDict["H_MM_fit1sub_DATA"] =     H_MM_fit1sub_DATA
    histDict["H_MM_pisub_DATA"] =     H_MM_pisub_DATA
    histDict["H_MM_nosub_DATA"] =     H_MM_nosub_DATA
    histDict["H_th_DATA"] =     H_th_DATA
    histDict["H_ph_DATA"] =     H_ph_DATA
    histDict["H_ph_q_DATA"] =     H_ph_q_DATA
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

    ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+'(')

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

    CMMfit2sub.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))    

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

    CMMfit1sub.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

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

    CMMpisub.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

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

    CMMsub.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # MM dummy plots    
    CMMdummy = TCanvas()

    histDict["H_MM_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_dummy_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_dummy_DATA"].SetFillColor(kBlack)  # Set fill color to black
    histDict["H_MM_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_dummy_DATA"].Draw("hist same")

    CMMdummy.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # MM rand dummy plots    
    CMMranddummy = TCanvas()

    histDict["H_MM_rand_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_rand_dummy_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_rand_dummy_DATA"].SetFillColor(kBlack)  # Set fill color to black    
    histDict["H_MM_rand_dummy_DATA"].SetLineColor(1)
    histDict["H_MM_rand_dummy_DATA"].Draw("hist same")

    CMMranddummy.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    
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
        c_pid.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    else:
        c_pid.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')

    if ParticleType == "kaon":
        
        ##
        # HGCer Hole Plots
        c_hgcervsMM = TCanvas()

        c_hgcervsMM.Divide(2,2)

        c_hgcervsMM.cd(1)
        P_hgcer_xAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(2)
        P_hgcer_nohole_xAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(3)
        P_hgcer_yAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_yAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(4)
        P_hgcer_nohole_yAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Draw("colz")    

        c_hgcervsMM.Draw()

        c_hgcervsMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

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

        c_hgcer_hole.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')      

    _print_rand_timer("rand_sub plotting {}".format(phi_setting), perf_counter() - stage_start)
    _print_rand_timer("rand_sub total {}".format(phi_setting), perf_counter() - total_start)
    return histDict
