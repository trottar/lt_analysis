#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-18 02:57:48 trottar"
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
from ROOT import TF1
import sys, os, math

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

################################################################################################################################################

def get_fit_histogram_in_range(fit_func, hist, mm_min, mm_max):
    """
    Create a histogram from fit_func matching 'hist', but zero outside [mm_min, mm_max].
    """
    h_fit = hist.Clone(hist.GetName() + "_fit_inrange")
    h_fit.Reset()  # Clear bin contents and errors

    for ibin in range(1, h_fit.GetNbinsX() + 1):
        x = h_fit.GetBinCenter(ibin)
        if mm_min <= x <= mm_max:
            h_fit.SetBinContent(ibin, fit_func.Eval(x))
            # Set error if needed; left as 0 by default
        else:
            h_fit.SetBinContent(ibin, 0.0)
            h_fit.SetBinError(ibin, 0.0)
    return h_fit

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(phi_setting, inpDict, hist):

    if no_bg_subtract:
        # Create a zero function over the fit range
        fit_func = TF1("fit_func_zero", "0", inpDict["mm_min"], inpDict["mm_max"])
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        # Return zero fit, visual fit, and zero background
        return fit_func, fit_vis, bg_par    
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]

    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    # ------------------------------------------------------------------
    # Dynamic side‑band background estimate (Chebyshev, order‑2)
    # ------------------------------------------------------------------
    #
    #  Signal window (Λ peak)   : 1.107 – 1.123 GeV/c²
    #  Side‑bands               : 1.070 – 1.095 GeV/c²  and  1.135 – 1.180 GeV/c²
    #
    #  The function returns:
    #     fit_func  – the TF1 that describes the background shape
    #     bg_par    – the **integrated** background counts under the Λ window
    #

    num_evts = hist.GetEntries()

    sig_lo, sig_hi  = mm_min, mm_max
    tail_gap        = 0.022

    sb_left         = (mm_min,  mm_max)
    sb_right_start  = sig_hi + tail_gap
    sb_right_end    = mm_max
    sb_right        = (max(sb_right_start, mm_min), sb_right_end)

    # --- Diagnostics ---
    print("Sideband left:", sb_left)
    print("Sideband right:", sb_right)

    # ---- build a copy that keeps only the side‑bands ----
    h_sb = hist.Clone(hist.GetName() + "_sb")

    # Helper to check if sideband has width
    def is_valid_sideband(sb):
        return sb[1] > sb[0]

    n_sb_bins = 0
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)
        in_sb_left = is_valid_sideband(sb_left) and (sb_left[0] <= x <= sb_left[1])
        in_sb_right = is_valid_sideband(sb_right) and (sb_right[0] <= x <= sb_right[1])
        if not (in_sb_left or in_sb_right):
            h_sb.SetBinContent(ib, 0.0)
            h_sb.SetBinError(ib, 0.0)
        else:
            n_sb_bins += 1

    print("Nonzero bins in h_sb after masking:", n_sb_bins)

    # ---- Handle empty sidebands robustly ----
    if n_sb_bins < 2:
        print("[WARNING] All sideband bins are empty after masking! Fit will be flat zero.")
        fit_func = TF1("fit_func_zero", "0", mm_min, mm_max)
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        fit_hist_inrange = get_fit_histogram_in_range(fit_func, hist, mm_min, mm_max)
        return fit_hist_inrange, fit_vis, bg_par

    # ---- Fit sidebands ----
    fit_func = TF1("fit_func", "pol1", mm_min, mm_max)
    h_sb.Fit(fit_func, "Q0")                # silent fit

    # ---- integrate background under the Λ window ----
    bin_w     = h_sb.GetBinWidth(1)
    bg_par    = fit_func.Integral(sig_lo, sig_hi) / bin_w
    bg_err    = abs(fit_func.GetParError(0)) * (sig_hi - sig_lo) / bin_w

    # --------------------------------------------------------------
    # Keep the ORIGINAL fit to draw on H_MM_nosub_DATA
    # --------------------------------------------------------------
    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    # --------------------------------------------------------------
    # Create a STAT‑MATCHED copy to subtract from H_MM_DATA
    # The calling code puts the scale in  inpDict["bg_stat_scale"].
    # --------------------------------------------------------------
    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip,
                                  fit_func.GetParameter(ip) * scale)
        bg_par *= scale
        bg_err *= scale

    fit_hist_inrange = get_fit_histogram_in_range(fit_func, hist, mm_min, mm_max)
    # Return (scaled_function, visual_function, background_counts)
    return fit_hist_inrange, fit_vis, bg_par
