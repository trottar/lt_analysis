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

def get_fit_histogram_padded(fit_func, hist, mm_min, mm_max, n_pad=0):
    """
    Clone *hist*, fill it with the fitted function, then keep the
    contents only in the window [mm_min, mm_max] extended by *n_pad*
    extra **bins** on each side.  Everything outside that window is
    forced to zero (contents *and* errors).

    Parameters
    ----------
    fit_func : ROOT.TF1
        The background-fit function you want to visualise / subtract.
    hist : ROOT.TH1
        The template histogram that defines the binning.
    mm_min, mm_max : float
        Nominal physics window (in the same units as the histogram X-axis).
    n_pad : int, optional
        How many whole bins to pad on each side.  Default is **0**
        so the function is strictly limited to [mm_min, mm_max] unless
        you explicitly ask for padding.
    """
    # --- determine bin numbers instead of working in x so that
    #     we cannot straddle a bin edge by accident ------------------------
    ax   = hist.GetXaxis()
    i_lo = ax.FindBin(mm_min)            # first bin INSIDE the signal region
    i_hi = ax.FindBin(mm_max)            # last  bin INSIDE the signal region

    # expand by the requested padding (but keep within histogram bounds)
    i_lo = max(1,                 i_lo - n_pad)
    i_hi = min(hist.GetNbinsX(),  i_hi + n_pad)

    # --- build the padded-fit histogram ----------------------------------
    h_fit = hist.Clone(hist.GetName() + "_bg_fit_pad")
    h_fit.Reset()

    for ib in range(1, h_fit.GetNbinsX() + 1):
        if i_lo <= ib <= i_hi:
            x = h_fit.GetBinCenter(ib)
            y = fit_func.Eval(x)
            h_fit.SetBinContent(ib, y)
            h_fit.SetBinError  (ib, 0.0)      # or propagate your own errors
        else:
            h_fit.SetBinContent(ib, 0.0)
            h_fit.SetBinError  (ib, 0.0)

    return h_fit

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(phi_setting, inpDict, hist):

    if no_bg_subtract:
        fit_func = TF1("fit_func_zero", "0", inpDict["mm_min"], inpDict["mm_max"])
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        return fit_func, fit_vis, bg_par    
    
    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    # --- Use physics-motivated wide sidebands for fitting, NOT just mm_min/mm_max ---
    sb_left = inpDict.get("sb_left", (1.070, 1.095))
    sb_right = inpDict.get("sb_right", (1.135, 1.180))
    sig_lo = max(1.107, mm_min)
    sig_hi = min(1.123, mm_max)

    print(f"Signal window: [{sig_lo:.3f}, {sig_hi:.3f}]")
    print(f"Sideband left: [{sb_left[0]:.3f}, {sb_left[1]:.3f}]")
    print(f"Sideband right: [{sb_right[0]:.3f}, {sb_right[1]:.3f}]")

    h_sb = hist.Clone(hist.GetName() + "_sb")

    def is_valid_sideband(sb):
        return sb[1] > sb[0]

    n_sb_bins = 0
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)
        in_sb_left  = is_valid_sideband(sb_left)  and (sb_left[0]  <= x <= sb_left[1])
        in_sb_right = is_valid_sideband(sb_right) and (sb_right[0] <= x <= sb_right[1])
        if not (in_sb_left or in_sb_right):
            h_sb.SetBinContent(ib, 0.0)
            h_sb.SetBinError(ib, 0.0)
        else:
            n_sb_bins += 1

    print("Nonzero bins in h_sb after masking:", n_sb_bins)

    if n_sb_bins < 2:
        print("[WARNING] All sideband bins are empty after masking! Fit will be flat zero.")
        fit_func = TF1("fit_func_zero", "0", mm_min, mm_max)
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        fit_hist_inrange = get_fit_histogram_in_range(fit_func, hist, mm_min, mm_max)
        return fit_hist_inrange, fit_vis, bg_par

    fit_min = sb_left[0]
    fit_max = sb_right[1]
    fit_func = TF1("fit_func", "pol1", fit_min, fit_max)
    h_sb.Fit(fit_func, "Q0")

    bin_w = h_sb.GetBinWidth(1)
    bg_par = fit_func.Integral(sig_lo, sig_hi) / bin_w
    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        bg_par *= scale

    fit_hist_inrange = get_fit_histogram_padded(fit_func, hist, mm_min, mm_max, n_pad=0)
    return fit_hist_inrange, fit_vis, bg_par
