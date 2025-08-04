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

def get_fit_histogram_padded(fit_func,
                             hist,
                             mm_min, mm_max,
                             n_pad: int = 0,
                             allow_negative: bool = False,
                             ref_binwidth: float | None = None):
    """
    Return a histogram that contains the background function *fit_func*
    inside [mm_min , mm_max] (± n_pad *bins*).  Outside that window the
    contents – and errors – are forced to zero.

    Parameters
    ----------
    fit_func      : ROOT.TF1      – function already fitted on the **wide** MM spectrum
    hist          : ROOT.TH1      – histogram whose axis you want to copy
    mm_min, mm_max: float         – MM window limits
    n_pad         : int           – extra *bins* to keep on each side (default 0)
    allow_negative: bool          – if False, clamp negative predictions to 0
    ref_binwidth  : float | None  – bin-width used to normalise *fit_func*.
                                    • Pass the wide-spectrum bin-width (recommended).
                                    • If None, fall back to hist’s own bin-width.
    """
    ax  = hist.GetXaxis()
    nb  = hist.GetNbinsX()

    # use the caller-supplied reference width or fall back to this axis
    if ref_binwidth is None:
        ref_binwidth = ax.GetBinWidth(1)

    # translate physics window → bin indices, include optional padding
    i_lo = max(1, ax.FindBin(mm_min) - n_pad)
    i_hi = min(nb, ax.FindBin(mm_max) + n_pad)

    # fresh empty clone (keeps binning, titles, Sumw2 array)
    h_bg = hist.Clone(hist.GetName() + "_bg_fit_pad")
    h_bg.Reset("ICES")                         # clean slate, keep Sumw2

    # fill background only where requested
    for ib in range(1, nb + 1):
        if i_lo <= ib <= i_hi:
            # exact integral of the density across the bin, then convert to counts
            x_lo = ax.GetBinLowEdge(ib)
            x_hi = x_lo + ax.GetBinWidth(ib)
            y_pred = fit_func.Integral(x_lo, x_hi) / ref_binwidth

            if not allow_negative and y_pred < 0.0:
                y_pred = 0.0

            h_bg.SetBinContent(ib, y_pred)
            h_bg.SetBinError  (ib, 0.0)       # param-error propagation goes here if desired
        else:
            h_bg.SetBinContent(ib, 0.0)
            h_bg.SetBinError  (ib, 0.0)

    return h_bg

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(phi_setting, inpDict, hist, hist_mm_cut=None):

    if hist_mm_cut is None:
         hist_mm_cut = hist

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

    bg_par = fit_func.Integral(sig_lo, sig_hi) / hist.GetBinWidth(1)
    bg_par = max(0.0, bg_par)  # physical prior
    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        bg_par *= scale

    fit_hist_inrange = get_fit_histogram_padded(
                       fit_func,
                       hist_mm_cut,                         # use the cut axis
                       mm_min, mm_max,
                       n_pad=0,
                       ref_binwidth=hist.GetXaxis().GetBinWidth(1))  # wide Δx

    return fit_hist_inrange, fit_vis, bg_par