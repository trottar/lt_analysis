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
                             hist_axis_ref,         # the cut histogram
                             mm_min, mm_max,
                             n_pad=0,
                             allow_negative=False,
                             ref_binwidth=None):    # NEW -----------------

    ax  = hist_axis_ref.GetXaxis()
    nb  = hist_axis_ref.GetNbinsX()
    if ref_binwidth is None:
        ref_binwidth = ax.GetBinWidth(1)            # fallback

    h_bg = hist_axis_ref.Clone(hist_axis_ref.GetName() + "_bg_fit_pad")
    h_bg.Reset("ICES")

    i_lo = max(1, ax.FindBin(mm_min) - n_pad)
    i_hi = min(nb, ax.FindBin(mm_max) + n_pad)

    for ib in range(1, nb + 1):
        if i_lo <= ib <= i_hi:
            x_lo = ax.GetBinLowEdge(ib)
            x_hi = x_lo + ax.GetBinWidth(ib)
            y    = fit_func.Integral(x_lo, x_hi) / ref_binwidth   # exact counts
            if not allow_negative and y < 0.0:
                y = 0.0
            h_bg.SetBinContent(ib, y)
            h_bg.SetBinError  (ib, 0.0)
        else:
            h_bg.SetBinContent(ib, 0.0)
            h_bg.SetBinError  (ib, 0.0)

    return h_bg

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(phi_setting,
           inpDict,
           hist_full,            #  H_MM_pisub_DATA  (wide axis, no MM cut)
           hist_cut=None):       #  H_MM_DATA        (narrow axis, MM cut)

    if hist_cut is None:                         # keep old calls valid
        hist_cut = hist_full

    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    # ---------------- side-band mask & wide fit (unchanged) ----------------
    sb_left  = inpDict.get("sb_left",  (1.070, 1.095))
    sb_right = inpDict.get("sb_right", (1.135, 1.180))

    h_sb = hist_full.Clone(hist_full.GetName() + "_sb")
    ...                                    # (masking code stays exactly as is)

    fit_func = ROOT.TF1("fit_func", "pol1", sb_left[0], sb_right[1])
    h_sb.Fit(fit_func, "Q0")               # wide fit on **hist_full**

    # ---------------- build background on the cut axis ---------------------
    ref_bw = hist_full.GetXaxis().GetBinWidth(1)   # ← original bin-width

    fit_hist_inrange = get_fit_histogram_padded(
                           fit_func,
                           hist_cut,               # axis to copy  (H_MM_DATA)
                           mm_min, mm_max,
                           n_pad=0,
                           allow_negative=False,
                           ref_binwidth=ref_bw)    # ← pass the wide width

    # return tuple shape unchanged so your “[0]” indexing still works
    fit_vis = fit_func.Clone(f"{hist_cut.GetName()}_bg_vis")
    bg_par  = fit_func.Integral(mm_min, mm_max) / ref_bw
    return fit_hist_inrange, fit_vis, bg_par
