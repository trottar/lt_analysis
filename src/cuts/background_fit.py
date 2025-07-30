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

def bg_fit(phi_setting, inpDict, hist):
    
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

    # ------------------------------------------------------------------
    #  Signal window  : 1.107 – 1.123 GeV/c²
    #  Rad‑tail guard : 22 MeV to the right of the peak (skip real Λ)
    #  Side‑bands     : 1.070 – 1.095  and  1.145 – 1.180 GeV/c²
    # ------------------------------------------------------------------
    sig_lo, sig_hi  = 1.107, 1.123
    tail_gap        = 0.022                       # 22 MeV guard band
    sb_left         = (1.070, 1.095)
    sb_right        = (sig_hi + tail_gap, 1.180)  # starts at 1.145

    # ---- build a copy that keeps only the side‑bands ----
    h_sb = hist.Clone(hist.GetName() + "_sb")
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)

        # kill: signal region  +  radiative tail guard
        if (sig_lo <= x <= sig_hi + tail_gap):
            h_sb.SetBinContent(ib, 0)
            h_sb.SetBinError  (ib, 0)

        # also kill anything left of sb_left or right of sb_right
        if not (sb_left[0] <= x <= sb_left[1] or
                sb_right[0] <= x <= sb_right[1]):
            h_sb.SetBinContent(ib, 0)
            h_sb.SetBinError  (ib, 0)

    # ---- first‑order Chebyshev (pol1) is now sufficient ----
    fit_func = TF1("fit_func", "pol1", mm_min, mm_max)
    h_sb.Fit(fit_func, "Q0")                # silent fit

    # ---- integrate background under the Λ window ----
    bin_w     = h_sb.GetBinWidth(1)
    bg_par    = fit_func.Integral(sig_lo, sig_hi) / bin_w
    bg_err    = abs(fit_func.GetParError(0)) * (sig_hi - sig_lo) / bin_w

    return fit_func, bg_par
