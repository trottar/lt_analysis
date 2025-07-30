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

    # ---- build a copy that keeps only side‑band bins ----
    sig_lo, sig_hi = 1.107, 1.123
    sb_left  = (1.070, 1.095)
    sb_right = (1.135, 1.180)

    h_sb = hist.Clone(hist.GetName() + "_sb")
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)
        if sig_lo <= x <= sig_hi:
            h_sb.SetBinContent(ib, 0)
            h_sb.SetBinError  (ib, 0)

    # ---- Chebyshev (pol2) fit to the side‑bands only ----
    fit_func = TF1("fit_func", "pol2", mm_min, mm_max)
    h_sb.Fit(fit_func, "Q0")         # quiet, no draw

    # ---- integrate background under the signal window ----
    bin_w   = h_sb.GetBinWidth(1)
    bg_par  = fit_func.Integral(sig_lo, sig_hi) / bin_w

    # simple error estimate: propagate the constant‑term uncertainty
    bg_err  = abs(fit_func.GetParError(0)) * (sig_hi - sig_lo) / bin_w

    # --------------------------------------------------------------
    # Scale the fitted background to the stats that remain after
    # random, dummy, pion subtraction and MM cuts (set in rand_sub.py)
    # --------------------------------------------------------------
    scale = inpDict["bg_stat_scale"]
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        bg_par *= scale
        bg_err *= scale

    # ---- done ----
    return fit_func, bg_par
