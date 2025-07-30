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
    # Dynamic background fit that never looks outside [mm_min, mm_max]
    # ------------------------------------------------------------------
    #
    #  User‑supplied limits (already defined in your script)
    #       mm_min = 1.100
    #       mm_max = 1.140
    #
    #  We carve three regions inside that window:
    #    1) Guard   : ±Δ around the Λ peak (excluded from the fit)
    #    2) Side L  : [mm_min,  sig_lo - gap)
    #    3) Side R  : (sig_hi + gap,  mm_max]
    #
    #  Only the sidebands are fitted with a linear (pol1) Chebyshev.
    #  If one side‑band is too narrow (< 3 bins) it is dropped automatically.
    # ------------------------------------------------------------------
    sig_lo, sig_hi = 1.107, 1.123          # Λ window (GeV/c²)
    gap            = 0.004                 # 4 MeV guard band
    sb_left        = (mm_min,  sig_lo - gap)
    sb_right       = (sig_hi + gap, mm_max)

    h_sb = hist.Clone(hist.GetName() + "_sb")
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)

        # Kill central region (= signal + guard gap)
        if sig_lo - gap <= x <= sig_hi + gap:
            h_sb.SetBinContent(ib, 0)
            h_sb.SetBinError  (ib, 0)
            continue

        # Kill anything *outside* the allowed side‑bands
        if not (sb_left[0] <= x <= sb_left[1] or
                sb_right[0] <= x <= sb_right[1]):
            h_sb.SetBinContent(ib, 0)
            h_sb.SetBinError  (ib, 0)

    # --- decide order dynamically: use pol0 if statistics are tiny ---
    nbins_left  = h_sb.GetXaxis().FindBin(sb_left [1]) - \
                  h_sb.GetXaxis().FindBin(sb_left [0]) + 1
    nbins_right = h_sb.GetXaxis().FindBin(sb_right[1]) - \
                  h_sb.GetXaxis().FindBin(sb_right[0]) + 1
    n_sidebins  = max(0, nbins_left) + max(0, nbins_right)

    poly_order  = 1 if n_sidebins >= 6 else 0        # need ≥6 bins for pol1
    fit_func    = TF1("fit_func", f"pol{poly_order}", mm_min, mm_max)

    h_sb.Fit(fit_func, "Q0")                         # silent fit

    # --- integrate background under the Λ window *inside* the same limits ---
    bin_w    = h_sb.GetBinWidth(1)
    bg_par   = fit_func.Integral(sig_lo, sig_hi) / bin_w
    # propagate only constant‑term error for a conservative estimate
    bg_err   = abs(fit_func.GetParError(0)) * (sig_hi - sig_lo) / bin_w

    return fit_func, bg_par

