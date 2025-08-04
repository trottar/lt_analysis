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

####################################################################
# added helper -----------------------------------------------------
def project_fit_onto(hist_src, hist_dst, fit_func,
                     allow_negative=False, n_pad=0):
    """
    Return a clone of *hist_dst* whose bin contents are filled with the
    values of *fit_func* **in the dst range only**.
    """
    h = hist_dst.Clone(hist_dst.GetName() + "_bg_from_full")
    h.Reset()

    ax_src = hist_src.GetXaxis()
    ax_dst = hist_dst.GetXaxis()
    nb_dst = ax_dst.GetNbins()

    # Source bin width is what TH1::Fit used
    bw_src = ax_src.GetBinWidth(1)

    for ib in range(1, nb_dst + 1):
        x = ax_dst.GetBinCenter(ib)
        y = fit_func.Eval(x)
        if not allow_negative and y < 0.0:
            y = 0.0
        h.SetBinContent(ib, y)
        h.SetBinError  (ib, 0.0)

    return h
####################################################################

def bg_fit(phi_setting, inpDict, full_hist):
    """
    full_hist : TH1 that still contains *both* sidebands
                (e.g. 0.70–1.50 GeV for H(e,e′K⁺) on hydrogen)

    Returns
    --------
    h_bg_cut  : TH1 – background shape but *only* inside
                        [mm_min, mm_max] (zeros elsewhere)
    bg_int    : float – background counts in the signal window
    bg_err    : float – propagated uncertainty on that integral
    """
    mm_min, mm_max = inpDict["mm_min"], inpDict["mm_max"]

    # -----------------------------------------------------------------
    # 1) Fit polynomial to the two sidebands on the full histogram
    # -----------------------------------------------------------------
    sb_left  = inpDict.get("sb_left",  (1.070, 1.095))
    sb_right = inpDict.get("sb_right", (1.135, 1.180))

    bg_func = ROOT.TF1("bg", "pol1", sb_left[0], sb_right[1])

    # Fit left band, then add the right band (equal bin weights)
    full_hist.Fit(bg_func, "Q0R",  "", sb_left[0],  sb_left[1])
    full_hist.Fit(bg_func, "Q0R+", "", sb_right[0], sb_right[1])

    # -----------------------------------------------------------------
    # 2) Project that function onto the Λ window only
    #    (bins outside the window → 0)
    # -----------------------------------------------------------------
    h_bg_cut = get_fit_histogram_padded(bg_func,
                                        full_hist,
                                        mm_min, mm_max,
                                        n_pad=0, allow_negative=False)

    # -----------------------------------------------------------------
    # 3) Background integral and error in the same window
    # -----------------------------------------------------------------
    sig_lo = max(1.107, mm_min)
    sig_hi = min(1.123, mm_max)

    pars  = np.array([bg_func.GetParameter(i) for i in range(bg_func.GetNpar())])
    cov   = np.zeros((bg_func.GetNpar(), bg_func.GetNpar()))
    bg_func.GetCovarianceMatrix(cov.flatten().tolist())

    bw      = full_hist.GetBinWidth(1)
    bg_int  =  bg_func.Integral(     sig_lo, sig_hi) / bw
    bg_err  =  bg_func.IntegralError(sig_lo, sig_hi,
                                     pars, cov.flatten().tolist()) / bw

    return h_bg_cut, bg_int, bg_err
