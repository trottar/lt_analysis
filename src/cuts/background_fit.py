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

# -------------------------------------------------------------------------
# Catalogue of background-shape models and their default side-bands
# -------------------------------------------------------------------------
BG_MODELS = {
    # --- ** GOOD ** linear background -----------
    "linear": {
        "func_expr": "pol1",          # Linear
        "n_par":      2,
        "sidebands":  {
            #"left":  (1.070, 1.095),
            #"right": (1.135, 1.180)
            "left":  (1.00, 1.06),
            "right": (1.20, 1.25), 
            },
    },

    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
    "cheb2": {
        "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
        "n_par":      3,
        "sidebands": {
            #"left":  (1.00, 1.06),
            #"right": (1.14, 1.20),
            "left":  (1.00, 1.06),
            "right": (1.20, 1.25), 
        }
    },

    # --- ** GOOD ** 3rd-order polynomial ---------------------------
    "pol3": {
        "func_expr": "pol3",    # a0 + a1·x + a2·x² + a3·x³
        "n_par":      4,
        "sidebands": {
            #"left":  (1.00, 1.06),
            #"right": (1.20, 1.22),
            "left":  (1.00, 1.06),
            "right": (1.20, 1.25),
        }
    },

    # --- Sigma peak 2nd-order Chebyshev ---------------------------
    "sigma_peak": {
        "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
        "n_par":      3,
        "sidebands": {
            #"left":  (1.00, 1.06),
            #"right": (1.14, 1.20),
            "left":  (1.165, 1.20),
            "right": (1.20, 1.21), 
        }
    },    

    # Next functions all look similar to Chebyshev order 2
    # --- Crystal-Ball (Gaussian + tail) *** NEEED TO FIX SIDEBANDS *** -----------------
    "crystalball": {
        "func_expr": "crystalball",  # TMath::CrystalBall(x; m, σ, α, n)
        "n_par":      5,
        "sidebands": {
            "left":  (1.00, 1.06),
            "right": (1.20, 1.22),
        }
    },

    # --- Landau distribution ----------------------------
    "landau": {
        "func_expr": "landau",      # p[0]·Landau(x; mpv, σ)
        "n_par": 3,
        "sidebands": {
            "left":  (1.00, 1.06),
            "right": (1.20, 1.22),
        }
    },

    # --- Exponential × quadratic polynomial -------------
    "exppol2": {
        "func_expr": "[0]*exp([1]*x)*(1+[2]*x+[3]*x*x)",  # norm·e^(slope·x)·(1 + c1·x + c2·x²)
        "n_par":      4,
        "sidebands": {
            "left":  (1.00, 1.06),
            "right": (1.20, 1.22),
        }
    },

}

################################################################################################################################################

def get_fit_histogram_padded(fit_func,
                             hist,
                             mm_min, mm_max,
                             n_pad=0,
                             allow_negative=False,
                             ref_binwidth=None):
    """
    Build a background histogram from *fit_func*.
    • Only bins inside [mm_min, mm_max] (± n_pad bins) are filled.
    • Everything outside that window is set to 0.
    • *ref_binwidth* must be the bin-width of the wide/full MM histogram
      that was used to fit *fit_func*.
      If you pass None, the routine falls back to hist’s own bin-width.
    """
    ax  = hist.GetXaxis()
    nb  = hist.GetNbinsX()

    # if caller did not supply a reference width, fall back to this axis
    if ref_binwidth is None:
        ref_binwidth = ax.GetBinWidth(1)

    # translate physics window → bin indices, add optional padding
    i_lo = max(1,  ax.FindBin(mm_min) - n_pad)
    i_hi = min(nb, ax.FindBin(mm_max) + n_pad)

    # make a clean clone for the background
    h_bg = hist.Clone(hist.GetName() + "_bg_fit_pad")
    h_bg.Reset("ICES")                     # wipe contents, keep Sumw2

    # fill background
    for ib in range(1, nb + 1):
        if i_lo <= ib <= i_hi:
            x_lo = ax.GetBinLowEdge(ib)
            x_hi = x_lo + ax.GetBinWidth(ib)
            # integral of density across the bin → counts
            y_pred = fit_func.Integral(x_lo, x_hi) / ref_binwidth
            if (not allow_negative) and (y_pred < 0.0):
                y_pred = 0.0
            h_bg.SetBinContent(ib, y_pred)
            h_bg.SetBinError  (ib, 0.0)
        else:
            h_bg.SetBinContent(ib, 0.0)
            h_bg.SetBinError  (ib, 0.0)

    return h_bg

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(
        phi_setting,
        inpDict,
        hist,
        hist_mm_cut=None,
        *,
        model_key="linear",   # ← just pick a key from BG_MODELS
        no_bg_subtract=False         # no background subtraction
):
    """
    Generic side-band fit and background subtraction.

    Parameters
    ----------
    model_key : str
        Which entry in BG_MODELS to use.
        Every model supplies the TF1 expression and its default side-bands.
    All other parameters are unchanged compared to the legacy version.
    """
    # ------------------------------- setup --------------------------------
    model = BG_MODELS[model_key]                 # raises KeyError if wrong
    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    # allow user to override SBs ad-hoc via inpDict; otherwise fall back to
    # the defaults stored in the model definition
    sb_left  = inpDict.get("sb_left",  model["sidebands"]["left"])
    sb_right = inpDict.get("sb_right", model["sidebands"]["right"])

    # signal window: keep the physics defaults unless over-ridden
    sig_lo = max(inpDict.get("sig_lo", 1.107), mm_min)
    sig_hi = min(inpDict.get("sig_hi", 1.123), mm_max)

    # ---------------------------------------------------------------------
    if hist_mm_cut is None:
        hist_mm_cut = hist

    if no_bg_subtract:
        fit_func = TF1("fit_func_zero", "0", mm_min, mm_max)
        return fit_func, fit_func.Clone(f"{hist.GetName()}_bg_vis"), 0.0

    # -------------------------- mask to side-bands ------------------------
    h_sb = hist.Clone(hist.GetName() + "_sb")
    h_sb.ResetStats()  # avoids ROOT complaints about negative stats

    def in_any_sb(x):
        return (sb_left[0]  <= x <= sb_left[1]) or \
               (sb_right[0] <= x <= sb_right[1])

    for ib in range(1, h_sb.GetNbinsX() + 1):
        if not in_any_sb(h_sb.GetBinCenter(ib)):
            h_sb.SetBinContent(ib, 0.0)
            h_sb.SetBinError(ib, 0.0)

    # ------------------------------ fit -----------------------------------
    fit_min, fit_max = sb_left[0], sb_right[1]
    fit_func = TF1("fit_func", model["func_expr"], fit_min, fit_max)
    h_sb.Fit(fit_func, "Q0")  # quiet, no UI

    # integral of background under the signal window, per-bin normalised
    bg_par = max(0.0, fit_func.Integral(sig_lo, sig_hi) / hist.GetBinWidth(1))

    # global scaling of stat errors if user supplied it
    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        bg_par *= scale

    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    # produce a histogram of the fit on the same binning as (hist_mm_cut)
    fit_hist_inrange = get_fit_histogram_padded(
        fit_func,
        hist_mm_cut,
        mm_min, mm_max,
        n_pad=0,
        ref_binwidth=hist.GetXaxis().GetBinWidth(1)
    )

    return fit_hist_inrange, fit_vis, bg_par
