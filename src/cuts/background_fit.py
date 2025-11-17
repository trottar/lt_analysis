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
    "Q3p0W2p32" : {

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {                
                #"left":  (1.08, 1.10), #Q4p4W2p74
                #"right": (1.20, 1.25), #Q4p4W2p74
                "left":  (1.06, 1.10),
                "right": (1.18, 1.21),                
            }
        },   

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.05, 1.10), #Q4p4W2p74
                #"right": (1.20, 1.22), #Q4p4W2p74              
                "left":  (1.06, 1.11),
                "right": (1.14, 1.20),
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.08, 1.10), #Q4p4W2p74
                #"right": (1.20, 1.25), #Q4p4W2p74                 
                "left":  (1.08, 1.10),
                "right": (1.18, 1.23),
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.05, 1.10), #Q4p4W2p74
                #"right": (1.20, 1.22), #Q4p4W2p74
                "left":  (1.08, 1.10),
                "right": (1.18, 1.20),                
            }
        },      

            # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Right_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.08, 1.10), #Q4p4W2p74
                #"right": (1.20, 1.25), #Q4p4W2p74 
                "left":  (1.06, 1.10),
                "right": (1.18, 1.21),                 
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
    },    
    "Q4p4W2p74" : {

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
                "right": (1.20, 1.25), 
            }
        },   

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.05, 1.10), 
                "right": (1.20, 1.22),
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
                "right": (1.20, 1.25), 
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.05, 1.10),
                "right": (1.20, 1.22),
            }
        },      

            # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Right_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
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
    },
    "testing" : {
        # --- ** GOOD ** linear background -----------
        "linear": {
            "func_expr": "pol1",          # Linear
            "n_par":      2,
            "sidebands":  {
                "left":  (1.00, 1.06),
                "right": (1.20, 1.25), 
                },
        },

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.05, 1.10),
                "right": (1.20, 1.22),
            }
        },  

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
                "right": (1.20, 1.25), 
            }
        },   

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_lowe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.05, 1.10), 
                "right": (1.20, 1.22),
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Center_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
                "right": (1.20, 1.25), 
            }
        },     

        # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Left_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.05, 1.10),
                "right": (1.20, 1.22),
            }
        },      

            # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
        "cheb2_Right_highe": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                "left":  (1.08, 1.10),
                "right": (1.20, 1.25), 
            }
        },   

        # --- ** GOOD ** 3rd-order polynomial ---------------------------
        "pol3": {
            "func_expr": "pol3",    # a0 + a1·x + a2·x² + a3·x³
            "n_par":      4,
            "sidebands": {
                "left":  (1.05, 1.10),
                "right": (1.20, 1.25),
            }
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
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

def is_good_background_shape(fit_func, x_min, x_max, *, n_grid=400):
    """
    Return True if the fitted background is acceptable on [x_min, x_max]:

      - Non-negative:        f(x) >= 0        for all x in [x_min, x_max]
      - Concave (not convex): d²f/dx² <= 0   for all x in [x_min, x_max]

    The second derivative is estimated numerically using finite differences
    on a uniform grid.
    """
    # Sample the fit on a dense grid
    xs = np.linspace(x_min, x_max, n_grid)
    ys = np.array([fit_func.Eval(x) for x in xs], dtype="float64")

    # Reject if the function evaluation itself is pathological
    if not np.all(np.isfinite(ys)):
        return False

    '''
    # Reject if it ever goes negative
    if np.any(ys < 0.0):
        return False
    '''
    # Numerical first and second derivatives
    dy = np.gradient(ys, xs)
    d2y = np.gradient(dy, xs)

    if not np.all(np.isfinite(d2y)):
        return False

    # Reject if the second derivative is positive anywhere (convex region)
    if np.any(d2y > 0.0):
        return False

    return True

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
        no_bg_subtract=False  # no background subtraction
):
    """
    Generic side-band fit and background subtraction.

    Parameters
    ----------
    model_key : str
        Which entry in BG_MODELS to use.
        Every model supplies the TF1 expression and its default side-bands.
    All other parameters are unchanged compared to the legacy version.

    Returns
    -------
    fit_hist_inrange : TH1
        Histogram of the fitted background, binned like `hist_mm_cut`, but
        only non-zero inside [mm_min, mm_max].
    fit_vis : TF1
        The fitted function, useful for plotting.
    bg_par : float
        Legacy quantity: background counts per bin in the signal window.
    f_sig : float
        NEW: signal fraction = max(0, 1 - N_bg / N_tot) in the MM signal window.
    """
    # ------------------------------- setup --------------------------------
    Q2 = inpDict["Q2"]
    W  = inpDict["W"]
    model  = BG_MODELS[f"Q{Q2}W{W}"][model_key]  # raises KeyError if wrong
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

    # Special case: explicitly no background subtraction
    if no_bg_subtract:
        fit_func = TF1("fit_func_zero", "0", mm_min, mm_max)
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0.0
        f_sig    = 1.0  # keep everything as "signal"
        return fit_func, fit_vis, bg_par, f_sig

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

    # ------------------- shape sanity check (new) -------------------------
    # Reject fits that:
    #   - are convex anywhere (d²f/dx² > 0) on [mm_min, mm_max], or
    #   - go negative anywhere on [mm_min, mm_max].
    #
    # If the fit fails this check, fall back to "no background subtraction"
    # for this histogram: zero background, f_sig = 1.
    if not is_good_background_shape(fit_func, mm_min, mm_max):
        # zero background function on the full MM range
        fit_func_zero = TF1("fit_func_zero_bad", "0", mm_min, mm_max)
        fit_vis = fit_func_zero.Clone(f"{hist.GetName()}_bg_vis")

        # build a zero-valued background histogram on the MM-cut binning
        fit_hist_inrange = get_fit_histogram_padded(
            fit_func_zero,
            hist_mm_cut,
            mm_min,
            mm_max,
            n_pad=0
        )

        bg_par = 0.0
        f_sig = 1.0  # treat everything as signal if the background is unphysical

        return fit_hist_inrange, fit_vis, bg_par, f_sig

    # ------------------- proceed with accepted fit ------------------------
    # integral of background under the signal window
    # (N_bg = expected background counts in [sig_lo, sig_hi])
    N_bg = max(0.0, fit_func.Integral(sig_lo, sig_hi))

    # legacy quantity: per-bin background in the window
    bw = hist.GetBinWidth(1)
    bg_par = 0.0
    if bw > 0.0:
        bg_par = N_bg / bw

    # global scaling of stat errors if user supplied it
    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        N_bg  *= scale
        bg_par *= scale

    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    # ---------------------- signal fraction in MM -------------------------
    ax_mm    = hist_mm_cut.GetXaxis()
    ib_sig_lo = ax_mm.FindBin(sig_lo)
    ib_sig_hi = ax_mm.FindBin(sig_hi)

    # total data counts in the MM signal window
    N_tot = hist_mm_cut.Integral(ib_sig_lo, ib_sig_hi)

    # default: if N_tot <= 0, just keep everything
    f_sig = 1.0
    if N_tot > 0.0:
        f_sig = max(0.0, min(1.0, 1.0 - N_bg / N_tot))

    # ------------------------ build background hist -----------------------
    fit_hist_inrange = get_fit_histogram_padded(
        fit_func,
        hist_mm_cut,
        mm_min, mm_max,
        n_pad=0,
        ref_binwidth=bw
    )

    # NOTE: now returning 4 values instead of 3
    return fit_hist_inrange, fit_vis, bg_par, f_sig
