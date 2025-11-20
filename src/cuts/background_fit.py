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

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.15) + [1]*(x-1.15)*(x-1.15)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.15, 1.18),
                        "right": (1.20, 1.30), 
                    }
                },   

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.15) + [1]*(x-1.15)*(x-1.15)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.15, 1.18),
                        "right": (1.20, 1.30), 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.15) + [1]*(x-1.15)*(x-1.15)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.15, 1.18),
                        "right": (1.20, 1.30), 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.15) + [1]*(x-1.15)*(x-1.15)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.15, 1.18),
                        "right": (1.20, 1.30), 
                    }
                },      

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.15) + [1]*(x-1.15)*(x-1.15)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.15, 1.18),
                        "right": (1.20, 1.30), 
                    }
                },   
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.05, 1.10),
                        "right": (1.15, 1.20),   
                    }
                },   

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74              
                        #"left":  (1.06, 1.11),
                        #"right": (1.14, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.15, 1.20),                 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74                 
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.23),
                        "left":  (1.05, 1.10),
                        "right": (1.15, 1.20),                 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.20),                
                        "left":  (1.05, 1.10),
                        "right": (1.15, 1.20),                 
                    }
                },      

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74 
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                 
                        "left":  (1.05, 1.10),
                        "right": (1.15, 1.20),                 
                    }
                },   
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

def is_good_background_shape(
        fit_func,
        x_min,
        x_max,
        *,
        neg_tol=0.25,
):
    """
    Sanity check on the *fit function itself*.

    The bin is accepted iff:
      - the function has a finite minimum and maximum on [x_min, x_max];
      - the maximum is > 0 (non-trivial background);
      - the minimum is >= -neg_tol.

    With neg_tol = 0.0 this is literally:
        "if it goes negative anywhere in y, it's a bad bin".
    """
    # ROOT does the extremum search
    f_min = float(fit_func.GetMinimum(x_min, x_max))
    f_max = float(fit_func.GetMaximum(x_min, x_max))

    # Non-finite or completely non-positive → reject
    if not (math.isfinite(f_min) and math.isfinite(f_max)):
        return False
    if f_max <= -neg_tol:
        return False

    # Reject if the minimum is below the allowed tolerance
    if f_min <= -neg_tol:
        return False

    return True

def move_edges_to_nonnegative(
        fit_func,
        sig_lo,
        sig_hi,
        *,
        neg_tol=0.01,
        max_iter=1000,
        step_frac=0.05,
        min_width=1e-4,
):
    """
    Adjust [sig_lo, sig_hi] so that f(sig_lo) and f(sig_hi) are >= -neg_tol.

    • Uses ONLY the current fit_func (no refitting).
    • On each iteration, any edge with f(edge) < -neg_tol is moved inward
      by step_frac * current_width.
    • Stops when both edges are >= -neg_tol, the window collapses, or
      max_iter is reached.

    Returns (new_lo, new_hi) on success, or (None, None) on failure.
    """
    orig_lo, orig_hi = float(sig_lo), float(sig_hi)

    for it in range(max_iter):
        width = sig_hi - sig_lo
        if width <= min_width:
            print(
                f"[move_edges] window collapsed at iter={it}: "
                f"orig=[{orig_lo:.5f},{orig_hi:.5f}], "
                f"final=[{sig_lo:.5f},{sig_hi:.5f}], width={width:.5g}"
            )
            return None, None

        f_lo = float(fit_func.Eval(sig_lo))
        f_hi = float(fit_func.Eval(sig_hi))

        print(
            f"[move_edges] iter={it} "
            f"window=[{sig_lo:.5f},{sig_hi:.5f}] width={width:.5g} "
            f"f_lo={f_lo:.5g} f_hi={f_hi:.5g}"
        )

        # If both edges are above threshold, we are done.
        if (f_lo >= -neg_tol) and (f_hi >= -neg_tol):
            print(
                f"[move_edges] success after {it} iters: "
                f"orig=[{orig_lo:.5f},{orig_hi:.5f}] -> "
                f"new=[{sig_lo:.5f},{sig_hi:.5f}]"
            )
            return sig_lo, sig_hi

        step = step_frac * width

        if f_lo < -neg_tol:
            new_lo = sig_lo + step
            print(
                f"  -> move lower edge: {sig_lo:.5f} -> {new_lo:.5f} "
                f"(f_lo={f_lo:.5g})"
            )
            sig_lo = new_lo

        if f_hi < -neg_tol:
            new_hi = sig_hi - step
            print(
                f"  -> move upper edge: {sig_hi:.5f} -> {new_hi:.5f} "
                f"(f_hi={f_hi:.5g})"
            )
            sig_hi = new_hi

    print(
        f"[move_edges] max_iter={max_iter} reached without both edges "
        f"above -neg_tol; final=[{sig_lo:.5f},{sig_hi:.5f}]"
    )
    return None, None

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
        fit_name=None
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
    if fit_name:
        model  = BG_MODELS[f"Q{Q2}W{W}"][fit_name][model_key] 
    else:
        model  = BG_MODELS[f"Q{Q2}W{W}"][model_key]

    mm_min = inpDict["mm_min"]  
    mm_max = inpDict["mm_max"]

    # allow user to override SBs ad-hoc via inpDict; otherwise fall back to
    # the defaults stored in the model definition
    sb_left  = inpDict.get("sb_left",  model["sidebands"]["left"])
    sb_right = inpDict.get("sb_right", model["sidebands"]["right"])

    # signal window: keep the physics defaults unless over-ridden
    sig_lo = max(inpDict.get("sig_lo", 1.05), mm_min)
    sig_hi = min(inpDict.get("sig_hi", 1.25), mm_max)

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

    # ------------------------------ fit once ------------------------------
    # Define the function over the FULL MM range, but only fit in sidebands.
    func_min, func_max = mm_min, mm_max
    fit_min,  fit_max  = sb_left[0], sb_right[1]

    fit_func = TF1("fit_func", model["func_expr"], func_min, func_max)
    h_sb.Fit(fit_func, "Q0", "", fit_min, fit_max)

    # ----------------------------------------------------------------------
    # Now enforce: the green background must not go negative in the
    # signal window. If it does, move negative edges inward until they
    # are above -neg_tol, then re-check the global shape once.
    # ----------------------------------------------------------------------
    neg_tol     = inpDict.get("bg_neg_tol", 0.01)
    max_refit   = inpDict.get("bg_max_refit", 1000)       # now = max edge moves
    shrink_frac = inpDict.get("bg_shrink_frac", 0.05)
    min_width   = inpDict.get("bg_min_sig_width", 1e-4)

    orig_sig_lo, orig_sig_hi = sig_lo, sig_hi
    shape_ok = True

    # 1) Check original window
    if not is_good_background_shape(fit_func, sig_lo, sig_hi, neg_tol=neg_tol):
        print(
            f"[bg_fit] initial shape bad for {hist.GetName()}: "
            f"window=[{sig_lo:.5f},{sig_hi:.5f}]"
        )

        # 2) Move any negative edges inward (no refit)
        new_lo, new_hi = move_edges_to_nonnegative(
            fit_func,
            sig_lo,
            sig_hi,
            neg_tol=neg_tol,
            max_iter=max_refit,
            step_frac=shrink_frac,
            min_width=min_width,
        )

        if new_lo is None:
            shape_ok = False
        else:
            sig_lo, sig_hi = new_lo, new_hi
            # 3) After edge adjustment, check shape once more
            if not is_good_background_shape(fit_func, sig_lo, sig_hi, neg_tol=neg_tol):
                print(
                    f"[bg_fit] shape still bad after edge adjustment for "
                    f"{hist.GetName()}: window=[{sig_lo:.5f},{sig_hi:.5f}]"
                )
                shape_ok = False
            elif (sig_lo != orig_sig_lo) or (sig_hi != orig_sig_hi):
                print(
                    f"[bg_fit] adjusted signal window for {hist.GetName()}: "
                    f"[{orig_sig_lo:.5f},{orig_sig_hi:.5f}] -> "
                    f"[{sig_lo:.5f},{sig_hi:.5f}]"
                )

    # 4) If shape is still bad, fall back to zero background and return
    if not shape_ok:
        f_min = float(fit_func.GetMinimum(orig_sig_lo, orig_sig_hi))
        f_max = float(fit_func.GetMaximum(orig_sig_lo, orig_sig_hi))

        hist_name = hist.GetName()
        parts = hist_name.split("_")

        try:
            tbin = int(parts[-2]) + 1
            phibin = int(parts[-1]) + 1
            print(
                f"Bad fit for: {hist_name}  "
                f"(tbin={tbin}, phibin={phibin})  "
                f"f_min={f_min:.6g}  f_max={f_max:.6g}"
            )
        except ValueError:
            try:
                tbin = int(parts[-1]) + 1
                print(
                    f"Bad fit for: {hist_name}  "
                    f"(tbin={tbin})  "
                    f"f_min={f_min:.6g}  f_max={f_max:.6g}"
                )
            except ValueError:
                print(
                    "ERROR!"
                    f" Bad fit for: {hist_name}  "
                    "\nClosing script..."
                )
                sys.exit(2)

        fit_func_zero = TF1("fit_func_zero_bad", "0", mm_min, mm_max)
        fit_vis = fit_func_zero.Clone(f"{hist.GetName()}_bg_vis")

        fit_hist_inrange = get_fit_histogram_padded(
            fit_func_zero,
            hist_mm_cut,
            mm_min,
            mm_max,
            n_pad=0
        )

        bg_par = 0.0
        f_sig = 1.0  # 0.0 remove bin, 1.0 dont apply subtraction

        return fit_hist_inrange, fit_vis, bg_par, f_sig

    # If we get here, fit_func is good on the (possibly adjusted) [sig_lo, sig_hi]
    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        
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
