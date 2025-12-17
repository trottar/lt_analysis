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
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21), 
                    }
                },   

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21), 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21), 
                    }
                },     

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21), 
                    }
                },      

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,                    
                    "sidebands": {                
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),                
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21), 
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
                        "right": (1.12, 1.17),   
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
                        "right": (1.12, 1.17),                 
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
                        "right": (1.12, 1.17),                 
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
                        "right": (1.12, 1.17),                 
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
                        "right": (1.12, 1.17),                 
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

import math

def is_good_background_shape(
    fit_func,
    x_min,
    x_max,
    *,
    neg_tol=0.25,
    pos_tol=0.0,                 # require f_max > pos_tol
    concavity_rel_tol=0.02,       # max deviation from a straight line, as a fraction of scale
    concavity_abs_tol=None,       # if set, overrides rel tol (same units as y)
    n_samples=64,
):
    """
    Sanity check on the *fit function itself*.

    Accept iff:
      - finite min/max on [x_min, x_max]
      - f_max > pos_tol
      - f_min >= -neg_tol
      - not "too concave": close to the secant (endpoint-to-endpoint line)
    """
    if not (x_max > x_min):
        return False

    # ROOT extremum search
    f_min = float(fit_func.GetMinimum(x_min, x_max))
    f_max = float(fit_func.GetMaximum(x_min, x_max))

    # Sample the function explicitly (including near the edges) to catch dips
    # that ROOT's extremum search can miss, especially right at the boundaries.
    span = x_max - x_min
    n_uniform = max(n_samples, 3)
    xs_uniform = [x_min + i * span / (n_uniform - 1) for i in range(n_uniform)]

    # Add extra edge-biased probes to enforce the negative tolerance near the ends
    edge_fracs = (0.001, 0.01)
    xs_minmax = xs_uniform + [
        x_min + frac * span for frac in edge_fracs
    ] + [
        x_max - frac * span for frac in edge_fracs
    ]
    xs_minmax = sorted(set(x for x in xs_minmax if x_min <= x <= x_max))

    y_at = {}
    for x in xs_minmax:
        y = float(fit_func.Eval(x))
        if not math.isfinite(y):
            return False
        y_at[x] = y

    sample_min = min(y_at.values())
    sample_max = max(y_at.values())
    f_min = min(f_min, sample_min)
    f_max = max(f_max, sample_max)

    if not (math.isfinite(f_min) and math.isfinite(f_max)):
        return False
    if f_max <= pos_tol:
        return False
    if f_min <= -neg_tol:
        return False

    # Concavity / curvature check: max deviation from the secant line
    if n_samples >= 3 and (concavity_rel_tol is not None or concavity_abs_tol is not None):
        xs = xs_uniform
        ys = [y_at[x] for x in xs]

        y0, y1 = ys[0], ys[-1]
        span = x_max - x_min
        max_dev = 0.0
        for x, y in zip(xs, ys):
            t = (x - x_min) / span
            y_lin = y0 + t * (y1 - y0)
            dev = abs(y - y_lin)
            if dev > max_dev:
                max_dev = dev

        if concavity_abs_tol is not None:
            if max_dev > concavity_abs_tol:
                return False
        else:
            scale = max(abs(f_max), abs(f_min), 1.0)
            if max_dev / scale > concavity_rel_tol:
                return False

    return True

def bg_integral_norm_and_err_from_cov(
        fit_func,
        fit_res,
        hist_ref,
        x_min,
        x_max,
        *,
        ref_binwidth,
        allow_negative=False
):
    """
    Returns (N_bg_norm, dN_bg_norm) where N_bg_norm matches the sum of bin-contents
    produced by get_fit_histogram_padded(..., n_pad=0) over [x_min, x_max].

    Error is propagated from the fit covariance:
        Var(N) = grad^T * Cov * grad
    where grad_i = dN/dp_i computed by finite differences.

    Robust to: null TFitResultPtr / missing covariance matrix.
    """
    ax = hist_ref.GetXaxis()
    nb = ax.GetNbins()
    bw = float(ref_binwidth) if (ref_binwidth is not None and ref_binwidth > 0.0) else float(ax.GetBinWidth(1))

    i_lo = max(1,  ax.FindBin(x_min))
    i_hi = min(nb, ax.FindBin(x_max))

    def eval_N():
        tot = 0.0
        for ib in range(i_lo, i_hi + 1):
            x_lo = ax.GetBinLowEdge(ib)
            x_hi = x_lo + ax.GetBinWidth(ib)
            y_pred = fit_func.Integral(x_lo, x_hi) / bw
            if (not allow_negative) and (y_pred < 0.0):
                y_pred = 0.0
            tot += y_pred
        return tot

    N0 = eval_N()

    npar = int(fit_func.GetNpar())

    # --- get covariance safely ---
    cov = None
    cov_ok = False

    # NOTE: in PyROOT, a null TFitResultPtr is *not* None but evaluates False
    if (fit_res is not None) and bool(fit_res):
        try:
            cov = fit_res.GetCovarianceMatrix()
            cov_ok = True
        except ReferenceError:
            cov_ok = False

        # If available, reject bad/empty covariance status
        if cov_ok and hasattr(fit_res, "CovMatrixStatus"):
            try:
                if int(fit_res.CovMatrixStatus()) <= 0:
                    cov_ok = False
            except Exception:
                pass

    # Fallback: diagonal covariance from TF1 parameter errors
    if not cov_ok:
        cov = ROOT.TMatrixDSym(npar)
        any_nonzero = False
        for i in range(npar):
            ei = float(fit_func.GetParError(i))
            if ei > 0.0:
                any_nonzero = True
            # set diagonal
            try:
                cov[i][i] = ei * ei
            except Exception:
                cov(i, i)  # touch to ensure exists
                cov[i][i] = ei * ei

        # If even par-errors are zero (no info), be conservative: 100% of N0
        if not any_nonzero:
            return N0, abs(N0)

    grads = np.zeros(npar, dtype=float)

    # finite-difference derivatives dN/dp_i
    for ip in range(npar):
        p0 = float(fit_func.GetParameter(ip))
        pe = float(fit_func.GetParError(ip))

        dp = 0.1 * pe if pe > 0.0 else 1e-4 * (abs(p0) + 1.0)
        dp = max(dp, 1e-8)

        fit_func.SetParameter(ip, p0 + dp)
        Np = eval_N()

        fit_func.SetParameter(ip, p0 - dp)
        Nm = eval_N()

        fit_func.SetParameter(ip, p0)  # restore
        grads[ip] = (Np - Nm) / (2.0 * dp)

    # Var(N) = g^T C g
    var = 0.0
    for i in range(npar):
        for j in range(npar):
            try:
                cij = float(cov[i][j])
            except Exception:
                cij = float(cov(i, j))
            var += grads[i] * cij * grads[j]

    if (not math.isfinite(var)) or (var < 0.0):
        var = 0.0

    return N0, math.sqrt(var)

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(
        phi_setting,
        inpDict,
        hist,
        hist_mm_cut=None,
        *,
        scaling=1.0,
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
        return fit_func, fit_vis, bg_par, f_sig, 0.0

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

    # ------------------------------ fit + shape rescue -------------------
    # We will:
    #   1) Fit the sidebands defined in BG_MODELS.
    #   2) Check the background shape ONLY between the BG_MODELS bounds:
    #          [sb_left[0], sb_right[1]]
    #      i.e. exactly the x–range where the green fit is drawn.
    #   3) If the fit goes negative at those bounds, pull those bounds inward
    #      and refit, until the function is above -neg_tol or the window
    #      collapses.
    #
    neg_tol     = inpDict.get("bg_neg_tol", 1e-4)      # how negative we tolerate
    max_refit   = inpDict.get("bg_max_refit", 2000)   # max refits per histogram
    shrink_frac = inpDict.get("bg_shrink_frac", 0.05)  # fraction of width to move per step
    min_width   = inpDict.get("bg_min_sig_width", 1e-10)  # minimum BG window width

    # Background SHAPE window = BG_MODELS bounds (this is what you care about)
    bg_lo = sb_left[0]
    bg_hi = sb_right[1]
    orig_bg_lo, orig_bg_hi = bg_lo, bg_hi

    fit_func = None
    fit_res = None  # keep TFitResultPtr for covariance

    for irefit in range(max_refit):
        # (Re)fit on the same sideband histogram each iteration
        fit_min, fit_max = sb_left[0], sb_right[1]
        fit_func = TF1("fit_func", model["func_expr"], fit_min, fit_max)
        fit_res_ptr = h_sb.Fit(fit_func, "SQ0")  # S: keep fit result (covariance) + quiet
        fit_res = None
        try:
            fit_res = fit_res_ptr.Get()  # nullptr if fit failed
            if not fit_res:
                fit_res = None
        except Exception:
            fit_res = None

        # Check the shape INSIDE THE BG_MODELS BOUNDS [bg_lo, bg_hi]
        if is_good_background_shape(fit_func, bg_lo, bg_hi, neg_tol=neg_tol):
            if irefit > 0:
                print(
                    f"[bg_fit] rescued background for {hist.GetName()}: "
                    f"BG window {orig_bg_lo:.3f}-{orig_bg_hi:.3f} "
                    f"-> {bg_lo:.3f}-{bg_hi:.3f} after {irefit} refits"
                )
            break  # good shape, keep this fit + BG window

        # Shape still bad ⇒ shrink the BG window and try again
        width = bg_hi - bg_lo
        if width <= min_width:
            print(
                f"[bg_fit] unable to rescue {hist.GetName()}: "
                f"BG window collapsed (width={width:.4g} <= {min_width:.4g})"
            )
            fit_func = None
            break

        f_lo = float(fit_func.Eval(bg_lo))
        f_hi = float(fit_func.Eval(bg_hi))
        step = shrink_frac * width

        print(
            f"[bg_fit] refit {irefit}: "
            f"BG window=[{bg_lo:.5f},{bg_hi:.5f}] width={width:.5g} "
            f"f_lo={f_lo:.5g} f_hi={f_hi:.5g} step={step:.5g}"
        )

        # If one edge is clearly more negative, pull that edge in.
        if (f_lo < -neg_tol) or (f_hi < -neg_tol):
            if f_lo <= f_hi:
                print(
                    f"  -> shrink lower BG edge: {bg_lo:.5f} -> {bg_lo + step:.5f} "
                    f"(f_lo={f_lo:.5g})"
                )
                bg_lo += step
            if f_hi < f_lo:
                print(
                    f"  -> shrink upper BG edge: {bg_hi:.5f} -> {bg_hi - step:.5f} "
                    f"(f_hi={f_hi:.5g})"
                )
                bg_hi -= step
        else:
            # Edges are OK but the interior in [bg_lo, bg_hi] is still bad:
            # shrink symmetrically.
            new_lo = bg_lo + 0.5 * step
            new_hi = bg_hi - 0.5 * step
            print(
                "  -> symmetric BG shrink: "
                f"{bg_lo:.5f}-{bg_hi:.5f} -> {new_lo:.5f}-{new_hi:.5f}"
            )
            bg_lo, bg_hi = new_lo, new_hi

    # After the loop, if we still have no acceptable fit, fall back to zero BG
    if fit_func is None or not is_good_background_shape(fit_func, bg_lo, bg_hi, neg_tol=neg_tol):
        # For debug, look at the shape on the ORIGINAL BG window if possible
        if fit_func is not None:
            f_min = float(fit_func.GetMinimum(orig_bg_lo, orig_bg_hi))
            f_max = float(fit_func.GetMaximum(orig_bg_lo, orig_bg_hi))
        else:
            f_min = float("nan")
            f_max = float("nan")

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

        # Zero background over the full MM range
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

        # Return NaN so calculate_yield flags this bin with large (100%) BG fractional uncertainty
        return fit_hist_inrange, fit_vis, bg_par, f_sig, float("nan")


    # If we get here, fit_func is good on the (possibly shrunken) BG window.
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
    if scaling != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scaling)
        N_bg  *= scaling
        bg_par *= scaling

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

    # --- background integral uncertainty (covariance-propagated) ---
    # matches how calculate_yield currently uses bg_hist.Integral()
    bw_for_norm = hist.GetBinWidth(1)
    if bw_for_norm <= 0.0:
        bw_for_norm = hist_mm_cut.GetXaxis().GetBinWidth(1)

    _, N_bg_norm_err = bg_integral_norm_and_err_from_cov(
        fit_func,
        fit_res,
        hist_mm_cut,
        mm_min,
        mm_max,
        ref_binwidth=bw_for_norm,
        allow_negative=False
    )

    # if you apply a global scaling to the fit, scale the uncertainty the same way
    if scaling != 1.0:
        N_bg_norm_err *= abs(scaling)

    return fit_hist_inrange, fit_vis, bg_par, f_sig, N_bg_norm_err
