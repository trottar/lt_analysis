#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-06-03 11:47:09 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import ROOT
from ROOT import TGraphErrors, TF1, TF2, TGraph2DErrors, TCanvas
from ROOT import TString, TNtuple, TMinuit
from array import array
import math
import ctypes
import os, sys

ParticleType = sys.argv[1]
POL = float(sys.argv[2])

if POL > 0:
    polID = 'pl'
else:
    polID = 'mn'

Q2 = sys.argv[3]
W = sys.argv[4]

LOEPS = float(sys.argv[5])
HIEPS = float(sys.argv[6])

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

outputpdf  = "{}/{}_lt_fit_Q{}W{}.pdf".format(OUTPATH, ParticleType, Q2, W)

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Constants
#pt_to_pt_systematic_error = 2.9 # Percent, just matching Bill's for now
pt_to_pt_systematic_error = 3.6 # In percent, matches PAC propsal projections (https://redmine.jlab.org/attachments/download/635/k12_proposal.pdf)
PI = math.pi

###############################################################################################################################################

# ------------------------------------------------------------------
TOL_AT_LIMIT = 1e-3          # how close to the hard wall counts as “at the limit”
PENALTY_LAMBDA = 1.0e3   # strength of χ² penalty once |p| exceeds its limit
# --- Adaptive-limit parameters -----------------------------
HIGH_PULL_FRAC  = 0.85   # when |σ| > 85 % of its limit, consider loosening
SCALE_UP_FACTOR = 1.30   # enlarge the limit by ×1.3 when it keeps getting hit
# --- Decorrelation parameters ------------------------------
CORR_THRESHOLD      = 0.92   # if |ρ(L,T)| > 0.92 ⇒ re-fit in Σ/Δ basis
MIN_CORR_IMPROVE    = 0.15   # require ≥ 0.15 reduction, else keep old fit
# --- Gaussian priors (LT & TT) -----------------------------
# Prior:  χ²_prior = ((p_i – μ_i)/σ_i)²
PRIOR_SIGMA_LT = 0.10   # nb; adjust to taste / external theory
PRIOR_SIGMA_TT = 0.10   # nb
ENFORCE_PRIOR  = True   # quick on/off switch
# === GLOBAL-FIT CONTROLS ===
USE_GLOBAL_FIT = True         # True → do one global fit; False → run per-bin fits (Steps 1–5)
PAR_MODEL      = "exp"        # "exp" for A·e^{B t}, or "poly" for A + B·t
N_PAR_MODEL    = 2            # number of parameters per σ_i(t) model
# === MINOS ERROR ANALYSIS ===
USE_MINOS = True    # True → run Minos on each fit to get asymmetric errors
# ===================================

# Per-parameter scaling factors (start at 1.0 = “base” limit)
LIMIT_SCALE = {
    'sigL' : 1.0,
    'sigT' : 1.0,
    'sigLT': 1.0,
    'sigTT': 1.0,
}

# Map TF1 parameter index -> key used in adapt_limits
IDX2KEY = {0: 'sigL', 1: 'sigT', 2: 'sigLT', 3: 'sigTT'}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def get_scaled_limit(key):
    """
    Return the current hard wall for 'sigL', 'sigT', 'sigLT', or 'sigTT',
    after applying the adaptive SCALE factor.
    """
    base_lim = adapt_limits(key)[1]
    return base_lim * LIMIT_SCALE[key]
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def refit_until_inside_limits(fit_func, graph, limit_map):
    """
    Iteratively fixes any parameter that sits on its hard limit and refits
    until all free parameters lie strictly inside their allowed range.
    `limit_map` must map parameter index -> (abs_limit_value).
    """
    while True:
        hit_any = False
        for idx, lim in limit_map.items():
            # Is this parameter already fixed?  If so, skip.
            if fit_func.IsFixed(idx):
                continue
            # Check distance to either +lim or –lim
            if abs(abs(fit_func.GetParameter(idx)) - lim) < TOL_AT_LIMIT:
                fit_func.FixParameter(idx, 0.0)
                hit_any = True
        if not hit_any:
            break
        graph.Fit(fit_func, "MRQ")   # silent refit
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def run_penalized_fit(graph, tf1, limit_map, lam=PENALTY_LAMBDA):
    """
    One pass of Minuit2 that minimises
        χ²_eff = χ²_data  +  lam * Σ max(0, |p_i| - L_i)².
    It re-uses the current TF1 as the model, updates its parameters
    in-place, and returns nothing (tf1 now holds the penalised minimum).
    """
    from ROOT import Math, TVectorD, TArrayD
    npar     = tf1.GetNpar()
    # -- build a callable χ² --
    def chi2(par_vec):
        # 1) regular χ² from the data points
        chi2_data = 0.0
        x  = ctypes.c_double()
        y  = ctypes.c_double()
        ex = ctypes.c_double()
        ey = ctypes.c_double()
        # copy parameters into the TF1 so EvalPar uses them
        for i in range(npar):
            tf1.SetParameter(i, par_vec[i])
        for i in range(graph.GetN()):
            graph.GetPoint(i, x, y)
            ey = graph.GetErrorY(i)
            if ey <= 0: ey = 1e-9
            model = tf1.Eval(x.value)
            chi2_data += ((y.value - model)/ey)**2
        # 2) quadratic soft wall
        soft_wall = 0.0
        for i, L_i in limit_map.items():
            excess = abs(par_vec[i]) - L_i
            if excess > 0.0:
                soft_wall += excess*excess
        # 3) optional Gaussian priors
        prior_pen = 0.0
        if ENFORCE_PRIOR:
            prior_pen += (par_vec[2] / PRIOR_SIGMA_LT)**2
            prior_pen += (par_vec[3] / PRIOR_SIGMA_TT)**2
        return chi2_data + lam*soft_wall + prior_pen

    # -- wrap in a ROOT.Math.Functor so Minuit2 will accept it --
    minim = Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    # directly give the Python chi2() and its dimensionality
    minim.SetFunction(chi2, npar)

    # initial guesses = current TF1 parameters
    for i in range(npar):
        minim.SetVariable(i,
                          tf1.GetParName(i),
                          tf1.GetParameter(i),
                          tf1.GetParError(i) if tf1.GetParError(i) > 0 else 0.1)
    minim.Minimize()

    # update the TF1 with the penalised best-fit values
    best = minim.X()
    for i in range(npar):
        tf1.SetParameter(i, best[i])
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def adapt_limits_after_fit(tf1):
    """
    Inspect the final parameters in tf1 and, if a value keeps pressing
    near its wall, enlarge that wall for future t-bins.
    """
    for idx, key in IDX2KEY.items():
        cur_lim  = get_scaled_limit(key)
        cur_val  = abs(tf1.GetParameter(idx))
        if cur_val > HIGH_PULL_FRAC * cur_lim:
            LIMIT_SCALE[key] *= SCALE_UP_FACTOR
            # Optional: print a one-line diagnostic
            print(f"[adaptive-limit] Enlarged {key} to {get_scaled_limit(key):.3g}")
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def lt_correlation(tf2):
    """
    Return the correlation coefficient ρ(σ_L, σ_T) for the current TF2 fit.
    Assumes parameter index 0 ↔ σ_T and 1 ↔ σ_L (the order used in fff2).
    """
    cov   = tf2.GetCovarianceMatrix()   # TMatrixDSym
    var_T = cov(0,0)
    var_L = cov(1,1)
    cov_TL= cov(0,1)
    return cov_TL / math.sqrt(var_T * var_L) if var_T*var_L > 0 else 0.0
# ------------------------------------------------------------------

# ------------------------------------------------------------------
def refit_in_sum_diff_basis(graph, tf2, limit_map):
    """
    Perform a second fit in the rotated basis:
        Σ = σ_L + σ_T        (param [0])
        Δ = σ_L – σ_T        (param [1])
        σ_LT, σ_TT unchanged ([2], [3])
    If the correlation improves by at least MIN_CORR_IMPROVE, keep the
    rotated result; otherwise revert to the original parameter set.
    """

    # -- Build a new TF2 with the Σ/Δ parameterisation -------------
    expr = ("0.5*((1+y)*[0] + (y-1)*[1])"                # Σ/Δ core
            " + sqrt(2*y*(1+y))*cos(x*0.017453)*[2]"     # σ_LT term
            " + y*cos(2*x*0.017453)*[3]")                # σ_TT term
    f_rot = ROOT.TF2("f_rot", expr, 0, 360, 0.0, 1.0)

    # -- Seed parameters from current fit --------------------------
    sigT = tf2.GetParameter(0)
    sigL = tf2.GetParameter(1)
    f_rot.SetParameters(sigL + sigT,          # Σ
                        sigL - sigT,          # Δ
                        tf2.GetParameter(2),  # σ_LT
                        tf2.GetParameter(3))  # σ_TT

    # Re-use *errors* to guide step sizes
    for i in range(4):
        f_rot.SetParError(i, tf2.GetParError(i))

    # -- Perform a quick fit with the same penalty limits ----------
    refit_until_inside_limits(f_rot, graph, limit_map)

    # -- Evaluate whether correlation improved ---------------------
    rho_before = abs(lt_correlation(tf2))
    rho_after  = abs(lt_correlation(f_rot))

    if rho_before - rho_after >= MIN_CORR_IMPROVE:
        # Accept rotated parameters – convert back to L & T
        sigSum = f_rot.GetParameter(0)
        sigDiff= f_rot.GetParameter(1)
        tf2.SetParameter(0, 0.5*(sigSum - sigDiff))   # σ_T
        tf2.SetParameter(1, 0.5*(sigSum + sigDiff))   # σ_L
        tf2.SetParameter(2, f_rot.GetParameter(2))    # σ_LT
        tf2.SetParameter(3, f_rot.GetParameter(3))    # σ_TT
        print(f"[decorrelate] ρ(L,T) improved {rho_before:.2f} → {rho_after:.2f}")
    else:
        print(f"[decorrelate] kept original fit (ρ {rho_before:.2f} → {rho_after:.2f})")
# ------------------------------------------------------------------

# === GLOBAL FIT FUNCTION ===
def run_global_fit(global_pts):
    """
    Simultaneous χ² fit to ALL (t, eps, φ, σ_unsep, σ_err) in global_pts.
    Models each σ_i(t) with PAR_MODEL (“exp” or “poly”), 
    using 2 parameters [A_i, B_i] for i = {L, T, LT, TT}.
    Returns [A_L, B_L, A_T, B_T, A_LT, B_LT, A_TT, B_TT].
    """
    import math
    from ROOT import Math
    npar = 4 * N_PAR_MODEL

    # model f(t; A, B)
    if PAR_MODEL == "exp":
        f = lambda A,B,t: A * math.exp(B*t)
    else:
        f = lambda A,B,t: A + B*t

    # χ² over all points
    def chi2(params):
        χ2 = 0.0
        # unpack in order [A_L,B_L, A_T,B_T, A_LT,B_LT, A_TT,B_TT]
        for t_val, eps, phi_deg, σ, σ_err in global_pts:
            AL, BL, AT, BT, ALT, BLT, ATT, BTT = params
            sigL  = f(AL, BL, t_val)
            sigT  = f(AT, BT, t_val)
            sigLT = f(ALT,BLT,t_val)
            sigTT = f(ATT,BTT,t_val)
            φ     = math.radians(phi_deg)
            model = (eps*sigL + sigT
                     + math.sqrt(2*eps*(1+eps))*sigLT*math.cos(φ)
                     + eps*sigTT*math.cos(2*φ))
            χ2   += ((σ - model)/σ_err)**2
        return χ2

    functor = Math.Functor(chi2, npar)
    minim   = Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    minim.SetFunction(functor)
    # initialize parameters
    for i in range(npar):
        minim.SetVariable(i, f"p{i}", 0.05, 0.01)
    minim.Minimize()
    # === STEP-7 GLOBAL MINOS CALLS ===
    if USE_MINOS:
        from ctypes import c_double
        # prepare containers
        lo_errs = [c_double(0) for _ in range(4*N_PAR_MODEL)]
        hi_errs = [c_double(0) for _ in range(4*N_PAR_MODEL)]
        for ipar in range(4 * N_PAR_MODEL):
            minim.Minos(ipar, lo_errs[ipar], hi_errs[ipar])
        # now lo_errs[ipar].value and hi_errs[ipar].value hold the asymmetric errors
        # you can return them alongside the best-fit parameters if you like
    # ===================================

    return list(minim.X())
# ===================================

###############################################################################################################################################
# -------------------------  DYNAMIC LIMITS  -------------------------------------------------
# Edit PARAM_LIMITS below to adjust allowed ranges per fit-step (0-based index).
# Provide either:
#   • one tuple  -> same limits every pass
#   • list/tuple of tuples -> individual limits per step; last entry reused if
#     the sequence is shorter than the maximum number of steps (7 here).
#
#'''
PARAM_LIMITS = {
    'sigT' : [(-5.0, 1000.0)]*7,   # parameter 0
    'sigL' : [(-5.0, 1000.0)]*7,   # parameter 1
    'sigLT': [(-5.0,  5.0)]*7,  # parameter 2
    'sigTT': [(-5.0,  5.0)]*7   # parameter 3
}
'''
PARAM_LIMITS = {
    'sigT' : [(-10.0, 1000.0)]*7,   # parameter 0
    'sigL' : [(-10.0, 1000.0)]*7,   # parameter 1
    'sigLT': [(-10.0, 100.0)]*7,  # parameter 2
    'sigTT': [(-10.0, 100.0)]*7   # parameter 3
}
'''

def adapt_limits(param_name_or_idx, step=0):
    """
    Return (low, high) limits for a given parameter and step index.

    param_name_or_idx : int 0-3 OR 'sigT'|'sigL'|'sigLT'|'sigTT'
    step              : 0-based fit pass counter
    """
    name_map = {0:'sigT', 1:'sigL', 2:'sigLT', 3:'sigTT'}
    pname = name_map.get(param_name_or_idx, param_name_or_idx)
    limdef = PARAM_LIMITS.get(pname, (0.0, 1.0e6))

    # If user supplied per-step list:
    if isinstance(limdef[0], (list, tuple)):
        if step < len(limdef):
            return limdef[step]
        return limdef[-1]
    # Single tuple provided:
    return limdef
# --------------------------------------------------------------------------------------------
###############################################################################################################################################

# Import separated xsects models
from lt_active import LT_sep_x_lo_fun_wrapper, LT_sep_x_lo_fun_unsep_wrapper, LT_sep_x_hi_fun_wrapper, LT_sep_x_hi_fun_unsep_wrapper

###############################################################################################################################################

def single_setting(q2_set, w_set, fn_lo, fn_hi):

    # Set epsilon for lt_active script
    #set_val(LOEPS, HIEPS)

    eps_diff = HIEPS-LOEPS
    
    sig_L_g  = TGraphErrors()
    sig_T_g  = TGraphErrors()
    sig_LT_g = TGraphErrors()
    sig_TT_g = TGraphErrors()

    sig_lo = TGraphErrors()
    sig_hi = TGraphErrors()
    sig_diff_g = TGraphErrors()
    
    nlo = TNtuple("nlo", "nlo", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nlo.ReadFile(fn_lo)
    
    nhi = TNtuple("nhi", "nhi", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nhi.ReadFile(fn_hi)

    q2_list = []
    w_list = []
    theta_list = []
    t_list = []
    lo_eps_list = []
    hi_eps_list = []

    for evt in nlo:
        if evt.t not in t_list:
            q2_list.append(evt.Q2)
            w_list.append(evt.w)
            theta_list.append(evt.theta)
            t_list.append(evt.t)
            lo_eps_list.append(evt.eps)

    tmp_t_list = []
    for evt in nhi:
        if evt.t not in tmp_t_list:
            tmp_t_list.append(evt.t)
            hi_eps_list.append(evt.eps)
        
    t_bin_num = len(t_list)

    lo_cross_sec = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec = np.zeros(t_bin_num, dtype=float)
    lo_cross_sec_err = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec_err = np.zeros(t_bin_num, dtype=float)

    # === STEP-6 GLOBAL ACCUMULATION ===
    if USE_GLOBAL_FIT:
        global_pts = []   # will collect tuples: (t_bin, eps, phi_deg, σ_unsep, σ_err)
    # ================================

    for i in range(0, t_bin_num):    
        
        print("\n/*--------------------------------------------------*/")
        print(" Starting t-bin {0} (t={1:.4f})...".format(i+1, float(t_list[i])))
        print("\n/*--------------------------------------------------*/\n\n")

        tcut = "t=={0} && x!=0.0".format(float(t_list[i]))
        print(tcut)
        
        lo_eps = lo_eps_list[i]
        hi_eps = hi_eps_list[i]

        nlo.Draw("x:phi:dx", tcut, "goff")

        glo_tmp = TGraphErrors()
        for j in range(nlo.GetSelectedRows()):
            glo_tmp.SetPoint(j, nlo.GetV2()[j], nlo.GetV1()[j])
            glo_tmp.SetPointError(j, 0, nlo.GetV3()[j])

        LT_sep_x_lo_fun = LT_sep_x_lo_fun_wrapper(lo_eps)
        flo = TF1("lo_eps_fit", LT_sep_x_lo_fun, 0, 360, 4)
        LT_sep_x_lo_fun_unsep = LT_sep_x_lo_fun_unsep_wrapper(lo_eps)
        flo_unsep = TF1("lo_eps_unsep", LT_sep_x_lo_fun_unsep, 0, 2*PI, 4)
        
        glo = glo_tmp.Clone("glo")
        ave_sig_lo = glo.GetMean(2)
        err_sig_lo = glo.GetRMS(2)

        sig_lo.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_lo.GetYaxis().SetTitle("#bar{#it{#sigma}_{Low}}")
        sig_lo.SetTitle("t = {:.3f}".format(t_list[i]))
        sig_lo.SetPoint(sig_lo.GetN(), float(t_list[i]), ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)
        
        nhi.Draw("x:phi:dx", tcut, "goff")

        ghi_tmp = TGraphErrors()
        for j in range(nhi.GetSelectedRows()):
            ghi_tmp.SetPoint(j, nhi.GetV2()[j], nhi.GetV1()[j])
            ghi_tmp.SetPointError(j, 0, nhi.GetV3()[j])

        LT_sep_x_hi_fun = LT_sep_x_hi_fun_wrapper(hi_eps)
        fhi = TF1("hi_eps_fit", LT_sep_x_hi_fun, 0, 360, 4)
        LT_sep_x_hi_fun_unsep = LT_sep_x_hi_fun_unsep_wrapper(hi_eps)
        fhi_unsep = TF1("hi_eps_unsep", LT_sep_x_hi_fun_unsep, 0, 2*PI, 4)
            
        ghi = ghi_tmp.Clone("ghi")
        ave_sig_hi = ghi.GetMean(2)
        err_sig_hi = ghi.GetRMS(2)

        sig_hi.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_hi.GetYaxis().SetTitle("#bar{#it{#sigma}_{High}}")
        sig_hi.SetTitle("t = {:.3f}".format(t_list[i]))
        sig_hi.SetPoint(sig_hi.GetN(), float(t_list[i]), ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)
        
        g_plot_err = TGraph2DErrors()
        g_xx, g_yy, g_yy_err = ctypes.c_double(0),ctypes.c_double(0),ctypes.c_double(0)

        for ii in range(glo.GetN()):
            glo.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((glo.GetErrorY(ii) / g_yy.value)**2 + (pt_to_pt_systematic_error/100)**2) * g_yy.value
            lo_cross_sec_err[i] += 1 / (g_yy_err**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, lo_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0,
                                     math.sqrt((glo.GetErrorY(ii))**2 + (pt_to_pt_systematic_error/100)**2))

        for ii in range(ghi.GetN()):
            ghi.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((ghi.GetErrorY(ii) / g_yy.value)**2 + (pt_to_pt_systematic_error/100)**2) * g_yy.value
            hi_cross_sec_err[i] += 1 / (g_yy_err**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, hi_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0,
                                     math.sqrt((ghi.GetErrorY(ii))**2 + (pt_to_pt_systematic_error/100)**2))

        try:
            lo_cross_sec_err[i] = 1/math.sqrt(lo_cross_sec_err[i])            
            hi_cross_sec_err[i] = 1/math.sqrt(hi_cross_sec_err[i])
        except ZeroDivisionError:
            lo_cross_sec_err[i] = -1000
            hi_cross_sec_err[i] = -1000
            
        g_plot_err.SetFillColor(29)
        g_plot_err.SetMarkerSize(0.8)
        g_plot_err.SetMarkerStyle(20)
        g_plot_err.SetMarkerColor(ROOT.kRed)
        g_plot_err.SetLineColor(ROOT.kBlue-3)
        g_plot_err.SetLineWidth(2)
        
        # === STEP-6 COLLECT & SKIP PER-BIN (fixed for TGraph2DErrors) ===
        if USE_GLOBAL_FIT:
            for idx_pt, eps_val in enumerate((lo_eps, hi_eps)):
                # prepare three C-doubles: x (φ), y (ε), z (σ)
                x_cd = ctypes.c_double()
                y_cd = ctypes.c_double()
                z_cd = ctypes.c_double()

                # note: four args here
                g_plot_err.GetPoint(idx_pt, x_cd, y_cd, z_cd)

                # pull the z-error (σ uncertainty)
                err_z = g_plot_err.GetErrorZ(idx_pt)

                # append (t, ε, φ, σ, σ_err)
                global_pts.append((
                    float(t_list[i]),  # t-bin center
                    eps_val,           # ε (lo or hi)
                    x_cd.value,        # φ (deg)
                    z_cd.value,        # σ_unsep
                    err_z              # σ uncertainty
                ))
        # ================================================================

        fff2 = TF2("fff2",
                "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]",
                0, 360, 0.0, 1.0)

        '''
        fff2 = TF2("fff2",
                "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]",
                0, 360, LOEPS-0.1, HIEPS+0.1)
        '''

        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()
        sigLT_change = TGraphErrors()
        sigTT_change = TGraphErrors()

        # ---------------- FIT SEQUENCE ------------------
        fit_step = 0  # counter for adapt_limits

        # --- Fit 1: L & T ---
        #fff2.SetParameter(0, 1)
        low, high = adapt_limits(0, fit_step)
        fff2.SetParLimits(0, low, high)

        #fff2.SetParameter(1, 1)
        low, high = adapt_limits(1, fit_step)
        fff2.SetParLimits(1, low, high)

        fff2.FixParameter(2, 0.0)
        fff2.FixParameter(3, 0.0)

        g_plot_err.Fit(fff2, "MRQ")

        sigL_change.SetTitle("t = {:.3f}".format(t_list[i]))
        sigL_change.GetXaxis().SetTitle("Fit Step")
        sigL_change.GetYaxis().SetTitle("#it{#sigma}_{L}")

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))

        sigT_change.SetTitle("t = {:.3f}".format(t_list[i]))
        sigT_change.GetXaxis().SetTitle("Fit Step")
        sigT_change.GetYaxis().SetTitle("#it{#sigma}_{T}")
        
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 2: LT ---
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(3, fff2.GetParameter(3))

        fff2.ReleaseParameter(2)
        fff2.SetParameter(2, 0.0)
        low, high = adapt_limits(2, fit_step)
        fff2.SetParLimits(2, low, high)

        g_plot_err.Fit(fff2, "MRQ")

        # -----------------------------------------------------------
        # Soft-constraint refinement
        # -----------------------------------------------------------
        limit_map = {
            0: get_scaled_limit('sigL'),
            1: get_scaled_limit('sigT'),
            2: get_scaled_limit('sigLT'),
            3: get_scaled_limit('sigTT'),
        }
        run_penalized_fit(g_plot_err, fff2, limit_map)
        print(f"[priors] σ_LT = {fff2.GetParameter(2):+.3f} nb   "
          f"σ_TT = {fff2.GetParameter(3):+.3f} nb")
        # -----------------------------------------------------------

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 3: L & T again ---
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)
        fff2.SetParameter(0, fff2.GetParameter(0))
        fff2.SetParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(2, fff2.GetParameter(2))
        fff2.FixParameter(3, fff2.GetParameter(3))

        g_plot_err.Fit(fff2, "MRQ")

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 4: TT ---
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(2, fff2.GetParameter(2))

        fff2.ReleaseParameter(3)
        fff2.SetParameter(3, 0.0)
        low, high = adapt_limits(3, fit_step)
        fff2.SetParLimits(3, low, high)

        g_plot_err.Fit(fff2, "MRQ")

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 5: L & T once more ---
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)
        fff2.SetParameter(0, fff2.GetParameter(0))
        fff2.SetParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(2, fff2.GetParameter(2))
        fff2.FixParameter(3, fff2.GetParameter(3))

        g_plot_err.Fit(fff2, "MRQ")

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 6: LT & TT together ---
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))

        fff2.ReleaseParameter(2)
        fff2.ReleaseParameter(3)
        fff2.SetParameter(2, fff2.GetParameter(2))
        fff2.SetParameter(3, fff2.GetParameter(3))

        g_plot_err.Fit("fff2", "MRQ")
        
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1

        # --- Fit 7: All free ---
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)
        fff2.ReleaseParameter(2)
        fff2.ReleaseParameter(3)

        low, high = adapt_limits(0, fit_step)
        fff2.SetParLimits(0, low, high)
        low, high = adapt_limits(1, fit_step)
        fff2.SetParLimits(1, low, high)
        low, high = adapt_limits(2, fit_step)
        fff2.SetParLimits(2, low, high)
        low, high = adapt_limits(3, fit_step)
        fff2.SetParLimits(3, low, high)

        # --- “All-free” fit, now capturing the result ---
        fitres = g_plot_err.Fit(fff2, "MRQS")   # “S” → store and return TFitResultPtr

        # === PER-BIN MINOS CALLS ===
        if USE_MINOS:
            # param order in fff2: 0→σ_T, 1→σ_L, 2→σ_LT, 3→σ_TT
            # run Minos and grab lower/upper for each
            sig_t_lo, sig_t_hi   = fitres.LowerError(0), fitres.UpperError(0)
            sig_l_lo, sig_l_hi   = fitres.LowerError(1), fitres.UpperError(1)
            sig_lt_lo, sig_lt_hi = fitres.LowerError(2), fitres.UpperError(2)
            sig_tt_lo, sig_tt_hi = fitres.LowerError(3), fitres.UpperError(3)
        # ===================================

        # -----------------------------------------------------------
        # Iterative limit enforcement for all four cross-section terms
        # -----------------------------------------------------------
        limit_map = {
            0: get_scaled_limit('sigL'),
            1: get_scaled_limit('sigT'),
            2: get_scaled_limit('sigLT'),
            3: get_scaled_limit('sigTT'),
        }
        refit_until_inside_limits(fff2, g_plot_err, limit_map)
        # -----------------------------------------------------------

        # -----------------------------------------------------------
        # Record pulls and adapt future hard walls
        # -----------------------------------------------------------
        adapt_limits_after_fit(fff2)
        # -----------------------------------------------------------

        # -----------------------------------------------------------
        # Re-fit in Σ/Δ basis if L & T are highly correlated
        # -----------------------------------------------------------
        if abs(lt_correlation(fff2)) > CORR_THRESHOLD:
            refit_in_sum_diff_basis(g_plot_err, fff2, limit_map)

        c1 =  TCanvas()

        c1.Update()
                
        # Print plots for c1 canvases
        g_plot_err.Draw("perr")
        g_plot_err.SetTitle("t = {:.3f}".format(t_list[i]))
        g_plot_err.GetXaxis().SetTitle("#it{#phi} [degree]")
        g_plot_err.GetXaxis().CenterTitle()
        g_plot_err.GetXaxis().SetTitleOffset(1.5)
        g_plot_err.GetYaxis().SetTitle("Epsilon")
        g_plot_err.GetYaxis().CenterTitle()
        g_plot_err.GetYaxis().SetTitleOffset(1.5)
        g_plot_err.GetZaxis().SetTitle("#it{#sigma}")
        g_plot_err.GetZaxis().CenterTitle()
        g_plot_err.GetZaxis().SetTitleOffset(1.5)
        if i == 0:
            c1.Print(outputpdf+'(')
        else:
            c1.Print(outputpdf)
        c1.Clear()

        c2 =  TCanvas()
                
        c2.Update()
        c2.cd()
        
        glo.SetMarkerStyle(5)
        glo.GetXaxis().SetLimits(0, 360)

        ghi.SetMarkerColor(2)
        ghi.SetLineColor(2)
        ghi.SetMarkerStyle(4)

        g = ROOT.TMultiGraph()
        g.Add(glo)
        g.Add(ghi)

        g.Draw("AP")

        g.GetYaxis().SetTitle("Unseparated Cross Section [nb/GeV^{2}]")
        g.GetYaxis().CenterTitle()
        g.GetYaxis().SetTitleOffset(1.2)
        g.GetYaxis().SetTitleSize(0.04)

        g.GetXaxis().SetTitle("#it{#phi} [degree]")
        g.GetXaxis().CenterTitle()
        g.GetXaxis().SetLimits(0, 360)
        g.GetXaxis().SetTitleSize(0.04)
        
        c2.Update()

        flo.FixParameter(0, fff2.GetParameter(0))
        flo.FixParameter(1, fff2.GetParameter(1))
        flo.FixParameter(2, fff2.GetParameter(2))
        flo.FixParameter(3, fff2.GetParameter(3))

        flo_unsep.FixParameter(0, fff2.GetParameter(0))
        flo_unsep.FixParameter(1, fff2.GetParameter(1))
        flo_unsep.FixParameter(2, fff2.GetParameter(2))
        flo_unsep.FixParameter(3, fff2.GetParameter(3))

        fhi.FixParameter(0, fff2.GetParameter(0))
        fhi.FixParameter(1, fff2.GetParameter(1))
        fhi.FixParameter(2, fff2.GetParameter(2))
        fhi.FixParameter(3, fff2.GetParameter(3))

        fhi_unsep.FixParameter(0, fff2.GetParameter(0))
        fhi_unsep.FixParameter(1, fff2.GetParameter(1))
        fhi_unsep.FixParameter(2, fff2.GetParameter(2))
        fhi_unsep.FixParameter(3, fff2.GetParameter(3))

        glo.Fit(flo, "MRQ")
        ghi.Fit(fhi, "MRQ")
        
        flo.SetLineColor(1)
        fhi.SetLineColor(2)
        flo.SetLineWidth(2)
        fhi.SetLineWidth(2)
        fhi.SetLineStyle(2)

        ghi.SetLineColor(2)

        flo.Draw("same")
        fhi.Draw("same")

        leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
        leg.SetFillColor(0)
        leg.SetMargin(0.4)
        leg.SetTextSize(0.04)
        leg.AddEntry(glo, "Low #it{#font[120]{e}}", "pl")
        leg.AddEntry(ghi, "High #it{#font[120]{e}}", "pl")
        leg.Draw()

        flo_status = flo.GetNDF()
        fhi_status = fhi.GetNDF()        
        
        fit_status = ROOT.TLatex()
        fit_status.SetTextSize(0.04)
        fit_status.DrawTextNDC(0.15, 0.85, "t={:.3f}, Q2={:.1f}, W={:.2f}".format(
            t_list[i], float(q2_set.replace("p",".")), float(w_set.replace("p","."))))

        if ghi.GetMaximum() > glo.GetMaximum():
            glo.SetMaximum(ghi.GetMaximum() * 1.1)

        if ghi.GetMinimum() < glo.GetMinimum():
            glo.SetMinimum(ghi.GetMinimum() * 0.9)

        sig_diff_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_diff_g.GetYaxis().SetTitle("#Delta#sigma_{ave}/#Delta#epsilon")

        try:
            sig_diff = (ave_sig_hi-ave_sig_lo)/eps_diff
            sig_diff_g.SetPoint(sig_diff_g.GetN(), float(t_list[i]), sig_diff)
            sig_diff_err = sig_diff*math.sqrt((err_sig_hi/ave_sig_hi)**2+(err_sig_lo/ave_sig_lo)**2)
            sig_diff_g.SetPointError(sig_diff_g.GetN()-1, 0, sig_diff_err)
        except ZeroDivisionError:
            sig_diff = 0.0
            sig_diff_g.SetPoint(sig_diff_g.GetN(), float(t_list[i]), sig_diff)
            sig_diff_err = 0.0
            sig_diff_g.SetPointError(sig_diff_g.GetN()-1, 0, sig_diff_err)            
            
        sig_l, sig_t, sig_lt, sig_tt = (fff2.GetParameter(1), fff2.GetParameter(0),
                                        fff2.GetParameter(2), fff2.GetParameter(3))
        #sig_l_err, sig_t_err, sig_lt_err, sig_tt_err = (fff2.GetParError(1), fff2.GetParError(0),
        #                                                fff2.GetParError(2), fff2.GetParError(3))
        sig_t_err   = (sig_t_hi + sig_t_lo)/2
        sig_l_err   = (sig_l_hi + sig_l_lo)/2
        sig_lt_err  = (sig_lt_hi + sig_lt_lo)/2
        sig_tt_err  = (sig_tt_hi + sig_tt_lo)/2

        print("\nBin {}: Outputting...  ".format(i+1), "sig_l: ", sig_l, "sig_t: ", sig_t, \
              "sig_lt: ", sig_lt, "sig_tt: ", sig_tt, \
              "t: ", t_list[i], "theta: ", theta_list[i], "W: ", w_list[i], "Q2:", q2_list[i], \
              "eps_lo: ", lo_eps_list[i], "eps_hi: ", hi_eps_list[i])

        fn_sep = "{}/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(
            LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""))
        try:
            mode = 'w' if i == 0 else 'a'
            with open(fn_sep, mode) as f:
                f.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                    sig_l, sig_l_err, sig_t, sig_t_err, sig_lt, sig_lt_err, sig_tt, sig_tt_err,
                    fff2.GetChisquare(), t_list[i], w_list[i], q2_list[i], theta_list[i]))
        except IOError:
            print("Error writing to file {}.".format(fn_sep))
            
        del g_plot_err
        
        g_sig_l_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_l_total.GetYaxis().SetTitle("#it{#sigma}_{L} [nb/GeV^{2}]")
        
        g_sig_t_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_t_total.GetYaxis().SetTitle("#it{#sigma}_{T} [nb/GeV^{2}]")

        g_sig_lt_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_lt_total.GetYaxis().SetTitle("#it{#sigma}_{LT} [nb/GeV^{2}]")

        g_sig_tt_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_tt_total.GetYaxis().SetTitle("#it{#sigma}_{TT} [nb/GeV^{2}]")
        
        g_sig_l_total.SetPoint(g_sig_l_total.GetN(), float(t_list[i]), sig_l)
        g_sig_l_total.SetPointError(g_sig_l_total.GetN()-1, 0, sig_l_err)

        g_sig_t_total.SetPoint(g_sig_t_total.GetN(), float(t_list[i]), sig_t)
        g_sig_t_total.SetPointError(g_sig_t_total.GetN()-1, 0, sig_t_err)

        g_sig_lt_total.SetPoint(g_sig_lt_total.GetN(), float(t_list[i]), sig_lt)
        g_sig_lt_total.SetPointError(g_sig_lt_total.GetN()-1, 0, sig_lt_err)

        g_sig_tt_total.SetPoint(g_sig_tt_total.GetN(), float(t_list[i]), sig_tt)
        g_sig_tt_total.SetPointError(g_sig_tt_total.GetN()-1, 0, sig_tt_err)

        sig_L_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_L_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
        sig_L_g.SetTitle("t = {:.3f}".format(t_list[i]))
        
        sig_T_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_T_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
        sig_T_g.SetTitle("t = {:.3f}".format(t_list[i]))

        sig_LT_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_LT_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
        sig_LT_g.SetTitle("t = {:.3f}".format(t_list[i]))

        sig_TT_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_TT_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
        sig_TT_g.SetTitle("t = {:.3f}".format(t_list[i]))        
        
        sig_L_g.SetPoint(i, float(t_list[i]), sig_l)
        sig_T_g.SetPoint(i, float(t_list[i]), sig_t)
        sig_LT_g.SetPoint(i, float(t_list[i]), sig_lt)
        sig_TT_g.SetPoint(i, float(t_list[i]), sig_tt)

        sig_L_g.SetPointError(i, 0.0, sig_l_err)
        sig_T_g.SetPointError(i, 0.0, sig_t_err)
        sig_LT_g.SetPointError(i, 0.0, sig_lt_err)
        sig_TT_g.SetPointError(i, 0.0, sig_tt_err)

        c4 = TCanvas()

        sigL_change.Draw("a*")
        c4.Print(outputpdf)

        sigT_change.Draw("a*")
        c4.Print(outputpdf)    
        c4.Clear()

        c2.SetTopMargin(0.03)
        c2.SetRightMargin(0.03)
        c2.Print(outputpdf)
        c2.Clear()
            
        c3 = TCanvas()

        sig_L_g.Draw("a*")
        c3.Print(outputpdf)

        sig_T_g.Draw("a*")
        c3.Print(outputpdf)

        sig_LT_g.Draw("a*")
        c3.Print(outputpdf)

        sig_TT_g.Draw("a*")
        c3.Print(outputpdf)
        c3.Clear()

        c5 = TCanvas()

        sig_lo.Draw("a*")
        c5.Print(outputpdf)
        
        sig_hi.Draw("a*")
        c5.Print(outputpdf)

        sig_diff_g.Draw("a*")
        c5.Print(outputpdf)
        c5.Clear()
        
        lo_cross_sec[i] = flo_unsep.Integral(0, 2*PI) / (2*PI)
        hi_cross_sec[i] = fhi_unsep.Integral(0, 2*PI) / (2*PI)
        
        g_unsep_lo.SetPoint(g_unsep_lo.GetN(), LOEPS, lo_cross_sec[i])
        g_unsep_lo.SetPointError(g_unsep_lo.GetN()-1, 0, lo_cross_sec_err[i])        
        
        g_unsep_hi.SetPoint(g_unsep_hi.GetN(), HIEPS, hi_cross_sec[i])
        g_unsep_hi.SetPointError(g_unsep_hi.GetN()-1, 0, hi_cross_sec_err[i])                
        
        fn_unsep = "{}/src/{}/xsects/unsep_Q{}W{}.csv".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))
        try:
            mode = 'w' if i == 0 else 'a'
            with open(fn_unsep, mode) as f:
                f.write("{:.3f} {:.3f} {:.3f} {:.3f} {:.4f}\n".format(
                    lo_cross_sec[i], lo_cross_sec_err[i], hi_cross_sec[i], hi_cross_sec_err[i], t_list[i]))
        except IOError:
            print("Error writing to file {}.".format(fn_unsep))

        del c1, c2, c3, c4, c5

    # === EXECUTE GLOBAL FIT & FILL TOTALS ===
    if USE_GLOBAL_FIT:
        # 1) run the fit
        best = run_global_fit(global_pts)
        AL, BL, AT, BT, ALT, BLT, ATT, BTT = best

        # 2) (Re)initialize the TOTAL graphs so they're empty:
        g_sig_l_total.Reset()
        g_sig_t_total.Reset()
        g_sig_lt_total.Reset()
        g_sig_tt_total.Reset()

        # 3) evaluate at each t and fill
        for i_t, t_val in enumerate(t_list):
            if PAR_MODEL == "exp":
                sigL  = AL  * math.exp(BL  * t_val)
                sigT  = AT  * math.exp(BT  * t_val)
                sigLT = ALT * math.exp(BLT * t_val)
                sigTT = ATT * math.exp(BTT * t_val)
            else:
                sigL  = AL  + BL  * t_val
                sigT  = AT  + BT  * t_val
                sigLT = ALT + BLT * t_val
                sigTT = ATT + BTT * t_val

            # fill exactly where your per-bin code used to
            g_sig_l_total.SetPoint(g_sig_l_total.GetN(), float(t_val), sigL)
            g_sig_l_total.SetPointError(g_sig_l_total.GetN()-1, 0, sig_l_err)   # reuse per-bin err arrays if desired

            g_sig_t_total.SetPoint(g_sig_t_total.GetN(), float(t_val), sigT)
            g_sig_t_total.SetPointError(g_sig_t_total.GetN()-1, 0, sig_t_err)

            g_sig_lt_total.SetPoint(g_sig_lt_total.GetN(), float(t_val), sigLT)
            g_sig_lt_total.SetPointError(g_sig_lt_total.GetN()-1, 0, sig_lt_err)

            g_sig_tt_total.SetPoint(g_sig_tt_total.GetN(), float(t_val), sigTT)
            g_sig_tt_total.SetPointError(g_sig_tt_total.GetN()-1, 0, sig_tt_err)
    # ================================================    

    return t_list

fn_lo =  "{}/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(
            LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100)
fn_hi =  "{}/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(
            LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100)

g_sig_l_total = TGraphErrors()
g_sig_t_total = TGraphErrors()
g_sig_lt_total = TGraphErrors()
g_sig_tt_total = TGraphErrors()

g_unsep_lo = TGraphErrors()
g_unsep_hi = TGraphErrors()

t_list = single_setting(Q2, W, fn_lo, fn_hi) # Main function that performs fitting

f_lin_l = TF1("f_lin_l", "[0]+[1]*x", 0.0, 1.0)
f_lin_t = TF1("f_lin_t", "[0]+[1]*x", 0.0, 1.0)

num_events = g_unsep_lo.GetN()

for i in range(num_events):

    c = ROOT.TCanvas("c", "c", 800, 600)
    
    g_unsep_mult = ROOT.TMultiGraph()
    
    x_lo, y_lo = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_lo.GetPoint(i, x_lo, y_lo)
    x_err_lo = g_unsep_lo.GetErrorX(i)
    y_err_lo = g_unsep_lo.GetErrorY(i)
    g_lo_event = ROOT.TGraphErrors(1, array('d', [x_lo.value]), array('d', [y_lo.value]),
                                   array('d', [x_err_lo]), array('d', [y_err_lo]))
    g_unsep_mult.Add(g_lo_event)
    
    x_hi, y_hi = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_hi.GetPoint(i, x_hi, y_hi)
    x_err_hi = g_unsep_hi.GetErrorX(i)
    y_err_hi = g_unsep_hi.GetErrorY(i)
    g_hi_event = ROOT.TGraphErrors(1, array('d', [x_hi.value]), array('d', [y_hi.value]),
                                   array('d', [x_err_hi]), array('d', [y_err_hi]))
    g_unsep_mult.Add(g_hi_event)

    g_lo_event.SetMarkerStyle(5)
    g_hi_event.SetMarkerStyle(5)
    
    g_unsep_mult.Draw("AP")
    
    g_unsep_mult.GetYaxis().SetTitle("Unseparated Cross Section [nb/GeV^{2}]")
    g_unsep_mult.GetYaxis().SetTitleOffset(1.2)
    g_unsep_mult.GetXaxis().SetTitle("#epsilon")
    g_unsep_mult.GetXaxis().SetTitleOffset(1.2)
    
    f_lin = ROOT.TF1("f_lin", "[0]*x + [1]", 0, 1)
    g_unsep_mult.Fit(f_lin, "MRQ")
        
    f_lin.SetLineColor(2)
    f_lin.SetLineWidth(2)    
    f_lin.Draw("same")
    
    leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
    leg.SetFillColor(0)
    leg.SetMargin(0.4)
    leg.AddEntry(g_hi_event, "t = {:.3f}".format(t_list[i]), "")
    leg.AddEntry(g_lo_event, "Low #it{#font[120]{e}}", "p")
    leg.AddEntry(g_hi_event, "High #it{#font[120]{e}}", "p")
    leg.Draw()
    
    c.Print(outputpdf)
    c.Clear()

    del c
    
f_exp_l = TF1("f_exp_l", "[0]*exp(-[1]*x)", 0.0, 2.0)
f_exp_t = TF1("f_exp_t", "[0]*exp(-[1]*x)", 0.0, 2.0)

c_total_l_t = TCanvas()

g_sig_mult = ROOT.TMultiGraph()
g_sig_mult.Add(g_sig_l_total)
g_sig_mult.Add(g_sig_t_total)

g_sig_mult.Draw("AP")

g_sig_mult.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right) [nb/GeV^{2}]")
g_sig_mult.GetYaxis().SetTitleOffset(1.2)

g_sig_mult.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
g_sig_mult.GetXaxis().SetTitleOffset(1.2)

g_sig_l_total.Fit(f_exp_l, "MRQ")
g_sig_t_total.Fit(f_exp_t, "MRQ")

g_sig_l_total.SetLineColor(1)
g_sig_l_total.SetMarkerStyle(5)
g_sig_t_total.SetLineColor(2)
g_sig_t_total.SetMarkerColor(2)
g_sig_t_total.SetMarkerStyle(4)

f_exp_l.SetLineColor(1)
f_exp_l.SetLineWidth(2)
f_exp_t.SetLineColor(2)
f_exp_t.SetLineWidth(2)
f_exp_t.SetLineStyle(2)

f_exp_l.Draw("same")
f_exp_t.Draw("same")

leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
leg.SetFillColor(0)
leg.SetMargin(0.4)
leg.AddEntry(g_sig_l_total, "#it{#sigma}_{L}", "p")
leg.AddEntry(g_sig_t_total, "#it{#sigma}_{T}", "p")
leg.Draw()

c_total_l_t.Print(outputpdf)
c_total_l_t.Clear()

ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetStatW(0.15)
ROOT.gStyle.SetStatH(0.15)

f_exp = TF1("f_exp", "[0]*exp(-[1]*x)", 0.0, 2.0)

c_total = TCanvas()

g_sig_l_total.Draw("A*")
g_sig_l_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_t_total.SetMarkerColor(1)
g_sig_t_total.SetLineColor(1)
g_sig_t_total.Draw("A*")
g_sig_t_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_lt_total.Draw("A*")
g_sig_lt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_tt_total.Draw("A*")
g_sig_tt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf+')')
c_total.Clear()
