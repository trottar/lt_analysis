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
# ---------------------------  DYNAMIC LIMITS  ---------------------------------
#  PARAM_LIMITS encodes the *physical* boundaries that each cross-section term
#  may occupy and lets you tighten / loosen them at any individual fit-pass.
#
#  ▸ 3 fit-passes → index 0,1,2.  Fewer tuples than passes?  The last tuple
#    is reused for all remaining passes.
#
#  ▸ σT  , σL   (nb)        :  must be ≥ 0 by definition; an upper cap of
#                              1 000 nb is safely above any Hall-C kaon data.
#
#  ▸ ρ_LT, ρ_TT (dimension-less ratios) obey positivity:
#        |σ_LT| ≤ √(σT σL)  →  |ρ_LT| ≤ 1
#        |σ_TT| ≤    σT     →  |ρ_TT| ≤ 1
#    so the natural, model-independent range is [−1 … +1].
#
#  If you ever run a special kinematic slice (e.g. t near the kaon pole) and
#  want to *pre-shrink* the search volume at a given pass, just replace the
#  relevant tuple, e.g.
#        "sigL": [(0.001, 1e3), (0.001, 500), (0.001, 200)]
#  to narrow σL after the first pass.
# ------------------------------------------------------------------------------
PARAM_LIMITS = {
    "sigT" : [(0.001, 1e3)]*3,   # σ_T  : transverse
    "sigL" : [(0.001, 1e3)]*3,   # σ_L  : longitudinal
    "rhoLT": [(-1.0, 1.0)]*3,    # ρ_LT : σ_LT / √(σT σL)
    "rhoTT": [(-1.0, 1.0)]*3,     # ρ_TT : σ_TT / σT
    "Norm" : [(0.10, 10.0)]*3          # ← new
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------
# Initial starting values for the fit (nb)
# ------------------------------------------------------------------
SEED_SIGT = 10.0      # seed for σ_T  — typical transverse scale
SEED_SIGL = 100.0      # seed for σ_L  — ~50 % of σ_T so pole-dominance is reachable
# ------------------------------------------------------------------

# ------------------------------------------------------------------
#  ε-lever-arm quality guard
#  -------------------------------------------------
#  For two-epsilon LT separation the design matrix is
#         A = [[1, ε_lo],
#              [1, ε_hi]] .
#  Its condition number  κ ≈ 1/|ε_hi – ε_lo|.
#  If κ grows much beyond ~20 the σ_L extraction becomes
#  numerically unstable and MINUIT tends to peg σ_L at its
#  bound.  We use COND_MAX as the cut-off: when κ > COND_MAX
#  we apply the soft floor to σ_L (see Pass-2 code).
COND_MAX   = 20.0

# ---------------------------------------------------------------
def reset_limits_from_table(func, idx, key, stage):
    """
    Reset parameter *idx* of TF2 *func* to the lower/upper limits
    stored in PARAM_LIMITS[key][stage], and give it a non-zero step.
    """
    lo, hi = PARAM_LIMITS[key][stage]
    func.SetParLimits(idx, lo, hi)

    # give MINUIT a first step = 5 % of the allowed range,
    # or a small absolute step if the range is tiny
    step = 0.05 * (hi - lo) if hi > lo else 0.02
    func.SetParError(idx, step)

# ---------------------------------------------------------------
def get_limits(fcn, idx):
    """
    Return (lo, hi) limits for parameter *idx* of TF1/TF2 *fcn*.
    Works with both old (3-arg) and new (tuple) PyROOT signatures.
    """
    try:                                        # modern PyROOT (tuple)
        return fcn.GetParLimits(idx)
    except TypeError:                           # old signature → use refs
        lo_ref = ctypes.c_double(0.0)
        hi_ref = ctypes.c_double(0.0)
        fcn.GetParLimits(idx, lo_ref, hi_ref)
        return lo_ref.value, hi_ref.value

# --------------------------------------------------------------------------------------------
def penalty(f):
    sigT, sigL, rhoLT, rhoTT = (f.GetParameter(i) for i in range(4))
    p = 0.0
    # soft quadratic walls once the limits are exceeded
    if abs(rhoLT) > math.sqrt(max(sigT*sigL,0)):
        diff = abs(rhoLT) - math.sqrt(max(sigT*sigL,0))
        p += (diff / rhoLT_err)**2          # scale by current error
    if abs(rhoTT) > sigT:
        diff = abs(rhoTT) - sigT
        p += (diff / rhoTT_err)**2
    return p

# ---------------------------------------------------------------
# Positivity guard + auto-refit (robust version)
# ---------------------------------------------------------------
def check_sigma_positive(fcn, graph,
                         phi_set=(0, 90, 180, 270),
                         eps_set=(LOEPS, HIEPS),
                         shrink_factor=0.90):
    """
    Ensure σ(φ,ε) ≥ 0 for φ-sample × ε-set.
    If violated, shrink ρTT limits by `shrink_factor` and refit once.
    """

    def _is_positive():
        return min(fcn.Eval(phi, eps) for phi in phi_set for eps in eps_set) >= 0

    # early exit if already fine
    if _is_positive():
        return True

    print("WARNING: negative σ detected → tightening ρTT limits and refitting")

    # -----------------------------------------------------------
    # Fetch current limits on parameter 3 (ρTT)  –  two ways
    # -----------------------------------------------------------
    try:
        lo, hi = fcn.GetParLimits(3)        # modern PyROOT returns tuple
    except TypeError:
        # fallback for older signatures that need C-style refs
        lo_ref = ctypes.c_double(0.0)
        hi_ref = ctypes.c_double(0.0)
        fcn.GetParLimits(3, lo_ref, hi_ref)
        lo, hi = lo_ref.value, hi_ref.value

    # tighten limits and refit
    fcn.SetParLimits(3, lo * shrink_factor, hi * shrink_factor)
    graph.Fit(fcn, "WMRQ")                    # quiet, no redraw

    # final check
    if not _is_positive():
        print("WARNING: σ still negative after shrink; "
              "consider excluding this t-bin.")
        return False

    return True

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

        '''
        fff2 = TF2("fff2",
                   "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]",
                   0, 360, 0.0, 1.0)

        fff2 = TF2("fff2",
                   "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]",
                   0, 360, LOEPS-0.1, HIEPS+0.1)
        '''
        # ------------------------------------------------------------------
        # Re-parameterised version enforcing |ρ| ≤ 1 
        # ------------------------------------------------------------------
        fff2 = ROOT.TF2("fff2",
            "[4]*( "                              
            "  [0] + y*[1] "
            "+ sqrt(2*y*(1.+y))*cos(x*0.017453)*[2]*sqrt([0]*[1]) "
            "+ y*cos(2*x*0.017453)*[3]*[0] "
            ")" ,
            0, 360, 0.0, 1.0)
        # ------------------------------------------------------------------
        
        for k in range(4):
            fff2.ReleaseParameter(k)

        # ---------------------------------------------------------------
        par_keys  = ["sigT", "sigL", "rhoLT", "rhoTT", "Norm"]
        current_i = 0

        for idx, key in enumerate(par_keys):
            fff2.SetParName(idx, key)
            if key in PARAM_LIMITS:
                lo, hi = PARAM_LIMITS[key][current_i]
                fff2.SetParLimits(idx, lo, hi)
            else:
                raise KeyError(f"{key} not found in PARAM_LIMITS")

            # --- give MINUIT a sensible first step ---------------------
            if key.startswith("rho"):
                fff2.SetParError(idx, 0.02)          # ρLT, ρTT  → fixed 0.02
            elif key == "Norm":
                fff2.SetParError(idx, 0.05)          # Norm      → 5 % step
            else:
                step = 0.05 * (hi - lo) if hi > lo else 0.1      # σT, σL
                fff2.SetParError(idx, step)
        # ---------------------------------------------------------------

        # SEED the parameters (otherwise they all start at zero)
        fff2.SetParameters( SEED_SIGT,   # σ_T
                            SEED_SIGL,   # σ_L
                            0.0,         # ρ_LT
                            0.0,         # ρ_TT
                            1.0)                # Norm

        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()
        sigLT_change = TGraphErrors()
        sigTT_change = TGraphErrors()

        # ---------------- FIT SEQUENCE ------------------
        fit_step = 0  # counter for adapt_limits

        # --- Fit 1: T ---
        fff2.FixParameter(1, 0.0)   # σL
        fff2.FixParameter(2, 0.0)   # ρLT
        fff2.FixParameter(3, 0.0)   # ρTT
        fff2.FixParameter(4, 1.0)   # Norm
        g_plot_err.Fit(fff2, "WMRQ")       # quiet, no redraw
        check_sigma_positive(fff2, g_plot_err)

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

        # --- Fit 2: L ---
        fff2.ReleaseParameter(1)    # σL now floats
        reset_limits_from_table(fff2, 1, "sigL", stage=1)
        g_plot_err.Fit(fff2, "WMRQ")
        check_sigma_positive(fff2, g_plot_err)

        # ---------- soft floor on σ_L when ε-lever arm is weak -------------
        eps_diff   = abs(HIEPS - LOEPS)
        cond_num   = math.sqrt(1+LOEPS**2)*math.sqrt(1+HIEPS**2) / max(eps_diff, 1e-6)

        if cond_num > COND_MAX:
            # matrix is ill-conditioned → apply soft floor to σ_L
            sigL     = fff2.GetParameter(1)
            sigL_err = fff2.GetParError(1)
            floor    = max(0.25*sigL_err, 1e-3)
            if sigL < floor:
                fff2.SetParameter(1, floor)

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1    

        # --- Fit 3: σ_L , ρ_LT , ρ_TT all float together --------------------------
        fff2.FixParameter(0, fff2.GetParameter(0))   # σT fixed
        fff2.FixParameter(1, fff2.GetParameter(1))   # σL fixed
        stage_idx = 2            # third-pass entry in PARAM_LIMITS

        # --- Grab a crude statistical error scale for this t-bin -------------
        # RMS of the Y-errors in the graph is a quick, stable proxy
        stat_err_estimate = math.sqrt(
            sum((g_plot_err.GetErrorY(i))**2 for i in range(g_plot_err.GetN()))
            / max(1, g_plot_err.GetN())
        )
        # --------------------------------------------------------------------------

        for p_idx, p_key in ((2,"rhoLT"), (3,"rhoTT"), (4, "Norm")):
            fff2.ReleaseParameter(p_idx)
            lo_lim, hi_lim = PARAM_LIMITS[p_key][stage_idx]
            fff2.SetParLimits(p_idx, lo_lim, hi_lim)

            # --- Seeding & step size ----------------------------
            if p_key.startswith("rho"):
                fff2.SetParameter(p_idx, 0.0)
                step = max(0.1, 0.6*stat_err_estimate)
                fff2.SetParError(p_idx, step)
            else:
                fff2.SetParError(p_idx, 0.05*(hi_lim - lo_lim))
            # ---------------------------------------------------------------------

        g_plot_err.Fit(fff2, "WMRQ")
        check_sigma_positive(fff2, g_plot_err)

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1    

        # --- Report reduced χ² ---
        chi2     = fff2.GetChisquare() + penalty(fff2)
        ndf      = max(1, fff2.GetNDF())   # avoid divide-by-zero
        red_chi2 = chi2 / ndf
        print(f"Reduced χ²: {red_chi2:.2f}")
    
        print(f"  **Norm = {fff2.GetParameter(4):.3f}**")

        phi_extrema = np.array([0, 180, 360], float)
        for ang in phi_extrema:
            model_val_lo = fff2.Eval(ang, lo_eps)  # low-ε curve
            model_val_hi = fff2.Eval(ang, hi_eps)  # high-ε curve
            print(f"φ={ang:3.0f}°  low model={model_val_lo:7.2f} nb   "
                f"high model={model_val_hi:7.2f} nb")
        
        # -----------------------  remainder of original code  -----------------------
        # (all canvases, output files, plots, integration, etc. unchanged)
        # ---------------------------------------------------------------------------

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

        glo.Fit(flo, "WMRQ")
        ghi.Fit(fhi, "WMRQ")
        
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

        # ---------------------------------------------------------------
        # Central values -------------------------------------------------
        sig_t   = fff2.GetParameter(0)
        sig_l   = fff2.GetParameter(1)
        rho_lt  = fff2.GetParameter(2)
        rho_tt  = fff2.GetParameter(3)

        sig_lt  = rho_lt * math.sqrt(sig_t * sig_l)
        sig_tt  = rho_tt * sig_t

        # One-sigma errors ----------------------------------------------
        sig_t_err   = fff2.GetParError(0)
        sig_l_err   = fff2.GetParError(1)
        rho_lt_err  = fff2.GetParError(2)
        rho_tt_err  = fff2.GetParError(3)

        # ---------------------------------------------------------------
        # Error propagation (fully guarded) -----------------------------
        _eps = 1e-6                           # numerical floor

        safe_sig_t  = max(abs(sig_t),  _eps)
        safe_sig_l  = max(abs(sig_l),  _eps)
        safe_rho_lt = max(abs(rho_lt), _eps)

        sig_lt_err = abs(sig_lt) * math.sqrt(
            (rho_lt_err / safe_rho_lt)**2
            + (sig_t_err / (2.0 * safe_sig_t))**2
            + (sig_l_err / (2.0 * safe_sig_l))**2
        )

        sig_tt_err = math.hypot( safe_sig_t * rho_tt_err,
                                rho_tt      * sig_t_err )
        # ---------------------------------------------------------------

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
    g_unsep_mult.Fit(f_lin, "WMRQ")
        
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

g_sig_l_total.Fit(f_exp_l, "WMRQ")
g_sig_t_total.Fit(f_exp_t, "WMRQ")

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
g_sig_l_total.Fit(f_exp, "WMRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_t_total.SetMarkerColor(1)
g_sig_t_total.SetLineColor(1)
g_sig_t_total.Draw("A*")
g_sig_t_total.Fit(f_exp, "WMRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_lt_total.Draw("A*")
g_sig_lt_total.Fit(f_exp, "WMRQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_tt_total.Draw("A*")
g_sig_tt_total.Fit(f_exp, "WMRQ")
c_total.Print(outputpdf+')')
c_total.Clear()
