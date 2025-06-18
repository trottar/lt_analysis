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
if polID == "pl":
    mtar = 0.938272046 # GeV/c^2, mass of the target (proton)  
else:
    mtar = 0.9395654133 # GeV/c^2, mass of the target (neutron)
w_set = float(W.replace("p",".")) # W value

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
#    "sigT" : [(0.001, 1e3)]*3,   # σ_T  : transverse
#    "sigL" : [(0.001, 1e3)]*3,   # σ_L  : longitudinal
    "sigT" : [(-1e3, 1e3)]*3,   # σ_T  : transverse
    "sigL" : [(-1e3, 1e3)]*3,   # σ_L  : longitudinal
    "sigLT": [(-50.0, 50.0)]*3,    # σ_LT
    "sigTT": [(-50.0, 50.0)]*3,     # σ_TT
    "rhoLT": [(-1.0, 1.0)]*3,    # ρ_LT : σ_LT / √(σT σL)
    "rhoTT": [(-1.0, 1.0)]*3     # ρ_TT : σ_TT / σT
}
# ------------------------------------------------------------------------------

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

    print("WARNING: negative σ detected → tightening ρTT and ρLT limits and refitting")

    # -----------------------------------------------------------
    # Fetch current limits on parameter 3 (ρTT)  –  two ways
    # -----------------------------------------------------------
    # tighten both ρ_LT (idx=2) and ρ_TT (idx=3) and refit
    for idx in (2, 3):
        try:
            lo_i, hi_i = fcn.GetParLimits(idx)
        except TypeError:
            lo_ref = ctypes.c_double(0.0)
            hi_ref = ctypes.c_double(0.0)
            fcn.GetParLimits(idx, lo_ref, hi_ref)
            lo_i, hi_i = lo_ref.value, hi_ref.value
        fcn.SetParLimits(idx, lo_i * shrink_factor, hi_i * shrink_factor)
    graph.Fit(fcn, "MRQ")                    # quiet, no redraw

    # final check
    if not _is_positive():
        print("WARNING: σ still negative after shrink; "
              "consider excluding this t-bin.")
        return False

    return True

###############################################################################################################################################

# Import separated xsects models
from lt_active import LT_sep_x_fun_wrapper, LT_sep_x_fun_unsep_wrapper

###############################################################################################################################################

def single_setting(q2_set, w_set, fn_lo, fn_hi):

    w_set_num = float(w_set.replace("p",".")) # W value

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

        xsect_scalefac = 1.0 # scale factor

        # ——— Low graph ———
        print("Drawing low data into internal arrays…")
        nlo_count = nlo.Draw("phi:x_mod:dx", tcut, "goff")
        nlo_rows  = nlo.GetSelectedRows()

        # Build the low TGraphErrors
        glo_tmp = ROOT.TGraphErrors()
        for idx in range(nlo_rows):
            phi_val  = nlo.GetV1()[idx]   # φ-values            
            x_val    = nlo.GetV2()[idx]   # x-values
            dx_error = nlo.GetV3()[idx]   # dx → use as y-error

            print(f"  [LO] Point {idx}: φ={phi_val:.4f}, x={x_val:.4f}, dx_error={dx_error:.4f}")
            glo_tmp.SetPoint(idx, phi_val, x_val)
            glo_tmp.SetPointError(idx, 0.0, dx_error)

        LT_sep_x_lo_fun = LT_sep_x_fun_wrapper(lo_eps)
        flo = TF1("lo_eps_fit", LT_sep_x_lo_fun, 0, 360, 4)
        LT_sep_x_lo_fun_unsep = LT_sep_x_fun_unsep_wrapper(lo_eps)
        flo_unsep = TF1("lo_eps_unsep", LT_sep_x_lo_fun_unsep, 0, 2*PI, 4)
        
        glo = glo_tmp.Clone("glo")
        ave_sig_lo = glo.GetMean(2)
        err_sig_lo = glo.GetRMS(2)/math.sqrt(glo.GetN())
        print(f"Epsilon: {lo_eps:.4f} (low)")
        print(f"Average low σ: {ave_sig_lo:.4f} ± {err_sig_lo:.4f}")        

        print("#"*25)

        sig_lo.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_lo.GetYaxis().SetTitle("#bar{#it{#sigma}_{Low}}")
        sig_lo.SetTitle("t = {:.3f}".format(t_list[i]))
        sig_lo.SetPoint(sig_lo.GetN(), float(t_list[i]), ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)
        
        # ——— High graph ———
        print("Drawing high data into internal arrays…")
        nhi_count = nhi.Draw("phi:x:dx", tcut, "goff")
        nhi_rows  = nhi.GetSelectedRows()

        # Build the high TGraphErrors
        ghi_tmp = ROOT.TGraphErrors()
        for idx in range(nhi_rows):
            phi_val  = nhi.GetV1()[idx]   # φ-values            
            x_val    = nhi.GetV2()[idx]   # x-values
            dx_error = nhi.GetV3()[idx]   # dx → use as y-error

            print(f"  [HI] Point {idx}: φ={phi_val:.4f}, x={x_val:.4f}, dx_error={dx_error:.4f}")
            ghi_tmp.SetPoint(idx, phi_val, x_val)
            ghi_tmp.SetPointError(idx, 0.0, dx_error)

        LT_sep_x_hi_fun = LT_sep_x_fun_wrapper(hi_eps)
        fhi = TF1("hi_eps_fit", LT_sep_x_hi_fun, 0, 360, 4)
        LT_sep_x_hi_fun_unsep = LT_sep_x_fun_unsep_wrapper(hi_eps)
        fhi_unsep = TF1("hi_eps_unsep", LT_sep_x_hi_fun_unsep, 0, 2*PI, 4)
            
        ghi = ghi_tmp.Clone("ghi")
        ave_sig_hi = ghi.GetMean(2)
        err_sig_hi = ghi.GetRMS(2)/math.sqrt(ghi.GetN())
        print(f"Epsilon: {hi_eps:.4f} (high)")
        print(f"Average high σ: {ave_sig_hi:.4f} ± {err_sig_hi:.4f}")        

        print("#"*25)   

        sig_hi.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_hi.GetYaxis().SetTitle("#bar{#it{#sigma}_{High}}")
        sig_hi.SetTitle("t = {:.3f}".format(t_list[i]))
        sig_hi.SetPoint(sig_hi.GetN(), float(t_list[i]), ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)
        
        g_plot_err = TGraph2DErrors()
        g_xx, g_yy, g_yy_err = ctypes.c_double(0),ctypes.c_double(0),ctypes.c_double(0)

        # Fractional point-to-point systematic uncertainty (same in both loops)
        syst_frac = pt_to_pt_systematic_error / 100.0

        # Loop over low-ε points (mapping: φ→x, ε_lo→y-axis coordinate for TF2)
        for ii in range(glo.GetN()):
            # Read φ into phi_val, σ into sigma_val
            glo.GetPoint(ii, g_xx, g_yy)
            phi_val   = g_xx.value
            sigma_val = g_yy.value

            # Statistical + systematic uncertainty on σ
            stat_err   = glo.GetErrorY(ii)
            syst_err   = syst_frac * sigma_val
            total_err  = math.sqrt(stat_err**2 + syst_err**2)

            # Accumulate inverse-variance weight
            lo_cross_sec_err[i] += 1.0 / (total_err**2)

            # Insert into TGraph2DErrors: x=φ, y=ε_lo, z=σ
            g_plot_err.SetPoint(g_plot_err.GetN(),
                                phi_val,
                                lo_eps,
                                sigma_val)
            g_plot_err.SetPointError(g_plot_err.GetN() - 1,
                                     0.0,    # no φ error
                                     0.0,    # no ε error
                                     total_err)
            
        # Loop over high-ε points (mapping: φ→x, ε_hi→y-axis coordinate for TF2)
        for ii in range(ghi.GetN()):
            ghi.GetPoint(ii, g_xx, g_yy)
            phi_val   = g_xx.value
            sigma_val = g_yy.value

            stat_err  = ghi.GetErrorY(ii)
            syst_err  = syst_frac * sigma_val
            total_err = math.sqrt(stat_err**2 + syst_err**2)

            hi_cross_sec_err[i] += 1.0 / (total_err**2)

            g_plot_err.SetPoint(g_plot_err.GetN(),
                                phi_val,
                                hi_eps,
                                sigma_val)
            g_plot_err.SetPointError(g_plot_err.GetN() - 1,
                                     0.0,
                                     0.0,
                                     total_err)

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

        # ------------------------------------------------------------------
        # Re-parameterised version enforcing |ρ| ≤ 1 
        # ------------------------------------------------------------------

        fff2_normfactor = 1.0 # scale factor for the fit function

        w_dep = 1/((w_list[i]**2) - (mtar**2))**(0.85*(w_set_num**2) - 5.97*w_set_num + 12.68)
        fff2_normfactor_wdep =  (1/w_dep) # change W dependence
        fff2_normfactor_qdep = 1e-1 * np.exp(-q2_list[i])
        #fff2_normfactor = fff2_normfactor_qdep

        a = 1.0
        b = 1.0
        c = 1.0
        d = 1.0

        # Re-parameterised LT/TT enforcing |ρ|≤1:
        # fff2_normfactor, a, b, c, d, PI must be in scope here
        '''
        fff2 = TF2(
            "fff2",
            (
                f"{fff2_normfactor} * ("
                f"{a} * [0]"                                             # σ_T
                f"+ {b} * y*[1]"                                        # ε·σ_L
                f"+ {c} * sqrt(2*y*(1.+y))*cos(x*({PI}/180))*[2]*sqrt([0]*[1])"  # ρ_LT·√(σₜσₗ)
                f"+ {d} * y*cos(2*x*({PI}/180))*[3]*[0]"               # ρ_TT·σₜ
                f")"
            ),
            0, 360,
            LOEPS-0.1, HIEPS+0.1
        )      
        '''
        fff2 = TF2(
            "fff2",
            (
                f"{fff2_normfactor} * ("
                f"{a} * [0]"                                             # σ_T
                f"+ {b} * y*[1]"                                        # ε·σ_L
                f"+ {c} * sqrt(2*y*(1.+y))*cos(x*({PI}/180))*[2]"  # ρ_LT·√(σₜσₗ)
                f"+ {d} * y*cos(2*x*({PI}/180))*[3]"               # ρ_TT·σₜ
                f")"
            ),
            0, 360,
            LOEPS-0.1, HIEPS+0.1
        )             

        for k in range(4):
            fff2.ReleaseParameter(k)

        # ---------------------------------------------------------------
        #par_keys  = ["sigT", "sigL", "rhoLT", "rhoTT"]
        par_keys  = ["sigT", "sigL", "sigLT", "sigTT"]
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
                fff2.SetParError(idx, 0.02)            # ±0.02 for ρ’s
            else:
                step = 0.05 * (hi - lo) if hi > lo else 0.1
                fff2.SetParError(idx, step)            # 5 % of range for σT, σL
        # ---------------------------------------------------------------

        # — Dynamic seeds based on data averages —
        # Equations 1.13 and 1.14 from Bill's thesis
        SEED_SIGT = ((HIEPS * ave_sig_lo) - (LOEPS * ave_sig_hi)) / eps_diff
        SEED_SIGL = (ave_sig_hi - ave_sig_lo) / eps_diff
        print(f"SEED_SIGT = {SEED_SIGT}")
        print(f"SEED_SIGL = {SEED_SIGL}")
        fff2.SetParameters(
            SEED_SIGT,      # σ_T
            SEED_SIGL,      # σ_L
            0.0,         # ρ_LT
            0.0          # ρ_TT
        )
        
        # — Give Minuit a finite “kick size” on each parameter —
        fff2.SetParError(0, max(1.0, 0.1 * SEED_SIGT))     # σ_T step ≃10% of its seed (but at least 1)
        fff2.SetParError(1, max(0.1, 0.1 * abs(SEED_SIGL)))# σ_L step
        fff2.SetParError(2, 0.5)                        # ρ_LT step
        fff2.SetParError(3, 0.5)                        # ρ_TT step        

        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()
        sigLT_change = TGraphErrors()
        sigTT_change = TGraphErrors()

        # ---------------- FIT SEQUENCE ------------------
        fit_step = 0  # counter for adapt_limits

        # --- Fit 1: T ---
        fff2.FixParameter(1, SEED_SIGL)   # σL
        fff2.FixParameter(2, 0.0)   # ρLT
        fff2.FixParameter(3, 0.0)   # ρTT
        # — Apply limits for all parameters in stage 0 —
        #for idx, name in enumerate(["sigT","sigL","rhoLT","rhoTT"]):
        for idx, name in enumerate(["sigT","sigL","sigLT","sigTT"]):
            reset_limits_from_table(fff2, idx, name, stage=0)
        # — Give Minuit a finite “kick size” on each parameter —
        fff2.SetParError(0, max(1.0, 0.1 * SEED_SIGT))     # σ_T step ≃10% of its seed (but at least 1)
        fff2.SetParError(1, max(0.1, 0.1 * abs(SEED_SIGL)))# σ_L step
        fff2.SetParError(2, 0.5)                        # ρ_LT step
        fff2.SetParError(3, 0.5)                        # ρ_TT step  
        g_plot_err.Fit(fff2, "SEWQ")       # quiet, no redraw
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

        # --- Fit 2: L (fix T) ---
        fff2.FixParameter(0, fff2.GetParameter(0))  # σT now fixed
        fff2.ReleaseParameter(1)    # σL now floats
        # — Apply limits for all parameters in stage 1 —
        #for idx, name in enumerate(["sigT","sigL","rhoLT","rhoTT"]):
        for idx, name in enumerate(["sigT","sigL","sigLT","sigTT"]):
            reset_limits_from_table(fff2, idx, name, stage=1)
        # — Give Minuit a finite “kick size” on each parameter —
        fff2.SetParError(0, max(1.0, 0.1 * SEED_SIGT))     # σ_T step ≃10% of its seed (but at least 1)
        fff2.SetParError(1, max(0.1, 0.1 * abs(SEED_SIGL)))# σ_L step
        fff2.SetParError(2, 0.5)                        # ρ_LT step
        fff2.SetParError(3, 0.5)                        # ρ_TT step         
        g_plot_err.Fit(fff2, "SEWQ")
        check_sigma_positive(fff2, g_plot_err)

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1    

        # --- Fit 3: ρ_LT , ρ_TT --------------------------
        fff2.FixParameter(0, fff2.GetParameter(0))  # σT now fixed
        fff2.FixParameter(1, fff2.GetParameter(1))  # σL now fixed
        fff2.ReleaseParameter(2)    # ρ_LT now floats
        fff2.ReleaseParameter(3)    # ρ_TT now floats
        # — Apply limits for all parameters in stage 2 —
        #for idx, name in enumerate(["sigT","sigL","rhoLT","rhoTT"]):
        for idx, name in enumerate(["sigT","sigL","sigLT","sigTT"]):
            reset_limits_from_table(fff2, idx, name, stage=2)
        # — Give Minuit a finite “kick size” on each parameter —
        fff2.SetParError(0, max(1.0, 0.1 * SEED_SIGT))     # σ_T step ≃10% of its seed (but at least 1)
        fff2.SetParError(1, max(0.1, 0.1 * abs(SEED_SIGL)))# σ_L step
        fff2.SetParError(2, 0.5)                        # ρ_LT step
        fff2.SetParError(3, 0.5)                        # ρ_TT step         
        g_plot_err.Fit(fff2, "SEWQ")
        check_sigma_positive(fff2, g_plot_err)

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1         

        # --- Fit 4: ALL --------------------------
        fff2.ReleaseParameter(0)    # σL now floats
        fff2.ReleaseParameter(1)    # σL now floats
        # — Apply limits for all parameters in stage 2 —
        #for idx, name in enumerate(["sigT","sigL","rhoLT","rhoTT"]):
        for idx, name in enumerate(["sigT","sigL","sigLT","sigTT"]):
            reset_limits_from_table(fff2, idx, name, stage=2)
        # — Give Minuit a finite “kick size” on each parameter —
        fff2.SetParError(0, max(1.0, 0.1 * SEED_SIGT))     # σ_T step ≃10% of its seed (but at least 1)
        fff2.SetParError(1, max(0.1, 0.1 * abs(SEED_SIGL)))# σ_L step
        fff2.SetParError(2, 0.5)                        # ρ_LT step
        fff2.SetParError(3, 0.5)                        # ρ_TT step         
        g_plot_err.Fit(fff2, "SEWQ")
        check_sigma_positive(fff2, g_plot_err)     

        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN()+1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN()-1, 0, fff2.GetParError(1))
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN()+1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN()-1, 0, fff2.GetParError(0))

        fit_step += 1    

        # --- Report reduced χ² ---
        chi2     = fff2.GetChisquare()
        ndf      = max(1, fff2.GetNDF())   # avoid divide-by-zero
        red_chi2 = chi2 / ndf
        
        # -----------------------  remainder of original code  -----------------------
        # (all canvases, output files, plots, integration, etc. unchanged)
        # ---------------------------------------------------------------------------

        c1 =  TCanvas()

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

        glo.Fit(flo, "SEWQ")
        ghi.Fit(fhi, "SEWQ")
        
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

        '''
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
        '''
        sig_t   = fff2.GetParameter(0)
        sig_l   = fff2.GetParameter(1)        
        sig_lt  = fff2.GetParameter(2)
        sig_tt  = fff2.GetParameter(3)
        sig_t_err   = fff2.GetParError(0)
        sig_l_err   = fff2.GetParError(1)        
        sig_lt_err  = fff2.GetParError(2)
        sig_tt_err  = fff2.GetParError(3)
 
        print(f"\n=== Bin {i+1} Summary ===")
        print(f"  t = {t_list[i]:.3f} GeV²   θ = {theta_list[i]:.1f}°   W = {w_list[i]:.3f} GeV   Q² = {q2_list[i]:.3f} GeV²")
        print(f"  ε_lo = {lo_eps_list[i]:.3f}   ε_hi = {hi_eps_list[i]:.3f}\n")
        #print(f"  ρ_LT = {rho_lt:.3f} ± {rho_lt_err:.3f}")
        #print(f"  ρ_TT = {rho_tt:.3f} ± {rho_tt_err:.3f}")
        print(f"  Reduced χ² = {red_chi2:.2f} (NDF = {ndf})")        
        print(f"  σ_T  = {sig_t:.3f} ± {sig_t_err:.3f}")
        print(f"  σ_L  = {sig_l:.3f} ± {sig_l_err:.3f}")
        print(f"  σ_LT = {sig_lt:.3f} ± {sig_lt_err:.3f}")
        print(f"  σ_TT = {sig_tt:.3f} ± {sig_tt_err:.3f}")
        print("=== End of Bin Summary ===\n")

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
    g_unsep_mult.Fit(f_lin, "SEWQ")
        
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

g_sig_l_total.Fit(f_exp_l, "SEWQ")
g_sig_t_total.Fit(f_exp_t, "SEWQ")

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
g_sig_l_total.Fit(f_exp, "SEWQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_t_total.SetMarkerColor(1)
g_sig_t_total.SetLineColor(1)
g_sig_t_total.Draw("A*")
g_sig_t_total.Fit(f_exp, "SEWQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_lt_total.Draw("A*")
g_sig_lt_total.Fit(f_exp, "SEWQ")
c_total.Print(outputpdf)
c_total.Clear()

g_sig_tt_total.Draw("A*")
g_sig_tt_total.Fit(f_exp, "SEWQ")
c_total.Print(outputpdf+')')
c_total.Clear()
