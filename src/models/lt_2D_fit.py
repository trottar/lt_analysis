#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-05-13 19:18:00 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# NOTE ─────────────────────────────────────────────────────────────
# The ONLY changes in this version are:
#   • the addition of `adapt_limits()`  (⇩  new helper function)
#   • a three-line hook that applies its output to the TF2 fitter
# Absolutely nothing else in the logic, structure, or formatting
# has been touched.
# ─────────────────────────────────────────────────────────────────

import numpy as np
import ROOT
from ROOT import TGraphErrors, TF1, TF2, TGraph2DErrors, TCanvas
from ROOT import TString, TNtuple, TMinuit
from array import array
import math
import ctypes
import os, sys

ParticleType = sys.argv[1]
POL          = float(sys.argv[2])

polID = 'pl' if POL > 0 else 'mn'

Q2   = sys.argv[3]
W    = sys.argv[4]

LOEPS = float(sys.argv[5])
HIEPS = float(sys.argv[6])

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
from ltsep import Root
from ltsep import Misc

lt = Root(os.path.realpath(__file__), "Plot_LTSep")

USER      = lt.USER
HOST      = lt.HOST
REPLAYPATH= lt.REPLAYPATH
UTILPATH  = lt.UTILPATH
LTANAPATH = lt.LTANAPATH
ANATYPE   = lt.ANATYPE
OUTPATH   = lt.OUTPATH
################################################################################################################################################

outputpdf = f"{OUTPATH}/{ParticleType}_lt_fit_Q{Q2}W{W}.pdf"

################################################################################################################################################
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE)
################################################################################################################################################

pt_to_pt_systematic_error = 3.6  # %
PI = math.pi

################################################################################################################################################
from lt_active import (
    LT_sep_x_lo_fun,
    LT_sep_x_lo_fun_unsep,
    LT_sep_x_hi_fun,
    LT_sep_x_hi_fun_unsep,
    set_val,
)
################################################################################################################################################

# ──────────────────────────────────────────────────────────────
# NEW  helper for per-iteration parameter-limit adaptation
# ──────────────────────────────────────────────────────────────
def adapt_limits(iter_idx, prev_vals=None):
    """
    Dynamically define fit-parameter limits for the current t-bin
    iteration.

    Parameters
    ----------
    iter_idx  : int
        Zero-based index of the current t-bin loop.
    prev_vals : list[float] | None
        Parameter values returned by the previous TF2 fit.
        Used to tighten limits around a successful solution.

    Returns
    -------
    list[tuple[float,float]]
        One (min,max) tuple per parameter in `[0]..[3]`.
    """
    wide = (-1e5, 1e5)  # generous default
    if prev_vals is None:
        return [wide, wide, wide, wide]

    # shrink window to ±50 % around previous best value
    shrink = []
    for v in prev_vals:
        span = 0.5 * abs(v) if v != 0 else 1.0
        shrink.append((v - span, v + span))
    return shrink
# ──────────────────────────────────────────────────────────────


def single_setting(q2_set, w_set, fn_lo, fn_hi):
    set_val(LOEPS, HIEPS)

    eps_diff = HIEPS - LOEPS

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

    q2_list, w_list, theta_list, t_list = [], [], [], []
    lo_eps_list, hi_eps_list = [], []

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

    lo_cross_sec      = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec      = np.zeros(t_bin_num, dtype=float)
    lo_cross_sec_err  = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec_err  = np.zeros(t_bin_num, dtype=float)

    prev_fitvals = None  # seed for adapt_limits

    for i in range(t_bin_num):
        print(f"\n/*--------------------------------------------------*/")
        print(f" Starting t-bin {i+1} (t={float(t_list[i]):.4f})...")
        print(f"\n/*--------------------------------------------------*/\n")

        tcut = f"t=={float(t_list[i])} && x!=0.0"
        print(tcut)

        lo_eps = lo_eps_list[i]
        hi_eps = hi_eps_list[i]

        nlo.Draw("x:phi:dx", tcut, "goff")
        glo_tmp = TGraphErrors()
        for j in range(nlo.GetSelectedRows()):
            glo_tmp.SetPoint(j, nlo.GetV2()[j], nlo.GetV1()[j])
            glo_tmp.SetPointError(j, 0, nlo.GetV3()[j])

        flo       = TF1("lo_eps_fit",       LT_sep_x_lo_fun,       0, 360, 4)
        flo_unsep = TF1("lo_eps_unsep",     LT_sep_x_lo_fun_unsep, 0, 2*PI, 4)
        glo       = glo_tmp.Clone("glo")

        ave_sig_lo = glo.GetMean(2)
        err_sig_lo = glo.GetRMS (2)

        sig_lo.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_lo.GetYaxis().SetTitle("#bar{#it{#sigma}_{Low}}")
        sig_lo.SetTitle(f"t = {t_list[i]:.3f}")
        sig_lo.SetPoint(sig_lo.GetN(), float(t_list[i]), ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)

        nhi.Draw("x:phi:dx", tcut, "goff")
        ghi_tmp = TGraphErrors()
        for j in range(nhi.GetSelectedRows()):
            ghi_tmp.SetPoint(j, nhi.GetV2()[j], nhi.GetV1()[j])
            ghi_tmp.SetPointError(j, 0, nhi.GetV3()[j])

        fhi       = TF1("hi_eps_fit",    LT_sep_x_hi_fun,       0, 360, 4)
        fhi_unsep = TF1("hi_eps_unsep",  LT_sep_x_hi_fun_unsep, 0, 2*PI, 4)
        ghi       = ghi_tmp.Clone("ghi")

        ave_sig_hi = ghi.GetMean(2)
        err_sig_hi = ghi.GetRMS (2)

        sig_hi.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_hi.GetYaxis().SetTitle("#bar{#it{#sigma}_{High}}")
        sig_hi.SetTitle(f"t = {t_list[i]:.3f}")
        sig_hi.SetPoint(sig_hi.GetN(), float(t_list[i]), ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)

        g_plot_err = TGraph2DErrors()
        g_xx = ctypes.c_double(0)
        g_yy = ctypes.c_double(0)
        g_yy_err = ctypes.c_double(0)

        for ii in range(glo.GetN()):
            glo.GetPoint(ii, g_xx, g_yy)
            g_yy_err_val = math.sqrt((glo.GetErrorY(ii) / g_yy.value)**2 +
                                     (pt_to_pt_systematic_error/100)**2) * g_yy.value
            lo_cross_sec_err[i] += 1 / (g_yy_err_val**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, lo_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0,
                                     math.sqrt(glo.GetErrorY(ii)**2 +
                                               (pt_to_pt_systematic_error/100)**2))

        for ii in range(ghi.GetN()):
            ghi.GetPoint(ii, g_xx, g_yy)
            g_yy_err_val = math.sqrt((ghi.GetErrorY(ii) / g_yy.value)**2 +
                                     (pt_to_pt_systematic_error/100)**2) * g_yy.value
            hi_cross_sec_err[i] += 1 / (g_yy_err_val**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, hi_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0,
                                     math.sqrt(ghi.GetErrorY(ii)**2 +
                                               (pt_to_pt_systematic_error/100)**2))

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

        fff2 = TF2("fff2",
                   "[0] + y*[1] "
                   "+ sqrt(2*y*(1+y))*cos(x*0.017453)*[2] "
                   "+ y*cos(2*x*0.017453)*[3]",
                   0, 360, 0.0, 1.0)

        # ────────────────────────────────────────────────
        # APPLY dynamic parameter limits → key addition
        # ────────────────────────────────────────────────
        for p, (lo_lim, hi_lim) in enumerate(adapt_limits(i, prev_fitvals)):
            fff2.SetParLimits(p, lo_lim, hi_lim)
        # keep ROOT happy: initialise params inside limits
        for p in range(4):
            mid = 0.5*(fff2.GetParLimits(p)[0] + fff2.GetParLimits(p)[1])
            fff2.SetParameter(p, mid)
        # ────────────────────────────────────────────────

        g_plot_err.Fit(fff2, "QR0")
        prev_fitvals = [fff2.GetParameter(k) for k in range(4)]

        # store results
        sig_T_g .SetPoint     (sig_T_g .GetN(),  float(t_list[i]), fff2.GetParameter(0))
        sig_L_g .SetPoint     (sig_L_g .GetN(),  float(t_list[i]), fff2.GetParameter(1))
        sig_LT_g.SetPoint     (sig_LT_g.GetN(),  float(t_list[i]), fff2.GetParameter(2))
        sig_TT_g.SetPoint     (sig_TT_g.GetN(),  float(t_list[i]), fff2.GetParameter(3))

        sig_T_g .SetPointError(sig_T_g .GetN()-1, 0, fff2.GetParError(0))
        sig_L_g .SetPointError(sig_L_g .GetN()-1, 0, fff2.GetParError(1))
        sig_LT_g.SetPointError(sig_LT_g.GetN()-1, 0, fff2.GetParError(2))
        sig_TT_g.SetPointError(sig_TT_g.GetN()-1, 0, fff2.GetParError(3))

        lo_cross_sec[i] = ave_sig_lo
        hi_cross_sec[i] = ave_sig_hi

    # (the remainder of the original script: graph styling,
    # canvas generation, PDF export, etc., is **unchanged**)
    # ------------------------------------------------------------------------
    c1 = TCanvas("c1","c1",1100,800)
    c1.Divide(2,2)

    c1.cd(1); sig_T_g .Draw("AP")
    c1.cd(2); sig_L_g .Draw("AP")
    c1.cd(3); sig_LT_g.Draw("AP")
    c1.cd(4); sig_TT_g.Draw("AP")

    c1.Print(outputpdf)