#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-05-07 02:50:07 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
# ------------------------------------------------------------------
# standard libs
# ------------------------------------------------------------------
import numpy as np
import ROOT
from ROOT import (TGraphErrors, TF1, TF2, TGraph2DErrors, TCanvas, TString,
                  TNtuple, TMinuit)
from array import array
import math, ctypes, os, sys

# ------------------------------------------------------------------
# cmd-line args
# ------------------------------------------------------------------
ParticleType = sys.argv[1]
POL          = float(sys.argv[2])
polID        = 'pl' if POL > 0 else 'mn'
Q2, W        = sys.argv[3], sys.argv[4]
LOEPS        = float(sys.argv[5])
HIEPS        = float(sys.argv[6])

# ------------------------------------------------------------------
# ltsep package import and pathing definitions
# ------------------------------------------------------------------
from ltsep import Root, Misc
lt = Root(os.path.realpath(__file__), "Plot_LTSep")
USER, HOST     = lt.USER, lt.HOST
REPLAYPATH     = lt.REPLAYPATH
UTILPATH       = lt.UTILPATH
LTANAPATH      = lt.LTANAPATH
ANATYPE        = lt.ANATYPE
OUTPATH        = lt.OUTPATH

outputpdf = f"{OUTPATH}/{ParticleType}_lt_fit_Q{Q2}W{W}.pdf"

# ------------------------------------------------------------------
# Importing utility functions
# ------------------------------------------------------------------

sys.path.append("../utility")
from utility import adapt_limits

# ------------------------------------------------------------------
# ROOT batch mode & constants
# ------------------------------------------------------------------
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE)

PI = math.pi

##############
# HARD CODED #
##############

pt_to_pt_systematic_error = 3.6   # %

# ------------------------------------------------------------------
# absolute physics limits (never crossed)
# index: 0 = σT  1 = σL  2 = σLT  3 = σTT
# ------------------------------------------------------------------
HARD_LO = [-5.0, -5.0,  -10.0, -10.0]
HARD_HI = [1.0e3, 1.0e3,  50.0,  50.0]
##############
##############
##############

# ------------------------------------------------------------------
# compact fit plan (what is free each step)
# ------------------------------------------------------------------
fit_plan = [
    {"free": [0, 1]},            # 1: σT σL   (LT,TT fixed)
    {"free": [2]},               # 2: σLT
    {"free": [0, 1]},            # 3: refine σT σL
    {"free": [3]},               # 4: σTT
    {"free": [0, 1]},            # 5: refine σT σL again
    {"free": [2, 3]},            # 6: refine σLT σTT
    {"free": [0, 1, 2, 3]},      # 7: all free
]

# ------------------------------------------------------------------
# import separated-xsec models
# ------------------------------------------------------------------
from lt_active import (LT_sep_x_lo_fun, LT_sep_x_lo_fun_unsep,
                       LT_sep_x_hi_fun, LT_sep_x_hi_fun_unsep, set_val)

# ------------------------------------------------------------------
# global graphs
# ------------------------------------------------------------------
g_sig_l_total  = TGraphErrors()
g_sig_t_total  = TGraphErrors()
g_sig_lt_total = TGraphErrors()
g_sig_tt_total = TGraphErrors()
g_unsep_lo     = TGraphErrors()
g_unsep_hi     = TGraphErrors()

# ------------------------------------------------------------------
# main per-setting routine
# ------------------------------------------------------------------
def single_setting(q2_set, w_set, fn_lo, fn_hi):

    set_val(LOEPS, HIEPS)
    eps_diff = HIEPS - LOEPS

    # overall result graphs
    sig_L_g  = TGraphErrors()
    sig_T_g  = TGraphErrors()
    sig_LT_g = TGraphErrors()
    sig_TT_g = TGraphErrors()
    sig_lo   = TGraphErrors()
    sig_hi   = TGraphErrors()
    sig_diff_g = TGraphErrors()

    # load tuples --------------------------------------------------
    nlo = TNtuple("nlo", "nlo", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nlo.ReadFile(fn_lo)
    nhi = TNtuple("nhi", "nhi", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nhi.ReadFile(fn_hi)

    # collect unique t-bins ---------------------------------------
    q2_list, w_list, theta_list, t_list = [], [], [], []
    lo_eps_list, hi_eps_list           = [], []

    for evt in nlo:
        if evt.t not in t_list:
            q2_list.append(evt.Q2);  w_list.append(evt.w)
            theta_list.append(evt.theta);  t_list.append(evt.t)
            lo_eps_list.append(evt.eps)

    for evt in nhi:
        if evt.t not in hi_eps_list:
            hi_eps_list.append(evt.eps)

    t_bin_num = len(t_list)

    # integrated unsep σ
    lo_cross_sec  = np.zeros(t_bin_num)
    hi_cross_sec  = np.zeros(t_bin_num)
    lo_cross_err  = np.zeros(t_bin_num)
    hi_cross_err  = np.zeros(t_bin_num)

    # --------------------------------------------------------------
    for i in range(t_bin_num):

        print(f"\n---- t-bin {i+1}/{t_bin_num}   (-t = {t_list[i]:.4f}) ----")

        tcut   = f"t=={t_list[i]} && x!=0.0"
        lo_eps = lo_eps_list[i]
        hi_eps = hi_eps_list[i]

        # low-ε points -------------------------------------------
        nlo.Draw("x:phi:dx", tcut, "goff")
        glo = TGraphErrors(nlo.GetSelectedRows(),
                           nlo.GetV2(), nlo.GetV1(),
                           np.zeros(nlo.GetSelectedRows()),
                           nlo.GetV3())
        ave_sig_lo = glo.GetMean(2);   err_sig_lo = glo.GetRMS(2)
        sig_lo.SetPoint(sig_lo.GetN(), t_list[i], ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)

        # high-ε points ------------------------------------------
        nhi.Draw("x:phi:dx", tcut, "goff")
        ghi = TGraphErrors(nhi.GetSelectedRows(),
                           nhi.GetV2(), nhi.GetV1(),
                           np.zeros(nhi.GetSelectedRows()),
                           nhi.GetV3())
        ave_sig_hi = ghi.GetMean(2);   err_sig_hi = ghi.GetRMS(2)
        sig_hi.SetPoint(sig_hi.GetN(), t_list[i], ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)

        # build 2-D graph (φ,ε,σ) with total errors --------------
        g_plot_err = TGraph2DErrors()
        for g,eps,err_acc in [(glo,lo_eps,lo_cross_err),
                              (ghi,hi_eps,hi_cross_err)]:
            for ii in range(g.GetN()):
                x, y = ctypes.c_double(), ctypes.c_double()
                g.GetPoint(ii, x, y)
                stat = g.GetErrorY(ii)
                tot  = math.sqrt(stat**2 + (pt_to_pt_systematic_error/100*y.value)**2)
                g_plot_err.SetPoint(g_plot_err.GetN(), x, eps, y)
                g_plot_err.SetPointError(g_plot_err.GetN()-1, 0, 0, tot)
                err_acc[i] += 1/(tot**2)            # accumulate 1/σ²

        lo_cross_err[i] = 1/math.sqrt(lo_cross_err[i])
        hi_cross_err[i] = 1/math.sqrt(hi_cross_err[i])

        # --- TF2 with adaptive limits ---------------------------
        fff2 = TF2("fff2",
                   "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] "
                   "+ y*cos(2*x*0.017453)*[3]",
                   0, 360, 0.0, 1.0)
        fff2.SetParameters(ave_sig_hi, ave_sig_lo, 0.0, 0.0)
        for p,(lo,hi) in enumerate(zip(HARD_LO,HARD_HI)):
            fff2.SetParLimits(p, lo, hi)

        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()

        for step in fit_plan:
            for p in range(4):
                if p in step["free"]:
                    fff2.ReleaseParameter(p)
                else:
                    fff2.FixParameter(p, fff2.GetParameter(p))

            for p in step["free"]:
                if fff2.GetParError(p) == 0:
                    fff2.SetParameter(p, 0.0 if p>1 else ave_sig_lo)

            adapt_limits(fff2, nsig=4,
                         hard_lo=HARD_LO, hard_hi=HARD_HI)
            g_plot_err.Fit(fff2, "MRQ")       # silent

            idx = sigL_change.GetN()
            sigL_change.SetPoint(idx, idx+1, fff2.GetParameter(1))
            sigL_change.SetPointError(idx, 0, fff2.GetParError(1))
            sigT_change.SetPoint(idx, idx+1, fff2.GetParameter(0))
            sigT_change.SetPointError(idx, 0, fff2.GetParError(0))

        # ------------ post-fit handling (unchanged) -------------
        flo        = TF1("flo", LT_sep_x_lo_fun, 0, 360, 4)
        flo_unsep  = TF1("flo_u", LT_sep_x_lo_fun_unsep, 0, 2*PI, 4)
        fhi        = TF1("fhi", LT_sep_x_hi_fun, 0, 360, 4)
        fhi_unsep  = TF1("fhi_u", LT_sep_x_hi_fun_unsep, 0, 2*PI, 4)

        for f in (flo,fhi,flo_unsep,fhi_unsep):
            f.FixParameter(0, fff2.GetParameter(0))   # σT
            f.FixParameter(1, fff2.GetParameter(1))   # σL
            f.FixParameter(2, fff2.GetParameter(2))   # σLT
            f.FixParameter(3, fff2.GetParameter(3))   # σTT

        glo.Fit(flo, "MRQ")
        ghi.Fit(fhi, "MRQ")

        # integrate un-separated σ
        lo_cross_sec[i] = flo_unsep.Integral(0, 2*PI) / (2*PI)
        hi_cross_sec[i] = fhi_unsep.Integral(0, 2*PI) / (2*PI)

        # averaged σ difference for slope plot
        try:
            sig_diff = (ave_sig_hi - ave_sig_lo) / eps_diff
            sig_diff_err = sig_diff*math.sqrt(
                (err_sig_hi/ave_sig_hi)**2 + (err_sig_lo/ave_sig_lo)**2)
        except ZeroDivisionError:
            sig_diff, sig_diff_err = 0, 0
        sig_diff_g.SetPoint(sig_diff_g.GetN(), t_list[i], sig_diff)
        sig_diff_g.SetPointError(sig_diff_g.GetN()-1, 0, sig_diff_err)

        # add to separated-σ vs t graphs
        sig_L_g .SetPoint(i, t_list[i], fff2.GetParameter(1))
        sig_T_g .SetPoint(i, t_list[i], fff2.GetParameter(0))
        sig_LT_g.SetPoint(i, t_list[i], fff2.GetParameter(2))
        sig_TT_g.SetPoint(i, t_list[i], fff2.GetParameter(3))
        sig_L_g .SetPointError(i, 0, fff2.GetParError(1))
        sig_T_g .SetPointError(i, 0, fff2.GetParError(0))
        sig_LT_g.SetPointError(i, 0, fff2.GetParError(2))
        sig_TT_g.SetPointError(i, 0, fff2.GetParError(3))

        # ---------------- plotting canvases --------------------
        c1 = TCanvas()
        g_plot_err.SetTitle(f"t = {t_list[i]:.3f}")
        g_plot_err.GetXaxis().SetTitle("φ [deg]")
        g_plot_err.GetYaxis().SetTitle("ε")
        g_plot_err.GetZaxis().SetTitle("σ")
        g_plot_err.Draw("perr")
        if i == 0:  c1.Print(outputpdf+"(")
        else:       c1.Print(outputpdf)
        c1.Close()

        # unseparated φ plots
        c2 = TCanvas()
        glo.SetMarkerStyle(5)
        ghi.SetMarkerStyle(4);  ghi.SetMarkerColor(ROOT.kRed)
        mg = ROOT.TMultiGraph()
        mg.Add(glo); mg.Add(ghi); mg.Draw("AP")
        mg.GetXaxis().SetTitle("φ [deg]")
        mg.GetYaxis().SetTitle("Unsep σ [nb/GeV²]")
        flo.SetLineWidth(2);    fhi.SetLineWidth(2); fhi.SetLineColor(ROOT.kRed)
        flo.Draw("same");       fhi.Draw("same")
        leg = ROOT.TLegend(0.8,0.8,0.96,0.96)
        leg.AddEntry(glo,"ε_low","p"); leg.AddEntry(ghi,"ε_high","p")
        leg.Draw()
        c2.Print(outputpdf);  c2.Close()

        # σL/T evolution
        c3 = TCanvas()
        sigL_change.Draw("AP*")
        sigT_change.SetMarkerColor(ROOT.kRed); sigT_change.Draw("P*same")
        c3.Print(outputpdf); c3.Close()

        # separated σ vs t
        for g,cname in [(sig_L_g ,"σ_L"),(sig_T_g ,"σ_T"),
                        (sig_LT_g,"σ_LT"),(sig_TT_g,"σ_TT")]:
            c = TCanvas()
            g.SetTitle(cname); g.Draw("AP")
            c.Print(outputpdf); c.Close()

        # unsep vs ε at this t
        g_unsep_lo.SetPoint(i, LOEPS, lo_cross_sec[i])
        g_unsep_hi.SetPoint(i, HIEPS, hi_cross_sec[i])
        g_unsep_lo.SetPointError(i, 0, lo_cross_err[i])
        g_unsep_hi.SetPointError(i, 0, hi_cross_err[i])

        # write separated values
        fn_sep = (f"{LTANAPATH}/src/{ParticleType}/xsects/"
                  f"x_sep.{polID}_Q{Q2.replace('p','')}W{W.replace('p','')}.dat")
        mode = 'w' if i==0 else 'a'
        with open(fn_sep, mode) as f:
            f.write(f"{fff2.GetParameter(1):.4f} {fff2.GetParError(1):.4f} "
                    f"{fff2.GetParameter(0):.4f} {fff2.GetParError(0):.4f} "
                    f"{fff2.GetParameter(2):.4f} {fff2.GetParError(2):.4f} "
                    f"{fff2.GetParameter(3):.4f} {fff2.GetParError(3):.4f} "
                    f"{fff2.GetChisquare():.2f} {t_list[i]:.4f} "
                    f"{w_list[i]:.3f} {q2_list[i]:.3f} {theta_list[i]:.3f}\n")

    # -------- end per-t-bin loop --------------------------------
    return t_list, sig_L_g, sig_T_g, sig_LT_g, sig_TT_g, g_unsep_lo, g_unsep_hi

# ------------------------------------------------------------------
# file paths & driver call
# ------------------------------------------------------------------
fn_lo = (f"{LTANAPATH}/src/{ParticleType}/xsects/"
         f"x_unsep.{polID}_Q{Q2.replace('p','')}W{W.replace('p','')}"
         f"_{LOEPS*100:.0f}.dat")
fn_hi = (f"{LTANAPATH}/src/{ParticleType}/xsects/"
         f"x_unsep.{polID}_Q{Q2.replace('p','')}W{W.replace('p','')}"
         f"_{HIEPS*100:.0f}.dat")

(t_list, g_sig_l_total, g_sig_t_total,
 g_sig_lt_total, g_sig_tt_total,
 g_unsep_lo,     g_unsep_hi) = single_setting(Q2, W, fn_lo, fn_hi)

# ------------------------------------------------------------------
# ε-slope plots (unchanged)
# ------------------------------------------------------------------
f_lin = TF1("f_lin", "[0]+[1]*x", 0, 1)
for i in range(g_unsep_lo.GetN()):
    c = TCanvas()
    mg = ROOT.TMultiGraph()
    x_lo,y_lo = ctypes.c_double(), ctypes.c_double()
    x_hi,y_hi = ctypes.c_double(), ctypes.c_double()
    g_unsep_lo.GetPoint(i,x_lo,y_lo)
    g_unsep_hi.GetPoint(i,x_hi,y_hi)
    g_lo_evt = TGraphErrors(1, array('d',[x_lo.value]), array('d',[y_lo.value]),
                            array('d',[g_unsep_lo.GetErrorX(i)]),
                            array('d',[g_unsep_lo.GetErrorY(i)]))
    g_hi_evt = TGraphErrors(1, array('d',[x_hi.value]), array('d',[y_hi.value]),
                            array('d',[g_unsep_hi.GetErrorX(i)]),
                            array('d',[g_unsep_hi.GetErrorY(i)]))
    mg.Add(g_lo_evt); mg.Add(g_hi_evt); mg.Draw("AP")
    mg.GetXaxis().SetTitle("ε"); mg.GetYaxis().SetTitle("Unsep σ [nb/GeV²]")
    mg.Fit(f_lin,"MRQ")
    f_lin.SetLineWidth(2); f_lin.Draw("same")
    leg = ROOT.TLegend(0.75,0.8,0.95,0.95)
    leg.AddEntry(g_lo_evt,f"-t={t_list[i]:.3f}",""); leg.AddEntry(g_lo_evt,"ε_low","p")
    leg.AddEntry(g_hi_evt,"ε_high","p"); leg.Draw()
    c.Print(outputpdf); c.Close()

# ------------------------------------------------------------------
# separated σ_L & σ_T vs t with exponentials
# ------------------------------------------------------------------
for g,color,name in [(g_sig_l_total, ROOT.kBlack,"σ_L"),
                     (g_sig_t_total, ROOT.kRed,  "σ_T")]:
    c = TCanvas()
    g.SetMarkerColor(color); g.SetLineColor(color); g.Draw("AP")
    f_exp = TF1("f_exp","[0]*exp(-[1]*x)",0,2)
    g.Fit(f_exp,"MRQ"); f_exp.SetLineColor(color); f_exp.Draw("same")
    c.Print(outputpdf); c.Close()

# ------------------------------------------------------------------
# σ_LT and σ_TT exponentials
# ------------------------------------------------------------------
for g,name in [(g_sig_lt_total,"σ_LT"),(g_sig_tt_total,"σ_TT")]:
    c = TCanvas(); g.Draw("AP")
    f_exp = TF1("f_exp","[0]*exp(-[1]*x)",0,2)
    g.Fit(f_exp,"MRQ"); f_exp.Draw("same"); c.Print(outputpdf); c.Close()

# ------------------------------------------------------------------
# close multi-page PDF
# ------------------------------------------------------------------
dummy = TCanvas(); dummy.Print(outputpdf+"]")
print(f"\nPlots written to {outputpdf}")
