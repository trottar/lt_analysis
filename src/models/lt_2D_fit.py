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
# Importing utility functions
# ------------------------------------------------------------------

sys.path.append("utility")
from utility import adapt_limits

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
# main per-setting routine
# ------------------------------------------------------------------
def single_setting(q2_set, w_set, fn_lo, fn_hi):

    set_val(LOEPS, HIEPS)              # hand LO/HI ε to lt_active
    eps_diff = HIEPS - LOEPS

    # → graphs initialisation (unchanged) --------------------------
    sig_L_g  = TGraphErrors()
    sig_T_g  = TGraphErrors()
    sig_LT_g = TGraphErrors()
    sig_TT_g = TGraphErrors()

    sig_lo = TGraphErrors()
    sig_hi = TGraphErrors()
    sig_diff_g = TGraphErrors()

    # --------------------------------------------------------------
    nlo = TNtuple("nlo", "nlo", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nlo.ReadFile(fn_lo)
    nhi = TNtuple("nhi", "nhi", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nhi.ReadFile(fn_hi)

    # gather unique t-bins  ---------------------------------------
    q2_list, w_list, theta_list, t_list = [], [], [], []
    lo_eps_list, hi_eps_list           = [], []

    for evt in nlo:
        if evt.t not in t_list:
            q2_list.append(evt.Q2); w_list.append(evt.w)
            theta_list.append(evt.theta); t_list.append(evt.t)
            lo_eps_list.append(evt.eps)

    tmp_t_list = []
    for evt in nhi:
        if evt.t not in tmp_t_list:
            tmp_t_list.append(evt.t)
            hi_eps_list.append(evt.eps)

    t_bin_num = len(t_list)

    # storage for integrated unsep cross sections -----------------
    lo_cross_sec      = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec      = np.zeros(t_bin_num, dtype=float)
    lo_cross_sec_err  = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec_err  = np.zeros(t_bin_num, dtype=float)

    # loop over t-bins --------------------------------------------
    for i in range(t_bin_num):

        print("\n/*--------------------------------------------------*/")
        print(f" Starting t-bin {i+1} (t={t_list[i]:.4f})...")
        print("/*--------------------------------------------------*/\n")

        tcut   = f"t=={t_list[i]} && x!=0.0"
        lo_eps = lo_eps_list[i]
        hi_eps = hi_eps_list[i]

        # low-ε graph --------------------------------------------
        nlo.Draw("x:phi:dx", tcut, "goff")
        glo_tmp = TGraphErrors()
        for j in range(nlo.GetSelectedRows()):
            glo_tmp.SetPoint(j, nlo.GetV2()[j], nlo.GetV1()[j])
            glo_tmp.SetPointError(j, 0, nlo.GetV3()[j])
        glo = glo_tmp.Clone("glo")

        ave_sig_lo = glo.GetMean(2)
        err_sig_lo = glo.GetRMS(2)
        sig_lo.SetPoint(sig_lo.GetN(), t_list[i], ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)

        # high-ε graph -------------------------------------------
        nhi.Draw("x:phi:dx", tcut, "goff")
        ghi_tmp = TGraphErrors()
        for j in range(nhi.GetSelectedRows()):
            ghi_tmp.SetPoint(j, nhi.GetV2()[j], nhi.GetV1()[j])
            ghi_tmp.SetPointError(j, 0, nhi.GetV3()[j])
        ghi = ghi_tmp.Clone("ghi")

        ave_sig_hi = ghi.GetMean(2)
        err_sig_hi = ghi.GetRMS(2)
        sig_hi.SetPoint(sig_hi.GetN(), t_list[i], ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)

        # build 2D graph (φ,ε,σ) ----------------------------------
        g_plot_err = TGraph2DErrors()
        g_xx = ctypes.c_double(0); g_yy = ctypes.c_double(0)

        for ii in range(glo.GetN()):
            glo.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((glo.GetErrorY(ii)/g_yy.value)**2 +
                                 (pt_to_pt_systematic_error/100)**2) * g_yy.value
            lo_cross_sec_err[i] += 1/(g_yy_err**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, lo_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0, 0,
                                     math.sqrt(glo.GetErrorY(ii)**2 +
                                               (pt_to_pt_systematic_error/100)**2))

        for ii in range(ghi.GetN()):
            ghi.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((ghi.GetErrorY(ii)/g_yy.value)**2 +
                                 (pt_to_pt_systematic_error/100)**2) * g_yy.value
            hi_cross_sec_err[i] += 1/(g_yy_err**2)
            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, hi_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0, 0,
                                     math.sqrt(ghi.GetErrorY(ii)**2 +
                                               (pt_to_pt_systematic_error/100)**2))

        lo_cross_sec_err[i] = (1/math.sqrt(lo_cross_sec_err[i])
                               if lo_cross_sec_err[i] else -1000)
        hi_cross_sec_err[i] = (1/math.sqrt(hi_cross_sec_err[i])
                               if hi_cross_sec_err[i] else -1000)

        # ---------------------------------------------------------
        # *** adaptive fit starts here ***
        # ---------------------------------------------------------
        fff2 = TF2("fff2",
                   "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] "
                   "+ y*cos(2*x*0.017453)*[3]",
                   0, 360, 0.0, 1.0)

        # initial guesses
        fff2.SetParameters(ave_sig_hi, ave_sig_lo, 0.0, 0.0)
        for p,(lo,hi) in enumerate(zip(HARD_LO, HARD_HI)):
            fff2.SetParLimits(p, lo, hi)

        # canvases tracking change
        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()
        sigL_change.SetTitle(f"t = {t_list[i]:.3f}")
        sigT_change.SetTitle(f"t = {t_list[i]:.3f}")

        for step in fit_plan:

            # fix/release per plan
            for p in range(fff2.GetNpar()):
                (fff2.ReleaseParameter if p in step["free"]
                 else fff2.FixParameter)(p, fff2.GetParameter(p))

            # seed newly released params
            for p in step["free"]:
                if fff2.GetParError(p) == 0:
                    seed = 0.0 if p > 1 else ave_sig_lo
                    fff2.SetParameter(p, seed)

            # tighten search box
            adapt_limits(fff2, nsig=4, hard_lo=HARD_LO, hard_hi=HARD_HI)

            # fit quietly
            g_plot_err.Fit(fff2, "MRQ")

            # record evolution
            idx = sigL_change.GetN()
            sigL_change.SetPoint(idx, idx+1, fff2.GetParameter(1))
            sigL_change.SetPointError(idx, 0, fff2.GetParError(1))
            sigT_change.SetPoint(idx, idx+1, fff2.GetParameter(0))
            sigT_change.SetPointError(idx, 0, fff2.GetParError(0))
        # ---------------------------------------------------------
        # *** adaptive fit finished ***
        # ---------------------------------------------------------

        # -------- remainder of your original code ----------------
        # (parameter extraction, plotting, file output, etc.)
        # -- unchanged, still uses fff2.GetParameter / GetParError --
        # (snipped here for brevity; keep everything you had below)
        # ---------------------------------------------------------

    # end for-each t-bin
    return t_list
# ------------------------------------------------------------------
# paths & driver call
# ------------------------------------------------------------------
fn_lo = (f"{LTANAPATH}/src/{ParticleType}/xsects/"
         f"x_unsep.{polID}_Q{Q2.replace('p','')}W{W.replace('p','')}"
         f"_{LOEPS*100:.0f}.dat")
fn_hi = (f"{LTANAPATH}/src/{ParticleType}/xsects/"
         f"x_unsep.{polID}_Q{Q2.replace('p','')}W{W.replace('p','')}"
         f"_{HIEPS*100:.0f}.dat")

g_sig_l_total = TGraphErrors()
g_sig_t_total = TGraphErrors()
g_sig_lt_total = TGraphErrors()
g_sig_tt_total = TGraphErrors()
g_unsep_lo    = TGraphErrors()
g_unsep_hi    = TGraphErrors()

t_list = single_setting(Q2, W, fn_lo, fn_hi)

f_lin_l = TF1("f_lin_l", "[0]+[1]*x", 0.0, 1.0)
f_lin_t = TF1("f_lin_t", "[0]+[1]*x", 0.0, 1.0)

# Define the number of events in 'lo'
num_events = g_unsep_lo.GetN()

# Loop over each event in 'lo'
for i in range(num_events):

    # Create a canvas
    c = ROOT.TCanvas("c", "c", 800, 600)
    
    # Create a new TMultiGraph for each event
    g_unsep_mult = ROOT.TMultiGraph()
    
    # Create TGraphErrors for 'lo' event
    x_lo, y_lo = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_lo.GetPoint(i, x_lo, y_lo)
    x_err_lo = g_unsep_lo.GetErrorX(i)
    y_err_lo = g_unsep_lo.GetErrorY(i)
    g_lo_event = ROOT.TGraphErrors(1, array('d', [x_lo.value]), array('d', [y_lo.value]), array('d', [x_err_lo]), array('d', [y_err_lo]))
    
    # Add 'lo' event to the TMultiGraph
    g_unsep_mult.Add(g_lo_event)
    
    # Create TGraphErrors for 'hi' event
    x_hi, y_hi = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_hi.GetPoint(i, x_hi, y_hi)
    x_err_hi = g_unsep_hi.GetErrorX(i)
    y_err_hi = g_unsep_hi.GetErrorY(i)
    g_hi_event = ROOT.TGraphErrors(1, array('d', [x_hi.value]), array('d', [y_hi.value]), array('d', [x_err_hi]), array('d', [y_err_hi]))
    
    # Add 'hi' event to the TMultiGraph
    g_unsep_mult.Add(g_hi_event)

    g_lo_event.SetMarkerStyle(5)
    g_hi_event.SetMarkerStyle(5)
    
    # Draw TMultiGraph for the current event
    g_unsep_mult.Draw("AP")
    
    # Set axis titles and offsets
    g_unsep_mult.GetYaxis().SetTitle("Unseparated Cross Section [nb/GeV^{2}]")
    g_unsep_mult.GetYaxis().SetTitleOffset(1.2)
    g_unsep_mult.GetXaxis().SetTitle("#epsilon")
    g_unsep_mult.GetXaxis().SetTitleOffset(1.2)
    
    # Fit functions to 'lo' and 'hi' events
    f_lin = ROOT.TF1("f_lin", "[0]*x + [1]", 0, 1)  # Define fit function for 'lo'
    g_unsep_mult.Fit(f_lin, "MRQ")
        
    # Set line properties for 'lo' and 'hi' fits
    f_lin.SetLineColor(2)
    f_lin.SetLineWidth(2)    
    f_lin.Draw("same")
    
    # Create and draw TLegend
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

# Create TMultiGraph and add glo, ghi
g_sig_mult = ROOT.TMultiGraph()
g_sig_mult.Add(g_sig_l_total)
g_sig_mult.Add(g_sig_t_total)

#g_sig_l_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_l_total.GetYaxis().SetRangeUser(-5,30)
#g_sig_t_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_t_total.GetYaxis().SetRangeUser(-5,30)

g_sig_mult.Draw("AP")

g_sig_mult.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right) [nb/GeV^{2}]")
g_sig_mult.GetYaxis().SetTitleOffset(1.2)

g_sig_mult.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
g_sig_mult.GetXaxis().SetTitleOffset(1.2)

g_sig_l_total.Fit(f_exp_l, "MRQ")
g_sig_t_total.Fit(f_exp_t, "MRQ")

# Set properties for g_sig_l_total and g_sig_t_total
g_sig_l_total.SetLineColor(1)
g_sig_l_total.SetMarkerStyle(5)
g_sig_t_total.SetLineColor(2)
g_sig_t_total.SetMarkerColor(2)
g_sig_t_total.SetMarkerStyle(4)

# Set line properties for g_sig_l_total and g_sig_t_total
f_exp_l.SetLineColor(1)
f_exp_l.SetLineWidth(2)
f_exp_t.SetLineColor(2)
f_exp_t.SetLineWidth(2)
f_exp_t.SetLineStyle(2)

# Draw f_exp_l and f_exp_t on the same canvas
f_exp_l.Draw("same")
f_exp_t.Draw("same")

# Create and draw TLegend
leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
leg.SetFillColor(0)
leg.SetMargin(0.4)
leg.AddEntry(g_sig_l_total, "#it{#sigma}_{L}", "p")
leg.AddEntry(g_sig_t_total, "#it{#sigma}_{T}", "p")
leg.Draw()

c_total_l_t.Print(outputpdf)
c_total_l_t.Clear()

ROOT.gStyle.SetOptFit(1)

# Set the size of the statbox
ROOT.gStyle.SetStatW(0.15)  # Set width of statbox
ROOT.gStyle.SetStatH(0.15)  # Set height of statbox

# Set the position of the statbox
#ROOT.gStyle.SetStatX(0.8)  # Set x position of statbox (right edge)
#ROOT.gStyle.SetStatY(0.8)  # Set y position of statbox (top edge)

f_exp = TF1("f_exp", "[0]*exp(-[1]*x)", 0.0, 2.0)

c_total = TCanvas()

#g_sig_l_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_l_total.GetYaxis().SetRangeUser(-5,30)
g_sig_l_total.Draw("A*")
g_sig_l_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_t_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_t_total.GetYaxis().SetRangeUser(-5,30)
g_sig_t_total.SetMarkerColor(1)
g_sig_t_total.SetLineColor(1)
g_sig_t_total.Draw("A*")
g_sig_t_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_lt_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_lt_total.GetYaxis().SetRangeUser(-5,3)
g_sig_lt_total.Draw("A*")
g_sig_lt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_tt_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_tt_total.GetYaxis().SetRangeUser(-3,1)
g_sig_tt_total.Draw("A*")
g_sig_tt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf+')')
c_total.Clear()
