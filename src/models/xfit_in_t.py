#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-27 06:02:39 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TF1, TGraph, TGraphErrors, TCanvas
import math

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

################################################################################################################################################

# Define constants
pi = 3.14159
mtar_gev = 0.93827231
mpipl=0.139570
mkpl=0.493677

hi_bound =  0.7;
lo_bound = -0.1;

def x_fit_in_t():
    lt_ratio_file = open("x_sep/LT_ratio.txt", "w")

    single_setting("160", lt_ratio_file)
    single_setting("245", lt_ratio_file)

    lt_ratio_file.close()

def single_setting(closest_date, q2_set):
    prv_par_vec = []
    g_vec = []
    w_vec = []
    q2_vec = []
    th_vec = []
    logq2_vec = []
    par_vec = []
    par_err_vec = []
    par_chi2_vec = []

    t0, t1, t2, t3 = 0, 0, 0, 0
    l0, l1, l2, l3 = 0, 0, 0, 0
    lt0, lt1, lt2, lt3 = 0, 0, 0, 0
    tt0, tt1, tt2, tt3 = 0, 0, 0, 0

    fit_status = TText()

    fn_sep = "x_sep.pl_" + q2_set

    nsep = TNtuple("nsep", "nsep", "t:t_e:l:l_e:lt:lt_e:tt:tt_e:chi:t:t_min:w:q2")
    nsep.ReadFile(fn_sep)

    prv_par_vec = []
    para_file_in = CACHEPATH + closest_date + "par.pl_" + q2_set
    try:
        with open(para_file_in, 'r') as para_file_in:
            for line in para_file_in:
                data = line.split()
                par, par_err, indx, chi2 = map(float, data)
                print("!!!!  {}".format(par))
                prv_par_vec.append(par)
    except FileNotFoundError:
        print("File {} not found.".format(para_file_in))

    t0, t1, t2, t3, l0, l1, l2, l3, lt0, lt1, lt2, lt3, tt0, tt1, tt2, tt3 = prv_par_vec[:16]
    

    ave_file_in = LTANAPATH + "src/averages/avek." + q2_set + ".dat"
    with open(ave_file_in, 'r') as ave_file:
        for line in ave_file:
            w, w_e, q2, q2_e, tt, tt_e, thetacm, it = map(float, line.strip().split())

            g = 1 / ((w**2) - (m_p**2))**2
            g_vec.append(g)
            w_vec.append(w)
            q2_vec.append(q2)
            th_vec.append(thetacm)
            logq2_vec.append(log(q2))

    g_sigt_prv = TGraph()
    g_sigl_prv = TGraph()
    g_siglt_prv = TGraph()
    g_sigtt_prv = TGraph()

    g_sigt_fit = TGraphErrors()
    g_sigl_fit = TGraphErrors()
    g_siglt_fit = TGraphErrors()
    g_sigtt_fit = TGraphErrors()

    g_sigt_fit_tot = TGraph()
    g_sigl_fit_tot = TGraph()
    g_siglt_fit_tot = TGraph()
    g_sigtt_fit_tot = TGraph()

    c1 = TCanvas("c1", "c1", 800, 800)
    c1.Divide(2, 2)

    c2 = TCanvas("c2", "c2", 800, 800)
    c2.Divide(2, 2)

    n1.Draw("u:u_min", "", "goff")
    u_umin_map = TGraph(n1.GetSelectedRows(), n1.GetV1(), n1.GetV2())

    u_list = u_umin_map.GetX()
    u_min_list = u_umin_map.GetY()

    ########
    # SigT #
    ########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig T")
    
    c1.cd(1).SetLeftMargin(0.12)
    n1.Draw("t:u:t_e", "", "goff")
    g_sigt = TGraphErrors(n1.GetSelectedRows(), n1.GetV2(), n1.GetV1(), [0] * n1.GetSelectedRows(), n1.GetV3())
    
    f_sigT_pre = TF1("sig_T_pre", fun_Sig_T, 0, 0.5, 2)
    f_sigT_pre.SetParameters(t0, t1)

    for i in range(len(w_vec)):
        sigt_X_pre = 0.0

        q2_term = t2 / logq2_vec[i] + t3 * g_sigt.GetX()[i] / logq2_vec[i]
        q2_dep = sqrt(q2_vec[i])

        sigt_X_pre = (f_sigT_pre.Eval(g_sigt.GetX()[i]) / q2_dep + q2_term) * g_vec[i]
        g_sigt_prv.SetPoint(i, g_sigt.GetX()[i], sigt_X_pre)

        sigt_X_fit = ((g_sigt.GetY()[i]) / g_vec[i] - q2_term) * q2_dep
        sigt_X_fit_err = g_sigt.GetEY()[i] / g_vec[i] * q2_dep

        g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)
    
    g_max = g_sigt.GetYaxis().GetXmax()
    gp_max = max(g_sigt_prv.GetN(), key=lambda i: g_sigt_prv.GetY()[i])
    g_min = g_sigt.GetYaxis().GetXmin()
    gp_min = min(g_sigt_prv.GetN(), key=lambda i: g_sigt_prv.GetY()[i])

    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_sigt.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_sigt.SetMinimum(gp_min - difff)

    g_sigt.SetTitle("Sig T")

    g_sigt.SetMarkerStyle(5)
    g_sigt.Draw("AP")

    g_sigt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt.GetXaxis().CenterTitle()
    g_sigt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{T} [#mub/GeV^{2}]")
    g_sigt.GetYaxis().SetTitleOffset(1.5)
    g_sigt.GetYaxis().SetTitleSize(0.035)
    g_sigt.GetYaxis().CenterTitle()

    g_sigt_prv.SetMarkerColor(4)
    g_sigt_prv.SetMarkerStyle(21)
    g_sigt_prv.Draw("P")
        
    c2.cd(1)
    g_sigt_fit.SetTitle("Sigma T Model Fit")
    g_sigt_fit.Draw("A*")

    f_sigT = ROOT.TF1("sig_T", fun_Sig_T, 0, 0.5, 2)
    f_sigT.SetParameters(t0, t1)

    fit_t_result = g_sigt_fit.Fit(f_sigT, "S")

    for i in range(len(w_vec)):
        sigt_X = 0.0

        q2_term = t2 / logq2_vec[i] + t3 * g_sigt.GetX()[i] / logq2_vec[i]

        q2_dep = sqrt(q2_vec[i])

        sigt_X = (f_sigT.Eval(g_sigt.GetX()[i]) / q2_dep + q2_term) * g_vec[i]

        print(f_sigT.Eval(g_sigt.GetX()[i]), sigt_X)

        g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)
    
    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + fit_t_result.Status())

    c1.cd(1)

    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    t0 = f_sigT.GetParameter(0)
    t1 = f_sigT.GetParameter(1)
    
    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + fit_t_result.Status())

    c1.cd(1)

    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    t0 = f_sigT.GetParameter(0)
    t1 = f_sigT.GetParameter(1)

    ########
    # SigL #
    ########
    
    print("/*--------------------------------------------------*/")
    print("Fit for Sig L")

    c1.cd(2).SetLeftMargin(0.12)
    n1.Draw("l:u:l_e", "", "goff")

    f_sigL_pre = TF1("sig_L", fun_Sig_L, 0, 0.5, 2)
    f_sigL_pre.SetParameters(l0, l1)
    
    g_sigl = TGraphErrors(n1.GetSelectedRows(), n1.GetV2(), n1.GetV1(), 0, n1.GetV3())

    for i in range(len(w_vec)):
        sigl_X_pre = 0.0
        q2_term = l2 / q2_vec[i] + l3 * g_sigl.GetX()[i] / q2_vec[i]
        q2_dep = q2_vec[i] * q2_vec[i]

        sigl_X_pre = (f_sigL_pre.Eval(g_sigl.GetX()[i]) / q2_dep + q2_term) * g_vec[i]
        g_sigl_prv.SetPoint(i, g_sigl.GetX()[i], sigl_X_pre)

        sigl_X_fit = (g_sigl.GetY()[i] / g_vec[i] - q2_term) * q2_dep
        sigl_X_fit_err = g_sigl.GetEY()[i] / g_vec[i] * q2_dep

        g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

    g_max = g_sigl.GetYaxis().GetXmax()
    gp_max = max(g_sigl_prv.GetN(), max(g_sigl_prv.GetY()))

    g_min = g_sigl.GetYaxis().GetXmin()
    gp_min = min(g_sigl_prv.GetN(), min(g_sigl_prv.GetY()))

    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_sigl.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_sigl.SetMinimum(gp_min - difff)
        
    g_sigl.SetTitle("Sig L")

    g_sigl.SetMarkerStyle(5)
    g_sigl.Draw("AP")

    g_sigl.GetXaxis().SetTitle("#it{-u} [GeV^{2}]")
    g_sigl.GetXaxis().CenterTitle()
    g_sigl.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{L} [#mub/GeV^{2}]")
    g_sigl.GetYaxis().SetTitleOffset(1.5)
    g_sigl.GetYaxis().SetTitleSize(0.035)
    g_sigl.GetYaxis().CenterTitle()

    g_sigl_prv.SetMarkerColor(4)
    g_sigl_prv.SetMarkerStyle(21)
    g_sigl_prv.Draw("P")

    c2.cd(2)
    g_sigl_fit.SetTitle("Sigma L Model Fit")
    g_sigl_fit.Draw("A*")

    
