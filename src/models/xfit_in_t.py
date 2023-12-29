#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-29 12:47:02 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TGraph, TGraphErrors, TCanvas
from ROOT import TF1, TFitResultPtr
import math
import os

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
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Define constants
PI = ROOT.TMath.Pi()
m_p = 0.93827231

hi_bound =  0.7;
lo_bound = -0.1;

# Function for SigT
def fun_Sig_T(x, par):
    xx = x[0]
    f = par[0] + par[1]*xx
    return f

# Function for SigL
def fun_Sig_L(x, par):
    xx = x[0]
    f = par[0] + par[1]*xx
    return f

# Function for SigLT
def fun_Sig_LT(x, par):
    xx = x[0]
    f = par[0] + par[1]*xx
    return f

# Function for SigTT
def fun_Sig_TT(x, par):
    xx = x[0]
    f = par[0] + par[1]*xx
    return f

def x_fit_in_t(ParticleType, polID, closest_date, Q2):

    single_setting(ParticleType, polID, closest_date, Q2)

def single_setting(ParticleType, polID, dir_iter, q2_set):

    outputpdf  = OUTPATH + "/" + ParticleType + "_xfit_in_t.pdf"
    
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

    fn_sep = "{}/src/{}/xsects/x_sep.{}_{}.dat".format(LTANAPATH, ParticleType, polID, q2_set.replace("p",""), float(LOEPS)*100)
    nsep = TNtuple("nsep", "nsep", "sigt:sigt_e:sigl:sigl_e:siglt:siglt_e:sigtt:sigtt_e:chi:t:t_min:w:q2")
    nsep.ReadFile(fn_sep)

    prv_par_vec = []
    para_file_in =  "{}/{}/{}/{}/parameters/par.{}_{}.dat".format(CACHEPATH, USER, ParticleType, dir_iter, polID, q2_set.replace("p",""))
    print("Reading {}...".format(para_file_in))
    try:
        with open(para_file_in, 'r') as f:
            for line in f:
                data = line.split()
                par, par_err, indx, chi2 = map(float, data)
                print("  {}".format(par))
                prv_par_vec.append(par)
    except FileNotFoundError:
        print("File {} not found.".format(para_file_in))
        
    #t0, t1, t2, t3, l0, l1, l2, l3, lt0, lt1, lt2, lt3, tt0, tt1, tt2, tt3 = prv_par_vec[:16]
    t0, t1, t2, l0, l1, l2, lt0, lt1, lt2, tt0, tt1, tt2 = prv_par_vec[:12] # RLT too many parameters
    
    ave_file_in = "{}/src/{}/averages/avek.{}.dat".format(LTANAPATH, ParticleType, q2_set.replace("p",""))
    with open(ave_file_in, 'r') as f:
        for line in f:
            w, w_e, q2, q2_e, tt, tt_e, thetacm, it = map(float, line.strip().split())

            g = 1 / ((w**2) - (m_p**2))**2
            g_vec.append(g)
            w_vec.append(w)
            q2_vec.append(q2)
            th_vec.append(thetacm)
            logq2_vec.append(math.log(q2))

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

    nsep.Draw("t:t_min", "", "goff")
    t_tmin_map = TGraph(nsep.GetSelectedRows(), nsep.GetV1(), nsep.GetV2())

    t_list = t_tmin_map.GetX()
    t_min_list = t_tmin_map.GetY()

    ########
    # SigT #
    ########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig T")
    
    c1.cd(1).SetLeftMargin(0.12)
    nsep.Draw("sigt:t:sigt_e", "", "goff")
    g_sigt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), [0] * nsep.GetSelectedRows(), nsep.GetV3())
    
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
    g_sigt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [#mub/GeV^{2}]")
    g_sigt.GetYaxis().SetTitleOffset(1.5)
    g_sigt.GetYaxis().SetTitleSize(0.035)
    g_sigt.GetYaxis().CenterTitle()

    g_sigt_prv.SetMarkerColor(4)
    g_sigt_prv.SetMarkerStyle(21)
    g_sigt_prv.Draw("P")
        
    c2.cd(1)
    g_sigt_fit.SetTitle("Sigma T Model Fit")
    g_sigt_fit.Draw("A*")

    f_sigT = TF1("sig_T", fun_Sig_T, 0, 0.5, 2)
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
    nsep.Draw("sigl:t:sigl_e", "", "goff")

    f_sigL_pre = TF1("sig_L", fun_Sig_L, 0, 0.5, 2)
    f_sigL_pre.SetParameters(l0, l1)
    
    g_sigl = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), 0, nsep.GetV3())

    for i in range(len(w_vec)):
        
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

    g_sigl.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigl.GetXaxis().CenterTitle()
    g_sigl.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [#mub/GeV^{2}]")
    g_sigl.GetYaxis().SetTitleOffset(1.5)
    g_sigl.GetYaxis().SetTitleSize(0.035)
    g_sigl.GetYaxis().CenterTitle()

    g_sigl_prv.SetMarkerColor(4)
    g_sigl_prv.SetMarkerStyle(21)
    g_sigl_prv.Draw("P")

    c2.cd(2)
    g_sigl_fit.SetTitle("Sigma L Model Fit")
    g_sigl_fit.Draw("A*")

    f_sigL = TF1("sig_L", fun_Sig_L, 0, 0.5, 2)
    f_sigL.SetParameters(l0, l1)
    g_sigl_fit.Fit(f_sigL)

    for i in range(len(w_vec)):

        sigl_X = 0.0
        q2_term = l2 * q2_vec[i] + l3 * g_sigl.GetX()[i] * q2_vec[i]

        q2_dep = q2_vec[i] * q2_vec[i]

        sigl_X = (f_sigL.Eval(g_sigl.GetX()[i]) / q2_dep + q2_term) * g_vec[i]
        g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)
    
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)
    c1.cd(2)

    g_sigl_fit_tot.SetMarkerStyle(26)
    g_sigl_fit_tot.SetMarkerColor(2)
    g_sigl_fit_tot.SetLineColor(2)
    g_sigl_fit_tot.Draw("LP")

    l0 = f_sigL.GetParameter(0)
    l1 = f_sigL.GetParameter(1)

    par_vec.append(l0)
    par_vec.append(l1)
    par_vec.append(l2)
    par_vec.append(l3)

    par_err_vec.append(f_sigL.GetParError(0))
    par_err_vec.append(f_sigL.GetParError(1))
    par_err_vec.append(0)
    par_err_vec.append(0)

    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    
    ########
    # SigLT #
    ########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig LT")

    c1.cd(3).SetLeftMargin(0.12)
    nsep.Draw("siglt:t:siglt_e", "", "goff")

    f_sigLT_pre = TF1("sig_LT", fun_Sig_LT, 0, 0.5, 2)
    f_sigLT_pre.SetParameters(lt0, lt1)
    
    g_siglt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), ROOT.nullptr, nsep.GetV3())

    for i in range(len(w_vec)):

        siglt_X_pre = 0.0
        #q2_term = lt2 / logq2_vec[i] + lt3 * g_siglt.GetX()[i] / logq2_vec[i]
        q2_term = 0.0 # RLT Too many parameters

        q2_dep = q2_vec[i]

        siglt_X_pre = (f_sigLT_pre.Eval(g_siglt.GetX()[i]) / q2_dep + q2_term) * g_vec[i] * math.sin(th_vec[i] * PI / 180)
        g_siglt_prv.SetPoint(i, g_sigl.GetX()[i], siglt_X_pre)

        siglt_X_fit, siglt_X_fit_err = 0.0, 1.0

        if th_vec[i] != 180:
            siglt_X_fit = (g_siglt.GetY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180) - q2_term) * q2_dep
            siglt_X_fit_err = g_siglt.GetEY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180) * q2_dep

        g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

    g_max = g_siglt.GetYaxis().GetXmax()
    gp_max = max(g_siglt_prv.GetN(), g_siglt_prv.GetY())
    g_min = g_siglt.GetYaxis().GetXmin()
    gp_min = min(g_siglt_prv.GetN(), g_siglt_prv.GetY())

    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_siglt.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_siglt.SetMinimum(gp_min - difff)
        
    g_siglt.SetTitle("Sig LT")

    g_siglt.SetMarkerStyle(5)
    g_siglt.Draw("AP")

    g_siglt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_siglt.GetXaxis().CenterTitle()
    g_siglt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [#mub/GeV^{2}]")
    g_siglt.GetYaxis().SetTitleOffset(1.5)
    g_siglt.GetYaxis().SetTitleSize(0.035)
    g_siglt.GetYaxis().CenterTitle()

    g_siglt_prv.SetMarkerColor(4)
    g_siglt_prv.SetMarkerStyle(21)
    g_siglt_prv.Draw("P")

    c2.cd(3)
    g_siglt_fit.SetTitle("Sigma LT Model Fit")
    g_siglt_fit.Draw("A*")

    f_sigLT = TF1("sig_LT", fun_Sig_LT, 0, 0.5, 2)
    f_sigLT.SetParameters(lt0, lt1)

    g_siglt_fit.Fit(f_sigLT)
    
    for i in range(len(w_vec)):
        siglt_X = 0.0
        #q2_term = lt2 / logq2_vec[i] + lt3 * g_siglt.GetX()[i] / logq2_vec[i]
        q2_term = 0.0 # RLT Too many parameters

        q2_dep = q2_vec[i]

        if th_vec[i] != 180:
            siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i]) / q2_dep + q2_term) * g_vec[i] * sin(th_vec[i] * PI / 180)

        g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)

    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)

    c1.cd(3)

    g_siglt_fit_tot.SetMarkerStyle(26)
    g_siglt_fit_tot.SetMarkerColor(2)
    g_siglt_fit_tot.SetLineColor(2)
    g_siglt_fit_tot.Draw("LP")

    lt0 = f_sigLT.GetParameter(0)
    lt1 = f_sigLT.GetParameter(1)
        
    par_vec.append(lt0)
    par_vec.append(lt1)
    par_vec.append(lt2)
    par_vec.append(lt3)

    par_err_vec.append(f_sigLT.GetParError(0))
    par_err_vec.append(f_sigLT.GetParError(1))
    par_err_vec.append(0.0)
    par_err_vec.append(0.0)

    par_chi2_vec.append(f_sigLT.GetChisquare())
    par_chi2_vec.append(f_sigLT.GetChisquare())
    par_chi2_vec.append(f_sigLT.GetChisquare())
    par_chi2_vec.append(f_sigLT.GetChisquare())

    ########
    # SigTT #
    ########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig TT")

    c1.cd(4).SetLeftMargin(0.12)
    nsep.Draw("sigtt:t:sigtt_e", "", "goff")
    
    f_sigTT_pre = TF1("sig_TT", fun_Sig_TT, 0, 0.5, 2)
    f_sigTT_pre.SetParameters(tt0, tt1)
    
    g_sigtt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), [0]*nsep.GetSelectedRows(), nsep.GetV3())

    for i in range(len(w_vec)):
        sigtt_X_pre = 0.0
        #q2_term = tt2 / logq2_vec[i] + tt3 * g_sigtt.GetX()[i] / logq2_vec[i]
        q2_term = 0.0 # RLT Too many parameters
        q2_dep = q2_vec[i]

        sigtt_X_pre = (f_sigTT_pre.Eval(g_sigtt.GetX()[i]) / q2_dep + q2_term) * g_vec[i] * \
                      sin(th_vec[i] * PI / 180) * sin(th_vec[i] * PI / 180)

        g_sigtt_prv.SetPoint(i, nsep.GetV2()[i], sigtt_X_pre)

        sigtt_X_fit, sigtt_X_fit_err = 0.0, 1.0

        if th_vec[i] != 180:
            sigtt_X_fit = (g_sigtt.GetY()[i] / g_vec[i] / sin(th_vec[i] * PI / 180) / sin(th_vec[i] * PI / 180) - q2_term) * q2_dep
            sigtt_X_fit_err = g_sigtt.GetEY()[i] / g_vec[i] / sin(th_vec[i] * PI / 180) / sin(th_vec[i] * PI / 180) * q2_dep

        g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)
    
    g_max = g_sigtt.GetYaxis().GetXmax()
    gp_max = ROOT.TMath.MaxElement(g_sigtt_prv.GetN(), g_sigtt_prv.GetY())

    g_min = g_sigtt.GetYaxis().GetXmin()
    gp_min = ROOT.TMath.MinElement(g_sigtt_prv.GetN(), g_sigtt_prv.GetY())

    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_sigtt.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_sigtt.SetMinimum(gp_min - difff)

    g_sigtt.SetTitle("Sig TT")

    g_sigtt.SetMarkerStyle(5)
    g_sigtt.Draw("AP")

    g_sigtt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigtt.GetXaxis().CenterTitle()
    g_sigtt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [#mub/GeV^{2}]")
    g_sigtt.GetYaxis().SetTitleOffset(1.5)
    g_sigtt.GetYaxis().SetTitleSize(0.035)
    g_sigtt.GetYaxis().CenterTitle()

    g_sigtt_prv.SetMarkerColor(4)
    g_sigtt_prv.SetMarkerStyle(21)
    g_sigtt_prv.Draw("P")

    c2.cd(4)

    g_sigtt_fit.SetTitle("Sigma TT Model Fit")
    g_sigtt_fit.Draw("A*")

    f_sigTT = TF1("sig_TT", fun_Sig_TT, 0, 0.5, 2)
    f_sigTT.SetParameters(tt0, tt1)
    g_sigtt_fit.Fit(f_sigTT)
        
    for i in range(len(w_vec)):
        sigtt_X = 0.0
        #q2_term = tt2 / logq2_vec[i] + tt3 * g_sigtt.GetX()[i] / logq2_vec[i]
        q2_term = 0.0 # RLT Too many parameters
        q2_dep = q2_vec[i]

        sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i]) / q2_dep + q2_term) * g_vec[i] * \
                  sin(th_vec[i] * PI / 180) * sin(th_vec[i] * PI / 180)

        g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)
        
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)

    c1.cd(4)

    g_sigtt_fit_tot.SetMarkerStyle(26)
    g_sigtt_fit_tot.SetMarkerColor(2)
    g_sigtt_fit_tot.SetLineColor(2)
    g_sigtt_fit_tot.Draw("LP")

    c1.Print(outputpdf+'(')
    c2.Print(outputpdf+')')
    
    with open("{}/src/{}/parameters/par.{}_{}.dat".format(LTANAPATH, ParticleType, polID, q2_set), 'w') as f:
        sys.stdout.write(fixed)
        f.write(fixed)
        f.write("{:.5f}\n".format(par_value))

        for i in range(len(par_vec)):
            sys.stdout.write("{:>10}   {:>15.4f} {:>14} {:>10}\n".format(par_vec[i], par_err_vec[i], par_chi2_vec[i], i))
            f.write("{:>12}   {:>15.4f} {:>12} {:>5}\n".format(par_vec[i], par_err_vec[i], par_chi2_vec[i], i))
        
    l_t_ratio = 0.0
    l_t_ratio_err = 0.0

    lt_ratio = TGraphErrors()

    for i in range(len(w_vec)):
        l_t_ratio = f_sigL.Eval(g_sigt.GetX()[i]) / f_sigT.Eval(g_sigt.GetX()[i])

        l_t_ratio_err = (0.03 / g_sigl_fit_tot.GetY()[i])**2 + (0.03 / g_sigt_fit_tot.GetY()[i])**2
        l_t_ratio_err = pow(l_t_ratio_err, 0.5) * l_t_ratio
