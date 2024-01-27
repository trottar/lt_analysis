#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-27 16:00:53 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TGraph, TGraphErrors, TCanvas
from ROOT import TF2, TFitResultPtr
import math
import os, sys

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
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Define constants
PI = ROOT.TMath.Pi()
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

def x_fit_in_t(ParticleType, pol_str, closest_date, Q2, W):

    single_setting(ParticleType, pol_str, closest_date, Q2, W)

def single_setting(ParticleType, pol_str, dir_iter, q2_set, w_set):

    hi_bound =  0.7;
    lo_bound = -0.1;

    tav = (0.1112 + 0.0066*math.log(float(q2_set.replace("p","."))))*float(q2_set.replace("p","."))

    # Function for SigT
    def fun_Sig_T(x, par):
        tt = abs(x[0])
        qq = abs(x[1])
        ftav = (abs(tt)-tav)/tav
        print("Calculating params for func_SigT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
        return f

    # Function for SigL
    def fun_Sig_L(x, par):
        tt = abs(x[0])
        qq = abs(x[1])
        print("Calculating params for func_SigL...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
        return f

    # Function for SigLT
    # thetacm term is defined on function calling
    def fun_Sig_LT(x, par):
        tt = abs(x[0])
        qq = abs(x[1])
        print("Calculating params for func_SigLT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))
        return f

    # Function for SigTT
    # thetacm term is defined on function calling
    def fun_Sig_TT(x, par):
        tt = abs(x[0])
        qq = abs(x[1])
        if pol_str == "pl":
            f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
        print("Calculating params for func_SigTT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]*qq*math.exp(-qq))*f_tt
        return f
    
    outputpdf  = "{}/{}_xfit_in_t_Q{}W{}.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
    
    prv_par_vec = []
    g_vec = []
    w_vec = []
    q2_vec = []
    th_vec = []
    par_vec = []
    par_err_vec = []
    par_chi2_vec = []

    t0, t1, t2, t3 = 0, 0, 0, 0
    l0, l1, l2, l3 = 0, 0, 0, 0
    lt0, lt1, lt2, lt3 = 0, 0, 0, 0
    tt0, tt1, tt2, tt3 = 0, 0, 0, 0

    fit_status = TText()

    fn_sep = "{}/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(LTANAPATH, ParticleType, pol_str, q2_set.replace("p",""), w_set.replace("p",""))
    nsep = TNtuple("nsep", "nsep", "sigt:sigt_e:sigl:sigl_e:siglt:siglt_e:sigtt:sigtt_e:chi:t:t_min:w:q2")
    nsep.ReadFile(fn_sep)

    print("Reading {}...".format(fn_sep))
    for entry in nsep:
        print("sigt: {}, sigt_e: {}, sigl: {}, sigl_e: {}, siglt: {}, siglt_e: {}, sigtt: {}, sigtt_e: {}, chi: {}, t: {}, t_min: {}, w: {}, q2: {}".format(
            entry.sigt, entry.sigt_e, entry.sigl, entry.sigl_e, entry.siglt, entry.siglt_e, entry.sigtt, entry.sigtt_e, entry.chi, entry.t, entry.t_min, entry.w, entry.q2
        ))

    prv_par_vec = []
    para_file_in =  "{}/{}/{}/{}/parameters/par.{}_Q{}W{}.dat".format(CACHEPATH, USER, ParticleType, dir_iter, pol_str, q2_set.replace("p",""), w_set.replace("p",""))
    print("Reading {}...".format(para_file_in))
    try:
        with open(para_file_in, 'r') as f:
            for line in f:
                data = line.split()
                par, par_err, indx, chi2 = map(float, data)
                print("  {} {} {} {}".format(par, par_err, indx, chi2))
                prv_par_vec.append(par)
    except FileNotFoundError:
        print("File {} not found.".format(para_file_in))
        
    t0, t1, t2, t3, l0, l1, l2, l3, lt0, lt1, lt2, lt3, tt0, tt1, tt2, tt3 = prv_par_vec[:16]
    
    ave_file_in = "{}/src/{}/averages/avek.Q{}W{}.dat".format(LTANAPATH, ParticleType, q2_set.replace("p",""), w_set.replace("p",""))
    with open(ave_file_in, 'r') as f:
        for line in f:
            w, w_e, q2, q2_e, tt, tt_e, thetacm, it = map(float, line.strip().split())

            if pol_str == "pl":
                g = 1 / ((w**2) - (m_p**2))**2
            else:
                g = 1 / ((w**2) - (m_n**2))**2
            g_vec.append(g)
            w_vec.append(w)
            q2_vec.append(q2)
            th_vec.append(thetacm)

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
    #t_tmin_map = TGraph(nsep.GetSelectedRows(), nsep.GetV1(), nsep.GetV2())
    t_tmin_map = TGraph()
    # Use SetPoint to fill the graph
    for i in range(nsep.GetSelectedRows()):
        t_tmin_map.SetPoint(i, nsep.GetV1()[i], nsep.GetV2()[i])

    t_list = t_tmin_map.GetX()
    t_min_list = t_tmin_map.GetY()

    ########
    # SigT #
    ########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig T")
    
    c1.cd(1).SetLeftMargin(0.12)
    nsep.Draw("sigt:t:sigt_e", "", "goff")
    #g_sigt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), [0] * nsep.GetSelectedRows(), nsep.GetV3())
    g_sigt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

    f_sigT_pre = TF2("sig_T_pre", fun_Sig_T, 0, 0.5, 0.0, 2.0, 4)
    f_sigT_pre.SetParameters(t0, t1, t2, t3)

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue

        sigt_X_pre = (f_sigT_pre.Eval(g_sigt.GetX()[i], q2_vec[i])) * g_vec[i]
        g_sigt_prv.SetPoint(i, g_sigt.GetX()[i], sigt_X_pre)

        sigt_X_fit = (g_sigt.GetY()[i]) / g_vec[i]
        sigt_X_fit_err = g_sigt.GetEY()[i] / g_vec[i]

        g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)
    
    g_max = g_sigt.GetYaxis().GetXmax()
    gp_max = max(range(g_sigt_prv.GetN()), key=lambda i: g_sigt_prv.GetY()[i])
    g_min = g_sigt.GetYaxis().GetXmin()
    gp_min = min(range(g_sigt_prv.GetN()), key=lambda i: g_sigt_prv.GetY()[i])

    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_sigt.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_sigt.SetMinimum(gp_min - difff)

    g_sigt.SetTitle("Sig T")

    g_sigt.SetMarkerStyle(5)
    g_sigt.Draw("AP")

    g_sigt.SetMaximum(hi_bound)
    g_sigt.SetMinimum(lo_bound)
    
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

    f_sigT = TF2("sig_T", fun_Sig_T, 0, 0.5, 0.0, 2.0, 4)
    f_sigT.SetParameters(t0, t1, t2, t3)

    fit_t_result = g_sigt_fit.Fit(f_sigT, "S")

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue        

        sigt_X = (f_sigT.Eval(g_sigt.GetX()[i], q2_vec[i])) * g_vec[i]
        g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

    gMinuit = ROOT.TMinuit()
        
    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)

    c1.cd(1)

    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    t0 = f_sigT.GetParameter(0)
    t1 = f_sigT.GetParameter(1)
    t2 = f_sigT.GetParameter(2)
    t3 = f_sigT.GetParameter(3)
    
    par_vec.append(t0)
    par_vec.append(t1)
    par_vec.append(t2)
    par_vec.append(t3)

    par_err_vec.append(f_sigT.GetParError(0))
    par_err_vec.append(f_sigT.GetParError(1))
    par_err_vec.append(f_sigT.GetParError(2))
    par_err_vec.append(f_sigT.GetParError(3))

    par_chi2_vec.append(f_sigT.GetChisquare())
    par_chi2_vec.append(f_sigT.GetChisquare())
    par_chi2_vec.append(f_sigT.GetChisquare())
    par_chi2_vec.append(f_sigT.GetChisquare())
    
    ########
    # SigL #
    ########
    
    print("/*--------------------------------------------------*/")
    print("Fit for Sig L")

    c1.cd(2).SetLeftMargin(0.12)
    nsep.Draw("sigl:t:sigl_e", "", "goff")

    f_sigL_pre = TF2("sig_L", fun_Sig_L, 0, 0.5, 0.0, 2.0, 4)
    f_sigL_pre.SetParameters(l0, l1, l2, l3)
    
    #g_sigl = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), 0, nsep.GetV3())
    g_sigl = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigl.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigl.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue
        
        sigl_X_pre = (f_sigL_pre.Eval(g_sigl.GetX()[i], q2_vec[i])) * g_vec[i]
        g_sigl_prv.SetPoint(i, g_sigl.GetX()[i], sigl_X_pre)

        sigl_X_fit = g_sigl.GetY()[i] / g_vec[i]
        sigl_X_fit_err = g_sigl.GetEY()[i] / g_vec[i]

        g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

    g_max = g_sigl.GetYaxis().GetXmax()
    gp_max = max(range(g_sigl_prv.GetN()), key=lambda i: g_sigl_prv.GetY()[i])
    g_min = g_sigl.GetYaxis().GetXmin()
    gp_min = min(range(g_sigl_prv.GetN()), key=lambda i: g_sigl_prv.GetY()[i])
    
    difff = (g_max - g_min) / 5

    if g_max < gp_max:
        g_sigl.SetMaximum(gp_max + difff)

    if g_min > gp_min:
        g_sigl.SetMinimum(gp_min - difff)
        
    g_sigl.SetTitle("Sig L")

    g_sigl.SetMarkerStyle(5)
    g_sigl.Draw("AP")

    g_sigl.SetMaximum(hi_bound)
    g_sigl.SetMinimum(lo_bound)

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

    f_sigL = TF2("sig_L", fun_Sig_L, 0, 0.5, 0.0, 2.0, 4)
    f_sigL.SetParameters(l0, l1, l2, l3)
    g_sigl_fit.Fit(f_sigL)

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue
        
        sigl_X = (f_sigL.Eval(g_sigl.GetX()[i], q2_vec[i])) * g_vec[i]
        g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)

    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)
    
    c1.cd(2)

    g_sigl_fit_tot.SetMarkerStyle(26)
    g_sigl_fit_tot.SetMarkerColor(2)
    g_sigl_fit_tot.SetLineColor(2)
    g_sigl_fit_tot.Draw("LP")

    l0 = f_sigL.GetParameter(0)
    l1 = f_sigL.GetParameter(1)
    l2 = f_sigL.GetParameter(2)
    l3 = f_sigL.GetParameter(3)

    par_vec.append(l0)
    par_vec.append(l1)
    par_vec.append(l2)
    par_vec.append(l3)

    par_err_vec.append(f_sigL.GetParError(0))
    par_err_vec.append(f_sigL.GetParError(1))
    par_err_vec.append(f_sigL.GetParError(2))
    par_err_vec.append(f_sigL.GetParError(3))

    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    par_chi2_vec.append(f_sigL.GetChisquare())
    
    #########
    # SigLT #
    #########

    print("/*--------------------------------------------------*/")
    print("Fit for Sig LT")

    c1.cd(3).SetLeftMargin(0.12)
    nsep.Draw("siglt:t:siglt_e", "", "goff")

    f_sigLT_pre = TF2("sig_LT", fun_Sig_LT, 0, 0.5, 0.0, 2.0, 4)
    f_sigLT_pre.SetParameters(lt0, lt1, lt2, lt3)
    
    #g_siglt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), ROOT.nullptr, nsep.GetV3())
    g_siglt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_siglt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_siglt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue
        
        siglt_X_pre = (f_sigLT_pre.Eval(g_siglt.GetX()[i], q2_vec[i]) * math.sin(th_vec[i] * PI / 180)) * g_vec[i]
        g_siglt_prv.SetPoint(i, g_sigl.GetX()[i], siglt_X_pre)

        if th_vec[i] != 180:
            siglt_X_fit = g_siglt.GetY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180)
            siglt_X_fit_err = g_siglt.GetEY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180)
        else:
            siglt_X_fit = 0.0
            siglt_X_fit_err = g_siglt.GetEY()[i]

        g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

    g_max = g_siglt.GetYaxis().GetXmax()
    gp_max = max(range(g_siglt_prv.GetN()), key=lambda i: g_siglt_prv.GetY()[i])
    g_min = g_siglt.GetYaxis().GetXmin()
    gp_min = min(range(g_siglt_prv.GetN()), key=lambda i: g_siglt_prv.GetY()[i])
    
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

    g_siglt.SetMaximum(hi_bound)
    g_siglt.SetMinimum(lo_bound)
    
    g_siglt_prv.SetMarkerColor(4)
    g_siglt_prv.SetMarkerStyle(21)
    g_siglt_prv.Draw("P")

    c2.cd(3)
    g_siglt_fit.SetTitle("Sigma LT Model Fit")
    g_siglt_fit.Draw("A*")

    f_sigLT = TF2("sig_LT", fun_Sig_LT, 0, 0.5, 0.0, 2.0, 4)
    f_sigLT.SetParameters(lt0, lt1, lt2, lt3)

    g_siglt_fit.Fit(f_sigLT)
    
    for i in range(0, len(w_vec)-1):        

        if q2_vec[i] == 0.0:
            continue
        
        siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i], q2_vec[i]) * math.sin(th_vec[i] * PI / 180)) * g_vec[i]
        g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)


    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)

    c1.cd(3)

    g_siglt_fit_tot.SetMarkerStyle(26)
    g_siglt_fit_tot.SetMarkerColor(2)
    g_siglt_fit_tot.SetLineColor(2)
    g_siglt_fit_tot.Draw("LP")

    lt0 = f_sigLT.GetParameter(0)
    lt1 = f_sigLT.GetParameter(1)
    lt2 = f_sigLT.GetParameter(2)
    lt3 = f_sigLT.GetParameter(3)
        
    par_vec.append(lt0)
    par_vec.append(lt1)
    par_vec.append(lt2)
    par_vec.append(lt3)

    par_err_vec.append(f_sigLT.GetParError(0))
    par_err_vec.append(f_sigLT.GetParError(1))
    par_err_vec.append(f_sigLT.GetParError(2))
    par_err_vec.append(f_sigLT.GetParError(3))

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
    
    f_sigTT_pre = TF2("sig_TT", fun_Sig_TT, 0, 0.5, 0.0, 2.0, 4)
    f_sigTT_pre.SetParameters(tt0, tt1, tt2, tt3)
    
    #g_sigtt = TGraphErrors(nsep.GetSelectedRows(), nsep.GetV2(), nsep.GetV1(), [0]*nsep.GetSelectedRows(), nsep.GetV3())
    g_sigtt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigtt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigtt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue
        
        sigtt_X_pre = (f_sigTT_pre.Eval(g_sigtt.GetX()[i], q2_vec[i]) * math.sin(th_vec[i] * PI / 180)**2) * g_vec[i]

        g_sigtt_prv.SetPoint(i, nsep.GetV2()[i], sigtt_X_pre)

        if th_vec[i] != 180:
            sigtt_X_fit = g_sigtt.GetY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180) / math.sin(th_vec[i] * PI / 180)
            sigtt_X_fit_err = g_sigtt.GetEY()[i] / g_vec[i] / math.sin(th_vec[i] * PI / 180) / math.sin(th_vec[i] * PI / 180)
        else:
            sigtt_X_fit = 0.0
            sigtt_X_fit_err = g_sigtt.GetEY()

        g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)
    
    g_max = g_sigtt.GetYaxis().GetXmax()
    gp_max = max(range(g_sigtt_prv.GetN()), key=lambda i: g_sigtt_prv.GetY()[i])
    g_min = g_sigtt.GetYaxis().GetXmin()
    gp_min = min(range(g_sigtt_prv.GetN()), key=lambda i: g_sigtt_prv.GetY()[i])

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

    g_sigtt.SetMaximum(hi_bound)
    g_sigtt.SetMinimum(lo_bound)
    
    g_sigtt_prv.SetMarkerColor(4)
    g_sigtt_prv.SetMarkerStyle(21)
    g_sigtt_prv.Draw("P")

    c2.cd(4)

    g_sigtt_fit.SetTitle("Sigma TT Model Fit")
    g_sigtt_fit.Draw("A*")

    f_sigTT = TF2("sig_TT", fun_Sig_TT, 0, 0.5, 0.0, 2.0, 4)
    f_sigTT.SetParameters(tt0, tt1, tt2, tt3)
    g_sigtt_fit.Fit(f_sigTT)
        
    for i in range(0, len(w_vec)-1):

        if q2_vec[i] == 0.0:
            continue
        
        sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i], q2_vec[i]) * math.sin(th_vec[i] * PI / 180)**2) * g_vec[i]
        g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)
        
    fit_status.DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit.fCstatu)

    c1.cd(4)

    g_sigtt_fit_tot.SetMarkerStyle(26)
    g_sigtt_fit_tot.SetMarkerColor(2)
    g_sigtt_fit_tot.SetLineColor(2)
    g_sigtt_fit_tot.Draw("LP")

    tt0 = f_sigTT.GetParameter(0)
    tt1 = f_sigTT.GetParameter(1)
    tt2 = f_sigTT.GetParameter(2)
    tt3 = f_sigTT.GetParameter(3)
    
    par_vec.append(tt0)
    par_vec.append(tt1)
    par_vec.append(tt2)
    par_vec.append(tt3)

    par_err_vec.append(f_sigTT.GetParError(0))
    par_err_vec.append(f_sigTT.GetParError(1))
    par_err_vec.append(f_sigTT.GetParError(2))
    par_err_vec.append(f_sigTT.GetParError(3))

    par_chi2_vec.append(f_sigTT.GetChisquare())
    par_chi2_vec.append(f_sigTT.GetChisquare())
    par_chi2_vec.append(f_sigTT.GetChisquare())
    par_chi2_vec.append(f_sigTT.GetChisquare())
    
    c1.Print(outputpdf+'(')
    c2.Print(outputpdf+')')

    para_file_out = "{}/src/{}/parameters/par.{}_{}.dat".format(LTANAPATH, ParticleType, pol_str, q2_set.replace("p",""))
    print("\nWriting {}...".format(para_file_out))
    with open(para_file_out, 'w') as f:
        format_specifier = "{:>13.5E} {:>13.5E} {:>3} {:>12.1f}"

        for i in range(len(par_vec)):
            sys.stdout.write(format_specifier.format(par_vec[i], par_err_vec[i], par_chi2_vec[i], i) + "\n")
            f.write(format_specifier.format(par_vec[i], par_err_vec[i], i, par_chi2_vec[i]) + "\n")
            print("  {}".format(par_vec[i], par_err_vec[i], i, par_chi2_vec[i]))
