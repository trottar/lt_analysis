#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-03 16:48:20 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
from collections import defaultdict
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TLatex
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

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

################################################################################################################################################

def plot_data_vs_simc(t_bins, phi_bins, histlist, phisetlist, inpDict):

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]
    data_charge_right = inpDict["data_charge_right"] 
    data_charge_left = inpDict["data_charge_left"] 
    data_charge_center = inpDict["data_charge_center"] 
    dummy_charge_right = inpDict["dummy_charge_right"] 
    dummy_charge_left = inpDict["dummy_charge_left"] 
    dummy_charge_center = inpDict["dummy_charge_center"] 
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"] 
    efficiency_table = inpDict["efficiency_table"] 
    ParticleType = inpDict["ParticleType"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    
    eff_plt = TCanvas()
    G_eff_plt = TMultiGraph()
    l_eff_plt = TLegend(0.115,0.35,0.33,0.5)

    eff_plt.SetGrid()

    for i,hist in enumerate(histlist):
        hist["G_data_eff"].SetMarkerStyle(21)
        hist["G_data_eff"].SetMarkerSize(1)
        hist["G_data_eff"].SetMarkerColor(i+1)
        G_eff_plt.Add(hist["G_data_eff"])

    G_eff_plt.Draw("AP")

    G_eff_plt.SetTitle(" ;Run Numbers; Total Efficiency")

    i=0
    for i,hist in enumerate(histlist):
        while i <= G_eff_plt.GetXaxis().GetXmax():
            bin_ix = G_eff_plt.GetXaxis().FindBin(i)
            if str(i) in hist["runNums"]: 
                G_eff_plt.GetXaxis().SetBinLabel(bin_ix,"%d" % i)
            i+=1

    G_eff_plt.GetYaxis().SetTitleOffset(1.5)
    G_eff_plt.GetXaxis().SetTitleOffset(1.5)
    G_eff_plt.GetXaxis().SetLabelSize(0.04)

    for i,hist in enumerate(histlist):
        l_eff_plt.AddEntry(hist["G_data_eff"],hist["phi_setting"])

    l_eff_plt.Draw()

    eff_plt.Print(outputpdf + '(')

    # Plot histograms
    c_pid = TCanvas()

    c_pid.Divide(2,3)

    c_pid.cd(1)
    gPad.SetLogy()

    for i,hist in enumerate(histlist):
        hist["H_cal_etottracknorm_DATA"].SetLineColor(i+1)
        hist["H_cal_etottracknorm_DATA"].Draw("same, E1")

    c_pid.cd(2)
    gPad.SetLogy()

    for i,hist in enumerate(histlist):
        hist["H_cer_npeSum_DATA"].SetLineColor(i+1)
        hist["H_cer_npeSum_DATA"].Draw("same, E1")

    c_pid.cd(3)
    gPad.SetLogy()
    for i,hist in enumerate(histlist):
        hist["P_cal_etottracknorm_DATA"].SetLineColor(i+1)
        hist["P_cal_etottracknorm_DATA"].Draw("same, E1")

    c_pid.cd(4)
    gPad.SetLogy()
    for i,hist in enumerate(histlist):
        hist["P_hgcer_npeSum_DATA"].SetLineColor(i+1)
        hist["P_hgcer_npeSum_DATA"].Draw("same, E1")

    c_pid.cd(5)
    gPad.SetLogy()
    for i,hist in enumerate(histlist):
        hist["P_aero_npeSum_DATA"].SetLineColor(i+1)
        hist["P_aero_npeSum_DATA"].Draw("same, E1")

    c_pid.Draw()

    c_pid.Print(outputpdf)

    ct = TCanvas()
    
    gPad.SetLogy()
    
    for i,hist in enumerate(histlist):
        hist["H_ct_DATA"].SetLineColor(i+1)
        hist["H_ct_DATA"].Draw("same, E1")

    ct.Print(outputpdf)


    CQ2 = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_Q2_DATA"].SetLineColor(i+1)
        hist["H_Q2_DATA"].Draw("same, E1")
        hist["H_Q2_SIMC"].SetLineColor(40)
        hist["H_Q2_SIMC"].SetLineStyle(10-i)
        hist["H_Q2_SIMC"].Draw("same, E1")    

    CQ2.Print(outputpdf)

    CW = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_W_DATA"].SetLineColor(i+1)
        hist["H_W_DATA"].Draw("same, E1")
        hist["H_W_SIMC"].SetLineColor(40)
        hist["H_W_SIMC"].SetLineStyle(10-i)
        hist["H_W_SIMC"].Draw("same, E1")        

    CW.Print(outputpdf)

    Ct = TCanvas()
    l_t = TLegend(0.115,0.45,0.33,0.95)
    l_t.SetTextSize(0.0235)

    binmax = []
    for i,hist in enumerate(histlist):
        hist["H_t_DATA"].SetLineColor(i+1)
        l_t.AddEntry(hist["H_t_DATA"],hist["phi_setting"])
        hist["H_t_DATA"].Draw("same, E1")
        hist["H_t_SIMC"].SetLineColor(40)
        hist["H_t_SIMC"].SetLineStyle(10-i)
        hist["H_t_SIMC"].Draw("same, E1")
        binmax.append(hist["H_t_DATA"].GetMaximum())
    binmax = max(binmax)

    t_bins = np.sort(t_bins)
    tBin_line = TLine()
    events_between = []
    for j in range(0, len(t_bins)-1):
        events_between_tmp = []
        for i,hist in enumerate(histlist):
            b = t_bins[j]
            # Find the bins corresponding to the given bin centers
            bin1 = hist["H_t_DATA"].FindBin(t_bins[j])
            bin2 = hist["H_t_DATA"].FindBin(t_bins[j+1])
            # Get the content of the bins and calculate the number of events between them
            events_between_tmp.append(sum(hist["H_t_DATA"].GetBinContent(k) for k in range(bin1, bin2+1)))
        events_between.append(sum(events_between_tmp))

    for j in range(0, len(t_bins)-1):
        b = t_bins[j]
        if j == 0:
            tBin_line.SetLineColor(5)
        else:
            tBin_line.SetLineColor(4)
        l_t.AddEntry(tBin_line,"Evts in {:.2f}-{:.2f}: {:.0f}".format(t_bins[j],t_bins[j+1], events_between[j]))     
        tBin_line.SetLineWidth(5)
        tBin_line.DrawLine(b,0,b,binmax)
        
    b = t_bins[j+1]
    tBin_line.SetLineColor(5)
    tBin_line.SetLineWidth(4)
    tBin_line.DrawLine(b,0,b,binmax)
        
    l_t.Draw()    

    Ct.Print(outputpdf)

    Cepsilon = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_epsilon_DATA"].SetLineColor(i+1)
        hist["H_epsilon_DATA"].Draw("same, E1")
        hist["H_epsilon_SIMC"].SetLineColor(40)
        hist["H_epsilon_SIMC"].SetLineStyle(10-i)
        hist["H_epsilon_SIMC"].Draw("same, E1")

    Cepsilon.Print(outputpdf)

    CMM = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_MM_DATA"].SetLineColor(i+1)
        hist["H_MM_DATA"].Draw("same, E1")
        hist["H_MM_SIMC"].SetLineColor(40)
        hist["H_MM_SIMC"].SetLineStyle(10-i)
        hist["H_MM_SIMC"].Draw("same, E1")

    CMM.Print(outputpdf)

    xfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssxfp_DATA"].SetLineColor(i+1)
        hist["H_ssxfp_DATA"].Draw("same, E1")
        hist["H_ssxfp_SIMC"].SetLineColor(40)
        hist["H_ssxfp_SIMC"].SetLineStyle(10-i)
        hist["H_ssxfp_SIMC"].Draw("same, E1")

    xfp.Print(outputpdf)

    yfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssyfp_DATA"].SetLineColor(i+1)
        hist["H_ssyfp_DATA"].Draw("same, E1")
        hist["H_ssyfp_SIMC"].SetLineColor(40)
        hist["H_ssyfp_SIMC"].SetLineStyle(10-i)
        hist["H_ssyfp_SIMC"].Draw("same, E1")

    yfp.Print(outputpdf)

    xpfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssxpfp_DATA"].SetLineColor(i+1)
        hist["H_ssxpfp_DATA"].Draw("same, E1")
        hist["H_ssxpfp_SIMC"].SetLineColor(40)
        hist["H_ssxpfp_SIMC"].SetLineStyle(10-i)
        hist["H_ssxpfp_SIMC"].Draw("same, E1")

    xpfp.Print(outputpdf)

    ypfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssxpfp_DATA"].SetLineColor(i+1)
        hist["H_ssxpfp_DATA"].Draw("same, E1")
        hist["H_ssxpfp_SIMC"].SetLineColor(40)
        hist["H_ssxpfp_SIMC"].SetLineStyle(10-i)
        hist["H_ssxpfp_SIMC"].Draw("same, E1")

    ypfp.Print(outputpdf)

    hxfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsxfp_DATA"].SetLineColor(i+1)
        hist["H_hsxfp_DATA"].Draw("same, E1")
        hist["H_hsxfp_SIMC"].SetLineColor(40)
        hist["H_hsxfp_SIMC"].SetLineStyle(10-i)
        hist["H_hsxfp_SIMC"].Draw("same, E1")

    hxfp.Print(outputpdf)

    hyfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsyfp_DATA"].SetLineColor(i+1)
        hist["H_hsyfp_DATA"].Draw("same, E1")
        hist["H_hsyfp_SIMC"].SetLineColor(40)
        hist["H_hsyfp_SIMC"].SetLineStyle(10-i)
        hist["H_hsyfp_SIMC"].Draw("same, E1")

    hyfp.Print(outputpdf)

    hxpfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsxpfp_DATA"].SetLineColor(i+1)
        hist["H_hsxpfp_DATA"].Draw("same, E1")
        hist["H_hsxpfp_SIMC"].SetLineColor(40)
        hist["H_hsxpfp_SIMC"].SetLineStyle(10-i)
        hist["H_hsxpfp_SIMC"].Draw("same, E1")

    hxpfp.Print(outputpdf)

    hypfp = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsypfp_DATA"].SetLineColor(i+1)
        hist["H_hsypfp_DATA"].Draw("same, E1")
        hist["H_hsypfp_SIMC"].SetLineColor(40)
        hist["H_hsypfp_SIMC"].SetLineStyle(10-i)
        hist["H_hsypfp_SIMC"].Draw("same, E1")

    hypfp.Print(outputpdf)

    xptar = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssxptar_DATA"].SetLineColor(i+1)
        hist["H_ssxptar_DATA"].Draw("same, E1")
        hist["H_ssxptar_SIMC"].SetLineColor(40)
        hist["H_ssxptar_SIMC"].SetLineStyle(10-i)
        hist["H_ssxptar_SIMC"].Draw("same, E1")

    xptar.Print(outputpdf)

    yptar = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssyptar_DATA"].SetLineColor(i+1)
        hist["H_ssyptar_DATA"].Draw("same, E1")
        hist["H_ssyptar_SIMC"].SetLineColor(40)
        hist["H_ssyptar_SIMC"].SetLineStyle(10-i)
        hist["H_ssyptar_SIMC"].Draw("same, E1")

    yptar.Print(outputpdf)

    hxptar = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsxptar_DATA"].SetLineColor(i+1)
        hist["H_hsxptar_DATA"].Draw("same, E1")
        hist["H_hsxptar_SIMC"].SetLineColor(40)
        hist["H_hsxptar_SIMC"].SetLineStyle(10-i)
        hist["H_hsxptar_SIMC"].Draw("same, E1")

    hxptar.Print(outputpdf)

    hyptar = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsyptar_DATA"].SetLineColor(i+1)
        hist["H_hsyptar_DATA"].Draw("same, E1")
        hist["H_hsyptar_SIMC"].SetLineColor(40)
        hist["H_hsyptar_SIMC"].SetLineStyle(10-i)
        hist["H_hsyptar_SIMC"].Draw("same, E1")

    hyptar.Print(outputpdf)

    Delta = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ssdelta_DATA"].SetLineColor(i+1)
        hist["H_ssdelta_DATA"].Draw("same, E1")
        hist["H_ssdelta_SIMC"].SetLineColor(40)
        hist["H_ssdelta_SIMC"].SetLineStyle(10-i)
        hist["H_ssdelta_SIMC"].Draw("same, E1")

    Delta.Print(outputpdf)

    hDelta = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_hsdelta_DATA"].SetLineColor(i+1)
        hist["H_hsdelta_DATA"].Draw("same, E1")
        hist["H_hsdelta_SIMC"].SetLineColor(40)
        hist["H_hsdelta_SIMC"].SetLineStyle(10-i)
        hist["H_hsdelta_SIMC"].Draw("same, E1")

    hDelta.Print(outputpdf)

    Cph_q = TCanvas()

    binmax = []
    for i,hist in enumerate(histlist):
        hist["H_ph_q_DATA"].SetLineColor(i+1)
        l_t.AddEntry(hist["H_ph_q_DATA"],hist["phi_setting"])
        hist["H_ph_q_DATA"].Draw("same, E1")
        hist["H_ph_q_SIMC"].SetLineColor(40)
        hist["H_ph_q_SIMC"].SetLineStyle(10-i)
        hist["H_ph_q_SIMC"].Draw("same, E1")
        binmax.append(hist["H_ph_q_DATA"].GetMaximum())
    binmax = max(binmax)

    binned_phi_tmp = []
    for val in phi_bins:
        binned_phi_tmp.append(((float(val)/180))*math.pi)
    phiBin_line = TLine()
    for i,b in enumerate(binned_phi_tmp):
        b = float(b)
        phiBin_line.SetLineColor(4)
        phiBin_line.SetLineWidth(4)
        phiBin_line.DrawLine(b,0,b,binmax)
        l_t.AddEntry(phiBin_line,"Bin Edge %s" % i )
        l_t.AddEntry(phiBin_line,"BinCenter = %.2f" % b)

    Cph_q.Print(outputpdf)

    Cth_q = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_th_q_DATA"].SetLineColor(i+1)
        hist["H_th_q_DATA"].Draw("same, E1")
        hist["H_th_q_SIMC"].SetLineColor(40)
        hist["H_th_q_SIMC"].SetLineStyle(10-i)
        hist["H_th_q_SIMC"].Draw("same, E1")

    Cth_q.Print(outputpdf)

    Cph_recoil = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_ph_recoil_DATA"].SetLineColor(i+1)
        hist["H_ph_recoil_DATA"].Draw("same, E1")

    Cph_recoil.Print(outputpdf)

    Cth_recoil = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_th_recoil_DATA"].SetLineColor(i+1)
        hist["H_th_recoil_DATA"].Draw("same, E1")

    Cth_recoil.Print(outputpdf)

    Cpmiss = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_pmiss_DATA"].SetLineColor(i+1)
        hist["H_pmiss_DATA"].Draw("same, E1")
        hist["H_pmiss_SIMC"].SetLineColor(40)
        hist["H_pmiss_SIMC"].SetLineStyle(10-i)
        hist["H_pmiss_SIMC"].Draw("same, E1")

    Cpmiss.Print(outputpdf)

    Cemiss = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_emiss_DATA"].SetLineColor(i+1)
        hist["H_emiss_DATA"].Draw("same, E1")
        hist["H_emiss_SIMC"].SetLineColor(40)
        hist["H_emiss_SIMC"].SetLineStyle(10-i)
        hist["H_emiss_SIMC"].Draw("same, E1")

    Cemiss.Print(outputpdf)

    Cpmiss_x = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_pmx_DATA"].SetLineColor(i+1)
        hist["H_pmx_DATA"].Draw("same, E1")
        hist["H_pmx_SIMC"].SetLineColor(40)
        hist["H_pmx_SIMC"].SetLineStyle(10-i)
        hist["H_pmx_SIMC"].Draw("same, E1")

    Cpmiss_x.Print(outputpdf)

    Cpmiss_y = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_pmy_DATA"].SetLineColor(i+1)
        hist["H_pmy_DATA"].Draw("same, E1")
        hist["H_pmy_SIMC"].SetLineColor(40)
        hist["H_pmy_SIMC"].SetLineStyle(10-i)
        hist["H_pmy_SIMC"].Draw("same, E1")

    Cpmiss_y.Print(outputpdf)

    Cpmiss_z = TCanvas()

    for i,hist in enumerate(histlist):
        hist["H_pmz_DATA"].SetLineColor(i+1)
        hist["H_pmz_DATA"].Draw("same, E1")
        hist["H_pmz_SIMC"].SetLineColor(40)
        hist["H_pmz_SIMC"].SetLineStyle(10-i)
        hist["H_pmz_SIMC"].Draw("same, E1")

    Cpmiss_z.Print(outputpdf)

    Cmmct = TCanvas()

    Cmmct.Divide(2,2)

    gPad.SetLogy()
    
    for i,hist in enumerate(histlist):
        Cmmct.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_CoinTime_DATA"].SetMinimum(1)

        hist["MM_vs_CoinTime_DATA"].SetLineColor(i+1)
        hist["MM_vs_CoinTime_DATA"].Draw("same, COLZ")
        hist["MM_vs_CoinTime_DATA"].SetTitle(phisetlist[i])

    Cmmct.Print(outputpdf)

    Cctbeta = TCanvas()

    Cctbeta.Divide(2,2)

    for i,hist in enumerate(histlist):
        Cctbeta.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["CoinTime_vs_beta_DATA"].SetMinimum(1)

        hist["CoinTime_vs_beta_DATA"].SetLineColor(i+1)
        hist["CoinTime_vs_beta_DATA"].Draw("same, COLZ")
        hist["CoinTime_vs_beta_DATA"].SetTitle(phisetlist[i])

    Cctbeta.Print(outputpdf)

    Cmmbeta = TCanvas()

    Cmmbeta.Divide(2,2)

    for i,hist in enumerate(histlist):
        Cmmbeta.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_beta_DATA"].SetMinimum(1)

        hist["MM_vs_beta_DATA"].SetLineColor(i+1)
        hist["MM_vs_beta_DATA"].Draw("same, COLZ")
        hist["MM_vs_beta_DATA"].SetTitle(phisetlist[i])

    Cmmbeta.Print(outputpdf)

    Cqw = TCanvas()

    Cqw.Divide(2,2)

    for i,hist in enumerate(histlist):
        Cqw.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["Q2_vs_W_DATA"].SetMinimum(1)

        hist["Q2_vs_W_DATA"].SetLineColor(i+1)
        hist["Q2_vs_W_DATA"].Draw("same, COLZ")
        hist["Q2_vs_W_DATA"].SetTitle(phisetlist[i])

    Cqw.Print(outputpdf)

    Cpht_simc = TCanvas()

    # Create a list to store all polar plots
    polar_plots = []

    for i, hist in enumerate(histlist):
        polar_plot = TGraphPolar(hist["polar_phiq_vs_t_SIMC"].GetN(), hist["polar_phiq_vs_t_SIMC"].GetX(), hist["polar_phiq_vs_t_SIMC"].GetY())
        polar_plot.SetMarkerColor(i+1)
        polar_plot.SetMarkerSize(0.5)
        polar_plot.SetMarkerStyle(20)
        polar_plots.append(polar_plot)  # Store the plot in the list
        polar_plot.Draw("AP same")

    # Set titles and axes for the last plot
    polar_plots[-1].GetXaxis().SetName("#Phi")
    polar_plots[-1].GetYaxis().SetName("-t")
    polar_plots[-1].SetTitle("SIMC")
    
    Cpht_simc.Print(outputpdf)
    
    Cpht = TCanvas()

    # Create a list to store all polar plots
    polar_plots = []

    for i, hist in enumerate(histlist):
        polar_plot = TGraphPolar(hist["polar_phiq_vs_t_DATA"].GetN(), hist["polar_phiq_vs_t_DATA"].GetX(), hist["polar_phiq_vs_t_DATA"].GetY())
        polar_plot.SetMarkerColor(i+1)
        polar_plot.SetMarkerSize(0.5)
        polar_plot.SetMarkerStyle(20)
        polar_plots.append(polar_plot)  # Store the plot in the list
        polar_plot.Draw("AP same")

    # Set titles and axes for the last plot
    polar_plots[-1].GetXaxis().SetName("#Phi")
    polar_plots[-1].GetYaxis().SetName("-t")
    polar_plots[-1].SetTitle("DATA")
    
    Cpht.Print(outputpdf)

    cut_summary_lst = ""
    for i,hist in enumerate(histlist):
        texlist = []
        Ctext = TCanvas()
        for j,line in enumerate(hist["pid_text"]):
            if j == 0:
                tex = TLatex(0.8,0.+(0.95-(0.3)),"{}".format(hist["phi_setting"]))
                tex.SetTextSize(0.03)
                tex.SetTextColor(i+1)
                texlist.append(tex)
                cut_summary_lst += "\n"+hist["phi_setting"]
            tex = TLatex(0.,0.+(0.95-(0.3+(0.05*j/2))),"{}".format(line))
            tex.SetTextSize(0.03)
            tex.SetTextColor(i+1)
            texlist.append(tex)
            cut_summary_lst += "\n"+line
        j = len(hist["pid_text"])
        tex = TLatex(0.,0.+(0.95-(0.3+(0.05*(j+1)/2))),"t_range = ({}-{})".format(tmin,tmax))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)
        cut_summary_lst += "\n"+"t_range = ({}-{})".format(tmin,tmax)
        tex = TLatex(0.,0.+(0.95-(0.3+(0.05*(j+2)/2))),"t_bins-> {}".format(t_bins))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)
        cut_summary_lst += "\n"+"t_bins-> {}".format(t_bins)
        tex = TLatex(0.,0.+(0.95-(0.3+(0.05*(j+3)/2))),"phi_bins-> {}".format(phi_bins))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)
        cut_summary_lst += "\n"+"phi_bins-> {}".format(phi_bins)
        cut_summary_lst += "\n\n"+"Diamond Cut Paramters"
        for p in [1,2,3,4]:
            tex = TLatex(0.,0.+(0.95-(0.3+(0.05*(j+3+p)/2))),"a{} = {}, b{} = {}".format(p,inpDict["a%i" % p],p,inpDict["b%i" % p]))
            tex.SetTextSize(0.03)
            tex.SetTextColor(i+1)
            texlist.append(tex)
            cut_summary_lst += "\n"+"a{} = {}, b{} = {}".format(p,inpDict["a%i" % p],p,inpDict["b%i" % p])

        for j, tex in enumerate(texlist):
            tex.Draw()

        if i == len(histlist)-1:
            Ctext.Print(outputpdf+')')
        else:
            Ctext.Print(outputpdf)
            
    return cut_summary_lst
