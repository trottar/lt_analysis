#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-13 09:52:17 trottar"
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

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import TH1D_to_TH2D, create_polar_plot

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
# Disable statistics box by default
ROOT.gStyle.SetOptStat(0)
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
    data_charge_err_right = inpDict["data_charge_err_right"] 
    data_charge_err_left = inpDict["data_charge_err_left"] 
    data_charge_err_center = inpDict["data_charge_err_center"] 
    dummy_charge_err_right = inpDict["dummy_charge_err_right"] 
    dummy_charge_err_left = inpDict["dummy_charge_err_left"] 
    dummy_charge_err_center = inpDict["dummy_charge_err_center"]     
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"]
    InData_error_efficiency_right = inpDict["InData_error_efficiency_right"] 
    InData_error_efficiency_left = inpDict["InData_error_efficiency_left"] 
    InData_error_efficiency_center = inpDict["InData_error_efficiency_center"]         
    efficiency_table = inpDict["efficiency_table"] 
    ParticleType = inpDict["ParticleType"]

    # Step 5 uses the precomputed histograms from rand_sub.py/compare_simc.py.
    # These overlays intentionally happen before Step 6 yield extraction, so they
    # do not include the per-(t,phi) background treatment used to form the final
    # data/SIMC ratios.
    # Create an empty list to store copied histograms
    histlist_copy = []

    # Check if the value of the dictionary is a TObject that can use .Clone()
    for hist in histlist:
        hist_copy = {}
        for key, val in hist.items():
            if hasattr(val, 'Clone') and callable(getattr(val, 'Clone')):
                # Clone the TObject if it has the 'Clone' method
                hist_copy[key] = val.Clone()
            else:
                # Otherwise, just copy the value
                hist_copy[key] = val
        # Append the copied histogram dictionary to the new list
        histlist_copy.append(hist_copy)

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    
    # Histograms arriving here are already physics normalized upstream in
    # rand_sub.py / compare_simc.py. This plotting step should never
    # renormalize or area-match SIMC to data.
        
    # Creating clone of efficiency plots because otherwise things
    # crash due to pointer issues with TMultiGraph's Add() function
    data_eff_dict = {}
    for i,hist in enumerate(histlist_copy):
        data_eff_dict[hist["phi_setting"]] = hist["G_data_eff"].Clone()
    
    eff_plt = TCanvas()
    G_eff_plt = TMultiGraph()
    l_eff_plt = TLegend(0.115,0.35,0.33,0.5)

    eff_plt.SetGrid()

    for i,hist in enumerate(histlist_copy):
        data_eff_dict[hist["phi_setting"]].SetMarkerStyle(21)
        data_eff_dict[hist["phi_setting"]].SetMarkerSize(1)
        data_eff_dict[hist["phi_setting"]].SetMarkerColor(i+1)
        G_eff_plt.Add(data_eff_dict[hist["phi_setting"]])
        
    G_eff_plt.Draw("AP")

    G_eff_plt.SetTitle(" ;Run Numbers; Total Efficiency")

    i=0
    for i,hist in enumerate(histlist_copy):
        while i <= G_eff_plt.GetXaxis().GetXmax():
            bin_ix = G_eff_plt.GetXaxis().FindBin(i)
            if str(i) in hist["runNums"]: 
                G_eff_plt.GetXaxis().SetBinLabel(bin_ix,"%d" % i)
            i+=1

    G_eff_plt.GetYaxis().SetTitleOffset(1.5)
    G_eff_plt.GetXaxis().SetTitleOffset(1.5)
    G_eff_plt.GetXaxis().SetLabelSize(0.03)

    for i,hist in enumerate(histlist_copy):
        l_eff_plt.AddEntry(data_eff_dict[hist["phi_setting"]],hist["phi_setting"])

    l_eff_plt.Draw()

    eff_plt.Print(outputpdf + '(')

    # Plot histograms
    c_pid = TCanvas()

    c_pid.Divide(2,3)

    c_pid.cd(1)
    gPad.SetLogy()

    for i,hist in enumerate(histlist_copy):
        hist["H_cal_etottracknorm_DATA"].SetLineColor(i+1)
        hist["H_cal_etottracknorm_DATA"].Draw("same, E1")

    c_pid.cd(2)
    gPad.SetLogy()

    for i,hist in enumerate(histlist_copy):
        hist["H_cer_npeSum_DATA"].SetLineColor(i+1)
        hist["H_cer_npeSum_DATA"].Draw("same, E1")

    c_pid.cd(3)
    gPad.SetLogy()
    for i,hist in enumerate(histlist_copy):
        hist["P_cal_etottracknorm_DATA"].SetLineColor(i+1)
        hist["P_cal_etottracknorm_DATA"].Draw("same, E1")

    c_pid.cd(4)
    gPad.SetLogy()
    for i,hist in enumerate(histlist_copy):
        hist["P_hgcer_npeSum_DATA"].SetLineColor(i+1)
        hist["P_hgcer_npeSum_DATA"].Draw("same, E1")

    c_pid.cd(5)
    gPad.SetLogy()
    for i,hist in enumerate(histlist_copy):
        hist["P_aero_npeSum_DATA"].SetLineColor(i+1)
        hist["P_aero_npeSum_DATA"].Draw("same, E1")

    c_pid.Draw()

    c_pid.Print(outputpdf)

    ct = TCanvas()
    
    gPad.SetLogy()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_ct_DATA"].SetLineColor(i+1)
        hist["H_ct_DATA"].Draw("same, E1")

    ct.Print(outputpdf)
    
    CQ2 = TCanvas()
    l_Q2 = TLegend(0.1, 0.75, 0.35, 0.95)
    
    for i,hist in enumerate(histlist_copy):
        hist["H_Q2_DATA"].SetLineColor(i+1)
        hist["H_Q2_DATA"].Draw("same, E1")
        hist["H_Q2_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_Q2_SIMC"].Draw("same, E1")
        hist["H_Q2_SIMC"].Draw("same, HIST")
        l_Q2.AddEntry(hist["H_Q2_DATA"],hist["phi_setting"]+" Data")
        l_Q2.AddEntry(hist["H_Q2_SIMC"],hist["phi_setting"]+" Simc")

    l_Q2.Draw()
    CQ2.Print(outputpdf)

    CW = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_W_DATA"].SetLineColor(i+1)
        hist["H_W_DATA"].Draw("same, E1")
        hist["H_W_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_W_SIMC"].Draw("same, E1")
        hist["H_W_SIMC"].Draw("same, HIST")

    CW.Print(outputpdf)

    Ct = TCanvas()
    l_t = TLegend(0.1, 0.75, 0.35, 0.95)
    l_t.SetTextSize(0.0135)
    
    binmax = []
    for i,hist in enumerate(histlist_copy):
        hist["H_t_DATA"].SetLineColor(i+1)
        hist["H_t_DATA"].Draw("same, E1")
        hist["H_t_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_t_SIMC"].Draw("same, E1")
        hist["H_t_SIMC"].Draw("same, HIST")
        binmax.append(hist["H_t_DATA"].GetMaximum())
    binmax = max(binmax)

    t_bins = np.sort(t_bins)
    tBin_line = TLine()
    yields_between = []
    for j in range(0, len(t_bins)-1):
        yields_between_tmp = []
        for i,hist in enumerate(histlist_copy):
            b = t_bins[j]
            # Find the bins corresponding to the given bin centers
            bin1 = hist["H_t_DATA"].FindBin(t_bins[j])
            bin2 = hist["H_t_DATA"].FindBin(t_bins[j+1])
            # Sum the already physics-normalized yield in the displayed t-range.
            yields_between_tmp.append(sum(hist["H_t_DATA"].GetBinContent(k) for k in range(bin1, bin2+1)))
        yields_between.append(sum(yields_between_tmp))

        
    for j in range(0, len(t_bins)-1):
        b = t_bins[j]
        if j == 0:
            tBin_line.SetLineColor(7)
        else:
            tBin_line.SetLineColor(7)
        l_t.AddEntry(tBin_line,"Yield in {:.2f}-{:.2f}: {:.3e}".format(t_bins[j],t_bins[j+1], yields_between[j]))
        tBin_line.SetLineWidth(4)
        tBin_line.DrawLine(b,0,b,binmax)
        
    b = t_bins[j+1]
    tBin_line.SetLineColor(7)
    tBin_line.SetLineWidth(4)
    tBin_line.DrawLine(b,0,b,binmax)

    l_t.Draw()    

    Ct.Print(outputpdf)

    Cepsilon = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_epsilon_DATA"].SetLineColor(i+1)
        hist["H_epsilon_DATA"].Draw("same, E1")
        hist["H_epsilon_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_epsilon_SIMC"].Draw("same, E1")
        hist["H_epsilon_SIMC"].Draw("same, HIST")

    Cepsilon.Print(outputpdf)

    CMM = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_MM_DATA"].SetLineColor(i+1)
        hist["H_MM_DATA"].Draw("same, E1")
        hist["H_MM_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_MM_SIMC"].Draw("same, E1")
        hist["H_MM_SIMC"].Draw("same, HIST")

    CMM.Print(outputpdf)

    xfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssxfp_DATA"].SetLineColor(i+1)
        hist["H_ssxfp_DATA"].Draw("same, E1")
        hist["H_ssxfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssxfp_SIMC"].Draw("same, E1")
        hist["H_ssxfp_SIMC"].Draw("same, HIST")
        
    xfp.Print(outputpdf)

    yfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssyfp_DATA"].SetLineColor(i+1)
        hist["H_ssyfp_DATA"].Draw("same, E1")
        hist["H_ssyfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssyfp_SIMC"].Draw("same, E1")
        hist["H_ssyfp_SIMC"].Draw("same, HIST")

    yfp.Print(outputpdf)

    xpfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssxpfp_DATA"].SetLineColor(i+1)
        hist["H_ssxpfp_DATA"].Draw("same, E1")
        hist["H_ssxpfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssxpfp_SIMC"].Draw("same, E1")
        hist["H_ssxpfp_SIMC"].Draw("same, HIST")

    xpfp.Print(outputpdf)

    ypfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssypfp_DATA"].SetLineColor(i+1)
        hist["H_ssypfp_DATA"].Draw("same, E1")
        hist["H_ssypfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssypfp_SIMC"].Draw("same, E1")
        hist["H_ssypfp_SIMC"].Draw("same, HIST")

    ypfp.Print(outputpdf)

    hxfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsxfp_DATA"].SetLineColor(i+1)
        hist["H_hsxfp_DATA"].Draw("same, E1")
        hist["H_hsxfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsxfp_SIMC"].Draw("same, E1")
        hist["H_hsxfp_SIMC"].Draw("same, HIST")

    hxfp.Print(outputpdf)

    hyfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsyfp_DATA"].SetLineColor(i+1)
        hist["H_hsyfp_DATA"].Draw("same, E1")
        hist["H_hsyfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsyfp_SIMC"].Draw("same, E1")
        hist["H_hsyfp_SIMC"].Draw("same, HIST")

    hyfp.Print(outputpdf)

    hxpfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsxpfp_DATA"].SetLineColor(i+1)
        hist["H_hsxpfp_DATA"].Draw("same, E1")
        hist["H_hsxpfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsxpfp_SIMC"].Draw("same, E1")
        hist["H_hsxpfp_SIMC"].Draw("same, HIST")

    hxpfp.Print(outputpdf)

    hypfp = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsypfp_DATA"].SetLineColor(i+1)
        hist["H_hsypfp_DATA"].Draw("same, E1")
        hist["H_hsypfp_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsypfp_SIMC"].Draw("same, E1")
        hist["H_hsypfp_SIMC"].Draw("same, HIST")

    hypfp.Print(outputpdf)

    xptar = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssxptar_DATA"].SetLineColor(i+1)
        hist["H_ssxptar_DATA"].Draw("same, E1")
        hist["H_ssxptar_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssxptar_SIMC"].Draw("same, E1")
        hist["H_ssxptar_SIMC"].Draw("same, HIST")

    xptar.Print(outputpdf)

    yptar = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssyptar_DATA"].SetLineColor(i+1)
        hist["H_ssyptar_DATA"].Draw("same, E1")
        hist["H_ssyptar_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssyptar_SIMC"].Draw("same, E1")
        hist["H_ssyptar_SIMC"].Draw("same, HIST")

    yptar.Print(outputpdf)

    hxptar = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsxptar_DATA"].SetLineColor(i+1)
        hist["H_hsxptar_DATA"].Draw("same, E1")
        hist["H_hsxptar_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsxptar_SIMC"].Draw("same, E1")
        hist["H_hsxptar_SIMC"].Draw("same, HIST")

    hxptar.Print(outputpdf)

    hyptar = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsyptar_DATA"].SetLineColor(i+1)
        hist["H_hsyptar_DATA"].Draw("same, E1")
        hist["H_hsyptar_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsyptar_SIMC"].Draw("same, E1")
        hist["H_hsyptar_SIMC"].Draw("same, HIST")

    hyptar.Print(outputpdf)

    Delta = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ssdelta_DATA"].SetLineColor(i+1)
        hist["H_ssdelta_DATA"].Draw("same, E1")
        hist["H_ssdelta_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ssdelta_SIMC"].Draw("same, E1")
        hist["H_ssdelta_SIMC"].Draw("same, HIST")

    Delta.Print(outputpdf)

    hDelta = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_hsdelta_DATA"].SetLineColor(i+1)
        hist["H_hsdelta_DATA"].Draw("same, E1")
        hist["H_hsdelta_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_hsdelta_SIMC"].Draw("same, E1")
        hist["H_hsdelta_SIMC"].Draw("same, HIST")

    hDelta.Print(outputpdf)

    C_ssxptar_ssyptar = []        
    for i,hist in enumerate(histlist_copy):
        C_ssxptar_ssyptar.append(TCanvas())
        C_ssxptar_ssyptar[i].Divide(1,2)
        C_ssxptar_ssyptar[i].cd(1)
        h_ssyptar_ssxptar_data = TH1D_to_TH2D(hist["H_ssxptar_DATA"], hist["H_ssyptar_DATA"], h2d_name="h_ssyptar_ssxptar_data", title=f"{phisetlist[i]} Data;SHMS xptar;SHMS yptar", n_bins=25)
        h_ssyptar_ssxptar_data.Draw("COLZ")
        C_ssxptar_ssyptar[i].cd(2)
        h_ssyptar_ssxptar_simc = TH1D_to_TH2D(hist["H_ssxptar_SIMC"], hist["H_ssyptar_SIMC"], h2d_name="h_ssyptar_ssxptar_simc", title=f"{phisetlist[i]} Simc;SHMS xptar;SHMS yptar", n_bins=25, z_min=1e-12, z_max=h_ssyptar_ssxptar_data.GetMaximum())
        h_ssyptar_ssxptar_simc.Draw("COLZ")
        C_ssxptar_ssyptar[i].Print(outputpdf)

    C_hsxptar_hsyptar = []        
    for i,hist in enumerate(histlist_copy):
        C_hsxptar_hsyptar.append(TCanvas())
        C_hsxptar_hsyptar[i].Divide(1,2)
        C_hsxptar_hsyptar[i].cd(1)
        h_hsyptar_hsxptar_data = TH1D_to_TH2D(hist["H_hsxptar_DATA"], hist["H_hsyptar_DATA"], h2d_name="h_hsyptar_hsxptar_data", title=f"{phisetlist[i]} Data;HMS xptar;HMS yptar", n_bins=25)
        h_hsyptar_hsxptar_data.Draw("COLZ")
        C_hsxptar_hsyptar[i].cd(2)
        h_hsyptar_hsxptar_simc = TH1D_to_TH2D(hist["H_hsxptar_SIMC"], hist["H_hsyptar_SIMC"], h2d_name="h_hsyptar_hsxptar_simc", title=f"{phisetlist[i]} Simc;HMS xptar;HMS yptar", n_bins=25, z_min=1e-12, z_max=h_hsyptar_hsxptar_data.GetMaximum())
        h_hsyptar_hsxptar_simc.Draw("COLZ")
        C_hsxptar_hsyptar[i].Print(outputpdf)        

    C_ssdelta_hsdelta = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_hsdelta.append(TCanvas())
        C_ssdelta_hsdelta[i].Divide(1,2)
        C_ssdelta_hsdelta[i].cd(1)
        h_hsdelta_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_hsdelta_DATA"], h2d_name="h_hsdelta_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;HMS Delta", n_bins=25)
        h_hsdelta_ssdelta_data.Draw("COLZ")
        C_ssdelta_hsdelta[i].cd(2)
        h_hsdelta_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_hsdelta_SIMC"], h2d_name="h_hsdelta_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;HMS Delta", n_bins=25, z_min=1e-12, z_max=h_hsdelta_ssdelta_data.GetMaximum())
        h_hsdelta_ssdelta_simc.Draw("COLZ")
        C_ssdelta_hsdelta[i].Print(outputpdf)
        
    C_ssdelta_ssxptar = []    
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_ssxptar.append(TCanvas())
        C_ssdelta_ssxptar[i].Divide(1,2)
        C_ssdelta_ssxptar[i].cd(1)
        h_ssxptar_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_ssxptar_DATA"], h2d_name="h_ssxptar_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS Delta;SHMS xptar", n_bins=25)
        h_ssxptar_ssdelta_data.Draw("COLZ")
        C_ssdelta_ssxptar[i].cd(2)
        h_ssxptar_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_ssxptar_SIMC"], h2d_name="h_ssxptar_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS Delta;SHMS xptar", n_bins=25, z_min=1e-12, z_max=h_ssxptar_ssdelta_data.GetMaximum())
        h_ssxptar_ssdelta_simc.Draw("COLZ")
        C_ssdelta_ssxptar[i].Print(outputpdf)

    C_hsdelta_hsxptar = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_hsxptar.append(TCanvas())
        C_hsdelta_hsxptar[i].Divide(1,2)
        C_hsdelta_hsxptar[i].cd(1)
        h_hsxptar_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_hsxptar_DATA"], h2d_name="h_hsxptar_hsdelta_data", title=f"{phisetlist[i]} Data;HMS Delta;HMS xptar", n_bins=25)
        h_hsxptar_hsdelta_data.Draw("COLZ")
        C_hsdelta_hsxptar[i].cd(2)
        h_hsxptar_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_hsxptar_SIMC"], h2d_name="h_hsxptar_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS Delta;HMS xptar", n_bins=25, z_min=1e-12, z_max=h_hsxptar_hsdelta_data.GetMaximum())
        h_hsxptar_hsdelta_simc.Draw("COLZ")
        C_hsdelta_hsxptar[i].Print(outputpdf)

    C_ssdelta_ssyptar = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_ssyptar.append(TCanvas())
        C_ssdelta_ssyptar[i].Divide(1,2)
        C_ssdelta_ssyptar[i].cd(1)
        h_ssyptar_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_ssyptar_DATA"], h2d_name="h_ssyptar_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS Delta;SHMS yptar", n_bins=25)
        h_ssyptar_ssdelta_data.Draw("COLZ")
        C_ssdelta_ssyptar[i].cd(2)
        h_ssyptar_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_ssyptar_SIMC"], h2d_name="h_ssyptar_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS Delta;SHMS yptar", n_bins=25, z_min=1e-12, z_max=h_ssyptar_ssdelta_data.GetMaximum())
        h_ssyptar_ssdelta_simc.Draw("COLZ")
        C_ssdelta_ssyptar[i].Print(outputpdf)

    C_hsdelta_hsyptar = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_hsyptar.append(TCanvas())
        C_hsdelta_hsyptar[i].Divide(1,2)
        C_hsdelta_hsyptar[i].cd(1)
        h_hsyptar_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_hsyptar_DATA"], h2d_name="h_hsyptar_hsdelta_data", title=f"{phisetlist[i]} Data;HMS Delta;HMS yptar", n_bins=25)
        h_hsyptar_hsdelta_data.Draw("COLZ")
        C_hsdelta_hsyptar[i].cd(2)
        h_hsyptar_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_hsyptar_SIMC"], h2d_name="h_hsyptar_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS Delta;HMS yptar", n_bins=25, z_min=1e-12, z_max=h_hsyptar_hsdelta_data.GetMaximum())
        h_hsyptar_hsdelta_simc.Draw("COLZ")
        C_hsdelta_hsyptar[i].Print(outputpdf)

    C_ssdelta_Q2 = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_Q2.append(TCanvas())
        C_ssdelta_Q2[i].Divide(1,2)
        C_ssdelta_Q2[i].cd(1)
        h_Q2_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_Q2_DATA"], h2d_name="h_Q2_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;Q2", n_bins=25)
        h_Q2_ssdelta_data.Draw("COLZ")
        C_ssdelta_Q2[i].cd(2)
        h_Q2_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_Q2_SIMC"], h2d_name="h_Q2_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;Q2", n_bins=25, z_min=1e-12, z_max=h_Q2_ssdelta_data.GetMaximum())
        h_Q2_ssdelta_simc.Draw("COLZ")
        C_ssdelta_Q2[i].Print(outputpdf)

    C_hsdelta_Q2 = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_Q2.append(TCanvas())
        C_hsdelta_Q2[i].Divide(1,2)
        C_hsdelta_Q2[i].cd(1)
        h_Q2_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_Q2_DATA"], h2d_name="h_Q2_hsdelta_data", title=f"{phisetlist[i]} Data;HMS delta;Q2", n_bins=25)
        h_Q2_hsdelta_data.Draw("COLZ")
        C_hsdelta_Q2[i].cd(2)
        h_Q2_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_Q2_SIMC"], h2d_name="h_Q2_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS delta;Q2", n_bins=25, z_min=1e-12, z_max=h_Q2_hsdelta_data.GetMaximum())
        h_Q2_hsdelta_simc.Draw("COLZ")
        C_hsdelta_Q2[i].Print(outputpdf)

    C_ssdelta_MM = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_MM.append(TCanvas())
        C_ssdelta_MM[i].Divide(1,2)
        C_ssdelta_MM[i].cd(1)
        h_MM_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_MM_DATA"], h2d_name="h_MM_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;MM", n_bins=25)
        h_MM_ssdelta_data.Draw("COLZ")
        C_ssdelta_MM[i].cd(2)
        h_MM_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_MM_SIMC"], h2d_name="h_MM_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;MM", n_bins=25, z_min=1e-12, z_max=h_MM_ssdelta_data.GetMaximum())
        h_MM_ssdelta_simc.Draw("COLZ")
        C_ssdelta_MM[i].Print(outputpdf)

    C_hsdelta_MM = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_MM.append(TCanvas())
        C_hsdelta_MM[i].Divide(1,2)
        C_hsdelta_MM[i].cd(1)
        h_MM_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_MM_DATA"], h2d_name="h_MM_hsdelta_data", title=f"{phisetlist[i]} Data;HMS delta;MM", n_bins=25)
        h_MM_hsdelta_data.Draw("COLZ")
        C_hsdelta_MM[i].cd(2)
        h_MM_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_MM_SIMC"], h2d_name="h_MM_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS delta;MM", n_bins=25, z_min=1e-12, z_max=h_MM_hsdelta_data.GetMaximum())
        h_MM_hsdelta_simc.Draw("COLZ")
        C_hsdelta_MM[i].Print(outputpdf)        
        
    C_ssdelta_W = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_W.append(TCanvas())
        C_ssdelta_W[i].Divide(1,2)
        C_ssdelta_W[i].cd(1)
        h_W_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_W_DATA"], h2d_name="h_W_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;W", n_bins=25)
        h_W_ssdelta_data.Draw("COLZ")
        C_ssdelta_W[i].cd(2)
        h_W_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_W_SIMC"], h2d_name="h_W_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;W", n_bins=25, z_min=1e-12, z_max=h_W_ssdelta_data.GetMaximum())
        h_W_ssdelta_simc.Draw("COLZ")
        C_ssdelta_W[i].Print(outputpdf)

    C_hsdelta_W = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_W.append(TCanvas())
        C_hsdelta_W[i].Divide(1,2)
        C_hsdelta_W[i].cd(1)
        h_W_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_W_DATA"], h2d_name="h_W_hsdelta_data", title=f"{phisetlist[i]} Data;HMS delta;W", n_bins=25)
        h_W_hsdelta_data.Draw("COLZ")
        C_hsdelta_W[i].cd(2)
        h_W_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_W_SIMC"], h2d_name="h_W_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS delta;W", n_bins=25, z_min=1e-12, z_max=h_W_hsdelta_data.GetMaximum())
        h_W_hsdelta_simc.Draw("COLZ")
        C_hsdelta_W[i].Print(outputpdf)
        
    C_ssdelta_t = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_t.append(TCanvas())
        C_ssdelta_t[i].Divide(1,2)
        C_ssdelta_t[i].cd(1)
        h_t_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_t_DATA"], h2d_name="h_t_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;|-t|", n_bins=25)
        h_t_ssdelta_data.Draw("COLZ")
        C_ssdelta_t[i].cd(2)
        h_t_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_t_SIMC"], h2d_name="h_t_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;|-t|", n_bins=25, z_min=1e-12, z_max=h_t_ssdelta_data.GetMaximum())
        h_t_ssdelta_simc.Draw("COLZ")
        C_ssdelta_t[i].Print(outputpdf)
        
    C_hsdelta_t = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_t.append(TCanvas())
        C_hsdelta_t[i].Divide(1,2)
        C_hsdelta_t[i].cd(1)
        h_t_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_t_DATA"], h2d_name="h_t_hsdelta_data", title=f"{phisetlist[i]} Data;HMS delta;|-t|", n_bins=25)
        h_t_hsdelta_data.Draw("COLZ")
        C_hsdelta_t[i].cd(2)
        h_t_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_t_SIMC"], h2d_name="h_t_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS delta;|-t|", n_bins=25, z_min=1e-12, z_max=h_t_hsdelta_data.GetMaximum())
        h_t_hsdelta_simc.Draw("COLZ")
        C_hsdelta_t[i].Print(outputpdf)

    C_ssdelta_ph_q = []        
    for i,hist in enumerate(histlist_copy):
        C_ssdelta_ph_q.append(TCanvas())
        C_ssdelta_ph_q[i].Divide(1,2)
        C_ssdelta_ph_q[i].cd(1)
        h_ph_q_ssdelta_data = TH1D_to_TH2D(hist["H_ssdelta_DATA"], hist["H_ph_q_DATA"], h2d_name="h_ph_q_ssdelta_data", title=f"{phisetlist[i]} Data;SHMS delta;ph_q", n_bins=25)
        h_ph_q_ssdelta_data.Draw("COLZ")
        C_ssdelta_ph_q[i].cd(2)
        h_ph_q_ssdelta_simc = TH1D_to_TH2D(hist["H_ssdelta_SIMC"], hist["H_ph_q_SIMC"], h2d_name="h_ph_q_ssdelta_simc", title=f"{phisetlist[i]} Simc;SHMS delta;ph_q", n_bins=25, z_min=1e-12, z_max=h_ph_q_ssdelta_data.GetMaximum())
        h_ph_q_ssdelta_simc.Draw("COLZ")
        C_ssdelta_ph_q[i].Print(outputpdf)

    C_hsdelta_ph_q = []        
    for i,hist in enumerate(histlist_copy):
        C_hsdelta_ph_q.append(TCanvas())
        C_hsdelta_ph_q[i].Divide(1,2)
        C_hsdelta_ph_q[i].cd(1)
        h_ph_q_hsdelta_data = TH1D_to_TH2D(hist["H_hsdelta_DATA"], hist["H_ph_q_DATA"], h2d_name="h_ph_q_hsdelta_data", title=f"{phisetlist[i]} Data;HMS delta;ph_q", n_bins=25)
        h_ph_q_hsdelta_data.Draw("COLZ")
        C_hsdelta_ph_q[i].cd(2)
        h_ph_q_hsdelta_simc = TH1D_to_TH2D(hist["H_hsdelta_SIMC"], hist["H_ph_q_SIMC"], h2d_name="h_ph_q_hsdelta_simc", title=f"{phisetlist[i]} Simc;HMS delta;ph_q", n_bins=25, z_min=1e-12, z_max=h_ph_q_hsdelta_data.GetMaximum())
        h_ph_q_hsdelta_simc.Draw("COLZ")
        C_hsdelta_ph_q[i].Print(outputpdf)
        
    Cph_q = TCanvas()

    binmax = []
    for i,hist in enumerate(histlist_copy):
        hist["H_ph_q_DATA"].SetLineColor(i+1)
        l_t.AddEntry(hist["H_ph_q_DATA"],hist["phi_setting"])
        hist["H_ph_q_DATA"].Draw("same, E1")
        hist["H_ph_q_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ph_q_SIMC"].Draw("same, E1")
        hist["H_ph_q_SIMC"].Draw("same, HIST")
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

    for i,hist in enumerate(histlist_copy):
        hist["H_th_q_DATA"].SetLineColor(i+1)
        hist["H_th_q_DATA"].Draw("same, E1")
        hist["H_th_q_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_th_q_SIMC"].Draw("same, E1")
        hist["H_th_q_SIMC"].Draw("same, HIST")

    Cth_q.Print(outputpdf)

    Cph_recoil = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_ph_recoil_DATA"].SetLineColor(i+1)
        hist["H_ph_recoil_DATA"].Draw("same, E1")

    Cph_recoil.Print(outputpdf)

    Cth_recoil = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_th_recoil_DATA"].SetLineColor(i+1)
        hist["H_th_recoil_DATA"].Draw("same, E1")

    Cth_recoil.Print(outputpdf)

    Cpmiss = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_pmiss_DATA"].SetLineColor(i+1)
        hist["H_pmiss_DATA"].Draw("same, E1")
        hist["H_pmiss_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_pmiss_SIMC"].Draw("same, E1")
        hist["H_pmiss_SIMC"].Draw("same, HIST")

    Cpmiss.Print(outputpdf)

    Cemiss = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_emiss_DATA"].SetLineColor(i+1)
        hist["H_emiss_DATA"].Draw("same, E1")
        hist["H_emiss_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_emiss_SIMC"].Draw("same, E1")
        hist["H_emiss_SIMC"].Draw("same, HIST")

    Cemiss.Print(outputpdf)

    Cpmiss_x = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_pmx_DATA"].SetLineColor(i+1)
        hist["H_pmx_DATA"].Draw("same, E1")
        hist["H_pmx_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_pmx_SIMC"].Draw("same, E1")
        hist["H_pmx_SIMC"].Draw("same, HIST")

    Cpmiss_x.Print(outputpdf)

    Cpmiss_y = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_pmy_DATA"].SetLineColor(i+1)
        hist["H_pmy_DATA"].Draw("same, E1")
        hist["H_pmy_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_pmy_SIMC"].Draw("same, E1")
        hist["H_pmy_SIMC"].Draw("same, HIST")

    Cpmiss_y.Print(outputpdf)

    Cpmiss_z = TCanvas()

    for i,hist in enumerate(histlist_copy):
        hist["H_pmz_DATA"].SetLineColor(i+1)
        hist["H_pmz_DATA"].Draw("same, E1")
        hist["H_pmz_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_pmz_SIMC"].Draw("same, E1")
        hist["H_pmz_SIMC"].Draw("same, HIST")

    Cpmiss_z.Print(outputpdf)

    Cmmct = TCanvas()

    Cmmct.Divide(2,2)

    gPad.SetLogy()
    
    for i,hist in enumerate(histlist_copy):
        Cmmct.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_CoinTime_DATA"].SetMinimum(1)
        hist["MM_vs_CoinTime_DATA"].SetLineColor(i+1)
        hist["MM_vs_CoinTime_DATA"].Draw("same, COLZ")
        hist["MM_vs_CoinTime_DATA"].SetTitle(phisetlist[i])

    Cmmct.Print(outputpdf)

    Cctbeta = TCanvas()

    Cctbeta.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cctbeta.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["CoinTime_vs_beta_DATA"].SetMinimum(1)
        hist["CoinTime_vs_beta_DATA"].SetLineColor(i+1)
        hist["CoinTime_vs_beta_DATA"].Draw("same, COLZ")
        hist["CoinTime_vs_beta_DATA"].SetTitle(phisetlist[i])

    Cctbeta.Print(outputpdf)

    Cmmbeta = TCanvas()

    Cmmbeta.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cmmbeta.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_beta_DATA"].SetMinimum(1)
        hist["MM_vs_beta_DATA"].SetLineColor(i+1)
        hist["MM_vs_beta_DATA"].Draw("same, COLZ")
        hist["MM_vs_beta_DATA"].SetTitle(phisetlist[i])

    Cmmbeta.Print(outputpdf)

    CmmH_cer = TCanvas()

    CmmH_cer.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        CmmH_cer.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_H_cer_DATA"].SetMinimum(1)
        hist["MM_vs_H_cer_DATA"].SetLineColor(i+1)
        hist["MM_vs_H_cer_DATA"].Draw("same, COLZ")
        hist["MM_vs_H_cer_DATA"].SetTitle(phisetlist[i])

    CmmH_cer.Print(outputpdf)

    CmmH_cal = TCanvas()

    CmmH_cal.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        CmmH_cal.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_H_cal_DATA"].SetMinimum(1)
        hist["MM_vs_H_cal_DATA"].SetLineColor(i+1)
        hist["MM_vs_H_cal_DATA"].Draw("same, COLZ")
        hist["MM_vs_H_cal_DATA"].SetTitle(phisetlist[i])

    CmmH_cal.Print(outputpdf)    

    CmmP_cal = TCanvas()

    CmmP_cal.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        CmmP_cal.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_P_cal_DATA"].SetMinimum(1)
        hist["MM_vs_P_cal_DATA"].SetLineColor(i+1)
        hist["MM_vs_P_cal_DATA"].Draw("same, COLZ")
        hist["MM_vs_P_cal_DATA"].SetTitle(phisetlist[i])

    CmmP_cal.Print(outputpdf)
    
    CmmP_hgcer = TCanvas()

    CmmP_hgcer.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        CmmP_hgcer.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_P_hgcer_DATA"].SetMinimum(1)
        hist["MM_vs_P_hgcer_DATA"].SetLineColor(i+1)
        hist["MM_vs_P_hgcer_DATA"].Draw("same, COLZ")
        hist["MM_vs_P_hgcer_DATA"].SetTitle(phisetlist[i])

    CmmP_hgcer.Print(outputpdf)

    CmmP_aero = TCanvas()

    CmmP_aero.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        CmmP_aero.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["MM_vs_P_aero_DATA"].SetMinimum(1)
        hist["MM_vs_P_aero_DATA"].SetLineColor(i+1)
        hist["MM_vs_P_aero_DATA"].Draw("same, COLZ")
        hist["MM_vs_P_aero_DATA"].SetTitle(phisetlist[i])

    CmmP_aero.Print(outputpdf)
    
    Cqw = TCanvas()

    Cqw.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cqw.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["Q2_vs_W_DATA"].SetMinimum(1)
        hist["Q2_vs_W_DATA"].SetLineColor(i+1)
        hist["Q2_vs_W_DATA"].Draw("same, COLZ")
        hist["Q2_vs_W_DATA"].SetTitle(phisetlist[i])

    Cqw.Print(outputpdf)

    Cq2t = TCanvas()

    Cq2t.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cq2t.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["Q2_vs_t_DATA"].SetMinimum(1)
        hist["Q2_vs_t_DATA"].SetLineColor(i+1)
        hist["Q2_vs_t_DATA"].Draw("same, COLZ")
        hist["Q2_vs_t_DATA"].SetTitle(phisetlist[i])

    Cq2t.Print(outputpdf)    
    
    Cwt = TCanvas()

    Cwt.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cwt.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["W_vs_t_DATA"].SetMinimum(1)
        hist["W_vs_t_DATA"].SetLineColor(i+1)
        hist["W_vs_t_DATA"].Draw("same, COLZ")
        hist["W_vs_t_DATA"].SetTitle(phisetlist[i])

    Cwt.Print(outputpdf)    

    Cepst = TCanvas()

    Cepst.Divide(2,2)

    for i,hist in enumerate(histlist_copy):
        Cepst.cd(i+1)

        # Set color representing zero events to transparent (alpha = 0)
        hist["EPS_vs_t_DATA"].SetMinimum(1)
        hist["EPS_vs_t_DATA"].SetLineColor(i+1)
        hist["EPS_vs_t_DATA"].Draw("same, COLZ")
        hist["EPS_vs_t_DATA"].SetTitle(phisetlist[i])

    Cepst.Print(outputpdf)

    CMMt = TCanvas("CMMt", "MM_vs_t_Data_Simc", 1200, max(400, 400 * len(histlist_copy)))

    CMMt.Divide(2, max(1, len(histlist_copy)))

    for i,hist in enumerate(histlist_copy):
        data_pad = (2 * i) + 1
        simc_pad = data_pad + 1

        CMMt.cd(data_pad)
        hist["MM_vs_t_DATA"].SetMinimum(1)
        hist["MM_vs_t_DATA"].SetLineColor(i+1)
        hist["MM_vs_t_DATA"].SetTitle(phisetlist[i] + " Data")
        hist["MM_vs_t_DATA"].Draw("COLZ")

        CMMt.cd(simc_pad)
        hist["MM_vs_t_SIMC"].SetMinimum(1)
        hist["MM_vs_t_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["MM_vs_t_SIMC"].SetTitle(phisetlist[i] + " Simc")
        hist["MM_vs_t_SIMC"].Draw("COLZ")

    CMMt.Print(outputpdf)
    
    # For SIMC plots
    Cpht_simc = ROOT.TCanvas()
    polar_plots_simc = []

    for i, hist in enumerate(histlist_copy):
        polar_plot = create_polar_plot(hist["polar_phiq_vs_t_SIMC"], marker_color=i+1)
        polar_plots_simc.append(polar_plot)

        if i == 0:
            polar_plot.Draw("AP")
        else:
            polar_plot.Draw("P same")

        # Set titles and axes for the last plot
        polar_plots_simc[-1].GetXaxis().SetName("#Phi")
        polar_plots_simc[-1].GetYaxis().SetName("-t")
        polar_plots_simc[-1].SetTitle("") # SIMC            
            
    Cpht_simc.Print(outputpdf)

    # For DATA plots
    Cpht = ROOT.TCanvas()
    polar_plots_data = []

    for i, hist in enumerate(histlist_copy):
        polar_plot = create_polar_plot(hist["polar_phiq_vs_t_DATA"], marker_color=i+1)
        polar_plots_data.append(polar_plot)

        if i == 0:
            polar_plot.Draw("AP")
        else:
            polar_plot.Draw("P same")

        # Set titles and axes for the last plot
        polar_plots_data[-1].GetXaxis().SetName("#Phi")
        polar_plots_data[-1].GetYaxis().SetName("-t")
        polar_plots_data[-1].SetTitle("") # DATA            
            
    Cpht.Print(outputpdf)

    cut_summary_lst = ""
    for i,hist in enumerate(histlist_copy):
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
        cut_summary_lst += "\n\nCut Polygon Parameters"
        cut_mode = inpDict.get("cut_mode", "unknown")
        tex = TLatex(0., 0.+(0.95-(0.3+(0.05*(j+4)/2))), "cut_mode = {}".format(cut_mode))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)
        cut_summary_lst += "\ncut_mode = {}".format(cut_mode)

        poly_points = inpDict.get("poly_points", [])
        tex = TLatex(0., 0.+(0.95-(0.3+(0.05*(j+5)/2))), "vertices = {}".format(len(poly_points)))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)
        cut_summary_lst += "\nvertices = {}".format(len(poly_points))

        for p_idx, pval in enumerate(poly_points, start=1):
            tex = TLatex(
                0.,
                0.+(0.95-(0.3+(0.05*(j+5+p_idx)/2))),
                "v{} = ({:.6f}, {:.6f})".format(p_idx, float(pval[0]), float(pval[1])),
            )
            tex.SetTextSize(0.03)
            tex.SetTextColor(i+1)
            texlist.append(tex)
            cut_summary_lst += "\nv{} = ({:.6f}, {:.6f})".format(p_idx, float(pval[0]), float(pval[1]))

        for j, tex in enumerate(texlist):
            tex.Draw()

        if i == len(histlist_copy)-1:
            Ctext.Print(outputpdf+')')
        else:
            Ctext.Print(outputpdf)

    mm_shift_summary = inpDict.get("mm_shift_summary", {})
    if mm_shift_summary:
        cut_summary_lst += "\n\nMissing Mass Shifts"
        for phiset in phisetlist:
            shift_info = mm_shift_summary.get(phiset)
            if not shift_info:
                continue
            cut_summary_lst += "\n{} MM_shift = {:+.6f}".format(
                phiset,
                float(shift_info["shift"]),
            )

    t_shift_summary = inpDict.get("t_shift_summary", {})
    if t_shift_summary:
        cut_summary_lst += "\n\n-t Shifts"
        for phiset in phisetlist:
            shift_info = t_shift_summary.get(phiset)
            if not shift_info:
                continue
            cut_summary_lst += "\n{} t_shift = {:+.6f}".format(
                phiset,
                float(shift_info["shift"]),
            )
            
    return cut_summary_lst
