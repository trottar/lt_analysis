#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-22 18:32:26 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
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
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
################################################################################################################################################

def plot_iteration(histlist, phisetlist, inpDict):

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    ParticleType = inpDict["ParticleType"]
    closest_date = inpDict["closest_date"]
    formatted_date = inpDict["formatted_date"]

    # Define diamond cut parameters
    a1 = inpDict["a1"]
    b1 = inpDict["b1"]
    a2 = inpDict["a2"]
    b2 = inpDict["b2"]
    a3 = inpDict["a3"]
    b3 = inpDict["b3"]
    a4 = inpDict["a4"]
    b4 = inpDict["b4"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################


    # Create an empty list to store copied histograms
    histlist_copy = []

    # Copy root tree and phi_setting strings to cloned dictionary
    for hist in histlist:
        hist_copy = {}
        hist_copy["InFile_SIMC"] = hist["InFile_SIMC"]
        hist_copy["phi_setting"] = hist["phi_setting"]
        # Append the copied histogram dictionary to the new list
        histlist_copy.append(hist_copy)
    
    for hist in histlist_copy:

        TBRANCH_SIMC  = hist["InFile_SIMC"].Get("h10")

        hist["H_Weight_SIMC"] = TH1D("H_Weight_SIMC", "Simc Weight", 100, 0, 1e-5)
        hist["H_hsdelta_SIMC"] = TH1D("H_hsdelta_SIMC","HMS Delta", 100, -20.0, 20.0)
        hist["H_hsxptar_SIMC"] = TH1D("H_hsxptar_SIMC","HMS xptar", 100, -0.1, 0.1)
        hist["H_hsyptar_SIMC"] = TH1D("H_hsyptar_SIMC","HMS yptar", 100, -0.1, 0.1)
        hist["H_ssxfp_SIMC"] = TH1D("H_ssxfp_SIMC","SHMS xfp", 100, -25.0, 25.0)
        hist["H_ssyfp_SIMC"] = TH1D("H_ssyfp_SIMC","SHMS yfp", 100, -25.0, 25.0)
        hist["H_ssxpfp_SIMC"] = TH1D("H_ssxpfp_SIMC","SHMS xpfp", 100, -0.09, 0.09)
        hist["H_ssypfp_SIMC"] = TH1D("H_ssypfp_SIMC","SHMS ypfp", 100, -0.05, 0.04)
        hist["H_hsxfp_SIMC"] = TH1D("H_hsxfp_SIMC","HMS xfp", 100, -40.0, 40.0)
        hist["H_hsyfp_SIMC"] = TH1D("H_hsyfp_SIMC","HMS yfp", 100, -20.0, 20.0)
        hist["H_hsxpfp_SIMC"] = TH1D("H_hsxpfp_SIMC","HMS xpfp", 100, -0.09, 0.05)
        hist["H_hsypfp_SIMC"] = TH1D("H_hsypfp_SIMC","HMS ypfp", 100, -0.05, 0.04)
        hist["H_ssdelta_SIMC"] = TH1D("H_ssdelta_SIMC","SHMS delta", 100, -20.0, 20.0)
        hist["H_ssxptar_SIMC"] = TH1D("H_ssxptar_SIMC","SHMS xptar", 100, -0.1, 0.1)
        hist["H_ssyptar_SIMC"] = TH1D("H_ssyptar_SIMC","SHMS yptar", 100, -0.04, 0.04)
        hist["H_q_SIMC"] = TH1D("H_q_SIMC","q", 100, 0.0, 10.0)
        hist["H_Q2_SIMC"] = TH1D("H_Q2_SIMC","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist["H_W_SIMC"] = TH1D("H_W_SIMC","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist["H_t_SIMC"] = TH1D("H_t_SIMC","-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist["H_epsilon_SIMC"] = TH1D("H_epsilon_SIMC","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist["H_MM_SIMC"] = TH1D("H_MM_SIMC","MM_{K}", 100, 0.7, 1.5)
        hist["H_th_SIMC"] = TH1D("H_th_SIMC","X' tar", 100, -0.1, 0.1)
        hist["H_ph_SIMC"] = TH1D("H_ph_SIMC","Y' tar", 100, -0.1, 0.1)
        hist["H_ph_q_SIMC"] = TH1D("H_ph_q_SIMC","Phi Detected (ph_xq)", 100, 0.0, 2*math.pi)
        hist["H_th_q_SIMC"] = TH1D("H_th_q_SIMC","Theta Detected (th_xq)", 100, -0.2, 0.2)
        hist["H_ph_recoil_SIMC"] = TH1D("H_ph_recoil_SIMC","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        hist["H_th_recoil_SIMC"] = TH1D("H_th_recoil_SIMC","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        hist["H_pmiss_SIMC"] = TH1D("H_pmiss_SIMC","pmiss", 100, 0.0, 10.0)
        hist["H_emiss_SIMC"] = TH1D("H_emiss_SIMC","emiss", 100, 0.0, 10.0)
        hist["H_pmx_SIMC"] = TH1D("H_pmx_SIMC","pmx", 100, -10.0, 10.0)
        hist["H_pmy_SIMC"] = TH1D("H_pmy_SIMC","pmy ", 100, -10.0, 10.0)
        hist["H_pmz_SIMC"] = TH1D("H_pmz_SIMC","pmz", 100, -10.0, 10.0)

        hist["H_Weight_SIMC_OLD"] = TH1D("H_Weight_SIMC_OLD", "Simc Weight", 100, 0, 1e-5)
        hist["H_hsdelta_SIMC_OLD"] = TH1D("H_hsdelta_SIMC_OLD","HMS Delta", 100, -20.0, 20.0)
        hist["H_hsxptar_SIMC_OLD"] = TH1D("H_hsxptar_SIMC_OLD","HMS xptar", 100, -0.1, 0.1)
        hist["H_hsyptar_SIMC_OLD"] = TH1D("H_hsyptar_SIMC_OLD","HMS yptar", 100, -0.1, 0.1)
        hist["H_ssxfp_SIMC_OLD"] = TH1D("H_ssxfp_SIMC_OLD","SHMS xfp", 100, -25.0, 25.0)
        hist["H_ssyfp_SIMC_OLD"] = TH1D("H_ssyfp_SIMC_OLD","SHMS yfp", 100, -25.0, 25.0)
        hist["H_ssxpfp_SIMC_OLD"] = TH1D("H_ssxpfp_SIMC_OLD","SHMS xpfp", 100, -0.09, 0.09)
        hist["H_ssypfp_SIMC_OLD"] = TH1D("H_ssypfp_SIMC_OLD","SHMS ypfp", 100, -0.05, 0.04)
        hist["H_hsxfp_SIMC_OLD"] = TH1D("H_hsxfp_SIMC_OLD","HMS xfp", 100, -40.0, 40.0)
        hist["H_hsyfp_SIMC_OLD"] = TH1D("H_hsyfp_SIMC_OLD","HMS yfp", 100, -20.0, 20.0)
        hist["H_hsxpfp_SIMC_OLD"] = TH1D("H_hsxpfp_SIMC_OLD","HMS xpfp", 100, -0.09, 0.05)
        hist["H_hsypfp_SIMC_OLD"] = TH1D("H_hsypfp_SIMC_OLD","HMS ypfp", 100, -0.05, 0.04)
        hist["H_ssdelta_SIMC_OLD"] = TH1D("H_ssdelta_SIMC_OLD","SHMS delta", 100, -20.0, 20.0)
        hist["H_ssxptar_SIMC_OLD"] = TH1D("H_ssxptar_SIMC_OLD","SHMS xptar", 100, -0.1, 0.1)
        hist["H_ssyptar_SIMC_OLD"] = TH1D("H_ssyptar_SIMC_OLD","SHMS yptar", 100, -0.04, 0.04)
        hist["H_q_SIMC_OLD"] = TH1D("H_q_SIMC_OLD","q", 100, 0.0, 10.0)
        hist["H_Q2_SIMC_OLD"] = TH1D("H_Q2_SIMC_OLD","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        hist["H_W_SIMC_OLD"] = TH1D("H_W_SIMC_OLD","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        hist["H_t_SIMC_OLD"] = TH1D("H_t_SIMC_OLD","-t", 100, inpDict["tmin"], inpDict["tmax"])
        hist["H_epsilon_SIMC_OLD"] = TH1D("H_epsilon_SIMC_OLD","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        hist["H_MM_SIMC_OLD"] = TH1D("H_MM_SIMC_OLD","MM_{K}", 100, 0.7, 1.5)
        hist["H_th_SIMC_OLD"] = TH1D("H_th_SIMC_OLD","X' tar", 100, -0.1, 0.1)
        hist["H_ph_SIMC_OLD"] = TH1D("H_ph_SIMC_OLD","Y' tar", 100, -0.1, 0.1)
        hist["H_ph_q_SIMC_OLD"] = TH1D("H_ph_q_SIMC_OLD","Phi Detected (ph_xq)", 100, 0.0, 2*math.pi)
        hist["H_th_q_SIMC_OLD"] = TH1D("H_th_q_SIMC_OLD","Theta Detected (th_xq)", 100, -0.2, 0.2)
        hist["H_ph_recoil_SIMC_OLD"] = TH1D("H_ph_recoil_SIMC_OLD","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        hist["H_th_recoil_SIMC_OLD"] = TH1D("H_th_recoil_SIMC_OLD","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        hist["H_pmiss_SIMC_OLD"] = TH1D("H_pmiss_SIMC_OLD","pmiss", 100, 0.0, 10.0)
        hist["H_emiss_SIMC_OLD"] = TH1D("H_emiss_SIMC_OLD","emiss", 100, 0.0, 10.0)
        hist["H_pmx_SIMC_OLD"] = TH1D("H_pmx_SIMC_OLD","pmx", 100, -10.0, 10.0)
        hist["H_pmy_SIMC_OLD"] = TH1D("H_pmy_SIMC_OLD","pmy ", 100, -10.0, 10.0)
        hist["H_pmz_SIMC_OLD"] = TH1D("H_pmz_SIMC_OLD","pmz", 100, -10.0, 10.0)        
        
        print("\n\nReweighing {} histograms for comparison...".format(hist["phi_setting"].lower()))
        for i,evt in enumerate(TBRANCH_SIMC):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

            # Define the acceptance cuts  
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

            Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

            #........................................

            #Fill SIMC events
            if(HMS_Acceptance & SHMS_Acceptance & Diamond):

                hist["H_Weight_SIMC"].Fill(evt.iter_weight)

                hist["H_ssxfp_SIMC"].Fill(evt.ssxfp, evt.iter_weight)
                hist["H_ssyfp_SIMC"].Fill(evt.ssyfp, evt.iter_weight)
                hist["H_ssxpfp_SIMC"].Fill(evt.ssxpfp, evt.iter_weight)
                hist["H_ssypfp_SIMC"].Fill(evt.ssypfp, evt.iter_weight)
                hist["H_hsxfp_SIMC"].Fill(evt.hsxfp, evt.iter_weight)
                hist["H_hsyfp_SIMC"].Fill(evt.hsyfp, evt.iter_weight)
                hist["H_hsxpfp_SIMC"].Fill(evt.hsxpfp, evt.iter_weight)
                hist["H_hsypfp_SIMC"].Fill(evt.hsypfp, evt.iter_weight)
                hist["H_ssdelta_SIMC"].Fill(evt.ssdelta, evt.iter_weight) 
                hist["H_hsdelta_SIMC"].Fill(evt.hsdelta, evt.iter_weight)	
                hist["H_ssxptar_SIMC"].Fill(evt.ssxptar, evt.iter_weight)
                hist["H_ssyptar_SIMC"].Fill(evt.ssyptar, evt.iter_weight)
                hist["H_hsxptar_SIMC"].Fill(evt.hsxptar, evt.iter_weight)	
                hist["H_hsyptar_SIMC"].Fill(evt.hsyptar, evt.iter_weight)

                hist["H_ph_q_SIMC"].Fill((evt.phipq+math.pi), evt.iter_weight)
                hist["H_th_q_SIMC"].Fill(evt.thetapq, evt.iter_weight)

                hist["H_pmiss_SIMC"].Fill(evt.Pm, evt.iter_weight)	
                hist["H_emiss_SIMC"].Fill(evt.Em, evt.iter_weight)	
                #hist["H_pmx_SIMC"].Fill(evt.Pmx, evt.iter_weight)
                #hist["H_pmy_SIMC"].Fill(evt.Pmy, evt.iter_weight)
                #hist["H_pmz_SIMC"].Fill(evt.Pmz, evt.iter_weight)
                hist["H_Q2_SIMC"].Fill(evt.Q2, evt.iter_weight)
                hist["H_W_SIMC"].Fill(evt.W, evt.iter_weight)
                hist["H_t_SIMC"].Fill(-evt.t, evt.iter_weight)
                hist["H_epsilon_SIMC"].Fill(evt.epsilon, evt.iter_weight)
                #hist["H_MM_SIMC"].Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.iter_weight)
                hist["H_MM_SIMC"].Fill(evt.missmass, evt.iter_weight)

                hist["H_Weight_SIMC_OLD"].Fill(evt.Weight)

                hist["H_ssxfp_SIMC_OLD"].Fill(evt.ssxfp, evt.Weight)
                hist["H_ssyfp_SIMC_OLD"].Fill(evt.ssyfp, evt.Weight)
                hist["H_ssxpfp_SIMC_OLD"].Fill(evt.ssxpfp, evt.Weight)
                hist["H_ssypfp_SIMC_OLD"].Fill(evt.ssypfp, evt.Weight)
                hist["H_hsxfp_SIMC_OLD"].Fill(evt.hsxfp, evt.Weight)
                hist["H_hsyfp_SIMC_OLD"].Fill(evt.hsyfp, evt.Weight)
                hist["H_hsxpfp_SIMC_OLD"].Fill(evt.hsxpfp, evt.Weight)
                hist["H_hsypfp_SIMC_OLD"].Fill(evt.hsypfp, evt.Weight)
                hist["H_ssdelta_SIMC_OLD"].Fill(evt.ssdelta, evt.Weight) 
                hist["H_hsdelta_SIMC_OLD"].Fill(evt.hsdelta, evt.Weight)	
                hist["H_ssxptar_SIMC_OLD"].Fill(evt.ssxptar, evt.Weight)
                hist["H_ssyptar_SIMC_OLD"].Fill(evt.ssyptar, evt.Weight)
                hist["H_hsxptar_SIMC_OLD"].Fill(evt.hsxptar, evt.Weight)	
                hist["H_hsyptar_SIMC_OLD"].Fill(evt.hsyptar, evt.Weight)

                hist["H_ph_q_SIMC_OLD"].Fill((evt.phipq+math.pi), evt.Weight)
                hist["H_th_q_SIMC_OLD"].Fill(evt.thetapq, evt.Weight)

                hist["H_pmiss_SIMC_OLD"].Fill(evt.Pm, evt.Weight)	
                hist["H_emiss_SIMC_OLD"].Fill(evt.Em, evt.Weight)	
                #hist["H_pmx_SIMC_OLD"].Fill(evt.Pmx, evt.Weight)
                #hist["H_pmy_SIMC_OLD"].Fill(evt.Pmy, evt.Weight)
                #hist["H_pmz_SIMC_OLD"].Fill(evt.Pmz, evt.Weight)
                hist["H_Q2_SIMC_OLD"].Fill(evt.Q2, evt.Weight)
                hist["H_W_SIMC_OLD"].Fill(evt.W, evt.Weight)
                hist["H_t_SIMC_OLD"].Fill(-evt.t, evt.Weight)
                hist["H_epsilon_SIMC_OLD"].Fill(evt.epsilon, evt.Weight)
                #hist["H_MM_SIMC_OLD"].Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
                hist["H_MM_SIMC_OLD"].Fill(evt.missmass, evt.Weight)
                
    CWeight = TCanvas()
    CWeight.Divide(2,2)
    l_Weight = TLegend(0.1, 0.75, 0.35, 0.95)
    
    for i,hist in enumerate(histlist_copy):
        CWeight.cd(i+1)
        hist["H_Weight_SIMC_OLD"].SetLineColor(i+1)
        hist["H_Weight_SIMC_OLD"].Draw("same, HIST")
        hist["H_Weight_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_Weight_SIMC"].Draw("same, HIST")
        l_Weight.AddEntry(hist["H_Weight_SIMC_OLD"],hist["phi_setting"]+" Simc Old")
        l_Weight.AddEntry(hist["H_Weight_SIMC"],hist["phi_setting"]+" Simc New")
    l_Weight.Draw()
        
    CWeight.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date))+'(')    
    CQ2 = TCanvas()
    CQ2.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        CQ2.cd(i+1)
        hist["H_Q2_SIMC_OLD"].SetLineColor(i+1)
        hist["H_Q2_SIMC_OLD"].Draw("same, HIST")
        hist["H_Q2_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_Q2_SIMC"].Draw("same, HIST")

    CQ2.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    CW = TCanvas()
    CW.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        CW.cd(i+1)
        hist["H_W_SIMC_OLD"].SetLineColor(i+1)
        hist["H_W_SIMC_OLD"].Draw("same, HIST")
        hist["H_W_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_W_SIMC"].Draw("same, HIST")

    CW.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    Ct = TCanvas()
    Ct.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        Ct.cd(i+1)
        hist["H_t_SIMC_OLD"].SetLineColor(i+1)
        hist["H_t_SIMC_OLD"].Draw("same, HIST")
        hist["H_t_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_t_SIMC"].Draw("same, HIST")

    Ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    Cepsilon = TCanvas()
    Cepsilon.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        Cepsilon.cd(i+1)
        hist["H_epsilon_SIMC_OLD"].SetLineColor(i+1)
        hist["H_epsilon_SIMC_OLD"].Draw("same, HIST")
        hist["H_epsilon_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_epsilon_SIMC"].Draw("same, HIST")

    Cepsilon.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    CMM = TCanvas()
    CMM.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        CMM.cd(i+1)
        hist["H_MM_SIMC_OLD"].SetLineColor(i+1)
        hist["H_MM_SIMC_OLD"].Draw("same, HIST")
        hist["H_MM_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_MM_SIMC"].Draw("same, HIST")

    CMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))
    
    Cph_q = TCanvas()
    Cph_q.Divide(2,2)
    
    for i,hist in enumerate(histlist_copy):
        Cph_q.cd(i+1)
        hist["H_ph_q_SIMC_OLD"].SetLineColor(i+1)
        hist["H_ph_q_SIMC_OLD"].Draw("same, HIST")
        hist["H_ph_q_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ph_q_SIMC"].Draw("same, HIST")

    Cph_q.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date))+')')
