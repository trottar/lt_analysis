#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-23 19:10:14 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
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

def compare_simc(rootFileSimc, hist, inpDict):
    
    phi_setting = hist["phi_setting"]
    
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
    InData_error_efficiency_right = inpDict["InData_error_efficiency_right"]
    InData_error_efficiency_left = inpDict["InData_error_efficiency_left"]
    InData_error_efficiency_center = inpDict["InData_error_efficiency_center"]
    efficiency_table = inpDict["efficiency_table"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

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
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, simc=True)
    
    ################################################################################################################################################
    # Define simc root file trees of interest

    # Names don't match so need to do some string rearrangement
    if not os.path.isfile(rootFileSimc):
        print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
        return histDict

    # Opening new simc root file with new iteration of weight
    InFile_SIMC = TFile.Open(rootFileSimc, "OPEN")

    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    ###############################################################################################################################################

    # Grabs simc number of events and normalizaton factor
    simc_hist = rootFileSimc.replace('_iter.root','.hist')
    f_simc = open(simc_hist)
    for line in f_simc:
        #print(line)
        if "Ncontribute" in line:
            val = line.split("=")
            simc_nevents = int(val[1])
        if "normfac" in line:
            val = line.split("=")
            simc_normfactor = float(val[1])
    if 'simc_nevents' and 'simc_normfactor' not in locals():
        print("\n\nERROR: Invalid simc hist file %s\n\n" % simc_hist)
        sys.exit(1)
    f_simc.close()    

    ################################################################################################################################################
    # Plot definitions

    H_Weight_SIMC = TH1D("H_Weight_SIMC", "Simc Weight", 100, 0, 1e-5)
    H_iWeight_SIMC = TH1D("H_iWeight_SIMC", "Simc iWeight", 100, 0, 1e-5)
    H_hsdelta_SIMC  = TH1D("H_hsdelta_SIMC","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_SIMC  = TH1D("H_hsxptar_SIMC","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_SIMC  = TH1D("H_hsyptar_SIMC","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_SIMC    = TH1D("H_ssxfp_SIMC","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_SIMC    = TH1D("H_ssyfp_SIMC","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_SIMC   = TH1D("H_ssxpfp_SIMC","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_SIMC   = TH1D("H_ssypfp_SIMC","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_SIMC    = TH1D("H_hsxfp_SIMC","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_SIMC    = TH1D("H_hsyfp_SIMC","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_SIMC   = TH1D("H_hsxpfp_SIMC","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_SIMC   = TH1D("H_hsypfp_SIMC","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_SIMC  = TH1D("H_ssdelta_SIMC","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_SIMC  = TH1D("H_ssxptar_SIMC","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_SIMC  = TH1D("H_ssyptar_SIMC","SHMS yptar", 100, -0.04, 0.04)
    H_q_SIMC        = TH1D("H_q_SIMC","q", 100, 0.0, 10.0)
    H_Q2_SIMC       = TH1D("H_Q2_SIMC","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_SIMC  = TH1D("H_W_SIMC","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_SIMC       = TH1D("H_t_SIMC","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_SIMC  = TH1D("H_epsilon_SIMC","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_SIMC  = TH1D("H_MM_SIMC","MM_{K}", 100, 0.7, 1.5)
    H_MM_unweighted_SIMC  = TH1D("H_MM_unweighted_SIMC","MM_unweighted_{K}", 100, 0.7, 1.5)    
    H_th_SIMC  = TH1D("H_th_SIMC","X' tar", 100, -0.1, 0.1)
    H_ph_SIMC  = TH1D("H_ph_SIMC","Y' tar", 100, -0.1, 0.1)
    H_ph_q_SIMC  = TH1D("H_ph_q_SIMC","Phi Detected (ph_xq)", 100, 0.0, 2*math.pi)
    H_th_q_SIMC  = TH1D("H_th_q_SIMC","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_SIMC  = TH1D("H_ph_recoil_SIMC","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_SIMC  = TH1D("H_th_recoil_SIMC","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_SIMC  = TH1D("H_pmiss_SIMC","pmiss", 100, 0.0, 10.0)
    H_emiss_SIMC  = TH1D("H_emiss_SIMC","emiss", 100, 0.0, 10.0)
    H_pmx_SIMC  = TH1D("H_pmx_SIMC","pmx", 100, -10.0, 10.0)
    H_pmy_SIMC  = TH1D("H_pmy_SIMC","pmy ", 100, -10.0, 10.0)
    H_pmz_SIMC  = TH1D("H_pmz_SIMC","pmz", 100, -10.0, 10.0)

    polar_phiq_vs_t_SIMC = TGraphPolar()
    polar_phiq_vs_t_SIMC.SetName("polar_phiq_vs_t_SIMC")
    
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_SIMC = TH2D("P_hgcer_xAtCer_vs_yAtCer_SIMC", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    
    ################################################################################################################################################
    # Fill data histograms for various trees called above

    print("\nGrabbing %s simc..." % phi_setting)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

      # Define the acceptance cuts  
      SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
      HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

      Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

      if ParticleType == "kaon":
          
          ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.phgcer_x_det, evt.phgcer_y_det)
          NOHOLECUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond

      else:

          ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond
          
      if(NOHOLECUTS):
          # HGCer hole comparison            
          P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC.Fill(evt.phgcer_x_det,evt.phgcer_y_det)
          
      #Fill SIMC events
      if(ALLCUTS):

          # HGCer hole comparison
          P_hgcer_xAtCer_vs_yAtCer_SIMC.Fill(evt.phgcer_x_det,evt.phgcer_y_det)
      
          polar_phiq_vs_t_SIMC.SetPoint(polar_phiq_vs_t_SIMC.GetN(), (evt.phipq+math.pi)*(180/math.pi), -evt.t)
          
          H_Weight_SIMC.Fill(evt.Weight)
          H_iWeight_SIMC.Fill(evt.iter_weight)

          H_ssxfp_SIMC.Fill(evt.ssxfp, evt.iter_weight)
          H_ssyfp_SIMC.Fill(evt.ssyfp, evt.iter_weight)
          H_ssxpfp_SIMC.Fill(evt.ssxpfp, evt.iter_weight)
          H_ssypfp_SIMC.Fill(evt.ssypfp, evt.iter_weight)
          H_hsxfp_SIMC.Fill(evt.hsxfp, evt.iter_weight)
          H_hsyfp_SIMC.Fill(evt.hsyfp, evt.iter_weight)
          H_hsxpfp_SIMC.Fill(evt.hsxpfp, evt.iter_weight)
          H_hsypfp_SIMC.Fill(evt.hsypfp, evt.iter_weight)
          H_ssdelta_SIMC.Fill(evt.ssdelta, evt.iter_weight) 
          H_hsdelta_SIMC.Fill(evt.hsdelta, evt.iter_weight)	
          H_ssxptar_SIMC.Fill(evt.ssxptar, evt.iter_weight)
          H_ssyptar_SIMC.Fill(evt.ssyptar, evt.iter_weight)
          H_hsxptar_SIMC.Fill(evt.hsxptar, evt.iter_weight)	
          H_hsyptar_SIMC.Fill(evt.hsyptar, evt.iter_weight)

          H_ph_q_SIMC.Fill((evt.phipq+math.pi), evt.iter_weight)
          H_th_q_SIMC.Fill(evt.thetapq, evt.iter_weight)

          H_pmiss_SIMC.Fill(evt.Pm, evt.iter_weight)	
          H_emiss_SIMC.Fill(evt.Em, evt.iter_weight)	
          #H_pmx_SIMC.Fill(evt.Pmx, evt.iter_weight)
          #H_pmy_SIMC.Fill(evt.Pmy, evt.iter_weight)
          #H_pmz_SIMC.Fill(evt.Pmz, evt.iter_weight)
          H_Q2_SIMC.Fill(evt.Q2, evt.iter_weight)
          H_W_SIMC.Fill(evt.W, evt.iter_weight)
          H_t_SIMC.Fill(-evt.t, evt.iter_weight)
          H_epsilon_SIMC.Fill(evt.epsilon, evt.iter_weight)
          #H_MM_SIMC.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.iter_weight)
          H_MM_SIMC.Fill(evt.missmass, evt.iter_weight)
          H_MM_unweighted_SIMC.Fill(evt.missmass)
              
    ################################################################################################################################################    

    histDict["InFile_SIMC"] = InFile_SIMC
    histDict["simc_normfactor"] = simc_normfactor
    histDict["simc_nevents"] = simc_nevents
    histDict["H_Weight_SIMC"] =     H_Weight_SIMC
    histDict["H_iWeight_SIMC"] =     H_iWeight_SIMC
    histDict["H_hsdelta_SIMC"] =     H_hsdelta_SIMC
    histDict["H_hsxptar_SIMC"] =     H_hsxptar_SIMC
    histDict["H_hsyptar_SIMC"] =     H_hsyptar_SIMC
    histDict["H_ssxfp_SIMC"] =     H_ssxfp_SIMC  
    histDict["H_ssyfp_SIMC"] =     H_ssyfp_SIMC  
    histDict["H_ssxpfp_SIMC"] =     H_ssxpfp_SIMC 
    histDict["H_ssypfp_SIMC"] =     H_ssypfp_SIMC 
    histDict["H_hsxfp_SIMC"] =     H_hsxfp_SIMC  
    histDict["H_hsyfp_SIMC"] =     H_hsyfp_SIMC  
    histDict["H_hsxpfp_SIMC"] =     H_hsxpfp_SIMC 
    histDict["H_hsypfp_SIMC"] =     H_hsypfp_SIMC 
    histDict["H_ssdelta_SIMC"] =     H_ssdelta_SIMC
    histDict["H_ssxptar_SIMC"] =     H_ssxptar_SIMC
    histDict["H_ssyptar_SIMC"] =     H_ssyptar_SIMC
    histDict["H_q_SIMC"] =     H_q_SIMC      
    histDict["H_Q2_SIMC"] =     H_Q2_SIMC     
    histDict["H_t_SIMC"] =     H_t_SIMC     
    histDict["H_epsilon_SIMC"] =     H_epsilon_SIMC
    histDict["H_MM_SIMC"] =     H_MM_SIMC
    histDict["H_th_SIMC"] =     H_th_SIMC
    histDict["H_ph_SIMC"] =     H_ph_SIMC
    histDict["H_ph_q_SIMC"] =     H_ph_q_SIMC
    histDict["H_th_q_SIMC"] =     H_th_q_SIMC
    histDict["H_ph_recoil_SIMC"] =     H_ph_recoil_SIMC
    histDict["H_th_recoil_SIMC"] =     H_th_recoil_SIMC
    histDict["H_pmiss_SIMC"] =     H_pmiss_SIMC
    histDict["H_emiss_SIMC"] =     H_emiss_SIMC
    histDict["H_pmx_SIMC"] =     H_pmx_SIMC
    histDict["H_pmy_SIMC"] =     H_pmy_SIMC
    histDict["H_pmz_SIMC"] =     H_pmz_SIMC
    histDict["H_W_SIMC"] =     H_W_SIMC
    histDict["polar_phiq_vs_t_SIMC"] = polar_phiq_vs_t_SIMC
    histDict["NumEvts_MM_SIMC"] = H_MM_SIMC.Integral()
    histDict["NumEvts_MM_unweighted_SIMC"] = H_MM_unweighted_SIMC.Integral()
    
    ################################################################################################################################################

    ###
    # Q2 plots    
    CQ2 = TCanvas()

    histDict["H_Q2_SIMC"].SetLineColor(1)
    histDict["H_Q2_SIMC"].Draw("same, E1")

    CQ2.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType))+'(')
    
    ###
    # W plots    
    CW = TCanvas()

    histDict["H_W_SIMC"].SetLineColor(1)
    histDict["H_W_SIMC"].Draw("same, E1")

    CW.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType)))

    ###
    # MM plots    
    CMM = TCanvas()

    histDict["H_MM_SIMC"].SetLineColor(1)
    histDict["H_MM_SIMC"].Draw("same, E1")

    CMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType)))

    ###
    # t-Phi plots        
    Cpht_data = TCanvas()

    # Create a list to store all polar plots
    polar_plots = []

    polar_plot = TGraphPolar(histDict["polar_phiq_vs_t_SIMC"].GetN(), histDict["polar_phiq_vs_t_SIMC"].GetX(), histDict["polar_phiq_vs_t_SIMC"].GetY())
    polar_plot.SetMarkerColor(1)
    polar_plot.SetMarkerSize(0.5)
    polar_plot.SetMarkerStyle(20)
    polar_plots.append(polar_plot)  # Store the plot in the list
    polar_plot.Draw("AP same")

    # Set titles and axes for the last plot
    polar_plots[-1].GetXaxis().SetName("#Phi")
    polar_plots[-1].GetYaxis().SetName("-t")
    polar_plots[-1].SetTitle("")
    
    Cpht_data.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType)))

    ###
    # t plots            
    Ct = TCanvas()
    l_t = TLegend(0.115,0.45,0.33,0.95)
    l_t.SetTextSize(0.0135)

    histDict["H_t_SIMC"].SetLineColor(1)
    l_t.AddEntry(histDict["H_t_SIMC"],phi_setting)
    histDict["H_t_SIMC"].Draw("same, E1")

    Ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType)))

    ###
    # phi plots            
    Cphi = TCanvas()
    l_phi = TLegend(0.115,0.45,0.33,0.95)
    l_phi.SetTextSize(0.0135)
    
    histDict["H_ph_q_SIMC"].SetLineColor(1)
    l_phi.AddEntry(histDict["H_ph_q_SIMC"],phi_setting)
    histDict["H_ph_q_SIMC"].Draw("same, E1")    

    Cphi.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType)))
    
    ##
    # HGCer Hole Plots
    c_hgcer_hole = TCanvas()

    c_hgcer_hole.Divide(2,2)

    c_hgcer_hole.cd(1)
    P_hgcer_xAtCer_vs_yAtCer_SIMC.SetMinimum(1e-6) # Remove color of empty bins
    P_hgcer_xAtCer_vs_yAtCer_SIMC.Draw("colz")

    c_hgcer_hole.cd(2)
    P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC.SetMinimum(1e-6) # Remove color of empty bins
    P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC.Draw("colz")

    c_hgcer_hole.cd(3)
    P_hgcer_xAtCer_vs_yAtCer_SIMC.SetMinimum(1e-6) # Remove color of empty bins
    P_hgcer_xAtCer_vs_yAtCer_SIMC.Draw("colz")
    hgcer_cutg.SetLineColor(7)
    hgcer_cutg.Draw("same")

    c_hgcer_hole.cd(4)
    P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC.SetMinimum(1e-6) # Remove color of empty bins
    P_hgcer_nohole_xAtCer_vs_yAtCer_SIMC.Draw("colz")
    hgcer_cutg.SetLineColor(7)
    hgcer_cutg.Draw("same")
    
    c_hgcer_hole.Draw()    
    
    c_hgcer_hole.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(phi_setting,ParticleType))+')')
    
    return histDict
