#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-29 15:16:40 trottar"
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

def compare_simc(hist, inpDict):

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
    efficiency_table = inpDict["efficiency_table"] 
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
    # Define simc root file trees of interest

    # Names don't match so need to do some string rearrangement
    InSIMCFilename = "Prod_Coin_{}.root".format(kinematics[0]+phi_setting.lower()+"_"+kinematics[1])
    rootFileSimc = OUTPATH+"/"+InSIMCFilename
    if not os.path.isfile(rootFileSimc):
        print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
        return histDict

    InFile_SIMC = TFile.Open(rootFileSimc, "OPEN")

    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    ###############################################################################################################################################

    # Grabs simc number of events and normalizaton factor
    simc_hist = rootFileSimc.replace('.root','.hist')
    f_simc = open(simc_hist)
    for line in f_simc:
        #print(line)
        if "Ngen" in line:
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

    H_Weight_SIMC = TH1D("H_Weight_SIMC", "Simc Weight", 500, 0, 1e-8)
    H_hsdelta_SIMC  = TH1D("H_hsdelta_SIMC","HMS Delta", 500, -20.0, 20.0)
    H_hsxptar_SIMC  = TH1D("H_hsxptar_SIMC","HMS xptar", 500, -0.1, 0.1)
    H_hsyptar_SIMC  = TH1D("H_hsyptar_SIMC","HMS yptar", 500, -0.1, 0.1)
    H_ssxfp_SIMC    = TH1D("H_ssxfp_SIMC","SHMS xfp", 500, -25.0, 25.0)
    H_ssyfp_SIMC    = TH1D("H_ssyfp_SIMC","SHMS yfp", 500, -25.0, 25.0)
    H_ssxpfp_SIMC   = TH1D("H_ssxpfp_SIMC","SHMS xpfp", 500, -0.09, 0.09)
    H_ssypfp_SIMC   = TH1D("H_ssypfp_SIMC","SHMS ypfp", 500, -0.05, 0.04)
    H_hsxfp_SIMC    = TH1D("H_hsxfp_SIMC","HMS xfp", 500, -40.0, 40.0)
    H_hsyfp_SIMC    = TH1D("H_hsyfp_SIMC","HMS yfp", 500, -20.0, 20.0)
    H_hsxpfp_SIMC   = TH1D("H_hsxpfp_SIMC","HMS xpfp", 500, -0.09, 0.05)
    H_hsypfp_SIMC   = TH1D("H_hsypfp_SIMC","HMS ypfp", 500, -0.05, 0.04)
    H_ssdelta_SIMC  = TH1D("H_ssdelta_SIMC","SHMS delta", 500, -20.0, 20.0)
    H_ssxptar_SIMC  = TH1D("H_ssxptar_SIMC","SHMS xptar", 500, -0.1, 0.1)
    H_ssyptar_SIMC  = TH1D("H_ssyptar_SIMC","SHMS yptar", 500, -0.04, 0.04)
    H_q_SIMC        = TH1D("H_q_SIMC","q", 500, 0.0, 10.0)
    H_Q2_SIMC       = TH1D("H_Q2_SIMC","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_SIMC  = TH1D("H_W_SIMC","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_SIMC       = TH1D("H_t_SIMC","-t", 500, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_SIMC  = TH1D("H_epsilon_SIMC","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_SIMC  = TH1D("H_MM_SIMC","MM_{K}", 500, 0.0, 1.5)
    H_th_SIMC  = TH1D("H_th_SIMC","X' tar", 500, -0.1, 0.1)
    H_ph_SIMC  = TH1D("H_ph_SIMC","Y' tar", 500, -0.1, 0.1)
    H_ph_q_SIMC  = TH1D("H_ph_q_SIMC","Phi Detected (ph_xq)", 500, 0.0, 2*math.pi)
    H_th_q_SIMC  = TH1D("H_th_q_SIMC","Theta Detected (th_xq)", 500, -0.2, 0.2)
    H_ph_recoil_SIMC  = TH1D("H_ph_recoil_SIMC","Phi Recoil (ph_bq)", 500, -10.0, 10.0)
    H_th_recoil_SIMC  = TH1D("H_th_recoil_SIMC","Theta Recoil (th_bq)", 500, -10.0, 10.0)
    H_pmiss_SIMC  = TH1D("H_pmiss_SIMC","pmiss", 500, 0.0, 10.0)
    H_emiss_SIMC  = TH1D("H_emiss_SIMC","emiss", 500, 0.0, 10.0)
    H_pmx_SIMC  = TH1D("H_pmx_SIMC","pmx", 500, -10.0, 10.0)
    H_pmy_SIMC  = TH1D("H_pmy_SIMC","pmy ", 500, -10.0, 10.0)
    H_pmz_SIMC  = TH1D("H_pmz_SIMC","pmz", 500, -10.0, 10.0)

    polar_phiq_vs_t_SIMC = TGraphPolar()
    polar_phiq_vs_t_SIMC.SetName("polar_phiq_vs_t_SIMC")
    
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

      #........................................

      #Fill SIMC events
      if(HMS_Acceptance & SHMS_Acceptance & Diamond):

          # SIMC goes from 0 to 2pi so no need for +pi
          polar_phiq_vs_t_SIMC.SetPoint(polar_phiq_vs_t_SIMC.GetN(), (evt.phipq)*(180/math.pi), -evt.t)
          
          H_Weight_SIMC.Fill(evt.Weight)

          H_ssxfp_SIMC.Fill(evt.ssxfp, evt.Weight)
          H_ssyfp_SIMC.Fill(evt.ssyfp, evt.Weight)
          H_ssxpfp_SIMC.Fill(evt.ssxpfp, evt.Weight)
          H_ssypfp_SIMC.Fill(evt.ssypfp, evt.Weight)
          H_hsxfp_SIMC.Fill(evt.hsxfp, evt.Weight)
          H_hsyfp_SIMC.Fill(evt.hsyfp, evt.Weight)
          H_hsxpfp_SIMC.Fill(evt.hsxpfp, evt.Weight)
          H_hsypfp_SIMC.Fill(evt.hsypfp, evt.Weight)
          H_ssdelta_SIMC.Fill(evt.ssdelta, evt.Weight) 
          H_hsdelta_SIMC.Fill(evt.hsdelta, evt.Weight)	
          H_ssxptar_SIMC.Fill(evt.ssxptar, evt.Weight)
          H_ssyptar_SIMC.Fill(evt.ssyptar, evt.Weight)
          H_hsxptar_SIMC.Fill(evt.hsxptar, evt.Weight)	
          H_hsyptar_SIMC.Fill(evt.hsyptar, evt.Weight)

          # SIMC goes from 0 to 2pi, but data -pi to pi. So matching data with -pi
          H_ph_q_SIMC.Fill(evt.phipq, evt.Weight)
          H_th_q_SIMC.Fill(evt.thetapq, evt.Weight)

          H_pmiss_SIMC.Fill(evt.Pm, evt.Weight)	
          H_emiss_SIMC.Fill(evt.Em, evt.Weight)	
          #H_pmx_SIMC.Fill(evt.Pmx, evt.Weight)
          #H_pmy_SIMC.Fill(evt.Pmy, evt.Weight)
          #H_pmz_SIMC.Fill(evt.Pmz, evt.Weight)
          H_Q2_SIMC.Fill(evt.Q2, evt.Weight)
          H_W_SIMC.Fill(evt.W, evt.Weight)
          H_t_SIMC.Fill(-evt.t, evt.Weight)
          H_epsilon_SIMC.Fill(evt.epsilon, evt.Weight)
          H_MM_SIMC.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)

    ################################################################################################################################################
    # Normalize simc by normfactor/nevents

    normfac_simc = (simc_normfactor)/(simc_nevents)
              
    ################################################################################################################################################    

    histDict["InFile_SIMC"] = InFile_SIMC
    histDict["normfac_simc"] = normfac_simc
    histDict["H_Weight_SIMC"] =     H_Weight_SIMC
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
          
    ################################################################################################################################################

    #################
    # HARD CODED
    #################    
    #H_ssxfp_SIMC.Scale(10)
    #H_ssyfp_SIMC.Scale(10)
    #H_ssxpfp_SIMC.Scale(10)
    #H_ssypfp_SIMC.Scale(10)
    #H_hsxfp_SIMC.Scale(10)
    #H_hsyfp_SIMC.Scale(10)
    #H_hsxpfp_SIMC.Scale(10)
    #H_hsypfp_SIMC.Scale(10)
    #H_ssxptar_SIMC.Scale(10)
    #H_ssyptar_SIMC.Scale(10)
    #H_hsxptar_SIMC.Scale(10)
    #H_hsyptar_SIMC.Scale(10)
    #H_ssdelta_SIMC.Scale(10)
    #H_hsdelta_SIMC.Scale(10)
    #H_Q2_SIMC.Scale(10)
    #H_t_SIMC.Scale(10)
    #H_epsilon_SIMC.Scale(10)
    #H_MM_SIMC.Scale(10)
    #H_ph_q_SIMC.Scale(10)
    #H_th_q_SIMC.Scale(10)
    #H_ph_recoil_SIMC.Scale(10)
    #H_th_recoil_SIMC.Scale(10)
    #H_pmiss_SIMC.Scale(10)
    #H_emiss_SIMC.Scale(10)
    #H_W_SIMC.Scale(10)
    #################
    #################
    #################    
    
    return histDict
