#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-14 23:45:48 trottar"
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
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
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

################################################################################################################################################

def particle_subtraction(Q2, W, EPSSET, ParticleType, SubtractedParticle, scale_factor=1.0):

    ################################################################################################################################################
    # Define simc root file trees of interest

    # Names don't match so need to do some string rearrangement
    #InSIMCFilename = "{}_Sub_Q{}W{}_{}e.root".format(SubtractedParticle, Q2, W, EPSSET)
    InSIMCFilename = "Prod_Coin_{}.root".format(kinematics[0]+phi_setting.lower()+"_"+kinematics[1])
    rootFileSimc = OUTPATH+"/"+InSIMCFilename
    if not os.path.isfile(rootFileSimc):
        print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
        return histDict

    InFile_SIMC = ROOT.TFile.Open(rootFileSimc, "OPEN")

    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    ################################################################################################################################################
    # Plot definitions

    H_Weight_SUB_SIMC = TH1D("H_Weight_SUB_SIMC", "Sub_Simc Weight", 100, 0, 1e-5)
    H_hsdelta_SUB_SIMC  = TH1D("H_hsdelta_SUB_SIMC","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_SUB_SIMC  = TH1D("H_hsxptar_SUB_SIMC","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_SUB_SIMC  = TH1D("H_hsyptar_SUB_SIMC","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_SUB_SIMC    = TH1D("H_ssxfp_SUB_SIMC","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_SUB_SIMC    = TH1D("H_ssyfp_SUB_SIMC","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_SUB_SIMC   = TH1D("H_ssxpfp_SUB_SIMC","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_SUB_SIMC   = TH1D("H_ssypfp_SUB_SIMC","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_SUB_SIMC    = TH1D("H_hsxfp_SUB_SIMC","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_SUB_SIMC    = TH1D("H_hsyfp_SUB_SIMC","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_SUB_SIMC   = TH1D("H_hsxpfp_SUB_SIMC","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_SUB_SIMC   = TH1D("H_hsypfp_SUB_SIMC","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_SUB_SIMC  = TH1D("H_ssdelta_SUB_SIMC","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_SUB_SIMC  = TH1D("H_ssxptar_SUB_SIMC","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_SUB_SIMC  = TH1D("H_ssyptar_SUB_SIMC","SHMS yptar", 100, -0.04, 0.04)
    H_q_SUB_SIMC        = TH1D("H_q_SUB_SIMC","q", 100, 0.0, 10.0)
    H_Q2_SUB_SIMC       = TH1D("H_Q2_SUB_SIMC","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_SUB_SIMC  = TH1D("H_W_SUB_SIMC","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_SUB_SIMC       = TH1D("H_t_SUB_SIMC","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_SUB_SIMC  = TH1D("H_epsilon_SUB_SIMC","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_SUB_SIMC  = TH1D("H_MM_SUB_SIMC","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
    H_MM_unweighted_SUB_SIMC  = TH1D("H_MM_unweighted_SUB_SIMC","MM_unweighted_{}".format(SubtractedParticle), 100, 0.7, 1.5)
    H_th_SUB_SIMC  = TH1D("H_th_SUB_SIMC","X' tar", 100, -0.1, 0.1)
    H_ph_SUB_SIMC  = TH1D("H_ph_SUB_SIMC","Y' tar", 100, -0.1, 0.1)
    H_ph_q_SUB_SIMC  = TH1D("H_ph_q_SUB_SIMC","Phi Detected (ph_xq)", 100, 0.0, 2*math.pi)
    H_th_q_SUB_SIMC  = TH1D("H_th_q_SUB_SIMC","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_SUB_SIMC  = TH1D("H_ph_recoil_SUB_SIMC","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_SUB_SIMC  = TH1D("H_th_recoil_SUB_SIMC","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_SUB_SIMC  = TH1D("H_pmiss_SUB_SIMC","pmiss", 100, 0.0, 10.0)
    H_emiss_SUB_SIMC  = TH1D("H_emiss_SUB_SIMC","emiss", 100, 0.0, 10.0)
    H_pmx_SUB_SIMC  = TH1D("H_pmx_SUB_SIMC","pmx", 100, -10.0, 10.0)
    H_pmy_SUB_SIMC  = TH1D("H_pmy_SUB_SIMC","pmy ", 100, -10.0, 10.0)
    H_pmz_SUB_SIMC  = TH1D("H_pmz_SUB_SIMC","pmz", 100, -10.0, 10.0)

    ################################################################################################################################################    
    # Fill data histograms for various trees called above

    print("\nGrabbing %s subtraction simc..." % SubtractedParticle)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

      # Define the acceptance cuts  
      SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
      HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
      
      Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

      if ParticleType == "kaon":
          
          ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.phgcer_x_det, evt.phgcer_y_det)

      else:

          ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond
                    
      #Fill SIMC events
      if(ALLCUTS):
          
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

          # SIMC goes from 0 to 2pi so no need for +pi/2
          H_ph_q_SIMC.Fill((evt.phipq+math.pi), evt.Weight)
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
          #H_MM_SIMC.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          H_MM_SIMC.Fill(evt.missmass, evt.Weight)
          H_MM_unweighted_SIMC.Fill(evt.missmass)
    
    H_Weight_SUB_SIMC.Scale(scale_factor)
    H_hsdelta_SUB_SIMC.Scale(scale_factor)
    H_hsxptar_SUB_SIMC.Scale(scale_factor)
    H_hsyptar_SUB_SIMC.Scale(scale_factor)
    H_ssxfp_SUB_SIMC.Scale(scale_factor)
    H_ssyfp_SUB_SIMC.Scale(scale_factor)
    H_ssxpfp_SUB_SIMC.Scale(scale_factor)
    H_ssypfp_SUB_SIMC.Scale(scale_factor)
    H_hsxfp_SUB_SIMC.Scale(scale_factor)
    H_hsyfp_SUB_SIMC.Scale(scale_factor)
    H_hsxpfp_SUB_SIMC.Scale(scale_factor)
    H_hsypfp_SUB_SIMC.Scale(scale_factor)
    H_ssdelta_SUB_SIMC.Scale(scale_factor)
    H_ssxptar_SUB_SIMC.Scale(scale_factor)
    H_ssyptar_SUB_SIMC.Scale(scale_factor)
    H_q_SUB_SIMC.Scale(scale_factor)
    H_Q2_SUB_SIMC.Scale(scale_factor)
    H_W_SUB_SIMC.Scale(scale_factor)
    H_t_SUB_SIMC.Scale(scale_factor)
    H_epsilon_SUB_SIMC.Scale(scale_factor)
    H_MM_SUB_SIMC.Scale(scale_factor)
    H_MM_unweighted_SUB_SIMC.Scale(scale_factor)
    H_th_SUB_SIMC.Scale(scale_factor)
    H_ph_SUB_SIMC.Scale(scale_factor)
    H_ph_q_SUB_SIMC.Scale(scale_factor)
    H_th_q_SUB_SIMC.Scale(scale_factor)
    H_ph_recoil_SUB_SIMC.Scale(scale_factor)
    H_th_recoil_SUB_SIMC.Scale(scale_factor)
    H_pmiss_SUB_SIMC.Scale(scale_factor)
    H_emiss_SUB_SIMC.Scale(scale_factor)
    H_pmx_SUB_SIMC.Scale(scale_factor)
    H_pmy_SUB_SIMC.Scale(scale_factor)
    H_pmz_SUB_SIMC.Scale(scale_factor)
        
    return H_MM_SUB_SIMC
