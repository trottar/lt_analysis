#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-15 00:08:42 trottar"
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
from ROOT import TH1D
import sys, os, math

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

def particle_subtraction(inpDict, phi_setting, SubtractedParticle, scale_factor=1.0):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
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
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
    # Define simc root file trees of interest

    # Names don't match so need to do some string rearrangement
    #InSIMCFilename = "{}_Sub_Q{}W{}{}_{}e.root".format(SubtractedParticle, Q2, W, phi_setting.lower(), EPSSET)
    InSIMCFilename = "Prod_Coin_Q{}W{}{}_{}e.root".format(Q2, W, phi_setting.lower(), EPSSET)
    rootFileSimc = OUTPATH+"/"+InSIMCFilename
    if not os.path.isfile(rootFileSimc):
        print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
        sys.exit(2)

    InFile_SIMC = ROOT.TFile.Open(rootFileSimc, "OPEN")

    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    ################################################################################################################################################
    # Plot definitions

    H_MM_SUB_SIMC  = TH1D("H_MM_SUB_SIMC","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)

    ################################################################################################################################################    
    # Fill data histograms for various trees called above

    print("\nGrabbing {} {} subtraction simc...".format(phi_setting, SubtractedParticle))
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
          
          #H_MM_SUB_SIMC.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          H_MM_SUB_SIMC.Fill(evt.missmass, evt.Weight)    

    H_MM_SUB_SIMC.Scale(scale_factor)
        
    return H_MM_SUB_SIMC
