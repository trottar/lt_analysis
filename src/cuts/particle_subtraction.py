#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-01 16:39:59 trottar"
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
from ROOT import TH1D, TCutG, TFile
import sys, os, math

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

################################################################################################################################################

def particle_subtraction(subDict, inpDict, SubtractedParticle, hgcer_cutg=None, scale_factor=1.0):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = OUTPATH + "/" + "{}".format(SubtractedParticle) + "_" + InDATAFilename + "_%s.root" % (phi_setting)
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        histDict.update({ "phi_setting" : phi_setting})
        sys.exit(2)

    InFile_DATA = TFile.Open(rootFileData, "OPEN")

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(SubtractedParticle.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_RF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = OUTPATH + "/" + "{}".format(SubtractedParticle) + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        return histDict

    InFile_DUMMY = TFile.Open(rootFileDummy, "OPEN")  

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_RF".format(SubtractedParticle.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_RF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################

    H_MM_SUB_DATA = subDict["H_MM_SUB_DATA"]
    H_MM_SUB_DUMMY = subDict["H_MM_SUB_DUMMY"]
    H_MM_SUB_RAND = subDict["H_MM_SUB_RAND"]
    H_MM_SUB_DUMMY_RAND = subDict["H_MM_SUB_DUMMY_RAND"]
    
    ##############
    # HARD CODED #
    ##############

    mm_min = 0.7
    mm_max = 1.15
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

    for c0, p in zip(c0_list, h_momentum_list):
        if p == 0.889:
            c0_dict["Q2p1W2p95_lowe"] = c0 # Proper value 0.888
        elif p == 0.968:
            c0_dict["Q0p5W2p40_lowe"] = c0
            c0_dict["Q3p0W3p14_lowe"] = c0 # Proper value 1.821
            c0_dict["Q5p5W3p02_lowe"] = c0 # Proper value 0.962
        elif p == 2.185:
            c0_dict["Q0p5W2p40_highe"] = c0 # Proper value 2.066
            c0_dict["Q3p0W2p32_lowe"] = c0
        elif p == 2.328:
            c0_dict["Q4p4W2p74_lowe"] = c0
        elif p == 3.266:
            c0_dict["Q5p5W3p02_highe"] = c0            
        elif p == 4.2:
            c0_dict["Q3p0W3p14_highe"] = c0 # Proper value 4.204
        elif p == 4.712:
            c0_dict["Q4p4W2p74_highe"] = c0            
        elif p == 5.292:
            c0_dict["Q2p1W2p95_highe"] = c0
        elif p == 6.59:
            c0_dict["Q3p0W2p32_highe"] = c0
            
    ##############
    ##############        
    ##############

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
        else:
            ALLCUTS = apply_data_cuts(evt)
            
        if(ALLCUTS):
    
          #H_MM_SUB_DATA.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          #H_MM_SUB_DATA.Fill(evt.MM, evt.Weight)
          H_MM_SUB_DATA.Fill(evt.MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
        else:
            ALLCUTS = apply_data_cuts(evt)
            
        if(ALLCUTS):
    
          #H_MM_SUB_DUMMY.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          #H_MM_SUB_DUMMY.Fill(evt.MM, evt.Weight)
          H_MM_SUB_DUMMY.Fill(evt.MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} rand...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
        else:
            ALLCUTS = apply_data_cuts(evt)
            
        if(ALLCUTS):
    
          #H_MM_SUB_RAND.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          #H_MM_SUB_RAND.Fill(evt.MM, evt.Weight)
          H_MM_SUB_RAND.Fill(evt.MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} dummy_rand...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)
        else:
            ALLCUTS = apply_data_cuts(evt)
            
        if(ALLCUTS):
    
          #H_MM_SUB_DUMMY_RAND.Fill(np.sqrt(abs(pow(evt.Em, 2) - pow(evt.Pm, 2))), evt.Weight)
          #H_MM_SUB_DUMMY_RAND.Fill(evt.MM, evt.Weight)
          H_MM_SUB_DUMMY_RAND.Fill(evt.MM)


    H_MM_SUB_RAND.Scale(1/nWindows)
    H_MM_SUB_DUMMY_RAND.Scale(1/nWindows)
    
    H_MM_SUB_DATA.Add(H_MM_SUB_RAND,-1)
    H_MM_SUB_DUMMY.Add(H_MM_SUB_DUMMY_RAND,-1)

    H_MM_SUB_DATA.Add(H_MM_SUB_DUMMY,-1)
    
    H_MM_SUB_DATA.Scale(scale_factor)
