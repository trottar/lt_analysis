#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-19 02:25:54 trottar"
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
from scipy.integrate import quad
import matplotlib.pyplot as plt
from collections import defaultdict
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
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
from utility import remove_bad_bins, get_centroid

##################################################################################################################################################

def process_hist_data(tree_data, tree_dummy, t_bins, nWindows, phi_setting, inpDict):

    processed_dict = {}

    OutFilename = inpDict["OutFilename"]    
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"]
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, set_val
    set_val(inpDict) # Set global variables for optimization

    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
        
    TBRANCH_DATA  = tree_data.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
    TBRANCH_RAND  = tree_data.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))
    
    TBRANCH_DUMMY  = tree_dummy.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
    TBRANCH_DUMMY_RAND  = tree_dummy.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))

    ##############
    # HARD CODED #
    ##############
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

    if ParticleType == "kaon":
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
    else:
        c0_dict["Q0p4W2p20_lowe"] = 0.0
        c0_dict["Q0p4W2p20_highe"] = 0.0
            
    ##############
    ##############        
    ##############

    hist_bin_dict = {}

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        from particle_subtraction import particle_subtraction_ave
        SubtractedParticle = "pion"
        subDict = {}

    # Fit background and subtract
    from background_fit import bg_fit
        
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):

        hist_bin_dict["H_Q2_DATA_{}".format(j)]       = TH1D("H_Q2_DATA_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DATA_{}".format(j)]  = TH1D("H_W_DATA_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DATA_{}".format(j)]       = TH1D("H_t_DATA_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DATA_{}".format(j)]  = TH1D("H_epsilon_DATA_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DATA_{}".format(j)]       = TH1D("H_MM_DATA_{}".format(j),"MM", 1000, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)]       = TH1D("H_MM_nosub_DATA_{}".format(j),"MM", 1000, 0.7, 1.5)

        hist_bin_dict["H_Q2_RAND_{}".format(j)]       = TH1D("H_Q2_RAND_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_RAND_{}".format(j)]  = TH1D("H_W_RAND_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_RAND_{}".format(j)]       = TH1D("H_t_RAND_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_RAND_{}".format(j)]  = TH1D("H_epsilon_RAND_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_RAND_{}".format(j)]       = TH1D("H_MM_RAND_{}".format(j),"MM", 1000, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_nosub_RAND_{}".format(j)]       = TH1D("H_MM_nosub_RAND_{}".format(j),"MM", 1000, 0.7, 1.5)

        hist_bin_dict["H_Q2_DUMMY_{}".format(j)]       = TH1D("H_Q2_DUMMY_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DUMMY_{}".format(j)]  = TH1D("H_W_DUMMY_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DUMMY_{}".format(j)]       = TH1D("H_t_DUMMY_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DUMMY_{}".format(j)]  = TH1D("H_epsilon_DUMMY_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DUMMY_{}".format(j)]       = TH1D("H_MM_DUMMY_{}".format(j),"MM", 1000, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)]       = TH1D("H_MM_nosub_DUMMY_{}".format(j),"MM", 1000, 0.7, 1.5)

        hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)]       = TH1D("H_Q2_DUMMY_RAND_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)]  = TH1D("H_W_DUMMY_RAND_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)]       = TH1D("H_t_DUMMY_RAND_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)]  = TH1D("H_epsilon_DUMMY_RAND_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
        hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_DUMMY_RAND_{}".format(j),"MM", 1000, inpDict["mm_min"], inpDict["mm_max"])
        hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)]       = TH1D("H_MM_nosub_DUMMY_RAND_{}".format(j),"MM", 1000, 0.7, 1.5)

        # Pion subtraction by scaling simc to peak size
        if ParticleType == "kaon":
            
            subDict["H_Q2_SUB_DATA_{}".format(j)]       = TH1D("H_Q2_SUB_DATA_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DATA_{}".format(j)]  = TH1D("H_W_SUB_DATA_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DATA_{}".format(j)]       = TH1D("H_t_SUB_DATA_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DATA_{}".format(j)]  = TH1D("H_epsilon_SUB_DATA_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DATA_{}".format(j)]  = TH1D("H_MM_SUB_DATA_{}".format(j),"MM_{}".format(SubtractedParticle), 1000, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DATA_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DATA_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 1000, 0.7, 1.5)

            subDict["H_Q2_SUB_RAND_{}".format(j)]       = TH1D("H_Q2_SUB_RAND_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_RAND_{}".format(j)]  = TH1D("H_W_SUB_RAND_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_RAND_{}".format(j)]       = TH1D("H_t_SUB_RAND_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_RAND_{}".format(j)]  = TH1D("H_epsilon_SUB_RAND_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_RAND_{}".format(j)]  = TH1D("H_MM_SUB_RAND_{}".format(j),"MM_{}".format(SubtractedParticle), 1000, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_RAND_{}".format(j)]  = TH1D("H_MM_nosub_SUB_RAND_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 1000, 0.7, 1.5)

            subDict["H_Q2_SUB_DUMMY_{}".format(j)]       = TH1D("H_Q2_SUB_DUMMY_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DUMMY_{}".format(j)]  = TH1D("H_W_SUB_DUMMY_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DUMMY_{}".format(j)]       = TH1D("H_t_SUB_DUMMY_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DUMMY_{}".format(j)]  = TH1D("H_epsilon_SUB_DUMMY_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DUMMY_{}".format(j)]  = TH1D("H_MM_SUB_DUMMY_{}".format(j),"MM_{}".format(SubtractedParticle), 1000, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DUMMY_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DUMMY_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 1000, 0.7, 1.5)

            subDict["H_Q2_SUB_DUMMY_RAND_{}".format(j)]       = TH1D("H_Q2_SUB_DUMMY_RAND_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
            subDict["H_W_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_W_SUB_DUMMY_RAND_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
            subDict["H_t_SUB_DUMMY_RAND_{}".format(j)]       = TH1D("H_t_SUB_DUMMY_RAND_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
            subDict["H_epsilon_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_epsilon_SUB_DUMMY_RAND_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])
            subDict["H_MM_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_MM_SUB_DUMMY_RAND_{}".format(j),"MM_{}".format(SubtractedParticle), 1000, inpDict["mm_min"], inpDict["mm_max"])
            subDict["H_MM_nosub_SUB_DUMMY_RAND_{}".format(j)]  = TH1D("H_MM_nosub_SUB_DUMMY_RAND_{}".format(j),"MM_nosub_{}".format(SubtractedParticle), 1000, 0.7, 1.5)

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        particle_subtraction_ave(t_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg)
            
    print("\nBinning data...")
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = math.sqrt(abs(evt.emiss**2-evt.pmiss**2))

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].Fill(adj_MM)            

        if(ALLCUTS):            

            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_bin_dict["H_t_DATA_{}".format(j)].Fill(-evt.MandelT)
                    hist_bin_dict["H_Q2_DATA_{}".format(j)].Fill(evt.Q2)
                    hist_bin_dict["H_W_DATA_{}".format(j)].Fill(evt.W)                        
                    hist_bin_dict["H_epsilon_DATA_{}".format(j)].Fill(evt.epsilon)
                    hist_bin_dict["H_MM_DATA_{}".format(j)].Fill(adj_MM)

    print("\nBinning dummy...")
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = math.sqrt(abs(evt.emiss**2-evt.pmiss**2))

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:                
                    hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)].Fill(adj_MM)            

        if(ALLCUTS):

            # Loop through bins in t_dummy and identify events in specified bins
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_bin_dict["H_t_DUMMY_{}".format(j)].Fill(-evt.MandelT)
                    hist_bin_dict["H_Q2_DUMMY_{}".format(j)].Fill(evt.Q2)
                    hist_bin_dict["H_W_DUMMY_{}".format(j)].Fill(evt.W)                        
                    hist_bin_dict["H_epsilon_DUMMY_{}".format(j)].Fill(evt.epsilon)
                    hist_bin_dict["H_MM_DUMMY_{}".format(j)].Fill(adj_MM)                    
                    
    print("\nBinning rand...")
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = math.sqrt(abs(evt.emiss**2-evt.pmiss**2))

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:                
                    hist_bin_dict["H_MM_nosub_RAND_{}".format(j)].Fill(adj_MM)            

        if(ALLCUTS):

            # Loop through bins in t_rand and identify events in specified bins
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_bin_dict["H_t_RAND_{}".format(j)].Fill(-evt.MandelT)
                    hist_bin_dict["H_Q2_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_bin_dict["H_W_RAND_{}".format(j)].Fill(evt.W)                        
                    hist_bin_dict["H_epsilon_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_bin_dict["H_MM_RAND_{}".format(j)].Fill(adj_MM)
                    
    print("\nBinning dummy_rand...")
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = math.sqrt(abs(evt.emiss**2-evt.pmiss**2))

        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:                
                    hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Fill(adj_MM)            

        if(ALLCUTS):

            # Loop through bins in t_dummy_rand and identify events in specified bins
            for j in range(len(t_bins)-1):                
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)].Fill(-evt.MandelT)
                    hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)].Fill(evt.W)                        
                    hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)].Fill(adj_MM)
                    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
                    
        hist_bin_dict["H_Q2_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_W_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_bin_dict["H_t_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_epsilon_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_nosub_RAND_{}".format(j)].Scale(1/nWindows)        

        hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(hist_bin_dict["H_Q2_RAND_{}".format(j)],-1)
        hist_bin_dict["H_W_DATA_{}".format(j)].Add(hist_bin_dict["H_W_RAND_{}".format(j)],-1)
        hist_bin_dict["H_t_DATA_{}".format(j)].Add(hist_bin_dict["H_t_RAND_{}".format(j)],-1)
        hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(hist_bin_dict["H_epsilon_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_bin_dict["H_MM_nosub_RAND_{}".format(j)],-1)        

        hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)        

        hist_bin_dict["H_Q2_DUMMY_{}".format(j)].Add(hist_bin_dict["H_Q2_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_W_DUMMY_{}".format(j)].Add(hist_bin_dict["H_W_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_t_DUMMY_{}".format(j)].Add(hist_bin_dict["H_t_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_epsilon_DUMMY_{}".format(j)].Add(hist_bin_dict["H_epsilon_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_DUMMY_RAND_{}".format(j)],-1)
        hist_bin_dict["H_MM_nosub_DUMMY_{}".format(j)].Add(hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)],-1)        

        # Pion subtraction by scaling simc to peak size
        if ParticleType == "kaon":
            
            try:
                ##############
                # HARD CODED #
                ##############
                scale_factor = hist_bin_dict["H_MM_nosub_DATA_{}".format(j)] \
                               .Integral(\
                                         hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].FindBin(0.89), \
                                         hist_bin_dict["H_MM_nosub_DATA_{}".format(j)].FindBin(0.94)) \
                                         /subDict["H_MM_nosub_SUB_DATA_{}".format(j)].\
                                         Integral(\
                                                  subDict["H_MM_nosub_SUB_DATA_{}".format(j)].FindBin(0.89), \
                                                  subDict["H_MM_nosub_SUB_DATA_{}".format(j)].FindBin(0.94))
                ##############
                ##############
                ##############                
            except ZeroDivisionError:
                scale_factor = 0.0

            # Scale pion to subtraction proper peak
            subDict["H_Q2_SUB_DATA_{}".format(j)].Scale(scale_factor)
            subDict["H_W_SUB_DATA_{}".format(j)].Scale(scale_factor)
            subDict["H_t_SUB_DATA_{}".format(j)].Scale(scale_factor)
            subDict["H_epsilon_SUB_DATA_{}".format(j)].Scale(scale_factor)
            subDict["H_MM_SUB_DATA_{}".format(j)].Scale(scale_factor)
                
            hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(subDict["H_Q2_SUB_DATA_{}".format(j)],-1)
            hist_bin_dict["H_W_DATA_{}".format(j)].Add(subDict["H_W_SUB_DATA_{}".format(j)],-1)
            hist_bin_dict["H_t_DATA_{}".format(j)].Add(subDict["H_t_SUB_DATA_{}".format(j)],-1)
            hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(subDict["H_epsilon_SUB_DATA_{}".format(j)],-1)
            hist_bin_dict["H_MM_DATA_{}".format(j)].Add(subDict["H_MM_SUB_DATA_{}".format(j)],-1)

        # Fit background and subtract
        background_data_fit = bg_fit(phi_setting, inpDict, hist_bin_dict["H_MM_nosub_DATA_{}".format(j)])
        hist_bin_dict["H_Q2_DATA_{}".format(j)].Add(background_data_fit[0], -1)
        hist_bin_dict["H_W_DATA_{}".format(j)].Add(background_data_fit[0], -1)
        hist_bin_dict["H_t_DATA_{}".format(j)].Add(background_data_fit[0], -1)
        hist_bin_dict["H_epsilon_DATA_{}".format(j)].Add(background_data_fit[0], -1)
        hist_bin_dict["H_MM_DATA_{}".format(j)].Add(background_data_fit[0], -1)
            
        processed_dict["t_bin{}".format(j+1)] = {
            "H_Q2_DATA" : remove_bad_bins(hist_bin_dict["H_Q2_DATA_{}".format(j)]),
            "H_W_DATA" : remove_bad_bins(hist_bin_dict["H_W_DATA_{}".format(j)]),
            "H_t_DATA" : remove_bad_bins(hist_bin_dict["H_t_DATA_{}".format(j)]),
            "H_epsilon_DATA" : remove_bad_bins(hist_bin_dict["H_epsilon_DATA_{}".format(j)]),
            "H_MM_DATA" : remove_bad_bins(hist_bin_dict["H_MM_DATA_{}".format(j)]),
            "H_Q2_DUMMY" : remove_bad_bins(hist_bin_dict["H_Q2_DUMMY_{}".format(j)]),
            "H_W_DUMMY" : remove_bad_bins(hist_bin_dict["H_W_DUMMY_{}".format(j)]),
            "H_t_DUMMY" : remove_bad_bins(hist_bin_dict["H_t_DUMMY_{}".format(j)]),
            "H_epsilon_DUMMY" : remove_bad_bins(hist_bin_dict["H_epsilon_DUMMY_{}".format(j)]),
            "H_MM_DUMMY" : remove_bad_bins(hist_bin_dict["H_MM_DUMMY_{}".format(j)]),
        }

        # Sort dictionary keys alphabetically
        processed_dict["t_bin{}".format(j+1)] = {key : processed_dict["t_bin{}".format(j+1)][key] \
                                                 for key in sorted(processed_dict["t_bin{}".format(j+1)].keys())}

        # Include Stat box
        ROOT.gStyle.SetOptStat(1)
        for i, (key,val) in enumerate(processed_dict["t_bin{}".format(j+1)].items()):
            canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
            centroid = get_centroid(val, val.GetXaxis().GetXmin(), val.GetXaxis().GetXmax()) # [centroid, centroid_error]
            val.Draw()
            val.SetTitle("{} {}".format(phi_setting, val.GetName()))
            # Create a TLatex object to add text to the plot
            text = ROOT.TLatex()
            text.SetNDC();
            text.SetTextSize(0.04);
            text.SetTextAlign(22); # Centered alignment
            text.SetTextColor(ROOT.kRed); # Set text color to red            
            # Add the centroid value to the plot
            text.DrawLatex(0.7, 0.65, "Centroid: {:.2f}+/-{:.3f}".format(centroid[0], centroid[1]))
            if i==0 and j==0:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType))+'(')
            elif i==len(processed_dict["t_bin{}".format(j+1)].items())-1 and j==len(t_bins)-2:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType))+')')
            else:
                canvas.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_averages_data_".format(phi_setting, ParticleType)))
            del canvas
            
    return processed_dict

def bin_data(kinematic_types, tree_data, tree_dummy, t_bins, nWindows, phi_setting, inpDict):

    processed_dict = process_hist_data(tree_data, tree_dummy, t_bins, nWindows, phi_setting, inpDict)
    
    binned_dict = {}

    for kin_type in kinematic_types:

        # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
        binned_t_data = []
        binned_hist_data = []
        binned_hist_dummy = []
        
        # Loop through bins in t_data and identify events in specified bins
        for j in range(len(t_bins)-1):

            H_Q2_DATA = processed_dict["t_bin{}".format(j+1)]["H_Q2_DATA"]
            H_W_DATA = processed_dict["t_bin{}".format(j+1)]["H_W_DATA"]
            H_t_DATA = processed_dict["t_bin{}".format(j+1)]["H_t_DATA"]
            H_epsilon_DATA = processed_dict["t_bin{}".format(j+1)]["H_epsilon_DATA"]

            H_Q2_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_Q2_DUMMY"]
            H_W_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_W_DUMMY"]
            H_t_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_t_DUMMY"]
            H_epsilon_DUMMY = processed_dict["t_bin{}".format(j+1)]["H_epsilon_DUMMY"]

            # Initialize lists for tmp_binned_t_data, tmp_binned_hist_data, and tmp_binned_hist_dummy
            tmp_binned_t_data = []
            tmp_binned_hist_data = []
            tmp_binned_hist_dummy = []

            tmp_hist_data = [[],[]]
            for i in range(1, H_t_DATA.GetNbinsX() + 1):
                tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                
            tmp_binned_t_data.append(tmp_hist_data)

            if kin_type == "t":
                tmp_hist_data = [[],[]]                
                for i in range(1, H_t_DATA.GetNbinsX() + 1):        
                    tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)
            if kin_type == "Q2":
                tmp_hist_data = [[],[]]                
                for i in range(1, H_Q2_DATA.GetNbinsX() + 1):        
                    tmp_hist_data[0].append(H_Q2_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_Q2_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)
            if kin_type == "W":
                tmp_hist_data = [[],[]]                
                for i in range(1, H_W_DATA.GetNbinsX() + 1):
                    tmp_hist_data[0].append(H_W_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_W_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)        
            if kin_type == "epsilon":
                tmp_hist_data = [[],[]]                
                for i in range(1, H_epsilon_DATA.GetNbinsX() + 1):
                    tmp_hist_data[0].append(H_epsilon_DATA.GetBinCenter(i))
                    tmp_hist_data[1].append(H_epsilon_DATA.GetBinContent(i))                    
                tmp_binned_hist_data.append(tmp_hist_data)

            if kin_type == "t":
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_t_DUMMY.GetNbinsX() + 1):        
                    tmp_hist_dummy[0].append(H_t_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_t_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)
            if kin_type == "Q2":
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_Q2_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_Q2_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_Q2_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)
            if kin_type == "W":
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_W_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_W_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_W_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)        
            if kin_type == "epsilon":
                tmp_hist_dummy = [[],[]]                
                for i in range(1, H_epsilon_DUMMY.GetNbinsX() + 1):
                    tmp_hist_dummy[0].append(H_epsilon_DUMMY.GetBinCenter(i))
                    tmp_hist_dummy[1].append(H_epsilon_DUMMY.GetBinContent(i))                    
                tmp_binned_hist_dummy.append(tmp_hist_dummy)

            binned_t_data.append(tmp_binned_t_data[0]) # Save a list of hists where each one is a t-bin
            binned_hist_data.append(tmp_binned_hist_data[0])
            binned_hist_dummy.append(tmp_binned_hist_dummy[0])

            if j+1 == len(t_bins)-1:
                binned_dict[kin_type] = {
                    "binned_t_data" : binned_t_data,
                    "binned_hist_data" : binned_hist_data,
                    "binned_hist_dummy" : binned_hist_dummy
                }
        
    return binned_dict
    
def calculate_ave_data(kinematic_types, hist, t_bins, phi_bins, inpDict):

    tree_data, tree_dummy = hist["InFile_DATA"], hist["InFile_DUMMY"]
    nWindows, phi_setting = hist["nWindows"], hist["phi_setting"]
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_dict = bin_data(kinematic_types, tree_data, tree_dummy, t_bins, nWindows, phi_setting, inpDict)

    group_dict = {}
    
    for kin_type in kinematic_types:
        binned_t_data = binned_dict[kin_type]["binned_t_data"]
        binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
        binned_hist_dummy = binned_dict[kin_type]["binned_hist_dummy"]

        ave_hist = []
        ave_err_hist = []
        binned_sub_data = [[],[]]
        i=0 # iter
        print("-"*25)
        # Subtract binned_hist_dummy from binned_hist_data element-wise
        for data, dummy in zip(binned_hist_data, binned_hist_dummy):
            bin_val_data, hist_val_data = data
            bin_val_dummy, hist_val_dummy = dummy
            sub_val = np.subtract(hist_val_data, hist_val_dummy)
            try:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(sub_val * bin_val_data)
                total_count = np.sum(sub_val)
                average = weighted_sum / total_count
                if math.isnan(average) or math.isinf(average):
                    print("Empty binning for data {} (t-bin={})... ".format(kin_type, i+1))
                    #sys.exit(2)
                    average = 0.0
                ave_hist.append(average)
                # Calculate the standard deviation of the data points within the bin
                std_dev = np.std(bin_val_data)
                # Determine the number of data points within the bin
                n = len(bin_val_data)
                # Calculate the standard error of the mean (SEM) for the bin
                sem = std_dev / np.sqrt(n)
                # Append the uncertainty (SEM) to the list
                ave_err_hist.append(sem)
                #print("Weighted Sum:",weighted_sum)
                #print("Total Count:",total_count)
                #print("Average for t-bin {}:".format(i+1),average)
                binned_sub_data[0].append(bin_val_data)
                binned_sub_data[1].append(sub_val)
            except ZeroDivisionError:
                print("Empty binning for data {} (t-bin={})... ".format(kin_type, i+1))
                sys.exit(2)
            i+=1

        # Print statements to check sizes
        #print("Size of binned_t_data:", len(binned_t_data))
        #print("Size of binned_hist_data:", len(binned_hist_data))
        #print("Size of binned_hist_dummy:", len(binned_hist_dummy))
        #print("Size of binned_sub_data:", len(binned_sub_data[1]))
        #print("Size of ave_hist:", len(ave_hist))
        #print("Size of t_bins:", len(t_bins))
        #print("Size of phi_bins:", len(phi_bins), "\n")

        dict_lst = []
        for j in range(len(t_bins) - 1):
            tbin_index = j
            for k in range(len(phi_bins) - 1):
                phibin_index = k
                hist_val = [binned_sub_data[0][j], binned_sub_data[1][j]]
                ave_val = ave_hist[j]
                ave_err_val = ave_err_hist[j]
                print("Data average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, j+1, k+1, ave_val, ave_err_val))
                dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))
                
        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_ave".format(kin_type) : tup[2],
                "{}_ave_err".format(kin_type) : tup[3],
            }

        group_dict[kin_type] = groups
            
    return group_dict

def ave_per_bin_data(histlist, inpDict):

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data averages for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        binned_dict = calculate_ave_data(kinematic_types, hist, t_bins, phi_bins, inpDict)
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
                
    return {"binned_DATA" : aveDict}

##################################################################################################################################################

def process_hist_simc(tree_simc, t_bins, inpDict, iteration):

    processed_dict = {}

    OutFilename = inpDict["OutFilename"]    
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]    
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]    

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_simc_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
        
    TBRANCH_SIMC  = tree_simc.Get("h10")
    
    hist_bin_dict = {}
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):

        hist_bin_dict["H_Q2_SIMC_{}".format(j)]       = TH1D("H_Q2_SIMC_{}".format(j),"Q2", 1000, inpDict["Q2min"], inpDict["Q2max"])
        hist_bin_dict["H_W_SIMC_{}".format(j)]  = TH1D("H_W_SIMC_{}".format(j),"W ", 1000, inpDict["Wmin"], inpDict["Wmax"])
        hist_bin_dict["H_t_SIMC_{}".format(j)]       = TH1D("H_t_SIMC_{}".format(j),"-t", 1000, inpDict["tmin"], inpDict["tmax"])
        hist_bin_dict["H_epsilon_SIMC_{}".format(j)]  = TH1D("H_epsilon_SIMC_{}".format(j),"epsilon", 1000, inpDict["Epsmin"], inpDict["Epsmax"])

    print("\nBinning simc...")
    for i,evt in enumerate(TBRANCH_SIMC):

        # Progress bar
        Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

        if ParticleType == "kaon":
            ALLCUTS =  apply_simc_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.phgcer_x_det, evt.phgcer_y_det)
        else:
            ALLCUTS = apply_simc_cuts(evt, mm_min, mm_max)

        #Fill SIMC events
        if(ALLCUTS):

            # Loop through bins in t_simc and identify events in specified bins
            for j in range(len(t_bins)-1):            
                if t_bins[j] <= -evt.t <= t_bins[j+1]:
                    if iteration:
                        hist_bin_dict["H_t_SIMC_{}".format(j)].Fill(-evt.t, evt.iter_weight)
                        hist_bin_dict["H_Q2_SIMC_{}".format(j)].Fill(evt.Q2, evt.iter_weight)
                        hist_bin_dict["H_W_SIMC_{}".format(j)].Fill(evt.W, evt.iter_weight)
                        hist_bin_dict["H_epsilon_SIMC_{}".format(j)].Fill(evt.epsilon, evt.iter_weight)
                    else:
                        hist_bin_dict["H_t_SIMC_{}".format(j)].Fill(-evt.t, evt.Weight)
                        hist_bin_dict["H_Q2_SIMC_{}".format(j)].Fill(evt.Q2, evt.Weight)
                        hist_bin_dict["H_W_SIMC_{}".format(j)].Fill(evt.W, evt.Weight)
                        hist_bin_dict["H_epsilon_SIMC_{}".format(j)].Fill(evt.epsilon, evt.Weight)                    

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
                        
        processed_dict["t_bin{}".format(j+1)] = {
            "H_Q2_SIMC" : hist_bin_dict["H_Q2_SIMC_{}".format(j)],
            "H_W_SIMC" : hist_bin_dict["H_W_SIMC_{}".format(j)],
            "H_t_SIMC" : hist_bin_dict["H_t_SIMC_{}".format(j)],
            "H_epsilon_SIMC" : hist_bin_dict["H_epsilon_SIMC_{}".format(j)],
        }
        
    return processed_dict                    
        
def bin_simc(kinematic_types, tree_simc, t_bins, inpDict, iteration):

    processed_dict = process_hist_simc(tree_simc, t_bins, inpDict, iteration)
    
    binned_dict = {}
    
    for kin_type in kinematic_types:
        
        # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
        binned_t_simc = []
        binned_hist_simc = []
    
        # Loop through bins in t_simc and identify events in specified bins
        for j in range(len(t_bins)-1):

            H_Q2_SIMC = processed_dict["t_bin{}".format(j+1)]["H_Q2_SIMC"]
            H_W_SIMC = processed_dict["t_bin{}".format(j+1)]["H_W_SIMC"]
            H_t_SIMC = processed_dict["t_bin{}".format(j+1)]["H_t_SIMC"]
            H_epsilon_SIMC = processed_dict["t_bin{}".format(j+1)]["H_epsilon_SIMC"]
            
            # Initialize lists for tmp_binned_t_simc, tmp_binned_hist_simc, and tmp_binned_hist_dummy
            tmp_binned_t_simc = []
            tmp_binned_hist_simc = []

            tmp_hist_simc = [[],[]]
            for i in range(1, H_t_SIMC.GetNbinsX() + 1):
                tmp_hist_simc[0].append(H_t_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_t_SIMC.GetBinContent(i))                
            tmp_binned_t_simc.append(tmp_hist_simc)

            if kin_type == "t":
                tmp_binned_hist_simc.append(tmp_hist_simc)
            if kin_type == "Q2":
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_Q2_SIMC.GetNbinsX() + 1):        
                    tmp_hist_simc[0].append(H_Q2_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_Q2_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)
            if kin_type == "W":
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_W_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc[0].append(H_W_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_W_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)        
            if kin_type == "epsilon":
                tmp_hist_simc = [[],[]]                
                for i in range(1, H_epsilon_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc[0].append(H_epsilon_SIMC.GetBinCenter(i))
                    tmp_hist_simc[1].append(H_epsilon_SIMC.GetBinContent(i))                    
                tmp_binned_hist_simc.append(tmp_hist_simc)

            binned_t_simc.append(tmp_binned_t_simc[0]) # Save a list of hists where each one is a t-bin
            binned_hist_simc.append(tmp_binned_hist_simc[0])

            if j+1 == len(t_bins)-1:
                binned_dict[kin_type] = {
                    "binned_t_simc" : binned_t_simc,
                    "binned_hist_simc" : binned_hist_simc
                }
        
    return binned_dict                

def calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict, iteration):

    tree_simc = hist["InFile_SIMC"]
    
    # Initialize lists for binned_t_data, binned_hist_data
    binned_dict = bin_simc(kinematic_types, tree_simc, t_bins, inpDict, iteration)

    group_dict = {}
    
    for kin_type in kinematic_types:
        binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
        binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]
    
        ave_hist = []
        ave_err_hist = []
        binned_sub_simc = [[],[]]
        i=0 # iter
        print("-"*25)
        for simc in binned_hist_simc:
            bin_val_simc, hist_val_simc = simc
            sub_val = np.array(hist_val_simc) # No dummy subtraction for simc
            try:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(sub_val * bin_val_simc)
                total_count = np.sum(sub_val)
                average = weighted_sum / total_count
                if math.isnan(average) or math.isinf(average):
                    print("Empty binning for simc {} (t-bin={})... ".format(kin_type, i+1))
                    #sys.exit(2)
                    average = 0.0
                ave_hist.append(average)
                # Calculate the standard deviation of the simc points within the bin
                std_dev = np.std(bin_val_simc)
                # Determine the number of simc points within the bin
                n = len(bin_val_simc)
                # Calculate the standard error of the mean (SEM) for the bin
                sem = std_dev / np.sqrt(n)
                # Append the uncertainty (SEM) to the list
                ave_err_hist.append(sem)                
                #print("Weighted Sum:",weighted_sum)
                #print("Total Count:",total_count)
                #print("Average for t-bin {}:".format(i+1),average)
                binned_sub_simc[0].append(bin_val_simc)
                binned_sub_simc[1].append(sub_val)
            except ZeroDivisionError:
                print("Empty binning for simc {} (t-bin={})... ".format(kin_type, i+1))
                sys.exit(2)
            i+=1

        # Print statements to check sizes
        #print("Size of binned_t_simc:", len(binned_t_simc))
        #print("Size of binned_hist_simc:", len(binned_hist_simc))
        #print("Size of binned_sub_simc:", len(binned_sub_simc[1]))
        #print("Size of ave_hist:", len(ave_hist))
        #print("Size of t_bins:", len(t_bins))
        #print("Size of phi_bins:", len(phi_bins), "\n")

        dict_lst = []
        for j in range(len(t_bins) - 1):
            tbin_index = j
            for k in range(len(phi_bins) - 1):
                phibin_index = k
                hist_val = [binned_sub_simc[0][j], binned_sub_simc[1][j]]
                ave_val = ave_hist[j]
                ave_err_val = ave_err_hist[j]
                print("Simc average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, j+1, k+1, ave_val, ave_err_val))
                dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))

        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_ave".format(kin_type) : tup[2],
                "{}_ave_err".format(kin_type) : tup[3],
            }                    
    
        group_dict[kin_type] = groups
            
    return group_dict

def ave_per_bin_simc(histlist, inpDict, iteration=False):

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding simc averages for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        binned_dict = calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict, iteration)
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
        
    return {"binned_SIMC" : aveDict}

##################################################################################################################################################

def grab_ave_data(histlist, inpDict):

    NumtBins = inpDict["NumtBins"]
    W = inpDict["W"]
    Q2 = inpDict["Q2"]
    EPSVAL = float(inpDict["EPSVAL"] )    
    ParticleType = inpDict["ParticleType"]

    POL = float(inpDict["POL"])
        
    if POL > 0:
        polID = 'pl'
    else:
        polID = 'mn'
    
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]

    pThetaValRight = inpDict["pThetaValRight"]
    EbeamValRight = inpDict["EbeamValRight"]
    pThetaValLeft = inpDict["pThetaValLeft"]
    EbeamValLeft = inpDict["EbeamValLeft"]
    pThetaValCenter = inpDict["pThetaValCenter"]
    EbeamValCenter = inpDict["EbeamValCenter"]
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    aveDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }

    # Define thpq vector relative to middle setting
    if hist["phi_setting"] == "Right":
        runNums = np.array([int(x) for x in runNumRight.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,hist["phi_setting"],ParticleType,run)
            if os.path.exists(pid_log):
                thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                ebeam_right = float(EbeamValRight[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_-{}.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))                
                break
            else:
                continue

    if hist["phi_setting"] == "Left":
        runNums = np.array([int(x) for x in runNumLeft.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,hist["phi_setting"],ParticleType,run)
            if os.path.exists(pid_log):
                thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                ebeam_left = float(EbeamValLeft[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))                
                break
            else:
                continue

    if hist["phi_setting"] == "Center":
        runNums = np.array([int(x) for x in runNumCenter.split(' ')])
        for i, run in enumerate(runNums):
            pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,hist["phi_setting"],ParticleType,run)
            if os.path.exists(pid_log):
                thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                ebeam_center = float(EbeamValCenter[i])
                f_kindata = '{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat'.format(ParticleType, polID, Q2.replace("p",""), \
                                                                                 W.replace("p",""), float(EPSVAL)*100, int(thpq_center*1000))                
                break
            else:
                continue
    
    #f_avek = '{}/averages/avek.Q{}W{}.dat'.format(ParticleType, Q2.replace("p",""), W.replace("p",""))
    
    # List of kinematic types
    kinematic_types = ["Q2", "W", "t", "epsilon"]

    # Loop through histlist and update aveDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data averages for {}...".format(hist["phi_setting"]))
        print("\nIteration, therefore grabbing data from {}...".format(f_kindata))        
        print("-"*25)
        print("-"*25)
        aveDict[hist["phi_setting"]] = {}
        group_dict = {}
        for kin_type in kinematic_types:
            with open(f_kindata, 'r') as f:
                lines = f.readlines()
            dict_lst = []
            for j, line in enumerate(lines): # j is t-bin
                for k in range(len(phi_bins) - 1):
                    line_lst = line.split(" ") # aveQ2, errQ2, aveW, errW, avett, errtt
                    if kin_type == "Q2":
                        ave_val = float(line_lst[0])
                        ave_err_val = float(line_lst[1])                        
                    if kin_type == "W":
                        ave_val = float(line_lst[2])
                        ave_err_val = float(line_lst[3])
                    if kin_type == "t":
                        ave_val = float(line_lst[4])
                        ave_err_val = float(line_lst[5])
                    tbin_index = j
                    phibin_index = k
                    print("Data average {} for t-bin {} phi-bin {}: {:.3f} +/- {:.3e}".format(kin_type, tbin_index+1, phibin_index+1, ave_val, ave_err_val))
                    dict_lst.append((tbin_index, phibin_index, ave_val, ave_err_val))

            # Group the tuples by the first two elements using defaultdict
            groups = defaultdict(list)
            for tup in dict_lst:
                key = (tup[0], tup[1])
                groups[key] = {
                    "{}_ave".format(kin_type) : tup[2],
                    "{}_ave_err".format(kin_type) : tup[3],                    
                }
                
            group_dict[kin_type] = groups
        
        binned_dict = group_dict
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
                
    return {"binned_DATA" : aveDict}
