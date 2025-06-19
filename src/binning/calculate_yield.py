#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 15:20:09 trottar"
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
from utility import is_hist, remove_bad_bins, fit_gaussian

##################################################################################################################################################

def process_hist_data(tree_data, tree_dummy, normfac_data, normfac_dummy, t_bins, phi_bins, nWindows, phi_setting, inpDict):

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

    # Pion subtraction by scaling pion background to peak size
    if ParticleType == "kaon":
        from particle_subtraction import particle_subtraction_yield
        SubtractedParticle = "pion"
        subDict = {}

    # Fit background and subtract
    from background_fit import bg_fit
        
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_DATA_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DATA_{}_{}".format(j, k),"MM", 200, 0.7, 1.5)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)]       = TH1D("H_t_DATA_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_RAND_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_RAND_{}_{}".format(j, k),"MM", 200, 0.7, 1.5)
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)]       = TH1D("H_t_RAND_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DUMMY_{}_{}".format(j, k),"MM", 200, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_RAND_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k),"MM", 200, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_RAND_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])

            # Pion subtraction by scaling simc to peak size
            if ParticleType == "kaon":

                subDict["H_t_SUB_DATA_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DATA_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DATA_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DATA_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DATA_{}_{}".format(j, k),"MM_nosub_{}".format(SubtractedParticle), 200, 0.7, 1.5)

                subDict["H_t_SUB_RAND_{}_{}".format(j, k)]       = TH1D("H_t_SUB_RAND_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_RAND_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_RAND_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_RAND_{}_{}".format(j, k),"MM_nosub_{}".format(SubtractedParticle), 200, 0.7, 1.5)

                subDict["H_t_SUB_DUMMY_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DUMMY_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DUMMY_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DUMMY_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, 0.7, 1.5)

                subDict["H_t_SUB_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_t_SUB_DUMMY_RAND_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])
                subDict["H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k)]  = TH1D("H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, inpDict["mm_min"], inpDict["mm_max"])
                subDict["H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k)]  \
                    = TH1D("H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k),"MM_{}".format(SubtractedParticle), 200, 0.7, 1.5)
                
    print("\nBinning data...")
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM

        ##############
        ##############        
        ##############
        
        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        

        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)
            
        if(NOMMCUTS):
            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Fill(adj_MM)
        
        if(ALLCUTS):

            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Fill(adj_MM)
                            MM_offset_DATA = evt.MM_shift-evt.MM

    print("\nBinning dummy...")
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM

        ##############
        ##############        
        ##############        
        
        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        

        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])                
                            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)

        if(ALLCUTS):                

            # Loop through bins in t_dummy and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)
                            
    print("\nBinning rand...")
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM

        ##############
        ##############        
        ##############
                
        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        

        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])                
                            hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Fill(adj_MM)

        if(ALLCUTS):                

            # Loop through bins in t_rand and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)].Fill(adj_MM)
                            
    print("\nBinning dummy_rand...")
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM

        ##############
        ##############        
        ##############        
        
       # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        

        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])                
                            hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)

        if(ALLCUTS):                

            # Loop through bins in t_dummy_rand and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(phi_shift)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)

    # Pion subtraction by scaling pion background to peak size
    if ParticleType == "kaon":
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        subDict["MM_offset_DATA"] = MM_offset_DATA
        particle_subtraction_yield(t_bins, phi_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg)        
        
    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_φ = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_φ)] for _ in range(n_t)]

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
                            
            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_nosub_RAND_{}_{}".format(j, k)],-1)            
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_RAND_{}_{}".format(j, k)],-1)

            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)],-1)            
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)],-1)        

            # Pion subtraction by scaling pion background to peak size
            if ParticleType == "kaon":
                
                try:
                    ##############
                    # HARD CODED #
                    ##############
                    pi_mm_min = 0.88 + MM_offset_DATA
                    pi_mm_max = 0.92 + MM_offset_DATA
                    # Fit amplitudes: pion background (from raw DATA) and kaon (from SUB_DATA)
                    kaon_amp = fit_gaussian(
                        hist_bin_dict[f"H_MM_nosub_DATA_{j}_{k}"],
                        pi_mm_min, pi_mm_max,
                        show_fit=False
                    )[2]

                    pion_background_amp = fit_gaussian(
                        subDict[f"H_MM_nosub_SUB_DATA_{j}_{k}"],
                        pi_mm_min, pi_mm_max,
                        show_fit=False
                    )[2]
                    scale_factor = (kaon_amp / pion_background_amp) * 0.85
                    ##############
                    ##############
                    ##############
                except ZeroDivisionError:
                    scale_factor = 0.0
                '''
                if scale_factor > 10.0:
                    print("\n\nWARNING: Pion scaling factor too large, likely no pion peak. Setting to zero....")
                    scale_factor = 0.0
                '''  
                arr_scale_factor[j][k] = scale_factor
                print(f"tbin {j} phibin {k} | Scale factor: {scale_factor:.3e}")

                subDict["H_t_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                subDict["H_MM_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)

                subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)].Scale(scale_factor)
                subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)].Scale(normfac_data)
  
                hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(subDict["H_t_SUB_DATA_{}_{}".format(j, k)],-1)
                hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(subDict["H_MM_SUB_DATA_{}_{}".format(j, k)],-1)

            # Fit background and subtract
            background_fit = bg_fit(phi_setting, inpDict, hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)])
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(background_fit[0], -1)
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(background_fit[0], -1)

            # Normalize for yields
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Scale(normfac_data)
            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Scale(normfac_dummy)          
            
            # Dummy subtraction
            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)], -1)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)], -1)             

    # Checks for first plots and calls +'(' to Print
    canvas_iter = 0
    total_plots = (len(t_bins)-1) * (len(phi_bins)-1) * len(list(["H_MM_DATA_{}_{}".format(j, k), "H_t_DATA_{}_{}".format(j, k), "H_MM_DUMMY_{}_{}".format(j, k), "H_t_DUMMY_{}_{}".format(j, k)]))-1 # '-1' to remove t-phi bin edges
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1): 
            
            processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)] = {
                "H_MM_DATA" : hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                "H_t_DATA" : hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                "H_MM_SUB_DATA" : subDict["H_MM_SUB_DATA_{}_{}".format(j, k)],
                "H_t_SUB_DATA" : subDict["H_t_SUB_DATA_{}_{}".format(j, k)],
                "scale_factor" : arr_scale_factor[j][k],
            }

            # Sort dictionary keys alphabetically
            processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)] = {key : processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)][key] \
                                                                  for key in sorted(processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)].keys())}
            
            # Clone dictionary
            cloned_dict = {}

            # Iterate through processed_dict to clone all objects
            for key, sub_dict in processed_dict.items():
                cloned_dict[key] = {}
                for sub_key, obj in sub_dict.items():
                    if is_hist(obj):
                        # Clone the object and assign it to the new dictionary
                        cloned_dict[key][sub_key] = obj.Clone()

            # Optionally sort the keys in cloned_dict if needed
            for key in cloned_dict.keys():
                cloned_dict[key] = {sub_key: cloned_dict[key][sub_key]
                                    for sub_key in sorted(cloned_dict[key].keys())}    

            # Include Stat box
            ROOT.gStyle.SetOptStat(1)
            for i, (key,val) in enumerate(cloned_dict["t_bin{}phi_bin{}".format(j+1,k+1)].items()):

                # Track the absolute first and last plots across all iterations
                is_absolute_first = (canvas_iter == 0)
                is_absolute_last = (canvas_iter == total_plots)

                #print("Processing plot: {}, Canvas iter: {}".format(key, canvas_iter))

                if is_hist(val):
                    
                    if "MM_DATA" in key:
                        # Create a new canvas for each plot
                        canvas2 = ROOT.TCanvas("canvas2_{}".format(canvas_iter), "Canvas", 800, 600)
                        
                        hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].SetLineColor(1)
                        hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Draw()
                        if ParticleType == "kaon":
                            subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)].SetLineColor(2)
                            subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)].Draw("same, E1")
                        background_fit[0].SetLineColor(3)
                        background_fit[0].Draw("same")
                        hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].SetTitle(hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].GetName())
                        
                        # Ensure correct PDF opening and closing
                        pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                        if is_absolute_first:
                            print("(")
                            canvas2.Print(pdf_name + '(')
                        elif is_absolute_last:
                            print(")")
                            canvas2.Print(pdf_name + ')')
                        else:
                            canvas2.Print(pdf_name)
                            
                        # Close the canvas2 to free up memory
                        canvas2.Close()

                    # Create a new canvas for each plot
                    canvas = ROOT.TCanvas("canvas_{}".format(canvas_iter), "Canvas", 800, 600)
                    val.Draw()
                    val.SetTitle(val.GetName())
                    
                    if "MM_DATA" in key:
                        
                        # Create a TLatex object to add text to the plot
                        text = ROOT.TLatex()
                        text.SetNDC()
                        text.SetTextSize(0.04)
                        text.SetTextAlign(22) # Centered alignment
                        text.SetTextColor(ROOT.kBlack)

                        # Add the number of mesons to the plot
                        text.DrawLatex(0.7, 0.65, "{} Yield: {:.3e}".format(ParticleType.capitalize(), val.Integral()))
                    
                    # Ensure correct PDF opening and closing
                    pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_data_".format(phi_setting, ParticleType))

                    if is_absolute_first:
                        print("(")
                        canvas.Print(pdf_name + '(')
                    elif is_absolute_last:
                        print(")")
                        canvas.Print(pdf_name + ')')
                    else:
                        canvas.Print(pdf_name)

                    # Close the canvas to free up memory
                    canvas.Close()
                        
                    # Increment canvas iterator AFTER printing
                    canvas_iter += 1
            
    return processed_dict

def bin_data(kin_type, tree_data, tree_dummy, normfac_data, normfac_dummy, t_bins, phi_bins, nWindows, phi_setting, inpDict):

    processed_dict = process_hist_data(tree_data, tree_dummy, normfac_data, normfac_dummy, t_bins, phi_bins, nWindows, phi_setting, inpDict)
    
    binned_dict = {}

    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_sub
    binned_t_data = []
    binned_hist_data = []
    binned_hist_sub = []
    mm_hist_data = []
    mm_hist_sub = []

    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_φ = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_φ)] for _ in range(n_t)]

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_DATA"]
            H_t_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_DATA"]

            H_MM_SUB_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_SUB_DATA"]
            H_t_SUB_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_SUB_DATA"]

            arr_scale_factor[j][k] = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["scale_factor"]
            print(f"2 | tbin {j} phibin {k} | Scale factor: {arr_scale_factor[j][k]:.3e}")

            mm_hist_data.append(H_MM_DATA.Clone())
            mm_hist_sub.append(H_MM_SUB_DATA.Clone())

            # Initialize lists for tmp_binned_t_data, tmp_binned_hist_data, and tmp_binned_hist_sub
            tmp_binned_t_data = []
            tmp_binned_hist_data = []
            tmp_binned_hist_sub = []

            tmp_hist_data = [[],[]]
            for i in range(1, H_t_DATA.GetNbinsX() + 1):
                tmp_hist_data[0].append(H_t_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_t_DATA.GetBinContent(i))                
            tmp_binned_t_data.append(tmp_hist_data)

            tmp_hist_data = [[],[]]                
            for i in range(1, H_MM_DATA.GetNbinsX() + 1):        
                tmp_hist_data[0].append(H_MM_DATA.GetBinCenter(i))
                tmp_hist_data[1].append(H_MM_DATA.GetBinContent(i))
            tmp_binned_hist_data.append(tmp_hist_data)

            tmp_hist_sub = [[],[]]                
            for i in range(1, H_MM_SUB_DATA.GetNbinsX() + 1):
                tmp_hist_sub[0].append(H_MM_SUB_DATA.GetBinCenter(i))
                tmp_hist_sub[1].append(H_MM_SUB_DATA.GetBinContent(i))                    
            tmp_binned_hist_sub.append(tmp_hist_sub)

            binned_t_data.append(tmp_binned_t_data[0]) # Save a list of hists where each one is a t-bin
            binned_hist_data.append(tmp_binned_hist_data[0])
            binned_hist_sub.append(tmp_binned_hist_sub[0])

            if j+1 == len(t_bins)-1:
                binned_dict[kin_type] = {
                    "binned_t_data" : binned_t_data,
                    "binned_hist_data" : binned_hist_data,
                    "binned_hist_sub" : binned_hist_sub,
                    "mm_hist_data" : mm_hist_data,
                    "mm_hist_sub" : mm_hist_sub,
                    "scale_factor" : arr_scale_factor,
                }
        
    return binned_dict

def calculate_yield_data(kin_type, hist, t_bins, phi_bins, inpDict):

    tree_data, tree_dummy = hist["InFile_DATA"], hist["InFile_DUMMY"]
    normfac_data, normfac_dummy = hist["normfac_data"], hist["normfac_dummy"]
    nWindows, phi_setting = hist["nWindows"], hist["phi_setting"]

    # Grab the setting by setting normalized error
    data_charge_err = inpDict["data_charge_err_{}".format(hist["phi_setting"].lower())]
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_dict = bin_data(kin_type, tree_data, tree_dummy, normfac_data, normfac_dummy, t_bins, phi_bins, nWindows, phi_setting, inpDict)

    binned_t_data = binned_dict[kin_type]["binned_t_data"]
    binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
    binned_hist_sub = binned_dict[kin_type]["binned_hist_sub"]
    mm_hist_data = binned_dict[kin_type]["mm_hist_data"]
    mm_hist_sub = binned_dict[kin_type]["mm_hist_sub"]

    # Initialize list saving scaled pion values    
    n_t = len(t_bins) - 1
    n_φ = len(phi_bins) - 1
    arr_scale_factor = [[0.0 for _ in range(n_φ)] for _ in range(n_t)]

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            arr_scale_factor[j][k] = binned_dict[kin_type]["scale_factor"]

    nφ = len(phi_bins) - 1
    yield_hist = []
    yield_err_hist = []
    binned_sub_data = [[],[]]
    i=0 # iter
    print("-"*25)
    # Subtract binned_hist_dummy from binned_hist_data element-wise
    for data, sub in zip(binned_hist_data, binned_hist_sub):
        j = i // nφ
        k = i %  nφ
        bin_val_data, hist_val_data = data
        bin_val_sub, hist_val_sub = sub
        bin_width_data = np.mean(np.diff(bin_val_data))
        arr_data = np.array(hist_val_data)
        bin_width_sub = np.mean(np.diff(bin_val_sub))
        arr_sub = np.array(hist_val_sub)
        try:
            yld = np.sum(arr_data)/bin_width_data
            print(f"{i} | DATA Yield: {yld:.3e} =  NumEvts: {np.sum(arr_data):.3e} / BinWidth: {bin_width_data:.3e}")
            # Calculate experimental yield error (relative error)
            # Divide by norm factor to cancel out since we need raw counts            
            yld_data_err = np.sqrt(data_charge_err**2+(1/np.sqrt(np.sum(arr_data/normfac_data)))**2)
            if math.isnan(yld_data_err) or math.isinf(yld_data_err):
                yld_data_err = 0.0
            yld_sub_err = np.sqrt(data_charge_err**2+(1/np.sqrt(np.sum(arr_sub/normfac_data)))**2)
            if math.isnan(yld_sub_err) or math.isinf(yld_sub_err):
                yld_sub_err = 0.0            
            # Convert to absolute error (required for average_ratio.f)
            yld_err = (yld_data_err**2 + (arr_scale_factor[j][k] * yld_sub_err)**2) * yld
            print(f"    | DATA Yield Error: {yld_data_err:.3e} = {np.sum(arr_data/normfac_data):.3e}")
            print(f"    | SUB Yield Error: {yld_sub_err:.3e} = {np.sum(arr_sub/normfac_data):.3e}, SCALE: {arr_scale_factor[j][k]:.3e}")            
        except ZeroDivisionError:
            yld = 0.0
            yld_err = 0.0
        yield_hist.append(yld)
        yield_err_hist.append(yld_err)
        binned_sub_data[0].append(bin_val_data)
        binned_sub_data[1].append(arr_data)
        i+=1

    # Print statements to check sizes
    #print("Size of binned_t_data:", len(binned_t_data))
    #print("Size of binned_phi_data:", len(binned_phi_data))
    #print("Size of binned_hist_data:", len(binned_hist_data))
    #print("Size of binned_hist_dummy:", len(binned_hist_dummy))
    #print("Size of binned_sub_data:", len(binned_sub_data[1]))
    #print("Size of yield_hist:", len(yield_hist))
    #print("Size of t_bins:", len(t_bins)-1)
    #print("Size of phi_bins:", len(phi_bins)-1, "\n")

    i = 0
    dict_lst = []
    for j in range(len(t_bins) - 1):
        tbin_index = j
        for k in range(len(phi_bins) - 1):
            phibin_index = k
            yield_val = yield_hist[i]
            yield_err_val = yield_err_hist[i]
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
            i+=1

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}".format(kin_type) : tup[2],
            "{}_err".format(kin_type) : tup[3],
        }            
            
    return groups

def find_yield_data(histlist, inpDict):
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_data("yield", hist, t_bins, phi_bins, inpDict)

    return {"binned_DATA" : yieldDict}

##################################################################################################################################################

def process_hist_simc(tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration):

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
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
        
    ################################################################################################################################################
        
    TBRANCH_SIMC  = tree_simc.Get("h10")
    
    hist_bin_dict = {}
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)]       = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])
            hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)]       = TH1D("H_t_SIMC_{}_{}".format(j, k),"-t", 200, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)] = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 200, inpDict["mm_min"], inpDict["mm_max"])

    print("\nBinning simc...")
    for i,evt in enumerate(TBRANCH_SIMC):

        # Progress bar
        Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        # Check if variable shift branch exists
        try:
            adj_missmass = evt.missmass_shift
        except AttributeError:
            adj_missmass = evt.missmass

        ##############
        ##############        
        ##############        
        
        if ParticleType == "kaon":          
            ALLCUTS =  apply_simc_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.phgcer_x_det, evt.phgcer_y_det)          
        else:
            ALLCUTS = apply_simc_cuts(evt, mm_min, mm_max)

        #Fill SIMC events
        if(ALLCUTS):
            
            # Phase shift to right setting
            #phi_shift = (evt.phipq+math.pi)
            phi_shift = (evt.phipq)            
            
            # Loop through bins in t_simc and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.t <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            if iteration:                                
                                hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Fill(-evt.t, evt.iter_weight)
                                hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Fill(adj_missmass, evt.iter_weight)
                            else:
                                hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Fill(-evt.t, evt.Weight)
                                hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Fill(adj_missmass, evt.Weight)
                            hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)].Fill(adj_missmass)

    # Checks for first plots and calls +'(' to Print
    canvas_iter = 0
    total_plots = (len(t_bins)-1) * (len(phi_bins)-1) * len(list(["H_MM_SIMC", "H_t_SIMC"]))-1 # '-1' to remove t-phi bin edges and NumEvts_bin_MM_SIMC_unweighted

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            # Normalize for yields
            hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Scale(normfac_simc)
            
            processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)] = {
                "H_MM_SIMC" : hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)],
                "H_t_SIMC" : hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)],
                "NumEvts_bin_MM_SIMC_unweighted" : hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)].Integral(),
            }

            # Sort dictionary keys alphabetically
            processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)] = {key : processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)][key] \
                                                                  for key in sorted(processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)].keys())}

            # Clone dictionary
            cloned_dict = {}

            # Iterate through processed_dict to clone all objects
            for key, sub_dict in processed_dict.items():
                cloned_dict[key] = {}
                for sub_key, obj in sub_dict.items():
                    if is_hist(obj):
                        # Clone the object and assign it to the new dictionary
                        cloned_dict[key][sub_key] = obj.Clone()

            # Optionally sort the keys in cloned_dict if needed
            for key in cloned_dict.keys():
                cloned_dict[key] = {sub_key: cloned_dict[key][sub_key]
                                    for sub_key in sorted(cloned_dict[key].keys())}    
            
            # Include Stat box
            ROOT.gStyle.SetOptStat(1)
            for i, (key,val) in enumerate(cloned_dict["t_bin{}phi_bin{}".format(j+1,k+1)].items()):
                
                # Create a new canvas for each plot
                canvas = ROOT.TCanvas("canvas_{}".format(canvas_iter), "Canvas", 800, 600)

                # Track the absolute first and last plots across all iterations
                is_absolute_first = (canvas_iter == 0)
                is_absolute_last = (canvas_iter == total_plots)

                #print("Processing plot: {}, Canvas iter: {}".format(key, canvas_iter))

                if is_hist(val):
                    
                    if "MM_SIMC" in key:
                        hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].SetLineColor(1)
                        hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Draw()
                        hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].SetTitle(hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].GetName())

                        # Create a TLatex object to add text to the plot
                        text = ROOT.TLatex()
                        text.SetNDC()
                        text.SetTextSize(0.04)
                        text.SetTextAlign(22) # Centered alignment
                        text.SetTextColor(ROOT.kBlack)

                        # Add the number of mesons to the plot
                        text.DrawLatex(0.7, 0.65, "{} Yield: {:.3e}".format(ParticleType.capitalize(), val.Integral()))

                    else:
                        val.Draw()
                        val.SetTitle(val.GetName())

                    # Ensure correct PDF opening and closing
                    pdf_name = outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_yield_simc_".format(phi_setting, ParticleType))

                    if is_absolute_first:
                        print("(")
                        canvas.Print(pdf_name + '(')
                    elif is_absolute_last:
                        print(")")
                        canvas.Print(pdf_name + ')')
                    else:
                        canvas.Print(pdf_name)

                    # Increment canvas iterator AFTER printing
                    canvas_iter += 1

                    # Close the canvas to free up memory
                    canvas.Close()
                
    return processed_dict

def bin_simc(kin_type, tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration):

    processed_dict = process_hist_simc(tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration)
    
    binned_dict = {}

    # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
    binned_t_simc = []
    binned_hist_simc = []
    mm_hist_simc = []
    
    binned_unweighted_NumEvts_simc = []

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_SIMC"]
            H_t_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_SIMC"]
            NumEvts_bin_MM_SIMC_unweighted = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["NumEvts_bin_MM_SIMC_unweighted"]

            mm_hist_simc.append(H_MM_SIMC.Clone())            
            
            # Initialize lists for tmp_binned_t_simc, tmp_binned_hist_simc, and tmp_binned_hist_dummy
            tmp_binned_t_simc = []
            tmp_binned_hist_simc = []

            tmp_hist_simc = [[],[]]
            for i in range(1, H_t_SIMC.GetNbinsX() + 1):
                tmp_hist_simc[0].append(H_t_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_t_SIMC.GetBinContent(i))                
            tmp_binned_t_simc.append(tmp_hist_simc)

            tmp_hist_simc = [[],[]]                
            for i in range(1, H_MM_SIMC.GetNbinsX() + 1):        
                tmp_hist_simc[0].append(H_MM_SIMC.GetBinCenter(i))
                tmp_hist_simc[1].append(H_MM_SIMC.GetBinContent(i))                    
            tmp_binned_hist_simc.append(tmp_hist_simc)

            binned_t_simc.append(tmp_binned_t_simc[0]) # Save a list of hists where each one is a t-bin
            binned_hist_simc.append(tmp_binned_hist_simc[0])
            
            binned_unweighted_NumEvts_simc.append(NumEvts_bin_MM_SIMC_unweighted)

            if j+1 == len(t_bins)-1:
                binned_dict[kin_type] = {
                    "binned_t_simc" : binned_t_simc,
                    "binned_hist_simc" : binned_hist_simc,
                    "mm_hist_simc" : mm_hist_simc,
                    "binned_unweighted_NumEvts_simc" : binned_unweighted_NumEvts_simc
                }
        
    return binned_dict

def calculate_yield_simc(kin_type, hist, t_bins, phi_bins, inpDict, iteration):

    tree_simc, normfac_simc = hist["InFile_SIMC"], hist["normfac_simc"]
    phi_setting = hist["phi_setting"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Initialize lists for binned_t_data, binned_hist_data
    binned_dict = bin_simc(kin_type, tree_simc, normfac_simc, t_bins, phi_bins, phi_setting, inpDict, iteration)

    binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
    binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]
    mm_hist_simc = binned_dict[kin_type]["mm_hist_simc"]
    
    binned_unweighted_NumEvts_simc = binned_dict[kin_type]["binned_unweighted_NumEvts_simc"]

    yield_hist = []
    yield_err_hist = []
    binned_sub_simc = [[],[]]
    i=0 # iter
    print("-"*25)
    for simc in binned_hist_simc:
        bin_val_simc, hist_val_simc = simc
        bin_width_simc = np.mean(np.diff(bin_val_simc))
        arr_simc = np.array(hist_val_simc)
        try:
            yld = np.sum(arr_simc)/bin_width_simc
            # Calculate simc yield error (relative error)
            # No norm_fac, shouldn't normalize non-weighted distribution
            yld_err = (1/np.sqrt(binned_unweighted_NumEvts_simc[i]))
            # Convert to absolute error (required for average_ratio.f)
            yld_err = yld_err*yld
        except ZeroDivisionError:
            yld = 0.0
            yld_err = 0.0
        if yld < 0.0:
            yld = 0.0
            yld_err = 0.0
        if math.isnan(yld) or math.isnan(yld_err):
            yld = 0.0
            yld_err = 0.0
        if math.isinf(yld) or math.isinf(yld_err):
            yld = 0.0
            yld_err = 0.0            
        yield_hist.append(yld)
        yield_err_hist.append(yld_err)
        binned_sub_simc[0].append(bin_val_simc)
        binned_sub_simc[1].append(arr_simc)
        i+=1

    # Print statements to check sizes
    #print("Size of binned_t_simc:", len(binned_t_simc))
    #print("Size of binned_phi_simc:", len(binned_phi_simc))
    #print("Size of binned_hist_simc:", len(binned_hist_simc))
    #print("Size of binned_sub_simc:", len(binned_sub_simc[1]))
    #print("Size of yield_hist:", len(yield_hist))
    #print("Size of t_bins:", len(t_bins)-1)
    #print("Size of phi_bins:", len(phi_bins)-1, "\n")

    i = 0
    dict_lst = []
    for j in range(len(t_bins) - 1):
        tbin_index = j
        for k in range(len(phi_bins) - 1):
            phibin_index = k
            yield_val = yield_hist[i]
            yield_err_val = yield_err_hist[i]
            print("Simc yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
            i+=1

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}".format(kin_type) : tup[2],
            "{}_err".format(kin_type) : tup[3],
        }            
            
    return groups

def find_yield_simc(histlist, inpDict, iteration=False):
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding simc yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_simc("yield", hist, t_bins, phi_bins, inpDict, iteration)
            
    return {"binned_SIMC" : yieldDict}

##################################################################################################################################################

def grab_yield_data(histlist, phisetlist, inpDict):

    OutFilename = inpDict["OutFilename"]
    
    Ws = inpDict["W"]
    Qs = inpDict["Q2"]
    Q2 = float(Qs.replace("p","."))
    W = float(Ws.replace("p","."))
    EPSVAL = float(inpDict["EPSVAL"] )
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]    
    ParticleType = inpDict["ParticleType"]
    pThetaValRight = inpDict["pThetaValRight"]
    EbeamValRight = inpDict["EbeamValRight"]
    pThetaValLeft = inpDict["pThetaValLeft"]
    EbeamValLeft = inpDict["EbeamValLeft"]
    pThetaValCenter = inpDict["pThetaValCenter"]
    EbeamValCenter = inpDict["EbeamValCenter"]    
    POL = float(inpDict["POL"])
    
    if POL > 0:
        polID = 'pl'
    else:
        polID = 'mn'
            
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # Loop through histlist and update yieldDict
    for hist in histlist:

        phiset = hist["phi_setting"]
        
        if phiset == "Right":
            runNums = np.array([int(x) for x in runNumRight.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                    ebeam_right = float(EbeamValRight[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                    Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))

        if phiset == "Left":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                    ebeam_left = float(EbeamValLeft[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                    Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))

        if phiset == "Center":
            runNums = np.array([int(x) for x in runNumCenter.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phiset}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                    ebeam_center = float(EbeamValCenter[i])
                    break
                else:
                    continue            
            f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                      Ws.replace("p",""), float(EPSVAL)*100)
            
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("\nIteration, therefore grabbing data from {}...".format(f_yield))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        with open(f_yield, 'r') as f:
            lines = f.readlines()
        dict_lst = []
        for line in lines:
            line_lst = line.split(" ") # yield, yield_err, phibin, tbin
            yield_val = float(line_lst[0])
            yield_err_val = float(line_lst[1])
            phibin_index = int(line_lst[2])-1
            tbin_index = int(line_lst[3])-1
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(tbin_index+1, phibin_index+1, yield_val, yield_err_val))
            dict_lst.append((tbin_index, phibin_index, yield_val, yield_err_val))
        
        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}".format("yield") : tup[2],
                "{}_err".format("yield") : tup[3],
            }

        yieldDict[hist["phi_setting"]]["yield"] = groups
        
    return {"binned_DATA" : yieldDict}
