#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-19 15:23:49 trottar"
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

lt=Root(os.path.realpath(__file__),"Plot_Prod")

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
from utility import get_histogram

##################################################################################################################################################

def process_hist_data(tree_data, tree_dummy, t_bins, phi_bins, nWindows, inpDict):

    processed_dict = {}
    
    ParticleType = inpDict["ParticleType"]

    Q2 = inpDict["Q2"]
    W = inpDict["W"]
    EPSSET = inpDict["EPSSET"]
    
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
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
    
    TBRANCH_DATA  = tree_data.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
    TBRANCH_RAND  = tree_data.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))
    
    TBRANCH_DUMMY  = tree_dummy.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
    TBRANCH_DUMMY_RAND  = tree_dummy.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))

    ##############
    # HARD CODED #
    ##############

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

    hist_bin_dict = {}
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)]       = TH1D("H_MM_DATA_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)]       = TH1D("H_t_DATA_{}_{}".format(j, k),"-t", 500, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_RAND_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)]       = TH1D("H_t_RAND_{}_{}".format(j, k),"-t", 500, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_{}_{}".format(j, k),"-t", 500, inpDict["tmin"], inpDict["tmax"])

            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_MM_DUMMY_RAND_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)]       = TH1D("H_t_DUMMY_RAND_{}_{}".format(j, k),"-t", 500, inpDict["tmin"], inpDict["tmax"])

    print("\nBinning data...")
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

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):                

            # Loop through bins in t_data and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (evt.ph_q+math.pi)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(evt.ph_q+math.pi)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Fill(evt.MM)

    print("\nBinning dummy...")
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

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):                

            # Loop through bins in t_dummy and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (evt.ph_q+math.pi)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(evt.ph_q+math.pi)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Fill(evt.MM)

    print("\nBinning rand...")
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

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):                

            # Loop through bins in t_rand and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (evt.ph_q+math.pi)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(evt.ph_q+math.pi)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)].Fill(evt.MM)

    print("\nBinning dummy_rand...")
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

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):                

            # Loop through bins in t_dummy_rand and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (evt.ph_q+math.pi)*(180 / math.pi) <= phi_bins[k+1]:
                            #print(phi_bins[k]," <= ",(evt.ph_q+math.pi)*(180 / math.pi)," <= ",phi_bins[k+1])
                            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Fill(evt.MM)
                            
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):                            
                            
            hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_t_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_RAND_{}_{}".format(j, k)],-1)

            hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)

            hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)].Add(hist_bin_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)],-1)

            processed_dict["t_bin{}phi_bin{}".format(j+1,k+1)] = {
                "H_MM_DATA" : hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)],
                "H_t_DATA" : hist_bin_dict["H_t_DATA_{}_{}".format(j, k)],
                "H_MM_DUMMY" : hist_bin_dict["H_MM_DUMMY_{}_{}".format(j, k)],
                "H_t_DUMMY" : hist_bin_dict["H_t_DUMMY_{}_{}".format(j, k)],
            }
    
    return processed_dict

def bin_data(kin_type, tree_data, tree_dummy, t_bins, phi_bins, nWindows, inpDict):

    processed_dict = process_hist_data(tree_data, tree_dummy, t_bins, phi_bins, nWindows, inpDict)
    
    binned_dict = {}

    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_t_data = []
    binned_hist_data = []
    binned_hist_dummy = []

    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_DATA"]
            H_t_DATA = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_DATA"]

            H_MM_DUMMY = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_DUMMY"]
            H_t_DUMMY = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_DUMMY"]

            # Initialize lists for tmp_binned_t_data, tmp_binned_hist_data, and tmp_binned_hist_dummy
            tmp_binned_t_data = []
            tmp_binned_hist_data = []
            tmp_binned_hist_dummy = []

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

            tmp_hist_dummy = [[],[]]                
            for i in range(1, H_MM_DUMMY.GetNbinsX() + 1):
                tmp_hist_dummy[0].append(H_MM_DUMMY.GetBinCenter(i))
                tmp_hist_dummy[1].append(H_MM_DUMMY.GetBinContent(i))                    
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

def calculate_yield_data(kin_type, hist, t_bins, phi_bins, inpDict):

    tree_data, tree_dummy = hist["InFile_DATA"], hist["InFile_DUMMY"]
    nWindows, normfac_data, normfac_dummy = hist["nWindows"], hist["normfac_data"], hist["normfac_dummy"]

    # Total number of events selected for setting
    NumEvts_MM_DATA = hist["NumEvts_MM_DATA"]

    # Grab the setting by setting normalized error
    data_charge_err = inpDict["data_charge_err_{}".format(hist["phi_setting"].lower())] 
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_dict = bin_data(kin_type, tree_data, tree_dummy, t_bins, phi_bins, nWindows, inpDict)

    binned_t_data = binned_dict[kin_type]["binned_t_data"]
    binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
    binned_hist_dummy = binned_dict[kin_type]["binned_hist_dummy"]

    yield_hist = []
    yield_err_hist = []
    binned_sub_data = [[],[]]
    i=0 # iter
    print("-"*25)
    # Subtract binned_hist_dummy from binned_hist_data element-wise
    for data, dummy in zip(binned_hist_data, binned_hist_dummy):
        bin_val_data, hist_val_data = data
        bin_val_dummy, hist_val_dummy = dummy
        # Scale the lists before subtraction
        #print("Y_data = {:.5e}*{:.5e}".format(np.sum(hist_val_data), normfac_data))
        scaled_hist_val_data = [val * normfac_data for val in hist_val_data]
        scaled_hist_val_dummy = [val * normfac_dummy for val in hist_val_dummy]
        sub_val = np.subtract(scaled_hist_val_data, scaled_hist_val_dummy)
        total_count = np.sum(sub_val)
        try:
            yld = total_count # Normalization applied above
            # Calculate experimental yield error (relative error)
            yld_err = np.sqrt(data_charge_err**2+(1/np.sqrt(sum(hist_val_data)))**2)
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
        binned_sub_data[0].append(bin_val_data)
        binned_sub_data[1].append(sub_val)
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
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, (yield_err_val)*yield_val))
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

def process_hist_simc(tree_simc, t_bins, phi_bins, inpDict, iter_file=""):

    processed_dict = {}
    
    Q2 = inpDict["Q2"]
    W = inpDict["W"]
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
        sys.path.append("cuts")
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET, simc=True)
        
    ################################################################################################################################################

    if iter_file != "":
        iteration = True
        tree_simc = iter_file
    else:
        iteration = False
        
    TBRANCH_SIMC  = tree_simc.Get("h10")

    hist_bin_dict = {}
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)]       = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)
            hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)]       = TH1D("H_t_SIMC_{}_{}".format(j, k),"-t", 500, inpDict["tmin"], inpDict["tmax"])
            hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)] = TH1D("H_MM_SIMC_{}_{}".format(j, k),"MM", 500, 0.7, 1.5)

    print("\nBinning simc...")
    for i,evt in enumerate(TBRANCH_SIMC):

        # Progress bar
        Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

        # Define the acceptance cuts  
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
        HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond and not hgcer_cutg.IsInside(evt.phgcer_y_det, evt.phgcer_x_det)
            NOHOLECUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond

        else:

            ALLCUTS =  HMS_Acceptance and SHMS_Acceptance and Diamond                

        #Fill SIMC events
        if(ALLCUTS):

            # Loop through bins in t_simc and identify events in specified bins
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):            
                    if t_bins[j] <= -evt.t <= t_bins[j+1]:
                        if phi_bins[k] <= (evt.phipq+math.pi)*(180 / math.pi) <= phi_bins[k+1]:
                            if iteration:
                                hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Fill(-evt.t, evt.iter_weight)
                                hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Fill(evt.missmass, evt.iter_weight)
                            else:
                                hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)].Fill(-evt.t, evt.Weight)
                                hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)].Fill(evt.missmass, evt.Weight)
                            hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)].Fill(evt.missmass)

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
                                
            processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)] = {
                "H_MM_SIMC" : hist_bin_dict["H_MM_SIMC_{}_{}".format(j, k)],
                "H_t_SIMC" : hist_bin_dict["H_t_SIMC_{}_{}".format(j, k)],
                "NumEvts_bin_MM_SIMC_unweighted" : hist_bin_dict["H_MM_SIMC_unweighted_{}_{}".format(j, k)].Integral(),
            }
    
    return processed_dict

def bin_simc(kin_type, tree_simc, t_bins, phi_bins, inpDict, iter_file=""):

    processed_dict = process_hist_simc(tree_simc, t_bins, phi_bins, inpDict, iter_file="")
    
    binned_dict = {}

    # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
    binned_t_simc = []
    binned_hist_simc = []

    binned_unweighted_NumEvts_simc = []

    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):

            H_MM_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_MM_SIMC"]
            H_t_SIMC = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["H_t_SIMC"]
            NumEvts_bin_MM_SIMC_unweighted = processed_dict["t_bin{}phi_bin{}".format(j+1, k+1)]["NumEvts_bin_MM_SIMC_unweighted"]

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
                    "binned_unweighted_NumEvts_simc" : binned_unweighted_NumEvts_simc
                }
        
    return binned_dict

def calculate_yield_simc(kin_type, hist, t_bins, phi_bins, inpDict, iter_file=""):

    tree_simc, normfac_simc = hist["InFile_SIMC"], hist["normfac_simc"]

    # Total number of events selected for setting
    NumEvts_MM_SIMC = hist["NumEvts_MM_SIMC"]    
    
    # Initialize lists for binned_t_data, binned_hist_data
    binned_dict = bin_simc(kin_type, tree_simc, t_bins, phi_bins, inpDict, iter_file="")

    binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
    binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]

    binned_unweighted_NumEvts_simc = binned_dict[kin_type]["binned_unweighted_NumEvts_simc"]

    yield_hist = []
    yield_err_hist = []
    binned_sub_simc = [[],[]]
    i=0 # iter
    print("-"*25)
    for simc in binned_hist_simc:
        bin_val_simc, hist_val_simc = simc
        #print("Y_simc = {:.5e}*{:.5e}".format(np.sum(hist_val_simc), normfac_simc))
        sub_val = np.array(hist_val_simc) # No dummy subtraction for simc, duh
        total_count = np.sum(sub_val)
        try:
            yld = total_count*normfac_simc
            # Calculate simc yield error (relative error)
            #print("!!!!!!!!!!!!!",(1/np.sqrt(binned_unweighted_NumEvts_simc[i]))*normfac_simc)
            #yld_err = (1/np.sqrt(binned_unweighted_NumEvts_simc[i]))*normfac_simc/100 # FIX (/100 is temporary)
            # FIX: Removed norm_fac because shouldn't normalize non-weighted distribution
            yld_err = (1/np.sqrt(binned_unweighted_NumEvts_simc[i]))
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
        binned_sub_simc[1].append(sub_val)
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
            print("Simc yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(j+1, k+1, yield_val, (yield_err_val)*yield_val))
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

def find_yield_simc(histlist, inpDict, iter_file=""):
    
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
        yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_simc("yield", hist, t_bins, phi_bins, inpDict, iter_file="")
            
    return {"binned_SIMC" : yieldDict}

##################################################################################################################################################

def grab_yield_data(histlist, phisetlist, inpDict):

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
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
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
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
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
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
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
            if yield_val > 0.0:
                yield_err_val = float(line_lst[1])/yield_val
            else:
                yield_err_val = 0.0
            phibin_index = int(line_lst[2])
            tbin_index = int(line_lst[3])
            print("Data yield for t-bin {} phi-bin {}: {:.3e} +/- {:.3e}".format(tbin_index, phibin_index, yield_val, (yield_err_val)*yield_val))
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
