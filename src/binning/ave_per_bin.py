#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-17 19:26:59 trottar"
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

def bin_data(kinematic_types, tree_data, tree_dummy, t_bins, nWindows, inpDict):

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
    
    binned_dict = {}
    
    TBRANCH_DATA  = tree_data.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
    TBRANCH_RAND  = tree_data.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))
    
    TBRANCH_DUMMY  = tree_dummy.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
    TBRANCH_DUMMY_RAND  = tree_dummy.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))

    for kin_type in kinematic_types:

        # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
        binned_t_data = []
        binned_hist_data = []
        binned_hist_dummy = []
        
        # Loop through bins in t_data and identify events in specified bins
        for j in range(len(t_bins)-1):

            H_Q2_DATA       = TH1D("H_Q2_DATA","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
            H_W_DATA  = TH1D("H_W_DATA","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
            H_t_DATA       = TH1D("H_t_DATA","-t", 500, inpDict["tmin"], inpDict["tmax"])
            H_epsilon_DATA  = TH1D("H_epsilon_DATA","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])

            H_Q2_RAND       = TH1D("H_Q2_RAND","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
            H_W_RAND  = TH1D("H_W_RAND","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
            H_t_RAND       = TH1D("H_t_RAND","-t", 500, inpDict["tmin"], inpDict["tmax"])
            H_epsilon_RAND  = TH1D("H_epsilon_RAND","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])

            H_Q2_DUMMY       = TH1D("H_Q2_DUMMY","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
            H_W_DUMMY  = TH1D("H_W_DUMMY","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
            H_t_DUMMY       = TH1D("H_t_DUMMY","-t", 500, inpDict["tmin"], inpDict["tmax"])
            H_epsilon_DUMMY  = TH1D("H_epsilon_DUMMY","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])

            H_Q2_DUMMY_RAND       = TH1D("H_Q2_DUMMY_RAND","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
            H_W_DUMMY_RAND  = TH1D("H_W_DUMMY_RAND","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
            H_t_DUMMY_RAND       = TH1D("H_t_DUMMY_RAND","-t", 500, inpDict["tmin"], inpDict["tmax"])
            H_epsilon_DUMMY_RAND  = TH1D("H_epsilon_DUMMY_RAND","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])

            print("\nBinning t-bin {} data...".format(j+1))
            for i,evt in enumerate(TBRANCH_DATA):

                # Progress bar
                Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

                #CUTs Definations 
                SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
                SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

                HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
                HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

                Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

                if ParticleType == "kaon":

                    # Defined HGCer Geometric cuts
                    cutg = TCutG("cutg",21)
                    cutg.SetVarX("P_hgcer_yAtCer")
                    cutg.SetVarY("P_hgcer_xAtCer")
                    cutg.SetPoint(0,-25,2)
                    cutg.SetPoint(1,-2,2)
                    cutg.SetPoint(2,-1,2.5)
                    cutg.SetPoint(3,0,3)
                    cutg.SetPoint(4,1,3)
                    cutg.SetPoint(5,2,3.3)
                    cutg.SetPoint(6,3,3.0)
                    cutg.SetPoint(7,4,2.5)
                    cutg.SetPoint(8,5,2)
                    cutg.SetPoint(9,25,2)
                    cutg.SetPoint(10,25,0.5)
                    cutg.SetPoint(11,5,0.5)
                    cutg.SetPoint(12,4,1)
                    cutg.SetPoint(13,3,-1)
                    cutg.SetPoint(14,2,-2)
                    cutg.SetPoint(15,1,-2.3)
                    cutg.SetPoint(16,0,-1.5)
                    cutg.SetPoint(17,-1,-1)
                    cutg.SetPoint(18,-2,0.5)
                    cutg.SetPoint(19,-25,0.5)
                    cutg.SetPoint(20,-25,2)

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

                else:

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

                if(ALLCUTS):

                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        H_t_DATA.SetBinContent(j+1, -evt.MandelT)
                        H_Q2_DATA.SetBinContent(j+1, evt.Q2)
                        H_W_DATA.SetBinContent(j+1, evt.W)                        
                        H_epsilon_DATA.SetBinContent(j+1, evt.epsilon)

            print("\nBinning t-bin {} rand...".format(j+1))
            for i,evt in enumerate(TBRANCH_RAND):

                # Progress bar
                Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

                #CUTs Definations 
                SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
                SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

                HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
                HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

                Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

                if ParticleType == "kaon":

                    # Defined HGCer Geometric cuts
                    cutg = TCutG("cutg",21)
                    cutg.SetVarX("P_hgcer_yAtCer")
                    cutg.SetVarY("P_hgcer_xAtCer")
                    cutg.SetPoint(0,-25,2)
                    cutg.SetPoint(1,-2,2)
                    cutg.SetPoint(2,-1,2.5)
                    cutg.SetPoint(3,0,3)
                    cutg.SetPoint(4,1,3)
                    cutg.SetPoint(5,2,3.3)
                    cutg.SetPoint(6,3,3.0)
                    cutg.SetPoint(7,4,2.5)
                    cutg.SetPoint(8,5,2)
                    cutg.SetPoint(9,25,2)
                    cutg.SetPoint(10,25,0.5)
                    cutg.SetPoint(11,5,0.5)
                    cutg.SetPoint(12,4,1)
                    cutg.SetPoint(13,3,-1)
                    cutg.SetPoint(14,2,-2)
                    cutg.SetPoint(15,1,-2.3)
                    cutg.SetPoint(16,0,-1.5)
                    cutg.SetPoint(17,-1,-1)
                    cutg.SetPoint(18,-2,0.5)
                    cutg.SetPoint(19,-25,0.5)
                    cutg.SetPoint(20,-25,2)

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

                else:

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

                if(ALLCUTS):

                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        H_t_RAND.SetBinContent(j+1, -evt.MandelT)
                        H_Q2_RAND.SetBinContent(j+1, evt.Q2)
                        H_W_RAND.SetBinContent(j+1, evt.W)                        
                        H_epsilon_RAND.SetBinContent(j+1, evt.epsilon)

            print("\nBinning t-bin {} dummy...".format(j+1))
            for i,evt in enumerate(TBRANCH_DUMMY):

                # Progress bar
                Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

                #CUTs Definations 
                SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
                SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

                HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
                HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

                Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

                if ParticleType == "kaon":

                    # Defined HGCer Geometric cuts
                    cutg = TCutG("cutg",21)
                    cutg.SetVarX("P_hgcer_yAtCer")
                    cutg.SetVarY("P_hgcer_xAtCer")
                    cutg.SetPoint(0,-25,2)
                    cutg.SetPoint(1,-2,2)
                    cutg.SetPoint(2,-1,2.5)
                    cutg.SetPoint(3,0,3)
                    cutg.SetPoint(4,1,3)
                    cutg.SetPoint(5,2,3.3)
                    cutg.SetPoint(6,3,3.0)
                    cutg.SetPoint(7,4,2.5)
                    cutg.SetPoint(8,5,2)
                    cutg.SetPoint(9,25,2)
                    cutg.SetPoint(10,25,0.5)
                    cutg.SetPoint(11,5,0.5)
                    cutg.SetPoint(12,4,1)
                    cutg.SetPoint(13,3,-1)
                    cutg.SetPoint(14,2,-2)
                    cutg.SetPoint(15,1,-2.3)
                    cutg.SetPoint(16,0,-1.5)
                    cutg.SetPoint(17,-1,-1)
                    cutg.SetPoint(18,-2,0.5)
                    cutg.SetPoint(19,-25,0.5)
                    cutg.SetPoint(20,-25,2)

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

                else:

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

                if(ALLCUTS):

                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        H_t_DUMMY.SetBinContent(j+1, -evt.MandelT)
                        H_Q2_DUMMY.SetBinContent(j+1, evt.Q2)
                        H_W_DUMMY.SetBinContent(j+1, evt.W)                        
                        H_epsilon_DUMMY.SetBinContent(j+1, evt.epsilon)

            print("\nBinning t-bin {} dummy_rand...".format(j+1))
            for i,evt in enumerate(TBRANCH_DUMMY_RAND):

                # Progress bar
                Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

                #CUTs Definations 
                SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
                SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

                HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
                HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

                Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

                if ParticleType == "kaon":

                    # Defined HGCer Geometric cuts
                    cutg = TCutG("cutg",21)
                    cutg.SetVarX("P_hgcer_yAtCer")
                    cutg.SetVarY("P_hgcer_xAtCer")
                    cutg.SetPoint(0,-25,2)
                    cutg.SetPoint(1,-2,2)
                    cutg.SetPoint(2,-1,2.5)
                    cutg.SetPoint(3,0,3)
                    cutg.SetPoint(4,1,3)
                    cutg.SetPoint(5,2,3.3)
                    cutg.SetPoint(6,3,3.0)
                    cutg.SetPoint(7,4,2.5)
                    cutg.SetPoint(8,5,2)
                    cutg.SetPoint(9,25,2)
                    cutg.SetPoint(10,25,0.5)
                    cutg.SetPoint(11,5,0.5)
                    cutg.SetPoint(12,4,1)
                    cutg.SetPoint(13,3,-1)
                    cutg.SetPoint(14,2,-2)
                    cutg.SetPoint(15,1,-2.3)
                    cutg.SetPoint(16,0,-1.5)
                    cutg.SetPoint(17,-1,-1)
                    cutg.SetPoint(18,-2,0.5)
                    cutg.SetPoint(19,-25,0.5)
                    cutg.SetPoint(20,-25,2)

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

                else:

                    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

                if(ALLCUTS):

                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        H_t_DUMMY_RAND.SetBinContent(j+1, -evt.MandelT)
                        H_Q2_DUMMY_RAND.SetBinContent(j+1, evt.Q2)
                        H_W_DUMMY_RAND.SetBinContent(j+1, evt.W)                        
                        H_epsilon_DUMMY_RAND.SetBinContent(j+1, evt.epsilon)

            H_Q2_RAND.Scale(1/nWindows)
            H_W_RAND.Scale(1/nWindows)    
            H_t_RAND.Scale(1/nWindows)
            H_epsilon_RAND.Scale(1/nWindows)

            H_Q2_DATA.Add(H_Q2_RAND,-1)
            H_W_DATA.Add(H_W_RAND,-1)
            H_t_DATA.Add(H_t_RAND,-1)
            H_epsilon_DATA.Add(H_epsilon_RAND,-1)    

            H_Q2_DUMMY_RAND.Scale(1/nWindows)
            H_W_DUMMY_RAND.Scale(1/nWindows)    
            H_t_DUMMY_RAND.Scale(1/nWindows)
            H_epsilon_DUMMY_RAND.Scale(1/nWindows)

            H_Q2_DUMMY.Add(H_Q2_DUMMY_RAND,-1)
            H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
            H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
            H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)

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

            tmp_hist_dummy = [[],[]]                
            for i in range(1, H_t_DUMMY.GetNbinsX() + 1):
                tmp_hist_dummy[0].append(H_t_DUMMY.GetBinCenter(i))
                tmp_hist_dummy[1].append(H_t_DUMMY.GetBinContent(i))                
            tmp_binned_t_dummy.append(tmp_hist_dummy)

            if kin_type == "t":
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

    tree_data, tree_dummy, nWindows = hist["InFile_DATA"], hist["InFile_DUMMY"], hist["nWindows"]
    
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_dict = bin_data(kinematic_types, tree_data, tree_dummy, t_bins, nWindows, inpDict)

    group_dict = {}
    
    for kin_type in kinematic_types:
        binned_t_data = binned_dict[kin_type]["binned_t_data"]
        binned_hist_data = binned_dict[kin_type]["binned_hist_data"]
        if kin_type == "Q2":
            print("!!!!!!!!!!!",binned_hist_data)
        binned_hist_dummy = binned_dict[kin_type]["binned_hist_dummy"]

        ave_hist = []
        binned_sub_data = [[],[]]
        i=0 # iter
        print("-"*25)
        # Subtract binned_hist_dummy from binned_hist_data element-wise
        for data, dummy in zip(binned_hist_data, binned_hist_dummy):
            bin_val_data, hist_val_data = data
            bin_val_dummy, hist_val_dummy = dummy
            sub_val = np.subtract(hist_val_data, hist_val_dummy)
            if sub_val.size != 0:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(sub_val * bin_val_data)
                total_count = np.sum(sub_val)
                average = weighted_sum / total_count            
                ave_hist.append(average)
                print("Weighted Sum:",weighted_sum)
                print("Total Count:",total_count)
                print("Average for t-bin {}:".format(i+1),average)
                binned_sub_data[0].append(bin_val_data)
                binned_sub_data[1].append(sub_val)
            else:
                ave_hist.append(0)
                #print("Weighted Sum: N/A")
                #print("Total Count: N/A")
                #print("Average for t-bin {}: 0.0".format(i+1))
                binned_sub_data[0].append(bin_val_data)
                binned_sub_data[1].append([0]*len(bin_val_data))
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
                print("Average {} for t-bin {} phi-bin {}: {:.3f}".format(kin_type, j+1, k+1, ave_val))
                dict_lst.append((tbin_index, phibin_index, hist_val, ave_val))
                
        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_arr".format(kin_type) : tup[2],
                "{}_ave".format(kin_type) : tup[3],
            }

        group_dict[kin_type] = groups
            
    return group_dict

def bin_simc(kinematic_types, tree_simc, t_bins, inpDict):

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
    
    # Initialize lists for binned_t_simc, binned_hist_simc, and binned_hist_dummy
    binned_t_simc = []
    binned_hist_simc = []

    H_Q2_SIMC       = TH1D("H_Q2_SIMC","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_SIMC  = TH1D("H_W_SIMC","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_SIMC       = TH1D("H_t_SIMC","-t", 500, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_SIMC  = TH1D("H_epsilon_SIMC","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    
    TBRANCH_SIMC  = tree_simc.Get("h10")
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        print("\nBinning t-bin {} simc...".format(j+1))
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
                
                if t_bins[j] <= -evt.t <= t_bins[j+1]:
                    H_t_SIMC.SetBinContent(j+1, -evt.t)
                    H_Q2_SIMC.SetBinContent(j+1, evt.Q2)
                    H_W_SIMC.SetBinContent(j+1, evt.W)
                    H_epsilon_SIMC.SetBinContent(j+1, evt.epsilon)

    binned_dict = {}
                    
    binned_t_simc = []
    binned_hist_simc = []

    for kin_type in kinematic_types:
        for j in range(len(t_bins)-1):
            tmp_hist_simc = []            
            for i in range(1, H_t_SIMC.GetNbinsX() + 1):
                tmp_hist_simc.append(H_t_SIMC.GetBinCenter(i))
            binned_t_simc.append([tmp_hist_simc, len(tmp_hist_simc)])

            if kin_type == "t":
                binned_hist_simc.append([tmp_hist_simc, len(tmp_hist_simc)])
            if kin_type == "Q2":
                tmp_hist_simc = []                
                for i in range(1, H_Q2_SIMC.GetNbinsX() + 1):        
                    tmp_hist_simc.append(H_Q2_SIMC.GetBinCenter(i))
                binned_hist_simc.append([tmp_hist_simc, len(tmp_hist_simc)])
            if kin_type == "W":
                tmp_hist_simc = []                
                for i in range(1, H_W_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc.append(H_W_SIMC.GetBinCenter(i))
                binned_hist_simc.append([tmp_hist_simc, len(tmp_hist_simc)])        
            if kin_type == "epsilon":
                tmp_hist_simc = []                
                for i in range(1, H_epsilon_SIMC.GetNbinsX() + 1):
                    tmp_hist_simc.append(H_epsilon_SIMC.GetBinCenter(i))
                binned_hist_simc.append([tmp_hist_simc, len(tmp_hist_simc)])

        binned_dict[kin_type] = {
            "binned_t_data" : binned_t_data,
            "binned_hist_data" : binned_hist_data,
            "binned_hist_dummy" : binned_hist_dummy
        }
        
    return binned_dict                

def calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict):

    tree_simc = hist["InFile_SIMC"]
    
    # Initialize lists for binned_t_data, binned_hist_data
    binned_dict = bin_simc(kinematic_types, tree_simc, t_bins, inpDict)

    group_dict = {}
    
    for kin_type in kinematic_types:
        binned_t_simc = binned_dict[kin_type]["binned_t_simc"]
        binned_hist_simc = binned_dict[kin_type]["binned_hist_simc"]
    
        ave_hist = []
        binned_sub_simc = [[],[]]
        i=0 # iter
        print("-"*25)
        for simc in binned_hist_simc:
            bin_val_simc, hist_val_simc = simc
            sub_val = np.array(hist_val_simc) # No dummy subtraction for simc
            if sub_val.size != 0:
                # Calculate the weighted sum of frequencies and divide by the total count
                weighted_sum = np.sum(sub_val * bin_val_simc)
                total_count = np.sum(sub_val)
                average = weighted_sum / total_count            
                ave_hist.append(average)
                #print("Weighted Sum:",weighted_sum)
                #print("Total Count:",total_count)
                #print("Average for t-bin {}:".format(i+1),average)
                binned_sub_simc[0].append(bin_val_simc)
                binned_sub_simc[1].append(sub_val)
            else:
                ave_hist.append(0)
                #print("Weighted Sum: N/A")
                #print("Total Count: N/A")
                #print("Average for t-bin {}: 0.0".format(i+1))
                binned_sub_simc[0].append(bin_val_simc)
                binned_sub_simc[1].append([0]*len(bin_val_simc))
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
                print("Average {} for t-bin {} phi-bin {}: {:.3f}".format(kin_type, j+1, k+1, ave_val))            
                dict_lst.append((tbin_index, phibin_index, hist_val, ave_val))

        # Group the tuples by the first two elements using defaultdict
        groups = defaultdict(list)
        for tup in dict_lst:
            key = (tup[0], tup[1])
            groups[key] = {
                "{}_arr".format(kin_type) : tup[2],
                "{}_ave".format(kin_type) : tup[3],
            }                    
    
        group_dict[kin_type] = groups
            
    return group_dict

##################################################################################################################################################

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

def ave_per_bin_simc(histlist, inpDict):

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
        binned_dict = calculate_ave_simc(kinematic_types, hist, t_bins, phi_bins, inpDict)
        for kin_type in kinematic_types:
            aveDict[hist["phi_setting"]][kin_type] = binned_dict[kin_type]
        
    sys.exit(1)
    return {"binned_SIMC" : aveDict}
