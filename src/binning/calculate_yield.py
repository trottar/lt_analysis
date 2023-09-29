#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-29 15:39:43 trottar"
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
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

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

def calculate_yield_data(kin_type, hist_data, hist_dummy, t_data, t_bins, phi_data, phi_bins, normfac_data):
    # Initialize lists
    binned_t_data = []
    binned_phi_data = []
    binned_hist_data = []
    binned_hist_dummy = []
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            tmp_t_data = [[],[]]
            tmp_phi_data = [[],[]]
            tmp_hist_data = [[],[]]
            tmp_hist_dummy = [[],[]]
            for tbin_index in range(1, t_data.GetNbinsX() + 1):
                tbin_center = t_data.GetBinCenter(tbin_index)
                if t_bins[j] <= tbin_center <= t_bins[j+1]:
                    if hist_data.GetBinContent(tbin_index) > 0:
                        for phibin_index in range(1, phi_data.GetNbinsX() + 1):
                            phibin_center = (phi_data.GetBinCenter(phibin_index))*(180 / math.pi)
                            if phi_bins[k] <= phibin_center <= phi_bins[k+1]:
                                if hist_data.GetBinContent(phibin_index) > 0:
                                    #print("Checking if t: {} <= {} <= {:.3f}".format(t_bins[j], tbin_center, t_bins[j+1]))
                                    #print("t-bin {}, Hist bin {:.3f} Passed with content {:.3f}".format(j+1, hist_data.GetBinCenter(tbin_index), hist_data.GetBinContent(tbin_index)))
                                    #print("Checking if phi: {} <= {:.0f} <= {}".format(phi_bins[k], phibin_center, phi_bins[k+1]))
                                    #print("phi-bin {}, Hist bin {:.0f} Passed with content {:.0f}".format(k+1, hist_data.GetBinCenter(phibin_index), hist_data.GetBinContent(phibin_index)))
                                    tmp_t_data[0].append(t_data.GetBinCenter(tbin_index))
                                    tmp_t_data[1].append(t_data.GetBinContent(tbin_index))
                                    tmp_phi_data[0].append(phi_data.GetBinCenter(phibin_index))
                                    tmp_phi_data[1].append(phi_data.GetBinContent(phibin_index))                                    
                                    tmp_hist_data[0].append(hist_data.GetBinCenter(phibin_index))
                                    tmp_hist_data[1].append(hist_data.GetBinContent(phibin_index))
                                    tmp_hist_dummy[0].append(hist_dummy.GetBinCenter(phibin_index))
                                    tmp_hist_dummy[1].append(hist_dummy.GetBinContent(phibin_index))
            binned_t_data.append(tmp_t_data)
            binned_phi_data.append(tmp_phi_data)
            binned_hist_data.append(tmp_hist_data)
            binned_hist_dummy.append(tmp_hist_dummy)

    yield_hist = []
    binned_sub_data = [[],[]]
    i=0 # iter
    print("-"*25)
    # Subtract binned_hist_dummy from binned_hist_data element-wise
    for data, dummy in zip(binned_hist_data, binned_hist_dummy):
        bin_val_data, hist_val_data = data
        bin_val_dummy, hist_val_dummy = dummy
        sub_val = np.subtract(hist_val_data, hist_val_dummy)
        total_count = np.sum(sub_val)
        yld = total_count*normfac_data
        yield_hist.append(yld)
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
            hist_val = [binned_sub_data[0][i], binned_sub_data[1][i]]
            yield_val = yield_hist[i]
            print("Yield for t-bin {} phi-bin {}: {:.3f}".format(j+1, k+1, yield_val))
            dict_lst.append((tbin_index, phibin_index, hist_val, yield_val))
            i+=1

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}_arr".format(kin_type) : tup[2],
            "{}".format(kin_type) : tup[3],
        }            
            
    return groups

def calculate_yield_simc(kin_type, hist_simc, t_simc, t_bins, phi_simc, phi_bins, normfac_simc):
    # Initialize lists
    binned_t_simc = []
    binned_phi_simc = []
    binned_hist_simc = []
    
    # Loop through bins in t_simc and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            tmp_t_simc = [[],[]]
            tmp_phi_simc = [[],[]]
            tmp_hist_simc = [[],[]]
            for tbin_index in range(1, t_simc.GetNbinsX() + 1):
                tbin_center = t_simc.GetBinCenter(tbin_index)
                if t_bins[j] <= tbin_center <= t_bins[j+1]:
                    if hist_simc.GetBinContent(tbin_index) > 0:
                        for phibin_index in range(1, phi_simc.GetNbinsX() + 1):
                            phibin_center = (phi_simc.GetBinCenter(phibin_index))*(180 / math.pi)
                            if phi_bins[k] <= phibin_center <= phi_bins[k+1]:
                                if hist_simc.GetBinContent(phibin_index) > 0:
                                    #print("Checking if t: {} <= {} <= {}".format(t_bins[j], tbin_center, t_bins[j+1]))
                                    #print("t-bin {}, Hist bin {} Passed with content {}".format(j+1, hist_simc.GetBinCenter(tbin_index), hist_simc.GetBinContent(tbin_index)))
                                    #print("Checking if phi: {} <= {} <= {}".format(phi_bins[k], phibin_center, phi_bins[k+1]))
                                    #print("phi-bin {}, Hist bin {} Passed with content {}".format(k+1, hist_simc.GetBinCenter(phibin_index), hist_simc.GetBinContent(phibin_index)))
                                    tmp_t_simc[0].append(t_simc.GetBinCenter(tbin_index))
                                    tmp_t_simc[1].append(t_simc.GetBinContent(tbin_index))
                                    tmp_phi_simc[0].append(phi_simc.GetBinCenter(phibin_index))
                                    tmp_phi_simc[1].append(phi_simc.GetBinContent(phibin_index))
                                    tmp_hist_simc[0].append(hist_simc.GetBinCenter(phibin_index))
                                    tmp_hist_simc[1].append(hist_simc.GetBinContent(phibin_index))
            binned_t_simc.append(tmp_t_simc)
            binned_phi_simc.append(tmp_phi_simc)
            binned_hist_simc.append(tmp_hist_simc)

    yield_hist = []
    binned_sub_simc = [[],[]]
    i=0 # iter
    print("-"*25)
    for simc in binned_hist_simc:
        bin_val_simc, hist_val_simc = simc
        sub_val = np.array(hist_val_simc) # No dummy subtraction for simc
        total_count = np.sum(sub_val)
        yld = total_count*normfac_simc
        yield_hist.append(yld)
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
            hist_val = [binned_sub_simc[0][i], binned_sub_simc[1][i]]
            yield_val = yield_hist[i]
            print("Yield for t-bin {} phi-bin {}: {:.3f}".format(j+1, k+1, yield_val))
            dict_lst.append((tbin_index, phibin_index, hist_val, yield_val))
            i+=1
            
    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}_arr".format(kin_type) : tup[2],
            "{}".format(kin_type) : tup[3],
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
        
    # List of kinematic types
    kinematic_types = ["MM"]

    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        for kin_type in kinematic_types:
            yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_data("yield", hist["H_{}_DATA".format(kin_type)], hist["H_{}_DUMMY".format(kin_type)], hist["H_t_DATA"], t_bins, hist["H_ph_q_DATA"], phi_bins, hist["normfac_data"])
            
    return {"binned_DATA" : yieldDict}

def find_yield_simc(histlist, inpDict):
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    yieldDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }
        
    # List of kinematic types
    kinematic_types = ["MM"]

    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding simc yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        yieldDict[hist["phi_setting"]] = {}
        for kin_type in kinematic_types:
            yieldDict[hist["phi_setting"]]["yield"] = calculate_yield_simc("yield", hist["H_{}_SIMC".format(kin_type)], hist["H_t_SIMC"], t_bins, hist["H_ph_q_SIMC"], phi_bins, hist["normfac_simc"])
            
    return {"binned_SIMC" : yieldDict}
