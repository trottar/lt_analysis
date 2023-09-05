#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-04 20:05:17 trottar"
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
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
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

def calculate_yield_data(kin_type, hist_data, hist_dummy, t_data, t_bins, phi_data, phi_bins, eff_charge):
    # Initialize lists for binned_t_data, binned_hist_data, and binned_hist_dummy
    binned_t_data = []
    binned_phi_data = []
    binned_hist_data = []
    binned_hist_dummy = []
    
    t_bins = np.append(t_bins, 0.0) # Needed to fully finish loop over bins
    phi_bins = np.append(phi_bins, 0.0) # Needed to fully finish loop over bins
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
                            phibin_center = (phi_data.GetBinCenter(phibin_index)+math.pi)*(180 / math.pi)
                            if phi_bins[k] <= phibin_center <= phi_bins[k+1]:
                                if hist_data.GetBinContent(phibin_index) > 0:
                                    print("Checking if t: {} <= {} <= {}".format(t_bins[j], tbin_center, t_bins[j+1]))
                                    print("t-bin {}, Hist bin {} Passed with content {}".format(j, hist_data.GetBinCenter(tbin_index), hist_data.GetBinContent(tbin_index)))
                                    print("Checking if phi: {} <= {} <= {}".format(phi_bins[k], phibin_center, phi_bins[k+1]))
                                    print("phi-bin {}, Hist bin {} Passed with content {}".format(k, hist_data.GetBinCenter(phibin_index), hist_data.GetBinContent(phibin_index)))
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
    print("Finding yield {} per t/phi-bin...".format(kin_type))
    print("-"*25)
    # Subtract binned_hist_dummy from binned_hist_data element-wise
    for data, dummy in zip(binned_hist_data, binned_hist_dummy):
        bin_val_data, hist_val_data = data
        bin_val_dummy, hist_val_dummy = dummy
        sub_val = np.subtract(hist_val_data, hist_val_dummy)
        if sub_val.size != 0:
            total_count = np.sum(sub_val)
            yld = total_count/eff_charge
            yield_hist.append(yld)
            print("Total Count:",total_count)
            print("Effective Charge:",eff_charge)
            #print("Yield for t-bin {} phi-bin {}:".format(binned_t_data[0][i],i),yld)
            binned_sub_data[0].append(bin_val_data)
            binned_sub_data[1].append(sub_val)
        else:
            yield_hist.append(0)
            print("Total Count: N/A")
            print("Effective Charge:",eff_charge)
            #print("Yield for t-bin {} phi-bin {}: 0.0".format(binned_t_data[0][i],i))
            binned_sub_data[0].append(bin_val_data)
            binned_sub_data[1].append([0]*len(bin_val_data))
        i+=1
        print("-"*25)
    
    # Print statements to check sizes
    print("Size of binned_t_data:", len(binned_t_data))
    print("Size of binned_phi_data:", len(binned_phi_data))
    print("Size of binned_hist_data:", len(binned_hist_data))
    print("Size of binned_hist_dummy:", len(binned_hist_dummy))
    print("Size of binned_sub_data:", len(binned_sub_data[1]))
    print("Size of yield_hist:", len(yield_hist))
    print("Size of t_bins:", len(t_bins))
    print("Size of phi_bins:", len(phi_bins), "\n")

    dict_lst = []
    for j in range(len(t_bins) - 1):
        tbin_index = j
        for k in range(len(phi_bins) - 1):
            phibin_index = k
            hist_val = [binned_sub_data[0][j], binned_sub_data[1][j]]
            yield_val = yield_hist[j]
            #print("----------------------",(tbin_index, phibin_index, len(hist_val), yield_val))
            dict_lst.append((tbin_index, phibin_index, hist_val, yield_val))

    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "{}_arr".format(kin_type) : tup[2],
            "{}_yield".format(kin_type) : tup[3],
        }            
            
    return groups

def find_yield_data(histlist, inpDict):
    
    for hist in histlist:
        eff_charge = hist["normfac_data"]
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
            yieldDict[hist["phi_setting"]][kin_type] = calculate_yield_data(kin_type, hist["H_{}_DATA".format(kin_type)], hist["H_{}_DUMMY".format(kin_type)], hist["H_t_DATA"], t_bins, hist["H_ph_q_DATA"], phi_bins, eff_charge)
                
    return {"binned_DATA" : yieldDict}
