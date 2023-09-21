#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-21 16:46:41 trottar"
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

def calculate_ratio(kin_type, phisetlist, yieldDict):

    dict_lst = []
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(yieldDict["binned_DATA"][phiset][kin_type])
        simc_key_tuples = list(yieldDict["binned_SIMC"][phiset][kin_type])
        for data_key_tuple,simc_key_tuple in zip(data_key_tuples,simc_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict[kin_type][data_key_tuple]["{}".format(kin_type)]), data_nested_dict[kin_type][data_key_tuple][kin_type]))
            # Fill histogram
            yield_data = data_nested_dict[kin_type][data_key_tuple]["{}".format(kin_type)]
            yield_simc = simc_nested_dict[kin_type][simc_key_tuple]["{}".format(kin_type)]
            try:
                ratio = yield_data/yield_simc
            except ZeroDivisionError:
                ratio = 0.0
            if math.isnan(ratio):
                ratio = 0.0
            if math.isinf(ratio):
                ratio = 0.0                
            print("Ratio for t-bin {} phi-bin {}: {:.3f}".format(i+1, j+1, ratio))
            dict_lst.append((i, j, ratio))
    
    # Group the tuples by the first two elements using defaultdict
    groups = defaultdict(list)
    for tup in dict_lst:
        key = (tup[0], tup[1])
        groups[key] = {
            "ratio".format(kin_type) : tup[2],
        }            
            
    return groups

def find_ratio(histlist, inpDict, phisetlist, yieldDict):

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

    ratioDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins
    }

    # List of kinematic types
    kinematic_types = ["yield"]    

    # Loop through histlist and update yieldDict
    for hist in histlist:
        print("\n\n")
        print("-"*25)
        print("-"*25)
        print("Finding data yields for {}...".format(hist["phi_setting"]))
        print("-"*25)
        print("-"*25)
        ratioDict[hist["phi_setting"]] = {}
        for kin_type in kinematic_types:
            ratioDict[hist["phi_setting"]]["ratio"] = calculate_ratio(kin_type, phisetlist, yieldDict)
            
    return {"binned" : ratioDict}
