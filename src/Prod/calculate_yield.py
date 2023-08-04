#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-04 10:55:11 trottar"
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
# Importing utility functions

from utility import weight_bins

##################################################################################################################################################

def calculate_yield(histlist, inpDict):

    # Initialize NumPy arrays before the loop
    t = np.array([])
    phi_deg = np.array([])
    Q2 = np.array([])
    W = np.array([])
    MM = np.array([])

    for hist in histlist:

        # Convert to NumPy arrays
        t = np.append(t, weight_bins(hist["H_t_DATA"]))
        phi_deg = np.append(phi_deg, [(phi + math.pi)*(180 / math.pi) for phi in weight_bins(hist["H_ph_q_DATA"])])
        Q2 = np.append(Q2, weight_bins(hist["H_Q2_DATA"]))
        W = np.append(W, weight_bins(hist["H_W_DATA"]))
        MM = np.append(MM, weight_bins(hist["H_MM_DATA"]))

    print("@@@@@@@@@@@@@@@@@--> t", t)
    print("@@@@@@@@@@@@@@@@@--> phi_deg", phi_deg)
    print("@@@@@@@@@@@@@@@@@--> Q2", Q2)
    print("@@@@@@@@@@@@@@@@@--> W", W)
    print("@@@@@@@@@@@@@@@@@--> MM", MM)

    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]
    
    for j in range(len(t_bins) - 1):
        tbin_indices = np.where((t_bins[j] <= t) & (t < t_bins[j + 1]))[0]
        if len(tbin_indices) > 0:
            tbin_index = j
            Q2_val = Q2[tbin_indices]
            W_val = W[tbin_indices]
            t_val = t[tbin_indices]
            # Append tbin_index, Q2, W, and t to aver_lst
            for k in range(len(phi_bins) - 1):
                phibin_indices = np.where((phi_bins[k] <= phi_deg) & (phi_deg < phi_bins[k + 1]))[0]
                if len(phibin_indices) > 0:
                    phibin_index = k
                    # t binning
                    #MM_val = MM[tbin_indices]
                    # t+phi binning
                    # Combine tbin_indices and phibin_indices using logical AND
                    combined_indices = np.intersect1d(tbin_indices, phibin_indices)
                    MM_val = MM[combined_indices]
                    print("________________",tbin_index, phibin_index, (MM), (Q2), (W), (t),"________________")
                    print("----------------",tbin_index, phibin_index, (MM_val), (Q2_val), (W_val), (t_val),"----------------\n")
