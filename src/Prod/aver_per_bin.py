#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-11 11:34:02 trottar"
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

from utility import calculate_aver_simc, convert_TH1F_to_numpy

##################################################################################################################################################

def calculate_aver_data(hist_data, hist_dummy, t_data, t_dummy, t_bins):

    # Ensure that the input histograms are properly initialized
    if not hist_data or not hist_dummy or not t_data or not t_dummy:
        print("Error: Input histograms are not properly initialized.")
        return []

    # Bin t_data and t_dummy in t_bins
    binned_t_data = t_data.Rebin(len(t_bins)-1, "binned_t_data", array('d', t_bins))
    binned_t_dummy = t_dummy.Rebin(len(t_bins)-1, "binned_t_dummy", array('d', t_bins))

    print("Binned t data contents:", [binned_t_data.GetBinContent(i) for i in range(1, binned_t_data.GetNbinsX()+1)])
    print("Binned t dummy contents:", [binned_t_dummy.GetBinContent(i) for i in range(1, binned_t_dummy.GetNbinsX()+1)])
    
    # Get the bin numbers of t_data and t_dummy
    bin_numbers_t_data = []
    bin_numbers_t_dummy = []
    for bin_idx in range(1, len(t_bins)):
        bin_center = binned_t_data.GetBinCenter(bin_idx)
        bin_number = hist_data.FindBin(bin_center)
        bin_numbers_t_data.append(bin_number)

        bin_center_dummy = binned_t_dummy.GetBinCenter(bin_idx)
        bin_number_dummy = hist_dummy.FindBin(bin_center_dummy)
        bin_numbers_t_dummy.append(bin_number_dummy)

    # Bin hist_data and hist_dummy using the calculated bin numbers
    binned_hist_data = hist_data.Rebin(len(bin_numbers_t_data) - 1, "binned_hist_data", array('i', bin_numbers_t_data))
    binned_hist_dummy = hist_dummy.Rebin(len(bin_numbers_t_dummy) - 1, "binned_hist_dummy", array('i', bin_numbers_t_dummy))

    # Debugging step: Print histogram contents and properties
    print("Binned hist data contents:", [binned_hist_data.GetBinContent(i) for i in range(1, binned_hist_data.GetNbinsX()+1)])
    print("Binned hist dummy contents:", [binned_hist_dummy.GetBinContent(i) for i in range(1, binned_hist_dummy.GetNbinsX()+1)])
    
    # Create subtracted_hist by cloning binned_hist_data
    subtracted_hist = binned_hist_data.Clone("subtracted_hist")
    subtracted_hist.Add(binned_hist_dummy, -1)
    
    # Calculate the average per bin of the subtracted bins
    num_bins = subtracted_hist.GetNbinsX()
    averaged_values = []
    for bin_idx in range(1, num_bins+1):
        bin_content = subtracted_hist.GetBinContent(bin_idx)
        bin_width = subtracted_hist.GetBinWidth(bin_idx)
        average_value = bin_content / bin_width
        averaged_values.append(average_value)
    
    return averaged_values

##################################################################################################################################################

def aver_per_bin(histlist, inpDict):

    # Create empty histograms
    empty_hist = ROOT.TH1F()
    
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]

        # Assign histograms for Q2
        if hist["phi_setting"] == "Center":
            Q2_Center_DATA = hist["H_Q2_DATA"]
        else:
            Q2_Center_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            Q2_Left_DATA = hist["H_Q2_DATA"]
        else:
            Q2_Left_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            Q2_Right_DATA = hist["H_Q2_DATA"]
        else:
            Q2_Right_DATA = empty_hist.Clone()

        # Assign histograms for W
        if hist["phi_setting"] == "Center":
            W_Center_DATA = hist["H_W_DATA"]
        else:
            W_Center_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            W_Left_DATA = hist["H_W_DATA"]
        else:
            W_Left_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            W_Right_DATA = hist["H_W_DATA"]
        else:
            W_Right_DATA = empty_hist.Clone()

        # Assign histograms for t
        if hist["phi_setting"] == "Center":
            t_Center_DATA = hist["H_t_DATA"]
        else:
            t_Center_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            t_Left_DATA = hist["H_t_DATA"]
        else:
            t_Left_DATA = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            t_Right_DATA = hist["H_t_DATA"]
        else:
            t_Right_DATA = empty_hist.Clone()

        # Assign histograms for Q2
        if hist["phi_setting"] == "Center":
            Q2_Center_DUMMY = hist["H_Q2_DUMMY"]
        else:
            Q2_Center_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            Q2_Left_DUMMY = hist["H_Q2_DUMMY"]
        else:
            Q2_Left_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            Q2_Right_DUMMY = hist["H_Q2_DUMMY"]
        else:
            Q2_Right_DUMMY = empty_hist.Clone()

        # Assign histograms for W
        if hist["phi_setting"] == "Center":
            W_Center_DUMMY = hist["H_W_DUMMY"]
        else:
            W_Center_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            W_Left_DUMMY = hist["H_W_DUMMY"]
        else:
            W_Left_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            W_Right_DUMMY = hist["H_W_DUMMY"]
        else:
            W_Right_DUMMY = empty_hist.Clone()

        # Assign histograms for t
        if hist["phi_setting"] == "Center":
            t_Center_DUMMY = hist["H_t_DUMMY"]
        else:
            t_Center_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            t_Left_DUMMY = hist["H_t_DUMMY"]
        else:
            t_Left_DUMMY = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            t_Right_DUMMY = hist["H_t_DUMMY"]
        else:
            t_Right_DUMMY = empty_hist.Clone()

        # Assign histograms for Q2
        if hist["phi_setting"] == "Center":
            Q2_Center_SIMC = hist["H_Q2_SIMC"]
        else:
            Q2_Center_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            Q2_Left_SIMC = hist["H_Q2_SIMC"]
        else:
            Q2_Left_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            Q2_Right_SIMC = hist["H_Q2_SIMC"]
        else:
            Q2_Right_SIMC = empty_hist.Clone()

        # Assign histograms for W
        if hist["phi_setting"] == "Center":
            W_Center_SIMC = hist["H_W_SIMC"]
        else:
            W_Center_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            W_Left_SIMC = hist["H_W_SIMC"]
        else:
            W_Left_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            W_Right_SIMC = hist["H_W_SIMC"]
        else:
            W_Right_SIMC = empty_hist.Clone()

        # Assign histograms for t
        if hist["phi_setting"] == "Center":
            t_Center_SIMC = hist["H_t_SIMC"]
        else:
            t_Center_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Left":
            t_Left_SIMC = hist["H_t_SIMC"]
        else:
            t_Left_SIMC = empty_hist.Clone()
        if hist["phi_setting"] == "Right":
            t_Right_SIMC = hist["H_t_SIMC"]
        else:
            t_Right_SIMC = empty_hist.Clone()            
            
    # Combine histograms for Q2_data
    Q2_data = ROOT.TH1F("Q2_data", "Combined Q2_data Histogram")
    for bin in range(1, Q2_Center_DATA.GetNbinsX() + 1):
        combined_content = Q2_Center_DATA.GetBinContent(bin) + Q2_Left_DATA.GetBinContent(bin) + Q2_Right_DATA.GetBinContent(bin)
        Q2_data.SetBinContent(bin, combined_content)

    # Combine histograms for W_data
    W_data = ROOT.TH1F("W_data", "Combined W_data Histogram")
    for bin in range(1, W_Center_DATA.GetNbinsX() + 1):
        combined_content = W_Center_DATA.GetBinContent(bin) + W_Left_DATA.GetBinContent(bin) + W_Right_DATA.GetBinContent(bin)
        W_data.SetBinContent(bin, combined_content)

    # Combine histograms for t_data
    t_data = ROOT.TH1F("t_data", "Combined t_data Histogram")
    for bin in range(1, t_Center_DATA.GetNbinsX() + 1):
        combined_content = t_Center_DATA.GetBinContent(bin) + t_Left_DATA.GetBinContent(bin) + t_Right_DATA.GetBinContent(bin)
        t_data.SetBinContent(bin, combined_content)

    # Combine histograms for Q2_dummy
    Q2_dummy = ROOT.TH1F("Q2_dummy")
    for bin in range(1, Q2_Center_DUMMY.GetNbinsX() + 1):
        combined_content = Q2_Center_DUMMY.GetBinContent(bin) + Q2_Left_DUMMY.GetBinContent(bin) + Q2_Right_DUMMY.GetBinContent(bin)
        Q2_dummy.SetBinContent(bin, combined_content)

    # Combine histograms for W_dummy
    W_dummy = ROOT.TH1F("W_dummy", "Combined W_dummy Histogram")
    for bin in range(1, W_Center_DUMMY.GetNbinsX() + 1):
        combined_content = W_Center_DUMMY.GetBinContent(bin) + W_Left_DUMMY.GetBinContent(bin) + W_Right_DUMMY.GetBinContent(bin)
        W_dummy.SetBinContent(bin, combined_content)

    # Combine histograms for t_dummy
    t_dummy = ROOT.TH1F("t_dummy", "Combined t_dummy Histogram")
    for bin in range(1, t_Center_DUMMY.GetNbinsX() + 1):
        combined_content = t_Center_DUMMY.GetBinContent(bin) + t_Left_DUMMY.GetBinContent(bin) + t_Right_DUMMY.GetBinContent(bin)
        t_dummy.SetBinContent(bin, combined_content)

    print("@@@@@@@@@@@@@@@@@@",Q2_data.GetNbinsX(), Q2_dummy.GetNbinsX(), t_data.GetNbinsX(), t_dummy.GetNbinsX(), t_bins)
        
    Q2_aver_data = calculate_aver_data(Q2_data, Q2_dummy, t_data, t_dummy, t_bins)
    W_aver_data = calculate_aver_data(W_data, W_dummy, t_data, t_dummy, t_bins)
    t_aver_data = calculate_aver_data(t_data, t_dummy, t_data, t_dummy, t_bins)

    # Combine histograms for Q2_simc
    Q2_simc = ROOT.TH1F("Q2_simc", "Combined Q2_simc Histogram")
    for bin in range(1, Q2_Center_SIMC.GetNbinsX() + 1):
        combined_content = Q2_Center_SIMC.GetBinContent(bin) + Q2_Left_SIMC.GetBinContent(bin) + Q2_Right_SIMC.GetBinContent(bin)
        Q2_simc.SetBinContent(bin, combined_content)

    # Combine histograms for W_simc
    W_simc = ROOT.TH1F("W_simc", "Combined W_simc Histogram")
    for bin in range(1, W_Center_SIMC.GetNbinsX() + 1):
        combined_content = W_Center_SIMC.GetBinContent(bin) + W_Left_SIMC.GetBinContent(bin) + W_Right_SIMC.GetBinContent(bin)
        W_simc.SetBinContent(bin, combined_content)

    # Combine histograms for t_simc
    t_simc = ROOT.TH1F("t_simc", "Combined t_simc Histogram")
    for bin in range(1, t_Center_SIMC.GetNbinsX() + 1):
        combined_content = t_Center_SIMC.GetBinContent(bin) + t_Left_SIMC.GetBinContent(bin) + t_Right_SIMC.GetBinContent(bin)
        t_simc.SetBinContent(bin, combined_content)
    
    Q2_aver_simc = calculate_aver_simc(Q2_simc, t_bins)
    W_aver_simc = calculate_aver_simc(W_simc, t_bins)
    t_aver_simc = calculate_aver_simc(t_simc, t_bins)
    
    averDict = {
        "t_bins" : t_bins,
        "phi_bins" : phi_bins,                       
        "Q2_aver_data" : Q2_aver_data,
        "W_aver_data" : W_aver_data,
        "t_aver_data" : t_aver_data,
        "Q2_aver_simc" : Q2_aver_simc,
        "W_aver_simc" : W_aver_simc,
        "t_aver_simc" : t_aver_simc,        
    }

    return averDict
