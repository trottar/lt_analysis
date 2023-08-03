#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-03 12:53:36 trottar"
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

################################################################################################################################################

foutname = OUTPATH + "/" + inpDict["ParticleType"] + "_" + OutFilename + ".root"
fouttxt  = OUTPATH + "/" + inpDict["ParticleType"] + "_" + OutFilename + ".txt"
outputpdf  = OUTPATH + "/" + inpDict["ParticleType"] + "_" + OutFilename + ".pdf"

################################################################################################################################################

# Convert TH1F to NumPy array
def hist_to_numpy(histogram, data):
    
    # Get the number of bins in the histogram
    n_bins = histogram.GetNbinsX()

    # Calculate the integral for each bin and store it along with the bin edges
    bin_edges = []
    bin_integrals = []

    for i in range(1, n_bins + 1):
        bin_low_edge = histogram.GetXaxis().GetBinLowEdge(i)
        bin_integral = histogram.Integral(1, i)
        bin_edges.append(bin_low_edge)
        bin_integrals.append(bin_integral)

    # The last bin edge is the upper edge of the last bin
    bin_edges.append(histogram.GetXaxis().GetBinUpEdge(n_bins))

    print("^^^^^^^^^^^^^^^^^",np.average(bin_edges),"^^^^^^^^^^^^^^^^^")

    # Calculate the total integral of the histogram (integral up to the last bin)
    total_integral = histogram.Integral()

    # Calculate the weights for each bin based on their integrals
    bin_weights = [integral / total_integral for integral in bin_integrals]

    # Weight the bin edges by the bin weights
    weighted_bin_edges = [edge * weight for edge, weight in zip(bin_edges, bin_weights)]
    
    print("¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬",np.average(weighted_bin_edges),"¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬")
    
    return weighted_bin_edges

def find_bins(histlist, inpDict):

    ################################################################################################################################################
    # Define root file trees of interest

    # Initialize NumPy arrays
    H_t_Right = np.array([])
    H_t_Left = np.array([])
    H_t_Center = np.array([])

    H_phi_Right = np.array([])
    H_phi_Left = np.array([])
    H_phi_Center = np.array([])
    
    for i,hist in enumerate(histlist):
        
        t = hist_to_numpy(hist["H_t_DATA"])
        phi_deg = [(phi + math.pi)*(180 / math.pi) for phi in hist_to_numpy(hist["H_ph_q_DATA"])]
        
        if hist["phi_setting"] == 'Right':
            H_t_Right = np.append(H_t_Right, t)
            H_phi_Right = np.append(H_phi_Right, phi_deg)

        elif hist["phi_setting"] == 'Left':
            H_t_Left = np.append(H_t_Left, t)
            H_phi_Left = np.append(H_phi_Left, phi_deg)

        elif hist["phi_setting"] == 'Center':
            H_t_Center = np.append(H_t_Center, t)
            H_phi_Center = np.append(H_phi_Center, phi_deg)

    ################################################################################################################################################

    # Concatenate the H_t arrays for Right, Left, and Center
    H_t_BinTest = np.concatenate((H_t_Right, H_t_Left, H_t_Center))

    # Concatenate the H_phi arrays for Right, Left, and Center
    H_phi_BinTest = np.concatenate((H_phi_Right, H_phi_Left, H_phi_Center))
        
    return [find_phibins(H_phi_BinTest), find_tbins(H_t_BinTest)]

def find_phibins(H_phi_BinTest):

    print("\nFinding phi bins...")
    phi_arr = np.linspace(0.0, 360.0, NumPhiBins+1)

    n, bins, patches = plt.hist(H_phi_BinTest, phi_arr)

    return [n,bins]

def find_tbins(H_t_BinTest):

                
    ################################################################################################################################################

    def histedges_equalN(x, nbin):
        # Grab number of events in array
        npt = len(x)
        # One-dimensional linear interpolation for monotonically increasing sample points.
        # Returns the one-dimensional piecewise linear interpolant to a function with given
        # discrete data points (xp, fp), evaluated at x.
        #
        # np.interp(x, xp, fp)
        # x -> np.linspace(0, npt, nbin + 1) : The x-coordinates at which to evaluate the interpolated values
        # In this case, this is an array of evenly spaced t-bins
        # xp -> np.arange(npt) : The x-coordinates of the data points
        # In this case, this returns evenly spaced values within a given interval
        # yp -> np.sort(x) : he y-coordinates of the data points
        # In this case, this returns a sorted copy of the array
        return np.interp(np.linspace(0, npt, nbin + 1),np.arange(npt),np.sort(x))

    print("\nFinding t bins...")
    # Histogram takes the array data set and the bins as input
    # The bins are determined by a linear interpolation (see function above)
    # This returns the binned data with equal number of events per bin
    # Returns...
    # n -> The values of the histogram bins
    # bins -> The edges of the bins
    # patches -> Container of individual artists used to create the histogram or list of
    # such containers if there are multiple input datasets.
    n, bins, patches = plt.hist(H_t_BinTest, histedges_equalN(H_t_BinTest, NumtBins))
    
    # Write t_bin_interval for lt_analysis scripts
    lines = []
    #with open("{}/src/t_bin_interval_{}_{:.0f}".format(LTANAPATH,Q2.replace("p",""),float(EPSVAL)*100), "w") as file:
    with open("{}/src/t_bin_interval".format(LTANAPATH), "w") as file:
        file.write("{}\t{}\t{}\n".format(Q2.replace("p","."),NumtBins,NumPhiBins))
        for i,t in enumerate(bins):
            lines.append("\t{:.2f}".format(float(t)))
        file.writelines(lines)
    
    return [n,bins]

