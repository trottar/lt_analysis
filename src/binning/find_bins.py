#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-19 14:33:48 trottar"
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

sys.path.append("utility")
from utility import weight_bins

################################################################################################################################################

def find_bins(histlist, inpDict):

    ################################################################################################################################################
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"]
    tmax = inpDict["tmax"]
    
    ################################################################################################################################################

    foutname = OUTPATH + "/" + inpDict["ParticleType"] + "_" + inpDict["OutFilename"] + ".root"
    fouttxt  = OUTPATH + "/" + inpDict["ParticleType"] + "_" + inpDict["OutFilename"] + ".txt"
    outputpdf  = OUTPATH + "/" + inpDict["ParticleType"] + "_" + inpDict["OutFilename"] + ".pdf"
    
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
        
        t = weight_bins(hist["H_t_DATA"])
        phi_deg = [(phi + math.pi)*(180 / math.pi) for phi in weight_bins(hist["H_ph_q_DATA"])]
        
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
    H_t_BinTest = np.append(H_t_BinTest, tmin)
    H_t_BinTest = np.append(H_t_BinTest, tmax)

    # Concatenate the H_phi arrays for Right, Left, and Center
    H_phi_BinTest = np.concatenate((H_phi_Right, H_phi_Left, H_phi_Center))

    def find_phibins(H_phi_BinTest):

        print("\nFinding phi bins...")
        phi_arr = np.linspace(0.0, 360.0, inpDict["NumPhiBins"]+1)

        n, bins, patches = plt.hist(H_phi_BinTest, phi_arr)

        bin_centers = (bins[:-1] + bins[1:]) / 2

        print("phi_bins = ", bin_centers)
        
        # Write phibin_interval for lt_analysis scripts
        lines = []
        with open("{}/src/{}/phi_bin_interval".format(LTANAPATH, ParticleType), "w") as file:
            file.write("{}\t{}\t{}\n".format(inpDict["Q2"].replace("p","."),inpDict["NumtBins"],inpDict["NumPhiBins"]))
            for i,phi in enumerate(bin_centers):
                lines.append("\t{}".format(float(phi)))
            file.writelines(lines)

        return [n,bin_centers]

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
            # x -> np.linspace(0, npt, nbin) : The x-coordinates at which to evaluate the interpolated values
            # In this case, this is an array of evenly spaced t-bins
            # xp -> np.arange(npt) : The x-coordinates of the data points
            # In this case, this returns evenly spaced values within a given interval
            # yp -> np.sort(x) : the y-coordinates of the data points
            # In this case, this returns a sorted copy of the array
            npt = len(x) 
            indices = np.linspace(0, npt-1, nbin+1).astype(int)  # Corrected line
            sorted_x = np.sort(x)
            equalN_values = sorted_x[indices]
            return np.interp(np.linspace(0, npt, nbin+2), indices, equalN_values)

        print("\nFinding t bins...")
        # Histogram takes the array data set and the bins as input
        # The bins are determined by a linear interpolation (see function above)
        # This returns the binned data with equal number of events per bin
        # Returns...
        # n -> The values of the histogram bins
        # bins -> The edges of the bins
        # patches -> Container of individual artists used to create the histogram or list of
        # such containers if there are multiple input datasets.
        n, bins, patches = plt.hist(H_t_BinTest, histedges_equalN(H_t_BinTest, inpDict["NumtBins"]))

        bin_centers = bins[1:-1]
        
        #bin_centers = (bins[:-1] + bins[1:]) / 2
        
        print("t_bins = ", bin_centers)
        
        # Write t_bin_interval for lt_analysis scripts
        lines = []
        with open("{}/src/{}/t_bin_interval".format(LTANAPATH, ParticleType), "w") as file:
            file.write("{}\t{}\t{}\n".format(inpDict["Q2"].replace("p","."),inpDict["NumtBins"],inpDict["NumPhiBins"]))
            for i,t in enumerate(bin_centers):
                lines.append("\t{:.2f}".format(float(t)))
            file.writelines(lines)

        return [n,bin_centers]
    
    return [find_phibins(H_phi_BinTest), find_tbins(H_t_BinTest)]
