#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-12 12:29:38 trottar"
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
from utility import flatten_hist

################################################################################################################################################

def find_bins(histlist, inpDict):

    ################################################################################################################################################
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"]
    tmax = inpDict["tmax"]
        
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
        
        t = flatten_hist(hist["H_t_DATA"])
        phi_deg = [(phi)*(180 / math.pi) for phi in flatten_hist(hist["H_ph_q_DATA"])]
        
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

    # Apply proper boundaries for t
    # Find indices of min and max values
    min_index = np.argmin(H_t_BinTest)
    max_index = np.argmax(H_t_BinTest)
    # Replace min and max values with tmin and tmax
    H_t_BinTest[min_index] = tmin
    H_t_BinTest[max_index] = tmax
    
    # Concatenate the H_phi arrays for Right, Left, and Center
    H_phi_BinTest = np.concatenate((H_phi_Right, H_phi_Left, H_phi_Center))
    
    def find_phibins(H_phi_BinTest):

        print("\nFinding phi bins...")
        bin_edges = np.linspace(-180.0, 180.0, inpDict["NumPhiBins"]+1)
        n, bins = np.histogram(H_phi_BinTest, bin_edges)

        ##############
        # HARD CODED #
        ##############
        # Set minimum threhold number of events per bin
        # Aim for >100 events
        bad_bins_threshold = 25
        #bad_bins_threshold = 50
        #bad_bins_threshold = 75
        ##############
        ##############
        ##############

        num_phi_bins = inpDict["NumPhiBins"]
        bin_edges = np.linspace(-180.0, 180.0, num_phi_bins + 1)

        while True:
            n, bins = np.histogram(H_phi_BinTest, bin_edges)

            # Check for bad bins
            bad_bins = np.where(n < bad_bins_threshold)[0]

            # If no bad bins, we are done
            if len(bad_bins) == 0:
                break

            # Adjust the number of bins by reducing it
            num_phi_bins -= 1
            if num_phi_bins <= 4:  # Ensure there are at least 4 bins
                print("ERROR: Only {} phi-bins achieve a minimum of {} events!\nTry decreasing the minimum number of events per bin.".format(num_phi_bins, bad_bins_threshold))
                sys.exit(2)

            bin_edges = np.linspace(-180.0, 180.0, num_phi_bins + 1)

        # Redefine number of phi-bins
        if num_phi_bins != inpDict["NumPhiBins"]:
            
            print("Number of phi-bins changed from {} to: {} (threshold = {})".format(inpDict["NumPhiBins"], num_phi_bins, bad_bins_threshold))
            inpDict["NumPhiBins"] = num_phi_bins
            n, bins = np.histogram(H_phi_BinTest, bin_edges)
        
        for i,val in enumerate(n):
            print("Bin {} from {:.1f} to {:.1f} has {} events".format(i+1, bins[i], bins[i+1], n[i]))
        
        bin_centers = (bins[:-1] + bins[1:]) / 2

        print("phi_bins = ", bins)
        
        return bins

    def find_tbins(H_t_BinTest):

        ##############
        # HARD CODED #
        ##############
        # Set minimum threhold number of events per bin
        # Aim for >1000 events
        bad_bins_threshold = 200 # Q2=5.5, W=3.02 
        #bad_bins_threshold = 500
        #bad_bins_threshold = 1000
        #bad_bins_threshold = 2000
        ##############
        ##############
        ##############

        # Tolerance is determined by the number of sigfigs for |-t|
        def adjust_bins(x, nbin, tolerance=1e-3, max_iterations=10):
            # Account for bin range
            nbin += 1
            npt = len(x)  # Total number of data points
            n_per_bin = npt // nbin  # Calculate the number of events per bin
            remainder = npt % nbin  # Calculate remainder for uneven division

            # Initialize bin edges with the minimum and maximum data points
            bin_edges = []
            #bin_edges.extend(np.linspace(tmin, tmax, num=nbin+1))
            bin_edges.extend(np.linspace(tmin, tmax, num=nbin))

            # Calculate the number of events in each bin
            counts, _ = np.histogram(x, bins=bin_edges)
            # Check if any bin has fewer events than the threshold
            while np.any(counts < bad_bins_threshold):
                
                # Adjust bin edges to ensure each bin has at least bad_bins_threshold events
                for i in range(1, len(bin_edges) - 1):
                    
                    if counts[i - 1] < bad_bins_threshold:
                        # Increase bin edge to meet minimum events
                        bin_edges[i] += tolerance / 2
                    elif counts[i] < bad_bins_threshold:
                        # Decrease bin edge to meet minimum events
                        bin_edges[i] -= tolerance / 2

                    # Perform iterations to adjust bin edges
                    for _ in range(max_iterations):
                        # Calculate the number of events in each bin
                        counts, _ = np.histogram(x, bins=bin_edges)

                        # Calculate the cumulative sum of counts
                        cum_counts = np.cumsum(counts)

                        # Calculate the difference between target and actual counts per bin
                        diff_counts = n_per_bin - counts[:-1]

                        # Find the bins where the difference exceeds the tolerance
                        exceed_tolerance = np.abs(diff_counts) > tolerance

                        # If all differences are within tolerance, break the loop
                        if not np.any(exceed_tolerance):
                            break

                        # Adjust bin edges based on the differences
                        for j, exceed in enumerate(exceed_tolerance):
                            if exceed:
                                if diff_counts[j] > 0:
                                    # Increase bin edge
                                    bin_edges[j + 1] += tolerance / 2
                                else:
                                    # Decrease bin edge
                                    bin_edges[j + 1] -= tolerance / 2
                        
                # Decrease the number of bins by one
                nbin -= 1
                if nbin <= 1:
                    print("ERROR: Only {} t-bin achieved a minimum of {} events!\nTry decreasing minimum number of events per bin.".format(nbin, bad_bins_threshold))
                    sys.exit(2)

                # Recalculate bin edges and counts with the new number of bins
                bin_edges = np.linspace(tmin, tmax, num=nbin)
                counts, _ = np.histogram(x, bins=bin_edges)
                
            return np.array(bin_edges)
        
        print("\nFinding t bins...")
        # Histogram takes the array data set and the bins as input
        # The bins are determined by an iterative algorithm (see function above)
        #print("H_t_BinTest: ", H_t_BinTest, type(H_t_BinTest))
        try:
            bin_edges = adjust_bins(H_t_BinTest, inpDict["NumtBins"], max_iterations=10)
            n, bins = np.histogram(H_t_BinTest, bin_edges)
        except ValueError:
            print("ERROR: Unavoidable empty bins. Tighten t-range or adjust number of t-bins...")
            sys.exit(2)

        # Check there are good t-bins
        if np.size(n) == 0:
            print("\n\nt-binning Failed: no valid bins avaliable!\nIncrease initial number of t-bins or change t-range...")
            sys.exit(2)

        # Redefine number of t-bins
        if len(n) != inpDict["NumtBins"]:
            print("Number of t-bins changed from {} to: {} (threhold = {})".format(inpDict["NumtBins"], len(n), bad_bins_threshold))
            inpDict["NumtBins"] = len(n)
            if len(n) <= 1:
                print("ERROR: Only {} t-bin achieved a minimum of {} events!\nTry decreasing minimum number of events per bin.".format(len(n), bad_bins_threshold))
                sys.exit(2)

        for i,val in enumerate(n):
            print("Bin {} from {:.3f} to {:.3f} has {} events".format(i+1, bins[i], bins[i+1], n[i]))

        bin_centers = (bins[:-1] + bins[1:]) / 2
        
        print("t_bins = ", bins)
        
        return bins

    new_phi_bin = find_phibins(H_phi_BinTest)
    new_t_bin = find_tbins(H_t_BinTest)
    
    # Write phibin_interval for lt_analysis scripts
    lines = []
    with open("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, inpDict["Q2"].replace("p",""), inpDict["W"].replace("p","")), "w") as file:
        file.write("{}\t{}\t{}\t{}\n".format(inpDict["Q2"].replace("p","."),inpDict["W"].replace("p","."),inpDict["NumtBins"],inpDict["NumPhiBins"]))
        # Loop over phi-bin edges
        for i,phi in enumerate(new_phi_bin):
            lines.append("\t{}".format(float(phi)))
        file.writelines(lines)

    # Write t_bin_interval for lt_analysis scripts
    lines = []
    with open("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, inpDict["Q2"].replace("p",""), inpDict["W"].replace("p","")), "w") as file:
        file.write("{}\t{}\t{}\t{}\n".format(inpDict["Q2"].replace("p","."),inpDict["W"].replace("p","."),inpDict["NumtBins"],inpDict["NumPhiBins"]))
        # Loop over t-bin edges
        for i,t in enumerate(new_t_bin):
            lines.append("\t{:.2f}".format(float(t)))
        file.writelines(lines)
            
def check_bins(histlist, inpDict):

    ################################################################################################################################################
    
    ParticleType = inpDict["ParticleType"]

    tmin = inpDict["tmin"]
    tmax = inpDict["tmax"]

    t_bins = histlist[0]["t_bins"]
    phi_bins = histlist[0]["phi_bins"]

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
        
        t = flatten_hist(hist["H_t_DATA"])
        phi_deg = [(phi)*(180 / math.pi) for phi in flatten_hist(hist["H_ph_q_DATA"])]
        
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

    # Apply proper boundaries for t
    # Find indices of min and max values
    min_index = np.argmin(H_t_BinTest)
    max_index = np.argmax(H_t_BinTest)
    # Replace min and max values with tmin and tmax
    H_t_BinTest[min_index] = tmin
    H_t_BinTest[max_index] = tmax

    # Concatenate the H_phi arrays for Right, Left, and Center
    H_phi_BinTest = np.concatenate((H_phi_Right, H_phi_Left, H_phi_Center))
    
    def find_phibins(H_phi_BinTest):

        print("\nFinding phi bins...")
        n, bins = np.histogram(H_phi_BinTest, phi_bins)

        # Redefine number of t-bins
        if len(phi_bins) != inpDict["NumPhiBins"]:
            print("Number of phi-bins changed from {} to: {}".format(inpDict["NumPhiBins"], len(n)))
            inpDict["NumPhiBins"] = len(n)

        for i,val in enumerate(n):
            print("Bin {} from {:.1f} to {:.1f} has {} events".format(i+1, bins[i], bins[i+1], n[i]))
        
        print("phi_bins = ", bins)
        
    def find_tbins(H_t_BinTest):
        
        print("\nFinding t bins...")
        n, bins = np.histogram(H_t_BinTest, t_bins)

        # Redefine number of t-bins
        if len(t_bins) != inpDict["NumtBins"]:
            print("Number of t-bins changed from {} to: {}".format(inpDict["NumtBins"], len(n)))
            inpDict["NumtBins"] = len(n)
        
        for i,val in enumerate(n):
            print("Bin {} from {:.3f} to {:.3f} has {} events".format(i+1, bins[i], bins[i+1], n[i]))
        
        print("t_bins = ", bins)

    new_phi_bin = find_phibins(H_phi_BinTest)
    new_t_bin = find_tbins(H_t_BinTest)
