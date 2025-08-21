#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-05-22 10:48:36 trottar"
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
# Import package for progress bar
from ltsep import Misc

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

    # Create an empty list to store copied histograms
    histlist_copy = []

    # Check if the value of the dictionary is a TObject that can use .Clone()
    for hist in histlist:
        hist_copy = {}
        for key, val in hist.items():
            if hasattr(val, 'Clone') and callable(getattr(val, 'Clone')):
                # Clone the TObject if it has the 'Clone' method
                hist_copy[key] = val.Clone()
            else:
                # Otherwise, just copy the value
                hist_copy[key] = val
        # Append the copied histogram dictionary to the new list
        histlist_copy.append(hist_copy)
    
    for i,hist in enumerate(histlist_copy):

        # Unscale data to get raw number of events
        hist["H_t_DATA"].Scale(1.0 / inpDict["normfac_data"])
        hist["H_ph_q_DATA"].Scale(1.0 / inpDict["normfac_data"])
        
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
        #bad_bins_threshold = 25
        #bad_bins_threshold = 50
        #bad_bins_threshold = 75
        bad_bins_threshold = 100
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
        #bad_bins_threshold = 200 # Q2=5.5, W=3.02 
        #bad_bins_threshold = 500
        bad_bins_threshold = 1000
        #bad_bins_threshold = 2500 # Q2=4.4, W=2.74
        ##############
        ##############
        ##############
        
        def adjust_bins(x, nbin, tolerance=1e-3, max_iterations=10, edge_bias=2.0):
            """
            Adaptive binning function with edge-biased bin movement.

            Parameters:
            - x: Input data array
            - nbin: Initial number of bins
            - tolerance: Tolerance for bin edge adjustments
            - max_iterations: Maximum iterations for bin edge adjustment
            - edge_bias: Controls the non-linear scaling of bin edge movements
                         Higher values make edge bins move more aggressively
                         Default is 2.0, meaning quadratic scaling

            Returns:
            - Adjusted bin edges array
            """
            # Account for bin range
            nbin += 1
            npt = len(x)  # Total number of data points
            n_per_bin = npt // nbin  # Calculate the number of events per bin
            tmin, tmax = np.min(x), np.max(x)  # Min and max of input data

            # Initialize bin edges 
            bin_edges = np.linspace(tmin, tmax, num=nbin)

            # Calculate the number of events in each bin
            counts, _ = np.histogram(x, bins=bin_edges)

            def calculate_edge_scaling(index, total_bins):
                """
                Calculate a non-linear scaling factor for bin edge movement.
                Bins closer to the edges move more aggressively.

                Parameters:
                - index: Current bin edge index
                - total_bins: Total number of bins
                - edge_bias: Controls the non-linear scaling

                Returns:
                - Scaling factor for bin edge movement
                """
                # Normalize index to range [0, 1]
                normalized_index = index / (total_bins - 1)

                # Use non-linear scaling - quadratic or higher order based on edge_bias
                # Closer to edges (0 or 1), scaling is more aggressive
                scaling = (1 - abs(2 * normalized_index - 1)) ** edge_bias

                return scaling

            iteration = 0
            while np.any(counts < bad_bins_threshold) and iteration < max_iterations:
                # Create a copy of bin edges to modify
                new_bin_edges = bin_edges.copy()

                # Adjust bin edges
                for i in range(1, len(bin_edges) - 1):
                    # Calculate scaling factor for this bin edge
                    edge_scaling = calculate_edge_scaling(i, len(bin_edges))

                    if counts[i - 1] < bad_bins_threshold:
                        # Increase bin edge with scaled movement
                        new_bin_edges[i] += tolerance * edge_scaling
                    elif counts[i] < bad_bins_threshold:
                        # Decrease bin edge with scaled movement
                        new_bin_edges[i] -= tolerance * edge_scaling

                # Update bin edges
                bin_edges = new_bin_edges

                # Recalculate counts
                counts, _ = np.histogram(x, bins=bin_edges)

                iteration += 1

                Misc.progressBar(iteration, max_iterations-1, bar_length=25)
                
                # Progress tracking (optional)
                #print(f"Iteration {iteration}: Bin counts = {counts}")

            # Handle case where minimum events cannot be achieved
            if np.any(counts < bad_bins_threshold):
                print(f"WARNING: Could not achieve minimum {bad_bins_threshold} events in all bins after {max_iterations} iterations.")
                print("Returned bin edges may not meet the event count requirement.")

            return bin_edges
        
        print("\nFinding t bins...")
        # Histogram takes the array data set and the bins as input
        # The bins are determined by an iterative algorithm (see function above)
        #print("H_t_BinTest: ", H_t_BinTest, type(H_t_BinTest))
        try:
            bin_edges = adjust_bins(H_t_BinTest, inpDict["NumtBins"], max_iterations=50000, edge_bias=2.0)
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

    # Create an empty list to store copied histograms
    histlist_copy = []

    # Check if the value of the dictionary is a TObject that can use .Clone()
    for hist in histlist:
        hist_copy = {}
        for key, val in hist.items():
            if hasattr(val, 'Clone') and callable(getattr(val, 'Clone')):
                # Clone the TObject if it has the 'Clone' method
                hist_copy[key] = val.Clone()
            else:
                # Otherwise, just copy the value
                hist_copy[key] = val
        # Append the copied histogram dictionary to the new list
        histlist_copy.append(hist_copy)

    t_bins = histlist_copy[0]["t_bins"]
    phi_bins = histlist_copy[0]["phi_bins"]

    ################################################################################################################################################
    # Define root file trees of interest

    # Initialize NumPy arrays
    H_t_Right = np.array([])
    H_t_Left = np.array([])
    H_t_Center = np.array([])

    H_phi_Right = np.array([])
    H_phi_Left = np.array([])
    H_phi_Center = np.array([])
    
    for i,hist in enumerate(histlist_copy):

        # Unscale data to get raw number of events
        hist["H_t_DATA"].Scale(1.0 / inpDict["normfac_data"])
        hist["H_ph_q_DATA"].Scale(1.0 / inpDict["normfac_data"])
        
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
