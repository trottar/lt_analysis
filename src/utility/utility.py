#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-16 12:33:24 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from array import array
import numpy as np
import os, subprocess

################################################################################################################################################

def show_pdf_with_evince(file_path):
    try:
        process = subprocess.Popen(['evince', file_path])
        process.wait()  # Pauses the script until Evince is closed
    except FileNotFoundError:
        print("Evince not found. Please make sure it is installed.")
    except Exception as e:
        print("An error occurred: {}".format(e))

################################################################################################################################################

def write_to_file(f_out,line,write_mode='a'):
    # Open a file in append mode
    with open(f_out, write_mode) as f:
        # Write the value of the variable to the file
        f.write(line)        

################################################################################################################################################

# Create a new dir specified by argument
def create_dir(dir_name):
    # Check if dir exists
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

################################################################################################################################################

# Convert TH1D to NumPy array
def weight_bins(histogram):
    
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

    # Calculate the total integral of the histogram (integral up to the last bin)
    total_integral = histogram.Integral()

    # Calculate the weights for each bin based on their integrals
    bin_weights = [integral / total_integral for integral in bin_integrals]

    # Weight the bin edges by the bin weights
    weighted_bin_edges = [edge * weight for edge, weight in zip(bin_edges, bin_weights)]
    
    return weighted_bin_edges

################################################################################################################################################

def match_to_bin(data):
    # Initialize a dictionary to store the matches
    match_dict = {}

    # Iterate through the given data
    for item in data:
        match = item[0][0]  # Extract the match value (e.g., 0.33)
        value = item[1][0]  # Extract the corresponding value

        # Check if the match exists in the dictionary
        if match in match_dict:
            match_dict[match].append(value)
        else:
            match_dict[match] = [value]

    # Convert the dictionary to the desired format
    return [[match, np.array(values)] for match, values in match_dict.items()]

################################################################################################################################################

# Function to check if an object is of type TH1D or TH2D
def is_hist(obj):
    return isinstance(obj, (ROOT.TH1D, ROOT.TH2D))

################################################################################################################################################

# Save histograms to root file
def hist_to_root(hist, file_name, tree_name):
    # Check if the ROOT file already exists
    root_file = ROOT.TFile.Open(file_name, "UPDATE")

    # If the file exists, check if the tree already exists
    if root_file and tree_name in [key.GetName() for key in root_file.GetListOfKeys()]:
        # If the tree already exists, open it
        tree = root_file.Get(tree_name)
    else:
        # If the tree does not exist, create a new one
        tree = ROOT.TTree(tree_name, "{} Histograms".format(tree_name.capitalize()))

    # Clone the histogram to an object (since TH1D/TH2D cannot be directly stored in a TTree)
    cloned_hist = hist.Clone()

    # Create a pointer to the cloned histogram
    hist_address = ROOT.AddressOf(cloned_hist)

    # Add the cloned histogram as a branch to the tree
    tree.Branch(hist.GetName(), hist_address)

    # Fill the tree with the cloned histogram
    tree.Fill()

    # Write the tree and close the file
    root_file.Write()
    root_file.Close()

################################################################################################################################################    
