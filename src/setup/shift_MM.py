#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-24 15:31:30 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TH1F, kRed
import os, sys

##################################################################################################################################################

# Input params - run number, particle type, and max number of events
runNum = sys.argv[1]
ParticleType = sys.argv[2]
ROOTPrefix = sys.argv[3]
MaxEvent = "-1"

###############################################################################################################################################
# ltsep package import and pathing definitions

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
CACHEPATH=lt.CACHEPATH

#################################################################################################################################################################

filename = f"{OUTPATH}/{ParticleType}_{runNum}_{MaxEvent}_Raw_Data.root"

# Check if the file exists
if not os.path.exists(filename):
    print(f"Error: File '{filename}' does not exist.")
    sys.exit(1)  # Exit the script with an error code

trees = [f"Uncut_{ParticleType.capitalize()}_Events", \
         f"Cut_{ParticleType.capitalize()}_Events_all_noRF", f"Cut_{ParticleType.capitalize()}_Events_prompt_noRF", f"Cut_{ParticleType.capitalize()}_Events_rand_noRF", \
         f"Cut_{ParticleType.capitalize()}_Events_all_RF", f"Cut_{ParticleType.capitalize()}_Events_prompt_RF", f"Cut_{ParticleType.capitalize()}_Events_rand_RF"]

if ParticleType == "kaon":
    MM_str = "lambda"
    MM_true = 1.1156  # Lambda peak
else:
    MM_str = "proton"
    MM_true = 0.938  # Proton peak

MM_min = MM_true-0.05
MM_max = MM_true+0.05
    
print(f"\n\nShifting missing mass to {MM_str} peak of {MM_true} for file {filename}...")
    
# Open ROOT file
file = TFile.Open(filename)

# Define a function for fitting a Gaussian with dynamically determined FWHM range
def fit_gaussian(hist, x_min, x_max):

    print(hist.GetName(),"-"*25)
    
    # Find the corresponding bin numbers
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)
    
    # Find the maximum value within the specified range
    max_bin = bin_min
    max_value = hist.GetBinContent(max_bin)
    for i in range(bin_min, bin_max):
        if hist.GetBinContent(i) > max_value:
            max_bin = i
            max_value = hist.GetBinContent(i)

    # Print the results
    print("max_bin", max_bin)
    print("max_value", max_value)
    print("bin_center",hist.GetBinCenter(max_bin))
    
    half_max = max_value*0.75
    
    # Find left and right bins closest to half-max value
    left_bin = max_bin
    right_bin = max_bin
    while hist.GetBinContent(left_bin) > half_max and left_bin > 1:
        left_bin -= 1
    while hist.GetBinContent(right_bin) > half_max and right_bin < hist.GetNbinsX():
        right_bin += 1

    #min_range = hist.GetBinCenter(max_bin-100)
    #max_range = hist.GetBinCenter(max_bin+100)

    min_range = hist.GetBinCenter(left_bin)
    print("min_range",min_range)
    max_range = hist.GetBinCenter(right_bin)
    print("max_range",max_range)
    print("-"*25)

    hist.Fit("gaus", "Q", "", min_range, max_range)
    fit_func = hist.GetFunction('gaus')
    
    fit_func.SetLineColor(kRed)
    
    mean = fit_func.GetParameter(1)
    mean_err = fit_func.GetParError(1)
    
    return [mean, mean_err]

# Function to apply mass shift
def shift_mass(tree, branch_name, shift):
    new_values = []
    for event in tree:
        original_mass = getattr(event, branch_name)
        shifted_mass = original_mass + shift
        new_values.append(shifted_mass)
    return new_values

# Fit the peak for the reference tree
reference_tree_name = f"Cut_{ParticleType.capitalize()}_Events_prompt_noRF"
print(f"\nFitting MM_{ParticleType[0].upper()} of tree {reference_tree_name}...")
reference_tree = file.Get(reference_tree_name)

# Define histogram for fitting
hist = TH1F("hist", f"MM_{ParticleType[0].upper()}", 200, 0.7, 1.5)  # Adjust binning/range as needed
reference_tree.Draw("MM>>hist")

# Fit the histogram to find the peak
peak_position = fit_gaussian(hist, MM_min, MM_max)[0] # Grab Mean for shift

# Calculate the shift
shift = MM_true - peak_position
print(f"Calculated shift: {shift:.4f}")

# Apply the shift to all trees
for tree_name in trees:
    tree = file.Get(tree_name)
    if not tree:
        print(f"Tree {tree_name} not found!")
        continue

    shifted_masses = shift_mass(tree, "MM", shift)
    print(f"Applied shift to {tree_name}")

# Clean up
file.Close()
