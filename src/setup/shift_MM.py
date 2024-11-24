#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-24 15:20:52 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import ROOT
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

trees = [f"Uncut_{ParticleType.upper()}_Events", \
         f"Cut_{ParticleType.upper()}_Events_all_noRF", f"Cut_{ParticleType.upper()}_Events_prompt_noRF", f"Cut_{ParticleType.upper()}_Events_rand_noRF", \
         f"Cut_{ParticleType.upper()}_Events_all_RF", f"Cut_{ParticleType.upper()}_Events_prompt_RF", f"Cut_{ParticleType.upper()}_Events_rand_RF"]

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
file = ROOT.TFile.Open(filename)

# Function to apply mass shift
def shift_mass(tree, branch_name, shift):
    new_values = []
    for event in tree:
        original_mass = getattr(event, branch_name)
        shifted_mass = original_mass + shift
        new_values.append(shifted_mass)
    return new_values

# Fit the peak for the reference tree
reference_tree_name = f"Cut_{ParticleType.upper()}_Events_prompt_noRF"
print(f"Fitting MM_{ParticleType[0].upper()} peak for {reference_tree_name}...")
reference_tree = file.Get(reference_tree_name)

# Define histogram for fitting
hist = ROOT.TH1F("hist", f"MM_{ParticleType[0].upper()}", 200, 0.7, 1.5)  # Adjust binning/range as needed
reference_tree.Draw("MM>>hist")

# Fit the histogram to find the peak
gaus = ROOT.TF1("gaus", "gaus", MM_min, MM_max)  # Adjust range if necessary
hist.Fit(gaus, "Q")
peak_position = gaus.GetParameter(1)  # Extract the mean of the Gaussian fit

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
