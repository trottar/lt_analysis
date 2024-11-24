#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-24 16:10:23 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TH1F, TCanvas, TLine, kRed, kGreen
from ROOT import gStyle
from array import array
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

lt = Root(os.path.realpath(__file__), "Plot_LTSep")

# Add this to all files for more dynamic pathing
USER = lt.USER  # Grab user info for file finding
HOST = lt.HOST
REPLAYPATH = lt.REPLAYPATH
UTILPATH = lt.UTILPATH
LTANAPATH = lt.LTANAPATH
ANATYPE = lt.ANATYPE
OUTPATH = lt.OUTPATH
CACHEPATH = lt.CACHEPATH

#################################################################################################################################################################

filename = f"{OUTPATH}/{ParticleType}_{runNum}_{MaxEvent}_Raw_Data.root"

# Check if the file exists
if not os.path.exists(filename):
    print(f"Error: File '{filename}' does not exist.")
    sys.exit(1)  # Exit the script with an error code

trees = [f"Uncut_{ParticleType.capitalize()}_Events",
         f"Cut_{ParticleType.capitalize()}_Events_all_noRF", f"Cut_{ParticleType.capitalize()}_Events_prompt_noRF", f"Cut_{ParticleType.capitalize()}_Events_rand_noRF",
         f"Cut_{ParticleType.capitalize()}_Events_all_RF", f"Cut_{ParticleType.capitalize()}_Events_prompt_RF", f"Cut_{ParticleType.capitalize()}_Events_rand_RF"]

if ParticleType == "kaon":
    MM_str = "lambda"
    MM_true = 1.1156  # Lambda peak
else:
    MM_str = "proton"
    MM_true = 0.938  # Proton peak

MM_min = MM_true - 0.05
MM_max = MM_true + 0.05

print(f"\n\nShifting missing mass to {MM_str} peak of {MM_true} for file {filename}...")

# Suppress canvas splash
ROOT.gROOT.SetBatch(True)

# Open ROOT file
file = TFile.Open(filename)

# Define a function for fitting a Gaussian with dynamically determined FWHM range
def fit_gaussian(hist, x_min, x_max):

    print("-" * 25)

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

    half_max = max_value * 0.75

    # Find left and right bins closest to half-max value
    left_bin = max_bin
    right_bin = max_bin
    while hist.GetBinContent(left_bin) > half_max and left_bin > 1:
        left_bin -= 1
    while hist.GetBinContent(right_bin) > half_max and right_bin < hist.GetNbinsX():
        right_bin += 1

    min_range = hist.GetBinCenter(left_bin)
    max_range = hist.GetBinCenter(right_bin)

    print(f"min_range: {min_range:.4f}")
    print(f"max_range: {max_range:.4f}")
    print("-" * 25)

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
peak_position = fit_gaussian(hist, MM_min, MM_max)[0]  # Grab Mean for shift

# Calculate the shift
shift = MM_true - peak_position
print(f"Calculated shift: {shift:.4f}\n\n")

# Create a canvas for plotting
canvas = TCanvas("canvas", "Fit and Shift", 800, 600)

# Save plots to PDF
pdf_filename = f"{OUTPATH}/{ParticleType}_{runNum}_MM_Shift.pdf"
canvas.Print(f"{pdf_filename}[")

# Disable the stats box globally
gStyle.SetOptStat(0)

# Draw the histogram and fit
hist.Draw("E1")  # Draw the histogram

# Draw the vertical line at MM_true
line = TLine(MM_true, 0, MM_true, hist.GetMaximum())
line.SetLineColor(kGreen)
line.SetLineStyle(2)
line.SetLineWidth(3)
line.Draw("same")

# Add a breakdown of MM_true, fit mean, and shift value
text = ROOT.TPaveText(0.6, 0.7, 0.9, 0.9, "NDC")  # Define the text box in normalized device coordinates
text.SetFillColor(0)  # Transparent background
text.SetBorderSize(1)  # Thin border
text.AddText(f"MM_true: {MM_true:.4f}")
text.AddText(f"Fit Mean: {peak_position:.4f}")
text.AddText(f"Shift: {shift:.4f}")
text.SetTextSize(0.03)
text.SetTextAlign(13)  # Align text to the top-left within the box
text.Draw("same")

# Save the plot to the PDF
canvas.Print(pdf_filename)

# Enable the stats box globally
gStyle.SetOptStat(1)

# Apply the shift to all trees
for tree_name in trees:
    tree = file.Get(tree_name)
    if not tree:
        print(f"Tree {tree_name} not found!")
        continue

    shifted_masses = shift_mass(tree, "MM", shift)

    # Plot shifted histogram
    hist_shifted = TH1F(f"hist_shifted_{tree_name}", f"MM_{ParticleType[0].upper()} Shifted", 200, 0.7, 1.5)
    for mass in shifted_masses:
        hist_shifted.Fill(mass)

    hist_shifted.Draw()
    line.Draw("same")
    canvas.Print(pdf_filename)

canvas.Print(f"{pdf_filename}]")

# Clean up
file.Close()

# Function to apply mass shift and update the tree in place
def apply_shift_to_tree(tree, branch_name, shift):
    # Get the original branch
    branch = tree.GetBranch(branch_name)
    if not branch:
        print(f"Branch '{branch_name}' not found!")
        return

    # Create a temporary array to hold the shifted values
    new_mass = array('f', [0.0])
    
    # Set the branch address to the new array
    branch.SetAddress(new_mass)

    # Loop over the tree and apply the shift
    for event in tree:
        original_mass = getattr(event, branch_name)
        shifted_mass = original_mass + shift
        new_mass[0] = shifted_mass
        branch.Fill()  # Update the tree with the shifted mass

# Open the ROOT file in UPDATE mode
file = TFile.Open(filename, "UPDATE")

# Apply the shift to all trees
for tree_name in trees:
    tree = file.Get(tree_name)
    if not tree:
        print(f"Tree {tree_name} not found!")
        continue

    # Apply the shift and update the tree with the new shifted masses
    apply_shift_to_tree(tree, "MM", shift)
    print(f"Applied shift to {tree_name}")

# Write the changes to the file and close it
file.Write()
file.Close()
