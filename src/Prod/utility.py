#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-04 10:34:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import subprocess

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

# Convert TH1F to NumPy array
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

    print("^^^^^^^^^^^^^^^^^",np.average(bin_edges),"^^^^^^^^^^^^^^^^^")

    # Calculate the total integral of the histogram (integral up to the last bin)
    total_integral = histogram.Integral()

    # Calculate the weights for each bin based on their integrals
    bin_weights = [integral / total_integral for integral in bin_integrals]

    # Weight the bin edges by the bin weights
    weighted_bin_edges = [edge * weight for edge, weight in zip(bin_edges, bin_weights)]
    
    print("¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬",np.average(weighted_bin_edges),"¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬")
    
    return weighted_bin_edges
