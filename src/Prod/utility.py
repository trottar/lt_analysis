#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-10 17:53:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
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

    # Calculate the total integral of the histogram (integral up to the last bin)
    total_integral = histogram.Integral()

    # Calculate the weights for each bin based on their integrals
    bin_weights = [integral / total_integral for integral in bin_integrals]

    # Weight the bin edges by the bin weights
    weighted_bin_edges = [edge * weight for edge, weight in zip(bin_edges, bin_weights)]
    
    return weighted_bin_edges

################################################################################################################################################

def calculate_aver_data(H_histogram_DATA, H_histogram_DUMMY, t_bins, phi_bins):
    
    # Create histograms to store the sum and count for each t-phi bin
    sum_histogram = ROOT.TH2F("sum_histogram", "Sum of histogram Values", len(t_bins) - 1, t_bins, len(phi_bins) - 1, phi_bins)
    count_histogram = ROOT.TH2F("count_histogram", "Number of Entries", len(t_bins) - 1, t_bins, len(phi_bins) - 1, phi_bins)

    # Loop over histogram histogram bins and fill the sum and count histograms
    for t_bin in range(1, len(t_bins)):
        for phi_bin in range(1, len(phi_bins)):
            t_low = t_bins[t_bin - 1]
            t_high = t_bins[t_bin]
            phi_low = phi_bins[phi_bin - 1]
            phi_high = phi_bins[phi_bin]

            hist_sum = 0.0
            hist_count = 0

            for bin_x in range(1, H_histogram_DATA.GetNbinsX() + 1):
                t_value = H_histogram_DATA.GetXaxis().GetBinCenter(bin_x)
                phi_value = H_histogram_DATA.GetYaxis().GetBinCenter(bin_y)

                if t_low <= t_value < t_high and phi_low <= phi_value < phi_high:
                    hist_sum += (H_histogram_DATA.GetBinContent(bin_x, bin_y) - H_histogram_DUMMY.GetBinContent(bin_x, bin_y))
                    hist_count += 1

            sum_histogram.SetBinContent(t_bin, phi_bin, hist_sum)
            count_histogram.SetBinContent(t_bin, phi_bin, hist_count)

    # Calculate the average histogram value within each t-phi bin
    histogram_aver = sum_histogram.Clone("histogram_aver")
    histogram_aver.Divide(count_histogram)

    return histogram_aver

################################################################################################################################################

def calculate_aver_simc(H_histogram_SIMC, t_bins, phi_bins):
    
    # Create histograms to store the sum and count for each t-phi bin
    sum_histogram = ROOT.TH2F("sum_histogram", "Sum of histogram Values", len(t_bins) - 1, t_bins, len(phi_bins) - 1, phi_bins)
    count_histogram = ROOT.TH2F("count_histogram", "Number of Entries", len(t_bins) - 1, t_bins, len(phi_bins) - 1, phi_bins)

    # Loop over histogram histogram bins and fill the sum and count histograms
    for t_bin in range(1, len(t_bins)):
        for phi_bin in range(1, len(phi_bins)):
            t_low = t_bins[t_bin - 1]
            t_high = t_bins[t_bin]
            phi_low = phi_bins[phi_bin - 1]
            phi_high = phi_bins[phi_bin]

            hist_sum = 0.0
            hist_count = 0

            for bin_x in range(1, H_histogram_SIMC.GetNbinsX() + 1):
                t_value = H_histogram_SIMC.GetXaxis().GetBinCenter(bin_x)
                phi_value = H_histogram_SIMC.GetYaxis().GetBinCenter(bin_y)

                if t_low <= t_value < t_high and phi_low <= phi_value < phi_high:
                    hist_sum += (H_histogram_SIMC.GetBinContent(bin_x, bin_y))
                    hist_count += 1

            sum_histogram.SetBinContent(t_bin, phi_bin, hist_sum)
            count_histogram.SetBinContent(t_bin, phi_bin, hist_count)

    # Calculate the average histogram value within each t-phi bin
    histogram_aver = sum_histogram.Clone("histogram_aver")
    histogram_aver.Divide(count_histogram)

    return histogram_aver
