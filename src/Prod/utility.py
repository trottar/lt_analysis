#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-10 21:54:50 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from array import array
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

def convert_TH1F_to_numpy(histogram):

    # Get the bin contents as a 1D NumPy array
    bin_contents = np.zeros((histogram.GetNbinsX()))
    for i in range(histogram.GetNbinsX()):
        bin_contents[i] = histogram.GetBinContent(i+1)

    return bin_contents

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

def calculate_aver_data2(hist_data, hist_dummy, t_bins):
    """
    Process histograms hist_data and hist_dummy using provided t_bins and phi_bins.

    Parameters:
    hist_data (TH1F): Histogram containing data
    hist_dummy (TH1F): Histogram containing dummy data
    t_bins (list): List of bin edges for t

    Returns:
    average_hist_data (TH1F): Histogram containing average data after processing
    """
    
    # Create histograms for storing processed data
    num_bins = len(t_bins) - 1
    bin_edges = array('d', t_bins)
    average_hist_data = ROOT.TH1F("average_hist_data", "Average Data", num_bins, bin_edges)

    for t_bin in range(1, num_bins+1):
        # Get the lower and upper edges of the current t_bin
        t_bin_low = t_bins[t_bin - 1]
        t_bin_up = t_bins[t_bin]

        # Find events in hist_data and hist_dummy within the current t_bin
        events_data = hist_data.Integral(hist_data.FindBin(t_bin_low), hist_data.FindBin(t_bin_up))
        events_dummy = hist_dummy.Integral(hist_dummy.FindBin(t_bin_low), hist_dummy.FindBin(t_bin_up))
        events_dummy_sub = events_data - events_dummy
        print("---------------",events_data, events_dummy, events_dummy_sub)
        
        # Calculate the average value for the current t_bin
        bin_width = t_bin_up - t_bin_low
        average_value = events_dummy_sub / bin_width
        print("_______________",average_value)
        
        # Fill the histogram with the calculated average value
        #average_hist_data.Fill((t_bin_low + t_bin_up) / 2, average_value)
        average_hist_data.Fill(average_value)

    return convert_TH1F_to_numpy(average_hist_data)  # Return the processed histogram as a numpy array

def calculate_aver_data(hist_data, hist_dummy, t_data, t_dummy, t_bins):

    # Bin t_data and t_dummy in t_bins
    binned_t_data = t_data.Rebin(len(t_bins)-1, "binned_t_data", array('d', t_bins))
    binned_t_dummy = t_dummy.Rebin(len(t_bins)-1, "binned_t_dummy", array('d', t_bins))
    
    # Bin hist_data and hist_dummy using the binned t_data and t_dummy
    binned_hist_data = hist_data.Rebin(len(t_bins)-1, "binned_hist_data", array('d', t_bins))
    binned_hist_dummy = hist_dummy.Rebin(len(t_bins)-1, "binned_hist_dummy", array('d', t_bins))
    
    # Subtract hist_dummy from hist_data per bin
    subtracted_hist = binned_hist_data.Clone("subtracted_hist")
    subtracted_hist.Add(binned_hist_dummy, -1)
    
    # Calculate the average per bin of the subtracted bins
    num_bins = subtracted_hist.GetNbinsX()
    averaged_values = []
    for bin_idx in range(1, num_bins+1):
        bin_content = subtracted_hist.GetBinContent(bin_idx)
        bin_width = subtracted_hist.GetBinWidth(bin_idx)
        average_value = bin_content / bin_width
        averaged_values.append(average_value)
    
    return averaged_values

################################################################################################################################################

def calculate_aver_simc(hist_data, t_bins):
    """
    Process histogram hist_data using provided t_bins and phi_bins.

    Parameters:
    hist_data (TH1F): Histogram containing data
    t_bins (list): List of bin edges for t

    Returns:
    average_hist_data (TH1F): Histogram containing average data after processing
    """

    # Create histogram for storing processed data
    average_hist_data = ROOT.TH1F("average_hist_data", "Average Data", len(t_bins)-1, array('d', t_bins))

    for t_bin in range(1, hist_data.GetNbinsX()+1):
        
        # Find events in hist_data within bins of t
        events_data = hist_data.GetBinContent(t_bin)

        # Calculate the average hist_data value per t bin
        bin_width = average_hist_data.GetXaxis().GetBinWidth(t_bin)
        average_value = events_data / bin_width
        average_hist_data.SetBinContent(t_bin, average_value)

        
    return convert_TH1F_to_numpy(average_hist_data)
