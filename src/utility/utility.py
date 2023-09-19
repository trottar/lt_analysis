#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-19 15:39:26 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from array import array
import numpy as np
from datetime import datetime
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

# Flatten TH1D to NumPy array
def flatten_hist(histogram):

    # Extract the bin contents and edges
    contents = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX()+1)]
    edges = [histogram.GetXaxis().GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX()+2)]

    # Flatten the histogram
    flattened_histogram = [edge for edge, count in zip(edges, contents) for _ in range(int(count))]

    return flattened_histogram

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

# Function to check if an object is of a ROOT type
def is_hist(obj):
    return isinstance(obj, (ROOT.TH1D, ROOT.TH2D, ROOT.TGraphPolar))

################################################################################################################################################

# Function to check if an object is of a ROOT type
def is_root_obj(obj):
    return isinstance(obj, (ROOT.TH1D, ROOT.TH2D, ROOT.TGraphErrors, ROOT.TGraphPolar, ROOT.TFile, ROOT.TMultiGraph))

################################################################################################################################################

# Save histograms to root file
def hist_to_root(hist, file_name, directory_name):
    # Check if the ROOT file already exists
    root_file = ROOT.TFile.Open(file_name, "UPDATE")

    # Split the directory names
    directories = directory_name.split('/')

    # Create or navigate through the nested directories
    current_dir = root_file    
    for directory in directories:
        # Check if the directory exists
        dir_exists = bool(current_dir.GetDirectory(directory))
        if not dir_exists:
            current_dir.mkdir(directory)
        current_dir.cd(directory)
        current_dir = ROOT.gDirectory  # Update the current directory

    # Check if the histogram already exists in the file
    existing_hist = current_dir.Get(hist.GetName())
    if existing_hist:
        current_dir.Delete(hist.GetName() + ";*")  # Delete existing histogram
        
    #print("Saving {} to {}".format(hist.GetName(), file_name))
        
    # Clone the histogram since we're storing it in a directory
    cloned_hist = hist.Clone()
    cloned_hist.Write()

################################################################################################################################################

# Save histograms to root file
def hist_to_root(hist, file_name, directory_name):
    # Check if the ROOT file already exists
    root_file = ROOT.TFile.Open(file_name, "UPDATE")

    # Split the directory names
    directories = directory_name.split('/')

    # Create or navigate through the nested directories
    current_dir = root_file    
    for directory in directories:
        # Check if the directory exists
        dir_exists = bool(current_dir.GetDirectory(directory))
        if not dir_exists:
            current_dir.mkdir(directory)
        current_dir.cd(directory)
        current_dir = ROOT.gDirectory  # Update the current directory

    # Check if the histogram already exists in the file
    existing_hist = current_dir.Get(hist.GetName())
    if existing_hist:
        current_dir.Delete(hist.GetName() + ";*")  # Delete existing histogram
        
    #print("Saving {} to {}".format(hist.GetName(), file_name))
        
    # Clone the histogram since we're storing it in a directory
    cloned_hist = hist.Clone()
    cloned_hist.Write()

################################################################################################################################################    

# Used to check if object is non-serializable and converts to a list
# so that it can be saved in json file
def custom_encoder(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError("Type not serializable")

################################################################################################################################################

# Find closest date relative to current date, based off of the iteration list (see main.py for more info)
def last_iter(file_name, current_date):

    # Read formatted dates from the file
    formatted_dates = []

    with open(file_name, "r") as file:
        for line in file:
            formatted_dates.append(line.strip())

    # Function to convert formatted date to datetime object for comparison
    def convert_to_datetime(date_str):
        # Format of dates, see run_Prod_Analysis.sh
        date_format = "H%HM%MS%S_%Y%B%d"
        return datetime.strptime(date_str, date_format)

    # Find the closest date
    closest_date = min(formatted_dates, key=lambda date: abs((convert_to_datetime(date) - convert_to_datetime(current_date)).total_seconds()))

    return closest_date
    
################################################################################################################################################

# Save histograms to root file
def get_histogram(root_file, directory_name, histogram_name):

    if not root_file:
        print("Error: Unable to open file {}.".format(root_file))
        return None

    # Split the directory names
    directories = directory_name.split('/')

    # Initialize current_dir to the root of the file
    current_dir = root_file
    
    for directory in directories:
        print("Checking directory:", directory)  # Debug statement
        # Check if the directory exists
        dir_exists = bool(current_dir.GetDirectory(directory))
        if not dir_exists:
            print("Error: Unable to find directory {}.".format(directory))
            root_file.Close()
            return None
        current_dir.cd(directory)
        current_dir = ROOT.gDirectory
        histograms_in_dir = current_dir.GetListOfKeys()
        print("Histograms in directory:", [histogram.GetName() for histogram in histograms_in_dir])  # Debug statement

    # Get the histogram
    histogram = current_dir.Get(histogram_name)

    if not histogram:
        print("Error: Unable to find histogram {}.".format(histogram_name))
        root_file.Close()
        return None

    # Check the number of entries in the histogram
    print("Number of entries in {}: {}".format(histogram.GetName(),histogram.GetEntries()))  # Debug statement
    
    # Clone the histogram to avoid ownership issues
    cloned_histogram = histogram.Clone()

    # Close the ROOT file
    #root_file.Close()

    # Check the number of entries in the cloned_histogram
    print("Number of entries in the {} cloned_histogram: {}".format(cloned_histogram.GetName(),cloned_histogram.GetEntries()))  # Debug statement

    return cloned_histogram

################################################################################################################################################

# Save histograms to root file
def hist_in_dir(root_file, directory_name):

    histDict = {}
    
    if not root_file or root_file.IsZombie():
        print("Error: Unable to open file {}.".format(root_file))
        return None

    # Split the directory names
    directories = directory_name.split('/')

    # Initialize current_dir to the root of the file
    current_dir = root_file
    
    for directory in directories:
        # Check if the directory exists
        dir_exists = bool(current_dir.GetDirectory(directory))
        if not dir_exists:
            print("Error: Unable to find directory {}.".format(directory))
            root_file.Close()
            return None
        current_dir.cd(directory)
        current_dir = ROOT.gDirectory
        histograms_in_dir = current_dir.GetListOfKeys()
        for hist_key in histograms_in_dir:

            # Get the histogram
            histogram = current_dir.Get(hist_key.GetName())

            if not histogram:
                print("Error: Unable to find histogram {}.".format(hist_key.GetName()))
                root_file.Close()
                return {}

            # Clone the histogram to avoid ownership issues
            cloned_histogram = histogram.Clone()

            histDict["{}".format(histogram.GetName())] = cloned_histogram

    return histDict


################################################################################################################################################

# Run fortran script from python with variable argument
def run_fortran(fort_script, inp_val=""):

    # Define the command to compile and run the Fortran script with input
    command = 'gfortran {} -o output && ./output {}'.format(fort_script, inp_val)

    # Execute the command and capture the output
    completed_process = subprocess.Popen(command, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Return the standard output and standard error
    return completed_process.communicate()

################################################################################################################################################
