#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-08-27 16:29:05 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TGraph, TGraphErrors, TCanvas
from ROOT import TF1, TFitResultPtr
from ROOT import Math
import csv
from array import array
import numpy as np
from datetime import datetime
import shutil
import signal
import random
import math
import sys, os, subprocess

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

################################################################################################################################################

# Checks if run number if found in analysed root files
def check_runs_in_effcharge(run, ParticleType, OUTPATH):

    if run != 0:
        root_file_path = "%s/%s_%s_-1_Raw_Data.root" % (OUTPATH, ParticleType, run)
        if not os.path.exists(root_file_path):
            return False
        else:
            return True
    else:
        return False

################################################################################################################################################
    
# Checks if run number if found in analysed root files
def check_runs_in_main(OUTPATH, phiset, inpDict):

    ParticleType = inpDict["ParticleType"]
    runs = inpDict["runNum{}".format(phiset)].split(' ')
    efficiencies = inpDict["InData_efficiency_{}".format(phiset.lower())].split(' ')
    for run, eff in zip(runs, efficiencies):
        print("!!!!!!",run, eff)
        if int(run) != 0:
            root_file_path = "%s/%s_%s_%s_Raw_Data.root" % (OUTPATH, ParticleType, run, -1)
            if not os.path.exists(root_file_path):
                print("\n\nRun number {} not found in {}\n\t Removing...".format(run, root_file_path))
                runs.remove(run)
                efficiencies.remove(eff)
                inpDict["runNum{}".format(phiset)] = ' '.join(map(str, runs))
                inpDict["InData_efficiency_{}".format(phiset.lower())] = ' '.join(map(str, efficiencies))
        else:
            print("No {} phi setting found...".format(phiset.lower()))
            
    return inpDict    

################################################################################################################################################

def data_to_csv(file_path, column_name, new_value, run_number):
    # Check if file exists
    file_exists = os.path.isfile(file_path)
    
    if file_exists:
        # Read existing data
        with open(file_path, 'r', newline='') as file:
            reader = csv.DictReader(file)
            data = list(reader)
            fieldnames = reader.fieldnames
    else:
        # Initialize empty data and fieldnames
        data = []
        fieldnames = ['Run Number']
    
    # Add new column name if it doesn't exist
    if column_name not in fieldnames and column_name != 'Run Number':
        fieldnames.append(column_name)
    
    # Find or create the row for this run number
    row = next((row for row in data if row['Run Number'] == str(run_number)), None)
    if row is None:
        row = {fn: '' for fn in fieldnames}
        row['Run Number'] = str(run_number)
        data.append(row)
    
    # Update the value in the correct column
    row[column_name] = new_value
    
    # Write updated data back to file
    with open(file_path, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

################################################################################################################################################

def show_pdf_with_evince(file_path):

    def signal_handler(sig, frame):
        print("\n\n\tCTRL+C detected...Exiting analysis!")
        sys.exit(2)
        
    # Set up the signal handler
    signal.signal(signal.SIGINT, signal_handler)
    
    try:
        while True:            
            user_input = input("\nDo you want to open {}? (y/n): ".format(file_path))
            
            if user_input.lower() == 'y':
                process = subprocess.Popen(['evince', file_path])
                print("\n\n\nPress CTRL+C to exit or close window to continue...")
                process.wait()  # Pauses the script until Evince is closed

                break
            elif user_input.lower() == 'c':
                print("File closed...")
                break
            elif user_input.lower() == 'n':
                break
            elif user_input.lower() == 'q':
                print("Quitting...")
                sys.exit(2)
            else:
                print("Invalid input. Please enter 'y' to open or 'n'/'c' to continue or 'q' to quit.")
                continue
    
    except FileNotFoundError:
        print("Evince not found. Please make sure it is installed.")
    except Exception as e:
        print("An error occurred: {}".format(e))
    except KeyboardInterrupt:
        # Handle Ctrl+C gracefully
        print("\n\n\tCTRL+C detected...Exiting analysis!")
        process.terminate()
        sys.exit(2)        
        
################################################################################################################################################

def write_to_file(f_out,line,write_mode='a'):
    # Open a file in append mode
    with open(f_out, write_mode) as f:
        # Write the value of the variable to the file
        #print("Writing {} to {}...".format(line, f_out))
        f.write(line)

################################################################################################################################################

def replace_line(file_path, line_number, new_line):
    # Read the file and store its lines in a list
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Check if the line number is valid
    if 1 <= line_number <= len(lines):
        # Replace the specified line
        lines[line_number - 1] = new_line

        # Write the modified lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)
    else:
        print("ERROR: Invalid line number {} in file {}".format(line_number, file_path))
        sys.exit(2)

################################################################################################################################################

# Create a new dir specified by argument
def create_dir(dir_name):
    # Check if dir exists
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

################################################################################################################################################        
        
# Create a new file specified by argument        
def create_file(file_path):
    
    # Check if the file already exists
    if not os.path.exists(file_path):
        with open(file_path, 'w') as f:
            pass  # Creates an empty file
        
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

def plot1DAs2D(h1, h2, h2d_name="h2d", title="2D Histogram;X axis;Y axis"):
    # Get the number of bins for each histogram
    n_bins_x = h1.GetNbinsX()
    n_bins_y = h2.GetNbinsX()

    # Get the bin edges directly from the TH1Ds
    x_bins = np.array([h1.GetXaxis().GetBinLowEdge(i) for i in range(1, n_bins_x + 2)])
    y_bins = np.array([h2.GetXaxis().GetBinLowEdge(i) for i in range(1, n_bins_y + 2)])

    # Create a 2D histogram
    h2d = ROOT.TH2D(h2d_name, title, n_bins_x, x_bins, n_bins_y, y_bins)

    # Fill the 2D histogram
    for i in range(1, n_bins_x + 1):
        for j in range(1, n_bins_y + 1):
            x_value = h1.GetBinCenter(i)
            y_value = h2.GetBinCenter(j)
            
            # Get bin contents and errors
            z_value1 = h1.GetBinContent(i)
            z_value2 = h2.GetBinContent(j)
            z_error1 = h1.GetBinError(i)
            z_error2 = h2.GetBinError(j)
            
            # Calculate the new bin content and error
            z_value = z_value1 * z_value2

            # Calculate the new bin error, avoiding divide by zero
            if z_value1 != 0 and z_value2 != 0:
                rel_error1 = z_error1 / z_value1
                rel_error2 = z_error2 / z_value2
                z_error = z_value * np.sqrt(rel_error1**2 + rel_error2**2)
            elif z_value1 == 0 and z_value2 == 0:
                z_error = 0
            elif z_value1 == 0:
                z_error = z_value2 * z_error1
            else:  # z_value2 == 0
                z_error = z_value1 * z_error2            
            
            # Fill the bin
            bin_val = h2d.Fill(x_value, y_value, z_value)
            
            # Set the bin error
            h2d.SetBinError(bin_val, z_error)

    return h2d

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

    # Sort the dictionary keys in numerical order
    sorted_match_dict = sorted(match_dict.keys())
    # Create a new dictionary with sorted keys
    sorted_match_dict = {key: match_dict[key] for key in sorted_match_dict}

    # Convert the dictionary to the desired format
    return [[match, np.array(values)] for match, values in sorted_match_dict.items()]

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
        date_format = "%Y%B%d_H%HM%MS%S"
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
        #print("Checking directory:", directory)  # Debug statement
        # Check if the directory exists
        dir_exists = bool(current_dir.GetDirectory(directory))
        if not dir_exists:
            print("Error: Unable to find directory {}.".format(directory))
            root_file.Close()
            return None
        current_dir.cd(directory)
        current_dir = ROOT.gDirectory
        histograms_in_dir = current_dir.GetListOfKeys()
        #print("Histograms in directory:", [histogram.GetName() for histogram in histograms_in_dir])  # Debug statement

    # Get the histogram
    histogram = current_dir.Get(histogram_name)

    if not histogram:
        print("Error: Unable to find histogram {}.".format(histogram_name))
        root_file.Close()
        return None

    # Check the number of entries in the histogram
    #print("Number of entries in {}: {}".format(histogram.GetName(),histogram.GetEntries()))  # Debug statement
    
    # Clone the histogram to avoid ownership issues
    cloned_histogram = histogram.Clone()

    # Close the ROOT file
    #root_file.Close()

    # Check the number of entries in the cloned_histogram
    #print("Number of entries in the {} cloned_histogram: {}".format(cloned_histogram.GetName(),cloned_histogram.GetEntries()))  # Debug statement

    return cloned_histogram

################################################################################################################################################

# Save histograms to root file
def hist_in_dir(root_file, directory_name):

    histDict = {}
    
    if not root_file or root_file.IsZombie():
        print("Error: Unable to open file {}.".format(root_file))
        return {}

    # Initialize current_dir to the root of the file
    current_dir = root_file
    
    # Check if the directory exists
    dir_exists = bool(current_dir.GetDirectory(directory_name))
    if not dir_exists:
        print("Error: Unable to find directory {}.".format(directory_name))
        return {}
    
    current_dir.cd(directory_name)
    current_dir = ROOT.gDirectory
    histograms_in_dir = current_dir.GetListOfKeys()
    
    for hist_key in histograms_in_dir:
        # Get the TObject associated with the key
        obj = hist_key.ReadObj()

        # Get the histogram
        histogram = current_dir.Get(hist_key.GetName())

        if not histogram:
            print("Error: Unable to find histogram {}.".format(hist_key.GetName()))
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

def set_dynamic_axis_ranges(inp_str, histlist, range_factor="Default", hist_type = "DATA"):
    
    x_axis_range_factor = 0.005
    
    # Check if a custom range factor is provided
    if range_factor != "Default":
        x_axis_range_factor = float(range_factor)
        try:
            # Attempt to convert range_factor to float
            x_axis_range_factor = float(range_factor)
        except ValueError as e:
            # Handle the case where conversion fails with a specific error message
            print("Error: Unable to convert to float. {}".format(e))

    min_values = []
    max_values = []

    # Iterate through the histograms in the list
    for i, hist in enumerate(histlist):
        # Access the histogram using the specified input string and histogram type
        histogram = hist["H_{}_{}".format(inp_str, hist_type)]

        # Get the number of bins
        num_bins = histogram.GetNbinsX()

        # Initialize variables to track non-empty bins
        non_empty_bins = []

        # Find non-empty bins along the x-axis
        for bin_idx in range(1, num_bins + 1):
            bin_content = histogram.GetBinContent(bin_idx)
            if bin_content > 0:
                non_empty_bins.append(bin_idx)

        if not non_empty_bins:
            # Handle the case where all bins are empty
            print("Warning: All bins are empty for histogram {}".format(i))
            continue

        # Set x-axis range dynamically based on non-empty bins
        x_axis_min = histogram.GetBinLowEdge(min(non_empty_bins)) - histogram.GetBinLowEdge(min(non_empty_bins)) * x_axis_range_factor
        x_axis_max = histogram.GetBinLowEdge(max(non_empty_bins) + 1) + histogram.GetBinLowEdge(max(non_empty_bins) + 1) * x_axis_range_factor

        min_values.append(x_axis_min)
        max_values.append(x_axis_max)

    # Calculate the average minimum and maximum values
    avg_min = np.average(min_values)
    avg_max = np.average(max_values)

    for i, hist in enumerate(histlist):    
        # Set the x-axis range for the current histogram
        hist["H_{}_{}".format(inp_str, hist_type)].SetBins(200, avg_min, avg_max)
    
    return avg_min, avg_max
    
################################################################################################################################################    

def notify_email(email_address):

    # Email notification
    email_notify = "swif2 notify LTSep_$USER -when done -email {}".format(email_address)
    #email_notify = "swif2 notify LTSep_$USER -when done"
    
    # Run the bash command
    subprocess.call(email_notify, shell=True)

################################################################################################################################################        

# Remove negative bins
def remove_negative_bins(histogram):
    nbins = histogram.GetNbinsX()
    for i in range(1, nbins + 1):
        content = histogram.GetBinContent(i)
        if content < 0.0:
            histogram.SetBinContent(i, 0)

    return histogram
            
################################################################################################################################################        

# Define a function for fitting a Gaussian with dynamically determined FWHM range
def get_centroid(hist, x_min, x_max):
    
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
    
    half_max = max_value*0.75
    
    # Find left and right bins closest to half-max value
    left_bin = max_bin
    right_bin = max_bin
    while hist.GetBinContent(left_bin) > half_max and left_bin > 1:
        left_bin -= 1
    while hist.GetBinContent(right_bin) > half_max and right_bin < hist.GetNbinsX():
        right_bin += 1

    min_range = hist.GetBinCenter(left_bin)
    max_range = hist.GetBinCenter(right_bin)

    hist.Fit("gaus", "Q", "", min_range, max_range)
    fit_func = hist.GetFunction('gaus')

    try:
        fit_func.SetLineColor(ROOT.kRed)
        mean = fit_func.GetParameter(1)
        mean_err = fit_func.GetParError(1)
    except ReferenceError:
        pass
        mean = -1000.0
        mean_err = -1000.0
        
    return [mean, mean_err]

################################################################################################################################################

def adaptive_cooling(initial_temp, iteration, max_iterations):
    return initial_temp * (1 - iteration / max_iterations)

def simulated_annealing(param, temperature, perturbation_factor=0.1):
    # Perturbation factor determines the maximum percentage change
    max_perturbation = abs(param) * perturbation_factor
    perturbation = random.uniform(-max_perturbation, max_perturbation) * temperature
    return param + perturbation

def acceptance_probability(old_cost, new_cost, temperature):
    # Calculate the probability of accepting a worse solution
    if new_cost < old_cost:
        return 1.0
    elif temperature == 0:
        return 0.0
    else:
        return math.exp((old_cost - new_cost) / temperature)
    
def adjust_params(params, adjustment_factor=0.1):
    return params + np.random.uniform(-adjustment_factor, adjustment_factor, size=len(params)) * params

################################################################################################################################################

# Create a PyROOT callable object
class PyFunc:
    def __call__(self, par):
        return chi2_func(par)

# RedHat7    
#minimizer = Math.Factory.CreateMinimizer("Minuit2", "Migrad")
# Alma9
minimizer = Math.Factory.CreateMinimizer("Minuit", "Migrad")
minimizer.SetMaxFunctionCalls(1000000)
minimizer.SetMaxIterations(100000)
minimizer.SetTolerance(0.001)
minimizer.SetPrintLevel(0)

def local_search(params, inp_func, num_params):

    if num_params+1 > 2:
        # Create a wrapper function that can be called by the minimizer
        def chi2_func(par):
            for i in range(num_params+1):
                inp_func.SetParameter(i, par[i])
            return inp_func.GetChisquare()

        py_func = PyFunc()

        # Create the functor
        func = Math.Functor(py_func, num_params+1)  # num_params+1 is the number of parameters
        minimizer.SetFunction(func)

        # Set initial values and step sizes
        for i, param in enumerate(params):
            minimizer.SetVariable(i, "p{}".format(i), param, 0.01 * abs(param))

        # Perform the minimization
        minimizer.Minimize()

        # Get the improved parameters
        improved_params = [minimizer.X()[i] for i in range(num_params+1)]

        minimizer.Delete()
        func.Delete()

        return improved_params

    else:

        # Create a wrapper function that can be called by the minimizer
        def chi2_func(par):
            inp_func.SetParameter(1, par)
            return inp_func.GetChisquare()

        py_func = PyFunc()

        # Create the functor
        func = Math.Functor(py_func, 1)  # 1 is the number of parameters
        minimizer.SetFunction(func)

        # Set initial values and step sizes
        minimizer.SetVariable(1, "p1", param, 0.01 * abs(param))

        # Perform the minimization
        minimizer.Minimize()

        # Get the improved parameters
        improved_params = minimizer.X()

        minimizer.Delete()
        func.Delete()
        
        return improved_params            
            
################################################################################################################################################
