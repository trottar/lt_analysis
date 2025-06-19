#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-02-12 19:15:55 trottar"
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
import scipy.stats as stats
import csv
from array import array
import numpy as np
from datetime import datetime
import select
import shutil
import signal
import random
import math
import re
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
def request_yn_response(string=""):
    """
    Request a yes/no response from user with customizable output message.
    
    Args:
        string: Optional custom message to display instead of default yes/no confirmation
        
    Returns:
        bool: True for 'y', False for 'n'
    """
    while True:
        response = input("Please enter 'y' or 'n': ").lower()
        if response == 'y':
            if string == "":
                print("You selected Yes.")
            else:
                print(string)
            return True
        elif response == 'n':
            if string == "":
                print("You selected No.")
            else:
                print(string)
            return False
        else:
            print("Invalid input. Please enter 'y' for yes or 'n' for no.")

################################################################################################################################################

def process_lines(string, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Identify lines containing string
    string_lines = [i for i, line in enumerate(lines) if string in line]

    # If more than one line with string exists, mark all but the last for removal
    if len(string_lines) > 1:
        to_remove = set(string_lines[:-1])  # All except the last string line
    else:
        to_remove = set()

    # Write back the filtered lines
    with open(file_path, 'w') as file:
        for i, line in enumerate(lines):
            if i not in to_remove:
                file.write(line)

################################################################################################################################################

def run_bash_script(filename, *args):
    """
    Runs a bash script with given parameters and outputs real-time logs.
    
    :param filename: The bash script file to run.
    :param args: Additional arguments to pass to the bash script.
    """
    # Start the subprocess with args unpacked
    process = subprocess.Popen(
        ['bash', filename, *map(str, args)],  # Ensure all args are converted to strings
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1  # Line-buffered
    )

    error_detected = False

    # Monitor both stdout and stderr
    streams = [process.stdout, process.stderr]

    while True:
        readable, _, _ = select.select(streams, [], [])

        for stream in readable:
            line = stream.readline()
            if line:
                if stream is process.stdout:
                    sys.stdout.write(line)  # Print stdout to terminal
                    if "ERROR:" in line:
                        error_detected = True
                elif stream is process.stderr:
                    sys.stderr.write(line)  # Print stderr to terminal
                    if "ERROR:" in line:
                        error_detected = True

        # Break when both streams are exhausted and process ends
        if process.poll() is not None and all(stream.closed or stream.read() == '' for stream in streams):
            break

    # Ensure process completes
    process.wait()

    # Exit with error code if an error is detected
    if error_detected:
        sys.exit(2)
                
################################################################################################################################################
                
def open_root_file(root_file, option="OPEN"):
    """
    Safely open a ROOT file with error handling and cache retrieval.
    
    Args:
        root_file: Path to ROOT file
        option: File opening mode (default: "OPEN")
        
    Returns:
        TFile: Opened ROOT file object
        
    Raises:
        SystemExit: If file not found and cannot be retrieved
    """
    try:
        opened_root_file = TFile.Open(root_file, option)
        return opened_root_file
    except OSError:
        # Handle cached files separately
        if "cache" in root_file:
            print(f"ERROR: {root_file} not found. It may have been removed from cache, would you like to request retrieval from silo? (y/n)")
            if request_yn_response(string="Running jcache get..."):
                # Attempt to retrieve file from cache
                subprocess.call(f"jcache get {root_file}.replace('/lustre/expphy','')", shell=True)
                sys.exit(2)
        else:
            print(f"ERROR: {root_file} not found.")
            sys.exit(2)

################################################################################################################################################
            
def check_runs_in_effcharge(run, ParticleType, OUTPATH):
    """
    Check if a specific run number exists in analyzed ROOT files.
    
    Args:
        run: Run number to check
        ParticleType: Type of particle data
        OUTPATH: Directory path for ROOT files
        
    Returns:
        bool: True if run exists, False otherwise
    """
    if run != 0:
        root_file_path = "%s/%s_%s_-1_Raw_Data.root" % (OUTPATH, ParticleType, run)
        return os.path.exists(root_file_path)
    return False

################################################################################################################################################

def data_to_csv(file_path, column_name, new_value, run_number):
    """
    Update or create a CSV file with run data, adding or modifying columns as needed.
    
    Args:
        file_path: Path to CSV file
        column_name: Name of column to update
        new_value: Value to insert
        run_number: Run number identifier
    """
    # Check if file exists
    file_exists = os.path.isfile(file_path)
    
    if file_exists:
        # Check if file has columns
        with open(file_path, 'r', newline='') as file:
            reader = csv.DictReader(file)
            fieldnames = reader.fieldnames
            if not fieldnames:  # If no columns are found
                os.remove(file_path)  # Delete the file
                file_exists = False  # Proceed as if file doesn't exist
            else:
                data = list(reader)
    else:
        data = []
        fieldnames = ['Run Number']

    if not file_exists:
        # Initialize new file structure
        data = []
        fieldnames = ['Run Number']

    # Add new column if it doesn't exist
    if column_name not in fieldnames and column_name != 'Run Number':
        fieldnames.append(column_name)
    
    # Find or create row for this run number
    row = next((row for row in data if row['Run Number'] == str(run_number)), None)
    if row is None:
        row = {fn: '' for fn in fieldnames}
        row['Run Number'] = str(run_number)
        data.append(row)
    
    # Update value and write back to file
    row[column_name] = new_value
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

# Iterate over all bins and apply the condition to zero out bins with less than 5 events
def apply_bin_threshold(histogram, threshold):
    n_bins_x = histogram.GetNbinsX()
    n_bins_y = histogram.GetNbinsY()

    for i in range(1, n_bins_x + 1):  # Bins start from 1 to n_binsX
        for j in range(1, n_bins_y + 1):  # Bins start from 1 to n_binsY
            if histogram.GetBinContent(i, j) < threshold:
                histogram.SetBinContent(i, j, 0)  # Zero out the bin content if less than 5 events
        
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

def TH1D_to_TH2D(h1, h2, h2d_name="h2d", title="2D Histogram;X axis;Y axis", 
                 n_bins=None, x_min=None, x_max=None, y_min=None, y_max=None,
                 z_min=None, z_max=None):
    
    # Get the number of bins and ranges for each histogram
    orig_n_bins_x = h1.GetNbinsX()
    orig_n_bins_y = h2.GetNbinsX()
    
    # Determine the number of bins to use
    if n_bins is None:
        n_bins = min(orig_n_bins_x, orig_n_bins_y)
    else:
        # Check that n_bins does not exceed the original number of bins
        max_allowed_bins = min(orig_n_bins_x, orig_n_bins_y)
        if n_bins > max_allowed_bins:
            print(f"Warning: Requested {n_bins} bins, but input histograms have {orig_n_bins_x} and {orig_n_bins_y} bins.")
            print(f"Using {max_allowed_bins} bins instead.")
            n_bins = max_allowed_bins

    # Use original ranges if custom ranges are not specified
    x_min = x_min if x_min is not None else h1.GetXaxis().GetXmin()
    x_max = x_max if x_max is not None else h1.GetXaxis().GetXmax()
    y_min = y_min if y_min is not None else h2.GetXaxis().GetXmin()
    y_max = y_max if y_max is not None else h2.GetXaxis().GetXmax()
    
    # Create a 2D histogram with specified or original binning
    h2d = ROOT.TH2D(h2d_name, title, n_bins, x_min, x_max, n_bins, y_min, y_max)
    
    # Calculate bin widths for input and output histograms
    input_bin_width_x = (h1.GetXaxis().GetXmax() - h1.GetXaxis().GetXmin()) / orig_n_bins_x
    input_bin_width_y = (h2.GetXaxis().GetXmax() - h2.GetXaxis().GetXmin()) / orig_n_bins_y
    output_bin_width_x = (x_max - x_min) / n_bins
    output_bin_width_y = (y_max - y_min) / n_bins
    
    # Fill the 2D histogram
    for i in range(1, orig_n_bins_x + 1):
        for j in range(1, orig_n_bins_y + 1):
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
            
            # Scale the content and error based on the change in bin widths
            scale_factor_x = input_bin_width_x / output_bin_width_x
            scale_factor_y = input_bin_width_y / output_bin_width_y
            z_value *= scale_factor_x * scale_factor_y
            z_error *= scale_factor_x * scale_factor_y
            
            # Fill the bin
            bin_val = h2d.Fill(x_value, y_value, z_value)
            
            # Set the bin error
            h2d.SetBinError(bin_val, z_error)

    # Set the z-axis range if specified
    if z_min is not None and z_max is not None:
        h2d.GetZaxis().SetRangeUser(z_min, z_max)            
    
    return h2d

################################################################################################################################################

# Create a polar plot from a histogram with phi range from -pi to pi.
def create_polar_plot(hist, title="", r_title='|-t|', phi_title='#Phi', marker_color=1, marker_size=0.5, marker_style=20):    
    
    # Extract data from histogram
    n_points = hist.GetN()
    phi_values = hist.GetX()
    r_values = hist.GetY()

    # Convert phi from [-pi, pi] to [0, 2pi]
    phi_converted = []
    for i in range(n_points):
        phi = phi_values[i]
        if phi < 0:
            phi += 2 * math.pi
        phi_converted.append(phi*(180/math.pi))

    # Create the polar plot
    polar_plot = ROOT.TGraphPolar(n_points, np.array(phi_converted, dtype='float64'), np.array(r_values, dtype='float64'))
    polar_plot.SetMarkerColor(marker_color)
    polar_plot.SetMarkerSize(marker_size)
    polar_plot.SetMarkerStyle(marker_style)
    
    # Set titles and axes
    polar_plot.SetTitle(title)
    polar_plot.GetXaxis().SetName(phi_title)
    polar_plot.GetYaxis().SetName(r_title)
    
    return polar_plot

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
    root_file = open_root_file(file_name, "UPDATE")

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
    root_file = open_root_file(file_name, "UPDATE")

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
        sys.exit(2)

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

def remove_bad_bins(histogram, error_max=10.00):
    nbins = histogram.GetNbinsX()
    for i in range(1, nbins + 1):
        content = histogram.GetBinContent(i)
        error = histogram.GetBinError(i)
        # Set bins with negative content, or large uncertainty to zero
        '''
        if content < 0.0:
            histogram.SetBinContent(i, 0)
        if error > error_max * content:
            histogram.SetBinContent(i, 0)
            histogram.SetBinError(i, error_max)
        '''
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

def adaptive_cooling(initial_temp, iteration, max_iterations, cooling_rate=0.99, scaling=50):
    frac = iteration / float(max_iterations)
    exponent = scaling * frac
    return initial_temp * (cooling_rate ** exponent)

def sanitize_params(params, clip_min=-1e4, clip_max=1e4):
    """Clip parameters to a reasonable range to avoid huge values."""
    return [max(min(p, clip_max), clip_min) for p in params]

def simulated_annealing(param, temperature, perturbation_factor=0.1, min_scale=1e-6, clip_min=-1e4, clip_max=1e4):
    """Perturb a parameter using a factor of its scale (with a minimum scale to avoid zero) and then clip the result."""
    scale = abs(param) if abs(param) > min_scale else min_scale
    max_perturbation = scale * perturbation_factor
    perturbation = random.uniform(-max_perturbation, max_perturbation) * temperature
    new_param = param + perturbation
    # Clip new_param to keep it in a reasonable range:
    return max(min(new_param, clip_max), clip_min)

def calculate_cost(f_sig, g_sig, current_params, num_events, num_params, lambda_reg=0.01):
    """
    Calculate cost (modified reduced chi-square) with adaptive regularization.
    Includes parameter clipping to avoid huge l2_reg values.
    """
    # Clip parameters before cost calculation
    current_params = sanitize_params(current_params, clip_min=-1e4, clip_max=1e4)
    
    # Compute l2 regularization using the sanitized parameters.
    l2_reg = np.sum(np.square(current_params))
    
    lambda_min = 1e-6
    lambda_max = 100.0

    def calculate_residuals():
        residuals = []
        for i in range(num_events):
            observed = g_sig.GetY()[i]
            expected = f_sig.Eval(g_sig.GetX()[i])
            error = g_sig.GetEY()[i]
            residual = (observed - expected) / (error if error != 0 else 1.0)
            residuals.append(residual)
        return np.array(residuals)

    residuals = calculate_residuals()
    
    # If any residual is non-finite, return a very large cost.
    if not np.all(np.isfinite(residuals)):
        print("Non-finite residual detected. Parameters:", current_params)
        return 1e12, lambda_reg

    def compute_cost(lambda_val, residuals):
        mse = np.mean(np.square(residuals))
        complexity_penalty = 0.1 * num_params / num_events
        if num_events <= num_params:
            cost = (mse + lambda_val * l2_reg) / (num_events + complexity_penalty)
        else:
            chi_square = f_sig.GetChisquare()
            nu = f_sig.GetNDF()
            # Safeguard against division by very small nu:
            if nu < 1e-6:
                nu = 1e-6
            if lambda_reg == 0:
                cost = (chi_square) / nu
            else:
                cost = (chi_square + lambda_val * l2_reg) / nu
        return cost

    if num_events <= num_params:
        lambda_values = np.logspace(np.log10(lambda_min), np.log10(lambda_max), 20)
        costs = [compute_cost(lam, residuals) for lam in lambda_values]
        best_index = np.argmin(costs)
        best_lambda = lambda_values[best_index]
        best_cost = costs[best_index]
    else:
        best_lambda = lambda_reg
        best_cost = compute_cost(best_lambda, residuals)
    return best_cost, best_lambda

def acceptance_probability(old_cost, new_cost, temperature):
    """
    Return 1 if the new cost is lower.
    Otherwise, return exp(-(new_cost - old_cost)/temperature)
    to decide probabilistically.
    """
    if new_cost < old_cost:
        return 1.0
    else:
        try:
            return math.exp(-(new_cost - old_cost) / temperature)
        except OverflowError:
            return 0.0

def adjust_params(params, adjustment_factor=0.1):
    """Adjust parameters randomly by a percentage of their value."""
    return params + np.random.uniform(-adjustment_factor, adjustment_factor, size=len(params)) * params

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

def adaptive_regularization(cost_history, lambda_reg, min_improvement=1e-4, max_lambda=1.0, min_lambda=1e-6):
    """
    Dynamically adapt regularization strength based on cost history.
    """
    if len(cost_history) < 2:
        return lambda_reg
    try:
        improvement = (cost_history[-2] - cost_history[-1]) / (abs(cost_history[-2]) + 1e-10)
    except Exception:
        improvement = 0.0
    if improvement > min_improvement:
        lambda_reg *= 0.9
    elif improvement < 0:
        lambda_reg *= 1.1
    return max(min_lambda, min(lambda_reg, max_lambda))

def get_central_value(lst):
    """Return the median value from a list."""
    n = len(lst)
    if n % 2 == 1:
        return lst[n // 2]
    else:
        mid1, mid2 = n // 2 - 1, n // 2
        return (lst[mid1] + lst[mid2]) / 2

def calculate_information_criteria(n_samples, n_parameters, chi2):
    """Calculate AIC and BIC for model selection."""
    log_likelihood = -chi2 / 2  # approximate
    aic = 2 * n_parameters - 2 * log_likelihood
    bic = n_parameters * np.log(n_samples) - 2 * log_likelihood
    return aic, bic

def select_valid_parameter(param, lower_bound, upper_bound):
    """Ensure the parameter stays within the given bounds."""
    return max(lower_bound, min(param, upper_bound))

def compute_p_value(chi2, ndf):
    """Compute the p-value for a given chi-square and degrees of freedom."""
    return 1 - stats.chi2.cdf(chi2, ndf)

def log_iteration(iteration, cost, params, acceptance, temperature, logfile):
    """Log iteration details to a file for later analysis."""
    with open(logfile, 'a') as f:
        f.write(f"{iteration}\t{cost:.6f}\t{' '.join(f'{p:.6f}' for p in params)}\t{acceptance:.6f}\t{temperature:.6f}\n")

################################################################################################################################################

def check_chi_squared_values(par_chi2_vec, chi2_threshold, fit_params, equations):
    """
    Check chi-squared values where every 4 elements are identical.
    Only prints warning once per unique chi-squared value above threshold.
    
    Args:
        par_chi2_vec (list): List of chi-squared values where every 4 elements are identical
        chi2_threshold (float): Threshold value for acceptable chi-squared
        fit_params (dict): Dictionary of fit parameters
        equations: Function to find parameter wrapper
    
    Returns:
        bool: True if any chi-squared values exceed threshold, False otherwise
    """
    # Verify input length is multiple of 4
    if len(par_chi2_vec) % 4 != 0:
        raise ValueError("Input vector length must be multiple of 4")
        
    unique_bad_chi2_count = 0
    checked_values = set()  # Track unique chi-squared values we've already warned about

    bad_chi2_indices = []
    
    # Only need to check every 4th value since they repeat
    for i in range(0, len(par_chi2_vec), 4):
        chi2 = par_chi2_vec[i]
        
        # Skip if we've already warned about this chi2 value
        if chi2 in checked_values:
            continue
            
        if chi2 > chi2_threshold:
            # Get corresponding signal name and parameters
            sig_name, initial_params = list(fit_params.items())[i // 4]
            # Extract equation string from wrapper
            _, _, equation_str = find_params_wrapper(equations)(sig_name, initial_params)
            
            print(
                f"\nWARNING: Reduced Chi-Squared of {chi2:.5f} found, "
                f"which is above the threshold of {chi2_threshold}.\n\t"
                f"Increase fit iterations or adjust functional form of...{equation_str}"
            )
            
            checked_values.add(chi2)
            bad_chi2_indices.append(i // 4)
            unique_bad_chi2_count += 1
            
            if unique_bad_chi2_count >= 4:  # Limit to 4 unique warnings
                break
                
    return (unique_bad_chi2_count > 0, bad_chi2_indices)
    
################################################################################################################################################

def load_equations(filename='variables.inp'):
    equations = {}
    full_path = f"{LTANAPATH}/src/models/{filename}"
    
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"The file {full_path} does not exist.")
    
    with open(full_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line and not line.startswith('#'):  # Ignore empty lines and comments
                try:
                    key, value = line.split('=', 1)
                    equations[key.strip()] = value.strip()
                except ValueError:
                    print(f"Warning: Invalid format in line {line_num}: {line}")
    
    #print(f"Loaded {len(equations)} equations from {full_path}")
    return equations

################################################################################################################################################

# Function to extract values from the file
def extract_values(filename):
    values = []

    # Regular expression to match the pattern 'key=value' or 'key = value'
    pattern = re.compile(r'(\w+)\s*=\s*([^\s#]+)')

    # Open the file and search for matching patterns
    with open(filename, 'r') as file:
        for line in file:
            # Ignore lines that start with '#'
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            # Search for key=value pairs
            match = pattern.search(line)
            if match:
                # Append the key and value as a tuple (key, value) to the list
                key, value = match.groups()
                values.append((key, value))
    
    return values

################################################################################################################################################

def prepare_equations(equations, sig_type):
    tiny_offset = 1e-15  # Define a tiny offset to avoid division by zero

    if sig_type == "sig_L":
        eq_lst = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_T', 'sig_LT', 'sig_TT', 'wfactor')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par1, par2, par3, par4):\n"
        func_str += "        par1 = par1 if par1 > 1e-15 else par1 + tiny_offset\n"
        func_str += "        par2 = par2 if par2 > 1e-15 else par2 + tiny_offset\n"
        func_str += "        par3 = par3 if par3 > 1e-15 else par3 + tiny_offset\n"
        func_str += "        par4 = par4 if par4 > 1e-15 else par4 + tiny_offset\n"
    elif sig_type == "sig_T":
        eq_lst = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_LT', 'sig_TT', 'wfactor')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par5, par6, par7, par8):\n"
        func_str += "        par6 = par6 if par6 > 1e-15 else par5 + tiny_offset\n"
        func_str += "        par6 = par6 if par6 > 1e-15 else par6 + tiny_offset\n"
        func_str += "        par7 = par7 if par7 > 1e-15 else par7 + tiny_offset\n"
        func_str += "        par8 = par8 if par8 > 1e-15 else par8 + tiny_offset\n"
    elif sig_type == "sig_LT":
        eq_lst = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_TT', 'wfactor')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par9, par10, par11, par12):\n"
        func_str += "        par9 = par9 if par9 > 1e-15 else par9 + tiny_offset\n"
        func_str += "        par10 = par10 if par10 > 1e-15 else par10 + tiny_offset\n"
        func_str += "        par11 = par11 if par11 > 1e-15 else par11 + tiny_offset\n"
        func_str += "        par12 = par12 if par12 > 1e-15 else par12 + tiny_offset\n"
    elif sig_type == "sig_TT":
        eq_lst = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_LT', 'wfactor')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par13, par14, par15, par16):\n"
        func_str += "        par13 = par13 if par13 > 1e-15 else par13 + tiny_offset\n"
        func_str += "        par14 = par14 if par14 > 1e-15 else par14 + tiny_offset\n"
        func_str += "        par15 = par15 if par15 > 1e-15 else par15 + tiny_offset\n"
        func_str += "        par16 = par16 if par16 > 1e-15 else par16 + tiny_offset\n"
    elif sig_type == "wfactor":
        eq_lst = [f"{k} = {v}" for k, v in equations.items() if k in ('wfactor', 'mtar', 'mpipl', 'mkpl')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt):\n"
        # Need to initialize since isn't an argument
        func_str += "        theta_cm = 0.0\n"        
    else:
        print(f"ERROR: Unrecognized sig_type '{sig_type}'!")
        sys.exit(2)

    matches = list(filter(lambda e: sig_type in e, eq_lst))
    if not matches:
        print(f"ERROR: Issue with function {sig_type}! Check input model file...")
        sys.exit(2)

    # Add checks to avoid zero values
    func_str += "        qq = qq if qq > 1e-15 else qq + tiny_offset\n"
    func_str += "        ww = ww if ww > 1e-15 else ww + tiny_offset\n"
    func_str += "        tt = tt if tt > 1e-15 else tt + tiny_offset\n"
    func_str += "        theta_cm = theta_cm if theta_cm > 1e-15 else theta_cm + tiny_offset\n"

    # Build function body with equations
    func_str += "        " + "\n        ".join(eq_lst) + "\n"
    func_str += f"        return {sig_type}\n"

    exec_globals = {'__builtins__': None, 'math': math, 'tiny_offset': tiny_offset}
    exec(func_str, exec_globals)
    return exec_globals[f'{sig_type}_optimized']

##################################################################################################################################################

def find_params_wrapper(equations):
    def tmp_func(sig_type, param_vals, eqns=equations):
        return find_params(eqns, sig_type, param_vals)
    return tmp_func
        
def find_params(equations, sig_type, param_vals):
    new_param_lst = []
    num_params = 0.0
    if sig_type == "L":
        eq_str = '\n'.join([f"{k} = {v}" for k, v in equations.items() if k not in ('sig_T', 'sig_LT', 'sig_TT', 'wfactor')])
        num_params = eq_str.count('par')
        for i in range(num_params):
            new_param_lst.append(param_vals[i])
        return num_params, new_param_lst, eq_str.split(f"sig_{sig_type}")[1].replace("=", f"sig_{sig_type} = ").strip()
    if sig_type == "T":
        eq_str = '\n'.join([f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_LT', 'sig_TT', 'wfactor')])
        num_params = eq_str.count('par')
        for i in range(num_params):
            new_param_lst.append(param_vals[i])
        return num_params, new_param_lst, eq_str.split(f"sig_{sig_type}")[1].replace("=", f"sig_{sig_type} = ").strip()
    if sig_type == "LT":
        eq_str = '\n'.join([f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_TT', 'wfactor')])
        num_params = eq_str.count('par')
        for i in range(num_params):
            new_param_lst.append(param_vals[i])
        return num_params, new_param_lst, eq_str.split(f"sig_{sig_type}")[1].replace("=", f"sig_{sig_type} = ").strip()
    if sig_type == "TT":
        eq_str = '\n'.join([f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_LT', 'wfactor')])
        num_params = eq_str.count('par')
        for i in range(num_params):
            new_param_lst.append(param_vals[i])
        return num_params, new_param_lst, eq_str.split(f"sig_{sig_type}")[1].replace("=", f"sig_{sig_type} = ").strip()
    if sig_type == "wfactor":
        eq_str = '\n'.join([f"{k} = {v}" for k, v in equations.items() if k in ('wfactor', 'mtar', 'mpipl', 'mkpl')])
        return num_params, new_param_lst, eq_str.split(f"wfactor")[1].replace("=", f"wfactor = ").strip()
    else:
        print("ERROR: Invalid function request!")
        sys.exit(2)
        
##################################################################################################################################################

def select_valid_parameter(sig_name, elements):
    
    sig_dict = {
        "L" : [ f"p{v}" for v in range(1,5)],
        "T" : [ f"p{v}" for v in range(5,9)],
        "LT" : [ f"p{v}" for v in range(9,13)],
        "TT" : [ f"p{v}" for v in range(13,17)]
    }
    
    valid_params = [ f"{s}" for s, e in zip(sig_dict[sig_name], elements) if e != 0.0]
    
    while True:
        # Prompt user for input
        user_input = input(f"\n\nPlease enter parameter to fit for this iteration ({', '.join(map(str, valid_params))}):")
        # Check if input is within the valid range and the element is not zero
        if user_input in valid_params:
            print(f"\n\nParameter {user_input} selected for Sig {sig_name}...")
            return [e for s, e in zip(sig_dict[sig_name], elements) if s == user_input]
        else:
            print(f"ERROR: Invalid parameter! Please select one of the following...{', '.join(map(str, valid_params))}")

##################################################################################################################################################            

# Define a function for fitting a Gaussian with dynamically determined FWHM range
def fit_gaussian(hist_original, x_min, x_max, show_fit=True):

    hist = hist_original.Clone()
    
    #print("-" * 25)

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

    #print(f"min_range: {min_range:.4f}")
    #print(f"max_range: {max_range:.4f}")
    #print("-" * 25)
        
    if show_fit:
        hist.Fit("gaus", "Q", "", min_range, max_range)
    else:
        hist.Fit("gaus", "Q0", "", min_range, max_range)
    fit_func = hist.GetFunction('gaus')

    fit_func.SetLineColor(ROOT.kRed)

    mean = fit_func.GetParameter(1)
    mean_err = fit_func.GetParError(1)

    integral = fit_func.Integral(min_range, max_range)

    del hist

    return [mean, mean_err, integral]

##################################################################################################################################################

def adapt_limits(f, nsig=5, hard_lo=None, hard_hi=None, min_width=1e-3):
    """
    Shrink/grow parameter limits around last best value ± nsig·σ,
    but never exceed the hard envelopes hard_lo / hard_hi.
    """
    npar = f.GetNpar()
    hard_lo = hard_lo or [-np.inf]*npar
    hard_hi = hard_hi or [ np.inf]*npar

    for i in range(npar):
        p  = f.GetParameter(i)
        dp = f.GetParError(i) or 1e-6
        lo = max(hard_lo[i], p - nsig*dp)
        hi = min(hard_hi[i], p + nsig*dp)
        if hi - lo < min_width:
            lo, hi = p - min_width/2, p + min_width/2
        f.SetParLimits(i, lo, hi)

##################################################################################################################################################
