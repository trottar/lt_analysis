#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-28 14:45:50 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import os, sys, subprocess

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

inp_dir = sys.argv[1]

def generate_file_batches(source_dir, batch_size=100):
    """
    Generate batches of file paths to avoid exceeding command-line length limits.

    :param source_dir: Local source directory path.
    :param batch_size: Number of files per batch.
    :return: A generator yielding batches of file paths.
    """
    file_paths = []
    for root, _, files in os.walk(source_dir):
        for file in files:
            file_paths.append(os.path.join(root, file))
            if len(file_paths) >= batch_size:
                yield file_paths
                file_paths = []
    if file_paths:
        yield file_paths

def run_jput_in_batches(source_dir, dest_prefix, batch_size=100):
    """
    Run jput command in smaller batches to avoid exceeding command-line length limits.

    :param source_dir: Local source directory path.
    :param dest_prefix: Target directory prefix for the stub.
    :param batch_size: Number of files per batch.
    """
    file_batches = list(generate_file_batches(source_dir, batch_size))
    print(len(file_batches))
    for i, batch in enumerate(file_batches):
        # Progress bar
        Misc.progressBar(i, batch_size, bar_length=25)
        command = ["jput", "-r", source_dir, dest_prefix] + batch
        print("!!!!!",batch)
        try:
            #print(" ".join(command))
            #subprocess.run(command, check=True)
            subprocess.run(" ".join(command), check=True)
            print(f"Batch with {len(batch)} files uploaded successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing jput command for batch: {e}")
        except FileNotFoundError:
            print("jput command not found. Make sure it is installed and in your PATH.")
            break

# Parameters
source_directory = f"/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT/cache_transfer/kaon/Q3p0W3p14/{inp_dir}"  # Source directory
destination_prefix = f"/mss/hallc/kaonlt/trottar/kaon/Q3p0W3p14/{inp_dir}"  # Destination prefix

# Run jput command in batches
run_jput_in_batches(source_directory, destination_prefix, batch_size=5)
