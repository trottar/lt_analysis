#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-12 14:52:19 trottar"
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
