#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-03 12:13:19 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import subprocess

################################################################################################################################################

def show_pdf_with_evince(file_path):
    try:
        subprocess.Popen(['evince', file_path])
    except FileNotFoundError:
        print("Evince not found. Please make sure it is installed.")
    except Exception as e:
        print("An error occurred: {}".format(e))
