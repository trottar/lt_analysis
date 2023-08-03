#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-03 11:49:04 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import subprocess

################################################################################################################################################

def show_pdf_with_evince(pdf_file_path):
    try:
        subprocess.run(['evince', pdf_file_path])
    except FileNotFoundError:
        print("Evince not found. Please make sure it is installed.")
    except Exception as e:
        print("An error occurred: {}".format(e))
