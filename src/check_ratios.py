#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-11 03:05:43 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import os, sys

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

# Read command-line arguments
inp_pid = sys.argv[1]
inp_pol = sys.argv[2]
if inp_pol == "+1":
    inp_pol = "pl"
else:
    inp_pol = "mn"
inp_Q2 = sys.argv[3]
inp_W = sys.argv[4]
inp_loeps = sys.argv[5]
inp_hieps = sys.argv[6]

aver_lo_file = '{}/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(LTANAPATH, inp_pid, inp_pol, inp_Q2.replace(".",""), inp_W.replace(".",""), float(inp_loeps)*100)
aver_hi_file = '{}/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(LTANAPATH, inp_pid, inp_pol, inp_Q2.replace(".",""), inp_W.replace(".",""), float(inp_hieps)*100)

for f in [aver_hi_file, aver_lo_file]:

    print(f"\nChecking {f}...")

    # Open the file, read the lines, and filter them
    with open(f, "r") as infile:
        lines = infile.readlines()
    print(lines)

    # Filter out lines that contain '*'
    lines = [line for line in lines if '*' not in line]

    num_lines = len(lines)
    
    # Remove the last line if it's empty or just a newline
    lines = [line.replace("\n", "") for line in lines if len(line) != num_lines]
    
    # Open the file again, this time in write mode to overwrite the content
    with open(f, "w") as outfile:
        outfile.writelines(lines)
