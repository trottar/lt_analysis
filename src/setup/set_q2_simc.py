#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-07 16:16:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import os, sys

Q2 = sys.argv[1].replace("p",".")

##############
# HARD CODED #
##############
if "2.1" in Q2:
    # True value of this Q2
    Q2 = 2.115
##############
##############
##############

    
###############################################################################################################################################
# ltsep package import and pathing definitions

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
SIMCPATH=lt.SIMCPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

###############################################################################################################################################

file_path = "{}/physics_iterate.f".format(SIMCPATH)

print("\n\nUpdating {} with proper Q2...".format(file_path))
# Read the file and store its contents
with open(file_path, 'r') as file:
    lines = file.readlines()

# Find and replace the specific variable definition
for i, line in enumerate(lines):
    if 'q2_set=' in line:
        print('''
        Changing {} to q2_set={}
        '''.format(line.strip(),Q2))
        # Preserve the spaces before 'q2_set'
        prefix_spaces = line[:line.find('q2_set=')]
        lines[i] = '{}q2_set={}\n'.format(prefix_spaces,Q2)

# Write the modified content back to the file
with open(file_path, 'w') as file:
    file.writelines(lines)
