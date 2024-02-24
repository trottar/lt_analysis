#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-24 15:26:08 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import os, sys

W = sys.argv[1].replace("p",".")
    
###############################################################################################################################################
# ltsep package import and pathing definitions

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
SIMCPATH=lt.SIMCPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

###############################################################################################################################################

file_path = "{}/physics_iterate.f".format(SIMCPATH)

print("\n\nUpdating {} with proper W...".format(file_path))
# Read the file and store its contents
with open(file_path, 'r') as file:
    lines = file.readlines()

# Find and replace the specific variable definition
for i, line in enumerate(lines):
    if 'w_set=' in line:
        print('''
        Changing {} to w_set={}
        '''.format(line.strip(),W))
        # Preserve the spaces before 'w_set'
        prefix_spaces = line[:line.find('w_set=')]
        lines[i] = '{}w_set={}\n'.format(prefix_spaces,W)

# Write the modified content back to the file
with open(file_path, 'w') as file:
    file.writelines(lines)
