#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 03:45:59 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import os, sys

Q2 = sys.argv[1]
W = sys.argv[2]
ParticleType = sys.argv[3]
POL = sys.argv[4]

if POL == "+1":
    pol_str = "pl"
else:
    pol_str = "mn"

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
# Define paths for the Fortran file and the test file
file_path = f"{LTANAPATH}/src/models/xmodel_{ParticleType}_{pol_str}.f"
test_file_path = f"{LTANAPATH}/src/models/Q{Q2}W{W}.model"

print("\n\nUpdating {} with proper sig_L...".format(file_path))

# Step 1: Read sig_L from test.txt
with open(test_file_path, 'r') as test_file:
    sigl_str = test_file.read().strip()  # Assuming the test file contains just the sig_L value

# Step 2: Read the Fortran file and store its contents
with open(file_path, 'r') as file:
    lines = file.readlines()

# Step 3: Find and replace the specific variable definition
for i, line in enumerate(lines):
    if 'sig_L=' in line:
        print(f'''
        Changing {line.strip()} to sig_L={sigl_str}
        ''')
        # Preserve the spaces before 'sig_L'
        prefix_spaces = line[:line.find('sig_L=')]
        # Replace the line with the new sig_L value
        for par in range(1,16):
            lines[i] = f'{prefix_spaces}sig_L={sigl_str}\n'.replace("p{par}", "par({par})")

# Step 4: Write the modified content back to the Fortran file
with open(file_path, 'w') as file:
    file.writelines(lines)

print(f"Updated sig_L in {file_path}")
