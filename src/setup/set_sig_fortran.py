#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 07:42:11 trottar"
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

##################################################################################################################################################
# Importing utility functions

sys.path.append(f"{LTANAPATH}/src/utility")
from utility import extract_values

###############################################################################################################################################
# Define paths for the Fortran file and the test file
file_path = f"{LTANAPATH}/src/models/xmodel_{ParticleType}_{pol_str}.f"
test_file_path = f"{LTANAPATH}/src/models/Q{Q2}W{W}.model"

# Define max line length for Fortran (fixed-format typically has a limit of 72 characters)
max_fortran_line_length = 65

# Grab all equations in models definition file
sig_var = [var[0] for var in extract_values(test_file_path)]
for sig_val in sig_var:
    print("-"*50)
    print(f"Updating {file_path} with proper {sig_val}...")

    # Step 1: Read and extract {sig_val} from test.txt
    sigl_str = None  # Initialize variable to hold the {sig_val} value
    with open(test_file_path, 'r') as test_file:
        for line in test_file:
            if '#' not in line:
                if (f'{sig_val}=' in line) or (f'{sig_val} =' in line):  # Look for the line that defines {sig_val}
                    sigl_str = line.split('=')[1].strip()  # Extract the value after '=' and strip extra spaces
                    # Convert math to fortran syntax
                    sigl_str = sigl_str.replace(f"math.pi",f"3.14159")
                    sigl_str = sigl_str.replace(f"math.exp",f"exp")
                    sigl_str = sigl_str.replace(f"math.log",f"log")
                    sigl_str = sigl_str.replace(f"math.fabs",f"abs")
                    sigl_str = sigl_str.replace(f"math.sin",f"sin")
                    # Update parameter names to match fortran vectors
                    if "sig_L" == sig_val:
                        for par in range(1,5):
                            sigl_str = sigl_str.replace(f"p{par:d}",f"par({par:d})")
                    if "sig_T" == sig_val:
                        for par in range(5,9):
                            sigl_str = sigl_str.replace(f"p{par:d}",f"par({par:d})")
                    if "sig_LT" == sig_val:
                        for par in range(9,13):
                            sigl_str = sigl_str.replace(f"p{par:d}",f"par({par:d})")
                    if "sig_TT" == sig_val:
                        for par in range(13,17):
                            sigl_str = sigl_str.replace(f"p{par:d}",f"par({par:d})")
                    break  # No need to search further once {sig_val} is found

    if sigl_str is None:
        print(f"{sig_val} not found in {test_file_path}!")
    else:
        # Step 2: Read the Fortran file and store its contents
        with open(file_path, 'r') as file:
            lines = file.readlines()

        SigSet = False
        # Step 3: Find and replace the specific variable definition in the Fortran file
        i = 0
        while i < len(lines):
            line = lines[i]
            if '#' not in line:
                if f'{sig_val}=' in line:
                    if "(" in line.split('=')[0]:
                        sigl_str = sigl_str + ')'
                    if not SigSet:
                        print(f"\tChanging {line.strip()} to {sig_val}={sigl_str}")
                        # Preserve the spaces before '{sig_val}'
                        prefix_spaces = line[:line.find(f'{sig_val}=')]
                        # Construct the new line with {sig_val}, ensuring we don't exceed the Fortran line length
                        new_line = f'{prefix_spaces}{sig_val}={sigl_str}\n'
                        # Check if the new line exceeds the maximum Fortran line length
                        if len(new_line) > max_fortran_line_length:
                            # Split the line into multiple lines with proper continuation
                            first_line = new_line[:max_fortran_line_length].rstrip() + "\n"
                            continuation_lines = []
                            remaining_line = new_line[max_fortran_line_length:]
                            # Add continuation lines, placing & (or digit) in column 6
                            while len(remaining_line) > max_fortran_line_length - 6:
                                continuation_lines.append(f"     >     {remaining_line[:max_fortran_line_length - 6]}\n")
                                remaining_line = remaining_line[max_fortran_line_length - 6:]
                            # Add the last continuation line
                            continuation_lines.append(f"     >     {remaining_line}")
                            # Replace the original line with the properly formatted multi-line version
                            lines[i] = first_line + ''.join(continuation_lines)
                        else:
                            # No need for continuation lines if the length is okay
                            lines[i] = new_line

                        # Delete the next immediate line
                        if i + 1 < len(lines):
                            del lines[i + 1]

                        # Add a blank line
                        lines.insert(i + 1, '\n')

                        # Once found, no need to replace anything else
                        SigSet = True
            i += 1
    
        # Step 4: Write the modified content back to the Fortran file
        with open(file_path, 'w') as file:
            file.writelines(lines)

        #print(f"Updated {sig_val} in {file_path}")
