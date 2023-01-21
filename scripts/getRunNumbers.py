#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-01-21 12:57:05 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os

inp_f = sys.argv[1]

if not os.path.exists(f_list):
    print("{} not found...")
else:
    # Read in run numbers from run list and save as a string where each run is separated by a space (format for bash input)
    f_data = ''
    with open(inp_f, 'r') as f:
        for line in f:
            f_data += str(line)+' '
    print(f_data)
