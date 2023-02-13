#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-13 15:52:06 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys
import os
import ROOT

root_path = sys.argv[1]
input_file_name = sys.argv[2]
input_tree_name = sys.argv[3]
output_file_name = sys.argv[4]
string_run_nums = sys.argv[5]

###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Misc


arr_run_nums = [int(x) for x in string_run_nums.split()]

chain = ROOT.TChain(input_tree_name)

for i,n in enumerate(arr_run_nums):
    # Progress bar
    Misc.progressBar(i, len(arr_run_nums),bar_length=25)
    filepath = root_path + str(n) + input_file_name + ".root"
    tempfile = ROOT.TFile.Open(filepath)
    if tempfile == None or not tempfile.IsOpen() or tempfile.TestBit(ROOT.TFile.kRecovered):
        print("File {} not found or not opened or corrupted. Skipping this file.".format(filepath))
        continue
    #print("Adding {}...".format(filepath))
    chain.Add(filepath)

outfile = ROOT.TFile(root_path + output_file_name + ".root", "RECREATE")
if not outfile.IsOpen():
    print("Output file {} cannot be opened. Exiting the function.".format(outfile.GetName()))

chain.Merge(outfile.GetName())

print("\n\tTree {} added to {}.root".format(input_tree_name,output_file_name))

outfile.Close()

