#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-04-24 11:02:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os
import ROOT

root_path = sys.argv[1]
input_file_name = sys.argv[2]
input_tree_names = sys.argv[3]
output_file_name = sys.argv[4]
string_run_nums = sys.argv[5]
particle = sys.argv[6]

###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import progres bar
from ltsep import Misc

###############################################################################################################################################

outfile = ROOT.TFile(root_path + output_file_name + ".root", "RECREATE")
if not outfile.IsOpen():
    print("ERROR: Output file {} cannot be opened. Exiting the function.".format(outfile.GetName()))
    sys.exit(1)

arr_run_nums = [int(x) for x in string_run_nums.split()]

for tree in input_tree_names.split():

    chain = ROOT.TChain(tree)

    for i,n in enumerate(arr_run_nums):
        # Progress bar
        if len(arr_run_nums) > 1:
            Misc.progressBar(i, len(arr_run_nums)-1,bar_length=25)
        else:
            Misc.progressBar(len(arr_run_nums), len(arr_run_nums),bar_length=25)
        filepath = root_path + particle + "_" + str(n) + input_file_name + ".root"
        if not os.path.isfile(filepath):
            print("WARNING: File {} not found.".format(filepath))
            continue
        tempfile = ROOT.TFile.Open(filepath)
        if tempfile == None or not tempfile.IsOpen() or tempfile.TestBit(ROOT.TFile.kRecovered):
            print("WARNING: File {} not found or not opened or corrupted.".format(filepath))
            continue
        #print("Adding {}...".format(filepath))
        chain.Add(filepath)
        
    outfile.cd()
    chain.Write(tree, ROOT.TObject.kWriteDelete)
    
    print("\n\tTree {} added to {}.root".format(tree,output_file_name))
    
outfile.Close()

