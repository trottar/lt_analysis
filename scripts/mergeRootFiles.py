#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-13 14:46:10 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys
import os
import ROOT

def mergeRootFiles(root_path, input_file_name, input_tree_name, output_file_name, string_run_nums):
    arr_run_nums = [int(x) for x in string_run_nums.split()]

    chain = ROOT.TChain(input_tree_name)

    for n in arr_run_nums:
        filepath = root_path + str(n) + input_file_name + ".root"
        tempfile = ROOT.TFile.Open(filepath)
        if tempfile == None or not tempfile.IsOpen() or tempfile.TestBit(ROOT.TFile.kRecovered):
            print(f"File {filepath} not found or not opened or corrupted. Skipping this file.")
            continue
        chain.Add(filepath)

    outfile = ROOT.TFile(root_path + output_file_name + ".root", "RECREATE")
    if not outfile.IsOpen():
        print(f"Output file {outfile.GetName()} cannot be opened. Exiting the function.")
        return 0

    chain.Merge(outfile.GetName())

    outfile.Close()

    return 0
