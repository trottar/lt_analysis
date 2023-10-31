#! /usr/bin/python

#
# Description: Adapted from fortran code wt28_3.f
# ================================================================
# Time-stamp: "2023-10-31 14:22:26 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
from array import array
import sys, math, os, subprocess

# Import the dynamic script
import importlib.util

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import run_fortran

##################################################################################################################################################
# Importing param model for weight iteration

sys.path.append("models")
from param_active import iterWeight    

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

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
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def iter_weight(param_file, simc_root, inpDict, phi_setting):

    
    formatted_date  = inpDict["formatted_date"]
    ParticleType  = inpDict["ParticleType"]
    Q2 = inpDict["Q2"].replace("p","")
    POL = inpDict["POL"]
    
    if int(POL) == 1:
        pol_str = "pl"
    elif int(POL) == -1:
        pol_str = "mn"
    else:
        print("ERROR: Invalid polarity...must be +1 or -1")
        sys.exit(2)
    
    # Define diamond cut parameters
    a1 = inpDict["a1"]
    b1 = inpDict["b1"]
    a2 = inpDict["a2"]
    b2 = inpDict["b2"]
    a3 = inpDict["a3"]
    b3 = inpDict["b3"]
    a4 = inpDict["a4"]
    b4 = inpDict["b4"]
    
    param_arr = []
    with open(param_file, 'r') as f:
        for i, line in enumerate(f):
            columns = line.split()
            param_arr.append(str(columns[0]))    

    if not os.path.isfile(simc_root):
        print("\n\nERROR: No simc file found called {}\n\n".format(simc_root))

    InFile_SIMC = TFile.Open(simc_root, "OPEN")
    TBRANCH_SIMC  = InFile_SIMC.Get("h10")
    Weight_SIMC  = TBRANCH_SIMC.GetBranch("Weight")
    sig_SIMC  = TBRANCH_SIMC.GetBranch("sigcm")

    # Create a new ROOT file for writing
    new_InFile_SIMC = ROOT.TFile.Open(simc_root.replace(".root","_new.root"), "RECREATE")

    # Clone the TTree from the original file
    new_TBRANCH_SIMC = TBRANCH_SIMC.CloneTree(-1, "fast")

    # Get the Weight branch from the new tree
    new_Weight_SIMC = new_TBRANCH_SIMC.GetBranch("Weight")

    # Get the sig branch from the new tree
    new_sig_SIMC = new_TBRANCH_SIMC.GetBranch("sigcm")    
    
    # Create a new branch with the updated values
    iweight = array('d', [0])  # Assuming 'd' is the data type, change if needed
    new_weight_branch = new_TBRANCH_SIMC.Branch("iter_weight", iweight, "iter_weight/D")  # 'D' for double, change if needed
    isig = array('d', [0])  # Assuming 'd' is the data type, change if needed
    new_sig_branch = new_TBRANCH_SIMC.Branch("iter_sig", isig, "iter_sig/D")  # 'D' for double, change if needed
    
    ################################################################################################################################################
    # Run over simc root branch to determine new weight

    print("\nRecalculating weight for %s simc..." % phi_setting)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

      TBRANCH_SIMC.GetEntry(i)

      # phipq == phicm
      #inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetacm, evt.phipq, evt.sigcm, evt.Weight)+' '.join(param_arr)
      inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetapq, evt.phipq, evt.sigcm, evt.Weight)+' '.join(param_arr)
      #print("-"*25,"\n",i,"\n",inp_param)

      iter_lst = iterWeight(inp_param)

      print(iter_lst)

      # Set the value of iweight
      iweight[0] = iter_lst[0]
      new_Weight_SIMC.SetAddress(iweight)

      # Fill the new branch with the new value for this entry
      new_weight_branch.Fill()
      
      # Set the value of isig
      isig[0] = iter_lst[1]
      new_sig_SIMC.SetAddress(isig)
      
      # Fill the new branch with the new value for this entry
      new_sig_branch.Fill()
          
    new_TBRANCH_SIMC.Write()
    
    new_InFile_SIMC.Close()
    InFile_SIMC.Close()
