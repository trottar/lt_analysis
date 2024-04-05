#! /usr/bin/python

#
# Description: Adapted from fortran code wt28_3.f
# ================================================================
# Time-stamp: "2024-04-05 19:11:47 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
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

def iter_weight(param_file, simc_root, inpDict, phi_setting):
    
    formatted_date  = inpDict["formatted_date"]
    iter_num = inpDict["iter_num"]
    ParticleType  = inpDict["ParticleType"]
    Q2 = inpDict["Q2"].replace("p",".")
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


        
    if iter_num > 1:

        import root_numpy as rnp
        
        InFile_SIMC = TFile.Open(simc_root, "READ")
        TBRANCH_SIMC  = InFile_SIMC.Get("h10")

        Weight_root = rnp.tree2array(TBRANCH_SIMC, branches=["Weight"])
        Weight_array = np.array(Weight_root) # Save as array
        sigcm_root = rnp.tree2array(TBRANCH_SIMC, branches=["sigcm"])
        sigcm_array = np.array(sigcm_root) # Save as array

        print("!!!!!!", Weight_root)
        print("!!!!!!", sigcm_root)
        sys.exit(2)
        
        # Create a new ROOT file for writing
        new_InFile_SIMC = TFile.Open(simc_root.replace("iter_{}".format(iter_num-1),"iter_{}".format(iter_num)), "RECREATE")
        #new_TBRANCH_SIMC = ROOT.TTree("h10", "Iteration {}".format(iter_num))
        new_TBRANCH_SIMC = rnp.array2tree([], tree_name="h10")

        # Create a new branch with the updated values
        iter_weight = ROOT.Double(0)  # Assuming iter branch is of type float
        new_iter_weight = new_TBRANCH_SIMC.Branch("iter_weight", iter_weight, "iter_weight/D")
        iter_sig = ROOT.Double(0)  # Assuming iter branch is of type float
        new_iter_sig = new_TBRANCH_SIMC.Branch("iter_sig", iter_sig, "iter_sig/D")
        
    else:

        InFile_SIMC = TFile.Open(simc_root, "OPEN")
        TBRANCH_SIMC  = InFile_SIMC.Get("h10")
        
        Weight_SIMC  = TBRANCH_SIMC.GetBranch("Weight")
        sig_SIMC  = TBRANCH_SIMC.GetBranch("sigcm")        
        # Create a new ROOT file for writing
        new_InFile_SIMC = TFile.Open(simc_root.replace(".root","_iter_{}.root".format(iter_num)), "RECREATE")
        
        # Clone the TTree from the original file
        new_TBRANCH_SIMC = TBRANCH_SIMC.CloneTree(-1, "fast")

        # Get the Weight branch from the new tree
        new_Weight_SIMC = new_TBRANCH_SIMC.GetBranch("Weight")
        # Get the sig branch from the new tree
        new_sig_SIMC = new_TBRANCH_SIMC.GetBranch("sigcm")
                
        # Create a new branch with the updated values
        iweight = array('f', [0])  # Assuming 'f' is the data type, change if needed
        new_weight_branch = new_TBRANCH_SIMC.Branch("iter_weight", iweight, "iter_weight/F")  # 'f' for float, change if needed
        isig = array('f', [0])  # Assuming 'f' is the data type, change if needed
        new_sig_branch = new_TBRANCH_SIMC.Branch("iter_sig", isig, "iter_sig/F")  # 'f' for float, change if needed
    
    ################################################################################################################################################
    # Run over simc root branch to determine new weight

    print("\nRecalculating weight for %s simc..." % phi_setting)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

      TBRANCH_SIMC.GetEntry(i)

      if iter_num > 1:

          #evt.Weight = evt.iter_weight
          #evt.sigcm = evt.iter_sig

          # Note: ti is used instead of t, ti = main%t which matches its calculation in simc
          #       while t is calculated in recon_hcana (but should be invariant?? Not sure the issue)
          #       This goes for Q2i, Wi, and phiqpi as well
          #inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #print("-"*25,"\n",i,"\n",inp_param)
          inp_param = '{} {} {} {} {} {} {} {} {} '\
                      .format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.iter_sig, evt.iter_weight)+' '.join(param_arr)
          #.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #inp_param = '{} {} {} {} {} {} {} {} {} '\
              #.format(Q2, evt.Q2, evt.W, evt.t, evt.epscm, evt.thetacm, evt.phipq, evt.sigcm, evt.Weight)+' '.join(param_arr)

          iter_lst = iterWeight(inp_param)
          
          iter_weight = iter_lst[0]
          new_iter_weight.Fill()
          iter_sig = iter_lst[1]
          new_iter_sig.Fill()
          
      else:

          # Note: ti is used instead of t, ti = main%t which matches its calculation in simc
          #       while t is calculated in recon_hcana (but should be invariant?? Not sure the issue)
          #       This goes for Q2i, Wi, and phiqpi as well
          #inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #print("-"*25,"\n",i,"\n",inp_param)
          inp_param = '{} {} {} {} {} {} {} {} {} '\
                      .format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #inp_param = '{} {} {} {} {} {} {} {} {} '\
              #.format(Q2, evt.Q2, evt.W, evt.t, evt.epscm, evt.thetacm, evt.phipq, evt.sigcm, evt.Weight)+' '.join(param_arr)

          iter_lst = iterWeight(inp_param)
          
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
