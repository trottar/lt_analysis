#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-19 19:33:18 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
##################################################################################################################################################
# Import relevant packages

import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
from ROOT import TF1
import sys, os, math

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

##############
# HARD CODED #
##############

bg_dict ={
    # Q2=3p0, W=2p32
    "Q3p0W2p32_highe" : 50,
    "Q3p0W2p32_lowe" : 50,
}

##############
##############
##############

################################################################################################################################################

def bg_fit(inpDict, hist):
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    bg_factor = bg_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]    

    fit_func = TF1("fit_func", "[0]", mm_min, mm_max)
    
    fit_func.FixParameter(0, bg_factor) 

    hist.Fit("fit_func", "Q")
    
    # Get the fitted constant value and its uncertainties
    bg_par = fit_func.GetParameter(0)
    bg_err = fit_func.GetParError(0)

    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    hist.Draw()
    
    hist.GetFunction("fit_func").Delete()  # Delete any previous fit function
    
    return fit_func, bg_par
