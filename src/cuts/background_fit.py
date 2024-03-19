#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-19 16:56:51 trottar"
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

# Define a constant background function
def constant(x, par):
    return par[0]

def bg_fit(inpDict, hist):

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Fit the background function to the histogram
    fit_func = TF1("fit_func", constant, mm_min, mm_max, 1)
    fit_func.SetParameter(0, 50)

    hist.Fit("fit_func", "Q")
    
    # Get the fitted constant value and its uncertainties
    bg_par = fit_func.GetParameter(0)
    bg_err = fit_func.GetParError(0)

    return fit_func, bg_par
