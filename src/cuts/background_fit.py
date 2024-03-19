#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-19 16:37:54 trottar"
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

def bg_fit(inpDict, hist):

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    # Fit the background function to the histogram
    fit_func = TF1("fit_func", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", mm_min, mm_max, 1)
    fit_func.SetParameters(40, 0.0, 0.0)  # Set initial parameters for the Gaussian
    
    hist.Fit("fit_func", "MRQ")

    # Get the fitted parameters and their uncertainties
    amplitude = fit_func.GetParameter(0)
    mean = fit_func.GetParameter(1)
    sigma = fit_func.GetParameter(2)
    
    # Get the fitted constant value and its uncertainties
    bg_par = amplitude

    return fit_func, bg_par
