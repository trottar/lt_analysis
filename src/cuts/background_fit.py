#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-18 02:57:48 trottar"
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
    # Q2=4p4, W=2p74
    "Q4p4W2p74Right_highe" : 40, # Background value divided by number of events, check number of events of MM
    "Q4p4W2p74Left_highe" : 50,
    "Q4p4W2p74Center_highe" : 50,
    "Q4p4W2p74Left_lowe" : 150,
    "Q4p4W2p74Center_lowe" : 100,    
}

##############
##############
##############

################################################################################################################################################

def bg_fit(phi_setting, inpDict, hist):
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]

    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    norm_factor_data = inpDict["normfac_data"]

    num_evts = hist.GetEntries()

    norm_tot_evts = num_evts/inpDict["bg_tot_num_evts_{}".format(phi_setting)]

    bg_factor = bg_dict["Q{}W{}{}_{}e".format(Q2, W, phi_setting, EPSSET)]*norm_tot_evts
    # No background fit
    #bg_factor = 0.0

    fit_func = TF1("fit_func", "[0]", mm_min, mm_max)
    
    fit_func.FixParameter(0, bg_factor) 

    hist.Fit("fit_func", "Q")
    
    # Get the fitted constant value and its uncertainties
    bg_par = fit_func.GetParameter(0)
    bg_err = fit_func.GetParError(0)

    if num_evts == 0:
        return fit_func, bg_par
    else:
        # Removes function from histogram
        hist.GetFunction("fit_func").Delete()
    
        return fit_func, bg_par
