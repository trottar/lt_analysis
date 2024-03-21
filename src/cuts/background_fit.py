#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-21 17:34:15 trottar"
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
    # t-range = 0.1-1.0 (no MM cut)
    #"Q3p0W2p32Right_highe" : 15/682, # Background value divided by number of events
    #"Q3p0W2p32Left_highe" : 10/1570,
    #"Q3p0W2p32Center_highe" : 15/2829,
    #"Q3p0W2p32Left_lowe" : 50/3130,
    #"Q3p0W2p32Center_lowe" : 50/3301,
    # t-range = 0.1-1.0 (MM cut)
    "Q3p0W2p32Right_highe" : 50/517, # Background value divided by number of events
    "Q3p0W2p32Left_highe" : 1/979,
    "Q3p0W2p32Center_highe" : 50/1706,
    "Q3p0W2p32Left_lowe" : 50/2508,
    "Q3p0W2p32Center_lowe" : 50/2975,    
    # t-range = 0.45-1.0 (MM cut)
    #"Q3p0W2p32Right_highe" : 15/484, # Background value divided by number of events
    #"Q3p0W2p32Left_highe" : 10/978,
    #"Q3p0W2p32Center_highe" : 15/1594,
    #"Q3p0W2p32Left_lowe" : 50/2508,
    #"Q3p0W2p32Center_lowe" : 50/2734,    
    # Q2=3p0, W=3p14
    "Q3p0W3p14Right_highe" : 15/682, # Background value divided by number of events
    "Q3p0W3p14Left_highe" : 50/1570,
    "Q3p0W3p14Center_highe" : 50/2829,
    "Q3p0W3p14Left_lowe" : 50/3130,
    "Q3p0W3p14Center_lowe" : 50/3301,    
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

    num_evts = hist.GetEntries()

    bg_factor = bg_dict["Q{}W{}{}_{}e".format(Q2, W, phi_setting, EPSSET)]*num_evts

    fit_func = TF1("fit_func", "[0]", mm_min, mm_max)
    
    fit_func.FixParameter(0, bg_factor) 

    hist.Fit("fit_func", "Q")
    
    # Get the fitted constant value and its uncertainties
    bg_par = fit_func.GetParameter(0)
    bg_err = fit_func.GetParError(0)

    #canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    #hist.Draw()
    #canvas.SaveAs("{}/{}_Q{}W{}_{}e.png".format(LTANAPATH, hist.GetName(), Q2, W, EPSSET))

    if num_evts == 0:
        return fit_func, bg_par
    else:
        # Removes function from histogram
        hist.GetFunction("fit_func").Delete()
    
        return fit_func, bg_par
