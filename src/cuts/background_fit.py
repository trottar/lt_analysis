#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-19 18:50:53 trottar"
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
    "Q2p1W2p95Right_highe" : 2.75e-2,
    "Q2p1W2p95Left_highe" : 2.00e-2,
    "Q2p1W2p95Center_highe" : 1.25e-2,
    "Q2p1W2p95Left_lowe" : 0.5e-2,
    "Q2p1W2p95Center_lowe" : 0.5e-2,
    # Q2=3p0, W=2p32
    "Q3p0W2p32Right_highe" : 50,
    "Q3p0W2p32Left_highe" : 50,
    "Q3p0W2p32Center_highe" : 50,
    "Q3p0W2p32Left_lowe" : 50,
    "Q3p0W2p32Center_lowe" : 50,
    # Q2=3p0, W=3p14
    "Q3p0W3p14Right_highe" : 1.25e-2,
    "Q3p0W3p14Left_highe" : 2.75e-2,
    "Q3p0W3p14Center_highe" : 1.00e-2,
    "Q3p0W3p14Left_lowe" : 1.00e-2,
    "Q3p0W3p14Center_lowe" : 0.75e-2,
    # Q2=4p4, W=2p74
    "Q4p4W2p74Right_highe" : 1.25e-2,
    "Q4p4W2p74Left_highe" : 1.25e-2,
    "Q4p4W2p74Center_highe" : 0.75e-2,
    "Q4p4W2p74Left_lowe" : 0.50e-2,
    "Q4p4W2p74Center_lowe" : 0.25e-2,
    # Q2=5p5, W=3p02
    "Q5p5W3p02Right_highe" : 1.75e-2,
    "Q5p5W3p02Left_highe" : 3.25e-2,
    "Q5p5W3p02Center_highe" : 1.25e-2,
    "Q5p5W3p02Left_lowe" : 0.75e-2,
    "Q5p5W3p02Center_lowe" : 0.75e-2,
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

    phi_setting = subDict["phi_setting"]
    
    bg_factor = bg_dict["Q{}W{}{}_{}e".format(Q2,W,phi_setting,EPSSET)]    

    fit_func = TF1("fit_func", "[0]", mm_min, mm_max)
    
    fit_func.FixParameter(0, bg_factor) 

    fit_func.SetLineColorAlpha(ROOT.kWhite, 0)  # Transparent color

    hist.Fit("fit_func", "Q")
    
    # Get the fitted constant value and its uncertainties
    bg_par = fit_func.GetParameter(0)
    bg_err = fit_func.GetParError(0)

    return fit_func, bg_par
