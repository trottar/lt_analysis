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
    # Q2=2p1, W=2p95
    # t-range = 0.150-0.250 (MM cut)
    "Q2p1W2p95Right_highe" : 150, # Background value divided by number of events, check number of events of MM
    "Q2p1W2p95Left_highe" : 50,
    "Q2p1W2p95Center_highe" : 300,
    "Q2p1W2p95Left_lowe" : 25,
    "Q2p1W2p95Center_lowe" : 100,
    # Q2=3p0, W=2p32
    # t-range = 0.1-1.0 (no MM cut)
    #"Q3p0W2p32Right_highe" : 15/682, # Background value divided by number of events
    #"Q3p0W2p32Left_highe" : 10/1570,
    #"Q3p0W2p32Center_highe" : 15/2829,
    #"Q3p0W2p32Left_lowe" : 50/3130,
    #"Q3p0W2p32Center_lowe" : 50/3301,
    # t-range = 0.1-1.0 (MM cut)
    #"Q3p0W2p32Right_highe" : 25/517, # Background value divided by number of events
    #"Q3p0W2p32Left_highe" : 25/979,
    #"Q3p0W2p32Center_highe" : 60/1706,
    #"Q3p0W2p32Left_lowe" : 50/2508,
    #"Q3p0W2p32Center_lowe" : 50/2975,    
    # t-range = 0.45-1.0 (MM cut)
    #"Q3p0W2p32Right_highe" : 15/484, # Background value divided by number of events
    #"Q3p0W2p32Left_highe" : 10/978,
    #"Q3p0W2p32Center_highe" : 15/1594,
    #"Q3p0W2p32Left_lowe" : 50/2508,
    #"Q3p0W2p32Center_lowe" : 50/2734,    
    # Q2=3p0, W=3p14
    # t-range = 0.1-1.0 (MM cut)
    #"Q3p0W3p14Right_highe" : 950/8632, # Background value divided by number of events, check number of events of MM
    #"Q3p0W3p14Left_highe" : 700/5656,
    #"Q3p0W3p14Center_highe" : 2000/21704,
    #"Q3p0W3p14Left_lowe" : 850/10763,
    #"Q3p0W3p14Center_lowe" : 1500/12573,
    # t-range = 0.1-0.3 (MM cut)
    #"Q3p0W3p14Right_highe" : 250/4572,
    #"Q3p0W3p14Left_highe" : 100/879,
    #"Q3p0W3p14Center_highe" : 850/19952,
    #"Q3p0W3p14Left_lowe" : 100/827,
    #"Q3p0W3p14Center_lowe" : 850/9258,
    # t-range = 0.1-1.0 (MM cut)
    #"Q3p0W3p14Right_highe" : 950, # Background value divided by number of events, check number of events of MM
    #"Q3p0W3p14Left_highe" : 700,
    #"Q3p0W3p14Center_highe" : 2000,
    #"Q3p0W3p14Left_lowe" : 850,
    #"Q3p0W3p14Center_lowe" : 1500,
    # t-range = 0.1-0.3 (MM cut)
    #"Q3p0W3p14Right_highe" : 250,
    #"Q3p0W3p14Left_highe" : 100,
    #"Q3p0W3p14Center_highe" : 700,
    #"Q3p0W3p14Left_lowe" : 100,
    #"Q3p0W3p14Center_lowe" : 500,
    # t-range = 0.1-0.3 (MM cut, Center only) bg fit test
    #"Q3p0W3p14Center_highe" : 500,
    #"Q3p0W3p14Center_lowe" : 500,
    # t-range = 0.1-1.0 (MM cut), bg=500
    #"Q3p0W3p14Right_highe" : 500, # Background value divided by number of events, check number of events of MM
    #"Q3p0W3p14Left_highe" : 500,
    #"Q3p0W3p14Center_highe" : 500,
    #"Q3p0W3p14Left_lowe" : 500,
    #"Q3p0W3p14Center_lowe" : 500,
    # t-range = 0.1-1.0 (MM cut), bg=500
    #"Q3p0W3p14Right_highe" : 500, # Background value divided by number of events, check number of events of MM
    #"Q3p0W3p14Left_highe" : 500,
    #"Q3p0W3p14Center_highe" : 500,
    #"Q3p0W3p14Left_lowe" : 500,
    #"Q3p0W3p14Center_lowe" : 500,    
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
    bg_par = fit_func.GetParameter(0) * norm_factor_data
    bg_err = fit_func.GetParError(0)

    if num_evts == 0:
        return fit_func, bg_par
    else:
        # Removes function from histogram
        hist.GetFunction("fit_func").Delete()
    
        return fit_func, bg_par
