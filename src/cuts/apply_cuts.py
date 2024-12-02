#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-02 05:56:55 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import sys, os, math

###############################################################################################################################################
# Establish global variables so that they're not called each iteration of loop

# First, define empty variables
W = ""
Q2 = ""
EPSSET = ""
tmin = ""
tmax = ""
a1 = ""
b1 = ""
a2 = ""
b2 = ""
a3 = ""
b3 = ""
a4 = ""
b4 = ""
c0_dict = {}

# Then, set global variables which is called with arguments
def set_val(inpDict):
    
    global W, Q2, EPSSET, ParticleType
    global tmin, tmax
    global a1, b1, a2, b2, a3, b3, a4, b4
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    tmin = float(inpDict["tmin"] )
    tmax = float(inpDict["tmax"] )
    
    # Define diamond cut parameters
    a1 = float(inpDict["a1"])
    b1 = float(inpDict["b1"])
    a2 = float(inpDict["a2"])
    b2 = float(inpDict["b2"])
    a3 = float(inpDict["a3"])
    b3 = float(inpDict["b3"])
    a4 = float(inpDict["a4"])
    b4 = float(inpDict["b4"])
    
    ##############
    # HARD CODED #
    ##############

    global c0_dict
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    if ParticleType == "kaon":
        for c0, p in zip(c0_list, h_momentum_list):
            if p == 0.889:
                c0_dict["Q2p1W2p95_lowe"] = c0 # Proper value 0.888
            elif p == 0.968:
                c0_dict["Q0p5W2p40_lowe"] = c0
                c0_dict["Q3p0W3p14_lowe"] = c0 # Proper value 1.821
                c0_dict["Q5p5W3p02_lowe"] = c0 # Proper value 0.962
            elif p == 2.185:
                c0_dict["Q0p5W2p40_highe"] = c0 # Proper value 2.066
                c0_dict["Q3p0W2p32_lowe"] = c0
            elif p == 2.328:
                c0_dict["Q4p4W2p74_lowe"] = c0
            elif p == 3.266:
                c0_dict["Q5p5W3p02_highe"] = c0            
            elif p == 4.2:
                c0_dict["Q3p0W3p14_highe"] = c0 # Proper value 4.204
            elif p == 4.712:
                c0_dict["Q4p4W2p74_highe"] = c0            
            elif p == 5.292:
                c0_dict["Q2p1W2p95_highe"] = c0
            elif p == 6.59:
                c0_dict["Q3p0W2p32_highe"] = c0
    else:
        c0_dict["Q0p4W2p20_lowe"] = 0.0
        c0_dict["Q0p4W2p20_highe"] = 0.0
            
    ##############
    ##############        
    ##############
    
###############################################################################################################################################

def apply_data_cuts(evt, mm_min=0.7, mm_max=1.5):

    ##############
    # HARD CODED #
    ##############

    adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

    adj_MM = evt.MM_shift
    
    ##############
    ##############        
    ##############
    
    #CUTs Definations 
    SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
    SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

    HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
    HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

    t_RANGE =  (tmin<-evt.MandelT) & (-evt.MandelT<tmax)

    MMCUT =  (mm_min<adj_MM) & (adj_MM<mm_max)
    
    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and t_RANGE and MMCUT

    return ALLCUTS

###############################################################################################################################################

# Subtraction cuts
def apply_data_sub_cuts(evt):

    ##############
    # HARD CODED #
    ##############

    adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp
    
    ##############
    ##############        
    ##############
    
    #CUTs Definations 
    SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1)
    SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

    HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
    HMS_Acceptance = (adj_hsdelta>=-8.0) & (adj_hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

    t_RANGE =  (tmin<-evt.MandelT) & (-evt.MandelT<tmax)

    # No MM cut    
    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and t_RANGE

    return ALLCUTS

################################################################################################################################################

def apply_simc_cuts(evt, mm_min=0.7, mm_max=1.5):

    ##############
    # HARD CODED #
    ##############

    adj_missmass = evt.missmass_shift

    ##############
    ##############        
    ##############        
    
    # Define the acceptance cuts  
    SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
    HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
      
    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

    t_RANGE =  (tmin<-evt.t) & (-evt.t<tmax)

    MMCUT =  (mm_min<adj_missmass) & (adj_missmass<mm_max)
      
    ALLCUTS = HMS_Acceptance and SHMS_Acceptance and Diamond and t_RANGE and MMCUT
    
    return ALLCUTS
