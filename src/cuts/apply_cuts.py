#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-29 13:52:55 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

###############################################################################################################################################
# Establish global variables so that they're not called each iteration of loop

# First, define empty strings
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

# Then, set global variables which is called with arguments
def set_val(inpDict):
    
    global W, Q2, EPSSET
    global tmin, tmax
    global a1, b1, a2, b2, a3, b3
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSSET = inpDict["EPSSET"]

    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    
    # Define diamond cut parameters
    a1 = inpDict["a1"]
    b1 = inpDict["b1"]
    a2 = inpDict["a2"]
    b2 = inpDict["b2"]
    a3 = inpDict["a3"]
    b3 = inpDict["b3"]
    a4 = inpDict["a4"]
    b4 = inpDict["b4"]
    
###############################################################################################################################################
    
def apply_data_cuts(evt, mm_min=0.7, mm_max=1.5):
    
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

    MMCUT =  (mm_min<evt.MM) & (evt.MM<mm_max)
    
    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and t_RANGE and MMCUT

    return ALLCUTS

################################################################################################################################################

def apply_simc_cuts(evt, mm_min=0.7, mm_max=1.5):

    # Define the acceptance cuts  
    SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
    HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
      
    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

    t_RANGE =  (tmin<-evt.t) & (-evt.t<tmax)

    MMCUT =  (mm_min<evt.missmass) & (evt.missmass<mm_max)
      
    ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and t_RANGE and MMCUT
    
    return ALLCUTS
