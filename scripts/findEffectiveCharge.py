#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-14 17:25:46 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import sys, os

################################################################################################################################################
'''
User Inputs
'''
efficiency_table = sys.argv[1]
RUNLIST = sys.argv[2].split(" ")

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH

################################################################################################################################################
# Grab and calculate efficiency 

from getEfficiencyValue import getEfficiencyValue

################################################################################################################################################

effective_charge = ""
effective_charge_uncern = ""
tot_efficiency = ""
tot_efficiency_uncern = ""
ebeam_val = ""
pTheta_val = ""

for runNum in RUNLIST:

    tot_efficiency += " " + getEfficiencyValue(runNum,efficiency_table,"efficiency")
    tot_efficiency_uncern += " " + "1"
    
    charge  = getEfficiencyValue(runNum,efficiency_table,"bcm")
    print(charge)
    # Need to convert to int value for bash to interpret correctly
    effective_charge += " " + str(int(1000*(float(charge/1000)*float(tot_efficiency))))
    effective_charge_uncern += " " + "1"
    
    ebeam_val += " " + getEfficiencyValue(runNum,efficiency_table,"ebeam")
    pTheta_val += " " + getEfficiencyValue(runNum,efficiency_table,"pTheta")    


BashInput=("({}) ({}) ({}) ({}) ({}) ({})".format(effective_charge, effective_charge_uncern, tot_efficiency, tot_efficiency_uncern, ebeam_val, pTheta_val))
print(BashInput)
