#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-24 14:29:11 trottar"
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
ParticleType = sys.argv[3]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
LTANAPATH=lt.LTANAPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import check_runs_in_effcharge

################################################################################################################################################
# Grab and calculate efficiency 

from getEfficiencyValue import getEfficiencyValue

################################################################################################################################################

efficiency_lst = ""
efficiency_err_lst = ""
eff_charge_lst = ""
eff_charge_err_lst = ""
ebeam_val_lst = ""
pTheta_val_lst = ""

for runNum in RUNLIST:

    # Check if run number exists in analysed root files
    if check_runs_in_effcharge(runNum, ParticleType, "{}/OUTPUT/Analysis/{}LT".format(LTANAPATH, ANATYPE)):
    
        efficiency = getEfficiencyValue(runNum,efficiency_table,"efficiency")
        efficiency_err = getEfficiencyValue(runNum,efficiency_table,"efficiency_err")
        eff_charge = getEfficiencyValue(runNum,efficiency_table,"eff_charge")
        eff_charge_err = getEfficiencyValue(runNum,efficiency_table,"eff_charge_err")

        # Run by run list of effective charge
        eff_charge_lst += " " + str(float(eff_charge))
        eff_charge_err_lst += " " + str(float(eff_charge_err))

        # Calculated run by run total efficiency and efficiency error
        efficiency_lst += " " + str(efficiency)
        efficiency_err_lst += " " + str(efficiency_err)

        ebeam_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,"ebeam"))
        pTheta_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,"pTheta"))
        
    else:
        RUNLIST.remove(runNum)

# Convert back to string        
RUNLIST = ' '.join(map(str, RUNLIST))        
        
BashInput=("{}\n{}\n{}\n{}\n{}\n{}\n{}".format(eff_charge_lst, eff_charge_err_lst, efficiency_lst, efficiency_err_lst, pTheta_val_lst, ebeam_val_lst, RUNLIST))
print(BashInput)
