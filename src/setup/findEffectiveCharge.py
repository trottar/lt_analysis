#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-08 16:40:13 trottar"
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
kinematics = sys.argv[4]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
OUTPATH=lt.OUTPATH
ANATYPE=lt.ANATYPE
LTANAPATH=lt.LTANAPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import check_runs_in_effcharge

##################################################################################################################################################

if "None" in OUTPATH:
    OUTPATH = OUTPATH.replace("None", f"{ANATYPE}LT")

OutFilename = f"table_{ParticleType}_{kinematics}"
foutcsv = OUTPATH + "/" + OutFilename + ".csv"

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
    
    if ParticleType == "heep":

            efficiency = getEfficiencyValue(runNum,efficiency_table,foutcsv,"efficiency")
            efficiency_err = getEfficiencyValue(runNum,efficiency_table,foutcsv,"efficiency_err")
            eff_charge = getEfficiencyValue(runNum,efficiency_table,foutcsv,"eff_charge")
            eff_charge_err = getEfficiencyValue(runNum,efficiency_table,foutcsv,"eff_charge_err")

            # Run by run list of effective charge
            eff_charge_lst += " " + str(float(eff_charge))
            eff_charge_err_lst += " " + str(float(eff_charge_err))

            # Calculated run by run total efficiency and efficiency error
            efficiency_lst += " " + str(efficiency)
            efficiency_err_lst += " " + str(efficiency_err)

            ebeam_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,foutcsv,"ebeam"))
            pTheta_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,foutcsv,"pTheta"))
        
    else:
    
        # Check if run number exists in analysed root files
        if check_runs_in_effcharge(runNum, ParticleType, "{}/OUTPUT/Analysis/{}LT".format(LTANAPATH, ANATYPE)):

            efficiency = getEfficiencyValue(runNum,efficiency_table,foutcsv,"efficiency")
            efficiency_err = getEfficiencyValue(runNum,efficiency_table,foutcsv,"efficiency_err")
            eff_charge = getEfficiencyValue(runNum,efficiency_table,foutcsv,"eff_charge")
            eff_charge_err = getEfficiencyValue(runNum,efficiency_table,foutcsv,"eff_charge_err")

            # Run by run list of effective charge
            eff_charge_lst += " " + str(float(eff_charge))
            eff_charge_err_lst += " " + str(float(eff_charge_err))

            # Calculated run by run total efficiency and efficiency error
            efficiency_lst += " " + str(efficiency)
            efficiency_err_lst += " " + str(efficiency_err)

            ebeam_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,foutcsv,"ebeam"))
            pTheta_val_lst += " " + str(getEfficiencyValue(runNum,efficiency_table,foutcsv,"pTheta"))

        else:
            RUNLIST.remove(runNum)

# Convert back to string        
RUNLIST = ' '.join(map(str, RUNLIST))        
        
BashInput=("{}\n{}\n{}\n{}\n{}\n{}\n{}".format(eff_charge_lst, eff_charge_err_lst, efficiency_lst, efficiency_err_lst, pTheta_val_lst, ebeam_val_lst, RUNLIST))
print(BashInput)
