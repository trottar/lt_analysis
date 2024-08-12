#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-08-12 15:45:07 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import sys,os

###############################################################################################################################################

'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH

def getEfficiencyValue(runNum,efficiency_table,foutcsv,table_val):

    if table_val == "efficiency":
        ################################################################################################################################################
        # Grab and calculate efficiency 

        from getDataTable import calculate_efficiency

        tot_efficiency = calculate_efficiency(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return tot_efficiency

    if table_val == "efficiency_err":
        ################################################################################################################################################
        # Grab and calculate efficiency_err 

        from getDataTable import calculate_efficiency_err

        tot_efficiency_err = calculate_efficiency_err(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return tot_efficiency_err

    if table_val == "eff_charge":
        ################################################################################################################################################
        # Grab and calculate eff_charge 

        from getDataTable import calculate_eff_charge

        eff_charge = calculate_eff_charge(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return eff_charge

    if table_val == "eff_charge_err":
        ################################################################################################################################################
        # Grab and calculate eff_charge_err 

        from getDataTable import calculate_eff_charge_err

        eff_charge_err = calculate_eff_charge_err(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return eff_charge_err
    
    if table_val == "bcm":
        ################################################################################################################################################
        # Grab beam energy

        from getDataTable import get_bcm

        bcm_val = get_bcm(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return bcm_val
    
    if table_val == "ebeam":
        ################################################################################################################################################
        # Grab beam energy

        from getDataTable import get_ebeam

        ebeam_val = get_ebeam(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return ebeam_val

    if table_val == "pTheta":
        ################################################################################################################################################
        # Grab pTheta

        from getDataTable import get_pTheta

        pTheta_val = get_pTheta(runNum,efficiency_table,foutcsv)

        ################################################################################################################################################

        return pTheta_val

