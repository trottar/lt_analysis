#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-14 17:23:18 trottar"
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

def getEfficiencyValue(runNum,efficiency_table,table_val):

    if table_val == "efficiency":
        ################################################################################################################################################
        # Grab and calculate efficiency 

        from getDataTable import calculate_efficiency

        tot_efficiency = calculate_efficiency(runNum,efficiency_table)

        ################################################################################################################################################

        return tot_efficiency

    if table_val == "bcm":
        ################################################################################################################################################
        # Grab beam energy

        from getDataTable import get_bcm

        bcm_val = get_bcm(runNum,efficiency_table)

        ################################################################################################################################################

        return bcm_val
    
    if table_val == "ebeam":
        ################################################################################################################################################
        # Grab beam energy

        from getDataTable import get_ebeam

        ebeam_val = get_ebeam(runNum,efficiency_table)

        ################################################################################################################################################

        return ebeam_val

    if table_val == "pTheta":
        ################################################################################################################################################
        # Grab pTheta

        from getDataTable import get_pTheta

        pTheta_val = get_pTheta(runNum,efficiency_table)

        ################################################################################################################################################

        return pTheta_val

