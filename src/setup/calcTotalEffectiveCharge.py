#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-10 16:05:32 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import numpy and np
import sys, os

################################################################################################################################################
'''
User Inputs
'''
runs_effective_charge = sys.argv[1].split(" ")
runs_effective_charge_uncern = sys.argv[2].split(" ")

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH

################################################################################################################################################

# Sum of all the effective charge per run
tot_effective_charge = np.sum(runs_effective_charge)

# Normalized uncertainty (converted to %)
tot_effective_charge_uncern = (np.sum(runs_effective_charge_uncern**2)/tot_effective_charge**2)*100
        
BashInput=("{}\n{}".format(tot_effective_charge, tot_effective_charge_uncern))
print(BashInput)
