#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-23 17:51:42 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import numpy as np
import sys, os

################################################################################################################################################
'''
User Inputs
'''
try:
    runs_eff_charge = [float(q) for q in sys.argv[1].split(" ")]
    runs_eff_charge_err = [float(err) for err in sys.argv[2].split(" ")]
except ValueError:
    print("ERROR: No valid run numbers found for this setting, check efficiency table...")
    sys.exit(2)

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
tot_eff_charge = np.sum(runs_eff_charge)

# Normalized uncertainty (relative error)
tot_eff_charge_err = np.sqrt((np.sum([err**2 for err in runs_eff_charge_err])/tot_eff_charge**2))

BashInput=("{:.3f}\n{:.5f}".format(tot_eff_charge, tot_eff_charge_err))
print(BashInput)
