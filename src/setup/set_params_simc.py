#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-01 04:21:02 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import shutil
import os, sys

Q2 = sys.argv[1].replace("p", "")
W = sys.argv[2].replace("p", "")
    
###############################################################################################################################################
# ltsep package import and pathing definitions

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
SIMCPATH=lt.SIMCPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

###############################################################################################################################################

lt_param_file = "{}/src/models/par_pl_Q{}W{}".format(LTANAPATH, Q2, W)
simc_param_file = "{}/par.pl".format(SIMCPATH)

print("\n\nUpdating {} with {}...".format(simc_param_file, lt_param_file))
# Copy and overwrite the destination file with the contents of the source file
shutil.copy2(lt_param_file, simc_param_file)
