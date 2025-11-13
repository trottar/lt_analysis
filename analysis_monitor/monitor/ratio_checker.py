#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 14:32:45 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
from pathlib import Path
import numpy as np
import os

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def check_ratio(inpDict, iter_dir, RATIO_THRESHOLD_SPREAD):

    # Unpack input dictionary
    Q2_str = inpDict["Q2"].replace('p','')
    W_str = inpDict["W"].replace('p','')
    LOEPS_str = float(inpDict["LOEPS"].replace('p',''))*100
    HIEPS_str = float(inpDict["HIEPS"].replace('p',''))*100

    base = Path(iter_dir) / "averages"
    f_lo = base / f"aver.pl_Q{Q2_str}W{W_str}_{LOEPS_str:.0f}.dat"
    f_hi = base / f"aver.pl_Q{Q2_str}W{W_str}_{HIEPS_str:.0f}.dat"

    # ratio, dratio, phi_bin, t_bin
    lo = np.loadtxt(f_lo, comments="#", ndmin=2)
    hi = np.loadtxt(f_hi, comments="#", ndmin=2)

    lo_ratio, lo_dratio = lo[:, 0], lo[:, 1]
    hi_ratio, hi_dratio = hi[:, 0], hi[:, 1]

    # Guard against nonpositive errors
    lo_mask = lo_dratio > 0
    hi_mask = hi_dratio > 0
    if not (lo_mask.all() and hi_mask.all()):
        raise ValueError("Found nonpositive dratio values.")

    # 1) Point-by-point band check: |r-1| <= threshold
    lo_points_pass = np.all(np.abs(lo_ratio - 1.0) <= RATIO_THRESHOLD_SPREAD)
    hi_points_pass = np.all(np.abs(hi_ratio - 1.0) <= RATIO_THRESHOLD_SPREAD)

    # 2) Error-weighted average checks
    def wmean(r, dr):
        w = 1.0 / (dr ** 2)
        return np.sum(w * r) / np.sum(w)

    lo_mean = wmean(lo_ratio, lo_dratio)
    hi_mean = wmean(hi_ratio, hi_dratio)

    lo_avg_pass = abs(lo_mean - 1.0) <= RATIO_THRESHOLD_SPREAD
    hi_avg_pass = abs(hi_mean - 1.0) <= RATIO_THRESHOLD_SPREAD

    #overall = lo_points_pass and hi_points_pass and lo_avg_pass and hi_avg_pass
    overall = True

    print(f"LOEPS points: {'PASS' if lo_points_pass else 'FAIL'}")
    print(f"HIEPS points: {'PASS' if hi_points_pass else 'FAIL'}")
    print(f"LOEPS ⟨r⟩_w = {lo_mean:.6f} -> {'PASS' if lo_avg_pass else 'FAIL'}")
    print(f"HIEPS ⟨r⟩_w = {hi_mean:.6f} -> {'PASS' if hi_avg_pass else 'FAIL'}")
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}")

    CONTINUE = overall

    return CONTINUE