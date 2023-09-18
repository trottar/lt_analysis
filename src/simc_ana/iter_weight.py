#! /usr/bin/python

#
# Description: Adapted from fortran code wt28_3.f
# ================================================================
# Time-stamp: "2023-09-18 00:31:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

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
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def iter_weight(param_file):
    '''
    # Fortran script converted to python
    
    # Define constants
    pi = 3.14159
    mtar_gev = 0.93827231
    
    # Define parameters
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    
    if abs(q2_set - 245) < 1:
        p1 =  0.25961E+02 
        p2 = -0.10000E+02 
        p3 = -0.15838E+02 
        p4 =  0.00000E+00 
        p5 =  0.46859E+02 
        p6 = -0.30000E+02 
        p7 = -0.33572E+01 
        p8 =  0.00000E+00 
        p9 =  0.10000E+04 
        p10 = -0.28000E+02 
        p11 =  0.35000E+01 
        p12 = -0.67276E+02 
    else:
        print('wtn: q2 error', q2_set)
        return

    # Parameterization based upon Fpi-1 pi+ IT25, 12.04.18
    # Revised for IT21, 12.11.06
    q2_gev = float(q2_set) / 100.0
    tav = (0.0735 + 0.028 * math.log(q2_gev)) * q2_gev
    ftav = (ti - tav) / tav
    ft = ti / ((ti + 0.139570**2)**2)

    nsigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (ti - 0.2))
    nsigt = p5 + p6 * math.log(q2_gev) + (p7 + p8 * math.log(q2_gev)) * ftav

    nsiglt = (p9 * math.exp(p10 * ti) + p11 / ti) * math.sin(thetacmi)
    nsigtt = (p12 * q2_gev * math.exp(-q2_gev)) * ft * math.sin(thetacmi)**2

    nsig219 = (nsigt + epsiloni * nsigl + epsiloni * math.cos(2.0 * phicmi) * nsigtt
              + math.sqrt(2.0 * epsiloni * (1.0 + epsiloni)) * math.cos(phicmi) * nsiglt) / 1.0

    wfactor = 1.0 / (Wcmi**2 - mtar_gev**2)**2
    nsig = nsig219 * wfactor

    nsig = nsig / (2.0 * pi * 1.0E+06)  # dsig/dtdphicm in microbarns/MeV**2/rad

    wtn = Weight * nsig / dsigdt

    wtn_limit = 0.20
    if wtn < wtn_limit and wtn > 0.0:
        continue
    else:
        wtn = 0.0

    return
    '''

    paramDict = {}
    with open(param_file, 'r') as f:
        for i, line in enumerate(f):
            columns = line.split()
            paramDict["p{}".format(i)] = str(columns)

    for key,val in enumerate(paramDict.items()):
        print("{} = {}".format(key,val))
