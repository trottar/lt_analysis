#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 02:25:51 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations

##################################################################################################################################################

# Define constants
pi = 3.14159
mtar_gev = 0.93827231
mpipl=0.139570
mkpl=0.493677

def iterWeight(arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, q2_sim, w_sim, t_sim, eps_sim, thetacm_sim, phicm_sim, sigcm_sim, wt_sim, *params = args
    
    # Load equations
    equations = load_equations(f"Q{q2_set.replace('.','p')}W{w_set.replace('.','p')}.model")

    q2_gev = q2_sim # Already GeV
    t_gev = t_sim  # Already GeV, issue here!!! t_sim makes no sense
    s = w_sim**2
    s_gev = s # Already GeV

    ##############
    # HARD CODED #
    ##############
    if q2_set == 2.1:
        q2_set = 2.115
    ##############
    ##############
    ##############

    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        try:
            local_vars[key] = eval(equation, {"__builtins__": None}, local_vars)
        except OverflowError:
            local_vars[key] = -1000.0

    sigl, sigt, siglt, sigtt, wfactor = [local_vars[key] for key in ['sigl', 'sigt', 'siglt', 'sigtt', 'wfactor']]
    
    sigl = sigl*wfactor
    sigt = sigt*wfactor
    sigtt = sigtt*wfactor
    siglt = siglt*wfactor

    sig = (sigt + eps_sim * sigl + eps_sim * math.cos(2. * phicm_sim) * sigtt +
             math.sqrt(2.0 * eps_sim * (1. + eps_sim)) * math.cos(phicm_sim) * siglt)
    
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    #sig = sig / 2.0 / pi  # dsig/dtdphicm in microbarns/GeV**2/rad

    wtn = wt_sim * sig / sigcm_sim

    #print("sig",sig)
    #print("sigcm",sigcm_sim)
    #print("wtn",wtn)
    #print("wt_sim",wt_sim)
    
    return [float(wtn),float(sig)]
