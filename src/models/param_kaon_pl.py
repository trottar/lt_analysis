#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 02:50:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math, sys
import logging

logging.basicConfig(level=logging.DEBUG)

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations

##################################################################################################################################################

def iterWeight(arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, qq, ww, tt, eps, theta_cm, phi_cm, sig_prev_iter, wtt, *params = args

    ss = ww**2
    
    # Load equations
    equations = load_equations(f"Q{str(q2_set).replace('.','p')}W{str(w_set).replace('.','p')}.model")
    logging.debug(f"Loaded equations: {equations}")

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
            logging.debug(f"Evaluating equation for {key}: {equation}")
            local_vars[key] = eval(equation, {"__builtins__": None, "math": math}, local_vars)
            logging.debug(f"Result for {key}: {local_vars[key]}")
        except OverflowError:
            logging.warning(f"OverflowError for {key}, setting to -1000.0")
            local_vars[key] = -1000.0
        except Exception as e:
            logging.error(f"Error evaluating equation for {key}: {equation}")
            logging.error(f"Error message: {str(e)}")
            logging.error(f"Local variables: {local_vars}")
            raise
        
    sigl, sigt, siglt, sigtt, wfactor = [local_vars[key] for key in ['sigl', 'sigt', 'siglt', 'sigtt', 'wfactor']]
    
    sigl = sigl*wfactor
    sigt = sigt*wfactor
    sigtt = sigtt*wfactor
    siglt = siglt*wfactor

    sig = (sigt + eps * sigl + eps * math.cos(2. * phi_cm) * sigtt +
             math.sqrt(2.0 * eps * (1. + eps)) * math.cos(phi_cm) * siglt)
    
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    #sig = sig / 2.0 / pi  # dsig/dtdphicm in microbarns/GeV**2/rad

    wtn = wtt * sig / sig_prev_iter

    #print("sig",sig)
    #print("sigcm",sig_prev_iter)
    #print("wtn",wtn)
    #print("wtt",wtt)
    
    return [float(wtn),float(sig)]
