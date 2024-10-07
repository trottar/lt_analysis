#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 07:03:22 trottar"
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

###############################################################################################################################################
# Need to grab polarity Q2 and W string values from xfit script

# Check output equations
DEBUG=False

# First, define empty strings
pol_str = ""
Q2 = ""
W = ""
equations = ""

# Then, set global variables which is called with arguments defined in xfit script
def set_val(inp_pol_str, inp_Q2, inp_W):
    global pol_str, Q2, W, equations
    pol_str = inp_pol_str
    Q2 = inp_Q2
    W = inp_W
    # Load equations
    equations = load_equations(f"Q{Q2}W{W}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")
        
###############################################################################################################################################

def iterWeight(arg_str):
    
    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, qq, ww, tt, eps, theta_cm, phi_cm, sig_prev_iter, weight_prev_iter, *params = args
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params
        
    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        try:
            if DEBUG:
                logging.debug(f"Evaluating equation for {key}: {equation}")
            local_vars[key] = eval(equation, {"__builtins__": None, "math": math}, local_vars)
            if DEBUG:
                logging.debug(f"Result for {key}: {local_vars[key]}")
        except OverflowError:
            logging.warning(f"OverflowError for {key}, setting to -1000.0")
            local_vars[key] = -1000.0
        except Exception as e:
            logging.error(f"Error evaluating equation for {key}: {equation}")
            logging.error(f"Error message: {str(e)}")
            logging.error(f"Local variables: {local_vars}")
            raise
        
    sig_L, sig_T, sig_LT, sig_TT, wfactor = [local_vars[key] for key in ['sig_L', 'sig_T', 'sig_LT', 'sig_TT', 'wfactor']]
    
    sig_L = sig_L*wfactor
    sig_T = sig_T*wfactor
    sig_TT = sig_TT*wfactor
    sig_LT = sig_LT*wfactor

    sig = (sig_T + eps * sig_L + eps * math.cos(2. * phi_cm) * sig_TT +
             math.sqrt(2.0 * eps * (1. + eps)) * math.cos(phi_cm) * sig_LT)
    
    sig = sig / 2.0 / math.pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    #sig = sig / 2.0 / math.pi  # dsig/dtdphicm in microbarns/GeV**2/rad

    wtn = weight_prev_iter * (sig / sig_prev_iter)

    #print("sig",sig)
    #print("sigcm",sig_prev_iter)
    #print("wtn",wtn)
    #print("weight_prev_iter",weight_prev_iter)
    
    return [float(wtn),float(sig)]
