#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 07:47:10 trottar"
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

sys.path.append("../utility")
from utility import load_equations

##################################################################################################################################################

# Check output equations
DEBUG=False

def import_model(inp_model, arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, qq, ww, tt, eps, theta_cm, phi_cm, sig_prev_iter, weight_prev_iter, *params = args
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params
    # Load equations
    equations = load_equations(f"Q{str(q2_set).replace('.','p')}W{str(w_set).replace('.','p')}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")
    
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

    modelDict = {
        "sigL": sig_L,
        "sigT": sig_T,
        "sigLT": sig_LT,
        "sigTT": sig_TT,
    }

    sig_sep = modelDict[inp_model]

    sig_sep = sig_sep/2.0/math.pi
    
    # Convert from ub/GeV**2 to nb/GeV**2
    sig_sep = sig_sep*1e3
    
    print("Model {} = {:.4e}".format(inp_model, sig_sep))
    
    return sig_sep
