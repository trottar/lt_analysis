#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 02:30:26 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math, sys

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations

##################################################################################################################################################

def import_model(inp_model, arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    Q2, W, theta_cm, tt, qq, ww, *params = args
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params

    # Load equations
    equations = load_equations(f"Q{Q2.replace('.','p')}W{W.replace('.','p')}.model")
    
    # Evaluate equations defined in src/models
    local_vars = locals()
    for key, equation in equations.items():
        try:
            local_vars[key] = eval(equation, {"__builtins__": None, "math": math}, local_vars)
        except (ValueError, OverflowError):
            local_vars[key] = -1000.0

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
