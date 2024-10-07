#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 10:16:21 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math, sys

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import load_equations, prepare_equations

##################################################################################################################################################

# Check output equations
DEBUG=False

def import_model(inp_model, arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, theta_cm, tt, qq, ww, *params = args
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params
    # Load equations
    equations = load_equations(f"Q{str(q2_set).replace('.','p')}W{str(w_set).replace('.','p')}.model")

    fun_Sig_L_optimized = prepare_equations(equations, 'sig_L')
    fun_Sig_T_optimized = prepare_equations(equations, 'sig_T')
    fun_Sig_LT_optimized = prepare_equations(equations, 'sig_LT')
    fun_Sig_TT_optimized = prepare_equations(equations, 'sig_TT')
    fun_wfactor_optimized = prepare_equations(equations, 'wfactor')

    sig_L = fun_Sig_L_optimized(q2_set, w_set, qq, ww, tt, p1, p2, p3, p4)
    sig_T = fun_Sig_T_optimized(q2_set, w_set, qq, ww, tt, p5, p6, p7, p8)
    sig_LT = fun_Sig_LT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p9, p10, p11, p12)
    sig_TT = fun_Sig_TT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p13, p14, p15, p16)
    wfactor = fun_wfactor_optimized(q2_set, w_set, qq, ww, tt)
        
    modelDict = {
        "sigL": sig_L,
        "sigT": sig_T,
        "sigLT": sig_LT,
        "sigTT": sig_TT,
    }

    sig_sep = modelDict[inp_model]

    sig_sep = sig_sep*wfactor
    
    sig_sep = sig_sep/2.0/math.pi
    
    # Convert from ub/GeV**2 to nb/GeV**2
    sig_sep = sig_sep*1e3
    
    print("Model {} = {:.4e}".format(inp_model, sig_sep))
    
    return sig_sep
