#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 10:13:22 trottar"
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
from utility import load_equations

##################################################################################################################################################

def prepare_equations(equations, sig_type):
    if sig_type == "sig_L":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_T', 'sig_LT', 'sig_TT')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, p1, p2, p3, p4):\n"
    if sig_type == "sig_T":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_LT', 'sig_TT')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, p5, p6, p7, p8):\n"
    if sig_type == "sig_LT":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_TT')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p9, p10, p11, p12):\n"
    if sig_type == "sig_TT":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_LT')]
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p13, p14, p15, p16):\n"
    if sig_type == "wfactor":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k in 'wfactor']
        func_str = f"def {sig_type}_optimized(q2_set, w_set, qq, ww, tt):\n"        
    
    func_str += "    " + "\n    ".join(eq_list) + "\n"
    func_str += f"    return {sig_type}"

    print("!!!!!!!!!!",func_str)
    
    exec_globals = {'__builtins__': None, 'math': math}
    exec(func_str, exec_globals)
    return exec_globals[f'{sig_type}_optimized']

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
