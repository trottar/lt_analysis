#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-30 14:12:23 trottar"
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
    par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12, par13, par14, par15, par16 = params

    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180
    
    # Load equations from model input file of given setting
    equations = load_equations(f"Q{str(q2_set).replace('.','p')}W{str(w_set).replace('.','p')}.model")

    # Grab functional forms from model input file
    fun_Sig_L_optimized = prepare_equations(equations, 'sig_L')
    fun_Sig_T_optimized = prepare_equations(equations, 'sig_T')
    fun_Sig_LT_optimized = prepare_equations(equations, 'sig_LT')
    fun_Sig_TT_optimized = prepare_equations(equations, 'sig_TT')
    fun_wfactor_optimized = prepare_equations(equations, 'wfactor')

    # Calculate SigL, SigT, SigLT, SigTT
    sig_L = fun_Sig_L_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par1, par2, par3, par4)
    sig_T = fun_Sig_T_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par5, par6, par7, par8)
    sig_LT = fun_Sig_LT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par9, par10, par11, par12)
    sig_TT = fun_Sig_TT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par13, par14, par15, par16)
    wfactor = fun_wfactor_optimized(q2_set, w_set, qq, ww, tt)
        
    modelDict = {
        "sigL": sig_L,
        "sigT": sig_T,
        "sigLT": sig_LT,
        "sigTT": sig_TT,
    }

    sig_sep = modelDict[inp_model]

    sig_sep = sig_sep*wfactor
    
    #sig_sep = sig_sep/2.0/math.pi

    sig_sep = sig_sep / 2.0 / math.pi / 1e6
        
    print("Model {} = {:.4e}".format(inp_model, sig_sep))
    
    return sig_sep
