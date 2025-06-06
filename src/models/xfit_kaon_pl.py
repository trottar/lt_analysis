
#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-03-09 12:24:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math
import sys

################################################################################################################################################
# Importing utility functions
sys.path.append("utility")
from utility import load_equations, prepare_equations

################################################################################################################################################

# Global variables
pol_str = ""
Q2 = ""
W = ""
equations = {}
fun_Sig_L_optimized = None
fun_Sig_T_optimized = None
fun_Sig_LT_optimized = None
fun_Sig_TT_optimized = None

def set_val(inp_pol_str, inp_Q2, inp_W):
    global pol_str, Q2, W, equations, fun_Sig_L_optimized, fun_Sig_T_optimized, fun_Sig_LT_optimized, fun_Sig_TT_optimized
    pol_str = inp_pol_str
    Q2 = inp_Q2
    W = inp_W

    # Load equations from model input file of given setting
    equations = load_equations(f"Q{Q2}W{W}.model")

    # Grab functional forms from model input file
    fun_Sig_L_optimized = prepare_equations(equations, 'sig_L')
    fun_Sig_T_optimized = prepare_equations(equations, 'sig_T')
    fun_Sig_LT_optimized = prepare_equations(equations, 'sig_LT')
    fun_Sig_TT_optimized = prepare_equations(equations, 'sig_TT')

################################################################################################################################################

def fun_Sig_L_wrapper(wfactor, q2, w, theta):
    def tmp_func(x, par, g=wfactor, qq=q2, ww=w, theta_cm=theta):
        return fun_Sig_L(g, qq, ww, theta_cm, x, par)
    return tmp_func

def fun_Sig_L(g, qq, ww, theta_cm, x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    par1, par2, par3, par4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180    
    # Calculate SigL
    return fun_Sig_L_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par1, par2, par3, par4)

def fun_Sig_T_wrapper(wfactor, q2, w, theta):
    def tmp_func(x, par, g=wfactor, qq=q2, ww=w, theta_cm=theta):
        return fun_Sig_T(g, qq, ww, theta_cm, x, par)
    return tmp_func

def fun_Sig_T(g, qq, ww, theta_cm, x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    par1, par2, par3, par4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180    
    # Calculate SigT
    return fun_Sig_T_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par1, par2, par3, par4)

def fun_Sig_LT_wrapper(wfactor, q2, w, theta):
    def tmp_func(x, par, g=wfactor, qq=q2, ww=w, theta_cm=theta):
        return fun_Sig_LT(g, qq, ww, theta_cm, x, par)
    return tmp_func

def fun_Sig_LT(g, qq, ww, theta_cm, x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    par9, par10, par11, par12 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180
    # Calculate SigLT
    return fun_Sig_LT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par9, par10, par11, par12)

def fun_Sig_TT_wrapper(wfactor, q2, w, theta):
    def tmp_func(x, par, g=wfactor, qq=q2, ww=w, theta_cm=theta):
        return fun_Sig_TT(g, qq, ww, theta_cm, x, par)
    return tmp_func

def fun_Sig_TT(g, qq, ww, theta_cm, x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    par13, par14, par15, par16 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180    
    # Calculate SigTT
    return fun_Sig_TT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par13, par14, par15, par16)
