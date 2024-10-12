#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-12 03:24:53 trottar"
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
    
def fun_Sig_L(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    p1, p2, p3, p4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Calculate SigL
    return fun_Sig_L_optimized(q2_set, w_set, qq, ww, tt, p1, p2, p3, p4)

def fun_Sig_T(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    p5, p6, p7, p8 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Calculate SigT
    return fun_Sig_T_optimized(q2_set, w_set, qq, ww, tt, p5, p6, p7, p8)

def fun_Sig_LT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    # Sine term called separately so setting to 1.0
    theta_cm = math.pi/2
    p9, p10, p11, p12 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Calculate SigLT
    return fun_Sig_LT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p9, p10, p11, p12)

def fun_Sig_TT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    # Sine term called separately so setting to 1.0
    theta_cm = math.pi/2
    p13, p14, p15, p16 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    # Calculate SigTT
    return fun_Sig_TT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p13, p14, p15, p16)
