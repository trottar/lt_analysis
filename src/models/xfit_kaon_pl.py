#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 09:29:49 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math
import sys
import logging

logging.basicConfig(level=logging.DEBUG)

# Importing utility functions
sys.path.append("utility")
from utility import load_equations

DEBUG = True

# Global variables
pol_str = ""
Q2 = ""
W = ""
equations = {}
fun_Sig_L_optimized = None
fun_Sig_T_optimized = None
fun_Sig_LT_optimized = None
fun_Sig_TT_optimized = None

def prepare_equations(equations, sig_type):
    if sig_type == "sig_L":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_T', 'sig_LT', 'sig_TT')]
    if sig_type == "sig_T":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_LT', 'sig_TT')]
    if sig_type == "sig_LT":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_TT')]
    if sig_type == "sig_TT":
        eq_list = [f"{k} = {v}" for k, v in equations.items() if k not in ('sig_L', 'sig_T', 'sig_LT')]
    eq_list.append(sig_type)
    print("!!!!!!!!!!",eq_list)
    
    func_str = f"def {sig_type}_optimized(tt, qq, ww, p1, p2, p3, p4):\n"
    func_str += "    " + "\n    ".join(eq_list) + "\n"
    func_str += f"    return {sig_type}"
    
    exec_globals = {'__builtins__': None, 'math': math}
    exec(func_str, exec_globals)
    return exec_globals[f'{sig_type}_optimized']

def set_val(inp_pol_str, inp_Q2, inp_W):
    global pol_str, Q2, W, equations, fun_Sig_L_optimized, fun_Sig_T_optimized, fun_Sig_LT_optimized, fun_Sig_TT_optimized
    pol_str = inp_pol_str
    Q2 = inp_Q2
    W = inp_W
    equations = load_equations(f"Q{Q2}W{W}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")

    fun_Sig_L_optimized = prepare_equations(equations, 'sig_L')
    fun_Sig_T_optimized = prepare_equations(equations, 'sig_T')
    fun_Sig_LT_optimized = prepare_equations(equations, 'sig_LT')
    fun_Sig_TT_optimized = prepare_equations(equations, 'sig_TT')

def fun_Sig_L(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    p1, p2, p3, p4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    return fun_Sig_L_optimized(tt, q2_set, w_set, p1, p2, p3, p4)

def fun_Sig_T(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    p5, p6, p7, p8 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    return fun_Sig_T_optimized(tt, q2_set, w_set, p1, p2, p3, p4)

def fun_Sig_LT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    p9, p10, p11, p12 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    return fun_Sig_LT_optimized(tt, q2_set, w_set, p1, p2, p3, p4)

def fun_Sig_TT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    p13, p14, p15, p16 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    return fun_Sig_TT_optimized(tt, q2_set, w_set, p1, p2, p3, p4)
