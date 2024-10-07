#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 06:15:20 trottar"
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
# Need to grab polarity and Q2 string values from xfit script

# Check output equations
DEBUG=True

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

# Function for SigL
def fun_Sig_L(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = float(q2_set)

    try:
        p1 = par[0]
    except:
        p1 = 0.0
    try:
        p2 = par[1]
    except:
        p2 = 0.0        
    try:
        p3 = par[2]
    except:
        p3 = 0.0
    try:
        p4 = par[3]
    except:
        p4 = 0.0
        
    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        if (key != 'sig_T') and (key != 'sig_LT') and (key != 'sig_TT'):
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

    f = [local_vars[key] for key in ['sig_L']]

    return f

###############################################################################################################################################
    
# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = float(q2_set)
    
    try:
        p5 = par[0]
    except:
        p5 = 0.0
    try:
        p6 = par[1]
    except:
        p6 = 0.0        
    try:
        p7 = par[2]
    except:
        p7 = 0.0
    try:
        p8 = par[3]
    except:
        p8 = 0.0

    print("!!!!!!!!!",p5,p6,p7,p8)
        
    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        if (key != 'sig_L') and (key != 'sig_LT') and (key != 'sig_TT'):
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

    f = [local_vars[key] for key in ['sig_T']]

    return f

###############################################################################################################################################

# Function for SigLT
# thetacm term is defined on function calling
def fun_Sig_LT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = float(q2_set)
    
    try:
        p9 = par[0]
    except:
        p9 = 0.0
    try:
        p10 = par[1]
    except:
        p10 = 0.0        
    try:
        p11 = par[2]
    except:
        p11 = 0.0
    try:
        p12 = par[3]
    except:
        p12 = 0.0

    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        if (key != 'sig_L') and (key != 'sig_T') and (key != 'sig_TT'):
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

    f = [local_vars[key] for key in ['sig_LT']]

    return f

###############################################################################################################################################

# Function for SigTT
# thetacm term is defined on function calling
def fun_Sig_TT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = float(q2_set)

    try:
        p13 = par[0]
    except:
        p13 = 0.0
    try:
        p14 = par[1]
    except:
        p14 = 0.0        
    try:
        p15 = par[2]
    except:
        p15 = 0.0
    try:
        p16 = par[3]
    except:
        p16 = 0.0

    # Evaluate equations
    local_vars = locals()
    for key, equation in equations.items():
        if (key != 'sig_L') and (key != 'sig_T') and (key != 'sig_LT'):
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

    f = [local_vars[key] for key in ['sig_TT']]

    return f

###############################################################################################################################################
