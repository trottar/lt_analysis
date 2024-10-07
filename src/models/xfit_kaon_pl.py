#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 05:44:02 trottar"
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

# First, define empty strings
pol_str = ""
q2_set = ""
equations = ""

# Then, set global variables which is called with arguments defined in xfit script
def set_val(inp_pol_str, inp_q2_set):
    global pol_str, q2_set, equations
    pol_str = inp_pol_str
    q2_set = inp_q2_set
    # Load equations
    equations = load_equations(f"Q{str(q2_set).replace('.','p')}W{str(w_set).replace('.','p')}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")
        
###############################################################################################################################################

# Function for SigL
def fun_Sig_L(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    
    p1, p2, p3, p4 = par
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

    f = [local_vars[key] for key in ['sig_L']]

    return f

###############################################################################################################################################
    
# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    
    p5, p6, p7, p8 = par
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

    f = [local_vars[key] for key in ['sig_T']]

    return f

###############################################################################################################################################

# Function for SigLT
# thetacm term is defined on function calling
def fun_Sig_LT(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    
    p9, p10, p11, p12 = par
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

    f = [local_vars[key] for key in ['sig_LT']]

    return f

###############################################################################################################################################

# Function for SigTT
# thetacm term is defined on function calling
def fun_Sig_TT(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    
    p13, p14, p15, p16 = par
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

    f = [local_vars[key] for key in ['sig_TT']]

    return f

###############################################################################################################################################
