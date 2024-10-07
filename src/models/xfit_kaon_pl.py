#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-07 09:08:09 trottar"
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

# Prepare equations for numexpr
def prepare_equation(eq):
    # Replace Python math functions with numpy equivalents
    for func in dir(math):
        if callable(getattr(math, func)) and func in eq:
            eq = eq.replace(f"math.{func}", f"np.{func}")
    return eq


###############################################################################################################################################
# Need to grab polarity Q2 and W string values from xfit script

# Check output equations
DEBUG=False

# First, define empty strings
pol_str = ""
Q2 = ""
W = ""
equations = ""
fun_Sig_L_optimized = ""

# Then, set global variables which is called with arguments defined in xfit script
def set_val(inp_pol_str, inp_Q2, inp_W):
    global pol_str, Q2, W, equations, fun_Sig_L_optimized
    pol_str = inp_pol_str
    Q2 = inp_Q2
    W = inp_W
    # Load equations
    equations = load_equations(f"Q{Q2}W{W}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")
        
    prepared_equations = {k: prepare_equation(v) for k, v in equations.items() if k not in ('sig_T', 'sig_LT', 'sig_TT')}
        
    def create_sig_L_function():
        # Combine all equations into a single expression
        combined_eq = '; '.join([f"{k} = {v}" for k, v in prepared_equations.items()])
        combined_eq += f"; sig_L"  # Add final computation

        def sig_L_vectorized(tt, p1, p2, p3, p4):
            return ne.evaluate(combined_eq, {
                'tt': tt, 'qq': q2_set, 'ww': w_set,
                'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4
            })

        return sig_L_vectorized
        
    # Create the optimized function
    fun_Sig_L_optimized = create_sig_L_function()
        

###############################################################################################################################################

# Function for SigL
def fun_Sig_L(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    
    # Use a list comprehension to get parameters, defaulting to 0.0 if not available
    p1, p2, p3, p4 = [par[i] if i < len(par) else 0.0 for i in range(4)]

    f = fun_Sig_L_optimized(tt, p1, p2, p3, p4)

    return f

###############################################################################################################################################

# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    
    # Use a list comprehension to get parameters, defaulting to 0.0 if not available
    p5, p6, p7, p8 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    
    # Create a dictionary for local variables
    local_vars = {'q2_set':q2_set, 'w_set':w_set, 'tt': tt, 'qq': qq, 'ww': ww, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4}
    
    # Add math functions to local_vars
    #math_functions = {name: getattr(math, name) for name in dir(math) if callable(getattr(math, name))}
    #local_vars.update(math_functions)
    
    # Evaluate equations
    for key, equation in equations.items():
        if key not in ('sig_L', 'sig_LT', 'sig_TT'):
            try:
                if DEBUG:
                    logging.debug(f"Evaluating equation for {key}: {equation}")
                #local_vars[key] = eval(equation, {"__builtins__": None}, local_vars)
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
    
    return local_vars['sig_T']

###############################################################################################################################################

# Function for SigLT
def fun_Sig_LT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    
    # Use a list comprehension to get parameters, defaulting to 0.0 if not available
    p1, p2, p3, p4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    
    # Create a dictionary for local variables
    local_vars = {'q2_set':q2_set, 'w_set':w_set, 'tt': tt, 'qq': qq, 'ww': ww, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4}
    
    # Add math functions to local_vars
    #math_functions = {name: getattr(math, name) for name in dir(math) if callable(getattr(math, name))}
    #local_vars.update(math_functions)
    
    # Evaluate equations
    for key, equation in equations.items():
        if key not in ('sig_L', 'sig_T', 'sig_TT'):
            try:
                if DEBUG:
                    logging.debug(f"Evaluating equation for {key}: {equation}")
                #local_vars[key] = eval(equation, {"__builtins__": None}, local_vars)
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
    
    return local_vars['sig_LT']

###############################################################################################################################################

# Function for SigTT
def fun_Sig_TT(x, par):
    tt = abs(x[0])
    q2_set = float(Q2.replace("p","."))
    w_set = float(W.replace("p","."))
    qq = q2_set
    ww = w_set
    
    # Use a list comprehension to get parameters, defaulting to 0.0 if not available
    p1, p2, p3, p4 = [par[i] if i < len(par) else 0.0 for i in range(4)]
    
    # Create a dictionary for local variables
    local_vars = {'q2_set':q2_set, 'w_set':w_set, 'tt': tt, 'qq': qq, 'ww': ww, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4}
    
    # Add math functions to local_vars
    #math_functions = {name: getattr(math, name) for name in dir(math) if callable(getattr(math, name))}
    #local_vars.update(math_functions)
    
    # Evaluate equations
    for key, equation in equations.items():
        if key not in ('sig_L', 'sig_T', 'sig_LT'):
            try:
                if DEBUG:
                    logging.debug(f"Evaluating equation for {key}: {equation}")
                #local_vars[key] = eval(equation, {"__builtins__": None}, local_vars)
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
    
    return local_vars['sig_TT']

###############################################################################################################################################
