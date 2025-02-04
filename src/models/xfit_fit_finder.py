#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-02-03 20:19:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import numpy as np
from ROOT import TGraph, TGraphErrors, TF1, TCanvas, TText, TLatex, TLegend, kRed, kBlue, kGreen, kMagenta, kBlack
import sys, math, time, random

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import (adaptive_regularization, calculate_cost, adaptive_cooling,
                     simulated_annealing, acceptance_probability, adjust_params, 
                     local_search, select_valid_parameter, get_central_value, 
                     calculate_information_criteria)

##################################################################################################################################################

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

###############################################################################################################################################
# Import separated xsects models

from xfit_active import fun_Sig_L_wrapper, fun_Sig_T_wrapper, fun_Sig_LT_wrapper, fun_Sig_TT_wrapper

##################################################################################################################################################

###############################################################################
# Helper function to do one run of simulated annealing + local search
###############################################################################
def global_fit_attempt(fun_wrapper, g_sig_data, num_params,
                       x_min, x_max,
                       max_iterations=200, initial_temperature=1.0,
                       param_bound=1e4, local_search_interval=25):
    """
    Performs one run of simulated annealing + optional local search
    from a random initial guess in [-param_bound, param_bound].
    Returns best_params, best_cost, plus optional debug info.
    """

    # Create the TF1 for the user-defined function
    f_sig = TF1("model_func", fun_wrapper, x_min, x_max, num_params)

    # # For demonstration: Name the parameters p0, p1, etc.
    for i in range(num_params):
        f_sig.SetParName(i, f"p{i}")

    # We want to feed the function the best_params each iteration:
    # But inside cost calc, we do: f_sig.SetParameter(i, param_val).
    # We'll do that before each cost evaluation.

    num_events = g_sig_data.GetN()

    # Random initial guess
    current_params = [random.uniform(-param_bound, param_bound)
                      for _ in range(num_params)]
    current_params = sanitize_params(current_params)

    # Evaluate cost using our utility function
    for i in range(num_params):
        f_sig.SetParameter(i, current_params[i])
    best_cost, lambda_reg = calculate_cost(f_sig, g_sig_data,
                                           current_params, 
                                           num_events, num_params)
    best_params = current_params[:]

    temperature = initial_temperature
    cost_history = [best_cost]

    for iteration in range(max_iterations):
        # --- Simulated Annealing Step ---
        # Perturb each param
        trial_params = [simulated_annealing(p, temperature)
                        for p in best_params]

        # Evaluate cost
        for i in range(num_params):
            f_sig.SetParameter(i, trial_params[i])
        new_cost, new_lambda = calculate_cost(f_sig, g_sig_data,
                                              trial_params,
                                              num_events, num_params,
                                              lambda_reg)

        cost_history.append(new_cost)
        # Update adaptive regularization
        lambda_reg = adaptive_regularization(cost_history, new_lambda)

        # Accept or reject
        acc_prob = acceptance_probability(best_cost, new_cost, temperature)
        if acc_prob > random.random():
            best_params = trial_params[:]
            best_cost = new_cost

        # --- Optional Local Search every local_search_interval steps ---
        if (iteration+1) % local_search_interval == 0:
            # local search starts from best_params
            for i in range(num_params):
                f_sig.SetParameter(i, best_params[i])
            refined_params = local_search(best_params, f_sig, num_params)
            
            # Evaluate cost after local search
            for i in range(num_params):
                f_sig.SetParameter(i, refined_params[i])
            refined_cost, _ = calculate_cost(f_sig, g_sig_data,
                                             refined_params,
                                             num_events, num_params,
                                             lambda_reg)
            if refined_cost < best_cost:
                best_params = refined_params[:]
                best_cost = refined_cost

        # Update temperature
        temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

    return best_params, best_cost

###############################################################################
# Main parameterize function that does the multi-start approach
###############################################################################
def parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec,
                 prv_par_vec, prv_err_vec, prv_chi2_vec,
                 fixed_params, outputpdf, full_optimization=True):

    # Extract from inpDict:
    num_optimizations = inpDict.get("num_optimizations", 10)
    max_iterations    = inpDict.get("max_iterations", 200)
    fit_params        = inpDict.get("fit_params", {})
    x_min             = inpDict.get("tmin_range", 0.0)
    x_max             = inpDict.get("tmax_range", 1.0)

    # Example: the graph with data
    # In reality, you’d build from your TTree or hist/histos:
    g_sig_data = TGraphErrors()
    # Fill in your data points here...
    # e.g. g_sig_data.SetPoint(0, x_val, y_val)
    #      g_sig_data.SetPointError(0, 0, y_err)

    # For demonstration, let's assume we have 5 points:
    for i in range(5):
        xval = 0.2 * i
        yval = math.sin(xval) + 0.1*random.random()
        g_sig_data.SetPoint(i, xval, yval)
        g_sig_data.SetPointError(i, 0.0, 0.02)

    # Let’s suppose we handle signals "L" and "T"
    for sig_name, val in fit_params.items():
        if sig_name in fixed_params:
            # Skip if user wants it fixed
            continue

        # Suppose val holds: (num_params, initial_guess, eqn_str),
        # or just the number of parameters, etc.
        # Example: we define 2 parameters for a test.
        num_params, initial_guess, eqn_str = val

        print(f"\n=== Global Minimization for {sig_name} with up to {num_params} params ===")

        # Build the function wrapper
        if sig_name == "L":
            # Use the wrapper from above
            fun_wrapper = fun_Sig_L_wrapper(g=1.0, Q2=2.5, W=2.0, th=45.0)
        elif sig_name == "T":
            fun_wrapper = fun_Sig_T_wrapper(g=1.0, Q2=2.5, W=2.0, th=45.0)
        else:
            raise ValueError("Unknown signal name. Provide the correct function wrapper.")

        best_global_params = None
        best_global_cost   = float("inf")

        # Multi-Start / Multi-run approach
        for run_idx in range(num_optimizations):
            params_candidate, cost_candidate = global_fit_attempt(
                fun_wrapper, g_sig_data,
                num_params=num_params,
                x_min=x_min, x_max=x_max,
                max_iterations=max_iterations,
                initial_temperature=2.0,  # maybe higher for exploration
                param_bound=1e4,         # domain for random initialization
                local_search_interval=25
            )
            if cost_candidate < best_global_cost:
                best_global_cost   = cost_candidate
                best_global_params = params_candidate[:]

        print(f"Best overall cost for {sig_name} = {best_global_cost:.5f}")
        print("Best params:", best_global_params)
        print("---------------------------------------------------")

    # Return or store your best-fit parameters as desired
    return
