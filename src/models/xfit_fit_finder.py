#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-22 03:22:18 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import random
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TGraph, TGraphErrors, TCanvas
from ROOT import TF1, TFitResultPtr
import numpy as np
import math
import time
import gc
import os, sys

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import adaptive_cooling, simulated_annealing, acceptance_probability, adjust_params, local_search, select_valid_parameter, within_tolerance, get_central_value

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

def find_fit(inpDict, par_vec, par_err_vec, par_chi2_vec):

    # Create lists to store graph objects outside the loop
    graphs_sig_fit = []
    graphs_sig_p0 = []
    graphs_sig_p1 = []
    graphs_sig_p2 = []
    graphs_sig_p3 = []
    graphs_sig_converge = []
    graphs_sig_temp = []
    graphs_sig_accept = []

    c2 = TCanvas("c2", "c2", 800, 800)
    c2.Divide(2, 2)

    # Create ROOT canvases for additional parameter convergence plots
    c3 = TCanvas("c3", "Parameter Convergence", 800, 800)
    c3.Divide(2, 2)
    c4 = TCanvas("c4", "Fit Convergence", 800, 800)
    c4.Divide(2, 2)
    c5 = TCanvas("c5", "Temperature", 800, 800)
    c5.Divide(2, 2)
    c6 = TCanvas("c6", "Acceptance Probability", 800, 800)
    c6.Divide(2, 2)

    q2_set = inpDict["q2_set"]
    w_set = inpDict["w_set"]
    nsep, g_vec, w_vec, q2_vec, th_vec = inpDict["objects"]
    max_iterations = inpDict["max_iterations"]
    num_optimizations = inpDict["num_optimizations"]
    sine_exp_LT = inpDict["sine_exp_LT"]
    sine_exp_TT = inpDict["sine_exp_TT"]
    tmin_range = inpDict["tmin_range"]
    tmax_range = inpDict["tmax_range"]
    Q2min_range = inpDict["Q2min_range"]
    Q2max_range = inpDict["Q2max_range"]
    iter_num = inpDict["iter_num"]            
    outputpdf = inpDict["outputpdf"]
    
    err_sets = inpDict["err_sets"]
    chi2_sets = inpDict["chi2_sets"]
    fit_params = inpDict["fit_params"]

    # Using central bin value to determine best fit
    q2_center_val = get_central_value(q2_vec)
    w_center_val = get_central_value(w_vec)
    th_center_val = get_central_value(th_vec)
    fun_Sig_L = fun_Sig_L_wrapper(q2_center_val, w_center_val)
    fun_Sig_T = fun_Sig_T_wrapper(q2_center_val, w_center_val)
    fun_Sig_LT = fun_Sig_LT_wrapper(q2_center_val, w_center_val) # Sine terms defined in xfit_in_t.py
    fun_Sig_TT = fun_Sig_TT_wrapper(q2_center_val, w_center_val) # Sine terms defined in xfit_in_t.py

    '''
    fun_Sig_L = fun_Sig_L_wrapper(q2_set, w_set)
    fun_Sig_T = fun_Sig_T_wrapper(q2_set, w_set)
    fun_Sig_LT = fun_Sig_LT_wrapper(q2_set, w_set) # Sine terms defined in xfit_in_t.py
    fun_Sig_TT = fun_Sig_TT_wrapper(q2_set, w_set) # Sine terms defined in xfit_in_t.py
    '''
    
    '''
    # Build the final dictionary excluding good fits from previous iteration (within tolerance of 1e-3)
    sig_fit_dict = {
        key: {"params": fit_params[key]}
        for key, values in chi2_sets.items()
        if not within_tolerance(values)
    }

    # Safely handle the case where sig_fit_dict is empty
    if not sig_fit_dict:
        print("ERROR: All fits are within 1e-3 of untity!")
        sys.exit(2)

    # Create a boolean dictionary to indicate which keys are still present in sig_fit_dict
    presence_dict = {key: key in sig_fit_dict for key in fit_params}
    '''
    
    num_events = nsep.GetEntries()    
    
    for it, (key, val) in enumerate(fit_params.items()):

        sig_name = key
        # Grab parameters used by functional forms
        num_params, initial_params, equation_str = inpDict["initial_params"](sig_name, val)

        # Checks initial parameters and replaces zeros to avoid errors
        #initial_params = [v if abs(v) > 0.0 else max_iterations for v in initial_params]
        initial_params = [v if abs(v) > 0.0 else 1.0 for v in initial_params]

        # String list of initial parameters
        param_str = ', '.join(str(param) for param in initial_params)

        if num_events <= num_params:
            print(f"\n\nWARNING: The number of parameters ({num_params}) for Sig {sig_name} is greater than or equal to the number of data points ({num_events})! Using adaptive regularization methods for determining quality of fit...")
            fit_convergence_type = "Adapt. Reg." # Adaptive Regularization            
        else:
            fit_convergence_type = "Red. Chi-Square"
            
        if num_params == 1:

            # 1 param
            #######
            # Sig #
            #######

            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Paramters: ({param_str})")
            print(f"{equation_str}")
            print("/*--------------------------------------------------*/")

            best_overall_params = None
            best_overall_cost = float('inf')
            total_iteration = 0
            max_param_value = 1e4

            # Regularization strength (used when num_events > num_params)
            # Initialize adaptive regularization parameters
            lambda_reg = 0.01  # Initial regularization strength
            lambda_increase = 1.05  # Factor to increase lambda (+ 5%)
            lambda_decrease = 0.95  # Factor to decrease lambda (- 5%)
            lambda_min = 1e-6  # Minimum lambda value
            lambda_max = 1.0  # Maximum lambda value
            cost_history = []
            tolerance = 1e-6  # Early stopping tolerance
            alpha = 0.1  # Small constant to penalize model complexity            
            
            # Store the parameter values and chi-square values for each iteration
            params_sig_history = {'p0': []}

            # Create TGraphs for parameter convergence
            graph_sig_p0 = TGraph()
            graph_sig_chi2 = TGraph()
            graph_sig_temp = TGraph()
            graph_sig_accept = TGraph()
            graphs_sig_p0.append(graph_sig_p0)
            graphs_sig_p1.append(0.0)
            graphs_sig_p2.append(0.0)
            graphs_sig_p3.append(0.0)
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)

            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Record the start time
            start_time = time.time()

            for start in range(num_optimizations):
                print("\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                iteration = 0

                initial_temperature = 1.0
                temperature = initial_temperature
                unchanged_iterations = 0
                max_unchanged_iterations = 5

                # Initialize adaptive parameter limits
                par_sig_0 = initial_params[0]
                par_sig_err_0 = 0.0

                # Track the best solution
                best_params = par_sig_0
                best_cost = float('inf')
                previous_params = best_params
                best_errors = par_sig_err_0

                # Check for local minima
                local_minima = []
                local_iterations = 0
                tabu_list = set()

                # Local search
                local_search_interval = 25

                while iteration <= max_iterations:

                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
                    sys.stdout.flush()

                    try:
                        # Perturb parameters
                        current_params = simulated_annealing(par_sig_0, temperature)

                        # Insert tabu list check here
                        if current_params not in tabu_list:
                            tabu_list.add(current_params)
                        else:
                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        g_sig = TGraphErrors()
                        for i in range(nsep.GetSelectedRows()):
                            g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                            g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                        for i in range(len(w_vec)):
                            if sig_name == "LT":
                                sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                            elif sig_name == "TT":
                                sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                            else:
                                sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                            graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                        if sig_name == "L":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
                        elif sig_name == "T":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
                        elif sig_name == "LT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
                        elif sig_name == "TT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
                        f_sig.SetParNames("p0")
                        f_sig.SetParameter(0, current_params)
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)

                        r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                        #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                        f_sig_status = f_sig.GetNDF() != 0

                        params_sig_history['p0'].append(current_params)

                        if num_events > num_params:
                            # Calculate the cost (reduced chi-square value) for the current parameters                            
                            current_cost = f_sig.GetChisquare()/(num_events-num_params) # Divided by DoF for red. chi-squared
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, current_cost, temperature)
                        else:
                            # Define lambda values to try (in log space for better exploration)
                            lambda_values = np.logspace(np.log10(lambda_min), np.log10(lambda_max), 10)
                            best_cost_iteration = float('inf')  # Initialize to a large value
                            best_lambda = lambda_reg  # Start with initial lambda                            
                            # Adaptive search for the best lambda
                            for lambda_try in lambda_values:
                                residuals = []  # Store residuals for this lambda value
                                # Loop through each data point (event)
                                for i in range(num_events):
                                    observed = g_sig.GetY()[i]  # Observed value
                                    expected = f_sig.Eval(g_sig.GetX()[i])  # Model prediction
                                    # Compute residuals (normalized by error if available)
                                    if g_sig.GetEY()[i] != 0:
                                        residual = (observed - expected) / g_sig.GetEY()[i]
                                    else:
                                        residual = observed - expected
                                    residuals.append(residual)
                                # Calculate the mean squared error (MSE)
                                mse = np.mean(np.square(residuals))
                                # Compute L2 regularization term (sum of parameter squares)
                                l2_reg = sum(p**2 for p in current_params)
                                # Generalized reduced chi-squared formula
                                current_cost_try = (mse + lambda_try * l2_reg) / (num_events + alpha * num_params)
                                # Check if this lambda gives the best cost so far
                                if current_cost_try < best_cost_iteration:
                                    best_cost_iteration = current_cost_try
                                    best_lambda = lambda_try
                                    best_params = current_params.copy()  # Save parameters for this lambda
                                # Optional: Early stopping if the improvement is small
                                if abs(current_cost_try - best_cost_iteration) < tolerance:
                                    break
                            # Use the best lambda and parameters found
                            lambda_reg = best_lambda
                            current_params = best_params  # Update model parameters
                            current_cost = best_cost_iteration  # Final cost for this iteration
                            # Store cost history for analysis
                            cost_history.append(current_cost)
                            # Calculate acceptance probability for simulated annealing (if applicable)
                            accept_prob = acceptance_probability(current_cost, best_cost_iteration, temperature)
                            
                        current_params = f_sig.GetParameter(0)
                        current_errors = f_sig.GetParError(0)

                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, current_params)
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                        # If the new cost is better or accepted by the acceptance probability, update the best parameters
                        if accept_prob > random.random():
                            best_params = current_params
                            best_cost = current_cost
                            best_errors = current_errors

                        if iteration % local_search_interval == 0:
                            current_params = local_search(current_params, f_sig, num_params)
                            par_sig_0 = current_params

                        # Check if current parameters haven't changed for the past N iterations
                        if len(params_sig_history['p0']) >= max_unchanged_iterations:
                            if np.allclose(round(params_sig_history['p0'][-2], 3), round(params_sig_history['p0'][-1], 3), atol=5.0):
                                unchanged_iterations += 1
                            else:
                                unchanged_iterations = 0

                        # Adjust the cooling rate if parameters haven't changed for N iterations
                        if unchanged_iterations >= max_unchanged_iterations:
                            if not any(np.allclose([current_params], minima, atol=5.0) for minima in local_minima):                    
                                local_minima.append([
                                    current_params
                                ])
                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        previous_params = current_params                

                        # Update parameters with the best found so far
                        par_sig_0 = best_params
                        par_sig_err_0 = best_errors

                        # Update the temperature
                        temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                        # Check if current_params are close to any local minimum
                        if any(np.allclose([current_params], minima, atol=5.0) for minima in local_minima):
                            #print("WARNING: Parameters p0={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params))

                            current_params = adjust_params(best_params)
                            par_sig_0 = current_params
                            par_sig_err_0 = 0.0

                    except (TypeError or ZeroDivisionError) as e:
                        #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                        # Adjust parameter limits within a random number
                        par_sig_0 = initial_params
                        par_sig_err_0 = 0.0

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0                

                # After the while loop, check if this run found a better solution
                if abs(best_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = best_cost
                    best_overall_params = best_params
                    best_overall_errors = best_errors
                
            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            best_overall_params = [best_overall_params]
            best_overall_errors = [best_overall_errors]

            try:            
                while len(best_overall_params) < 4:
                    best_overall_params.append(0.0)
                    best_overall_errors.append(0.0)
            except TypeError:
                print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                sys.exit(2)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[1])
            par_err_vec.append(best_overall_errors[2])
            par_err_vec.append(best_overall_errors[3])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_fit = TGraphErrors()
            
            graphs_sig_fit.append(g_sig_fit)
            
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                if sig_name == "LT":
                    sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                elif sig_name == "TT":
                    sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                else:
                    sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            graphs_sig_fit[it].GetYaxis().SetTitleOffset(1.5)
            graphs_sig_fit[it].GetYaxis().SetTitleSize(0.035)
            graphs_sig_fit[it].GetYaxis().CenterTitle()
            
            # Set axis limits to ensure everything is shown
            x_min = min(graphs_sig_fit[it].GetX())
            x_max = max(graphs_sig_fit[it].GetX())
            y_min = min(graphs_sig_fit[it].GetY())
            y_max = max(graphs_sig_fit[it].GetY())

            # You can also set a margin to ensure all points are visible
            margin = 0.1
            graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

            if sig_name == "L":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
            f_sig.SetParNames("p0")
            f_sig.FixParameter(0, best_overall_params[0])
            
            # Evaluate the fit function at several points to determine its range
            n_points = 100  # Number of points to evaluate the fit function
            fit_y_values = [f_sig.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)

            # Extend the y-axis range to include the fit function range
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)

            # Set a margin to ensure all points are visible
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

            r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
            f_sig.Draw("same")
            
            #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
            f_sig_status = f_sig.GetNDF() != 0
            f_sig_status_message = "Fit Successful" if f_sig_status else "Fit Failed"

            fit_status = TText()
            fit_status.SetTextSize(0.04)
            fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sig_status_message))
            
            # Calculate the minimum and maximum values from the graphs
            min_sig_y = float('inf')
            max_sig_y = float('-inf')

            # Update min_sig_y and max_sig_y based on each graph's values
            for graph in [graphs_sig_p0[it]]:
                n_points = graph.GetN()
                for i in range(n_points):
                    y = graph.GetY()[i]
                    if y < min_sig_y:
                        min_sig_y = y
                    if y > max_sig_y:
                        max_sig_y = y

            # Scale the y-axis
            graphs_sig_p0[it].SetMinimum(min_sig_y * 0.9)
            graphs_sig_p0[it].SetMaximum(max_sig_y * 1.1)
            c2.Update()

            # Plot parameter convergence
            c3.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_p0[it].SetTitle(f"Sig {sig_name} Parameter Convergence;Optimization Run;Parameter")
            graphs_sig_p0[it].SetLineColor(ROOT.kRed)
            graphs_sig_p0[it].Draw("ALP")
            c3.Update()
            
            # Plot chi-square convergence
            c4.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
            graphs_sig_converge[it].SetLineColor(ROOT.kBlack)
            graphs_sig_converge[it].Draw("ALP")
            converge_status = TText()
            converge_status.SetTextSize(0.04)
            converge_status.DrawTextNDC(0.35, 0.85, f"Best cost: {best_overall_cost:.3f}")
            c4.Update()
            
            # Plot temperature
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability
            c6.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(ROOT.kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()
            
            print("\n")    

        elif num_params == 2:

            # 2 params
            #######
            # Sig #
            #######

            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Paramters: ({param_str})")
            print(f"{equation_str}")            
            print("/*--------------------------------------------------*/")

            best_overall_params = None
            best_overall_cost = float('inf')
            total_iteration = 0
            max_param_value = 1e4

            # Regularization strength (used when num_events > num_params)
            # Initialize adaptive regularization parameters
            lambda_reg = 0.01  # Initial regularization strength
            lambda_increase = 1.05  # Factor to increase lambda (+ 5%)
            lambda_decrease = 0.95  # Factor to decrease lambda (- 5%)
            lambda_min = 1e-6  # Minimum lambda value
            lambda_max = 1.0  # Maximum lambda value
            cost_history = []
            tolerance = 1e-6  # Early stopping tolerance
            alpha = 0.1  # Small constant to penalize model complexity            
            
            # Store the parameter values and chi-square values for each iteration
            params_sig_history = {'p0': [], 'p1': []}

            # Create TGraphs for parameter convergence
            graph_sig_p0 = TGraph()
            graph_sig_p1 = TGraph()
            graph_sig_chi2 = TGraph()
            graph_sig_temp = TGraph()
            graph_sig_accept = TGraph()
            graphs_sig_p0.append(graph_sig_p0)
            graphs_sig_p1.append(graph_sig_p1)
            graphs_sig_p2.append(0.0)
            graphs_sig_p3.append(0.0)            
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)

            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Record the start time
            start_time = time.time()

            for start in range(num_optimizations):
                print("\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                iteration = 0

                initial_temperature = 1.0
                temperature = initial_temperature
                unchanged_iterations = 0
                max_unchanged_iterations = 5

                # Initialize adaptive parameter limits
                par_sig_0 = initial_params[0]
                par_sig_1 = initial_params[1]
                par_sig_err_0 = 0.0
                par_sig_err_1 = 0.0

                # Track the best solution
                best_params = [par_sig_0, par_sig_1]
                best_cost = float('inf')
                best_errors = [par_sig_err_0, par_sig_err_1]
                previous_params = best_params[:]

                # Check for local minima
                local_minima = []
                tabu_list = set()

                # Local search
                local_search_interval = 25

                while iteration <= max_iterations:

                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
                    sys.stdout.flush()

                    try:
                        # Perturb parameters
                        current_params = [
                            simulated_annealing(par_sig_0, temperature),
                            simulated_annealing(par_sig_1, temperature)
                        ]

                        # Insert tabu list check here
                        if tuple(current_params) not in tabu_list:
                            tabu_list.add(tuple(current_params))
                        else:
                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        g_sig = TGraphErrors()
                        for i in range(nsep.GetSelectedRows()):
                            g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                            g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                        for i in range(len(w_vec)):
                            if sig_name == "LT":
                                sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                            elif sig_name == "TT":
                                sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                            else:
                                sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                            graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                        if sig_name == "L":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
                        elif sig_name == "T":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
                        elif sig_name == "LT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
                        elif sig_name == "TT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
                        f_sig.SetParNames("p0", "p1")
                        f_sig.SetParameter(0, current_params[0])
                        f_sig.SetParameter(1, current_params[1])
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)

                        r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                        #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                        f_sig_status = f_sig.GetNDF() != 0

                        params_sig_history['p0'].append(current_params[0])
                        params_sig_history['p1'].append(current_params[1])

                        if num_events > num_params:
                            # Calculate the cost (reduced chi-square value) for the current parameters                            
                            current_cost = f_sig.GetChisquare()/(num_events-num_params) # Divided by DoF for red. chi-squared
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, current_cost, temperature)
                        else:
                            # Define lambda values to try (in log space for better exploration)
                            lambda_values = np.logspace(np.log10(lambda_min), np.log10(lambda_max), 10)
                            best_cost_iteration = float('inf')  # Initialize to a large value
                            best_lambda = lambda_reg  # Start with initial lambda                            
                            # Adaptive search for the best lambda
                            for lambda_try in lambda_values:
                                residuals = []  # Store residuals for this lambda value
                                # Loop through each data point (event)
                                for i in range(num_events):
                                    observed = g_sig.GetY()[i]  # Observed value
                                    expected = f_sig.Eval(g_sig.GetX()[i])  # Model prediction
                                    # Compute residuals (normalized by error if available)
                                    if g_sig.GetEY()[i] != 0:
                                        residual = (observed - expected) / g_sig.GetEY()[i]
                                    else:
                                        residual = observed - expected
                                    residuals.append(residual)
                                # Calculate the mean squared error (MSE)
                                mse = np.mean(np.square(residuals))
                                # Compute L2 regularization term (sum of parameter squares)
                                l2_reg = sum(p**2 for p in current_params)
                                # Generalized reduced chi-squared formula
                                current_cost_try = (mse + lambda_try * l2_reg) / (num_events + alpha * num_params)
                                # Check if this lambda gives the best cost so far
                                if current_cost_try < best_cost_iteration:
                                    best_cost_iteration = current_cost_try
                                    best_lambda = lambda_try
                                    best_params = current_params.copy()  # Save parameters for this lambda
                                # Optional: Early stopping if the improvement is small
                                if abs(current_cost_try - best_cost_iteration) < tolerance:
                                    break
                            # Use the best lambda and parameters found
                            lambda_reg = best_lambda
                            current_params = best_params  # Update model parameters
                            current_cost = best_cost_iteration  # Final cost for this iteration
                            # Store cost history for analysis
                            cost_history.append(current_cost)
                            # Calculate acceptance probability for simulated annealing (if applicable)
                            accept_prob = acceptance_probability(current_cost, best_cost_iteration, temperature)
                            
                        current_params = [
                            f_sig.GetParameter(0),
                            f_sig.GetParameter(1)
                        ]

                        current_errors = [
                            f_sig.GetParError(0),
                            f_sig.GetParError(1)
                        ]

                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, current_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, current_params[1])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                        # If the new cost is better or accepted by the acceptance probability, update the best parameters
                        if accept_prob > random.random():
                            best_params = current_params
                            best_cost = current_cost
                            best_errors = current_errors

                        if iteration % local_search_interval == 0:
                            current_params = local_search(current_params, f_sig, num_params)
                            par_sig_0, par_sig_1 = current_params

                        # Check if current parameters haven't changed for the past N iterations
                        if len(params_sig_history['p0']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p1']) >= max_unchanged_iterations:
                            if np.allclose(round(params_sig_history['p0'][-2], 3), round(params_sig_history['p0'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p1'][-2], 3), round(params_sig_history['p1'][-1], 3), atol=5.0):
                                unchanged_iterations += 1        
                            else:
                                unchanged_iterations = 0

                        # Adjust the cooling rate if parameters haven't changed for N iterations
                        if unchanged_iterations >= max_unchanged_iterations:
                            if not any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):                    
                                local_minima.append([
                                    current_params[0],
                                    current_params[1]
                                ])

                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        previous_params = current_params[:]

                        # Update parameters with the best found so far
                        par_sig_0, par_sig_1 = best_params
                        par_sig_err_0, par_sig_err_1 = best_errors

                        # Update the temperature
                        temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                        # Check if current_params are close to any local minimum
                        if any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):
                            #print("WARNING: Parameters p0={:.3e}, p1={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1]))

                            current_params = adjust_params(best_params)
                            par_sig_0, par_sig_1 = current_params
                            par_sig_err_0, par_sig_err_1 = [0.0 for _ in range(num_params)]

                    except (TypeError or ZeroDivisionError) as e:
                        #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                        # Adjust parameter limits within a random number
                        par_sig_0, par_sig_1 = initial_params
                        par_sig_err_0, par_sig_err_1 = [0.0 for _ in range(num_params)]

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                # After the while loop, check if this run found a better solution
                if abs(best_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = best_cost
                    best_overall_params = best_params[:]
                    best_overall_errors = best_errors[:]
                
            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            try:            
                while len(best_overall_params) < 4:
                    best_overall_params.append(0.0)
                    best_overall_errors.append(0.0)
            except TypeError:
                print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                sys.exit(2)
                
            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[1])
            par_err_vec.append(best_overall_errors[2])
            par_err_vec.append(best_overall_errors[3])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_fit = TGraphErrors()

            graphs_sig_fit.append(g_sig_fit)
            
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                if sig_name == "LT":
                    sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                elif sig_name == "TT":
                    sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                else:
                    sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            graphs_sig_fit[it].GetYaxis().SetTitleOffset(1.5)
            graphs_sig_fit[it].GetYaxis().SetTitleSize(0.035)
            graphs_sig_fit[it].GetYaxis().CenterTitle()
            
            # Set axis limits to ensure everything is shown
            x_min = min(graphs_sig_fit[it].GetX())
            x_max = max(graphs_sig_fit[it].GetX())
            y_min = min(graphs_sig_fit[it].GetY())
            y_max = max(graphs_sig_fit[it].GetY())

            # You can also set a margin to ensure all points are visible
            margin = 0.1
            graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

            if sig_name == "L":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)    
            f_sig.SetParNames("p0", "p1")
            f_sig.FixParameter(0, best_overall_params[0])
            f_sig.FixParameter(1, best_overall_params[1])
        
            # Evaluate the fit function at several points to determine its range
            n_points = 100  # Number of points to evaluate the fit function
            fit_y_values = [f_sig.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)

            # Extend the y-axis range to include the fit function range
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)

            # Set a margin to ensure all points are visible
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

            r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
            f_sig.Draw("same")
            
            #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
            f_sig_status = f_sig.GetNDF() != 0
            f_sig_status_message = "Fit Successful" if f_sig_status else "Fit Failed"

            fit_status = TText()
            fit_status.SetTextSize(0.04)
            fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sig_status_message))
            
            # Calculate the minimum and maximum values from the graphs
            min_sig_y = float('inf')
            max_sig_y = float('-inf')

            # Update min_sig_y and max_sig_y based on each graph's values
            for graph in [graphs_sig_p0[it], graphs_sig_p1[it]]:
                n_points = graph.GetN()
                for i in range(n_points):
                    y = graph.GetY()[i]
                    if y < min_sig_y:
                        min_sig_y = y
                    if y > max_sig_y:
                        max_sig_y = y

            # Scale the y-axis
            graphs_sig_p0[it].SetMinimum(min_sig_y * 0.9)
            graphs_sig_p0[it].SetMaximum(max_sig_y * 1.1)
            c2.Update()

            # Plot parameter convergence
            c3.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_p0[it].SetTitle(f"Sig {sig_name} Parameter Convergence;Optimization Run;Parameter")
            graphs_sig_p0[it].SetLineColor(ROOT.kRed)
            graphs_sig_p1[it].SetLineColor(ROOT.kBlue)
            graphs_sig_p0[it].Draw("ALP")
            graphs_sig_p1[it].Draw("LP SAME")
            c3.Update()
            
            # Plot chi-square convergence
            c4.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
            graphs_sig_converge[it].SetLineColor(ROOT.kBlack)
            graphs_sig_converge[it].Draw("ALP")
            converge_status = TText()
            converge_status.SetTextSize(0.04)
            converge_status.DrawTextNDC(0.35, 0.85, f"Best cost: {best_overall_cost:.3f}")
            c4.Update()
            
            # Plot temperature
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability
            c6.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(ROOT.kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()
            
            print("\n")    

        elif num_params == 3:

            # 3 params
            #######
            # Sig #
            #######

            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Paramters: ({param_str})")
            print(f"{equation_str}")            
            print("/*--------------------------------------------------*/")

            best_overall_params = None
            best_overall_cost = float('inf')
            total_iteration = 0
            max_param_value = 1e4

            # Regularization strength (used when num_events > num_params)
            # Initialize adaptive regularization parameters
            lambda_reg = 0.01  # Initial regularization strength
            lambda_increase = 1.05  # Factor to increase lambda (+ 5%)
            lambda_decrease = 0.95  # Factor to decrease lambda (- 5%)
            lambda_min = 1e-6  # Minimum lambda value
            lambda_max = 1.0  # Maximum lambda value
            cost_history = []
            tolerance = 1e-6  # Early stopping tolerance
            alpha = 0.1  # Small constant to penalize model complexity            
            
            # Store the parameter values and chi-square values for each iteration
            params_sig_history = {'p0': [], 'p1': [], 'p2': []}

            # Create TGraphs for parameter convergence
            graph_sig_p0 = TGraph()
            graph_sig_p1 = TGraph()
            graph_sig_p2 = TGraph()
            graph_sig_chi2 = TGraph()
            graph_sig_temp = TGraph()
            graph_sig_accept = TGraph()
            graphs_sig_p0.append(graph_sig_p0)
            graphs_sig_p1.append(graph_sig_p1)
            graphs_sig_p2.append(graph_sig_p2)
            graphs_sig_p3.append(0.0)
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)
            
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Record the start time
            start_time = time.time()

            for start in range(num_optimizations):
                print("\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                iteration = 0

                initial_temperature = 1.0
                temperature = initial_temperature
                unchanged_iterations = 0
                max_unchanged_iterations = 5

                # Initialize adaptive parameter limits
                par_sig_0 = initial_params[0]
                par_sig_1 = initial_params[1]
                par_sig_2 = initial_params[2]
                par_sig_err_0 = 0.0
                par_sig_err_1 = 0.0
                par_sig_err_2 = 0.0

                # Track the best solution
                best_params = [par_sig_0, par_sig_1, par_sig_2]
                best_cost = float('inf')
                best_errors = [par_sig_err_0, par_sig_err_1, par_sig_err_2]
                previous_params = best_params[:]

                # Check for local minima
                local_minima = []
                tabu_list = set()

                # Local search
                local_search_interval = 25

                while iteration <= max_iterations:

                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
                    sys.stdout.flush()

                    try:

                        # Perturb parameters
                        current_params = [
                            simulated_annealing(par_sig_0, temperature),
                            simulated_annealing(par_sig_1, temperature),
                            simulated_annealing(par_sig_2, temperature)
                        ]

                        # Insert tabu list check here
                        if tuple(current_params) not in tabu_list:
                            tabu_list.add(tuple(current_params))
                        else:
                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        g_sig = TGraphErrors()
                        for i in range(nsep.GetSelectedRows()):
                            g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                            g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                        for i in range(len(w_vec)):
                            if sig_name == "LT":
                                sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                            elif sig_name == "TT":
                                sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                            else:
                                sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                            graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                        if sig_name == "L":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
                        elif sig_name == "T":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
                        elif sig_name == "LT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
                        elif sig_name == "TT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
                        f_sig.SetParNames("p0", "p1", "p2")
                        f_sig.SetParameter(0, current_params[0])
                        f_sig.SetParameter(1, current_params[1])
                        f_sig.SetParameter(2, current_params[2])
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)
                        f_sig.SetParLimits(2, -max_param_value, max_param_value)

                        r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                        #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                        f_sig_status = f_sig.GetNDF() != 0

                        params_sig_history['p0'].append(current_params[0])
                        params_sig_history['p1'].append(current_params[1])
                        params_sig_history['p2'].append(current_params[2])

                        if num_events > num_params:
                            # Calculate the cost (reduced chi-square value) for the current parameters                            
                            current_cost = f_sig.GetChisquare()/(num_events-num_params) # Divided by DoF for red. chi-squared
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, current_cost, temperature)
                        else:
                            # Define lambda values to try (in log space for better exploration)
                            lambda_values = np.logspace(np.log10(lambda_min), np.log10(lambda_max), 10)
                            best_cost_iteration = float('inf')  # Initialize to a large value
                            best_lambda = lambda_reg  # Start with initial lambda                            
                            # Adaptive search for the best lambda
                            for lambda_try in lambda_values:
                                residuals = []  # Store residuals for this lambda value
                                # Loop through each data point (event)
                                for i in range(num_events):
                                    observed = g_sig.GetY()[i]  # Observed value
                                    expected = f_sig.Eval(g_sig.GetX()[i])  # Model prediction
                                    # Compute residuals (normalized by error if available)
                                    if g_sig.GetEY()[i] != 0:
                                        residual = (observed - expected) / g_sig.GetEY()[i]
                                    else:
                                        residual = observed - expected
                                    residuals.append(residual)
                                # Calculate the mean squared error (MSE)
                                mse = np.mean(np.square(residuals))
                                # Compute L2 regularization term (sum of parameter squares)
                                l2_reg = sum(p**2 for p in current_params)
                                # Generalized reduced chi-squared formula
                                current_cost_try = (mse + lambda_try * l2_reg) / (num_events + alpha * num_params)
                                # Check if this lambda gives the best cost so far
                                if current_cost_try < best_cost_iteration:
                                    best_cost_iteration = current_cost_try
                                    best_lambda = lambda_try
                                    best_params = current_params.copy()  # Save parameters for this lambda
                                # Optional: Early stopping if the improvement is small
                                if abs(current_cost_try - best_cost_iteration) < tolerance:
                                    break
                            # Use the best lambda and parameters found
                            lambda_reg = best_lambda
                            current_params = best_params  # Update model parameters
                            current_cost = best_cost_iteration  # Final cost for this iteration
                            # Store cost history for analysis
                            cost_history.append(current_cost)
                            # Calculate acceptance probability for simulated annealing (if applicable)
                            accept_prob = acceptance_probability(current_cost, best_cost_iteration, temperature)
                            
                        current_params = [
                            f_sig.GetParameter(0),
                            f_sig.GetParameter(1),
                            f_sig.GetParameter(2)
                        ]

                        current_errors = [
                            f_sig.GetParError(0),
                            f_sig.GetParError(1),
                            f_sig.GetParError(2)
                        ]

                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, current_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, current_params[1])
                        graphs_sig_p2[it].SetPoint(total_iteration, total_iteration, current_params[2])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                        # If the new cost is better or accepted by the acceptance probability, update the best parameters
                        if accept_prob > random.random():
                            best_params = current_params
                            best_cost = current_cost
                            best_errors = current_errors

                        if iteration % local_search_interval == 0:
                            current_params = local_search(current_params, f_sig, num_params)
                            par_sig_0, par_sig_1, par_sig_2 = current_params

                        # Check if current parameters haven't changed for the past N iterations
                        if len(params_sig_history['p0']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p1']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p2']) >= max_unchanged_iterations:
                            if np.allclose(round(params_sig_history['p0'][-2], 3), round(params_sig_history['p0'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p1'][-2], 3), round(params_sig_history['p1'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p2'][-2], 3), round(params_sig_history['p2'][-1], 3), atol=5.0):
                                unchanged_iterations += 1
                            else:
                                unchanged_iterations = 0

                        # Adjust the cooling rate if parameters haven't changed for N iterations
                        if unchanged_iterations >= max_unchanged_iterations:
                            if not any(np.allclose([current_params[0], current_params[1], current_params[2]], minima, atol=5.0) for minima in local_minima):
                                local_minima.append([
                                    current_params[0],
                                    current_params[1],
                                    current_params[2]
                                ])

                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        previous_params = current_params[:]

                        # Update parameters with the best found so far
                        par_sig_0, par_sig_1, par_sig_2 = best_params
                        par_sig_err_0, par_sig_err_1, par_sig_err_2 = best_errors

                        # Update the temperature
                        temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                        # Check if current_params are close to any local minimum
                        if any(np.allclose([current_params[0], current_params[1], current_params[2]], minima, atol=5.0) for minima in local_minima):
                            #print("WARNING: Parameters p0={:.3e}, p1={:.3e}, p2={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2]))

                            current_params = adjust_params(best_params)
                            par_sig_0, par_sig_1, par_sig_2 = current_params
                            par_sig_err_0, par_sig_err_1, par_sig_err_2 = [0.0 for _ in range(num_params)]

                    except (TypeError or ZeroDivisionError) as e:
                        #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                        # Adjust parameter limits within a random number
                        par_sig_0, par_sig_1, par_sig_2 = initial_params
                        par_sig_err_0, par_sig_err_1, par_sig_err_2 = [0.0 for _ in range(num_params)]

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                # After the while loop, check if this run found a better solution
                if abs(best_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = best_cost
                    best_overall_params = best_params[:]
                    best_overall_errors = best_errors[:]
                
            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            try:            
                while len(best_overall_params) < 4:
                    best_overall_params.append(0.0)
                    best_overall_errors.append(0.0)
            except TypeError:
                print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                sys.exit(2)
                
            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[1])
            par_err_vec.append(best_overall_errors[2])
            par_err_vec.append(best_overall_errors[3])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_fit = TGraphErrors()

            graphs_sig_fit.append(g_sig_fit)
                        
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                if sig_name == "LT":
                    sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                elif sig_name == "TT":
                    sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                else:
                    sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            graphs_sig_fit[it].GetYaxis().SetTitleOffset(1.5)
            graphs_sig_fit[it].GetYaxis().SetTitleSize(0.035)
            graphs_sig_fit[it].GetYaxis().CenterTitle()

            # Set axis limits to ensure everything is shown
            x_min = min(graphs_sig_fit[it].GetX())
            x_max = max(graphs_sig_fit[it].GetX())
            y_min = min(graphs_sig_fit[it].GetY())
            y_max = max(graphs_sig_fit[it].GetY())

            # You can also set a margin to ensure all points are visible
            margin = 0.1
            graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

            if sig_name == "L":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
            f_sig.SetParNames("p0", "p1", "p2")
            f_sig.FixParameter(0, best_overall_params[0])
            f_sig.FixParameter(1, best_overall_params[1])
            f_sig.FixParameter(2, best_overall_params[2])
        
            # Evaluate the fit function at several points to determine its range
            n_points = 100  # Number of points to evaluate the fit function
            fit_y_values = [f_sig.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)

            # Extend the y-axis range to include the fit function range
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)

            # Set a margin to ensure all points are visible
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

            r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
            f_sig.Draw("same")
            
            #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
            f_sig_status = f_sig.GetNDF() != 0
            f_sig_status_message = "Fit Successful" if f_sig_status else "Fit Failed"

            fit_status = TText()
            fit_status.SetTextSize(0.04)
            fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sig_status_message))
            
            # Calculate the minimum and maximum values from the graphs
            min_sig_y = float('inf')
            max_sig_y = float('-inf')

            # Update min_sig_y and max_sig_y based on each graph's values
            for graph in [graphs_sig_p0[it], graphs_sig_p1[it], graphs_sig_p2[it]]:
                n_points = graph.GetN()
                for i in range(n_points):
                    y = graph.GetY()[i]
                    if y < min_sig_y:
                        min_sig_y = y
                    if y > max_sig_y:
                        max_sig_y = y

            # Scale the y-axis
            graphs_sig_p0[it].SetMinimum(min_sig_y * 0.9)
            graphs_sig_p0[it].SetMaximum(max_sig_y * 1.1)
            c2.Update()            

            # Plot parameter convergence
            c3.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_p0[it].SetTitle(f"Sig {sig_name} Parameter Convergence;Optimization Run;Parameter")
            graphs_sig_p0[it].SetLineColor(ROOT.kRed)
            graphs_sig_p1[it].SetLineColor(ROOT.kBlue)
            graphs_sig_p2[it].SetLineColor(ROOT.kGreen)
            graphs_sig_p0[it].Draw("ALP")
            graphs_sig_p1[it].Draw("LP SAME")
            graphs_sig_p2[it].Draw("LP SAME")
            c3.Update()
            
            # Plot chi-square convergence
            c4.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
            graphs_sig_converge[it].SetLineColor(ROOT.kBlack)
            graphs_sig_converge[it].Draw("ALP")
            converge_status = TText()
            converge_status.SetTextSize(0.04)
            converge_status.DrawTextNDC(0.35, 0.85, f"Best cost: {best_overall_cost:.3f}")
            c4.Update()
            
            # Plot temperature
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability
            c6.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(ROOT.kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()
            
            print("\n")    

        elif num_params == 4:

            # 4 params
            #######
            # Sig #
            #######

            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Paramters: ({param_str})")
            print(f"{equation_str}")            
            print("/*--------------------------------------------------*/")    

            best_overall_params = None
            best_overall_cost = float('inf')
            total_iteration = 0
            max_param_value = 1e4

            # Regularization strength (used when num_events > num_params)
            # Initialize adaptive regularization parameters
            lambda_reg = 0.01  # Initial regularization strength
            lambda_increase = 1.05  # Factor to increase lambda (+ 5%)
            lambda_decrease = 0.95  # Factor to decrease lambda (- 5%)
            lambda_min = 1e-6  # Minimum lambda value
            lambda_max = 1.0  # Maximum lambda value
            cost_history = []
            tolerance = 1e-6  # Early stopping tolerance
            alpha = 0.1  # Small constant to penalize model complexity            
            
            # Store the parameter values and chi-square values for each iteration
            params_sig_history = {'p0': [], 'p1': [], 'p2': [], 'p3': []}

            # Create TGraphs for parameter convergence
            graph_sig_p0 = TGraph()
            graph_sig_p1 = TGraph()
            graph_sig_p2 = TGraph()
            graph_sig_p3 = TGraph()
            graph_sig_chi2 = TGraph()
            graph_sig_temp = TGraph()
            graph_sig_accept = TGraph()
            graphs_sig_p0.append(graph_sig_p0)
            graphs_sig_p1.append(graph_sig_p1)
            graphs_sig_p2.append(graph_sig_p2)
            graphs_sig_p3.append(graph_sig_p3)
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)
            
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Record the start time
            start_time = time.time()

            for start in range(num_optimizations):
                print("\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                iteration = 0

                initial_temperature = 1.0
                temperature = initial_temperature
                unchanged_iterations = 0
                max_unchanged_iterations = 5

                # Initialize adaptive parameter limits
                par_sig_0 = initial_params[0]
                par_sig_1 = initial_params[1]
                par_sig_2 = initial_params[2]
                par_sig_3 = initial_params[3]
                par_sig_err_0 = 0.0
                par_sig_err_1 = 0.0
                par_sig_err_2 = 0.0
                par_sig_err_3 = 0.0

                # Track the best solution
                best_params = [par_sig_0, par_sig_1, par_sig_2, par_sig_3]
                best_cost = float('inf')
                best_errors = [par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3]
                previous_params = best_params[:]

                # Check for local minima
                local_minima = []
                local_iterations = 0
                tabu_list = set()

                # Local search
                local_search_interval = 25

                while iteration <= max_iterations:

                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
                    sys.stdout.flush()

                    try:
                        # Perturb parameters
                        current_params = [
                            simulated_annealing(par_sig_0, temperature),
                            simulated_annealing(par_sig_1, temperature),
                            simulated_annealing(par_sig_2, temperature),
                            simulated_annealing(par_sig_3, temperature),
                        ]

                        # Insert tabu list check here
                        if tuple(current_params) not in tabu_list:
                            tabu_list.add(tuple(current_params))
                        else:
                            # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        g_sig = TGraphErrors()
                        for i in range(nsep.GetSelectedRows()):
                            g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                            g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                        for i in range(len(w_vec)):
                            if sig_name == "LT":
                                sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                            elif sig_name == "TT":
                                sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                                sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                            else:
                                sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                                sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                            graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                        if sig_name == "L":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
                        elif sig_name == "T":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
                        elif sig_name == "LT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
                        elif sig_name == "TT":
                            #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                            f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
                        f_sig.SetParNames("p0", "p1", "p2", "p3")
                        f_sig.SetParameter(0, current_params[0])
                        f_sig.SetParameter(1, current_params[1])
                        f_sig.SetParameter(2, current_params[2])
                        f_sig.SetParameter(3, current_params[3])
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)
                        f_sig.SetParLimits(2, -max_param_value, max_param_value)
                        f_sig.SetParLimits(3, -max_param_value, max_param_value)                

                        r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                        #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                        f_sig_status = f_sig.GetNDF() != 0

                        params_sig_history['p0'].append(current_params[0])
                        params_sig_history['p1'].append(current_params[1])
                        params_sig_history['p2'].append(current_params[2])
                        params_sig_history['p3'].append(current_params[3])

                        if num_events > num_params:
                            # Calculate the cost (reduced chi-square value) for the current parameters                            
                            current_cost = f_sig.GetChisquare()/(num_events-num_params) # Divided by DoF for red. chi-squared
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, current_cost, temperature)
                        else:
                            # Define lambda values to try (in log space for better exploration)
                            lambda_values = np.logspace(np.log10(lambda_min), np.log10(lambda_max), 10)
                            best_cost_iteration = float('inf')  # Initialize to a large value
                            best_lambda = lambda_reg  # Start with initial lambda                            
                            # Adaptive search for the best lambda
                            for lambda_try in lambda_values:
                                residuals = []  # Store residuals for this lambda value
                                # Loop through each data point (event)
                                for i in range(num_events):
                                    observed = g_sig.GetY()[i]  # Observed value
                                    expected = f_sig.Eval(g_sig.GetX()[i])  # Model prediction
                                    # Compute residuals (normalized by error if available)
                                    if g_sig.GetEY()[i] != 0:
                                        residual = (observed - expected) / g_sig.GetEY()[i]
                                    else:
                                        residual = observed - expected
                                    residuals.append(residual)
                                # Calculate the mean squared error (MSE)
                                mse = np.mean(np.square(residuals))
                                # Compute L2 regularization term (sum of parameter squares)
                                l2_reg = sum(p**2 for p in current_params)
                                # Generalized reduced chi-squared formula
                                current_cost_try = (mse + lambda_try * l2_reg) / (num_events + alpha * num_params)
                                # Check if this lambda gives the best cost so far
                                if current_cost_try < best_cost_iteration:
                                    best_cost_iteration = current_cost_try
                                    best_lambda = lambda_try
                                    best_params = current_params.copy()  # Save parameters for this lambda
                                # Optional: Early stopping if the improvement is small
                                if abs(current_cost_try - best_cost_iteration) < tolerance:
                                    break
                            # Use the best lambda and parameters found
                            lambda_reg = best_lambda
                            current_params = best_params  # Update model parameters
                            current_cost = best_cost_iteration  # Final cost for this iteration
                            # Store cost history for analysis
                            cost_history.append(current_cost)
                            # Calculate acceptance probability for simulated annealing (if applicable)
                            accept_prob = acceptance_probability(current_cost, best_cost_iteration, temperature)
                            
                        current_params = [
                            f_sig.GetParameter(0),
                            f_sig.GetParameter(1),
                            f_sig.GetParameter(2),
                            f_sig.GetParameter(3)
                        ]

                        current_errors = [
                            f_sig.GetParError(0),
                            f_sig.GetParError(1),
                            f_sig.GetParError(2),
                            f_sig.GetParError(3)
                        ]

                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, current_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, current_params[1])
                        graphs_sig_p2[it].SetPoint(total_iteration, total_iteration, current_params[2])
                        graphs_sig_p3[it].SetPoint(total_iteration, total_iteration, current_params[3])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                        # If the new cost is better or accepted by the acceptance probability, update the best parameters
                        if accept_prob > random.random():
                            best_params = current_params
                            best_cost = current_cost
                            best_errors = current_errors

                        if iteration % local_search_interval == 0:
                            current_params = local_search(current_params, f_sig, num_params)
                            par_sig_0, par_sig_1, par_sig_2, par_sig_3 = current_params

                        # Check if current parameters haven't changed for the past N iterations
                        if len(params_sig_history['p0']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p1']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p2']) >= max_unchanged_iterations  and \
                           len(params_sig_history['p3']) >= max_unchanged_iterations:
                            if np.allclose(round(params_sig_history['p0'][-2], 3), round(params_sig_history['p0'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p1'][-2], 3), round(params_sig_history['p1'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p2'][-2], 3), round(params_sig_history['p2'][-1], 3), atol=5.0) and \
                               np.allclose(round(params_sig_history['p3'][-2], 3), round(params_sig_history['p3'][-1], 3), atol=5.0):
                                unchanged_iterations += 1
                            else:
                                unchanged_iterations = 0

                        # Adjust the cooling rate if parameters haven't changed for N iterations
                        if unchanged_iterations >= max_unchanged_iterations:
                            if not any(np.allclose([current_params[0], current_params[1], current_params[2], current_params[3]], minima, atol=5.0) for minima in local_minima):                    
                                local_minima.append([
                                    current_params[0],
                                    current_params[1],
                                    current_params[2],
                                    current_params[3]
                                ])
                                # Restart from initial parameters
                            current_params = initial_params
                            temperature = initial_temperature
                            unchanged_iterations = 0

                        previous_params = current_params[:]                

                        # Update parameters with the best found so far
                        par_sig_0, par_sig_1, par_sig_2, par_sig_3 = best_params
                        par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3 = best_errors

                        # Update the temperature
                        temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0

                        # Check if current_params are close to any local minimum
                        if any(np.allclose([current_params[0], current_params[1], current_params[2], current_params[3]], minima, atol=5.0) for minima in local_minima):
                            #print("WARNING: Parameters p0={:.3e}, p1={:.3e}, p2={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2]))

                            current_params = adjust_params(best_params)
                            par_sig_0, par_sig_1, par_sig_2, par_sig_3 = current_params
                            par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3 = [0.0 for _ in range(num_params)]

                    except (TypeError or ZeroDivisionError) as e:
                        #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                        # Adjust parameter limits within a random number
                        par_sig_0, par_sig_1, par_sig_2, par_sig_3 = initial_params
                        par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3 = [0.0 for _ in range(num_params)]

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0                

                # After the while loop, check if this run found a better solution
                if abs(best_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = best_cost
                    best_overall_params = best_params[:]
                    best_overall_errors = best_errors[:]
                
            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            try:            
                while len(best_overall_params) < 4:
                    best_overall_params.append(0.0)
                    best_overall_errors.append(0.0)
            except TypeError:
                print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                sys.exit(2)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[1])
            par_err_vec.append(best_overall_errors[2])
            par_err_vec.append(best_overall_errors[3])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_fit = TGraphErrors()

            graphs_sig_fit.append(g_sig_fit)
                        
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                if sig_name == "LT":
                    sig_X_fit = (g_sig.GetY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (math.sin(th_vec[i] * PI / 180)**sine_exp_LT * (g_vec[i]))
                if sig_name == "TT":
                    sig_X_fit = (g_sig.GetY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i]))
                    sig_X_fit_err = (g_sig.GetEY()[i] * math.sin(th_vec[i] * PI / 180)**sine_exp_TT * (g_vec[i])
                else:
                    sig_X_fit = (g_sig.GetY()[i]) / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i]) / (g_vec[i])

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s } [nb/GeV^{2}]" % sig_name)
            graphs_sig_fit[it].GetYaxis().SetTitleOffset(1.5)
            graphs_sig_fit[it].GetYaxis().SetTitleSize(0.035)
            graphs_sig_fit[it].GetYaxis().CenterTitle()

            # Set axis limits to ensure everything is shown
            x_min = min(graphs_sig_fit[it].GetX())
            x_max = max(graphs_sig_fit[it].GetX())
            y_min = min(graphs_sig_fit[it].GetY())
            y_max = max(graphs_sig_fit[it].GetY())

            # You can also set a margin to ensure all points are visible
            margin = 0.1
            graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

            if sig_name == "L":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)    
                f_sig.SetParNames("p0", "p1", "p2", "p3")
            f_sig.FixParameter(0, best_overall_params[0])
            f_sig.FixParameter(1, best_overall_params[1])
            f_sig.FixParameter(2, best_overall_params[2])
            f_sig.FixParameter(3, best_overall_params[3])
        
            # Evaluate the fit function at several points to determine its range
            n_points = 100  # Number of points to evaluate the fit function
            fit_y_values = [f_sig.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)

            # Extend the y-axis range to include the fit function range
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)

            # Set a margin to ensure all points are visible
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

            r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
            f_sig.Draw("same")
            
            #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
            f_sig_status = f_sig.GetNDF() != 0
            f_sig_status_message = "Fit Successful" if f_sig_status else "Fit Failed"

            fit_status = TText()
            fit_status.SetTextSize(0.04)
            fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sig_status_message))
            
            # Calculate the minimum and maximum values from the graphs
            min_sig_y = float('inf')
            max_sig_y = float('-inf')

            # Update min_sig_y and max_sig_y based on each graph's values
            for graph in [graphs_sig_p0[it], graphs_sig_p1[it], graphs_sig_p2[it], graphs_sig_p3[it]]:
                n_points = graph.GetN()
                for i in range(n_points):
                    y = graph.GetY()[i]
                    if y < min_sig_y:
                        min_sig_y = y
                    if y > max_sig_y:
                        max_sig_y = y

            # Scale the y-axis
            graphs_sig_p0[it].SetMinimum(min_sig_y * 0.9)
            graphs_sig_p0[it].SetMaximum(max_sig_y * 1.1)
            c2.Update()            

            # Plot parameter convergence
            c3.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_p0[it].SetTitle(f"Sig {sig_name} Parameter Convergence;Optimization Run;Parameter")
            graphs_sig_p0[it].SetLineColor(ROOT.kRed)
            graphs_sig_p1[it].SetLineColor(ROOT.kBlue)
            graphs_sig_p2[it].SetLineColor(ROOT.kGreen)
            graphs_sig_p3[it].SetLineColor(ROOT.kPink)
            graphs_sig_p0[it].Draw("ALP")
            graphs_sig_p1[it].Draw("LP SAME")
            graphs_sig_p2[it].Draw("LP SAME")
            graphs_sig_p3[it].Draw("LP SAME")
            c3.Update()
            
            # Plot chi-square convergence
            c4.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
            graphs_sig_converge[it].SetLineColor(ROOT.kBlack)
            graphs_sig_converge[it].Draw("ALP")
            converge_status = TText()
            converge_status.SetTextSize(0.04)
            converge_status.DrawTextNDC(0.35, 0.85, f"Best cost: {best_overall_cost:.3f}")
            c4.Update()
            
            # Plot temperature
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability
            c6.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(ROOT.kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()
            
            print("\n")

        c2.Update()
        c3.Update()
        c4.Update()
        c5.Update()
        c6.Update()
    
    c2.Print(outputpdf+'(')
    c3.Print(outputpdf)
    c4.Print(outputpdf)
    c5.Print(outputpdf)
    c6.Print(outputpdf+')')
