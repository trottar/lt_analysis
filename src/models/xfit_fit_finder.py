#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-11 13:51:38 trottar"
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
from utility import adaptive_regularization, calculate_cost, adaptive_cooling, simulated_annealing, acceptance_probability, adjust_params, local_search, select_valid_parameter, get_central_value

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

def parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization=True):

    # Create lists to store graph objects outside the loop
    graphs_sig_fit = []
    graphs_sig_p0 = []
    graphs_sig_p1 = []
    graphs_sig_p2 = []
    graphs_sig_p3 = []
    graphs_sig_chi2 = []
    graphs_sig_temp = []
    graphs_sig_accept = []
    graphs_sig_converge = []
    
    c2 = TCanvas("c2", "c2", 800, 800)
    c2.Divide(2, 2)

    # Create ROOT canvases for additional parameter convergence plots
    c3 = TCanvas("c3", "Parameter Convergence", 800, 800)
    c3.Divide(2, 2)
    c4 = TCanvas("c4", "Red. Chi-Square Convergence", 800, 800)
    c4.Divide(2, 2)
    c5 = TCanvas("c5", "Temperature", 800, 800)
    c5.Divide(2, 2)
    c6 = TCanvas("c6", "Acceptance Probability", 800, 800)
    c6.Divide(2, 2)

    q2_set = inpDict["q2_set"]
    w_set = inpDict["w_set"]
    nsep, t_vec, g_vec, w_vec, q2_vec, th_vec = inpDict["objects"]
    max_iterations = inpDict["max_iterations"]
    num_optimizations = inpDict["num_optimizations"]
    initial_param_bounds = inpDict["initial_param_bounds"]
    tmin_range = inpDict["tmin_range"]
    tmax_range = inpDict["tmax_range"]
    Q2min_range = inpDict["Q2min_range"]
    Q2max_range = inpDict["Q2max_range"]
    iter_num = inpDict["iter_num"]            
    fit_params = inpDict["fit_params"]
    chi2_threshold = inpDict["chi2_threshold"]
    
    # Using central bin value to determine best fit, which should have the best statistics
    q2_center_val = get_central_value(q2_vec)
    w_center_val = get_central_value(w_vec)
    g_center_val = get_central_value(g_vec)
    th_center_val = get_central_value(th_vec)
    #print(f"\n\nDetermining best fit off the central bin values...\n Q2={q2_center_val:.3f}, W={w_center_val:.3f}, theta={th_center_val:.3f}")
    #fun_Sig_L = fun_Sig_L_wrapper(q2_center_val, w_center_val)
    #fun_Sig_T = fun_Sig_T_wrapper(q2_center_val, w_center_val)
    #fun_Sig_LT = fun_Sig_LT_wrapper(q2_center_val, w_center_val, th_center_val)
    #fun_Sig_TT = fun_Sig_TT_wrapper(q2_center_val, w_center_val, th_center_val)

    num_events = nsep.GetEntries()

    for it, (key, val) in enumerate(fit_params.items()):

        # Don't find new fits if debugging
        if key not in fixed_params:
            
            # Find optimized fits
            sig_name = key
            # Grab parameters used by functional forms
            num_params, initial_params, equation_str = inpDict["initial_params"](sig_name, val)

            # Checks initial parameters and replaces zeros to avoid errors
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
                best_overall_bin = None
                total_iteration = 0
                max_param_bounds = initial_param_bounds

                # Regularization strength (used when num_events > num_params)
                # Initialize adaptive regularization parameters
                lambda_reg = 0.01  # Initial regularization strength
                cost_history = []

                # Parameter range offset (% of param value)
                param_offset_0 = 0.1
                param_offset_1 = 0.1                
                
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
                    print("\n\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                    set_optimization = full_optimization
                    
                    for b in range(len(w_vec)):

                        print(f"Determining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")

                        iteration = 0

                        initial_temperature = 1.0
                        temperature = initial_temperature
                        unchanged_iterations = 0
                        max_unchanged_iterations = 5

                        # Initialize adaptive parameter limits
                        par_sig_0 = initial_params[0]
                        par_sig_err_0 = 0.0

                        # Track the best solution
                        best_params = [par_sig_0]
                        best_cost = float('inf')
                        best_bin = None
                        best_errors = [par_sig_err_0]
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

                                current_params = [simulated_annealing(par_sig_0, temperature)]

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
                                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                                    graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                                    graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                                if sig_name == "L":
                                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                                elif sig_name == "T":
                                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                                elif sig_name == "LT":
                                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                                elif sig_name == "TT":
                                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
                                f_sig.SetParNames("p0")
                                f_sig.SetParameter(0, current_params[0])
                                if set_optimization:
                                    f_sig.SetParLimits(0, -max_param_bounds, max_param_bounds)
                                else:
                                    f_sig.SetParLimits(0, current_params[0]-param_offset_0*current_params[0], current_params[0]+param_offset_0*current_params[0])

                                r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                                #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                                f_sig_status = f_sig.GetNDF() != 0

                                params_sig_history['p0'].append(current_params[0])

                                # Calculate cost with consistent regularization
                                current_cost, lambda_reg = calculate_cost(
                                    f_sig, g_sig, current_params,
                                    num_events, num_params, lambda_reg
                                )
                                # Store cost for history
                                cost_history.append(current_cost)            
                                # Adapt regularization strength based on history
                                if len(cost_history) >= 2:
                                    lambda_reg = adaptive_regularization(cost_history, lambda_reg)
                                # Update acceptance probability for simulated annealing
                                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                                current_params = [f_sig.GetParameter(0)]
                                current_errors = [f_sig.GetParError(0)]
                                current_bin = b

                                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                                if accept_prob > random.random():
                                    best_params = current_params
                                    best_cost = current_cost
                                    best_bin = current_bin
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
                                    if not any(np.allclose([current_params[0]], minima, atol=5.0) for minima in local_minima):                    
                                        local_minima.append([
                                            current_params[0]
                                        ])
                                    # Restart from initial parameters
                                    current_params = initial_params
                                    temperature = initial_temperature
                                    unchanged_iterations = 0

                                previous_params = current_params[:]

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
                                    par_sig_err_0 = [0.0]

                            except (TypeError or ZeroDivisionError) as e:
                                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                                # Generate safer parameter values within reasonable bounds
                                recovery_params = [
                                    p + random.uniform(-0.1 * abs(p), 0.1 * abs(p)) 
                                    for p in (best_params if best_params != [float('inf')] * len(initial_params) else initial_params)
                                ]

                                # Ensure parameters stay within bounds
                                recovery_params = [
                                    max(min(p, max_param_bounds), -max_param_bounds) 
                                    for p in recovery_params
                                ]

                                # Reset function parameters
                                for i, param in enumerate(recovery_params):
                                    f_sig.SetParameter(i, param)

                                # Don't update best_cost to inf, keep previous best
                                current_params = recovery_params
                                current_cost = best_cost * 1.1 if math.isfinite(best_cost) else 1000.0

                                # Increase temperature slightly to encourage exploration
                                temperature = min(temperature * 1.2, initial_temperature)

                                max_param_bounds = max_param_bounds * random.random()
                                iteration += 1
                                total_iteration += 1 if iteration % max_iterations == 0 else 0                
                        
                            # After the while loop, check if this run found a better solution
                            if abs(best_cost - 1) < abs(best_overall_cost - 1):
                                best_overall_cost = best_cost
                                best_overall_bin = best_bin
                                best_overall_params = best_params[:]
                                best_overall_errors = best_errors[:]
                                if best_overall_cost < chi2_threshold:
                                    set_optimization = False                                    
                                    
                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, best_overall_params[0])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(best_overall_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                        print(f"\nBest Cost: {best_overall_cost:.3f}")
                    
                try:
                    print(f"\n\nBest overall solution: {best_overall_params}")
                    print(f"Best overall cost: {best_overall_cost:.5f}")
                    print(f"Best overall bin: t={t_vec[best_overall_bin]:.3f}, Q2={q2_vec[best_overall_bin]:.3f}, W={w_vec[best_overall_bin]:.3f}, theta={th_vec[best_overall_bin]:.3f}")
                except TypeError:
                    print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                    sys.exit(2)
                    
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
                    
                for j in range(4):
                    par_vec[4*it+j] = best_overall_params[j]
                    par_err_vec[4*it+j] = best_overall_errors[j]
                    par_chi2_vec[4*it+j] = best_overall_cost

                g_sig_fit = TGraphErrors()

                graphs_sig_fit.append(g_sig_fit)

                g_sig = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
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
                # Create a TLatex object
                latex = ROOT.TLatex()
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
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
                best_overall_bin = None
                total_iteration = 0
                max_param_bounds = initial_param_bounds

                # Regularization strength (used when num_events > num_params)
                # Initialize adaptive regularization parameters
                lambda_reg = 0.01  # Initial regularization strength
                cost_history = []

                # Parameter range offset (% of param value)
                param_offset_0 = 0.1
                param_offset_1 = 0.1                                

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
                    print("\n\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))

                    set_optimization = full_optimization
                    
                    for b in range(len(w_vec)):

                        print(f"Determining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")

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
                        best_bin = None
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
                                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                                    graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                                    graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                                if sig_name == "L":
                                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                                elif sig_name == "T":
                                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                                elif sig_name == "LT":
                                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                                elif sig_name == "TT":
                                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
                                f_sig.SetParNames("p0", "p1")
                                f_sig.SetParameter(0, current_params[0])
                                f_sig.SetParameter(1, current_params[1])
                                if set_optimization:
                                    f_sig.SetParLimits(0, -max_param_bounds, max_param_bounds)
                                    f_sig.SetParLimits(1, -max_param_bounds, max_param_bounds)
                                else:
                                    f_sig.SetParLimits(0, current_params[0]-param_offset_0*current_params[0], current_params[0]+param_offset_0*current_params[0])
                                    f_sig.SetParLimits(1, current_params[1]-param_offset_1*current_params[1], current_params[1]+param_offset_1*current_params[1])
                                
                                r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                                #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                                f_sig_status = f_sig.GetNDF() != 0

                                params_sig_history['p0'].append(current_params[0])
                                params_sig_history['p1'].append(current_params[1])

                                # Calculate cost with consistent regularization
                                current_cost, lambda_reg = calculate_cost(
                                    f_sig, g_sig, current_params,
                                    num_events, num_params, lambda_reg
                                )
                                # Store cost for history
                                cost_history.append(current_cost)            
                                # Adapt regularization strength based on history
                                if len(cost_history) >= 2:
                                    lambda_reg = adaptive_regularization(cost_history, lambda_reg)            
                                # Update acceptance probability for simulated annealing
                                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                                current_params = [
                                    f_sig.GetParameter(0),
                                    f_sig.GetParameter(1)
                                ]

                                current_errors = [
                                    f_sig.GetParError(0),
                                    f_sig.GetParError(1)
                                ]

                                current_bin = b

                                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                                if accept_prob > random.random():
                                    best_params = current_params
                                    best_cost = current_cost
                                    best_bin = current_bin
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
                                # Generate safer parameter values within reasonable bounds
                                recovery_params = [
                                    p + random.uniform(-0.1 * abs(p), 0.1 * abs(p)) 
                                    for p in (best_params if best_params != [float('inf')] * len(initial_params) else initial_params)
                                ]

                                # Ensure parameters stay within bounds
                                recovery_params = [
                                    max(min(p, max_param_bounds), -max_param_bounds) 
                                    for p in recovery_params
                                ]

                                # Reset function parameters
                                for i, param in enumerate(recovery_params):
                                    f_sig.SetParameter(i, param)

                                # Don't update best_cost to inf, keep previous best
                                current_params = recovery_params
                                current_cost = best_cost * 1.1 if math.isfinite(best_cost) else 1000.0

                                # Increase temperature slightly to encourage exploration
                                temperature = min(temperature * 1.2, initial_temperature)

                                max_param_bounds = max_param_bounds * random.random()
                                iteration += 1
                                total_iteration += 1 if iteration % max_iterations == 0 else 0
                        
                            # After the while loop, check if this run found a better solution
                            if abs(best_cost - 1) < abs(best_overall_cost - 1):
                                best_overall_cost = best_cost
                                best_overall_bin = best_bin
                                best_overall_params = best_params[:]
                                best_overall_errors = best_errors[:]
                                if best_overall_cost < chi2_threshold:
                                    set_optimization = False                                    
                                    
                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, best_overall_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, best_overall_params[1])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(best_overall_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                        print(f"\nBest Cost: {best_overall_cost:.3f}")
                    
                try:
                    print(f"\n\nBest overall solution: {best_overall_params}")
                    print(f"Best overall cost: {best_overall_cost:.5f}")
                    print(f"Best overall bin: t={t_vec[best_overall_bin]:.3f}, Q2={q2_vec[best_overall_bin]:.3f}, W={w_vec[best_overall_bin]:.3f}, theta={th_vec[best_overall_bin]:.3f}")
                except TypeError:
                    print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                    sys.exit(2)
                    
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
                    
                for j in range(4):
                    par_vec[4*it+j] = best_overall_params[j]
                    par_err_vec[4*it+j] = best_overall_errors[j]
                    par_chi2_vec[4*it+j] = best_overall_cost

                g_sig_fit = TGraphErrors()

                graphs_sig_fit.append(g_sig_fit)

                g_sig = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
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
                # Create a TLatex object
                latex = ROOT.TLatex()
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
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
                best_overall_bin = None
                total_iteration = 0
                max_param_bounds = initial_param_bounds

                # Regularization strength (used when num_events > num_params)
                # Initialize adaptive regularization parameters
                lambda_reg = 0.01  # Initial regularization strength
                cost_history = []

                # Parameter range offset (% of param value)
                param_offset_0 = 0.1
                param_offset_1 = 0.1
                param_offset_2 = 0.1

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
                    print("\n\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                    set_optimization = full_optimization
                    
                    for b in range(len(w_vec)):

                        print(f"Determining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")

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
                        best_bin = None
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
                                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                                    graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                                    graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                                if sig_name == "L":
                                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                                elif sig_name == "T":
                                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                                elif sig_name == "LT":
                                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                                elif sig_name == "TT":
                                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
                                f_sig.SetParNames("p0", "p1", "p2")
                                f_sig.SetParameter(0, current_params[0])
                                f_sig.SetParameter(1, current_params[1])
                                f_sig.SetParameter(2, current_params[2])
                                if set_optimization:
                                    f_sig.SetParLimits(0, -max_param_bounds, max_param_bounds)
                                    f_sig.SetParLimits(1, -max_param_bounds, max_param_bounds)                                
                                    f_sig.SetParLimits(2, -max_param_bounds, max_param_bounds)
                                else:
                                    f_sig.SetParLimits(0, current_params[0]-param_offset_0*current_params[0], current_params[0]+param_offset_0*current_params[0])
                                    f_sig.SetParLimits(1, current_params[1]-param_offset_1*current_params[1], current_params[1]+param_offset_1*current_params[1])
                                    f_sig.SetParLimits(2, current_params[2]-param_offset_2*current_params[2], current_params[2]+param_offset_2*current_params[2])

                                r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                                #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                                f_sig_status = f_sig.GetNDF() != 0

                                params_sig_history['p0'].append(current_params[0])
                                params_sig_history['p1'].append(current_params[1])
                                params_sig_history['p2'].append(current_params[2])

                                # Calculate cost with consistent regularization
                                current_cost, lambda_reg = calculate_cost(
                                    f_sig, g_sig, current_params,
                                    num_events, num_params, lambda_reg
                                )
                                # Store cost for history
                                cost_history.append(current_cost)            
                                # Adapt regularization strength based on history
                                if len(cost_history) >= 2:
                                    lambda_reg = adaptive_regularization(cost_history, lambda_reg)            
                                # Update acceptance probability for simulated annealing
                                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

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

                                current_bin = b

                                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                                if accept_prob > random.random():
                                    best_params = current_params
                                    best_cost = current_cost
                                    best_bin = current_bin
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
                                # Generate safer parameter values within reasonable bounds
                                recovery_params = [
                                    p + random.uniform(-0.1 * abs(p), 0.1 * abs(p)) 
                                    for p in (best_params if best_params != [float('inf')] * len(initial_params) else initial_params)
                                ]

                                # Ensure parameters stay within bounds
                                recovery_params = [
                                    max(min(p, max_param_bounds), -max_param_bounds) 
                                    for p in recovery_params
                                ]

                                # Reset function parameters
                                for i, param in enumerate(recovery_params):
                                    f_sig.SetParameter(i, param)

                                # Don't update best_cost to inf, keep previous best
                                current_params = recovery_params
                                current_cost = best_cost * 1.1 if math.isfinite(best_cost) else 1000.0

                                # Increase temperature slightly to encourage exploration
                                temperature = min(temperature * 1.2, initial_temperature)

                                max_param_bounds = max_param_bounds * random.random()
                                iteration += 1
                                total_iteration += 1 if iteration % max_iterations == 0 else 0
                        
                            # After the while loop, check if this run found a better solution
                            if abs(best_cost - 1) < abs(best_overall_cost - 1):
                                best_overall_cost = best_cost
                                best_overall_bin = best_bin
                                best_overall_params = best_params[:]
                                best_overall_errors = best_errors[:]
                                if best_overall_cost < chi2_threshold:
                                    set_optimization = False                                    
                                    
                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, best_overall_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, best_overall_params[1])
                        graphs_sig_p2[it].SetPoint(total_iteration, total_iteration, best_overall_params[2])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(best_overall_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                        print(f"\nBest Cost: {best_overall_cost:.3f}")
                    
                try:
                    print(f"\n\nBest overall solution: {best_overall_params}")
                    print(f"Best overall cost: {best_overall_cost:.5f}")
                    print(f"Best overall bin: t={t_vec[best_overall_bin]:.3f}, Q2={q2_vec[best_overall_bin]:.3f}, W={w_vec[best_overall_bin]:.3f}, theta={th_vec[best_overall_bin]:.3f}")
                except TypeError:
                    print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                    sys.exit(2)
                    
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
                    
                for j in range(4):
                    par_vec[4*it+j] = best_overall_params[j]
                    par_err_vec[4*it+j] = best_overall_errors[j]
                    par_chi2_vec[4*it+j] = best_overall_cost

                g_sig_fit = TGraphErrors()

                graphs_sig_fit.append(g_sig_fit)

                g_sig = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
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
                # Create a TLatex object
                latex = ROOT.TLatex()
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
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
                best_overall_bin = None
                total_iteration = 0
                max_param_bounds = initial_param_bounds

                # Regularization strength (used when num_events > num_params)
                # Initialize adaptive regularization parameters
                lambda_reg = 0.01  # Initial regularization strength
                cost_history = []

                # Parameter range offset (% of param value)
                param_offset_0 = 0.1
                param_offset_1 = 0.1                                
                param_offset_2 = 0.1
                param_offset_3 = 0.1
                
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
                    print("\n\nStarting optimization run {0}/{1}".format(start + 1, num_optimizations))    

                    set_optimization = full_optimization
                    
                    for b in range(len(w_vec)):

                        print(f"Determining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")

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
                        best_bin = None
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
                                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                                    graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                                    graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                                if sig_name == "L":
                                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                                elif sig_name == "T":
                                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                                elif sig_name == "LT":
                                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                                elif sig_name == "TT":
                                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
                                f_sig.SetParNames("p0", "p1", "p2", "p3")
                                f_sig.SetParameter(0, current_params[0])
                                f_sig.SetParameter(1, current_params[1])
                                f_sig.SetParameter(2, current_params[2])
                                f_sig.SetParameter(3, current_params[3])
                                if set_optimization:
                                    f_sig.SetParLimits(0, -max_param_bounds, max_param_bounds)                                
                                    f_sig.SetParLimits(2, -max_param_bounds, max_param_bounds)
                                    f_sig.SetParLimits(3, -max_param_bounds, max_param_bounds)
                                    f_sig.SetParLimits(1, -max_param_bounds, max_param_bounds)
                                else:
                                    f_sig.SetParLimits(0, current_params[0]-param_offset_0*current_params[0], current_params[0]+param_offset_0*current_params[0])
                                    f_sig.SetParLimits(1, current_params[1]-param_offset_1*current_params[1], current_params[1]+param_offset_1*current_params[1])
                                    f_sig.SetParLimits(2, current_params[2]-param_offset_2*current_params[2], current_params[2]+param_offset_2*current_params[2])
                                    f_sig.SetParLimits(3, current_params[3]-param_offset_3*current_params[3], current_params[3]+param_offset_3*current_params[3])                                
                                
                                r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                                #f_sig_status = (r_sig_fit.Status() == 0 and r_sig_fit.IsValid())
                                f_sig_status = f_sig.GetNDF() != 0

                                params_sig_history['p0'].append(current_params[0])
                                params_sig_history['p1'].append(current_params[1])
                                params_sig_history['p2'].append(current_params[2])
                                params_sig_history['p3'].append(current_params[3])

                                # Calculate cost with consistent regularization
                                current_cost, lambda_reg = calculate_cost(
                                    f_sig, g_sig, current_params,
                                    num_events, num_params, lambda_reg
                                )
                                # Store cost for history
                                cost_history.append(current_cost)            
                                # Adapt regularization strength based on history
                                if len(cost_history) >= 2:
                                    lambda_reg = adaptive_regularization(cost_history, lambda_reg)
                                # Update acceptance probability for simulated annealing
                                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

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

                                current_bin = b

                                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                                if accept_prob > random.random():
                                    best_params = current_params
                                    best_cost = current_cost
                                    best_bin = current_bin
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
                                # Generate safer parameter values within reasonable bounds
                                recovery_params = [
                                    p + random.uniform(-0.1 * abs(p), 0.1 * abs(p)) 
                                    for p in (best_params if best_params != [float('inf')] * len(initial_params) else initial_params)
                                ]

                                # Ensure parameters stay within bounds
                                recovery_params = [
                                    max(min(p, max_param_bounds), -max_param_bounds) 
                                    for p in recovery_params
                                ]

                                # Reset function parameters
                                for i, param in enumerate(recovery_params):
                                    f_sig.SetParameter(i, param)

                                # Don't update best_cost to inf, keep previous best
                                current_params = recovery_params
                                current_cost = best_cost * 1.1 if math.isfinite(best_cost) else 1000.0

                                # Increase temperature slightly to encourage exploration
                                temperature = min(temperature * 1.2, initial_temperature)

                                max_param_bounds = max_param_bounds * random.random()
                                iteration += 1
                                total_iteration += 1 if iteration % max_iterations == 0 else 0                
                        
                            # After the while loop, check if this run found a better solution
                            if abs(best_cost - 1) < abs(best_overall_cost - 1):
                                best_overall_cost = best_cost
                                best_overall_bin = best_bin
                                best_overall_params = best_params[:]
                                best_overall_errors = best_errors[:]
                                if best_overall_cost < chi2_threshold:
                                    set_optimization = False
                                    
                        # Update ROOT TGraphs for plotting
                        graphs_sig_p0[it].SetPoint(total_iteration, total_iteration, best_overall_params[0])
                        graphs_sig_p1[it].SetPoint(total_iteration, total_iteration, best_overall_params[1])
                        graphs_sig_p2[it].SetPoint(total_iteration, total_iteration, best_overall_params[2])
                        graphs_sig_p3[it].SetPoint(total_iteration, total_iteration, best_overall_params[3])
                        graphs_sig_converge[it].SetPoint(total_iteration, total_iteration, round(best_overall_cost, 4))
                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)
                        graphs_sig_accept[it].SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                        print(f"\nBest Cost: {best_overall_cost:.3f}")
                    
                try:
                    print(f"\n\nBest overall solution: {best_overall_params}")
                    print(f"Best overall cost: {best_overall_cost:.5f}")
                    print(f"Best overall bin: t={t_vec[best_overall_bin]:.3f}, Q2={q2_vec[best_overall_bin]:.3f}, W={w_vec[best_overall_bin]:.3f}, theta={th_vec[best_overall_bin]:.3f}")
                except TypeError:
                    print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                    sys.exit(2)
                    
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
                    
                for j in range(4):
                    par_vec[4*it+j] = best_overall_params[j]
                    par_err_vec[4*it+j] = best_overall_errors[j]
                    par_chi2_vec[4*it+j] = best_overall_cost

                g_sig_fit = TGraphErrors()

                graphs_sig_fit.append(g_sig_fit)

                g_sig = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                    sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 3.0, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 3.0, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 3.0, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    #f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 3.0, num_params)
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
                # Create a TLatex object
                latex = ROOT.TLatex()
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
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
            
        ## plot_fit
        else:

            sig_name = key
            # Grab parameters used by functional forms
            num_params, initial_params, equation_str = inpDict["initial_params"](sig_name, val)

            # Checks initial parameters and replaces zeros to avoid errors
            initial_params = [v if abs(v) > 0.0 else 1.0 for v in initial_params]

            # String list of initial parameters
            param_str = ', '.join(str(param) for param in initial_params)

            fit_convergence_type = "Fixed"

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

                nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")

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
                
                best_overall_cost = float('inf')  # Initialize with infinity
                best_overall_params = []
                best_overall_bin = 0  # Initialize best bin tracker
                
                for b in range(len(w_vec)):
                    
                    print(f"\n\nDetermining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    
                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    g_sig = TGraphErrors()
                    for i in range(nsep.GetSelectedRows()):
                        g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                        g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                    for i in range(len(w_vec)):
                        sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                        sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                        graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                        graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                    if sig_name == "L":
                        fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    elif sig_name == "T":
                        fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    elif sig_name == "LT":
                        fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    elif sig_name == "TT":
                        fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    f_sig.SetParNames("p0")
                    f_sig.FixParameter(0, par_vec[4*it])

                    r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                    # Retrieve the chi-squared and degrees of freedom
                    chi2 = f_sig.GetChisquare()  # Get the chi-squared value
                    ndf = f_sig.GetNDF()         # Get the number of degrees of freedom
                    red_chi2 = chi2 / ndf    # Calculate reduced chi-squared

                    print(f"\tCost: {red_chi2:.3f}")
                    
                    # After the while loop, check if this run found a better solution
                    if abs(red_chi2 - 1) < abs(best_overall_cost - 1):
                        best_overall_cost = red_chi2
                        best_overall_bin = b
                        for j in range(4):
                            best_overall_params.append(par_vec[4*it+j])
                            
                print(f"\n\nBest overall solution: {best_overall_params}")
                print(f"Best overall cost: {best_overall_cost:.5f}")
                    
                for j in range(4):
                    par_chi2_vec[4*it+j] = best_overall_cost

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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig.SetParNames("p0")
                f_sig.FixParameter(0, par_vec[4*it])

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

                # Create a TLatex object
                latex = ROOT.TLatex()# Set the font size and alignment (optional but recommended)
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates# Format the text
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
                c2.Update()
                    
                print("\n")    

            elif num_params == 2:

                # 2 param
                #######
                # Sig #
                #######

                print("\n/*--------------------------------------------------*/")
                print(f"Fit for Sig {sig_name} ({num_params} parameters)")
                print(f"Initial Paramters: ({param_str})")
                print(f"{equation_str}")
                print("/*--------------------------------------------------*/")

                nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")

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
                
                best_overall_cost = float('inf')  # Initialize with infinity
                best_overall_params = []
                best_overall_bin = 0  # Initialize best bin tracker
                
                for b in range(len(w_vec)):

                    print(f"\n\nDetermining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    
                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    g_sig = TGraphErrors()
                    for i in range(nsep.GetSelectedRows()):
                        g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                        g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                    for i in range(len(w_vec)):
                        sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                        sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                        graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                        graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                    if sig_name == "L":
                        fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    elif sig_name == "T":
                        fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    elif sig_name == "LT":
                        fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    elif sig_name == "TT":
                        fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    f_sig.SetParNames("p0", "p1")
                    f_sig.FixParameter(0, par_vec[4*it])
                    f_sig.FixParameter(1, par_vec[4*it+1])

                    r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                    # Retrieve the chi-squared and degrees of freedom
                    chi2 = f_sig.GetChisquare()  # Get the chi-squared value
                    ndf = f_sig.GetNDF()         # Get the number of degrees of freedom
                    red_chi2 = chi2 / ndf    # Calculate reduced chi-squared
                    
                    print(f"\tCost: {red_chi2:.3f}")
                    
                    # After the while loop, check if this run found a better solution
                    if abs(red_chi2 - 1) < abs(best_overall_cost - 1):
                        best_overall_cost = red_chi2
                        best_overall_bin = b
                        for j in range(4):
                            best_overall_params.append(par_vec[4*it+j])
                            
                print(f"\n\nBest overall solution: {best_overall_params}")
                print(f"Best overall cost: {best_overall_cost:.5f}")
                    
                for j in range(4):
                    par_chi2_vec[4*it+j] = best_overall_cost

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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig.SetParNames("p0", "p1")
                f_sig.FixParameter(0, par_vec[4*it])
                f_sig.FixParameter(1, par_vec[4*it+1])

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

                # Create a TLatex object
                latex = ROOT.TLatex()# Set the font size and alignment (optional but recommended)
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates# Format the text
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
                c2.Update()
                    
                print("\n")    

            elif num_params == 3:

                # 3 param
                #######
                # Sig #
                #######

                print("\n/*--------------------------------------------------*/")
                print(f"Fit for Sig {sig_name} ({num_params} parameters)")
                print(f"Initial Paramters: ({param_str})")
                print(f"{equation_str}")
                print("/*--------------------------------------------------*/")

                nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")

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
                
                best_overall_cost = float('inf')  # Initialize with infinity
                best_overall_params = []
                best_overall_bin = 0  # Initialize best bin tracker
                
                for b in range(len(w_vec)):
                    
                    print(f"\n\nDetermining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    
                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    g_sig = TGraphErrors()
                    for i in range(nsep.GetSelectedRows()):
                        g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                        g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                    for i in range(len(w_vec)):
                        sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                        sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                        graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                        graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                    if sig_name == "L":
                        fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    elif sig_name == "T":
                        fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    elif sig_name == "LT":
                        fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    elif sig_name == "TT":
                        fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    f_sig.SetParNames("p0", "p1", "p2")
                    f_sig.FixParameter(0, par_vec[4*it])
                    f_sig.FixParameter(1, par_vec[4*it+1])
                    f_sig.FixParameter(2, par_vec[4*it+2])

                    r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                    # Retrieve the chi-squared and degrees of freedom
                    chi2 = f_sig.GetChisquare()  # Get the chi-squared value
                    ndf = f_sig.GetNDF()         # Get the number of degrees of freedom
                    red_chi2 = chi2 / ndf    # Calculate reduced chi-squared

                    print(f"\tCost: {red_chi2:.3f}")
                    
                    # After the while loop, check if this run found a better solution
                    if abs(red_chi2 - 1) < abs(best_overall_cost - 1):
                        best_overall_cost = red_chi2
                        best_overall_bin = b
                        for j in range(4):
                            best_overall_params.append(par_vec[4*it+j])
                            
                print(f"\n\nBest overall solution: {best_overall_params}")
                print(f"Best overall cost: {best_overall_cost:.5f}")
                    
                for j in range(4):
                    par_chi2_vec[4*it+j] = best_overall_cost

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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig.SetParNames("p0", "p1", "p2")
                f_sig.FixParameter(0, par_vec[4*it])
                f_sig.FixParameter(1, par_vec[4*it+1]) 
                f_sig.FixParameter(2, par_vec[4*it+2])

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

                # Create a TLatex object
                latex = ROOT.TLatex()# Set the font size and alignment (optional but recommended)
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates# Format the text
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
                c2.Update()
                    
                print("\n")    

            elif num_params == 4:

                # 4 param
                #######
                # Sig #
                #######

                print("\n/*--------------------------------------------------*/")
                print(f"Fit for Sig {sig_name} ({num_params} parameters)")
                print(f"Initial Paramters: ({param_str})")
                print(f"{equation_str}")
                print("/*--------------------------------------------------*/")

                nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")

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
                
                best_overall_cost = float('inf')  # Initialize with infinity
                best_overall_params = []
                best_overall_bin = 0  # Initialize best bin tracker
                
                for b in range(len(w_vec)):
                    
                    print(f"\n\nDetermining best fit off the bin values...\n t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    
                    g_sig_fit = TGraphErrors()

                    graphs_sig_fit.append(g_sig_fit)

                    g_sig = TGraphErrors()
                    for i in range(nsep.GetSelectedRows()):
                        g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                        g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                    for i in range(len(w_vec)):
                        sig_X_fit = (g_sig.GetY()[i])# / (g_vec[i])
                        sig_X_fit_err = (g_sig.GetEY()[i])# / (g_vec[i])
                        graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                        graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

                    if sig_name == "L":
                        fun_Sig_L = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                    elif sig_name == "T":
                        fun_Sig_T = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                    elif sig_name == "LT":
                        fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                    elif sig_name == "TT":
                        fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                        f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                    f_sig.SetParNames("p0", "p1", "p2", "p3")
                    f_sig.FixParameter(0, par_vec[4*it])
                    f_sig.FixParameter(1, par_vec[4*it+1])
                    f_sig.FixParameter(2, par_vec[4*it+2])
                    f_sig.FixParameter(3, par_vec[4*it+3])                    

                    r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")

                    # Retrieve the chi-squared and degrees of freedom
                    chi2 = f_sig.GetChisquare()  # Get the chi-squared value
                    ndf = f_sig.GetNDF()         # Get the number of degrees of freedom
                    red_chi2 = chi2 / ndf    # Calculate reduced chi-squared

                    print(f"\tCost: {red_chi2:.3f}")
                    
                    # After the while loop, check if this run found a better solution
                    if abs(red_chi2 - 1) < abs(best_overall_cost - 1):
                        best_overall_cost = red_chi2
                        best_overall_bin = b
                        for j in range(4):
                            best_overall_params.append(par_vec[4*it+j])
                            
                print(f"\n\nBest overall solution: {best_overall_params}")
                print(f"Best overall cost: {best_overall_cost:.5f}")
                    
                for j in range(4):
                    par_chi2_vec[4*it+j] = best_overall_cost

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

                # Set a margin to ensure all points are visible
                margin = 0.1
                graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
                graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

                if sig_name == "L":
                    fun_Sig_L = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                elif sig_name == "T":
                    fun_Sig_T = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                elif sig_name == "LT":
                    fun_Sig_LT = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                elif sig_name == "TT":
                    fun_Sig_TT = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig.SetParNames("p0", "p1", "p2", "p3")
                f_sig.FixParameter(0, par_vec[4*it])
                f_sig.FixParameter(1, par_vec[4*it+1]) 
                f_sig.FixParameter(2, par_vec[4*it+2])
                f_sig.FixParameter(3, par_vec[4*it+3])                

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

                # Create a TLatex object
                latex = ROOT.TLatex()# Set the font size and alignment (optional but recommended)
                latex.SetTextSize(0.04)  # Adjust size as needed
                latex.SetNDC(True)       # Enable normalized device coordinates# Format the text
                best_cost_text = f"Best #chi^{{2}}: {best_overall_cost:.3f}"
                latex.DrawLatex(0.35, 0.85, best_cost_text)                
                c2.Update()
                    
                print("\n")    

            c2.Update()

    c2.Print(outputpdf+'(')
    c3.Print(outputpdf)
    c4.Print(outputpdf)
    c5.Print(outputpdf)
    c6.Print(outputpdf+')')

    print(f"\n\nFits saved to {outputpdf}...")
