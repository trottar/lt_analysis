#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-13 14:28:44 trottar"
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
from utility import adaptive_cooling, simulated_annealing, acceptance_probability, adjust_params, local_search, select_valid_parameter

##################################################################################################################################################

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

###############################################################################################################################################
# Import separated xsects models

from xfit_active import fun_Sig_L, fun_Sig_T, fun_Sig_LT, fun_Sig_TT

##################################################################################################################################################

def find_fit(sig_fit_dict, inp_dict, par_vec, par_err_vec, par_chi2_vec):

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
    
    nsep, g_vec, w_vec, q2_vec, th_vec = inp_dict["objects"]
    max_iterations = inp_dict["max_iterations"]
    num_optimizations = inp_dict["num_optimizations"]
    tmin_range = inp_dict["tmin_range"]
    tmax_range = inp_dict["tmax_range"]
    Q2min_range = inp_dict["Q2min_range"]
    Q2max_range = inp_dict["Q2max_range"]
    iter_num = inp_dict["iter_num"]            
    outputpdf = inp_dict["outputpdf"]

    num_events = nsep.GetEntries()    
    
    for it, (key, val) in enumerate(sig_fit_dict.items()):

        sig_name = key
        # Grab parameters used by functional forms
        num_params, initial_params = inp_dict["initial_params"](sig_name, val["params"])

        #
        initial_params = [v for v in initial_params if v > 0.0 else 1e-5]
        #prev_params = [v for v in val["params"] if v != 0.0]
        #if num_params > len(prev_params):            
        #if num_params < len(prev_params):
            

        if num_events <= num_params:
            print(f"\n\nWARNING: The number of parameters ({num_params}) for Sig {sig_name} is greater than or equal to the number of data points ({num_events})! Using adaptive regression methods for determining quality of fit...")
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

                    g_sig_prv = TGraph()
                    g_sig_fit = TGraphErrors()
                    g_sig_fit_tot = TGraph()

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
                            sig_X_fit = g_sig.GetY()[i]
                            sig_X_fit_err = g_sig.GetEY()[i]

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
                        #f_sig.SetParLimits(0, current_params - abs(current_params * par_sig_0), current_params + abs(current_params * par_sig_0))
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(0, -5, 5)
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)

                        g_q2_sig_fit = TGraphErrors()
                        for i in range(len(w_vec)):
                            g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                            sig_X = (f_sig.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                            g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                            residuals = []
                            for i in range(num_events):
                                observed = g_sig.GetY()[i]
                                expected = f_sig.Eval(g_sig.GetX()[i])
                                residual = (observed - expected) / g_sig.GetEY()[i] if g_sig.GetEY()[i] != 0 else (observed - expected)
                                residuals.append(residual)
                            # Mean Squared Error (MSE)
                            mse = np.mean(np.square(residuals))
                            # L2 regularization term
                            l2_reg = sum(p**2 for p in current_params)
                            # Regularized cost function
                            current_cost = mse + lambda_reg * l2_reg
                            # Effective degrees of freedom
                            effective_dof = max(num_events - num_params, 1)
                            # Adjusted cost (similar to reduced chi-squared)
                            adjusted_cost = current_cost / effective_dof
                            # Store cost history
                            cost_history.append(adjusted_cost)
                            # Adaptive regularization
                            if len(cost_history) > 1:
                                if adjusted_cost < cost_history[-2]:
                                    # If cost is decreasing, slightly decrease regularization
                                    lambda_reg = max(lambda_reg * lambda_decrease, lambda_min)
                                else:
                                    # If cost is increasing or stagnant, increase regularization
                                    lambda_reg = min(lambda_reg * lambda_increase, lambda_max)
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, adjusted_cost, temperature)
                            
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

            while len(best_overall_params) < 4:
                best_overall_params.append(0.0)
                best_overall_errors.append(0.0)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_prv = TGraph()
            g_sig_fit = TGraphErrors()
            g_sig_fit_tot = TGraph()

            graphs_sig_fit.append(graphs_sig_fit[it])            

            if sig_name == "L":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
            f_sig_pre.SetParNames("p0")
            f_sig_pre.FixParameter(0, best_overall_params[0])

            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                sig_X_pre = (f_sig_pre.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                g_sig_prv.SetPoint(i, g_sig.GetX()[i], sig_X_pre)

                sig_X_fit = g_sig.GetY()[i]
                sig_X_fit_err = g_sig.GetEY()[i]

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            g_sig.SetTitle(f"Sig {sig_name}")
            g_sig.SetMarkerStyle(5)
            g_sig.Draw("AP")
            g_sig.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            g_sig.GetXaxis().CenterTitle()
            g_sig.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            g_sig.GetYaxis().SetTitleOffset(1.5)
            g_sig.GetYaxis().SetTitleSize(0.035)
            g_sig.GetYaxis().CenterTitle()

            g_sig_prv.SetMarkerColor(4)
            g_sig_prv.SetMarkerStyle(25)
            g_sig_prv.Draw("P")

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            c2.Update()
            
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

            g_q2_sig_fit = TGraphErrors()
            for i in range(len(w_vec)):
                g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                sig_X = (f_sig.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
            
            # Plot temperature convergence
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability convergence
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

                    g_sig_prv = TGraph()
                    g_sig_fit = TGraphErrors()
                    g_sig_fit_tot = TGraph()
                    
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
                            sig_X_fit = g_sig.GetY()[i]
                            sig_X_fit_err = g_sig.GetEY()[i]

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
                        #f_sig.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sig_0), current_params[0] + abs(current_params[0] * par_sig_0))
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sig_1), current_params[1] + abs(current_params[1] * par_sig_1))
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)

                        g_q2_sig_fit = TGraphErrors()
                        for i in range(len(w_vec)):
                            g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                            sig_X = (f_sig.Eval(g_sig.GetX()[i])) * (g_vec[i])
                            g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                            residuals = []
                            for i in range(num_events):
                                observed = g_sig.GetY()[i]
                                expected = f_sig.Eval(g_sig.GetX()[i])
                                residual = (observed - expected) / g_sig.GetEY()[i] if g_sig.GetEY()[i] != 0 else (observed - expected)
                                residuals.append(residual)
                            # Mean Squared Error (MSE)
                            mse = np.mean(np.square(residuals))
                            # L2 regularization term
                            l2_reg = sum(p**2 for p in current_params)
                            # Regularized cost function
                            current_cost = mse + lambda_reg * l2_reg
                            # Effective degrees of freedom
                            effective_dof = max(num_events - num_params, 1)
                            # Adjusted cost (similar to reduced chi-squared)
                            adjusted_cost = current_cost / effective_dof
                            # Store cost history
                            cost_history.append(adjusted_cost)
                            # Adaptive regularization
                            if len(cost_history) > 1:
                                if adjusted_cost < cost_history[-2]:
                                    # If cost is decreasing, slightly decrease regularization
                                    lambda_reg = max(lambda_reg * lambda_decrease, lambda_min)
                                else:
                                    # If cost is increasing or stagnant, increase regularization
                                    lambda_reg = min(lambda_reg * lambda_increase, lambda_max)
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, adjusted_cost, temperature)
                            
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

                        graphs_sig_temp[it].SetPoint(total_iteration, total_iteration, temperature)

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
                    best_overall_errors = best_errors

            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            while len(best_overall_params) < 4:
                best_overall_params.append(0.0)
                best_overall_errors.append(0.0)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_prv = TGraph()
            g_sig_fit = TGraphErrors()
            g_sig_fit_tot = TGraph()

            graphs_sig_fit.append(g_sig_fit)            

            if sig_name == "L":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)    
            f_sig_pre.SetParNames("p0", "p1")
            f_sig_pre.FixParameter(0, best_overall_params[0])
            f_sig_pre.FixParameter(1, best_overall_params[1])

            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                sig_X_pre = (f_sig_pre.Eval(g_sig.GetX()[i])) * (g_vec[i])
                g_sig_prv.SetPoint(i, g_sig.GetX()[i], sig_X_pre)

                sig_X_fit = g_sig.GetY()[i]
                sig_X_fit_err = g_sig.GetEY()[i]

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            g_sig.SetTitle(f"Sig {sig_name}")
            g_sig.SetMarkerStyle(5)
            g_sig.Draw("AP")
            g_sig.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            g_sig.GetXaxis().CenterTitle()
            g_sig.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            g_sig.GetYaxis().SetTitleOffset(1.5)
            g_sig.GetYaxis().SetTitleSize(0.035)
            g_sig.GetYaxis().CenterTitle()

            g_sig_prv.SetMarkerColor(4)
            g_sig_prv.SetMarkerStyle(25)
            g_sig_prv.Draw("P")

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            c2.Update()
            
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

            g_q2_sig_fit = TGraphErrors()
            for i in range(len(w_vec)):
                g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                sig_X = (f_sig.Eval(g_sig.GetX()[i])) * (g_vec[i])
                g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                par_sig_1 = initial_params[0]
                par_sig_2 = initial_params[0]
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

                    g_sig_prv = TGraph()
                    g_sig_fit = TGraphErrors()
                    g_sig_fit_tot = TGraph()    

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
                            sig_X_fit = g_sig.GetY()[i]
                            sig_X_fit_err = g_sig.GetEY()[i]

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
                        #f_sig.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sig_0), current_params[0] + abs(current_params[0] * par_sig_0))
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sig_1), current_params[1] + abs(current_params[1] * par_sig_1))
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(2, current_params[2] - abs(current_params[2] * par_sig_2), current_params[2] + abs(current_params[2] * par_sig_2))
                        f_sig.SetParLimits(2, -max_param_value, max_param_value)

                        g_q2_sig_fit = TGraphErrors()
                        for i in range(len(w_vec)):
                            g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                            sig_X = (f_sig.Eval(g_sig.GetX()[i])) * (g_vec[i])
                            g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                            residuals = []
                            for i in range(num_events):
                                observed = g_sig.GetY()[i]
                                expected = f_sig.Eval(g_sig.GetX()[i])
                                residual = (observed - expected) / g_sig.GetEY()[i] if g_sig.GetEY()[i] != 0 else (observed - expected)
                                residuals.append(residual)
                            # Mean Squared Error (MSE)
                            mse = np.mean(np.square(residuals))
                            # L2 regularization term
                            l2_reg = sum(p**2 for p in current_params)
                            # Regularized cost function
                            current_cost = mse + lambda_reg * l2_reg
                            # Effective degrees of freedom
                            effective_dof = max(num_events - num_params, 1)
                            # Adjusted cost (similar to reduced chi-squared)
                            adjusted_cost = current_cost / effective_dof
                            # Store cost history
                            cost_history.append(adjusted_cost)
                            # Adaptive regularization
                            if len(cost_history) > 1:
                                if adjusted_cost < cost_history[-2]:
                                    # If cost is decreasing, slightly decrease regularization
                                    lambda_reg = max(lambda_reg * lambda_decrease, lambda_min)
                                else:
                                    # If cost is increasing or stagnant, increase regularization
                                    lambda_reg = min(lambda_reg * lambda_increase, lambda_max)
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, adjusted_cost, temperature)
                            
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
                    best_overall_errors = best_errors

            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            while len(best_overall_params) < 4:
                best_overall_params.append(0.0)
                best_overall_errors.append(0.0)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_prv = TGraph()
            g_sig_fit = TGraphErrors()
            g_sig_fit_tot = TGraph()        

            graphs_sig_fit.append(g_sig_fit)
            
            if sig_name == "L":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)
            f_sig_pre.SetParNames("p0", "p1", "p2")
            f_sig_pre.FixParameter(0, best_overall_params[0])
            f_sig_pre.FixParameter(1, best_overall_params[1])
            f_sig_pre.FixParameter(2, best_overall_params[2])

            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                sig_X_pre = (f_sig_pre.Eval(g_sig.GetX()[i])) * (g_vec[i])
                g_sig_prv.SetPoint(i, g_sig.GetX()[i], sig_X_pre)

                sig_X_fit = g_sig.GetY()[i]
                sig_X_fit_err = g_sig.GetEY()[i]

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            g_sig.SetTitle(f"Sig {sig_name}")
            g_sig.SetMarkerStyle(5)
            g_sig.Draw("AP")
            g_sig.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            g_sig.GetXaxis().CenterTitle()
            g_sig.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            g_sig.GetYaxis().SetTitleOffset(1.5)
            g_sig.GetYaxis().SetTitleSize(0.035)
            g_sig.GetYaxis().CenterTitle()

            g_sig_prv.SetMarkerColor(4)
            g_sig_prv.SetMarkerStyle(25)
            g_sig_prv.Draw("P")

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            c2.Update()
            
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

            g_q2_sig_fit = TGraphErrors()
            for i in range(len(w_vec)):
                g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                sig_X = (f_sig.Eval(g_sig.GetX()[i])) * (g_vec[i])
                g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                previous_params = best_params[:]
                best_errors = [par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3]

                # Check for local minima
                local_minima = []
                local_iterations = 0
                tabu_list = set()

                # Local search
                local_search_interval = 25

                while iteration <= max_iterations:

                    g_sig_prv = TGraph()
                    g_sig_fit = TGraphErrors()
                    g_sig_fit_tot = TGraph()    

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
                            sig_X_fit = g_sig.GetY()[i]
                            sig_X_fit_err = g_sig.GetEY()[i]

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
                        #f_sig.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sig_0), current_params[0] + abs(current_params[0] * par_sig_0))
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sig_1), current_params[1] + abs(current_params[1] * par_sig_1))
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(2, current_params[2] - abs(current_params[2] * par_sig_2), current_params[2] + abs(current_params[2] * par_sig_2))
                        f_sig.SetParLimits(2, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(3, current_params[3] - abs(current_params[3] * par_sig_3), current_params[3] + abs(current_params[3] * par_sig_3))
                        f_sig.SetParLimits(3, -max_param_value, max_param_value)                
                        #f_sig.SetParLimits(0, -5, 5)
                        f_sig.SetParLimits(0, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(1, -5, 5)
                        f_sig.SetParLimits(1, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(2, -5, 5)
                        f_sig.SetParLimits(2, -max_param_value, max_param_value)
                        #f_sig.SetParLimits(3, -5, 5)
                        f_sig.SetParLimits(3, -max_param_value, max_param_value)                

                        g_q2_sig_fit = TGraphErrors()
                        for i in range(len(w_vec)):
                            g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                            g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                            sig_X = (f_sig.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                            g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
                            residuals = []
                            for i in range(num_events):
                                observed = g_sig.GetY()[i]
                                expected = f_sig.Eval(g_sig.GetX()[i])
                                residual = (observed - expected) / g_sig.GetEY()[i] if g_sig.GetEY()[i] != 0 else (observed - expected)
                                residuals.append(residual)
                            # Mean Squared Error (MSE)
                            mse = np.mean(np.square(residuals))
                            # L2 regularization term
                            l2_reg = sum(p**2 for p in current_params)
                            # Regularized cost function
                            current_cost = mse + lambda_reg * l2_reg
                            # Effective degrees of freedom
                            effective_dof = max(num_events - num_params, 1)
                            # Adjusted cost (similar to reduced chi-squared)
                            adjusted_cost = current_cost / effective_dof
                            # Store cost history
                            cost_history.append(adjusted_cost)
                            # Adaptive regularization
                            if len(cost_history) > 1:
                                if adjusted_cost < cost_history[-2]:
                                    # If cost is decreasing, slightly decrease regularization
                                    lambda_reg = max(lambda_reg * lambda_decrease, lambda_min)
                                else:
                                    # If cost is increasing or stagnant, increase regularization
                                    lambda_reg = min(lambda_reg * lambda_increase, lambda_max)
                            # Acceptance probability
                            accept_prob = acceptance_probability(best_cost, adjusted_cost, temperature)
                            
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

                        # If the new cost is better or accepted by the acceptance probability, update the best parameters
                        if accept_prob > random.random():
                            best_params = current_params
                            best_cost = current_cost

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

                        par_sig_0, par_sig_1, par_sig_2, par_sig_3 = initial_params
                        par_sig_err_0, par_sig_err_1, par_sig_err_2, par_sig_err_3 = [0.0 for _ in range(num_params)]

                        iteration += 1
                        total_iteration += 1 if iteration % max_iterations == 0 else 0                

                # After the while loop, check if this run found a better solution
                if abs(best_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = best_cost
                    best_overall_params = best_params[:]
                    best_overall_errors = best_errors

            print("\nBest overall solution: {0}".format(best_overall_params))
            print("Best overall cost: {0}".format(best_overall_cost))

            # Record the end time
            end_time = time.time()
            # Calculate the total duration
            total_duration = end_time - start_time
            print("The loop took {:.2f} seconds.".format(total_duration))

            while len(best_overall_params) < 4:
                best_overall_params.append(0.0)
                best_overall_errors.append(0.0)

            par_vec.append(best_overall_params[0])
            par_vec.append(best_overall_params[1])
            par_vec.append(best_overall_params[2])
            par_vec.append(best_overall_params[3])

            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])
            par_err_vec.append(best_overall_errors[0])

            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)
            par_chi2_vec.append(best_overall_cost)

            g_sig_prv = TGraph()
            g_sig_fit = TGraphErrors()
            g_sig_fit_tot = TGraph()    

            graphs_sig_fit.append(g_sig_fit)
            
            if sig_name == "L":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_L, 0.0, 2.0, num_params)
            elif sig_name == "T":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_T, 0.0, 2.0, num_params)
            elif sig_name == "LT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_LT, 0.0, 2.0, num_params)
            elif sig_name == "TT":
                #f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, tmin_range, tmax_range, num_params)
                f_sig_pre = TF1(f"sig_{sig_name}", fun_Sig_TT, 0.0, 2.0, num_params)    
            f_sig_pre.SetParNames("p0", "p1", "p2", "p3")
            f_sig_pre.FixParameter(0, best_overall_params[0])
            f_sig_pre.FixParameter(1, best_overall_params[1])
            f_sig_pre.FixParameter(2, best_overall_params[2])
            f_sig_pre.FixParameter(3, best_overall_params[3])

            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])

            for i in range(len(w_vec)):
                sig_X_pre = (f_sig_pre.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                g_sig_prv.SetPoint(i, g_sig.GetX()[i], sig_X_pre)

                sig_X_fit = g_sig.GetY()[i]
                sig_X_fit_err = g_sig.GetEY()[i]

                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)

            g_sig.SetTitle(f"Sig {sig_name}")
            g_sig.SetMarkerStyle(5)
            g_sig.Draw("AP")
            g_sig.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            g_sig.GetXaxis().CenterTitle()
            g_sig.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{%s} [nb/GeV^{2}]" % sig_name)
            g_sig.GetYaxis().SetTitleOffset(1.5)
            g_sig.GetYaxis().SetTitleSize(0.035)
            g_sig.GetYaxis().CenterTitle()

            g_sig_prv.SetMarkerColor(4)
            g_sig_prv.SetMarkerStyle(25)
            g_sig_prv.Draw("P")

            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            c2.Update()
            
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

            g_q2_sig_fit = TGraphErrors()
            for i in range(len(w_vec)):
                g_q2_sig_fit.SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                g_q2_sig_fit.SetPointError(i, 0.0, sig_X_fit_err)
                sig_X = (f_sig.Eval(g_sig.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                g_sig_fit_tot.SetPoint(i, g_sig.GetX()[i], sig_X)

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
            
            # Plot temperature convergence
            c5.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graphs_sig_temp[it].SetLineColor(ROOT.kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()
            
            # Plot acceptance probability convergence
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
