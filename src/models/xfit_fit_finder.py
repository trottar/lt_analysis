#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-02-02 22:47:35 trottar"
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
                     calculate_information_criteria, compute_p_value, log_iteration)

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

def parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec,
                 prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params,
                 outputpdf, full_optimization=True):
    # Create lists to store graph objects for each fit
    graphs_sig_fit       = []  # TGraphErrors for data and final fit curve
    graphs_sig_params    = []  # List of lists: one TGraph per parameter (for convergence)
    graphs_sig_chi2      = []  # Reduced chi-square evolution
    graphs_sig_temp      = []  # Temperature evolution
    graphs_sig_accept    = []  # Acceptance probability evolution
    graphs_sig_converge  = []  # Alternate chi2 convergence
    graphs_sig_residuals = []  # Residual evolution
    graphs_sig_ic_aic    = []  # AIC evolution
    graphs_sig_ic_bic    = []  # BIC evolution

    # Create ROOT canvases for plotting
    c2 = TCanvas("c2", "Data & Final Fit", 800, 800)
    c2.Divide(2,2)
    c3 = TCanvas("c3", "Parameter Convergence", 800, 800)
    c3.Divide(2,2)
    c4 = TCanvas("c4", "Red. Chi-Square Convergence", 800, 800)
    c4.Divide(2,2)
    c5 = TCanvas("c5", "Temperature", 800, 800)
    c5.Divide(2,2)
    c6 = TCanvas("c6", "Acceptance Probability", 800, 800)
    c6.Divide(2,2)
    c7 = TCanvas("c7", "Residuals", 800,600)
    c7.Divide(2,2)
    c8 = TCanvas("c8", "Information Criteria", 800,600)
    c8.Divide(2,2)

    # Unpack input dictionary
    q2_set = inpDict["q2_set"]
    w_set = inpDict["w_set"]
    nsep, t_vec, g_vec, w_vec, q2_vec, th_vec = inpDict["objects"]
    max_iterations      = inpDict["max_iterations"]
    num_optimizations   = inpDict["num_optimizations"]
    initial_param_bounds= inpDict["initial_param_bounds"]
    tmin_range          = inpDict["tmin_range"]
    tmax_range          = inpDict["tmax_range"]
    Q2min_range         = inpDict["Q2min_range"]
    Q2max_range         = inpDict["Q2max_range"]
    iter_num            = inpDict["iter_num"]
    fit_params          = inpDict["fit_params"]
    chi2_threshold      = inpDict["chi2_threshold"]
    xfit_log            = inpDict["xfit_log"]
    
    num_events = nsep.GetEntries()

    # Determine central values for guidance
    q2_center_val = get_central_value(q2_vec)
    w_center_val = get_central_value(w_vec)
    g_center_val = get_central_value(g_vec)
    th_center_val = get_central_value(th_vec)

    # Loop over each signal in fit_params
    for it, (sig_name, val) in enumerate(fit_params.items()):
        # Get the number of parameters, initial parameters, and the equation string.
        num_params, initial_params, equation_str = inpDict["initial_params"](sig_name, val)
        # Replace any zero parameters (to avoid perturbation issues)
        initial_params = [v if abs(v) > 0.0 else 1.0 for v in initial_params]
        param_str = ', '.join(str(p) for p in initial_params)

        # Create dynamic parameter convergence graphs (one TGraph per parameter)
        param_graphs = []
        for i in range(num_params):
            param_graphs.append(TGraph())
        graphs_sig_params.append(param_graphs)

        # Create diagnostic graphs for this signal
        graph_sig_fit       = TGraphErrors()
        graph_sig_chi2      = TGraph()
        graph_sig_temp      = TGraph()
        graph_sig_accept    = TGraph()
        graph_sig_converge  = TGraph()
        graph_sig_residuals = TGraph()
        graph_sig_aic       = TGraph()
        graph_sig_bic       = TGraph()

        graphs_sig_fit.append(graph_sig_fit)
        graphs_sig_chi2.append(graph_sig_chi2)
        graphs_sig_temp.append(graph_sig_temp)
        graphs_sig_accept.append(graph_sig_accept)
        graphs_sig_converge.append(graph_sig_converge)
        graphs_sig_residuals.append(graph_sig_residuals)
        graphs_sig_ic_aic.append(graph_sig_aic)
        graphs_sig_ic_bic.append(graph_sig_bic)

        # Branch based on whether parameters are to be re-fitted or fixed.
        if sig_name not in fixed_params:
            # --------------------- Optimization Branch ---------------------
            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Parameters: ({param_str})")
            print(f"{equation_str}")
            print("/*--------------------------------------------------*/")
            
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Initialize adaptive regularization and solution tracking.
            lambda_reg = 0.01
            best_overall_cost = float('inf')
            best_overall_params = current_params = initial_params[:]  # working copy
            best_overall_bin = None
            total_iteration = 0
            cost_history = []
            local_search_interval = 25
            current_errors = [0.0] * num_params

            # For demonstration, loop over a chosen bin (e.g. bin index 2)
            for b in [2]:
                print(f"\nOptimizing Sig {sig_name} for bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                iteration = 0
                initial_temperature = 1.0
                temperature = initial_temperature

                while iteration <= max_iterations:
                    # Perturb parameters using simulated annealing
                    new_params = [simulated_annealing(p, temperature) for p in current_params]
                    current_params = new_params

                    # Build data graph from nsep
                    g_sig = TGraphErrors()
                    for i in range(nsep.GetSelectedRows()):
                        g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                        g_sig.SetPointError(i, 0, nsep.GetV3()[i])

                    # Select model function based on sig_name
                    if sig_name == "L":
                        fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                    elif sig_name == "T":
                        fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                    elif sig_name == "LT":
                        fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                    elif sig_name == "TT":
                        fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])

                    # Create the fit function with dynamic parameter number
                    f_sig = TF1(f"sig_{sig_name}", fun_Sig, tmin_range, tmax_range, num_params)
                    f_sig.SetParNames(*[f"p{i}" for i in range(num_params)])
                    for i in range(num_params):
                        f_sig.SetParameter(i, current_params[i])
                        if full_optimization:
                            f_sig.SetParLimits(i, -initial_param_bounds, initial_param_bounds)
                        else:
                            offset = 0.1 * abs(current_params[i])
                            f_sig.SetParLimits(i, current_params[i]-offset, current_params[i]+offset)

                    # Fit the function to the data graph in quiet mode ("SQ")
                    fit_result = graph_sig_fit.Fit(f_sig, "SQ")

                    # Update parameter convergence graphs.
                    for i in range(num_params):
                        graphs_sig_params[it][i].SetPoint(total_iteration, total_iteration, current_params[i])

                    # Compute cost with adaptive regularization.
                    current_cost, lambda_reg = calculate_cost(f_sig, g_sig, current_params, num_events, num_params, lambda_reg)
                    cost_history.append(current_cost)
                    accept_prob = acceptance_probability(best_overall_cost, current_cost, temperature)

                    # ------------------ Diagnostic Logging ------------------
                    log_iteration(total_iteration, current_cost, current_params, accept_prob, temperature, xfit_log)
                    # ----------------------------------------------------------

                    # Update diagnostic graphs.
                    graph_sig_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                    graph_sig_temp.SetPoint(total_iteration, total_iteration, round(temperature, 4))
                    graph_sig_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                    # Compute average residual for this iteration.
                    residuals = []
                    for i in range(g_sig.GetN()):
                        x = g_sig.GetX()[i]
                        y_data = g_sig.GetY()[i]
                        y_err  = g_sig.GetEY()[i]
                        y_fit  = f_sig.Eval(x)
                        residuals.append(abs((y_data - y_fit) / (y_err if y_err != 0 else 1.0)))
                    avg_residual = np.mean(residuals) if residuals else 0.0
                    graph_sig_residuals.SetPoint(total_iteration, total_iteration, round(avg_residual, 4))

                    # Update information criteria graphs.
                    n_samples = len(w_vec)
                    ic_aic, ic_bic = calculate_information_criteria(n_samples, num_params, current_cost)
                    graph_sig_aic.SetPoint(total_iteration, total_iteration, round(ic_aic, 4))
                    graph_sig_bic.SetPoint(total_iteration, total_iteration, round(ic_bic, 4))

                    # Accept new state based on probability.
                    if accept_prob > random.random():
                        best_overall_params = [f_sig.GetParameter(i) for i in range(num_params)]
                        best_cost = current_cost
                        best_overall_bin = b
                        best_overall_errors = [f_sig.GetParError(i) for i in range(num_params)]

                    # Periodically perform a local search to refine parameters.
                    if iteration % local_search_interval == 0:
                        current_params = local_search(current_params, f_sig, num_params)

                    previous_params = current_params[:]
                    temperature = adaptive_cooling(initial_temperature, iteration)
                    iteration += 1
                    total_iteration += 1

                # End bin loop.
            end_time = time.time()
            total_duration = end_time - start_time
            print(f"\nBest overall cost for Sig {sig_name}: {best_overall_cost:.3f}")
            print(f"Best overall parameters: {best_overall_params}")
            print(f"Time taken: {total_duration:.2f} seconds")

            # Update output parameter vectors dynamically.
            for j in range(num_params):
                par_vec[num_params * it + j] = best_overall_params[j]
                par_err_vec[num_params * it + j] = best_overall_errors[j]
                par_chi2_vec[num_params * it + j] = best_overall_cost

            # Prepare final data graph for the best bin.
            g_sig_final = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig_final.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig_final.SetPointError(i, 0, nsep.GetV3()[i])

            # Rebuild the final fit function using best_overall_bin.
            if sig_name == "L":
                fun_Sig_final = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin],
                                                  w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "T":
                fun_Sig_final = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin],
                                                  w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "LT":
                fun_Sig_final = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin],
                                                   w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "TT":
                fun_Sig_final = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin],
                                                   w_vec[best_overall_bin], th_vec[best_overall_bin])
            f_sig_final = TF1(f"sig_{sig_name}_final", fun_Sig_final, tmin_range, tmax_range, num_params)
            f_sig_final.SetParNames(*[f"p{i}" for i in range(num_params)])
            for i in range(num_params):
                f_sig_final.FixParameter(i, best_overall_params[i])

            # Draw the final fit on canvas c2.
            c2.cd(it+1)
            g_sig_final.Draw("AP")
            f_sig_final.Draw("same")
            graph_sig_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graph_sig_fit.GetXaxis().CenterTitle()
            graph_sig_fit.GetYaxis().SetTitle(f"#left(#frac{{#it{{d#sigma}}}}{{#it{{dt}}}}#right)_{{{sig_name}}} [nb/GeV^2]")
            graph_sig_fit.GetYaxis().SetTitleOffset(1.5)
            graph_sig_fit.GetYaxis().SetTitleSize(0.035)
            graph_sig_fit.GetYaxis().CenterTitle()
            x_min = min(graph_sig_fit.GetX())
            x_max = max(graph_sig_fit.GetX())
            y_min = min(graph_sig_fit.GetY())
            y_max = max(graph_sig_fit.GetY())
            margin = 0.1
            graph_sig_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graph_sig_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            n_points = 100
            fit_y_values = [f_sig_final.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)
            margin = 0.1 * (y_max - y_min)
            graph_sig_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            r_sig_fit = graph_sig_fit.Fit(f_sig_final, "SQ")
            f_sig_final.Draw("same")

            # Compute and display chi-square p-value for the optimization branch.
            latex = TLatex()
            latex.SetTextSize(0.04)
            latex.SetNDC(True)
            ndf = f_sig_final.GetNDF()
            chi2_val = f_sig_final.GetChisquare()
            p_value = compute_p_value(chi2_val, ndf)
            print(f"Chi-square p-value (optimization branch): {p_value:.4f}")
            latex.DrawLatex(0.35, 0.85, f"p-value: {p_value:.4f}")
            c2.Update()

            # ------------------ Plotting Diagnostics for Optimization Branch ------------------
            # Parameter Convergence on Canvas c3
            c3.cd(it+1)
            for i in range(num_params):
                param_graph = graphs_sig_params[it][i]
                if i == 0:
                    param_graph.SetLineColor(kRed)
                elif i == 1:
                    param_graph.SetLineColor(kBlue)
                elif i == 2:
                    param_graph.SetLineColor(kGreen)
                elif i == 3:
                    param_graph.SetLineColor(kMagenta)
                if i == 0:
                    param_graph.Draw("ALP")
                else:
                    param_graph.Draw("LP SAME")
            c3.Update()

            # Reduced Chi-square Convergence on Canvas c4
            c4.cd(it+1)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} Chi-square Convergence")
            graphs_sig_converge[it].SetLineColor(kBlack)
            graphs_sig_converge[it].Draw("ALP")
            c4.Update()

            # Temperature Evolution on Canvas c5
            c5.cd(it+1)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Evolution")
            graphs_sig_temp[it].SetLineColor(kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()

            # Acceptance Probability on Canvas c6
            c6.cd(it+1)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()

            # Residual Evolution on Canvas c7
            c7.cd(it+1)
            graphs_sig_residuals[it].SetTitle(f"Sig {sig_name} Residual Evolution")
            graphs_sig_residuals[it].SetLineColor(kBlack)
            graphs_sig_residuals[it].Draw("ALP")
            c7.Update()

            # Information Criteria on Canvas c8
            c8.cd(it+1)
            graphs_sig_ic_aic[it].SetTitle(f"Sig {sig_name} Information Criteria (AIC)")
            graphs_sig_ic_aic[it].SetLineColor(kRed)
            graphs_sig_ic_aic[it].Draw("ALP")
            graphs_sig_ic_bic[it].SetTitle(f"Sig {sig_name} Information Criteria (BIC)")
            graphs_sig_ic_bic[it].SetLineColor(kGreen)
            graphs_sig_ic_bic[it].Draw("LP SAME")
            c8.Update()
            # ---------------------------------------------------------------------
            print("\n")
        else:
            # --------------------- Fixed Parameter Branch ---------------------
            print("\n/*--------------------------------------------------*/")
            print(f"Fixed Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Parameters: ({param_str})")
            print(f"{equation_str}")
            print("/*--------------------------------------------------*/")
            
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            
            # Do NOT re-create and append diagnostic graphs here—use those created above.
            lambda_reg = 0.01
            best_overall_cost = float('inf')
            best_overall_params = []
            best_overall_bin = 0
            
            for b in [2]:
                print(f"\nDetermining best fit for fixed Sig {sig_name} using bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                g_sig_fit = TGraphErrors()
                graphs_sig_fit.append(g_sig_fit)
                
                # Build data graph from nsep
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
                    fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "T":
                    fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "LT":
                    fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "TT":
                    fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                
                f_sig = TF1(f"sig_{sig_name}", fun_Sig, tmin_range, tmax_range, num_params)
                f_sig.SetParNames(*[f"p{i}" for i in range(num_params)])
                for i in range(num_params):
                    f_sig.FixParameter(i, par_vec[num_params * it + i])
                
                r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
                
                current_cost, lambda_reg = calculate_cost(f_sig, g_sig, par_vec[num_params * it : num_params * (it+1)],
                                                          num_events, num_params, lambda_reg)
                print(f"\tCost: {current_cost:.3f}")
                if abs(current_cost - 1) < abs(best_overall_cost - 1):
                    best_overall_cost = current_cost
                    best_overall_bin = b
                    best_overall_params = [par_vec[num_params * it + j] for j in range(num_params)]
            
            print(f"\nBest overall solution (fixed): {best_overall_params}")
            print(f"Best overall cost (fixed): {best_overall_cost:.5f}")
            for j in range(num_params):
                par_chi2_vec[num_params * it + j] = best_overall_cost

            c2.cd(it+1)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit (Fixed Parameters)")
            graphs_sig_fit[it].Draw("A*")
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle(f"#left(#frac{{#it{{d#sigma}}}}{{#it{{dt}}}}#right)_{{{sig_name}}} [nb/GeV^2]")
            graphs_sig_fit[it].GetYaxis().SetTitleOffset(1.5)
            graphs_sig_fit[it].GetYaxis().SetTitleSize(0.035)
            graphs_sig_fit[it].GetYaxis().CenterTitle()
            x_min = min(graphs_sig_fit[it].GetX())
            x_max = max(graphs_sig_fit[it].GetX())
            y_min = min(graphs_sig_fit[it].GetY())
            y_max = max(graphs_sig_fit[it].GetY())
            margin = 0.1
            graphs_sig_fit[it].GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            n_points = 100
            fit_y_values = [f_sig.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            r_sig_fit = graphs_sig_fit[it].Fit(f_sig, "SQ")
            f_sig.Draw("same")
            latex = TLatex()
            latex.SetTextSize(0.04)
            latex.SetNDC(True)
            ndf = f_sig.GetNDF()
            chi2_val = f_sig.GetChisquare()
            p_value = compute_p_value(chi2_val, ndf)
            print(f"Chi-square p-value (fixed branch): {p_value:.4f}")
            latex.DrawLatex(0.35, 0.85, f"p-value: {p_value:.4f}")
            c2.Update()

            # ------------------ Plotting Diagnostics for Fixed Branch ------------------
            # Parameter Convergence on Canvas c3
            c3.cd(it+1)
            for i in range(num_params):
                param_graph = graphs_sig_params[it][i]
                if i == 0:
                    param_graph.SetLineColor(kRed)
                elif i == 1:
                    param_graph.SetLineColor(kBlue)
                elif i == 2:
                    param_graph.SetLineColor(kGreen)
                elif i == 3:
                    param_graph.SetLineColor(kMagenta)
                if i == 0:
                    param_graph.Draw("ALP")
                else:
                    param_graph.Draw("LP SAME")
            c3.Update()

            # Reduced Chi-square Convergence on Canvas c4
            c4.cd(it+1)
            graphs_sig_converge[it].SetTitle(f"Sig {sig_name} Chi-square Convergence")
            graphs_sig_converge[it].SetLineColor(kBlack)
            graphs_sig_converge[it].Draw("ALP")
            c4.Update()

            # Temperature Evolution on Canvas c5
            c5.cd(it+1)
            graphs_sig_temp[it].SetTitle(f"Sig {sig_name} Temperature Evolution")
            graphs_sig_temp[it].SetLineColor(kBlack)
            graphs_sig_temp[it].Draw("ALP")
            c5.Update()

            # Acceptance Probability on Canvas c6
            c6.cd(it+1)
            graphs_sig_accept[it].SetTitle(f"Sig {sig_name} Acceptance Probability")
            graphs_sig_accept[it].SetLineColor(kBlack)
            graphs_sig_accept[it].Draw("ALP")
            c6.Update()

            # Residual Evolution on Canvas c7
            c7.cd(it+1)
            graphs_sig_residuals[it].SetTitle(f"Sig {sig_name} Residual Evolution")
            graphs_sig_residuals[it].SetLineColor(kBlack)
            graphs_sig_residuals[it].Draw("ALP")
            c7.Update()

            # Information Criteria on Canvas c8
            c8.cd(it+1)
            graphs_sig_ic_aic[it].SetTitle(f"Sig {sig_name} Information Criteria (AIC)")
            graphs_sig_ic_aic[it].SetLineColor(kRed)
            graphs_sig_ic_aic[it].Draw("ALP")
            graphs_sig_ic_bic[it].SetTitle(f"Sig {sig_name} Information Criteria (BIC)")
            graphs_sig_ic_bic[it].SetLineColor(kGreen)
            graphs_sig_ic_bic[it].Draw("LP SAME")
            c8.Update()
            # ---------------------------------------------------------------------
            print("\n")
        # End of each signal's branch.
    # Print canvases to the output PDF file.
    c2.Print(outputpdf + '(')
    c3.Print(outputpdf)
    c4.Print(outputpdf)
    c5.Print(outputpdf)
    c6.Print(outputpdf)
    c7.Print(outputpdf)
    c8.Print(outputpdf + ')')

    print(f"\n\nFits saved to {outputpdf}...")
