#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-02-04 06:30:22 trottar"
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
from utility import (
    adaptive_regularization, calculate_cost, adaptive_cooling,
    simulated_annealing, acceptance_probability, adjust_params, 
    local_search, select_valid_parameter, get_central_value, 
    calculate_information_criteria, sanitize_params
)

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

def parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization=True, debug=False):
    """
    'parameterize' function including:
      - Simple ±initial_param_bounds random init (no expansions).
      - Simulated annealing step each iteration (on all parameters).
      - Adaptive cooling of temperature each iteration.
      - Local search every 15 iterations.
      - Stall check => re-random if no improvement in cost for many accepted attempts.
      - Multiple runs => mitigate local minima.
      - Unconditional best cost updates => no infinite.
      - TGraph updates each iteration => full diagnostic.
      - Debug prints (if debug=True) every 100 iterations => param/cost/temperature/stall info.
      - Data lines from nsep remain in same order, no reordering.
    """

    # -----------------------------------------------------------------------------
    # 1. Create lists to store global graphs (unchanged structure)
    # -----------------------------------------------------------------------------
    graphs_sig_fit      = []
    graphs_sig_params_all = []
    graphs_sig_converge = []
    graphs_sig_temp     = []
    graphs_sig_accept   = []
    graphs_sig_residuals= []
    graphs_sig_ic_aic   = []
    graphs_sig_ic_bic   = []    
    funcs_sig = []
    fits_sig = []
    
    c2 = TCanvas("c2", "c2", 800, 800)
    c2.Divide(2, 2)
    c3 = TCanvas("c3", "Parameter Convergence", 800, 800)
    c3.Divide(2, 2)
    c4 = TCanvas("c4", "Red. Chi-Square Convergence", 800, 800)
    c4.Divide(2, 2)
    c5 = TCanvas("c5", "Temperature", 800, 800)
    c5.Divide(2, 2)
    c6 = TCanvas("c6", "Acceptance Probability", 800, 600)
    c6.Divide(2, 2)
    c7 = TCanvas("c7", "Residuals", 800, 600)
    c7.Divide(2, 2)
    c8 = TCanvas("c8", "Information Criteria", 800, 600)
    c8.Divide(2, 2)

    # -----------------------------------------------------------------------------
    # 2. Unpack input objects/settings
    # -----------------------------------------------------------------------------
    q2_set, w_set = inpDict["q2_set"], inpDict["w_set"]
    nsep, t_vec, g_vec, w_vec, q2_vec, th_vec = inpDict["objects"]
    max_iterations     = inpDict["max_iterations"]
    num_optimizations  = inpDict["num_optimizations"]
    initial_param_bounds = inpDict["initial_param_bounds"]
    tmin_range, tmax_range = inpDict["tmin_range"], inpDict["tmax_range"]
    Q2min_range, Q2max_range = inpDict["Q2min_range"], inpDict["Q2max_range"]
    iter_num    = inpDict["iter_num"]
    fit_params  = inpDict["fit_params"]
    chi2_threshold = inpDict["chi2_threshold"]

    q2_center_val = get_central_value(q2_vec)
    w_center_val  = get_central_value(w_vec)
    g_center_val  = get_central_value(g_vec)
    th_center_val = get_central_value(th_vec)

    num_events = nsep.GetEntries()
    colors = [kRed, kBlue, kGreen, kMagenta]

    # -----------------------------------------------------------------------------
    # 3. Main loop over each fit in fit_params
    # -----------------------------------------------------------------------------
    for it, (sig_name, val) in enumerate(fit_params.items()):
        if sig_name not in fixed_params:

            num_params, init_params, equation_str = inpDict["initial_params"](sig_name, val)
            # Avoid zero init param
            init_params = [v if abs(v) > 0.0 else 1.0 for v in init_params]
            param_str = ', '.join(str(p) for p in init_params)

            if num_events <= num_params:
                print(f"\nWARNING: For Sig {sig_name}, #params={num_params} >= #data={num_events}. Using adapt. reg.")
                fit_convergence_type = "Adapt. Reg."
            else:
                fit_convergence_type = "Red. Chi-Square"

            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Parameters: ({param_str})")
            print(equation_str)
            print("/*--------------------------------------------------*/")

            # Track best across all runs
            best_overall_params   = None
            best_overall_errors   = None
            best_overall_cost     = float('inf')
            best_overall_bin      = None
            best_overall_temp     = float('inf')
            best_overall_prob     = 1.0
            best_overall_residual= float('inf')
            best_overall_ic_aic  = float('inf')
            best_overall_ic_bic  = float('inf')

            total_iteration = 0
            lambda_reg = 0.01
            cost_history = []

            param_offsets = [0.1 for _ in range(num_params)]
            params_sig_history = [[] for _ in range(num_params)]
            graph_sig_params = [TGraph() for _ in range(num_params)]
            graphs_sig_params_all.append(graph_sig_params)
            graph_sig_chi2   = TGraph()
            graph_sig_temp   = TGraph()
            graph_sig_accept = TGraph()
            graph_sig_residuals = TGraph()
            graph_sig_aic    = TGraph()
            graph_sig_bic    = TGraph()
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)
            graphs_sig_residuals.append(graph_sig_residuals)
            graphs_sig_ic_aic.append(graph_sig_aic)
            graphs_sig_ic_bic.append(graph_sig_bic)

            # Draw data
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")
            start_time = time.time()

            # -----------------------------------------------------------------------------
            # 4. Outer loop => multiple runs => helps avoid local minima
            # -----------------------------------------------------------------------------
            for run_idx in range(num_optimizations):
                print(f"\nStarting optimization run {run_idx+1}/{num_optimizations}")
                set_optimization   = full_optimization

                # For example: bin=2
                for b in [2]:
                    print(f"Determining best fit for bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, "
                          f"W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    iteration = 0
                    stagnation_count = 0
                    initial_temperature = 1.0
                    temperature         = initial_temperature

                    # Random initialization in ± initial_param_bounds
                    current_params = [
                        random.uniform(-initial_param_bounds, initial_param_bounds)
                        for _ in range(num_params)
                    ]
                    current_errors = [0.0]*num_params
                    best_params = list(current_params)
                    best_errors = list(current_errors)
                    best_cost   = float('inf')
                    accept_prob = 0.0
                    residual    = float('inf')

                    while iteration <= max_iterations:
                        # Debug prints every 100 iterations
                        if debug and iteration % 100 == 0:
                            print(f"[DEBUG run={run_idx+1}] iteration={iteration}, "
                                  f"best_cost={best_cost:.3f}, temp={temperature:.3f}, "
                                  f"stagnation={stagnation_count}, params={current_params}")

                        # -----------------------------------------------------------------
                        # Simulated annealing step => perturb each param by some fraction
                        # using the current 'temperature'
                        # -----------------------------------------------------------------
                        current_params = [simulated_annealing(p, temperature)
                                          for p in current_params]

                        g_sig_fit = TGraphErrors()
                        graphs_sig_fit.append(g_sig_fit)
                        sys.stdout.write(f"\rSearching for best parameters...({iteration}/{max_iterations})")
                        sys.stdout.flush()

                        try:
                            # Data lines from nsep => build g_sig => fill graphs_sig_fit[it]
                            g_sig = TGraphErrors()
                            for i_data in range(nsep.GetSelectedRows()):
                                x_val = nsep.GetV2()[i_data]
                                y_val = nsep.GetV1()[i_data]
                                y_err = nsep.GetV3()[i_data]
                                g_sig.SetPoint(i_data, x_val, y_val)
                                g_sig.SetPointError(i_data, 0, y_err)

                            for i_pt in range(len(w_vec)):
                                sig_X_fit = g_sig.GetY()[i_pt]# / (g_vec[i_pt])
                                sig_X_fit_err = g_sig.GetEY()[i_pt]# / (g_vec[i_pt])
                                graphs_sig_fit[it].SetPoint(i_pt, g_sig.GetX()[i_pt], sig_X_fit)
                                graphs_sig_fit[it].SetPointError(i_pt, 0, sig_X_fit_err)

                            # Pick function
                            if sig_name == "L":
                                fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                            elif sig_name == "T":
                                fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                            elif sig_name == "LT":
                                fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                            elif sig_name == "TT":
                                fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                            else:
                                raise ValueError("Unknown signal name")
                            funcs_sig.append(fun_Sig)

                            f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
                            fits_sig.append(f_sig)
                            fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
                            for i_par in range(num_params):
                                if abs(current_params[i_par]) > abs(initial_param_bounds) or abs(current_params[i_par]) < 1e-10:
                                    print("!!!!!!!!", abs(current_params[i_par])))
                                    current_params[i_par] = 0.0
                                fits_sig[it].SetParameter(i_par, current_params[i_par])
                                if set_optimization:
                                    fits_sig[it].SetParLimits(i_par, -initial_param_bounds, initial_param_bounds)
                                else:
                                    off = param_offsets[i_par]
                                    fits_sig[it].SetParLimits(
                                        i_par,
                                        current_params[i_par] - off*abs(current_params[i_par]),
                                        current_params[i_par] + off*abs(current_params[i_par])
                                    )

                            # Fit
                            r_sig_fit = graphs_sig_fit[it].Fit(fits_sig[it], "SQ")

                            # Evaluate cost
                            current_cost, lambda_reg = calculate_cost(
                                fits_sig[it], g_sig, current_params,
                                num_events, num_params, lambda_reg
                            )

                            # Simple residual from last data point
                            residual = 0.0
                            for i_pt2 in range(num_events):
                                x_pt = g_sig.GetX()[i_pt2]
                                y_data = g_sig.GetY()[i_pt2]
                                y_err  = g_sig.GetEY()[i_pt2]
                                y_fit  = fits_sig[it].Eval(x_pt)
                                if y_err != 0:
                                    residual = (y_data - y_fit)/y_err
                                else:
                                    residual = (y_data - y_fit)

                            cost_history.append(current_cost)
                            if len(cost_history) >= 2:
                                lambda_reg = adaptive_regularization(cost_history, lambda_reg)

                            accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                            # Update current params from the fit
                            current_params = [fits_sig[it].GetParameter(i_par) for i_par in range(num_params)]
                            current_errors = [fits_sig[it].GetParError(i_par) for i_par in range(num_params)]
                            
                            # Accept or not
                            if accept_prob > random.random():
                                best_params = list(current_params)
                                best_cost   = current_cost
                                best_errors = list(current_errors)
                                if debug and iteration % 50 == 0:
                                    print(f"[DEBUG] ACCEPT => cost={best_cost:.3f}, iteration={iteration}")
                                if current_cost > best_overall_cost:
                                    stagnation_count += 1
                            else:
                                stagnation_count += 1

                            # Local search every 15 iterations
                            if iteration % 15 == 0:
                                current_params = local_search(current_params, fits_sig[it], num_params)

                            # If stalling
                            if stagnation_count > 5:
                                # re-random
                                current_params = [
                                    random.uniform(-initial_param_bounds, initial_param_bounds)
                                    for _ in range(num_params)
                                ]
                                stagnation_count = 0
                                if debug:
                                    print(f"[DEBUG] Stalling => re-random => params={current_params}")

                            # Sync best
                            current_params = list(best_params)
                            # Adjust temperature with your adaptive_cooling
                            temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)
                            iteration += 1

                        except (TypeError, ZeroDivisionError, OverflowError) as e:
                            # On error => re-random
                            if debug:
                                print(f"[DEBUG] Exception => {str(e)}, re-randomizing.")
                            initial_param_bounds = random.uniform(0.0, initial_param_bounds)
                            current_params = [
                                random.uniform(-initial_param_bounds, initial_param_bounds)
                                for _ in range(num_params)
                            ]                        
                            iteration += 1
                            continue

                        # Unconditionally update best_overall_cost
                        if best_cost < best_overall_cost:
                            best_overall_cost    = best_cost
                            best_overall_bin     = b
                            best_overall_params  = best_params[:]
                            best_overall_errors  = best_errors[:]
                            best_overall_temp    = temperature
                            best_overall_prob    = accept_prob
                            best_overall_residual= residual

                        # TGraph updates
                        for i_par in range(num_params):
                            graph_sig_params[i_par].SetPoint(total_iteration, total_iteration, best_params[i_par])
                        graph_sig_chi2.SetPoint(total_iteration, total_iteration, round(best_cost, 4))
                        graph_sig_temp.SetPoint(total_iteration, total_iteration, round(temperature, 4))
                        graph_sig_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                        graph_sig_residuals.SetPoint(total_iteration, total_iteration, round(residual, 4))
                        graph_sig_aic.SetPoint(total_iteration, total_iteration, 0)
                        graph_sig_bic.SetPoint(total_iteration, total_iteration, 0)

                        total_iteration += 1

                    # end while iteration <= max_iterations
                    print(f"\nBest Cost: {best_overall_cost:.3f}")

            # end for run_idx
            try:
                if best_overall_params is not None:
                    print(f"\nBest overall solution: {best_overall_params}")
                    print(f"Best overall cost: {best_overall_cost:.5f}")
                    print(f"Best overall bin: t={t_vec[best_overall_bin]:.3f}, "
                          f"Q2={q2_vec[best_overall_bin]:.3f}, "
                          f"W={w_vec[best_overall_bin]:.3f}, "
                          f"theta={th_vec[best_overall_bin]:.3f}")
                else:
                    print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                    sys.exit(2)                    
            except TypeError:
                print(f"ERROR: Fit failed! Check {equation_str} in input model file...")
                sys.exit(2)

            end_time = time.time()
            print("The loop took {:.2f} seconds.".format(end_time - start_time))
            
            # Make sure best_overall_params has the expected length
            while len(best_overall_params) < num_params:
                best_overall_params.append(0.0)
                best_overall_errors.append(0.0)
            for j in range(num_params):
                par_vec[4*it + j]     = best_overall_params[j]
                par_err_vec[4*it + j] = best_overall_errors[j]
                par_chi2_vec[4*it + j]  = best_overall_cost

            # Plot the final model fit
            g_sig_fit = TGraphErrors()
            graphs_sig_fit.append(g_sig_fit)
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])
            for i in range(len(w_vec)):
                sig_X_fit = g_sig.GetY()[i]# / (g_vec[i])
                sig_X_fit_err = g_sig.GetEY()[i]# / (g_vec[i])
                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)
            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle(f"#left(#frac{{#it{{d#sigma}}}}{{#it{{dt}}}}#right)_{{{sig_name}}} [nb/GeV^{2}]")
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
            # Now create and draw a final TF1 with the best parameters fixed
            if sig_name == "L":
                fun_Sig = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "T":
                fun_Sig = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "LT":
                fun_Sig = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "TT":
                fun_Sig = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            funcs_sig.append(fun_Sig)                
            f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
            fits_sig.append(f_sig)
            fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
            for i in range(num_params):
                fits_sig[it].FixParameter(i, best_overall_params[i])
            n_points = 100
            fit_y_values = [fits_sig[it].Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            r_sig_fit = graphs_sig_fit[it].Fit(fits_sig[it], "SQ")
            fits_sig[it].Draw("same")
            f_sig_status = "Fit Successful" if fits_sig[it].GetNDF() != 0 else "Fit Failed"
            fit_status = TText()
            fit_status.SetTextSize(0.04)
            fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: " + f_sig_status)

            # Plot the parameter convergence on canvas c3 (each parameter drawn with its own color)
            c3.cd(it+1).SetLeftMargin(0.12)
            for i in range(num_params):
                if i == 0:
                    graph = graph_sig_params[i]
                    graph.SetTitle(f"Sig {sig_name} Parameter Convergence;Optimization Run;Parameter")
                    graph.SetLineColor(colors[i])
                    graph.Draw("ALP")
                else:
                    graph = graph_sig_params[i]
                    graph.SetLineColor(colors[i])
                    graph.Draw("LP SAME")
            c3.Update()

            # Plot chi2 convergence
            c4.cd(it+1).SetLeftMargin(0.12)
            graph_sig_chi2.SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
            graph_sig_chi2.SetLineColor(1)
            graph_sig_chi2.Draw("ALP")
            latex = TLatex()
            latex.SetTextSize(0.04)
            latex.SetNDC(True)
            latex.DrawLatex(0.35, 0.85, f"Best #chi^{{2}}: {best_overall_cost:.3f}")
            c4.Update()

            # Plot temperature, acceptance probability, residuals, and information criteria on their canvases
            c5.cd(it+1).SetLeftMargin(0.12)
            graph_sig_temp.SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
            graph_sig_temp.SetLineColor(1)
            graph_sig_temp.Draw("ALP")
            c5.Update()

            c6.cd(it+1).SetLeftMargin(0.12)
            graph_sig_accept.SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
            graph_sig_accept.SetLineColor(1)
            graph_sig_accept.Draw("ALP")
            c6.Update()

            c7.cd(it+1).SetLeftMargin(0.12)
            graph_sig_residuals.SetTitle(f"Sig {sig_name} Residuals Evolution")
            graph_sig_residuals.SetLineColor(1)
            graph_sig_residuals.Draw("ALP")
            c7.Update()

            c8.cd(it+1).SetLeftMargin(0.12)
            graph_sig_aic.SetTitle(f"Sig {sig_name} Information Criteria Evolution")
            graph_sig_aic.SetLineColor(kRed)
            graph_sig_aic.Draw("ALP")
            graph_sig_bic.SetLineColor(kGreen)
            graph_sig_bic.Draw("same")
            leg = TLegend(0.7, 0.7, 0.9, 0.9)
            leg.AddEntry(graph_sig_aic, "AIC", "lp")
            leg.AddEntry(graph_sig_bic, "BIC", "lp")
            leg.Draw()
            c8.Update()
            print("\n")
        else:

            num_params, init_params, equation_str = inpDict["initial_params"](sig_name, val)
            param_str = ', '.join(str(p) for p in init_params)
            print("\n/*--------------------------------------------------*/")
            print(f"Fit for Sig {sig_name} ({num_params} parameters)")
            print(f"Initial Parameters: ({param_str})")
            print(equation_str)
            print("/*--------------------------------------------------*/")

            # Track best across all runs
            best_overall_params   = None
            best_overall_errors   = None
            best_overall_cost     = float('inf')
            best_overall_bin      = None
            best_overall_temp     = float('inf')
            best_overall_prob     = 1.0
            best_overall_residual= float('inf')
            best_overall_ic_aic  = float('inf')
            best_overall_ic_bic  = float('inf')
            
            total_iteration = 0
            lambda_reg = 0.0
            cost_history = []
            
            graph_sig_params = [TGraph() for _ in range(num_params)]
            graphs_sig_params_all.append(graph_sig_params)
            graph_sig_chi2   = TGraph()
            graph_sig_temp   = TGraph()
            graph_sig_accept = TGraph()
            graph_sig_residuals = TGraph()
            graph_sig_aic    = TGraph()
            graph_sig_bic    = TGraph()
            graphs_sig_converge.append(graph_sig_chi2)
            graphs_sig_temp.append(graph_sig_temp)
            graphs_sig_accept.append(graph_sig_accept)
            graphs_sig_residuals.append(graph_sig_residuals)
            graphs_sig_ic_aic.append(graph_sig_aic)
            graphs_sig_ic_bic.append(graph_sig_bic)

            # Draw data
            nsep.Draw(f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e", "", "goff")

            for b in [2]:
                print(f"Determining best fit for bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, "
                      f"W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                g_sig_fit = TGraphErrors()
                graphs_sig_fit.append(g_sig_fit)
                g_sig = TGraphErrors()
                for i_data in range(nsep.GetSelectedRows()):
                    x_val = nsep.GetV2()[i_data]
                    y_val = nsep.GetV1()[i_data]
                    y_err = nsep.GetV3()[i_data]
                    g_sig.SetPoint(i_data, x_val, y_val)
                    g_sig.SetPointError(i_data, 0, y_err)
                for i_pt in range(len(w_vec)):
                    sig_X_fit = g_sig.GetY()[i_pt]# / (g_vec[i_pt])
                    sig_X_fit_err = g_sig.GetEY()[i_pt]# / (g_vec[i_pt])
                    graphs_sig_fit[it].SetPoint(i_pt, g_sig.GetX()[i_pt], sig_X_fit)
                    graphs_sig_fit[it].SetPointError(i_pt, 0, sig_X_fit_err)
                if sig_name == "L":
                    fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "T":
                    fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "LT":
                    fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                elif sig_name == "TT":
                    fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], th_vec[b])
                funcs_sig.append(fun_Sig)                    
                f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
                fits_sig.append(f_sig)
                fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
                for i in range(num_params):
                    fits_sig[it].FixParameter(i, par_vec[4*it + i])
                r_sig_fit = graphs_sig_fit[it].Fit(fits_sig[it], "SQ")
                current_cost, lambda_reg = calculate_cost(
                    fits_sig[it], g_sig, par_vec[4*it:4*(it+1)],
                    num_events, num_params, lambda_reg
                )
                print(f"\tCost: {current_cost:.3f}")
                if current_cost < best_overall_cost:
                    best_overall_cost = current_cost
                    best_overall_bin = b
                    best_overall_params = [par_vec[4*it + j] for j in range(num_params)]
            print(f"\nBest overall solution: {best_overall_params}")
            print(f"Best overall cost: {best_overall_cost:.5f}")
            
            for j in range(num_params):
                par_chi2_vec[4*it + j] = best_overall_cost
                
            # Plot the final model fit
            g_sig_fit = TGraphErrors()
            graphs_sig_fit.append(g_sig_fit)
            g_sig = TGraphErrors()
            for i in range(nsep.GetSelectedRows()):
                g_sig.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                g_sig.SetPointError(i, 0, nsep.GetV3()[i])
            for i in range(len(w_vec)):
                sig_X_fit = g_sig.GetY()[i]# / (g_vec[i])
                sig_X_fit_err = g_sig.GetEY()[i]# / (g_vec[i])
                graphs_sig_fit[it].SetPoint(i, g_sig.GetX()[i], sig_X_fit)
                graphs_sig_fit[it].SetPointError(i, 0, sig_X_fit_err)                
            c2.cd(it+1).SetLeftMargin(0.12)
            graphs_sig_fit[it].SetTitle(f"Sigma {sig_name} Model Fit")
            graphs_sig_fit[it].Draw("A*")
            graphs_sig_fit[it].GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
            graphs_sig_fit[it].GetXaxis().CenterTitle()
            graphs_sig_fit[it].GetYaxis().SetTitle(f"#left(#frac{{#it{{d#sigma}}}}{{#it{{dt}}}}#right)_{{{sig_name}}} [nb/GeV^{2}]")
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
            if sig_name == "L":
                fun_Sig = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "T":
                fun_Sig = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "LT":
                fun_Sig = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            elif sig_name == "TT":
                fun_Sig = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], th_vec[best_overall_bin])
            funcs_sig.append(fun_Sig)
            f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
            fits_sig.append(f_sig)
            fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
            for i in range(num_params):
                fits_sig[it].FixParameter(i, par_vec[4*it + i])
            n_points = 100
            fit_y_values = [fits_sig[it].Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
            fit_y_min = min(fit_y_values)
            fit_y_max = max(fit_y_values)
            y_min = min(y_min, fit_y_min)
            y_max = max(y_max, fit_y_max)
            margin = 0.1 * (y_max - y_min)
            graphs_sig_fit[it].GetYaxis().SetRangeUser(y_min - margin, y_max + margin)
            r_sig_fit = graphs_sig_fit[it].Fit(fits_sig[it], "SQ")
            fits_sig[it].Draw("same")
            latex = TLatex()
            latex.SetTextSize(0.04)
            latex.SetNDC(True)
            latex.DrawLatex(0.35, 0.85, f"Best #chi^{{2}}: {best_overall_cost:.3f}")
            c2.Update()
            print("\n")
        c2.Update()
        c3.Update()
        c4.Update()
        c5.Update()
        c6.Update()
        c7.Update()
        c8.Update()
        # Print all canvases to the output PDF
        c2.Print(outputpdf+'(')
        c3.Print(outputpdf)
        c4.Print(outputpdf)
        c5.Print(outputpdf)
        c6.Print(outputpdf)
        c7.Print(outputpdf)
        c8.Print(outputpdf+')')
    print(f"\nFits saved to {outputpdf}...")

