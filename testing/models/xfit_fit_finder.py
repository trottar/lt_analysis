#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-05-23 01:45:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import numpy as np
from ROOT import TGraph, TGraphErrors, TF1, TCanvas, TText, TLatex, TLegend, kRed, kBlue, kGreen, kMagenta, kBlack
import sys, math, time, random, re

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

def compute_iteration_ic(red_chi2, num_events, num_params):
    ndf = num_events - num_params
    if ndf <= 0 or not math.isfinite(red_chi2):
        return 0.0, 0.0
    return calculate_information_criteria(num_events, num_params, red_chi2 * ndf)

def build_signal_graph(nsep, sig_name, selection_expr=""):
    draw_expr = f"sig{sig_name.lower()}:t:sig{sig_name.lower()}_e"
    nsep.Draw(draw_expr, selection_expr, "goff")
    n_points = nsep.GetSelectedRows()
    graph = TGraphErrors()
    for i_point in range(n_points):
        x_val = nsep.GetV2()[i_point]
        y_val = nsep.GetV1()[i_point]
        y_err = nsep.GetV3()[i_point]
        graph.SetPoint(i_point, x_val, y_val)
        graph.SetPointError(i_point, 0, y_err)
    return graph, n_points

def infer_positive_parameter_indices(equation_str):
    ordered_tokens = []
    for token in re.findall(r"par\d+", equation_str):
        if token not in ordered_tokens:
            ordered_tokens.append(token)
    token_to_index = {token: idx for idx, token in enumerate(ordered_tokens)}
    positive_indices = set()
    for token in re.findall(r"math\.exp\(\s*-\s*(par\d+)\s*\*", equation_str):
        if token in token_to_index:
            positive_indices.add(token_to_index[token])
    return positive_indices

def build_parameter_bounds(init_params, equation_str, initial_param_bounds):
    positive_indices = infer_positive_parameter_indices(equation_str)
    lower_bounds = []
    upper_bounds = []
    bound_cap = max(5.0, min(abs(initial_param_bounds), 1.0e4))

    for i_par, init_value in enumerate(init_params):
        center = init_value if math.isfinite(init_value) else 0.0
        span = min(bound_cap, max(5.0, 10.0 * max(abs(center), 1.0)))
        if i_par in positive_indices:
            lower_bound = 0.0
            upper_bound = min(bound_cap, max(center + span, span))
        else:
            lower_bound = max(-bound_cap, center - span)
            upper_bound = min(bound_cap, center + span)
        if lower_bound >= upper_bound:
            lower_bound = max(-bound_cap, center - 1.0)
            upper_bound = min(bound_cap, center + 1.0)
        lower_bounds.append(lower_bound)
        upper_bounds.append(upper_bound)

    return lower_bounds, upper_bounds

def format_metric(value, precision=3, sci_threshold=1.0e-3, large_threshold=1.0e4):
    if not math.isfinite(value):
        return str(value)
    abs_value = abs(value)
    if abs_value != 0.0 and (abs_value < sci_threshold or abs_value >= large_threshold):
        return f"{value:.6e}"
    return f"{value:.{precision}f}"

def set_graph_vertical_range(graph, pad_fraction=0.1, min_pad=1.0e-6):
    n_points = graph.GetN()
    if n_points <= 0:
        return
    y_values = [graph.GetY()[i] for i in range(n_points) if math.isfinite(graph.GetY()[i])]
    if not y_values:
        return
    y_min = min(y_values)
    y_max = max(y_values)
    span = y_max - y_min
    if span <= 0.0:
        pad = max(abs(y_min) * pad_fraction, min_pad)
    else:
        pad = max(span * pad_fraction, min_pad)
    graph.SetMinimum(y_min - pad)
    graph.SetMaximum(y_max + pad)

def parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec, fixed_params, outputpdf, full_optimization=True, debug=False):
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
    graphs_sig_fit      = [None] * 4  # one graph per signal index (L/T/LT/TT)
    graphs_sig_params_all = []
    graphs_sig_converge = []
    graphs_sig_temp     = []
    graphs_sig_accept   = []
    graphs_sig_residuals= []
    graphs_sig_ic_aic   = []
    graphs_sig_ic_bic   = []    
    funcs_sig = [None]*4 # L, T, LT, TT
    fits_sig = [None]*4 # L, T, LT, TT
    
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
    c9 = TCanvas("c9", "Search Diagnostics", 900, 900)
    c9.Divide(2, 2)
    c10 = TCanvas("c10", "Search Control", 900, 900)
    c10.Divide(2, 2)

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

    params_used = []

    # -----------------------------------------------------------------------------
    # 3. Main loop over each fit in fit_params
    # -----------------------------------------------------------------------------
    for it, (sig_name, val) in enumerate(fit_params.items()):

        # Track best across all runs
        best_overall_params   = None
        best_overall_errors   = None
        best_overall_cost     = float('inf')
        best_overall_chi2     = float('inf')
        best_overall_chi_square = float('inf')
        best_overall_ndf      = 0
        best_overall_lambda   = float('inf')
        best_overall_bin      = None
        best_overall_run_idx  = None
        best_overall_temp     = float('inf')
        best_overall_prob     = 1.0
        best_overall_residual= float('inf')
        best_overall_rms_pull = float('inf')
        best_overall_restarts = 0
        best_overall_ic_aic  = float('inf')
        best_overall_ic_bic  = float('inf')

        # Keep graph indexing aligned with signal indexing.
        graphs_sig_fit[it] = TGraphErrors()

        for b in range(len(t_vec)):
            if sig_name not in fixed_params:

                num_params, init_params, equation_str = inpDict["initial_params"](sig_name, val)
                # Avoid zero init param
                init_params = [v if abs(v) > 0.0 else 1e-15 for v in init_params]
                param_str = ', '.join(str(p) for p in init_params)

                print("\n/*--------------------------------------------------*/")
                print(f"Fit for Sig {sig_name} ({num_params} parameters)")
                print(f"Initial Parameters: ({param_str})")
                print(equation_str)
                print("/*--------------------------------------------------*/")

                if sig_name == "L":
                    fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "T":
                    fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "LT":
                    fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "TT":
                    fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                else:
                    raise ValueError("Unknown signal name") 
                funcs_sig[it] = fun_Sig                           
                f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
                fits_sig[it] = f_sig

                param_offsets = [0.1 for _ in range(num_params)]
                param_lower_bounds, param_upper_bounds = build_parameter_bounds(
                    init_params, equation_str, initial_param_bounds
                )
                params_sig_history = [[] for _ in range(num_params)]
                graph_sig_params = [TGraph() for _ in range(num_params)]
                graphs_sig_params_all.append(graph_sig_params)
                graph_sig_chi2   = TGraph()
                graph_sig_temp   = TGraph()
                graph_sig_accept = TGraph()
                graph_sig_residuals = TGraph()
                graph_sig_aic    = TGraph()
                graph_sig_bic    = TGraph()
                graph_sig_cost_current = TGraph()
                graph_sig_cost_best = TGraph()
                graph_sig_chi2_current = TGraph()
                graph_sig_chi2_best = TGraph()
                graph_sig_lambda = TGraph()
                graph_sig_stagnation = TGraph()
                graph_sig_restarts = TGraph()
                graphs_sig_converge.append(graph_sig_chi2)
                graphs_sig_temp.append(graph_sig_temp)
                graphs_sig_accept.append(graph_sig_accept)
                graphs_sig_residuals.append(graph_sig_residuals)
                graphs_sig_ic_aic.append(graph_sig_aic)
                graphs_sig_ic_bic.append(graph_sig_bic)

                g_sig, fit_num_events = build_signal_graph(nsep, sig_name)
                if fit_num_events == 0:
                    print(f"WARNING: No data points selected for Sig {sig_name}. Skipping.")
                    continue

                if fit_num_events <= num_params:
                    print(f"\nWARNING: For Sig {sig_name}, #params={num_params} >= #data={fit_num_events}. Using adapt. reg.")
                    fit_convergence_type = "Adapt. Reg."
                else:
                    fit_convergence_type = "Red. Chi-Square"

                total_iteration = 0
                restart_patience = max(250, 100 * num_params)
                local_search_interval = max(100, 40 * num_params)
                start_time = time.time()

                for i_pt in range(fit_num_events):
                    sig_X_fit = g_sig.GetY()[i_pt]
                    sig_X_fit_err = g_sig.GetEY()[i_pt]
                    graphs_sig_fit[it].SetPoint(i_pt, g_sig.GetX()[i_pt], sig_X_fit)
                    graphs_sig_fit[it].SetPointError(i_pt, 0, sig_X_fit_err)

                # -----------------------------------------------------------------------------
                # 4. Outer loop => multiple runs => helps avoid local minima
                # -----------------------------------------------------------------------------
                for run_idx in range(num_optimizations):
                    print(f"\nStarting optimization run {run_idx+1}/{num_optimizations}")
                    set_optimization   = full_optimization
                    lambda_reg = 0.01
                    cost_history = []

                    print(f"Determining best fit for bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, "
                          f"W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")
                    iteration = 0
                    stagnation_count = 0
                    initial_temperature = 1.0
                    temperature         = initial_temperature

                    def sample_random_params():
                        return [
                            random.uniform(low, high)
                            for low, high in zip(param_lower_bounds, param_upper_bounds)
                        ]

                    if run_idx == 0:
                        current_params = [
                            select_valid_parameter(value, low, high)
                            for value, low, high in zip(init_params, param_lower_bounds, param_upper_bounds)
                        ]
                    else:
                        current_params = sample_random_params()
                    current_errors = [0.0]*num_params
                    current_cost   = float('inf')
                    current_chi2   = float('inf')
                    current_chi_square = float('inf')
                    current_ndf    = 0
                    current_residual = float('inf')
                    current_rms_pull = float('inf')
                    current_lambda = float('inf')
                    best_params = list(current_params)
                    best_errors = list(current_errors)
                    best_cost   = float('inf')
                    best_chi2   = float('inf')
                    best_chi_square = float('inf')
                    best_ndf    = 0
                    best_residual = float('inf')
                    best_rms_pull = float('inf')
                    best_lambda = float('inf')
                    accept_prob = 0.0
                    restart_count = 0

                    def evaluate_candidate(seed_params, lambda_value):
                        trial_params = list(seed_params)
                        for i_par in range(num_params):
                            trial_params[i_par] = select_valid_parameter(
                                trial_params[i_par], param_lower_bounds[i_par], param_upper_bounds[i_par]
                            )
                            if abs(trial_params[i_par]) < 1e-15:
                                trial_params[i_par] = 0.0
                            fits_sig[it].SetParameter(i_par, trial_params[i_par])
                            if set_optimization:
                                fits_sig[it].SetParLimits(i_par, param_lower_bounds[i_par], param_upper_bounds[i_par])
                            else:
                                off = param_offsets[i_par]
                                local_half_width = max(off * max(abs(trial_params[i_par]), 1.0), 1e-6)
                                local_low = max(param_lower_bounds[i_par], trial_params[i_par] - local_half_width)
                                local_high = min(param_upper_bounds[i_par], trial_params[i_par] + local_half_width)
                                if local_low >= local_high:
                                    local_low = param_lower_bounds[i_par]
                                    local_high = param_upper_bounds[i_par]
                                fits_sig[it].SetParLimits(
                                    i_par,
                                    local_low,
                                    local_high
                                )

                        graphs_sig_fit[it].Fit(fits_sig[it], "SQ")
                        fitted_params = [
                            select_valid_parameter(
                                fits_sig[it].GetParameter(i_par),
                                param_lower_bounds[i_par],
                                param_upper_bounds[i_par]
                            )
                            for i_par in range(num_params)
                        ]
                        fitted_errors = [fits_sig[it].GetParError(i_par) for i_par in range(num_params)]
                        if not all(math.isfinite(value) for value in fitted_params + fitted_errors):
                            raise OverflowError("Non-finite fit parameter encountered")
                        for i_par in range(num_params):
                            fits_sig[it].SetParameter(i_par, fitted_params[i_par])
                        fitted_cost, next_lambda, fitted_chi2 = calculate_cost(
                            fits_sig[it], g_sig, fitted_params, fit_num_events, num_params, lambda_value
                        )
                        fitted_chi_square = fits_sig[it].GetChisquare()
                        fitted_ndf = fits_sig[it].GetNDF()

                        residuals = []
                        for i_pt2 in range(fit_num_events):
                            x_pt = g_sig.GetX()[i_pt2]
                            y_data = g_sig.GetY()[i_pt2]
                            y_err = g_sig.GetEY()[i_pt2]
                            y_fit = fits_sig[it].Eval(x_pt)
                            if not math.isfinite(y_fit):
                                residuals.append(float('inf'))
                                continue
                            if y_err != 0:
                                residual_value = (y_data - y_fit)/y_err
                            else:
                                residual_value = (y_data - y_fit)
                            residuals.append(residual_value)

                        max_abs_pull = max(abs(res) for res in residuals) if residuals else float('inf')
                        rms_pull = math.sqrt(sum(res*res for res in residuals) / len(residuals)) if residuals else float('inf')

                        return (
                            fitted_params, fitted_errors, fitted_cost, next_lambda, fitted_chi2,
                            fitted_chi_square, fitted_ndf, max_abs_pull, rms_pull
                        )

                    while iteration < max_iterations:
                        # Debug prints every 100 iterations
                        if debug and iteration % 100 == 0:
                            print(f"[DEBUG run={run_idx+1}] iteration={iteration}, "
                                  f"current_cost={current_cost:.3f}, run_best={best_cost:.3f}, temp={temperature:.3f}, "
                                  f"stagnation={stagnation_count}, params={current_params}")

                        # -----------------------------------------------------------------
                        # Simulated annealing step => perturb each param by some fraction
                        # using the current 'temperature'
                        # -----------------------------------------------------------------
                        proposal_params = [
                            simulated_annealing(p, temperature)
                            for p in current_params
                        ]

                        sys.stdout.write(
                            f"\rSearching for best parameters...({iteration}/{max_iterations}) "
                            f"T={temperature:.4f}"
                        )
                        sys.stdout.flush()

                        try:
                            fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
                            candidate_params, candidate_errors, candidate_cost, candidate_lambda, candidate_chi2, candidate_chi_square, candidate_ndf, candidate_residual, candidate_rms_pull = evaluate_candidate(
                                proposal_params, lambda_reg
                            )

                            cost_history.append(candidate_cost)
                            if len(cost_history) >= 2:
                                lambda_reg = adaptive_regularization(cost_history, lambda_reg)

                            accept_prob = acceptance_probability(current_cost, candidate_cost, temperature)

                            if accept_prob > random.random():
                                current_params = list(candidate_params)
                                current_errors = list(candidate_errors)
                                current_cost = candidate_cost
                                current_chi2 = candidate_chi2
                                current_chi_square = candidate_chi_square
                                current_ndf = candidate_ndf
                                current_residual = candidate_residual
                                current_rms_pull = candidate_rms_pull
                                current_lambda = candidate_lambda

                                if current_cost < best_cost:
                                    best_params = list(current_params)
                                    best_errors = list(current_errors)
                                    best_cost = current_cost
                                    best_chi2 = current_chi2
                                    best_chi_square = current_chi_square
                                    best_ndf = current_ndf
                                    best_residual = current_residual
                                    best_rms_pull = current_rms_pull
                                    best_lambda = current_lambda
                                    stagnation_count = 0
                                else:
                                    stagnation_count += 1

                                if debug and iteration % 50 == 0:
                                    print(f"[DEBUG] ACCEPT => cost={current_cost:.3f}, chi2/ndf={current_chi2:.3f}, iteration={iteration}")
                            else:
                                stagnation_count += 1

                            if (
                                stagnation_count > 0
                                and stagnation_count % local_search_interval == 0
                                and math.isfinite(current_cost)
                                and math.isfinite(best_cost)
                            ):
                                local_params, local_cost, local_chi2 = local_search(
                                    current_params, fits_sig[it], g_sig,
                                    fit_num_events, num_params, lambda_reg,
                                    lower_bounds=param_lower_bounds,
                                    upper_bounds=param_upper_bounds
                                )
                                if local_cost < current_cost:
                                    local_params, local_errors, local_cost, local_lambda, local_chi2, local_chi_square, local_ndf, local_residual, local_rms_pull = evaluate_candidate(
                                        local_params, lambda_reg
                                    )
                                    current_params = list(local_params)
                                    current_errors = list(local_errors)
                                    current_cost = local_cost
                                    current_chi2 = local_chi2
                                    current_chi_square = local_chi_square
                                    current_ndf = local_ndf
                                    current_residual = local_residual
                                    current_rms_pull = local_rms_pull
                                    current_lambda = local_lambda
                                    if current_cost < best_cost:
                                        best_params = list(current_params)
                                        best_errors = list(current_errors)
                                        best_cost = current_cost
                                        best_chi2 = current_chi2
                                        best_chi_square = current_chi_square
                                        best_ndf = current_ndf
                                        best_residual = current_residual
                                        best_rms_pull = current_rms_pull
                                        best_lambda = current_lambda
                                    stagnation_count = 0

                            if stagnation_count >= restart_patience:
                                current_params = sample_random_params()
                                current_errors = [0.0] * num_params
                                current_cost = float('inf')
                                current_chi2 = float('inf')
                                current_chi_square = float('inf')
                                current_ndf = 0
                                current_residual = float('inf')
                                current_rms_pull = float('inf')
                                current_lambda = float('inf')
                                restart_count += 1
                                stagnation_count = 0
                                if debug:
                                    print(f"[DEBUG] Stalling => re-random => params={current_params}")

                            temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)
                            iteration += 1

                        except (TypeError, ZeroDivisionError, OverflowError) as e:
                            if debug:
                                print(f"[DEBUG] Exception => {str(e)}, re-randomizing.")
                            current_params = sample_random_params()
                            current_errors = [0.0] * num_params
                            current_cost = float('inf')
                            current_chi2 = float('inf')
                            current_chi_square = float('inf')
                            current_ndf = 0
                            current_residual = float('inf')
                            current_rms_pull = float('inf')
                            current_lambda = float('inf')
                            restart_count += 1
                            stagnation_count = 0
                            iteration += 1
                            continue

                        if best_cost < best_overall_cost:
                            best_overall_cost    = best_cost
                            best_overall_chi2    = best_chi2
                            best_overall_chi_square = best_chi_square
                            best_overall_ndf     = best_ndf
                            best_overall_lambda  = best_lambda
                            best_overall_bin     = b
                            best_overall_run_idx = run_idx + 1
                            best_overall_params  = best_params[:]
                            best_overall_errors  = best_errors[:]
                            best_overall_temp    = temperature
                            best_overall_prob    = accept_prob
                            best_overall_residual= best_residual
                            best_overall_rms_pull = best_rms_pull
                            best_overall_restarts = restart_count

                        current_aic, current_bic = compute_iteration_ic(current_chi2, fit_num_events, num_params)
                        best_aic, best_bic = compute_iteration_ic(best_chi2, fit_num_events, num_params)

                        plotted_params = best_params if math.isfinite(best_cost) else current_params
                        plotted_best_cost = best_cost if math.isfinite(best_cost) else current_cost
                        plotted_best_chi2 = best_chi2 if math.isfinite(best_chi2) else current_chi2

                        for i_par in range(num_params):
                            graph_sig_params[i_par].SetPoint(total_iteration, total_iteration, plotted_params[i_par])
                        graph_sig_chi2.SetPoint(total_iteration, total_iteration, plotted_best_chi2)
                        graph_sig_temp.SetPoint(total_iteration, total_iteration, temperature)
                        graph_sig_accept.SetPoint(total_iteration, total_iteration, accept_prob)
                        graph_sig_residuals.SetPoint(total_iteration, total_iteration, current_residual)
                        graph_sig_aic.SetPoint(total_iteration, total_iteration, best_aic)
                        graph_sig_bic.SetPoint(total_iteration, total_iteration, best_bic)
                        graph_sig_cost_current.SetPoint(total_iteration, total_iteration, current_cost)
                        graph_sig_cost_best.SetPoint(total_iteration, total_iteration, plotted_best_cost)
                        graph_sig_chi2_current.SetPoint(total_iteration, total_iteration, current_chi2)
                        graph_sig_chi2_best.SetPoint(total_iteration, total_iteration, plotted_best_chi2)
                        graph_sig_lambda.SetPoint(total_iteration, total_iteration, lambda_reg)
                        graph_sig_stagnation.SetPoint(total_iteration, total_iteration, stagnation_count)
                        graph_sig_restarts.SetPoint(total_iteration, total_iteration, restart_count)

                        total_iteration += 1

                    # end while iteration <= max_iterations
                    print(f"\nRun Best Objective Cost: {format_metric(best_cost, precision=3)}")
                    print(f"Run Best chi2/ndf: {format_metric(best_chi2, precision=3)}")
                    print(
                        f"Run Diagnostics: chi2={best_chi_square:.6e}, ndf={best_ndf}, "
                        f"max|pull|={best_residual:.6e}, rms_pull={best_rms_pull:.6e}, "
                        f"lambda={best_lambda:.6e}, restarts={restart_count}, "
                        f"npts={fit_num_events}, npar={num_params}"
                    )

                # end for run_idx
                try:
                    if best_overall_params is not None:
                        print(f"\nBest overall solution: {best_overall_params}")
                        print(f"Best overall objective cost: {format_metric(best_overall_cost, precision=5)}")
                        print(f"Best overall chi2/ndf: {format_metric(best_overall_chi2, precision=5)}")
                        print(
                            f"Best overall diagnostics: chi2={best_overall_chi_square:.6e}, "
                            f"ndf={best_overall_ndf}, max|pull|={best_overall_residual:.6e}, "
                            f"rms_pull={best_overall_rms_pull:.6e}, lambda={best_overall_lambda:.6e}, "
                            f"restarts={best_overall_restarts}, run={best_overall_run_idx}, "
                            f"npts={fit_num_events}, npar={num_params}"
                        )
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
                    par_chi2_vec[4*it + j]  = best_overall_chi2

                # Plot the final model fit
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
                fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
                for i in range(num_params):
                    fits_sig[it].FixParameter(i, best_overall_params[i])
                try:
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
                except OverflowError:
                    print("WARNING: OverflowError during final fit drawing. Skipping fit overlay. Likely due to no t-dependence in model.")                
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
                set_graph_vertical_range(graph_sig_chi2)
                set_graph_vertical_range(graph_sig_chi2_current)
                graph_sig_chi2.SetTitle(f"Sig {sig_name} {fit_convergence_type} Convergence;Optimization Run;{fit_convergence_type}")
                graph_sig_chi2.SetLineColor(kBlue)
                graph_sig_chi2.Draw("ALP")
                graph_sig_chi2_current.SetLineColor(kRed)
                graph_sig_chi2_current.Draw("LP SAME")
                latex = TLatex()
                latex.SetTextSize(0.04)
                latex.SetNDC(True)
                latex.DrawLatex(0.35, 0.85, f"Best #chi^{{2}}/ndf: {best_overall_chi2:.3f}")
                leg_chi2 = TLegend(0.6, 0.7, 0.9, 0.9)
                leg_chi2.AddEntry(graph_sig_chi2, "Run Best", "lp")
                leg_chi2.AddEntry(graph_sig_chi2_current, "Current", "lp")
                leg_chi2.Draw()
                c4.Update()

                # Plot temperature, acceptance probability, residuals, and information criteria on their canvases
                c5.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_temp)
                graph_sig_temp.SetTitle(f"Sig {sig_name} Temperature Convergence;Optimization Run;Temperature")
                graph_sig_temp.SetLineColor(1)
                graph_sig_temp.Draw("ALP")
                c5.Update()

                c6.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_accept)
                graph_sig_accept.SetTitle(f"Sig {sig_name} Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
                graph_sig_accept.SetLineColor(1)
                graph_sig_accept.Draw("ALP")
                c6.Update()

                c7.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_residuals)
                graph_sig_residuals.SetTitle(f"Sig {sig_name} Residuals Evolution")
                graph_sig_residuals.SetLineColor(1)
                graph_sig_residuals.Draw("ALP")
                c7.Update()

                c8.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_aic)
                set_graph_vertical_range(graph_sig_bic)
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

                c9.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_cost_current)
                set_graph_vertical_range(graph_sig_cost_best)
                graph_sig_cost_current.SetTitle(f"Sig {sig_name} Objective Cost;Optimization Run;Objective Cost")
                graph_sig_cost_current.SetLineColor(kBlack)
                graph_sig_cost_current.Draw("ALP")
                graph_sig_cost_best.SetLineColor(kBlue)
                graph_sig_cost_best.Draw("LP SAME")
                leg_cost = TLegend(0.6, 0.7, 0.9, 0.9)
                leg_cost.AddEntry(graph_sig_cost_current, "Current", "lp")
                leg_cost.AddEntry(graph_sig_cost_best, "Run Best", "lp")
                leg_cost.Draw()
                c9.Update()

                c10.cd(it+1).SetLeftMargin(0.12)
                set_graph_vertical_range(graph_sig_stagnation)
                set_graph_vertical_range(graph_sig_restarts)
                set_graph_vertical_range(graph_sig_lambda)
                graph_sig_stagnation.SetTitle(f"Sig {sig_name} Search Control;Optimization Run;Control Value")
                graph_sig_stagnation.SetLineColor(kBlack)
                graph_sig_stagnation.Draw("ALP")
                graph_sig_restarts.SetLineColor(kRed)
                graph_sig_restarts.Draw("LP SAME")
                graph_sig_lambda.SetLineColor(kMagenta)
                graph_sig_lambda.Draw("LP SAME")
                leg_control = TLegend(0.55, 0.65, 0.9, 0.9)
                leg_control.AddEntry(graph_sig_stagnation, "Stagnation", "lp")
                leg_control.AddEntry(graph_sig_restarts, "Restarts", "lp")
                leg_control.AddEntry(graph_sig_lambda, "Lambda", "lp")
                leg_control.Draw()
                c10.Update()
                print("\n")
            else:

                num_params, init_params, equation_str = inpDict["initial_params"](sig_name, val)
                param_str = ', '.join(str(p) for p in init_params)

                total_iteration = 0
                lambda_reg = 0.0

                print("\n/*--------------------------------------------------*/")
                print(f"Fit for Sig {sig_name} ({num_params} parameters)")
                print(f"Initial Parameters: ({param_str})")
                print(equation_str)
                print("/*--------------------------------------------------*/")                

                if sig_name == "L":
                    fun_Sig = fun_Sig_L_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "T":
                    fun_Sig = fun_Sig_T_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "LT":
                    fun_Sig = fun_Sig_LT_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                elif sig_name == "TT":
                    fun_Sig = fun_Sig_TT_wrapper(g_vec[b], q2_vec[b], w_vec[b], math.radians(th_vec[b]))
                else:
                    raise ValueError("Unknown signal name")
                funcs_sig[it] = fun_Sig                
                f_sig = TF1(f"sig_{sig_name}", funcs_sig[it], tmin_range, tmax_range, num_params)
                fits_sig[it] = f_sig

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

                print(f"Determining best fit for bin: t={t_vec[b]:.3f}, Q2={q2_vec[b]:.3f}, "
                      f"W={w_vec[b]:.3f}, theta={th_vec[b]:.3f}")

                g_sig, fit_num_events = build_signal_graph(nsep, sig_name)
                if fit_num_events == 0:
                    print(f"WARNING: No data points selected for Sig {sig_name}. Skipping.")
                    continue

                for i_pt in range(fit_num_events):
                    sig_X_fit = g_sig.GetY()[i_pt]# / (g_vec[i_pt])# / 1e3
                    sig_X_fit_err = g_sig.GetEY()[i_pt]# / (g_vec[i_pt])# / 1e3
                    graphs_sig_fit[it].SetPoint(i_pt, g_sig.GetX()[i_pt], sig_X_fit)
                    graphs_sig_fit[it].SetPointError(i_pt, 0, sig_X_fit_err)
                fits_sig[it].SetParNames(*[f"p{4*it + i}" for i in range(num_params)])
                for i in range(num_params):
                    fits_sig[it].FixParameter(i, par_vec[4*it + i])
                r_sig_fit = graphs_sig_fit[it].Fit(fits_sig[it], "SQ")
                current_cost, lambda_reg, current_chi2 = calculate_cost(
                    fits_sig[it], g_sig, par_vec[4*it:4*(it+1)],
                    fit_num_events, num_params, lambda_reg
                )
                if current_cost < best_overall_cost:
                    best_overall_cost = current_cost
                    best_overall_chi2 = current_chi2
                    best_overall_bin = b
                    best_overall_params = [par_vec[4*it + j] for j in range(num_params)]
                print(f"\nBest overall solution: {best_overall_params}")
                print(f"Best overall objective cost: {best_overall_cost:.5f}")
                print(f"Best overall chi2/ndf: {best_overall_chi2:.5f}")

                for j in range(num_params):
                    par_chi2_vec[4*it + j] = best_overall_chi2

                # Plot the final model fit
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
                    fun_Sig = fun_Sig_L_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], math.radians(th_vec[best_overall_bin]))
                elif sig_name == "T":
                    fun_Sig = fun_Sig_T_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], math.radians(th_vec[best_overall_bin]))
                elif sig_name == "LT":
                    fun_Sig = fun_Sig_LT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], math.radians(th_vec[best_overall_bin]))
                elif sig_name == "TT":
                    fun_Sig = fun_Sig_TT_wrapper(g_vec[best_overall_bin], q2_vec[best_overall_bin], w_vec[best_overall_bin], math.radians(th_vec[best_overall_bin]))
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
                latex.DrawLatex(0.35, 0.85, f"Best #chi^{{2}}/ndf: {best_overall_chi2:.3f}")
                c2.Update()
                print("\n")
            c2.Update()
            c3.Update()
            c4.Update()
            c5.Update()
            c6.Update()
            c7.Update()
            c8.Update()
            c9.Update()
            c10.Update()
        params_used.append(param_str)
        # Print all canvases to the output PDF
        c2.Print(outputpdf+'(')
        c3.Print(outputpdf)
        c4.Print(outputpdf)
        c5.Print(outputpdf)
        c6.Print(outputpdf)
        c7.Print(outputpdf)
        c8.Print(outputpdf)
        c9.Print(outputpdf)
        c10.Print(outputpdf+')')
    print(f"\nFits saved to {outputpdf}...")

    return params_used
