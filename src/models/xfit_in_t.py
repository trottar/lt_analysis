#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-07-19 22:33:48 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import random
import ROOT
from ROOT import TFile, TNtuple, TText
from ROOT import TGraph, TGraphErrors, TCanvas
from ROOT import TF1, TFitResultPtr
from multiprocessing import Process, Queue
import numpy as np
import math
import time
import gc
import os, sys

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import adaptive_cooling, simulated_annealing, acceptance_probability, adjust_params, local_search

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
##################################################################################################################################################

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

DEBUG=False
#DEBUG=True

###############################################################################################################################################
# Import separated xsects models

from xfit_active import fun_Sig_L, fun_Sig_T, fun_Sig_LT, fun_Sig_TT, set_val

###############################################################################################################################################

def x_fit_in_t(ParticleType, pol_str, closest_date, Q2, W, inpDict):

    tmin_range = inpDict["tmin"]
    tmax_range = inpDict["tmax"]
    Q2min_range = inpDict["Q2min"]
    Q2max_range = inpDict["Q2max"]

    iter_num = inpDict["iter_num"]
    
    ##############
    # HARD CODED #
    ##############
    # Maximum iterations before ending loop
    max_iterations = 100
    #max_iterations = 500
    #max_iterations = 1000
    #max_iterations = 10000
    ##############
    ##############
    ##############
    
    single_setting(ParticleType, pol_str, closest_date, Q2, W, tmin_range, tmax_range, Q2min_range, Q2max_range, iter_num, max_iterations)
    
def single_setting(ParticleType, pol_str, dir_iter, q2_set, w_set, tmin_range, tmax_range, Q2min_range, Q2max_range, iter_num, max_iterations):

    # Set pol_str, q2_set for xfit_active script
    set_val(pol_str, q2_set)
    
    outputpdf  = "{}/{}_xfit_in_t_Q{}W{}.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
    
    prv_par_vec = []
    g_vec = []
    w_vec = []
    q2_vec = []
    th_vec = []
    par_vec = []
    par_err_vec = []
    par_chi2_vec = []

    l0, l1, l2, l3 = 0, 0, 0, 0
    t0, t1, t2, t3 = 0, 0, 0, 0
    lt0, lt1, lt2, lt3 = 0, 0, 0, 0
    tt0, tt1, tt2, tt3 = 0, 0, 0, 0

    fn_sep = "{}/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(LTANAPATH, ParticleType, pol_str, q2_set.replace("p",""), w_set.replace("p",""))
    nsep = TNtuple("nsep", "nsep", "sigl:sigl_e:sigt:sigt_e:siglt:siglt_e:sigtt:sigtt_e:chi:t:w:q2:thetacm")
    nsep.ReadFile(fn_sep)

    '''
    print("Reading {}...".format(fn_sep))
    for entry in nsep:
        print("sigl: {}, sigl_e: {}, sigt: {}, sigt_e: {}, siglt: {}, siglt_e: {}, sigtt: {}, sigtt_e: {}, chi: {}, t: {}, w: {}, q2: {}, thetacm: {}".format(
            entry.sigl, entry.sigl_e, entry.sigt, entry.sigt_e, entry.siglt, entry.siglt_e, entry.sigtt, entry.sigtt_e, entry.chi, entry.t, entry.w, entry.q2, entry.thetacm
        ))
    '''

    prv_par_vec = []
    para_file_in =  "{}/{}/{}/Q{}W{}/{}/parameters/par.{}_Q{}W{}.dat".format(CACHEPATH, USER, ParticleType, q2_set, w_set, dir_iter, \
                                                                             pol_str, q2_set.replace("p",""), w_set.replace("p",""))
    print("Reading {}...".format(para_file_in))
    try:
        with open(para_file_in, 'r') as f:
            for line in f:
                data = line.split()
                par, par_err, indx, chi2 = map(float, data)
                print("  {} {} {} {}".format(par, par_err, indx, chi2))
                prv_par_vec.append(par)
    except FileNotFoundError:
        print("File {} not found.".format(para_file_in))

    #if iter_num > 2:
    l0, l1, l2, l3, t0, t1, t2, t3, lt0, lt1, lt2, lt3, tt0, tt1, tt2, tt3 = prv_par_vec[:16]
    
    ave_file_in = "{}/src/{}/averages/avek.Q{}W{}.dat".format(LTANAPATH, ParticleType, q2_set.replace("p",""), w_set.replace("p",""))
    with open(ave_file_in, 'r') as f:
        for line in f:
            w, w_e, q2, q2_e, tt, tt_e, thetacm, it = map(float, line.strip().split())

            if pol_str == "pl":
                g = (1 / ((w**2) - (m_p**2))**2)
            else:
                g = (1 / ((w**2) - (m_n**2))**2)
            g_vec.append(g)
            w_vec.append(w)
            q2_vec.append(q2)
            th_vec.append(thetacm)

    c1 = TCanvas("c1", "c1", 800, 800)
    c1.Divide(2, 2)

    c2 = TCanvas("c2", "c2", 800, 800)
    c2.Divide(2, 2)

    # Create ROOT canvases for additional parameter convergence plots
    c3 = TCanvas("c3", "Parameter Convergence", 800, 800)
    c3.Divide(2, 2)
    c4 = TCanvas("c4", "Chi-Square Convergence", 800, 800)
    c4.Divide(2, 2)
    c5 = TCanvas("c5", "Temperature", 800, 800)
    c5.Divide(2, 2)
    c6 = TCanvas("c6", "Acceptance Probability", 800, 800)
    c6.Divide(2, 2)


    # 2 params
    ########
    # SigL #
    ########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig L")
    print("/*--------------------------------------------------*/")

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigL_history = {'p1': [], 'p2': []}

    # Create TGraphs for parameter convergence
    graph_sigL_p1 = TGraph()
    graph_sigL_p2 = TGraph()
    graph_sigL_chi2 = TGraph()
    graph_sigL_temp = TGraph()
    graph_sigL_accept = TGraph()

    c1.cd(1).SetLeftMargin(0.12)
    nsep.Draw("sigl:t:sigl_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigl_0 = l0
        par_sigl_1 = l1
        par_sigl_err_0 = 0.0
        par_sigl_err_1 = 0.0

        # Track the best solution
        best_params = [par_sigl_0, par_sigl_1]
        best_cost = float('inf')
        best_errors = [par_sigl_err_0, par_sigl_err_1]
        previous_params = best_params[:]

        # Check for local minima
        local_minima = []
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigl_prv = TGraph()
            g_sigl_fit = TGraphErrors()
            g_sigl_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigl_0, temperature),
                    simulated_annealing(par_sigl_1, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [l0, l1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigl = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigl.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigl.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigl_X_fit = g_sigl.GetY()[i]
                    sigl_X_fit_err = g_sigl.GetEY()[i]

                    g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
                    g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

                f_sigL = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 2)
                f_sigL.SetParNames("p1", "p2")
                f_sigL.SetParameter(0, current_params[0])
                f_sigL.SetParameter(1, current_params[1])
                #f_sigL.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigl_0), current_params[0] + abs(current_params[0] * par_sigl_0))
                #f_sigL.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigl_1), current_params[1] + abs(current_params[1] * par_sigl_1))

                g_q2_sigl_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
                    g_q2_sigl_fit.SetPointError(i, 0.0, sigl_X_fit_err)
                    sigl_X = (f_sigL.Eval(g_sigl.GetX()[i])) * (g_vec[i])
                    g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)

                r_sigl_fit = g_sigl_fit.Fit(f_sigL, "SQ")

                #f_sigL_status = (r_sigl_fit.Status() == 0 and r_sigl_fit.IsValid())
                f_sigL_status = f_sigL.GetNDF() != 0

                params_sigL_history['p1'].append(current_params[0])
                params_sigL_history['p2'].append(current_params[1])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigL.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigL.GetParameter(0),
                    f_sigL.GetParameter(1)
                ]

                current_errors = [
                    f_sigL.GetParError(0),
                    f_sigL.GetParError(1)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigL_p1.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigL_p2.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigL_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigL_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigL, 2)
                    par_sigl_0, par_sigl_1 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigL_history['p1']) >= max_unchanged_iterations  and \
                   len(params_sigL_history['p2']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigL_history['p1'][-2], 3), round(params_sigL_history['p1'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigL_history['p2'][-2], 3), round(params_sigL_history['p2'][-1], 3), atol=5.0):
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
                    current_params = [l0, l1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]

                # Update parameters with the best found so far
                par_sigl_0, par_sigl_1 = best_params
                par_sigl_err_0, par_sigl_err_1 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                graph_sigL_temp.SetPoint(total_iteration, total_iteration, temperature)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p1={:.3e}, p2={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1]))

                    current_params = adjust_params(best_params)
                    par_sigl_0, par_sigl_1 = current_params
                    par_sigl_err_0, par_sigl_err_1 = [0.0 for _ in range(2)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                # Adjust parameter limits within a random number
                par_sigl_0, par_sigl_1 = [l0, l1]
                par_sigl_err_0, par_sigl_err_1 = [0.0 for _ in range(2)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigl_prv = TGraph()
    g_sigl_fit = TGraphErrors()
    g_sigl_fit_tot = TGraph()        

    f_sigL_pre = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 2)
    f_sigL_pre.SetParNames("p1", "p2")
    f_sigL_pre.FixParameter(0, best_overall_params[0])
    f_sigL_pre.FixParameter(1, best_overall_params[1])

    g_sigl = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigl.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigl.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigl_X_pre = (f_sigL_pre.Eval(g_sigl.GetX()[i])) * (g_vec[i])
        g_sigl_prv.SetPoint(i, g_sigl.GetX()[i], sigl_X_pre)

        sigl_X_fit = g_sigl.GetY()[i]
        sigl_X_fit_err = g_sigl.GetEY()[i]

        g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

    g_sigl.SetTitle("Sig L")
    g_sigl.SetMarkerStyle(5)
    g_sigl.Draw("AP")
    g_sigl.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigl.GetXaxis().CenterTitle()
    g_sigl.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
    g_sigl.GetYaxis().SetTitleOffset(1.5)
    g_sigl.GetYaxis().SetTitleSize(0.035)
    g_sigl.GetYaxis().CenterTitle()

    g_sigl_prv.SetMarkerColor(4)
    g_sigl_prv.SetMarkerStyle(25)
    g_sigl_prv.Draw("P")

    c2.cd(1).SetLeftMargin(0.12)
    g_sigl_fit.SetTitle("Sigma L Model Fit")
    g_sigl_fit.Draw("A*")

    g_sigl_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigl_fit.GetXaxis().CenterTitle()
    g_sigl_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
    g_sigl_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigl_fit.GetYaxis().SetTitleSize(0.035)
    g_sigl_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigl_fit.GetX())
    x_max = max(g_sigl_fit.GetX())
    y_min = min(g_sigl_fit.GetY())
    y_max = max(g_sigl_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigl_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigl_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigL = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 2)
    f_sigL.SetParNames("p1", "p2")
    f_sigL.FixParameter(0, best_overall_params[0])
    f_sigL.FixParameter(1, best_overall_params[1])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigL.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigl_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigl_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_q2_sigl_fit.SetPointError(i, 0.0, sigl_X_fit_err)
        sigl_X = (f_sigL.Eval(g_sigl.GetX()[i])) * (g_vec[i])
        g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)

    r_sigl_fit = g_sigl_fit.Fit(f_sigL, "SQ")
    f_sigL.Draw("same")

    #f_sigL_status = (r_sigl_fit.Status() == 0 and r_sigl_fit.IsValid())
    f_sigL_status = f_sigL.GetNDF() != 0
    f_sigL_status_message = "Fit Successful" if f_sigL_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigL_status_message))

    c1.cd(1)
    g_sigl_fit_tot.SetMarkerStyle(26)
    g_sigl_fit_tot.SetMarkerColor(2)
    g_sigl_fit_tot.SetLineColor(2)
    g_sigl_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigL_y = float('inf')
    max_sigL_y = float('-inf')

    # Update min_sigL_y and max_sigL_y based on each graph's values
    for graph in [graph_sigL_p1, graph_sigL_p2]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigL_y:
                min_sigL_y = y
            if y > max_sigL_y:
                max_sigL_y = y

    # Scale the y-axis
    graph_sigL_p1.SetMinimum(min_sigL_y * 0.9)
    graph_sigL_p1.SetMaximum(max_sigL_y * 1.1)    

    # Plot parameter convergence
    c3.cd(1).SetLeftMargin(0.12)
    graph_sigL_p1.SetTitle("Sig L Parameter Convergence;Optimization Run;Parameter")
    graph_sigL_p1.SetLineColor(ROOT.kRed)
    graph_sigL_p2.SetLineColor(ROOT.kBlue)
    graph_sigL_p1.Draw("ALP")
    graph_sigL_p2.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(1).SetLeftMargin(0.12)
    graph_sigL_chi2.SetTitle("Sig L Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigL_chi2.SetLineColor(ROOT.kBlack)
    graph_sigL_chi2.Draw("ALP")

    # Plot temperature
    c5.cd(1).SetLeftMargin(0.12)
    graph_sigL_temp.SetTitle("Sig L Temperature Convergence;Optimization Run;Temperature")
    graph_sigL_temp.SetLineColor(ROOT.kBlack)
    graph_sigL_temp.Draw("ALP")

    # Plot acceptance probability
    c6.cd(1).SetLeftMargin(0.12)
    graph_sigL_accept.SetTitle("Sig L Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigL_accept.SetLineColor(ROOT.kBlack)
    graph_sigL_accept.Draw("ALP")

    print("\n")    

    '''   
    # 3 params
    ########
    # SigL #
    ########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig L")
    print("/*--------------------------------------------------*/")

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigL_history = {'p1': [], 'p2': [], 'p3': []}

    # Create TGraphs for parameter convergence
    graph_sigL_p1 = TGraph()
    graph_sigL_p2 = TGraph()
    graph_sigL_p3 = TGraph()
    graph_sigL_chi2 = TGraph()
    graph_sigL_temp = TGraph()
    graph_sigL_accept = TGraph()

    c1.cd(1).SetLeftMargin(0.12)
    nsep.Draw("sigl:t:sigl_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigl_0 = l0
        par_sigl_1 = l1
        par_sigl_2 = l2
        par_sigl_err_0 = 0.0
        par_sigl_err_1 = 0.0
        par_sigl_err_2 = 0.0

        # Track the best solution
        best_params = [par_sigl_0, par_sigl_1, par_sigl_2]
        best_cost = float('inf')
        best_errors = [par_sigl_err_0, par_sigl_err_1, par_sigl_err_2]
        previous_params = best_params[:]

        # Check for local minima
        local_minima = []
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigl_prv = TGraph()
            g_sigl_fit = TGraphErrors()
            g_sigl_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:

                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigl_0, temperature),
                    simulated_annealing(par_sigl_1, temperature),
                    simulated_annealing(par_sigl_2, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [l0, l1, l2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigl = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigl.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigl.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigl_X_fit = g_sigl.GetY()[i]
                    sigl_X_fit_err = g_sigl.GetEY()[i]

                    g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
                    g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

                f_sigL = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 3)
                f_sigL.SetParNames("p1", "p2", "p3")
                f_sigL.SetParameter(0, current_params[0])
                f_sigL.SetParameter(1, current_params[1])
                f_sigL.SetParameter(2, current_params[2])
                #f_sigL.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigl_0), current_params[0] + abs(current_params[0] * par_sigl_0))
                #f_sigL.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigl_1), current_params[1] + abs(current_params[1] * par_sigl_1))
                #f_sigL.SetParLimits(2, current_params[2] - abs(current_params[2] * par_sigl_2), current_params[2] + abs(current_params[2] * par_sigl_2))

                g_q2_sigl_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
                    g_q2_sigl_fit.SetPointError(i, 0.0, sigl_X_fit_err)
                    sigl_X = (f_sigL.Eval(g_sigl.GetX()[i])) * (g_vec[i])
                    g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)

                r_sigl_fit = g_sigl_fit.Fit(f_sigL, "SQ")

                #f_sigL_status = (r_sigl_fit.Status() == 0 and r_sigl_fit.IsValid())
                f_sigL_status = f_sigL.GetNDF() != 0

                params_sigL_history['p1'].append(current_params[0])
                params_sigL_history['p2'].append(current_params[1])
                params_sigL_history['p3'].append(current_params[2])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigL.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigL.GetParameter(0),
                    f_sigL.GetParameter(1),
                    f_sigL.GetParameter(2)
                ]

                current_errors = [
                    f_sigL.GetParError(0),
                    f_sigL.GetParError(1),
                    f_sigL.GetParError(2)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigL_p1.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigL_p2.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigL_p3.SetPoint(total_iteration, total_iteration, current_params[2])
                graph_sigL_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigL_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigL_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigL, 3)
                    par_sigl_0, par_sigl_1, par_sigl_2 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigL_history['p1']) >= max_unchanged_iterations  and \
                   len(params_sigL_history['p2']) >= max_unchanged_iterations  and \
                   len(params_sigL_history['p3']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigL_history['p1'][-2], 3), round(params_sigL_history['p1'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigL_history['p2'][-2], 3), round(params_sigL_history['p2'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigL_history['p3'][-2], 3), round(params_sigL_history['p3'][-1], 3), atol=5.0):
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
                    current_params = [l0, l1, l2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]

                # Update parameters with the best found so far
                par_sigl_0, par_sigl_1, par_sigl_2 = best_params
                par_sigl_err_0, par_sigl_err_1, par_sigl_err_2 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1], current_params[2]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p1={:.3e}, p2={:.3e}, p3={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2]))

                    current_params = adjust_params(best_params)
                    par_sigl_0, par_sigl_1, par_sigl_2 = current_params
                    par_sigl_err_0, par_sigl_err_1, par_sigl_err_2 = [0.0 for _ in range(3)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                # Adjust parameter limits within a random number
                par_sigl_0, par_sigl_1, par_sigl_2 = [l0, l1, l2]
                par_sigl_err_0, par_sigl_err_1, par_sigl_err_2 = [0.0 for _ in range(3)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigl_prv = TGraph()
    g_sigl_fit = TGraphErrors()
    g_sigl_fit_tot = TGraph()        

    f_sigL_pre = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 3)
    f_sigL_pre.SetParNames("p1", "p2", "p3")
    f_sigL_pre.FixParameter(0, best_overall_params[0])
    f_sigL_pre.FixParameter(1, best_overall_params[1])
    f_sigL_pre.FixParameter(2, best_overall_params[2])

    g_sigl = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigl.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigl.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigl_X_pre = (f_sigL_pre.Eval(g_sigl.GetX()[i])) * (g_vec[i])
        g_sigl_prv.SetPoint(i, g_sigl.GetX()[i], sigl_X_pre)

        sigl_X_fit = g_sigl.GetY()[i]
        sigl_X_fit_err = g_sigl.GetEY()[i]

        g_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_sigl_fit.SetPointError(i, 0, sigl_X_fit_err)

    g_sigl.SetTitle("Sig L")
    g_sigl.SetMarkerStyle(5)
    g_sigl.Draw("AP")
    g_sigl.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigl.GetXaxis().CenterTitle()
    g_sigl.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
    g_sigl.GetYaxis().SetTitleOffset(1.5)
    g_sigl.GetYaxis().SetTitleSize(0.035)
    g_sigl.GetYaxis().CenterTitle()

    g_sigl_prv.SetMarkerColor(4)
    g_sigl_prv.SetMarkerStyle(25)
    g_sigl_prv.Draw("P")

    c2.cd(1).SetLeftMargin(0.12)
    g_sigl_fit.SetTitle("Sigma L Model Fit")
    g_sigl_fit.Draw("A*")

    g_sigl_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigl_fit.GetXaxis().CenterTitle()
    g_sigl_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
    g_sigl_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigl_fit.GetYaxis().SetTitleSize(0.035)
    g_sigl_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigl_fit.GetX())
    x_max = max(g_sigl_fit.GetX())
    y_min = min(g_sigl_fit.GetY())
    y_max = max(g_sigl_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigl_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigl_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigL = TF1("sig_L", fun_Sig_L, tmin_range, tmax_range, 3)
    f_sigL.SetParNames("p1", "p2", "p3")
    f_sigL.FixParameter(0, best_overall_params[0])
    f_sigL.FixParameter(1, best_overall_params[1])
    f_sigL.FixParameter(2, best_overall_params[2])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigL.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigl_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigl_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigl_fit.SetPoint(i, g_sigl.GetX()[i], sigl_X_fit)
        g_q2_sigl_fit.SetPointError(i, 0.0, sigl_X_fit_err)
        sigl_X = (f_sigL.Eval(g_sigl.GetX()[i])) * (g_vec[i])
        g_sigl_fit_tot.SetPoint(i, g_sigl.GetX()[i], sigl_X)

    r_sigl_fit = g_sigl_fit.Fit(f_sigL, "SQ")
    f_sigL.Draw("same")

    #f_sigL_status = (r_sigl_fit.Status() == 0 and r_sigl_fit.IsValid())
    f_sigL_status = f_sigL.GetNDF() != 0
    f_sigL_status_message = "Fit Successful" if f_sigL_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigL_status_message))

    c1.cd(1)
    g_sigl_fit_tot.SetMarkerStyle(26)
    g_sigl_fit_tot.SetMarkerColor(2)
    g_sigl_fit_tot.SetLineColor(2)
    g_sigl_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigL_y = float('inf')
    max_sigL_y = float('-inf')

    # Update min_sigL_y and max_sigL_y based on each graph's values
    for graph in [graph_sigL_p1, graph_sigL_p2, graph_sigL_p3]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigL_y:
                min_sigL_y = y
            if y > max_sigL_y:
                max_sigL_y = y

    # Scale the y-axis
    graph_sigL_p1.SetMinimum(min_sigL_y * 0.9)
    graph_sigL_p1.SetMaximum(max_sigL_y * 1.1)    

    # Plot parameter convergence
    c3.cd(1).SetLeftMargin(0.12)
    graph_sigL_p1.SetTitle("Sig L Parameter Convergence;Optimization Run;Parameter")
    graph_sigL_p1.SetLineColor(ROOT.kRed)
    graph_sigL_p2.SetLineColor(ROOT.kBlue)
    graph_sigL_p3.SetLineColor(ROOT.kGreen)
    graph_sigL_p1.Draw("ALP")
    graph_sigL_p2.Draw("LP SAME")
    graph_sigL_p3.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(1).SetLeftMargin(0.12)
    graph_sigL_chi2.SetTitle("Sig L Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigL_chi2.SetLineColor(ROOT.kBlack)
    graph_sigL_chi2.Draw("ALP")

    # Plot temperature
    c5.cd(1).SetLeftMargin(0.12)
    graph_sigL_temp.SetTitle("Sig L Temperature Convergence;Optimization Run;Temperature")
    graph_sigL_temp.SetLineColor(ROOT.kBlack)
    graph_sigL_temp.Draw("ALP")

    # Plot acceptance probability
    c6.cd(1).SetLeftMargin(0.12)
    graph_sigL_accept.SetTitle("Sig L Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigL_accept.SetLineColor(ROOT.kBlack)
    graph_sigL_accept.Draw("ALP")

    print("\n")    
    '''
    
    '''
    # 2 params
    ########
    # SigT #
    ########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig T")
    print("/*--------------------------------------------------*/")

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigT_history = {'p5': [], 'p6': []}

    # Create TGraphs for parameter convergence
    graph_sigT_p5 = TGraph()
    graph_sigT_p6 = TGraph()
    graph_sigT_chi2 = TGraph()
    graph_sigT_temp = TGraph()
    graph_sigT_accept = TGraph()

    c1.cd(2).SetLeftMargin(0.12)
    nsep.Draw("sigt:t:sigt_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigt_0 = t0
        par_sigt_1 = t1
        par_sigt_err_0 = 0.0
        par_sigt_err_1 = 0.0

        # Track the best solution
        best_params = [par_sigt_0, par_sigt_1]
        best_cost = float('inf')
        previous_params = best_params[:]
        best_errors = [par_sigt_err_0, par_sigt_err_1]

        # Check for local minima
        local_minima = []
        local_iterations = 0
        tabu_list = set()    

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigt_prv = TGraph()
            g_sigt_fit = TGraphErrors()
            g_sigt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigt_0, temperature),
                    simulated_annealing(par_sigt_1, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [t0, t1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigt_X_fit = g_sigt.GetY()[i]
                    sigt_X_fit_err = g_sigt.GetEY()[i]

                    g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

                f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 2)
                f_sigT.SetParNames("p5", "p6")
                f_sigT.SetParameter(0, current_params[0])
                f_sigT.SetParameter(1, current_params[1])
                #f_sigT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigt_0), current_params[0] + abs(current_params[0] * par_sigt_0))
                #f_sigT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigt_1), current_params[1] + abs(current_params[1] * par_sigt_1))

                g_q2_sigt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
                    sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
                    g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

                r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")

                #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
                f_sigT_status = f_sigT.GetNDF() != 0

                params_sigT_history['p5'].append(current_params[0])
                params_sigT_history['p6'].append(current_params[1])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigT.GetParameter(0),
                    f_sigT.GetParameter(1)
                ]

                current_errors = [
                    f_sigT.GetParError(0),
                    f_sigT.GetParError(1)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigT_p5.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigT_p6.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigT, 2)
                    par_sigt_0, par_sigt_1 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigT_history['p5']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p6']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigT_history['p5'][-2], 3), round(params_sigT_history['p5'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p6'][-2], 3), round(params_sigT_history['p6'][-1], 3), atol=5.0):
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
                    current_params = [t0, t1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]                

                # Update parameters with the best found so far
                par_sigt_0, par_sigt_1 = best_params
                par_sigt_err_0, par_sigt_err_1 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p5={:.3e}, p6={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1]))

                    current_params = adjust_params(best_params)
                    par_sigt_0, par_sigt_1 = current_params
                    par_sigt_err_0, par_sigt_err_1 = [0.0 for _ in range(2)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                par_sigt_0, par_sigt_1 = [t0, t1]
                par_sigt_err_0, par_sigt_err_1 = [0.0 for _ in range(2)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigt_prv = TGraph()
    g_sigt_fit = TGraphErrors()
    g_sigt_fit_tot = TGraph()    

    f_sigT_pre = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 2)
    f_sigT_pre.SetParNames("p5", "p6")
    f_sigT_pre.FixParameter(0, best_overall_params[0])
    f_sigT_pre.FixParameter(1, best_overall_params[1])

    g_sigt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigt_X_pre = (f_sigT_pre.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_prv.SetPoint(i, g_sigt.GetX()[i], sigt_X_pre)

        sigt_X_fit = g_sigt.GetY()[i]
        sigt_X_fit_err = g_sigt.GetEY()[i]

        g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

    g_sigt.SetTitle("Sig T")
    g_sigt.SetMarkerStyle(5)
    g_sigt.Draw("AP")
    g_sigt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt.GetXaxis().CenterTitle()
    g_sigt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt.GetYaxis().SetTitleOffset(1.5)
    g_sigt.GetYaxis().SetTitleSize(0.035)
    g_sigt.GetYaxis().CenterTitle()

    g_sigt_prv.SetMarkerColor(4)
    g_sigt_prv.SetMarkerStyle(25)
    g_sigt_prv.Draw("P")

    c2.cd(2).SetLeftMargin(0.12)
    g_sigt_fit.SetTitle("Sigma T Model Fit")
    g_sigt_fit.Draw("A*")

    g_sigt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt_fit.GetXaxis().CenterTitle()
    g_sigt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigt_fit.GetYaxis().SetTitleSize(0.035)
    g_sigt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigt_fit.GetX())
    x_max = max(g_sigt_fit.GetX())
    y_min = min(g_sigt_fit.GetY())
    y_max = max(g_sigt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 2)
    f_sigT.SetParNames("p5", "p6")
    f_sigT.FixParameter(0, best_overall_params[0])
    f_sigT.FixParameter(1, best_overall_params[1])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
        sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

    r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")
    f_sigT.Draw("same")

    #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
    f_sigT_status = f_sigT.GetNDF() != 0
    f_sigT_status_message = "Fit Successful" if f_sigT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigT_status_message))

    c1.cd(2)
    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigT_y = float('inf')
    max_sigT_y = float('-inf')

    # Update min_sigT_y and max_sigT_y based on each graph's values
    for graph in [graph_sigT_p5, graph_sigT_p6]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigT_y:
                min_sigT_y = y
            if y > max_sigT_y:
                max_sigT_y = y

    # Scale the y-axis
    graph_sigT_p5.SetMinimum(min_sigT_y * 0.9)
    graph_sigT_p6.SetMaximum(max_sigT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(2).SetLeftMargin(0.12)
    graph_sigT_p5.SetTitle("Sig T Parameter Convergence;Optimization Run;Parameter")
    graph_sigT_p5.SetLineColor(ROOT.kRed)
    graph_sigT_p6.SetLineColor(ROOT.kBlue)
    graph_sigT_p5.Draw("ALP")
    graph_sigT_p6.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(2).SetLeftMargin(0.12)
    graph_sigT_chi2.SetTitle("Sig T Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigT_chi2.Draw("ALP")

    # Plot temperature convergence
    c5.cd(2).SetLeftMargin(0.12)
    graph_sigT_temp.SetTitle("Sig T Temperature Convergence;Optimization Run;Temperature")
    graph_sigT_temp.SetLineColor(ROOT.kBlack)
    graph_sigT_temp.Draw("ALP")

    # Plot acceptance probability convergence
    c6.cd(2).SetLeftMargin(0.12)
    graph_sigT_accept.SetTitle("Sig T Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigT_accept.SetLineColor(ROOT.kBlack)
    graph_sigT_accept.Draw("ALP")

    print("\n")    
    '''

    '''
    # 3 params
    ########
    # SigT #
    ########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig T")
    print("/*--------------------------------------------------*/")

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigT_history = {'p5': [], 'p6': [], 'p7': []}

    # Create TGraphs for parameter convergence
    graph_sigT_p5 = TGraph()
    graph_sigT_p6 = TGraph()
    graph_sigT_p7 = TGraph()
    graph_sigT_chi2 = TGraph()
    graph_sigT_temp = TGraph()
    graph_sigT_accept = TGraph()

    c1.cd(2).SetLeftMargin(0.12)
    nsep.Draw("sigt:t:sigt_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigt_0 = t0
        par_sigt_1 = t1
        par_sigt_2 = t2
        par_sigt_err_0 = 0.0
        par_sigt_err_1 = 0.0
        par_sigt_err_2 = 0.0

        # Track the best solution
        best_params = [par_sigt_0, par_sigt_1, par_sigt_2]
        best_cost = float('inf')
        best_errors = [par_sigt_err_0, par_sigt_err_1, par_sigt_err_2]
        previous_params = best_params[:]

        # Check for local minima
        local_minima = []
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigt_prv = TGraph()
            g_sigt_fit = TGraphErrors()
            g_sigt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigt_0, temperature),
                    simulated_annealing(par_sigt_1, temperature),
                    simulated_annealing(par_sigt_2, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [t0, t1, t2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigt_X_fit = g_sigt.GetY()[i]
                    sigt_X_fit_err = g_sigt.GetEY()[i]

                    g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

                f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 3)
                f_sigT.SetParNames("p5", "p6", "p7")
                f_sigT.SetParameter(0, current_params[0])
                f_sigT.SetParameter(1, current_params[1])
                f_sigT.SetParameter(2, current_params[2])
                #f_sigT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigt_0), current_params[0] + abs(current_params[0] * par_sigt_0))
                #f_sigT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigt_1), current_params[1] + abs(current_params[1] * par_sigt_1))
                #f_sigT.SetParLimits(2, current_params[2] - abs(current_params[2] * par_sigt_2), current_params[2] + abs(current_params[2] * par_sigt_2))

                g_q2_sigt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
                    sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
                    g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

                r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")

                #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
                f_sigT_status = f_sigT.GetNDF() != 0

                params_sigT_history['p5'].append(current_params[0])
                params_sigT_history['p6'].append(current_params[1])
                params_sigT_history['p7'].append(current_params[2])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigT.GetParameter(0),
                    f_sigT.GetParameter(1),
                    f_sigT.GetParameter(2)
                ]

                current_errors = [
                    f_sigT.GetParError(0),
                    f_sigT.GetParError(1),
                    f_sigT.GetParError(2)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigT_p5.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigT_p6.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigT_p7.SetPoint(total_iteration, total_iteration, current_params[2])
                graph_sigT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigT, 3)
                    par_sigt_0, par_sigt_1, par_sigt_2 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigT_history['p5']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p6']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p7']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigT_history['p5'][-2], 3), round(params_sigT_history['p5'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p6'][-2], 3), round(params_sigT_history['p6'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p7'][-2], 3), round(params_sigT_history['p7'][-1], 3), atol=5.0):
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
                    current_params = [t0, t1, t2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]

                # Update parameter limits
                par_sigt_0, par_sigt_1, par_sigt_2 = best_params
                par_sigt_err_0, par_sigt_err_1, par_sigt_err_2 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1], current_params[2]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p5={:.3e}, p6={:.3e}, p7={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2]))

                    current_params = adjust_params(best_params)
                    par_sigt_0, par_sigt_1, par_sigt_2 = current_params
                    par_sigt_err_0, par_sigt_err_1, par_sigt_err_2 = [0.0 for _ in range(3)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                # Adjust parameter limits within a random number
                par_sigt_0, par_sigt_1, par_sigt_2 = [t0, t1, t2]
                par_sigt_err_0, par_sigt_err_1, par_sigt_err_2 = [0.0 for _ in range(3)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigt_prv = TGraph()
    g_sigt_fit = TGraphErrors()
    g_sigt_fit_tot = TGraph()        

    f_sigT_pre = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 3)
    f_sigT_pre.SetParNames("p5", "p6", "p7")
    f_sigT_pre.FixParameter(0, best_overall_params[0])
    f_sigT_pre.FixParameter(1, best_overall_params[1])
    f_sigT_pre.FixParameter(2, best_overall_params[2])

    g_sigt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigt_X_pre = (f_sigT_pre.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_prv.SetPoint(i, g_sigt.GetX()[i], sigt_X_pre)

        sigt_X_fit = g_sigt.GetY()[i]
        sigt_X_fit_err = g_sigt.GetEY()[i]

        g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

    g_sigt.SetTitle("Sig T")
    g_sigt.SetMarkerStyle(5)
    g_sigt.Draw("AP")
    g_sigt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt.GetXaxis().CenterTitle()
    g_sigt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt.GetYaxis().SetTitleOffset(1.5)
    g_sigt.GetYaxis().SetTitleSize(0.035)
    g_sigt.GetYaxis().CenterTitle()

    g_sigt_prv.SetMarkerColor(4)
    g_sigt_prv.SetMarkerStyle(25)
    g_sigt_prv.Draw("P")

    c2.cd(2).SetLeftMargin(0.12)
    g_sigt_fit.SetTitle("Sigma T Model Fit")
    g_sigt_fit.Draw("A*")

    g_sigt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt_fit.GetXaxis().CenterTitle()
    g_sigt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigt_fit.GetYaxis().SetTitleSize(0.035)
    g_sigt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigt_fit.GetX())
    x_max = max(g_sigt_fit.GetX())
    y_min = min(g_sigt_fit.GetY())
    y_max = max(g_sigt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 3)
    f_sigT.SetParNames("p5", "p6", "p7")
    f_sigT.FixParameter(0, best_overall_params[0])
    f_sigT.FixParameter(1, best_overall_params[1])
    f_sigT.FixParameter(2, best_overall_params[2])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
        sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

    r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")
    f_sigT.Draw("same")

    #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
    f_sigT_status = f_sigT.GetNDF() != 0
    f_sigT_status_message = "Fit Successful" if f_sigT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigT_status_message))

    c1.cd(2)
    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigT_y = float('inf')
    max_sigT_y = float('-inf')

    # Update min_sigT_y and max_sigT_y based on each graph's values
    for graph in [graph_sigT_p5, graph_sigT_p6, graph_sigT_p7]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigT_y:
                min_sigT_y = y
            if y > max_sigT_y:
                max_sigT_y = y

    # Scale the y-axis
    graph_sigT_p5.SetMinimum(min_sigT_y * 0.9)
    graph_sigT_p5.SetMaximum(max_sigT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(2).SetLeftMargin(0.12)
    graph_sigT_p5.SetTitle("Sig T Parameter Convergence;Optimization Run;Parameter")
    graph_sigT_p5.SetLineColor(ROOT.kRed)
    graph_sigT_p6.SetLineColor(ROOT.kBlue)
    graph_sigT_p7.SetLineColor(ROOT.kGreen)
    graph_sigT_p5.Draw("ALP")
    graph_sigT_p6.Draw("LP SAME")
    graph_sigT_p7.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(2).SetLeftMargin(0.12)
    graph_sigT_chi2.SetTitle("Sig T Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigT_chi2.Draw("ALP")

    # Plot temperature
    c5.cd(2).SetLeftMargin(0.12)
    graph_sigT_temp.SetTitle("Sig T Temperature Convergence;Optimization Run;Temperature")
    graph_sigT_temp.SetLineColor(ROOT.kBlack)
    graph_sigT_temp.Draw("ALP")

    # Plot acceptance probability
    c6.cd(2).SetLeftMargin(0.12)
    graph_sigT_accept.SetTitle("Sig T Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigT_accept.SetLineColor(ROOT.kBlack)
    graph_sigT_accept.Draw("ALP")

    print("\n")    
    '''

    # 4 params
    ########
    # SigT #
    ########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig T")
    print("/*--------------------------------------------------*/")

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigT_history = {'p5': [], 'p6': [], 'p7': [], 'p8': []}

    # Create TGraphs for parameter convergence
    graph_sigT_p5 = TGraph()
    graph_sigT_p6 = TGraph()
    graph_sigT_p7 = TGraph()
    graph_sigT_p8 = TGraph()
    graph_sigT_chi2 = TGraph()
    graph_sigT_temp = TGraph()
    graph_sigT_accept = TGraph()

    c1.cd(2).SetLeftMargin(0.12)
    nsep.Draw("sigt:t:sigt_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigt_0 = t0
        par_sigt_1 = t1
        par_sigt_2 = t2
        par_sigt_3 = t3
        par_sigt_err_0 = 0.0
        par_sigt_err_1 = 0.0
        par_sigt_err_2 = 0.0
        par_sigt_err_3 = 0.0

        # Track the best solution
        best_params = [par_sigt_0, par_sigt_1, par_sigt_2, par_sigt_3]
        best_cost = float('inf')
        best_errors = [par_sigt_err_0, par_sigt_err_1, par_sigt_err_2, par_sigt_err_3]
        previous_params = best_params[:]

        # Check for local minima
        local_minima = []
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigt_prv = TGraph()
            g_sigt_fit = TGraphErrors()
            g_sigt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigt_0, temperature),
                    simulated_annealing(par_sigt_1, temperature),
                    simulated_annealing(par_sigt_2, temperature),
                    simulated_annealing(par_sigt_3, temperature)                    
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [t0, t1, t2, t3]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigt_X_fit = g_sigt.GetY()[i]
                    sigt_X_fit_err = g_sigt.GetEY()[i]

                    g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

                f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 4)
                f_sigT.SetParNames("p5", "p6", "p7", "p8")
                f_sigT.SetParameter(0, current_params[0])
                f_sigT.SetParameter(1, current_params[1])
                f_sigT.SetParameter(2, current_params[2])
                f_sigT.SetParameter(2, current_params[3])
                #f_sigT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigt_0), current_params[0] + abs(current_params[0] * par_sigt_0))
                #f_sigT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigt_1), current_params[1] + abs(current_params[1] * par_sigt_1))
                #f_sigT.SetParLimits(2, current_params[2] - abs(current_params[2] * par_sigt_2), current_params[2] + abs(current_params[2] * par_sigt_2))
                #f_sigT.SetParLimits(3, current_params[3] - abs(current_params[3] * par_sigt_3), current_params[3] + abs(current_params[3] * par_sigt_3))

                g_q2_sigt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
                    g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
                    sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
                    g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

                r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")

                #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
                f_sigT_status = f_sigT.GetNDF() != 0

                params_sigT_history['p5'].append(current_params[0])
                params_sigT_history['p6'].append(current_params[1])
                params_sigT_history['p7'].append(current_params[2])
                params_sigT_history['p8'].append(current_params[3])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigT.GetParameter(0),
                    f_sigT.GetParameter(1),
                    f_sigT.GetParameter(2),
                    f_sigT.GetParameter(3)
                ]

                current_errors = [
                    f_sigT.GetParError(0),
                    f_sigT.GetParError(1),
                    f_sigT.GetParError(2),
                    f_sigT.GetParError(3)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigT_p5.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigT_p6.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigT_p7.SetPoint(total_iteration, total_iteration, current_params[2])
                graph_sigT_p8.SetPoint(total_iteration, total_iteration, current_params[3])
                graph_sigT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigT, 4)
                    par_sigt_0, par_sigt_1, par_sigt_2, par_sigt_3 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigT_history['p5']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p6']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p7']) >= max_unchanged_iterations  and \
                   len(params_sigT_history['p8']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigT_history['p5'][-2], 3), round(params_sigT_history['p5'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p6'][-2], 3), round(params_sigT_history['p6'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p7'][-2], 3), round(params_sigT_history['p7'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigT_history['p8'][-2], 3), round(params_sigT_history['p8'][-1], 3), atol=5.0):
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
                    current_params = [t0, t1, t2, t3]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]

                # Update parameter limits
                par_sigt_0, par_sigt_1, par_sigt_2, par_sigt_3 = best_params
                par_sigt_err_0, par_sigt_err_1, par_sigt_err_2, par_sigt_err_3 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1], current_params[2], current_params[3]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p5={:.3e}, p6={:.3e}, p7={:.3e}, p8={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2], current_params[3]))

                    current_params = adjust_params(best_params)
                    par_sigt_0, par_sigt_1, par_sigt_2, par_sigt_3 = current_params
                    par_sigt_err_0, par_sigt_err_1, par_sigt_err_2, par_sigt_err_3 = [0.0 for _ in range(4)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                # Adjust parameter limits within a random number
                par_sigt_0, par_sigt_1, par_sigt_2, par_sigt_3 = [t0, t1, t2, t3]
                par_sigt_err_0, par_sigt_err_1, par_sigt_err_2, par_sigt_err_3 = [0.0 for _ in range(4)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigt_prv = TGraph()
    g_sigt_fit = TGraphErrors()
    g_sigt_fit_tot = TGraph()        

    f_sigT_pre = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 4)
    f_sigT_pre.SetParNames("p5", "p6", "p7", "p8")
    f_sigT_pre.FixParameter(0, best_overall_params[0])
    f_sigT_pre.FixParameter(1, best_overall_params[1])
    f_sigT_pre.FixParameter(2, best_overall_params[2])
    f_sigT_pre.FixParameter(3, best_overall_params[3])

    g_sigt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigt_X_pre = (f_sigT_pre.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_prv.SetPoint(i, g_sigt.GetX()[i], sigt_X_pre)

        sigt_X_fit = g_sigt.GetY()[i]
        sigt_X_fit_err = g_sigt.GetEY()[i]

        g_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_sigt_fit.SetPointError(i, 0, sigt_X_fit_err)

    g_sigt.SetTitle("Sig T")
    g_sigt.SetMarkerStyle(5)
    g_sigt.Draw("AP")
    g_sigt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt.GetXaxis().CenterTitle()
    g_sigt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt.GetYaxis().SetTitleOffset(1.5)
    g_sigt.GetYaxis().SetTitleSize(0.035)
    g_sigt.GetYaxis().CenterTitle()

    g_sigt_prv.SetMarkerColor(4)
    g_sigt_prv.SetMarkerStyle(25)
    g_sigt_prv.Draw("P")

    c2.cd(2).SetLeftMargin(0.12)
    g_sigt_fit.SetTitle("Sigma T Model Fit")
    g_sigt_fit.Draw("A*")

    g_sigt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigt_fit.GetXaxis().CenterTitle()
    g_sigt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
    g_sigt_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigt_fit.GetYaxis().SetTitleSize(0.035)
    g_sigt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigt_fit.GetX())
    x_max = max(g_sigt_fit.GetX())
    y_min = min(g_sigt_fit.GetY())
    y_max = max(g_sigt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigT = TF1("sig_T", fun_Sig_T, tmin_range, tmax_range, 4)
    f_sigT.SetParNames("p5", "p6", "p7")
    f_sigT.FixParameter(0, best_overall_params[0])
    f_sigT.FixParameter(1, best_overall_params[1])
    f_sigT.FixParameter(2, best_overall_params[2])
    f_sigT.FixParameter(3, best_overall_params[3])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigt_fit.SetPoint(i, g_sigt.GetX()[i], sigt_X_fit)
        g_q2_sigt_fit.SetPointError(i, 0.0, sigt_X_fit_err)
        sigt_X = (f_sigT.Eval(g_sigt.GetX()[i])) * (g_vec[i])
        g_sigt_fit_tot.SetPoint(i, g_sigt.GetX()[i], sigt_X)

    r_sigt_fit = g_sigt_fit.Fit(f_sigT, "SQ")
    f_sigT.Draw("same")

    #f_sigT_status = (r_sigt_fit.Status() == 0 and r_sigt_fit.IsValid())
    f_sigT_status = f_sigT.GetNDF() != 0
    f_sigT_status_message = "Fit Successful" if f_sigT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigT_status_message))

    c1.cd(2)
    g_sigt_fit_tot.SetMarkerStyle(26)
    g_sigt_fit_tot.SetMarkerColor(2)
    g_sigt_fit_tot.SetLineColor(2)
    g_sigt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigT_y = float('inf')
    max_sigT_y = float('-inf')

    # Update min_sigT_y and max_sigT_y based on each graph's values
    for graph in [graph_sigT_p5, graph_sigT_p6, graph_sigT_p7, graph_sigT_p8]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigT_y:
                min_sigT_y = y
            if y > max_sigT_y:
                max_sigT_y = y

    # Scale the y-axis
    graph_sigT_p5.SetMinimum(min_sigT_y * 0.9)
    graph_sigT_p5.SetMaximum(max_sigT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(2).SetLeftMargin(0.12)
    graph_sigT_p5.SetTitle("Sig T Parameter Convergence;Optimization Run;Parameter")
    graph_sigT_p5.SetLineColor(ROOT.kRed)
    graph_sigT_p6.SetLineColor(ROOT.kBlue)
    graph_sigT_p7.SetLineColor(ROOT.kGreen)
    graph_sigT_p5.Draw("ALP")
    graph_sigT_p6.Draw("LP SAME")
    graph_sigT_p7.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(2).SetLeftMargin(0.12)
    graph_sigT_chi2.SetTitle("Sig T Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigT_chi2.Draw("ALP")

    # Plot temperature
    c5.cd(2).SetLeftMargin(0.12)
    graph_sigT_temp.SetTitle("Sig T Temperature Convergence;Optimization Run;Temperature")
    graph_sigT_temp.SetLineColor(ROOT.kBlack)
    graph_sigT_temp.Draw("ALP")

    # Plot acceptance probability
    c6.cd(2).SetLeftMargin(0.12)
    graph_sigT_accept.SetTitle("Sig T Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigT_accept.SetLineColor(ROOT.kBlack)
    graph_sigT_accept.Draw("ALP")

    print("\n")    
    
    # 2 params
    #########
    # SigLT #
    #########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig LT")
    print("/*--------------------------------------------------*/")    

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigLT_history = {'p9': [], 'p10': []}

    # Create TGraphs for parameter convergence
    graph_sigLT_p9 = TGraph()
    graph_sigLT_p10 = TGraph()
    graph_sigLT_chi2 = TGraph()
    graph_sigLT_temp = TGraph()
    graph_sigLT_accept = TGraph()

    c1.cd(3).SetLeftMargin(0.12)
    nsep.Draw("siglt:t:siglt_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_siglt_0 = lt0
        par_siglt_1 = lt1
        par_siglt_err_0 = 0.0
        par_siglt_err_1 = 0.0
        par_siglt_err_2 = 0.0

        # Track the best solution
        best_params = [par_siglt_0, par_siglt_1]
        best_cost = float('inf')
        previous_params = best_params[:]
        best_errors = [par_siglt_err_0, par_siglt_err_1]

        # Check for local minima
        local_minima = []
        local_iterations = 0
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_siglt_prv = TGraph()
            g_siglt_fit = TGraphErrors()
            g_siglt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_siglt_0, temperature),
                    simulated_annealing(par_siglt_1, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [lt0, lt1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_siglt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_siglt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_siglt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    siglt_X_fit = g_siglt.GetY()[i]
                    siglt_X_fit_err = g_siglt.GetEY()[i]

                    g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
                    g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

                f_sigLT = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 2)
                f_sigLT.SetParNames("p9", "p10")
                f_sigLT.SetParameter(0, current_params[0])
                f_sigLT.SetParameter(1, current_params[1])
                #f_sigLT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_siglt_0), current_params[0] + abs(current_params[0] * par_siglt_0))
                #f_sigLT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_siglt_1), current_params[1] + abs(current_params[1] * par_siglt_1))
                #f_sigLT.SetParLimits(0, -5, 5)
                #f_sigLT.SetParLimits(1, -5, 5)

                g_q2_siglt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
                    g_q2_siglt_fit.SetPointError(i, 0.0, siglt_X_fit_err)
                    siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
                    g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)

                r_siglt_fit = g_siglt_fit.Fit(f_sigLT, "SQ")

                #f_sigLT_status = (r_siglt_fit.Status() == 0 and r_siglt_fit.IsValid())
                f_sigLT_status = f_sigLT.GetNDF() != 0

                params_sigLT_history['p9'].append(current_params[0])
                params_sigLT_history['p10'].append(current_params[1])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigLT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigLT.GetParameter(0),
                    f_sigLT.GetParameter(1)
                ]

                current_errors = [
                    f_sigLT.GetParError(0),
                    f_sigLT.GetParError(1)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigLT_p9.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigLT_p10.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigLT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigLT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigLT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

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
                    current_params = local_search(current_params, f_sigLT, 3)
                    par_siglt_0, par_siglt_1 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigLT_history['p9']) >= max_unchanged_iterations  and \
                   len(params_sigLT_history['p10']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigLT_history['p9'][-2], 3), round(params_sigLT_history['p9'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigLT_history['p10'][-2], 3), round(params_sigLT_history['p10'][-1], 3), atol=5.0):
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
                    current_params = [lt0, lt1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]                

                # Update parameters with the best found so far
                par_siglt_0, par_siglt_1 = best_params
                par_siglt_err_0, par_siglt_err_1 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p9={:.3e}, p10={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1]))

                    current_params = adjust_params(best_params)
                    par_siglt_0, par_siglt_1 = current_params
                    par_siglt_err_0, par_siglt_err_1 = [0.0 for _ in range(2)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                par_siglt_0, par_siglt_1 = [lt0, lt1]
                par_siglt_err_0, par_siglt_err_1 = [0.0 for _ in range(2)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0                

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_siglt_prv = TGraph()
    g_siglt_fit = TGraphErrors()
    g_siglt_fit_tot = TGraph()    

    f_sigLT_pre = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 2)
    f_sigLT_pre.SetParNames("p9", "p10")
    f_sigLT_pre.FixParameter(0, best_overall_params[0])
    f_sigLT_pre.FixParameter(1, best_overall_params[1])

    g_siglt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_siglt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_siglt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        siglt_X_pre = (f_sigLT_pre.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
        g_siglt_prv.SetPoint(i, g_siglt.GetX()[i], siglt_X_pre)

        siglt_X_fit = g_siglt.GetY()[i]
        siglt_X_fit_err = g_siglt.GetEY()[i]

        g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

    g_siglt.SetTitle("Sig LT")
    g_siglt.SetMarkerStyle(5)
    g_siglt.Draw("AP")
    g_siglt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_siglt.GetXaxis().CenterTitle()
    g_siglt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
    g_siglt.GetYaxis().SetTitleOffset(1.5)
    g_siglt.GetYaxis().SetTitleSize(0.035)
    g_siglt.GetYaxis().CenterTitle()

    g_siglt_prv.SetMarkerColor(4)
    g_siglt_prv.SetMarkerStyle(25)
    g_siglt_prv.Draw("P")

    c2.cd(3).SetLeftMargin(0.12)
    g_siglt_fit.SetTitle("Sigma LT Model Fit")
    g_siglt_fit.Draw("A*")

    g_siglt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_siglt_fit.GetXaxis().CenterTitle()
    g_siglt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
    g_siglt_fit.GetYaxis().SetTitleOffset(1.5)
    g_siglt_fit.GetYaxis().SetTitleSize(0.035)
    g_siglt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_siglt_fit.GetX())
    x_max = max(g_siglt_fit.GetX())
    y_min = min(g_siglt_fit.GetY())
    y_max = max(g_siglt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_siglt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_siglt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigLT = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 2)
    f_sigLT.SetParNames("p9", "p10")
    f_sigLT.FixParameter(0, best_overall_params[0])
    f_sigLT.FixParameter(1, best_overall_params[1])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigLT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_siglt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_siglt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_q2_siglt_fit.SetPointError(i, 0.0, siglt_X_fit_err)
        siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
        g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)

    r_siglt_fit = g_siglt_fit.Fit(f_sigLT, "SQ")
    f_sigLT.Draw("same")

    #f_sigLT_status = (r_siglt_fit.Status() == 0 and r_siglt_fit.IsValid())
    f_sigLT_status = f_sigLT.GetNDF() != 0
    f_sigLT_status_message = "Fit Successful" if f_sigLT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigLT_status_message))

    c1.cd(3)
    g_siglt_fit_tot.SetMarkerStyle(26)
    g_siglt_fit_tot.SetMarkerColor(2)
    g_siglt_fit_tot.SetLineColor(2)
    g_siglt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigLT_y = float('inf')
    max_sigLT_y = float('-inf')

    # Update min_sigLT_y and max_sigLT_y based on each graph's values
    for graph in [graph_sigLT_p9, graph_sigLT_p10]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigLT_y:
                min_sigLT_y = y
            if y > max_sigLT_y:
                max_sigLT_y = y

    # Scale the y-axis
    graph_sigLT_p9.SetMinimum(min_sigLT_y * 0.9)
    graph_sigLT_p9.SetMaximum(max_sigLT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(3).SetLeftMargin(0.12)
    graph_sigLT_p9.SetTitle("Sig LT Parameter Convergence;Optimization Run;Parameter")
    graph_sigLT_p9.SetLineColor(ROOT.kRed)
    graph_sigLT_p10.SetLineColor(ROOT.kBlue)
    graph_sigLT_p9.Draw("ALP")
    graph_sigLT_p10.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(3).SetLeftMargin(0.12)
    graph_sigLT_chi2.SetTitle("Sig LT Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigLT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigLT_chi2.Draw("ALP")

    # Plot temperature convergence
    c5.cd(3).SetLeftMargin(0.12)
    graph_sigLT_temp.SetTitle("Sig LT Temperature Convergence;Optimization Run;Temperature")
    graph_sigLT_temp.SetLineColor(ROOT.kBlack)
    graph_sigLT_temp.Draw("ALP")

    # Plot acceptance probability convergence
    c6.cd(3).SetLeftMargin(0.12)
    graph_sigLT_accept.SetTitle("Sig LT Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigLT_accept.SetLineColor(ROOT.kBlack)
    graph_sigLT_accept.Draw("ALP")

    print("\n")    

    '''    
    # 3 params
    #########
    # SigLT #
    #########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig LT")
    print("/*--------------------------------------------------*/")    

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0

    # Store the parameter values and chi-square values for each iteration
    params_sigLT_history = {'p9': [], 'p10': [], 'p11': []}

    # Create TGraphs for parameter convergence
    graph_sigLT_p9 = TGraph()
    graph_sigLT_p10 = TGraph()
    graph_sigLT_p11 = TGraph()
    graph_sigLT_chi2 = TGraph()
    graph_sigLT_temp = TGraph()
    graph_sigLT_accept = TGraph()

    c1.cd(3).SetLeftMargin(0.12)
    nsep.Draw("siglt:t:siglt_e", "", "goff")

    # Record the start time
    start_time = time.time()

    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0

        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_siglt_0 = lt0
        par_siglt_1 = lt1
        par_siglt_2 = lt2
        par_siglt_err_0 = 0.0
        par_siglt_err_1 = 0.0
        par_siglt_err_2 = 0.0

        # Track the best solution
        best_params = [par_siglt_0, par_siglt_1, par_siglt_2]
        best_cost = float('inf')
        previous_params = best_params[:]
        best_errors = [par_siglt_err_0, par_siglt_err_1, par_siglt_err_2]

        # Check for local minima
        local_minima = []
        local_iterations = 0
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_siglt_prv = TGraph()
            g_siglt_fit = TGraphErrors()
            g_siglt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_siglt_0, temperature),
                    simulated_annealing(par_siglt_1, temperature),
                    simulated_annealing(par_siglt_2, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [lt0, lt1, lt2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_siglt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_siglt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_siglt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    siglt_X_fit = g_siglt.GetY()[i]
                    siglt_X_fit_err = g_siglt.GetEY()[i]

                    g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
                    g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

                f_sigLT = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 3)
                f_sigLT.SetParNames("p9", "p10", "p11")
                f_sigLT.SetParameter(0, current_params[0])
                f_sigLT.SetParameter(1, current_params[1])
                f_sigLT.SetParameter(2, current_params[2])
                #f_sigLT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_siglt_0), current_params[0] + abs(current_params[0] * par_siglt_0))
                #f_sigLT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_siglt_1), current_params[1] + abs(current_params[1] * par_siglt_1))
                #f_sigLT.SetParLimits(2, current_params[2] - abs(current_params[2] * par_siglt_2), current_params[2] + abs(current_params[2] * par_siglt_2))
                #f_sigLT.SetParLimits(0, -5, 5)
                #f_sigLT.SetParLimits(1, -5, 5)
                #f_sigLT.SetParLimits(2, -5, 5)

                g_q2_siglt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
                    g_q2_siglt_fit.SetPointError(i, 0.0, siglt_X_fit_err)
                    siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
                    g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)

                r_siglt_fit = g_siglt_fit.Fit(f_sigLT, "SQ")

                #f_sigLT_status = (r_siglt_fit.Status() == 0 and r_siglt_fit.IsValid())
                f_sigLT_status = f_sigLT.GetNDF() != 0

                params_sigLT_history['p9'].append(current_params[0])
                params_sigLT_history['p10'].append(current_params[1])
                params_sigLT_history['p11'].append(current_params[2])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigLT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigLT.GetParameter(0),
                    f_sigLT.GetParameter(1),
                    f_sigLT.GetParameter(2)
                ]

                current_errors = [
                    f_sigLT.GetParError(0),
                    f_sigLT.GetParError(1),
                    f_sigLT.GetParError(2)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigLT_p9.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigLT_p10.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigLT_p11.SetPoint(total_iteration, total_iteration, current_params[2])
                graph_sigLT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigLT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigLT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

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
                    current_params = local_search(current_params, f_sigLT, 3)
                    par_siglt_0, par_siglt_1, par_siglt_2 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigLT_history['p9']) >= max_unchanged_iterations  and \
                   len(params_sigLT_history['p10']) >= max_unchanged_iterations  and \
                   len(params_sigLT_history['p11']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigLT_history['p9'][-2], 3), round(params_sigLT_history['p9'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigLT_history['p10'][-2], 3), round(params_sigLT_history['p10'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigLT_history['p11'][-2], 3), round(params_sigLT_history['p11'][-1], 3), atol=5.0):
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
                    current_params = [lt0, lt1, lt2]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]                

                # Update parameters with the best found so far
                par_siglt_0, par_siglt_1, par_siglt_2 = best_params
                par_siglt_err_0, par_siglt_err_1, par_siglt_err_2 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0

                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1], current_params[2]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p9={:.3e}, p10={:.3e}, p11={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1], current_params[2]))

                    current_params = adjust_params(best_params)
                    par_siglt_0, par_siglt_1, par_siglt_2 = current_params
                    par_siglt_err_0, par_siglt_err_1, par_siglt_err_2 = [0.0 for _ in range(3)]

            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                par_siglt_0, par_siglt_1, par_siglt_2 = [lt0, lt1, lt2]
                par_siglt_err_0, par_siglt_err_1, par_siglt_err_2 = [0.0 for _ in range(3)]

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0                

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_siglt_prv = TGraph()
    g_siglt_fit = TGraphErrors()
    g_siglt_fit_tot = TGraph()    

    f_sigLT_pre = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 3)
    f_sigLT_pre.SetParNames("p9", "p10", "p11")
    f_sigLT_pre.FixParameter(0, best_overall_params[0])
    f_sigLT_pre.FixParameter(1, best_overall_params[1])
    f_sigLT_pre.FixParameter(2, best_overall_params[2])                

    g_siglt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_siglt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_siglt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        siglt_X_pre = (f_sigLT_pre.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
        g_siglt_prv.SetPoint(i, g_siglt.GetX()[i], siglt_X_pre)

        siglt_X_fit = g_siglt.GetY()[i]
        siglt_X_fit_err = g_siglt.GetEY()[i]

        g_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_siglt_fit.SetPointError(i, 0, siglt_X_fit_err)

    g_siglt.SetTitle("Sig LT")
    g_siglt.SetMarkerStyle(5)
    g_siglt.Draw("AP")
    g_siglt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_siglt.GetXaxis().CenterTitle()
    g_siglt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
    g_siglt.GetYaxis().SetTitleOffset(1.5)
    g_siglt.GetYaxis().SetTitleSize(0.035)
    g_siglt.GetYaxis().CenterTitle()

    g_siglt_prv.SetMarkerColor(4)
    g_siglt_prv.SetMarkerStyle(25)
    g_siglt_prv.Draw("P")

    c2.cd(3).SetLeftMargin(0.12)
    g_siglt_fit.SetTitle("Sigma LT Model Fit")
    g_siglt_fit.Draw("A*")

    g_siglt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_siglt_fit.GetXaxis().CenterTitle()
    g_siglt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
    g_siglt_fit.GetYaxis().SetTitleOffset(1.5)
    g_siglt_fit.GetYaxis().SetTitleSize(0.035)
    g_siglt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_siglt_fit.GetX())
    x_max = max(g_siglt_fit.GetX())
    y_min = min(g_siglt_fit.GetY())
    y_max = max(g_siglt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_siglt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_siglt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigLT = TF1("sig_LT", fun_Sig_LT, tmin_range, tmax_range, 3)
    f_sigLT.SetParNames("p9", "p10", "p11")
    f_sigLT.FixParameter(0, best_overall_params[0])
    f_sigLT.FixParameter(1, best_overall_params[1])
    f_sigLT.FixParameter(2, best_overall_params[2])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigLT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_siglt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_siglt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_siglt_fit.SetPoint(i, g_siglt.GetX()[i], siglt_X_fit)
        g_q2_siglt_fit.SetPointError(i, 0.0, siglt_X_fit_err)
        siglt_X = (f_sigLT.Eval(g_siglt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)) * (g_vec[i])
        g_siglt_fit_tot.SetPoint(i, g_siglt.GetX()[i], siglt_X)

    r_siglt_fit = g_siglt_fit.Fit(f_sigLT, "SQ")
    f_sigLT.Draw("same")

    #f_sigLT_status = (r_siglt_fit.Status() == 0 and r_siglt_fit.IsValid())
    f_sigLT_status = f_sigLT.GetNDF() != 0
    f_sigLT_status_message = "Fit Successful" if f_sigLT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigLT_status_message))

    c1.cd(3)
    g_siglt_fit_tot.SetMarkerStyle(26)
    g_siglt_fit_tot.SetMarkerColor(2)
    g_siglt_fit_tot.SetLineColor(2)
    g_siglt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigLT_y = float('inf')
    max_sigLT_y = float('-inf')

    # Update min_sigLT_y and max_sigLT_y based on each graph's values
    for graph in [graph_sigLT_p9, graph_sigLT_p10, graph_sigLT_p11]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigLT_y:
                min_sigLT_y = y
            if y > max_sigLT_y:
                max_sigLT_y = y

    # Scale the y-axis
    graph_sigLT_p9.SetMinimum(min_sigLT_y * 0.9)
    graph_sigLT_p9.SetMaximum(max_sigLT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(3).SetLeftMargin(0.12)
    graph_sigLT_p9.SetTitle("Sig LT Parameter Convergence;Optimization Run;Parameter")
    graph_sigLT_p9.SetLineColor(ROOT.kRed)
    graph_sigLT_p10.SetLineColor(ROOT.kBlue)
    graph_sigLT_p11.SetLineColor(ROOT.kGreen)
    graph_sigLT_p9.Draw("ALP")
    graph_sigLT_p10.Draw("LP SAME")
    graph_sigLT_p11.Draw("LP SAME")

    # Plot chi-square convergence
    c4.cd(3).SetLeftMargin(0.12)
    graph_sigLT_chi2.SetTitle("Sig LT Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigLT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigLT_chi2.Draw("ALP")

    # Plot temperature convergence
    c5.cd(3).SetLeftMargin(0.12)
    graph_sigLT_temp.SetTitle("Sig LT Temperature Convergence;Optimization Run;Temperature")
    graph_sigLT_temp.SetLineColor(ROOT.kBlack)
    graph_sigLT_temp.Draw("ALP")

    # Plot acceptance probability convergence
    c6.cd(3).SetLeftMargin(0.12)
    graph_sigLT_accept.SetTitle("Sig LT Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigLT_accept.SetLineColor(ROOT.kBlack)
    graph_sigLT_accept.Draw("ALP")

    print("\n")    
    '''

    '''
    # 1 param
    #########
    # SigTT #
    #########
    
    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig TT")
    print("/*--------------------------------------------------*/")
    
    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0
    
    # Store the parameter values and chi-square values for each iteration
    params_sigTT_history = {'p13': []}

    # Create TGraphs for parameter convergence
    graph_sigTT_p13 = TGraph()
    graph_sigTT_chi2 = TGraph()
    graph_sigTT_temp = TGraph()
    graph_sigTT_accept = TGraph()
    
    c1.cd(4).SetLeftMargin(0.12)
    nsep.Draw("sigtt:t:sigtt_e", "", "goff")

    # Record the start time
    start_time = time.time()
    
    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0
        
        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigtt_0 = tt0
        par_sigtt_err_0 = 0.0

        # Track the best solution
        best_params = par_sigtt_0
        best_cost = float('inf')
        previous_params = best_params
        best_errors = par_sigtt_err_0
        
        # Check for local minima
        local_minima = []
        local_iterations = 0
        tabu_list = set()

        # Local search
        local_search_interval = 25

        while iteration <= max_iterations:

            g_sigtt_prv = TGraph()
            g_sigtt_fit = TGraphErrors()
            g_sigtt_fit_tot = TGraph()    
            
            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = simulated_annealing(par_sigtt_0, temperature)

                # Insert tabu list check here
                if current_params not in tabu_list:
                    tabu_list.add(current_params)
                else:
                    # Restart from initial parameters
                    current_params = tt0
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigtt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigtt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigtt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigtt_X_fit = g_sigtt.GetY()[i]
                    sigtt_X_fit_err = g_sigtt.GetEY()[i]

                    g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
                    g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)

                f_sigTT = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
                f_sigTT.SetParNames("p13")
                f_sigTT.SetParameter(0, current_params)
                #f_sigTT.SetParLimits(0, current_params - abs(current_params * par_sigtt_0), current_params + abs(current_params * par_sigtt_0))
                #f_sigTT.SetParLimits(0, -5, 5)

                g_q2_sigtt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
                    g_q2_sigtt_fit.SetPointError(i, 0.0, sigtt_X_fit_err)
                    sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                    g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)

                r_sigtt_fit = g_sigtt_fit.Fit(f_sigTT, "SQ")

                #f_sigTT_status = (r_sigtt_fit.Status() == 0 and r_sigtt_fit.IsValid())
                f_sigTT_status = f_sigTT.GetNDF() != 0

                params_sigTT_history['p13'].append(current_params)

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigTT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)
                
                current_params = f_sigTT.GetParameter(0)

                current_errors = f_sigTT.GetParError(0)

                # Update ROOT TGraphs for plotting
                graph_sigTT_p13.SetPoint(total_iteration, total_iteration, current_params)
                graph_sigTT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigTT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigTT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))
                
                # If the new cost is better or accepted by the acceptance probability, update the best parameters
                if accept_prob > random.random():
                    best_params = current_params
                    best_cost = current_cost
                    best_errors = current_errors

                if iteration % local_search_interval == 0:
                    current_params = local_search(current_params, f_sigTT, 1)
                    par_sigtt_0 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigTT_history['p13']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigTT_history['p13'][-2], 3), round(params_sigTT_history['p13'][-1], 3), atol=5.0):
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
                    current_params = tt0
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params                

                # Update parameters with the best found so far
                par_sigtt_0 = best_params
                par_sigtt_err_0 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0
                
                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p13={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params))

                    current_params = adjust_params(best_params)
                    par_sigtt_0 = current_params
                    par_sigtt_err_0 = 0.0
                    
            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))

                # Adjust parameter limits within a random number
                par_sigtt_0 = tt0
                par_sigtt_err_0 = 0.0

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0                

        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)

    g_sigtt_prv = TGraph()
    g_sigtt_fit = TGraphErrors()
    g_sigtt_fit_tot = TGraph()    

    f_sigTT_pre = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
    f_sigTT_pre.SetParNames("p13")
    f_sigTT_pre.FixParameter(0, best_overall_params[0])

    g_sigtt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigtt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigtt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigtt_X_pre = (f_sigTT_pre.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
        g_sigtt_prv.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_pre)

        sigtt_X_fit = g_sigtt.GetY()[i]
        sigtt_X_fit_err = g_sigtt.GetEY()[i]

        g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)

    g_sigtt.SetTitle("Sig TT")
    g_sigtt.SetMarkerStyle(5)
    g_sigtt.Draw("AP")
    g_sigtt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigtt.GetXaxis().CenterTitle()
    g_sigtt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
    g_sigtt.GetYaxis().SetTitleOffset(1.5)
    g_sigtt.GetYaxis().SetTitleSize(0.035)
    g_sigtt.GetYaxis().CenterTitle()

    g_sigtt_prv.SetMarkerColor(4)
    g_sigtt_prv.SetMarkerStyle(25)
    g_sigtt_prv.Draw("P")

    c2.cd(4).SetLeftMargin(0.12)
    g_sigtt_fit.SetTitle("Sigma TT Model Fit")
    g_sigtt_fit.Draw("A*")

    g_sigtt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigtt_fit.GetXaxis().CenterTitle()
    g_sigtt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
    g_sigtt_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigtt_fit.GetYaxis().SetTitleSize(0.035)
    g_sigtt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigtt_fit.GetX())
    x_max = max(g_sigtt_fit.GetX())
    y_min = min(g_sigtt_fit.GetY())
    y_max = max(g_sigtt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigtt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigtt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigTT = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
    f_sigTT.SetParNames("p13")
    f_sigTT.FixParameter(0, best_overall_params[0])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigTT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigtt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigtt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_q2_sigtt_fit.SetPointError(i, 0.0, sigtt_X_fit_err)
        sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
        g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)

    r_sigtt_fit = g_sigtt_fit.Fit(f_sigTT, "SQ")
    f_sigTT.Draw("same")

    #f_sigTT_status = (r_sigtt_fit.Status() == 0 and r_sigtt_fit.IsValid())
    f_sigTT_status = f_sigTT.GetNDF() != 0
    f_sigTT_status_message = "Fit Successful" if f_sigTT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigTT_status_message))

    c1.cd(4)
    g_sigtt_fit_tot.SetMarkerStyle(26)
    g_sigtt_fit_tot.SetMarkerColor(2)
    g_sigtt_fit_tot.SetLineColor(2)
    g_sigtt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigTT_y = float('inf')
    max_sigTT_y = float('-inf')

    # Update min_sigTT_y and max_sigTT_y based on each graph's values
    for graph in [graph_sigTT_p13]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigTT_y:
                min_sigTT_y = y
            if y > max_sigTT_y:
                max_sigTT_y = y

    # Scale the y-axis
    graph_sigTT_p13.SetMinimum(min_sigTT_y * 0.9)
    graph_sigTT_p13.SetMaximum(max_sigTT_y * 1.1)    

    # Plot parameter convergence
    c3.cd(4).SetLeftMargin(0.12)
    graph_sigTT_p13.SetTitle("Sig TT Parameter Convergence;Optimization Run;Parameter")
    graph_sigTT_p13.SetLineColor(ROOT.kRed)
    graph_sigTT_p13.Draw("ALP")
    
    # Plot chi-square convergence
    c4.cd(4).SetLeftMargin(0.12)
    graph_sigTT_chi2.SetTitle("Sig TT Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigTT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigTT_chi2.Draw("ALP")
    
    # Plot temperature convergence
    c5.cd(4).SetLeftMargin(0.12)
    graph_sigTT_temp.SetTitle("Sig TT Temperature Convergence;Optimization Run;Temperature")
    graph_sigTT_temp.SetLineColor(ROOT.kBlack)
    graph_sigTT_temp.Draw("ALP")
    
    # Plot acceptance probability convergence
    c6.cd(4).SetLeftMargin(0.12)
    graph_sigTT_accept.SetTitle("Sig TT Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigTT_accept.SetLineColor(ROOT.kBlack)
    graph_sigTT_accept.Draw("ALP")
    
    print("\n")    
    '''    

    # 2 params
    #########
    # SigTT #
    #########

    print("\n/*--------------------------------------------------*/")
    print("Fit for Sig TT")
    print("/*--------------------------------------------------*/")    

    num_starts = 10  # Number of times to restart the algorithm
    best_overall_params = None
    best_overall_cost = float('inf')
    total_iteration = 0
    
    # Store the parameter values and chi-square values for each iteration
    params_sigTT_history = {'p13': [], 'p14': []}

    # Create TGraphs for parameter convergence
    graph_sigTT_p13 = TGraph()
    graph_sigTT_p14 = TGraph()
    graph_sigTT_chi2 = TGraph()
    graph_sigTT_temp = TGraph()
    graph_sigTT_accept = TGraph()
    
    c1.cd(4).SetLeftMargin(0.12)
    nsep.Draw("sigtt:t:sigtt_e", "", "goff")
    
    # Record the start time
    start_time = time.time()
        
    for start in range(num_starts):
        print("\nStarting optimization run {0}/{1}".format(start + 1, num_starts))    

        iteration = 0
    
        initial_temperature = 1.0
        temperature = initial_temperature
        unchanged_iterations = 0
        max_unchanged_iterations = 5

        # Initialize adaptive parameter limits
        par_sigtt_0 = tt0
        par_sigtt_1 = tt1
        par_sigtt_err_0 = 0.0
        par_sigtt_err_1 = 0.0
        par_sigtt_err_2 = 0.0

        # Track the best solution
        best_params = [par_sigtt_0, par_sigtt_1]
        best_cost = float('inf')
        previous_params = best_params[:]
        best_errors = [par_sigtt_err_0, par_sigtt_err_1]
        
        # Check for local minima
        local_minima = []
        local_iterations = 0
        tabu_list = set()

        # Local search
        local_search_interval = 25
        
        while iteration <= max_iterations:
            
            g_sigtt_prv = TGraph()
            g_sigtt_fit = TGraphErrors()
            g_sigtt_fit_tot = TGraph()    

            sys.stdout.write(" \rSearching for best parameters...({0}/{1})\r{2}".format(iteration, max_iterations, ''))
            sys.stdout.flush()

            try:
                # Perturb parameters
                current_params = [
                    simulated_annealing(par_sigtt_0, temperature),
                    simulated_annealing(par_sigtt_1, temperature)
                ]

                # Insert tabu list check here
                if tuple(current_params) not in tabu_list:
                    tabu_list.add(tuple(current_params))
                else:
                    # Restart from initial parameters
                    current_params = [tt0, tt1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                g_sigtt = TGraphErrors()
                for i in range(nsep.GetSelectedRows()):
                    g_sigtt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
                    g_sigtt.SetPointError(i, 0, nsep.GetV3()[i])

                for i in range(len(w_vec)):
                    sigtt_X_fit = g_sigtt.GetY()[i]
                    sigtt_X_fit_err = g_sigtt.GetEY()[i]

                    g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
                    g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)

                f_sigTT = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
                f_sigTT.SetParNames("p13", "p14")
                f_sigTT.SetParameter(0, current_params[0])
                f_sigTT.SetParameter(1, current_params[1])
                #f_sigTT.SetParLimits(0, current_params[0] - abs(current_params[0] * par_sigtt_0), current_params[0] + abs(current_params[0] * par_sigtt_0))
                #f_sigTT.SetParLimits(1, current_params[1] - abs(current_params[1] * par_sigtt_1), current_params[1] + abs(current_params[1] * par_sigtt_1))
                #f_sigTT.SetParLimits(0, -5, 5)
                #f_sigTT.SetParLimits(1, -5, 5)

                g_q2_sigtt_fit = TGraphErrors()
                for i in range(len(w_vec)):
                    g_q2_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
                    g_q2_sigtt_fit.SetPointError(i, 0.0, sigtt_X_fit_err)
                    sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
                    g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)

                r_sigtt_fit = g_sigtt_fit.Fit(f_sigTT, "SQ")

                #f_sigTT_status = (r_sigtt_fit.Status() == 0 and r_sigtt_fit.IsValid())
                f_sigTT_status = f_sigTT.GetNDF() != 0

                params_sigTT_history['p13'].append(current_params[0])
                params_sigTT_history['p14'].append(current_params[1])

                # Calculate the cost (chi-square value) for the current parameters
                current_cost = f_sigTT.GetChisquare()

                # Acceptance probability
                accept_prob = acceptance_probability(best_cost, current_cost, temperature)

                current_params = [
                    f_sigTT.GetParameter(0),
                    f_sigTT.GetParameter(1)
                ]

                current_errors = [
                    f_sigTT.GetParError(0),
                    f_sigTT.GetParError(1)
                ]

                # Update ROOT TGraphs for plotting
                graph_sigTT_p13.SetPoint(total_iteration, total_iteration, current_params[0])
                graph_sigTT_p14.SetPoint(total_iteration, total_iteration, current_params[1])
                graph_sigTT_chi2.SetPoint(total_iteration, total_iteration, round(current_cost, 4))
                graph_sigTT_temp.SetPoint(total_iteration, total_iteration, temperature)
                graph_sigTT_accept.SetPoint(total_iteration, total_iteration, round(accept_prob, 4))

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
                    current_params = local_search(current_params, f_sigTT, 3)
                    par_sigtt_0, par_sigtt_1 = current_params

                # Check if current parameters haven't changed for the past N iterations
                if len(params_sigTT_history['p13']) >= max_unchanged_iterations  and \
                   len(params_sigTT_history['p14']) >= max_unchanged_iterations:
                    if np.allclose(round(params_sigTT_history['p13'][-2], 3), round(params_sigTT_history['p13'][-1], 3), atol=5.0) and \
                       np.allclose(round(params_sigTT_history['p14'][-2], 3), round(params_sigTT_history['p14'][-1], 3), atol=5.0):
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
                    current_params = [tt0, tt1]
                    temperature = initial_temperature
                    unchanged_iterations = 0

                previous_params = current_params[:]                

                # Update parameters with the best found so far
                par_sigtt_0, par_sigtt_1 = best_params
                par_sigtt_err_0, par_sigtt_err_1 = best_errors

                # Update the temperature
                temperature = adaptive_cooling(initial_temperature, iteration, max_iterations)

                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0
                
                # Check if current_params are close to any local minimum
                if any(np.allclose([current_params[0], current_params[1]], minima, atol=5.0) for minima in local_minima):
                    #print("WARNING: Parameters p13={:.3e}, p14={:.3e} are a local minima. Adjusting parameter limits and retrying...".format(current_params[0], current_params[1]))

                    current_params = adjust_params(best_params)
                    par_sigtt_0, par_sigtt_1 = current_params
                    par_sigtt_err_0, par_sigtt_err_1 = [0.0 for _ in range(2)]
                    
            except (TypeError or ZeroDivisionError) as e:
                #print("WARNING: {}, Adjusting parameter limits and retrying...".format(e))
                
                par_sigtt_0, par_sigtt_1 = [tt0, tt1]
                par_sigtt_err_0, par_sigtt_err_1 = [0.0 for _ in range(2)]
                
                iteration += 1
                total_iteration += 1 if iteration % max_iterations == 0 else 0                
                
        # After the while loop, check if this run found a better solution
        if best_cost < best_overall_cost:
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

    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
    par_chi2_vec.append(best_cost)
            
    g_sigtt_prv = TGraph()
    g_sigtt_fit = TGraphErrors()
    g_sigtt_fit_tot = TGraph()    

    f_sigTT_pre = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
    f_sigTT_pre.SetParNames("p13", "p14")
    f_sigTT_pre.FixParameter(0, best_overall_params[0])
    f_sigTT_pre.FixParameter(1, best_overall_params[1])

    g_sigtt = TGraphErrors()
    for i in range(nsep.GetSelectedRows()):
        g_sigtt.SetPoint(i, nsep.GetV2()[i], nsep.GetV1()[i])
        g_sigtt.SetPointError(i, 0, nsep.GetV3()[i])

    for i in range(len(w_vec)):
        sigtt_X_pre = (f_sigTT_pre.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
        g_sigtt_prv.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_pre)

        sigtt_X_fit = g_sigtt.GetY()[i]
        sigtt_X_fit_err = g_sigtt.GetEY()[i]

        g_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_sigtt_fit.SetPointError(i, 0, sigtt_X_fit_err)

    g_sigtt.SetTitle("Sig TT")
    g_sigtt.SetMarkerStyle(5)
    g_sigtt.Draw("AP")
    g_sigtt.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigtt.GetXaxis().CenterTitle()
    g_sigtt.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
    g_sigtt.GetYaxis().SetTitleOffset(1.5)
    g_sigtt.GetYaxis().SetTitleSize(0.035)
    g_sigtt.GetYaxis().CenterTitle()

    g_sigtt_prv.SetMarkerColor(4)
    g_sigtt_prv.SetMarkerStyle(25)
    g_sigtt_prv.Draw("P")

    c2.cd(4).SetLeftMargin(0.12)
    g_sigtt_fit.SetTitle("Sigma TT Model Fit")
    g_sigtt_fit.Draw("A*")

    g_sigtt_fit.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
    g_sigtt_fit.GetXaxis().CenterTitle()
    g_sigtt_fit.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
    g_sigtt_fit.GetYaxis().SetTitleOffset(1.5)
    g_sigtt_fit.GetYaxis().SetTitleSize(0.035)
    g_sigtt_fit.GetYaxis().CenterTitle()

    # Set axis limits to ensure everything is shown
    x_min = min(g_sigtt_fit.GetX())
    x_max = max(g_sigtt_fit.GetX())
    y_min = min(g_sigtt_fit.GetY())
    y_max = max(g_sigtt_fit.GetY())

    # You can also set a margin to ensure all points are visible
    margin = 0.1
    g_sigtt_fit.GetXaxis().SetRangeUser(x_min - margin, x_max + margin)
    g_sigtt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)            

    f_sigTT = TF1("sig_TT", fun_Sig_TT, tmin_range, tmax_range, 2)
    f_sigTT.SetParNames("p13", "p14")
    f_sigTT.FixParameter(0, best_overall_params[0])
    f_sigTT.FixParameter(1, best_overall_params[1])

    # Evaluate the fit function at several points to determine its range
    n_points = 100  # Number of points to evaluate the fit function
    fit_y_values = [f_sigTT.Eval(x) for x in np.linspace(tmin_range, tmax_range, n_points)]
    fit_y_min = min(fit_y_values)
    fit_y_max = max(fit_y_values)

    # Extend the y-axis range to include the fit function range
    y_min = min(y_min, fit_y_min)
    y_max = max(y_max, fit_y_max)

    # Set a margin to ensure all points are visible
    margin = 0.1 * (y_max - y_min)
    g_sigtt_fit.GetYaxis().SetRangeUser(y_min - margin, y_max + margin)

    g_q2_sigtt_fit = TGraphErrors()
    for i in range(len(w_vec)):
        g_q2_sigtt_fit.SetPoint(i, g_sigtt.GetX()[i], sigtt_X_fit)
        g_q2_sigtt_fit.SetPointError(i, 0.0, sigtt_X_fit_err)
        sigtt_X = (f_sigTT.Eval(g_sigtt.GetX()[i]) * math.sin(th_vec[i] * PI / 180)**2) * (g_vec[i])
        g_sigtt_fit_tot.SetPoint(i, g_sigtt.GetX()[i], sigtt_X)

    r_sigtt_fit = g_sigtt_fit.Fit(f_sigTT, "SQ")
    f_sigTT.Draw("same")

    #f_sigTT_status = (r_sigtt_fit.Status() == 0 and r_sigtt_fit.IsValid())
    f_sigTT_status = f_sigTT.GetNDF() != 0
    f_sigTT_status_message = "Fit Successful" if f_sigTT_status else "Fit Failed"

    fit_status = TText()
    fit_status.SetTextSize(0.04)
    fit_status.DrawTextNDC(0.35, 0.85, " Fit Status: {}".format(f_sigTT_status_message))

    c1.cd(4)
    g_sigtt_fit_tot.SetMarkerStyle(26)
    g_sigtt_fit_tot.SetMarkerColor(2)
    g_sigtt_fit_tot.SetLineColor(2)
    g_sigtt_fit_tot.Draw("LP")

    # Calculate the minimum and maximum values from the graphs
    min_sigTT_y = float('inf')
    max_sigTT_y = float('-inf')

    # Update min_sigTT_y and max_sigTT_y based on each graph's values
    for graph in [graph_sigTT_p13, graph_sigTT_p14]:
        n_points = graph.GetN()
        for i in range(n_points):
            y = graph.GetY()[i]
            if y < min_sigTT_y:
                min_sigTT_y = y
            if y > max_sigTT_y:
                max_sigTT_y = y

    # Scale the y-axis
    graph_sigTT_p13.SetMinimum(min_sigTT_y * 0.9)
    graph_sigTT_p13.SetMaximum(max_sigTT_y * 1.1)    
    
    # Plot parameter convergence
    c3.cd(4).SetLeftMargin(0.12)
    graph_sigTT_p13.SetTitle("Sig TT Parameter Convergence;Optimization Run;Parameter")
    graph_sigTT_p13.SetLineColor(ROOT.kRed)
    graph_sigTT_p14.SetLineColor(ROOT.kBlue)
    graph_sigTT_p13.Draw("ALP")
    graph_sigTT_p14.Draw("LP SAME")
    
    # Plot chi-square convergence
    c4.cd(4).SetLeftMargin(0.12)
    graph_sigTT_chi2.SetTitle("Sig TT Chi-Square Convergence;Optimization Run;Chi-Square")
    graph_sigTT_chi2.SetLineColor(ROOT.kBlack)
    graph_sigTT_chi2.Draw("ALP")
    
    # Plot temperature convergence
    c5.cd(4).SetLeftMargin(0.12)
    graph_sigTT_temp.SetTitle("Sig TT Temperature Convergence;Optimization Run;Temperature")
    graph_sigTT_temp.SetLineColor(ROOT.kBlack)
    graph_sigTT_temp.Draw("ALP")
    
    # Plot acceptance probability convergence
    c6.cd(4).SetLeftMargin(0.12)
    graph_sigTT_accept.SetTitle("Sig TT Acceptance Probability Convergence;Optimization Run;Acceptance Probability")
    graph_sigTT_accept.SetLineColor(ROOT.kBlack)
    graph_sigTT_accept.Draw("ALP")

    print("\n")    
    
    c1.Print(outputpdf+'(')
    c2.Print(outputpdf)
    c3.Print(outputpdf)
    c4.Print(outputpdf)
    c5.Print(outputpdf)
    c6.Print(outputpdf+')')
    
    for i, (old, new) in enumerate(zip(prv_par_vec, par_vec)):
        if old != new:
            print("par{} changed from {:.3f} to {:.3f}".format(i+1, old, new))
            
    # Don't write to new parameter file if debugging
    if not DEBUG:
        para_file_out = "{}/src/{}/parameters/par.{}_Q{}W{}.dat".format(LTANAPATH, ParticleType, pol_str, q2_set.replace("p",""), w_set.replace("p",""))
        print("\nWriting {}...".format(para_file_out))
        with open(para_file_out, 'w') as f:
            for i in range(len(par_vec)):
                f.write("{:13.5e} {:13.5e} {:3d} {:12.1f}\n".format(par_vec[i], par_err_vec[i], i+1, par_chi2_vec[i]))
                print("  {:.3f} {:.3f} {:.1f} {:.1f}".format(par_vec[i], par_err_vec[i], i+1, par_chi2_vec[i]))
    else:
        print("\n\nDEBUG ENABLED: No changes to previous iteration...")
        for i,par in enumerate(prv_par_vec):
            print("par{} = {:.3f}".format(i+1, par))

    # Testing
    #sys.exit(2)            
