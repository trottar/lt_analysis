#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-04-17 16:18:27 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TNtuple
from ROOT import TCanvas
import numpy as np
import math
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

TEMP_CACHEPATH=f"{OUTPATH}/cache_transfer"

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations, prepare_equations, find_params_wrapper, check_chi_squared_values, request_yn_response

##################################################################################################################################################
# Import fit finder function

from xfit_fit_finder import parameterize

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen

###############################################################################################################################################
# Import separated xsects models

from xfit_active import set_val

###############################################################################################################################################

# dir_iter = closest_date
def x_fit_in_t(ParticleType, pol_str, dir_iter, q2_set, w_set, inpDict, output_file_lst, skip_optimization=False):

    tmin_range = inpDict["tmin"]
    tmax_range = inpDict["tmax"]
    Q2min_range = inpDict["Q2min"]
    Q2max_range = inpDict["Q2max"]

    iter_num = inpDict["iter_num"]
        
    ##############
    # HARD CODED #
    ##############

    # True - Set range of parameter search from +/- max_param_bounds
    # False - Set range of parameter search from current_params +/- off*abs(current_params) (off = 0.1, 10% param value)
    full_optimization = True

    #skip_optimization = True # Set to True to skip optimization and use fixed parameters

    # Fixed separated xsect parameterization
    if skip_optimization:
        fixed_params = ["L", "T", "LT", "TT"] # Skip optimization
    else:
        #fixed_params = ["L", "T", "LT", "TT"] # Skip optimization
        #fixed_params = ["L", "LT"]
        #fixed_params = ["L", "LT", "TT"]
        #fixed_params = ["TT"]
        fixed_params = [] # Update all
    
    # Maximum iterations before ending loop (should always aim for >10000)
    #max_iterations = 1000
    #max_iterations = 10000
    #max_iterations = 15000
    #max_iterations = 25000
    max_iterations = 50000
    #max_iterations = 100000

    # Number of times to run the algorithm
    #num_optimizations = 1
    num_optimizations = 3
    #num_optimizations = 5
    #num_optimizations = 10

    # Initial max/min bounds of finding parameter values (only used for iter=1)
    #initial_param_bounds = 1e2
    #initial_param_bounds = 1e4
    initial_param_bounds = 1e6

    # Threshold value of red. chi2
    #chi2_threshold = 1.0
    chi2_threshold = 3.0
    #chi2_threshold = 5.0
    #chi2_threshold = 10.0
    #chi2_threshold = 30.0
    #chi2_threshold = 600.0

    # Number of rechecks to assure chi2_threshold not reached (set to zero to ignore chi2_threshold)
    max_checks = 0
    #max_checks = 1
    
    ##############
    ##############
    ##############
    
    # Set pol_str, q2_set, w_set for xfit_active script
    set_val(pol_str, q2_set, w_set)        
    
    prv_par_vec = []
    prv_err_vec = []
    prv_chi2_vec = []
    t_vec = []
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
    para_file_in =  "{}/{}/Q{}W{}/{}/parameters/par.{}_Q{}W{}.dat".format(TEMP_CACHEPATH, ParticleType, q2_set, w_set, dir_iter, \
                                                                            pol_str, q2_set.replace("p",""), w_set.replace("p",""))
    
    print("Reading {}...".format(para_file_in))
    try:
        with open(para_file_in, 'r') as f:
            for line in f:
                data = line.split()
                par, par_err, indx, chi2 = map(float, data)
                print("  {} {} {} {}".format(par, par_err, indx, chi2))
                prv_par_vec.append(par)
                prv_err_vec.append(par_err)
                prv_chi2_vec.append(chi2)
    except FileNotFoundError:
        print("File {} not found.".format(para_file_in))

    # Read in parameters, errors, chi2 of previous iteration
    l0, l1, l2, l3, t0, t1, t2, t3, lt0, lt1, lt2, lt3, tt0, tt1, tt2, tt3 = prv_par_vec[:16]
    errl0, errl1, errl2, errl3, errt0, errt1, errt2, errt3, errlt0, errlt1, errlt2, errlt3, errtt0, errtt1, errtt2, errtt3 = prv_err_vec[:16]
    chi2l0, chi2l1, chi2l2, chi2l3, chi2t0, chi2t1, chi2t2, chi2t3, chi2lt0, chi2lt1, chi2lt2, chi2lt3, chi2tt0, chi2tt1, chi2tt2, chi2tt3 = prv_chi2_vec[:16]

    # Load equations from model of given setting
    equations = load_equations(f"Q{q2_set}W{w_set}.model")
    
    ave_file_in = "{}/src/{}/averages/avek.Q{}W{}.dat".format(LTANAPATH, ParticleType, q2_set.replace("p",""), w_set.replace("p",""))
    # Redefine strings for retrieving equation defintions
    q2_set = float(q2_set.replace("p","."))
    w_set = float(w_set.replace("p","."))
    with open(ave_file_in, 'r') as f:
        for line in f:
            ww, ww_e, qq, qq_e, tt, tt_e, theta_cm, it = map(float, line.strip().split())

            # Grab functional form from model input file
            fun_wfactor_optimized = prepare_equations(equations, 'wfactor')

            # Calculate wfactor
            g = fun_wfactor_optimized(q2_set, w_set, qq, ww, tt)

            t_vec.append(tt)
            g_vec.append(g)
            w_vec.append(ww)
            q2_vec.append(qq)
            th_vec.append(theta_cm)
            
    # Revert changes for rest of script
    q2_set = str(q2_set).replace(".","p")
    w_set = str(w_set).replace(".","p")

    # Find fits for L, T, LT, TT
    fit_params = {
        "L": [l0, l1, l2, l3],
        "T": [t0, t1, t2, t3],
        "LT": [lt0, lt1, lt2, lt3],
        "TT": [tt0, tt1, tt2, tt3],
    }
    
    inp_dict = {
        "q2_set" : q2_set.replace("p","."),
        "w_set" : w_set.replace("p","."),
        "objects" : [nsep, t_vec, g_vec, w_vec, q2_vec, th_vec],
        "max_iterations" : max_iterations,
        "num_optimizations" : num_optimizations,
        "initial_param_bounds" : initial_param_bounds,
        "initial_params" : find_params_wrapper(equations),
        "tmin_range" : tmin_range,
        "tmax_range" : tmax_range,
        "Q2min_range" : Q2min_range,
        "Q2max_range" : Q2max_range,
        "iter_num" : iter_num,
        "fit_params" : fit_params,
        "chi2_threshold" : chi2_threshold,
        "xfit_log" : "{}/{}_xfit_in_t_Q{}W{}.log".format(OUTPATH, ParticleType, q2_set, w_set)
    }

    par_vec = prv_par_vec
    par_err_vec = prv_err_vec
    par_chi2_vec = prv_chi2_vec

    if skip_optimization:

        # Define output file name
        outputpdf  = "{}/{}_lt_fit_in_t_Q{}W{}.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
        output_file_lst.append(outputpdf)

        parameterize(inp_dict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization) #, True)
        
    else:
            
        # Define output file name
        outputpdf  = "{}/{}_xfit_in_t_Q{}W{}_Start.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
        output_file_lst.append(outputpdf)

        parameterize(inp_dict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization) #, True)
        bad_chi2_bool, bad_chi2_indices = check_chi_squared_values(par_chi2_vec, chi2_threshold, fit_params, equations)

        # Store initial values
        best_par_vec = par_vec.copy()
        best_err_vec = par_err_vec.copy()
        best_chi2_vec = par_chi2_vec.copy()

        i = 0    
        while bad_chi2_bool and i < max_checks:
            fixed_params = ["L", "T", "LT", "TT"] # Check and rerun any with bad chi2
            fixed_params = [x for i, x in enumerate(fixed_params) if i not in bad_chi2_indices] # Rerun any settings with bad chi2
            print(f"\n\nChi2 above threshold of {chi2_threshold}! Check ({i} / {max_checks})...")

            # Define output file name
            outputpdf  = "{}/{}_xfit_in_t_Q{}W{}_{}.pdf".format(OUTPATH, ParticleType, q2_set, w_set, i)
            output_file_lst.append(outputpdf)

            parameterize(inp_dict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization)
            bad_chi2_bool, bad_chi2_indices = check_chi_squared_values(par_chi2_vec, chi2_threshold, fit_params, equations)

            # Update best values for each group of 4 elements
            for j in range(0, len(par_chi2_vec), 4):
                if np.abs(np.mean(par_chi2_vec[j:j+4]) - 1) < np.abs(np.mean(best_chi2_vec[j:j+4]) - 1):
                    print(f"\n\nNew set of best chi2 values found for sig_{list(fit_params.keys())[j // 4]} of {par_chi2_vec[j]:.1f}...")
                    best_par_vec[j:j+4] = par_vec[j:j+4].copy()
                    best_err_vec[j:j+4] = par_err_vec[j:j+4].copy()
                    best_chi2_vec[j:j+4] = par_chi2_vec[j:j+4].copy()
            i += 1

        par_vec = best_par_vec
        par_err_vec = best_err_vec
        par_chi2_vec = best_chi2_vec

        prv_par_vec = par_vec
        prv_err_vec = par_err_vec
        prv_chi2_vec = par_chi2_vec

        # Define output file name
        outputpdf  = "{}/{}_xfit_in_t_Q{}W{}_Final.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
        output_file_lst.append(outputpdf)

        # Update plots with best chi2
        fixed_params = ["L", "T", "LT", "TT"] # Using best found chi2 from above for all
        parameterize(inp_dict, par_vec, par_err_vec, par_chi2_vec, prv_par_vec, prv_err_vec, prv_chi2_vec, fixed_params, outputpdf, full_optimization)
    
        # Check if parameter values changed and print changes to terminal
        for i, (old, new) in enumerate(zip(prv_par_vec, par_vec)):
            if old != new:
                print("par{} changed from {:.3e} to {:.3e}".format(i+1, old, new))

        para_file_out = "{}/src/{}/parameters/par.{}_Q{}W{}.dat".format(LTANAPATH, ParticleType, pol_str, q2_set.replace("p",""), w_set.replace("p",""))
        print("\nWriting {}...".format(para_file_out))
        with open(para_file_out, 'w') as f:
            for i in range(len(par_vec)):
                f.write("{:13.5e} {:13.5e} {:3d} {:12.1f}\n".format(par_vec[i], par_err_vec[i], i+1, par_chi2_vec[i]))
                print("  {:.3e} {:.3e} {:.1e} {:.1e}".format(par_vec[i], par_err_vec[i], i+1, par_chi2_vec[i]))

    '''
    print("\n\nWould you like to continue with the analysis?\n")
    if not request_yn_response():
        print("-"*25)
        print("Exiting script...")
        print("-"*25)
        sys.exit(2)
    '''
    print("\n\n")
    
