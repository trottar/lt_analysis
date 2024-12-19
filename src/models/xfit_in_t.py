#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-19 08:51:24 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TNtuple
from ROOT import TCanvas
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
from utility import load_equations, prepare_equations, find_params_wrapper, check_chi_squared_values

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
def x_fit_in_t(ParticleType, pol_str, dir_iter, q2_set, w_set, inpDict):

    tmin_range = inpDict["tmin"]
    tmax_range = inpDict["tmax"]
    Q2min_range = inpDict["Q2min"]
    Q2max_range = inpDict["Q2max"]

    iter_num = inpDict["iter_num"]

    if iter_num > 1:
        DEBUG=False
    else:
        DEBUG=True
    
    #DEBUG=False
    
    ##############
    # HARD CODED #
    ##############
    # Maximum iterations before ending loop
    #max_iterations = 100
    max_iterations = 1000
    #max_iterations = 2000
    #max_iterations = 5000
    #max_iterations = 10000

    # Number of times to run the algorithm
    #num_optimizations = 5
    num_optimizations = 10
    #num_optimizations = 50
    #num_optimizations = 1000

    # Initial max/min bounds of finding parameter values
    initial_param_bounds = 1e4

    # Threshold on how bad red. chi2 can be
    #chi2_threshold = 1.0
    #chi2_threshold = 6.0
    #chi2_threshold = 10.0
    chi2_threshold = 600.0
    ##############
    ##############
    ##############

    # Set pol_str, q2_set, w_set for xfit_active script
    set_val(pol_str, q2_set, w_set)
    
    outputpdf  = "{}/{}_xfit_in_t_Q{}W{}.pdf".format(OUTPATH, ParticleType, q2_set, w_set)
    
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
    q2_set = q2_set.replace("p",".")
    w_set = w_set.replace("p",".")
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
    q2_set = q2_set.replace(".","p")
    w_set = w_set.replace(".","p")

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
        "outputpdf" : outputpdf,
        "fit_params" : fit_params,
    }

    par_vec = prv_par_vec
    par_err_vec = prv_err_vec
    par_chi2_vec = prv_chi2_vec

    fixed_params = ["L", "T", "LT", "TT"]
    #fixed_params = ["L", "T", "LT"]
    
    parameterize(inpDict, par_vec, par_err_vec, par_chi2_vec, fixed_params)
    
    if check_chi_squared_values(par_chi2_vec, chi2_threshold, fit_params, equations):
        sys.exit(2)

    # Check if parameter values changed and print changes to terminal
    for i, (old, new) in enumerate(zip(prv_par_vec, par_vec)):
        if old != new:
            print("par{} changed from {:.3e} to {:.3e}".format(i+1, old, new))
            
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
