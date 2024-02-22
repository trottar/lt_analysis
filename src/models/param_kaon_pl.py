#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-21 19:16:20 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

def iterWeight(arg_str):

    # Define constants
    pi = 3.14159
    mtar_gev = 0.93827231
    mpipl=0.139570
    mkpl=0.493677

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, q2_sim, w_sim, t_sim, eps_sim, thetacm_sim, phicm_sim, sigcm_sim, wt_sim, *params = args
    
    q2_gev = q2_sim # Already GeV
    t_gev = t_sim  # Already GeV, issue here!!! t_sim makes no sense
    s = w_sim**2
    s_gev = s # Already GeV

    ##############
    # HARD CODED #
    ##############
    if q2_set == 2.1:
        q2_set = 2.115
    ##############
    ##############
    ##############

    # Calculate tav, ftav, ft
    #tav = (0.0735 + 0.028 * math.log(q2_set)) * q2_set
    # RLT (10/8/2023): Testing new tav parameterization
    tav=(0.1112 + 0.0066*math.log(q2_set))*q2_set
    ftav = (abs(t_gev) - tav) / tav
    ft = abs(t_gev) / (abs(t_gev) + mkpl**2)**2
    
    # Calculate sigl, sigt, siglt, sigtt, sig219, sig
    # RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for the
    #                 xfit_in_t.py script to work. LT/TT are zeros
    #p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = params
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params
    
    #sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev) - 0.2))
    # RLT (10/12/2023): Removed 0.2 to keep things as simple as possible for initial start parameterization
    # RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the extreme slope at high t
    #sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev)))
    sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev)+0.2))
    # RLT (2/15/2024): Removing t dependence from sigT because it seems
    #                  to be driving poor sep xsects results
    # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
    # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
    #sigt = p5 + p6 * math.log(q2_gev) + (p7 + p8 * math.log(q2_gev)) * ftav
    #sigt = p5 + p6 * math.log(q2_gev)
    #sigt = p5 * math.log(q2_gev) + p6 / (q2_gev**2)
    sigt = p5 / (1 + p6*q2_gev)
    siglt = (p9 * math.exp(p10 * abs(t_gev)) + p11 / abs(t_gev)) * math.sin(thetacm_sim)
    # RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for the
    #                 xfit_in_t.py script to work. LT/TT are zeros
    #                 Therefore param 12 was also changed to 13    
    sigtt = (p13 * q2_gev * math.exp(-q2_gev)) * ft * math.sin(thetacm_sim)**2

    # RLT (9/25/2023): There are two tav parameterizations in here.
    #                  I am only using the one above, for now.    
    # tav = (-0.178 + 0.315 * math.log(q2_gev)) * q2_gev

    sig219 = (sigt + eps_sim * sigl + eps_sim * math.cos(2. * phicm_sim) * sigtt +
             math.sqrt(2.0 * eps_sim * (1. + eps_sim)) * math.cos(phicm_sim) * siglt)

    wfactor = 1.0 / (s_gev - mtar_gev**2)**2
    sig = sig219*wfactor
    sigl = sigl*wfactor
    sigt = sigt*wfactor
    sigtt = sigtt*wfactor
    siglt = siglt*wfactor
    
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    #sig = sig / 2.0 / pi  # dsig/dtdphicm in microbarns/GeV**2/rad

    wtn = wt_sim * sig / sigcm_sim

    #print("sig",sig)
    #print("sigcm",sigcm_sim)
    #print("wtn",wtn)
    #print("wt_sim",wt_sim)
    
    return [float(wtn),float(sig)]
