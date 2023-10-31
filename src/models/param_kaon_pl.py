#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-31 14:21:32 trottar"
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

    q2_gev = q2_set # Already GeV
    t_gev = t_sim  # Already GeV
    s = w_sim**2
    s_gev = s # Already GeV

    # Calculate tav, ftav, ft
    #tav = (0.0735 + 0.028 * math.log(q2_gev)) * q2_gev
    # RLT (10/8/2023): Testing new tav parameterization
    tav=(0.1112 + 0.0066*math.log(q2_set))*q2_set
    ftav = (abs(t_gev) - tav) / tav
    ft = abs(t_gev) / (abs(t_gev) + mkpl**2)**2

    # Calculate sigl, sigt, siglt, sigtt, sig219, sig
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = params
    #sigl = (p1 + p2 * math.log(q2_sim)) * math.exp((p3 + p4 * math.log(q2_sim)) * (abs(t_gev) - 0.2))
    # RLT (10/12/2023): Removed 0.2 to keep things as simple as possible for initial start parameterization
    sigl = (p1 + p2 * math.log(q2_sim)) * math.exp((p3 + p4 * math.log(q2_sim)) * (abs(t_gev)))
    sigt = p5 + p6 * math.log(q2_sim) + (p7 + p8 * math.log(q2_sim)) * ftav
    siglt = (p9 * math.exp(p10 * abs(t_gev)) + p11 / abs(t_gev)) * math.sin(thetacm_sim)
    sigtt = (p12 * q2_sim * math.exp(-q2_sim)) * ft * math.sin(thetacm_sim)**2

    # RLT (9/25/2023): There are two tav parameterizations in here.
    #                  I am only using the one above, for now.    
    # tav = (-0.178 + 0.315 * math.log(q2_sim)) * q2_sim

    sig219 = (sigt + eps_sim * sigl + eps_sim * math.cos(2. * phicm_sim) * sigtt +
             math.sqrt(2.0 * eps_sim * (1. + eps_sim)) * math.cos(phicm_sim) * siglt) / 1.0

    wfactor = 1.0 / (s_gev - mtar_gev**2)**2
    sig = sig219*wfactor
    sigl = sigl*wfactor
    sigt = sigt*wfactor
    sigtt = sigtt*wfactor
    siglt = siglt*wfactor    
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad

    wtn = wt_sim * sig / sigcm_sim

    print("sig",sig)
    print("sigcm",sigcm_sim)
    print("wtn",wtn)
    print("wt_sim",wt_sim)
    
    return [float(wtn),float(sig)]
