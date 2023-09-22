#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-21 20:51:22 trottar"
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

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, q2_sim, w_sim, t_sim, eps_sim, thetacm_sim, phicm_sim, sigcm_sim, wt_sim, *params = args

    q2_gev = q2_set / 1e6
    t_gev = t_sim
    s = w_sim**2
    s_gev = s / 1e6

    # Calculate tav, ftav, ft
    tav = (0.0735 + 0.028 * math.log(q2_gev)) * q2_gev
    ftav = (abs(t_gev) - tav) / tav
    ft = abs(t_gev) / (abs(t_gev) + 0.139570**2)**2

    # Calculate sigl, sigt, siglt, sigtt, sig219, sig
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = params
    sigl = (p1 + p2 * math.log(q2_sim)) * math.exp((p3 + p4 * math.log(q2_sim)) * (abs(t_gev) - 0.2))
    sigt = p5 + p6 * math.log(q2_sim) + (p7 + p8 * math.log(q2_sim)) * ftav
    siglt = (p9 * math.exp(p10 * abs(t_gev)) + p11 / abs(t_gev)) * math.sin(thetacm_sim)
    sigtt = (p12 * q2_sim * math.exp(-q2_sim)) * ft * math.sin(thetacm_sim)**2

    tav = (-0.178 + 0.315 * math.log(q2_sim)) * q2_sim

    sig219 = (sigt + eps_sim * sigl + eps_sim * math.cos(2. * phicm_sim) * sigtt +
             math.sqrt(2.0 * eps_sim * (1. + eps_sim)) * math.cos(phicm_sim) * siglt) / 1.0

    wfactor = 1.0 / (s_gev - mtar_gev**2)**2
    sig = sig219 * wfactor
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad

    wtn = wt_sim * sig / sigcm_sim

    if wtn > 0.0:
        pass
    else:
        wtn = 0.0

    return float(wtn)
