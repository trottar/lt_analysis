#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-05 16:55:38 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

def import_model(inp_model, arg_str):

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    Q2, theta_cm, tt, qq, ww, *params = args

    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params

    g = 1 / ((ww**2) - (m_p**2))**2
    
    # Function for SigL
    def sig_L(*par):
        print("\n\nCalculating function for sigL...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
        f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
        print("sigL = {:.4e}".format(f))
        return f

    # Function for SigT
    def sig_T(*par):
        tav = (0.1112 + 0.0066*math.log(Q2))*Q2
        ftav = (abs(tt)-tav)/tav
        print("\n\nCalculating function for sigT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
        f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
        print("sigT = {:.4e}".format(f))
        return f

    # Function for SigLT
    # thetacm term is defined on function calling
    def sig_LT(*par):
        print("\n\nCalculating function for sigLT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
        f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))*math.sin(theta_cm*PI/180)
        print("sigLT = {:.4e}".format(f))
        return f

    # Function for SigTT
    # thetacm term is defined on function calling
    def sig_TT(*par):
        f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
        print("\n\nCalculating function for sigTT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
        f = (par[0]*qq*math.exp(-qq))*f_tt*(math.sin(theta_cm*PI/180)**2)
        print("sigTT = {:.4e}".format(f))
        return f

    modelDict = {
        "sigL" : sig_L(p1, p2, p3, p4),
        "sigT" : sig_T(p5, p6, p7, p8),
        "sigLT" : sig_LT(p9, p10, p11, p12),
        "sigTT" : sig_TT(p13, p14, p15, p16),
    }

    sig = modelDict[inp_model]
    sig = sig*g
    
    return sig/2./PI/1e6
