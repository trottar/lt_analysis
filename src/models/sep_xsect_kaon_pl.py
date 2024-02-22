#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-21 19:19:29 trottar"
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
    
    # Function for SigL
    def sig_L(*par):
        if inp_model == "sigL":
            print("Calculating function for sigL...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            # RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the extreme slope at high t
            #f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
            f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)+0.2))
            return f

    # Function for SigT
    def sig_T(*par):
        if inp_model == "sigT":
            print("Calculating function for sigT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            tav = (0.1112 + 0.0066*math.log(Q2))*Q2
            ftav = (abs(tt)-tav)/tav
            # RLT (2/15/2024): Removing t dependence from sigT because it seems
            #                  to be driving poor sep xsects results
            # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
            # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)            
            #f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
            #f = par[0]+par[1]*math.log(qq)
            #f = par[0]*math.log(qq)+par[1]/(qq**2)
            f = par[0] / (1 + par[1]*qq)
            return f

    # Function for SigLT
    # thetacm term is defined on function calling
    def sig_LT(*par):
        if inp_model == "sigLT":
            print("Calculating function for sigLT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))*math.sin(theta_cm)
            return f

    # Function for SigTT
    # thetacm term is defined on function calling
    def sig_TT(*par):
        if inp_model == "sigTT":
            print("Calculating function for sigTT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
            f = (par[0]*qq*math.exp(-qq))*f_tt*(math.sin(theta_cm)**2)
            return f

    modelDict = {
        "sigL" : sig_L(p1, p2, p3, p4),
        "sigT" : sig_T(p5, p6, p7, p8),
        "sigLT" : sig_LT(p9, p10, p11, p12),
        "sigTT" : sig_TT(p13, p14, p15, p16),
    }

    sig_sep = modelDict[inp_model]

    # Apply weight factor
    g = 1 / ((ww**2) - (m_p**2))**2
    sig_sep = sig_sep*g
    
    print("Model {} = {:.4e}".format(inp_model, sig_sep))
    
    return sig_sep
