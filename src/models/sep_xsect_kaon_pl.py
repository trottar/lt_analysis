#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-05 15:51:03 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

# Define constants
PI = ROOT.TMath.Pi()
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677


def import_model(inp_model, Q2, theta_cm, tt, qq, ww, par):
        
    # Function for SigL
    def sig_L(tt, qq, ww, par):
        g = 1 / ((ww**2) - (m_p**2))**2
        print("Calculating function for sigL...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
        return f

    # Function for SigT
    def sig_T(Q2, tt, qq, ww, par):
        g = 1 / ((ww**2) - (m_p**2))**2
        tav = (0.1112 + 0.0066*math.log(float(Q2.replace("p","."))))*float(Q2.replace("p","."))
        ftav = (abs(tt)-tav)/tav
        print("Calculating function for sigT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
        return f

    # Function for SigLT
    # thetacm term is defined on function calling
    def sig_LT(theta_cm, tt, qq, ww, par):
        g = 1 / ((ww**2) - (m_p**2))**2
        print("Calculating function for sigLT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))*math.sin(theta_cm)
        return f

    # Function for SigTT
    # thetacm term is defined on function calling
    def sig_TT(theta_cm, tt, qq, ww, par):
        g = 1 / ((ww**2) - (m_p**2))**2
        f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
        print("Calculating function for sigTT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        f = (par[0]*qq*math.exp(-qq))*f_tt*(math.sin(theta_cm)**2)
        return f

    modelDict = {
        "sigL" : sig_L(tt, qq, ww, par),
        "sigT" : sig_T(Q2, tt, qq, ww, par),
        "sigLT" : sig_LT(theta_cm, tt, qq, ww, par),
        "sigTT" : sig_TT(theta_cm, tt, qq, ww, par),
    }

    modelDict[inp_model]
