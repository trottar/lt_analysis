#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-21 19:17:08 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

###############################################################################################################################################

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

###############################################################################################################################################
# Need to grab polarity and Q2 string values from xfit script

# First, define empty strings
pol_str = ""
q2_set = ""

# Then, set global variables which is called with arguments defined in xfit script
def set_val(inp_pol_str, inp_q2_set):
    global pol_str, q2_set
    pol_str = inp_pol_str
    q2_set = inp_q2_set

###############################################################################################################################################

# Function for SigL
def fun_Sig_L(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    #print("Calculating function for func_SigL...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    # RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the extreme slope at high t
    #f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
    f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)+0.2))
    return f

###############################################################################################################################################
    
# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    tav = (0.1112 + 0.0066*math.log(float(q2_set.replace("p","."))))*float(q2_set.replace("p","."))
    ftav = (abs(tt)-tav)/tav
    #print("Calculating function for func_SigT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    # RLT (2/15/2024): Removing t dependence from sigT because it seems
    #                  to be driving poor sep xsects results
    # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
    # RLT (2/21/2024): Reintroducing t-dependence
    # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
    #f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
    #f = par[0]+par[1]*math.log(qq)
    #f = par[0]*math.log(qq)+par[1]/(qq**2)
    f = par[0] / (1 + par[1]*qq)
    return f

###############################################################################################################################################

# Function for SigLT
# thetacm term is defined on function calling
def fun_Sig_LT(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    #print("Calculating function for func_SigLT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))
    return f

###############################################################################################################################################

# Function for SigTT
# thetacm term is defined on function calling
def fun_Sig_TT(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    if pol_str == "pl":
        f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
    #print("Calculating function for func_SigTT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    f = (par[0]*qq*math.exp(-qq))*f_tt
    return f

###############################################################################################################################################
