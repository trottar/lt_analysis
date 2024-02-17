#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-17 00:17:00 trottar"
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

###############################################################################################################################################

# Function for SigL
def fun_Sig_L(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    #print("Calculating function for func_SigL...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
    return f

###############################################################################################################################################

# HERE!!!! Need to incorporate q2_set

q2_set = 0.0

def q2_set_val(inp_q2_set):
    global q2_set
    q2_set = inp_q2_set

# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    qq = abs(x[1])
    tav = (0.1112 + 0.0066*math.log(float(q2_set.replace("p","."))))*float(q2_set.replace("p","."))
    ftav = (abs(tt)-tav)/tav
    #print("Calculating function for func_SigT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
    # RLT (2/15/2024): Removing t dependence from sigT because it seems
    #                  to be driving poor sep xsects results    
    #f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
    f = par[0]+par[1]*math.log(qq)
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
