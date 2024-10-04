#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-04 17:51:06 trottar"
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
    qq = float(q2_set.replace("p","."))
    #qq = abs(x[1])
    try:
        #print("Calculating function for func_SigL...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        # RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the extreme slope at high t
        # RLT (3/09/2024): Removing +0.2 term for better parameterization of Q2=3.0, W=2.32
        #
        #f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
        #f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)+0.2))
        # RLT (4/23/2024): Marco's thesis functional forms
        #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + par[2]*qq))
        # RLT (6/04/2024): Testing simplier exp form for L+T
        ##
        ##f = (par[0]+par[1]*math.log(qq)) * math.exp(par[2] * (abs(tt)))
        #f = (par[0] * ((abs(tt)/qq)-1)) * math.exp(par[1] * (abs(tt)))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        f_tt = abs(tt) / (abs(tt) + mkpl**2)**2 # pole term
        Qdep_L=qq/(1.0+(1.77*qq)+0.12*(qq**2))
        ##f=(par[0]*Qdep_L*f_tt)*math.exp(-par[1]*(abs(tt)))
        f=(par[0]*Qdep_L*f_tt)*math.exp(-par[1]*(abs(tt)))
        
    except OverflowError:
        f = -1000.0
        #print("WARNING: Overflowerror on sigL, setting to zero for this event...")        
    return f

###############################################################################################################################################
    
# Function for SigT
def fun_Sig_T(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    #qq = abs(x[1])
    tav = (0.1112 + 0.0066*math.log(float(q2_set.replace("p","."))))*float(q2_set.replace("p","."))
    ftav = (abs(tt)-tav)/tav
    f_tt = abs(tt) / (abs(tt) + mkpl**2)**2 # pole term
    try:
        #print("Calculating function for func_SigT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        # RLT (2/15/2024): Removing t dependence from sigT because it seems
        #                  to be driving poor sep xsects results
        # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
        # RLT (2/21/2024): Reintroducing t-dependence
        # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
        #
        #f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
        #f = par[0]+par[1]*math.log(qq)
        #f = par[0]*math.log(qq)+par[1]/(qq**2)
        #f = par[0] / (1 + par[1]*qq)
        # RLT (4/20/2024): Adding in t-dependence
        #f = (par[0] / (1 + par[1]*qq)) * ftav
        #f = (par[0] / (1 + par[1]*qq)) * abs(tt)
        # RLT (4/23/2024): Exponential t-dependence
        #f = (par[0] / (1 + par[1]*qq)) * math.exp(par[2]*abs(tt))
        # RLT (4/23/2024): Marco's thesis functional forms
        #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + par[2]*qq))
        # RLT (6/04/2024): Testing simplier exp form for L+T
        ##
        ##f = (par[0] * ((abs(tt)/qq)-1)) * math.exp(par[1] * (abs(tt)))
        #f = (par[0]+par[1]*math.log(qq)) * math.exp(par[2] * (abs(tt)))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        #f=(par[0]/qq)*math.exp(-par[1]*(qq**2))
        Qdep_T=(math.exp(-qq**2))/qq
        ##f=par[0]*(par[1]+math.exp(-par[2]*(abs(tt))))*(Qdep_T**par[3])
        f=(par[0]*math.exp(-par[1]*(abs(tt)))+par[2]*(abs(tt)))*(Qdep_T**par[3])
        
    except OverflowError:
        f = -1000.0
        #print("WARNING: Overflowerror on sigT, setting to zero for this event...")        
    return f

###############################################################################################################################################

# Function for SigLT
# thetacm term is defined on function calling
def fun_Sig_LT(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    #qq = abs(x[1])
    f_tt = abs(tt) / (abs(tt) + mkpl**2)**2 # pole term
    try:
        #print("Calculating function for func_SigLT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        ##
        f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))
        #f = (par[0]+par[2]/abs(tt))
        # RLT (4/23/2024): Marco's thesis functional forms
        #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + (qq**2)*par[2]))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        ###f=(par[0]/(1+qq))*math.exp(-par[1]*(abs(tt)))
        #f=(par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))
        
    except OverflowError:
        f = -1000.0
        #print("WARNING: Overflowerror on sigLT, setting to zero for this event...")        
    return f

###############################################################################################################################################

# Function for SigTT
# thetacm term is defined on function calling
def fun_Sig_TT(x, par):
    tt = abs(x[0])
    qq = float(q2_set.replace("p","."))
    #qq = abs(x[1])
    f_tt = abs(tt) / (abs(tt) + mkpl**2)**2 # pole term
    if pol_str == "pl":
        f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
    try:
        #print("Calculating function for func_SigTT...\nQ2={:.1e}, t={:.3e}\npar=({:.2e}, {:.2e}, {:.2e}, {:.2e})\n\n".format(qq, tt, *par))
        ##
        ##f = (par[0]*qq*math.exp(-qq))*f_tt
        # RLT (4/23/2024): Marco's thesis functional forms
        #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + (qq**2)*par[2]))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        #f=(-par[0]/(1+qq))*math.exp(-par[1]*(abs(tt)))
        ###f=(par[0]/(1+qq))*f_tt*math.exp(-par[1]*(qq))
        f = ((-par[0]*abs(tt)+par[1])*(abs(tt)**(qq/par[2]))-par[3]*qq)
        
    except OverflowError:
        f = -1000.0
        #print("WARNING: Overflowerror on sigTT, setting to zero for this event...")
    return f

###############################################################################################################################################
