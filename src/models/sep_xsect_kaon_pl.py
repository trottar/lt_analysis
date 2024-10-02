#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-10-01 22:54:15 trottar"
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
            # RLT (3/09/2024): Removing +0.2 term for better parameterization of Q2=3.0, W=2.32
            #f = (par[0]+par[1]*math.log(qq)) * math.exp((par[2]+par[3]*math.log(qq)) * (abs(tt)))
            try:
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
                ##                
                # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
                #                  that incorporates Q2-dep based of pi FF
                ft = abs(tt) / (abs(tt) + mkpl**2)**2 # pole term
                Qdep_L=qq/(1.0+(1.77*qq)+0.12*(qq**2))
                ##f=(par[0]*Qdep_L*ft)*math.exp(-par[1]*(abs(tt)))
                f=(par[0]*Qdep_L*ft)*math.exp(-par[1]*(abs(tt)))
                
            except ValueError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigL, setting to zero for this event...")
            except OverflowError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigL, setting to zero for this event...")
            return f

    # Function for SigT
    def sig_T(*par):
        if inp_model == "sigT":
            print("Calculating function for sigT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            tav = (0.1112 + 0.0066*math.log(Q2))*Q2
            ftav = (abs(tt)-tav)/tav
            f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
            try:
                # RLT (2/15/2024): Removing t dependence from sigT because it seems
                #                  to be driving poor sep xsects results
                # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
                # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
                #
                #f = par[0]+par[1]*math.log(qq)+(par[2]+par[3]*math.log(qq)) * ftav
                #f = par[0]+par[1]*math.log(qq)
                #f = par[0]*math.log(qq)+par[1]/(qq**2)
                #f = par[0] / (1 + par[1]*qq)
                # RLT (4/20/2024): Adding in t-dependence
                #f = (par[0] / (1 + par[1]*qq)) * ftav
                #f = (par[0] / (1 + par[1]*qq)) * abs(tt)
                # RLT (4/20/2024): Exponential t-dependence
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
                
            except ValueError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigT, setting to zero for this event...")
            except OverflowError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigT, setting to zero for this event...")
            return f

    # Function for SigLT
    # thetacm term is defined on function calling
    def sig_LT(*par):
        if inp_model == "sigLT":
            print("Calculating function for sigLT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            try:
                f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
                ##
                f = (par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))*math.sin(theta_cm)
                #f = (par[0]+par[2]/abs(tt))*math.sin(theta_cm)
                # RLT (4/23/2024): Marco's thesis functional forms
                #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + (qq**2)*par[2]))
                ##
                # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
                #                  that incorporates Q2-dep based of pi FF
                ###f=(par[0]/(1+qq))*math.sin(theta_cm)*math.exp(-par[1]*(abs(tt)))
                #f=(par[0]*math.exp(par[1]*abs(tt))+par[2]/abs(tt))*math.sin(theta_cm)
                
            except ValueError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigLT, setting to zero for this event...")
            except OverflowError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigLT, setting to zero for this event...")
            return f

    # Function for SigTT
    # thetacm term is defined on function calling
    def sig_TT(*par):
        if inp_model == "sigTT":
            print("Calculating function for sigTT...\nQ2={:.4e}, t={:.4e}\npar=({:.4e}, {:.4e}, {:.4e}, {:.4e})".format(qq, tt, *par))
            try:
                #  RLT (7/11/2024): Moved below for Q2dep func form
                f_tt=abs(tt)/(abs(tt)+mkpl**2)**2 # pole factor
                ##
                f = (par[0]*qq*math.exp(-qq))*f_tt*(math.sin(theta_cm)**2)
                # RLT (4/23/2024): Marco's thesis functional forms
                #f = par[0] * math.exp(-par[1]*abs(tt)) * (1.0 / (1 + (qq**2)*par[2]))
                ##
                # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
                #                  that incorporates Q2-dep based of pi FF
                ##f=(-par[0]/(1+qq))*(math.sin(theta_cm)**2)*math.exp(-par[1]*(abs(tt)))
                ###f=(par[0]/(1+qq))*(math.sin(theta_cm)**2)*f_tt*math.exp(-par[1]*(qq))
                ####f = ((-par[0]*abs(tt)+par[1])*(abs(tt)**(qq/par[2]))-par[3]*qq)*(math.sin(theta_cm)**2)
                
            except ValueError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigTT, setting to zero for this event...")
            except OverflowError:
                f = -1000.0
                #print("WARNING: Overflowerror on sigTT, setting to zero for this event...")                
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
    ##g = 1 / ((ww**2) - (m_p**2))**2.25
    sig_sep = sig_sep*g

    sig_sep = sig_sep/2.0/math.pi
    
    # Convert from ub/GeV**2 to nb/GeV**2
    sig_sep = sig_sep*1e3
    
    print("Model {} = {:.4e}".format(inp_model, sig_sep))
    
    return sig_sep
