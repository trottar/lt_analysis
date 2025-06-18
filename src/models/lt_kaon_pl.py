#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-06-03 11:38:12 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math
import ROOT

###############################################################################################################################################

# Define constants
PI = math.pi

###############################################################################################################################################    
    
# Low epsilon drawing function
def LT_sep_x_lo_fun_wrapper(inp_eps):
    def LT_sep_x_lo_fun(x, par):
        eps = inp_eps
        xx = x[0]
        #  ρ_LT term = ρₗₜ · √(σ_T·σ_L)  ;  ρ_TT term = ρₜₜ · σ_T
        xs = ( par[0]
              + eps * par[1]
              + ROOT.TMath.Sqrt(2*eps*(1+eps))
                * par[2]
                * ROOT.TMath.Sqrt(par[0] * par[1])
                * ROOT.TMath.Cos(xx)
              + eps
                * par[3]
                * par[0]
                * ROOT.TMath.Cos(2*xx) )
        return xs
    return LT_sep_x_lo_fun
    
###############################################################################################################################################

# High epsilon drawing function
def LT_sep_x_hi_fun_wrapper(inp_eps):
    def LT_sep_x_hi_fun(x, par):
        eps = inp_eps
        xx = x[0]
        #  ρ_LT term = ρₗₜ · √(σ_T·σ_L)  ;  ρ_TT term = ρₜₜ · σ_T
        xs = ( par[0]
              + eps * par[1]
              + ROOT.TMath.Sqrt(2*eps*(1+eps))
                * par[2]
                * ROOT.TMath.Sqrt(par[0] * par[1])
                * ROOT.TMath.Cos(xx)
              + eps
                * par[3]
                * par[0]
                * ROOT.TMath.Cos(2*xx) )
        return xs
    return LT_sep_x_hi_fun

###############################################################################################################################################

# Low epsilon calculating unseparated cross section
def LT_sep_x_lo_fun_unsep_wrapper(inp_eps):
    def LT_sep_x_lo_fun_unsep(x, par):
        eps = inp_eps
        xx = x[0]
        #  ρ_LT term = ρₗₜ · √(σ_T·σ_L)  ;  ρ_TT term = ρₜₜ · σ_T
        xs = ( par[0]
              + eps * par[1]
              + ROOT.TMath.Sqrt(2*eps*(1+eps))
                * par[2]
                * ROOT.TMath.Sqrt(par[0] * par[1])
                * ROOT.TMath.Cos(xx)
              + eps
                * par[3]
                * par[0]
                * ROOT.TMath.Cos(2*xx) )
        return xs
    return LT_sep_x_lo_fun_unsep

###############################################################################################################################################

# High epsilon calculating unseparated cross section
def LT_sep_x_hi_fun_unsep_wrapper(inp_eps):
    def LT_sep_x_hi_fun_unsep(x, par):
        eps = inp_eps
        xx = x[0]
        #  ρ_LT term = ρₗₜ · √(σ_T·σ_L)  ;  ρ_TT term = ρₜₜ · σ_T
        xs = ( par[0]
              + eps * par[1]
              + ROOT.TMath.Sqrt(2*eps*(1+eps))
                * par[2]
                * ROOT.TMath.Sqrt(par[0] * par[1])
                * ROOT.TMath.Cos(xx)
              + eps
                * par[3]
                * par[0]
                * ROOT.TMath.Cos(2*xx) )
        return xs
    return LT_sep_x_hi_fun_unsep
