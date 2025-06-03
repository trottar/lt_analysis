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
        xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx * PI / 180) + eps * par[3] * ROOT.TMath.Cos(2 * xx * PI / 180)
        return xs
    return LT_sep_x_lo_fun
    
###############################################################################################################################################

# High epsilon drawing function
def LT_sep_x_hi_fun_wrapper(inp_eps):
    def LT_sep_x_hi_fun(x, par):
        eps = inp_eps
        xx = x[0]
        xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx * PI / 180) + eps * par[3] * ROOT.TMath.Cos(2 * xx * PI / 180)
        return xs
    return LT_sep_x_hi_fun

###############################################################################################################################################

# Low epsilon calculating unseparated cross section
def LT_sep_x_lo_fun_unsep_wrapper(inp_eps):
    def LT_sep_x_lo_fun_unsep(x, par):
        eps = inp_eps
        xx = x[0]
        xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx) + eps * par[3] * ROOT.TMath.Cos(2 * xx)
        return xs
    return LT_sep_x_lo_fun_unsep

###############################################################################################################################################

# High epsilon calculating unseparated cross section
def LT_sep_x_hi_fun_unsep_wrapper(inp_eps):
    def LT_sep_x_hi_fun_unsep(x, par):
        eps = inp_eps
        xx = x[0]
        xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx) + eps * par[3] * ROOT.TMath.Cos(2 * xx)
        return xs
    return LT_sep_x_hi_fun_unsep
