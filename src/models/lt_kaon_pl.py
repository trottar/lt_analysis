#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-17 16:28:52 trottar"
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
# Need to grab epsilon string values from lt_2D_fit script

# First, define empty strings
LOEPS = ""
HIEPS = ""

# Then, set global variables which is called with arguments defined in lt_2D_fit script
def set_val(inp_LOEPS, inp_HIEPS):
    global LOEPS, HIEPS
    LOEPS = inp_LOEPS
    HIEPS = inp_HIEPS

###############################################################################################################################################    
    
# Low epsilon drawing function
def LT_sep_x_lo_fun(x, par):
    eps = float(LOEPS)
    xx = x[0]
    xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx * PI / 180) + eps * par[3] * ROOT.TMath.Cos(2 * xx * PI / 180)
    return xs

###############################################################################################################################################

# High epsilon drawing function
def LT_sep_x_hi_fun(x, par):
    eps = float(HIEPS)
    xx = x[0]
    xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx * PI / 180) + eps * par[3] * ROOT.TMath.Cos(2 * xx * PI / 180)
    return xs

###############################################################################################################################################

# Low epsilon calculating unseparated cross section
def LT_sep_x_lo_fun_unsep(x, par):
    eps = float(LOEPS)
    xx = x[0]
    xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx) + eps * par[3] * ROOT.TMath.Cos(2 * xx)
    return xs

###############################################################################################################################################

# High epsilon calculating unseparated cross section
def LT_sep_x_hi_fun_unsep(x, par):
    eps = float(HIEPS)
    xx = x[0]
    xs = par[0] + eps * par[1] + ROOT.TMath.Sqrt(2 * eps * (1 + eps)) * par[2] * ROOT.TMath.Cos(xx) + eps * par[3] * ROOT.TMath.Cos(2 * xx)
    return xs
