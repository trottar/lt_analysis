#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-03 11:17:59 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
from collections import defaultdict
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=25:
    print("!!!!! ERROR !!!!!\n Expected 25 arguments\n Usage is with - KIN W Q2 EPSVAL OutDATAFilename OutDUMMYFilename OutFullAnalysisFilename tmin tmax NumtBins NumPhiBins runNumRight runNumLeft runNumCenter data_charge_right data_charge_left data_charge_center dummy_charge_right dummy_charge_left dummy_charge_center InData_efficiency_right InData_efficiency_left InData_efficiency_center efficiency_table ParticleType\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    

DEBUG = False # Flag for no cut plots

# Input params
kinematics = sys.argv[1].split("_")
W = sys.argv[2]
Q2 = sys.argv[3]
EPSVAL = sys.argv[4]
InDATAFilename = sys.argv[5]
InDUMMYFilename = sys.argv[6]
OutFilename = sys.argv[7]
tmin = float(sys.argv[8])
tmax = float(sys.argv[9])
NumtBins = int(sys.argv[10])
NumPhiBins = int(sys.argv[11])
runNumRight = sys.argv[12]
runNumLeft = sys.argv[13]
runNumCenter = sys.argv[14]
data_charge_right = int(sys.argv[15])/1000 # Convert from uC to C
data_charge_left = int(sys.argv[16])/1000 # Convert from uC to C
data_charge_center = int(sys.argv[17])/1000 # Convert from uC to C
dummy_charge_right = int(sys.argv[18])/1000 # Convert from uC to C
dummy_charge_left = int(sys.argv[19])/1000 # Convert from uC to C
dummy_charge_center = int(sys.argv[20])/1000 # Convert from uC to C
InData_efficiency_right = sys.argv[21]
InData_efficiency_left = sys.argv[22]
InData_efficiency_center = sys.argv[23]
efficiency_table = sys.argv[24]
ParticleType = sys.argv[25]

inpDict = {
    "kinematics" : kinematics,
    "W" : W,
    "Q2" : Q2,
    "EPSVAL" : EPSVAL,
    "InDATAFilename" : InDATAFilename,
    "InDUMMYFilename" : InDUMMYFilename,
    "OutFilename" : OutFilename,
    "tmin" : tmin,
    "tmax" : tmax,
    "NumtBins" : NumtBins,
    "NumPhiBins" : NumPhiBins,
    "runNumRight" : runNumRight,
    "runNumLeft" : runNumLeft,
    "runNumCenter" : runNumCenter,
    "data_charge_right" : data_charge_right,
    "data_charge_left" : data_charge_left,
    "data_charge_center" : data_charge_center,
    "dummy_charge_right" : dummy_charge_right,
    "dummy_charge_left" : dummy_charge_left,
    "dummy_charge_center" : dummy_charge_center,
    "InData_efficiency_right" : InData_efficiency_right,
    "InData_efficiency_left" : InData_efficiency_left,
    "InData_efficiency_center" : InData_efficiency_center,
    "efficiency_table" : efficiency_table,
    "ParticleType" : ParticleType,
}

###############################################################################################################################################
# ltsep package import and pathing definitions

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

################################################################################################################################################    

##############################
# Step 1 of the lt_analysis: #
##############################
'''
* All analysis cuts per run (i.e. PID, acceptance, timing) are applied. 
* This should have been completed before this script using...

> lt_analysis/scripts/Prod/analyze_prod.py

* This is script is called when using the '-a' flag in...

> lt_analysis/run_Prod_Analysis.sh
'''

##############################
# Step 2 of the lt_analysis: #
##############################
'''
* Diamond cuts are drawn on the center setting.

* Make sure to check that high and low eps overlap in...

> lt_analysis/OUTPUT/Analysis/<ANATYPE>LT/<KIN>_<SETTING>_Diamond_Cut.pdf

** TODO: Add individual setting diamond plots
'''

#Importing diamond cut script
from diamond import DiamondPlot

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

Q2min = Q2Val - (2/7)*Q2Val # Minimum value of Q2 on the Q2 vs W plot
Q2max = Q2Val + (2/7)*Q2Val # Maximum value of Q2 on the Q2 vs W plot
Wmin = WVal - (2/7)*WVal # min y-range for Q2vsW plot
Wmax = WVal + (2/7)*WVal # max y-range for Q2vsW plot

phisetlist = ["Center","Left","Right"]
for phiset in phisetlist:
    # Call diamond cut script and append paramters to dictionary
    inpDict.update(DiamondPlot(ParticleType, Q2Val, Q2min, Q2max, WVal, Wmin, Wmax, phiset, tmin, tmax))

##############################
# Step 3 of the lt_analysis: #
##############################
'''
Apply random subtraction to data and dummy
'''

from rand_sub import rand_sub

# Call histogram function above to define dictonaries for right, left, center settings
# Put these all into an array so that if we are missing a setting it is easier to remove
# Plus it makes the code below less repetitive
phisetlist = ["Center","Left","Right"]
histlist = []
for phiset in phisetlist:
    histlist.append(rand_sub(phiset,inpDict))

print("\n\n")

settingList = []
for i,hist in enumerate(histlist):    
    if not bool(hist): # If hist is empty
        histlist.remove(hist)
    else:
        settingList.append(hist["phi_setting"])
    
##############################
# Step 4 of the lt_analysis: #
##############################
'''
* Combine all settings and choose t/phi bins for low eps.

* These bins will also be used of high eps, so check high eps as well.
'''

#from find_bins import find_bins

##############################
# Step 5 of the lt_analysis: #
##############################
'''
* Compare SIMC to data/dummy setting by setting.

* SIMC weight may need to be recalculated.

* If so, use script...

> lt_analysis/src/SIMC/??????????????.f
'''

#from compare_simc import compare_simc

##############################
# Step 6 of the lt_analysis: #
##############################
'''
* Combine settings again at high and low eps for each 
  Q2,W,tbin,phibin,eps for data, dummy, and SIMC.

* Dummy is subtracted from data bin by bin.
* The yield is calculated using the effective charge from data and 
  the normfactor/nevents.

* The data and SIMC yields are compared and the R value per bin is obtained.
'''

#from calculate_yield import calculate_yield

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data and SIMC.

* Find the mean data values of W,Q2,theta,eps for each t bin of high and low epsilon
  for both data and SIMC.
'''
