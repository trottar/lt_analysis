#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-03 10:35:21 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages

###############################################################################################################################################

'''
ltsep package import and pathing definitions
'''

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

# Define various inputs 
kinematics = inpDict["kinematics"] 
W = inpDict["W"] 
Q2 = inpDict["Q2"] 
EPSVAL = inpDict["EPSVAL"] 
InDATAFilename = inpDict["InDATAFilename"] 
InDUMMYFilename = inpDict["InDUMMYFilename"] 
OutFilename = inpDict["OutFilename"] 
tmin = inpDict["tmin"] 
tmax = inpDict["tmax"] 
NumtBins = inpDict["NumtBins"] 
NumPhiBins = inpDict["NumPhiBins"] 
runNumRight = inpDict["runNumRight"] 
runNumLeft = inpDict["runNumLeft"] 
runNumCenter = inpDict["runNumCenter"]
data_charge_right = inpDict["data_charge_right"] 
data_charge_left = inpDict["data_charge_left"] 
data_charge_center = inpDict["data_charge_center"] 
dummy_charge_right = inpDict["dummy_charge_right"] 
dummy_charge_left = inpDict["dummy_charge_left"] 
dummy_charge_center = inpDict["dummy_charge_center"] 
InData_efficiency_right = inpDict["InData_efficiency_right"] 
InData_efficiency_left = inpDict["InData_efficiency_left"] 
InData_efficiency_center = inpDict["InData_efficiency_center"] 
efficiency_table = inpDict["efficiency_table"] 
ParticleType = inpDict["ParticleType"]

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
'''

#Importing diamond cut script
from diamond import DiamondPlot

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

Q2min = Q2Val - (2/7)*Q2Val # Minimum value of Q2 on the Q2 vs W plot
Q2max = Q2Val + (2/7)*Q2Val # Maximum value of Q2 on the Q2 vs W plot
Wmin = WVal - (2/7)*WVal # min y-range for Q2vsW plot
Wmax = WVal + (2/7)*WVal # max y-range for Q2vsW plot

# Call diamond cut script
paramDict = DiamondPlot(ParticleType, Q2Val, Q2min, Q2max, WVal, Wmin, Wmax, phi_setting, tmin, tmax)

# Define diamond cut parameters
a1 = paramDict["a1"]
b1 = paramDict["b1"]
a2 = paramDict["a2"]
b2 = paramDict["b2"]
a3 = paramDict["a3"]
b3 = paramDict["b3"]
a4 = paramDict["a4"]
b4 = paramDict["b4"]

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
    histlist.append(defineHists(phiset,inpDict))

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

from find_bins import find_bins

##############################
# Step 5 of the lt_analysis: #
##############################
'''
* Compare SIMC to data/dummy setting by setting.

* SIMC weight may need to be recalculated.

* If so, use script...

> lt_analysis/src/SIMC/??????????????.f
'''

from compare_simc import compare_simc

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

from calculate_yield import calculate_yield

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data and SIMC.

* Find the mean data values of W,Q2,theta,eps for each t bin of high and low epsilon
  for both data and SIMC.
'''
