#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-10-01 08:44:00 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import ROOT
import shutil
import os, sys

##################################################################################################################################################
kinematics = sys.argv[1].split("_")
W = sys.argv[2]
Q2 = sys.argv[3]
tmin = float(sys.argv[4])
tmax = float(sys.argv[5])
ParticleType = sys.argv[6]
POL = sys.argv[7]
OutFilename = sys.argv[8]
formatted_date = sys.argv[9]
inp_debug =  sys.argv[10]

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import show_pdf_with_evince, create_dir

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=10:
    print("!!!!! ERROR !!!!!\n Expected 12 arguments\n Usage is with - KIN W Q2 TMIN TMAX ParticleType POL OutFilename formatted_date inp_debug\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    


###############################################################################################################################################
# ltsep package import and pathing definitions

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

TEMP_CACHEPATH=f"{OUTPATH}/testing_env"

output_file_lst = []

# Append csv of efficiencies
output_file_lst.append(TEMP_CACHEPATH + "/" + f"table_{ParticleType}_{kinematics}" + ".csv")

###############################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

if inp_debug == "False":
    DEBUG = False # Flag for no plot splash
else:
    DEBUG = True # Flag for plot splash

if int(POL) == 1:
    pol_str = "pl"
elif int(POL) == -1:
    pol_str = "mn"
else:
    print("ERROR: Invalid polarity...must be +1 or -1")
    sys.exit(2)    

closest_date = formatted_date
print("\n\n")
print("="*50)
iter_num = 1 # Default iteration number
print("\n\tIteration number", iter_num)
print("="*50)
print("\n\n")    

# Create a new directory for each iteration in cache
new_dir = "{}/{}/Q{}W{}/{}".format(TEMP_CACHEPATH, ParticleType.lower(), Q2, W, closest_date)
create_dir(new_dir)

TEMP_CACHEPATH = new_dir

inpDict = {
    "kinematics" : kinematics,
    "W" : W,
    "Q2" : Q2,
    "OutFilename" : OutFilename,
    "tmin" : tmin,
    "tmax" : tmax,   
    "ParticleType" : ParticleType,
    "POL" : POL,
    "formatted_date" : formatted_date,
    "TEMP_CACHEPATH" : TEMP_CACHEPATH
}

# Add closest and formatted dates to inpDict (used in plot comparison)
inpDict["closest_date"] = closest_date
inpDict["formatted_date"] = formatted_date
inpDict["iter_num"] = iter_num

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

##############
# HARD CODED #
##############
# May need to adjust these for diamond plots to work
if Q2Val == 2.1:
    inpDict["Q2min"] = 2.115 - (1/4)*2.115
    inpDict["Q2max"] = 2.115 + (1/4)*2.115
    inpDict["Wmin"] = WVal - (1/10)*WVal
    inpDict["Wmax"] = WVal + (1/10)*WVal
elif Q2Val == 0.4: # Q2=0.38, pion
    inpDict["Q2min"] = 0.385 - (1/4)*0.385
    inpDict["Q2max"] = 0.385 + (1/4)*0.385
    inpDict["Wmin"] = WVal - (1/9)*WVal
    inpDict["Wmax"] = WVal + (1/9)*WVal
else:
    inpDict["Q2min"] = Q2Val - (1/4)*Q2Val
    inpDict["Q2max"] = Q2Val + (1/4)*Q2Val
    inpDict["Wmin"] = WVal - (1/10)*WVal
    inpDict["Wmax"] = WVal + (1/10)*WVal
##############
##############

##############
# Save input model
output_file_lst.append('functions/Q{}W{}.model'.format(Q2, W))

# Save python script that contain separated xsect models for xfit script
py_xfit = 'models/xfit_{}_{}.py'.format(ParticleType, pol_str)
output_file_lst.append(py_xfit)

# Active scripts to make file selection dynamic
# Needs to be done this way because of fortran compiler limitations
py_xfit_active = 'models/xfit_active.py'
# Copying content of used models to actively used files
print("\nCopying {} to {}".format(LTANAPATH+"/testing/"+py_xfit, LTANAPATH+"/testing/"+py_xfit_active))
shutil.copy(LTANAPATH+"/testing/"+py_xfit, LTANAPATH+"/testing/"+py_xfit_active)

# Run weight iteration script for optimizing parameterization
sys.path.append("models")
from xfit_in_t import x_fit_in_t
x_fit_in_t(ParticleType, pol_str, closest_date, Q2, W, inpDict, output_file_lst)
if DEBUG:
    show_pdf_with_evince(TEMP_CACHEPATH+"/{}_xfit_in_t_Q{}W{}*.pdf".format(ParticleType, Q2, W))

# ***Parameter file for new iteration!***
# ***These parameters are newly generated for this iteration above. See README for more info on procedure!***
new_param_file = '{}/testing/parameters/new_par.{}_Q{}W{}.dat'.format(LTANAPATH, pol_str, Q2.replace("p",""), W.replace("p",""))
print("\nCopying {} to {}".format(new_param_file, new_param_file.replace(LTANAPATH+"/testing/parameters", TEMP_CACHEPATH)))
shutil.copy(new_param_file, new_param_file.replace(LTANAPATH+"/testing/parameters", TEMP_CACHEPATH))