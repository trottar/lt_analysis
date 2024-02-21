#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-21 13:27:58 trottar"
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
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce
import csv
import json
import shutil

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import show_pdf_with_evince, create_dir, is_root_obj, is_hist, hist_to_root, custom_encoder, set_dynamic_axis_ranges, notify_email

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=45:
    print("!!!!! ERROR !!!!!\n Expected 45 arguments\n Usage is with - KIN W Q2 LOEPS HIEPS OutDATAFilename OutDUMMYFilename OutFullAnalysisFilename tmin tmax NumtBins NumPhiBins runNumRight runNumLeft runNumCenter data_charge_right data_charge_left data_charge_center dummy_charge_right dummy_charge_left dummy_charge_center data_charge_err_right data_charge_err_left data_charge_err_center dummy_charge_err_right dummy_charge_err_left dummy_charge_err_center InData_efficiency_right InData_efficiency_left InData_efficiency_center InData_error_efficiency_right InData_error_efficiency_left InData_error_efficiency_center efficiency_table ParticleType EPSSET pThetaValRight pThetaValLeft pThetaValCenter EbeamValRight EbeamValLeft EbeamValCenter POL formatted_date inp_debug\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    

# Input params
kinematics = sys.argv[1].split("_")
W = sys.argv[2]
Q2 = sys.argv[3]
LOEPS = sys.argv[4]
HIEPS = sys.argv[5]
InDATAFilename = sys.argv[6]
InDUMMYFilename = sys.argv[7]
OutFilename = sys.argv[8]
tmin = float(sys.argv[9])
tmax = float(sys.argv[10])
NumtBins = int(sys.argv[11])
NumPhiBins = int(sys.argv[12])
runNumRight = sys.argv[13]
runNumLeft = sys.argv[14]
runNumCenter = sys.argv[15]
data_charge_right = float(sys.argv[16])
data_charge_left = float(sys.argv[17])
data_charge_center = float(sys.argv[18])
dummy_charge_right = float(sys.argv[19])
dummy_charge_left = float(sys.argv[20])
dummy_charge_center = float(sys.argv[21])
data_charge_err_right = float(sys.argv[22])
data_charge_err_left = float(sys.argv[23])
data_charge_err_center = float(sys.argv[24])
dummy_charge_err_right = float(sys.argv[25])
dummy_charge_err_left = float(sys.argv[26])
dummy_charge_err_center = float(sys.argv[27])
InData_efficiency_right = sys.argv[28]
InData_efficiency_left = sys.argv[29]
InData_efficiency_center = sys.argv[30]
InData_error_efficiency_right = sys.argv[31]
InData_error_efficiency_left = sys.argv[32]
InData_error_efficiency_center = sys.argv[33]
efficiency_table = sys.argv[34]
ParticleType = sys.argv[35]
EPSSET = sys.argv[36]
pThetaValRight = list(sys.argv[37].split(" "))
pThetaValLeft = list(sys.argv[38].split(" "))
pThetaValCenter = list(sys.argv[39].split(" "))
EbeamValRight = list(sys.argv[40].split(" "))
EbeamValLeft = list(sys.argv[41].split(" "))
EbeamValCenter = list(sys.argv[42].split(" "))
POL = sys.argv[43]
formatted_date = sys.argv[44]
inp_debug =  sys.argv[45]

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

if EPSSET == "low":
    EPSVAL = LOEPS
else:
    EPSVAL = HIEPS
    
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
    "data_charge_err_right" : data_charge_err_right,
    "data_charge_err_left" : data_charge_err_left,
    "data_charge_err_center" : data_charge_err_center,
    "dummy_charge_err_right" : dummy_charge_err_right,
    "dummy_charge_err_left" : dummy_charge_err_left,
    "dummy_charge_err_center" : dummy_charge_err_center,    
    "InData_efficiency_right" : InData_efficiency_right,
    "InData_efficiency_left" : InData_efficiency_left,
    "InData_efficiency_center" : InData_efficiency_center,
    "InData_error_efficiency_right" : InData_error_efficiency_right,
    "InData_error_efficiency_left" : InData_error_efficiency_left,
    "InData_error_efficiency_center" : InData_error_efficiency_center,    
    "efficiency_table" : efficiency_table,
    "ParticleType" : ParticleType,
    "EPSSET" : EPSSET,
    "pThetaValRight" : pThetaValRight,
    "pThetaValLeft" : pThetaValLeft,
    "pThetaValCenter" : pThetaValCenter,
    "EbeamValRight" : EbeamValRight,
    "EbeamValLeft" : EbeamValLeft,
    "EbeamValCenter" : EbeamValCenter,
    "POL" : POL,
    "formatted_date" : formatted_date,
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
CACHEPATH=lt.CACHEPATH

foutroot = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
foutjson  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".json"
outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

output_file_lst = []

###############################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Removes this file to reset iteration count (see below for more details)
f_path = "{}/{}_Q{}W{}_iter.dat".format(LTANAPATH,ParticleType,Q2,W)
# Check if the file exists
if os.path.exists(f_path):
    os.remove(f_path)

# Create a new directory for each iteration in cache
new_dir = "{}/{}/{}/Q{}W{}".format(CACHEPATH, USER, ParticleType.lower(), Q2, W)
create_dir(new_dir)
new_dir = "{}/{}/{}/Q{}W{}/{}".format(CACHEPATH, USER, ParticleType.lower(), Q2, W, formatted_date)
create_dir(new_dir)
    
##############################
# Step 1 of the lt_analysis: # DONE
##############################
'''
* All analysis cuts per run (i.e. PID, acceptance, timing) are applied. 
* This should have been completed before this script using...

> lt_analysis/init/analyze_prod.py

* This is script is called when using the '-a' flag in...

> lt_analysis/run_Prod_Analysis.sh
'''

##############################
# Step 2 of the lt_analysis: # DONE
##############################
'''
* Diamond cuts are drawn on the CENTER setting.

* Make sure to check that high and low eps overlap in...

> lt_analysis/OUTPUT/Analysis/<ANATYPE>LT/<ParticleType>_<KIN>_Center_Diamond_Cut.pdf

'''

# Importing diamond cut script
sys.path.append("cuts")
from diamond import DiamondPlot

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

##############
# HARD CODED #
##############
# May need to adjust these for diamond plots to work
# 2/7 seems to work for most Q2 of KaonLT 2018-19
if Q2Val == 2.1:
    inpDict["Q2min"] = Q2Val - (2/9)*Q2Val
    inpDict["Q2max"] = Q2Val + (2/9)*Q2Val
    inpDict["Wmin"] = WVal - (2/9)*WVal
    inpDict["Wmax"] = WVal + (2/9)*WVal
elif Q2Val == 3.0:
    inpDict["Q2min"] = Q2Val - (3/7)*Q2Val
    inpDict["Q2max"] = Q2Val + (3/7)*Q2Val
    inpDict["Wmin"] = WVal - (3/7)*WVal
    inpDict["Wmax"] = WVal + (3/7)*WVal
elif Q2Val == 5.5:
    inpDict["Q2min"] = Q2Val - (3/8)*Q2Val
    inpDict["Q2max"] = Q2Val + (3/8)*Q2Val
    inpDict["Wmin"] = WVal - (3/8)*WVal
    inpDict["Wmax"] = WVal + (3/8)*WVal
else:
    inpDict["Q2min"] = Q2Val - (2/7)*Q2Val
    inpDict["Q2max"] = Q2Val + (2/7)*Q2Val
    inpDict["Wmin"] = WVal - (2/7)*WVal
    inpDict["Wmax"] = WVal + (2/7)*WVal
##############
##############
##############
# Default starting values no need to change
inpDict["Epsmin"] = 0.0
inpDict["Epsmax"] = 1.0

phisetlist = ["Center","Left","Right"]

for phiset in phisetlist:
    # Call diamond cut script and append paramters to dictionary
    inpDict.update(DiamondPlot(ParticleType, Q2Val, inpDict["Q2min"], inpDict["Q2max"], WVal, inpDict["Wmin"], inpDict["Wmax"], phiset, tmin, tmax, inpDict))
 
print("\n\nDiamond cut parameters: ")
for p in [1,2,3,4]:
    if inpDict["a%i" % p] == 0.0 or inpDict["b%i" % p] == 0.0:
        print("ERROR: Invalid diamond cut paramters")
        sys.exit(2)
    else:
        print("a{} = {}, b{} = {}".format(p,inpDict["a%i" % p],p,inpDict["b%i" % p]))    
    
if DEBUG:
    # Show plot pdf for each setting
    for phiset in phisetlist:
        show_pdf_with_evince(OUTPATH+"/{}_{}_diamond_{}.pdf".format(phiset, ParticleType, 'Q'+Q2+'W'+W))
for phiset in phisetlist:
    output_file_lst.append(OUTPATH+"/{}_{}_diamond_{}.pdf".format(phiset, ParticleType, 'Q'+Q2+'W'+W))
        
##############################
# Step 3 of the lt_analysis: # DONE
##############################
'''
Apply random subtraction to data and dummy.
'''

sys.path.append("cuts")
from rand_sub import rand_sub

# Call histogram function above to define dictonaries for right, left, center settings
# Put these all into an array so that if we are missing a setting it is easier to remove
# Plus it makes the code below less repetitive
histlist = []
for phiset in phisetlist:
    histlist.append(rand_sub(phiset,inpDict))
    
print("\n\n")

settingList = []
for i,hist in enumerate(histlist):
    if len(hist.keys()) <= 1: # If hist is empty (length of oen for phi setting check)
        print("No {} setting found...".format(hist["phi_setting"]))
        phisetlist.remove(hist["phi_setting"])
        histlist.remove(hist)
    else:
        settingList.append(hist["phi_setting"])

if DEBUG:
    # Show plot pdf for each setting
    for hist in histlist:        
        show_pdf_with_evince(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))        
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))

##############################
# Step 4 of the lt_analysis: # Done
##############################
'''
* Combine all settings and choose t/phi bins for low eps.

* These bins will also be used of high eps, so check high eps as well.
'''

'''
##############
# HARD CODED #
##############
# Redefine boundaries of Q2, W and eps
# May need to adjust these for binning
# Too many zero bins can result in empty bins when t/phi binning
# Good method is to use std dev around central set value
if Q2 == "2p1" and W == "2p95":
    inpDict["Q2min"] = Q2Val - (0.24)*Q2Val
    inpDict["Q2max"] = Q2Val + (0.24)*Q2Val
    inpDict["Wmin"] = WVal - (0.06)*WVal
    inpDict["Wmax"] = WVal + (0.06)*WVal
    inpDict["Epsmin"] = float(EPSVAL) - (0.02)*float(EPSVAL)
    inpDict["Epsmax"] = float(EPSVAL) + (0.02)*float(EPSVAL)
if Q2 == "5p5" and W == "3p02":
    inpDict["Q2min"] = Q2Val - (0.24)*Q2Val
    inpDict["Q2max"] = Q2Val + (0.24)*Q2Val
    inpDict["Wmin"] = WVal - (0.06)*WVal
    inpDict["Wmax"] = WVal + (0.06)*WVal
    inpDict["Epsmin"] = float(EPSVAL) - (0.09)*float(EPSVAL)
    inpDict["Epsmax"] = float(EPSVAL) + (0.09)*float(EPSVAL)
##############
##############
##############
'''

##############
# HARD CODED #
##############
# Reset Q2, W, eps TH1F range dynamically
# Adjust range_factor argument for plot limits (Default=0.005)
#inpDict["Q2min"], inpDict["Q2max"] = set_dynamic_axis_ranges("Q2", histlist)
#inpDict["Wmin"], inpDict["Wmax"] = set_dynamic_axis_ranges("W", histlist)
#inpDict["Epsmin"], inpDict["Epsmax"] = set_dynamic_axis_ranges("epsilon", histlist)
##############
##############
##############

sys.path.append("binning")
from find_bins import find_bins, check_bins

if EPSSET == "low":
    bin_vals = find_bins(histlist, inpDict)

try:
    with open("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p","")), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            t_bins = all_lines[1].split("\t")
            del t_bins[0]
            t_bins = np.array([float(element) for element in t_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))
except IOError:
    print("Error reading {}...".format("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))    

try:
    with open("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p","")), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            phi_bins = all_lines[1].split("\t")
            del phi_bins[0]
            phi_bins = np.array([float(element) for element in phi_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))
except IOError:
    print("Error reading {}...".format("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))    
    
for hist in histlist:
    hist["t_bins"] = t_bins
    hist["phi_bins"] = phi_bins

if EPSSET == "high":
    check_bins(histlist, inpDict)
    
##############################
# Step 5 of the lt_analysis: #
##############################
'''
* Compare SIMC to data/dummy setting by setting.

* SIMC weight may need to be recalculated.

* If so...

** Need to put physics_iterate.f (found in src/simc_ana) in simc_gfortran (adjust reaction for specific needs)
** Once physics_iterate.f is in simc_gfortran, change object in Makefile from physics_<PID>.o to physics_iterate.o
   Note: This also needs initial parameterization. 
         par_pl (rename to par.pl in SIMC )is currently in there, Q2=2.45 of FPI2 analysis
** Run SIMC with this new cross section (the hard coded q2_set value in physics_iterate.f is dynamically updated using the set_ProdInput script). 
   There is no need to rerun SIMC from this point on unless model functional form changes.
** Run this script (i.e. src/main.py) for first iteration. This will generate json and root file with data info
   and current simc info. It will generate new parameters (located src/<PID>/parameters/par.<pol>_<Q2>.dat) from
   calc_xsect.f script.
** With these new parameters, you will be prompted if you would like to rerun script. Keep rerunning until ratio of
   data/simc yields start to approach unity. Once they approach unity, the correct parameterization has been found.
   If after many iterations, things don't improve, the functional form (i.e. defined in physics_iterate.f and in
   src/models/param*.f and src/models/xmodel*.f) may need to be changed. Make sure to change these across all files.
** Once unseparated cross sections are achieved with good parameterization, cross section separation can be performed.

'''

sys.path.append("simc_ana")    
from compare_simc import compare_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    hist.update(compare_simc(hist, inpDict))    

if DEBUG:
    # Show plot pdf for each setting
    for hist in histlist:        
        show_pdf_with_evince(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(hist["phi_setting"],ParticleType)))
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(hist["phi_setting"],ParticleType)))

sys.path.append("normalize")
from get_eff_charge import get_eff_charge

# Upate hist dictionary with effective charge
for hist in histlist:
    hist.update(get_eff_charge(hist, inpDict))
    
sys.path.append("plotting")
from data_vs_simc import plot_data_vs_simc

# Variable defines string of cuts applied during analysis
cut_summary_lst = plot_data_vs_simc(t_bins, phi_bins, histlist, phisetlist, inpDict)

if DEBUG:
    show_pdf_with_evince(outputpdf)
output_file_lst.append(outputpdf)

##############################
# Step 6 of the lt_analysis: #
##############################
'''
* Combine settings again at for each 
  Q2,W,tbin,phibin,eps for data, dummy, and SIMC.

* Find the mean data values of W,Q2,theta,eps for each t bin of high and low epsilon
  for both data and SIMC.

* Dummy is subtracted from data bin by bin.
* The yield is calculated using the effective charge from data and 
  normfactor/nevents is applied to normalize simc.

* The data and SIMC yields are compared and the R value per bin is obtained.

'''

sys.path.append("binning")
from calculate_yield import find_yield_data, find_yield_simc

yieldDict = {}
yieldDict.update(find_yield_data(histlist, inpDict))
yieldDict.update(find_yield_simc(histlist, inpDict))

sys.path.append("binning")
from calculate_ratio import find_ratio

ratioDict = {}
ratioDict.update(find_ratio(histlist, inpDict, yieldDict))

sys.path.append("binning")
from ave_per_bin import ave_per_bin_data, ave_per_bin_simc

aveDict = {}
aveDict.update(ave_per_bin_data(histlist, inpDict))
aveDict.update(ave_per_bin_simc(histlist, inpDict))

sys.path.append("plotting")
from binned import plot_binned

plot_binned(t_bins, phi_bins, histlist, phisetlist, inpDict, yieldDict, ratioDict, aveDict)
#notify_email(email_address="trotta@cua.edu")

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))    

# Save histograms to root file
# Check that root file doesnt already exist    
if not os.path.exists(foutroot):
#if os.path.exists(foutroot):
    for hist in histlist:
        print("\nSaving {} histograms to {}".format(hist["phi_setting"],foutroot))
        # Loop through all keggys,values of dictionary
        for i, (key, val) in enumerate(hist.items()):
            # Progress bar
            Misc.progressBar(i, len(hist.items())-1,bar_length=25)
            if "G_data_eff" in key:
                hist_to_root(val, foutroot, "{}/data".format(hist["phi_setting"]))
            if is_hist(val):
                if "ratio" in key:
                    continue
                if "DATA" in key:
                    if "yield" in key:
                        continue
                    elif "bin" in key:
                        continue
                    elif "totevts" in key:
                        continue
                    else:
                        hist_to_root(val, foutroot, "{}/data".format(hist["phi_setting"]))
                if "SIMC" in key:
                    if "yield" in key:
                        continue
                    elif "bin" in key:
                        continue
                    elif "totevts" in key:
                        continue
                    else:
                        hist_to_root(val, foutroot, "{}/simc".format(hist["phi_setting"]))
                if "DUMMY" in key:
                    hist_to_root(val, foutroot, "{}/dummy".format(hist["phi_setting"]))

    # Open the ROOT file
    root_file = TFile.Open(foutroot, "UPDATE")

    # Check if the file was opened successfully
    if root_file.IsOpen():
        # Close the file
        root_file.Close()
        print("\nThe root file {} has been successfully closed.".format(foutroot))
    else:
        print("\nError: Unable to close the root file {}.".format(foutroot))
output_file_lst.append(foutroot)

# Check that root file doesnt already exist
if not os.path.exists(foutjson):
#if os.path.exists(foutjson):
    # Create combined dictionary of all non-histogram information        
    combineDict = {}
    combineDict.update({"inpDict" : inpDict})
    tmp_lst = []
    for hist in histlist:
        print("\nSaving {} information to {}".format(hist["phi_setting"],foutjson))
        tmp_dict = {}
        for i, (key, val) in enumerate(hist.items()):
            # Progress bar
            Misc.progressBar(i, len(hist.items())-1,bar_length=25)
            if not is_root_obj(val):
                tmp_dict[key] = val
        tmp_lst.append(tmp_dict)
    combineDict.update({ "histlist" : tmp_lst})

    # Save combined dictionary to json file
    # Open the file in write mode and use json.dump() to save the dictionary to JSON
    with open(foutjson, 'w') as f_json:
        json.dump(combineDict, f_json, default=custom_encoder)
output_file_lst.append(foutjson)

from physics_lists import create_lists
create_lists(aveDict, yieldDict, histlist, inpDict, phisetlist, output_file_lst)

# Copy initial parameterization to specific particle type directory
shutil.copy('{}/src/models/par_{}_Q{}W{}'.format(LTANAPATH, pol_str, Q2.replace("p",""), W.replace("p","")), '{}/src/{}/parameters/par.{}_Q{}W{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")))

# ***Parameter file from last iteration!***
# ***These old parameters are needed for this iteration. See README for more info on procedure!***
old_param_file = '{}/src/{}/parameters/par.{}_Q{}W{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""))
cut_summary_lst += "\n\nUnsep Parameterization for {}...\n".format(formatted_date)
with open(old_param_file, 'r') as file:
    for line in file:
        cut_summary_lst += line
        
print("\n\n")
print("="*25)
print("{} Epsilon Summary...".format(EPSSET.capitalize()),cut_summary_lst)
print("="*25)
inpDict["cut_summary_lst"] = cut_summary_lst

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data

* Calculate the unseparated cross section
'''

if EPSSET == "high":
    
    # Save fortran scripts that contain iteration functional form of parameterization
    py_param = 'models/param_{}_{}.py'.format(ParticleType, pol_str)
    output_file_lst.append(py_param) 
    fort_xmodel = 'models/xmodel_{}_{}.f'.format(ParticleType, pol_str)
    output_file_lst.append(fort_xmodel)

    # Active scripts to make file selection dynamic
    # Needs to be done this way because of fortran compiler limitations
    py_param_active = 'models/param_active.py'
    fort_xmodel_active = 'models/xmodel_active.f'
    # Copying content of used models to actively used files
    print("Copying {} to {}".format(LTANAPATH+"/src/"+fort_xmodel, LTANAPATH+"/src/"+fort_xmodel_active))
    shutil.copy(LTANAPATH+"/src/"+fort_xmodel, LTANAPATH+"/src/"+fort_xmodel_active)
    print("Copying {} to {}".format(LTANAPATH+"/src/"+py_param, LTANAPATH+"/src/"+py_param_active))    
    shutil.copy(LTANAPATH+"/src/"+py_param, LTANAPATH+"/src/"+py_param_active)

    # Save python script that contain separated xsect models for lt script
    py_lt = 'models/lt_{}_{}.py'.format(ParticleType, pol_str)
    output_file_lst.append(py_lt)
    # Active scripts to make file selection dynamic
    # Needs to be done this way because of fortran compiler limitations
    py_lt_active = 'models/lt_active.py'
    # Copying content of used models to actively used files
    print("Copying {} to {}".format(LTANAPATH+"/src/"+py_lt, LTANAPATH+"/src/"+py_lt_active))
    shutil.copy(LTANAPATH+"/src/"+py_lt, LTANAPATH+"/src/"+py_lt_active)
    
    # run_xsect bash script calls average_kinematics.f to find error weighted average of data.
    # It then runs calc_xsect.f to find unseparated cross section as well as new set of parameters
    # if still iterating weights
    try:
        subprocess.call(['bash', '{}/run_xsect.sh'.format(LTANAPATH), Q2, W, ParticleType, POL, str(NumtBins), str(NumPhiBins)])
    except Exception as e:
        print("{}".format(e))
        sys.exit(2)

    if DEBUG:
        show_pdf_with_evince(OUTPATH+"/{}_lt_fit_Q{}W{}.pdf".format(ParticleType, Q2, W))
        show_pdf_with_evince(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))
    output_file_lst.append(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))
    output_file_lst.append(OUTPATH+"/{}_lt_fit_Q{}W{}.pdf".format(ParticleType, Q2, W))
    output_file_lst.append('models/lt_2D_fit.py')
    
    # Save new parameters and unsep values from current iteration
    # ***Old parameter file defined in step 7, the new parameter values are saved here!***
    # ***The old parameters, used for this iteration, are saved in the summary!***
    new_param_file = '{}/parameters/par.{}_Q{}W{}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""))
    output_file_lst.append(new_param_file) 
    sep_file = '{}/xsects/x_sep.{}_Q{}W{}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""))
    output_file_lst.append(sep_file)    
    unsep_lo_file = '{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100)
    output_file_lst.append(unsep_lo_file)
    unsep_hi_file = '{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100)
    output_file_lst.append(unsep_hi_file)        
    avek_file = '{}/averages/avek.Q{}W{}.dat'.format(ParticleType, Q2.replace("p",""), W.replace("p",""))
    output_file_lst.append(avek_file)
    aver_lo_file = '{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100)
    output_file_lst.append(aver_lo_file)
    aver_hi_file = '{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100)
    output_file_lst.append(aver_hi_file)    

##############################
# Step 8 of the lt_analysis: #
##############################
'''
* Save all information for this iteration
'''

# Create a new directory for each iteration in cache
new_dir = "{}/{}/{}/Q{}W{}/{}".format(CACHEPATH, USER, ParticleType.lower(), Q2, W, formatted_date)
create_dir(new_dir)

if EPSSET == "high":
    
    print("\n\n")

    # Grab simc root file
    for hist in histlist:
        for eps in ["highe","lowe"]:
            f_simc_root = OUTPATH+"/Prod_Coin_{}.root".format(kinematics[0]+hist["phi_setting"].lower()+"_"+eps)
            f_simc_hist = OUTPATH+"/Prod_Coin_{}.hist".format(kinematics[0]+hist["phi_setting"].lower()+"_"+eps)
            if os.path.exists(f_simc_root):
                output_file_lst.append(f_simc_root)
            if os.path.exists(f_simc_hist):
                output_file_lst.append(f_simc_hist)                

f_path = "{}/{}_Q{}W{}_iter.dat".format(LTANAPATH,ParticleType,Q2,W)
# Check if the file exists
if os.path.exists(f_path):
    # If it exists, update it with the string
    with open(f_path, 'a') as file:
        file.write('\n'+formatted_date)
else:
    # If not, create it and fill it with the string
    with open(f_path, 'x') as file:
        file.write(formatted_date)

# Get the total number of lines in the file
with open(f_path, 'r') as file:
    total_lines = len(file.readlines())

f_path_new = f_path.replace(LTANAPATH,new_dir).replace("iter","iter_{}".format(total_lines-1))
print("Copying {} to {}".format(f_path,f_path_new))
shutil.copy(f_path,f_path_new)

for f in output_file_lst:
    if OUTPATH in f:
        if ".pdf" in f:
            create_dir(new_dir+"/plots")
            f_new = f.replace(OUTPATH,new_dir+"/plots")
            print("Copying {} to {}".format(f,f_new))
            shutil.copy(f, f_new)
        if ".json" in f:
            create_dir(new_dir+"/json")
            f_new = f.replace(OUTPATH,new_dir+"/json")
            print("Copying {} to {}".format(f,f_new))
            shutil.copy(f, f_new)                
        if ".root" in f:
            create_dir(new_dir+"/root")
            f_new = f.replace(OUTPATH,new_dir+"/root")
            print("Copying {} to {}".format(f,f_new))
            shutil.copy(f, f_new)
        if ".hist" in f:
            create_dir(new_dir+"/root")
            f_new = f.replace(OUTPATH,new_dir+"/root")
            print("Copying {} to {}".format(f,f_new))
            shutil.copy(f, f_new)
    elif "{}/".format(ParticleType) in f:
        f_arr = f.split("/")
        f_tmp = f_arr.pop()
        for f_dir in f_arr:
            if "{}".format(ParticleType) not in f_dir:
                create_dir(new_dir+"/"+f_dir)
                f_new = new_dir+"/"+f_dir+"/"+f_tmp    
                print("Copying {} to {}".format(LTANAPATH+"/src/"+f,f_new))
                shutil.copy(LTANAPATH+"/src/"+f, f_new)
    else:
        f_new = new_dir
        print("Copying {} to {}".format(LTANAPATH+"/src/"+f,f_new))
        shutil.copy(LTANAPATH+"/src/"+f, f_new)
               
# Need summary for both high and low eps.
# All others should be saved once both are complete
with open(new_dir+'/{}_{}_summary_{}.txt'.format(ParticleType,OutFilename,formatted_date), 'w') as file:
    file.write(inpDict["cut_summary_lst"])

'''
for hist in histlist:
    key_str = ', '.join(hist.keys())
    print("{} keys: {}".format(hist["phi_setting"],key_str))
'''
