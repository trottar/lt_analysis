#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-31 00:20:12 trottar"
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
from utility import show_pdf_with_evince, create_dir, is_root_obj, is_hist, hist_to_root, custom_encoder

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=34:
    print("!!!!! ERROR !!!!!\n Expected 34 arguments\n Usage is with - KIN W Q2 EPSVAL OutDATAFilename OutDUMMYFilename OutFullAnalysisFilename tmin tmax NumtBins NumPhiBins runNumRight runNumLeft runNumCenter data_charge_right data_charge_left data_charge_center dummy_charge_right dummy_charge_left dummy_charge_center InData_efficiency_right InData_efficiency_left InData_efficiency_center efficiency_table ParticleType EPSSET pThetaValRight pThetaValLeft pThetaValCenter EbeamValRight EbeamValLeft EbeamValCenter POL formatted_date\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    

DEBUG = False # Flag for plots

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
EPSSET = sys.argv[26]
pThetaValRight = list(sys.argv[27].split(" "))
pThetaValLeft = list(sys.argv[28].split(" "))
pThetaValCenter = list(sys.argv[29].split(" "))
EbeamValRight = list(sys.argv[30].split(" "))
EbeamValLeft = list(sys.argv[31].split(" "))
EbeamValCenter = list(sys.argv[32].split(" "))
POL = sys.argv[33]
formatted_date = sys.argv[34]

if int(POL) == 1:
    pol_str = "pl"
elif int(POL) == -1:
    pol_str = "mn"
else:
    print("ERROR: Invalid polarity...must be +1 or -1")
    sys.exit(2)

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

#Importing diamond cut script
sys.path.append("cuts")
from diamond import DiamondPlot

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

##############
# HARD CODED #
##############
# May need to adjust these for diamond plots to work
# 2/7 seems to work for most Q2 of KaonLT 2018-19
if Q2Val == 3.0:
    inpDict["Q2min"] = Q2Val - (3/7)*Q2Val
    inpDict["Q2max"] = Q2Val + (3/7)*Q2Val
    inpDict["Wmin"] = WVal - (3/7)*WVal
    inpDict["Wmax"] = WVal + (3/7)*WVal
else:
    inpDict["Q2min"] = Q2Val - (2/7)*Q2Val
    inpDict["Q2max"] = Q2Val + (2/7)*Q2Val
    inpDict["Wmin"] = WVal - (2/7)*WVal
    inpDict["Wmax"] = WVal + (2/7)*WVal
##############
##############
##############

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
        show_pdf_with_evince(OUTPATH+"/{}_{}_{}_Diamond_Cut.pdf".format(ParticleType, 'Q'+Q2+'W'+W, phiset))
for phiset in phisetlist:
    output_file_lst.append(OUTPATH+"/{}_{}_{}_Diamond_Cut.pdf".format(ParticleType, 'Q'+Q2+'W'+W, phiset))
        
##############################
# Step 3 of the lt_analysis: # DONE
##############################
'''
Apply random subtraction to data and dummy.
'''

##############
# HARD CODED #
##############
# Redefine boundaries of Q2, W and eps
# May need to adjust these for binning
# Too many zero bins can result in empty bins when t/phi binning
# Good method is to use std dev around central set value
inpDict["Q2min"] = Q2Val - (0.24)*Q2Val
inpDict["Q2max"] = Q2Val + (0.24)*Q2Val
inpDict["Wmin"] = WVal - (0.06)*WVal
inpDict["Wmax"] = WVal + (0.06)*WVal
inpDict["Epsmin"] = float(EPSVAL) - (0.02)*float(EPSVAL)
inpDict["Epsmax"] = float(EPSVAL) + (0.02)*float(EPSVAL)
##############
##############
##############

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
        show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))

##############################
# Step 4 of the lt_analysis: # Done
##############################
'''
* Combine all settings and choose t/phi bins for low eps.

* These bins will also be used of high eps, so check high eps as well.
'''

sys.path.append("binning")
from find_bins import find_bins, check_bins

if EPSSET == "low":
    bin_vals = find_bins(histlist, inpDict)

try:
    with open("{}/src/{}/t_bin_interval".format(LTANAPATH, ParticleType), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            t_bins = all_lines[1].split("\t")
            del t_bins[0]
            t_bins = np.array([float(element) for element in t_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/t_bin_interval".format(LTANAPATH, ParticleType)))
except IOError:
    print("Error reading {}...".format("{}/src/{}/t_bin_interval".format(LTANAPATH, ParticleType)))    

try:
    with open("{}/src/{}/phi_bin_interval".format(LTANAPATH, ParticleType), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            phi_bins = all_lines[1].split("\t")
            del phi_bins[0]
            phi_bins = np.array([float(element) for element in phi_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/phi_bin_interval".format(LTANAPATH, ParticleType)))
except IOError:
    print("Error reading {}...".format("{}/src/{}/phi_bin_interval".format(LTANAPATH, ParticleType)))    
    
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

sys.path.append("normalize")
from get_eff_charge import get_eff_charge

# Upate hist dictionary with effective charge
for hist in histlist:
    hist.update(get_eff_charge(hist, inpDict))

sys.path.append("simc_ana")    
from compare_simc import compare_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    hist.update(compare_simc(hist, inpDict))
    
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
ratioDict.update(find_ratio(histlist, inpDict, phisetlist, yieldDict))

sys.path.append("binning")
from ave_per_bin import ave_per_bin_data, ave_per_bin_simc

aveDict = {}
aveDict.update(ave_per_bin_data(histlist, inpDict))
aveDict.update(ave_per_bin_simc(histlist, inpDict))

sys.path.append("plotting")
from binned import plot_binned

plot_binned(t_bins, phi_bins, histlist, phisetlist, inpDict, yieldDict, ratioDict, aveDict)

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))    

# Save histograms to root file
# Check that root file doesnt already exist    
if not os.path.exists(foutroot):
    for hist in histlist:
        print("\nSaving {} histograms to {}".format(hist["phi_setting"],foutroot))
        # Loop through all keggys,values of dictionary
        for i, (key, val) in enumerate(hist.items()):
            # Progress bar
            Misc.progressBar(i, len(hist.items())-1,bar_length=25)
            if is_hist(val):            
                if "ratio" in key:
                    hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))
                if "DATA" in key:
                    if "yield" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))                        
                    elif "bin" in key:
                        hist_to_root(val, foutroot, "{}/bins".format(hist["phi_setting"]))
                    elif "totevts" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))
                    else:
                        hist_to_root(val, foutroot, "{}/data".format(hist["phi_setting"]))
                if "SIMC" in key:
                    if "yield" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))                        
                    elif "bin" in key:
                        hist_to_root(val, foutroot, "{}/bins".format(hist["phi_setting"]))
                    elif "totevts" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))
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
create_lists(aveDict, ratioDict, histlist, inpDict, phisetlist, output_file_lst)

# Copy initial parameterization to specific particle type directory
shutil.copy(LTANAPATH+'/src/simc_ana/par_pl', '{}/src/{}/parameters/par.{}_{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p","")))

# ***Parameter file from last iteration!***
# ***These old parameters are needed for this iteration. See README for more info on procedure!***
old_param_file = '{}/src/{}/parameters/par.{}_{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""))
try:
    cut_summary_lst += "\n\nUnsep Parameterization for {}...".format(formatted_date)
    with open(old_param_file, 'r') as file:
        for line in file:
            cut_summary_lst += line
except FileNotFoundError:
    print('''
    \n\n
    File not found!
    Assuming first iteration!
    ''')
    
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

    # run_xsect bash script calls average_kinematics.f to find error weighted average of data.
    # It then runs calc_xsect.f to find unseparated cross section as well as new set of parameters
    # if still iterating weights
    try:
        subprocess.call(['bash', '{}/run_xsect.sh'.format(LTANAPATH), Q2, W, ParticleType, POL, str(NumtBins), str(NumPhiBins)])
    except Exception as e:
        print("{}".format(e))
        sys.exit(2)

    if DEBUG:
        show_pdf_with_evince(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))    
    output_file_lst.append(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))

    # Save new parameters and unsep values from current iteration
    # ***Old parameter file defined in step 7, the new parameter values are saved here!***
    # ***The old parameters, used for this iteration, are saved in the summary!***
    new_param_file = '{}/parameters/par.{}_{}.dat'.format(ParticleType, pol_str, Q2.replace("p",""))
    output_file_lst.append(new_param_file) 
    unsep_file = '{}/xsects/x_unsep.{}_{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(unsep_file)
    sep_file = '{}/xsects/x_sep.{}_{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(sep_file)

##############################
# Step 8 of the lt_analysis: #
##############################
'''
* Save all information for this iteration
'''

# Create a new directory for each iteration in cache
new_dir = CACHEPATH+"/"+USER+"/"+ParticleType.lower()+"/"+formatted_date
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
        
    # Grab low eps versions as well
    for f in output_file_lst:
        if OutFilename in f:
            f_lowe = f.replace("highe","lowe")
            if os.path.exists(f_lowe):
                output_file_lst.append(f_lowe)
                
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

    f_path_new = f_path.replace(LTANAPATH,new_dir).replace("iter","iter_{}".format(total_lines))
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
