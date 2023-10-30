#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-30 15:54:51 trottar"
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
from utility import show_pdf_with_evince, create_dir, is_root_obj, is_hist, hist_to_root, last_iter, get_histogram, hist_in_dir

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=11:
    print("!!!!! ERROR !!!!!\n Expected 11 arguments\n Usage is with - KIN W Q2 EPSVAL ParticleType EPSSET POL OutFilename formatted_date NumtBins NumPhiBins\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    

DEBUG = False # Flag for plots

# Input params
kinematics = sys.argv[1].split("_")
W = sys.argv[2]
Q2 = sys.argv[3]
EPSVAL = sys.argv[4]
ParticleType = sys.argv[5]
EPSSET = sys.argv[6]
POL = sys.argv[7]
OutFilename = sys.argv[8]
formatted_date = sys.argv[9]
NumtBins = sys.argv[10]
NumPhiBins = sys.argv[11]


if int(POL) == 1:
    pol_str = "pl"
elif int(POL) == -1:
    pol_str = "mn"
else:
    print("ERROR: Invalid polarity...must be +1 or -1")
    sys.exit(2)

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

################################################################################################################################################    
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

################################
# Step 1-4 of the lt_analysis: #
################################
'''
* Analysis already done for data, see main.py

* Here is just appending files to list so it can be copied to iteration in cache

* Need to read in previous iteration histograms, information, etc.

'''
    
# Find last iteration, based of closest date
f_iter = "{}/{}_Q{}W{}_iter.dat".format(LTANAPATH,ParticleType,Q2,W)
closest_date = last_iter(f_iter, formatted_date)
print("\n\nThe last iteration was ",closest_date)

# Save this as the directory to grab further information
prev_iter_dir = "{}/{}/{}/{}".format(CACHEPATH,USER,ParticleType.lower(),closest_date)

prev_iter_root = foutroot.replace(OUTPATH,prev_iter_dir+"/root")
prev_iter_json = foutjson.replace(OUTPATH,prev_iter_dir+"/json")

# Redefine dictionaries from old iteration information, see main.py
with open(prev_iter_json, 'r') as f:
    prev_iter_combineDict = json.load(f)

inpDict = prev_iter_combineDict["inpDict"]
histlist = prev_iter_combineDict["histlist"]

root_file = TFile.Open(prev_iter_root, "READ")
# Grab weight from previous iteration
for hist in histlist:
    hist.update(hist_in_dir(root_file, "{}/data".format(hist["phi_setting"])))
    hist.update(hist_in_dir(root_file, "{}/simc".format(hist["phi_setting"])))
    hist.update(hist_in_dir(root_file, "{}/dummy".format(hist["phi_setting"])))

# t/phi bins are the same for all settings
# so arbitrarily grabbing from first setting of list
t_bins = np.array(histlist[0]["t_bins"])
phi_bins = np.array(histlist[0]["phi_bins"])

print("\n\nt_bins = ", t_bins)
print("phi_bins = ", phi_bins)

phisetlist = []
for hist in histlist:
    phisetlist.append(hist["phi_setting"])

for phiset in phisetlist:
    output_file_lst.append(OUTPATH+"/{}_{}_{}_Diamond_Cut.pdf".format(ParticleType, 'Q'+Q2+'W'+W, phiset))
        
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))

'''
EXAMPLE: How to get histograms from previous iteration

# Open the ROOT file, must pass open root file so object exists here and in function
root_file = TFile.Open(prev_iter_root, "READ")
# Grab weight from previous iteration
iter_weight = get_histogram(root_file, "{}/simc".format(hist["phi_setting"]), "H_Weight_SIMC")
'''
    
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
** Run SIMC with this new cross section. There is no need to rerun SIMC from this point on unless model functional
   form changes.
** Run this script (i.e. src/main.py) for first iteration. This will generate json and root file with data info
   and current simc info. It will generate new parameters (located src/<PID>/parameters/par.<pol>_<Q2>.dat) from
   calc_xsect.f script.
** With these new parameters, you will be prompted if you would like to rerun script. Keep rerunning until ratio of
   data/simc yields start to approach unity. Once they approach unity, the correct parameterization has been found.
   If after many iterations, things don't improve, the functional form (i.e. defined in physics_iterate.f and in
   src/models/param*.f and src/models/xmodel*.f) may need to be changed. Make sure to change these across all files.
** Once unseparated cross sections are achieved with good parameterization, cross section separation can be performed.

'''

print("\n\n")

sys.path.append("normalize")
from get_eff_charge import get_eff_charge
# Upate hist dictionary with effective charge
for hist in histlist:
    hist.update(get_eff_charge(hist, inpDict))

# Create a new directory for each iteration in cache
# ***Moved up in procedure vs main.py since required for weight iteration***
new_dir = CACHEPATH+"/"+USER+"/"+ParticleType.lower()+"/"+formatted_date
create_dir(new_dir)
# ***Also must create new root directory in iter directory***
create_dir(new_dir+"/root")

# ***Parameter file from last iteration!***
# ***These old parameters are needed for this iteration. See README for more info on procedure!***
old_param_file = '{}/src/{}/parameters/par.{}_{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""))

# ***Moved from main.py location below because neede for weight iteration***
# Save fortran scripts that contain iteration functional form of parameterization
py_param = 'models/param_{}_{}.py'.format(ParticleType, pol_str)
output_file_lst.append(py_param) 

# Active scripts to make file selection dynamic
# Needs to be done this way because of fortran compiler limitations
py_param_active = 'models/param_active.py'
# Copying content of used models to actively used files
print("Copying {} to {}".format(LTANAPATH+"/src/"+py_param, LTANAPATH+"/src/"+py_param_active))    
shutil.copy(LTANAPATH+"/src/"+py_param, LTANAPATH+"/src/"+py_param_active)

sys.path.append("simc_ana")
from iter_weight import iter_weight
from compare_simc_iter import compare_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    # SIMC file with weight from last iteration
    old_simc_root = '{}/root/Prod_Coin_{}.root'.format(prev_iter_dir, kinematics[0]+hist["phi_setting"].lower()+"_"+kinematics[1])
    old_simc_hist = '{}/root/Prod_Coin_{}.hist'.format(prev_iter_dir, kinematics[0]+hist["phi_setting"].lower()+"_"+kinematics[1])
    new_simc_root = old_simc_root.replace(closest_date, formatted_date)
    new_simc_hist = old_simc_hist.replace(closest_date, formatted_date)
    # Make sure old simc root file exists
    if os.path.exists(old_simc_root):
        # Copy to new iteration so and then edit the weight
        print("Copying {} to {}".format(old_simc_root, new_simc_root))
        shutil.copy(old_simc_root,new_simc_root)
        shutil.copy(old_simc_hist,new_simc_hist)        
        # Make sure new simc root file exists
        if os.path.exists(new_simc_root):
            # Function to calculation new weight and apply it to simc root file 
            iter_weight(old_param_file, new_simc_root, inpDict, hist["phi_setting"])
            # Overwrite root file with updated weight
            os.rename(new_simc_root.replace(".root","_new.root"),new_simc_root)
            hist.update(compare_simc(new_simc_root, hist, inpDict))
        else:
            print("ERROR: {} not properly copied to {}".format(old_simc_root, new_simc_root))
            sys.exit(2)
            
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
        print("\nUpdating simc {} histograms in {}".format(hist["phi_setting"],foutroot))
        # Loop through all keys,values of dictionary
        for i, (key, val) in enumerate(hist.items()):
            # Progress bar
            Misc.progressBar(i, len(hist.items())-1,bar_length=25)
            if is_hist(val):
                if "ratio" in key:
                    hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))
                if "SIMC" in key:
                    if "yield" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))                        
                    elif "bin" in key:
                        hist_to_root(val, foutroot, "{}/bins".format(hist["phi_setting"]))
                    elif "totevts" in key:
                        hist_to_root(val, foutroot, "{}/yield".format(hist["phi_setting"]))
                    else:
                        hist_to_root(val, foutroot, "{}/simc".format(hist["phi_setting"]))

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

# Create combined dictionary of all non-histogram information        
combineDict = {}
combineDict.update({"inpDict" : inpDict})
tmp_lst = []
for hist in histlist:
    print("\nSaving {} information to {}".format(hist["phi_setting"],foutjson))
    for i, (key, val) in enumerate(hist.items()):
        # Progress bar
        Misc.progressBar(i, len(hist.items())-1,bar_length=25)
        if not is_root_obj(val):
            tmp_lst.append({key : val})
combineDict.update({ "histlist" : tmp_lst})

# Save combined dictionary to json file
# Check that root file doesnt already exist    
if not os.path.exists(foutjson):
    # Open the file in write mode and use json.dump() to save the dictionary to JSON
    with open(foutjson, 'w') as f_json:
        json.dump(combineDict, f_json, default=custom_encoder)
output_file_lst.append(foutjson)

from physics_lists import create_lists
create_lists(aveDict, ratioDict, histlist, inpDict, phisetlist, output_file_lst)

# Redefinition from above, but should be the same! This is just to stay consistent with main.py
# ***Parameter file from last iteration!***
# ***These old parameters are needed for this iteration. See README for more info on procedure!***
old_param_file = '{}/src/{}/parameters/par.{}_{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""))
try:
    cut_summary_lst += "\nUnsep Parameterization for {}...".format(formatted_date)
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

    # ***param python script is moved from main.py location to before weight iteration***
    
    # Save fortran scripts that contain iteration functional form of parameterization
    fort_xmodel = 'models/xmodel_{}_{}.f'.format(ParticleType, pol_str)
    output_file_lst.append(fort_xmodel)

    # Active scripts to make file selection dynamic
    # Needs to be done this way because of fortran compiler limitations
    fort_xmodel_active = 'models/xmodel_active.f'
    # Copying content of used models to actively used files
    print("Copying {} to {}".format(LTANAPATH+"/src/"+fort_xmodel, LTANAPATH+"/src/"+fort_xmodel_active))
    shutil.copy(LTANAPATH+"/src/"+fort_xmodel, LTANAPATH+"/src/"+fort_xmodel_active)

    # run_xsect bash script calls average_kinematics.f to find error weighted average of data.
    # It then runs calc_xsect.f to find unseparated cross section as well as new set of parameters
    # if still iterating weights
    try:
        subprocess.call(['bash', '{}/run_xsect.sh'.format(LTANAPATH), Q2, W, ParticleType, POL, str(NumtBins), str(NumPhiBins)])
    except Exception as e:
        print("1 ERROR: {}".format(e))
        sys.exit(2)

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

# ***Moved creation of iteration directory up from where it is in main.py. Now is near the new weight calculation***

if EPSSET == "high":
    
    print("\n\n")

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
