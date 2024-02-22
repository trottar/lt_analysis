#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-22 16:43:26 trottar"
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
from utility import show_pdf_with_evince, create_dir, is_root_obj, is_hist, hist_to_root, last_iter, get_histogram, hist_in_dir, custom_encoder, notify_email

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=13:
    print("!!!!! ERROR !!!!!\n Expected 13 arguments\n Usage is with - KIN W Q2 LOEPS HIEPS ParticleType EPSSET POL OutFilename formatted_date NumtBins NumPhiBins inp_debug\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################    

# Input params
kinematics = sys.argv[1].split("_")
W = sys.argv[2]
Q2 = sys.argv[3]
LOEPS = sys.argv[4]
HIEPS = sys.argv[5]
ParticleType = sys.argv[6]
EPSSET = sys.argv[7]
POL = sys.argv[8]
OutFilename = sys.argv[9]
formatted_date = sys.argv[10]
NumtBins = sys.argv[11]
NumPhiBins = sys.argv[12]
inp_debug =  sys.argv[13]

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
prev_iter_dir = "{}/{}/{}/Q{}W{}/{}".format(CACHEPATH,USER,ParticleType.lower(),Q2,W,closest_date)

# Copy all files from previous iteration to OUTPATH to assure consistency
# List all files in the source directory
files = os.listdir(prev_iter_dir+'/averages/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/averages/'+f, '{}/src/{}/averages/'.format(LTANAPATH, ParticleType)))
    shutil.copy(prev_iter_dir+'/averages/'+f, '{}/src/{}/averages/'.format(LTANAPATH, ParticleType))
files = os.listdir(prev_iter_dir+'/kindata/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/kindata/'+f, '{}/src/{}/kindata/'.format(LTANAPATH, ParticleType)))
    shutil.copy(prev_iter_dir+'/kindata/'+f, '{}/src/{}/kindata/'.format(LTANAPATH, ParticleType))
files = os.listdir(prev_iter_dir+'/parameters/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/parameters/'+f, '{}/src/{}/parameters/'.format(LTANAPATH, ParticleType)))
    shutil.copy(prev_iter_dir+'/parameters/'+f, '{}/src/{}/parameters/'.format(LTANAPATH, ParticleType))
files = os.listdir(prev_iter_dir+'/xsects/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/xsects/'+f, '{}/src/{}/xsects/'.format(LTANAPATH, ParticleType)))
    shutil.copy(prev_iter_dir+'/xsects/'+f, '{}/src/{}/xsects/'.format(LTANAPATH, ParticleType))
files = os.listdir(prev_iter_dir+'/yields/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/yields/'+f, '{}/src/{}/yields/'.format(LTANAPATH, ParticleType)))
    shutil.copy(prev_iter_dir+'/yields/'+f, '{}/src/{}/yields/'.format(LTANAPATH, ParticleType))
print("Copying {} to {}".format(prev_iter_dir+'/lt_2D_fit.py', '{}/src/models/lt_2D_fit.py'.format(LTANAPATH)))
shutil.copy(prev_iter_dir+'/lt_2D_fit.py', '{}/src/models/lt_2D_fit.py'.format(LTANAPATH))
print("Copying {} to {}".format(prev_iter_dir+'/lt_kaon_pl.py', '{}/src/models/lt_kaon_pl.py'.format(LTANAPATH)))
shutil.copy(prev_iter_dir+'/lt_kaon_pl.py', '{}/src/models/lt_kaon_pl.py'.format(LTANAPATH))
print("Copying {} to {}".format(prev_iter_dir+'/param_kaon_pl.py', '{}/src/models/param_kaon_pl.py'.format(LTANAPATH)))
shutil.copy(prev_iter_dir+'/param_kaon_pl.py', '{}/src/models/param_kaon_pl.py'.format(LTANAPATH))
print("Copying {} to {}".format(prev_iter_dir+'/xmodel_kaon_pl.f', '{}/src/models/xmodel_kaon_pl.f'.format(LTANAPATH)))
shutil.copy(prev_iter_dir+'/xmodel_kaon_pl.f', '{}/src/models/xmodel_kaon_pl.f'.format(LTANAPATH))
files = os.listdir(prev_iter_dir+'/root/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/root/'+f, OUTPATH))
    shutil.copy(prev_iter_dir+'/root/'+f, OUTPATH)
files = os.listdir(prev_iter_dir+'/json/')
for f in files:
    print("Copying {} to {}".format(prev_iter_dir+'/json/'+f, OUTPATH))
    shutil.copy(prev_iter_dir+'/json/'+f, OUTPATH)

prev_iter_root = foutroot.replace(OUTPATH,prev_iter_dir+"/root")
prev_iter_json = foutjson.replace(OUTPATH,prev_iter_dir+"/json")

# Redefine dictionaries from old iteration information, see main.py
with open(prev_iter_json, 'r') as f:
    prev_iter_combineDict = json.load(f)

inpDict = prev_iter_combineDict["inpDict"]
histlist = prev_iter_combineDict["histlist"]

# Add closest and formatted dates to inpDict (used in plot comparison)
inpDict["closest_date"] = closest_date
inpDict["formatted_date"] = formatted_date

if EPSSET == "low":
    # Save python script that contain separated xsect models for xfit script
    py_xfit = 'models/xfit_{}_{}.py'.format(ParticleType, pol_str)
    output_file_lst.append(py_xfit)

    # Active scripts to make file selection dynamic
    # Needs to be done this way because of fortran compiler limitations
    py_xfit_active = 'models/xfit_active.py'
    # Copying content of used models to actively used files
    print("Copying {} to {}".format(LTANAPATH+"/src/"+py_xfit, LTANAPATH+"/src/"+py_xfit_active))
    shutil.copy(LTANAPATH+"/src/"+py_xfit, LTANAPATH+"/src/"+py_xfit_active)
    
    # Run weight iteration script for optimizing parameterization
    sys.path.append("models")
    from xfit_in_t import x_fit_in_t
    x_fit_in_t(ParticleType, pol_str, closest_date, Q2, W, inpDict)
    if DEBUG:
        show_pdf_with_evince(OUTPATH+"/{}_xfit_in_t_Q{}W{}.pdf".format(ParticleType, Q2, W))
    #show_pdf_with_evince(OUTPATH+"/{}_xfit_in_t_Q{}W{}.pdf".format(ParticleType, Q2, W))
    output_file_lst.append(OUTPATH+"/{}_xfit_in_t_Q{}W{}.pdf".format(ParticleType, Q2, W))

# ***Parameter file for new iteration!***
# ***These parameters are newly generated for this iteration above. See README for more info on procedure!***
new_param_file = '{}/src/{}/parameters/par.{}_Q{}W{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""))
output_file_lst.append(new_param_file)

# Grab combined root files for data and dummy.
# Then save to dictionary
for hist in histlist:
    rootFileData = OUTPATH + "/" + "{}".format(ParticleType) + "_" + inpDict["InDATAFilename"] + "_%s.root" % (hist["phi_setting"])
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)
    InFile_DATA = TFile.Open(rootFileData, "OPEN")
    hist["InFile_DATA"]  = InFile_DATA
    
    rootFileDummy = OUTPATH + "/" + "{}".format(ParticleType) + "_" + inpDict["InDUMMYFilename"] + "_%s.root" % (hist["phi_setting"])
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)
    InFile_DUMMY = TFile.Open(rootFileDummy, "OPEN")
    hist["InFile_DUMMY"]  = InFile_DUMMY
    
prev_root_file = TFile.Open(prev_iter_root, "READ")
# Grab weight from previous iteration
for hist in histlist:
    hist.update(hist_in_dir(prev_root_file, "{}/data".format(hist["phi_setting"])))
    hist.update(hist_in_dir(prev_root_file, "{}/simc".format(hist["phi_setting"])))
    hist.update(hist_in_dir(prev_root_file, "{}/dummy".format(hist["phi_setting"])))
    
# t/phi bins are the same for all settings
# so arbitrarily grabbing from first setting of list
t_bins = np.array(histlist[0]["t_bins"])
phi_bins = np.array(histlist[0]["phi_bins"])

sys.path.append("binning")
from find_bins import find_bins, check_bins

#print("\n\nt_bins = ", t_bins)
#print("phi_bins = ", phi_bins)
check_bins(histlist, inpDict)

phisetlist = []
for hist in histlist:
    phisetlist.append(hist["phi_setting"])

for phiset in phisetlist:
    output_file_lst.append(OUTPATH+"/{}_{}_diamond_{}.pdf".format(phiset, ParticleType, 'Q'+Q2+'W'+W))
        
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(hist["phi_setting"],ParticleType)))

'''
EXAMPLE: How to get histograms from previous iteration

# Open the ROOT file, must pass open root file so object exists here and in function
prev_root_file = TFile.Open(prev_iter_root, "READ")
# Grab weight from previous iteration
iter_weight = get_histogram(prev_root_file, "{}/simc".format(hist["phi_setting"]), "H_Weight_SIMC")
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

# ***Removed effective charge calculation since it is grabbed from prev iter***

# Create a new directory for each iteration in cache
# ***Moved up in procedure vs main.py since required for weight iteration***
new_dir = "{}/{}/{}/Q{}W{}/{}".format(CACHEPATH, USER, ParticleType.lower(), Q2, W, formatted_date)
create_dir(new_dir)

# ***Moved from main.py location below because needed for weight iteration***
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
    # ***Create root directory here since it is used for weight iteration***
    create_dir(new_dir+"/root")
    # Make sure old simc root file exists
    if os.path.exists(old_simc_root):
        # Copy to new iteration so and then edit the weight
        print("Copying {} to {}".format(old_simc_root, new_simc_root))
        shutil.copy(old_simc_root,new_simc_root)
        shutil.copy(old_simc_hist,new_simc_hist)        
        # Make sure new simc root file exists
        if os.path.exists(new_simc_root):
            # Function to calculation new weight and apply it to simc root file 
            iter_weight(new_param_file, new_simc_root, inpDict, hist["phi_setting"])
            # Overwrite root file with updated weight
            os.rename(new_simc_root.replace(".root","_new.root"),new_simc_root)
            hist.update(compare_simc(new_simc_root, hist, inpDict))
        else:
            print("ERROR: {} not properly copied to {}".format(old_simc_root, new_simc_root))
            sys.exit(2)

if DEBUG:
    # Show plot pdf for each setting
    for hist in histlist:        
        show_pdf_with_evince(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(hist["phi_setting"],ParticleType)))
for hist in histlist:
    output_file_lst.append(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(hist["phi_setting"],ParticleType)))
            
sys.path.append("plotting")
from iter_check import plot_iteration

# Comparison plots of 0th to current iteration
plot_iteration(histlist, phisetlist, inpDict)

for hist in histlist:
    print("\n\n{} data total number of events: {:.3e}".format(hist["phi_setting"], hist["NumEvts_MM_DATA"]))
    print("{} dummy total number of events: {:.3e}".format(hist["phi_setting"], hist["NumEvts_MM_DUMMY"]))
    print("{} simc weighted total number of events: {:.3e}".format(hist["phi_setting"], hist["NumEvts_MM_SIMC"]))
    print("{} simc unweighted total number of events: {:.3e}".format(hist["phi_setting"], hist["NumEvts_MM_unweighted_SIMC"]))
    print("\n\n{} data normalization: {:.3e}".format(hist["phi_setting"], hist["normfac_data"]))
    print("{} dummy normalization: {:.3e}".format(hist["phi_setting"], hist["normfac_dummy"]))
    print("{} simc normalization: {:.3e}".format(hist["phi_setting"], hist["normfac_simc"]))

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))
#show_pdf_with_evince(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))    

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

# ***Grabbing data yield and average values from previous iteration rather than rebinning***
sys.path.append("binning")
from calculate_yield import grab_yield_data, find_yield_simc

yieldDict = {}
yieldDict.update(grab_yield_data(histlist, phisetlist, inpDict))
yieldDict.update(find_yield_simc(histlist, inpDict, iter_file=new_simc_root))

sys.path.append("binning")
from calculate_ratio import find_ratio

ratioDict = {}
ratioDict.update(find_ratio(histlist, inpDict, yieldDict))

sys.path.append("binning")
from ave_per_bin import grab_ave_data, ave_per_bin_simc

aveDict = {}
aveDict.update(grab_ave_data(histlist, inpDict))
aveDict.update(ave_per_bin_simc(histlist, inpDict, iter_file=new_simc_root))

sys.path.append("plotting")
from binned import plot_binned

plot_binned(t_bins, phi_bins, histlist, phisetlist, inpDict, yieldDict, ratioDict, aveDict)
#notify_email(email_address="trotta@cua.edu")

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))    

# Save histograms to root file
for hist in histlist:
    print("\nUpdating simc {} histograms in {}".format(hist["phi_setting"],foutroot))
    # Loop through all keys,values of dictionary
    for i, (key, val) in enumerate(hist.items()):
        # Progress bar
        Misc.progressBar(i, len(hist.items())-1,bar_length=25)
        if "G_data_eff" in key:
            hist_to_root(val, foutroot, "{}/data".format(hist["phi_setting"]))            
        if is_hist(val):
            if "ratio" in key:
                continue
            if "SIMC" in key:
                if "yield" in key:
                    continue
                elif "bin" in key:
                    continue
                elif "totevts" in key:
                    continue
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

# Redefinition from above, but should be the same! This is just to stay consistent with main.py
# ***Parameter files from last and this iteration!***
# FIX BELOW!!!
old_param_file = '{}/parameters/par.{}_Q{}W{}.dat'.format(prev_iter_dir, pol_str, Q2.replace("p",""), W.replace("p",""))

cut_summary_lst += "\nUnsep Parameterization for {}...".format(closest_date)
with open(old_param_file, 'r') as file:
    for line in file:
        cut_summary_lst += line
cut_summary_lst += "\nUnsep Parameterization for {}...".format(formatted_date)
with open(new_param_file, 'r') as file:
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
        print("1 ERROR: {}".format(e))
        sys.exit(2)

    if DEBUG:
        show_pdf_with_evince(OUTPATH+"/{}_lt_fit_Q{}W{}.pdf".format(ParticleType, Q2, W))
        show_pdf_with_evince(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))        
    output_file_lst.append(OUTPATH+"/{}_xsects_Q{}W{}.pdf".format(ParticleType, Q2, W))
    output_file_lst.append(OUTPATH+"/{}_lt_fit_Q{}W{}.pdf".format(ParticleType, Q2, W))
    
    # Save sep and unsep values from current iteration
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

# ***Moved creation of iteration directory up from where it is in main.py. Now is near the new weight calculation***
# ***Likewise for SIMC root/hist files***

if EPSSET == "high":
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
                print("Copying {} to {}".format(f,f_new))
                shutil.copy(f, f_new)
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
