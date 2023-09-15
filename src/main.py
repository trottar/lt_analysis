#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-14 23:27:25 trottar"
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
import csv
import shutil

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import show_pdf_with_evince, match_to_bin, create_dir

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

foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

output_file_lst = []

################################################################################################################################################    

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

** TODO: Add individual setting diamond plots
'''

#Importing diamond cut script
sys.path.append("cuts")
from diamond import DiamondPlot

Q2Val = float(Q2.replace("p","."))
WVal = float(W.replace("p","."))

inpDict["Q2min"] = Q2Val - (2/7)*Q2Val # Minimum value of Q2 on the Q2 vs W plot
inpDict["Q2max"] = Q2Val + (2/7)*Q2Val # Maximum value of Q2 on the Q2 vs W plot
inpDict["Wmin"] = WVal - (2/7)*WVal # min y-range for Q2vsW plot
inpDict["Wmax"] = WVal + (2/7)*WVal # max y-range for Q2vsW plot

phisetlist = ["Center","Left","Right"]
for phiset in phisetlist:
    # Call diamond cut script and append paramters to dictionary
    inpDict.update(DiamondPlot(ParticleType, Q2Val, inpDict["Q2min"], inpDict["Q2max"], WVal, inpDict["Wmin"], inpDict["Wmax"], phiset, tmin, tmax, inpDict))

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
from find_bins import find_bins

if EPSSET == "low":
    bin_vals = find_bins(histlist, inpDict)

try:
    with open("{}/out_data/{}/t_bin_interval".format(LTANAPATH, ParticleType), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            t_bins = all_lines[1].split("\t")
            del t_bins[0]
            t_bins = np.array([float(element) for element in t_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/out_data/{}/t_bin_interval".format(LTANAPATH, ParticleType)))
except IOError:
    print("Error reading {}...".format("{}/out_data/{}/t_bin_interval".format(LTANAPATH, ParticleType)))    

try:
    with open("{}/out_data/{}/phi_bin_interval".format(LTANAPATH, ParticleType), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            phi_bins = all_lines[1].split("\t")
            del phi_bins[0]
            phi_bins = np.array([float(element) for element in phi_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/out_data/{}/phi_bin_interval".format(LTANAPATH, ParticleType)))
except IOError:
    print("Error reading {}...".format("{}/out_data/{}/phi_bin_interval".format(LTANAPATH, ParticleType)))    
    
for hist in histlist:
    hist["t_bins"] = t_bins
    hist["phi_bins"] = phi_bins
    
##############################
# Step 5 of the lt_analysis: # DONE
##############################
'''
* Compare SIMC to data/dummy setting by setting.

* SIMC weight may need to be recalculated.

* If so, use script...

> lt_analysis/src/SIMC/??????????????.f

** TODO: Fix plots (e.g. polar) and find working simc weight script
'''

sys.path.append("normalize")
from get_eff_charge import get_eff_charge

# Upate hist dictionary with effective charge
for hist in histlist:
    hist.update(get_eff_charge(hist, inpDict))

sys.path.append("simc_ana")    
from compare_simc import compare_simc

sys.path.append("plotting")
from data_vs_simc import plot_data_vs_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    hist.update(compare_simc(hist, inpDict))

cut_summary_lst = plot_data_vs_simc(histlist, outputpdf)

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

histbinDict = {}

# t/phi binned histograms
histbinDict["H_phibins_DATA"] = ROOT.TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
histbinDict["H_tbins_DATA"] = ROOT.TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)

# Yield histograms
for phiset in phisetlist:
    histbinDict["H_yield_DATA_{}".format(phiset)] = ROOT.TH1D("H_yield_DATA_{}".format(phiset), "{} Data Yield".format(phiset), NumtBins*NumPhiBins, 0, 100.0)
    histbinDict["H_yield_SIMC_{}".format(phiset)] = ROOT.TH1D("H_yield_SIMC_{}".format(phiset), "{} Simc Yield".format(phiset), NumtBins*NumPhiBins, 0, 100.0)

# Ratio histogram
for phiset in phisetlist:
    histbinDict["H_ratio_{}".format(phiset)] = ROOT.TH1D("H_ratio_{}".format(phiset), "{} Ratio".format(phiset), NumtBins*NumPhiBins, -100.0, 100.0)
    
# Loop over each tuple key in the dictionary
for phiset in phisetlist:
    for k, data_key_tuple in enumerate(aveDict["binned_DATA"][phiset]['t']):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_Q2_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data Q2 (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Q2min"], inpDict["Q2max"])
        histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_W_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data W (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Wmin"], inpDict["Wmax"])
        histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_t_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data t (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["tmin"], inpDict["tmax"])

# Loop over each tuple key in the dictionary
for phiset in phisetlist:
    for k, simc_key_tuple in enumerate(aveDict["binned_SIMC"][phiset]['t']):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = aveDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin
        histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_Q2_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc Q2 (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Q2min"], inpDict["Q2max"])
        histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_W_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc W (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Wmin"], inpDict["Wmax"])
        histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_t_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc t (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["tmin"], inpDict["tmax"])

# Loop over each tuple key in the dictionary
for phiset in phisetlist:
    for k, data_key_tuple in enumerate(yieldDict["binned_DATA"][phiset]['yield']):
        # Access the nested dictionary using the tuple key
        data_nested_dict = yieldDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_totevts_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data Total Events (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, 0.0, 100.0)

# Loop over each tuple key in the dictionary
for phiset in phisetlist:
    for k, simc_key_tuple in enumerate(yieldDict["binned_SIMC"][phiset]['yield']):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = yieldDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin
        histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = ROOT.TH1D("H_totevts_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc Total Events (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, 0.0, 100.0)
        
C_Q2_tbin_DATA = TCanvas()
C_Q2_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["Q2"][data_key_tuple]["Q2_arr"]), data_nested_dict["Q2"][data_key_tuple]["Q2_ave"]))
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(data_nested_dict["Q2"][data_key_tuple]["Q2_arr"]):
            histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_Q2_tbin_DATA.cd(k+1)
        histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_Q2_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType))+'(')

C_Q2_tbin_SIMC = TCanvas()
C_Q2_tbin_SIMC.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
    for k, simc_key_tuple in enumerate(simc_key_tuples):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = aveDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin    
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(simc_nested_dict["Q2"][simc_key_tuple]["Q2_arr"]):
            histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_Q2_tbin_SIMC.cd(k+1)
        histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_Q2_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_W_tbin_DATA = TCanvas()
C_W_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin    
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(data_nested_dict["W"][data_key_tuple]["W_arr"]):
            histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_W_tbin_DATA.cd(k+1)
        histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_W_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_W_tbin_SIMC = TCanvas()
C_W_tbin_SIMC.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
    for k, simc_key_tuple in enumerate(simc_key_tuples):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = aveDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin    
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(simc_nested_dict["W"][simc_key_tuple]["W_arr"]):
            histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_W_tbin_SIMC.cd(k+1)
        histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_W_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_t_tbin_DATA = TCanvas()
C_t_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(data_nested_dict["t"][data_key_tuple]["t_arr"]):
            histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_t_tbin_DATA.cd(k+1)
        histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_t_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_t_tbin_SIMC = TCanvas()
C_t_tbin_SIMC.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
    for k, simc_key_tuple in enumerate(simc_key_tuples):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = aveDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin    
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(simc_nested_dict["t"][simc_key_tuple]["t_arr"]):
            histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_t_tbin_SIMC.cd(k+1)
        histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_t_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
    
C_t_bins_DATA = TCanvas()
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        # Fill histogram
        for val in aveDict["binned_DATA"]["t_bins"]:
            histbinDict["H_tbins_DATA"].Fill(float(val))
    histbinDict["H_tbins_DATA"].SetLineColor(it+1)
    histbinDict["H_tbins_DATA"].Draw("same, hist")
C_t_bins_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_phi_bins_DATA = TCanvas()
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        # Fill histogram
        for val in aveDict["binned_DATA"]["phi_bins"]:
            histbinDict["H_phibins_DATA"].Fill(float(val))
    histbinDict["H_phibins_DATA"].SetLineColor(it+1)            
    histbinDict["H_phibins_DATA"].Draw("same, hist")
C_phi_bins_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_totevts_DATA = TCanvas()
C_totevts_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = yieldDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["yield"][data_key_tuple]["yield_arr"]), data_nested_dict["yield"][data_key_tuple]["yield"]))
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(data_nested_dict["yield"][data_key_tuple]["yield_arr"]):
            histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_totevts_DATA.cd(k+1)
        histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_totevts_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_totevts_SIMC = TCanvas()
C_totevts_SIMC.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
    for k, simc_key_tuple in enumerate(simc_key_tuples):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = yieldDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield_arr"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
        # Fill histogram
        for (itt,jtt), val in np.ndenumerate(simc_nested_dict["yield"][simc_key_tuple]["yield_arr"]):
            histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
        C_totevts_SIMC.cd(k+1)
        histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
        histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
        histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
C_totevts_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

for it,phiset in enumerate(phisetlist):
    
    C_yieldvsphi_data_plt = TCanvas()
    C_yieldvsphi_data_plt.SetGrid()
    C_yieldvsphi_data_plt.Divide(1,NumtBins)
    l_yieldvsphi_data_plt = ROOT.TLegend(0.115,0.35,0.33,0.5)

    yield_data = []
    yield_simc = []
    phibins_data = []
    phibins_simc = []    
    data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
    simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
    for data_key_tuple,simc_key_tuple in zip(data_key_tuples,simc_key_tuples):
        tmp_yield_data = [[],[]]
        tmp_yield_simc = [[],[]]
        tmp_phibins_data = [[],[]]
        tmp_phibins_simc = [[],[]]
        # Access the nested dictionary using the tuple key
        data_nested_dict = yieldDict["binned_DATA"][phiset]
        simc_nested_dict = yieldDict["binned_SIMC"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        tmp_yield_data[0].append(yieldDict["binned_DATA"]["t_bins"][i])
        tmp_yield_data[1].append(data_nested_dict["yield"][data_key_tuple]["yield"])
        tmp_yield_simc[0].append(yieldDict["binned_SIMC"]["t_bins"][i])
        tmp_yield_simc[1].append(simc_nested_dict["yield"][simc_key_tuple]["yield"])
        tmp_phibins_data[0].append(yieldDict["binned_DATA"]["t_bins"][i])
        tmp_phibins_data[1].append(yieldDict["binned_DATA"]["phi_bins"][j])
        tmp_phibins_simc[0].append(yieldDict["binned_SIMC"]["t_bins"][i])
        tmp_phibins_simc[1].append(yieldDict["binned_SIMC"]["phi_bins"][j])
        yield_data.append(tmp_yield_data)
        yield_simc.append(tmp_yield_simc)
        phibins_data.append(tmp_phibins_data)
        phibins_simc.append(tmp_phibins_simc)

    # Match t-bins with list of yields
    yield_data = match_to_bin(yield_data)
    yield_simc = match_to_bin(yield_simc)
    phibins_data = match_to_bin(phibins_data)
    phibins_simc = match_to_bin(phibins_simc)

    multiDict = {}
    for i, val in enumerate(t_bins):

        multiDict["G_yieldvsphi_plt_{}".format(i)] = ROOT.TMultiGraph()

        G_yieldvsphi_data = ROOT.TGraphErrors(len(yield_data[i][1]),phibins_data[i][1],yield_data[i][1],np.array([0]*len(phibins_data[i][1])),np.array([0]*len(yield_data[i][1])))
        G_yieldvsphi_simc = ROOT.TGraphErrors(len(yield_simc[i][1]),phibins_simc[i][1],yield_simc[i][1],np.array([0]*len(phibins_simc[i][1])),np.array([0]*len(yield_simc[i][1])))

        G_yieldvsphi_data.SetMarkerStyle(21)
        G_yieldvsphi_data.SetMarkerSize(1)
        G_yieldvsphi_data.SetMarkerColor(1)
        multiDict["G_yieldvsphi_plt_{}".format(i)].Add(G_yieldvsphi_data)

        G_yieldvsphi_simc.SetMarkerStyle(21)
        G_yieldvsphi_simc.SetMarkerSize(1)
        G_yieldvsphi_simc.SetMarkerColor(2)
        multiDict["G_yieldvsphi_plt_{}".format(i)].Add(G_yieldvsphi_simc)

        C_yieldvsphi_data_plt.cd(i+1)

        multiDict["G_yieldvsphi_plt_{}".format(i)].Draw("AP")
        multiDict["G_yieldvsphi_plt_{}".format(i)].SetTitle("{} t = {};#phi; Yield".format(phiset,val))

        multiDict["G_yieldvsphi_plt_{}".format(i)].GetYaxis().SetTitleOffset(1.5)
        multiDict["G_yieldvsphi_plt_{}".format(i)].GetXaxis().SetTitleOffset(1.5)
        multiDict["G_yieldvsphi_plt_{}".format(i)].GetXaxis().SetLabelSize(0.04)

    l_yieldvsphi_data_plt.AddEntry(G_yieldvsphi_data,"Data")
    l_yieldvsphi_data_plt.AddEntry(G_yieldvsphi_simc,"Simc")
    l_yieldvsphi_data_plt.Draw()

    C_yieldvsphi_data_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

for it,phiset in enumerate(phisetlist):
    
    C_ratiovsphi_plt = TCanvas()
    C_ratiovsphi_plt.SetGrid()
    C_ratiovsphi_plt.Divide(1,NumtBins)

    ratio = []
    phibins = []
    key_tuples = list(ratioDict["binned"][phiset]['ratio'])
    for key_tuple in key_tuples:
        tmp_ratio = [[],[]]
        tmp_phibins = [[],[]]
        # Access the nested dictionary using the tuple key
        nested_dict = ratioDict["binned"][phiset]
        i = key_tuple[0] # t bin
        j = key_tuple[1] # phi bin
        tmp_ratio[0].append(ratioDict["binned"]["t_bins"][i])
        tmp_ratio[1].append(nested_dict["ratio"][key_tuple]["ratio"])
        tmp_phibins[0].append(ratioDict["binned"]["t_bins"][i])
        tmp_phibins[1].append(ratioDict["binned"]["phi_bins"][j])
        ratio.append(tmp_ratio)
        phibins.append(tmp_phibins)
        
    # Match t-bins with list of ratios
    ratio = match_to_bin(ratio)
    phibins = match_to_bin(phibins)

    multiDict = {}
    for i, val in enumerate(t_bins):

        multiDict["G_ratiovsphi_plt_{}".format(i)] = ROOT.TMultiGraph()

        G_ratiovsphi = ROOT.TGraphErrors(len(ratio[i][1]),phibins[i][1],ratio[i][1],np.array([0]*len(phibins[i][1])),np.array([0]*len(ratio[i][1])))

        G_ratiovsphi.SetMarkerStyle(21)
        G_ratiovsphi.SetMarkerSize(1)
        G_ratiovsphi.SetMarkerColor(1)
        multiDict["G_ratiovsphi_plt_{}".format(i)].Add(G_ratiovsphi)

        C_ratiovsphi_plt.cd(i+1)

        multiDict["G_ratiovsphi_plt_{}".format(i)].Draw("AP")
        multiDict["G_ratiovsphi_plt_{}".format(i)].SetTitle("{} t = {};#phi; Ratio".format(phiset,val))

        multiDict["G_ratiovsphi_plt_{}".format(i)].GetYaxis().SetTitleOffset(1.5)
        multiDict["G_ratiovsphi_plt_{}".format(i)].GetXaxis().SetTitleOffset(1.5)
        multiDict["G_ratiovsphi_plt_{}".format(i)].GetXaxis().SetLabelSize(0.04)

    C_ratiovsphi_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_yield_DATA = TCanvas()
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
    for k, data_key_tuple in enumerate(data_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = yieldDict["binned_DATA"][phiset]
        i = data_key_tuple[0] # t bin
        j = data_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["yield"][data_key_tuple]["yield"]), data_nested_dict["yield"][data_key_tuple]["yield"]))
        # Fill histogram
        val = data_nested_dict["yield"][data_key_tuple]["yield"]
        histbinDict["H_yield_DATA_{}".format(phiset)].Fill(val)
    histbinDict["H_yield_DATA_{}".format(phiset)].SetLineColor(it+1)
    histbinDict["H_yield_DATA_{}".format(phiset)].Draw("same, E1")
    histbinDict["H_yield_DATA_{}".format(phiset)].Draw("same, hist")
C_yield_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_yield_SIMC = TCanvas()
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
    for k, simc_key_tuple in enumerate(simc_key_tuples):
        # Access the nested dictionary using the tuple key
        simc_nested_dict = yieldDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
        # Fill histogram
        val = simc_nested_dict["yield"][simc_key_tuple]["yield"]
        histbinDict["H_yield_SIMC_{}".format(phiset)].Fill(val)
    histbinDict["H_yield_SIMC_{}".format(phiset)].SetLineColor(it+1)
    histbinDict["H_yield_SIMC_{}".format(phiset)].Draw("same, E1")
    histbinDict["H_yield_SIMC_{}".format(phiset)].Draw("same, hist")
C_yield_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_ratio = TCanvas()
# Loop over each tuple key in the dictionary
for it,phiset in enumerate(phisetlist):
    key_tuples = list(ratioDict["binned"][phiset]['ratio'])
    for k, key_tuple in enumerate(key_tuples):
        # Access the nested dictionary using the tuple key
        nested_dict = ratioDict["binned"][phiset]
        i = key_tuple[0] # t bin
        j = key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(nested_dict["ratio"][key_tuple]["ratio"]), nested_dict["ratio"][key_tuple]["ratio"]))
        # Fill histogram
        val = nested_dict["ratio"][key_tuple]["ratio"]
        histbinDict["H_ratio_{}".format(phiset)].Fill(val)
    histbinDict["H_ratio_{}".format(phiset)].SetLineColor(it+1)
    histbinDict["H_ratio_{}".format(phiset)].Draw("same, E1")
    histbinDict["H_ratio_{}".format(phiset)].Draw("same, hist")
C_ratio.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

C_yield_data_plt = TCanvas()
G_yield_data_plt = ROOT.TMultiGraph()
l_yield_data_plt = ROOT.TLegend(0.115,0.35,0.33,0.5)

C_yield_data_plt.SetGrid()

yield_data = np.array([])
yield_simc = np.array([])
setting = np.array([])
for it,phiset in enumerate(phisetlist):
    data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
    simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
    for data_key_tuple,simc_key_tuple in zip(data_key_tuples,simc_key_tuples):
        # Access the nested dictionary using the tuple key
        data_nested_dict = yieldDict["binned_DATA"][phiset]        
        simc_nested_dict = yieldDict["binned_SIMC"][phiset]
        i = simc_key_tuple[0] # t bin
        j = simc_key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
        # Fill histogram
        yield_data = np.append(yield_data, [data_nested_dict["yield"][data_key_tuple]["yield"]])        
        yield_simc = np.append(yield_simc, [simc_nested_dict["yield"][simc_key_tuple]["yield"]])
        if phiset == "Center": setting = np.append(setting,0)
        elif phiset == "Left": setting = np.append(setting,1)
        else: setting = np.append(setting,2)
        
G_yield_data = ROOT.TGraphErrors(len(yield_data),setting,yield_data,np.array([0]*len(setting)),np.array([0]*len(yield_data)))
G_yield_simc = ROOT.TGraphErrors(len(yield_simc),setting,yield_simc,np.array([0]*len(setting)),np.array([0]*len(yield_simc)))

G_yield_data.SetMarkerStyle(21)
G_yield_data.SetMarkerSize(1)
G_yield_data.SetMarkerColor(1)
G_yield_data_plt.Add(G_yield_data)

G_yield_simc.SetMarkerStyle(21)
G_yield_simc.SetMarkerSize(1)
G_yield_simc.SetMarkerColor(2)
G_yield_data_plt.Add(G_yield_simc)

G_yield_data_plt.Draw("AP")
G_yield_data_plt.SetTitle(" ;Setting; Yield")

for i,hist in enumerate(histlist):
    bin_ix = G_yield_data_plt.GetXaxis().FindBin(i)
    G_yield_data_plt.GetXaxis().SetBinLabel(bin_ix,hist["phi_setting"])

G_yield_data_plt.GetYaxis().SetTitleOffset(1.5)
G_yield_data_plt.GetXaxis().SetTitleOffset(1.5)
G_yield_data_plt.GetXaxis().SetLabelSize(0.04)

l_yield_data_plt.AddEntry(G_yield_data,"Data")
l_yield_data_plt.AddEntry(G_yield_simc,"Simc")
l_yield_data_plt.Draw()

C_yield_data_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
    
C_ratio_plt = TCanvas()
G_ratio_plt = ROOT.TMultiGraph()

C_ratio_plt.SetGrid()

ratio_data = np.array([])
setting = np.array([])
for it,phiset in enumerate(phisetlist):
    key_tuples = list(ratioDict["binned"][phiset]['ratio'])
    for k, key_tuple in enumerate(key_tuples):
        # Access the nested dictionary using the tuple key
        nested_dict = ratioDict["binned"][phiset]
        i = key_tuple[0] # t bin
        j = key_tuple[1] # phi bin
        #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(nested_dict["ratio"][key_tuple]["ratio"]), nested_dict["ratio"][key_tuple]["ratio"]))
        ratio_data = np.append(ratio_data, [nested_dict["ratio"][key_tuple]["ratio"]])
        if phiset == "Center": setting = np.append(setting,0)
        elif phiset == "Left": setting = np.append(setting,1)
        else: setting = np.append(setting,2)

G_ratio = ROOT.TGraphErrors(len(ratio_data),setting,ratio_data,np.array([0]*len(setting)),np.array([0]*len(ratio_data)))

G_ratio.SetMarkerStyle(21)
G_ratio.SetMarkerSize(1)
G_ratio.SetMarkerColor(1)
G_ratio_plt.Add(G_ratio)

G_ratio_plt.Draw("AP")

G_ratio_plt.SetTitle(" ;Setting; Ratio")

for i,hist in enumerate(histlist):
    bin_ix = G_ratio_plt.GetXaxis().FindBin(i)
    G_ratio_plt.GetXaxis().SetBinLabel(bin_ix,hist["phi_setting"])

G_ratio_plt.GetYaxis().SetTitleOffset(1.5)
G_ratio_plt.GetXaxis().SetTitleOffset(1.5)
G_ratio_plt.GetXaxis().SetLabelSize(0.04)

C_ratio_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType))+')')

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
output_file_lst.append(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))
    
# Run fortran script

from physics_lists import create_lists
create_lists(aveDict, ratioDict, histlist, inpDict, phisetlist, output_file_lst)

'''
*****************************************
* NEED TO ADD ROOT FILES FOR OTHERS USE *
*****************************************
'''

# ***Parameter file from last iteration!***
# ***These old parameters are needed for this iteration. See README for more info on procedure!***
old_param_file = '{}/out_data/{}/parameters/par.{}_{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""))
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
print("{} Cut Summary...".format(EPSSET.capitalize()),cut_summary_lst)
print("="*25)
inpDict["cut_summary_lst"] = cut_summary_lst

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data and SIMC.

* Calculate the unseparated cross section
'''

if EPSSET == "high":
    subprocess.call(['bash','{}/run_xsect.sh'.format(LTANAPATH), Q2, W])

    # Save new parameters and unsep values from current iteration
    # ***Old parameter file defined in step 7, the new parameter values are saved here!***
    # ***The old parameters, used for this iteration, are saved in the summary!***
    new_param_file = '{}/parameters/par.{}_{}.dat'.format(ParticleType, pol_str, Q2.replace("p",""))
    output_file_lst.append(new_param_file) 
    xsect_file = '{}/xsects/x_unsep.{}_{}_{:.0f}.dat'.format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(xsect_file)

    # Save fortran scripts that contain iteration functional form of parameterization
    fort_param = 'param_{}_{}.f'.format(ParticleType, pol_str)
    output_file_lst.append(fort_param) 
    fort_xmodel = 'xmodel_{}_{}.f'.format(ParticleType, pol_str)
    output_file_lst.append(fort_xmodel) 

##############################
# Step 8 of the lt_analysis: #
##############################
'''
* Save all information for this iteration
'''

if EPSSET == "high":
    
    print("\n\n")

    # Create a new directory for each iteration
    new_dir = CACHEPATH+"/"+USER+"/"+ParticleType.lower()+"/"+formatted_date
    create_dir(new_dir)

    for f in output_file_lst:
        if OUTPATH in f:
            if ".pdf" in f:
                create_dir(new_dir+"/plots")
                f_new = f.replace(OUTPATH,new_dir+"/plots")
                print("Copying {} to {}".format(f,f_new))
                shutil.copy(f, f_new)            
            if ".root" in f:
                create_dir(new_dir+"/rootfiles")
                f_new = f.replace(OUTPATH,new_dir+"/rootfiles")
                print("Copying {} to {}".format(f,f_new))
                shutil.copy(f, f_new)
        elif "{}/".format(ParticleType) in f:
            f_arr = f.split("/")
            f_tmp = f_arr.pop()
            for f_dir in f_arr:
                if "{}".format(ParticleType) not in f_dir:
                    create_dir(new_dir+"/"+f_dir)
                    f_new = new_dir+"/"+f_dir+"/"+f_tmp    
                    print("Copying {} to {}".format(LTANAPATH+"/out_data/"+f,f_new))
                    shutil.copy(LTANAPATH+"/out_data/"+f, f_new)
        else:
            f_new = new_dir
            print("Copying {} to {}".format(LTANAPATH+"/out_data/"+f,f_new))
            shutil.copy(LTANAPATH+"/out_data/"+f, f_new)

    with open(new_dir+'/{}_{}_summary_{}.txt'.format(ParticleType,OutFilename,formatted_date), 'w') as file:
        file.write(inpDict["cut_summary_lst"])
