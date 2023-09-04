#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-04 19:18:57 trottar"
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

##################################################################################################################################################
# Importing utility functions

from utility import show_pdf_with_evince

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=34:
    print("!!!!! ERROR !!!!!\n Expected 34 arguments\n Usage is with - KIN W Q2 EPSVAL OutDATAFilename OutDUMMYFilename OutFullAnalysisFilename tmin tmax NumtBins NumPhiBins runNumRight runNumLeft runNumCenter data_charge_right data_charge_left data_charge_center dummy_charge_right dummy_charge_left dummy_charge_center InData_efficiency_right InData_efficiency_left InData_efficiency_center efficiency_table ParticleType EPSSET pThetaValRight pThetaValLeft pThetaValCenter EbeamValRight EbeamValLeft EbeamValCenter POL KSet\n!!!!! ERROR !!!!!")
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
KSet = sys.argv[34]

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
    "KSet" : KSet,
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
# Step 1 of the lt_analysis: # DONE
##############################
'''
* All analysis cuts per run (i.e. PID, acceptance, timing) are applied. 
* This should have been completed before this script using...

> lt_analysis/scripts/Prod/analyze_prod.py

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
        
##############################
# Step 3 of the lt_analysis: # DONE
##############################
'''
Apply random subtraction to data and dummy.
'''

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

##############################
# Step 4 of the lt_analysis: # Done
##############################
'''
* Combine all settings and choose t/phi bins for low eps.

* These bins will also be used of high eps, so check high eps as well.
'''

from find_bins import find_bins

if EPSSET == "low":
    bin_vals = find_bins(histlist, inpDict)

try:
    with open("{}/src/t_bin_interval".format(LTANAPATH), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            t_bins = all_lines[1].split("\t")
            del t_bins[0]
            t_bins = np.array([float(element) for element in t_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/t_bin_interval".format(LTANAPATH)))
except IOError:
    print("Error reading {}...".format("{}/src/t_bin_interval".format(LTANAPATH)))    

try:
    with open("{}/src/phi_bin_interval".format(LTANAPATH), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            phi_bins = all_lines[1].split("\t")
            del phi_bins[0]
            phi_bins = np.array([float(element) for element in phi_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/phi_bin_interval".format(LTANAPATH)))
except IOError:
    print("Error reading {}...".format("{}/src/phi_bin_interval".format(LTANAPATH)))    
    
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

from get_eff_charge import get_eff_charge

# Upate hist dictionary with effective charge
for hist in histlist:
    hist.update(get_eff_charge(hist, inpDict))

from compare_simc import compare_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    hist.update(compare_simc(hist, inpDict))

eff_plt = TCanvas()
G_eff_plt = ROOT.TMultiGraph()
l_eff_plt = ROOT.TLegend(0.115,0.35,0.33,0.5)

eff_plt.SetGrid()

for i,hist in enumerate(histlist):
    hist["G_data_eff"].SetMarkerStyle(21)
    hist["G_data_eff"].SetMarkerSize(1)
    hist["G_data_eff"].SetMarkerColor(i+1)
    G_eff_plt.Add(hist["G_data_eff"])

G_eff_plt.Draw("AP")

G_eff_plt.SetTitle(" ;Run Numbers; Total Efficiency")

i=0
for i,hist in enumerate(histlist):
    while i <= G_eff_plt.GetXaxis().GetXmax():
        bin_ix = G_eff_plt.GetXaxis().FindBin(i)
        if str(i) in hist["runNums"]: 
            G_eff_plt.GetXaxis().SetBinLabel(bin_ix,"%d" % i)
        i+=1

G_eff_plt.GetYaxis().SetTitleOffset(1.5)
G_eff_plt.GetXaxis().SetTitleOffset(1.5)
G_eff_plt.GetXaxis().SetLabelSize(0.04)

for i,hist in enumerate(histlist):
    l_eff_plt.AddEntry(hist["G_data_eff"],hist["phi_setting"])

l_eff_plt.Draw()

eff_plt.Print(outputpdf + '(')
    
# Plot histograms
c_pid = TCanvas()

c_pid.Divide(2,3)

c_pid.cd(1)
gPad.SetLogy()

for i,hist in enumerate(histlist):
    hist["H_cal_etottracknorm_DATA"].SetLineColor(i+1)
    hist["H_cal_etottracknorm_DATA"].Draw("same, E1")

c_pid.cd(2)
gPad.SetLogy()

for i,hist in enumerate(histlist):
    hist["H_cer_npeSum_DATA"].SetLineColor(i+1)
    hist["H_cer_npeSum_DATA"].Draw("same, E1")

c_pid.cd(3)
gPad.SetLogy()
for i,hist in enumerate(histlist):
    hist["P_cal_etottracknorm_DATA"].SetLineColor(i+1)
    hist["P_cal_etottracknorm_DATA"].Draw("same, E1")

c_pid.cd(4)
gPad.SetLogy()
for i,hist in enumerate(histlist):
    hist["P_hgcer_npeSum_DATA"].SetLineColor(i+1)
    hist["P_hgcer_npeSum_DATA"].Draw("same, E1")

c_pid.cd(5)
gPad.SetLogy()
for i,hist in enumerate(histlist):
    hist["P_aero_npeSum_DATA"].SetLineColor(i+1)
    hist["P_aero_npeSum_DATA"].Draw("same, E1")
        
c_pid.Draw()

c_pid.Print(outputpdf)

ct = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ct_DATA"].SetLineColor(i+1)
    hist["H_ct_DATA"].Draw("same, E1")

ct.Print(outputpdf)


CQ2 = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_Q2_DATA"].SetLineColor(i+1)
    hist["H_Q2_DATA"].Draw("same, E1")
    hist["H_Q2_SIMC"].SetLineColor(40)
    hist["H_Q2_SIMC"].SetLineStyle(10-i)
    hist["H_Q2_SIMC"].Draw("same, E1")    

CQ2.Print(outputpdf)

CW = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_W_DATA"].SetLineColor(i+1)
    hist["H_W_DATA"].Draw("same, E1")
    hist["H_W_SIMC"].SetLineColor(40)
    hist["H_W_SIMC"].SetLineStyle(10-i)
    hist["H_W_SIMC"].Draw("same, E1")        
    
CW.Print(outputpdf)

Ct = TCanvas()
l_t = ROOT.TLegend(0.115,0.45,0.33,0.95)
l_t.SetTextSize(0.0235)

binmax = []
for i,hist in enumerate(histlist):
    hist["H_t_DATA"].SetLineColor(i+1)
    l_t.AddEntry(hist["H_t_DATA"],hist["phi_setting"])
    hist["H_t_DATA"].Draw("same, E1")
    hist["H_t_SIMC"].SetLineColor(40)
    hist["H_t_SIMC"].SetLineStyle(10-i)
    hist["H_t_SIMC"].Draw("same, E1")
    binmax.append(hist["H_t_DATA"].GetMaximum())
binmax = max(binmax)
    
tBin_line = TLine()
for i,b in enumerate(t_bins):
    b = float(b)
    tBin_line.SetLineColor(4)
    tBin_line.SetLineWidth(4)
    tBin_line.DrawLine(b,0,b,binmax)
    l_t.AddEntry(tBin_line,"Bin Edge %s" % i )
    l_t.AddEntry(tBin_line,"BinCenter = %.2f" % b)

l_t.Draw()    

Ct.Print(outputpdf)

Cepsilon = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_epsilon_DATA"].SetLineColor(i+1)
    hist["H_epsilon_DATA"].Draw("same, E1")
    hist["H_epsilon_SIMC"].SetLineColor(40)
    hist["H_epsilon_SIMC"].SetLineStyle(10-i)
    hist["H_epsilon_SIMC"].Draw("same, E1")
    
Cepsilon.Print(outputpdf)

CMM = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_MM_DATA"].SetLineColor(i+1)
    hist["H_MM_DATA"].Draw("same, E1")
    hist["H_MM_SIMC"].SetLineColor(40)
    hist["H_MM_SIMC"].SetLineStyle(10-i)
    hist["H_MM_SIMC"].Draw("same, E1")
    
CMM.Print(outputpdf)

xfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssxfp_DATA"].SetLineColor(i+1)
    hist["H_ssxfp_DATA"].Draw("same, E1")
    hist["H_ssxfp_SIMC"].SetLineColor(40)
    hist["H_ssxfp_SIMC"].SetLineStyle(10-i)
    hist["H_ssxfp_SIMC"].Draw("same, E1")
    
xfp.Print(outputpdf)

yfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssyfp_DATA"].SetLineColor(i+1)
    hist["H_ssyfp_DATA"].Draw("same, E1")
    hist["H_ssyfp_SIMC"].SetLineColor(40)
    hist["H_ssyfp_SIMC"].SetLineStyle(10-i)
    hist["H_ssyfp_SIMC"].Draw("same, E1")
    
yfp.Print(outputpdf)

xpfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssxpfp_DATA"].SetLineColor(i+1)
    hist["H_ssxpfp_DATA"].Draw("same, E1")
    hist["H_ssxpfp_SIMC"].SetLineColor(40)
    hist["H_ssxpfp_SIMC"].SetLineStyle(10-i)
    hist["H_ssxpfp_SIMC"].Draw("same, E1")
    
xpfp.Print(outputpdf)

ypfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssxpfp_DATA"].SetLineColor(i+1)
    hist["H_ssxpfp_DATA"].Draw("same, E1")
    hist["H_ssxpfp_SIMC"].SetLineColor(40)
    hist["H_ssxpfp_SIMC"].SetLineStyle(10-i)
    hist["H_ssxpfp_SIMC"].Draw("same, E1")
    
ypfp.Print(outputpdf)

hxfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsxfp_DATA"].SetLineColor(i+1)
    hist["H_hsxfp_DATA"].Draw("same, E1")
    hist["H_hsxfp_SIMC"].SetLineColor(40)
    hist["H_hsxfp_SIMC"].SetLineStyle(10-i)
    hist["H_hsxfp_SIMC"].Draw("same, E1")
    
hxfp.Print(outputpdf)

hyfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsyfp_DATA"].SetLineColor(i+1)
    hist["H_hsyfp_DATA"].Draw("same, E1")
    hist["H_hsyfp_SIMC"].SetLineColor(40)
    hist["H_hsyfp_SIMC"].SetLineStyle(10-i)
    hist["H_hsyfp_SIMC"].Draw("same, E1")
    
hyfp.Print(outputpdf)

hxpfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsxpfp_DATA"].SetLineColor(i+1)
    hist["H_hsxpfp_DATA"].Draw("same, E1")
    hist["H_hsxpfp_SIMC"].SetLineColor(40)
    hist["H_hsxpfp_SIMC"].SetLineStyle(10-i)
    hist["H_hsxpfp_SIMC"].Draw("same, E1")
    
hxpfp.Print(outputpdf)

hypfp = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsypfp_DATA"].SetLineColor(i+1)
    hist["H_hsypfp_DATA"].Draw("same, E1")
    hist["H_hsypfp_SIMC"].SetLineColor(40)
    hist["H_hsypfp_SIMC"].SetLineStyle(10-i)
    hist["H_hsypfp_SIMC"].Draw("same, E1")
    
hypfp.Print(outputpdf)

xptar = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssxptar_DATA"].SetLineColor(i+1)
    hist["H_ssxptar_DATA"].Draw("same, E1")
    hist["H_ssxptar_SIMC"].SetLineColor(40)
    hist["H_ssxptar_SIMC"].SetLineStyle(10-i)
    hist["H_ssxptar_SIMC"].Draw("same, E1")
    
xptar.Print(outputpdf)

yptar = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssyptar_DATA"].SetLineColor(i+1)
    hist["H_ssyptar_DATA"].Draw("same, E1")
    hist["H_ssyptar_SIMC"].SetLineColor(40)
    hist["H_ssyptar_SIMC"].SetLineStyle(10-i)
    hist["H_ssyptar_SIMC"].Draw("same, E1")
    
yptar.Print(outputpdf)

hxptar = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsxptar_DATA"].SetLineColor(i+1)
    hist["H_hsxptar_DATA"].Draw("same, E1")
    hist["H_hsxptar_SIMC"].SetLineColor(40)
    hist["H_hsxptar_SIMC"].SetLineStyle(10-i)
    hist["H_hsxptar_SIMC"].Draw("same, E1")
    
hxptar.Print(outputpdf)

hyptar = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsyptar_DATA"].SetLineColor(i+1)
    hist["H_hsyptar_DATA"].Draw("same, E1")
    hist["H_hsyptar_SIMC"].SetLineColor(40)
    hist["H_hsyptar_SIMC"].SetLineStyle(10-i)
    hist["H_hsyptar_SIMC"].Draw("same, E1")
    
hyptar.Print(outputpdf)

Delta = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ssdelta_DATA"].SetLineColor(i+1)
    hist["H_ssdelta_DATA"].Draw("same, E1")
    hist["H_ssdelta_SIMC"].SetLineColor(40)
    hist["H_ssdelta_SIMC"].SetLineStyle(10-i)
    hist["H_ssdelta_SIMC"].Draw("same, E1")
    
Delta.Print(outputpdf)

hDelta = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_hsdelta_DATA"].SetLineColor(i+1)
    hist["H_hsdelta_DATA"].Draw("same, E1")
    hist["H_hsdelta_SIMC"].SetLineColor(40)
    hist["H_hsdelta_SIMC"].SetLineStyle(10-i)
    hist["H_hsdelta_SIMC"].Draw("same, E1")
    
hDelta.Print(outputpdf)

Cph_q = TCanvas()

binmax = []
for i,hist in enumerate(histlist):
    hist["H_ph_q_DATA"].SetLineColor(i+1)
    l_t.AddEntry(hist["H_ph_q_DATA"],hist["phi_setting"])
    hist["H_ph_q_DATA"].Draw("same, E1")
    hist["H_ph_q_SIMC"].SetLineColor(40)
    hist["H_ph_q_SIMC"].SetLineStyle(10-i)
    hist["H_ph_q_SIMC"].Draw("same, E1")
    binmax.append(hist["H_ph_q_DATA"].GetMaximum())
binmax = max(binmax)

binned_phi_tmp = []
for val in phi_bins:
    binned_phi_tmp.append(((float(val)/180)-1)*math.pi)
phiBin_line = TLine()
for i,b in enumerate(binned_phi_tmp):
    b = float(b)
    phiBin_line.SetLineColor(4)
    phiBin_line.SetLineWidth(4)
    phiBin_line.DrawLine(b,0,b,binmax)
    l_t.AddEntry(phiBin_line,"Bin Edge %s" % i )
    l_t.AddEntry(phiBin_line,"BinCenter = %.2f" % b)
    
Cph_q.Print(outputpdf)

Cth_q = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_th_q_DATA"].SetLineColor(i+1)
    hist["H_th_q_DATA"].Draw("same, E1")
    hist["H_th_q_SIMC"].SetLineColor(40)
    hist["H_th_q_SIMC"].SetLineStyle(10-i)
    hist["H_th_q_SIMC"].Draw("same, E1")
    
Cth_q.Print(outputpdf)

Cph_recoil = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_ph_recoil_DATA"].SetLineColor(i+1)
    hist["H_ph_recoil_DATA"].Draw("same, E1")
    
Cph_recoil.Print(outputpdf)

Cth_recoil = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_th_recoil_DATA"].SetLineColor(i+1)
    hist["H_th_recoil_DATA"].Draw("same, E1")

Cth_recoil.Print(outputpdf)

Cpmiss = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmiss_DATA"].SetLineColor(i+1)
    hist["H_pmiss_DATA"].Draw("same, E1")
    hist["H_pmiss_SIMC"].SetLineColor(40)
    hist["H_pmiss_SIMC"].SetLineStyle(10-i)
    hist["H_pmiss_SIMC"].Draw("same, E1")
    
Cpmiss.Print(outputpdf)

Cemiss = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_emiss_DATA"].SetLineColor(i+1)
    hist["H_emiss_DATA"].Draw("same, E1")
    hist["H_emiss_SIMC"].SetLineColor(40)
    hist["H_emiss_SIMC"].SetLineStyle(10-i)
    hist["H_emiss_SIMC"].Draw("same, E1")
    
Cemiss.Print(outputpdf)

Cpmiss_x = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmx_DATA"].SetLineColor(i+1)
    hist["H_pmx_DATA"].Draw("same, E1")
    hist["H_pmx_SIMC"].SetLineColor(40)
    hist["H_pmx_SIMC"].SetLineStyle(10-i)
    hist["H_pmx_SIMC"].Draw("same, E1")
    
Cpmiss_x.Print(outputpdf)

Cpmiss_y = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmy_DATA"].SetLineColor(i+1)
    hist["H_pmy_DATA"].Draw("same, E1")
    hist["H_pmy_SIMC"].SetLineColor(40)
    hist["H_pmy_SIMC"].SetLineStyle(10-i)
    hist["H_pmy_SIMC"].Draw("same, E1")
    
Cpmiss_y.Print(outputpdf)

Cpmiss_z = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmz_DATA"].SetLineColor(i+1)
    hist["H_pmz_DATA"].Draw("same, E1")
    hist["H_pmz_SIMC"].SetLineColor(40)
    hist["H_pmz_SIMC"].SetLineStyle(10-i)
    hist["H_pmz_SIMC"].Draw("same, E1")
    
Cpmiss_z.Print(outputpdf)

Cmmct = TCanvas()

Cmmct.Divide(2,2)

for i,hist in enumerate(histlist):
    Cmmct.cd(i+1)
    
    # Set color representing zero events to transparent (alpha = 0)
    hist["MM_vs_CoinTime_DATA"].SetMinimum(1)

    hist["MM_vs_CoinTime_DATA"].SetLineColor(i+1)
    hist["MM_vs_CoinTime_DATA"].Draw("same, COLZ")
    hist["MM_vs_CoinTime_DATA"].SetTitle(phisetlist[i])

Cmmct.Print(outputpdf)

Cctbeta = TCanvas()

Cctbeta.Divide(2,2)

for i,hist in enumerate(histlist):
    Cctbeta.cd(i+1)

    # Set color representing zero events to transparent (alpha = 0)
    hist["CoinTime_vs_beta_DATA"].SetMinimum(1)
    
    hist["CoinTime_vs_beta_DATA"].SetLineColor(i+1)
    hist["CoinTime_vs_beta_DATA"].Draw("same, COLZ")
    hist["CoinTime_vs_beta_DATA"].SetTitle(phisetlist[i])

Cctbeta.Print(outputpdf)

Cmmbeta = TCanvas()

Cmmbeta.Divide(2,2)

for i,hist in enumerate(histlist):
    Cmmbeta.cd(i+1)

    # Set color representing zero events to transparent (alpha = 0)
    hist["MM_vs_beta_DATA"].SetMinimum(1)
    
    hist["MM_vs_beta_DATA"].SetLineColor(i+1)
    hist["MM_vs_beta_DATA"].Draw("same, COLZ")
    hist["MM_vs_beta_DATA"].SetTitle(phisetlist[i])

Cmmbeta.Print(outputpdf)

Cqw = TCanvas()

Cqw.Divide(2,2)

for i,hist in enumerate(histlist):
    Cqw.cd(i+1)

    # Set color representing zero events to transparent (alpha = 0)
    hist["Q2_vs_W_DATA"].SetMinimum(1)
    
    hist["Q2_vs_W_DATA"].SetLineColor(i+1)
    hist["Q2_vs_W_DATA"].Draw("same, COLZ")
    hist["Q2_vs_W_DATA"].SetTitle(phisetlist[i])

Cqw.Print(outputpdf)

Cpht = TCanvas()

# Create a list to store all polar plots
polar_plots = []

for i, hist in enumerate(histlist):
    polar_plot = TGraphPolar(hist["polar_phiq_vs_t_DATA"].GetN(), hist["polar_phiq_vs_t_DATA"].GetX(), hist["polar_phiq_vs_t_DATA"].GetY())
    polar_plot.SetMarkerColor(i+1)
    polar_plot.SetMarkerSize(0.5)
    polar_plot.SetMarkerStyle(20)
    polar_plots.append(polar_plot)  # Store the plot in the list
    polar_plot.Draw("AP same")

# Set titles and axes for the last plot
polar_plots[-1].GetXaxis().SetName("#Phi")
polar_plots[-1].GetYaxis().SetName("-t")

Cpht.Print(outputpdf)

for i,hist in enumerate(histlist):
    texlist = []
    Ctext = TCanvas()
    for j,line in enumerate(hist["pid_text"]):
        if j == 0:
            tex = TLatex(0.8,0.+(0.95-(0.3)),"{}".format(hist["phi_setting"]))
            tex.SetTextSize(0.03)
            tex.SetTextColor(i+1)
            texlist.append(tex)
        tex = TLatex(0.,0.+(0.95-(0.3+(0.05*j/2))),"{}".format(line))
        tex.SetTextSize(0.03)
        tex.SetTextColor(i+1)
        texlist.append(tex)

    for j, tex in enumerate(texlist):
        tex.Draw()
        
    if i == len(histlist)-1:
        Ctext.Print(outputpdf+')')
    else:
        Ctext.Print(outputpdf)

if DEBUG:
    show_pdf_with_evince(outputpdf)

##############################
# Step 6 of the lt_analysis: #
##############################
'''
* Combine settings again at for each 
  Q2,W,tbin,phibin,eps for data, dummy, and SIMC.

* Dummy is subtracted from data bin by bin.
* The yield is calculated using the effective charge from data and 
  normfactor/nevents is applied to normalize simc.

* The data and SIMC yields are compared and the R value per bin is obtained.
'''

from calculate_yield import find_yield_data

yieldDict = {}
yieldDict.update(yield_per_bin_data(histlist, inpDict))
print(yieldDict)

from ave_per_bin import ave_per_bin_data, ave_per_bin_simc

aveDict = {}
aveDict.update(ave_per_bin_data(histlist, inpDict))
aveDict.update(ave_per_bin_simc(histlist, inpDict))
#print(aveDict)

# t/phi binned histograms
H_phibins_DATA = ROOT.TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
H_tbins_DATA = ROOT.TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)

histbinDict = {}
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
            H_tbins_DATA.Fill(float(val))
    H_tbins_DATA.SetLineColor(it+1)
    H_tbins_DATA.Draw("same, hist")
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
            H_phibins_DATA.Fill(float(val))
    H_phibins_DATA.SetLineColor(it+1)            
    H_phibins_DATA.Draw("same, hist")
C_phi_bins_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType))+')')

if DEBUG:
    show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

# Run fortran script
from physics_lists import create_lists
create_lists(aveDict, inpDict)

'''
# Yield binned histograms
H_yield_DATA = ROOT.TH1D("H_yield_DATA", "Data Yield", NumtBins*NumPhiBins, 0, 1.0)
H_yield_SIMC = ROOT.TH1D("H_yield_SIMC", "Simc Yield", NumtBins*NumPhiBins, 0, 1.0)
H_ratio = ROOT.TH1D("H_ratio", "Ratio", NumtBins*NumPhiBins, 0, 1.0)

## !!!! Add yield vs phi variable

C_yield_DATA = TCanvas()
# Loop over each tuple key in the dictionary
for i, data_key_tuple in enumerate(aveDict["binned_DATA"]):
    # Access the nested dictionary using the tuple key
    data_nested_dict = aveDict["binned_DATA"][data_key_tuple]
    for hist in histlist:
        # Fill histogram
        H_yield_DATA.Fill(data_nested_dict["yield_data_{}".format(hist["phi_setting"])])
    H_yield_DATA.Draw("same")
    H_yield_DATA.SetLineColor(i+1)
C_yield_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType))+'(')

C_yield_SIMC = TCanvas()
# Loop over each tuple key in the dictionary
for i, simc_key_tuple in enumerate(aveDict["binned_SIMC"]):
    # Access the nested dictionary using the tuple key
    simc_nested_dict = aveDict["binned_SIMC"][simc_key_tuple]
    for hist in histlist:
        # Fill histogram
        H_yield_SIMC.Fill(simc_nested_dict["yield_simc_{}".format(hist["phi_setting"])])
    H_yield_SIMC.Draw("same")
    H_yield_SIMC.SetLineColor(i+1)
C_yield_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

C_ratio = TCanvas()
# Loop over each tuple key in the dictionary
for i, data_key_tuple in enumerate(aveDict["binned_DATA"]):
    # Access the nested dictionary using the tuple key
    data_nested_dict = aveDict["binned_DATA"][data_key_tuple]
    for hist in histlist:
        # Fill histogram
        H_ratio.Fill(data_nested_dict["ratio_{}".format(hist["phi_setting"])])
    H_ratio.Draw("same")
    H_ratio.SetLineColor(i+1)
C_ratio.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

C_yield_data_plt = TCanvas()
G_yield_data_plt = ROOT.TMultiGraph()
l_yield_data_plt = ROOT.TLegend(0.115,0.35,0.33,0.5)

C_yield_data_plt.SetGrid()

yield_data = np.array([])
yield_simc = np.array([])
setting = np.array([])
for hist in histlist:
    # Loop over each tuple key in the dictionary
    for i, data_key_tuple in enumerate(aveDict["binned_DATA"]):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][data_key_tuple]
        simc_nested_dict = aveDict["binned_SIMC"][simc_key_tuple]
        yield_data = np.append(yield_data, [data_nested_dict["yield_data_{}".format(hist["phi_setting"])]])
        yield_simc = np.append(yield_simc, [simc_nested_dict["yield_simc_{}".format(hist["phi_setting"])]])        
        if hist["phi_setting"] == "Center": setting = np.append(setting,0)
        elif hist["phi_setting"] == "Left": setting = np.append(setting,1)
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

i=0
for i,hist in enumerate(histlist):
    while i <= G_yield_data_plt.GetXaxis().GetXmax():
        bin_ix = G_yield_data_plt.GetXaxis().FindBin(i)
        if i == 0: 
            G_yield_data_plt.GetXaxis().SetBinLabel(bin_ix,"Center")
        elif i == 1:
            G_yield_data_plt.GetXaxis().SetBinLabel(bin_ix,"Left")
        else:
            G_yield_data_plt.GetXaxis().SetBinLabel(bin_ix,"Right")
        i+=1

G_yield_data_plt.GetYaxis().SetTitleOffset(1.5)
G_yield_data_plt.GetXaxis().SetTitleOffset(1.5)
G_yield_data_plt.GetXaxis().SetLabelSize(0.04)

l_yield_data_plt.AddEntry(G_yield_data,"Data")
l_yield_data_plt.AddEntry(G_yield_simc,"Simc")
l_yield_data_plt.Draw()

C_yield_data_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

C_ratio_plt = TCanvas()
G_ratio_plt = ROOT.TMultiGraph()

C_ratio_plt.SetGrid()

ratio_data = np.array([])
setting = np.array([])
for hist in histlist:
    # Loop over each tuple key in the dictionary
    for i, data_key_tuple in enumerate(aveDict["binned_DATA"]):
        # Access the nested dictionary using the tuple key
        data_nested_dict = aveDict["binned_DATA"][data_key_tuple]    
        ratio_data = np.append(ratio_data, [data_nested_dict["ratio_{}".format(hist["phi_setting"])]])
        if hist["phi_setting"] == "Center": setting = np.append(setting,0)
        elif hist["phi_setting"] == "Left": setting = np.append(setting,1)
        else: setting = np.append(setting,2)

G_ratio = ROOT.TGraphErrors(len(ratio_data),setting,ratio_data,np.array([0]*len(setting)),np.array([0]*len(ratio_data)))

G_ratio.SetMarkerStyle(21)
G_ratio.SetMarkerSize(1)
G_ratio.SetMarkerColor(1)
G_ratio_plt.Add(G_ratio)

G_ratio_plt.Draw("AP")

G_ratio_plt.SetTitle(" ;Setting; Ratio")

i=0
for i,hist in enumerate(histlist):
    while i <= G_ratio_plt.GetXaxis().GetXmax():
        bin_ix = G_ratio_plt.GetXaxis().FindBin(i)
        if i == 0: 
            G_ratio_plt.GetXaxis().SetBinLabel(bin_ix,"Center")
        elif i == 1:
            G_ratio_plt.GetXaxis().SetBinLabel(bin_ix,"Left")
        else:
            G_ratio_plt.GetXaxis().SetBinLabel(bin_ix,"Right")
        i+=1

G_ratio_plt.GetYaxis().SetTitleOffset(1.5)
G_ratio_plt.GetXaxis().SetTitleOffset(1.5)
G_ratio_plt.GetXaxis().SetLabelSize(0.04)

C_ratio_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType))+')')


#if DEBUG:
show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))
'''

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data and SIMC.

* Find the mean data values of W,Q2,theta,eps for each t bin of high and low epsilon
  for both data and SIMC.
'''
