#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-04 14:23:46 trottar"
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

if len(sys.argv)-1!=26:
    print("!!!!! ERROR !!!!!\n Expected 25 arguments\n Usage is with - KIN W Q2 EPSVAL OutDATAFilename OutDUMMYFilename OutFullAnalysisFilename tmin tmax NumtBins NumPhiBins runNumRight runNumLeft runNumCenter data_charge_right data_charge_left data_charge_center dummy_charge_right dummy_charge_left dummy_charge_center InData_efficiency_right InData_efficiency_left InData_efficiency_center efficiency_table ParticleType EPSSET\n!!!!! ERROR !!!!!")
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
EPSSET = sys.argv[26]

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
    "EPSSET" : EPSSET
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

> lt_analysis/OUTPUT/Analysis/<ANATYPE>LT/<KIN>_Center_Diamond_Cut.pdf

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
        show_pdf_with_evince(OUTPATH+"/%s_%s_Diamond_Cut.pdf" %(('Q'+Q2+'W'+W,phiset)))

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
    if not bool(hist): # If hist is empty
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
except FileNotFoundError:
    print("{} not found...".format("{}/src/phi_bin_interval".format(LTANAPATH)))
except IOError:
    print("Error reading {}...".format("{}/src/phi_bin_interval".format(LTANAPATH)))    
    
print(t_bins,"\n",phi_bins)

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

from compare_simc import compare_simc

# Upate hist dictionary with effective charge and simc histograms
for hist in histlist:
    hist.update(compare_simc(hist, inpDict))


### !!!!!!!!!!!! Need to edit t/phi plots, etc. to show bins and labels on plots    
    
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
    hist["MM_vs_CoinTime_DATA"].SetLineColor(i+1)
    hist["MM_vs_CoinTime_DATA"].Draw("same, COLZ")
    hist["MM_vs_CoinTime_DATA"].SetTitle(phisetlist[i])

Cmmct.Print(outputpdf)

Cctbeta = TCanvas()

Cctbeta.Divide(2,2)

for i,hist in enumerate(histlist):
    Cctbeta.cd(i+1)
    hist["CoinTime_vs_beta_DATA"].SetLineColor(i+1)
    hist["CoinTime_vs_beta_DATA"].Draw("same, COLZ")
    hist["CoinTime_vs_beta_DATA"].SetTitle(phisetlist[i])

Cctbeta.Print(outputpdf)

Cmmbeta = TCanvas()

Cmmbeta.Divide(2,2)

for i,hist in enumerate(histlist):
    Cmmbeta.cd(i+1)
    hist["MM_vs_beta_DATA"].SetLineColor(i+1)
    hist["MM_vs_beta_DATA"].Draw("same, COLZ")
    hist["MM_vs_beta_DATA"].SetTitle(phisetlist[i])

Cmmbeta.Print(outputpdf)

Cqw = TCanvas()

Cqw.Divide(2,2)

for i,hist in enumerate(histlist):
    Cqw.cd(i+1)
    hist["Q2_vs_W_DATA"].SetLineColor(i+1)
    hist["Q2_vs_W_DATA"].Draw("same, COLZ")
    hist["Q2_vs_W_DATA"].SetTitle(phisetlist[i])

Cqw.Print(outputpdf)

Cpht = TCanvas()

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Create a new TMultiGraph object
multi_graph = ROOT.TMultiGraph()

# Loop over each TGraphPolar object and add it to the TMultiGraph
for i, hist in enumerate(histlist):
    hist["polar_phiq_vs_t_DATA"].SetMarkerSize(2)
    hist["polar_phiq_vs_t_DATA"].SetMarkerColor(i+1)
    #hist["polar_phiq_vs_t_DATA"].SetMarkerStyle(ROOT.kFullCircle)
    #hist["polar_phiq_vs_t_DATA"].SetLineColor(i+1)
    multi_graph.Add(hist["polar_phiq_vs_t_DATA"])
    Cpht.Update()
    #hist["polar_phiq_vs_t_DATA"].GetPolargram().SetRangeRadial(0, 2.0)

multi_graph.Draw("COLZ")

# Customize the polar surface plot
multi_graph.GetXaxis().SetTitle("p_{T} (GeV/c)")
multi_graph.GetYaxis().SetTitle("#phi")
multi_graph.GetYaxis().SetTitleOffset(1.2)
#multi_graph.GetZaxis().SetTitle("Counts")
#multi_graph.GetZaxis().SetTitleOffset(1.5)

Cpht.Update()
    
'''
Cpht.Divide(2,2)

for i,hist in enumerate(histlist):
    Cpht.cd(i+1)
    hist["phiq_vs_t_DATA"].GetYaxis().SetRangeUser(tmin,tmax)
    hist["phiq_vs_t_DATA"].Draw("SURF2 POL")
    hist["phiq_vs_t_DATA"].SetTitle(phisetlist[i])
    
# Section for polar plotting
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
tvsphi_title = TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC")
tvsphi_title.AddText("-t vs #phi")
tvsphi_title.Draw()
Cpht.Update()
ptphizero = TPaveText(0.923951,0.513932,0.993778,0.574551,"NDC")
ptphizero.AddText("#phi = 0")
ptphizero.Draw()
Cpht.Update()
phihalfk = TLine(0,0,0,0.6)
phihalfk.SetLineColor(kBlack)
phihalfk.SetLineWidth(2)
phihalfk.Draw()
Cpht.Update()
ptphihalfk = TPaveText(0.417855,0.901876,0.486574,0.996358,"NDC")
ptphihalfk.AddText("#phi = #frac{K}{2}")
ptphihalfk.Draw()
Cpht.Update()
phik = TLine(0,0,-0.6,0)
phik.SetLineColor(kBlack)
phik.SetLineWidth(2)
phik.Draw()
Cpht.Update()
ptphik = TPaveText(0.0277092,0.514217,0.096428,0.572746,"NDC")
ptphik.AddText("#phi = K")
ptphik.Draw()
Cpht.Update()
phithreek = TLine(0,0,0,-0.6)
phithreek.SetLineColor(kBlack)
phithreek.SetLineWidth(2)
phithreek.Draw()
Cpht.Update()
ptphithreek = TPaveText(0.419517,0.00514928,0.487128,0.0996315,"NDC")
ptphithreek.AddText("#phi = #frac{3K}{2}")
ptphithreek.Draw()
Cpht.Update()
Arc = TArc()
for k in range(0, 10):
     Arc.SetFillStyle(0)
     Arc.SetLineWidth(2)
     # To change the arc radius we have to change number 0.6 in the lower line.
     Arc.DrawArc(0,0,0.6*(k+1)/(10),0.,360.,"same")
     Cpht.Update()
for i,(n,b) in enumerate(zip(tbinvals,tbinedges)):
     Arc.SetLineColor(3)
     Arc.SetLineWidth(2)
     # To change the arc radius we have to change number 0.6 in the lower line.
     Arc.DrawArc(0,0,0.6*b,0.,360.,"same")
     Cpht.Update()
tradius = TGaxis(0,0,0.6,0,tmin,tmax,10,"-+")
tradius.SetLineColor(2)
tradius.SetLabelColor(2)
tradius.Draw()
Cpht.Update()
'''

Cpht.Print(outputpdf)

Cphtsame = TCanvas()

for i,hist in enumerate(histlist):
    # set colors for the TGraphPolar object
    hist["polar_phiq_vs_t_DATA"].SetMarkerSize(2)
    hist["polar_phiq_vs_t_DATA"].SetMarkerColor(i+1)
    hist["polar_phiq_vs_t_DATA"].SetMarkerStyle(ROOT.kFullCircle)
    hist["polar_phiq_vs_t_DATA"].SetLineColor(i+1)
    hist["polar_phiq_vs_t_DATA"].Draw("AOP")
    Cphtsame.Update()
    hist["polar_phiq_vs_t_DATA"].GetPolargram().SetRangeRadial(0, 2.0)
    # Hide radial axis labels since redefined below
    hist["polar_phiq_vs_t_DATA"].GetPolargram().SetRadialLabelSize(0)
    Cphtsame.Update()

# Section for polar plotting
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
tvsphi_title = TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC")
tvsphi_title.AddText("-t vs #phi")
tvsphi_title.Draw()
phihalfk = TLine(0,0,0,tmax)
phihalfk.SetLineColor(kBlack)
phihalfk.SetLineWidth(2)
phihalfk.Draw()
phik = TLine(0,0,-tmax,0)
phik.SetLineColor(kBlack)
phik.SetLineWidth(2)
phik.Draw()
phithreek = TLine(0,0,0,-tmax)
phithreek.SetLineColor(kBlack)
phithreek.SetLineWidth(2)
phithreek.Draw()
Arc = TArc()
for k in range(0, 10):
    Arc.SetFillStyle(0)
    Arc.SetLineWidth(2)
    # To change the arc radius we have to change number tmax in the lower line.
    Arc.DrawArc(0,0,tmax*(k+1)/(10),0.,360.,"same")
for i,b in enumerate(t_bins):
    b = float(b)
    Arc.SetLineColor(9)
    Arc.SetLineWidth(2)
    # To change the arc radius we have to change number tmax in the lower line.
    Arc.DrawArc(0,0,tmax*b,0.,360.,"same")
tradius = TGaxis(0,0,tmax,0,tmin,tmax,10,"-+")
tradius.SetLineColor(9)
tradius.SetLabelColor(9)
tradius.Draw()
    
Cphtsame.Print(outputpdf)

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

from calculate_yield import calculate_yield

yieldDict = {}

for DataType in ["DATA","DUMMY","SIMC"]:
    yieldDict.update(calculate_yield(histlist, inpDict, DataType))

# Loop over each tuple key in the dictionary
for data_key_tuple,dummy_key_tuple in zip(yieldDict["binned_DATA"],yieldDict["binned_DUMMY"]):
    # Access the nested dictionary using the tuple key
    data_nested_dict = yieldDict["binned_DATA"][data_key_tuple]
    dummy_nested_dict = yieldDict["binned_DUMMY"][dummy_key_tuple]
    #print("\n\nData-> Tuple: {}, Nested Dictionary: {}".format(data_key_tuple,data_nested_dict))
    #print("Dummy-> Tuple: {}, Nested Dictionary: {}".format(dummy_key_tuple,dummy_nested_dict))
    for hist in histlist:
        print("{} Data-> Tuple: {}, Data yield: {}, ".format(hist["phi_setting"],data_key_tuple,data_nested_dict["nevents"]*hist["normfac_data"]))
        print("{} Dummy-> Tuple: {}, Dummy yield: {}, ".format(hist["phi_setting"],dummy_key_tuple,dummy_nested_dict["nevents"]*hist["normfac_dummy"]))
    # Calculate data yield per t/phi bin
    for hist in histlist:
        # Subtract dummy from data per t/phi bin and get data yield
        data_nested_dict["yield_{}".format(hist["phi_setting"])] = data_nested_dict["nevents"]*hist["normfac_data"] - dummy_nested_dict["nevents"]*hist["normfac_dummy"]
        
#print(yieldDict["binned_DATA"])

# t/phi binned histograms
H_phibins_DATA = ROOT.TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
H_tbins_DATA = ROOT.TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)
H_yield_DATA = ROOT.TH1D("H_yield_DATA", "Data Yield", NumtBins*NumPhiBins, 0, 1.0)

histbinDict = {}
# Loop over each tuple key in the dictionary
for data_key_tuple in yieldDict["binned_DATA"]:
    i = data_key_tuple[0] # t bin
    j = data_key_tuple[1] # phi bin
    histbinDict["H_Q2_tbin_DATA_{}_{}".format(i+1,j+1)] = ROOT.TH1D("H_Q2_tbin_DATA_{}_{}".format(i+1,j+1), "Q2 (t bin {}, phi bin {}".format(i+1,j+1), 500, inpDict["Q2min"], inpDict["Q2max"])
    histbinDict["H_W_tbin_DATA_{}_{}".format(i+1,j+1)] = ROOT.TH1D("H_W_tbin_DATA_{}_{}".format(i+1,j+1), "W (t bin {}, phi bin {}".format(i+1,j+1), 500, inpDict["Wmin"], inpDict["Wmax"])
    histbinDict["H_t_tbin_DATA_{}_{}".format(i+1,j+1)] = ROOT.TH1D("H_t_tbin_DATA_{}_{}".format(i+1,j+1), "t (t bin {}, phi bin {}".format(i+1,j+1), 500, inpDict["tmin"], inpDict["tmax"])
    histbinDict["H_MM_tbin_DATA_{}_{}".format(i+1,j+1)] = ROOT.TH1D("H_MM_tbin_DATA_{}_{}".format(i+1,j+1), "MM (t bin {}, phi bin {}".format(i+1,j+1), 500, 0.0, 1.5)   

C_Q2_tbin_DATA = TCanvas()
C_Q2_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for data_key_tuple in yieldDict["binned_DATA"]:
    i = data_key_tuple[0] # t bin
    j = data_key_tuple[1] # phi bin
    # Access the nested dictionary using the tuple key
    data_nested_dict = yieldDict["binned_DATA"][data_key_tuple]    
    # Fill histogram
    for val in data_nested_dict["Q2_arr"]:
        histbinDict["H_Q2_tbin_DATA_{}_{}".format(i+1,j+1)].Fill(val)
    C_Q2_tbin_DATA.cd(i+j+1)
    histbinDict["H_Q2_tbin_DATA_{}_{}".format(i+1,j+1)].Draw("same")
    histbinDict["H_Q2_tbin_DATA_{}_{}".format(i+1,j+1)].SetLineColor(i+1)
C_Q2_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType))+'(')

C_W_tbin_DATA = TCanvas()
C_W_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for data_key_tuple in yieldDict["binned_DATA"]:
    i = data_key_tuple[0] # t bin
    j = data_key_tuple[1] # phi bin
    # Access the nested dictionary using the tuple key
    data_nested_dict = yieldDict["binned_DATA"][data_key_tuple]    
    # Fill histogram
    for val in data_nested_dict["W_arr"]:
        histbinDict["H_W_tbin_DATA_{}_{}".format(i+1,j+1)].Fill(val)
    C_W_tbin_DATA.cd(i+j+1)
    histbinDict["H_W_tbin_DATA_{}_{}".format(i+1,j+1)].Draw("same")
    histbinDict["H_W_tbin_DATA_{}_{}".format(i+1,j+1)].SetLineColor(i+1)
C_W_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

C_t_tbin_DATA = TCanvas()
C_t_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for data_key_tuple in yieldDict["binned_DATA"]:
    i = data_key_tuple[0] # t bin
    j = data_key_tuple[1] # phi bin
    # Access the nested dictionary using the tuple key
    data_nested_dict = yieldDict["binned_DATA"][data_key_tuple]    
    # Fill histogram
    for val in data_nested_dict["t_arr"]:
        histbinDict["H_t_tbin_DATA_{}_{}".format(i+1,j+1)].Fill(val)
    C_t_tbin_DATA.cd(i+j+1)
    histbinDict["H_t_tbin_DATA_{}_{}".format(i+1,j+1)].Draw("same")
    histbinDict["H_t_tbin_DATA_{}_{}".format(i+1,j+1)].SetLineColor(i+1)
C_t_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

C_MM_tbin_DATA = TCanvas()
C_MM_tbin_DATA.Divide(NumtBins,NumPhiBins)
# Loop over each tuple key in the dictionary
for data_key_tuple in yieldDict["binned_DATA"]:
    i = data_key_tuple[0] # t bin
    j = data_key_tuple[1] # phi bin
    # Access the nested dictionary using the tuple key
    data_nested_dict = yieldDict["binned_DATA"][data_key_tuple]    
    # Fill histogram
    for val in data_nested_dict["MM_arr"]:
        histbinDict["H_MM_tbin_DATA_{}_{}".format(i+1,j+1)].Fill(val)
    C_MM_tbin_DATA.cd(i+j+1)
    histbinDict["H_MM_tbin_DATA_{}_{}".format(i+1,j+1)].Draw("same")
    histbinDict["H_MM_tbin_DATA_{}_{}".format(i+1,j+1)].SetLineColor(i+1)
C_MM_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType))+')')

#if DEBUG:
show_pdf_with_evince(outputpdf.replace("{}_".format(ParticleType),"{}_{}_yield_".format(hist["phi_setting"],ParticleType)))

##############################
# Step 7 of the lt_analysis: #
##############################
'''
* Calculate error weighted average of data and SIMC.

* Find the mean data values of W,Q2,theta,eps for each t bin of high and low epsilon
  for both data and SIMC.
'''
