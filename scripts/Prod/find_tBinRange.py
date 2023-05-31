#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-05-30 21:13:59 trottar"
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
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
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
'''
Importing dummy and random subtraction function
'''

from subtraction import defineHists

################################################################################################################################################

def bin_data(histlist):

    ################################################################################################################################################
    # Define root file trees of interest

    H_t_Right = []
    H_t_Left = []
    H_t_Center = []

    H_phi_Right = []
    H_phi_Left = []
    H_phi_Center = []
    
    for i,hist in enumerate(histlist):
        if hist["phi_setting"] == 'Right':
            InFile_RIGHT_DATA = hist["InFile_DATA"]
            TBRANCH_RIGHT_DATA  = InFile_RIGHT_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
            print("\nCreating right t-bin histogram...")
            # Grab t bin range
            H_list_Right = [(-evt.MandelT,(evt.ph_q+math.pi)*(180/math.pi)) for i,evt in enumerate(TBRANCH_RIGHT_DATA) if (tmin <= -evt.MandelT <= tmax)]
            H_t_Right = [t[0] for t in H_list_Right]
            H_phi_Right = [t[1] for t in H_list_Right]

        if hist["phi_setting"] == 'Left':
            InFile_LEFT_DATA = hist["InFile_DATA"]
            TBRANCH_LEFT_DATA  = InFile_LEFT_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
            print("\nCreating left t-bin histogram...")
            # Grab t bin range
            H_list_Left = [(-evt.MandelT,(evt.ph_q+math.pi)*(180/math.pi)) for i,evt in enumerate(TBRANCH_LEFT_DATA) if (tmin <= -evt.MandelT <= tmax)]
            H_t_Left = [t[0] for t in H_list_Left]
            H_phi_Left = [t[1] for t in H_list_Left]
            
        if hist["phi_setting"] == 'Center':
            InFile_CENTER_DATA = hist["InFile_DATA"]
            TBRANCH_CENTER_DATA  = InFile_CENTER_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))
            print("\nCreating center t-bin histogram...")
            # Grab t bin range
            H_list_Center = [(-evt.MandelT,(evt.ph_q+math.pi)*(180/math.pi)) for i,evt in enumerate(TBRANCH_CENTER_DATA) if (tmin <= -evt.MandelT <= tmax)]
            H_t_Center = [t[0] for t in H_list_Center]
            H_phi_Center = [t[1] for t in H_list_Center]
            
    ################################################################################################################################################

    H_t_BinTest = []
    H_phi_BinTest = []
    for val in settingList:
        if val == "Right":
            H_t_BinTest = np.concatenate((H_t_BinTest, H_t_Right))
            H_phi_BinTest = np.concatenate((H_phi_BinTest, H_phi_Right))
        if val == "Left":
            H_t_BinTest = np.concatenate((H_t_BinTest, H_t_Left))
            H_phi_BinTest = np.concatenate((H_phi_BinTest, H_phi_Left))
        if val == "Center":
            H_t_BinTest = np.concatenate((H_t_BinTest, H_t_Center))
            H_phi_BinTest = np.concatenate((H_phi_BinTest, H_phi_Center))
            
    return [find_phibins(H_phi_BinTest), find_tbins(H_t_BinTest)]


def find_phibins(H_phi_BinTest):

    print("\nFinding phi bins...")
    phi_arr = np.linspace(0.0, 360.0, NumPhiBins+1)

    n, bins, patches = plt.hist(H_phi_BinTest, phi_arr)

    return [n,bins]

def find_tbins(H_t_BinTest):

                
    ################################################################################################################################################

    def histedges_equalN(x, nbin):
        # Grab number of events in array
        npt = len(x)
        # One-dimensional linear interpolation for monotonically increasing sample points.
        # Returns the one-dimensional piecewise linear interpolant to a function with given
        # discrete data points (xp, fp), evaluated at x.
        #
        # np.interp(x, xp, fp)
        # x -> np.linspace(0, npt, nbin + 1) : The x-coordinates at which to evaluate the interpolated values
        # In this case, this is an array of evenly spaced t-bins
        # xp -> np.arange(npt) : The x-coordinates of the data points
        # In this case, this returns evenly spaced values within a given interval
        # yp -> np.sort(x) : he y-coordinates of the data points
        # In this case, this returns a sorted copy of the array
        return np.interp(np.linspace(0, npt, nbin + 1),np.arange(npt),np.sort(x))

    print("\nFinding t bins...")
    # Histogram takes the array data set and the bins as input
    # The bins are determined by a linear interpolation (see function above)
    # This returns the binned data with equal number of events per bin
    # Returns...
    # n -> The values of the histogram bins
    # bins -> The edges of the bins
    # patches -> Container of individual artists used to create the histogram or list of
    # such containers if there are multiple input datasets.
    n, bins, patches = plt.hist(H_t_BinTest, histedges_equalN(H_t_BinTest, NumtBins))
    
    # Write t_bin_interval for lt_analysis scripts
    lines = []
    #with open("{}/src/t_bin_interval_{}_{:.0f}".format(LTANAPATH,Q2.replace("p",""),float(EPSVAL)*100), "w") as file:
    with open("{}/src/t_bin_interval".format(LTANAPATH), "w") as file:
        file.write("{}\t{}\t{}\n".format(Q2.replace("p","."),NumtBins,NumPhiBins))
        for i,t in enumerate(bins):
            lines.append("\t{:.2f}".format(float(t)))
        file.writelines(lines)
    
    return [n,bins]
    
################################################################################################################################################

# Call histogram function above to define dictonaries for right, left, center settings
# Put these all into an array so that if we are missing a setting it is easier to remove
# Plus it makes the code below less repetitive
phisetlist = ["Center","Left","Right"]
histlist = []
for phiset in phisetlist:
    histlist.append(defineHists(phiset,inpDict))

print("$$$$$$$$$$$$$$$$$$$",type(histlist[0]))
    
print("\n\n")

settingList = []
for i,hist in enumerate(histlist):    
    if not bool(hist): # If hist is empty
        histlist.remove(hist)
    else:
        settingList.append(hist["phi_setting"])

for i,hist in enumerate(histlist):
    print("!!!!!!!!!!!!!!!!!!!!!!!!! MM", type(hist["H_MM_DATA"]), hist["H_MM_DATA"])
    print("!!!!!!!!!!!!!!!!!!!!!!!!! tbins", type(hist["H_tbins_DATA"]), hist["H_tbins_DATA"])
        
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

if ParticleType == "kaon":
    eff_plt.Print(outputpdf)
else:
    eff_plt.Print(outputpdf + '(')

c_bins = TCanvas()

binned_data = bin_data(histlist)

binned_phi = binned_data[0]
# binned_phi[0] is missing a value for the final bin
# so adding the first element allows the zip to include all bins
# this is okay because the number of events per bin should be the same
phibinvals = list(binned_phi[0])
phibinvals.append(binned_phi[0][0])

binned_t = binned_data[1]
# binned_t[0] is missing a value for the final bin
# so adding the first element allows the zip to include all bins
# this is okay because the number of events per bin should be the same
tbinvals = list(binned_t[0])
tbinvals.append(binned_t[0][0])

tbinedges = binned_t[1]
phibinedges = binned_phi[1]

for i,hist in enumerate(histlist):
    for j in range(NumtBins):
        for k in range(NumPhiBins):
            hist["H_tbins_DATA"].Fill((tbinedges[j]+tbinedges[j+1])/2)
            hist["H_phibins_DATA"].Fill((phibinedges[k]+phibinedges[k+1])/2)

c_bins.Divide(2,1)
        
for i,hist in enumerate(histlist):
    c_bins.cd(1)
    hist["H_phibins_DATA"].SetLineColor(i+1)
    hist["H_phibins_DATA"].Draw("same")
    
for i,hist in enumerate(histlist):
    c_bins.cd(2)
    hist["H_tbins_DATA"].SetLineColor(i+1)
    hist["H_tbins_DATA"].Draw("same")
    
c_bins.Print(outputpdf)
        
c_yield_data = TCanvas()
        
for i,hist in enumerate(histlist):

    InFile_DATA = hist["InFile_DATA"]
    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    mm_list = []
    aver_lst = []
    for evt in TBRANCH_DATA:
        for j in range(len(tbinedges) - 1):
            if tbinedges[j] <= -evt.MandelT < tbinedges[j+1]:
                tbin_index = j
            else:
                tbin_index = None
            if tbin_index != None:
                aver_lst.append((tbin_index, evt.Q2, evt.W, -evt.MandelT))
                for k in range(len(phibinedges) - 1):                    
                    if phibinedges[k] <= (evt.ph_q+math.pi)*(180/math.pi) < phibinedges[k+1]:
                        phibin_index = k
                    else:
                        phibin_index = None
                    if phibin_index != None:
                        mm_list.append((tbin_index, phibin_index, np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))))

    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in aver_lst:
        key = (t[0])
        if key in groups:
            groups[key].append((t[1], t[2], t[3]))
        else:
            groups[key] = [(t[1], t[2], t[3])]

    # Extract the desired values from each group
    Q2_aver = []
    W_aver = []
    t_aver = []
    for key, val in groups.items():
        Q2_tmp = []
        W_tmp = []
        t_tmp = []
        for tup in val:
            Q2_tmp.append(tup[0])
            W_tmp.append(tup[1])
            t_tmp.append(tup[2])
        Q2_aver.append((key, np.average(Q2_tmp)))
        W_aver.append((key, np.average(W_tmp)))
        t_aver.append((key, np.average(t_tmp)))

    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in mm_list:
        for j,a in enumerate(Q2_aver):
            if a[0] == t[0]:
                key = (t[0], t[1])
                if key in groups:
                    groups[key].append((t[2], Q2_aver[j][1], W_aver[j][1], t_aver[j][1]))
                else:
                    groups[key] = [(t[2], Q2_aver[j][1], W_aver[j][1], t_aver[j][1])]                    

    yieldValData = array('d', [0])
    Q2binValData = array('d', [0])
    WbinValData = array('d', [0])
    tbinValData = array('d', [0])
    
    tnum = array('d', [0])
    phinum = array('d', [0])
    tval = array('d', [0])
    phival = array('d', [0])
    
    hist["yieldTree"].Branch("yield_data", yieldValData, "yield_data/D")
    hist["yieldTree"].Branch("aver_Q2", Q2binValData, "aver_Q2/D")
    hist["yieldTree"].Branch("aver_W", WbinValData, "aver_W/D")
    hist["yieldTree"].Branch("aver_t", tbinValData, "aver_t/D")    
    hist["yieldTree"].Branch("tbins", tnum, "tbins/D")
    hist["yieldTree"].Branch("phibins", phinum, "phibins/D")
    hist["yieldTree"].Branch("tbincenter", tval, "tbincenter/D")
    hist["yieldTree"].Branch("phibincenter", phival, "phibincenter/D")

    tbinarr = []
    phibinarr = []
    # Extract the desired values from each group
    for key, val in groups.items():
        j = key[0]
        k = key[1]
        tbinarr.append(j)
        phibinarr.append(k)
        tnum[0] = j+1
        phinum[0] = k+1
        tval[0] = (tbinedges[j]+tbinedges[j+1])/2
        phival[0] = (phibinedges[k]+phibinedges[k+1])/2

        MM_tmp = []
        Q2_tmp = []
        W_tmp = []
        t_tmp = []
        for tup in val:
            MM_tmp.append(tup[0])
            Q2_tmp.append(tup[1])
            W_tmp.append(tup[2])
            t_tmp.append(tup[3])
        hist["H_yield_DATA"].Fill(integrate.simps(MM_tmp)*hist["normfac_data"])
        hist["yieldDictData"][key] = integrate.simps(MM_tmp)*hist["normfac_data"]
        yieldValData[0] = integrate.simps(MM_tmp)*hist["normfac_data"]
        Q2binValData[0] = Q2_tmp[0]
        WbinValData[0] = W_tmp[0]
        tbinValData[0] = t_tmp[0]
        hist["yieldTree"].Fill()

    hist["yieldTree"].ResetBranchAddresses()
    
    print("\n\n~~~~~~~~~~~~~~~",hist["yieldDictData"])
    print("~~~~~~~~~~~~~~~",hist["H_yield_DATA"])
    hist["H_yield_DATA"].SetLineColor(i+1)            
    hist["H_yield_DATA"].Draw("same")
        
c_yield_data.Print(outputpdf)

c_yield_simc = TCanvas()

for i,hist in enumerate(histlist):

    InFile_SIMC = hist["InFile_SIMC"]
    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    tmp_lst = []
    for evt in TBRANCH_SIMC:
        for j in range(len(tbinedges) - 1):
            if tbinedges[j] <= evt.t < tbinedges[j+1]:
                tbin_index = j
            else:
                tbin_index = None
            if tbin_index != None:
                for k in range(len(phibinedges) - 1):
                    if phibinedges[k] <= (evt.phipq)*(180/math.pi) < phibinedges[k+1]:
                        phibin_index = k
                    else:
                        phibin_index = None
                    if phibin_index != None:
                        tmp_lst.append((tbin_index, phibin_index, np.sqrt(pow(evt.Em, 2) - pow(evt.Pm, 2))*evt.Weight, evt.Q2, evt.W, evt.t))
            
    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in tmp_lst:
        for j,k in zip(tbinarr,phibinarr):
            if t[0] == j and t[1] == k:
                key = (t[0], t[1])
                if key in groups:
                    groups[key].append((t[2], t[3], t[4], t[5]))
                else:
                    groups[key] = [(t[2], t[3], t[4], t[5])]
            else:
                continue
            
    yieldValSimc = array('d', [0])
    
    hist["yieldTree"].Branch("yield_simc", yieldValSimc, "yield_simc/D")
    
    # Extract the desired values from each group
    for key, val in groups.items():
        MM_tmp = []
        Q2_tmp = []
        W_tmp = []
        t_tmp = []
        for tup in val:
            MM_tmp.append(tup[0])
            Q2_tmp.append(tup[1])
            W_tmp.append(tup[2])
            t_tmp.append(tup[3])
        hist["H_yield_SIMC"].Fill(integrate.simps(MM_tmp)*hist["normfac_simc"])
        hist["yieldDictSimc"][key] = integrate.simps(MM_tmp)*hist["normfac_simc"]
        yieldValSimc[0] = integrate.simps(MM_tmp)*hist["normfac_simc"]
        hist["yieldTree"].Fill()

    hist["yieldTree"].ResetBranchAddresses()
            
    print("\n\n~~~~~~~~~~~~~~~",hist["yieldDictSimc"])
    print("~~~~~~~~~~~~~~~",hist["H_yield_SIMC"])
    hist["H_yield_SIMC"].SetLineColor(i+1)            
    hist["H_yield_SIMC"].Draw("same")
        
c_yield_simc.Print(outputpdf)

c_Q2tbin = TCanvas()

c_Q2tbin.Divide(3, int(NumtBins/2))

for i,hist in enumerate(histlist):

    InFile_DATA = hist["InFile_DATA"]
    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    aver_lst = []
    for evt in TBRANCH_DATA:
        for j in range(len(tbinedges) - 1):
            if tbinedges[j] <= -evt.MandelT < tbinedges[j+1]:
                tbin_index = j
            else:
                tbin_index = None
            if tbin_index != None:
                aver_lst.append((tbin_index, evt.Q2))

    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in aver_lst:
        key = (t[0])
        if key in groups:
            groups[key].append((t[1]))
        else:
            groups[key] = [(t[1])]

    # Extract the desired values from each group
    for key, val in groups.items():
        for tup in val:
            hist["H_Q2_tbin_DATA_{}".format(key+1)].Fill(tup)
        c_Q2tbin.cd(key+1)
        hist["H_Q2_tbin_DATA_{}".format(key+1)].Draw("same")
        hist["H_Q2_tbin_DATA_{}".format(key+1)].SetLineColor(i+1)

c_Q2tbin.Print(outputpdf)

c_Wtbin = TCanvas()

c_Wtbin.Divide(3, int(NumtBins/2))

for i,hist in enumerate(histlist):

    InFile_DATA = hist["InFile_DATA"]
    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    aver_lst = []
    for evt in TBRANCH_DATA:
        for j in range(len(tbinedges) - 1):
            if tbinedges[j] <= -evt.MandelT < tbinedges[j+1]:
                tbin_index = j
            else:
                tbin_index = None
            if tbin_index != None:
                aver_lst.append((tbin_index, evt.W))

    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in aver_lst:
        key = (t[0])
        if key in groups:
            groups[key].append((t[1]))
        else:
            groups[key] = [(t[1])]

    # Extract the desired values from each group
    for key, val in groups.items():
        for tup in val:
            hist["H_W_tbin_DATA_{}".format(key+1)].Fill(tup)
        c_Wtbin.cd(key+1)
        hist["H_W_tbin_DATA_{}".format(key+1)].Draw("same")
        hist["H_W_tbin_DATA_{}".format(key+1)].SetLineColor(i+1)

c_Wtbin.Print(outputpdf)

c_ttbin = TCanvas()

c_ttbin.Divide(3, int(NumtBins/2))

for i,hist in enumerate(histlist):

    InFile_DATA = hist["InFile_DATA"]
    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    aver_lst = []
    for evt in TBRANCH_DATA:
        for j in range(len(tbinedges) - 1):
            if tbinedges[j] <= -evt.MandelT < tbinedges[j+1]:
                tbin_index = j
            else:
                tbin_index = None
            if tbin_index != None:
                aver_lst.append((tbin_index, -evt.MandelT))

    groups = {}
    # Group the tuples by the first two elements using a dictionary
    for t in aver_lst:
        key = (t[0])
        if key in groups:
            groups[key].append((t[1]))
        else:
            groups[key] = [(t[1])]

    # Extract the desired values from each group
    for key, val in groups.items():
        for tup in val:
            hist["H_t_tbin_DATA_{}".format(key+1)].Fill(tup)
        c_ttbin.cd(key+1)
        hist["H_t_tbin_DATA_{}".format(key+1)].Draw("same")
        hist["H_t_tbin_DATA_{}".format(key+1)].SetLineColor(i+1)

c_ttbin.Print(outputpdf)

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
for i,(n,b) in enumerate(zip(tbinvals,tbinedges)):
    tBin_line.SetLineColor(4)
    tBin_line.SetLineWidth(4)
    tBin_line.DrawLine(b,0,b,binmax)
    l_t.AddEntry(tBin_line,"Bin Edge %s" % i )
    l_t.AddEntry(tBin_line,"Evts = %.0f" % n)
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
for val in binned_phi[1]:
    binned_phi_tmp.append(((val/180)-1)*math.pi)
phiBin_line = TLine()
for i,(n,b) in enumerate(zip(phibinvals,binned_phi_tmp)):
    phiBin_line.SetLineColor(4)
    phiBin_line.SetLineWidth(4)
    phiBin_line.DrawLine(b,0,b,binmax)
    l_t.AddEntry(phiBin_line,"Bin Edge %s" % i )
    l_t.AddEntry(phiBin_line,"Evts = %.0f" % n)
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
    #hist["H_pmx_SIMC"].SetLineColor(40)
    #hist["H_pmx_SIMC"].SetLineStyle(10-i)
    #hist["H_pmx_SIMC"].Draw("same, E1")
    
Cpmiss_x.Print(outputpdf)

Cpmiss_y = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmy_DATA"].SetLineColor(i+1)
    hist["H_pmy_DATA"].Draw("same, E1")
    #hist["H_pmy_SIMC"].SetLineColor(40)
    #hist["H_pmy_SIMC"].SetLineStyle(10-i)
    #hist["H_pmy_SIMC"].Draw("same, E1")
    
Cpmiss_y.Print(outputpdf)

Cpmiss_z = TCanvas()

for i,hist in enumerate(histlist):
    hist["H_pmz_DATA"].SetLineColor(i+1)
    hist["H_pmz_DATA"].Draw("same, E1")
    #hist["H_pmz_SIMC"].SetLineColor(40)
    #hist["H_pmz_SIMC"].SetLineStyle(10-i)
    #hist["H_pmz_SIMC"].Draw("same, E1")
    
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
for i,(n,b) in enumerate(zip(tbinvals,tbinedges)):
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
        
#############################################################################################################################################
# Create new root file with trees representing cut simc and data used above. Good for those who see python as...problematic

outHistFile = ROOT.TFile.Open(foutname, "RECREATE")

for i,hist in enumerate(histlist):
    
    hist["yieldTree"].Write()
    
    if hist["phi_setting"] == "Right":
        d_Right_Data = outHistFile.mkdir("Right Data")
        d_Right_Simc = outHistFile.mkdir("Right Simc")
    if hist["phi_setting"] == "Left":
        d_Left_Data = outHistFile.mkdir("Left Data")
        d_Left_Simc = outHistFile.mkdir("Left Simc")
    if hist["phi_setting"] == "Center":
        d_Center_Data = outHistFile.mkdir("Center Data")
        d_Center_Simc = outHistFile.mkdir("Center Simc")
    
for i,hist in enumerate(histlist):
    if hist["phi_setting"] == "Right":
        d_Right_Data.cd()
    elif hist["phi_setting"] == "Left":
        d_Left_Data.cd()
    elif hist["phi_setting"] == "Center":
        d_Center_Data.cd()
    else:
        continue
    hist["H_hsdelta_DATA"].Write()
    hist["H_hsxptar_DATA"].Write()
    hist["H_hsyptar_DATA"].Write()
    hist["H_ssxfp_DATA"].Write()
    hist["H_ssyfp_DATA"].Write()
    hist["H_ssxpfp_DATA"].Write()
    hist["H_ssypfp_DATA"].Write()
    hist["H_hsxfp_DATA"].Write()
    hist["H_hsyfp_DATA"].Write()
    hist["H_hsxpfp_DATA"].Write()
    hist["H_hsypfp_DATA"].Write()
    hist["H_ssdelta_DATA"].Write()
    hist["H_ssxptar_DATA"].Write()
    hist["H_ssyptar_DATA"].Write()
    hist["H_q_DATA"].Write()
    hist["H_Q2_DATA"].Write()
    hist["H_W_DATA"].Write()
    hist["H_t_DATA"].Write()
    hist["H_epsilon_DATA"].Write()
    hist["H_MM_DATA"].Write()
    hist["H_th_DATA"].Write()
    hist["H_ph_DATA"].Write()
    hist["H_ph_q_DATA"].Write()
    hist["H_th_q_DATA"].Write()
    hist["H_ph_recoil_DATA"].Write()
    hist["H_th_recoil_DATA"].Write()
    hist["H_pmiss_DATA"].Write()
    hist["H_emiss_DATA"].Write()
    hist["H_pmx_DATA"].Write()
    hist["H_pmy_DATA"].Write()
    hist["H_pmz_DATA"].Write()
    hist["H_ct_DATA"].Write()
    for b in range(NumtBins):
        hist["H_Q2_tbin_DATA_{}".format(b+1)].Write()
        hist["H_W_tbin_DATA_{}".format(b+1)].Write()
        hist["H_t_tbin_DATA_{}".format(b+1)].Write()
    
for i,hist in enumerate(histlist):
    if hist["phi_setting"] == "Right":
        d_Right_Simc.cd()
    elif hist["phi_setting"] == "Left":
        d_Left_Simc.cd()
    elif hist["phi_setting"] == "Center":
        d_Center_Simc.cd()
    else:
        continue        
    hist["H_hsdelta_SIMC"].Write()
    hist["H_hsxptar_SIMC"].Write()
    hist["H_hsyptar_SIMC"].Write()
    hist["H_ssxfp_SIMC"].Write()
    hist["H_ssyfp_SIMC"].Write()
    hist["H_ssxpfp_SIMC"].Write()
    hist["H_ssypfp_SIMC"].Write()
    hist["H_hsxfp_SIMC"].Write()
    hist["H_hsyfp_SIMC"].Write()
    hist["H_hsxpfp_SIMC"].Write()
    hist["H_hsypfp_SIMC"].Write()
    hist["H_ssdelta_SIMC"].Write()
    hist["H_ssxptar_SIMC"].Write()
    hist["H_ssyptar_SIMC"].Write()
    hist["H_q_SIMC"].Write()
    hist["H_Q2_SIMC"].Write()
    hist["H_W_SIMC"].Write()
    hist["H_t_SIMC"].Write()
    hist["H_epsilon_SIMC"].Write()
    hist["H_MM_SIMC"].Write()
    hist["H_th_SIMC"].Write()
    hist["H_ph_SIMC"].Write()
    hist["H_ph_q_SIMC"].Write()
    hist["H_th_q_SIMC"].Write()
    hist["H_ph_recoil_SIMC"].Write()
    hist["H_th_recoil_SIMC"].Write()
    hist["H_pmiss_SIMC"].Write()
    hist["H_emiss_SIMC"].Write()
    hist["H_pmx_SIMC"].Write()
    hist["H_pmy_SIMC"].Write()
    hist["H_pmz_SIMC"].Write()

outHistFile.Close()

for i,hist in enumerate(histlist):
    hist["InFile_DATA"].Close()
    hist["InFile_DUMMY"].Close()
    hist["InFile_SIMC"].Close()
    if ParticleType == "kaon":
        hist["InFile_SUBPION_DATA"].Close()
        hist["InFile_SUBPION_DUMMY"].Close()
        hist["InFile_SUBPROTON_DATA"].Close()
        hist["InFile_SUBPROTON_DUMMY"].Close()

print ("Processing Complete")

