#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-12 16:54:24 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import root_pandas as rpd
import numpy as np
import ROOT
from ROOT import TCanvas, TGraph, TGraphErrors, TMultiGraph, TLegend
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce
import re
import sys, math, os, subprocess

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=10:
    print("!!!!! ERROR !!!!!\n Expected 10 arguments\n Usage is with - ParticleType POL Q2 W LOEPS HIEPS NumtBins NumPhiBins KIN OutUnsepxsectsFilename\n!!!!! ERROR !!!!!")
    sys.exit(1)

###############################################################################################################################################

ParticleType = sys.argv[1]
POL = sys.argv[2]

Q2 = sys.argv[3]
W = sys.argv[4]

LOEPS = sys.argv[5]
HIEPS = sys.argv[6]

NumtBins = int(sys.argv[7])
NumPhiBins = int(sys.argv[8])

kinematics = sys.argv[9]
OutFilename = sys.argv[10]

if int(POL) == 1:
    pol_str = "pl"
elif int(POL) == -1:
    pol_str = "mn"
else:
    print("ERROR: Invalid polarity...must be +1 or -1")
    sys.exit(2)
    
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

foutname = OUTPATH+"/" + OutFilename + ".root"
fouttxt  = OUTPATH+"/" + OutFilename + ".txt"
outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

################################################################################################################################################

def are_within_tolerance(num1, num2, tolerance=0.1):
    return abs(num1 - num2) <= tolerance

def file_to_df(f_name, columns):
    '''
    Read in file and convert to dataframe with custom column names
    '''

    lineskip=False
    
    # Open the file for reading
    with open(f_name, 'r') as file:
        lines = file.readlines()
        if '1.000000\n' in lines:
            lineskip=True

    if lineskip:
        df = pd.read_csv(f_name, header=None, sep=' ', skiprows=1, skipfooter=1, engine='python')
    else:
        df = pd.read_csv(f_name, header=None, sep=' ')    
    df.columns = columns
    return df

################################################################################################################################################

def fix_spacing(f_name):
    '''
    Fortran created files are bad with spacing. This fixes it.
    '''

    # Open the file for reading
    with open(f_name, 'r') as file:

        # Read the lines of the file and split based on whitespace
        lines = file.readlines()
        lines = [re.split(r'\s+', line.strip()) for line in lines]

        # Join the split lines with a single space
        lines = [' '.join(line) for line in lines]

        # Write the lines back to the file
        with open(f_name, 'w') as output:
            output.write('\n'.join(lines))

# Fix file spacing to work in pandas
fix_spacing(LTANAPATH+"/src/{}/averages/avek.{}.dat".format(ParticleType, Q2.replace("p","")))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_sep.{}_{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_sep.{}_{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100))
################################################################################################################################################
# Read in files and convert to dataframes

file_df_dict = {}

setting_file = LTANAPATH+"/src/{}/list.settings".format(ParticleType)
file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins'])

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

t_bin_centers = (t_bins[:-1] + t_bins[1:]) / 2

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

#phi_bin_centers = (phi_bins[:-1] + phi_bins[1:]) / 2
phi_bin_centers = phi_bins
    
for i,row in file_df_dict['setting_df'].iterrows():
    if row['Q2'] == float(Q2.replace("p",".")):
        file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/{}/beam/Eb_KLT.dat".format(ParticleType), ['ebeam', 'Q2', 'EPSVAL'])
        file_df_dict['avek_file'] = file_to_df(LTANAPATH+"/src/{}/averages/avek.{}.dat".format(ParticleType, Q2.replace("p","")) \
                                               , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])        
        
        if row['EPSVAL'] == float(LOEPS):
            file_df_dict['aver_loeps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin'])            
            if row['thpq'] < 0.0:
                file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 'tm', 'W', 'Q2'])

            file_df_dict['sep_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_sep.{}_{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), float(LOEPS)*100) \
                                                            , ['sigL', 'dsigL','sigT', 'dsigT','sigTT', 'dsigTT','sigLT', 'dsigLT', 'Q2', 'tm'])
        if row['EPSVAL'] == float(HIEPS):
            file_df_dict['aver_hieps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin'])            
            if row['thpq'] < 0.0:
                file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 'tm', 'W', 'Q2'])
            file_df_dict['sep_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_sep.{}_{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), float(HIEPS)*100) \
                                                            , ['sigL', 'dsigL','sigT', 'dsigT','sigTT', 'dsigTT','sigLT', 'dsigLT', 'Q2', 'tm'])

            
################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
################################################################################################################################################

C_ratio_phi = TCanvas()
C_ratio_phi.SetGrid()
C_ratio_phi.Divide(1,NumtBins)
l_ratio_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_ratio_phi_{}".format(k+1)] = TMultiGraph()
    
    G_ratio_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_loeps']['ratio'].tolist())[i])
            G_ratio_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['aver_loeps']['dratio'].tolist())[i])
            j+=1
    G_ratio_phi_loeps.SetMarkerStyle(21)
    G_ratio_phi_loeps.SetMarkerSize(1)
    G_ratio_phi_loeps.SetMarkerColor(1)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_loeps)

    G_ratio_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_hieps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_hieps']['ratio'].tolist())[i])
            G_ratio_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['aver_hieps']['dratio'].tolist())[i])
            j+=1
    G_ratio_phi_hieps.SetMarkerStyle(21)
    G_ratio_phi_hieps.SetMarkerSize(1)
    G_ratio_phi_hieps.SetMarkerColor(2)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_hieps)    
    
    C_ratio_phi.cd(k+1)
    
    multiDict["G_ratio_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_ratio_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; Ratio".format(t_bin_centers[k]))
    
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_ratio_phi.AddEntry(G_ratio_phi_loeps,"loeps")
l_ratio_phi.AddEntry(G_ratio_phi_hieps,"hieps")
l_ratio_phi.Draw()
C_ratio_phi.Print(outputpdf + '(')

C_Q2_t = TCanvas()
C_Q2_t.SetGrid()
l_Q2_t = TLegend(0.7, 0.6, 0.9, 0.9)


G_Q2_t = TMultiGraph()
for k in range(NumtBins):

    G_Q2_t_loeps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_Q2_t_loeps.SetPoint(j, t_bins[np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
            j+=1
    G_Q2_t_loeps.SetMarkerStyle(21)
    G_Q2_t_loeps.SetMarkerSize(1)
    G_Q2_t_loeps.SetMarkerColor(1)
    G_Q2_t.Add(G_Q2_t_loeps)

    G_Q2_t_hieps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_Q2_t_hieps.SetPoint(j, t_bins[np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i]], np.array(file_df_dict['unsep_file_hieps']['Q2'].tolist())[i])
            j+=1
    G_Q2_t_hieps.SetMarkerStyle(21)
    G_Q2_t_hieps.SetMarkerSize(1)
    G_Q2_t_hieps.SetMarkerColor(2)
    G_Q2_t.Add(G_Q2_t_hieps)    
    
    G_Q2_t.Draw('AP')
    G_Q2_t.SetTitle("; t; Q2")
    
    G_Q2_t.GetYaxis().SetTitleOffset(1.5)
    G_Q2_t.GetXaxis().SetTitleOffset(1.5)
    G_Q2_t.GetXaxis().SetLabelSize(0.04)
    
l_Q2_t.AddEntry(G_Q2_t_loeps,"loeps")
l_Q2_t.AddEntry(G_Q2_t_hieps,"hieps")
l_Q2_t.Draw()
C_Q2_t.Print(outputpdf)

C_Q2_phi = TCanvas()
C_Q2_phi.SetGrid()
C_Q2_phi.Divide(1,NumtBins)
l_Q2_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_Q2_phi_{}".format(k+1)] = TMultiGraph()
    
    G_Q2_phi_loeps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_Q2_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
            j+=1
    G_Q2_phi_loeps.SetMarkerStyle(21)
    G_Q2_phi_loeps.SetMarkerSize(1)
    G_Q2_phi_loeps.SetMarkerColor(1)
    multiDict["G_Q2_phi_{}".format(k+1)].Add(G_Q2_phi_loeps)

    G_Q2_phi_hieps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_Q2_phi_hieps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_hieps']['Q2'].tolist())[i])
            j+=1
    G_Q2_phi_hieps.SetMarkerStyle(21)
    G_Q2_phi_hieps.SetMarkerSize(1)
    G_Q2_phi_hieps.SetMarkerColor(2)
    multiDict["G_Q2_phi_{}".format(k+1)].Add(G_Q2_phi_hieps)    
    
    C_Q2_phi.cd(k+1)
    
    multiDict["G_Q2_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_Q2_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; Q2".format(t_bin_centers[k]))
    
    multiDict["G_Q2_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_Q2_phi.AddEntry(G_Q2_phi_loeps,"loeps")
l_Q2_phi.AddEntry(G_Q2_phi_hieps,"hieps")
l_Q2_phi.Draw()
C_Q2_phi.Print(outputpdf)

C_W_phi = TCanvas()
C_W_phi.SetGrid()
C_W_phi.Divide(1,NumtBins)
l_W_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_W_phi_{}".format(k+1)] = TMultiGraph()
    
    G_W_phi_loeps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_W_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['W'].tolist())[i])
            j+=1
    G_W_phi_loeps.SetMarkerStyle(21)
    G_W_phi_loeps.SetMarkerSize(1)
    G_W_phi_loeps.SetMarkerColor(1)
    multiDict["G_W_phi_{}".format(k+1)].Add(G_W_phi_loeps)

    G_W_phi_hieps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_W_phi_hieps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_hieps']['W'].tolist())[i])
            j+=1
    G_W_phi_hieps.SetMarkerStyle(21)
    G_W_phi_hieps.SetMarkerSize(1)
    G_W_phi_hieps.SetMarkerColor(2)
    multiDict["G_W_phi_{}".format(k+1)].Add(G_W_phi_hieps)    
    
    C_W_phi.cd(k+1)
    
    multiDict["G_W_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_W_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; W".format(t_bin_centers[k]))
    
    multiDict["G_W_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_W_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_W_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_W_phi.AddEntry(G_W_phi_loeps,"loeps")
l_W_phi.AddEntry(G_W_phi_hieps,"hieps")
l_W_phi.Draw()
C_W_phi.Print(outputpdf)

C_xreal_thcm = TCanvas()
C_xreal_thcm.SetGrid()
C_xreal_thcm.Divide(1,NumtBins)
l_xreal_thcm = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xreal_thcm_{}".format(k+1)] = TMultiGraph()
    
    G_xreal_thcm_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_xreal_thcm_loeps.SetPoint(j, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_thcm_loeps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
            j+=1
    G_xreal_thcm_loeps.SetMarkerStyle(21)
    G_xreal_thcm_loeps.SetMarkerSize(1)
    G_xreal_thcm_loeps.SetMarkerColor(1)
    multiDict["G_xreal_thcm_{}".format(k+1)].Add(G_xreal_thcm_loeps)

    G_xreal_thcm_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_xreal_thcm_hieps.SetPoint(j, np.array(file_df_dict['unsep_file_hieps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_thcm_hieps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
            j+=1
    G_xreal_thcm_hieps.SetMarkerStyle(21)
    G_xreal_thcm_hieps.SetMarkerSize(1)
    G_xreal_thcm_hieps.SetMarkerColor(2)
    multiDict["G_xreal_thcm_{}".format(k+1)].Add(G_xreal_thcm_hieps)    
    
    C_xreal_thcm.cd(k+1)
    
    multiDict["G_xreal_thcm_{}".format(k+1)].Draw('AP')
    multiDict["G_xreal_thcm_{}".format(k+1)].SetTitle("t = {:.2f} ; #theta_{{cm}}; xmod".format(t_bin_centers[k]))
    
    multiDict["G_xreal_thcm_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_thcm_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_thcm_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xreal_thcm.AddEntry(G_xreal_thcm_loeps,"loeps")
l_xreal_thcm.AddEntry(G_xreal_thcm_hieps,"hieps")
l_xreal_thcm.Draw()
C_xreal_thcm.Print(outputpdf)

C_xmod_thcm = TCanvas()
C_xmod_thcm.SetGrid()
C_xmod_thcm.Divide(1,NumtBins)
l_xmod_thcm = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xmod_thcm_{}".format(k+1)] = TMultiGraph()
    
    G_xmod_thcm_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_xmod_thcm_loeps.SetPoint(j, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
            G_xmod_thcm_loeps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
            j+=1
    G_xmod_thcm_loeps.SetMarkerStyle(21)
    G_xmod_thcm_loeps.SetMarkerSize(1)
    G_xmod_thcm_loeps.SetMarkerColor(1)
    multiDict["G_xmod_thcm_{}".format(k+1)].Add(G_xmod_thcm_loeps)

    G_xmod_thcm_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_xmod_thcm_hieps.SetPoint(j, np.array(file_df_dict['unsep_file_hieps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
            G_xmod_thcm_hieps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
            j+=1
    G_xmod_thcm_hieps.SetMarkerStyle(21)
    G_xmod_thcm_hieps.SetMarkerSize(1)
    G_xmod_thcm_hieps.SetMarkerColor(2)
    multiDict["G_xmod_thcm_{}".format(k+1)].Add(G_xmod_thcm_hieps)    
    
    C_xmod_thcm.cd(k+1)
    
    multiDict["G_xmod_thcm_{}".format(k+1)].Draw('AP')
    multiDict["G_xmod_thcm_{}".format(k+1)].SetTitle("t = {:.2f} ; #theta_{{cm}}; xmod".format(t_bin_centers[k]))
    
    multiDict["G_xmod_thcm_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_thcm_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_thcm_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xmod_thcm.AddEntry(G_xmod_thcm_loeps,"loeps")
l_xmod_thcm.AddEntry(G_xmod_thcm_hieps,"hieps")
l_xmod_thcm.Draw()
C_xmod_thcm.Print(outputpdf)

C_xreal_phi = TCanvas()
C_xreal_phi.SetGrid()
C_xreal_phi.Divide(1,NumtBins)
l_xreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)

'''

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xreal_phi_{}".format(k+1)] = TMultiGraph()
    
    G_xreal_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_xreal_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
            j+=1
    G_xreal_phi_loeps.SetMarkerStyle(21)
    G_xreal_phi_loeps.SetMarkerSize(1)
    G_xreal_phi_loeps.SetMarkerColor(1)
    multiDict["G_xreal_phi_{}".format(k+1)].Add(G_xreal_phi_loeps)
    
    C_xreal_phi.cd(k+1)
    
    multiDict["G_xreal_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xreal_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xreal".format(t_bin_centers[k]))
    
    multiDict["G_xreal_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xreal_phi.AddEntry(G_xreal_phi_loeps,"loeps")
l_xreal_phi.Draw()
C_xreal_phi.Print(outputpdf)

'''

C_xreal_phi = TCanvas()
C_xreal_phi.SetGrid()
C_xreal_phi.Divide(1,NumtBins)
l_xreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xreal_phi_{}".format(k+1)] = TMultiGraph()

    G_xreal_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_xreal_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
            j+=1
    G_xreal_phi_loeps.SetMarkerStyle(21)
    G_xreal_phi_loeps.SetMarkerSize(1)
    G_xreal_phi_loeps.SetMarkerColor(1)
    multiDict["G_xreal_phi_{}".format(k+1)].Add(G_xreal_phi_loeps)
    
    G_xreal_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_xreal_phi_hieps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
            j+=1
    G_xreal_phi_hieps.SetMarkerStyle(21)
    G_xreal_phi_hieps.SetMarkerSize(1)
    G_xreal_phi_hieps.SetMarkerColor(2)
    multiDict["G_xreal_phi_{}".format(k+1)].Add(G_xreal_phi_hieps)    
    
    C_xreal_phi.cd(k+1)
    
    multiDict["G_xreal_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xreal_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xreal".format(t_bin_centers[k]))
    
    multiDict["G_xreal_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xreal_phi.AddEntry(G_xreal_phi_hieps,"hieps")
l_xreal_phi.Draw()
C_xreal_phi.Print(outputpdf)

C_xmod_phi = TCanvas()
C_xmod_phi.SetGrid()
C_xmod_phi.Divide(1,NumtBins)
l_xmod_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xmod_phi_{}".format(k+1)] = TMultiGraph()
    
    G_xmod_phi_loeps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_xmod_phi_loeps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
            j+=1
    G_xmod_phi_loeps.SetMarkerStyle(21)
    G_xmod_phi_loeps.SetMarkerSize(1)
    G_xmod_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmod_phi_{}".format(k+1)].Add(G_xmod_phi_loeps)

    G_xmod_phi_hieps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_xmod_phi_hieps.SetPoint(j, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
            j+=1
    G_xmod_phi_hieps.SetMarkerStyle(21)
    G_xmod_phi_hieps.SetMarkerSize(1)
    G_xmod_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmod_phi_{}".format(k+1)].Add(G_xmod_phi_hieps)    
    
    C_xmod_phi.cd(k+1)
    
    multiDict["G_xmod_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xmod_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xmod".format(t_bin_centers[k]))
    
    multiDict["G_xmod_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xmod_phi.AddEntry(G_xmod_phi_loeps,"loeps")
l_xmod_phi.AddEntry(G_xmod_phi_hieps,"hieps")
l_xmod_phi.Draw()
C_xmod_phi.Print(outputpdf)

C_sigl_phi = TCanvas()
C_sigl_phi.Divide(1,NumtBins)
C_sigl_phi.SetGrid()
l_sigl_phi = TLegend(0.8,0.8,0.95,0.95)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_sigl_phi_{}".format(k+1)] = TMultiGraph()
    
    G_sigl_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_sigl_phi_loeps.SetPoint(j, phi_bin_centers[file_df_dict['aver_loeps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_loeps']['sigL'].tolist())[i])
            G_sigl_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['sep_file_loeps']['dsigL'].tolist())[i])
            j+=1
    G_sigl_phi_loeps.SetMarkerColor(1)
    G_sigl_phi_loeps.SetMarkerStyle(21)
    G_sigl_phi_loeps.SetMarkerSize(1)
    multiDict["G_sigl_phi_{}".format(k+1)].Add(G_sigl_phi_loeps)

    G_sigl_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_sigl_phi_hieps.SetPoint(j, phi_bin_centers[file_df_dict['aver_hieps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_hieps']['sigL'].tolist())[i])
            G_sigl_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['sep_file_hieps']['dsigL'].tolist())[i])
            j+=1
    G_sigl_phi_hieps.SetMarkerColor(2)
    G_sigl_phi_hieps.SetMarkerStyle(21)
    G_sigl_phi_hieps.SetMarkerSize(1)
    multiDict["G_sigl_phi_{}".format(k+1)].Add(G_sigl_phi_hieps)    

    C_sigl_phi.cd(k+1)

    multiDict["G_sigl_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_sigl_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; sigl".format(t_bin_centers[k]))
    
    multiDict["G_sigl_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_sigl_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_sigl_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)

l_sigl_phi.AddEntry(G_sigl_phi_loeps,"loeps")
l_sigl_phi.AddEntry(G_sigl_phi_hieps,"hieps")
l_sigl_phi.Draw()    
C_sigl_phi.Print(outputpdf)

C_sigt_phi = TCanvas()
C_sigt_phi.Divide(1,NumtBins)
C_sigt_phi.SetGrid()
l_sigt_phi = TLegend(0.8,0.8,0.95,0.95)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_sigt_phi_{}".format(k+1)] = TMultiGraph()
    
    G_sigt_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_sigt_phi_loeps.SetPoint(j, phi_bin_centers[file_df_dict['aver_loeps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_loeps']['sigT'].tolist())[i])
            G_sigt_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['sep_file_loeps']['dsigT'].tolist())[i])
            j+=1
    G_sigt_phi_loeps.SetMarkerColor(1)
    G_sigt_phi_loeps.SetMarkerStyle(21)
    G_sigt_phi_loeps.SetMarkerSize(1)
    multiDict["G_sigt_phi_{}".format(k+1)].Add(G_sigt_phi_loeps)

    G_sigt_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_sigt_phi_hieps.SetPoint(j, phi_bin_centers[file_df_dict['aver_hieps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_hieps']['sigT'].tolist())[i])
            G_sigt_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['sep_file_hieps']['dsigT'].tolist())[i])
            j+=1
    G_sigt_phi_hieps.SetMarkerColor(2)
    G_sigt_phi_hieps.SetMarkerStyle(21)
    G_sigt_phi_hieps.SetMarkerSize(1)
    multiDict["G_sigt_phi_{}".format(k+1)].Add(G_sigt_phi_hieps)    

    C_sigt_phi.cd(k+1)

    multiDict["G_sigt_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_sigt_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; sigt".format(t_bin_centers[k]))
    
    multiDict["G_sigt_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_sigt_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_sigt_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)

l_sigt_phi.AddEntry(G_sigt_phi_loeps,"loeps")
l_sigt_phi.AddEntry(G_sigt_phi_hieps,"hieps")
l_sigt_phi.Draw()    
C_sigt_phi.Print(outputpdf)

C_siglt_phi = TCanvas()
C_siglt_phi.Divide(1,NumtBins)
C_siglt_phi.SetGrid()
l_siglt_phi = TLegend(0.8,0.8,0.95,0.95)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_siglt_phi_{}".format(k+1)] = TMultiGraph()
    
    G_siglt_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_siglt_phi_loeps.SetPoint(j, phi_bin_centers[file_df_dict['aver_loeps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_loeps']['sigLT'].tolist())[i])
            G_siglt_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['sep_file_loeps']['dsigLT'].tolist())[i])
            j+=1
    G_siglt_phi_loeps.SetMarkerColor(1)
    G_siglt_phi_loeps.SetMarkerStyle(21)
    G_siglt_phi_loeps.SetMarkerSize(1)
    multiDict["G_siglt_phi_{}".format(k+1)].Add(G_siglt_phi_loeps)
    G_siglt_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_siglt_phi_loeps.SetPoint(j, np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigLT'].tolist())[i])
            G_siglt_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['sep_file_loeps']['dsigLT'].tolist())[i])
            j+=1
    G_siglt_phi_loeps.SetMarkerColor(1)
    G_siglt_phi_loeps.SetMarkerStyle(21)
    G_siglt_phi_loeps.SetMarkerSize(1)
    multiDict["G_siglt_phi_{}".format(k+1)].Add(G_siglt_phi_loeps)

    G_siglt_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_siglt_phi_hieps.SetPoint(j, phi_bin_centers[file_df_dict['aver_hieps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_hieps']['sigLT'].tolist())[i])
            G_siglt_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['sep_file_hieps']['dsigLT'].tolist())[i])
            j+=1
    G_siglt_phi_hieps.SetMarkerColor(2)
    G_siglt_phi_hieps.SetMarkerStyle(21)
    G_siglt_phi_hieps.SetMarkerSize(1)
    multiDict["G_siglt_phi_{}".format(k+1)].Add(G_siglt_phi_hieps)    

    C_siglt_phi.cd(k+1)

    multiDict["G_siglt_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_siglt_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; siglt".format(t_bin_centers[k]))
    
    multiDict["G_siglt_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_siglt_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_siglt_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)

l_siglt_phi.AddEntry(G_siglt_phi_loeps,"loeps")
l_siglt_phi.AddEntry(G_siglt_phi_hieps,"hieps")
l_siglt_phi.Draw()    
C_siglt_phi.Print(outputpdf)

C_sigtt_phi = TCanvas()
C_sigtt_phi.Divide(1,NumtBins)
C_sigtt_phi.SetGrid()
l_sigtt_phi = TLegend(0.8,0.8,0.95,0.95)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_sigtt_phi_{}".format(k+1)] = TMultiGraph()
    
    G_sigtt_phi_loeps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_sigtt_phi_loeps.SetPoint(j, np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigTT'].tolist())[i])
            G_sigtt_phi_loeps.SetPointError(j, 0, np.array(file_df_dict['sep_file_loeps']['dsigTT'].tolist())[i])
            j+=1
    G_sigtt_phi_loeps.SetMarkerColor(1)
    G_sigtt_phi_loeps.SetMarkerStyle(21)
    G_sigtt_phi_loeps.SetMarkerSize(1)
    multiDict["G_sigtt_phi_{}".format(k+1)].Add(G_sigtt_phi_loeps)

    G_sigtt_phi_hieps = TGraphErrors()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_sigtt_phi_hieps.SetPoint(j, phi_bin_centers[file_df_dict['aver_hieps']['phibin'].tolist()[i]], np.array(file_df_dict['sep_file_hieps']['sigTT'].tolist())[i])
            G_sigtt_phi_hieps.SetPointError(j, 0, np.array(file_df_dict['sep_file_hieps']['dsigTT'].tolist())[i])
            j+=1
    G_sigtt_phi_hieps.SetMarkerColor(2)
    G_sigtt_phi_hieps.SetMarkerStyle(21)
    G_sigtt_phi_hieps.SetMarkerSize(1)
    multiDict["G_sigtt_phi_{}".format(k+1)].Add(G_sigtt_phi_hieps)    

    C_sigtt_phi.cd(k+1)

    multiDict["G_sigtt_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_sigtt_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; sigtt".format(t_bin_centers[k]))
    
    multiDict["G_sigtt_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_sigtt_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_sigtt_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)

l_sigtt_phi.AddEntry(G_sigtt_phi_loeps,"loeps")
l_sigtt_phi.AddEntry(G_sigtt_phi_hieps,"hieps")
l_sigtt_phi.Draw()    
C_sigtt_phi.Print(outputpdf)

C_sigl_t = TCanvas()
C_sigl_t.SetGrid()

G_sigl_t = TGraphErrors()

l_sigl_t = TLegend(0.8,0.8,0.95,0.95)

G_sigl_t.SetTitle(" ; t; sigL")

for i in range(len(file_df_dict['sep_file_loeps']['sigL'].tolist())):
    G_sigl_t.SetPoint(i, np.array(file_df_dict['sep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigL'].tolist())[i])
    G_sigl_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_loeps']['dsigL'].tolist())[i])
    G_sigl_t.SetMarkerColor(1)

for i in range(len(file_df_dict['sep_file_hieps']['sigL'].tolist())):
    G_sigl_t.SetPoint(i, np.array(file_df_dict['sep_file_hieps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_hieps']['sigL'].tolist())[i])
    G_sigl_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_hieps']['dsigL'].tolist())[i])
    G_sigl_t.SetMarkerColor(2)
    
G_sigl_t.SetMarkerStyle(21)
G_sigl_t.SetMarkerSize(1)
    
G_sigl_t.Draw('AP')

C_sigl_t.Print(outputpdf)

C_sigt_t = TCanvas()
C_sigt_t.SetGrid()

G_sigt_t = TGraphErrors()

l_sigt_t = TLegend(0.8,0.8,0.95,0.95)

G_sigt_t.SetTitle(" ; t; sigt")

for i in range(len(file_df_dict['sep_file_loeps']['sigT'].tolist())):
    G_sigt_t.SetPoint(i, np.array(file_df_dict['sep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigT'].tolist())[i])
    G_sigt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_loeps']['dsigT'].tolist())[i])
    G_sigt_t.SetMarkerColor(1)

for i in range(len(file_df_dict['sep_file_hieps']['sigT'].tolist())):
    G_sigt_t.SetPoint(i, np.array(file_df_dict['sep_file_hieps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_hieps']['sigT'].tolist())[i])
    G_sigt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_hieps']['dsigT'].tolist())[i])
    G_sigt_t.SetMarkerColor(2)
    
G_sigt_t.SetMarkerStyle(21)
G_sigt_t.SetMarkerSize(1)
    
G_sigt_t.Draw('AP')

C_sigt_t.Print(outputpdf)

C_sigtt_t = TCanvas()
C_sigtt_t.SetGrid()

G_sigtt_t = TGraphErrors()

l_sigtt_t = TLegend(0.8,0.8,0.95,0.95)

G_sigtt_t.SetTitle(" ; t; sigtt")

for i in range(len(file_df_dict['sep_file_loeps']['sigTT'].tolist())):
    G_sigtt_t.SetPoint(i, np.array(file_df_dict['sep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigTT'].tolist())[i])
    G_sigtt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_loeps']['dsigTT'].tolist())[i])
    G_sigtt_t.SetMarkerColor(1)

for i in range(len(file_df_dict['sep_file_hieps']['sigTT'].tolist())):
    G_sigtt_t.SetPoint(i, np.array(file_df_dict['sep_file_hieps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_hieps']['sigTT'].tolist())[i])
    G_sigtt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_hieps']['dsigTT'].tolist())[i])
    G_sigtt_t.SetMarkerColor(2)
    
G_sigtt_t.SetMarkerStyle(21)
G_sigtt_t.SetMarkerSize(1)
    
G_sigtt_t.Draw('AP')

C_sigtt_t.Print(outputpdf)

C_siglt_t = TCanvas()
C_siglt_t.SetGrid()

G_siglt_t = TGraphErrors()

l_siglt_t = TLegend(0.8,0.8,0.95,0.95)

G_siglt_t.SetTitle(" ; t; siglt")

for i in range(len(file_df_dict['sep_file_loeps']['sigLT'].tolist())):
    G_siglt_t.SetPoint(i, np.array(file_df_dict['sep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_loeps']['sigLT'].tolist())[i])
    G_siglt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_loeps']['dsigLT'].tolist())[i])
    G_siglt_t.SetMarkerColor(1)

for i in range(len(file_df_dict['sep_file_hieps']['sigLT'].tolist())):
    G_siglt_t.SetPoint(i, np.array(file_df_dict['sep_file_hieps']['tm'].tolist())[i], np.array(file_df_dict['sep_file_hieps']['sigLT'].tolist())[i])
    G_siglt_t.SetPointError(i, 0, np.array(file_df_dict['sep_file_hieps']['dsigLT'].tolist())[i])
    G_siglt_t.SetMarkerColor(2)
    
G_siglt_t.SetMarkerStyle(21)
G_siglt_t.SetMarkerSize(1)
    
G_siglt_t.Draw('AP')

C_siglt_t.Print(outputpdf+')')
