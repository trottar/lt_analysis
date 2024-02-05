#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-05 13:44:57 trottar"
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
from ROOT import TCanvas, TGraph, TGraphErrors, TMultiGraph, TLegend, TLine
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
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
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
fix_spacing(LTANAPATH+"/src/{}/averages/avek.Q{}W{}.dat".format(ParticleType, Q2.replace("p",""), W.replace("p","")))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")))
################################################################################################################################################
# Read in files and convert to dataframes

file_df_dict = {}

setting_file = LTANAPATH+"/src/{}/list.settings".format(ParticleType)
file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'W', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins'])

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

t_bin_centers = (t_bins[:-1] + t_bins[1:]) / 2

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

#phi_bin_centers = (phi_bins[:-1] + phi_bins[1:]) / 2
phi_bin_centers = phi_bins
    
for i,row in file_df_dict['setting_df'].iterrows():
    if row['Q2'] == float(Q2.replace("p",".")):
        file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/{}/beam/Eb_KLT.dat".format(ParticleType), ['ebeam', 'Q2', 'W', 'EPSVAL'])
        file_df_dict['avek_file'] = file_to_df(LTANAPATH+"/src/{}/averages/avek.Q{}W{}.dat".format(ParticleType, Q2.replace("p",""), W.replace("p","")) \
                                               , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])        
        
        if row['EPSVAL'] == float(LOEPS):
            file_df_dict['aver_loeps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
            if row['thpq'] < 0.0:
                file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                       float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                      float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                        float(LOEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 't_min', 'W', 'Q2']).sort_values(by='t')

        if row['EPSVAL'] == float(HIEPS):
            file_df_dict['aver_hieps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
            if row['thpq'] < 0.0:
                file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                       float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                      float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                        float(HIEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 't_min', 'W', 'Q2']).sort_values(by='t')
        file_df_dict['sep_file'] = file_to_df( \
                                               LTANAPATH+"/src/{}/xsects/x_sep.{}_Q{}W{}.dat" \
                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")) \
                                               , ['sigT', 'dsigT', 'sigL', 'dsigL', 'sigLT', 'dsigLT', 'sigTT', 'dsigTT', 'chisq', 't', 'tm', 'W', 'Q2'])
            
################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
################################################################################################################################################

C_ratio_phi = TCanvas()
C_ratio_phi.SetGrid()
C_ratio_phi.Divide(1,NumtBins)
l_ratio_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    # Add a gray line at unity
    line_at_unity = TLine(0.0, 1.0, 360.0, 1.0)
    line_at_unity.SetLineColor(7)
    line_at_unity.SetLineStyle(2)   # Dashed line style
    line_at_unity.Draw()
    
    multiDict["G_ratio_phi_{}".format(k+1)] = TMultiGraph()
    
    G_ratio_phi_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_loeps.SetPoint(i, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_loeps']['ratio'].tolist())[i])
            G_ratio_phi_loeps.SetPointError(i, 0, np.array(file_df_dict['aver_loeps']['dratio'].tolist())[i])
    G_ratio_phi_loeps.SetMarkerStyle(21)
    G_ratio_phi_loeps.SetMarkerSize(1)
    G_ratio_phi_loeps.SetMarkerColor(1)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_loeps)

    G_ratio_phi_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_hieps.SetPoint(i, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_hieps']['ratio'].tolist())[i])
            G_ratio_phi_hieps.SetPointError(i, 0, np.array(file_df_dict['aver_hieps']['dratio'].tolist())[i])
    G_ratio_phi_hieps.SetMarkerStyle(23)
    G_ratio_phi_hieps.SetMarkerSize(1)
    G_ratio_phi_hieps.SetMarkerColor(2)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_hieps)    
    
    C_ratio_phi.cd(k+1)
    
    multiDict["G_ratio_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_ratio_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; Ratio".format(t_bin_centers[k]))
    
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetRangeUser(0.0, 2.0)    
    
l_ratio_phi.AddEntry(G_ratio_phi_loeps,"loeps")
l_ratio_phi.AddEntry(G_ratio_phi_hieps,"hieps")
l_ratio_phi.Draw()
C_ratio_phi.Print(outputpdf + '(')

multiDict = {}
for k in range(NumtBins):

    C_ratio_phi = TCanvas()
    C_ratio_phi.SetGrid()
    l_ratio_phi = TLegend(0.7, 0.6, 0.9, 0.9)

    multiDict["G_ratio_phi_{}".format(k+1)] = TMultiGraph()
    
    G_ratio_phi_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_loeps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_loeps.SetPoint(i, phi_bin_centers[np.array(file_df_dict['aver_loeps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_loeps']['ratio'].tolist())[i])
            G_ratio_phi_loeps.SetPointError(i, 0, np.array(file_df_dict['aver_loeps']['dratio'].tolist())[i])
    G_ratio_phi_loeps.SetMarkerStyle(21)
    G_ratio_phi_loeps.SetMarkerSize(1)
    G_ratio_phi_loeps.SetMarkerColor(1)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_loeps)

    G_ratio_phi_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['aver_hieps']['tbin'].tolist())[i] == (k+1):
            G_ratio_phi_hieps.SetPoint(i, phi_bin_centers[np.array(file_df_dict['aver_hieps']['phibin'].tolist())[i]], np.array(file_df_dict['aver_hieps']['ratio'].tolist())[i])
            G_ratio_phi_hieps.SetPointError(i, 0, np.array(file_df_dict['aver_hieps']['dratio'].tolist())[i])
    G_ratio_phi_hieps.SetMarkerStyle(23)
    G_ratio_phi_hieps.SetMarkerSize(1)
    G_ratio_phi_hieps.SetMarkerColor(2)
    multiDict["G_ratio_phi_{}".format(k+1)].Add(G_ratio_phi_hieps)
    
    multiDict["G_ratio_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_ratio_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; Ratio".format(t_bin_centers[k]))
    
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetRangeUser(0.0, 2.0)
    
    # Add a gray line at unity
    line_at_unity = TLine(0.0, 1.0, 360.0, 1.0)
    line_at_unity.SetLineColor(7)  # 7 corresponds to light blue
    line_at_unity.SetLineStyle(2)   # Dashed line style
    line_at_unity.Draw("same")
    
    l_ratio_phi.AddEntry(G_ratio_phi_loeps,"loeps")
    l_ratio_phi.AddEntry(G_ratio_phi_hieps,"hieps")
    #l_ratio_phi.Draw()
    C_ratio_phi.Print(outputpdf)

C_Q2_t = TCanvas()
C_Q2_t.SetGrid()
l_Q2_t = TLegend(0.7, 0.6, 0.9, 0.9)

G_Q2_t = TMultiGraph()
for k in range(NumtBins):

    G_Q2_t_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_Q2_t_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
    G_Q2_t_loeps.SetMarkerStyle(21)
    G_Q2_t_loeps.SetMarkerSize(1)
    G_Q2_t_loeps.SetMarkerColor(1)
    G_Q2_t.Add(G_Q2_t_loeps)

    G_Q2_t_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_Q2_t_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['Q2'].tolist())[i])
    G_Q2_t_hieps.SetMarkerStyle(23)
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
#l_Q2_t.Draw()
C_Q2_t.Print(outputpdf)

C_Q2_phi = TCanvas()
C_Q2_phi.SetGrid()
C_Q2_phi.Divide(1,NumtBins)
l_Q2_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_Q2_phi_{}".format(k+1)] = TMultiGraph()
    
    G_Q2_phi_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_Q2_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
    G_Q2_phi_loeps.SetMarkerStyle(21)
    G_Q2_phi_loeps.SetMarkerSize(1)
    G_Q2_phi_loeps.SetMarkerColor(1)
    multiDict["G_Q2_phi_{}".format(k+1)].Add(G_Q2_phi_loeps)

    G_Q2_phi_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_Q2_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['Q2'].tolist())[i])
    G_Q2_phi_hieps.SetMarkerStyle(23)
    G_Q2_phi_hieps.SetMarkerSize(1)
    G_Q2_phi_hieps.SetMarkerColor(2)
    multiDict["G_Q2_phi_{}".format(k+1)].Add(G_Q2_phi_hieps)    
    
    C_Q2_phi.cd(k+1)

    multiDict["G_Q2_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_Q2_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; Q2".format(t_bin_centers[k]))
    
    multiDict["G_Q2_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_Q2_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
l_Q2_phi.AddEntry(G_Q2_phi_loeps,"loeps")
l_Q2_phi.AddEntry(G_Q2_phi_hieps,"hieps")
#l_Q2_phi.Draw()
C_Q2_phi.Print(outputpdf)

C_W_phi = TCanvas()
C_W_phi.SetGrid()
C_W_phi.Divide(1,NumtBins)
l_W_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_W_phi_{}".format(k+1)] = TMultiGraph()
    
    G_W_phi_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_W_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['W'].tolist())[i])
    G_W_phi_loeps.SetMarkerStyle(21)
    G_W_phi_loeps.SetMarkerSize(1)
    G_W_phi_loeps.SetMarkerColor(1)
    multiDict["G_W_phi_{}".format(k+1)].Add(G_W_phi_loeps)

    G_W_phi_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_W_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['W'].tolist())[i])
    G_W_phi_hieps.SetMarkerStyle(23)
    G_W_phi_hieps.SetMarkerSize(1)
    G_W_phi_hieps.SetMarkerColor(2)
    multiDict["G_W_phi_{}".format(k+1)].Add(G_W_phi_hieps)    
    
    C_W_phi.cd(k+1)

    multiDict["G_W_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_W_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; W".format(t_bin_centers[k]))
    
    multiDict["G_W_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_W_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_W_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_W_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
l_W_phi.AddEntry(G_W_phi_loeps,"loeps")
l_W_phi.AddEntry(G_W_phi_hieps,"hieps")
#l_W_phi.Draw()
C_W_phi.Print(outputpdf)

C_xreal_thcm = TCanvas()
C_xreal_thcm.SetGrid()
C_xreal_thcm.Divide(1,NumtBins)
l_xreal_thcm = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xreal_thcm_{}".format(k+1)] = TMultiGraph()
    
    G_xreal_thcm_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xreal_thcm_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_thcm_loeps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    G_xreal_thcm_loeps.SetMarkerStyle(21)
    G_xreal_thcm_loeps.SetMarkerSize(1)
    G_xreal_thcm_loeps.SetMarkerColor(1)
    multiDict["G_xreal_thcm_{}".format(k+1)].Add(G_xreal_thcm_loeps)

    G_xreal_thcm_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xreal_thcm_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_thcm_hieps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
    G_xreal_thcm_hieps.SetMarkerStyle(23)
    G_xreal_thcm_hieps.SetMarkerSize(1)
    G_xreal_thcm_hieps.SetMarkerColor(2)
    multiDict["G_xreal_thcm_{}".format(k+1)].Add(G_xreal_thcm_hieps)    
    
    C_xreal_thcm.cd(k+1)
    
    multiDict["G_xreal_thcm_{}".format(k+1)].Draw('AP')
    multiDict["G_xreal_thcm_{}".format(k+1)].SetTitle("t = {:.2f} ; #theta_{{cm}}; xreal".format(t_bin_centers[k]))
    
    multiDict["G_xreal_thcm_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_thcm_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_thcm_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_xreal_thcm.AddEntry(G_xreal_thcm_loeps,"loeps")
l_xreal_thcm.AddEntry(G_xreal_thcm_hieps,"hieps")
#l_xreal_thcm.Draw()
C_xreal_thcm.Print(outputpdf)

C_xmod_thcm = TCanvas()
C_xmod_thcm.SetGrid()
C_xmod_thcm.Divide(1,NumtBins)
l_xmod_thcm = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xmod_thcm_{}".format(k+1)] = TMultiGraph()
    
    G_xmod_thcm_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xmod_thcm_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
            G_xmod_thcm_loeps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    G_xmod_thcm_loeps.SetMarkerStyle(21)
    G_xmod_thcm_loeps.SetMarkerSize(1)
    G_xmod_thcm_loeps.SetMarkerColor(1)
    multiDict["G_xmod_thcm_{}".format(k+1)].Add(G_xmod_thcm_loeps)

    G_xmod_thcm_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xmod_thcm_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
            G_xmod_thcm_hieps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
    G_xmod_thcm_hieps.SetMarkerStyle(23)
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
#l_xmod_thcm.Draw()
C_xmod_thcm.Print(outputpdf)

C_xreal_phi = TCanvas()
C_xreal_phi.SetGrid()
C_xreal_phi.Divide(1,NumtBins)
l_xreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)

C_xreal_phi = TCanvas()
C_xreal_phi.SetGrid()
C_xreal_phi.Divide(1,NumtBins)
l_xreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xreal_phi_{}".format(k+1)] = TMultiGraph()

    G_xreal_phi_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xreal_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_phi_loeps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    G_xreal_phi_loeps.SetMarkerStyle(21)
    G_xreal_phi_loeps.SetMarkerSize(1)
    G_xreal_phi_loeps.SetMarkerColor(1)
    multiDict["G_xreal_phi_{}".format(k+1)].Add(G_xreal_phi_loeps)
    
    G_xreal_phi_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xreal_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_phi_hieps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
    G_xreal_phi_hieps.SetMarkerStyle(23)
    G_xreal_phi_hieps.SetMarkerSize(1)
    G_xreal_phi_hieps.SetMarkerColor(2)
    multiDict["G_xreal_phi_{}".format(k+1)].Add(G_xreal_phi_hieps)    
    
    C_xreal_phi.cd(k+1)

    multiDict["G_xreal_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xreal_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xreal".format(t_bin_centers[k]))
    
    multiDict["G_xreal_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_xreal_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
l_xreal_phi.AddEntry(G_xreal_phi_hieps,"hieps")
#l_xreal_phi.Draw()
C_xreal_phi.Print(outputpdf)

C_xmod_phi = TCanvas()
C_xmod_phi.SetGrid()
C_xmod_phi.Divide(1,NumtBins)
l_xmod_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xmod_phi_{}".format(k+1)] = TMultiGraph()
    
    G_xmod_phi_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xmod_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    G_xmod_phi_loeps.SetMarkerStyle(21)
    G_xmod_phi_loeps.SetMarkerSize(1)
    G_xmod_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmod_phi_{}".format(k+1)].Add(G_xmod_phi_loeps)

    G_xmod_phi_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xmod_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
    G_xmod_phi_hieps.SetMarkerStyle(23)
    G_xmod_phi_hieps.SetMarkerSize(1)
    G_xmod_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmod_phi_{}".format(k+1)].Add(G_xmod_phi_hieps)    
    
    C_xmod_phi.cd(k+1)

    multiDict["G_xmod_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xmod_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xmod".format(t_bin_centers[k]))
    
    multiDict["G_xmod_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xmod_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_xmod_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
l_xmod_phi.AddEntry(G_xmod_phi_loeps,"loeps")
l_xmod_phi.AddEntry(G_xmod_phi_hieps,"hieps")
#l_xmod_phi.Draw()
C_xmod_phi.Print(outputpdf)


multiDict = {}
for k in range(NumtBins):

    C_xmodreal_phi = TCanvas()
    C_xmodreal_phi.SetGrid()
    l_xmodreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)    
    
    multiDict["G_xmodreal_phi_{}".format(k+1)] = TMultiGraph()

    G_xreal_phi_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xreal_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_phi_loeps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    G_xreal_phi_loeps.SetMarkerStyle(21)
    G_xreal_phi_loeps.SetMarkerSize(1)
    G_xreal_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xreal_phi_loeps)
    
    G_xreal_phi_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xreal_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_phi_hieps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
    G_xreal_phi_hieps.SetMarkerStyle(23)
    G_xreal_phi_hieps.SetMarkerSize(1)
    G_xreal_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xreal_phi_hieps)    
    
    G_xmod_phi_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xmod_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    G_xmod_phi_loeps.SetMarkerStyle(30)
    G_xmod_phi_loeps.SetMarkerSize(1)
    G_xmod_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xmod_phi_loeps)

    G_xmod_phi_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumtBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xmod_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
    G_xmod_phi_hieps.SetMarkerStyle(27)
    G_xmod_phi_hieps.SetMarkerSize(1)
    G_xmod_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xmod_phi_hieps)

    multiDict["G_xmodreal_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xmodreal_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xsect".format(t_bin_centers[k]))
    
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
    l_xmodreal_phi.AddEntry(G_xreal_phi_loeps,"loeps")
    l_xmodreal_phi.AddEntry(G_xreal_phi_hieps,"hieps")    
    l_xmodreal_phi.AddEntry(G_xmod_phi_loeps,"loeps model")
    l_xmodreal_phi.AddEntry(G_xmod_phi_hieps,"hieps model")
    #l_xmodreal_phi.Draw()
    C_xmodreal_phi.Print(outputpdf)

C_xmodreal_phi = TCanvas()
C_xmodreal_phi.SetGrid()
C_xmodreal_phi.Divide(1,NumtBins)
l_xmodreal_phi = TLegend(0.7, 0.6, 0.9, 0.9)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_xmodreal_phi_{}".format(k+1)] = TMultiGraph()

    G_xreal_phi_loeps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumPhiBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xreal_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
            G_xreal_phi_loeps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    G_xreal_phi_loeps.SetMarkerStyle(21)
    G_xreal_phi_loeps.SetMarkerSize(1)
    G_xreal_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xreal_phi_loeps)
    
    G_xreal_phi_hieps = TGraphErrors()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumPhiBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xreal_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_real'].tolist())[i])
            G_xreal_phi_hieps.SetPointError(i, 0, np.array(file_df_dict['unsep_file_hieps']['dx_real'].tolist())[i])
    G_xreal_phi_hieps.SetMarkerStyle(23)
    G_xreal_phi_hieps.SetMarkerSize(1)
    G_xreal_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xreal_phi_hieps)    
    
    G_xmod_phi_loeps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[k*NumPhiBins+int(i/NumPhiBins)] == np.array(file_df_dict['unsep_file_loeps']['t'].tolist())[i]:
            G_xmod_phi_loeps.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    G_xmod_phi_loeps.SetMarkerStyle(30)
    G_xmod_phi_loeps.SetMarkerSize(1)
    G_xmod_phi_loeps.SetMarkerColor(1)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xmod_phi_loeps)

    G_xmod_phi_hieps = TGraph()
    for i in range(NumtBins*NumPhiBins):
        if np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[k*NumPhiBins+int(i/NumPhiBins)] == np.array(file_df_dict['unsep_file_hieps']['t'].tolist())[i]:
            G_xmod_phi_hieps.SetPoint(i, np.array(file_df_dict['unsep_file_hieps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['x_mod'].tolist())[i])
    G_xmod_phi_hieps.SetMarkerStyle(27)
    G_xmod_phi_hieps.SetMarkerSize(1)
    G_xmod_phi_hieps.SetMarkerColor(2)
    multiDict["G_xmodreal_phi_{}".format(k+1)].Add(G_xmod_phi_hieps)
    
    C_xmodreal_phi.cd(k+1)

    multiDict["G_xmodreal_phi_{}".format(k+1)].Draw('AP')
    multiDict["G_xmodreal_phi_{}".format(k+1)].SetTitle("t = {:.2f} ; #phi; xsect".format(t_bin_centers[k]))
    
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    multiDict["G_xmodreal_phi_{}".format(k+1)].GetXaxis().SetRangeUser(0, 360)
    
l_xmodreal_phi.AddEntry(G_xreal_phi_loeps,"loeps")
l_xmodreal_phi.AddEntry(G_xreal_phi_hieps,"hieps")    
l_xmodreal_phi.AddEntry(G_xmod_phi_loeps,"loeps model")
l_xmodreal_phi.AddEntry(G_xmod_phi_hieps,"hieps model")
#l_xmodreal_phi.Draw()
C_xmodreal_phi.Print(outputpdf+')')
