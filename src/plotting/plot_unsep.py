#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-27 17:35:18 trottar"
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
l_ratio_phi = TLegend(0.115,0.35,0.33,0.5)

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
    multiDict["G_ratio_phi_{}".format(k+1)].SetTitle("t = {} ; #phi; Ratio".format(t_bin_centers[k]))
    
    multiDict["G_ratio_phi_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_ratio_phi_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_ratio_phi.AddEntry(G_ratio_phi_loeps,"loeps")
l_ratio_phi.AddEntry(G_ratio_phi_hieps,"hieps")
l_ratio_phi.Draw()
C_ratio_phi.Print(outputpdf + '(')

C_Q2_tbin = TCanvas()
C_Q2_tbin.SetGrid()
C_Q2_tbin.Divide(1,NumtBins)
l_Q2_tbin = TLegend(0.115,0.35,0.33,0.5)

multiDict = {}
for k in range(NumtBins):

    multiDict["G_Q2_tbin_{}".format(k+1)] = TMultiGraph()
    
    G_Q2_tbin_loeps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if "{:.1f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]) == "{:.1f}".format(t_bin_centers[k]):
            print("{:.1f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]), "==", "{:.1f}".format(t_bin_centers[k]))
            G_Q2_tbin_loeps.SetPoint(j, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
            j+=1
    G_Q2_tbin_loeps.SetMarkerStyle(21)
    G_Q2_tbin_loeps.SetMarkerSize(1)
    G_Q2_tbin_loeps.SetMarkerColor(1)
    multiDict["G_Q2_tbin_{}".format(k+1)].Add(G_Q2_tbin_loeps)

    G_Q2_tbin_hieps = TGraph()
    j=0
    for i in range(NumtBins*NumPhiBins):
        if "{:.1f}".format(np.array(file_df_dict['unsep_file_hieps']['tm'].tolist())[i]) == "{:.1f}".format(t_bin_centers[k]):
            print("{:.1f}".format(np.array(file_df_dict['unsep_file_hieps']['tm'].tolist())[i]), "==", "{:.1f}".format(t_bin_centers[k]))            
            G_Q2_tbin_hieps.SetPoint(j, np.array(file_df_dict['unsep_file_hieps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_hieps']['Q2'].tolist())[i])
            j+=1
    G_Q2_tbin_hieps.SetMarkerStyle(21)
    G_Q2_tbin_hieps.SetMarkerSize(1)
    G_Q2_tbin_hieps.SetMarkerColor(2)
    multiDict["G_Q2_tbin_{}".format(k+1)].Add(G_Q2_tbin_hieps)    
    
    C_Q2_tbin.cd(k+1)
    
    multiDict["G_Q2_tbin_{}".format(k+1)].Draw('AP')
    multiDict["G_Q2_tbin_{}".format(k+1)].SetTitle("t = {} ; #theta_{{cm}}; Q^2".format(t_bin_centers[k]))
    
    multiDict["G_Q2_tbin_{}".format(k+1)].GetYaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_tbin_{}".format(k+1)].GetXaxis().SetTitleOffset(1.5)
    multiDict["G_Q2_tbin_{}".format(k+1)].GetXaxis().SetLabelSize(0.04)
    
l_Q2_tbin.AddEntry(G_Q2_tbin_loeps,"loeps")
l_Q2_tbin.AddEntry(G_Q2_tbin_hieps,"hieps")
l_Q2_tbin.Draw()
C_Q2_tbin.Print(outputpdf)

C_Q2_tbin = TCanvas()
C_Q2_tbin.SetGrid()

G_Q2_tbin = TGraph()

l_Q2_tbin = TLegend(0.8,0.8,0.95,0.95)

G_Q2_tbin.SetTitle("eps = %s ; #theta_{cm}; Q^{2} Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['Q2'].tolist())):
    G_Q2_tbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
    l_Q2_tbin.AddEntry(G_Q2_tbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_Q2_tbin.SetMarkerColor(i)

G_Q2_tbin.SetMarkerStyle(21)
G_Q2_tbin.SetMarkerSize(1)
    
G_Q2_tbin.Draw('AP')
l_Q2_tbin.Draw()

C_Q2_tbin.Print(outputpdf)

C_W_tbin = TCanvas()
C_W_tbin.SetGrid()

G_W_tbin = TGraph()

l_W_tbin = TLegend(0.8,0.8,0.95,0.95)

G_W_tbin.SetTitle("eps = %s ; #theta_{cm}; W Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['W'].tolist())):
    G_W_tbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['W'].tolist())[i])
    l_W_tbin.AddEntry(G_W_tbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_W_tbin.SetMarkerColor(i)

G_W_tbin.SetMarkerStyle(21)
G_W_tbin.SetMarkerSize(1)
    
G_W_tbin.Draw('AP')
l_W_tbin.Draw()

C_W_tbin.Print(outputpdf)

C_t_tbin = TCanvas()
C_t_tbin.SetGrid()

G_t_tbin = TGraph()

l_t_tbin = TLegend(0.8,0.8,0.95,0.95)

G_t_tbin.SetTitle("eps = %s ; #theta_{cm}; t Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['tm'].tolist())):
    G_t_tbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i])
    l_t_tbin.AddEntry(G_t_tbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_t_tbin.SetMarkerColor(i)

G_t_tbin.SetMarkerStyle(21)
G_t_tbin.SetMarkerSize(1)
    
G_t_tbin.Draw('AP')
l_t_tbin.Draw()

C_t_tbin.Print(outputpdf)

C_Q2_phitbin = TCanvas()
C_Q2_phitbin.SetGrid()

G_Q2_phitbin = TGraph()

l_Q2_phitbin = TLegend(0.8,0.8,0.95,0.95)

G_Q2_phitbin.SetTitle("eps = %s ; #phi; Q^{2} Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['Q2'].tolist())):
    G_Q2_phitbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['Q2'].tolist())[i])
    l_Q2_phitbin.AddEntry(G_Q2_phitbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_Q2_phitbin.SetMarkerColor(i)

G_Q2_phitbin.SetMarkerStyle(21)
G_Q2_phitbin.SetMarkerSize(1)
    
G_Q2_phitbin.Draw('AP')
l_Q2_phitbin.Draw()

C_Q2_phitbin.Print(outputpdf)

C_W_phitbin = TCanvas()
C_W_phitbin.SetGrid()

G_W_phitbin = TGraph()

l_W_phitbin = TLegend(0.8,0.8,0.95,0.95)

G_W_phitbin.SetTitle("eps = %s ; #phi; W Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['W'].tolist())):
    G_W_phitbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['W'].tolist())[i])
    l_W_phitbin.AddEntry(G_W_phitbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_W_phitbin.SetMarkerColor(i)

G_W_phitbin.SetMarkerStyle(21)
G_W_phitbin.SetMarkerSize(1)
    
G_W_phitbin.Draw('AP')
l_W_phitbin.Draw()

C_W_phitbin.Print(outputpdf)

C_t_phitbin = TCanvas()
C_t_phitbin.SetGrid()

G_t_phitbin = TGraph()

l_t_phitbin = TLegend(0.8,0.8,0.95,0.95)

G_t_phitbin.SetTitle("eps = %s ; #phi; t Average" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['tm'].tolist())):
    G_t_phitbin.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i])
    l_t_phitbin.AddEntry(G_t_phitbin, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_t_phitbin.SetMarkerColor(i)

G_t_phitbin.SetMarkerStyle(21)
G_t_phitbin.SetMarkerSize(1)
    
G_t_phitbin.Draw('AP')
l_t_phitbin.Draw()

C_t_phitbin.Print(outputpdf)

C_xreal_th = TCanvas()
C_xreal_th.SetGrid()

G_xreal_th = TGraphErrors()

l_xreal_th = TLegend(0.8,0.8,0.95,0.95)

G_xreal_th.SetTitle("eps = %s ; #theta_{cm}; x_real" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_real'].tolist())):
    G_xreal_th.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
    G_xreal_th.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    l_xreal_th.AddEntry(G_xreal_th, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_xreal_th.SetMarkerColor(i)

G_xreal_th.SetMarkerStyle(21)
G_xreal_th.SetMarkerSize(1)
    
G_xreal_th.Draw('AP')
l_xreal_th.Draw()

C_xreal_th.Print(outputpdf)

C_xmod_th = TCanvas()
C_xmod_th.SetGrid()

G_xmod_th = TGraph()

l_xmod_th = TLegend(0.8,0.8,0.95,0.95)

G_xmod_th.SetTitle("eps = %s ; #theta_{cm}; x_mod" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_mod'].tolist())):
    G_xmod_th.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    l_xmod_th.AddEntry(G_xmod_th, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_xmod_th.SetMarkerColor(i)

G_xmod_th.SetMarkerStyle(21)
G_xmod_th.SetMarkerSize(1)
    
G_xmod_th.Draw('AP')
l_xmod_th.Draw()

C_xmod_th.Print(outputpdf)

C_xreal_phi = TCanvas()
C_xreal_phi.SetGrid()

G_xreal_phi = TGraphErrors()

l_xreal_phi = TLegend(0.8,0.8,0.95,0.95)

G_xreal_phi.SetTitle("eps = %s ; #phi; x_real" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_real'].tolist())):
    G_xreal_phi.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
    G_xreal_phi.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    l_xreal_phi.AddEntry(G_xreal_phi, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_xreal_phi.SetMarkerColor(i)

G_xreal_phi.SetMarkerStyle(21)
G_xreal_phi.SetMarkerSize(1)
    
G_xreal_phi.Draw('AP')
l_xreal_phi.Draw()

C_xreal_phi.Print(outputpdf)

C_xmod_phi = TCanvas()
C_xmod_phi.SetGrid()

G_xmod_phi = TGraph()

l_xmod_phi = TLegend(0.8,0.8,0.95,0.95)

G_xmod_phi.SetTitle("eps = %s ; #phi; x_mod" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_mod'].tolist())):
    G_xmod_phi.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    l_xmod_phi.AddEntry(G_xmod_phi, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_xmod_phi.SetMarkerColor(i)

G_xmod_phi.SetMarkerStyle(21)
G_xmod_phi.SetMarkerSize(1)
    
G_xmod_phi.Draw('AP')
l_xmod_phi.Draw()

C_xmod_phi.Print(outputpdf)

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

C_siglt_t.Print(outputpdf)

C_xreal_t = TCanvas()
C_xreal_t.SetGrid()

G_xreal_t = TGraphErrors()

l_xreal_t = TLegend(0.8,0.8,0.95,0.95)

G_xreal_t.SetTitle("eps = %s ; t; x_real" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_real'].tolist())):
    G_xreal_t.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_real'].tolist())[i])
    G_xreal_t.SetPointError(i, 0, np.array(file_df_dict['unsep_file_loeps']['dx_real'].tolist())[i])
    l_xreal_t.AddEntry(G_xreal_t, "t = {:.4f}".format(np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i]))
    G_xreal_t.SetMarkerColor(i)

G_xreal_t.SetMarkerStyle(21)
G_xreal_t.SetMarkerSize(1)
    
G_xreal_t.Draw('AP')
l_xreal_t.Draw()

C_xreal_t.Print(outputpdf)

C_xmod_t = TCanvas()
C_xmod_t.SetGrid()

G_xmod_t = TGraph()

G_xmod_t.SetTitle("eps = %s ; t; x_mod" % LOEPS)

for i in range(len(file_df_dict['unsep_file_loeps']['x_mod'].tolist())):
    G_xmod_t.SetPoint(i, np.array(file_df_dict['unsep_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['unsep_file_loeps']['x_mod'].tolist())[i])
    G_xmod_t.SetMarkerColor(i)

G_xmod_t.SetMarkerStyle(21)
G_xmod_t.SetMarkerSize(1)
    
G_xmod_t.Draw('AP')

C_xmod_t.Print(outputpdf + ')')
