#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-27 16:27:02 trottar"
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
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce
import re
import sys, math, os, subprocess

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=6:
    print("!!!!! ERROR !!!!!\n Expected 6 arguments\n Usage is with - Q2 W LOEPS HIEPS KIN OutUnsepxsectsFilename\n!!!!! ERROR !!!!!")
    sys.exit(1)

###############################################################################################################################################
    
Q2 = sys.argv[1]
W = sys.argv[2]
LOEPS = sys.argv[3]
HIEPS = sys.argv[4]

kinematics = sys.argv[5]
OutFilename = sys.argv[6]

particle = "kaon"

if particle == "kaon":
    PID = "k"
elif particle == "pion":
    PID = "pi"
    
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
        df = pd.read_csv(f_name, header=None, sep=' ', skiprows=1, skipfooter=1)
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
fix_spacing(LTANAPATH+"/src/averages/avek.{}.dat".format(Q2.replace("p","")))
fix_spacing(LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}".format(PID, Q2.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}".format(PID, Q2.replace("p",""), float(HIEPS)*100))
################################################################################################################################################
# Read in files and convert to dataframes

file_df_dict = {}

setting_file = LTANAPATH+"/src/list.settings"
file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins', 'Kset'])

for i,row in file_df_dict['setting_df'].iterrows():
    if row['Q2'] == float(Q2.replace("p",".")):
        file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/beam/Eb_KLT.dat", ['ebeam', 'Q2', 'EPSVAL'])
        file_df_dict['avek_file'] = file_to_df(LTANAPATH+"/src/averages/avek.{}.dat".format(Q2.replace("p","")) \
                                               , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])
        if row['EPSVAL'] == float(LOEPS):
            if row['thpq'] < 0.0:
                file_df_dict['aver_loeps_{}'.format('right')] = file_to_df( \
                                                                            LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                            .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                            , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_{}.dat" \
                                                                               .format(PID, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['aver_loeps_{}'.format('left')] = file_to_df( \
                                                                           LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                           .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                           , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(PID, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['aver_loeps_{}'.format('center')] = file_to_df( \
                                                                             LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                             .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                             , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_+0000.dat" \
                                                                                .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['xsects_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}" \
                                                            .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 'tm', 'um', 'um_min', 'W', 'Q2'])

        if row['EPSVAL'] == float(HIEPS):
            if row['thpq'] < 0.0:
                file_df_dict['aver_hieps_{}'.format('right')] = file_to_df( \
                                                                            LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                            .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                            , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_{}.dat" \
                                                                               .format(PID, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['aver_hieps_{}'.format('left')] = file_to_df( \
                                                                           LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                           .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                           , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(PID, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['aver_hieps_{}'.format('center')] = file_to_df( \
                                                                             LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                             .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                             , ['ratio', 'ratioerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/kindata/kindata.{}_{}_{:.0f}_+0000.dat" \
                                                                                .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['xsects_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}" \
                                                            .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 'tm', 'um', 'um_min', 'W', 'Q2'])

################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
################################################################################################################################################

G_Q2_tbin = ROOT.TGraphErrors()
C_Q2_tbin = TCanvas()
C_Q2_tbin.SetGrid()
l_Q2_tbin = ROOT.TLegend(0.115,0.35,0.33,0.5)

for i in range(len(file_df_dict['avek_file']['Q2'].tolist())):
    G_Q2_tbin.SetPoint(np.array(file_df_dict['avek_file']['th_pos'].tolist())[i], np.array(file_df_dict['avek_file']['Q2'].tolist())[i])
    G_Q2_tbin.SetPointErrors(np.array(len(file_df_dict['avek_file']['th_pos'].tolist())*[0])[i], np.array(file_df_dict['avek_file']['dQ2'].tolist())[i])
    G_Q2_tbin.SetMarkerColor(i+1)

G_Q2_tbin.SetMarkerStyle(21)
G_Q2_tbin.Draw('AP')

G_Q2_tbin.SetTitle(" ; #theta; Q^{2}")

for i,t in enumerate(file_df_dict['avek_file']['t'].tolist()):
    l_Q2_tbin.AddEntry(G_Q2_tbin, "t = {}".format(t))
l_Q2_tbin.Draw()

C_Q2_tbin.Print(outputpdf)
