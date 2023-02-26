#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-26 15:35:32 trottar"
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

    df = pd.read_csv(f_name, header=None, sep=' ')
    df.columns = columns
    return df

################################################################################################################################################
# Read in files and convert to dataframes

file_df_dict = {}

setting_file = LTANAPATH+"/src/list.settings"
file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins', 'Kset'])

for i,row in file_df_dict['setting_df'].iterrows():
    if row['Q2'] == Q2:
        file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/beam/Eb_KLT.dat", ['ebeam', 'Q2', 'EPSVAL'])        
        if row['EPSVAL'] == LOEPS:
            if row['thpq'] < 0.0:
                file_df_dict['aver_loeps_{}'.format('right')] = file_to_df( \
                                                                            LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                            .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                            , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/kindata.{}_{}_{:.0f}_-{}.dat" \
                                                                               .format(PID, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            if row['thpq'] > 0.0:
                file_df_dict['aver_loeps_{}'.format('left')] = file_to_df( \
                                                                           LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                           .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                           , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(PID, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            if row['thpq'] == 0.0:
                file_df_dict['aver_loeps_{}'.format('center')] = file_to_df( \
                                                                             LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                             .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                                             , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                                .format(PID, Q2.replace("p",""), float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                                , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            file_df_dict['xsects_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}" \
                                                            .format(PID, Q2.replace("p",""), float(LOEPS)*100) \
                                                            , [])

        if row['EPSVAL'] == HIEPS:
            if row['thpq'] < 0.0:
                file_df_dict['aver_hieps_{}'.format('right')] = file_to_df( \
                                                                            LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                            .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                            , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/kindata.{}_{}_{:.0f}_-{}.dat" \
                                                                               .format(PID, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            if row['thpq'] > 0.0:
                file_df_dict['aver_hieps_{}'.format('left')] = file_to_df( \
                                                                           LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                           .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                           , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                              .format(PID, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            if row['thpq'] == 0.0:
                file_df_dict['aver_hieps_{}'.format('center')] = file_to_df( \
                                                                             LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat" \
                                                                             .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                                             , ['relyield', 'relerr', 'phibin', 'tbin'])
                file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                                .format(PID, Q2.replace("p",""), float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                                , ['averQ2', 'averQ2err', 'averW', 'averWerr', 'avert', 'averterr'])
            file_df_dict['xsects_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}" \
                                                            .format(PID, Q2.replace("p",""), float(HIEPS)*100) \
                                                            , [])                    

avek_file = LTANAPATH+"/src/averages/aver.{}.dat".format(Q2.replace("p",""))
