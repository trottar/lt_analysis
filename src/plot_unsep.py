#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-26 14:24:03 trottar"
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

setting_file = LTANAPATH+"/src/list.settings"
setting_df = file_to_df(setting_file, ['POL', 'Q2', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins', 'Kset'])

thpq = setting_df['thpq']

print(thpq)

'''
beam_file = LTANAPATH+"/src/beam/Eb_KLT.dat"
aver_file_loeps = LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat".format(PID, Q2.replace("p",""), float(LOEPS)*100)
aver_file_hieps = LTANAPATH+"/src/averages/aver.{}_{}_{:.0f}.dat".format(PID, Q2.replace("p",""), float(HIEPS)*100)
avek_file = LTANAPATH+"/src/averages/aver.{}.dat".format(Q2.replace("p",""))
kindata_file_loeps_right = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_-{}.dat".format(PID, Q2.replace("p",""), float(LOEPS)*100, int(thpq_right*1000))
kindata_file_loeps_left = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat".format(PID, Q2.replace("p",""), float(LOEPS)*100, int(thpq_left*1000))
kindata_file_loeps_center = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat".format(PID, Q2.replace("p",""), float(LOEPS)*100, int(thpq_center*1000))
kindata_file_hieps_right = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_-{}.dat".format(PID, Q2.replace("p",""), float(HIEPS)*100, int(thpq_right*1000))
kindata_file_hieps_left = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat".format(PID, Q2.replace("p",""), float(HIEPS)*100, int(thpq_left*1000))
kindata_file_hieps_center = LTANAPATH+"/src/kindata.{}_{}_{:.0f}_+{}.dat".format(PID, Q2.replace("p",""), float(HIEPS)*100, int(thpq_center*1000))
xsects_file = LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}".format(PID, Q2.replace("p",""), float(LOEPS)*100)
'''
