#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-10-06 11:37:58 trottar"
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
import matplotlib.pyplot as plt
import sys, math, os, subprocess
import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

################################################################################################################################################
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

################################################################################################################################################

def rand_sub(phi_setting, inpDict):    

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]
    data_charge_right = inpDict["data_charge_right"] 
    data_charge_left = inpDict["data_charge_left"] 
    data_charge_center = inpDict["data_charge_center"] 
    dummy_charge_right = inpDict["dummy_charge_right"] 
    dummy_charge_left = inpDict["dummy_charge_left"] 
    dummy_charge_center = inpDict["dummy_charge_center"] 
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"] 
    efficiency_table = inpDict["efficiency_table"] 
    ParticleType = inpDict["ParticleType"]

    # Define diamond cut parameters
    a1 = inpDict["a1"]
    b1 = inpDict["b1"]
    a2 = inpDict["a2"]
    b2 = inpDict["b2"]
    a3 = inpDict["a3"]
    b3 = inpDict["b3"]
    a4 = inpDict["a4"]
    b4 = inpDict["b4"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = OUTPATH + "/" + "{}".format(ParticleType) + "_" + InDATAFilename + "_%s.root" % (phi_setting)
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        histDict.update({ "phi_setting" : phi_setting})
        return histDict

    InFile_DATA = TFile.Open(rootFileData, "OPEN")

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = OUTPATH + "/" + "{}".format(ParticleType) + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        return histDict

    InFile_DUMMY = TFile.Open(rootFileDummy, "OPEN")  

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize()))

    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize()))

    ################################################################################################################################################
    # Grabs PID cut string

    if phi_setting == "Right":
        runNums= runNumRight
        for run in runNumRight.split(' '):
            runNum = run
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,runNum)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue

        InData_efficiency = InData_efficiency_right
    if phi_setting == "Left":
        runNums= runNumLeft
        for run in runNumLeft.split(' '):
            runNum = run
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,runNum)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = InData_efficiency_left
    if phi_setting == "Center":
        runNums= runNumCenter
        for run in runNumCenter.split(' '):
            runNum = run
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,runNum)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = InData_efficiency_center

    if 'pid_text' in locals():
        print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
    else:
        print("ERROR: Invalid {} log file {}!".format(phi_setting.lower(),pid_log))
        pid_text = "\nNo {} cuts file found in logs...".format(phi_setting.lower())

    ###############################################################################################################################################
    # Grab windows for random subtraction

    # Section for grabing Prompt/Random selection parameters from PARAM file
    PARAMPATH = "%s/DB/PARAM" % UTILPATH
    print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, LTANAPATH))
    TimingCutFile = "%s/Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!
    TimingCutf = open(TimingCutFile)
    try:
        TimingCutFile
    except NameError:
        print("!!!!! ERRROR !!!!!\n One (or more) of the cut files not found!\n!!!!! ERRORR !!!!!")
        sys.exit(2)
    print("Reading timing cuts from %s" % TimingCutFile)
    PromptWindow = [0, 0]
    RandomWindows = [0, 0, 0, 0]
    linenum = 0 # Count line number we're on
    TempPar = -1 # To check later
    for line in TimingCutf: # Read all lines in the cut file
        linenum += 1 # Add one to line number at start of loop
        if(linenum > 1): # Skip first line
            line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
            line = line.rstrip()
            array = line.split(",") # Convert line into an array, anything after a comma is a new entry 
            if(int(runNum) in range (int(array[0]), int(array[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
                TempPar += 2 # If run number is in range, set to non -1 value
                BunchSpacing = float(array[2])
                CoinOffset = float(array[3]) # Coin offset value
                nSkip = float(array[4]) # Number of random windows skipped 
                nWindows = float(array[5]) # Total number of random windows
                PromptPeak = float(array[6]) # Pion CT prompt peak positon 
    TimingCutf.close() # After scanning all lines in file, close file

    if(TempPar == -1): # If value is still -1, run number provided din't match any ranges specified so exit 
        print("!!!!! ERROR !!!!!\n Run number specified does not fall within a set of runs for which cuts are defined in %s\n!!!!! ERROR !!!!!" % TimingCutFile)
        sys.exit(3)
    elif(TempPar > 1):
        print("!!! WARNING!!! Run number was found within the range of two (or more) line entries of %s !!! WARNING !!!" % TimingCutFile)
        print("The last matching entry will be treated as the input, you should ensure this is what you want")

    # From our values from the file, reconstruct our windows 
    PromptWindow[0] = PromptPeak - (BunchSpacing/2) - CoinOffset
    PromptWindow[1] = PromptPeak + (BunchSpacing/2) + CoinOffset
    RandomWindows[0] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing) - ((nWindows/2)*BunchSpacing)
    RandomWindows[1] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing)
    RandomWindows[2] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing)
    RandomWindows[3] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) + ((nWindows/2)*BunchSpacing)

    ################################################################################################################################################
    # Plot definitions

    H_hsdelta_DATA  = TH1D("H_hsdelta_DATA","HMS Delta", 500, -20.0, 20.0)
    H_hsxptar_DATA  = TH1D("H_hsxptar_DATA","HMS xptar", 500, -0.1, 0.1)
    H_hsyptar_DATA  = TH1D("H_hsyptar_DATA","HMS yptar", 500, -0.1, 0.1)
    H_ssxfp_DATA    = TH1D("H_ssxfp_DATA","SHMS xfp", 500, -25.0, 25.0)
    H_ssyfp_DATA    = TH1D("H_ssyfp_DATA","SHMS yfp", 500, -25.0, 25.0)
    H_ssxpfp_DATA   = TH1D("H_ssxpfp_DATA","SHMS xpfp", 500, -0.09, 0.09)
    H_ssypfp_DATA   = TH1D("H_ssypfp_DATA","SHMS ypfp", 500, -0.05, 0.04)
    H_hsxfp_DATA    = TH1D("H_hsxfp_DATA","HMS xfp", 500, -40.0, 40.0)
    H_hsyfp_DATA    = TH1D("H_hsyfp_DATA","HMS yfp", 500, -20.0, 20.0)
    H_hsxpfp_DATA   = TH1D("H_hsxpfp_DATA","HMS xpfp", 500, -0.09, 0.05)
    H_hsypfp_DATA   = TH1D("H_hsypfp_DATA","HMS ypfp", 500, -0.05, 0.04)
    H_ssdelta_DATA  = TH1D("H_ssdelta_DATA","SHMS delta", 500, -20.0, 20.0)
    H_ssxptar_DATA  = TH1D("H_ssxptar_DATA","SHMS xptar", 500, -0.1, 0.1)
    H_ssyptar_DATA  = TH1D("H_ssyptar_DATA","SHMS yptar", 500, -0.04, 0.04)
    H_q_DATA        = TH1D("H_q_DATA","q", 500, 0.0, 10.0)
    H_Q2_DATA       = TH1D("H_Q2_DATA","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DATA  = TH1D("H_W_DATA","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DATA       = TH1D("H_t_DATA","-t", 500, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DATA  = TH1D("H_epsilon_DATA","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DATA  = TH1D("H_MM_DATA","MM_{K}", 500, 0.0, 1.5)
    H_th_DATA  = TH1D("H_th_DATA","X' tar", 500, -0.1, 0.1)
    H_ph_DATA  = TH1D("H_ph_DATA","Y' tar", 500, -0.1, 0.1)
    H_ph_q_DATA  = TH1D("H_ph_q_DATA","Phi Detected (ph_xq)", 500, 0.0, 2*math.pi)
    H_th_q_DATA  = TH1D("H_th_q_DATA","Theta Detected (th_xq)", 500, -0.2, 0.2)
    H_ph_recoil_DATA  = TH1D("H_ph_recoil_DATA","Phi Recoil (ph_bq)", 500, -10.0, 10.0)
    H_th_recoil_DATA  = TH1D("H_th_recoil_DATA","Theta Recoil (th_bq)", 500, -10.0, 10.0)
    H_pmiss_DATA  = TH1D("H_pmiss_DATA","pmiss", 500, 0.0, 10.0)
    H_emiss_DATA  = TH1D("H_emiss_DATA","emiss", 500, 0.0, 10.0)
    H_pmx_DATA  = TH1D("H_pmx_DATA","pmx", 500, -10.0, 10.0)
    H_pmy_DATA  = TH1D("H_pmy_DATA","pmy ", 500, -10.0, 10.0)
    H_pmz_DATA  = TH1D("H_pmz_DATA","pmz", 500, -10.0, 10.0)
    H_ct_DATA = TH1D("H_ct_DATA", "Electron-{} CTime".format(ParticleType.capitalize()), 500, -50, 50)
    H_cal_etottracknorm_DATA = TH1D("H_cal_etottracknorm_DATA", "HMS Cal etottracknorm", 500, 0.2, 1.8)
    H_cer_npeSum_DATA = TH1D("H_cer_npeSum_DATA", "HMS Cer Npe Sum", 500, 0, 30)
    P_cal_etottracknorm_DATA = TH1D("P_cal_etottracknorm_DATA", "SHMS Cal etottracknorm", 500, 0, 1)
    P_hgcer_npeSum_DATA = TH1D("P_hgcer_npeSum_DATA", "SHMS HGCer Npe Sum", 500, 0, 50)
    P_aero_npeSum_DATA = TH1D("P_aero_npeSum_DATA", "SHMS Aero Npe Sum", 500, 0, 50)

    H_hsdelta_DUMMY  = TH1D("H_hsdelta_DUMMY","HMS Delta", 500, -20.0, 20.0)
    H_hsxptar_DUMMY  = TH1D("H_hsxptar_DUMMY","HMS xptar", 500, -0.1, 0.1)
    H_hsyptar_DUMMY  = TH1D("H_hsyptar_DUMMY","HMS yptar", 500, -0.1, 0.1)
    H_ssxfp_DUMMY    = TH1D("H_ssxfp_DUMMY","SHMS xfp", 500, -25.0, 25.0)
    H_ssyfp_DUMMY    = TH1D("H_ssyfp_DUMMY","SHMS yfp", 500, -25.0, 25.0)
    H_ssxpfp_DUMMY   = TH1D("H_ssxpfp_DUMMY","SHMS xpfp", 500, -0.09, 0.09)
    H_ssypfp_DUMMY   = TH1D("H_ssypfp_DUMMY","SHMS ypfp", 500, -0.05, 0.04)
    H_hsxfp_DUMMY    = TH1D("H_hsxfp_DUMMY","HMS xfp", 500, -40.0, 40.0)
    H_hsyfp_DUMMY    = TH1D("H_hsyfp_DUMMY","HMS yfp", 500, -20.0, 20.0)
    H_hsxpfp_DUMMY   = TH1D("H_hsxpfp_DUMMY","HMS xpfp", 500, -0.09, 0.05)
    H_hsypfp_DUMMY   = TH1D("H_hsypfp_DUMMY","HMS ypfp", 500, -0.05, 0.04)
    H_ssdelta_DUMMY  = TH1D("H_ssdelta_DUMMY","SHMS delta", 500, -20.0, 20.0)
    H_ssxptar_DUMMY  = TH1D("H_ssxptar_DUMMY","SHMS xptar", 500, -0.1, 0.1)
    H_ssyptar_DUMMY  = TH1D("H_ssyptar_DUMMY","SHMS yptar", 500, -0.04, 0.04)
    H_q_DUMMY        = TH1D("H_q_DUMMY","q", 500, 0.0, 10.0)
    H_Q2_DUMMY       = TH1D("H_Q2_DUMMY","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY  = TH1D("H_W_DUMMY","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY       = TH1D("H_t_DUMMY","-t", 500, inpDict["tmin"], inpDict["tmax"])  
    H_epsilon_DUMMY  = TH1D("H_epsilon_DUMMY","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY  = TH1D("H_MM_DUMMY","MM_{K}", 500, 0.0, 1.5)
    H_th_DUMMY  = TH1D("H_th_DUMMY","X' tar", 500, -0.1, 0.1)
    H_ph_DUMMY  = TH1D("H_ph_DUMMY","Y' tar", 500, -0.1, 0.1)
    H_ph_q_DUMMY  = TH1D("H_ph_q_DUMMY","Phi Detected (ph_xq)", 500, 0.0, 2*math.pi)
    H_th_q_DUMMY  = TH1D("H_th_q_DUMMY","Theta Detected (th_xq)", 500, -0.2, 0.2)
    H_ph_recoil_DUMMY  = TH1D("H_ph_recoil_DUMMY","Phi Recoil (ph_bq)", 500, -10.0, 10.0)
    H_th_recoil_DUMMY  = TH1D("H_th_recoil_DUMMY","Theta Recoil (th_bq)", 500, -10.0, 10.0)
    H_pmiss_DUMMY  = TH1D("H_pmiss_DUMMY","pmiss", 500, 0.0, 10.0)
    H_emiss_DUMMY  = TH1D("H_emiss_DUMMY","emiss", 500, 0.0, 10.0)
    H_pmx_DUMMY  = TH1D("H_pmx_DUMMY","pmx", 500, -10.0, 10.0)
    H_pmy_DUMMY  = TH1D("H_pmy_DUMMY","pmy ", 500, -10.0, 10.0)
    H_pmz_DUMMY  = TH1D("H_pmz_DUMMY","pmz", 500, -10.0, 10.0)
    H_ct_DUMMY = TH1D("H_ct_DUMMY", "Electron-{} CTime".format(ParticleType.capitalize()), 500, -50, 50)

    H_hsdelta_RAND  = TH1D("H_hsdelta_RAND","HMS Delta", 500, -20.0, 20.0)
    H_hsxptar_RAND  = TH1D("H_hsxptar_RAND","HMS xptar", 500, -0.1, 0.1)
    H_hsyptar_RAND  = TH1D("H_hsyptar_RAND","HMS yptar", 500, -0.1, 0.1)
    H_ssxfp_RAND    = TH1D("H_ssxfp_RAND","SHMS xfp", 500, -25.0, 25.0)
    H_ssyfp_RAND    = TH1D("H_ssyfp_RAND","SHMS yfp", 500, -25.0, 25.0)
    H_ssxpfp_RAND   = TH1D("H_ssxpfp_RAND","SHMS xpfp", 500, -0.09, 0.09)
    H_ssypfp_RAND   = TH1D("H_ssypfp_RAND","SHMS ypfp", 500, -0.05, 0.04)
    H_hsxfp_RAND    = TH1D("H_hsxfp_RAND","HMS xfp", 500, -40.0, 40.0)
    H_hsyfp_RAND    = TH1D("H_hsyfp_RAND","HMS yfp", 500, -20.0, 20.0)
    H_hsxpfp_RAND   = TH1D("H_hsxpfp_RAND","HMS xpfp", 500, -0.09, 0.05)
    H_hsypfp_RAND   = TH1D("H_hsypfp_RAND","HMS ypfp", 500, -0.05, 0.04)
    H_ssdelta_RAND  = TH1D("H_ssdelta_RAND","SHMS delta", 500, -20.0, 20.0)
    H_ssxptar_RAND  = TH1D("H_ssxptar_RAND","SHMS xptar", 500, -0.1, 0.1)
    H_ssyptar_RAND  = TH1D("H_ssyptar_RAND","SHMS yptar", 500, -0.04, 0.04)
    H_q_RAND        = TH1D("H_q_RAND","q", 500, 0.0, 10.0)
    H_Q2_RAND       = TH1D("H_Q2_RAND","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_RAND  = TH1D("H_W_RAND","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_RAND       = TH1D("H_t_RAND","-t", 500, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_RAND  = TH1D("H_epsilon_RAND","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_RAND  = TH1D("H_MM_RAND","MM_{K}", 500, 0.0, 1.5)
    H_th_RAND  = TH1D("H_th_RAND","X' tar", 500, -0.1, 0.1)
    H_ph_RAND  = TH1D("H_ph_RAND","Y' tar", 500, -0.1, 0.1)
    H_ph_q_RAND  = TH1D("H_ph_q_RAND","Phi Detected (ph_xq)", 500, 0.0, 2*math.pi)
    H_th_q_RAND  = TH1D("H_th_q_RAND","Theta Detected (th_xq)", 500, -0.2, 0.2)
    H_ph_recoil_RAND  = TH1D("H_ph_recoil_RAND","Phi Recoil (ph_bq)", 500, -10.0, 10.0)
    H_th_recoil_RAND  = TH1D("H_th_recoil_RAND","Theta Recoil (th_bq)", 500, -10.0, 10.0)
    H_pmiss_RAND  = TH1D("H_pmiss_RAND","pmiss", 500, 0.0, 10.0)
    H_emiss_RAND  = TH1D("H_emiss_RAND","emiss", 500, 0.0, 10.0)
    H_pmx_RAND  = TH1D("H_pmx_RAND","pmx", 500, -10.0, 10.0)
    H_pmy_RAND  = TH1D("H_pmy_RAND","pmy ", 500, -10.0, 10.0)
    H_pmz_RAND  = TH1D("H_pmz_RAND","pmz", 500, -10.0, 10.0)
    H_ct_RAND = TH1D("H_ct_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 500, -50, 50)

    H_hsdelta_DUMMY_RAND  = TH1D("H_hsdelta_DUMMY_RAND","HMS Delta", 500, -20.0, 20.0)
    H_hsxptar_DUMMY_RAND  = TH1D("H_hsxptar_DUMMY_RAND","HMS xptar", 500, -0.1, 0.1)
    H_hsyptar_DUMMY_RAND  = TH1D("H_hsyptar_DUMMY_RAND","HMS yptar", 500, -0.1, 0.1)
    H_ssxfp_DUMMY_RAND    = TH1D("H_ssxfp_DUMMY_RAND","SHMS xfp", 500, -25.0, 25.0)
    H_ssyfp_DUMMY_RAND    = TH1D("H_ssyfp_DUMMY_RAND","SHMS yfp", 500, -25.0, 25.0)
    H_ssxpfp_DUMMY_RAND   = TH1D("H_ssxpfp_DUMMY_RAND","SHMS xpfp", 500, -0.09, 0.09)
    H_ssypfp_DUMMY_RAND   = TH1D("H_ssypfp_DUMMY_RAND","SHMS ypfp", 500, -0.05, 0.04)
    H_hsxfp_DUMMY_RAND    = TH1D("H_hsxfp_DUMMY_RAND","HMS xfp", 500, -40.0, 40.0)
    H_hsyfp_DUMMY_RAND    = TH1D("H_hsyfp_DUMMY_RAND","HMS yfp", 500, -20.0, 20.0)
    H_hsxpfp_DUMMY_RAND   = TH1D("H_hsxpfp_DUMMY_RAND","HMS xpfp", 500, -0.09, 0.05)
    H_hsypfp_DUMMY_RAND   = TH1D("H_hsypfp_DUMMY_RAND","HMS ypfp", 500, -0.05, 0.04)
    H_ssdelta_DUMMY_RAND  = TH1D("H_ssdelta_DUMMY_RAND","SHMS delta", 500, -20.0, 20.0)
    H_ssxptar_DUMMY_RAND  = TH1D("H_ssxptar_DUMMY_RAND","SHMS xptar", 500, -0.1, 0.1)
    H_ssyptar_DUMMY_RAND  = TH1D("H_ssyptar_DUMMY_RAND","SHMS yptar", 500, -0.04, 0.04)
    H_q_DUMMY_RAND        = TH1D("H_q_DUMMY_RAND","q", 500, 0.0, 10.0)
    H_Q2_DUMMY_RAND       = TH1D("H_Q2_DUMMY_RAND","Q2", 500, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY_RAND  = TH1D("H_W_DUMMY_RAND","W ", 500, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY_RAND       = TH1D("H_t_DUMMY_RAND","-t", 500, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DUMMY_RAND  = TH1D("H_epsilon_DUMMY_RAND","epsilon", 500, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY_RAND  = TH1D("H_MM_DUMMY_RAND","MM_{K}", 500, 0.0, 1.5)
    H_th_DUMMY_RAND  = TH1D("H_th_DUMMY_RAND","X' tar", 500, -0.1, 0.1)
    H_ph_DUMMY_RAND  = TH1D("H_ph_DUMMY_RAND","Y' tar", 500, -0.1, 0.1)
    H_ph_q_DUMMY_RAND  = TH1D("H_ph_q_DUMMY_RAND","Phi Detected (ph_xq)", 500, 0.0, 2*math.pi)
    H_th_q_DUMMY_RAND  = TH1D("H_th_q_DUMMY_RAND","Theta Detected (th_xq)", 500, -0.2, 0.2)
    H_ph_recoil_DUMMY_RAND  = TH1D("H_ph_recoil_DUMMY_RAND","Phi Recoil (ph_bq)", 500, -10.0, 10.0)
    H_th_recoil_DUMMY_RAND  = TH1D("H_th_recoil_DUMMY_RAND","Theta Recoil (th_bq)", 500, -10.0, 10.0)
    H_pmiss_DUMMY_RAND  = TH1D("H_pmiss_DUMMY_RAND","pmiss", 500, 0.0, 10.0)
    H_emiss_DUMMY_RAND  = TH1D("H_emiss_DUMMY_RAND","emiss", 500, 0.0, 10.0)
    H_pmx_DUMMY_RAND  = TH1D("H_pmx_DUMMY_RAND","pmx", 500, -10.0, 10.0)
    H_pmy_DUMMY_RAND  = TH1D("H_pmy_DUMMY_RAND","pmy ", 500, -10.0, 10.0)
    H_pmz_DUMMY_RAND  = TH1D("H_pmz_DUMMY_RAND","pmz", 500, -10.0, 10.0)
    H_ct_DUMMY_RAND = TH1D("H_ct_DUMMY_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 500, -50, 50)

    ################################################################################################################################################
    # 2D histograms

    MM_vs_CoinTime_DATA = TH2D("MM_vs_CoinTime_DATA","Missing Mass vs CTime; MM; Coin_Time",500, 0, 2, 500, -50, 50)
    CoinTime_vs_beta_DATA = TH2D("CoinTime_vs_beta_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 500, -50, 50, 500, 0, 2)
    MM_vs_beta_DATA = TH2D("MM_vs_beta_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 500, 0, 2, 500, 0, 2)
    phiq_vs_t_DATA = TH2D("phiq_vs_t_DATA","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DATA = TGraphPolar()
    polar_phiq_vs_t_DATA.SetName("polar_phiq_vs_t_DATA")
    Q2_vs_W_DATA = TH2D("Q2_vs_W_DATA", "Q^{2} vs W; Q^{2}; W", 500, inpDict["Q2min"], inpDict["Q2max"], 500, inpDict["Wmin"], inpDict["Wmax"])

    MM_vs_CoinTime_DUMMY = TH2D("MM_vs_CoinTime_DUMMY","Missing Mass vs CTime; MM; Coin_Time",500, 0, 2, 500, -50, 50)
    CoinTime_vs_beta_DUMMY = TH2D("CoinTime_vs_beta_DUMMY", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 500, -50, 50, 500, 0, 2)
    MM_vs_beta_DUMMY = TH2D("MM_vs_beta_DUMMY", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 500, 0, 2, 500, 0, 2)
    phiq_vs_t_DUMMY = TH2D("phiq_vs_t_DUMMY","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DUMMY = TGraphPolar()
    polar_phiq_vs_t_DUMMY.SetName("polar_phiq_vs_t_DUMMY")
    Q2_vs_W_DUMMY = TH2D("Q2_vs_W_DUMMY", "Q^{2} vs W; Q^{2}; W", 500, inpDict["Q2min"], inpDict["Q2max"], 500, inpDict["Wmin"], inpDict["Wmax"])

    MM_vs_CoinTime_RAND = TH2D("MM_vs_CoinTime_RAND","Missing Mass vs CTime; MM; Coin_Time",500, 0, 2, 500, -50, 50)
    CoinTime_vs_beta_RAND = TH2D("CoinTime_vs_beta_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 500, -50, 50, 500, 0, 2)
    MM_vs_beta_RAND = TH2D("MM_vs_beta_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 500, 0, 2, 500, 0, 2)
    phiq_vs_t_RAND = TH2D("phiq_vs_t_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_RAND = TH2D("Q2_vs_W_RAND", "Q^{2} vs W; Q^{2}; W", 500, inpDict["Q2min"], inpDict["Q2max"], 500, inpDict["Wmin"], inpDict["Wmax"])


    MM_vs_CoinTime_DUMMY_RAND = TH2D("MM_vs_CoinTime_DUMMY_RAND","Missing Mass vs CTime; MM; Coin_Time",500, 0, 2, 500, -50, 50)
    CoinTime_vs_beta_DUMMY_RAND = TH2D("CoinTime_vs_beta_DUMMY_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 500, -50, 50, 500, 0, 2)
    MM_vs_beta_DUMMY_RAND = TH2D("MM_vs_beta_DUMMY_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 500, 0, 2, 500, 0, 2)
    phiq_vs_t_DUMMY_RAND = TH2D("phiq_vs_t_DUMMY_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_DUMMY_RAND = TH2D("Q2_vs_W_DUMMY_RAND", "Q^{2} vs W; Q^{2}; W", 500, inpDict["Q2min"], inpDict["Q2max"], 500, inpDict["Wmin"], inpDict["Wmax"])
    
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            # Defined Geomatrical cuts
            cutg = TCutG("cutg",21)
            cutg.SetVarX("P_hgcer_yAtCer")
            cutg.SetVarY("P_hgcer_xAtCer")
            cutg.SetPoint(0,-25,2)
            cutg.SetPoint(1,-2,2)
            cutg.SetPoint(2,-1,2.5)
            cutg.SetPoint(3,0,3)
            cutg.SetPoint(4,1,3)
            cutg.SetPoint(5,2,3.3)
            cutg.SetPoint(6,3,3.0)
            cutg.SetPoint(7,4,2.5)
            cutg.SetPoint(8,5,2)
            cutg.SetPoint(9,25,2)
            cutg.SetPoint(10,25,0.5)
            cutg.SetPoint(11,5,0.5)
            cutg.SetPoint(12,4,1)
            cutg.SetPoint(13,3,-1)
            cutg.SetPoint(14,2,-2)
            cutg.SetPoint(15,1,-2.3)
            cutg.SetPoint(16,0,-1.5)
            cutg.SetPoint(17,-1,-1)
            cutg.SetPoint(18,-2,0.5)
            cutg.SetPoint(19,-25,0.5)
            cutg.SetPoint(20,-25,2)

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):

          MM_vs_CoinTime_DATA.Fill(evt.MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DATA.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DATA.Fill(evt.MM,evt.P_gtr_beta)
          phiq_vs_t_DATA.Fill(evt.ph_q+math.pi, -evt.MandelT)
          Q2_vs_W_DATA.Fill(evt.Q2, evt.W)
          polar_phiq_vs_t_DATA.SetPoint(polar_phiq_vs_t_DATA.GetN(), (evt.ph_q+math.pi)*(180/math.pi), -evt.MandelT)
          
          H_ct_DATA.Fill(evt.CTime_ROC1)

          H_ssxfp_DATA.Fill(evt.ssxfp)
          H_ssyfp_DATA.Fill(evt.ssyfp)
          H_ssxpfp_DATA.Fill(evt.ssxpfp)
          H_ssypfp_DATA.Fill(evt.ssypfp)
          H_ssdelta_DATA.Fill(evt.ssdelta)
          H_ssxptar_DATA.Fill(evt.ssxptar)
          H_ssyptar_DATA.Fill(evt.ssyptar)

          H_hsxfp_DATA.Fill(evt.hsxfp)
          H_hsyfp_DATA.Fill(evt.hsyfp)
          H_hsxpfp_DATA.Fill(evt.hsxpfp)
          H_hsypfp_DATA.Fill(evt.hsypfp)
          H_hsdelta_DATA.Fill(evt.hsdelta)
          H_hsxptar_DATA.Fill(evt.hsxptar)	
          H_hsyptar_DATA.Fill(evt.hsyptar)

          H_ph_q_DATA.Fill(evt.ph_q+math.pi)
          H_th_q_DATA.Fill(evt.th_q)
          H_ph_recoil_DATA.Fill(evt.ph_recoil)
          H_th_recoil_DATA.Fill(evt.th_recoil)

          H_pmiss_DATA.Fill(evt.pmiss)	
          H_emiss_DATA.Fill(evt.emiss)	
          #H_emiss_DATA.Fill(evt.emiss_nuc)
          H_pmx_DATA.Fill(evt.pmx)
          H_pmy_DATA.Fill(evt.pmy)
          H_pmz_DATA.Fill(evt.pmz)
          H_Q2_DATA.Fill(evt.Q2)
          H_t_DATA.Fill(-evt.MandelT)
          H_W_DATA.Fill(evt.W)
          H_epsilon_DATA.Fill(evt.epsilon)
          H_MM_DATA.Fill(np.sqrt(abs(pow(evt.emiss, 2) - pow(evt.pmiss, 2))))
          #H_MM_DATA.Fill(pow(evt.MM, 2))  
          #H_MM_DATA.Fill(evt.Mrecoil)
          
          H_cal_etottracknorm_DATA.Fill(evt.H_cal_etottracknorm)
          H_cer_npeSum_DATA.Fill(evt.H_cer_npeSum)

          P_cal_etottracknorm_DATA.Fill(evt.P_cal_etottracknorm)
          P_hgcer_npeSum_DATA.Fill(evt.P_hgcer_npeSum)
          P_aero_npeSum_DATA.Fill(evt.P_aero_npeSum)          

    ################################################################################################################################################
    # Fill dummy histograms for various trees called above

    print("\nGrabbing {} {} dummy...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            # Defined Geomatrical cuts
            cutg = TCutG("cutg",21)
            cutg.SetVarX("P_hgcer_yAtCer")
            cutg.SetVarY("P_hgcer_xAtCer")
            cutg.SetPoint(0,-25,2)
            cutg.SetPoint(1,-2,2)
            cutg.SetPoint(2,-1,2.5)
            cutg.SetPoint(3,0,3)
            cutg.SetPoint(4,1,3)
            cutg.SetPoint(5,2,3.3)
            cutg.SetPoint(6,3,3.0)
            cutg.SetPoint(7,4,2.5)
            cutg.SetPoint(8,5,2)
            cutg.SetPoint(9,25,2)
            cutg.SetPoint(10,25,0.5)
            cutg.SetPoint(11,5,0.5)
            cutg.SetPoint(12,4,1)
            cutg.SetPoint(13,3,-1)
            cutg.SetPoint(14,2,-2)
            cutg.SetPoint(15,1,-2.3)
            cutg.SetPoint(16,0,-1.5)
            cutg.SetPoint(17,-1,-1)
            cutg.SetPoint(18,-2,0.5)
            cutg.SetPoint(19,-25,0.5)
            cutg.SetPoint(20,-25,2)

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):

          MM_vs_CoinTime_DUMMY.Fill(evt.MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DUMMY.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DUMMY.Fill(evt.MM,evt.P_gtr_beta)
          phiq_vs_t_DUMMY.Fill(evt.ph_q+math.pi, -evt.MandelT)
          Q2_vs_W_DUMMY.Fill(evt.Q2, evt.W)
          polar_phiq_vs_t_DUMMY.SetPoint(polar_phiq_vs_t_DUMMY.GetN(), (evt.ph_q+math.pi)*(180/math.pi), -evt.MandelT)            

          H_ct_DUMMY.Fill(evt.CTime_ROC1)

          H_ssxfp_DUMMY.Fill(evt.ssxfp)
          H_ssyfp_DUMMY.Fill(evt.ssyfp)
          H_ssxpfp_DUMMY.Fill(evt.ssxpfp)
          H_ssypfp_DUMMY.Fill(evt.ssypfp)
          H_ssdelta_DUMMY.Fill(evt.ssdelta)
          H_ssxptar_DUMMY.Fill(evt.ssxptar)
          H_ssyptar_DUMMY.Fill(evt.ssyptar)

          H_hsxfp_DUMMY.Fill(evt.hsxfp)
          H_hsyfp_DUMMY.Fill(evt.hsyfp)
          H_hsxpfp_DUMMY.Fill(evt.hsxpfp)
          H_hsypfp_DUMMY.Fill(evt.hsypfp)
          H_hsdelta_DUMMY.Fill(evt.hsdelta)
          H_hsxptar_DUMMY.Fill(evt.hsxptar)	
          H_hsyptar_DUMMY.Fill(evt.hsyptar)

          H_ph_q_DUMMY.Fill(evt.ph_q+math.pi)
          H_th_q_DUMMY.Fill(evt.th_q)
          H_ph_recoil_DUMMY.Fill(evt.ph_recoil)
          H_th_recoil_DUMMY.Fill(evt.th_recoil)

          H_pmiss_DUMMY.Fill(evt.pmiss)	
          H_emiss_DUMMY.Fill(evt.emiss)	
          #H_emiss_DUMMY.Fill(evt.emiss_nuc)
          H_pmx_DUMMY.Fill(evt.pmx)
          H_pmy_DUMMY.Fill(evt.pmy)
          H_pmz_DUMMY.Fill(evt.pmz)
          H_Q2_DUMMY.Fill(evt.Q2)
          H_t_DUMMY.Fill(-evt.MandelT)
          H_W_DUMMY.Fill(evt.W)
          H_epsilon_DUMMY.Fill(evt.epsilon)
          H_MM_DUMMY.Fill(np.sqrt(abs(pow(evt.emiss, 2) - pow(evt.pmiss, 2))))
          #H_MM_DUMMY.Fill(pow(evt.MM, 2))  
          #H_MM_DUMMY.Fill(evt.Mrecoil)

    ###################################################################################################################################################    
    # Fill random histograms for various trees called above

    print("\nGrabbing {} {} random data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            # Defined Geomatrical cuts
            cutg = TCutG("cutg",21)
            cutg.SetVarX("P_hgcer_yAtCer")
            cutg.SetVarY("P_hgcer_xAtCer")
            cutg.SetPoint(0,-25,2)
            cutg.SetPoint(1,-2,2)
            cutg.SetPoint(2,-1,2.5)
            cutg.SetPoint(3,0,3)
            cutg.SetPoint(4,1,3)
            cutg.SetPoint(5,2,3.3)
            cutg.SetPoint(6,3,3.0)
            cutg.SetPoint(7,4,2.5)
            cutg.SetPoint(8,5,2)
            cutg.SetPoint(9,25,2)
            cutg.SetPoint(10,25,0.5)
            cutg.SetPoint(11,5,0.5)
            cutg.SetPoint(12,4,1)
            cutg.SetPoint(13,3,-1)
            cutg.SetPoint(14,2,-2)
            cutg.SetPoint(15,1,-2.3)
            cutg.SetPoint(16,0,-1.5)
            cutg.SetPoint(17,-1,-1)
            cutg.SetPoint(18,-2,0.5)
            cutg.SetPoint(19,-25,0.5)
            cutg.SetPoint(20,-25,2)

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):

          MM_vs_CoinTime_RAND.Fill(evt.MM, evt.CTime_ROC1)
          CoinTime_vs_beta_RAND.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_RAND.Fill(evt.MM,evt.P_gtr_beta)
          phiq_vs_t_RAND.Fill(evt.ph_q+math.pi, -evt.MandelT)
          Q2_vs_W_RAND.Fill(evt.Q2, evt.W)

          H_ct_RAND.Fill(evt.CTime_ROC1)          
          
          H_ssxfp_RAND.Fill(evt.ssxfp)
          H_ssyfp_RAND.Fill(evt.ssyfp)
          H_ssxpfp_RAND.Fill(evt.ssxpfp)
          H_ssypfp_RAND.Fill(evt.ssypfp)
          H_ssdelta_RAND.Fill(evt.ssdelta)
          H_ssxptar_RAND.Fill(evt.ssxptar)
          H_ssyptar_RAND.Fill(evt.ssyptar)

          H_hsxfp_RAND.Fill(evt.hsxfp)
          H_hsyfp_RAND.Fill(evt.hsyfp)
          H_hsxpfp_RAND.Fill(evt.hsxpfp)
          H_hsypfp_RAND.Fill(evt.hsypfp)
          H_hsdelta_RAND.Fill(evt.hsdelta)
          H_hsxptar_RAND.Fill(evt.hsxptar)	
          H_hsyptar_RAND.Fill(evt.hsyptar)

          H_pmiss_RAND.Fill(evt.pmiss)	
          H_emiss_RAND.Fill(evt.emiss)	
          #H_emiss_RAND.Fill(evt.emiss_nuc)
          H_pmx_RAND.Fill(evt.pmx)
          H_pmy_RAND.Fill(evt.pmy)
          H_pmz_RAND.Fill(evt.pmz)
          H_Q2_RAND.Fill(evt.Q2)
          H_t_RAND.Fill(-evt.MandelT)
          H_W_RAND.Fill(evt.W)
          H_epsilon_RAND.Fill(evt.epsilon)
          H_MM_RAND.Fill(np.sqrt(abs(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))  )

    ###################################################################################################################################################    
    # Fill dummy random histograms for various trees called above

    print("\nGrabbing {} {} dummy random data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

        #CUTs Definations 
        SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
        SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

        HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
        HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

        Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)

        if ParticleType == "kaon":

            # Defined Geomatrical cuts
            cutg = TCutG("cutg",21)
            cutg.SetVarX("P_hgcer_yAtCer")
            cutg.SetVarY("P_hgcer_xAtCer")
            cutg.SetPoint(0,-25,2)
            cutg.SetPoint(1,-2,2)
            cutg.SetPoint(2,-1,2.5)
            cutg.SetPoint(3,0,3)
            cutg.SetPoint(4,1,3)
            cutg.SetPoint(5,2,3.3)
            cutg.SetPoint(6,3,3.0)
            cutg.SetPoint(7,4,2.5)
            cutg.SetPoint(8,5,2)
            cutg.SetPoint(9,25,2)
            cutg.SetPoint(10,25,0.5)
            cutg.SetPoint(11,5,0.5)
            cutg.SetPoint(12,4,1)
            cutg.SetPoint(13,3,-1)
            cutg.SetPoint(14,2,-2)
            cutg.SetPoint(15,1,-2.3)
            cutg.SetPoint(16,0,-1.5)
            cutg.SetPoint(17,-1,-1)
            cutg.SetPoint(18,-2,0.5)
            cutg.SetPoint(19,-25,0.5)
            cutg.SetPoint(20,-25,2)

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond and not cutg.IsInside(evt.P_hgcer_yAtCer, evt.P_hgcer_xAtCer)

        else:

            ALLCUTS = HMS_FixCut and HMS_Acceptance and SHMS_FixCut and SHMS_Acceptance and Diamond

        if(ALLCUTS):

          MM_vs_CoinTime_DUMMY_RAND.Fill(evt.MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DUMMY_RAND.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DUMMY_RAND.Fill(evt.MM,evt.P_gtr_beta)
          phiq_vs_t_DUMMY_RAND.Fill(evt.ph_q+math.pi, -evt.MandelT)
          Q2_vs_W_DUMMY_RAND.Fill(evt.Q2, evt.W)

          H_ct_DUMMY_RAND.Fill(evt.CTime_ROC1)                   

          H_ssxfp_DUMMY_RAND.Fill(evt.ssxfp)
          H_ssyfp_DUMMY_RAND.Fill(evt.ssyfp)
          H_ssxpfp_DUMMY_RAND.Fill(evt.ssxpfp)
          H_ssypfp_DUMMY_RAND.Fill(evt.ssypfp)
          H_ssdelta_DUMMY_RAND.Fill(evt.ssdelta)
          H_ssxptar_DUMMY_RAND.Fill(evt.ssxptar)
          H_ssyptar_DUMMY_RAND.Fill(evt.ssyptar)

          H_hsxfp_DUMMY_RAND.Fill(evt.hsxfp)
          H_hsyfp_DUMMY_RAND.Fill(evt.hsyfp)
          H_hsxpfp_DUMMY_RAND.Fill(evt.hsxpfp)
          H_hsypfp_DUMMY_RAND.Fill(evt.hsypfp)
          H_hsdelta_DUMMY_RAND.Fill(evt.hsdelta)
          H_hsxptar_DUMMY_RAND.Fill(evt.hsxptar)	
          H_hsyptar_DUMMY_RAND.Fill(evt.hsyptar)

          H_pmiss_DUMMY_RAND.Fill(evt.pmiss)	
          H_emiss_DUMMY_RAND.Fill(evt.emiss)	
          #H_emiss_DUMMY_RAND.Fill(evt.emiss_nuc)
          H_pmx_DUMMY_RAND.Fill(evt.pmx)
          H_pmy_DUMMY_RAND.Fill(evt.pmy)
          H_pmz_DUMMY_RAND.Fill(evt.pmz)
          H_Q2_DUMMY_RAND.Fill(evt.Q2)
          H_t_DUMMY_RAND.Fill(-evt.MandelT)
          H_W_DUMMY_RAND.Fill(evt.W)
          H_epsilon_DUMMY_RAND.Fill(evt.epsilon)
          H_MM_DUMMY_RAND.Fill(np.sqrt(abs(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))  )

    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge

    # Data Random subtraction window
    MM_vs_CoinTime_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_RAND.Scale(1/nWindows)
    MM_vs_beta_RAND.Scale(1/nWindows)
    phiq_vs_t_RAND.Scale(1/nWindows)
    Q2_vs_W_RAND.Scale(1/nWindows)
    H_ct_RAND.Scale(1/nWindows)
    H_ssxfp_RAND.Scale(1/nWindows)
    H_ssyfp_RAND.Scale(1/nWindows)
    H_ssxpfp_RAND.Scale(1/nWindows)
    H_ssypfp_RAND.Scale(1/nWindows)
    H_hsxfp_RAND.Scale(1/nWindows)
    H_hsyfp_RAND.Scale(1/nWindows)
    H_hsxpfp_RAND.Scale(1/nWindows)
    H_hsypfp_RAND.Scale(1/nWindows)
    H_ssxptar_RAND.Scale(1/nWindows)
    H_ssyptar_RAND.Scale(1/nWindows)
    H_hsxptar_RAND.Scale(1/nWindows)
    H_hsyptar_RAND.Scale(1/nWindows)
    H_ssdelta_RAND.Scale(1/nWindows)
    H_hsdelta_RAND.Scale(1/nWindows)
    H_Q2_RAND.Scale(1/nWindows)
    H_t_RAND.Scale(1/nWindows)
    H_epsilon_RAND.Scale(1/nWindows)
    H_MM_RAND.Scale(1/nWindows)
    H_pmiss_RAND.Scale(1/nWindows)
    H_emiss_RAND.Scale(1/nWindows)
    H_pmx_RAND.Scale(1/nWindows)
    H_pmy_RAND.Scale(1/nWindows)
    H_pmz_RAND.Scale(1/nWindows)
    H_W_RAND.Scale(1/nWindows)

    # Data Dummy_Random subtraction window
    MM_vs_CoinTime_DUMMY_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    phiq_vs_t_DUMMY_RAND.Scale(1/nWindows)
    Q2_vs_W_DUMMY_RAND.Scale(1/nWindows)
    H_ssxfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssyfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssxpfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssypfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsxfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsyfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsxpfp_DUMMY_RAND.Scale(1/nWindows)
    H_hsypfp_DUMMY_RAND.Scale(1/nWindows)
    H_ssxptar_DUMMY_RAND.Scale(1/nWindows)
    H_ssyptar_DUMMY_RAND.Scale(1/nWindows)
    H_hsxptar_DUMMY_RAND.Scale(1/nWindows)
    H_hsyptar_DUMMY_RAND.Scale(1/nWindows)
    H_ssdelta_DUMMY_RAND.Scale(1/nWindows)
    H_hsdelta_DUMMY_RAND.Scale(1/nWindows)
    H_Q2_DUMMY_RAND.Scale(1/nWindows)
    H_t_DUMMY_RAND.Scale(1/nWindows)
    H_epsilon_DUMMY_RAND.Scale(1/nWindows)
    H_MM_DUMMY_RAND.Scale(1/nWindows)
    H_pmiss_DUMMY_RAND.Scale(1/nWindows)
    H_emiss_DUMMY_RAND.Scale(1/nWindows)
    H_pmx_DUMMY_RAND.Scale(1/nWindows)
    H_pmy_DUMMY_RAND.Scale(1/nWindows)
    H_pmz_DUMMY_RAND.Scale(1/nWindows)
    H_W_DUMMY_RAND.Scale(1/nWindows)
    #H_ct_DUMMY_RAND.Scale(1/nWindows)

    ###
    # Data Random subtraction
    MM_vs_CoinTime_DATA.Add(MM_vs_CoinTime_RAND,-1)
    CoinTime_vs_beta_DATA.Add(CoinTime_vs_beta_RAND,-1)
    MM_vs_beta_DATA.Add(MM_vs_beta_RAND,-1)
    phiq_vs_t_DATA.Add(phiq_vs_t_RAND,-1)
    Q2_vs_W_DATA.Add(Q2_vs_W_RAND,-1)
    H_ssxfp_DATA.Add(H_ssxfp_RAND,-1)
    H_ssyfp_DATA.Add(H_ssyfp_RAND,-1)
    H_ssxpfp_DATA.Add(H_ssxpfp_RAND,-1)
    H_ssypfp_DATA.Add(H_ssypfp_RAND,-1)
    H_hsxfp_DATA.Add(H_hsxfp_RAND,-1)
    H_hsyfp_DATA.Add(H_hsyfp_RAND,-1)
    H_hsxpfp_DATA.Add(H_hsxpfp_RAND,-1)
    H_hsypfp_DATA.Add(H_hsypfp_RAND,-1)
    H_ssxptar_DATA.Add(H_ssxptar_RAND,-1)
    H_ssyptar_DATA.Add(H_ssyptar_RAND,-1)
    H_hsxptar_DATA.Add(H_hsxptar_RAND,-1)
    H_hsyptar_DATA.Add(H_hsyptar_RAND,-1)
    H_ssdelta_DATA.Add(H_ssdelta_RAND,-1)
    H_hsdelta_DATA.Add(H_hsdelta_RAND,-1)
    H_Q2_DATA.Add(H_Q2_RAND,-1)
    H_t_DATA.Add(H_t_RAND,-1)
    H_epsilon_DATA.Add(H_epsilon_RAND,-1)
    H_MM_DATA.Add(H_MM_RAND,-1)
    H_pmiss_DATA.Add(H_pmiss_RAND,-1)
    H_emiss_DATA.Add(H_emiss_RAND,-1)
    H_pmx_DATA.Add(H_pmx_RAND,-1)
    H_pmy_DATA.Add(H_pmy_RAND,-1)
    H_pmz_DATA.Add(H_pmz_RAND,-1)
    H_W_DATA.Add(H_W_RAND,-1)
    H_ct_DATA.Add(H_ct_RAND,-1)

    ###
    # Dummy Random subtraction
    MM_vs_CoinTime_DUMMY.Add(MM_vs_CoinTime_DUMMY_RAND,-1)
    CoinTime_vs_beta_DUMMY.Add(CoinTime_vs_beta_DUMMY_RAND,-1)
    MM_vs_beta_DUMMY.Add(MM_vs_beta_DUMMY_RAND,-1)
    phiq_vs_t_DUMMY.Add(phiq_vs_t_DUMMY_RAND,-1)
    Q2_vs_W_DUMMY.Add(Q2_vs_W_DUMMY_RAND,-1)
    H_ssxfp_DUMMY.Add(H_ssxfp_DUMMY_RAND,-1)
    H_ssyfp_DUMMY.Add(H_ssyfp_DUMMY_RAND,-1)
    H_ssxpfp_DUMMY.Add(H_ssxpfp_DUMMY_RAND,-1)
    H_ssypfp_DUMMY.Add(H_ssypfp_DUMMY_RAND,-1)
    H_hsxfp_DUMMY.Add(H_hsxfp_DUMMY_RAND,-1)
    H_hsyfp_DUMMY.Add(H_hsyfp_DUMMY_RAND,-1)
    H_hsxpfp_DUMMY.Add(H_hsxpfp_DUMMY_RAND,-1)
    H_hsypfp_DUMMY.Add(H_hsypfp_DUMMY_RAND,-1)
    H_ssxptar_DUMMY.Add(H_ssxptar_DUMMY_RAND,-1)
    H_ssyptar_DUMMY.Add(H_ssyptar_DUMMY_RAND,-1)
    H_hsxptar_DUMMY.Add(H_hsxptar_DUMMY_RAND,-1)
    H_hsyptar_DUMMY.Add(H_hsyptar_DUMMY_RAND,-1)
    H_ssdelta_DUMMY.Add(H_ssdelta_DUMMY_RAND,-1)
    H_hsdelta_DUMMY.Add(H_hsdelta_DUMMY_RAND,-1)
    H_Q2_DUMMY.Add(H_Q2_DUMMY_RAND,-1)
    H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
    H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)
    H_MM_DUMMY.Add(H_MM_DUMMY_RAND,-1)
    H_pmiss_DUMMY.Add(H_pmiss_DUMMY_RAND,-1)
    H_emiss_DUMMY.Add(H_emiss_DUMMY_RAND,-1)
    H_pmx_DUMMY.Add(H_pmx_DUMMY_RAND,-1)
    H_pmy_DUMMY.Add(H_pmy_DUMMY_RAND,-1)
    H_pmz_DUMMY.Add(H_pmz_DUMMY_RAND,-1)
    H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
    H_ct_DUMMY.Add(H_ct_DUMMY_RAND,-1)

    histDict["InFile_DATA"] = InFile_DATA
    histDict["InFile_DUMMY"] = InFile_DUMMY
    histDict["phi_setting"] = phi_setting
    histDict["pid_text"] = pid_text
    histDict["runNums"] = runNums.split(' ')
    histDict["H_hsdelta_DUMMY"] =     H_hsdelta_DUMMY
    histDict["H_hsxptar_DUMMY"] =     H_hsxptar_DUMMY
    histDict["H_hsyptar_DUMMY"] =     H_hsyptar_DUMMY
    histDict["H_ssxfp_DUMMY"] =     H_ssxfp_DUMMY  
    histDict["H_ssyfp_DUMMY"] =     H_ssyfp_DUMMY  
    histDict["H_ssxpfp_DUMMY"] =     H_ssxpfp_DUMMY 
    histDict["H_ssypfp_DUMMY"] =     H_ssypfp_DUMMY 
    histDict["H_hsxfp_DUMMY"] =     H_hsxfp_DUMMY  
    histDict["H_hsyfp_DUMMY"] =     H_hsyfp_DUMMY  
    histDict["H_hsxpfp_DUMMY"] =     H_hsxpfp_DUMMY 
    histDict["H_hsypfp_DUMMY"] =     H_hsypfp_DUMMY 
    histDict["H_ssdelta_DUMMY"] =     H_ssdelta_DUMMY
    histDict["H_ssxptar_DUMMY"] =     H_ssxptar_DUMMY
    histDict["H_ssyptar_DUMMY"] =     H_ssyptar_DUMMY
    histDict["H_q_DUMMY"] =     H_q_DUMMY      
    histDict["H_Q2_DUMMY"] =     H_Q2_DUMMY     
    histDict["H_t_DUMMY"] =     H_t_DUMMY     
    histDict["H_epsilon_DUMMY"] =     H_epsilon_DUMMY
    histDict["H_MM_DUMMY"] =     H_MM_DUMMY
    histDict["H_th_DUMMY"] =     H_th_DUMMY
    histDict["H_ph_DUMMY"] =     H_ph_DUMMY
    histDict["H_ph_q_DUMMY"] =     H_ph_q_DUMMY
    histDict["H_th_q_DUMMY"] =     H_th_q_DUMMY
    histDict["H_ph_recoil_DUMMY"] =     H_ph_recoil_DUMMY
    histDict["H_th_recoil_DUMMY"] =     H_th_recoil_DUMMY
    histDict["H_pmiss_DUMMY"] =     H_pmiss_DUMMY
    histDict["H_emiss_DUMMY"] =     H_emiss_DUMMY
    histDict["H_pmx_DUMMY"] =     H_pmx_DUMMY
    histDict["H_pmy_DUMMY"] =     H_pmy_DUMMY
    histDict["H_pmz_DUMMY"] =     H_pmz_DUMMY
    histDict["H_W_DUMMY"] =     H_W_DUMMY
    histDict["H_ct_DUMMY"] =     H_ct_DUMMY
    histDict["MM_vs_CoinTime_DUMMY"] = MM_vs_CoinTime_DUMMY
    histDict["CoinTime_vs_beta_DUMMY"] = CoinTime_vs_beta_DUMMY
    histDict["MM_vs_beta_DUMMY"] = MM_vs_beta_DUMMY
    histDict["phiq_vs_t_DUMMY"] = phiq_vs_t_DUMMY
    histDict["polar_phiq_vs_t_DUMMY"] = polar_phiq_vs_t_DUMMY    
    histDict["H_hsdelta_DATA"] =     H_hsdelta_DATA
    histDict["H_hsxptar_DATA"] =     H_hsxptar_DATA
    histDict["H_hsyptar_DATA"] =     H_hsyptar_DATA
    histDict["H_ssxfp_DATA"] =     H_ssxfp_DATA  
    histDict["H_ssyfp_DATA"] =     H_ssyfp_DATA  
    histDict["H_ssxpfp_DATA"] =     H_ssxpfp_DATA 
    histDict["H_ssypfp_DATA"] =     H_ssypfp_DATA 
    histDict["H_hsxfp_DATA"] =     H_hsxfp_DATA  
    histDict["H_hsyfp_DATA"] =     H_hsyfp_DATA  
    histDict["H_hsxpfp_DATA"] =     H_hsxpfp_DATA 
    histDict["H_hsypfp_DATA"] =     H_hsypfp_DATA 
    histDict["H_ssdelta_DATA"] =     H_ssdelta_DATA
    histDict["H_ssxptar_DATA"] =     H_ssxptar_DATA
    histDict["H_ssyptar_DATA"] =     H_ssyptar_DATA
    histDict["H_q_DATA"] =     H_q_DATA      
    histDict["H_Q2_DATA"] =     H_Q2_DATA     
    histDict["H_t_DATA"] =     H_t_DATA     
    histDict["H_epsilon_DATA"] =     H_epsilon_DATA
    histDict["H_MM_DATA"] =     H_MM_DATA
    histDict["H_th_DATA"] =     H_th_DATA
    histDict["H_ph_DATA"] =     H_ph_DATA
    histDict["H_ph_q_DATA"] =     H_ph_q_DATA
    histDict["H_th_q_DATA"] =     H_th_q_DATA
    histDict["H_ph_recoil_DATA"] =     H_ph_recoil_DATA
    histDict["H_th_recoil_DATA"] =     H_th_recoil_DATA
    histDict["H_pmiss_DATA"] =     H_pmiss_DATA
    histDict["H_emiss_DATA"] =     H_emiss_DATA
    histDict["H_pmx_DATA"] =     H_pmx_DATA
    histDict["H_pmy_DATA"] =     H_pmy_DATA
    histDict["H_pmz_DATA"] =     H_pmz_DATA
    histDict["H_W_DATA"] =     H_W_DATA
    histDict["H_ct_DATA"] =     H_ct_DATA
    histDict["H_cal_etottracknorm_DATA"] =     H_cal_etottracknorm_DATA
    histDict["H_cer_npeSum_DATA"] =     H_cer_npeSum_DATA
    histDict["P_cal_etottracknorm_DATA"] =     P_cal_etottracknorm_DATA
    histDict["P_hgcer_npeSum_DATA"] =     P_hgcer_npeSum_DATA
    histDict["P_aero_npeSum_DATA"] =     P_aero_npeSum_DATA
    histDict["MM_vs_CoinTime_DATA"] = MM_vs_CoinTime_DATA
    histDict["CoinTime_vs_beta_DATA"] = CoinTime_vs_beta_DATA
    histDict["MM_vs_beta_DATA"] = MM_vs_beta_DATA
    histDict["phiq_vs_t_DATA"] = phiq_vs_t_DATA
    histDict["polar_phiq_vs_t_DATA"] = polar_phiq_vs_t_DATA
    histDict["Q2_vs_W_DATA"] = Q2_vs_W_DATA

    ###
    # CT plots
    ct = TCanvas()
    l_ct = TLegend(0.115,0.65,0.33,0.95)
    l_ct.SetTextSize(0.0235)
    H_ct_DATA.SetLineColor(2)
    H_ct_DATA.Draw("same, HIST")
    l_ct.AddEntry(H_ct_DATA,"{}".format(ParticleType.capitalize()))
    l_ct.Draw()

    ct.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+'(')

    ###
    # PID Plots
    c_pid = TCanvas()

    c_pid.Divide(2,3)

    c_pid.cd(1)
    gPad.SetLogy()

    H_cal_etottracknorm_DATA.SetLineColor(2)
    H_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(2)
    gPad.SetLogy()

    H_cer_npeSum_DATA.SetLineColor(2)
    H_cer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(3)
    gPad.SetLogy()

    P_cal_etottracknorm_DATA.SetLineColor(2)
    P_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(4)
    gPad.SetLogy()

    P_hgcer_npeSum_DATA.SetLineColor(2)
    P_hgcer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(5)
    gPad.SetLogy()

    P_aero_npeSum_DATA.SetLineColor(2)
    P_aero_npeSum_DATA.Draw("same, HIST")

    c_pid.Draw()

    c_pid.Print(outputpdf.replace("{}_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')

    return histDict
