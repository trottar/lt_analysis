#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-02 05:32:19 trottar"
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

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import open_root_file, remove_bad_bins, create_polar_plot

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
# Disable statistics box by default
#ROOT.gStyle.SetOptStat(0)
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
    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]     
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
    data_charge_err_right = inpDict["data_charge_err_right"] 
    data_charge_err_left = inpDict["data_charge_err_left"] 
    data_charge_err_center = inpDict["data_charge_err_center"] 
    dummy_charge_err_right = inpDict["dummy_charge_err_right"] 
    dummy_charge_err_left = inpDict["dummy_charge_err_left"] 
    dummy_charge_err_center = inpDict["dummy_charge_err_center"]     
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"]
    InData_error_efficiency_right = inpDict["InData_error_efficiency_right"] 
    InData_error_efficiency_left = inpDict["InData_error_efficiency_left"] 
    InData_error_efficiency_center = inpDict["InData_error_efficiency_center"]    
    efficiency_table = inpDict["efficiency_table"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define HGCer hole cut for KaonLT 2018-19
    if ParticleType == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(Q2, W, EPSSET)
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        histDict.update({ "phi_setting" : phi_setting})
        return histDict

    InFile_DATA = open_root_file(rootFileData)

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{ParticleType}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        return histDict

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_noRF".format(ParticleType.capitalize()))

    ##############
    # HARD CODED #
    ##############
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

    if ParticleType == "kaon":
        for c0, p in zip(c0_list, h_momentum_list):
            if p == 0.889:
                c0_dict["Q2p1W2p95_lowe"] = c0 # Proper value 0.888
            elif p == 0.968:
                c0_dict["Q0p5W2p40_lowe"] = c0
                c0_dict["Q3p0W3p14_lowe"] = c0 # Proper value 1.821
                c0_dict["Q5p5W3p02_lowe"] = c0 # Proper value 0.962
            elif p == 2.185:
                c0_dict["Q0p5W2p40_highe"] = c0 # Proper value 2.066
                c0_dict["Q3p0W2p32_lowe"] = c0
            elif p == 2.328:
                c0_dict["Q4p4W2p74_lowe"] = c0
            elif p == 3.266:
                c0_dict["Q5p5W3p02_highe"] = c0            
            elif p == 4.2:
                c0_dict["Q3p0W3p14_highe"] = c0 # Proper value 4.204
            elif p == 4.712:
                c0_dict["Q4p4W2p74_highe"] = c0            
            elif p == 5.292:
                c0_dict["Q2p1W2p95_highe"] = c0
            elif p == 6.59:
                c0_dict["Q3p0W2p32_highe"] = c0
    else:
        c0_dict["Q0p4W2p20_lowe"] = 0.0
        c0_dict["Q0p4W2p20_highe"] = 0.0
        
    ##############
    ##############        
    ##############
    
    ################################################################################################################################################
    # Grabs PID cut string

    if phi_setting == "Right":
        runNums= runNumRight
        for run in runNumRight.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue

        InData_efficiency = InData_efficiency_right
    if phi_setting == "Left":
        runNums= runNumLeft
        for run in runNumLeft.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = InData_efficiency_left
    if phi_setting == "Center":
        runNums= runNumCenter
        for run in runNumCenter.split(' '):
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
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
            if(int(run) in range (int(array[0]), int(array[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
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

    H_hsdelta_DATA  = TH1D("H_hsdelta_DATA","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DATA  = TH1D("H_hsxptar_DATA","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DATA  = TH1D("H_hsyptar_DATA","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DATA    = TH1D("H_ssxfp_DATA","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DATA    = TH1D("H_ssyfp_DATA","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DATA   = TH1D("H_ssxpfp_DATA","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DATA   = TH1D("H_ssypfp_DATA","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DATA    = TH1D("H_hsxfp_DATA","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DATA    = TH1D("H_hsyfp_DATA","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DATA   = TH1D("H_hsxpfp_DATA","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DATA   = TH1D("H_hsypfp_DATA","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DATA  = TH1D("H_ssdelta_DATA","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DATA  = TH1D("H_ssxptar_DATA","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DATA  = TH1D("H_ssyptar_DATA","SHMS yptar", 100, -0.04, 0.04)
    H_q_DATA        = TH1D("H_q_DATA","q", 100, 0.0, 10.0)
    H_Q2_DATA       = TH1D("H_Q2_DATA","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DATA  = TH1D("H_W_DATA","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DATA       = TH1D("H_t_DATA","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DATA  = TH1D("H_epsilon_DATA","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DATA  = TH1D("H_MM_DATA","MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_nosub_DATA  = TH1D("H_MM_nosub_DATA","MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_th_DATA  = TH1D("H_th_DATA","X' tar", 100, -0.1, 0.1)
    H_ph_DATA  = TH1D("H_ph_DATA","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DATA  = TH1D("H_ph_q_DATA","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DATA  = TH1D("H_th_q_DATA","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_DATA  = TH1D("H_ph_recoil_DATA","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_DATA  = TH1D("H_th_recoil_DATA","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_DATA  = TH1D("H_pmiss_DATA","pmiss", 100, 0.0, 2.0)
    H_emiss_DATA  = TH1D("H_emiss_DATA","emiss", 100, 0.0, 2.0)
    H_pmx_DATA  = TH1D("H_pmx_DATA","pmx", 100, -10.0, 10.0)
    H_pmy_DATA  = TH1D("H_pmy_DATA","pmy ", 100, -10.0, 10.0)
    H_pmz_DATA  = TH1D("H_pmz_DATA","pmz", 100, -10.0, 10.0)
    H_ct_DATA = TH1D("H_ct_DATA", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)
    H_cal_etottracknorm_DATA = TH1D("H_cal_etottracknorm_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
    H_cer_npeSum_DATA = TH1D("H_cer_npeSum_DATA", "HMS Cer Npe Sum", 100, 0, 30)
    P_cal_etottracknorm_DATA = TH1D("P_cal_etottracknorm_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
    P_hgcer_npeSum_DATA = TH1D("P_hgcer_npeSum_DATA", "SHMS HGCer Npe Sum", 100, 0, 10)
    P_aero_npeSum_DATA = TH1D("P_aero_npeSum_DATA", "SHMS Aero Npe Sum", 100, 0, 30)

    H_hsdelta_DUMMY  = TH1D("H_hsdelta_DUMMY","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DUMMY  = TH1D("H_hsxptar_DUMMY","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DUMMY  = TH1D("H_hsyptar_DUMMY","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DUMMY    = TH1D("H_ssxfp_DUMMY","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DUMMY    = TH1D("H_ssyfp_DUMMY","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DUMMY   = TH1D("H_ssxpfp_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DUMMY   = TH1D("H_ssypfp_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DUMMY    = TH1D("H_hsxfp_DUMMY","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DUMMY    = TH1D("H_hsyfp_DUMMY","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DUMMY   = TH1D("H_hsxpfp_DUMMY","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DUMMY   = TH1D("H_hsypfp_DUMMY","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DUMMY  = TH1D("H_ssdelta_DUMMY","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DUMMY  = TH1D("H_ssxptar_DUMMY","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DUMMY  = TH1D("H_ssyptar_DUMMY","SHMS yptar", 100, -0.04, 0.04)
    H_q_DUMMY        = TH1D("H_q_DUMMY","q", 100, 0.0, 10.0)
    H_Q2_DUMMY       = TH1D("H_Q2_DUMMY","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY  = TH1D("H_W_DUMMY","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY       = TH1D("H_t_DUMMY","-t", 100, inpDict["tmin"], inpDict["tmax"])  
    H_epsilon_DUMMY  = TH1D("H_epsilon_DUMMY","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY  = TH1D("H_MM_DUMMY","MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_nosub_DUMMY  = TH1D("H_MM_nosub_DUMMY","MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_th_DUMMY  = TH1D("H_th_DUMMY","X' tar", 100, -0.1, 0.1)
    H_ph_DUMMY  = TH1D("H_ph_DUMMY","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DUMMY  = TH1D("H_ph_q_DUMMY","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DUMMY  = TH1D("H_th_q_DUMMY","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_DUMMY  = TH1D("H_ph_recoil_DUMMY","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_DUMMY  = TH1D("H_th_recoil_DUMMY","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_DUMMY  = TH1D("H_pmiss_DUMMY","pmiss", 100, 0.0, 2.0)
    H_emiss_DUMMY  = TH1D("H_emiss_DUMMY","emiss", 100, 0.0, 2.0)
    H_pmx_DUMMY  = TH1D("H_pmx_DUMMY","pmx", 100, -10.0, 10.0)
    H_pmy_DUMMY  = TH1D("H_pmy_DUMMY","pmy ", 100, -10.0, 10.0)
    H_pmz_DUMMY  = TH1D("H_pmz_DUMMY","pmz", 100, -10.0, 10.0)
    H_ct_DUMMY = TH1D("H_ct_DUMMY", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)

    H_hsdelta_RAND  = TH1D("H_hsdelta_RAND","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_RAND  = TH1D("H_hsxptar_RAND","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_RAND  = TH1D("H_hsyptar_RAND","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_RAND    = TH1D("H_ssxfp_RAND","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_RAND    = TH1D("H_ssyfp_RAND","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_RAND   = TH1D("H_ssxpfp_RAND","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_RAND   = TH1D("H_ssypfp_RAND","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_RAND    = TH1D("H_hsxfp_RAND","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_RAND    = TH1D("H_hsyfp_RAND","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_RAND   = TH1D("H_hsxpfp_RAND","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_RAND   = TH1D("H_hsypfp_RAND","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_RAND  = TH1D("H_ssdelta_RAND","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_RAND  = TH1D("H_ssxptar_RAND","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_RAND  = TH1D("H_ssyptar_RAND","SHMS yptar", 100, -0.04, 0.04)
    H_q_RAND        = TH1D("H_q_RAND","q", 100, 0.0, 10.0)
    H_Q2_RAND       = TH1D("H_Q2_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_RAND  = TH1D("H_W_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_RAND       = TH1D("H_t_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_RAND  = TH1D("H_epsilon_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_RAND  = TH1D("H_MM_RAND","MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_nosub_RAND  = TH1D("H_MM_nosub_RAND","MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_th_RAND  = TH1D("H_th_RAND","X' tar", 100, -0.1, 0.1)
    H_ph_RAND  = TH1D("H_ph_RAND","Y' tar", 100, -0.1, 0.1)
    H_ph_q_RAND  = TH1D("H_ph_q_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_RAND  = TH1D("H_th_q_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_RAND  = TH1D("H_ph_recoil_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_RAND  = TH1D("H_th_recoil_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_RAND  = TH1D("H_pmiss_RAND","pmiss", 100, 0.0, 2.0)
    H_emiss_RAND  = TH1D("H_emiss_RAND","emiss", 100, 0.0, 2.0)
    H_pmx_RAND  = TH1D("H_pmx_RAND","pmx", 100, -10.0, 10.0)
    H_pmy_RAND  = TH1D("H_pmy_RAND","pmy ", 100, -10.0, 10.0)
    H_pmz_RAND  = TH1D("H_pmz_RAND","pmz", 100, -10.0, 10.0)
    H_ct_RAND = TH1D("H_ct_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)

    H_hsdelta_DUMMY_RAND  = TH1D("H_hsdelta_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
    H_hsxptar_DUMMY_RAND  = TH1D("H_hsxptar_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
    H_hsyptar_DUMMY_RAND  = TH1D("H_hsyptar_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
    H_ssxfp_DUMMY_RAND    = TH1D("H_ssxfp_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
    H_ssyfp_DUMMY_RAND    = TH1D("H_ssyfp_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
    H_ssxpfp_DUMMY_RAND   = TH1D("H_ssxpfp_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
    H_ssypfp_DUMMY_RAND   = TH1D("H_ssypfp_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
    H_hsxfp_DUMMY_RAND    = TH1D("H_hsxfp_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
    H_hsyfp_DUMMY_RAND    = TH1D("H_hsyfp_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
    H_hsxpfp_DUMMY_RAND   = TH1D("H_hsxpfp_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
    H_hsypfp_DUMMY_RAND   = TH1D("H_hsypfp_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
    H_ssdelta_DUMMY_RAND  = TH1D("H_ssdelta_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
    H_ssxptar_DUMMY_RAND  = TH1D("H_ssxptar_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
    H_ssyptar_DUMMY_RAND  = TH1D("H_ssyptar_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
    H_q_DUMMY_RAND        = TH1D("H_q_DUMMY_RAND","q", 100, 0.0, 10.0)
    H_Q2_DUMMY_RAND       = TH1D("H_Q2_DUMMY_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
    H_W_DUMMY_RAND  = TH1D("H_W_DUMMY_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
    H_t_DUMMY_RAND       = TH1D("H_t_DUMMY_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
    H_epsilon_DUMMY_RAND  = TH1D("H_epsilon_DUMMY_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
    H_MM_DUMMY_RAND  = TH1D("H_MM_DUMMY_RAND","MM_{ParticleType[0].upper()}", 100, inpDict["mm_min"], inpDict["mm_max"])
    H_MM_nosub_DUMMY_RAND  = TH1D("H_MM_nosub_DUMMY_RAND","MM_nosub_{ParticleType[0].upper()}", 100, 0.7, 1.5)
    H_th_DUMMY_RAND  = TH1D("H_th_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
    H_ph_DUMMY_RAND  = TH1D("H_ph_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
    H_ph_q_DUMMY_RAND  = TH1D("H_ph_q_DUMMY_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
    H_th_q_DUMMY_RAND  = TH1D("H_th_q_DUMMY_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
    H_ph_recoil_DUMMY_RAND  = TH1D("H_ph_recoil_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
    H_th_recoil_DUMMY_RAND  = TH1D("H_th_recoil_DUMMY_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
    H_pmiss_DUMMY_RAND  = TH1D("H_pmiss_DUMMY_RAND","pmiss", 100, 0.0, 2.0)
    H_emiss_DUMMY_RAND  = TH1D("H_emiss_DUMMY_RAND","emiss", 100, 0.0, 2.0)
    H_pmx_DUMMY_RAND  = TH1D("H_pmx_DUMMY_RAND","pmx", 100, -10.0, 10.0)
    H_pmy_DUMMY_RAND  = TH1D("H_pmy_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
    H_pmz_DUMMY_RAND  = TH1D("H_pmz_DUMMY_RAND","pmz", 100, -10.0, 10.0)
    H_ct_DUMMY_RAND = TH1D("H_ct_DUMMY_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)

    ################################################################################################################################################
    # 2D histograms

    MM_vs_CoinTime_DATA = TH2D("MM_vs_CoinTime_DATA","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DATA = TH2D("CoinTime_vs_beta_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DATA = TH2D("MM_vs_beta_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DATA = TH2D("MM_vs_H_cer_DATA", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_DATA = TH2D("MM_vs_H_cal_DATA", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_DATA = TH2D("MM_vs_P_cal_DATA", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
    MM_vs_P_hgcer_DATA = TH2D("MM_vs_P_hgcer_DATA", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DATA = TH2D("MM_vs_P_aero_DATA", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    phiq_vs_t_DATA = TH2D("phiq_vs_t_DATA","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DATA = TGraphPolar()
    polar_phiq_vs_t_DATA.SetName("polar_phiq_vs_t_DATA")
    Q2_vs_W_DATA = TH2D("Q2_vs_W_DATA", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DATA = TH2D("Q2_vs_t_DATA", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DATA = TH2D("W_vs_t_DATA", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DATA = TH2D("EPS_vs_t_DATA", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DATA = TH2D("MM_vs_t_DATA", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DATA = TH2D("P_hgcer_xAtCer_vs_yAtCer_DATA", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DATA", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DATA = TH2D("P_hgcer_xAtCer_vs_MM_DATA", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DATA", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DATA = TH2D("P_hgcer_yAtCer_vs_MM_DATA", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DATA = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DATA", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_DUMMY = TH2D("MM_vs_CoinTime_DUMMY","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DUMMY = TH2D("CoinTime_vs_beta_DUMMY", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DUMMY = TH2D("MM_vs_beta_DUMMY", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DUMMY = TH2D("MM_vs_H_cer_DUMMY", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_P_cal_DUMMY = TH2D("MM_vs_P_cal_DUMMY", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)    
    MM_vs_H_cal_DUMMY = TH2D("MM_vs_H_cal_DUMMY", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
    MM_vs_P_hgcer_DUMMY = TH2D("MM_vs_P_hgcer_DUMMY", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DUMMY = TH2D("MM_vs_P_aero_DUMMY", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_DUMMY = TH2D("phiq_vs_t_DUMMY","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    polar_phiq_vs_t_DUMMY = TGraphPolar()
    polar_phiq_vs_t_DUMMY.SetName("polar_phiq_vs_t_DUMMY")
    Q2_vs_W_DUMMY = TH2D("Q2_vs_W_DUMMY", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DUMMY = TH2D("Q2_vs_t_DUMMY", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DUMMY = TH2D("W_vs_t_DUMMY", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DUMMY = TH2D("EPS_vs_t_DUMMY", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DUMMY = TH2D("MM_vs_t_DUMMY", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DUMMY = TH2D("P_hgcer_xAtCer_vs_yAtCer_DUMMY", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DUMMY = TH2D("P_hgcer_xAtCer_vs_MM_DUMMY", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DUMMY", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DUMMY = TH2D("P_hgcer_yAtCer_vs_MM_DUMMY", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DUMMY", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_RAND = TH2D("MM_vs_CoinTime_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_RAND = TH2D("CoinTime_vs_beta_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_RAND = TH2D("MM_vs_beta_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_RAND = TH2D("MM_vs_H_cer_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_RAND = TH2D("MM_vs_H_cal_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_RAND = TH2D("MM_vs_P_cal_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)    
    MM_vs_P_hgcer_RAND = TH2D("MM_vs_P_hgcer_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_RAND = TH2D("MM_vs_P_aero_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_RAND = TH2D("phiq_vs_t_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_RAND = TH2D("Q2_vs_W_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_RAND = TH2D("Q2_vs_t_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_RAND = TH2D("W_vs_t_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_RAND = TH2D("EPS_vs_t_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_RAND = TH2D("MM_vs_t_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_RAND = TH2D("P_hgcer_xAtCer_vs_yAtCer_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_RAND = TH2D("P_hgcer_xAtCer_vs_MM_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_MM_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_RAND = TH2D("P_hgcer_yAtCer_vs_MM_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_RAND = TH2D("P_hgcer_nohole_yAtCer_vs_MM_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)
    
    MM_vs_CoinTime_DUMMY_RAND = TH2D("MM_vs_CoinTime_DUMMY_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
    CoinTime_vs_beta_DUMMY_RAND = TH2D("CoinTime_vs_beta_DUMMY_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
    MM_vs_beta_DUMMY_RAND = TH2D("MM_vs_beta_DUMMY_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
    MM_vs_H_cer_DUMMY_RAND = TH2D("MM_vs_H_cer_DUMMY_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
    MM_vs_H_cal_DUMMY_RAND = TH2D("MM_vs_H_cal_DUMMY_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
    MM_vs_P_cal_DUMMY_RAND = TH2D("MM_vs_P_cal_DUMMY_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)    
    MM_vs_P_hgcer_DUMMY_RAND = TH2D("MM_vs_P_hgcer_DUMMY_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
    MM_vs_P_aero_DUMMY_RAND = TH2D("MM_vs_P_aero_DUMMY_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)    
    phiq_vs_t_DUMMY_RAND = TH2D("phiq_vs_t_DUMMY_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
    Q2_vs_W_DUMMY_RAND = TH2D("Q2_vs_W_DUMMY_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
    Q2_vs_t_DUMMY_RAND = TH2D("Q2_vs_t_DUMMY_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
    W_vs_t_DUMMY_RAND = TH2D("W_vs_t_DUMMY_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
    EPS_vs_t_DUMMY_RAND = TH2D("EPS_vs_t_DUMMY_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
    MM_vs_t_DUMMY_RAND = TH2D("MM_vs_t_DUMMY_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
    # HGCer hole comparison plots
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND = TH2D("P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_xAtCer_vs_MM_DUMMY_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_yAtCer_vs_MM_DUMMY_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND = TH2D("P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)        

    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":        
        
        from particle_subtraction import particle_subtraction_cuts
        SubtractedParticle = "pion"
        subDict = {}
        
        subDict["H_hsdelta_SUB_DATA"]  = TH1D("H_hsdelta_SUB_DATA","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DATA"]  = TH1D("H_hsxptar_SUB_DATA","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DATA"]  = TH1D("H_hsyptar_SUB_DATA","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DATA"]    = TH1D("H_ssxfp_SUB_DATA","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DATA"]    = TH1D("H_ssyfp_SUB_DATA","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DATA"]   = TH1D("H_ssxpfp_SUB_DATA","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DATA"]   = TH1D("H_ssypfp_SUB_DATA","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DATA"]    = TH1D("H_hsxfp_SUB_DATA","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DATA"]    = TH1D("H_hsyfp_SUB_DATA","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DATA"]   = TH1D("H_hsxpfp_SUB_DATA","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DATA"]   = TH1D("H_hsypfp_SUB_DATA","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DATA"]  = TH1D("H_ssdelta_SUB_DATA","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DATA"]  = TH1D("H_ssxptar_SUB_DATA","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DATA"]  = TH1D("H_ssyptar_SUB_DATA","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DATA"]        = TH1D("H_q_SUB_DATA","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DATA"]       = TH1D("H_Q2_SUB_DATA","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DATA"]  = TH1D("H_W_SUB_DATA","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DATA"]       = TH1D("H_t_SUB_DATA","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DATA"]  = TH1D("H_epsilon_SUB_DATA","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DATA"]  = TH1D("H_MM_SUB_DATA","MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_nosub_SUB_DATA"]  = TH1D("H_MM_nosub_SUB_DATA","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
        subDict["H_th_SUB_DATA"]  = TH1D("H_th_SUB_DATA","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DATA"]  = TH1D("H_ph_SUB_DATA","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DATA"]  = TH1D("H_ph_q_SUB_DATA","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DATA"]  = TH1D("H_th_q_SUB_DATA","Theta Detected (th_xq)", 100, -0.2, 0.2)
        subDict["H_ph_recoil_SUB_DATA"]  = TH1D("H_ph_recoil_SUB_DATA","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        subDict["H_th_recoil_SUB_DATA"]  = TH1D("H_th_recoil_SUB_DATA","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        subDict["H_pmiss_SUB_DATA"]  = TH1D("H_pmiss_SUB_DATA","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DATA"]  = TH1D("H_emiss_SUB_DATA","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DATA"]  = TH1D("H_pmx_SUB_DATA","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DATA"]  = TH1D("H_pmy_SUB_DATA","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DATA"]  = TH1D("H_pmz_SUB_DATA","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DATA"] = TH1D("H_ct_SUB_DATA", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DATA"] = TH1D("H_cal_etottracknorm_SUB_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DATA"] = TH1D("H_cer_npeSum_SUB_DATA", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DATA"] = TH1D("P_cal_etottracknorm_SUB_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DATA"] = TH1D("P_hgcer_npeSum_SUB_DATA", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DATA"] = TH1D("P_aero_npeSum_SUB_DATA", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["H_hsdelta_SUB_RAND"]  = TH1D("H_hsdelta_SUB_RAND","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_RAND"]  = TH1D("H_hsxptar_SUB_RAND","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_RAND"]  = TH1D("H_hsyptar_SUB_RAND","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_RAND"]    = TH1D("H_ssxfp_SUB_RAND","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_RAND"]    = TH1D("H_ssyfp_SUB_RAND","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_RAND"]   = TH1D("H_ssxpfp_SUB_RAND","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_RAND"]   = TH1D("H_ssypfp_SUB_RAND","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_RAND"]    = TH1D("H_hsxfp_SUB_RAND","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_RAND"]    = TH1D("H_hsyfp_SUB_RAND","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_RAND"]   = TH1D("H_hsxpfp_SUB_RAND","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_RAND"]   = TH1D("H_hsypfp_SUB_RAND","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_RAND"]  = TH1D("H_ssdelta_SUB_RAND","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_RAND"]  = TH1D("H_ssxptar_SUB_RAND","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_RAND"]  = TH1D("H_ssyptar_SUB_RAND","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_RAND"]        = TH1D("H_q_SUB_RAND","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_RAND"]       = TH1D("H_Q2_SUB_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_RAND"]  = TH1D("H_W_SUB_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_RAND"]       = TH1D("H_t_SUB_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_RAND"]  = TH1D("H_epsilon_SUB_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_RAND"]  = TH1D("H_MM_SUB_RAND","MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_nosub_SUB_RAND"]  = TH1D("H_MM_nosub_SUB_RAND","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
        subDict["H_th_SUB_RAND"]  = TH1D("H_th_SUB_RAND","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_RAND"]  = TH1D("H_ph_SUB_RAND","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_RAND"]  = TH1D("H_ph_q_SUB_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_RAND"]  = TH1D("H_th_q_SUB_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        subDict["H_ph_recoil_SUB_RAND"]  = TH1D("H_ph_recoil_SUB_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        subDict["H_th_recoil_SUB_RAND"]  = TH1D("H_th_recoil_SUB_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        subDict["H_pmiss_SUB_RAND"]  = TH1D("H_pmiss_SUB_RAND","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_RAND"]  = TH1D("H_emiss_SUB_RAND","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_RAND"]  = TH1D("H_pmx_SUB_RAND","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_RAND"]  = TH1D("H_pmy_SUB_RAND","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_RAND"]  = TH1D("H_pmz_SUB_RAND","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_RAND"] = TH1D("H_ct_SUB_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_RAND"] = TH1D("H_cal_etottracknorm_SUB_RAND", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_RAND"] = TH1D("H_cer_npeSum_SUB_RAND", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_RAND"] = TH1D("P_cal_etottracknorm_SUB_RAND", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_RAND"] = TH1D("P_hgcer_npeSum_SUB_RAND", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_RAND"] = TH1D("P_aero_npeSum_SUB_RAND", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["H_hsdelta_SUB_DUMMY"]  = TH1D("H_hsdelta_SUB_DUMMY","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DUMMY"]  = TH1D("H_hsxptar_SUB_DUMMY","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DUMMY"]  = TH1D("H_hsyptar_SUB_DUMMY","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DUMMY"]    = TH1D("H_ssxfp_SUB_DUMMY","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DUMMY"]    = TH1D("H_ssyfp_SUB_DUMMY","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DUMMY"]   = TH1D("H_ssxpfp_SUB_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DUMMY"]   = TH1D("H_ssypfp_SUB_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DUMMY"]    = TH1D("H_hsxfp_SUB_DUMMY","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DUMMY"]    = TH1D("H_hsyfp_SUB_DUMMY","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DUMMY"]   = TH1D("H_hsxpfp_SUB_DUMMY","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DUMMY"]   = TH1D("H_hsypfp_SUB_DUMMY","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DUMMY"]  = TH1D("H_ssdelta_SUB_DUMMY","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DUMMY"]  = TH1D("H_ssxptar_SUB_DUMMY","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DUMMY"]  = TH1D("H_ssyptar_SUB_DUMMY","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DUMMY"]        = TH1D("H_q_SUB_DUMMY","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DUMMY"]       = TH1D("H_Q2_SUB_DUMMY","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DUMMY"]  = TH1D("H_W_SUB_DUMMY","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DUMMY"]       = TH1D("H_t_SUB_DUMMY","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DUMMY"]  = TH1D("H_epsilon_SUB_DUMMY","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DUMMY"]  = TH1D("H_MM_SUB_DUMMY","MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_nosub_SUB_DUMMY"]  = TH1D("H_MM_nosub_SUB_DUMMY","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
        subDict["H_th_SUB_DUMMY"]  = TH1D("H_th_SUB_DUMMY","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DUMMY"]  = TH1D("H_ph_SUB_DUMMY","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DUMMY"]  = TH1D("H_ph_q_SUB_DUMMY","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DUMMY"]  = TH1D("H_th_q_SUB_DUMMY","Theta Detected (th_xq)", 100, -0.2, 0.2)
        subDict["H_ph_recoil_SUB_DUMMY"]  = TH1D("H_ph_recoil_SUB_DUMMY","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        subDict["H_th_recoil_SUB_DUMMY"]  = TH1D("H_th_recoil_SUB_DUMMY","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        subDict["H_pmiss_SUB_DUMMY"]  = TH1D("H_pmiss_SUB_DUMMY","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DUMMY"]  = TH1D("H_emiss_SUB_DUMMY","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DUMMY"]  = TH1D("H_pmx_SUB_DUMMY","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DUMMY"]  = TH1D("H_pmy_SUB_DUMMY","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DUMMY"]  = TH1D("H_pmz_SUB_DUMMY","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DUMMY"] = TH1D("H_ct_SUB_DUMMY", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DUMMY"] = TH1D("H_cal_etottracknorm_SUB_DUMMY", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DUMMY"] = TH1D("H_cer_npeSum_SUB_DUMMY", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DUMMY"] = TH1D("P_cal_etottracknorm_SUB_DUMMY", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DUMMY"] = TH1D("P_hgcer_npeSum_SUB_DUMMY", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DUMMY"] = TH1D("P_aero_npeSum_SUB_DUMMY", "SHMS Aero Npe Sum", 100, 0, 30)        

        subDict["H_hsdelta_SUB_DUMMY_RAND"]  = TH1D("H_hsdelta_SUB_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
        subDict["H_hsxptar_SUB_DUMMY_RAND"]  = TH1D("H_hsxptar_SUB_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
        subDict["H_hsyptar_SUB_DUMMY_RAND"]  = TH1D("H_hsyptar_SUB_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
        subDict["H_ssxfp_SUB_DUMMY_RAND"]    = TH1D("H_ssxfp_SUB_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
        subDict["H_ssyfp_SUB_DUMMY_RAND"]    = TH1D("H_ssyfp_SUB_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
        subDict["H_ssxpfp_SUB_DUMMY_RAND"]   = TH1D("H_ssxpfp_SUB_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
        subDict["H_ssypfp_SUB_DUMMY_RAND"]   = TH1D("H_ssypfp_SUB_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
        subDict["H_hsxfp_SUB_DUMMY_RAND"]    = TH1D("H_hsxfp_SUB_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
        subDict["H_hsyfp_SUB_DUMMY_RAND"]    = TH1D("H_hsyfp_SUB_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
        subDict["H_hsxpfp_SUB_DUMMY_RAND"]   = TH1D("H_hsxpfp_SUB_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
        subDict["H_hsypfp_SUB_DUMMY_RAND"]   = TH1D("H_hsypfp_SUB_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
        subDict["H_ssdelta_SUB_DUMMY_RAND"]  = TH1D("H_ssdelta_SUB_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
        subDict["H_ssxptar_SUB_DUMMY_RAND"]  = TH1D("H_ssxptar_SUB_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
        subDict["H_ssyptar_SUB_DUMMY_RAND"]  = TH1D("H_ssyptar_SUB_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
        subDict["H_q_SUB_DUMMY_RAND"]        = TH1D("H_q_SUB_DUMMY_RAND","q", 100, 0.0, 10.0)
        subDict["H_Q2_SUB_DUMMY_RAND"]       = TH1D("H_Q2_SUB_DUMMY_RAND","Q2", 100, inpDict["Q2min"], inpDict["Q2max"])
        subDict["H_W_SUB_DUMMY_RAND"]  = TH1D("H_W_SUB_DUMMY_RAND","W ", 100, inpDict["Wmin"], inpDict["Wmax"])
        subDict["H_t_SUB_DUMMY_RAND"]       = TH1D("H_t_SUB_DUMMY_RAND","-t", 100, inpDict["tmin"], inpDict["tmax"])
        subDict["H_epsilon_SUB_DUMMY_RAND"]  = TH1D("H_epsilon_SUB_DUMMY_RAND","epsilon", 100, inpDict["Epsmin"], inpDict["Epsmax"])
        subDict["H_MM_SUB_DUMMY_RAND"]  = TH1D("H_MM_SUB_DUMMY_RAND","MM_{}".format(SubtractedParticle), 100, inpDict["mm_min"], inpDict["mm_max"])
        subDict["H_MM_nosub_SUB_DUMMY_RAND"]  = TH1D("H_MM_nosub_SUB_DUMMY_RAND","MM_{}".format(SubtractedParticle), 100, 0.7, 1.5)
        subDict["H_th_SUB_DUMMY_RAND"]  = TH1D("H_th_SUB_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
        subDict["H_ph_SUB_DUMMY_RAND"]  = TH1D("H_ph_SUB_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
        subDict["H_ph_q_SUB_DUMMY_RAND"]  = TH1D("H_ph_q_SUB_DUMMY_RAND","Phi Detected (ph_xq)", 100, -math.pi, math.pi)
        subDict["H_th_q_SUB_DUMMY_RAND"]  = TH1D("H_th_q_SUB_DUMMY_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        subDict["H_ph_recoil_SUB_DUMMY_RAND"]  = TH1D("H_ph_recoil_SUB_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        subDict["H_th_recoil_SUB_DUMMY_RAND"]  = TH1D("H_th_recoil_SUB_DUMMY_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        subDict["H_pmiss_SUB_DUMMY_RAND"]  = TH1D("H_pmiss_SUB_DUMMY_RAND","pmiss", 100, 0.0, 2.0)
        subDict["H_emiss_SUB_DUMMY_RAND"]  = TH1D("H_emiss_SUB_DUMMY_RAND","emiss", 100, 0.0, 2.0)
        subDict["H_pmx_SUB_DUMMY_RAND"]  = TH1D("H_pmx_SUB_DUMMY_RAND","pmx", 100, -10.0, 10.0)
        subDict["H_pmy_SUB_DUMMY_RAND"]  = TH1D("H_pmy_SUB_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
        subDict["H_pmz_SUB_DUMMY_RAND"]  = TH1D("H_pmz_SUB_DUMMY_RAND","pmz", 100, -10.0, 10.0)
        subDict["H_ct_SUB_DUMMY_RAND"] = TH1D("H_ct_SUB_DUMMY_RAND", "Electron-{} CTime".format(ParticleType.capitalize()), 100, -50, 50)
        subDict["H_cal_etottracknorm_SUB_DUMMY_RAND"] = TH1D("H_cal_etottracknorm_SUB_DUMMY_RAND", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        subDict["H_cer_npeSum_SUB_DUMMY_RAND"] = TH1D("H_cer_npeSum_SUB_DUMMY_RAND", "HMS Cer Npe Sum", 100, 0, 30)
        subDict["P_cal_etottracknorm_SUB_DUMMY_RAND"] = TH1D("P_cal_etottracknorm_SUB_DUMMY_RAND", "SHMS Cal etottracknorm", 100, 0, 1)
        subDict["P_hgcer_npeSum_SUB_DUMMY_RAND"] = TH1D("P_hgcer_npeSum_SUB_DUMMY_RAND", "SHMS HGCer Npe Sum", 100, 0, 10)
        subDict["P_aero_npeSum_SUB_DUMMY_RAND"] = TH1D("P_aero_npeSum_SUB_DUMMY_RAND", "SHMS Aero Npe Sum", 100, 0, 30)

        subDict["MM_vs_CoinTime_SUB_DATA"] = TH2D("MM_vs_CoinTime_SUB_DATA","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DATA"] = TH2D("CoinTime_vs_beta_SUB_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DATA"] = TH2D("MM_vs_beta_SUB_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DATA"] = TH2D("MM_vs_H_cer_SUB_DATA", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DATA"] = TH2D("MM_vs_H_cal_SUB_DATA", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DATA"] = TH2D("MM_vs_P_cal_SUB_DATA", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DATA"] = TH2D("MM_vs_P_hgcer_SUB_DATA", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DATA"] = TH2D("MM_vs_P_aero_SUB_DATA", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DATA"] = TH2D("phiq_vs_t_SUB_DATA","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DATA"] = TH2D("Q2_vs_W_SUB_DATA", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DATA"] = TH2D("Q2_vs_t_SUB_DATA", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DATA"] = TH2D("W_vs_t_SUB_DATA", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DATA"] = TH2D("EPS_vs_t_SUB_DATA", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DATA"] = TH2D("MM_vs_t_SUB_DATA", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DATA", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DATA", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DATA", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_DUMMY"] = TH2D("MM_vs_CoinTime_SUB_DUMMY","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DUMMY"] = TH2D("CoinTime_vs_beta_SUB_DUMMY", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DUMMY"] = TH2D("MM_vs_beta_SUB_DUMMY", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DUMMY"] = TH2D("MM_vs_H_cer_SUB_DUMMY", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DUMMY"] = TH2D("MM_vs_H_cal_SUB_DUMMY", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DUMMY"] = TH2D("MM_vs_P_cal_SUB_DUMMY", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DUMMY"] = TH2D("MM_vs_P_hgcer_SUB_DUMMY", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DUMMY"] = TH2D("MM_vs_P_aero_SUB_DUMMY", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DUMMY"] = TH2D("phiq_vs_t_SUB_DUMMY","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DUMMY"] = TH2D("Q2_vs_W_SUB_DUMMY", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DUMMY"] = TH2D("Q2_vs_t_SUB_DUMMY", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DUMMY"] = TH2D("W_vs_t_SUB_DUMMY", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DUMMY"] = TH2D("EPS_vs_t_SUB_DUMMY", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DUMMY"] = TH2D("MM_vs_t_SUB_DUMMY", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DUMMY", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DUMMY", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_RAND"] = TH2D("MM_vs_CoinTime_SUB_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_RAND"] = TH2D("CoinTime_vs_beta_SUB_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_RAND"] = TH2D("MM_vs_beta_SUB_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_RAND"] = TH2D("MM_vs_H_cer_SUB_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_RAND"] = TH2D("MM_vs_H_cal_SUB_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_RAND"] = TH2D("MM_vs_P_cal_SUB_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_RAND"] = TH2D("MM_vs_P_hgcer_SUB_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_RAND"] = TH2D("MM_vs_P_aero_SUB_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_RAND"] = TH2D("phiq_vs_t_SUB_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_RAND"] = TH2D("Q2_vs_W_SUB_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_RAND"] = TH2D("Q2_vs_t_SUB_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_RAND"] = TH2D("W_vs_t_SUB_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_RAND"] = TH2D("EPS_vs_t_SUB_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_RAND"] = TH2D("MM_vs_t_SUB_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_RAND"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

        subDict["MM_vs_CoinTime_SUB_DUMMY_RAND"] = TH2D("MM_vs_CoinTime_SUB_DUMMY_RAND","Missing Mass vs CTime; MM; Coin_Time",100, inpDict["mm_min"], inpDict["mm_max"], 100, -50, 50)
        subDict["CoinTime_vs_beta_SUB_DUMMY_RAND"] = TH2D("CoinTime_vs_beta_SUB_DUMMY_RAND", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -10, 10, 100, 0, 2)
        subDict["MM_vs_beta_SUB_DUMMY_RAND"] = TH2D("MM_vs_beta_SUB_DUMMY_RAND", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 2)
        subDict["MM_vs_H_cer_SUB_DUMMY_RAND"] = TH2D("MM_vs_H_cer_SUB_DUMMY_RAND", "Missing Mass vs HMS Cerenkov; MM; HMS Cerenkov", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["MM_vs_H_cal_SUB_DUMMY_RAND"] = TH2D("MM_vs_H_cal_SUB_DUMMY_RAND", "Missing Mass vs HMS Cal eTrackNorm; MM; HMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0.2, 1.8)
        subDict["MM_vs_P_cal_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_cal_SUB_DUMMY_RAND", "Missing Mass vs SHMS Cal eTrackNorm; MM; SHMS Cal eTrackNorm", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 1)
        subDict["MM_vs_P_hgcer_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_hgcer_SUB_DUMMY_RAND", "Missing Mass vs SHMS HGCer; MM; SHMS HGCer", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 10)
        subDict["MM_vs_P_aero_SUB_DUMMY_RAND"] = TH2D("MM_vs_P_aero_SUB_DUMMY_RAND", "Missing Mass vs SHMS Aerogel; MM; SHMS Aerogel", 100, inpDict["mm_min"], inpDict["mm_max"], 100, 0, 30)
        subDict["phiq_vs_t_SUB_DUMMY_RAND"] = TH2D("phiq_vs_t_SUB_DUMMY_RAND","; #phi ;t", 12, -3.14, 3.14, 24, inpDict["tmin"], inpDict["tmax"])
        subDict["Q2_vs_W_SUB_DUMMY_RAND"] = TH2D("Q2_vs_W_SUB_DUMMY_RAND", "Q^{2} vs W; Q^{2}; W", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["Wmin"], inpDict["Wmax"])
        subDict["Q2_vs_t_SUB_DUMMY_RAND"] = TH2D("Q2_vs_t_SUB_DUMMY_RAND", "Q^{2} vs t; Q^{2}; t", 50, inpDict["Q2min"], inpDict["Q2max"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["W_vs_t_SUB_DUMMY_RAND"] = TH2D("W_vs_t_SUB_DUMMY_RAND", "W vs t; W; t", 50, inpDict["Wmin"], inpDict["Wmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["EPS_vs_t_SUB_DUMMY_RAND"] = TH2D("EPS_vs_t_SUB_DUMMY_RAND", "Epsilon vs t; Epsilon; t", 50, inpDict["Epsmin"], inpDict["Epsmax"], 50, inpDict["tmin"], inpDict["tmax"])
        subDict["MM_vs_t_SUB_DUMMY_RAND"] = TH2D("MM_vs_t_SUB_DUMMY_RAND", "Missing Mass vs t; MM; t", 100, inpDict["mm_min"], inpDict["mm_max"], 100, inpDict["tmin"], inpDict["tmax"])
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"] = TH2D("P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND", "X vs Y; X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND", "X vs Y (no hole cut); X; Y", 50, -30, 30, 50, -30, 30)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND", "X vs MM; X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND", "X vs MM (no hole cut); X; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND", "Y vs MM; Y; MM", 50, -30, 30, 50, 0, 2)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND"] = TH2D("P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND", "Y vs MM (no hole cut); Y; MM", 50, -30, 30, 50, 0, 2)

    # Fit background and subtract
    from background_fit import bg_fit
        
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = evt.MM_shift
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOHOLECUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            if(NOHOLECUTS):
                # HGCer hole comparison            
                P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
                P_hgcer_nohole_xAtCer_vs_MM_DATA.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
                P_hgcer_nohole_yAtCer_vs_MM_DATA.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)
            
        if(NOMMCUTS):
            H_MM_nosub_DATA.Fill(adj_MM)            
            
        if(ALLCUTS):

          # HGCer hole comparison
          P_hgcer_xAtCer_vs_yAtCer_DATA.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
          P_hgcer_xAtCer_vs_MM_DATA.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
          P_hgcer_yAtCer_vs_MM_DATA.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

          # Phase shift to fix polar plots
          #phi_shift = (evt.ph_q+math.pi)
          phi_shift = (evt.ph_q)          
          
          MM_vs_CoinTime_DATA.Fill(adj_MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DATA.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DATA.Fill(adj_MM,evt.P_gtr_beta)
          MM_vs_H_cer_DATA.Fill(adj_MM,evt.H_cer_npeSum)
          MM_vs_H_cal_DATA.Fill(adj_MM,evt.H_cal_etottracknorm)
          MM_vs_P_cal_DATA.Fill(adj_MM,evt.P_cal_etottracknorm)
          MM_vs_P_hgcer_DATA.Fill(adj_MM,evt.P_hgcer_npeSum)
          MM_vs_P_aero_DATA.Fill(adj_MM,evt.P_aero_npeSum)
          # SIMC goes from 0 to 2pi so no need for +pi
          phiq_vs_t_DATA.Fill(phi_shift, -evt.MandelT)
          Q2_vs_W_DATA.Fill(evt.Q2, evt.W)
          Q2_vs_t_DATA.Fill(evt.Q2, -evt.MandelT)
          W_vs_t_DATA.Fill(evt.W, -evt.MandelT)
          EPS_vs_t_DATA.Fill(evt.epsilon, -evt.MandelT)
          MM_vs_t_DATA.Fill(adj_MM, -evt.MandelT)
          polar_phiq_vs_t_DATA.SetPoint(polar_phiq_vs_t_DATA.GetN(), (phi_shift)*(180/math.pi), -evt.MandelT)
          
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
          H_hsdelta_DATA.Fill(adj_hsdelta)
          H_hsxptar_DATA.Fill(evt.hsxptar)	
          H_hsyptar_DATA.Fill(evt.hsyptar)

          # SIMC goes from 0 to 2pi so no need for +pi          
          H_ph_q_DATA.Fill((phi_shift))
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
          H_MM_DATA.Fill(adj_MM)
          #H_MM_DATA.Fill(pow(adj_MM, 2))
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

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = evt.MM_shift
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOHOLECUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            if(NOHOLECUTS):
                # HGCer hole comparison            
                P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
                P_hgcer_nohole_xAtCer_vs_MM_DUMMY.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
                P_hgcer_nohole_yAtCer_vs_MM_DUMMY.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)
            
        if(NOMMCUTS):
            H_MM_nosub_DUMMY.Fill(adj_MM)            
            
        if(ALLCUTS):

          # HGCer hole comparison
          P_hgcer_xAtCer_vs_yAtCer_DUMMY.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
          P_hgcer_xAtCer_vs_MM_DUMMY.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
          P_hgcer_yAtCer_vs_MM_DUMMY.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)


          # Phase shift to fix polar plots
          #phi_shift = (evt.ph_q+math.pi)
          phi_shift = (evt.ph_q)          
          
          MM_vs_CoinTime_DUMMY.Fill(adj_MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DUMMY.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DUMMY.Fill(adj_MM,evt.P_gtr_beta)
          MM_vs_H_cer_DUMMY.Fill(adj_MM,evt.H_cer_npeSum)
          MM_vs_H_cal_DUMMY.Fill(adj_MM,evt.H_cal_etottracknorm)
          MM_vs_P_cal_DUMMY.Fill(adj_MM,evt.P_cal_etottracknorm)
          MM_vs_P_hgcer_DUMMY.Fill(adj_MM,evt.P_hgcer_npeSum)
          MM_vs_P_aero_DUMMY.Fill(adj_MM,evt.P_aero_npeSum)          
          # SIMC goes from 0 to 2pi so no need for +pi
          phiq_vs_t_DUMMY.Fill(phi_shift, -evt.MandelT)
          Q2_vs_W_DUMMY.Fill(evt.Q2, evt.W)
          Q2_vs_t_DUMMY.Fill(evt.Q2, -evt.MandelT)
          W_vs_t_DUMMY.Fill(evt.W, -evt.MandelT)
          EPS_vs_t_DUMMY.Fill(evt.epsilon, -evt.MandelT)
          MM_vs_t_DUMMY.Fill(adj_MM, -evt.MandelT)
          polar_phiq_vs_t_DUMMY.SetPoint(polar_phiq_vs_t_DUMMY.GetN(), (phi_shift)*(180/math.pi), -evt.MandelT)

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
          H_hsdelta_DUMMY.Fill(adj_hsdelta)
          H_hsxptar_DUMMY.Fill(evt.hsxptar)	
          H_hsyptar_DUMMY.Fill(evt.hsyptar)

          # SIMC goes from 0 to 2pi so no need for +pi          
          H_ph_q_DUMMY.Fill((phi_shift))
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
          H_MM_DUMMY.Fill(adj_MM)
          #H_MM_DUMMY.Fill(pow(adj_MM, 2))  
          #H_MM_DUMMY.Fill(evt.Mrecoil)

    ###################################################################################################################################################    
    # Fill random histograms for various trees called above

    print("\nGrabbing {} {} random data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)
        
        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = evt.MM_shift
        
        ##############
        ##############        
        ##############

        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOHOLECUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            if(NOHOLECUTS):
                # HGCer hole comparison            
                P_hgcer_nohole_xAtCer_vs_yAtCer_RAND.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
                P_hgcer_nohole_xAtCer_vs_MM_RAND.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
                P_hgcer_nohole_yAtCer_vs_MM_RAND.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)
            
        if(NOMMCUTS):
            H_MM_nosub_RAND.Fill(adj_MM)
            
        if(ALLCUTS):

          # HGCer hole comparison
          P_hgcer_xAtCer_vs_yAtCer_RAND.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
          P_hgcer_xAtCer_vs_MM_RAND.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
          P_hgcer_yAtCer_vs_MM_RAND.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

          # Phase shift to fix polar plots
          #phi_shift = (evt.ph_q+math.pi)
          phi_shift = (evt.ph_q)          
          
          MM_vs_CoinTime_RAND.Fill(adj_MM, evt.CTime_ROC1)
          CoinTime_vs_beta_RAND.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_RAND.Fill(adj_MM,evt.P_gtr_beta)
          MM_vs_H_cer_RAND.Fill(adj_MM,evt.H_cer_npeSum)
          MM_vs_H_cal_RAND.Fill(adj_MM,evt.H_cal_etottracknorm)
          MM_vs_P_cal_RAND.Fill(adj_MM,evt.P_cal_etottracknorm)
          MM_vs_P_hgcer_RAND.Fill(adj_MM,evt.P_hgcer_npeSum)
          MM_vs_P_aero_RAND.Fill(adj_MM,evt.P_aero_npeSum)          
          # SIMC goes from 0 to 2pi so no need for +pi
          phiq_vs_t_RAND.Fill(phi_shift, -evt.MandelT)
          Q2_vs_W_RAND.Fill(evt.Q2, evt.W)
          Q2_vs_t_RAND.Fill(evt.Q2, -evt.MandelT)
          W_vs_t_RAND.Fill(evt.W, -evt.MandelT)
          EPS_vs_t_RAND.Fill(evt.epsilon, -evt.MandelT)
          MM_vs_t_RAND.Fill(adj_MM, -evt.MandelT)

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
          H_hsdelta_RAND.Fill(adj_hsdelta)
          H_hsxptar_RAND.Fill(evt.hsxptar)	
          H_hsyptar_RAND.Fill(evt.hsyptar)

          # SIMC goes from 0 to 2pi so no need for +pi          
          H_ph_q_RAND.Fill((phi_shift))
          H_th_q_RAND.Fill(evt.th_q)
          H_ph_recoil_RAND.Fill(evt.ph_recoil)
          H_th_recoil_RAND.Fill(evt.th_recoil)
          
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
          H_MM_RAND.Fill(adj_MM)

    ###################################################################################################################################################    
    # Fill dummy random histograms for various trees called above

    print("\nGrabbing {} {} dummy random data...".format(phi_setting,ParticleType))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        adj_MM = evt.MM_shift
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            NOHOLECUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) #and evt.P_hgcer_npeSum == 0.0
            if(NOHOLECUTS):
                # HGCer hole comparison            
                P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
                P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
                P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)
            
        if(NOMMCUTS):
            H_MM_nosub_DUMMY_RAND.Fill(adj_MM)
            
        if(ALLCUTS):

          # HGCer hole comparison
          P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
          P_hgcer_xAtCer_vs_MM_DUMMY_RAND.Fill(evt.P_hgcer_xAtCer,adj_MM, evt.P_hgcer_npeSum)
          P_hgcer_yAtCer_vs_MM_DUMMY_RAND.Fill(evt.P_hgcer_yAtCer,adj_MM, evt.P_hgcer_npeSum)

          # Phase shift to fix polar plots
          #phi_shift = (evt.ph_q+math.pi)
          phi_shift = (evt.ph_q)          
          
          MM_vs_CoinTime_DUMMY_RAND.Fill(adj_MM, evt.CTime_ROC1)
          CoinTime_vs_beta_DUMMY_RAND.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
          MM_vs_beta_DUMMY_RAND.Fill(adj_MM,evt.P_gtr_beta)
          MM_vs_H_cer_DUMMY_RAND.Fill(adj_MM,evt.H_cer_npeSum)
          MM_vs_H_cal_DUMMY_RAND.Fill(adj_MM,evt.H_cal_etottracknorm)
          MM_vs_P_cal_DUMMY_RAND.Fill(adj_MM,evt.P_cal_etottracknorm)
          MM_vs_P_hgcer_DUMMY_RAND.Fill(adj_MM,evt.P_hgcer_npeSum)
          MM_vs_P_aero_DUMMY_RAND.Fill(adj_MM,evt.P_aero_npeSum)          
          # SIMC goes from 0 to 2pi so no need for +pi
          phiq_vs_t_DUMMY_RAND.Fill(phi_shift, -evt.MandelT)
          Q2_vs_W_DUMMY_RAND.Fill(evt.Q2, evt.W)
          Q2_vs_t_DUMMY_RAND.Fill(evt.Q2, -evt.MandelT)
          W_vs_t_DUMMY_RAND.Fill(evt.W, -evt.MandelT)
          EPS_vs_t_DUMMY_RAND.Fill(evt.epsilon, -evt.MandelT)
          MM_vs_t_DUMMY_RAND.Fill(adj_MM, -evt.MandelT)
          
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
          H_hsdelta_DUMMY_RAND.Fill(adj_hsdelta)
          H_hsxptar_DUMMY_RAND.Fill(evt.hsxptar)	
          H_hsyptar_DUMMY_RAND.Fill(evt.hsyptar)

          # SIMC goes from 0 to 2pi so no need for +pi          
          H_ph_q_DUMMY_RAND.Fill((phi_shift))
          H_th_q_DUMMY_RAND.Fill(evt.th_q)
          H_ph_recoil_DUMMY_RAND.Fill(evt.ph_recoil)
          H_th_recoil_DUMMY_RAND.Fill(evt.th_recoil)
          
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
          H_MM_DUMMY_RAND.Fill(adj_MM)

    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge

    # Data Random subtraction window
    P_hgcer_xAtCer_vs_yAtCer_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_RAND.Scale(1/nWindows)
    P_hgcer_xAtCer_vs_MM_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_RAND.Scale(1/nWindows)
    P_hgcer_yAtCer_vs_MM_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_RAND.Scale(1/nWindows)        
    MM_vs_CoinTime_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_RAND.Scale(1/nWindows)
    MM_vs_beta_RAND.Scale(1/nWindows)
    MM_vs_H_cer_RAND.Scale(1/nWindows)
    MM_vs_H_cal_RAND.Scale(1/nWindows)
    MM_vs_P_cal_RAND.Scale(1/nWindows)
    MM_vs_P_hgcer_RAND.Scale(1/nWindows)
    MM_vs_P_aero_RAND.Scale(1/nWindows)
    phiq_vs_t_RAND.Scale(1/nWindows)
    Q2_vs_W_RAND.Scale(1/nWindows)
    Q2_vs_t_RAND.Scale(1/nWindows)
    W_vs_t_RAND.Scale(1/nWindows)
    EPS_vs_t_RAND.Scale(1/nWindows)
    MM_vs_t_RAND.Scale(1/nWindows)    
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
    H_ph_q_RAND.Scale(1/nWindows)
    H_th_q_RAND.Scale(1/nWindows)
    H_ph_recoil_RAND.Scale(1/nWindows)
    H_th_recoil_RAND.Scale(1/nWindows)
    H_Q2_RAND.Scale(1/nWindows)
    H_W_RAND.Scale(1/nWindows)    
    H_t_RAND.Scale(1/nWindows)
    H_epsilon_RAND.Scale(1/nWindows)
    H_MM_RAND.Scale(1/nWindows)
    H_MM_nosub_RAND.Scale(1/nWindows)
    H_pmiss_RAND.Scale(1/nWindows)
    H_emiss_RAND.Scale(1/nWindows)
    H_pmx_RAND.Scale(1/nWindows)
    H_pmy_RAND.Scale(1/nWindows)
    H_pmz_RAND.Scale(1/nWindows)

    # Data Dummy_Random subtraction window
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND.Scale(1/nWindows)
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND.Scale(1/nWindows)            
    MM_vs_CoinTime_DUMMY_RAND.Scale(1/nWindows)
    CoinTime_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_beta_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_H_cer_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_H_cal_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_P_cal_DUMMY_RAND.Scale(1/nWindows)    
    MM_vs_P_hgcer_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_P_aero_DUMMY_RAND.Scale(1/nWindows)    
    phiq_vs_t_DUMMY_RAND.Scale(1/nWindows)
    Q2_vs_W_DUMMY_RAND.Scale(1/nWindows)
    Q2_vs_t_DUMMY_RAND.Scale(1/nWindows)
    W_vs_t_DUMMY_RAND.Scale(1/nWindows)
    EPS_vs_t_DUMMY_RAND.Scale(1/nWindows)
    MM_vs_t_DUMMY_RAND.Scale(1/nWindows)
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
    H_ph_q_DUMMY_RAND.Scale(1/nWindows)
    H_th_q_DUMMY_RAND.Scale(1/nWindows)
    H_ph_recoil_DUMMY_RAND.Scale(1/nWindows)
    H_th_recoil_DUMMY_RAND.Scale(1/nWindows)    
    H_Q2_DUMMY_RAND.Scale(1/nWindows)
    H_W_DUMMY_RAND.Scale(1/nWindows)
    H_t_DUMMY_RAND.Scale(1/nWindows)
    H_epsilon_DUMMY_RAND.Scale(1/nWindows)
    H_MM_DUMMY_RAND.Scale(1/nWindows)
    H_MM_nosub_DUMMY_RAND.Scale(1/nWindows)
    H_pmiss_DUMMY_RAND.Scale(1/nWindows)
    H_emiss_DUMMY_RAND.Scale(1/nWindows)
    H_pmx_DUMMY_RAND.Scale(1/nWindows)
    H_pmy_DUMMY_RAND.Scale(1/nWindows)
    H_pmz_DUMMY_RAND.Scale(1/nWindows)
    #H_ct_DUMMY_RAND.Scale(1/nWindows)

    ###
    # Data Random subtraction
    P_hgcer_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_xAtCer_vs_yAtCer_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(P_hgcer_nohole_xAtCer_vs_yAtCer_RAND,-1)
    P_hgcer_xAtCer_vs_MM_DATA.Add(P_hgcer_xAtCer_vs_MM_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(P_hgcer_nohole_xAtCer_vs_MM_RAND,-1)
    P_hgcer_yAtCer_vs_MM_DATA.Add(P_hgcer_yAtCer_vs_MM_RAND,-1)
    if ParticleType == "kaon":    
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(P_hgcer_nohole_yAtCer_vs_MM_RAND,-1)        
    MM_vs_CoinTime_DATA.Add(MM_vs_CoinTime_RAND,-1)
    CoinTime_vs_beta_DATA.Add(CoinTime_vs_beta_RAND,-1)
    MM_vs_beta_DATA.Add(MM_vs_beta_RAND,-1)
    MM_vs_H_cer_DATA.Add(MM_vs_H_cer_RAND,-1)
    MM_vs_H_cal_DATA.Add(MM_vs_H_cal_RAND,-1)
    MM_vs_P_cal_DATA.Add(MM_vs_P_cal_RAND,-1)    
    MM_vs_P_hgcer_DATA.Add(MM_vs_P_hgcer_RAND,-1)
    MM_vs_P_aero_DATA.Add(MM_vs_P_aero_RAND,-1)
    phiq_vs_t_DATA.Add(phiq_vs_t_RAND,-1)
    Q2_vs_W_DATA.Add(Q2_vs_W_RAND,-1)
    Q2_vs_t_DATA.Add(Q2_vs_t_RAND,-1)
    W_vs_t_DATA.Add(W_vs_t_RAND,-1)
    EPS_vs_t_DATA.Add(EPS_vs_t_RAND,-1)
    MM_vs_t_DATA.Add(MM_vs_t_RAND,-1)    
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
    H_ph_q_DATA.Add(H_ph_q_RAND,-1)
    H_th_q_DATA.Add(H_th_q_RAND,-1)
    H_ph_recoil_DATA.Add(H_ph_recoil_RAND,-1)
    H_th_recoil_DATA.Add(H_th_recoil_RAND,-1)
    H_Q2_DATA.Add(H_Q2_RAND,-1)
    H_W_DATA.Add(H_W_RAND,-1)
    H_t_DATA.Add(H_t_RAND,-1)
    H_epsilon_DATA.Add(H_epsilon_RAND,-1)
    H_MM_DATA.Add(H_MM_RAND,-1)
    H_MM_nosub_DATA.Add(H_MM_nosub_RAND,-1)
    H_pmiss_DATA.Add(H_pmiss_RAND,-1)
    H_emiss_DATA.Add(H_emiss_RAND,-1)
    H_pmx_DATA.Add(H_pmx_RAND,-1)
    H_pmy_DATA.Add(H_pmy_RAND,-1)
    H_pmz_DATA.Add(H_pmz_RAND,-1)
    H_ct_DATA.Add(H_ct_RAND,-1)

    ###
    # Dummy Random subtraction
    P_hgcer_xAtCer_vs_yAtCer_DUMMY.Add(P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY.Add(P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND,-1)
    P_hgcer_xAtCer_vs_MM_DUMMY.Add(P_hgcer_xAtCer_vs_MM_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_xAtCer_vs_MM_DUMMY.Add(P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND,-1)
    P_hgcer_yAtCer_vs_MM_DUMMY.Add(P_hgcer_yAtCer_vs_MM_DUMMY_RAND,-1)
    if ParticleType == "kaon":
        P_hgcer_nohole_yAtCer_vs_MM_DUMMY.Add(P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND,-1)                
    MM_vs_CoinTime_DUMMY.Add(MM_vs_CoinTime_DUMMY_RAND,-1)
    CoinTime_vs_beta_DUMMY.Add(CoinTime_vs_beta_DUMMY_RAND,-1)
    MM_vs_beta_DUMMY.Add(MM_vs_beta_DUMMY_RAND,-1)
    MM_vs_H_cer_DUMMY.Add(MM_vs_H_cer_DUMMY_RAND,-1)
    MM_vs_H_cal_DUMMY.Add(MM_vs_H_cal_DUMMY_RAND,-1)
    MM_vs_P_cal_DATA.Add(MM_vs_P_cal_RAND,-1)    
    MM_vs_P_hgcer_DUMMY.Add(MM_vs_P_hgcer_DUMMY_RAND,-1)
    MM_vs_P_aero_DUMMY.Add(MM_vs_P_aero_DUMMY_RAND,-1)    
    phiq_vs_t_DUMMY.Add(phiq_vs_t_DUMMY_RAND,-1)
    Q2_vs_W_DUMMY.Add(Q2_vs_W_DUMMY_RAND,-1)
    Q2_vs_t_DUMMY.Add(Q2_vs_t_DUMMY_RAND,-1)
    W_vs_t_DUMMY.Add(W_vs_t_DUMMY_RAND,-1)
    EPS_vs_t_DUMMY.Add(EPS_vs_t_DUMMY_RAND,-1)
    MM_vs_t_DUMMY.Add(MM_vs_t_DUMMY_RAND,-1)
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
    H_ph_q_DUMMY.Add(H_ph_q_DUMMY_RAND,-1)
    H_th_q_DUMMY.Add(H_th_q_DUMMY_RAND,-1)
    H_ph_recoil_DUMMY.Add(H_ph_recoil_DUMMY_RAND,-1)
    H_th_recoil_DUMMY.Add(H_th_recoil_DUMMY_RAND,-1)    
    H_Q2_DUMMY.Add(H_Q2_DUMMY_RAND,-1)
    H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
    H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
    H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)
    H_MM_DUMMY.Add(H_MM_DUMMY_RAND,-1)
    H_MM_nosub_DUMMY.Add(H_MM_nosub_DUMMY_RAND,-1)
    H_pmiss_DUMMY.Add(H_pmiss_DUMMY_RAND,-1)
    H_emiss_DUMMY.Add(H_emiss_DUMMY_RAND,-1)
    H_pmx_DUMMY.Add(H_pmx_DUMMY_RAND,-1)
    H_pmy_DUMMY.Add(H_pmy_DUMMY_RAND,-1)
    H_pmz_DUMMY.Add(H_pmz_DUMMY_RAND,-1)
    H_ct_DUMMY.Add(H_ct_DUMMY_RAND,-1)
    
    # Pion subtraction by scaling simc to peak size
    if ParticleType == "kaon":
        subDict["nWindows"] = nWindows
        subDict["phi_setting"] = phi_setting
        particle_subtraction_cuts(subDict, inpDict, SubtractedParticle, hgcer_cutg)
        
        try:
        ##############
        # HARD CODED #
        ##############
            scale_factor = H_MM_nosub_DATA.Integral(H_MM_nosub_DATA.FindBin(0.89), H_MM_nosub_DATA.FindBin(0.94))/subDict["H_MM_nosub_SUB_DATA"]\
                                    .Integral(subDict["H_MM_nosub_SUB_DATA"].FindBin(0.89), subDict["H_MM_nosub_SUB_DATA"].FindBin(0.94))
        ##############
        ##############
        ##############
        except ZeroDivisionError:
            scale_factor = 0.0

        # Apply scale factor
        subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"].Scale(scale_factor)
        subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"].Scale(scale_factor)
        subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
        subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
        subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
        subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_CoinTime_SUB_DATA"].Scale(scale_factor)
        subDict["CoinTime_vs_beta_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_beta_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_H_cer_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_H_cal_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_P_cal_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_P_hgcer_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_P_aero_SUB_DATA"].Scale(scale_factor)
        subDict["phiq_vs_t_SUB_DATA"].Scale(scale_factor)
        subDict["Q2_vs_W_SUB_DATA"].Scale(scale_factor)
        subDict["Q2_vs_t_SUB_DATA"].Scale(scale_factor)
        subDict["W_vs_t_SUB_DATA"].Scale(scale_factor)
        subDict["EPS_vs_t_SUB_DATA"].Scale(scale_factor)
        subDict["MM_vs_t_SUB_DATA"].Scale(scale_factor)
        subDict["H_ct_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssxfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssyfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssxpfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssypfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsxfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsyfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsxpfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsypfp_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssxptar_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssyptar_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsxptar_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsyptar_SUB_DATA"].Scale(scale_factor)
        subDict["H_ssdelta_SUB_DATA"].Scale(scale_factor)
        subDict["H_hsdelta_SUB_DATA"].Scale(scale_factor)
        subDict["H_ph_q_SUB_DATA"].Scale(scale_factor)
        subDict["H_th_q_SUB_DATA"].Scale(scale_factor)
        subDict["H_ph_recoil_SUB_DATA"].Scale(scale_factor)
        subDict["H_th_recoil_SUB_DATA"].Scale(scale_factor)
        subDict["H_Q2_SUB_DATA"].Scale(scale_factor)
        subDict["H_W_SUB_DATA"].Scale(scale_factor)
        subDict["H_t_SUB_DATA"].Scale(scale_factor)
        subDict["H_epsilon_SUB_DATA"].Scale(scale_factor)
        subDict["H_MM_SUB_DATA"].Scale(scale_factor)
        subDict["H_MM_nosub_SUB_DATA"].Scale(scale_factor)
        subDict["H_pmiss_SUB_DATA"].Scale(scale_factor)
        subDict["H_emiss_SUB_DATA"].Scale(scale_factor)
        subDict["H_pmx_SUB_DATA"].Scale(scale_factor)
        subDict["H_pmy_SUB_DATA"].Scale(scale_factor)
        subDict["H_pmz_SUB_DATA"].Scale(scale_factor)
        histDict["H_MM_SUB_DATA"] = subDict["H_MM_SUB_DATA"]
        histDict["H_MM_nosub_SUB_DATA"] = subDict["H_MM_nosub_SUB_DATA"]

        # Subtract pion
        P_hgcer_xAtCer_vs_yAtCer_DATA.Add(subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"],-1)
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"],-1)
        P_hgcer_xAtCer_vs_MM_DATA.Add(subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"],-1)
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"],-1)
        P_hgcer_yAtCer_vs_MM_DATA.Add(subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"],-1)
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"],-1)        
        MM_vs_CoinTime_DATA.Add(subDict["MM_vs_CoinTime_SUB_DATA"],-1)
        CoinTime_vs_beta_DATA.Add(subDict["CoinTime_vs_beta_SUB_DATA"],-1)
        MM_vs_beta_DATA.Add(subDict["MM_vs_beta_SUB_DATA"],-1)
        MM_vs_H_cer_DATA.Add(subDict["MM_vs_H_cer_SUB_DATA"],-1)
        MM_vs_H_cal_DATA.Add(subDict["MM_vs_H_cal_SUB_DATA"],-1)
        MM_vs_P_cal_DATA.Add(subDict["MM_vs_P_cal_SUB_DATA"],-1)    
        MM_vs_P_hgcer_DATA.Add(subDict["MM_vs_P_hgcer_SUB_DATA"],-1)
        MM_vs_P_aero_DATA.Add(subDict["MM_vs_P_aero_SUB_DATA"],-1)
        phiq_vs_t_DATA.Add(subDict["phiq_vs_t_SUB_DATA"],-1)
        Q2_vs_W_DATA.Add(subDict["Q2_vs_W_SUB_DATA"],-1)
        Q2_vs_t_DATA.Add(subDict["Q2_vs_t_SUB_DATA"],-1)
        W_vs_t_DATA.Add(subDict["W_vs_t_SUB_DATA"],-1)
        EPS_vs_t_DATA.Add(subDict["EPS_vs_t_SUB_DATA"],-1)
        MM_vs_t_DATA.Add(subDict["MM_vs_t_SUB_DATA"],-1)    
        H_ssxfp_DATA.Add(subDict["H_ssxfp_SUB_DATA"],-1)
        H_ssyfp_DATA.Add(subDict["H_ssyfp_SUB_DATA"],-1)
        H_ssxpfp_DATA.Add(subDict["H_ssxpfp_SUB_DATA"],-1)
        H_ssypfp_DATA.Add(subDict["H_ssypfp_SUB_DATA"],-1)
        H_hsxfp_DATA.Add(subDict["H_hsxfp_SUB_DATA"],-1)
        H_hsyfp_DATA.Add(subDict["H_hsyfp_SUB_DATA"],-1)
        H_hsxpfp_DATA.Add(subDict["H_hsxpfp_SUB_DATA"],-1)
        H_hsypfp_DATA.Add(subDict["H_hsypfp_SUB_DATA"],-1)
        H_ssxptar_DATA.Add(subDict["H_ssxptar_SUB_DATA"],-1)
        H_ssyptar_DATA.Add(subDict["H_ssyptar_SUB_DATA"],-1)
        H_hsxptar_DATA.Add(subDict["H_hsxptar_SUB_DATA"],-1)
        H_hsyptar_DATA.Add(subDict["H_hsyptar_SUB_DATA"],-1)
        H_ssdelta_DATA.Add(subDict["H_ssdelta_SUB_DATA"],-1)
        H_hsdelta_DATA.Add(subDict["H_hsdelta_SUB_DATA"],-1)
        H_ph_q_DATA.Add(subDict["H_ph_q_SUB_DATA"],-1)
        H_th_q_DATA.Add(subDict["H_th_q_SUB_DATA"],-1)
        H_ph_recoil_DATA.Add(subDict["H_ph_recoil_SUB_DATA"],-1)
        H_th_recoil_DATA.Add(subDict["H_th_recoil_SUB_DATA"],-1)
        H_Q2_DATA.Add(subDict["H_Q2_SUB_DATA"],-1)
        H_W_DATA.Add(subDict["H_W_SUB_DATA"],-1)
        H_t_DATA.Add(subDict["H_t_SUB_DATA"],-1)
        H_epsilon_DATA.Add(subDict["H_epsilon_SUB_DATA"],-1)
        H_MM_DATA.Add(subDict["H_MM_SUB_DATA"],-1)
        H_pmiss_DATA.Add(subDict["H_pmiss_SUB_DATA"],-1)
        H_emiss_DATA.Add(subDict["H_emiss_SUB_DATA"],-1)
        H_pmx_DATA.Add(subDict["H_pmx_SUB_DATA"],-1)
        H_pmy_DATA.Add(subDict["H_pmy_SUB_DATA"],-1)
        H_pmz_DATA.Add(subDict["H_pmz_SUB_DATA"],-1)
        H_ct_DATA.Add(subDict["H_ct_SUB_DATA"],-1)

    # Fit background and subtract
    inpDict["bg_tot_num_evts_{}".format(phi_setting)] = H_MM_nosub_DATA.GetEntries()
    background_fit = bg_fit(phi_setting, inpDict, H_MM_nosub_DATA)
    # RLT (4/16/2023): Commented out because they return empty sometimes, probably a TH2D vs TH1D issue
    #P_hgcer_xAtCer_vs_yAtCer_DATA.Add(background_fit[0], -1)
    #P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Add(background_fit[0], -1)
    #P_hgcer_xAtCer_vs_MM_DATA.Add(background_fit[0], -1)
    #P_hgcer_nohole_xAtCer_vs_MM_DATA.Add(background_fit[0], -1)
    #P_hgcer_yAtCer_vs_MM_DATA.Add(background_fit[0], -1)
    #P_hgcer_nohole_yAtCer_vs_MM_DATA.Add(background_fit[0], -1)
    #MM_vs_CoinTime_DATA.Add(background_fit[0], -1)
    #CoinTime_vs_beta_DATA.Add(background_fit[0], -1)
    #MM_vs_beta_DATA.Add(background_fit[0], -1)
    #MM_vs_H_cer_DATA.Add(background_fit[0], -1)
    #MM_vs_H_cal_DATA.Add(background_fit[0], -1)
    #MM_vs_P_cal_DATA.Add(background_fit[0], -1)
    #MM_vs_P_hgcer_DATA.Add(background_fit[0], -1)
    #MM_vs_P_aero_DATA.Add(background_fit[0], -1)
    #phiq_vs_t_DATA.Add(background_fit[0], -1)
    #Q2_vs_W_DATA.Add(background_fit[0], -1)
    #Q2_vs_t_DATA.Add(background_fit[0], -1)
    #W_vs_t_DATA.Add(background_fit[0], -1)
    #EPS_vs_t_DATA.Add(background_fit[0], -1)
    #MM_vs_t_DATA.Add(background_fit[0], -1)
    H_ssxfp_DATA.Add(background_fit[0], -1)
    H_ssyfp_DATA.Add(background_fit[0], -1)
    H_ssxpfp_DATA.Add(background_fit[0], -1)
    H_ssypfp_DATA.Add(background_fit[0], -1)
    H_hsxfp_DATA.Add(background_fit[0], -1)
    H_hsyfp_DATA.Add(background_fit[0], -1)
    H_hsxpfp_DATA.Add(background_fit[0], -1)
    H_hsypfp_DATA.Add(background_fit[0], -1)
    H_ssxptar_DATA.Add(background_fit[0], -1)
    H_ssyptar_DATA.Add(background_fit[0], -1)
    H_hsxptar_DATA.Add(background_fit[0], -1)
    H_hsyptar_DATA.Add(background_fit[0], -1)
    H_ssdelta_DATA.Add(background_fit[0], -1)
    H_hsdelta_DATA.Add(background_fit[0], -1)
    H_ph_q_DATA.Add(background_fit[0], -1)
    H_th_q_DATA.Add(background_fit[0], -1)
    H_ph_recoil_DATA.Add(background_fit[0], -1)
    H_th_recoil_DATA.Add(background_fit[0], -1)
    H_Q2_DATA.Add(background_fit[0], -1)
    H_W_DATA.Add(background_fit[0], -1)
    H_t_DATA.Add(background_fit[0], -1)
    H_epsilon_DATA.Add(background_fit[0], -1)
    H_MM_DATA.Add(background_fit[0], -1)
    H_pmiss_DATA.Add(background_fit[0], -1)
    H_emiss_DATA.Add(background_fit[0], -1)
    H_pmx_DATA.Add(background_fit[0], -1)
    H_pmy_DATA.Add(background_fit[0], -1)
    H_pmz_DATA.Add(background_fit[0], -1)
    H_ct_DATA.Add(background_fit[0], -1)
        
    histDict["InFile_DATA"] = InFile_DATA
    histDict["InFile_DUMMY"] = InFile_DUMMY
    histDict["phi_setting"] = phi_setting
    histDict["pid_text"] = pid_text
    histDict["runNums"] = runNums.split(' ')
    histDict["nWindows"] = nWindows
    histDict["H_hsdelta_DUMMY"] =     remove_bad_bins(H_hsdelta_DUMMY)
    histDict["H_hsxptar_DUMMY"] =     remove_bad_bins(H_hsxptar_DUMMY)
    histDict["H_hsyptar_DUMMY"] =     remove_bad_bins(H_hsyptar_DUMMY)
    histDict["H_ssxfp_DUMMY"] =     remove_bad_bins(H_ssxfp_DUMMY)
    histDict["H_ssyfp_DUMMY"] =     remove_bad_bins(H_ssyfp_DUMMY)
    histDict["H_ssxpfp_DUMMY"] =     remove_bad_bins(H_ssxpfp_DUMMY)
    histDict["H_ssypfp_DUMMY"] =     remove_bad_bins(H_ssypfp_DUMMY)
    histDict["H_hsxfp_DUMMY"] =     remove_bad_bins(H_hsxfp_DUMMY)
    histDict["H_hsyfp_DUMMY"] =     remove_bad_bins(H_hsyfp_DUMMY)
    histDict["H_hsxpfp_DUMMY"] =     remove_bad_bins(H_hsxpfp_DUMMY)
    histDict["H_hsypfp_DUMMY"] =     remove_bad_bins(H_hsypfp_DUMMY)
    histDict["H_ssdelta_DUMMY"] =     remove_bad_bins(H_ssdelta_DUMMY)
    histDict["H_ssxptar_DUMMY"] =     remove_bad_bins(H_ssxptar_DUMMY)
    histDict["H_ssyptar_DUMMY"] =     remove_bad_bins(H_ssyptar_DUMMY)
    histDict["H_q_DUMMY"] =     remove_bad_bins(H_q_DUMMY)
    histDict["H_Q2_DUMMY"] =     remove_bad_bins(H_Q2_DUMMY)
    histDict["H_t_DUMMY"] =     remove_bad_bins(H_t_DUMMY)
    histDict["H_epsilon_DUMMY"] =     remove_bad_bins(H_epsilon_DUMMY)
    histDict["H_MM_DUMMY"] =     remove_bad_bins(H_MM_DUMMY)
    histDict["H_MM_nosub_DUMMY"] =     remove_bad_bins(H_MM_nosub_DUMMY)
    histDict["H_th_DUMMY"] =     remove_bad_bins(H_th_DUMMY)
    histDict["H_ph_DUMMY"] =     remove_bad_bins(H_ph_DUMMY)
    histDict["H_ph_q_DUMMY"] =     remove_bad_bins(H_ph_q_DUMMY)
    histDict["H_th_q_DUMMY"] =     remove_bad_bins(H_th_q_DUMMY)
    histDict["H_ph_recoil_DUMMY"] =     remove_bad_bins(H_ph_recoil_DUMMY)
    histDict["H_th_recoil_DUMMY"] =     remove_bad_bins(H_th_recoil_DUMMY)
    histDict["H_pmiss_DUMMY"] =     remove_bad_bins(H_pmiss_DUMMY)
    histDict["H_emiss_DUMMY"] =     remove_bad_bins(H_emiss_DUMMY)
    histDict["H_pmx_DUMMY"] =     remove_bad_bins(H_pmx_DUMMY)
    histDict["H_pmy_DUMMY"] =     remove_bad_bins(H_pmy_DUMMY)
    histDict["H_pmz_DUMMY"] =     remove_bad_bins(H_pmz_DUMMY)
    histDict["H_W_DUMMY"] =     remove_bad_bins(H_W_DUMMY)
    histDict["H_ct_DUMMY"] =     remove_bad_bins(H_ct_DUMMY)
    histDict["MM_vs_CoinTime_DUMMY"] = remove_bad_bins(MM_vs_CoinTime_DUMMY)
    histDict["CoinTime_vs_beta_DUMMY"] = remove_bad_bins(CoinTime_vs_beta_DUMMY)
    histDict["MM_vs_beta_DUMMY"] = remove_bad_bins(MM_vs_beta_DUMMY)
    histDict["phiq_vs_t_DUMMY"] = remove_bad_bins(phiq_vs_t_DUMMY)
    histDict["polar_phiq_vs_t_DUMMY"] = polar_phiq_vs_t_DUMMY
    histDict["H_hsdelta_DATA"] =     remove_bad_bins(H_hsdelta_DATA)
    histDict["H_hsxptar_DATA"] =     remove_bad_bins(H_hsxptar_DATA)
    histDict["H_hsyptar_DATA"] =     remove_bad_bins(H_hsyptar_DATA)
    histDict["H_ssxfp_DATA"] =     remove_bad_bins(H_ssxfp_DATA)
    histDict["H_ssyfp_DATA"] =     remove_bad_bins(H_ssyfp_DATA)
    histDict["H_ssxpfp_DATA"] =     remove_bad_bins(H_ssxpfp_DATA)
    histDict["H_ssypfp_DATA"] =     remove_bad_bins(H_ssypfp_DATA)
    histDict["H_hsxfp_DATA"] =     remove_bad_bins(H_hsxfp_DATA)
    histDict["H_hsyfp_DATA"] =     remove_bad_bins(H_hsyfp_DATA)
    histDict["H_hsxpfp_DATA"] =     remove_bad_bins(H_hsxpfp_DATA)
    histDict["H_hsypfp_DATA"] =     remove_bad_bins(H_hsypfp_DATA)
    histDict["H_ssdelta_DATA"] =     remove_bad_bins(H_ssdelta_DATA)
    histDict["H_ssxptar_DATA"] =     remove_bad_bins(H_ssxptar_DATA)
    histDict["H_ssyptar_DATA"] =     remove_bad_bins(H_ssyptar_DATA)
    histDict["H_q_DATA"] =     remove_bad_bins(H_q_DATA)
    histDict["H_Q2_DATA"] =     remove_bad_bins(H_Q2_DATA)
    histDict["H_t_DATA"] =     remove_bad_bins(H_t_DATA)
    histDict["H_epsilon_DATA"] =     remove_bad_bins(H_epsilon_DATA)
    histDict["H_MM_DATA"] =     remove_bad_bins(H_MM_DATA)
    histDict["H_MM_nosub_DATA"] =     remove_bad_bins(H_MM_nosub_DATA)
    histDict["H_th_DATA"] =     remove_bad_bins(H_th_DATA)
    histDict["H_ph_DATA"] =     remove_bad_bins(H_ph_DATA)
    histDict["H_ph_q_DATA"] =     remove_bad_bins(H_ph_q_DATA)
    histDict["H_th_q_DATA"] =     remove_bad_bins(H_th_q_DATA)
    histDict["H_ph_recoil_DATA"] =     remove_bad_bins(H_ph_recoil_DATA)
    histDict["H_th_recoil_DATA"] =     remove_bad_bins(H_th_recoil_DATA)
    histDict["H_pmiss_DATA"] =     remove_bad_bins(H_pmiss_DATA)
    histDict["H_emiss_DATA"] =     remove_bad_bins(H_emiss_DATA)
    histDict["H_pmx_DATA"] =     remove_bad_bins(H_pmx_DATA)
    histDict["H_pmy_DATA"] =     remove_bad_bins(H_pmy_DATA)
    histDict["H_pmz_DATA"] =     remove_bad_bins(H_pmz_DATA)
    histDict["H_W_DATA"] =     remove_bad_bins(H_W_DATA)
    histDict["H_ct_DATA"] =     remove_bad_bins(H_ct_DATA)
    histDict["H_cal_etottracknorm_DATA"] =     remove_bad_bins(H_cal_etottracknorm_DATA)
    histDict["H_cer_npeSum_DATA"] =     remove_bad_bins(H_cer_npeSum_DATA)
    histDict["P_cal_etottracknorm_DATA"] =     remove_bad_bins(P_cal_etottracknorm_DATA)
    histDict["P_hgcer_npeSum_DATA"] =     remove_bad_bins(P_hgcer_npeSum_DATA)
    histDict["P_aero_npeSum_DATA"] =     remove_bad_bins(P_aero_npeSum_DATA)
    histDict["MM_vs_CoinTime_DATA"] = remove_bad_bins(MM_vs_CoinTime_DATA)
    histDict["CoinTime_vs_beta_DATA"] = remove_bad_bins(CoinTime_vs_beta_DATA)
    histDict["MM_vs_beta_DATA"] = remove_bad_bins(MM_vs_beta_DATA)
    histDict["phiq_vs_t_DATA"] = remove_bad_bins(phiq_vs_t_DATA)
    histDict["polar_phiq_vs_t_DATA"] = polar_phiq_vs_t_DATA
    histDict["Q2_vs_W_DATA"] = remove_bad_bins(Q2_vs_W_DATA)
    histDict["Q2_vs_t_DATA"] = remove_bad_bins(Q2_vs_t_DATA)
    histDict["W_vs_t_DATA"] = remove_bad_bins(W_vs_t_DATA)
    histDict["EPS_vs_t_DATA"] = remove_bad_bins(EPS_vs_t_DATA)
    histDict["MM_vs_t_DATA"] = remove_bad_bins(MM_vs_t_DATA)
    histDict["MM_vs_H_cer_DATA"] = remove_bad_bins(MM_vs_H_cer_DATA)
    histDict["MM_vs_H_cal_DATA"] = remove_bad_bins(MM_vs_H_cal_DATA)
    histDict["MM_vs_P_cal_DATA"] = remove_bad_bins(MM_vs_P_cal_DATA)
    histDict["MM_vs_P_hgcer_DATA"] = remove_bad_bins(MM_vs_P_hgcer_DATA)
    histDict["MM_vs_P_aero_DATA"] = remove_bad_bins(MM_vs_P_aero_DATA)
    histDict["NumEvts_MM_DUMMY"] = int(H_MM_DUMMY.Integral())
    histDict["NumEvts_MM_DATA"] = int(H_MM_DATA.Integral())
    
    ###
    # CT plots
    ct = TCanvas()
    l_ct = TLegend(0.115,0.65,0.33,0.95)
    l_ct.SetTextSize(0.0235)
    H_ct_DATA.SetLineColor(1)
    H_ct_DATA.Draw("same, HIST")
    l_ct.AddEntry(H_ct_DATA,"{}".format(ParticleType.capitalize()))
    l_ct.Draw()

    ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+'(')

    ###
    # Q2 plots    
    CQ2 = TCanvas()

    histDict["H_Q2_DATA"].SetLineColor(1)
    histDict["H_Q2_DATA"].Draw("same, E1")

    CQ2.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # W plots    
    CW = TCanvas()

    histDict["H_W_DATA"].SetLineColor(1)
    histDict["H_W_DATA"].Draw("same, E1")

    CW.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))    
    
    ###
    # MM plots    
    CMM = TCanvas()

    histDict["H_MM_DATA"].SetLineColor(1)
    histDict["H_MM_DATA"].Draw("same, E1")

    CMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
        
    ###
    # MM sub plots    
    CMMsub = TCanvas()

    histDict["H_MM_nosub_DATA"].SetLineColor(1)
    histDict["H_MM_nosub_DATA"].SetFillStyle(3001)  # Set fill style to dots
    histDict["H_MM_nosub_DATA"].SetFillColor(kBlack)  # Set fill color to black
    #histDict["H_MM_nosub_DATA"].Draw("same, E1")
    histDict["H_MM_nosub_DATA"].Draw("hist same")
    if ParticleType == "kaon":        
        histDict["H_MM_nosub_SUB_DATA"].SetLineColor(2)
        histDict["H_MM_nosub_SUB_DATA"].Draw("same, E1")
    background_fit[0].SetLineColor(3)
    background_fit[0].Draw("same")

    CMMsub.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    
    ###
    # t-Phi plots        
    Cpht_data = ROOT.TCanvas()

    # Create the polar plot using the function
    polar_plot = create_polar_plot(histDict["polar_phiq_vs_t_DATA"])
    # Draw the plot
    polar_plot.Draw("AP")
    
    Cpht_data.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # t plots            
    Ct = TCanvas()
    l_t = TLegend(0.115,0.45,0.33,0.95)
    l_t.SetTextSize(0.0135)

    histDict["H_t_DATA"].SetLineColor(1)
    l_t.AddEntry(histDict["H_t_DATA"],histDict["phi_setting"])
    histDict["H_t_DATA"].Draw("same, E1")

    Ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

    ###
    # phi plots            
    Cphi = TCanvas()
    l_phi = TLegend(0.115,0.45,0.33,0.95)
    l_phi.SetTextSize(0.0135)
    histDict["H_ph_q_DATA"].SetLineColor(1)
    l_phi.AddEntry(histDict["H_ph_q_DATA"],histDict["phi_setting"])
    histDict["H_ph_q_DATA"].Draw("same, E1")    

    Cphi.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    
    ###
    # PID Plots
    c_pid = TCanvas()

    c_pid.Divide(2,3)

    c_pid.cd(1)
    gPad.SetLogy()

    H_cal_etottracknorm_DATA.SetLineColor(1)
    H_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(2)
    gPad.SetLogy()

    H_cer_npeSum_DATA.SetLineColor(1)
    H_cer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(3)
    gPad.SetLogy()

    P_cal_etottracknorm_DATA.SetLineColor(1)
    P_cal_etottracknorm_DATA.Draw("same, HIST")

    c_pid.cd(4)
    gPad.SetLogy()

    P_hgcer_npeSum_DATA.SetLineColor(1)
    P_hgcer_npeSum_DATA.Draw("same, HIST")

    c_pid.cd(5)
    gPad.SetLogy()

    P_aero_npeSum_DATA.SetLineColor(1)
    P_aero_npeSum_DATA.Draw("same, HIST")

    c_pid.Draw()

    if ParticleType == "kaon":
        c_pid.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))
    else:
        c_pid.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')

    if ParticleType == "kaon":
        
        ##
        # HGCer Hole Plots
        c_hgcervsMM = TCanvas()

        c_hgcervsMM.Divide(2,2)

        c_hgcervsMM.cd(1)
        P_hgcer_xAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(2)
        P_hgcer_nohole_xAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(3)
        P_hgcer_yAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_yAtCer_vs_MM_DATA.Draw("colz")

        c_hgcervsMM.cd(4)
        P_hgcer_nohole_yAtCer_vs_MM_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_yAtCer_vs_MM_DATA.Draw("colz")    

        c_hgcervsMM.Draw()

        c_hgcervsMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType)))

        ##
        # HGCer Hole Plots
        c_hgcer_hole = TCanvas()

        c_hgcer_hole.Divide(2,2)

        c_hgcer_hole.cd(1)
        P_hgcer_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_yAtCer_DATA.Draw("colz")

        c_hgcer_hole.cd(2)
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Draw("colz")

        c_hgcer_hole.cd(3)
        P_hgcer_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_xAtCer_vs_yAtCer_DATA.Draw("colz")
        hgcer_cutg.SetLineColor(7)
        hgcer_cutg.Draw("same")

        c_hgcer_hole.cd(4)
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.SetMinimum(1e-6) # Remove color of empty bins
        P_hgcer_nohole_xAtCer_vs_yAtCer_DATA.Draw("colz")
        hgcer_cutg.SetLineColor(7)
        hgcer_cutg.Draw("same")    

        c_hgcer_hole.Draw()

        c_hgcer_hole.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_rand_sub_".format(phi_setting,ParticleType))+')')      

    return histDict
