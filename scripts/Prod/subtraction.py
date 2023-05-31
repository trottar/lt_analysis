#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-05-30 21:37:43 trottar"
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
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
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

def defineHists(phi_setting, inpDict):    

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

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    '''
    Importing diamond cut script
    '''

    from DiamondPlot_all import DiamondPlot

    Q2Val = float(Q2.replace("p","."))
    WVal = float(W.replace("p","."))

    Q2min = Q2Val - (2/7)*Q2Val # Minimum value of Q2 on the Q2 vs W plot
    Q2max = Q2Val + (2/7)*Q2Val # Maximum value of Q2 on the Q2 vs W plot
    Wmin = WVal - (2/7)*WVal # min y-range for Q2vsW plot
    Wmax = WVal + (2/7)*WVal # max y-range for Q2vsW plot

    ################################################################################################################################################    

    ################################################################################################################################################
    # Call diamond cut script

    paramDict = DiamondPlot(ParticleType, Q2Val, Q2min, Q2max, WVal, Wmin, Wmax, phi_setting, tmin, tmax)

    a1 = paramDict["a1"]
    b1 = paramDict["b1"]
    a2 = paramDict["a2"]
    b2 = paramDict["b2"]
    a3 = paramDict["a3"]
    b3 = paramDict["b3"]
    a4 = paramDict["a4"]
    b4 = paramDict["b4"]

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Check particle type

    if ParticleType == "kaon":


        ################################################################################################################################################
        # Define simc root file trees of interest

        # Names don't match so need to do some string rearrangement
        InSIMCFilename = "Prod_Coin_{}.root".format(kinematics[0]+phi_setting.lower()+"_"+kinematics[1])
        rootFileSimc = OUTPATH+"/"+InSIMCFilename
        if not os.path.isfile(rootFileSimc):
            print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
            return histDict

        InFile_SIMC = ROOT.TFile.Open(rootFileSimc, "OPEN")

        TBRANCH_SIMC  = InFile_SIMC.Get("h10")

        ###############################################################################################################################################

        # Grabs simc number of events and normalizaton factor
        simc_hist = rootFileSimc.replace('.root','.hist')
        f_simc = open(simc_hist)
        for line in f_simc:
            #print(line)
            if "Ngen" in line:
                val = line.split("=")
                simc_nevents = int(val[1])
            if "normfac" in line:
                val = line.split("=")
                simc_normfactor = float(val[1])
        if 'simc_nevents' and 'simc_normfactor' not in locals():
            print("\n\nERROR: Invalid simc hist file %s\n\n" % simc_hist)
            sys.exit(1)
        f_simc.close()

        ################################################################################################################################################
        # Define data root file trees of interest

        rootFileData = OUTPATH + "/" + "kaon" + "_" + InDATAFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileData):
            print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
            return histDict

        InFile_DATA = ROOT.TFile.Open(rootFileData, "OPEN")

        TBRANCH_DATA  = InFile_DATA.Get("Cut_Kaon_Events_prompt_RF")

        TBRANCH_RAND  = InFile_DATA.Get("Cut_Kaon_Events_rand_RF")

        ################################################################################################################################################
        # Define dummy root file trees of interest

        rootFileDummy = OUTPATH + "/" + "kaon" + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileDummy):
            print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
            return histDict

        InFile_DUMMY = ROOT.TFile.Open(rootFileDummy, "OPEN")  

        TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_Kaon_Events_prompt_RF")

        TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_Kaon_Events_rand_RF")

        ################################################################################################################################################
        # Define data root file trees of interest

        rootFileSubPionData = OUTPATH + "/" + "pion" + "_" + InDATAFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileSubPionData):
            print("\n\nERROR: No data file found called {}\n\n".format(rootFileSubPionData))
            return histDict  

        InFile_SUBPION_DATA = ROOT.TFile.Open(rootFileSubPionData, "OPEN")

        TBRANCH_SUBPION_DATA  = InFile_SUBPION_DATA.Get("Cut_Pion_Events_prompt_noRF")

        TBRANCH_SUBPION_RAND  = InFile_SUBPION_DATA.Get("Cut_Pion_Events_rand_noRF")

        ################################################################################################################################################
        # Define dummy root file trees of interest

        rootFileSubPionDummy = OUTPATH + "/" + "pion" + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileSubPionDummy):
            print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileSubPionDummy))
            return histDict

        InFile_SUBPION_DUMMY = ROOT.TFile.Open(rootFileSubPionDummy, "OPEN")  

        TBRANCH_SUBPION_DUMMY  = InFile_SUBPION_DUMMY.Get("Cut_Pion_Events_prompt_noRF")

        TBRANCH_SUBPION_DUMMY_RAND  = InFile_SUBPION_DUMMY.Get("Cut_Pion_Events_rand_noRF")

        ################################################################################################################################################
        # Define data root file trees of interest

        rootFileSubProtonData = OUTPATH + "/" + "proton" + "_" + InDATAFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileSubProtonData):
            print("\n\nERROR: No data file found called {}\n\n".format(rootFileSubProtonData))
            return histDict  

        InFile_SUBPROTON_DATA = ROOT.TFile.Open(rootFileSubProtonData, "OPEN")

        TBRANCH_SUBPROTON_DATA  = InFile_SUBPROTON_DATA.Get("Cut_Proton_Events_prompt_noRF")

        TBRANCH_SUBPROTON_RAND  = InFile_SUBPROTON_DATA.Get("Cut_Proton_Events_rand_noRF")

        ################################################################################################################################################
        # Define dummy root file trees of interest

        rootFileSubProtonDummy = OUTPATH + "/" + "proton" + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileSubProtonDummy):
            print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileSubProtonDummy))
            return histDict

        InFile_SUBPROTON_DUMMY = ROOT.TFile.Open(rootFileSubProtonDummy, "OPEN")  

        TBRANCH_SUBPROTON_DUMMY  = InFile_SUBPROTON_DUMMY.Get("Cut_Proton_Events_prompt_noRF")

        TBRANCH_SUBPROTON_DUMMY_RAND  = InFile_SUBPROTON_DUMMY.Get("Cut_Proton_Events_rand_noRF")
        
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

        ################################################################################################################################################
        # Grab and calculate efficiency

        sys.path.append('../')
        from getDataTable import calculate_effError

        tot_effError_data = [calculate_effError(run,efficiency_table) for run in runNums.split(' ')]
        #print(InData_efficiency)
        #print(tot_effError_data)
        eff_errProp_data = sum(tot_effError_data) # Error propagation for addition

        print("\n\nTotal Data Efficiency Uncertainty =",eff_errProp_data)

        # Define total efficiency vs run number plots
        G_data_eff = ROOT.TGraphErrors(len(InData_efficiency.split(' ')), np.array([float(x) for x in runNums.split(' ')]),np.array([float(x) for x in InData_efficiency.split(' ')]),np.array([0]*len(tot_effError_data)),np.array(tot_effError_data)*np.array([float(x) for x in InData_efficiency.split(' ')]))

        ###############################################################################################################################################
        # Grab windows for random subtraction

        # Section for grabing Prompt/Random selection parameters from PARAM file
        PARAMPATH = "%s/DB/PARAM" % UTILPATH
        print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], LTANAPATH))
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

        H_hsdelta_SIMC  = ROOT.TH1D("H_hsdelta_SIMC","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SIMC  = ROOT.TH1D("H_hsxptar_SIMC","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SIMC  = ROOT.TH1D("H_hsyptar_SIMC","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SIMC    = ROOT.TH1D("H_ssxfp_SIMC","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SIMC    = ROOT.TH1D("H_ssyfp_SIMC","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SIMC   = ROOT.TH1D("H_ssxpfp_SIMC","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SIMC   = ROOT.TH1D("H_ssypfp_SIMC","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SIMC    = ROOT.TH1D("H_hsxfp_SIMC","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SIMC    = ROOT.TH1D("H_hsyfp_SIMC","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SIMC   = ROOT.TH1D("H_hsxpfp_SIMC","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SIMC   = ROOT.TH1D("H_hsypfp_SIMC","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SIMC  = ROOT.TH1D("H_ssdelta_SIMC","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SIMC  = ROOT.TH1D("H_ssxptar_SIMC","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SIMC  = ROOT.TH1D("H_ssyptar_SIMC","SHMS yptar", 100, -0.04, 0.04)
        H_q_SIMC        = ROOT.TH1D("H_q_SIMC","q", 100, 0.0, 10.0)
        H_Q2_SIMC       = ROOT.TH1D("H_Q2_SIMC","Q2", 100, Q2min, Q2max)
        H_W_SIMC  = ROOT.TH1D("H_W_SIMC","W ", 100, Wmin, Wmax)
        H_t_SIMC       = ROOT.TH1D("H_t_SIMC","-t", 100, -1.0, 1.5)  
        H_epsilon_SIMC  = ROOT.TH1D("H_epsilon_SIMC","epsilon", 100, 0., 1.0)
        H_MM_SIMC  = ROOT.TH1D("H_MM_SIMC","MM_{K}", 100, 0.0, 1.5)
        H_th_SIMC  = ROOT.TH1D("H_th_SIMC","X' tar", 100, -0.1, 0.1)
        H_ph_SIMC  = ROOT.TH1D("H_ph_SIMC","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SIMC  = ROOT.TH1D("H_ph_q_SIMC","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SIMC  = ROOT.TH1D("H_th_q_SIMC","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SIMC  = ROOT.TH1D("H_ph_recoil_SIMC","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SIMC  = ROOT.TH1D("H_th_recoil_SIMC","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SIMC  = ROOT.TH1D("H_pmiss_SIMC","pmiss", 100, 0.0, 10.0)
        H_emiss_SIMC  = ROOT.TH1D("H_emiss_SIMC","emiss", 100, 0.0, 10.0)
        H_pmx_SIMC  = ROOT.TH1D("H_pmx_SIMC","pmx", 100, -10.0, 10.0)
        H_pmy_SIMC  = ROOT.TH1D("H_pmy_SIMC","pmy ", 100, -10.0, 10.0)
        H_pmz_SIMC  = ROOT.TH1D("H_pmz_SIMC","pmz", 100, -10.0, 10.0)
        H_yield_SIMC = ROOT.TH1D("H_yield_SIMC", "Simc Yield", NumtBins*NumPhiBins, 0, 1.0)

        H_hsdelta_DATA  = ROOT.TH1D("H_hsdelta_DATA","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_DATA  = ROOT.TH1D("H_hsxptar_DATA","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_DATA  = ROOT.TH1D("H_hsyptar_DATA","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_DATA    = ROOT.TH1D("H_ssxfp_DATA","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_DATA    = ROOT.TH1D("H_ssyfp_DATA","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_DATA   = ROOT.TH1D("H_ssxpfp_DATA","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_DATA   = ROOT.TH1D("H_ssypfp_DATA","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_DATA    = ROOT.TH1D("H_hsxfp_DATA","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_DATA    = ROOT.TH1D("H_hsyfp_DATA","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_DATA   = ROOT.TH1D("H_hsxpfp_DATA","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_DATA   = ROOT.TH1D("H_hsypfp_DATA","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_DATA  = ROOT.TH1D("H_ssdelta_DATA","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_DATA  = ROOT.TH1D("H_ssxptar_DATA","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_DATA  = ROOT.TH1D("H_ssyptar_DATA","SHMS yptar", 100, -0.04, 0.04)
        H_q_DATA        = ROOT.TH1D("H_q_DATA","q", 100, 0.0, 10.0)
        H_Q2_DATA       = ROOT.TH1D("H_Q2_DATA","Q2", 100, Q2min, Q2max)
        H_W_DATA  = ROOT.TH1D("H_W_DATA","W ", 100, Wmin, Wmax)
        H_t_DATA       = ROOT.TH1D("H_t_DATA","-t", 100, -1.0, 1.5)  
        H_epsilon_DATA  = ROOT.TH1D("H_epsilon_DATA","epsilon", 100, 0., 1.0)
        H_MM_DATA  = ROOT.TH1D("H_MM_DATA","MM_{K}", 100, 0.0, 1.5)
        H_th_DATA  = ROOT.TH1D("H_th_DATA","X' tar", 100, -0.1, 0.1)
        H_ph_DATA  = ROOT.TH1D("H_ph_DATA","Y' tar", 100, -0.1, 0.1)
        H_ph_q_DATA  = ROOT.TH1D("H_ph_q_DATA","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_DATA  = ROOT.TH1D("H_th_q_DATA","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_DATA  = ROOT.TH1D("H_ph_recoil_DATA","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_DATA  = ROOT.TH1D("H_th_recoil_DATA","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_DATA  = ROOT.TH1D("H_pmiss_DATA","pmiss", 100, 0.0, 10.0)
        H_emiss_DATA  = ROOT.TH1D("H_emiss_DATA","emiss", 100, 0.0, 10.0)
        H_pmx_DATA  = ROOT.TH1D("H_pmx_DATA","pmx", 100, -10.0, 10.0)
        H_pmy_DATA  = ROOT.TH1D("H_pmy_DATA","pmy ", 100, -10.0, 10.0)
        H_pmz_DATA  = ROOT.TH1D("H_pmz_DATA","pmz", 100, -10.0, 10.0)
        H_ct_DATA = ROOT.TH1D("H_ct_DATA", "Electron-Kaon CTime", 100, -10, 10)
        H_cal_etottracknorm_DATA = ROOT.TH1D("H_cal_etottracknorm_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        H_cer_npeSum_DATA = ROOT.TH1D("H_cer_npeSum_DATA", "HMS Cer Npe Sum", 100, 0, 30)
        P_cal_etottracknorm_DATA = ROOT.TH1D("P_cal_etottracknorm_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
        P_hgcer_npeSum_DATA = ROOT.TH1D("P_hgcer_npeSum_DATA", "SHMS HGCer Npe Sum", 100, 0, 50)
        P_aero_npeSum_DATA = ROOT.TH1D("P_aero_npeSum_DATA", "SHMS Aero Npe Sum", 100, 0, 50)

        H_hsdelta_DUMMY  = ROOT.TH1D("H_hsdelta_DUMMY","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_DUMMY  = ROOT.TH1D("H_hsxptar_DUMMY","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_DUMMY  = ROOT.TH1D("H_hsyptar_DUMMY","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_DUMMY    = ROOT.TH1D("H_ssxfp_DUMMY","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_DUMMY    = ROOT.TH1D("H_ssyfp_DUMMY","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_DUMMY   = ROOT.TH1D("H_ssxpfp_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_DUMMY   = ROOT.TH1D("H_ssypfp_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_DUMMY    = ROOT.TH1D("H_hsxfp_DUMMY","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_DUMMY    = ROOT.TH1D("H_hsyfp_DUMMY","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_DUMMY   = ROOT.TH1D("H_hsxpfp_DUMMY","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_DUMMY   = ROOT.TH1D("H_hsypfp_DUMMY","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_DUMMY  = ROOT.TH1D("H_ssdelta_DUMMY","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_DUMMY  = ROOT.TH1D("H_ssxptar_DUMMY","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_DUMMY  = ROOT.TH1D("H_ssyptar_DUMMY","SHMS yptar", 100, -0.04, 0.04)
        H_q_DUMMY        = ROOT.TH1D("H_q_DUMMY","q", 100, 0.0, 10.0)
        H_Q2_DUMMY       = ROOT.TH1D("H_Q2_DUMMY","Q2", 100, Q2min, Q2max)
        H_W_DUMMY  = ROOT.TH1D("H_W_DUMMY","W ", 100, Wmin, Wmax)
        H_t_DUMMY       = ROOT.TH1D("H_t_DUMMY","-t", 100, -1.0, 1.5)  
        H_epsilon_DUMMY  = ROOT.TH1D("H_epsilon_DUMMY","epsilon", 100, 0., 1.0)
        H_MM_DUMMY  = ROOT.TH1D("H_MM_DUMMY","MM_{K}", 100, 0.0, 1.5)
        H_th_DUMMY  = ROOT.TH1D("H_th_DUMMY","X' tar", 100, -0.1, 0.1)
        H_ph_DUMMY  = ROOT.TH1D("H_ph_DUMMY","Y' tar", 100, -0.1, 0.1)
        H_ph_q_DUMMY  = ROOT.TH1D("H_ph_q_DUMMY","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_DUMMY  = ROOT.TH1D("H_th_q_DUMMY","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_DUMMY  = ROOT.TH1D("H_ph_recoil_DUMMY","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_DUMMY  = ROOT.TH1D("H_th_recoil_DUMMY","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_DUMMY  = ROOT.TH1D("H_pmiss_DUMMY","pmiss", 100, 0.0, 10.0)
        H_emiss_DUMMY  = ROOT.TH1D("H_emiss_DUMMY","emiss", 100, 0.0, 10.0)
        H_pmx_DUMMY  = ROOT.TH1D("H_pmx_DUMMY","pmx", 100, -10.0, 10.0)
        H_pmy_DUMMY  = ROOT.TH1D("H_pmy_DUMMY","pmy ", 100, -10.0, 10.0)
        H_pmz_DUMMY  = ROOT.TH1D("H_pmz_DUMMY","pmz", 100, -10.0, 10.0)
        H_ct_DUMMY = ROOT.TH1D("H_ct_DUMMY", "Electron-Kaon CTime", 100, -10, 10)

        H_hsdelta_RAND  = ROOT.TH1D("H_hsdelta_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_RAND  = ROOT.TH1D("H_hsxptar_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_RAND  = ROOT.TH1D("H_hsyptar_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_RAND    = ROOT.TH1D("H_ssxfp_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_RAND    = ROOT.TH1D("H_ssyfp_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_RAND   = ROOT.TH1D("H_ssxpfp_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_RAND   = ROOT.TH1D("H_ssypfp_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_RAND    = ROOT.TH1D("H_hsxfp_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_RAND    = ROOT.TH1D("H_hsyfp_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_RAND   = ROOT.TH1D("H_hsxpfp_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_RAND   = ROOT.TH1D("H_hsypfp_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_RAND  = ROOT.TH1D("H_ssdelta_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_RAND  = ROOT.TH1D("H_ssxptar_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_RAND  = ROOT.TH1D("H_ssyptar_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_RAND        = ROOT.TH1D("H_q_RAND","q", 100, 0.0, 10.0)
        H_Q2_RAND       = ROOT.TH1D("H_Q2_RAND","Q2", 100, Q2min, Q2max)
        H_W_RAND  = ROOT.TH1D("H_W_RAND","W ", 100, Wmin, Wmax)
        H_t_RAND       = ROOT.TH1D("H_t_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_RAND  = ROOT.TH1D("H_epsilon_RAND","epsilon", 100, 0., 1.0)
        H_MM_RAND  = ROOT.TH1D("H_MM_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_RAND  = ROOT.TH1D("H_th_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_RAND  = ROOT.TH1D("H_ph_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_RAND  = ROOT.TH1D("H_ph_q_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_RAND  = ROOT.TH1D("H_th_q_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_RAND  = ROOT.TH1D("H_ph_recoil_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_RAND  = ROOT.TH1D("H_th_recoil_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_RAND  = ROOT.TH1D("H_pmiss_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_RAND  = ROOT.TH1D("H_emiss_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_RAND  = ROOT.TH1D("H_pmx_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_RAND  = ROOT.TH1D("H_pmy_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_RAND  = ROOT.TH1D("H_pmz_RAND","pmz", 100, -10.0, 10.0)
        H_ct_RAND = ROOT.TH1D("H_ct_RAND", "Electron-Kaon CTime", 100, -10, 10)

        H_hsdelta_DUMMY_RAND  = ROOT.TH1D("H_hsdelta_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_DUMMY_RAND  = ROOT.TH1D("H_hsxptar_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_DUMMY_RAND  = ROOT.TH1D("H_hsyptar_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_DUMMY_RAND    = ROOT.TH1D("H_ssxfp_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_DUMMY_RAND    = ROOT.TH1D("H_ssyfp_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_DUMMY_RAND   = ROOT.TH1D("H_ssxpfp_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_DUMMY_RAND   = ROOT.TH1D("H_ssypfp_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_DUMMY_RAND    = ROOT.TH1D("H_hsxfp_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_DUMMY_RAND    = ROOT.TH1D("H_hsyfp_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_DUMMY_RAND   = ROOT.TH1D("H_hsxpfp_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_DUMMY_RAND   = ROOT.TH1D("H_hsypfp_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_DUMMY_RAND  = ROOT.TH1D("H_ssdelta_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_DUMMY_RAND  = ROOT.TH1D("H_ssxptar_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_DUMMY_RAND  = ROOT.TH1D("H_ssyptar_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_DUMMY_RAND        = ROOT.TH1D("H_q_DUMMY_RAND","q", 100, 0.0, 10.0)
        H_Q2_DUMMY_RAND       = ROOT.TH1D("H_Q2_DUMMY_RAND","Q2", 100, Q2min, Q2max)
        H_W_DUMMY_RAND  = ROOT.TH1D("H_W_DUMMY_RAND","W ", 100, Wmin, Wmax)
        H_t_DUMMY_RAND       = ROOT.TH1D("H_t_DUMMY_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_DUMMY_RAND  = ROOT.TH1D("H_epsilon_DUMMY_RAND","epsilon", 100, 0., 1.0)
        H_MM_DUMMY_RAND  = ROOT.TH1D("H_MM_DUMMY_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_DUMMY_RAND  = ROOT.TH1D("H_th_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_DUMMY_RAND  = ROOT.TH1D("H_ph_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_DUMMY_RAND  = ROOT.TH1D("H_ph_q_DUMMY_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_DUMMY_RAND  = ROOT.TH1D("H_th_q_DUMMY_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_DUMMY_RAND  = ROOT.TH1D("H_ph_recoil_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_DUMMY_RAND  = ROOT.TH1D("H_th_recoil_DUMMY_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_DUMMY_RAND  = ROOT.TH1D("H_pmiss_DUMMY_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_DUMMY_RAND  = ROOT.TH1D("H_emiss_DUMMY_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_DUMMY_RAND  = ROOT.TH1D("H_pmx_DUMMY_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_DUMMY_RAND  = ROOT.TH1D("H_pmy_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_DUMMY_RAND  = ROOT.TH1D("H_pmz_DUMMY_RAND","pmz", 100, -10.0, 10.0)
        H_ct_DUMMY_RAND = ROOT.TH1D("H_ct_DUMMY_RAND", "Electron-Kaon CTime", 100, -10, 10)

        H_hsdelta_SUBPION_DATA  = ROOT.TH1D("H_hsdelta_SUBPION_DATA","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPION_DATA  = ROOT.TH1D("H_hsxptar_SUBPION_DATA","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPION_DATA  = ROOT.TH1D("H_hsyptar_SUBPION_DATA","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPION_DATA    = ROOT.TH1D("H_ssxfp_SUBPION_DATA","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPION_DATA    = ROOT.TH1D("H_ssyfp_SUBPION_DATA","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPION_DATA   = ROOT.TH1D("H_ssxpfp_SUBPION_DATA","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPION_DATA   = ROOT.TH1D("H_ssypfp_SUBPION_DATA","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPION_DATA    = ROOT.TH1D("H_hsxfp_SUBPION_DATA","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPION_DATA    = ROOT.TH1D("H_hsyfp_SUBPION_DATA","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPION_DATA   = ROOT.TH1D("H_hsxpfp_SUBPION_DATA","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPION_DATA   = ROOT.TH1D("H_hsypfp_SUBPION_DATA","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPION_DATA  = ROOT.TH1D("H_ssdelta_SUBPION_DATA","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPION_DATA  = ROOT.TH1D("H_ssxptar_SUBPION_DATA","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPION_DATA  = ROOT.TH1D("H_ssyptar_SUBPION_DATA","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPION_DATA        = ROOT.TH1D("H_q_SUBPION_DATA","q", 100, 0.0, 10.0)
        H_Q2_SUBPION_DATA       = ROOT.TH1D("H_Q2_SUBPION_DATA","Q2", 100, Q2min, Q2max)
        H_W_SUBPION_DATA  = ROOT.TH1D("H_W_SUBPION_DATA","W ", 100, Wmin, Wmax)
        H_t_SUBPION_DATA       = ROOT.TH1D("H_t_SUBPION_DATA","-t", 100, -1.0, 1.5)  
        H_epsilon_SUBPION_DATA  = ROOT.TH1D("H_epsilon_SUBPION_DATA","epsilon", 100, 0., 1.0)
        H_MM_SUBPION_DATA  = ROOT.TH1D("H_MM_SUBPION_DATA","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPION_DATA  = ROOT.TH1D("H_th_SUBPION_DATA","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPION_DATA  = ROOT.TH1D("H_ph_SUBPION_DATA","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPION_DATA  = ROOT.TH1D("H_ph_q_SUBPION_DATA","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPION_DATA  = ROOT.TH1D("H_th_q_SUBPION_DATA","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPION_DATA  = ROOT.TH1D("H_ph_recoil_SUBPION_DATA","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPION_DATA  = ROOT.TH1D("H_th_recoil_SUBPION_DATA","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPION_DATA  = ROOT.TH1D("H_pmiss_SUBPION_DATA","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPION_DATA  = ROOT.TH1D("H_emiss_SUBPION_DATA","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPION_DATA  = ROOT.TH1D("H_pmx_SUBPION_DATA","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPION_DATA  = ROOT.TH1D("H_pmy_SUBPION_DATA","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPION_DATA  = ROOT.TH1D("H_pmz_SUBPION_DATA","pmz", 100, -10.0, 10.0)
        H_ct_epi_SUBPION_DATA = ROOT.TH1D("H_ct_epi_SUBPION_DATA", "Electron-Pion CTime", 100, -10, 10)
        H_cal_etottracknorm_SUBPION_DATA = ROOT.TH1D("H_cal_etottracknorm_SUBPION_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        H_cer_npeSum_SUBPION_DATA = ROOT.TH1D("H_cer_npeSum_SUBPION_DATA", "HMS Cer Npe Sum", 100, 0, 30)
        P_cal_etottracknorm_SUBPION_DATA = ROOT.TH1D("P_cal_etottracknorm_SUBPION_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
        P_hgcer_npeSum_SUBPION_DATA = ROOT.TH1D("P_hgcer_npeSum_SUBPION_DATA", "SHMS HGCer Npe Sum", 100, 0, 50)
        P_aero_npeSum_SUBPION_DATA = ROOT.TH1D("P_aero_npeSum_SUBPION_DATA", "SHMS Aero Npe Sum", 100, 0, 50)

        H_hsdelta_SUBPION_DUMMY  = ROOT.TH1D("H_hsdelta_SUBPION_DUMMY","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPION_DUMMY  = ROOT.TH1D("H_hsxptar_SUBPION_DUMMY","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPION_DUMMY  = ROOT.TH1D("H_hsyptar_SUBPION_DUMMY","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPION_DUMMY    = ROOT.TH1D("H_ssxfp_SUBPION_DUMMY","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPION_DUMMY    = ROOT.TH1D("H_ssyfp_SUBPION_DUMMY","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPION_DUMMY   = ROOT.TH1D("H_ssxpfp_SUBPION_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPION_DUMMY   = ROOT.TH1D("H_ssypfp_SUBPION_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPION_DUMMY    = ROOT.TH1D("H_hsxfp_SUBPION_DUMMY","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPION_DUMMY    = ROOT.TH1D("H_hsyfp_SUBPION_DUMMY","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPION_DUMMY   = ROOT.TH1D("H_hsxpfp_SUBPION_DUMMY","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPION_DUMMY   = ROOT.TH1D("H_hsypfp_SUBPION_DUMMY","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPION_DUMMY  = ROOT.TH1D("H_ssdelta_SUBPION_DUMMY","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPION_DUMMY  = ROOT.TH1D("H_ssxptar_SUBPION_DUMMY","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPION_DUMMY  = ROOT.TH1D("H_ssyptar_SUBPION_DUMMY","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPION_DUMMY        = ROOT.TH1D("H_q_SUBPION_DUMMY","q", 100, 0.0, 10.0)
        H_Q2_SUBPION_DUMMY       = ROOT.TH1D("H_Q2_SUBPION_DUMMY","Q2", 100, Q2min, Q2max)
        H_W_SUBPION_DUMMY  = ROOT.TH1D("H_W_SUBPION_DUMMY","W ", 100, Wmin, Wmax)
        H_t_SUBPION_DUMMY       = ROOT.TH1D("H_t_SUBPION_DUMMY","-t", 100, -1.0, 1.5)  
        H_epsilon_SUBPION_DUMMY  = ROOT.TH1D("H_epsilon_SUBPION_DUMMY","epsilon", 100, 0., 1.0)
        H_MM_SUBPION_DUMMY  = ROOT.TH1D("H_MM_SUBPION_DUMMY","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPION_DUMMY  = ROOT.TH1D("H_th_SUBPION_DUMMY","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPION_DUMMY  = ROOT.TH1D("H_ph_SUBPION_DUMMY","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPION_DUMMY  = ROOT.TH1D("H_ph_q_SUBPION_DUMMY","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPION_DUMMY  = ROOT.TH1D("H_th_q_SUBPION_DUMMY","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPION_DUMMY  = ROOT.TH1D("H_ph_recoil_SUBPION_DUMMY","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPION_DUMMY  = ROOT.TH1D("H_th_recoil_SUBPION_DUMMY","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPION_DUMMY  = ROOT.TH1D("H_pmiss_SUBPION_DUMMY","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPION_DUMMY  = ROOT.TH1D("H_emiss_SUBPION_DUMMY","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPION_DUMMY  = ROOT.TH1D("H_pmx_SUBPION_DUMMY","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPION_DUMMY  = ROOT.TH1D("H_pmy_SUBPION_DUMMY","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPION_DUMMY  = ROOT.TH1D("H_pmz_SUBPION_DUMMY","pmz", 100, -10.0, 10.0)
        H_ct_epi_SUBPION_DUMMY = ROOT.TH1D("H_ct_epi_SUBPION_DUMMY", "Electron-Pion CTime", 100, -10, 10)

        H_hsdelta_SUBPION_RAND  = ROOT.TH1D("H_hsdelta_SUBPION_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPION_RAND  = ROOT.TH1D("H_hsxptar_SUBPION_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPION_RAND  = ROOT.TH1D("H_hsyptar_SUBPION_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPION_RAND    = ROOT.TH1D("H_ssxfp_SUBPION_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPION_RAND    = ROOT.TH1D("H_ssyfp_SUBPION_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPION_RAND   = ROOT.TH1D("H_ssxpfp_SUBPION_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPION_RAND   = ROOT.TH1D("H_ssypfp_SUBPION_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPION_RAND    = ROOT.TH1D("H_hsxfp_SUBPION_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPION_RAND    = ROOT.TH1D("H_hsyfp_SUBPION_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPION_RAND   = ROOT.TH1D("H_hsxpfp_SUBPION_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPION_RAND   = ROOT.TH1D("H_hsypfp_SUBPION_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPION_RAND  = ROOT.TH1D("H_ssdelta_SUBPION_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPION_RAND  = ROOT.TH1D("H_ssxptar_SUBPION_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPION_RAND  = ROOT.TH1D("H_ssyptar_SUBPION_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPION_RAND        = ROOT.TH1D("H_q_SUBPION_RAND","q", 100, 0.0, 10.0)
        H_Q2_SUBPION_RAND       = ROOT.TH1D("H_Q2_SUBPION_RAND","Q2", 100, Q2min, Q2max)
        H_W_SUBPION_RAND  = ROOT.TH1D("H_W_SUBPION_RAND","W ", 100, Wmin, Wmax)
        H_t_SUBPION_RAND       = ROOT.TH1D("H_t_SUBPION_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_SUBPION_RAND  = ROOT.TH1D("H_epsilon_SUBPION_RAND","epsilon", 100, 0., 1.0)
        H_MM_SUBPION_RAND  = ROOT.TH1D("H_MM_SUBPION_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPION_RAND  = ROOT.TH1D("H_th_SUBPION_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPION_RAND  = ROOT.TH1D("H_ph_SUBPION_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPION_RAND  = ROOT.TH1D("H_ph_q_SUBPION_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPION_RAND  = ROOT.TH1D("H_th_q_SUBPION_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPION_RAND  = ROOT.TH1D("H_ph_recoil_SUBPION_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPION_RAND  = ROOT.TH1D("H_th_recoil_SUBPION_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPION_RAND  = ROOT.TH1D("H_pmiss_SUBPION_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPION_RAND  = ROOT.TH1D("H_emiss_SUBPION_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPION_RAND  = ROOT.TH1D("H_pmx_SUBPION_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPION_RAND  = ROOT.TH1D("H_pmy_SUBPION_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPION_RAND  = ROOT.TH1D("H_pmz_SUBPION_RAND","pmz", 100, -10.0, 10.0)
        H_ct_epi_SUBPION_RAND = ROOT.TH1D("H_ct_epi_SUBPION_RAND", "Electron-Pion CTime", 100, -10, 10)

        H_hsdelta_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_hsdelta_SUBPION_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_hsxptar_SUBPION_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_hsyptar_SUBPION_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPION_DUMMY_RAND    = ROOT.TH1D("H_ssxfp_SUBPION_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPION_DUMMY_RAND    = ROOT.TH1D("H_ssyfp_SUBPION_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPION_DUMMY_RAND   = ROOT.TH1D("H_ssxpfp_SUBPION_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPION_DUMMY_RAND   = ROOT.TH1D("H_ssypfp_SUBPION_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPION_DUMMY_RAND    = ROOT.TH1D("H_hsxfp_SUBPION_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPION_DUMMY_RAND    = ROOT.TH1D("H_hsyfp_SUBPION_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPION_DUMMY_RAND   = ROOT.TH1D("H_hsxpfp_SUBPION_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPION_DUMMY_RAND   = ROOT.TH1D("H_hsypfp_SUBPION_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ssdelta_SUBPION_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ssxptar_SUBPION_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ssyptar_SUBPION_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPION_DUMMY_RAND        = ROOT.TH1D("H_q_SUBPION_DUMMY_RAND","q", 100, 0.0, 10.0)
        H_Q2_SUBPION_DUMMY_RAND       = ROOT.TH1D("H_Q2_SUBPION_DUMMY_RAND","Q2", 100, Q2min, Q2max)
        H_W_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_W_SUBPION_DUMMY_RAND","W ", 100, Wmin, Wmax)
        H_t_SUBPION_DUMMY_RAND       = ROOT.TH1D("H_t_SUBPION_DUMMY_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_epsilon_SUBPION_DUMMY_RAND","epsilon", 100, 0., 1.0)
        H_MM_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_MM_SUBPION_DUMMY_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_th_SUBPION_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ph_SUBPION_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ph_q_SUBPION_DUMMY_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_th_q_SUBPION_DUMMY_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_ph_recoil_SUBPION_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_th_recoil_SUBPION_DUMMY_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_pmiss_SUBPION_DUMMY_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_emiss_SUBPION_DUMMY_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_pmx_SUBPION_DUMMY_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_pmy_SUBPION_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPION_DUMMY_RAND  = ROOT.TH1D("H_pmz_SUBPION_DUMMY_RAND","pmz", 100, -10.0, 10.0)
        H_ct_epi_SUBPION_DUMMY_RAND = ROOT.TH1D("H_ct_epi_SUBPION_DUMMY_RAND", "Electron-Pion CTime", 100, -10, 10)

        H_hsdelta_SUBPROTON_DATA  = ROOT.TH1D("H_hsdelta_SUBPROTON_DATA","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPROTON_DATA  = ROOT.TH1D("H_hsxptar_SUBPROTON_DATA","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPROTON_DATA  = ROOT.TH1D("H_hsyptar_SUBPROTON_DATA","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPROTON_DATA    = ROOT.TH1D("H_ssxfp_SUBPROTON_DATA","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPROTON_DATA    = ROOT.TH1D("H_ssyfp_SUBPROTON_DATA","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPROTON_DATA   = ROOT.TH1D("H_ssxpfp_SUBPROTON_DATA","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPROTON_DATA   = ROOT.TH1D("H_ssypfp_SUBPROTON_DATA","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPROTON_DATA    = ROOT.TH1D("H_hsxfp_SUBPROTON_DATA","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPROTON_DATA    = ROOT.TH1D("H_hsyfp_SUBPROTON_DATA","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPROTON_DATA   = ROOT.TH1D("H_hsxpfp_SUBPROTON_DATA","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPROTON_DATA   = ROOT.TH1D("H_hsypfp_SUBPROTON_DATA","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPROTON_DATA  = ROOT.TH1D("H_ssdelta_SUBPROTON_DATA","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPROTON_DATA  = ROOT.TH1D("H_ssxptar_SUBPROTON_DATA","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPROTON_DATA  = ROOT.TH1D("H_ssyptar_SUBPROTON_DATA","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPROTON_DATA        = ROOT.TH1D("H_q_SUBPROTON_DATA","q", 100, 0.0, 10.0)
        H_Q2_SUBPROTON_DATA       = ROOT.TH1D("H_Q2_SUBPROTON_DATA","Q2", 100, Q2min, Q2max)
        H_W_SUBPROTON_DATA  = ROOT.TH1D("H_W_SUBPROTON_DATA","W ", 100, Wmin, Wmax)
        H_t_SUBPROTON_DATA       = ROOT.TH1D("H_t_SUBPROTON_DATA","-t", 100, -1.0, 1.5)  
        H_epsilon_SUBPROTON_DATA  = ROOT.TH1D("H_epsilon_SUBPROTON_DATA","epsilon", 100, 0., 1.0)
        H_MM_SUBPROTON_DATA  = ROOT.TH1D("H_MM_SUBPROTON_DATA","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPROTON_DATA  = ROOT.TH1D("H_th_SUBPROTON_DATA","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPROTON_DATA  = ROOT.TH1D("H_ph_SUBPROTON_DATA","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPROTON_DATA  = ROOT.TH1D("H_ph_q_SUBPROTON_DATA","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPROTON_DATA  = ROOT.TH1D("H_th_q_SUBPROTON_DATA","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPROTON_DATA  = ROOT.TH1D("H_ph_recoil_SUBPROTON_DATA","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPROTON_DATA  = ROOT.TH1D("H_th_recoil_SUBPROTON_DATA","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPROTON_DATA  = ROOT.TH1D("H_pmiss_SUBPROTON_DATA","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPROTON_DATA  = ROOT.TH1D("H_emiss_SUBPROTON_DATA","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPROTON_DATA  = ROOT.TH1D("H_pmx_SUBPROTON_DATA","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPROTON_DATA  = ROOT.TH1D("H_pmy_SUBPROTON_DATA","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPROTON_DATA  = ROOT.TH1D("H_pmz_SUBPROTON_DATA","pmz", 100, -10.0, 10.0)
        H_ct_ep_SUBPROTON_DATA = ROOT.TH1D("H_ct_ep_SUBPROTON_DATA", "Electron-Proton CTime", 100, -10, 10)
        H_cal_etottracknorm_SUBPROTON_DATA = ROOT.TH1D("H_cal_etottracknorm_SUBPROTON_DATA", "HMS Cal etottracknorm", 100, 0.2, 1.8)
        H_cer_npeSum_SUBPROTON_DATA = ROOT.TH1D("H_cer_npeSum_SUBPROTON_DATA", "HMS Cer Npe Sum", 100, 0, 30)
        P_cal_etottracknorm_SUBPROTON_DATA = ROOT.TH1D("P_cal_etottracknorm_SUBPROTON_DATA", "SHMS Cal etottracknorm", 100, 0, 1)
        P_hgcer_npeSum_SUBPROTON_DATA = ROOT.TH1D("P_hgcer_npeSum_SUBPROTON_DATA", "SHMS HGCer Npe Sum", 100, 0, 50)
        P_aero_npeSum_SUBPROTON_DATA = ROOT.TH1D("P_aero_npeSum_SUBPROTON_DATA", "SHMS Aero Npe Sum", 100, 0, 50)

        H_hsdelta_SUBPROTON_DUMMY  = ROOT.TH1D("H_hsdelta_SUBPROTON_DUMMY","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPROTON_DUMMY  = ROOT.TH1D("H_hsxptar_SUBPROTON_DUMMY","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPROTON_DUMMY  = ROOT.TH1D("H_hsyptar_SUBPROTON_DUMMY","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPROTON_DUMMY    = ROOT.TH1D("H_ssxfp_SUBPROTON_DUMMY","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPROTON_DUMMY    = ROOT.TH1D("H_ssyfp_SUBPROTON_DUMMY","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPROTON_DUMMY   = ROOT.TH1D("H_ssxpfp_SUBPROTON_DUMMY","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPROTON_DUMMY   = ROOT.TH1D("H_ssypfp_SUBPROTON_DUMMY","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPROTON_DUMMY    = ROOT.TH1D("H_hsxfp_SUBPROTON_DUMMY","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPROTON_DUMMY    = ROOT.TH1D("H_hsyfp_SUBPROTON_DUMMY","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPROTON_DUMMY   = ROOT.TH1D("H_hsxpfp_SUBPROTON_DUMMY","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPROTON_DUMMY   = ROOT.TH1D("H_hsypfp_SUBPROTON_DUMMY","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPROTON_DUMMY  = ROOT.TH1D("H_ssdelta_SUBPROTON_DUMMY","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPROTON_DUMMY  = ROOT.TH1D("H_ssxptar_SUBPROTON_DUMMY","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPROTON_DUMMY  = ROOT.TH1D("H_ssyptar_SUBPROTON_DUMMY","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPROTON_DUMMY        = ROOT.TH1D("H_q_SUBPROTON_DUMMY","q", 100, 0.0, 10.0)
        H_Q2_SUBPROTON_DUMMY       = ROOT.TH1D("H_Q2_SUBPROTON_DUMMY","Q2", 100, Q2min, Q2max)
        H_W_SUBPROTON_DUMMY  = ROOT.TH1D("H_W_SUBPROTON_DUMMY","W ", 100, Wmin, Wmax)
        H_t_SUBPROTON_DUMMY       = ROOT.TH1D("H_t_SUBPROTON_DUMMY","-t", 100, -1.0, 1.5)  
        H_epsilon_SUBPROTON_DUMMY  = ROOT.TH1D("H_epsilon_SUBPROTON_DUMMY","epsilon", 100, 0., 1.0)
        H_MM_SUBPROTON_DUMMY  = ROOT.TH1D("H_MM_SUBPROTON_DUMMY","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPROTON_DUMMY  = ROOT.TH1D("H_th_SUBPROTON_DUMMY","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPROTON_DUMMY  = ROOT.TH1D("H_ph_SUBPROTON_DUMMY","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPROTON_DUMMY  = ROOT.TH1D("H_ph_q_SUBPROTON_DUMMY","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPROTON_DUMMY  = ROOT.TH1D("H_th_q_SUBPROTON_DUMMY","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPROTON_DUMMY  = ROOT.TH1D("H_ph_recoil_SUBPROTON_DUMMY","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPROTON_DUMMY  = ROOT.TH1D("H_th_recoil_SUBPROTON_DUMMY","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPROTON_DUMMY  = ROOT.TH1D("H_pmiss_SUBPROTON_DUMMY","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPROTON_DUMMY  = ROOT.TH1D("H_emiss_SUBPROTON_DUMMY","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPROTON_DUMMY  = ROOT.TH1D("H_pmx_SUBPROTON_DUMMY","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPROTON_DUMMY  = ROOT.TH1D("H_pmy_SUBPROTON_DUMMY","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPROTON_DUMMY  = ROOT.TH1D("H_pmz_SUBPROTON_DUMMY","pmz", 100, -10.0, 10.0)
        H_ct_ep_SUBPROTON_DUMMY = ROOT.TH1D("H_ct_ep_SUBPROTON_DUMMY", "Electron-Proton CTime", 100, -10, 10)

        H_hsdelta_SUBPROTON_RAND  = ROOT.TH1D("H_hsdelta_SUBPROTON_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPROTON_RAND  = ROOT.TH1D("H_hsxptar_SUBPROTON_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPROTON_RAND  = ROOT.TH1D("H_hsyptar_SUBPROTON_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPROTON_RAND    = ROOT.TH1D("H_ssxfp_SUBPROTON_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPROTON_RAND    = ROOT.TH1D("H_ssyfp_SUBPROTON_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPROTON_RAND   = ROOT.TH1D("H_ssxpfp_SUBPROTON_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPROTON_RAND   = ROOT.TH1D("H_ssypfp_SUBPROTON_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPROTON_RAND    = ROOT.TH1D("H_hsxfp_SUBPROTON_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPROTON_RAND    = ROOT.TH1D("H_hsyfp_SUBPROTON_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPROTON_RAND   = ROOT.TH1D("H_hsxpfp_SUBPROTON_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPROTON_RAND   = ROOT.TH1D("H_hsypfp_SUBPROTON_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPROTON_RAND  = ROOT.TH1D("H_ssdelta_SUBPROTON_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPROTON_RAND  = ROOT.TH1D("H_ssxptar_SUBPROTON_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPROTON_RAND  = ROOT.TH1D("H_ssyptar_SUBPROTON_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPROTON_RAND        = ROOT.TH1D("H_q_SUBPROTON_RAND","q", 100, 0.0, 10.0)
        H_Q2_SUBPROTON_RAND       = ROOT.TH1D("H_Q2_SUBPROTON_RAND","Q2", 100, Q2min, Q2max)
        H_W_SUBPROTON_RAND  = ROOT.TH1D("H_W_SUBPROTON_RAND","W ", 100, Wmin, Wmax)
        H_t_SUBPROTON_RAND       = ROOT.TH1D("H_t_SUBPROTON_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_SUBPROTON_RAND  = ROOT.TH1D("H_epsilon_SUBPROTON_RAND","epsilon", 100, 0., 1.0)
        H_MM_SUBPROTON_RAND  = ROOT.TH1D("H_MM_SUBPROTON_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPROTON_RAND  = ROOT.TH1D("H_th_SUBPROTON_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPROTON_RAND  = ROOT.TH1D("H_ph_SUBPROTON_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPROTON_RAND  = ROOT.TH1D("H_ph_q_SUBPROTON_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPROTON_RAND  = ROOT.TH1D("H_th_q_SUBPROTON_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPROTON_RAND  = ROOT.TH1D("H_ph_recoil_SUBPROTON_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPROTON_RAND  = ROOT.TH1D("H_th_recoil_SUBPROTON_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPROTON_RAND  = ROOT.TH1D("H_pmiss_SUBPROTON_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPROTON_RAND  = ROOT.TH1D("H_emiss_SUBPROTON_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPROTON_RAND  = ROOT.TH1D("H_pmx_SUBPROTON_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPROTON_RAND  = ROOT.TH1D("H_pmy_SUBPROTON_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPROTON_RAND  = ROOT.TH1D("H_pmz_SUBPROTON_RAND","pmz", 100, -10.0, 10.0)
        H_ct_ep_SUBPROTON_RAND = ROOT.TH1D("H_ct_ep_SUBPROTON_RAND", "Electron-Proton CTime", 100, -10, 10)

        H_hsdelta_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_hsdelta_SUBPROTON_DUMMY_RAND","HMS Delta", 100, -20.0, 20.0)
        H_hsxptar_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_hsxptar_SUBPROTON_DUMMY_RAND","HMS xptar", 100, -0.1, 0.1)
        H_hsyptar_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_hsyptar_SUBPROTON_DUMMY_RAND","HMS yptar", 100, -0.1, 0.1)
        H_ssxfp_SUBPROTON_DUMMY_RAND    = ROOT.TH1D("H_ssxfp_SUBPROTON_DUMMY_RAND","SHMS xfp", 100, -25.0, 25.0)
        H_ssyfp_SUBPROTON_DUMMY_RAND    = ROOT.TH1D("H_ssyfp_SUBPROTON_DUMMY_RAND","SHMS yfp", 100, -25.0, 25.0)
        H_ssxpfp_SUBPROTON_DUMMY_RAND   = ROOT.TH1D("H_ssxpfp_SUBPROTON_DUMMY_RAND","SHMS xpfp", 100, -0.09, 0.09)
        H_ssypfp_SUBPROTON_DUMMY_RAND   = ROOT.TH1D("H_ssypfp_SUBPROTON_DUMMY_RAND","SHMS ypfp", 100, -0.05, 0.04)
        H_hsxfp_SUBPROTON_DUMMY_RAND    = ROOT.TH1D("H_hsxfp_SUBPROTON_DUMMY_RAND","HMS xfp", 100, -40.0, 40.0)
        H_hsyfp_SUBPROTON_DUMMY_RAND    = ROOT.TH1D("H_hsyfp_SUBPROTON_DUMMY_RAND","HMS yfp", 100, -20.0, 20.0)
        H_hsxpfp_SUBPROTON_DUMMY_RAND   = ROOT.TH1D("H_hsxpfp_SUBPROTON_DUMMY_RAND","HMS xpfp", 100, -0.09, 0.05)
        H_hsypfp_SUBPROTON_DUMMY_RAND   = ROOT.TH1D("H_hsypfp_SUBPROTON_DUMMY_RAND","HMS ypfp", 100, -0.05, 0.04)
        H_ssdelta_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ssdelta_SUBPROTON_DUMMY_RAND","SHMS delta", 100, -20.0, 20.0)
        H_ssxptar_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ssxptar_SUBPROTON_DUMMY_RAND","SHMS xptar", 100, -0.1, 0.1)
        H_ssyptar_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ssyptar_SUBPROTON_DUMMY_RAND","SHMS yptar", 100, -0.04, 0.04)
        H_q_SUBPROTON_DUMMY_RAND        = ROOT.TH1D("H_q_SUBPROTON_DUMMY_RAND","q", 100, 0.0, 10.0)
        H_Q2_SUBPROTON_DUMMY_RAND       = ROOT.TH1D("H_Q2_SUBPROTON_DUMMY_RAND","Q2", 100, Q2min, Q2max)
        H_W_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_W_SUBPROTON_DUMMY_RAND","W ", 100, Wmin, Wmax)
        H_t_SUBPROTON_DUMMY_RAND       = ROOT.TH1D("H_t_SUBPROTON_DUMMY_RAND","-t", 100, -1.0, 1.5)
        H_epsilon_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_epsilon_SUBPROTON_DUMMY_RAND","epsilon", 100, 0., 1.0)
        H_MM_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_MM_SUBPROTON_DUMMY_RAND","MM_{K}", 100, 0.0, 1.5)
        H_th_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_th_SUBPROTON_DUMMY_RAND","X' tar", 100, -0.1, 0.1)
        H_ph_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ph_SUBPROTON_DUMMY_RAND","Y' tar", 100, -0.1, 0.1)
        H_ph_q_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ph_q_SUBPROTON_DUMMY_RAND","Phi Detected (ph_xq)", 100, -5.0, 5.0)
        H_th_q_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_th_q_SUBPROTON_DUMMY_RAND","Theta Detected (th_xq)", 100, -0.2, 0.2)
        H_ph_recoil_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_ph_recoil_SUBPROTON_DUMMY_RAND","Phi Recoil (ph_bq)", 100, -10.0, 10.0)
        H_th_recoil_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_th_recoil_SUBPROTON_DUMMY_RAND","Theta Recoil (th_bq)", 100, -10.0, 10.0)
        H_pmiss_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_pmiss_SUBPROTON_DUMMY_RAND","pmiss", 100, 0.0, 10.0)
        H_emiss_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_emiss_SUBPROTON_DUMMY_RAND","emiss", 100, 0.0, 10.0)
        H_pmx_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_pmx_SUBPROTON_DUMMY_RAND","pmx", 100, -10.0, 10.0)
        H_pmy_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_pmy_SUBPROTON_DUMMY_RAND","pmy ", 100, -10.0, 10.0)
        H_pmz_SUBPROTON_DUMMY_RAND  = ROOT.TH1D("H_pmz_SUBPROTON_DUMMY_RAND","pmz", 100, -10.0, 10.0)
        H_ct_ep_SUBPROTON_DUMMY_RAND = ROOT.TH1D("H_ct_ep_SUBPROTON_DUMMY_RAND", "Electron-Proton CTime", 100, -10, 10)
        
        ################################################################################################################################################
        # t/phi binned histograms

        H_phibins_DATA = ROOT.TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
        H_tbins_DATA = ROOT.TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)
        H_yield_DATA = ROOT.TH1D("H_yield_DATA", "Data Yield", NumtBins*NumPhiBins, 0, 1.0)

        tbinDict = {}
        for i in range(NumtBins):
            tbinDict["H_Q2_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_Q2_tbin_DATA_{}".format(i+1), "Q2 (t bin {})".format(i+1), 100, Q2min, Q2max)
            tbinDict["H_W_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_W_tbin_DATA_{}".format(i+1), "W (t bin {})".format(i+1), 100, Wmin, Wmax)
            tbinDict["H_t_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_t_tbin_DATA_{}".format(i+1), "t (t bin {})".format(i+1), 100, tmin, tmax)

        ################################################################################################################################################
        # 2D histograms

        MM_vs_CoinTime_DATA = ROOT.TH2D("MM_vs_CoinTime_DATA","Missing Mass vs CTime; MM; Coin_Time",100, 0, 2, 100, -2, 2)
        CoinTime_vs_beta_DATA = ROOT.TH2D("CoinTime_vs_beta_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0, 2)
        MM_vs_beta_DATA = ROOT.TH2D("MM_vs_beta_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, 0, 2, 100, 0, 2)
        phiq_vs_t_DATA = ROOT.TH2D("phiq_vs_t_DATA","; #phi ;t", 12, -3.14, 3.14, 24, tmin, tmax)
        polar_phiq_vs_t_DATA = ROOT.TGraphPolar()
        Q2_vs_W_DATA = ROOT.TH2D("Q2_vs_W_DATA", "Q^{2} vs W; Q^{2}; W", 100, Q2min, Q2max, 100, Wmin, Wmax)

        ################################################################################################################################################
        # Fill data histograms for various trees called above

        print("\nGrabbing %s simc..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SIMC):

          # Progress bar
          Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

          # Define the acceptance cuts  
          SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
          HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
          if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
              Diamond = True
          else:
              try:
                  Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
              except ZeroDivisionError:
                  Diamond = False

          #........................................

          #Fill SIMC events
          if(HMS_Acceptance & SHMS_Acceptance & Diamond):

              H_ssxfp_SIMC.Fill(evt.ssxfp, evt.Weight)
              H_ssyfp_SIMC.Fill(evt.ssyfp, evt.Weight)
              H_ssxpfp_SIMC.Fill(evt.ssxpfp, evt.Weight)
              H_ssypfp_SIMC.Fill(evt.ssypfp, evt.Weight)
              H_hsxfp_SIMC.Fill(evt.hsxfp, evt.Weight)
              H_hsyfp_SIMC.Fill(evt.hsyfp, evt.Weight)
              H_hsxpfp_SIMC.Fill(evt.hsxpfp, evt.Weight)
              H_hsypfp_SIMC.Fill(evt.hsypfp, evt.Weight)
              H_ssdelta_SIMC.Fill(evt.ssdelta, evt.Weight) 
              H_hsdelta_SIMC.Fill(evt.hsdelta, evt.Weight)	
              H_ssxptar_SIMC.Fill(evt.ssxptar, evt.Weight)
              H_ssyptar_SIMC.Fill(evt.ssyptar, evt.Weight)
              H_hsxptar_SIMC.Fill(evt.hsxptar, evt.Weight)	
              H_hsyptar_SIMC.Fill(evt.hsyptar, evt.Weight)

              H_ph_q_SIMC.Fill(evt.phipq, evt.Weight)
              H_th_q_SIMC.Fill(evt.thetapq, evt.Weight)

              H_pmiss_SIMC.Fill(evt.Pm, evt.Weight)	
              H_emiss_SIMC.Fill(evt.Em, evt.Weight)	
              #H_pmx_SIMC.Fill(evt.Pmx, evt.Weight)
              #H_pmy_SIMC.Fill(evt.Pmy, evt.Weight)
              #H_pmz_SIMC.Fill(evt.Pmz, evt.Weight)
              H_Q2_SIMC.Fill(evt.Q2, evt.Weight)
              H_W_SIMC.Fill(evt.W, evt.Weight)
              H_t_DATA.Fill(evt.t)
              H_epsilon_SIMC.Fill(evt.epsilon, evt.Weight)
              H_MM_SIMC.Fill(np.sqrt(pow(evt.Em, 2) - pow(evt.Pm, 2)), evt.Weight)


        ################################################################################################################################################
        # Fill histograms for various trees called above

        print("\nGrabbing %s kaon data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DATA):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            # Must be outside diamond cuts to avoid weird overflow errors
            polar_phiq_vs_t_DATA.SetPoint(polar_phiq_vs_t_DATA.GetN(), (evt.ph_q+math.pi)*(180/math.pi), abs(evt.MandelT))

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              MM_vs_CoinTime_DATA.Fill(evt.MM, evt.CTime_ROC1)
              CoinTime_vs_beta_DATA.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
              MM_vs_beta_DATA.Fill(evt.MM,evt.P_gtr_beta)
              phiq_vs_t_DATA.Fill(evt.ph_q, -evt.MandelT)
              #polar_phiq_vs_t_DATA.SetPoint(i, evt.ph_q, -evt.MandelT)          
              Q2_vs_W_DATA.Fill(evt.Q2, evt.W)

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

              H_ph_q_DATA.Fill(evt.ph_q)
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
              H_MM_DATA.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_DATA.Fill(pow(evt.MM, 2))  
              #H_MM_DATA.Fill(evt.Mrecoil)

              H_cal_etottracknorm_DATA.Fill(evt.H_cal_etottracknorm)
              H_cer_npeSum_DATA.Fill(evt.H_cer_npeSum)

              P_cal_etottracknorm_DATA.Fill(evt.P_cal_etottracknorm)
              P_hgcer_npeSum_DATA.Fill(evt.P_hgcer_npeSum)
              P_aero_npeSum_DATA.Fill(evt.P_aero_npeSum)          

        ################################################################################################################################################
        # Fill dummy histograms for various trees called above

        print("\nGrabbing %s kaon dummy..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DUMMY):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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

              H_ph_q_DUMMY.Fill(evt.ph_q)
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
              H_MM_DUMMY.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_DUMMY.Fill(pow(evt.MM, 2))  
              #H_MM_DUMMY.Fill(evt.Mrecoil)

        ###################################################################################################################################################    
        # Fill random histograms for various trees called above

        print("\nGrabbing %s kaon random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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
              H_MM_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ###################################################################################################################################################    
        # Fill dummy random histograms for various trees called above

        print("\nGrabbing %s kaon dummy random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DUMMY_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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
              H_MM_DUMMY_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ################################################################################################################################################
        # Fill histograms for various trees called above

        print("\nGrabbing %s pion data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPION_DATA):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPION_DATA.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ct_epi_SUBPION_DATA.Fill(evt.CTime_ROC1)

              H_ssxfp_SUBPION_DATA.Fill(evt.ssxfp)
              H_ssyfp_SUBPION_DATA.Fill(evt.ssyfp)
              H_ssxpfp_SUBPION_DATA.Fill(evt.ssxpfp)
              H_ssypfp_SUBPION_DATA.Fill(evt.ssypfp)
              H_ssdelta_SUBPION_DATA.Fill(evt.ssdelta)
              H_ssxptar_SUBPION_DATA.Fill(evt.ssxptar)
              H_ssyptar_SUBPION_DATA.Fill(evt.ssyptar)

              H_hsxfp_SUBPION_DATA.Fill(evt.hsxfp)
              H_hsyfp_SUBPION_DATA.Fill(evt.hsyfp)
              H_hsxpfp_SUBPION_DATA.Fill(evt.hsxpfp)
              H_hsypfp_SUBPION_DATA.Fill(evt.hsypfp)
              H_hsdelta_SUBPION_DATA.Fill(evt.hsdelta)
              H_hsxptar_SUBPION_DATA.Fill(evt.hsxptar)	
              H_hsyptar_SUBPION_DATA.Fill(evt.hsyptar)

              H_ph_q_SUBPION_DATA.Fill(evt.ph_q)
              H_th_q_SUBPION_DATA.Fill(evt.th_q)
              H_ph_recoil_SUBPION_DATA.Fill(evt.ph_recoil)
              H_th_recoil_SUBPION_DATA.Fill(evt.th_recoil)

              H_pmiss_SUBPION_DATA.Fill(evt.pmiss)	
              H_emiss_SUBPION_DATA.Fill(evt.emiss)	
              #H_emiss_SUBPION_DATA.Fill(evt.emiss_nuc)
              H_pmx_SUBPION_DATA.Fill(evt.pmx)
              H_pmy_SUBPION_DATA.Fill(evt.pmy)
              H_pmz_SUBPION_DATA.Fill(evt.pmz)
              H_Q2_SUBPION_DATA.Fill(evt.Q2)
              H_t_SUBPION_DATA.Fill(-evt.MandelT)
              H_W_SUBPION_DATA.Fill(evt.W)
              H_epsilon_SUBPION_DATA.Fill(evt.epsilon)
              H_MM_SUBPION_DATA.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_SUBPION_DATA.Fill(pow(evt.MM, 2))  
              #H_MM_SUBPION_DATA.Fill(evt.Mrecoil)

              H_cal_etottracknorm_SUBPION_DATA.Fill(evt.H_cal_etottracknorm)
              H_cer_npeSum_SUBPION_DATA.Fill(evt.H_cer_npeSum)

              P_cal_etottracknorm_SUBPION_DATA.Fill(evt.P_cal_etottracknorm)
              P_hgcer_npeSum_SUBPION_DATA.Fill(evt.P_hgcer_npeSum)
              P_aero_npeSum_SUBPION_DATA.Fill(evt.P_aero_npeSum)          

        ################################################################################################################################################
        # Fill dummy histograms for various trees called above

        print("\nGrabbing %s pion dummy..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPION_DUMMY):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPION_DUMMY.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ct_epi_SUBPION_DUMMY.Fill(evt.CTime_ROC1)

              H_ssxfp_SUBPION_DUMMY.Fill(evt.ssxfp)
              H_ssyfp_SUBPION_DUMMY.Fill(evt.ssyfp)
              H_ssxpfp_SUBPION_DUMMY.Fill(evt.ssxpfp)
              H_ssypfp_SUBPION_DUMMY.Fill(evt.ssypfp)
              H_ssdelta_SUBPION_DUMMY.Fill(evt.ssdelta)
              H_ssxptar_SUBPION_DUMMY.Fill(evt.ssxptar)
              H_ssyptar_SUBPION_DUMMY.Fill(evt.ssyptar)

              H_hsxfp_SUBPION_DUMMY.Fill(evt.hsxfp)
              H_hsyfp_SUBPION_DUMMY.Fill(evt.hsyfp)
              H_hsxpfp_SUBPION_DUMMY.Fill(evt.hsxpfp)
              H_hsypfp_SUBPION_DUMMY.Fill(evt.hsypfp)
              H_hsdelta_SUBPION_DUMMY.Fill(evt.hsdelta)
              H_hsxptar_SUBPION_DUMMY.Fill(evt.hsxptar)	
              H_hsyptar_SUBPION_DUMMY.Fill(evt.hsyptar)

              H_ph_q_SUBPION_DUMMY.Fill(evt.ph_q)
              H_th_q_SUBPION_DUMMY.Fill(evt.th_q)
              H_ph_recoil_SUBPION_DUMMY.Fill(evt.ph_recoil)
              H_th_recoil_SUBPION_DUMMY.Fill(evt.th_recoil)

              H_pmiss_SUBPION_DUMMY.Fill(evt.pmiss)	
              H_emiss_SUBPION_DUMMY.Fill(evt.emiss)	
              #H_emiss_SUBPION_DUMMY.Fill(evt.emiss_nuc)
              H_pmx_SUBPION_DUMMY.Fill(evt.pmx)
              H_pmy_SUBPION_DUMMY.Fill(evt.pmy)
              H_pmz_SUBPION_DUMMY.Fill(evt.pmz)
              H_Q2_SUBPION_DUMMY.Fill(evt.Q2)
              H_t_SUBPION_DUMMY.Fill(-evt.MandelT)
              H_W_SUBPION_DUMMY.Fill(evt.W)
              H_epsilon_SUBPION_DUMMY.Fill(evt.epsilon)
              H_MM_SUBPION_DUMMY.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_SUBPION_DUMMY.Fill(pow(evt.MM, 2))  
              #H_MM_SUBPION_DUMMY.Fill(evt.Mrecoil)

        ###################################################################################################################################################    
        # Fill random histograms for various trees called above

        print("\nGrabbing %s pion random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPION_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPION_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ssxfp_SUBPION_RAND.Fill(evt.ssxfp)
              H_ssyfp_SUBPION_RAND.Fill(evt.ssyfp)
              H_ssxpfp_SUBPION_RAND.Fill(evt.ssxpfp)
              H_ssypfp_SUBPION_RAND.Fill(evt.ssypfp)
              H_ssdelta_SUBPION_RAND.Fill(evt.ssdelta)
              H_ssxptar_SUBPION_RAND.Fill(evt.ssxptar)
              H_ssyptar_SUBPION_RAND.Fill(evt.ssyptar)

              H_hsxfp_SUBPION_RAND.Fill(evt.hsxfp)
              H_hsyfp_SUBPION_RAND.Fill(evt.hsyfp)
              H_hsxpfp_SUBPION_RAND.Fill(evt.hsxpfp)
              H_hsypfp_SUBPION_RAND.Fill(evt.hsypfp)
              H_hsdelta_SUBPION_RAND.Fill(evt.hsdelta)
              H_hsxptar_SUBPION_RAND.Fill(evt.hsxptar)	
              H_hsyptar_SUBPION_RAND.Fill(evt.hsyptar)

              H_pmiss_SUBPION_RAND.Fill(evt.pmiss)	
              H_emiss_SUBPION_RAND.Fill(evt.emiss)	
              #H_emiss_SUBPION_RAND.Fill(evt.emiss_nuc)
              H_pmx_SUBPION_RAND.Fill(evt.pmx)
              H_pmy_SUBPION_RAND.Fill(evt.pmy)
              H_pmz_SUBPION_RAND.Fill(evt.pmz)
              H_Q2_SUBPION_RAND.Fill(evt.Q2)
              H_t_SUBPION_RAND.Fill(-evt.MandelT)
              H_W_SUBPION_RAND.Fill(evt.W)
              H_epsilon_SUBPION_RAND.Fill(evt.epsilon)
              H_MM_SUBPION_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ###################################################################################################################################################    
        # Fill dummy random histograms for various trees called above

        print("\nGrabbing %s pion dummy random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPION_DUMMY_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPION_DUMMY_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ssxfp_SUBPION_DUMMY_RAND.Fill(evt.ssxfp)
              H_ssyfp_SUBPION_DUMMY_RAND.Fill(evt.ssyfp)
              H_ssxpfp_SUBPION_DUMMY_RAND.Fill(evt.ssxpfp)
              H_ssypfp_SUBPION_DUMMY_RAND.Fill(evt.ssypfp)
              H_ssdelta_SUBPION_DUMMY_RAND.Fill(evt.ssdelta)
              H_ssxptar_SUBPION_DUMMY_RAND.Fill(evt.ssxptar)
              H_ssyptar_SUBPION_DUMMY_RAND.Fill(evt.ssyptar)

              H_hsxfp_SUBPION_DUMMY_RAND.Fill(evt.hsxfp)
              H_hsyfp_SUBPION_DUMMY_RAND.Fill(evt.hsyfp)
              H_hsxpfp_SUBPION_DUMMY_RAND.Fill(evt.hsxpfp)
              H_hsypfp_SUBPION_DUMMY_RAND.Fill(evt.hsypfp)
              H_hsdelta_SUBPION_DUMMY_RAND.Fill(evt.hsdelta)
              H_hsxptar_SUBPION_DUMMY_RAND.Fill(evt.hsxptar)	
              H_hsyptar_SUBPION_DUMMY_RAND.Fill(evt.hsyptar)

              H_pmiss_SUBPION_DUMMY_RAND.Fill(evt.pmiss)	
              H_emiss_SUBPION_DUMMY_RAND.Fill(evt.emiss)	
              #H_emiss_SUBPION_DUMMY_RAND.Fill(evt.emiss_nuc)
              H_pmx_SUBPION_DUMMY_RAND.Fill(evt.pmx)
              H_pmy_SUBPION_DUMMY_RAND.Fill(evt.pmy)
              H_pmz_SUBPION_DUMMY_RAND.Fill(evt.pmz)
              H_Q2_SUBPION_DUMMY_RAND.Fill(evt.Q2)
              H_t_SUBPION_DUMMY_RAND.Fill(-evt.MandelT)
              H_W_SUBPION_DUMMY_RAND.Fill(evt.W)
              H_epsilon_SUBPION_DUMMY_RAND.Fill(evt.epsilon)
              H_MM_SUBPION_DUMMY_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ################################################################################################################################################
        # Fill histograms for various trees called above

        print("\nGrabbing %s proton data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPROTON_DATA):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPROTON_DATA.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ct_ep_SUBPROTON_DATA.Fill(evt.CTime_ROC1)

              H_ssxfp_SUBPROTON_DATA.Fill(evt.ssxfp)
              H_ssyfp_SUBPROTON_DATA.Fill(evt.ssyfp)
              H_ssxpfp_SUBPROTON_DATA.Fill(evt.ssxpfp)
              H_ssypfp_SUBPROTON_DATA.Fill(evt.ssypfp)
              H_ssdelta_SUBPROTON_DATA.Fill(evt.ssdelta)
              H_ssxptar_SUBPROTON_DATA.Fill(evt.ssxptar)
              H_ssyptar_SUBPROTON_DATA.Fill(evt.ssyptar)

              H_hsxfp_SUBPROTON_DATA.Fill(evt.hsxfp)
              H_hsyfp_SUBPROTON_DATA.Fill(evt.hsyfp)
              H_hsxpfp_SUBPROTON_DATA.Fill(evt.hsxpfp)
              H_hsypfp_SUBPROTON_DATA.Fill(evt.hsypfp)
              H_hsdelta_SUBPROTON_DATA.Fill(evt.hsdelta)
              H_hsxptar_SUBPROTON_DATA.Fill(evt.hsxptar)	
              H_hsyptar_SUBPROTON_DATA.Fill(evt.hsyptar)

              H_ph_q_SUBPROTON_DATA.Fill(evt.ph_q)
              H_th_q_SUBPROTON_DATA.Fill(evt.th_q)
              H_ph_recoil_SUBPROTON_DATA.Fill(evt.ph_recoil)
              H_th_recoil_SUBPROTON_DATA.Fill(evt.th_recoil)

              H_pmiss_SUBPROTON_DATA.Fill(evt.pmiss)	
              H_emiss_SUBPROTON_DATA.Fill(evt.emiss)	
              #H_emiss_SUBPROTON_DATA.Fill(evt.emiss_nuc)
              H_pmx_SUBPROTON_DATA.Fill(evt.pmx)
              H_pmy_SUBPROTON_DATA.Fill(evt.pmy)
              H_pmz_SUBPROTON_DATA.Fill(evt.pmz)
              H_Q2_SUBPROTON_DATA.Fill(evt.Q2)
              H_t_SUBPROTON_DATA.Fill(-evt.MandelT)
              H_W_SUBPROTON_DATA.Fill(evt.W)
              H_epsilon_SUBPROTON_DATA.Fill(evt.epsilon)
              H_MM_SUBPROTON_DATA.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_SUBPROTON_DATA.Fill(pow(evt.MM, 2))  
              #H_MM_SUBPROTON_DATA.Fill(evt.Mrecoil)

              H_cal_etottracknorm_SUBPROTON_DATA.Fill(evt.H_cal_etottracknorm)
              H_cer_npeSum_SUBPROTON_DATA.Fill(evt.H_cer_npeSum)

              P_cal_etottracknorm_SUBPROTON_DATA.Fill(evt.P_cal_etottracknorm)
              P_hgcer_npeSum_SUBPROTON_DATA.Fill(evt.P_hgcer_npeSum)
              P_aero_npeSum_SUBPROTON_DATA.Fill(evt.P_aero_npeSum)          

        ################################################################################################################################################
        # Fill dummy histograms for various trees called above

        print("\nGrabbing %s proton dummy..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPROTON_DUMMY):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPROTON_DUMMY.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ct_ep_SUBPROTON_DUMMY.Fill(evt.CTime_ROC1)

              H_ssxfp_SUBPROTON_DUMMY.Fill(evt.ssxfp)
              H_ssyfp_SUBPROTON_DUMMY.Fill(evt.ssyfp)
              H_ssxpfp_SUBPROTON_DUMMY.Fill(evt.ssxpfp)
              H_ssypfp_SUBPROTON_DUMMY.Fill(evt.ssypfp)
              H_ssdelta_SUBPROTON_DUMMY.Fill(evt.ssdelta)
              H_ssxptar_SUBPROTON_DUMMY.Fill(evt.ssxptar)
              H_ssyptar_SUBPROTON_DUMMY.Fill(evt.ssyptar)

              H_hsxfp_SUBPROTON_DUMMY.Fill(evt.hsxfp)
              H_hsyfp_SUBPROTON_DUMMY.Fill(evt.hsyfp)
              H_hsxpfp_SUBPROTON_DUMMY.Fill(evt.hsxpfp)
              H_hsypfp_SUBPROTON_DUMMY.Fill(evt.hsypfp)
              H_hsdelta_SUBPROTON_DUMMY.Fill(evt.hsdelta)
              H_hsxptar_SUBPROTON_DUMMY.Fill(evt.hsxptar)	
              H_hsyptar_SUBPROTON_DUMMY.Fill(evt.hsyptar)

              H_ph_q_SUBPROTON_DUMMY.Fill(evt.ph_q)
              H_th_q_SUBPROTON_DUMMY.Fill(evt.th_q)
              H_ph_recoil_SUBPROTON_DUMMY.Fill(evt.ph_recoil)
              H_th_recoil_SUBPROTON_DUMMY.Fill(evt.th_recoil)

              H_pmiss_SUBPROTON_DUMMY.Fill(evt.pmiss)	
              H_emiss_SUBPROTON_DUMMY.Fill(evt.emiss)	
              #H_emiss_SUBPROTON_DUMMY.Fill(evt.emiss_nuc)
              H_pmx_SUBPROTON_DUMMY.Fill(evt.pmx)
              H_pmy_SUBPROTON_DUMMY.Fill(evt.pmy)
              H_pmz_SUBPROTON_DUMMY.Fill(evt.pmz)
              H_Q2_SUBPROTON_DUMMY.Fill(evt.Q2)
              H_t_SUBPROTON_DUMMY.Fill(-evt.MandelT)
              H_W_SUBPROTON_DUMMY.Fill(evt.W)
              H_epsilon_SUBPROTON_DUMMY.Fill(evt.epsilon)
              H_MM_SUBPROTON_DUMMY.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_SUBPROTON_DUMMY.Fill(pow(evt.MM, 2))  
              #H_MM_SUBPROTON_DUMMY.Fill(evt.Mrecoil)

        ###################################################################################################################################################    
        # Fill random histograms for various trees called above

        print("\nGrabbing %s proton random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPROTON_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPROTON_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ssxfp_SUBPROTON_RAND.Fill(evt.ssxfp)
              H_ssyfp_SUBPROTON_RAND.Fill(evt.ssyfp)
              H_ssxpfp_SUBPROTON_RAND.Fill(evt.ssxpfp)
              H_ssypfp_SUBPROTON_RAND.Fill(evt.ssypfp)
              H_ssdelta_SUBPROTON_RAND.Fill(evt.ssdelta)
              H_ssxptar_SUBPROTON_RAND.Fill(evt.ssxptar)
              H_ssyptar_SUBPROTON_RAND.Fill(evt.ssyptar)

              H_hsxfp_SUBPROTON_RAND.Fill(evt.hsxfp)
              H_hsyfp_SUBPROTON_RAND.Fill(evt.hsyfp)
              H_hsxpfp_SUBPROTON_RAND.Fill(evt.hsxpfp)
              H_hsypfp_SUBPROTON_RAND.Fill(evt.hsypfp)
              H_hsdelta_SUBPROTON_RAND.Fill(evt.hsdelta)
              H_hsxptar_SUBPROTON_RAND.Fill(evt.hsxptar)	
              H_hsyptar_SUBPROTON_RAND.Fill(evt.hsyptar)

              H_pmiss_SUBPROTON_RAND.Fill(evt.pmiss)	
              H_emiss_SUBPROTON_RAND.Fill(evt.emiss)	
              #H_emiss_SUBPROTON_RAND.Fill(evt.emiss_nuc)
              H_pmx_SUBPROTON_RAND.Fill(evt.pmx)
              H_pmy_SUBPROTON_RAND.Fill(evt.pmy)
              H_pmz_SUBPROTON_RAND.Fill(evt.pmz)
              H_Q2_SUBPROTON_RAND.Fill(evt.Q2)
              H_t_SUBPROTON_RAND.Fill(-evt.MandelT)
              H_W_SUBPROTON_RAND.Fill(evt.W)
              H_epsilon_SUBPROTON_RAND.Fill(evt.epsilon)
              H_MM_SUBPROTON_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ###################################################################################################################################################    
        # Fill dummy random histograms for various trees called above

        print("\nGrabbing %s proton dummy random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SUBPROTON_DUMMY_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_SUBPROTON_DUMMY_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              H_ssxfp_SUBPROTON_DUMMY_RAND.Fill(evt.ssxfp)
              H_ssyfp_SUBPROTON_DUMMY_RAND.Fill(evt.ssyfp)
              H_ssxpfp_SUBPROTON_DUMMY_RAND.Fill(evt.ssxpfp)
              H_ssypfp_SUBPROTON_DUMMY_RAND.Fill(evt.ssypfp)
              H_ssdelta_SUBPROTON_DUMMY_RAND.Fill(evt.ssdelta)
              H_ssxptar_SUBPROTON_DUMMY_RAND.Fill(evt.ssxptar)
              H_ssyptar_SUBPROTON_DUMMY_RAND.Fill(evt.ssyptar)

              H_hsxfp_SUBPROTON_DUMMY_RAND.Fill(evt.hsxfp)
              H_hsyfp_SUBPROTON_DUMMY_RAND.Fill(evt.hsyfp)
              H_hsxpfp_SUBPROTON_DUMMY_RAND.Fill(evt.hsxpfp)
              H_hsypfp_SUBPROTON_DUMMY_RAND.Fill(evt.hsypfp)
              H_hsdelta_SUBPROTON_DUMMY_RAND.Fill(evt.hsdelta)
              H_hsxptar_SUBPROTON_DUMMY_RAND.Fill(evt.hsxptar)	
              H_hsyptar_SUBPROTON_DUMMY_RAND.Fill(evt.hsyptar)

              H_pmiss_SUBPROTON_DUMMY_RAND.Fill(evt.pmiss)	
              H_emiss_SUBPROTON_DUMMY_RAND.Fill(evt.emiss)	
              #H_emiss_SUBPROTON_DUMMY_RAND.Fill(evt.emiss_nuc)
              H_pmx_SUBPROTON_DUMMY_RAND.Fill(evt.pmx)
              H_pmy_SUBPROTON_DUMMY_RAND.Fill(evt.pmy)
              H_pmz_SUBPROTON_DUMMY_RAND.Fill(evt.pmz)
              H_Q2_SUBPROTON_DUMMY_RAND.Fill(evt.Q2)
              H_t_SUBPROTON_DUMMY_RAND.Fill(-evt.MandelT)
              H_W_SUBPROTON_DUMMY_RAND.Fill(evt.W)
              H_epsilon_SUBPROTON_DUMMY_RAND.Fill(evt.epsilon)
              H_MM_SUBPROTON_DUMMY_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )
              
        ################################################################################################################################################
        # Normalize simc by normfactor/nevents
        # Normalize dummy by effective charge and target correction
        # Normalize data by effective charge
        
        normfac_simc = (simc_normfactor)/(simc_nevents)
        H_ssxfp_SIMC.Scale(normfac_simc)
        H_ssyfp_SIMC.Scale(normfac_simc)
        H_ssxpfp_SIMC.Scale(normfac_simc)
        H_ssypfp_SIMC.Scale(normfac_simc)
        H_hsxfp_SIMC.Scale(normfac_simc)
        H_hsyfp_SIMC.Scale(normfac_simc)
        H_hsxpfp_SIMC.Scale(normfac_simc)
        H_hsypfp_SIMC.Scale(normfac_simc)
        H_ssxptar_SIMC.Scale(normfac_simc)
        H_ssyptar_SIMC.Scale(normfac_simc)
        H_hsxptar_SIMC.Scale(normfac_simc)
        H_hsyptar_SIMC.Scale(normfac_simc)
        H_ssdelta_SIMC.Scale(normfac_simc)
        H_hsdelta_SIMC.Scale(normfac_simc)
        H_Q2_SIMC.Scale(normfac_simc)
        H_t_SIMC.Scale(normfac_simc)
        H_epsilon_SIMC.Scale(normfac_simc)
        H_MM_SIMC.Scale(normfac_simc)
        H_ph_q_SIMC.Scale(normfac_simc)
        H_th_q_SIMC.Scale(normfac_simc)
        H_ph_recoil_SIMC.Scale(normfac_simc)
        H_th_recoil_SIMC.Scale(normfac_simc)
        H_pmiss_SIMC.Scale(normfac_simc)
        H_emiss_SIMC.Scale(normfac_simc)
        #H_pmx_SIMC.Scale(normfac_simc)
        #H_pmy_SIMC.Scale(normfac_simc)
        #H_pmz_SIMC.Scale(normfac_simc)
        H_W_SIMC.Scale(normfac_simc)

        dummy_target_corr = 4.8579
        if phi_setting == "Right":
            normfac_dummy = 1/(dummy_charge_right*dummy_target_corr)
            normfac_data = 1/(data_charge_right)
        if phi_setting == "Left":
            normfac_dummy = 1/(dummy_charge_left*dummy_target_corr)
            normfac_data = 1/(data_charge_left)
        if phi_setting == "Center":
            normfac_dummy = 1/(dummy_charge_center*dummy_target_corr)
            normfac_data = 1/(data_charge_center)

        H_ssxfp_DUMMY.Scale(normfac_dummy)
        H_ssyfp_DUMMY.Scale(normfac_dummy)
        H_ssxpfp_DUMMY.Scale(normfac_dummy)
        H_ssypfp_DUMMY.Scale(normfac_dummy)
        H_hsxfp_DUMMY.Scale(normfac_dummy)
        H_hsyfp_DUMMY.Scale(normfac_dummy)
        H_hsxpfp_DUMMY.Scale(normfac_dummy)
        H_hsypfp_DUMMY.Scale(normfac_dummy)
        H_ssxptar_DUMMY.Scale(normfac_dummy)
        H_ssyptar_DUMMY.Scale(normfac_dummy)
        H_hsxptar_DUMMY.Scale(normfac_dummy)
        H_hsyptar_DUMMY.Scale(normfac_dummy)
        H_ssdelta_DUMMY.Scale(normfac_dummy)
        H_hsdelta_DUMMY.Scale(normfac_dummy)
        H_Q2_DUMMY.Scale(normfac_dummy)
        H_t_DUMMY.Scale(normfac_dummy)
        H_epsilon_DUMMY.Scale(normfac_dummy)
        H_MM_DUMMY.Scale(normfac_dummy)
        H_ph_q_DUMMY.Scale(normfac_dummy)
        H_th_q_DUMMY.Scale(normfac_dummy)
        H_ph_recoil_DUMMY.Scale(normfac_dummy)
        H_th_recoil_DUMMY.Scale(normfac_dummy)
        H_pmiss_DUMMY.Scale(normfac_dummy)
        H_emiss_DUMMY.Scale(normfac_dummy)
        H_pmx_DUMMY.Scale(normfac_dummy)
        H_pmy_DUMMY.Scale(normfac_dummy)
        H_pmz_DUMMY.Scale(normfac_dummy)
        H_W_DUMMY.Scale(normfac_dummy)
        H_ct_DUMMY.Scale(normfac_dummy)

        H_ssxfp_DATA.Scale(normfac_data)
        H_ssyfp_DATA.Scale(normfac_data)
        H_ssxpfp_DATA.Scale(normfac_data)
        H_ssypfp_DATA.Scale(normfac_data)
        H_hsxfp_DATA.Scale(normfac_data)
        H_hsyfp_DATA.Scale(normfac_data)
        H_hsxpfp_DATA.Scale(normfac_data)
        H_hsypfp_DATA.Scale(normfac_data)
        H_ssxptar_DATA.Scale(normfac_data)
        H_ssyptar_DATA.Scale(normfac_data)
        H_hsxptar_DATA.Scale(normfac_data)
        H_hsyptar_DATA.Scale(normfac_data)
        H_ssdelta_DATA.Scale(normfac_data)
        H_hsdelta_DATA.Scale(normfac_data)
        H_Q2_DATA.Scale(normfac_data)
        H_t_DATA.Scale(normfac_data)
        H_epsilon_DATA.Scale(normfac_data)
        H_MM_DATA.Scale(normfac_data)
        H_ph_q_DATA.Scale(normfac_data)
        H_th_q_DATA.Scale(normfac_data)
        H_ph_recoil_DATA.Scale(normfac_data)
        H_th_recoil_DATA.Scale(normfac_data)
        H_pmiss_DATA.Scale(normfac_data)
        H_emiss_DATA.Scale(normfac_data)
        H_pmx_DATA.Scale(normfac_data)
        H_pmy_DATA.Scale(normfac_data)
        H_pmz_DATA.Scale(normfac_data)
        H_W_DATA.Scale(normfac_data)
        H_ct_DATA.Scale(normfac_data)

        # Data Random subtraction
        H_ssxfp_RAND.Scale(normfac_data/nWindows)
        H_ssyfp_RAND.Scale(normfac_data/nWindows)
        H_ssxpfp_RAND.Scale(normfac_data/nWindows)
        H_ssypfp_RAND.Scale(normfac_data/nWindows)
        H_hsxfp_RAND.Scale(normfac_data/nWindows)
        H_hsyfp_RAND.Scale(normfac_data/nWindows)
        H_hsxpfp_RAND.Scale(normfac_data/nWindows)
        H_hsypfp_RAND.Scale(normfac_data/nWindows)
        H_ssxptar_RAND.Scale(normfac_data/nWindows)
        H_ssyptar_RAND.Scale(normfac_data/nWindows)
        H_hsxptar_RAND.Scale(normfac_data/nWindows)
        H_hsyptar_RAND.Scale(normfac_data/nWindows)
        H_ssdelta_RAND.Scale(normfac_data/nWindows)
        H_hsdelta_RAND.Scale(normfac_data/nWindows)
        H_Q2_RAND.Scale(normfac_data/nWindows)
        H_t_RAND.Scale(normfac_data/nWindows)
        H_epsilon_RAND.Scale(normfac_data/nWindows)
        H_MM_RAND.Scale(normfac_data/nWindows)
        H_pmiss_RAND.Scale(normfac_data/nWindows)
        H_emiss_RAND.Scale(normfac_data/nWindows)
        H_pmx_RAND.Scale(normfac_data/nWindows)
        H_pmy_RAND.Scale(normfac_data/nWindows)
        H_pmz_RAND.Scale(normfac_data/nWindows)
        H_W_RAND.Scale(normfac_data/nWindows)
        #H_ct_RAND.Scale(normfac_data/nWindows)

        # Dummy Random subtraction
        H_ssxfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssyfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssxpfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssypfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsyfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxpfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsypfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssxptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssyptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsyptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssdelta_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsdelta_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_Q2_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_t_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_epsilon_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_MM_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmiss_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_emiss_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmx_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmy_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmz_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_W_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        #H_ct_DUMMY_RAND.Scale(normfac_dummy/nWindows)

        if phi_setting == "Right":
            normfac_subpion_dummy = 1/(300000)
            normfac_subpion_data = 1/(300000)
        if phi_setting == "Left":
            # 5p5, low
            normfac_subpion_dummy = 1/(300000)
            normfac_subpion_data = 1/(300000)
        if phi_setting == "Center":
            # 5p5, low
            normfac_subpion_dummy = 1/(500000)
            normfac_subpion_data = 1/(500000)

        H_ssxfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssyfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssxpfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssypfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsxfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsyfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsxpfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsypfp_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssxptar_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssyptar_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsxptar_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsyptar_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ssdelta_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_hsdelta_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_Q2_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_t_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_epsilon_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_MM_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ph_q_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_th_q_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ph_recoil_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_th_recoil_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_pmiss_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_emiss_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_pmx_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_pmy_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_pmz_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_W_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        H_ct_epi_SUBPION_DUMMY.Scale(normfac_subpion_dummy)
        
        H_ssxfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssyfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssxpfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssypfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsxfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsyfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsxpfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsypfp_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssxptar_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssyptar_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsxptar_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsyptar_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ssdelta_SUBPION_DATA.Scale(normfac_subpion_data)
        H_hsdelta_SUBPION_DATA.Scale(normfac_subpion_data)
        H_Q2_SUBPION_DATA.Scale(normfac_subpion_data)
        H_t_SUBPION_DATA.Scale(normfac_subpion_data)
        H_epsilon_SUBPION_DATA.Scale(normfac_subpion_data)
        H_MM_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ph_q_SUBPION_DATA.Scale(normfac_subpion_data)
        H_th_q_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ph_recoil_SUBPION_DATA.Scale(normfac_subpion_data)
        H_th_recoil_SUBPION_DATA.Scale(normfac_subpion_data)
        H_pmiss_SUBPION_DATA.Scale(normfac_subpion_data)
        H_emiss_SUBPION_DATA.Scale(normfac_subpion_data)
        H_pmx_SUBPION_DATA.Scale(normfac_subpion_data)
        H_pmy_SUBPION_DATA.Scale(normfac_subpion_data)
        H_pmz_SUBPION_DATA.Scale(normfac_subpion_data)
        H_W_SUBPION_DATA.Scale(normfac_subpion_data)
        H_ct_epi_SUBPION_DATA.Scale(normfac_subpion_data)

        # Data Random subtraction
        H_ssxfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssyfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssxpfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssypfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsxfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsyfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsxpfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsypfp_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssxptar_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssyptar_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsxptar_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsyptar_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_ssdelta_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_hsdelta_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_Q2_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_t_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_epsilon_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_MM_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_pmiss_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_emiss_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_pmx_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_pmy_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_pmz_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        H_W_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)
        #H_ct_epi_SUBPION_RAND.Scale(normfac_subpion_data/nWindows)

        # Dummy Random subtraction
        H_ssxfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssyfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssxpfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssypfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsxfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsyfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsxpfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsypfp_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssxptar_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssyptar_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsxptar_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsyptar_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_ssdelta_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_hsdelta_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_Q2_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_t_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_epsilon_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_MM_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_pmiss_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_emiss_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_pmx_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_pmy_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_pmz_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        H_W_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)
        #H_ct_epi_SUBPION_DUMMY_RAND.Scale(normfac_subpion_dummy/nWindows)

        if phi_setting == "Right":
            normfac_subproton_dummy = 1/(12000)
            normfac_subproton_data = 1/(12000)
        if phi_setting == "Left":
            # 5p5, low
            normfac_subproton_dummy = 1/(12000)
            normfac_subproton_data = 1/(12000)
        if phi_setting == "Center":
            # 5p5, low
            normfac_subproton_dummy = 1/(6500)
            normfac_subproton_data = 1/(6500)

        H_ssxfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssyfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssxpfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssypfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsxfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsyfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsxpfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsypfp_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssxptar_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssyptar_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsxptar_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsyptar_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ssdelta_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_hsdelta_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_Q2_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_t_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_epsilon_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_MM_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ph_q_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_th_q_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ph_recoil_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_th_recoil_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_pmiss_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_emiss_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_pmx_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_pmy_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_pmz_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_W_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)
        H_ct_ep_SUBPROTON_DUMMY.Scale(normfac_subproton_dummy)

        H_ssxfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssyfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssxpfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssypfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsxfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsyfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsxpfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsypfp_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssxptar_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssyptar_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsxptar_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsyptar_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ssdelta_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_hsdelta_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_Q2_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_t_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_epsilon_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_MM_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ph_q_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_th_q_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ph_recoil_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_th_recoil_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_pmiss_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_emiss_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_pmx_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_pmy_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_pmz_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_W_SUBPROTON_DATA.Scale(normfac_subproton_data)
        H_ct_ep_SUBPROTON_DATA.Scale(normfac_subproton_data)

        # Data Random subtraction
        H_ssxfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssyfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssxpfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssypfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsxfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsyfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsxpfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsypfp_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssxptar_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssyptar_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsxptar_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsyptar_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_ssdelta_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_hsdelta_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_Q2_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_t_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_epsilon_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_MM_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_pmiss_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_emiss_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_pmx_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_pmy_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_pmz_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        H_W_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)
        #H_ct_ep_SUBPROTON_RAND.Scale(normfac_subproton_data/nWindows)

        # Dummy Random subtraction
        H_ssxfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssyfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssxpfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssypfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsxfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsyfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsxpfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsypfp_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssxptar_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssyptar_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsxptar_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsyptar_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_ssdelta_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_hsdelta_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_Q2_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_t_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_epsilon_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_MM_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_pmiss_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_emiss_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_pmx_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_pmy_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_pmz_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        H_W_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        #H_ct_ep_SUBPROTON_DUMMY_RAND.Scale(normfac_subproton_dummy/nWindows)
        
        ###
        # Data Random subtraction
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

        ###
        # Data Random subtraction
        H_ssxfp_SUBPION_DATA.Add(H_ssxfp_SUBPION_RAND,-1)
        H_ssyfp_SUBPION_DATA.Add(H_ssyfp_SUBPION_RAND,-1)
        H_ssxpfp_SUBPION_DATA.Add(H_ssxpfp_SUBPION_RAND,-1)
        H_ssypfp_SUBPION_DATA.Add(H_ssypfp_SUBPION_RAND,-1)
        H_hsxfp_SUBPION_DATA.Add(H_hsxfp_SUBPION_RAND,-1)
        H_hsyfp_SUBPION_DATA.Add(H_hsyfp_SUBPION_RAND,-1)
        H_hsxpfp_SUBPION_DATA.Add(H_hsxpfp_SUBPION_RAND,-1)
        H_hsypfp_SUBPION_DATA.Add(H_hsypfp_SUBPION_RAND,-1)
        H_ssxptar_SUBPION_DATA.Add(H_ssxptar_SUBPION_RAND,-1)
        H_ssyptar_SUBPION_DATA.Add(H_ssyptar_SUBPION_RAND,-1)
        H_hsxptar_SUBPION_DATA.Add(H_hsxptar_SUBPION_RAND,-1)
        H_hsyptar_SUBPION_DATA.Add(H_hsyptar_SUBPION_RAND,-1)
        H_ssdelta_SUBPION_DATA.Add(H_ssdelta_SUBPION_RAND,-1)
        H_hsdelta_SUBPION_DATA.Add(H_hsdelta_SUBPION_RAND,-1)
        H_Q2_SUBPION_DATA.Add(H_Q2_SUBPION_RAND,-1)
        H_t_SUBPION_DATA.Add(H_t_SUBPION_RAND,-1)
        H_epsilon_SUBPION_DATA.Add(H_epsilon_SUBPION_RAND,-1)
        H_MM_SUBPION_DATA.Add(H_MM_SUBPION_RAND,-1)
        H_pmiss_SUBPION_DATA.Add(H_pmiss_SUBPION_RAND,-1)
        H_emiss_SUBPION_DATA.Add(H_emiss_SUBPION_RAND,-1)
        H_pmx_SUBPION_DATA.Add(H_pmx_SUBPION_RAND,-1)
        H_pmy_SUBPION_DATA.Add(H_pmy_SUBPION_RAND,-1)
        H_pmz_SUBPION_DATA.Add(H_pmz_SUBPION_RAND,-1)
        H_W_SUBPION_DATA.Add(H_W_SUBPION_RAND,-1)
        H_ct_epi_SUBPION_DATA.Add(H_ct_epi_SUBPION_RAND,-1)

        ###
        # Dummy Random subtraction
        H_ssxfp_SUBPION_DUMMY.Add(H_ssxfp_SUBPION_DUMMY_RAND,-1)
        H_ssyfp_SUBPION_DUMMY.Add(H_ssyfp_SUBPION_DUMMY_RAND,-1)
        H_ssxpfp_SUBPION_DUMMY.Add(H_ssxpfp_SUBPION_DUMMY_RAND,-1)
        H_ssypfp_SUBPION_DUMMY.Add(H_ssypfp_SUBPION_DUMMY_RAND,-1)
        H_hsxfp_SUBPION_DUMMY.Add(H_hsxfp_SUBPION_DUMMY_RAND,-1)
        H_hsyfp_SUBPION_DUMMY.Add(H_hsyfp_SUBPION_DUMMY_RAND,-1)
        H_hsxpfp_SUBPION_DUMMY.Add(H_hsxpfp_SUBPION_DUMMY_RAND,-1)
        H_hsypfp_SUBPION_DUMMY.Add(H_hsypfp_SUBPION_DUMMY_RAND,-1)
        H_ssxptar_SUBPION_DUMMY.Add(H_ssxptar_SUBPION_DUMMY_RAND,-1)
        H_ssyptar_SUBPION_DUMMY.Add(H_ssyptar_SUBPION_DUMMY_RAND,-1)
        H_hsxptar_SUBPION_DUMMY.Add(H_hsxptar_SUBPION_DUMMY_RAND,-1)
        H_hsyptar_SUBPION_DUMMY.Add(H_hsyptar_SUBPION_DUMMY_RAND,-1)
        H_ssdelta_SUBPION_DUMMY.Add(H_ssdelta_SUBPION_DUMMY_RAND,-1)
        H_hsdelta_SUBPION_DUMMY.Add(H_hsdelta_SUBPION_DUMMY_RAND,-1)
        H_Q2_SUBPION_DUMMY.Add(H_Q2_SUBPION_DUMMY_RAND,-1)
        H_t_SUBPION_DUMMY.Add(H_t_SUBPION_DUMMY_RAND,-1)
        H_epsilon_SUBPION_DUMMY.Add(H_epsilon_SUBPION_DUMMY_RAND,-1)
        H_MM_SUBPION_DUMMY.Add(H_MM_SUBPION_DUMMY_RAND,-1)
        H_pmiss_SUBPION_DUMMY.Add(H_pmiss_SUBPION_DUMMY_RAND,-1)
        H_emiss_SUBPION_DUMMY.Add(H_emiss_SUBPION_DUMMY_RAND,-1)
        H_pmx_SUBPION_DUMMY.Add(H_pmx_SUBPION_DUMMY_RAND,-1)
        H_pmy_SUBPION_DUMMY.Add(H_pmy_SUBPION_DUMMY_RAND,-1)
        H_pmz_SUBPION_DUMMY.Add(H_pmz_SUBPION_DUMMY_RAND,-1)
        H_W_SUBPION_DUMMY.Add(H_W_SUBPION_DUMMY_RAND,-1)
        H_ct_epi_SUBPION_DUMMY.Add(H_ct_epi_SUBPION_DUMMY_RAND,-1)

        ###
        # Data Random subtraction
        H_ssxfp_SUBPROTON_DATA.Add(H_ssxfp_SUBPROTON_RAND,-1)
        H_ssyfp_SUBPROTON_DATA.Add(H_ssyfp_SUBPROTON_RAND,-1)
        H_ssxpfp_SUBPROTON_DATA.Add(H_ssxpfp_SUBPROTON_RAND,-1)
        H_ssypfp_SUBPROTON_DATA.Add(H_ssypfp_SUBPROTON_RAND,-1)
        H_hsxfp_SUBPROTON_DATA.Add(H_hsxfp_SUBPROTON_RAND,-1)
        H_hsyfp_SUBPROTON_DATA.Add(H_hsyfp_SUBPROTON_RAND,-1)
        H_hsxpfp_SUBPROTON_DATA.Add(H_hsxpfp_SUBPROTON_RAND,-1)
        H_hsypfp_SUBPROTON_DATA.Add(H_hsypfp_SUBPROTON_RAND,-1)
        H_ssxptar_SUBPROTON_DATA.Add(H_ssxptar_SUBPROTON_RAND,-1)
        H_ssyptar_SUBPROTON_DATA.Add(H_ssyptar_SUBPROTON_RAND,-1)
        H_hsxptar_SUBPROTON_DATA.Add(H_hsxptar_SUBPROTON_RAND,-1)
        H_hsyptar_SUBPROTON_DATA.Add(H_hsyptar_SUBPROTON_RAND,-1)
        H_ssdelta_SUBPROTON_DATA.Add(H_ssdelta_SUBPROTON_RAND,-1)
        H_hsdelta_SUBPROTON_DATA.Add(H_hsdelta_SUBPROTON_RAND,-1)
        H_Q2_SUBPROTON_DATA.Add(H_Q2_SUBPROTON_RAND,-1)
        H_t_SUBPROTON_DATA.Add(H_t_SUBPROTON_RAND,-1)
        H_epsilon_SUBPROTON_DATA.Add(H_epsilon_SUBPROTON_RAND,-1)
        H_MM_SUBPROTON_DATA.Add(H_MM_SUBPROTON_RAND,-1)
        H_pmiss_SUBPROTON_DATA.Add(H_pmiss_SUBPROTON_RAND,-1)
        H_emiss_SUBPROTON_DATA.Add(H_emiss_SUBPROTON_RAND,-1)
        H_pmx_SUBPROTON_DATA.Add(H_pmx_SUBPROTON_RAND,-1)
        H_pmy_SUBPROTON_DATA.Add(H_pmy_SUBPROTON_RAND,-1)
        H_pmz_SUBPROTON_DATA.Add(H_pmz_SUBPROTON_RAND,-1)
        H_W_SUBPROTON_DATA.Add(H_W_SUBPROTON_RAND,-1)
        H_ct_ep_SUBPROTON_DATA.Add(H_ct_ep_SUBPROTON_RAND,-1)

        ###
        # Dummy Random subtraction
        H_ssxfp_SUBPROTON_DUMMY.Add(H_ssxfp_SUBPROTON_DUMMY_RAND,-1)
        H_ssyfp_SUBPROTON_DUMMY.Add(H_ssyfp_SUBPROTON_DUMMY_RAND,-1)
        H_ssxpfp_SUBPROTON_DUMMY.Add(H_ssxpfp_SUBPROTON_DUMMY_RAND,-1)
        H_ssypfp_SUBPROTON_DUMMY.Add(H_ssypfp_SUBPROTON_DUMMY_RAND,-1)
        H_hsxfp_SUBPROTON_DUMMY.Add(H_hsxfp_SUBPROTON_DUMMY_RAND,-1)
        H_hsyfp_SUBPROTON_DUMMY.Add(H_hsyfp_SUBPROTON_DUMMY_RAND,-1)
        H_hsxpfp_SUBPROTON_DUMMY.Add(H_hsxpfp_SUBPROTON_DUMMY_RAND,-1)
        H_hsypfp_SUBPROTON_DUMMY.Add(H_hsypfp_SUBPROTON_DUMMY_RAND,-1)
        H_ssxptar_SUBPROTON_DUMMY.Add(H_ssxptar_SUBPROTON_DUMMY_RAND,-1)
        H_ssyptar_SUBPROTON_DUMMY.Add(H_ssyptar_SUBPROTON_DUMMY_RAND,-1)
        H_hsxptar_SUBPROTON_DUMMY.Add(H_hsxptar_SUBPROTON_DUMMY_RAND,-1)
        H_hsyptar_SUBPROTON_DUMMY.Add(H_hsyptar_SUBPROTON_DUMMY_RAND,-1)
        H_ssdelta_SUBPROTON_DUMMY.Add(H_ssdelta_SUBPROTON_DUMMY_RAND,-1)
        H_hsdelta_SUBPROTON_DUMMY.Add(H_hsdelta_SUBPROTON_DUMMY_RAND,-1)
        H_Q2_SUBPROTON_DUMMY.Add(H_Q2_SUBPROTON_DUMMY_RAND,-1)
        H_t_SUBPROTON_DUMMY.Add(H_t_SUBPROTON_DUMMY_RAND,-1)
        H_epsilon_SUBPROTON_DUMMY.Add(H_epsilon_SUBPROTON_DUMMY_RAND,-1)
        H_MM_SUBPROTON_DUMMY.Add(H_MM_SUBPROTON_DUMMY_RAND,-1)
        H_pmiss_SUBPROTON_DUMMY.Add(H_pmiss_SUBPROTON_DUMMY_RAND,-1)
        H_emiss_SUBPROTON_DUMMY.Add(H_emiss_SUBPROTON_DUMMY_RAND,-1)
        H_pmx_SUBPROTON_DUMMY.Add(H_pmx_SUBPROTON_DUMMY_RAND,-1)
        H_pmy_SUBPROTON_DUMMY.Add(H_pmy_SUBPROTON_DUMMY_RAND,-1)
        H_pmz_SUBPROTON_DUMMY.Add(H_pmz_SUBPROTON_DUMMY_RAND,-1)
        H_W_SUBPROTON_DUMMY.Add(H_W_SUBPROTON_DUMMY_RAND,-1)
        H_ct_ep_SUBPROTON_DUMMY.Add(H_ct_ep_SUBPROTON_DUMMY_RAND,-1)
                
        ###
        # Dummy Subtraction
        H_ssxfp_DATA.Add(H_ssxfp_DUMMY,-1)
        H_ssyfp_DATA.Add(H_ssyfp_DUMMY,-1)
        H_ssxpfp_DATA.Add(H_ssxpfp_DUMMY,-1)
        H_ssypfp_DATA.Add(H_ssypfp_DUMMY,-1)
        H_hsxfp_DATA.Add(H_hsxfp_DUMMY,-1)
        H_hsyfp_DATA.Add(H_hsyfp_DUMMY,-1)
        H_hsxpfp_DATA.Add(H_hsxpfp_DUMMY,-1)
        H_hsypfp_DATA.Add(H_hsypfp_DUMMY,-1)
        H_ssxptar_DATA.Add(H_ssxptar_DUMMY,-1)
        H_ssyptar_DATA.Add(H_ssyptar_DUMMY,-1)
        H_hsxptar_DATA.Add(H_hsxptar_DUMMY,-1)
        H_hsyptar_DATA.Add(H_hsyptar_DUMMY,-1)
        H_ssdelta_DATA.Add(H_ssdelta_DUMMY,-1)
        H_hsdelta_DATA.Add(H_hsdelta_DUMMY,-1)
        H_Q2_DATA.Add(H_Q2_DUMMY,-1)
        H_t_DATA.Add(H_t_DUMMY,-1)
        H_epsilon_DATA.Add(H_epsilon_DUMMY,-1)
        H_MM_DATA.Add(H_MM_DUMMY,-1)
        H_pmiss_DATA.Add(H_pmiss_DUMMY,-1)
        H_emiss_DATA.Add(H_emiss_DUMMY,-1)
        H_pmx_DATA.Add(H_pmx_DUMMY,-1)
        H_pmy_DATA.Add(H_pmy_DUMMY,-1)
        H_pmz_DATA.Add(H_pmz_DUMMY,-1)
        H_W_DATA.Add(H_W_DUMMY,-1)
        H_ct_DATA.Add(H_ct_DUMMY,-1)

        ###
        # Dummy Subtraction
        H_ssxfp_SUBPION_DATA.Add(H_ssxfp_SUBPION_DUMMY,-1)
        H_ssyfp_SUBPION_DATA.Add(H_ssyfp_SUBPION_DUMMY,-1)
        H_ssxpfp_SUBPION_DATA.Add(H_ssxpfp_SUBPION_DUMMY,-1)
        H_ssypfp_SUBPION_DATA.Add(H_ssypfp_SUBPION_DUMMY,-1)
        H_hsxfp_SUBPION_DATA.Add(H_hsxfp_SUBPION_DUMMY,-1)
        H_hsyfp_SUBPION_DATA.Add(H_hsyfp_SUBPION_DUMMY,-1)
        H_hsxpfp_SUBPION_DATA.Add(H_hsxpfp_SUBPION_DUMMY,-1)
        H_hsypfp_SUBPION_DATA.Add(H_hsypfp_SUBPION_DUMMY,-1)
        H_ssxptar_SUBPION_DATA.Add(H_ssxptar_SUBPION_DUMMY,-1)
        H_ssyptar_SUBPION_DATA.Add(H_ssyptar_SUBPION_DUMMY,-1)
        H_hsxptar_SUBPION_DATA.Add(H_hsxptar_SUBPION_DUMMY,-1)
        H_hsyptar_SUBPION_DATA.Add(H_hsyptar_SUBPION_DUMMY,-1)
        H_ssdelta_SUBPION_DATA.Add(H_ssdelta_SUBPION_DUMMY,-1)
        H_hsdelta_SUBPION_DATA.Add(H_hsdelta_SUBPION_DUMMY,-1)
        H_Q2_SUBPION_DATA.Add(H_Q2_SUBPION_DUMMY,-1)
        H_t_SUBPION_DATA.Add(H_t_SUBPION_DUMMY,-1)
        H_epsilon_SUBPION_DATA.Add(H_epsilon_SUBPION_DUMMY,-1)
        H_MM_SUBPION_DATA.Add(H_MM_SUBPION_DUMMY,-1)
        H_pmiss_SUBPION_DATA.Add(H_pmiss_SUBPION_DUMMY,-1)
        H_emiss_SUBPION_DATA.Add(H_emiss_SUBPION_DUMMY,-1)
        H_pmx_SUBPION_DATA.Add(H_pmx_SUBPION_DUMMY,-1)
        H_pmy_SUBPION_DATA.Add(H_pmy_SUBPION_DUMMY,-1)
        H_pmz_SUBPION_DATA.Add(H_pmz_SUBPION_DUMMY,-1)
        H_W_SUBPION_DATA.Add(H_W_SUBPION_DUMMY,-1)
        H_ct_epi_SUBPION_DATA.Add(H_ct_epi_SUBPION_DUMMY,-1)
        
        ###
        # Dummy Subtraction
        H_ssxfp_SUBPROTON_DATA.Add(H_ssxfp_SUBPROTON_DUMMY,-1)
        H_ssyfp_SUBPROTON_DATA.Add(H_ssyfp_SUBPROTON_DUMMY,-1)
        H_ssxpfp_SUBPROTON_DATA.Add(H_ssxpfp_SUBPROTON_DUMMY,-1)
        H_ssypfp_SUBPROTON_DATA.Add(H_ssypfp_SUBPROTON_DUMMY,-1)
        H_hsxfp_SUBPROTON_DATA.Add(H_hsxfp_SUBPROTON_DUMMY,-1)
        H_hsyfp_SUBPROTON_DATA.Add(H_hsyfp_SUBPROTON_DUMMY,-1)
        H_hsxpfp_SUBPROTON_DATA.Add(H_hsxpfp_SUBPROTON_DUMMY,-1)
        H_hsypfp_SUBPROTON_DATA.Add(H_hsypfp_SUBPROTON_DUMMY,-1)
        H_ssxptar_SUBPROTON_DATA.Add(H_ssxptar_SUBPROTON_DUMMY,-1)
        H_ssyptar_SUBPROTON_DATA.Add(H_ssyptar_SUBPROTON_DUMMY,-1)
        H_hsxptar_SUBPROTON_DATA.Add(H_hsxptar_SUBPROTON_DUMMY,-1)
        H_hsyptar_SUBPROTON_DATA.Add(H_hsyptar_SUBPROTON_DUMMY,-1)
        H_ssdelta_SUBPROTON_DATA.Add(H_ssdelta_SUBPROTON_DUMMY,-1)
        H_hsdelta_SUBPROTON_DATA.Add(H_hsdelta_SUBPROTON_DUMMY,-1)
        H_Q2_SUBPROTON_DATA.Add(H_Q2_SUBPROTON_DUMMY,-1)
        H_t_SUBPROTON_DATA.Add(H_t_SUBPROTON_DUMMY,-1)
        H_epsilon_SUBPROTON_DATA.Add(H_epsilon_SUBPROTON_DUMMY,-1)
        H_MM_SUBPROTON_DATA.Add(H_MM_SUBPROTON_DUMMY,-1)
        H_pmiss_SUBPROTON_DATA.Add(H_pmiss_SUBPROTON_DUMMY,-1)
        H_emiss_SUBPROTON_DATA.Add(H_emiss_SUBPROTON_DUMMY,-1)
        H_pmx_SUBPROTON_DATA.Add(H_pmx_SUBPROTON_DUMMY,-1)
        H_pmy_SUBPROTON_DATA.Add(H_pmy_SUBPROTON_DUMMY,-1)
        H_pmz_SUBPROTON_DATA.Add(H_pmz_SUBPROTON_DUMMY,-1)
        H_W_SUBPROTON_DATA.Add(H_W_SUBPROTON_DUMMY,-1)
        H_ct_ep_SUBPROTON_DATA.Add(H_ct_ep_SUBPROTON_DUMMY,-1)

        H_MM_DATA_nosub = H_MM_DATA.Clone("H_MM_DATA_nosub")
        
        ###
        # Pion Subtraction
        H_ssxfp_DATA.Add(H_ssxfp_SUBPION_DATA,-1)
        H_ssyfp_DATA.Add(H_ssyfp_SUBPION_DATA,-1)
        H_ssxpfp_DATA.Add(H_ssxpfp_SUBPION_DATA,-1)
        H_ssypfp_DATA.Add(H_ssypfp_SUBPION_DATA,-1)
        H_hsxfp_DATA.Add(H_hsxfp_SUBPION_DATA,-1)
        H_hsyfp_DATA.Add(H_hsyfp_SUBPION_DATA,-1)
        H_hsxpfp_DATA.Add(H_hsxpfp_SUBPION_DATA,-1)
        H_hsypfp_DATA.Add(H_hsypfp_SUBPION_DATA,-1)
        H_ssxptar_DATA.Add(H_ssxptar_SUBPION_DATA,-1)
        H_ssyptar_DATA.Add(H_ssyptar_SUBPION_DATA,-1)
        H_hsxptar_DATA.Add(H_hsxptar_SUBPION_DATA,-1)
        H_hsyptar_DATA.Add(H_hsyptar_SUBPION_DATA,-1)
        H_ssdelta_DATA.Add(H_ssdelta_SUBPION_DATA,-1)
        H_hsdelta_DATA.Add(H_hsdelta_SUBPION_DATA,-1)
        H_Q2_DATA.Add(H_Q2_SUBPION_DATA,-1)
        H_t_DATA.Add(H_t_SUBPION_DATA,-1)
        H_epsilon_DATA.Add(H_epsilon_SUBPION_DATA,-1)
        H_MM_DATA.Add(H_MM_SUBPION_DATA,-1)
        H_pmiss_DATA.Add(H_pmiss_SUBPION_DATA,-1)
        H_emiss_DATA.Add(H_emiss_SUBPION_DATA,-1)
        H_pmx_DATA.Add(H_pmx_SUBPION_DATA,-1)
        H_pmy_DATA.Add(H_pmy_SUBPION_DATA,-1)
        H_pmz_DATA.Add(H_pmz_SUBPION_DATA,-1)
        H_W_DATA.Add(H_W_SUBPION_DATA,-1)
        H_ct_DATA.Add(H_ct_epi_SUBPION_DATA,-1)
        
        ###
        # Proton Subtraction
        H_ssxfp_DATA.Add(H_ssxfp_SUBPROTON_DATA,-1)
        H_ssyfp_DATA.Add(H_ssyfp_SUBPROTON_DATA,-1)
        H_ssxpfp_DATA.Add(H_ssxpfp_SUBPROTON_DATA,-1)
        H_ssypfp_DATA.Add(H_ssypfp_SUBPROTON_DATA,-1)
        H_hsxfp_DATA.Add(H_hsxfp_SUBPROTON_DATA,-1)
        H_hsyfp_DATA.Add(H_hsyfp_SUBPROTON_DATA,-1)
        H_hsxpfp_DATA.Add(H_hsxpfp_SUBPROTON_DATA,-1)
        H_hsypfp_DATA.Add(H_hsypfp_SUBPROTON_DATA,-1)
        H_ssxptar_DATA.Add(H_ssxptar_SUBPROTON_DATA,-1)
        H_ssyptar_DATA.Add(H_ssyptar_SUBPROTON_DATA,-1)
        H_hsxptar_DATA.Add(H_hsxptar_SUBPROTON_DATA,-1)
        H_hsyptar_DATA.Add(H_hsyptar_SUBPROTON_DATA,-1)
        H_ssdelta_DATA.Add(H_ssdelta_SUBPROTON_DATA,-1)
        H_hsdelta_DATA.Add(H_hsdelta_SUBPROTON_DATA,-1)
        H_Q2_DATA.Add(H_Q2_SUBPROTON_DATA,-1)
        H_t_DATA.Add(H_t_SUBPROTON_DATA,-1)
        H_epsilon_DATA.Add(H_epsilon_SUBPROTON_DATA,-1)
        H_MM_DATA.Add(H_MM_SUBPROTON_DATA,-1)
        H_pmiss_DATA.Add(H_pmiss_SUBPROTON_DATA,-1)
        H_emiss_DATA.Add(H_emiss_SUBPROTON_DATA,-1)
        H_pmx_DATA.Add(H_pmx_SUBPROTON_DATA,-1)
        H_pmy_DATA.Add(H_pmy_SUBPROTON_DATA,-1)
        H_pmz_DATA.Add(H_pmz_SUBPROTON_DATA,-1)
        H_W_DATA.Add(H_W_SUBPROTON_DATA,-1)
        H_ct_DATA.Add(H_ct_ep_SUBPROTON_DATA,-1)
        
        histDict["InFile_DATA"] = InFile_DATA
        histDict["InFile_DUMMY"] = InFile_DUMMY
        histDict["InFile_SIMC"] = InFile_SIMC
        histDict["InFile_SUBPION_DATA"] = InFile_SUBPION_DATA
        histDict["InFile_SUBPION_DUMMY"] = InFile_SUBPION_DUMMY
        histDict["InFile_SUBPROTON_DATA"] = InFile_SUBPROTON_DATA
        histDict["InFile_SUBPROTON_DUMMY"] = InFile_SUBPROTON_DUMMY
        histDict["phi_setting"] = phi_setting
        histDict["pid_text"] = pid_text
        histDict["runNums"] = runNums.split(' ')
        histDict["yieldTree"] = ROOT.TTree("{}".format(phi_setting), "{} Yields".format(phi_setting))
        histDict["InData_efficiency"] = InData_efficiency.split(' ')
        histDict["G_data_eff"] = G_data_eff
        histDict["normfac_data"] = normfac_data
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
        histDict["H_phibins_DATA"] = H_phibins_DATA
        histDict["H_tbins_DATA"] = H_tbins_DATA
        histDict["H_yield_DATA"] = H_yield_DATA
        histDict["normfac_simc"] = normfac_simc
        histDict["H_hsdelta_SIMC"] =     H_hsdelta_SIMC
        histDict["H_hsxptar_SIMC"] =     H_hsxptar_SIMC
        histDict["H_hsyptar_SIMC"] =     H_hsyptar_SIMC
        histDict["H_ssxfp_SIMC"] =     H_ssxfp_SIMC  
        histDict["H_ssyfp_SIMC"] =     H_ssyfp_SIMC  
        histDict["H_ssxpfp_SIMC"] =     H_ssxpfp_SIMC 
        histDict["H_ssypfp_SIMC"] =     H_ssypfp_SIMC 
        histDict["H_hsxfp_SIMC"] =     H_hsxfp_SIMC  
        histDict["H_hsyfp_SIMC"] =     H_hsyfp_SIMC  
        histDict["H_hsxpfp_SIMC"] =     H_hsxpfp_SIMC 
        histDict["H_hsypfp_SIMC"] =     H_hsypfp_SIMC 
        histDict["H_ssdelta_SIMC"] =     H_ssdelta_SIMC
        histDict["H_ssxptar_SIMC"] =     H_ssxptar_SIMC
        histDict["H_ssyptar_SIMC"] =     H_ssyptar_SIMC
        histDict["H_q_SIMC"] =     H_q_SIMC      
        histDict["H_Q2_SIMC"] =     H_Q2_SIMC     
        histDict["H_t_SIMC"] =     H_t_SIMC     
        histDict["H_epsilon_SIMC"] =     H_epsilon_SIMC
        histDict["H_MM_SIMC"] =     H_MM_SIMC
        histDict["H_th_SIMC"] =     H_th_SIMC
        histDict["H_ph_SIMC"] =     H_ph_SIMC
        histDict["H_ph_q_SIMC"] =     H_ph_q_SIMC
        histDict["H_th_q_SIMC"] =     H_th_q_SIMC
        histDict["H_ph_recoil_SIMC"] =     H_ph_recoil_SIMC
        histDict["H_th_recoil_SIMC"] =     H_th_recoil_SIMC
        histDict["H_pmiss_SIMC"] =     H_pmiss_SIMC
        histDict["H_emiss_SIMC"] =     H_emiss_SIMC
        histDict["H_pmx_SIMC"] =     H_pmx_SIMC
        histDict["H_pmy_SIMC"] =     H_pmy_SIMC
        histDict["H_pmz_SIMC"] =     H_pmz_SIMC
        histDict["H_W_SIMC"] =     H_W_SIMC
        histDict["H_yield_SIMC"] = H_yield_SIMC
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
        histDict["yieldDictData"] = {}
        histDict["yieldDictSimc"] = {}
        
        # Add t-binned histograms to dictionary
        histDict.update(tbinDict)

        for key,val in histDict.items():
            if val == H_MM_DATA:
                if isinstance(val, ROOT.TH1D):
                    print("##### key->", key, " val->", type(val), val.Print())
            if val == H_tbins_DATA:
                if isinstance(val, ROOT.TH1D):
                    print("##### key->", key, " val->", type(val), val.Print())


        ###
        # CT plots
        ct = TCanvas()
        l_ct = ROOT.TLegend(0.115,0.65,0.33,0.95)
        l_ct.SetTextSize(0.0235)
        H_ct_DATA.SetLineColor(2)
        H_ct_DATA.Draw("same, HIST")
        H_ct_epi_SUBPION_DATA.SetLineColor(3)
        H_ct_epi_SUBPION_DATA.Draw("same, HIST")
        H_ct_ep_SUBPROTON_DATA.SetLineColor(6)
        H_ct_ep_SUBPROTON_DATA.Draw("same, HIST")
        l_ct.AddEntry(H_ct_DATA,"Kaon")
        l_ct.AddEntry(H_ct_epi_SUBPION_DATA,"Pion")
        l_ct.AddEntry(H_ct_ep_SUBPROTON_DATA,"Proton")
        l_ct.Draw()

        ct.Print(outputpdf.replace("kaon_","{}_kaon_MM_subtract_".format(phi_setting))+'(')

                    
        ###
        # PID Plots
        c_pid = TCanvas()
        
        c_pid.Divide(2,3)

        c_pid.cd(1)
        gPad.SetLogy()

        H_cal_etottracknorm_DATA.SetLineColor(2)
        H_cal_etottracknorm_DATA.Draw("same, HIST")
        H_cal_etottracknorm_SUBPION_DATA.SetLineColor(3)
        H_cal_etottracknorm_SUBPION_DATA.Draw("same, HIST")
        H_cal_etottracknorm_SUBPROTON_DATA.SetLineColor(6)
        H_cal_etottracknorm_SUBPROTON_DATA.Draw("same, HIST")

        c_pid.cd(2)
        gPad.SetLogy()
    
        H_cer_npeSum_DATA.SetLineColor(2)
        H_cer_npeSum_DATA.Draw("same, HIST")
        H_cer_npeSum_SUBPION_DATA.SetLineColor(3)
        H_cer_npeSum_SUBPION_DATA.Draw("same, HIST")
        H_cer_npeSum_SUBPROTON_DATA.SetLineColor(6)
        H_cer_npeSum_SUBPROTON_DATA.Draw("same, HIST")

        c_pid.cd(3)
        gPad.SetLogy()
    
        P_cal_etottracknorm_DATA.SetLineColor(2)
        P_cal_etottracknorm_DATA.Draw("same, HIST")
        P_cal_etottracknorm_SUBPION_DATA.SetLineColor(3)
        P_cal_etottracknorm_SUBPION_DATA.Draw("same, HIST")
        P_cal_etottracknorm_SUBPROTON_DATA.SetLineColor(6)
        P_cal_etottracknorm_SUBPROTON_DATA.Draw("same, HIST")

        c_pid.cd(4)
        gPad.SetLogy()
    
        P_hgcer_npeSum_DATA.SetLineColor(2)
        P_hgcer_npeSum_DATA.Draw("same, HIST")
        P_hgcer_npeSum_SUBPION_DATA.SetLineColor(3)
        P_hgcer_npeSum_SUBPION_DATA.Draw("same, HIST")
        P_hgcer_npeSum_SUBPROTON_DATA.SetLineColor(6)
        P_hgcer_npeSum_SUBPROTON_DATA.Draw("same, HIST")

        c_pid.cd(5)
        gPad.SetLogy()

        P_aero_npeSum_DATA.SetLineColor(2)
        P_aero_npeSum_DATA.Draw("same, HIST")
        P_aero_npeSum_SUBPION_DATA.SetLineColor(3)
        P_aero_npeSum_SUBPION_DATA.Draw("same, HIST")
        P_aero_npeSum_SUBPROTON_DATA.SetLineColor(6)
        P_aero_npeSum_SUBPROTON_DATA.Draw("same, HIST")

        c_pid.Draw()

        c_pid.Print(outputpdf.replace("kaon_","{}_kaon_MM_subtract_".format(phi_setting)))

        ###
        # Plot MM for each particle type
        cmm = TCanvas()
        l_mm = ROOT.TLegend(0.115,0.65,0.33,0.95)
        l_mm.SetTextSize(0.0235)
        H_MM_DATA_nosub.SetLineColor(4)
        histDict["H_MM_DATA"].SetLineColor(2)
        H_MM_SUBPION_DATA.SetLineColor(3)
        H_MM_SUBPROTON_DATA.SetLineColor(6)
        histDict["H_MM_DATA"].SetLineWidth(3)
        H_MM_SUBPION_DATA.SetLineWidth(3)
        H_MM_SUBPROTON_DATA.SetLineWidth(3)
        l_mm.AddEntry(H_MM_DATA_nosub,"Kaon (no sub)")
        l_mm.AddEntry(histDict["H_MM_DATA"],"Kaon")
        l_mm.AddEntry(H_MM_SUBPION_DATA,"Pion")
        l_mm.AddEntry(H_MM_SUBPROTON_DATA,"Proton")
        H_MM_DATA_nosub.Draw("same, HIST")
        histDict["H_MM_DATA"].Draw("same, HIST")
        H_MM_SUBPION_DATA.Draw("same, HIST")
        H_MM_SUBPROTON_DATA.Draw("same, HIST")
        H_MM_DATA_nosub.SetFillColor(4)
        H_MM_DATA_nosub.SetTitle("%s M_{K}" % phi_setting) 
        l_mm.Draw()
        cmm.Draw()
        cmm.Print(outputpdf.replace("kaon_","{}_kaon_MM_subtract_".format(phi_setting))+')')
        cmm.Print(outputpdf+'(')
        
        print("@@@@@@@@@@@",histDict["H_tbins_DATA"])
        
        print("$$$$$$$$$$$$$$$$$$$",type(H_tbins_DATA))
        return histDict
        #return H_tbins_DATA
        
    else:    

        ################################################################################################################################################
        # Define simc root file trees of interest

        # Names don't match so need to do some string rearrangement
        InSIMCFilename = "Prod_Coin_{}.root".format(kinematics[0]+phi_setting.lower()+"_"+kinematics[1])
        rootFileSimc = OUTPATH+"/"+InSIMCFilename
        if not os.path.isfile(rootFileSimc):
            print("\n\nERROR: No simc file found called {}\n\n".format(rootFileSimc))
            return histDict

        InFile_SIMC = ROOT.TFile.Open(rootFileSimc, "OPEN")

        TBRANCH_SIMC  = InFile_SIMC.Get("h10")

        ###############################################################################################################################################

        # Grabs simc number of events and normalizaton factor
        simc_hist = rootFileSimc.replace('.root','.hist')
        f_simc = open(simc_hist)
        for line in f_simc:
            #print(line)
            if "Ngen" in line:
                val = line.split("=")
                simc_nevents = int(val[1])
            if "normfac" in line:
                val = line.split("=")
                simc_normfactor = float(val[1])
        if 'simc_nevents' and 'simc_normfactor' not in locals():
            print("\n\nERROR: Invalid simc hist file %s\n\n" % simc_hist)
            sys.exit(1)
        f_simc.close()    

        ################################################################################################################################################
        # Define data root file trees of interest

        rootFileData = OUTPATH + "/" + ParticleType + "_" + InDATAFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileData):
            print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
            return histDict  

        InFile_DATA = ROOT.TFile.Open(rootFileData, "OPEN")

        TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize())) # No RF unless replayed with pion/proton
        
        TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize())) # No RF unless replayed with pion/proton

        ################################################################################################################################################
        # Define dummy root file trees of interest

        rootFileDummy = OUTPATH + "/" + ParticleType + "_" + InDUMMYFilename + "_%s.root" % (phi_setting)
        if not os.path.isfile(rootFileDummy):
            print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
            return histDict

        InFile_DUMMY = ROOT.TFile.Open(rootFileDummy, "OPEN")  

        TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_RF".format(ParticleType.capitalize())) # No RF unless replayed with pion/proton

        TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_RF".format(ParticleType.capitalize())) # No RF unless replayed with pion/proton

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
                                if "coin_epi_cut_prompt_RF" in line:
                                    pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                    break
                                if "coin_ep_cut_prompt_RF" in line:
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
                                if "coin_epi_cut_prompt_RF" in line:
                                    pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                    break
                                if "coin_ep_cut_prompt_RF" in line:
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
                                if "coin_epi_cut_prompt_RF" in line:
                                    pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                    break
                                if "coin_ep_cut_prompt_RF" in line:
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

        ################################################################################################################################################
        # Grab and calculate efficiency

        sys.path.append('../')
        from getDataTable import calculate_effError

        tot_effError_data = [calculate_effError(run,efficiency_table) for run in runNums.split(' ')]
        #print(InData_efficiency)
        #print(tot_effError_data)
        eff_errProp_data = sum(tot_effError_data) # Error propagation for addition

        print("\n\nTotal Data Efficiency Uncertainty =",eff_errProp_data)

        # Define total efficiency vs run number plots
        G_data_eff = ROOT.TGraphErrors(len(InData_efficiency.split(' ')), np.array([float(x) for x in runNums.split(' ')]),np.array([float(x) for x in InData_efficiency.split(' ')]),np.array([0]*len(tot_effError_data)),np.array(tot_effError_data)*np.array([float(x) for x in InData_efficiency.split(' ')]))

        ###############################################################################################################################################
        # Grab windows for random subtraction

        # Section for grabing Prompt/Random selection parameters from PARAM file
        PARAMPATH = "%s/DB/PARAM" % UTILPATH
        print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], LTANAPATH))
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

        H_hsdelta_SIMC  = ROOT.TH1D("H_hsdelta_SIMC","HMS Delta", 200, -20.0, 20.0)
        H_hsxptar_SIMC  = ROOT.TH1D("H_hsxptar_SIMC","HMS xptar", 200, -0.1, 0.1)
        H_hsyptar_SIMC  = ROOT.TH1D("H_hsyptar_SIMC","HMS yptar", 200, -0.1, 0.1)
        H_ssxfp_SIMC    = ROOT.TH1D("H_ssxfp_SIMC","SHMS xfp", 200, -25.0, 25.0)
        H_ssyfp_SIMC    = ROOT.TH1D("H_ssyfp_SIMC","SHMS yfp", 200, -25.0, 25.0)
        H_ssxpfp_SIMC   = ROOT.TH1D("H_ssxpfp_SIMC","SHMS xpfp", 200, -0.09, 0.09)
        H_ssypfp_SIMC   = ROOT.TH1D("H_ssypfp_SIMC","SHMS ypfp", 200, -0.05, 0.04)
        H_hsxfp_SIMC    = ROOT.TH1D("H_hsxfp_SIMC","HMS xfp", 200, -40.0, 40.0)
        H_hsyfp_SIMC    = ROOT.TH1D("H_hsyfp_SIMC","HMS yfp", 200, -20.0, 20.0)
        H_hsxpfp_SIMC   = ROOT.TH1D("H_hsxpfp_SIMC","HMS xpfp", 200, -0.09, 0.05)
        H_hsypfp_SIMC   = ROOT.TH1D("H_hsypfp_SIMC","HMS ypfp", 200, -0.05, 0.04)
        H_ssdelta_SIMC  = ROOT.TH1D("H_ssdelta_SIMC","SHMS delta", 200, -20.0, 20.0)
        H_ssxptar_SIMC  = ROOT.TH1D("H_ssxptar_SIMC","SHMS xptar", 200, -0.1, 0.1)
        H_ssyptar_SIMC  = ROOT.TH1D("H_ssyptar_SIMC","SHMS yptar", 200, -0.04, 0.04)
        H_q_SIMC        = ROOT.TH1D("H_q_SIMC","q", 200, 0.0, 10.0)
        H_Q2_SIMC       = ROOT.TH1D("H_Q2_SIMC","Q2", 200, Q2min, Q2max)
        H_W_SIMC  = ROOT.TH1D("H_W_SIMC","W ", 200, Wmin, Wmax)
        H_t_SIMC       = ROOT.TH1D("H_t_SIMC","-t", 200, -1.0, 1.5)  
        H_epsilon_SIMC  = ROOT.TH1D("H_epsilon_SIMC","epsilon", 200, 0., 1.0)
        H_MM_SIMC  = ROOT.TH1D("H_MM_SIMC","MM_{K}", 200, 0.0, 1.5)
        H_th_SIMC  = ROOT.TH1D("H_th_SIMC","X' tar", 200, -0.1, 0.1)
        H_ph_SIMC  = ROOT.TH1D("H_ph_SIMC","Y' tar", 200, -0.1, 0.1)
        H_ph_q_SIMC  = ROOT.TH1D("H_ph_q_SIMC","Phi Detected (ph_xq)", 200, -5.0, 5.0)
        H_th_q_SIMC  = ROOT.TH1D("H_th_q_SIMC","Theta Detected (th_xq)", 200, -0.2, 0.2)
        H_ph_recoil_SIMC  = ROOT.TH1D("H_ph_recoil_SIMC","Phi Recoil (ph_bq)", 200, -10.0, 10.0)
        H_th_recoil_SIMC  = ROOT.TH1D("H_th_recoil_SIMC","Theta Recoil (th_bq)", 200, -10.0, 10.0)
        H_pmiss_SIMC  = ROOT.TH1D("H_pmiss_SIMC","pmiss", 200, 0.0, 10.0)
        H_emiss_SIMC  = ROOT.TH1D("H_emiss_SIMC","emiss", 200, 0.0, 10.0)
        H_pmx_SIMC  = ROOT.TH1D("H_pmx_SIMC","pmx", 200, -10.0, 10.0)
        H_pmy_SIMC  = ROOT.TH1D("H_pmy_SIMC","pmy ", 200, -10.0, 10.0)
        H_pmz_SIMC  = ROOT.TH1D("H_pmz_SIMC","pmz", 200, -10.0, 10.0)
        H_yield_SIMC = ROOT.TH1D("H_yield_SIMC", "Simc Yield", NumtBins*NumPhiBins, 0, 1.0)

        H_hsdelta_DATA  = ROOT.TH1D("H_hsdelta_DATA","HMS Delta", 200, -20.0, 20.0)
        H_hsxptar_DATA  = ROOT.TH1D("H_hsxptar_DATA","HMS xptar", 200, -0.1, 0.1)
        H_hsyptar_DATA  = ROOT.TH1D("H_hsyptar_DATA","HMS yptar", 200, -0.1, 0.1)
        H_ssxfp_DATA    = ROOT.TH1D("H_ssxfp_DATA","SHMS xfp", 200, -25.0, 25.0)
        H_ssyfp_DATA    = ROOT.TH1D("H_ssyfp_DATA","SHMS yfp", 200, -25.0, 25.0)
        H_ssxpfp_DATA   = ROOT.TH1D("H_ssxpfp_DATA","SHMS xpfp", 200, -0.09, 0.09)
        H_ssypfp_DATA   = ROOT.TH1D("H_ssypfp_DATA","SHMS ypfp", 200, -0.05, 0.04)
        H_hsxfp_DATA    = ROOT.TH1D("H_hsxfp_DATA","HMS xfp", 200, -40.0, 40.0)
        H_hsyfp_DATA    = ROOT.TH1D("H_hsyfp_DATA","HMS yfp", 200, -20.0, 20.0)
        H_hsxpfp_DATA   = ROOT.TH1D("H_hsxpfp_DATA","HMS xpfp", 200, -0.09, 0.05)
        H_hsypfp_DATA   = ROOT.TH1D("H_hsypfp_DATA","HMS ypfp", 200, -0.05, 0.04)
        H_ssdelta_DATA  = ROOT.TH1D("H_ssdelta_DATA","SHMS delta", 200, -20.0, 20.0)
        H_ssxptar_DATA  = ROOT.TH1D("H_ssxptar_DATA","SHMS xptar", 200, -0.1, 0.1)
        H_ssyptar_DATA  = ROOT.TH1D("H_ssyptar_DATA","SHMS yptar", 200, -0.04, 0.04)
        H_q_DATA        = ROOT.TH1D("H_q_DATA","q", 200, 0.0, 10.0)
        H_Q2_DATA       = ROOT.TH1D("H_Q2_DATA","Q2", 200, Q2min, Q2max)
        H_W_DATA  = ROOT.TH1D("H_W_DATA","W ", 200, Wmin, Wmax)
        H_t_DATA       = ROOT.TH1D("H_t_DATA","-t", 200, -1.0, 1.5)  
        H_epsilon_DATA  = ROOT.TH1D("H_epsilon_DATA","epsilon", 200, 0., 1.0)
        H_MM_DATA  = ROOT.TH1D("H_MM_DATA","MM_{K}", 200, 0.0, 1.5)
        H_th_DATA  = ROOT.TH1D("H_th_DATA","X' tar", 200, -0.1, 0.1)
        H_ph_DATA  = ROOT.TH1D("H_ph_DATA","Y' tar", 200, -0.1, 0.1)
        H_ph_q_DATA  = ROOT.TH1D("H_ph_q_DATA","Phi Detected (ph_xq)", 200, -5.0, 5.0)
        H_th_q_DATA  = ROOT.TH1D("H_th_q_DATA","Theta Detected (th_xq)", 200, -0.2, 0.2)
        H_ph_recoil_DATA  = ROOT.TH1D("H_ph_recoil_DATA","Phi Recoil (ph_bq)", 200, -10.0, 10.0)
        H_th_recoil_DATA  = ROOT.TH1D("H_th_recoil_DATA","Theta Recoil (th_bq)", 200, -10.0, 10.0)
        H_pmiss_DATA  = ROOT.TH1D("H_pmiss_DATA","pmiss", 200, 0.0, 10.0)
        H_emiss_DATA  = ROOT.TH1D("H_emiss_DATA","emiss", 200, 0.0, 10.0)
        H_pmx_DATA  = ROOT.TH1D("H_pmx_DATA","pmx", 200, -10.0, 10.0)
        H_pmy_DATA  = ROOT.TH1D("H_pmy_DATA","pmy ", 200, -10.0, 10.0)
        H_pmz_DATA  = ROOT.TH1D("H_pmz_DATA","pmz", 200, -10.0, 10.0)
        if ParticleType == "pion":
            H_ct_DATA = ROOT.TH1D("H_ct_DATA", "Electron-Pion CTime", 200, -10, 10)
        if ParticleType == "proton":
            H_ct_DATA = ROOT.TH1D("H_ct_DATA", "Electron-Proton CTime", 200, -10, 10)            
        H_cal_etottracknorm_DATA = ROOT.TH1D("H_cal_etottracknorm_DATA", "HMS Cal etottracknorm", 200, 0.2, 1.8)
        H_cer_npeSum_DATA = ROOT.TH1D("H_cer_npeSum_DATA", "HMS Cer Npe Sum", 200, 0, 30)
        P_cal_etottracknorm_DATA = ROOT.TH1D("P_cal_etottracknorm_DATA", "SHMS Cal etottracknorm", 200, 0, 1)
        P_hgcer_npeSum_DATA = ROOT.TH1D("P_hgcer_npeSum_DATA", "SHMS HGCer Npe Sum", 200, 0, 50)
        P_aero_npeSum_DATA = ROOT.TH1D("P_aero_npeSum_DATA", "SHMS Aero Npe Sum", 200, 0, 50)

        H_hsdelta_DUMMY  = ROOT.TH1D("H_hsdelta_DUMMY","HMS Delta", 200, -20.0, 20.0)
        H_hsxptar_DUMMY  = ROOT.TH1D("H_hsxptar_DUMMY","HMS xptar", 200, -0.1, 0.1)
        H_hsyptar_DUMMY  = ROOT.TH1D("H_hsyptar_DUMMY","HMS yptar", 200, -0.1, 0.1)
        H_ssxfp_DUMMY    = ROOT.TH1D("H_ssxfp_DUMMY","SHMS xfp", 200, -25.0, 25.0)
        H_ssyfp_DUMMY    = ROOT.TH1D("H_ssyfp_DUMMY","SHMS yfp", 200, -25.0, 25.0)
        H_ssxpfp_DUMMY   = ROOT.TH1D("H_ssxpfp_DUMMY","SHMS xpfp", 200, -0.09, 0.09)
        H_ssypfp_DUMMY   = ROOT.TH1D("H_ssypfp_DUMMY","SHMS ypfp", 200, -0.05, 0.04)
        H_hsxfp_DUMMY    = ROOT.TH1D("H_hsxfp_DUMMY","HMS xfp", 200, -40.0, 40.0)
        H_hsyfp_DUMMY    = ROOT.TH1D("H_hsyfp_DUMMY","HMS yfp", 200, -20.0, 20.0)
        H_hsxpfp_DUMMY   = ROOT.TH1D("H_hsxpfp_DUMMY","HMS xpfp", 200, -0.09, 0.05)
        H_hsypfp_DUMMY   = ROOT.TH1D("H_hsypfp_DUMMY","HMS ypfp", 200, -0.05, 0.04)
        H_ssdelta_DUMMY  = ROOT.TH1D("H_ssdelta_DUMMY","SHMS delta", 200, -20.0, 20.0)
        H_ssxptar_DUMMY  = ROOT.TH1D("H_ssxptar_DUMMY","SHMS xptar", 200, -0.1, 0.1)
        H_ssyptar_DUMMY  = ROOT.TH1D("H_ssyptar_DUMMY","SHMS yptar", 200, -0.04, 0.04)
        H_q_DUMMY        = ROOT.TH1D("H_q_DUMMY","q", 200, 0.0, 10.0)
        H_Q2_DUMMY       = ROOT.TH1D("H_Q2_DUMMY","Q2", 200, Q2min, Q2max)
        H_W_DUMMY  = ROOT.TH1D("H_W_DUMMY","W ", 200, Wmin, Wmax)
        H_t_DUMMY       = ROOT.TH1D("H_t_DUMMY","-t", 200, -1.0, 1.5)  
        H_epsilon_DUMMY  = ROOT.TH1D("H_epsilon_DUMMY","epsilon", 200, 0., 1.0)
        H_MM_DUMMY  = ROOT.TH1D("H_MM_DUMMY","MM_{K}", 200, 0.0, 1.5)
        H_th_DUMMY  = ROOT.TH1D("H_th_DUMMY","X' tar", 200, -0.1, 0.1)
        H_ph_DUMMY  = ROOT.TH1D("H_ph_DUMMY","Y' tar", 200, -0.1, 0.1)
        H_ph_q_DUMMY  = ROOT.TH1D("H_ph_q_DUMMY","Phi Detected (ph_xq)", 200, -5.0, 5.0)
        H_th_q_DUMMY  = ROOT.TH1D("H_th_q_DUMMY","Theta Detected (th_xq)", 200, -0.2, 0.2)
        H_ph_recoil_DUMMY  = ROOT.TH1D("H_ph_recoil_DUMMY","Phi Recoil (ph_bq)", 200, -10.0, 10.0)
        H_th_recoil_DUMMY  = ROOT.TH1D("H_th_recoil_DUMMY","Theta Recoil (th_bq)", 200, -10.0, 10.0)
        H_pmiss_DUMMY  = ROOT.TH1D("H_pmiss_DUMMY","pmiss", 200, 0.0, 10.0)
        H_emiss_DUMMY  = ROOT.TH1D("H_emiss_DUMMY","emiss", 200, 0.0, 10.0)
        H_pmx_DUMMY  = ROOT.TH1D("H_pmx_DUMMY","pmx", 200, -10.0, 10.0)
        H_pmy_DUMMY  = ROOT.TH1D("H_pmy_DUMMY","pmy ", 200, -10.0, 10.0)
        H_pmz_DUMMY  = ROOT.TH1D("H_pmz_DUMMY","pmz", 200, -10.0, 10.0)
        if ParticleType == "pion":
            H_ct_DUMMY = ROOT.TH1D("H_ct_DUMMY", "Electron-Pion CTime", 200, -10, 10)
        if ParticleType == "proton":
            H_ct_DUMMY = ROOT.TH1D("H_ct_DUMMY", "Electron-Proton CTime", 200, -10, 10)
            
        H_hsdelta_RAND  = ROOT.TH1D("H_hsdelta_RAND","HMS Delta", 200, -20.0, 20.0)
        H_hsxptar_RAND  = ROOT.TH1D("H_hsxptar_RAND","HMS xptar", 200, -0.1, 0.1)
        H_hsyptar_RAND  = ROOT.TH1D("H_hsyptar_RAND","HMS yptar", 200, -0.1, 0.1)
        H_ssxfp_RAND    = ROOT.TH1D("H_ssxfp_RAND","SHMS xfp", 200, -25.0, 25.0)
        H_ssyfp_RAND    = ROOT.TH1D("H_ssyfp_RAND","SHMS yfp", 200, -25.0, 25.0)
        H_ssxpfp_RAND   = ROOT.TH1D("H_ssxpfp_RAND","SHMS xpfp", 200, -0.09, 0.09)
        H_ssypfp_RAND   = ROOT.TH1D("H_ssypfp_RAND","SHMS ypfp", 200, -0.05, 0.04)
        H_hsxfp_RAND    = ROOT.TH1D("H_hsxfp_RAND","HMS xfp", 200, -40.0, 40.0)
        H_hsyfp_RAND    = ROOT.TH1D("H_hsyfp_RAND","HMS yfp", 200, -20.0, 20.0)
        H_hsxpfp_RAND   = ROOT.TH1D("H_hsxpfp_RAND","HMS xpfp", 200, -0.09, 0.05)
        H_hsypfp_RAND   = ROOT.TH1D("H_hsypfp_RAND","HMS ypfp", 200, -0.05, 0.04)
        H_ssdelta_RAND  = ROOT.TH1D("H_ssdelta_RAND","SHMS delta", 200, -20.0, 20.0)
        H_ssxptar_RAND  = ROOT.TH1D("H_ssxptar_RAND","SHMS xptar", 200, -0.1, 0.1)
        H_ssyptar_RAND  = ROOT.TH1D("H_ssyptar_RAND","SHMS yptar", 200, -0.04, 0.04)
        H_q_RAND        = ROOT.TH1D("H_q_RAND","q", 200, 0.0, 10.0)
        H_Q2_RAND       = ROOT.TH1D("H_Q2_RAND","Q2", 200, Q2min, Q2max)
        H_W_RAND  = ROOT.TH1D("H_W_RAND","W ", 200, Wmin, Wmax)
        H_t_RAND       = ROOT.TH1D("H_t_RAND","-t", 200, -1.0, 1.5)
        H_epsilon_RAND  = ROOT.TH1D("H_epsilon_RAND","epsilon", 200, 0., 1.0)
        H_MM_RAND  = ROOT.TH1D("H_MM_RAND","MM_{K}", 200, 0.0, 1.5)
        H_th_RAND  = ROOT.TH1D("H_th_RAND","X' tar", 200, -0.1, 0.1)
        H_ph_RAND  = ROOT.TH1D("H_ph_RAND","Y' tar", 200, -0.1, 0.1)
        H_ph_q_RAND  = ROOT.TH1D("H_ph_q_RAND","Phi Detected (ph_xq)", 200, -5.0, 5.0)
        H_th_q_RAND  = ROOT.TH1D("H_th_q_RAND","Theta Detected (th_xq)", 200, -0.2, 0.2)
        H_ph_recoil_RAND  = ROOT.TH1D("H_ph_recoil_RAND","Phi Recoil (ph_bq)", 200, -10.0, 10.0)
        H_th_recoil_RAND  = ROOT.TH1D("H_th_recoil_RAND","Theta Recoil (th_bq)", 200, -10.0, 10.0)
        H_pmiss_RAND  = ROOT.TH1D("H_pmiss_RAND","pmiss", 200, 0.0, 10.0)
        H_emiss_RAND  = ROOT.TH1D("H_emiss_RAND","emiss", 200, 0.0, 10.0)
        H_pmx_RAND  = ROOT.TH1D("H_pmx_RAND","pmx", 200, -10.0, 10.0)
        H_pmy_RAND  = ROOT.TH1D("H_pmy_RAND","pmy ", 200, -10.0, 10.0)
        H_pmz_RAND  = ROOT.TH1D("H_pmz_RAND","pmz", 200, -10.0, 10.0)
        if ParticleType == "pion":
            H_ct_RAND = ROOT.TH1D("H_ct_RAND", "Electron-Pion CTime", 200, -10, 10)
        if ParticleType == "proton":
            H_ct_RAND = ROOT.TH1D("H_ct_RAND", "Electron-Proton CTime", 200, -10, 10)            

        H_hsdelta_DUMMY_RAND  = ROOT.TH1D("H_hsdelta_DUMMY_RAND","HMS Delta", 200, -20.0, 20.0)
        H_hsxptar_DUMMY_RAND  = ROOT.TH1D("H_hsxptar_DUMMY_RAND","HMS xptar", 200, -0.1, 0.1)
        H_hsyptar_DUMMY_RAND  = ROOT.TH1D("H_hsyptar_DUMMY_RAND","HMS yptar", 200, -0.1, 0.1)
        H_ssxfp_DUMMY_RAND    = ROOT.TH1D("H_ssxfp_DUMMY_RAND","SHMS xfp", 200, -25.0, 25.0)
        H_ssyfp_DUMMY_RAND    = ROOT.TH1D("H_ssyfp_DUMMY_RAND","SHMS yfp", 200, -25.0, 25.0)
        H_ssxpfp_DUMMY_RAND   = ROOT.TH1D("H_ssxpfp_DUMMY_RAND","SHMS xpfp", 200, -0.09, 0.09)
        H_ssypfp_DUMMY_RAND   = ROOT.TH1D("H_ssypfp_DUMMY_RAND","SHMS ypfp", 200, -0.05, 0.04)
        H_hsxfp_DUMMY_RAND    = ROOT.TH1D("H_hsxfp_DUMMY_RAND","HMS xfp", 200, -40.0, 40.0)
        H_hsyfp_DUMMY_RAND    = ROOT.TH1D("H_hsyfp_DUMMY_RAND","HMS yfp", 200, -20.0, 20.0)
        H_hsxpfp_DUMMY_RAND   = ROOT.TH1D("H_hsxpfp_DUMMY_RAND","HMS xpfp", 200, -0.09, 0.05)
        H_hsypfp_DUMMY_RAND   = ROOT.TH1D("H_hsypfp_DUMMY_RAND","HMS ypfp", 200, -0.05, 0.04)
        H_ssdelta_DUMMY_RAND  = ROOT.TH1D("H_ssdelta_DUMMY_RAND","SHMS delta", 200, -20.0, 20.0)
        H_ssxptar_DUMMY_RAND  = ROOT.TH1D("H_ssxptar_DUMMY_RAND","SHMS xptar", 200, -0.1, 0.1)
        H_ssyptar_DUMMY_RAND  = ROOT.TH1D("H_ssyptar_DUMMY_RAND","SHMS yptar", 200, -0.04, 0.04)
        H_q_DUMMY_RAND        = ROOT.TH1D("H_q_DUMMY_RAND","q", 200, 0.0, 10.0)
        H_Q2_DUMMY_RAND       = ROOT.TH1D("H_Q2_DUMMY_RAND","Q2", 200, Q2min, Q2max)
        H_W_DUMMY_RAND  = ROOT.TH1D("H_W_DUMMY_RAND","W ", 200, Wmin, Wmax)
        H_t_DUMMY_RAND       = ROOT.TH1D("H_t_DUMMY_RAND","-t", 200, -1.0, 1.5)
        H_epsilon_DUMMY_RAND  = ROOT.TH1D("H_epsilon_DUMMY_RAND","epsilon", 200, 0., 1.0)
        H_MM_DUMMY_RAND  = ROOT.TH1D("H_MM_DUMMY_RAND","MM_{K}", 200, 0.0, 1.5)
        H_th_DUMMY_RAND  = ROOT.TH1D("H_th_DUMMY_RAND","X' tar", 200, -0.1, 0.1)
        H_ph_DUMMY_RAND  = ROOT.TH1D("H_ph_DUMMY_RAND","Y' tar", 200, -0.1, 0.1)
        H_ph_q_DUMMY_RAND  = ROOT.TH1D("H_ph_q_DUMMY_RAND","Phi Detected (ph_xq)", 200, -5.0, 5.0)
        H_th_q_DUMMY_RAND  = ROOT.TH1D("H_th_q_DUMMY_RAND","Theta Detected (th_xq)", 200, -0.2, 0.2)
        H_ph_recoil_DUMMY_RAND  = ROOT.TH1D("H_ph_recoil_DUMMY_RAND","Phi Recoil (ph_bq)", 200, -10.0, 10.0)
        H_th_recoil_DUMMY_RAND  = ROOT.TH1D("H_th_recoil_DUMMY_RAND","Theta Recoil (th_bq)", 200, -10.0, 10.0)
        H_pmiss_DUMMY_RAND  = ROOT.TH1D("H_pmiss_DUMMY_RAND","pmiss", 200, 0.0, 10.0)
        H_emiss_DUMMY_RAND  = ROOT.TH1D("H_emiss_DUMMY_RAND","emiss", 200, 0.0, 10.0)
        H_pmx_DUMMY_RAND  = ROOT.TH1D("H_pmx_DUMMY_RAND","pmx", 200, -10.0, 10.0)
        H_pmy_DUMMY_RAND  = ROOT.TH1D("H_pmy_DUMMY_RAND","pmy ", 200, -10.0, 10.0)
        H_pmz_DUMMY_RAND  = ROOT.TH1D("H_pmz_DUMMY_RAND","pmz", 200, -10.0, 10.0)
        if ParticleType == "pion":
            H_ct_DUMMY_RAND = ROOT.TH1D("H_ct_DUMMY_RAND", "Electron-Pion CTime", 200, -10, 10)
        if ParticleType == "proton":
            H_ct_DUMMY_RAND = ROOT.TH1D("H_ct_DUMMY_RAND", "Electron-Proton CTime", 200, -10, 10)            

        ################################################################################################################################################
        # t/phi binned histograms

        H_phibins_DATA = ROOT.TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
        H_tbins_DATA = ROOT.TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)
        H_yield_DATA = ROOT.TH1D("H_yield_DATA", "Data Yield", NumtBins*NumPhiBins, 0, 1.0)

        tbinDict = {}
        for i in range(NumtBins):
            tbinDict["H_Q2_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_Q2_tbin_DATA_{}".format(i+1), "Q2 (t bin {})".format(i+1), 200, Q2min, Q2max)
            tbinDict["H_W_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_W_tbin_DATA_{}".format(i+1), "W (t bin {})".format(i+1), 200, Wmin, Wmax)
            tbinDict["H_t_tbin_DATA_{}".format(i+1)] = ROOT.TH1D("H_t_tbin_DATA_{}".format(i+1), "t (t bin {})".format(i+1), 200, tmin, tmax)

        ################################################################################################################################################
        # 2D histograms

        MM_vs_CoinTime_DATA = ROOT.TH2D("MM_vs_CoinTime_DATA","Missing Mass vs CTime; MM; Coin_Time",100, 0, 2, 100, -2, 2)
        CoinTime_vs_beta_DATA = ROOT.TH2D("CoinTime_vs_beta_DATA", "CTime vs SHMS #beta; Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0, 2)
        MM_vs_beta_DATA = ROOT.TH2D("MM_vs_beta_DATA", "Missing Mass vs SHMS #beta; MM; SHMS_#beta", 100, 0, 2, 200, 0, 2)
        phiq_vs_t_DATA = ROOT.TH2D("phiq_vs_t_DATA","; #phi ;t", 12, -3.14, 3.14, 24, tmin, tmax)
        polar_phiq_vs_t_DATA = ROOT.TGraphPolar()
        Q2_vs_W_DATA = ROOT.TH2D("Q2_vs_W_DATA", "Q^{2} vs W; Q^{2}; W", 200, Q2min, Q2max, 200, Wmin, Wmax)

        ################################################################################################################################################
        # Fill data histograms for various trees called above

        print("\nGrabbing %s simc..." % phi_setting)
        for i,evt in enumerate(TBRANCH_SIMC):

          # Progress bar
          Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

          # Define the acceptance cuts  
          SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
          HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
          if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
              Diamond = True
          else:
              try:
                  Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
              except ZeroDivisionError:
                  Diamond = False

          #........................................

          #Fill SIMC events
          if(HMS_Acceptance & SHMS_Acceptance & Diamond):

              H_ssxfp_SIMC.Fill(evt.ssxfp, evt.Weight)
              H_ssyfp_SIMC.Fill(evt.ssyfp, evt.Weight)
              H_ssxpfp_SIMC.Fill(evt.ssxpfp, evt.Weight)
              H_ssypfp_SIMC.Fill(evt.ssypfp, evt.Weight)
              H_hsxfp_SIMC.Fill(evt.hsxfp, evt.Weight)
              H_hsyfp_SIMC.Fill(evt.hsyfp, evt.Weight)
              H_hsxpfp_SIMC.Fill(evt.hsxpfp, evt.Weight)
              H_hsypfp_SIMC.Fill(evt.hsypfp, evt.Weight)
              H_ssdelta_SIMC.Fill(evt.ssdelta, evt.Weight) 
              H_hsdelta_SIMC.Fill(evt.hsdelta, evt.Weight)	
              H_ssxptar_SIMC.Fill(evt.ssxptar, evt.Weight)
              H_ssyptar_SIMC.Fill(evt.ssyptar, evt.Weight)
              H_hsxptar_SIMC.Fill(evt.hsxptar, evt.Weight)	
              H_hsyptar_SIMC.Fill(evt.hsyptar, evt.Weight)

              H_ph_q_SIMC.Fill(evt.phipq, evt.Weight)
              H_th_q_SIMC.Fill(evt.thetapq, evt.Weight)

              H_pmiss_SIMC.Fill(evt.Pm, evt.Weight)	
              H_emiss_SIMC.Fill(evt.Em, evt.Weight)	
              #H_pmx_SIMC.Fill(evt.Pmx, evt.Weight)
              #H_pmy_SIMC.Fill(evt.Pmy, evt.Weight)
              #H_pmz_SIMC.Fill(evt.Pmz, evt.Weight)
              H_Q2_SIMC.Fill(evt.Q2, evt.Weight)
              H_W_SIMC.Fill(evt.W, evt.Weight)
              H_t_DATA.Fill(evt.t)
              H_epsilon_SIMC.Fill(evt.epsilon, evt.Weight)
              H_MM_SIMC.Fill(np.sqrt(pow(evt.Em, 2) - pow(evt.Pm, 2)), evt.Weight)


        ################################################################################################################################################
        # Fill histograms for various trees called above

        print("\nGrabbing %s data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DATA):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            # Must be outside diamond cuts to avoid weird overflow errors
            polar_phiq_vs_t_DATA.SetPoint(polar_phiq_vs_t_DATA.GetN(), (evt.ph_q+math.pi)*(180/math.pi), abs(evt.MandelT))

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

              MM_vs_CoinTime_DATA.Fill(evt.MM, evt.CTime_ROC1)
              CoinTime_vs_beta_DATA.Fill(evt.CTime_ROC1,evt.P_gtr_beta)
              MM_vs_beta_DATA.Fill(evt.MM,evt.P_gtr_beta)
              phiq_vs_t_DATA.Fill(evt.ph_q, -evt.MandelT)
              #polar_phiq_vs_t_DATA.SetPoint(i, evt.ph_q, -evt.MandelT)          
              Q2_vs_W_DATA.Fill(evt.Q2, evt.W)

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

              H_ph_q_DATA.Fill(evt.ph_q)
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
              H_MM_DATA.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_DATA.Fill(pow(evt.MM, 2))  
              #H_MM_DATA.Fill(evt.Mrecoil)

              H_cal_etottracknorm_DATA.Fill(evt.H_cal_etottracknorm)
              H_cer_npeSum_DATA.Fill(evt.H_cer_npeSum)

              P_cal_etottracknorm_DATA.Fill(evt.P_cal_etottracknorm)
              P_hgcer_npeSum_DATA.Fill(evt.P_hgcer_npeSum)
              P_aero_npeSum_DATA.Fill(evt.P_aero_npeSum)          

        ################################################################################################################################################
        # Fill dummy histograms for various trees called above

        print("\nGrabbing %s dummy..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DUMMY):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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

              H_ph_q_DUMMY.Fill(evt.ph_q)
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
              H_MM_DUMMY.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2)))
              #H_MM_DUMMY.Fill(pow(evt.MM, 2))  
              #H_MM_DUMMY.Fill(evt.Mrecoil)

        ###################################################################################################################################################    
        # Fill random histograms for various trees called above

        print("\nGrabbing %s random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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
              H_MM_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )

        ###################################################################################################################################################    
        # Fill dummy random histograms for various trees called above

        print("\nGrabbing %s dummy random data..." % phi_setting)
        for i,evt in enumerate(TBRANCH_DUMMY_RAND):

            # Progress bar
            Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)

            #CUTs Definations 
            SHMS_FixCut = (evt.P_hod_goodstarttime == 1) & (evt.P_dc_InsideDipoleExit == 1) # & P_hod_betanotrack > 0.5 & P_hod_betanotrack < 1.4
            SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)

            HMS_FixCut = (evt.H_hod_goodstarttime == 1) & (evt.H_dc_InsideDipoleExit == 1)
            HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
            if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
                Diamond = True
            else:
                try:
                    Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
                except ZeroDivisionError:
                    Diamond = False

            if(HMS_FixCut & HMS_Acceptance & SHMS_FixCut & SHMS_Acceptance & Diamond):

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
              H_MM_DUMMY_RAND.Fill(np.sqrt(pow(evt.emiss, 2) - pow(evt.pmiss, 2))  )


        ################################################################################################################################################
        # Normalize simc by normfactor/nevents
        # Normalize dummy by effective charge and target correction
        # Normalize data by effective charge

        normfac_simc = (simc_normfactor)/(simc_nevents)
        H_ssxfp_SIMC.Scale(normfac_simc)
        H_ssyfp_SIMC.Scale(normfac_simc)
        H_ssxpfp_SIMC.Scale(normfac_simc)
        H_ssypfp_SIMC.Scale(normfac_simc)
        H_hsxfp_SIMC.Scale(normfac_simc)
        H_hsyfp_SIMC.Scale(normfac_simc)
        H_hsxpfp_SIMC.Scale(normfac_simc)
        H_hsypfp_SIMC.Scale(normfac_simc)
        H_ssxptar_SIMC.Scale(normfac_simc)
        H_ssyptar_SIMC.Scale(normfac_simc)
        H_hsxptar_SIMC.Scale(normfac_simc)
        H_hsyptar_SIMC.Scale(normfac_simc)
        H_ssdelta_SIMC.Scale(normfac_simc)
        H_hsdelta_SIMC.Scale(normfac_simc)
        H_Q2_SIMC.Scale(normfac_simc)
        H_t_SIMC.Scale(normfac_simc)
        H_epsilon_SIMC.Scale(normfac_simc)
        H_MM_SIMC.Scale(normfac_simc)
        H_ph_q_SIMC.Scale(normfac_simc)
        H_th_q_SIMC.Scale(normfac_simc)
        H_ph_recoil_SIMC.Scale(normfac_simc)
        H_th_recoil_SIMC.Scale(normfac_simc)
        H_pmiss_SIMC.Scale(normfac_simc)
        H_emiss_SIMC.Scale(normfac_simc)
        #H_pmx_SIMC.Scale(normfac_simc)
        #H_pmy_SIMC.Scale(normfac_simc)
        #H_pmz_SIMC.Scale(normfac_simc)
        H_W_SIMC.Scale(normfac_simc)

        dummy_target_corr = 4.8579
        if phi_setting == "Right":
            normfac_dummy = 1/(dummy_charge_right*dummy_target_corr)
            normfac_data = 1/(data_charge_right)
        if phi_setting == "Left":
            normfac_dummy = 1/(dummy_charge_left*dummy_target_corr)
            normfac_data = 1/(data_charge_left)
        if phi_setting == "Center":
            normfac_dummy = 1/(dummy_charge_center*dummy_target_corr)
            normfac_data = 1/(data_charge_center)

        H_ssxfp_DUMMY.Scale(normfac_dummy)
        H_ssyfp_DUMMY.Scale(normfac_dummy)
        H_ssxpfp_DUMMY.Scale(normfac_dummy)
        H_ssypfp_DUMMY.Scale(normfac_dummy)
        H_hsxfp_DUMMY.Scale(normfac_dummy)
        H_hsyfp_DUMMY.Scale(normfac_dummy)
        H_hsxpfp_DUMMY.Scale(normfac_dummy)
        H_hsypfp_DUMMY.Scale(normfac_dummy)
        H_ssxptar_DUMMY.Scale(normfac_dummy)
        H_ssyptar_DUMMY.Scale(normfac_dummy)
        H_hsxptar_DUMMY.Scale(normfac_dummy)
        H_hsyptar_DUMMY.Scale(normfac_dummy)
        H_ssdelta_DUMMY.Scale(normfac_dummy)
        H_hsdelta_DUMMY.Scale(normfac_dummy)
        H_Q2_DUMMY.Scale(normfac_dummy)
        H_t_DUMMY.Scale(normfac_dummy)
        H_epsilon_DUMMY.Scale(normfac_dummy)
        H_MM_DUMMY.Scale(normfac_dummy)
        H_ph_q_DUMMY.Scale(normfac_dummy)
        H_th_q_DUMMY.Scale(normfac_dummy)
        H_ph_recoil_DUMMY.Scale(normfac_dummy)
        H_th_recoil_DUMMY.Scale(normfac_dummy)
        H_pmiss_DUMMY.Scale(normfac_dummy)
        H_emiss_DUMMY.Scale(normfac_dummy)
        H_pmx_DUMMY.Scale(normfac_dummy)
        H_pmy_DUMMY.Scale(normfac_dummy)
        H_pmz_DUMMY.Scale(normfac_dummy)
        H_W_DUMMY.Scale(normfac_dummy)
        H_ct_DUMMY.Scale(normfac_dummy)

        H_ssxfp_DATA.Scale(normfac_data)
        H_ssyfp_DATA.Scale(normfac_data)
        H_ssxpfp_DATA.Scale(normfac_data)
        H_ssypfp_DATA.Scale(normfac_data)
        H_hsxfp_DATA.Scale(normfac_data)
        H_hsyfp_DATA.Scale(normfac_data)
        H_hsxpfp_DATA.Scale(normfac_data)
        H_hsypfp_DATA.Scale(normfac_data)
        H_ssxptar_DATA.Scale(normfac_data)
        H_ssyptar_DATA.Scale(normfac_data)
        H_hsxptar_DATA.Scale(normfac_data)
        H_hsyptar_DATA.Scale(normfac_data)
        H_ssdelta_DATA.Scale(normfac_data)
        H_hsdelta_DATA.Scale(normfac_data)
        H_Q2_DATA.Scale(normfac_data)
        H_t_DATA.Scale(normfac_data)
        H_epsilon_DATA.Scale(normfac_data)
        H_MM_DATA.Scale(normfac_data)
        H_ph_q_DATA.Scale(normfac_data)
        H_th_q_DATA.Scale(normfac_data)
        H_ph_recoil_DATA.Scale(normfac_data)
        H_th_recoil_DATA.Scale(normfac_data)
        H_pmiss_DATA.Scale(normfac_data)
        H_emiss_DATA.Scale(normfac_data)
        H_pmx_DATA.Scale(normfac_data)
        H_pmy_DATA.Scale(normfac_data)
        H_pmz_DATA.Scale(normfac_data)
        H_W_DATA.Scale(normfac_data)
        H_ct_DATA.Scale(normfac_data)

        # Data Random subtraction
        H_ssxfp_RAND.Scale(normfac_data/nWindows)
        H_ssyfp_RAND.Scale(normfac_data/nWindows)
        H_ssxpfp_RAND.Scale(normfac_data/nWindows)
        H_ssypfp_RAND.Scale(normfac_data/nWindows)
        H_hsxfp_RAND.Scale(normfac_data/nWindows)
        H_hsyfp_RAND.Scale(normfac_data/nWindows)
        H_hsxpfp_RAND.Scale(normfac_data/nWindows)
        H_hsypfp_RAND.Scale(normfac_data/nWindows)
        H_ssxptar_RAND.Scale(normfac_data/nWindows)
        H_ssyptar_RAND.Scale(normfac_data/nWindows)
        H_hsxptar_RAND.Scale(normfac_data/nWindows)
        H_hsyptar_RAND.Scale(normfac_data/nWindows)
        H_ssdelta_RAND.Scale(normfac_data/nWindows)
        H_hsdelta_RAND.Scale(normfac_data/nWindows)
        H_Q2_RAND.Scale(normfac_data/nWindows)
        H_t_RAND.Scale(normfac_data/nWindows)
        H_epsilon_RAND.Scale(normfac_data/nWindows)
        H_MM_RAND.Scale(normfac_data/nWindows)
        H_pmiss_RAND.Scale(normfac_data/nWindows)
        H_emiss_RAND.Scale(normfac_data/nWindows)
        H_pmx_RAND.Scale(normfac_data/nWindows)
        H_pmy_RAND.Scale(normfac_data/nWindows)
        H_pmz_RAND.Scale(normfac_data/nWindows)
        H_W_RAND.Scale(normfac_data/nWindows)
        #H_ct_RAND.Scale(normfac_data/nWindows)

        # Dummy Random subtraction
        H_ssxfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssyfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssxpfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssypfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsyfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxpfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsypfp_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssxptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssyptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsxptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsyptar_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_ssdelta_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_hsdelta_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_Q2_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_t_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_epsilon_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_MM_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmiss_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_emiss_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmx_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmy_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_pmz_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        H_W_DUMMY_RAND.Scale(normfac_dummy/nWindows)
        #H_ct_DUMMY_RAND.Scale(normfac_dummy/nWindows)

        ###
        # Data Random subtraction
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

        ###
        # Dummy Subtraction
        H_ssxfp_DATA.Add(H_ssxfp_DUMMY,-1)
        H_ssyfp_DATA.Add(H_ssyfp_DUMMY,-1)
        H_ssxpfp_DATA.Add(H_ssxpfp_DUMMY,-1)
        H_ssypfp_DATA.Add(H_ssypfp_DUMMY,-1)
        H_hsxfp_DATA.Add(H_hsxfp_DUMMY,-1)
        H_hsyfp_DATA.Add(H_hsyfp_DUMMY,-1)
        H_hsxpfp_DATA.Add(H_hsxpfp_DUMMY,-1)
        H_hsypfp_DATA.Add(H_hsypfp_DUMMY,-1)
        H_ssxptar_DATA.Add(H_ssxptar_DUMMY,-1)
        H_ssyptar_DATA.Add(H_ssyptar_DUMMY,-1)
        H_hsxptar_DATA.Add(H_hsxptar_DUMMY,-1)
        H_hsyptar_DATA.Add(H_hsyptar_DUMMY,-1)
        H_ssdelta_DATA.Add(H_ssdelta_DUMMY,-1)
        H_hsdelta_DATA.Add(H_hsdelta_DUMMY,-1)
        H_Q2_DATA.Add(H_Q2_DUMMY,-1)
        H_t_DATA.Add(H_t_DUMMY,-1)
        H_epsilon_DATA.Add(H_epsilon_DUMMY,-1)
        H_MM_DATA.Add(H_MM_DUMMY,-1)
        H_pmiss_DATA.Add(H_pmiss_DUMMY,-1)
        H_emiss_DATA.Add(H_emiss_DUMMY,-1)
        H_pmx_DATA.Add(H_pmx_DUMMY,-1)
        H_pmy_DATA.Add(H_pmy_DUMMY,-1)
        H_pmz_DATA.Add(H_pmz_DUMMY,-1)
        H_W_DATA.Add(H_W_DUMMY,-1)
        H_ct_DATA.Add(H_ct_DUMMY,-1)

        histDict["InFile_DATA"] = InFile_DATA
        histDict["InFile_DUMMY"] = InFile_DUMMY
        histDict["InFile_SIMC"] = InFile_SIMC
        histDict["phi_setting"] = phi_setting
        histDict["pid_text"] = pid_text
        histDict["runNums"] = runNums.split(' ')
        histDict["yieldTree"] = ROOT.TTree("{}".format(phi_setting), "{} Yields".format(phi_setting))
        histDict["InData_efficiency"] = InData_efficiency.split(' ')
        histDict["G_data_eff"] = G_data_eff
        histDict["normfac_data"] = normfac_data
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
        histDict["H_phibins_DATA"] = H_phibins_DATA
        histDict["H_tbins_DATA"] = H_tbins_DATA
        histDict["H_yield_DATA"] = H_yield_DATA
        histDict["normfac_simc"] = normfac_simc
        histDict["H_hsdelta_SIMC"] =     H_hsdelta_SIMC
        histDict["H_hsxptar_SIMC"] =     H_hsxptar_SIMC
        histDict["H_hsyptar_SIMC"] =     H_hsyptar_SIMC
        histDict["H_ssxfp_SIMC"] =     H_ssxfp_SIMC  
        histDict["H_ssyfp_SIMC"] =     H_ssyfp_SIMC  
        histDict["H_ssxpfp_SIMC"] =     H_ssxpfp_SIMC 
        histDict["H_ssypfp_SIMC"] =     H_ssypfp_SIMC 
        histDict["H_hsxfp_SIMC"] =     H_hsxfp_SIMC  
        histDict["H_hsyfp_SIMC"] =     H_hsyfp_SIMC  
        histDict["H_hsxpfp_SIMC"] =     H_hsxpfp_SIMC 
        histDict["H_hsypfp_SIMC"] =     H_hsypfp_SIMC 
        histDict["H_ssdelta_SIMC"] =     H_ssdelta_SIMC
        histDict["H_ssxptar_SIMC"] =     H_ssxptar_SIMC
        histDict["H_ssyptar_SIMC"] =     H_ssyptar_SIMC
        histDict["H_q_SIMC"] =     H_q_SIMC      
        histDict["H_Q2_SIMC"] =     H_Q2_SIMC     
        histDict["H_t_SIMC"] =     H_t_SIMC     
        histDict["H_epsilon_SIMC"] =     H_epsilon_SIMC
        histDict["H_MM_SIMC"] =     H_MM_SIMC
        histDict["H_th_SIMC"] =     H_th_SIMC
        histDict["H_ph_SIMC"] =     H_ph_SIMC
        histDict["H_ph_q_SIMC"] =     H_ph_q_SIMC
        histDict["H_th_q_SIMC"] =     H_th_q_SIMC
        histDict["H_ph_recoil_SIMC"] =     H_ph_recoil_SIMC
        histDict["H_th_recoil_SIMC"] =     H_th_recoil_SIMC
        histDict["H_pmiss_SIMC"] =     H_pmiss_SIMC
        histDict["H_emiss_SIMC"] =     H_emiss_SIMC
        histDict["H_pmx_SIMC"] =     H_pmx_SIMC
        histDict["H_pmy_SIMC"] =     H_pmy_SIMC
        histDict["H_pmz_SIMC"] =     H_pmz_SIMC
        histDict["H_W_SIMC"] =     H_W_SIMC
        histDict["H_yield_SIMC"] = H_yield_SIMC
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
        histDict["yieldDictData"] = {}
        histDict["yieldDictSimc"] = {}

        # Add t-binned histograms to dictionary
        histDict.update(tbinDict)

        return histDict
