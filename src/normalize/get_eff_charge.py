#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-02 05:29:08 trottar"
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
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
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

################################################################################################################################################

def get_eff_charge(hist, inpDict, all_data=True):    

    phi_setting = hist["phi_setting"]
    if all_data:
        simc_normfactor = hist["simc_normfactor"]
        simc_nevents = hist["simc_nevents"]
    
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
    ParticleType = inpDict["ParticleType"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Grabs PID cut string

    if phi_setting == "Right":
        runNums = np.array([int(x) for x in runNumRight.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_right.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_right.split(' ')], dtype='float64')
        
    if phi_setting == "Left":
        runNums = np.array([int(x) for x in runNumLeft.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_left.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_left.split(' ')], dtype='float64')
        
    if phi_setting == "Center":
        runNums = np.array([int(x) for x in runNumCenter.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_center.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_center.split(' ')], dtype='float64')

    '''
    if 'pid_text' in locals():
        print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
    else:
        print("ERROR: Invalid {} log file {}!".format(phi_setting.lower(),pid_log))
        pid_text = "\nNo {} cuts file found in logs...".format(phi_setting.lower())
    '''
    print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
        
    ################################################################################################################################################
    # Plot calculate efficiency per run

    # Define total efficiency vs run number plots
    G_data_eff = TGraphErrors(len(InData_efficiency), \
                              np.array([float(x) for x in runNums], dtype='float64'), \
                              InData_efficiency, \
                              np.zeros(len(InData_error_efficiency)), \
                              InData_error_efficiency*InData_efficiency)
    G_data_eff.SetName("G_data_eff")    
    
    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge
    # Normalize simc by normfactor/nevents

    dummy_target_corr = 4.8579
    if phi_setting == "Right":
        normfac_dummy = 1/(dummy_charge_right*dummy_target_corr)
        normfac_data = 1/(data_charge_right)
        if all_data:
            normfac_simc = (simc_normfactor)/(simc_nevents)
    if phi_setting == "Left":
        normfac_dummy = 1/(dummy_charge_left*dummy_target_corr)
        normfac_data = 1/(data_charge_left)
        if all_data:
            normfac_simc = (simc_normfactor)/(simc_nevents)
    if phi_setting == "Center":
        normfac_dummy = 1/(dummy_charge_center*dummy_target_corr)
        normfac_data = 1/(data_charge_center)
        if all_data:
            normfac_simc = (simc_normfactor)/(simc_nevents)
        
    print("\n\n{} data normalization: {:.3e}".format(phi_setting, normfac_data))
    print("{} dummy normalization: {:.3e}".format(phi_setting, normfac_dummy))
    if all_data:
        print("{} simc normalization: {:.3e}".format(phi_setting, normfac_simc))

    ################################################################################################################################################        

    histDict["InData_efficiency"] = InData_efficiency
    histDict["InData_error_efficiency"] = InData_error_efficiency
    histDict["G_data_eff"] = G_data_eff
    histDict["normfac_data"] = normfac_data
    histDict["normfac_dummy"] = normfac_dummy
    if all_data:
        histDict["normfac_simc"] = normfac_simc

    inpDict["normfac_data"] = normfac_data
    inpDict["normfac_dummy"] = normfac_dummy
    if all_data:
        inpDict["normfac_simc"] = normfac_simc        
    
    return histDict

def find_events(hist, inpDict):    

    phi_setting = hist["phi_setting"]
    simc_normfactor = hist["simc_normfactor"]
    simc_nevents = hist["simc_nevents"]
    NumEvts_MM_DATA = hist["NumEvts_MM_DATA"]
    NumEvts_MM_DUMMY = hist["NumEvts_MM_DUMMY"]
    NumEvts_MM_SIMC = hist["NumEvts_MM_SIMC"]
    NumEvts_MM_unweighted_SIMC = hist["NumEvts_MM_unweighted_SIMC"]
    
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
    ParticleType = inpDict["ParticleType"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################
    # Define return dictionary of data
    histDict = {}

    ################################################################################################################################################
    # Grabs PID cut string

    if phi_setting == "Right":
        runNums = np.array([int(x) for x in runNumRight.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_right.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_right.split(' ')], dtype='float64')
        
    if phi_setting == "Left":
        runNums = np.array([int(x) for x in runNumLeft.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_left.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_left.split(' ')], dtype='float64')
        
    if phi_setting == "Center":
        runNums = np.array([int(x) for x in runNumCenter.split(' ')])
        for run in runNums:
            pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_noRF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_center.split(' ')], dtype='float64')
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_center.split(' ')], dtype='float64')

    '''
    if 'pid_text' in locals():
        print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
    else:
        print("ERROR: Invalid {} log file {}!".format(phi_setting.lower(),pid_log))
        pid_text = "\nNo {} cuts file found in logs...".format(phi_setting.lower())
    '''
    print('\n\n',phi_setting,'PID Cuts = ',pid_text,'\n\n')
        
    ################################################################################################################################################
    # Plot calculate efficiency per run

    # Define total efficiency vs run number plots
    G_data_eff = TGraphErrors(len(InData_efficiency), \
                              np.array([float(x) for x in runNums], dtype='float64'), \
                              InData_efficiency, \
                              np.zeros(len(InData_error_efficiency)), \
                              InData_error_efficiency*InData_efficiency)
    G_data_eff.SetName("G_data_eff")    
    
    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge
    # Normalize simc by normfactor/nevents

    dummy_target_corr = 4.8579
    if phi_setting == "Right":
        normfac_dummy = 1/(dummy_charge_right*dummy_target_corr)
        normfac_data = 1/(data_charge_right)
        normfac_simc = (simc_normfactor)/(simc_nevents)
    if phi_setting == "Left":
        normfac_dummy = 1/(dummy_charge_left*dummy_target_corr)
        normfac_data = 1/(data_charge_left)
        normfac_simc = (simc_normfactor)/(simc_nevents)
    if phi_setting == "Center":
        normfac_dummy = 1/(dummy_charge_center*dummy_target_corr)
        normfac_data = 1/(data_charge_center)
        normfac_simc = (simc_normfactor)/(simc_nevents)
        
    print("\n\n{} data total number of events: {:.3e}".format(phi_setting, NumEvts_MM_DATA))
    print("{} dummy total number of events: {:.3e}".format(phi_setting, NumEvts_MM_DUMMY))  
    print("{} simc weighted total number of events: {:.3e}".format(phi_setting, NumEvts_MM_SIMC))
    print("{} simc unweighted total number of events: {:.3e}".format(phi_setting, NumEvts_MM_unweighted_SIMC))
    print("\n\n{} data normalization: {:.3e}".format(phi_setting, normfac_data))
    print("{} dummy normalization: {:.3e}".format(phi_setting, normfac_dummy))
    print("{} simc normalization: {:.3e}".format(phi_setting, normfac_simc))

    ################################################################################################################################################        

    histDict["InData_efficiency"] = InData_efficiency
    histDict["InData_error_efficiency"] = InData_error_efficiency
    histDict["G_data_eff"] = G_data_eff
    histDict["normfac_data"] = normfac_data
    histDict["normfac_dummy"] = normfac_dummy
    histDict["normfac_simc"] = normfac_simc
    histDict["NumEvts_MM_DATA"] = NumEvts_MM_DATA    
    histDict["NumEvts_MM_DUMMY"] = NumEvts_MM_DUMMY
    histDict["NumEvts_MM_SIMC"] = NumEvts_MM_SIMC
    
    return histDict
