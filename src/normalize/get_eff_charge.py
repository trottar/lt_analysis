#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-29 22:51:49 trottar"
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

def get_eff_charge(hist, inpDict):    

    phi_setting = hist["phi_setting"]
    
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
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,run)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_epi_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_right.split(' ')])
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_right.split(' ')])
        
    if phi_setting == "Left":
        runNums = np.array([int(x) for x in runNumLeft.split(' ')])
        for run in runNums:
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,run)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_left.split(' ')])
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_left.split(' ')])
        
    if phi_setting == "Center":
        runNums = np.array([int(x) for x in runNumCenter.split(' ')])
        for run in runNums:
            pid_log = "%s/log/%s_Analysed_Prod_%s_%s.log" % (LTANAPATH,phi_setting,ParticleType,run)
            if os.path.exists(pid_log):
                    with open(pid_log, 'r') as f_log:
                        for line in f_log:
                            if "coin_ek_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                            
                            if "coin_epi_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break
                            if "coin_ep_cut_prompt_RF" in line:
                                pid_text = next(f_log).replace("[","").replace("]","").replace("{","").replace("}","").replace("'","").replace("&",",").split(",")
                                break                                
            else:
                print("WARNING: Run {} does not have a valid PID log!".format(run))
                continue
        InData_efficiency = np.array([float(x) for x in InData_efficiency_center.split(' ')])
        InData_error_efficiency = np.array([float(x) for x in InData_error_efficiency_center.split(' ')])

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
                              np.array([float(x) for x in runNums]), \
                              InData_efficiency, \
                              np.array([0]*len(InData_error_efficiency)), \
                              InData_error_efficiency*InData_efficiency)
    G_data_eff.SetName("G_data_eff")
    
    ################################################################################################################################################    

    ################################################################################################################################################
    # Normalize dummy by effective charge and target correction
    # Normalize data by effective charge
    # SIMC normalize is defined in src/simc_ana/compare_simc.py because of iteration version
    
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
          
    ################################################################################################################################################        

    histDict["InData_efficiency"] = InData_efficiency
    histDict["InData_error_efficiency"] = InData_error_efficiency
    histDict["G_data_eff"] = G_data_eff
    histDict["normfac_data"] = normfac_data
    histDict["normfac_dummy"] = normfac_dummy

    return histDict
