#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-08-21 23:24:01 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os, math
import ROOT
import numpy as np

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

###############################################################################################################################################

def create_lists(averDict, inpDict):

    ################################################################################################################################################

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

    if ParticleType == "kaon":
        PID = "k"
    elif ParticleType == "pion":
        PID = "pi"

    ################################################################################################################################################        

    phisetlist = ["Center","Left","Right"]
    for phiset in phisetlist:
        for k, data_key_tuple in enumerate(averDict["binned_DATA"][phiset]['t']):
            # Access the nested dictionary using the tuple key
            data_nested_dict = averDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            tbin = averDict["binned_DATA"]["t_bins"][i]
            phibin = averDict["binned_DATA"]["phi_bins"][j]
            if phiset == "Right":
                phibin_right_data = phibin
                tbin_right_data = tbin
                averQ2_right_data = data_nested_dict['Q2'][data_key_tuple]["Q2_arr"]
                averW_right_data = data_nested_dict['W'][data_key_tuple]["W_arr"]
                avert_right_data = data_nested_dict['t'][data_key_tuple]["t_arr"] 
            if phiset == "Left":
                phibin_left_data = phibin
                tbin_left_data = tbin
                averQ2_left_data = data_nested_dict['Q2'][data_key_tuple]["Q2_arr"]
                averW_left_data = data_nested_dict['W'][data_key_tuple]["W_arr"]
                avert_left_data = data_nested_dict['t'][data_key_tuple]["t_arr"]
            if phiset == "Center":
                phibin_center_data = phibin
                tbin_center_data = tbin
                averQ2_center_data = data_nested_dict['Q2'][data_key_tuple]["Q2_arr"]
                averW_center_data = data_nested_dict['W'][data_key_tuple]["W_arr"]
                avert_center_data = data_nested_dict['t'][data_key_tuple]["t_arr"]

    ################################################################################################################################################

    def write_to_file(f_out,line,write_mode='a'):
        # Open a file in append mode
        with open(f_out, write_mode) as f:
            # Write the value of the variable to the file
            f.write(line)

    ################################################################################################################################################

    # Define thpq vector relative to middle setting
    if float(runNumRight[0]) != 0:
        runNums= runNumRight
        for i, run in enumerate(runNumRight):
            runNum = run
            pid_log = "%s/log/Analysed_Prod_%s.log" % (LTANAPATH,runNum)
            if os.path.exists(pid_log):
                thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                #thpq_right = 3.000
                ebeam_right = float(EbeamValRight[i])
                break
            else:
                continue

    if float(runNumLeft[0]) != 0:
        runNums= runNumLeft
        for i, run in enumerate(runNumLeft):
            runNum = run
            pid_log = "%s/log/Analysed_Prod_%s.log" % (LTANAPATH,runNum)
            if os.path.exists(pid_log):
                thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                #thpq_left = 3.000
                ebeam_left = float(EbeamValLeft[i])
                break
            else:
                continue

    if float(runNumCenter[0]) != 0:
        runNums= runNumCenter
        for i, run in enumerate(runNumCenter):
            runNum = run
            pid_log = "%s/log/Analysed_Prod_%s.log" % (LTANAPATH,runNum)
            if os.path.exists(pid_log):
                thpq_center = 0.000
                ebeam_center = float(EbeamValCenter[i])
                break
            else:
                continue

    print(thpq_left)

    ################################################################################################################################################

    f_list_settings = '{}/src/beam/Eb_KLT.dat'.format(LTANAPATH)
    if not os.path.exists(f_list_settings):
        open(f_list_settings, "w").close()

    if float(runNumRight[0]) != 0:
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()    
            check_line = "{:.3f} {} {:.4f}\n".format(ebeam_right, Q2, EPSVAL)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)

    if float(runNumLeft[0]) != 0:
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()    
            check_line = "{:.3f} {} {:.4f}\n".format(ebeam_left, Q2, EPSVAL)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)

    if float(runNumCenter[0]) != 0:
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()        
            check_line = "{:.3f} {} {:.4f}\n".format(ebeam_center, Q2, EPSVAL)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)

    ################################################################################################################################################

    f_list_settings = '{}/src/list.settings'.format(LTANAPATH)
    if not os.path.exists(f_list_settings):
        open(f_list_settings, "w").close()

    if float(runNumRight[0]) != 0:    
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            check_line = "{} {} {:.4f} -{:.3f} {:.3f} {:.3f} {} {}\n".format(POL, Q2, EPSVAL, thpq_right, TMIN, TMAX, NumtBins, Kset)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)

    if float(runNumLeft[0]) != 0:
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()        
            check_line = "{} {} {:.4f} +{:.3f} {:.3f} {:.3f} {} {}\n".format(POL, Q2, EPSVAL, thpq_left, TMIN, TMAX, NumtBins, Kset)
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)

    if float(runNumCenter[0]) != 0:
        # First check if line exists
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()        
            check_line = "{} {} {:.4f} +{:.3f} {:.3f} {:.3f} {} {}\n".format(POL, Q2, EPSVAL, thpq_center, TMIN, TMAX, NumtBins, Kset)
            if check_line not in lines:
                write_to_file(f_list_settings,check_line)            

    ################################################################################################################################################

    if float(runNumRight[0]) != 0:
        f_list = '{}/src/kindata/kindata.{}_{}_{:.0f}_-{}.dat'.format(LTANAPATH, PID, Q2.replace(".",""), float(EPSVAL)*100, int(thpq_right*1000))

        if not os.path.exists(f_list):
            open(f_list, "w").close()    
        # Open a file in read mode
        with open(f_list, 'r') as f:
            lines = f.readlines()[1:-1]
            tbin_right_data.append(0.0)
            for i, Q2val in enumerate(averQ2_right_data):
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_right_data[i], 1.0, averW_right_data[i], 1.0, avert_right_data[i], 1.0)
                write_to_file(f_list,check_line)
            lines = f.readlines()
            first_line = "{:.6f}\n".format(1.0)
            last_line = "{:.6f}".format(1.0)
            lines.insert(0,first_line)
            lines.append(last_line)
            write_to_file(f_list,"".join(lines),write_mode='w')

    if float(runNumLeft[0]) != 0:
        f_list = '{}/src/kindata/kindata.{}_{}_{:.0f}_+{}.dat'.format(LTANAPATH, PID, Q2.replace(".",""), float(EPSVAL)*100, int(thpq_left*1000))

        if not os.path.exists(f_list):
            open(f_list, "w").close()    
        # Open a file in read mode
        with open(f_list, 'r') as f:
            lines = f.readlines()[1:-1]
            tbin_left_data.append(0.0)
            for i, Q2val in enumerate(averQ2_left_data):
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_left_data[i], 1.0, averW_left_data[i], 1.0, avert_left_data[i], 1.0)
                write_to_file(f_list,check_line)
            lines = f.readlines()
            first_line = "{:.6f}\n".format(1.0)
            last_line = "{:.6f}".format(1.0)
            lines.insert(0,first_line)
            lines.append(last_line)
            write_to_file(f_list,"".join(lines),write_mode='w')

    if float(runNumCenter[0]) != 0:
        f_list = '{}/src/kindata/kindata.{}_{}_{:.0f}_+0000.dat'.format(LTANAPATH, PID, Q2.replace(".",""), float(EPSVAL)*100)

        if not os.path.exists(f_list):
            open(f_list, "w").close()    
        # Open a file in read mode
        with open(f_list, 'r') as f:
            lines = f.readlines()[1:-1]
            tbin_center_data.append(0.0)
            for i, Q2val in enumerate(averQ2_center_data):
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_center_data[i], 1.0, averW_center_data[i], 1.0, avert_center_data[i], 1.0)
                write_to_file(f_list,check_line)
            lines = f.readlines()
            first_line = "{:.6f}\n".format(1.0)
            last_line = "{:.6f}".format(1.0)
            lines.insert(0,first_line)
            lines.append(last_line)
            write_to_file(f_list,"".join(lines),write_mode='w')

    ################################################################################################################################################

