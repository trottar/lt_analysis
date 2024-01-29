#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-29 15:22:28 trottar"
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
# Importing utility functions

sys.path.append("../utility")
from utility import write_to_file, replace_line

###############################################################################################################################################

def create_lists(aveDict, ratioDict, histlist, inpDict, phisetlist, output_file_lst):

    ################################################################################################################################################

    kinematics = inpDict["kinematics"] 
    Ws = inpDict["W"]
    Qs = inpDict["Q2"]
    Q2 = float(Qs.replace("p","."))
    W = float(Ws.replace("p","."))
    EPSVAL = float(inpDict["EPSVAL"] )
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
    pThetaValRight = inpDict["pThetaValRight"]
    EbeamValRight = inpDict["EbeamValRight"]
    pThetaValLeft = inpDict["pThetaValLeft"]
    EbeamValLeft = inpDict["EbeamValLeft"]
    pThetaValCenter = inpDict["pThetaValCenter"]
    EbeamValCenter = inpDict["EbeamValCenter"]
    POL = float(inpDict["POL"])
        
    if POL > 0:
        polID = 'pl'
    else:
        polID = 'mn'
        
    for hist in histlist:
        t_bins = hist["t_bins"]
        phi_bins = hist["phi_bins"]
        
    ################################################################################################################################################        

    averQ2_right_data = []
    averW_right_data = []
    avert_right_data = []
    averQ2_left_data = []
    averW_left_data = []
    avert_left_data = []
    averQ2_center_data = []
    averW_center_data = []
    avert_center_data = []

    ratio_right = []
    ratio_left = []
    ratio_center = []

    tbin_right = []
    tbin_left = []
    tbin_center = []

    phibin_right = []
    phibin_left = []
    phibin_center = []
    
    for phiset in phisetlist:
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for data_key_tuple in data_key_tuples:
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            if phiset == "Right":
                averQ2_right_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                averW_right_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avert_right_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
            if phiset == "Left":
                averQ2_left_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                averW_left_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avert_left_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
            if phiset == "Center":
                averQ2_center_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                averW_center_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avert_center_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
                
        key_tuples = list(ratioDict["binned"][phiset]['ratio'])
        for key_tuple in ratioDict["binned"][phiset]['ratio']:
            # Access the nested dictionary using the tuple key
            nested_dict = ratioDict["binned"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            if phiset == "Right":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == ratioDict["binned"]["t_bins"][i]:
                        tbin_right.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == ratioDict["binned"]["phi_bins"][j]:
                        phibin_right.append(k+1)                        
                ratio_right.append(nested_dict['ratio'][key_tuple]["ratio"])
            if phiset == "Left":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == ratioDict["binned"]["t_bins"][i]:
                        tbin_left.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == ratioDict["binned"]["phi_bins"][j]:
                        phibin_left.append(k+1)                        
                ratio_left.append(nested_dict['ratio'][key_tuple]["ratio"])
            if phiset == "Center":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == ratioDict["binned"]["t_bins"][i]:
                        tbin_center.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == ratioDict["binned"]["phi_bins"][j]:
                        phibin_center.append(k+1)                        
                ratio_center.append(nested_dict['ratio'][key_tuple]["ratio"])

    ################################################################################################################################################
    '''
    Need to sort data properly from least to greatest so that the fortran binning stays consistent (ie the lowest tbin corresponds to the first iteration)
    '''

    if float(runNumRight[0]) != 0:
        # Combine data from different lists into tuples
        data_right_tuples = list(zip(avert_right_data, averQ2_right_data, averW_right_data))
        # Sort based on avert data
        sorted_data_right_tuples = sorted(data_right_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avert_right_data, averQ2_right_data, averW_right_data = zip(*sorted_data_right_tuples[:len(avert_right_data)])

    if float(runNumLeft[0]) != 0:        
        # Combine data from different lists into tuples
        data_left_tuples = list(zip(avert_left_data, averQ2_left_data, averW_left_data))
        # Sort based on avert data
        sorted_data_left_tuples = sorted(data_left_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avert_left_data, averQ2_left_data, averW_left_data = zip(*sorted_data_left_tuples[:len(avert_left_data)])

    if float(runNumCenter[0]) != 0:        
        # Combine data from different lists into tuples
        data_center_tuples = list(zip(avert_center_data, averQ2_center_data, averW_center_data))
        # Sort based on avert data
        sorted_data_center_tuples = sorted(data_center_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avert_center_data, averQ2_center_data, averW_center_data = zip(*sorted_data_center_tuples[:len(avert_center_data)])

    if float(runNumRight[0]) != 0:
        # Combine data from different lists into tuples
        data_right_tuples = list(zip(tbin_right, phibin_right, ratio_right))
        # Sort based on tbin and phibin
        sorted_data_right_tuples = sorted(data_right_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_right, phibin_right, ratio_right = zip(*sorted_data_right_tuples[:len(tbin_right)])

    if float(runNumLeft[0]) != 0:        
        # Combine data from different lists into tuples
        data_left_tuples = list(zip(tbin_left, phibin_left, ratio_left))
        # Sort based on tbin and phibin
        sorted_data_left_tuples = sorted(data_left_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_left, phibin_left, ratio_left = zip(*sorted_data_left_tuples[:len(tbin_left)])

    if float(runNumCenter[0]) != 0:        
        # Combine data from different lists into tuples
        data_center_tuples = list(zip(tbin_center, phibin_center, ratio_center))
        # Sort based on tbin and phibin
        sorted_data_center_tuples = sorted(data_center_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_center, phibin_center, ratio_center = zip(*sorted_data_center_tuples[:len(tbin_center)])
        
    ################################################################################################################################################

    # Define thpq vector relative to middle setting
    for phiset in phisetlist:
        if phiset == "Right":
            runNums = np.array([int(x) for x in runNumRight.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                    ebeam_right = float(EbeamValRight[i])
                    break
                else:
                    continue
                
        if phiset == "Left":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                    ebeam_left = float(EbeamValLeft[i])
                    break
                else:
                    continue
                
        if phiset == "Center":
            runNums = np.array([int(x) for x in runNumCenter.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                    ebeam_center = float(EbeamValCenter[i])
                    break
                else:
                    continue


    ################################################################################################################################################

    '''    
    # Define thpq vector relative to middle setting
    for phiset in phisetlist:
        
        if phiset == "Right":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    if EbeamValRight[0] != EbeamValRight[i] or pThetaValRight[0] != pThetaValRight[i]:
                        print("Run {}".format(run))
                        print("{} | E_{} = {}".format(EbeamValRight[0], phiset, EbeamValRight[i]))
                        print("{} | theta_{} = {}".format(pThetaValRight[0], phiset, pThetaValRight[i]))

        if phiset == "Left":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    if EbeamValLeft[0] != EbeamValLeft[i] or pThetaValLeft[0] != pThetaValLeft[i]:
                        print("Run {}".format(run))
                        print("{} | E_{} = {}".format(EbeamValLeft[0], phiset, EbeamValLeft[i]))
                        print("{} | theta_{} = {}".format(pThetaValLeft[0], phiset, pThetaValLeft[i]))

        if phiset == "Center":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = "{}/log/{}_Analysed_Prod_{}_{}.log".format(LTANAPATH,phiset,ParticleType,run)
                if os.path.exists(pid_log):
                    if EbeamValCenter[0] != EbeamValCenter[i] or pThetaValCenter[0] != pThetaValCenter[i]:
                        print("Run {}".format(run))
                        print("{} | E_{} = {}".format(EbeamValCenter[0], phiset, EbeamValCenter[i]))
                        print("{} | theta_{} = {}".format(pThetaValCenter[0], phiset, pThetaValCenter[i]))
    '''
                        
    ################################################################################################################################################

    f_beam = '{}/src/{}/beam/Eb_KLT.dat'.format(LTANAPATH, ParticleType)
    output_file_lst.append(f_beam.split('/src/')[1])
    # Checks if file exists and creates if not
    if not os.path.exists(f_beam):
        open(f_beam, "w").close()

    if float(runNumRight[0]) != 0:
        # First check if line exists
        with open(f_beam, 'r') as f:
            lines = f.readlines()
            try:
                check_line = "{:.3f} {} {} {:.4f}\n".format(ebeam_right, Q2, W, EPSVAL)
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_beam,check_line)
            except ValueError:
                print("Error: Value not float...",type(ebeam_right), type(Q2), type(W), type(EPSVAL))

    if float(runNumLeft[0]) != 0:
        # First check if line exists
        with open(f_beam, 'r') as f:
            lines = f.readlines()
            try:
                check_line = "{:.3f} {} {} {:.4f}\n".format(ebeam_left, Q2, W, EPSVAL)
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_beam,check_line)
            except ValueError:
                print("Error: Value not float...",type(ebeam_left), type(Q2), type(W), type(EPSVAL))

    if float(runNumCenter[0]) != 0:
        # First check if line exists
        with open(f_beam, 'r') as f:
            lines = f.readlines()
            try:
                check_line = "{:.3f} {} {} {:.4f}\n".format(ebeam_center, Q2, W, EPSVAL)
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_beam,check_line)
            except ValueError:
                print("Error: Value not float...",type(ebeam_center), type(Q2), type(W), type(EPSVAL))
                
    ################################################################################################################################################

    f_list_settings = '{}/src/{}/list.settings'.format(LTANAPATH, ParticleType)
    output_file_lst.append(f_list_settings.split('/src/')[1])
    # Checks if file exists and creates if not
    if not os.path.exists(f_list_settings):
        open(f_list_settings, "w").close()

    '''
    if float(runNumRight[0]) != 0:    
        check_line = "{:d} {} {} {:.4f} -{:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_right, tmin, tmax, NumtBins)
        check_kin = ' '.join(check_line.split()[:4])
        inLine = False
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            if check_kin in line:
                    inLine = True
                    f.write(check_line)
                else:
                    f.write(line)                    
        # Append file if line not in file already
        if not inLine:
            write_to_file(f_list_settings,check_line)

    if float(runNumLeft[0]) != 0:    
        check_line = "{:d} {} {} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_left, tmin, tmax, NumtBins)
        check_kin = ' '.join(check_line.split()[:4])
        inLine = False
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
        # Overwrite current lines if already in file, this fixes the times bins or t-range is changed
        with open(f_list_settings, 'w') as f:
            for line in lines:
                if check_kin in line:
                    inLine = True
                    f.write(check_line)
                else:
                    f.write(line)                    
        # Append file if line not in file already
        if not inLine:
            write_to_file(f_list_settings,check_line)

    if float(runNumCenter[0]) != 0:    
        check_line = "{:d} {} {} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_center, tmin, tmax, NumtBins)
        check_kin = ' '.join(check_line.split()[:4])
        inLine = False
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
        # Overwrite current lines if already in file, this fixes the times bins or t-range is changed
        with open(f_list_settings, 'w') as f:
            for line in lines:
                if check_kin in line:
                    inLine = True
                    f.write(check_line)
                else:
                    f.write(line)                    
        # Append file if line not in file already
        if not inLine:
            write_to_file(f_list_settings,check_line)
    '''

    if float(runNumRight[0]) != 0:    
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} -{:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_right, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:4])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    if check_line in line:
                        replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)                        

    if float(runNumLeft[0]) != 0:    
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_left, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:4])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    if check_line in line:
                        print("\n\n!!!!!!!!!",check_line)
                        replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)                        

    if float(runNumCenter[0]) != 0:    
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_center, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:4])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    if check_line in line:
                        replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)                        
                
    ################################################################################################################################################

    if float(runNumRight[0]) != 0:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(averQ2_right_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_right_data[i], 0.001, averW_right_data[i], 0.001, avert_right_data[i], 0.001)
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)

    if float(runNumLeft[0]) != 0:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(averQ2_left_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_left_data[i], 0.001, averW_left_data[i], 0.001, avert_left_data[i], 0.001)
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)

    if float(runNumCenter[0]) != 0:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100)
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(averQ2_center_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(averQ2_center_data[i], 0.001, averW_center_data[i], 0.001, avert_center_data[i], 0.001)
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)
                
    ################################################################################################################################################

    f_aver = '{}/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(f_aver.split('/src/')[1])
    # Open the file in write mode, which creates a new empty file or overwrites the existing one
    open(f_aver, "w").close()

    if float(runNumRight[0]) != 0:        
        # Open a file in read mode
        with open(f_aver, 'r') as f:
            lines = f.readlines()
            for i, ratio in enumerate(ratio_right):
                check_line = "{:.4f} {:.4f} {} {}\n".format(ratio, 0.001, int(phibin_right[i]), int(tbin_right[i]))
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_aver,check_line)

    f_aver = '{}/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(f_aver.split('/src/')[1])
    # Open the file in write mode, which creates a new empty file or overwrites the existing one
    open(f_aver, "w").close()

    if float(runNumLeft[0]) != 0:        
        # Open a file in read mode
        with open(f_aver, 'r') as f:
            lines = f.readlines()
            for i, ratio in enumerate(ratio_left):
                check_line = "{:.4f} {:.4f} {} {}\n".format(ratio, 0.001, int(phibin_left[i]), int(tbin_left[i]))
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_aver,check_line)

    f_aver = '{}/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), Ws.replace("p",""), float(EPSVAL)*100)
    output_file_lst.append(f_aver.split('/src/')[1])
    # Open the file in write mode, which creates a new empty file or overwrites the existing one
    open(f_aver, "w").close()

    if float(runNumCenter[0]) != 0:        
        # Open a file in read mode
        with open(f_aver, 'r') as f:
            lines = f.readlines()
            for i, ratio in enumerate(ratio_center):
                check_line = "{:.4f} {:.4f} {} {}\n".format(ratio, 0.001, int(phibin_center[i]), int(tbin_center[i]))
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_aver,check_line)
                    
    # Check if the last character is a newline and remove if so
    with open(f_aver, "rb") as f:
        f.seek(-1, os.SEEK_END)
        last_char = f.read(1)
        if last_char == b'\n':
            # Remove the final line
            with open(f_aver, "rb+") as f:
                f.seek(-1, os.SEEK_END)
                f.truncate()                    
    ################################################################################################################################################
