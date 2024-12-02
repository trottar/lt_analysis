#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-02 14:33:09 trottar"
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

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

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

def create_lists(aveDict, yieldDict, histlist, inpDict, phisetlist, output_file_lst):

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

    aveQ2_right_data = []
    aveW_right_data = []
    avet_right_data = []
    aveQ2_left_data = []
    aveW_left_data = []
    avet_left_data = []
    aveQ2_center_data = []
    aveW_center_data = []
    avet_center_data = []

    aveQ2_err_right_data = []
    aveW_err_right_data = []
    avet_err_right_data = []
    aveQ2_err_left_data = []
    aveW_err_left_data = []
    avet_err_left_data = []
    aveQ2_err_center_data = []
    aveW_err_center_data = []
    avet_err_center_data = []
    
    yield_data_right = []
    yield_data_left = []
    yield_data_center = []

    yield_data_err_right = []
    yield_data_err_left = []
    yield_data_err_center = []

    yield_simc_right = []
    yield_simc_left = []
    yield_simc_center = []

    yield_simc_err_right = []
    yield_simc_err_left = []
    yield_simc_err_center = []
    
    tbin_data_right = []
    tbin_data_left = []
    tbin_data_center = []

    phibin_data_right = []
    phibin_data_left = []
    phibin_data_center = []

    tbin_simc_right = []
    tbin_simc_left = []
    tbin_simc_center = []

    phibin_simc_right = []
    phibin_simc_left = []
    phibin_simc_center = []
    
    for i,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        if any(value == 0 for value in data_key_tuples):
            zero_indices = [index for index, value in enumerate(data_key_tuples) if value == 0]
            print(f"ERROR: Zeros found for average kinematic values at t={data_key_tuples[zero_indices]}! Try different binning scheme or adjust t-range...")
            sys.exit(2)
        for data_key_tuple in data_key_tuples:
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            if phiset == "Right":
                aveQ2_right_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                aveW_right_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avet_right_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
                aveQ2_err_right_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave_err"])
                aveW_err_right_data.append(data_nested_dict['W'][data_key_tuple]["W_ave_err"])
                avet_err_right_data.append(data_nested_dict['t'][data_key_tuple]["t_ave_err"])                
            if phiset == "Left":
                aveQ2_left_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                aveW_left_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avet_left_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
                aveQ2_err_left_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave_err"])
                aveW_err_left_data.append(data_nested_dict['W'][data_key_tuple]["W_ave_err"])
                avet_err_left_data.append(data_nested_dict['t'][data_key_tuple]["t_ave_err"])                
            if phiset == "Center":
                aveQ2_center_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave"])
                aveW_center_data.append(data_nested_dict['W'][data_key_tuple]["W_ave"])
                avet_center_data.append(data_nested_dict['t'][data_key_tuple]["t_ave"])
                aveQ2_err_center_data.append(data_nested_dict['Q2'][data_key_tuple]["Q2_ave_err"])
                aveW_err_center_data.append(data_nested_dict['W'][data_key_tuple]["W_ave_err"])
                avet_err_center_data.append(data_nested_dict['t'][data_key_tuple]["t_ave_err"])                
                
        key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
        for key_tuple in key_tuples:
            # Access the nested dictionary using the tuple key
            nested_dict = yieldDict["binned_DATA"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            if phiset == "Right":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_DATA"]["t_bins"][i]:
                        tbin_data_right.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_DATA"]["phi_bins"][j]:
                        phibin_data_right.append(k+1)                        
                yield_data_right.append(nested_dict['yield'][key_tuple]["yield"])
                yield_data_err_right.append(nested_dict['yield'][key_tuple]["yield_err"])
            if phiset == "Left":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_DATA"]["t_bins"][i]:
                        tbin_data_left.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_DATA"]["phi_bins"][j]:
                        phibin_data_left.append(k+1)                        
                yield_data_left.append(nested_dict['yield'][key_tuple]["yield"])
                yield_data_err_left.append(nested_dict['yield'][key_tuple]["yield_err"])
            if phiset == "Center":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_DATA"]["t_bins"][i]:
                        tbin_data_center.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_DATA"]["phi_bins"][j]:
                        phibin_data_center.append(k+1)                        
                yield_data_center.append(nested_dict['yield'][key_tuple]["yield"])
                yield_data_err_center.append(nested_dict['yield'][key_tuple]["yield_err"])

        key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
        for key_tuple in key_tuples:
            # Access the nested dictionary using the tuple key
            nested_dict = yieldDict["binned_SIMC"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            if phiset == "Right":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_SIMC"]["t_bins"][i]:
                        tbin_simc_right.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_SIMC"]["phi_bins"][j]:
                        phibin_simc_right.append(k+1)                        
                yield_simc_right.append(nested_dict['yield'][key_tuple]["yield"])
                yield_simc_err_right.append(nested_dict['yield'][key_tuple]["yield_err"])
            if phiset == "Left":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_SIMC"]["t_bins"][i]:
                        tbin_simc_left.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_SIMC"]["phi_bins"][j]:
                        phibin_simc_left.append(k+1)                        
                yield_simc_left.append(nested_dict['yield'][key_tuple]["yield"])
                yield_simc_err_left.append(nested_dict['yield'][key_tuple]["yield_err"])
            if phiset == "Center":
                for k, t_bin in enumerate(t_bins):
                    if t_bin == yieldDict["binned_SIMC"]["t_bins"][i]:
                        tbin_simc_center.append(k+1)
                for k, phi_bin in enumerate(phi_bins):
                    if phi_bin == yieldDict["binned_SIMC"]["phi_bins"][j]:
                        phibin_simc_center.append(k+1)                        
                yield_simc_center.append(nested_dict['yield'][key_tuple]["yield"])
                yield_simc_err_center.append(nested_dict['yield'][key_tuple]["yield_err"])
                
    ################################################################################################################################################
    '''
    Need to sort data properly from least to greatest so that the fortran binning stays consistent (ie the lowest tbin corresponds to the first iteration)
    '''
    
    if "Right" in phisetlist:
        # Combine data from different lists into tuples
        data_right_tuples = list(zip(avet_right_data, aveQ2_right_data, aveW_right_data))
        # Sort based on avet data
        sorted_data_right_tuples = sorted(data_right_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_right_data, aveQ2_right_data, aveW_right_data = zip(*sorted_data_right_tuples[:len(avet_right_data)])

        # Combine data from different lists into tuples
        data_right_tuples = list(zip(avet_err_right_data, aveQ2_err_right_data, aveW_err_right_data))
        # Sort based on avet_err data
        sorted_data_right_tuples = sorted(data_right_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_err_right_data, aveQ2_err_right_data, aveW_err_right_data = zip(*sorted_data_right_tuples[:len(avet_err_right_data)])        
        
    if "Left" in phisetlist: 
        # Combine data from different lists into tuples
        data_left_tuples = list(zip(avet_left_data, aveQ2_left_data, aveW_left_data))
        # Sort based on avet data
        sorted_data_left_tuples = sorted(data_left_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_left_data, aveQ2_left_data, aveW_left_data = zip(*sorted_data_left_tuples[:len(avet_left_data)])

        # Combine data from different lists into tuples
        data_left_tuples = list(zip(avet_err_left_data, aveQ2_err_left_data, aveW_err_left_data))
        # Sort based on avet_err data
        sorted_data_left_tuples = sorted(data_left_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_err_left_data, aveQ2_err_left_data, aveW_err_left_data = zip(*sorted_data_left_tuples[:len(avet_err_left_data)])
        
    if "Center" in phisetlist: 
        # Combine data from different lists into tuples
        data_center_tuples = list(zip(avet_center_data, aveQ2_center_data, aveW_center_data))
        # Sort based on avet data
        sorted_data_center_tuples = sorted(data_center_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_center_data, aveQ2_center_data, aveW_center_data = zip(*sorted_data_center_tuples[:len(avet_center_data)])

        # Combine data from different lists into tuples
        data_center_tuples = list(zip(avet_err_center_data, aveQ2_err_center_data, aveW_err_center_data))
        # Sort based on avet_err data
        sorted_data_center_tuples = sorted(data_center_tuples, key=lambda x: x[0])
        # Extract sorted values back into separate lists
        avet_err_center_data, aveQ2_err_center_data, aveW_err_center_data = zip(*sorted_data_center_tuples[:len(avet_err_center_data)])
        
    if "Right" in phisetlist:
        # Combine data from different lists into tuples
        data_right_tuples = list(zip(tbin_data_right, phibin_data_right, yield_data_right, yield_data_err_right))
        # Sort based on tbin and phibin
        sorted_data_right_tuples = sorted(data_right_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_data_right, phibin_data_right, yield_data_right, yield_data_err_right = zip(*sorted_data_right_tuples[:len(tbin_data_right)])

    if "Left" in phisetlist:
        # Combine data from different lists into tuples
        data_left_tuples = list(zip(tbin_data_left, phibin_data_left, yield_data_left, yield_data_err_left))
        # Sort based on tbin and phibin
        sorted_data_left_tuples = sorted(data_left_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_data_left, phibin_data_left, yield_data_left, yield_data_err_left = zip(*sorted_data_left_tuples[:len(tbin_data_left)])

    if "Center" in phisetlist:
        # Combine data from different lists into tuples
        data_center_tuples = list(zip(tbin_data_center, phibin_data_center, yield_data_center, yield_data_err_center))
        # Sort based on tbin and phibin
        sorted_data_center_tuples = sorted(data_center_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_data_center, phibin_data_center, yield_data_center, yield_data_err_center = zip(*sorted_data_center_tuples[:len(tbin_data_center)])

    if "Right" in phisetlist:
        # Combine simc from different lists into tuples
        simc_right_tuples = list(zip(tbin_simc_right, phibin_simc_right, yield_simc_right, yield_simc_err_right))
        # Sort based on tbin and phibin
        sorted_simc_right_tuples = sorted(simc_right_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_simc_right, phibin_simc_right, yield_simc_right, yield_simc_err_right = zip(*sorted_simc_right_tuples[:len(tbin_simc_right)])

    if "Left" in phisetlist:
        # Combine simc from different lists into tuples
        simc_left_tuples = list(zip(tbin_simc_left, phibin_simc_left, yield_simc_left, yield_simc_err_left))
        # Sort based on tbin and phibin
        sorted_simc_left_tuples = sorted(simc_left_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_simc_left, phibin_simc_left, yield_simc_left, yield_simc_err_left = zip(*sorted_simc_left_tuples[:len(tbin_simc_left)])

    if "Center" in phisetlist:
        # Combine simc from different lists into tuples
        simc_center_tuples = list(zip(tbin_simc_center, phibin_simc_center, yield_simc_center, yield_simc_err_center))
        # Sort based on tbin and phibin
        sorted_simc_center_tuples = sorted(simc_center_tuples, key=lambda x: (x[0], x[1]))
        # Extract sorted values back into separate lists
        tbin_simc_center, phibin_simc_center, yield_simc_center, yield_simc_err_center = zip(*sorted_simc_center_tuples[:len(tbin_simc_center)])
            
    ################################################################################################################################################

    # Define thpq vector relative to middle setting
    for phiset in phisetlist:
        if phiset == "Right":
            runNums = np.array([int(x) for x in runNumRight.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_right = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))))
                    ebeam_right = float(EbeamValRight[i])
                    break
                else:
                    continue
                
        if phiset == "Left":
            runNums = np.array([int(x) for x in runNumLeft.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_left = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))))
                    ebeam_left = float(EbeamValLeft[i])
                    break
                else:
                    continue
                
        if phiset == "Center":
            runNums = np.array([int(x) for x in runNumCenter.split(' ')])
            for i, run in enumerate(runNums):
                pid_log = f"{LTANAPATH}/log/{phi_setting}_{ParticleType}_{run}_{OutFilename.replace('FullAnalysis_','')}.log"
                if os.path.exists(pid_log):
                    thpq_center = float("{:.3f}".format(abs(float(pThetaValCenter[i])-float(pThetaValCenter[i]))))
                    ebeam_center = float(EbeamValCenter[i])
                    break
                else:
                    continue
                        
    ################################################################################################################################################

    f_beam = '{}/src/{}/beam/Eb_KLT.dat'.format(LTANAPATH, ParticleType)
    output_file_lst.append(f_beam.split('/src/')[1])
    # Checks if file exists and creates if not
    if not os.path.exists(f_beam):
        open(f_beam, "w").close()

    if "Right" in phisetlist:
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

    if "Left" in phisetlist:
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

    if "Center" in phisetlist:
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

    if "Right" in phisetlist:
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} -{:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_right, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:5])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)                        

    if "Left" in phisetlist:
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_left, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:5])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)

    if "Center" in phisetlist:
        # Save lines in the file
        with open(f_list_settings, 'r') as f:
            lines = f.readlines()
            inLine = False
            check_line = "{:d} {:.1f} {:.2f} {:.4f} {:.3f} {:.3f} {:.3f} {}\n".format(int(POL), Q2, W, EPSVAL, thpq_center, tmin, tmax, NumtBins)
            check_kin = ' '.join(check_line.split()[:5])
            for i, line in enumerate(lines):
                if check_kin in line:
                    inLine = True
                    replace_line(f_list_settings, i+1, check_line)
            if not inLine:
                write_to_file(f_list_settings, check_line)                        
                
    ################################################################################################################################################

    if "Right" in phisetlist:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(aveQ2_right_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(aveQ2_right_data[i], aveQ2_err_right_data[i], \
                                                                                  aveW_right_data[i], aveW_err_right_data[i], \
                                                                                  avet_right_data[i], avet_err_right_data[i])
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)

    if "Left" in phisetlist:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(aveQ2_left_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(aveQ2_left_data[i], aveQ2_err_left_data[i], \
                                                                                  aveW_left_data[i], aveW_err_left_data[i], \
                                                                                  avet_left_data[i], avet_err_left_data[i])
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)

    if "Center" in phisetlist:
        f_kindata = '{}/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                  Ws.replace("p",""), float(EPSVAL)*100)
        output_file_lst.append(f_kindata.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_kindata, "w").close()

        processed_Q2vals = set()
        # Open a file in read mode
        with open(f_kindata, 'r') as f:
            for i, Q2val in enumerate(aveQ2_center_data):
                if Q2val in processed_Q2vals:
                    continue
                check_line = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(aveQ2_center_data[i], aveQ2_err_center_data[i], \
                                                                                  aveW_center_data[i], aveW_err_center_data[i], \
                                                                                  avet_center_data[i], avet_err_center_data[i])
                write_to_file(f_kindata,check_line)
                processed_Q2vals.add(Q2val)

    ################################################################################################################################################

    if "Right" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_data_right):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_data_right[i], yield_data_err_right[i], \
                                                                                  int(phibin_data_right[i]), int(tbin_data_right[i]))
                write_to_file(f_yield,check_line)

    if "Left" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_data_left):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_data_left[i], yield_data_err_left[i], \
                                                                                  int(phibin_data_left[i]), int(tbin_data_left[i]))
                write_to_file(f_yield,check_line)

    if "Center" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_data.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                  Ws.replace("p",""), float(EPSVAL)*100)
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_data_center):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_data_center[i], yield_data_err_center[i], \
                                                                                  int(phibin_data_center[i]), int(tbin_data_center[i]))
                write_to_file(f_yield,check_line)

    ################################################################################################################################################
    
    if "Right" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_simc.{}_Q{}W{}_{:.0f}_-{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_right*1000))
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_simc_right):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_simc_right[i], yield_simc_err_right[i], \
                                                                                  int(phibin_simc_right[i]), int(tbin_simc_right[i]))
                write_to_file(f_yield,check_line)

    if "Left" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_simc.{}_Q{}W{}_{:.0f}_+{}.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                Ws.replace("p",""), float(EPSVAL)*100, int(thpq_left*1000))
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_simc_left):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_simc_left[i], yield_simc_err_left[i], \
                                                                                  int(phibin_simc_left[i]), int(tbin_simc_left[i]))
                write_to_file(f_yield,check_line)

    if "Center" in phisetlist:
        f_yield = '{}/src/{}/yields/yield_simc.{}_Q{}W{}_{:.0f}_+0000.dat'.format(LTANAPATH, ParticleType, polID, Qs.replace("p",""), \
                                                                                  Ws.replace("p",""), float(EPSVAL)*100)
        output_file_lst.append(f_yield.split('/src/')[1])
        # Open the file in write mode, which creates a new empty file or overwrites the existing one
        open(f_yield, "w").close()

        # Open a file in read mode
        with open(f_yield, 'r') as f:
            for i, yieldval in enumerate(yield_simc_center):
                check_line = "{:.4f} {:.4f} {} {}\n".format(yield_simc_center[i], yield_simc_err_center[i], \
                                                                                  int(phibin_simc_center[i]), int(tbin_simc_center[i]))
                write_to_file(f_yield,check_line)
                                                                    
    ################################################################################################################################################
