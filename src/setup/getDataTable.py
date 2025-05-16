#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-08-12 15:49:29 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import numpy as np
from functools import reduce
import sys,os

###############################################################################################################################################

'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import data_to_csv

################################################################################################################################################
# Grab bcm value

def get_bcm(runNum,efficiency_table,foutcsv):
    
    inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s" % efficiency_table

    # Converts csv data to dataframe
    try:
        eff_data = pd.read_csv(inp_f)
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)
    #print(eff_data.keys())

    # Redefine table to be only the run number of interest
    eff_data = eff_data[eff_data['Run_Number'] == int(runNum)]
    #print(eff_data)
    
    return eff_data["BCM1_Beam_Cut_Charge"].iloc[0]

################################################################################################################################################
# Grab ebeam value

def get_ebeam(runNum,efficiency_table,foutcsv):
    
    inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s" % efficiency_table

    # Converts csv data to dataframe
    try:
        eff_data = pd.read_csv(inp_f)
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)
    #print(eff_data.keys())

    # Redefine table to be only the run number of interest
    eff_data = eff_data[eff_data['Run_Number'] == int(runNum)]
    #print(eff_data)

    return eff_data["Beam_Energy"].iloc[0]

################################################################################################################################################
# Grab pTheta value

def get_pTheta(runNum,efficiency_table,foutcsv):
    
    inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s" % efficiency_table

    # Converts csv data to dataframe
    try:
        eff_data = pd.read_csv(inp_f)
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)
    #print(eff_data.keys())

    # Redefine table to be only the run number of interest
    eff_data = eff_data[eff_data['Run_Number'] == int(runNum)]
    #print(eff_data)

    return eff_data["SHMS_Angle"].iloc[0]

################################################################################################################################################
# Define efficiencies

def get_efficiencies(runNum,efficiency_table,foutcsv):

    inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s" % efficiency_table

    if "heep" in efficiency_table:
        runType = "heep"
    else:
        runType = "production"

    # Converts csv data to dataframe
    try:
        eff_data = pd.read_csv(inp_f)
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)
    #print(eff_data.keys())

    ##############
    # HARD CODED #
    ##############
    
    # Find HMS Cherenkov efficiency and error based off golden singles run (6603) without contamination
    # See https://redmine.jlab.org/attachments/download/1758/KaonLT_Meeting_Feb_23.pdf for more info
    hcer_eff = 0.9730
    hcer_eff_error = 0.0006
    
    # Find HMS Calorimeter efficiency and error based off golden run () without contamination
    # See https://redmine.jlab.org/attachments/download/1758/KaonLT_Meeting_Feb_23.pdf for more info
    hcal_eff = 0.9961
    hcal_eff_error = 0.0005

    ##############
    ##############
    ##############
    
    # Redefine table to be only the run number of interest
    eff_data = eff_data[eff_data['Run_Number'] == int(runNum)]
    #print(eff_data)
    
    if len(eff_data) > 0:  # Check if eff_data is not empty
    
        if runType == "heep":
            # Define dictionary of efficiency values
            effDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time" : eff_data["Non_Scaler_EDTM_Live_Time"].iloc[0],
                # Hodo
                "HMS_Hodo_3_of_4_EFF" : eff_data["HMS_Hodo_3_of_4_EFF"].iloc[0],
                "SHMS_Hodo_3_of_4_EFF" : eff_data["SHMS_Hodo_3_of_4_EFF"].iloc[0],
                ##############
                # HARD CODED #
                ##############
                # HMS Cal, defined above...
                "HMS_Cal_ALL_Elec_Eff" : hcal_eff, 
                # HMS Cer, defined above...
                "HMS_Cer_ALL_Elec_Eff" : hcer_eff,
                ##############
                ##############
                ##############
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF" : eff_data["HMS_Elec_ALL_TRACK_EFF"].iloc[0],
                "SHMS_Prot_ALL_TRACK_EFF" : eff_data["SHMS_Prot_ALL_TRACK_EFF"].iloc[0],
                # Boiling Correction
                "BOIL_Eff" : eff_data["BOIL_Eff"].iloc[0],
            }
            # Define dictionary of efficiency uncertainty values
            efficiency_errDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time_ERROR" : eff_data["Non_Scaler_EDTM_Live_Time_ERROR"].iloc[0],
                ##############
                # HARD CODED #
                ##############
                # Hodo
                #"HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                #"SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                "HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Elec_ALL_TRACK_EFF_ERROR"].iloc[0], # FIX
                "SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Prot_ALL_TRACK_EFF_ERROR"].iloc[0], # FIX
                # HMS Cal, defined above...
                "HMS_Cal_ALL_Elec_Eff_ERROR" : hcal_eff_error,
                # HMS Cer, defined above...
                "HMS_Cer_ALL_Elec_Eff_ERROR" : hcer_eff_error,
                ##############
                ##############
                ##############
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF_ERROR" : eff_data["HMS_Elec_ALL_TRACK_EFF_ERROR"].iloc[0],
                "SHMS_Prot_ALL_TRACK_EFF_ERROR" : eff_data["SHMS_Prot_ALL_TRACK_EFF_ERROR"].iloc[0],
                # Boiling Correction
                "BOIL_Eff_ERROR" : eff_data["BOIL_Eff_ERROR"].iloc[0],
            }

        else:
            # Define dictionary of efficiency values
            effDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time" : eff_data["Non_Scaler_EDTM_Live_Time"].iloc[0],
                # Hodo
                "HMS_Hodo_3_of_4_EFF" : eff_data["HMS_Hodo_3_of_4_EFF"].iloc[0],
                "SHMS_Hodo_3_of_4_EFF" : eff_data["SHMS_Hodo_3_of_4_EFF"].iloc[0],
                ##############
                # HARD CODED #
                ##############
                # HMS Cal, defined above...
                "HMS_Cal_ALL_Elec_Eff" : hcal_eff, 
                # HMS Cer, defined above...
                "HMS_Cer_ALL_Elec_Eff" : hcer_eff,
                ##############
                ##############
                ##############
                # SHMS Aero
                "SHMS_Aero_ALL_Pion_Eff" : eff_data["SHMS_Aero_ALL_Pion_Eff"].iloc[0],
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF" : eff_data["HMS_Elec_ALL_TRACK_EFF"].iloc[0],
                "SHMS_Pion_ALL_TRACK_EFF" : eff_data["SHMS_Pion_ALL_TRACK_EFF"].iloc[0],
                # Boiling Correction
                "BOIL_Eff" : eff_data["BOIL_Eff"].iloc[0],
            }

            # Define dictionary of efficiency uncertainty values
            efficiency_errDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time_ERROR" : eff_data["Non_Scaler_EDTM_Live_Time_ERROR"].iloc[0],                
                ##############
                # HARD CODED #
                ##############                
                # Hodo
                #"HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                #"SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                "HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Elec_ALL_TRACK_EFF_ERROR"].iloc[0], # FIX
                "SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Pion_ALL_TRACK_EFF_ERROR"].iloc[0], # FIX
                # HMS Cal, defined above...
                "HMS_Cal_ALL_Elec_Eff_ERROR" : hcal_eff_error,
                # HMS Cer, defined above...
                "HMS_Cer_ALL_Elec_Eff_ERROR" : hcer_eff_error,
                ##############
                ##############
                ##############
                # SHMS Aero
                "SHMS_Aero_ALL_Pion_Eff_ERROR" : eff_data["SHMS_Aero_ALL_Pion_Eff_ERROR"].iloc[0],
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF_ERROR" : eff_data["HMS_Elec_ALL_TRACK_EFF_ERROR"].iloc[0],
                "SHMS_Pion_ALL_TRACK_EFF_ERROR" : eff_data["SHMS_Pion_ALL_TRACK_EFF_ERROR"].iloc[0],
                # Boiling Correction
                "BOIL_Eff_ERROR" : eff_data["BOIL_Eff_ERROR"].iloc[0],
            }        
    else:
        print("Error: DataFrame 'eff_data' is empty for {}.".format(runNum))

            
    return [effDict,efficiency_errDict]

def calculate_efficiency(runNum,efficiency_table,foutcsv):

    effDict = get_efficiencies(runNum,efficiency_table,foutcsv)[0] # Efficiency dictionary

    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!

    # Calculate run by run total efficiency
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))
    
    data_to_csv(foutcsv, "Run Number", runNum, runNum)
    data_to_csv(foutcsv, "Total Efficiency", tot_efficiency, runNum)
    
    return tot_efficiency

def calculate_efficiency_err(runNum,efficiency_table,foutcsv):

    runbyrun_efficiencies = get_efficiencies(runNum,efficiency_table,foutcsv)
    effDict = runbyrun_efficiencies[0] # Efficiency dictionary
    efficiency_errDict = runbyrun_efficiencies[1] # Efficiency errors dictionary
    
    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))

    # Calculate run by run total efficiency error
    # Error propagation by addition in quadrature
    d_eff = np.sqrt(sum((float(efferr)/float(eff))**2 for efferr,eff in zip(efficiency_errDict.values(),effDict.values())))
    # Error propagation by addition in quadrature
    tot_efficiency_err = np.sqrt((tot_efficiency**2)*(d_eff**2))

    data_to_csv(foutcsv, "Total Efficiency Error", tot_efficiency_err, runNum)
    
    return tot_efficiency_err

def calculate_eff_charge(runNum,efficiency_table,foutcsv):

    effDict = get_efficiencies(runNum,efficiency_table,foutcsv)[0] # Efficiency dictionary

    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!

    # Calculate run by run total efficiency
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))

    charge  = get_bcm(runNum,efficiency_table,foutcsv)

    eff_charge = tot_efficiency*charge

    data_to_csv(foutcsv, "Charge", charge, runNum)
    data_to_csv(foutcsv, "Effective Charge", eff_charge, runNum)
    
    return eff_charge

def calculate_eff_charge_err(runNum,efficiency_table,foutcsv):

    runbyrun_efficiencies = get_efficiencies(runNum,efficiency_table,foutcsv)
    effDict = runbyrun_efficiencies[0] # Efficiency dictionary
    efficiency_errDict = runbyrun_efficiencies[1] # Efficiency errors dictionary
    
    charge  = get_bcm(runNum,efficiency_table,foutcsv)
    
    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))

    eff_charge = tot_efficiency*charge
    
    # Calculate run by run total efficiency error
    # Error propagation by addition in quadrature
    d_eff = np.sqrt(sum((float(efferr)/float(eff))**2 for efferr,eff in zip(efficiency_errDict.values(),effDict.values())))
    # Dave Mack's write up on the charge error: dI/I(%) = 0.05 + 0.15/I (in uA) -> dQ/Q(%) = 0.05 + 0.15/Q
    # https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1034
    d_charge = 0.05 + 0.00015/eff_charge # in percent (0.15 uC converted to 0.00015 mC)
    # Error propagation by addition in quadrature (units of mC)
    eff_charge_err = np.sqrt((eff_charge**2)*(d_eff**2+d_charge**2))
    
    data_to_csv(foutcsv, "Effective Charge Error", eff_charge_err, runNum)

    for (key1, value1), (key2, value2) in zip(effDict.items(), efficiency_errDict.items()):    
        data_to_csv(foutcsv, key1, value1, runNum)
        data_to_csv(foutcsv, key2, value2, runNum)
    
    return eff_charge_err
