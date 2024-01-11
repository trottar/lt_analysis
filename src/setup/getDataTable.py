#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-11 16:08:25 trottar"
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

################################################################################################################################################
# Grab bcm value

def get_bcm(runNum,efficiency_table):
    
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

def get_ebeam(runNum,efficiency_table):
    
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

def get_pTheta(runNum,efficiency_table):
    
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

def get_efficiencies(runNum,efficiency_table):

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
                # HMS Cal
                "HMS_Cal_ALL_Elec_Eff" : eff_data["HMS_Cal_ALL_Elec_Eff"].iloc[0],            
                # HMS Cer
                "HMS_Cer_ALL_Elec_Eff" : eff_data["HMS_Cer_ALL_Elec_Eff"].iloc[0],
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF" : eff_data["HMS_Elec_ALL_TRACK_EFF"].iloc[0],
                "SHMS_Prot_ALL_TRACK_EFF" : eff_data["SHMS_Prot_ALL_TRACK_EFF"].iloc[0],
                # Boiling Correction
                "BOIL_Eff" : eff_data["BOIL_Eff"].iloc[0],
            }
            # Define dictionary of efficiency uncertainty values
            effErrorDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time_ERROR" : eff_data["Non_Scaler_EDTM_Live_Time_ERROR"].iloc[0],
                # Hodo
                #"HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                #"SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                "HMS_Hodo_3_of_4_EFF_ERROR" : 0.00001*eff_data["HMS_Hodo_3_of_4_EFF"].iloc[0], # FIX
                "SHMS_Hodo_3_of_4_EFF_ERROR" : 0.00001*eff_data["SHMS_Hodo_3_of_4_EFF"].iloc[0], # FIX
                # HMS Cal
                "HMS_Cal_ALL_Elec_Eff_ERROR" : eff_data["HMS_Cal_ALL_Elec_Eff_ERROR"].iloc[0],            
                # HMS Cer
                "HMS_Cer_ALL_Elec_Eff_ERROR" : eff_data["HMS_Cer_ALL_Elec_Eff_ERROR"].iloc[0],
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
                # HMS Cal
                "HMS_Cal_ALL_Elec_Eff" : eff_data["HMS_Cal_ALL_Elec_Eff"].iloc[0],            
                # HMS Cer
                "HMS_Cer_ALL_Elec_Eff" : eff_data["HMS_Cer_ALL_Elec_Eff"].iloc[0],
                # SHMS Aero
                "SHMS_Aero_ALL_Pion_Eff" : eff_data["SHMS_Aero_ALL_Pion_Eff"].iloc[0],
                # Tracking
                "HMS_Elec_ALL_TRACK_EFF" : eff_data["HMS_Elec_ALL_TRACK_EFF"].iloc[0],
                "SHMS_Pion_ALL_TRACK_EFF" : eff_data["SHMS_Pion_ALL_TRACK_EFF"].iloc[0],
                # Boiling Correction
                "BOIL_Eff" : eff_data["BOIL_Eff"].iloc[0],
            }

            # Define dictionary of efficiency uncertainty values
            effErrorDict ={
                # EDTM
                "Non_Scaler_EDTM_Live_Time_ERROR" : eff_data["Non_Scaler_EDTM_Live_Time_ERROR"].iloc[0],
                # Hodo
                #"HMS_Hodo_3_of_4_EFF_ERROR" : eff_data["HMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                #"SHMS_Hodo_3_of_4_EFF_ERROR" : eff_data["SHMS_Hodo_3_of_4_EFF_ERROR"].iloc[0],
                "HMS_Hodo_3_of_4_EFF_ERROR" : 0.00001*eff_data["HMS_Hodo_3_of_4_EFF"].iloc[0], # FIX
                "SHMS_Hodo_3_of_4_EFF_ERROR" : 0.00001*eff_data["SHMS_Hodo_3_of_4_EFF"].iloc[0], # FIX
                # HMS Cal
                "HMS_Cal_ALL_Elec_Eff_ERROR" : eff_data["HMS_Cal_ALL_Elec_Eff_ERROR"].iloc[0],            
                # HMS Cer
                "HMS_Cer_ALL_Elec_Eff_ERROR" : eff_data["HMS_Cer_ALL_Elec_Eff_ERROR"].iloc[0],
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

            
    return [effDict,effErrorDict]

def calculate_efficiency(runNum,efficiency_table):

    effDict = get_efficiencies(runNum,efficiency_table)[0] # Efficiency dictionary

    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))
    return tot_efficiency

def calculate_effError(runNum,efficiency_table):

    effDict = get_efficiencies(runNum,efficiency_table)[0] # Efficiency dictionary
    effErrorDict = get_efficiencies(runNum,efficiency_table)[1] # Efficiency errors dictionary

    # Calculate total efficiency. The reduce function pretty much iterates on
    # its arguments which in this case is a lambda function. This lambda function
    # takes x,y from the list (ie the list of efficiencies) and multiplies them.
    # This is all pythonic mumbo-jumbo for doing the product of everything in the
    # list. Enjoy!
    tot_efficiency = reduce(lambda x, y: x*y, list(effDict.values()))

    # Calculate total efficiency error. 
    tot_effError = tot_efficiency*np.sqrt(sum((float(efferr)/float(eff))*(float(efferr)/float(eff)) for efferr,eff in zip(effErrorDict.values(),effDict.values()))) # Error propagation for multiplication, using quadratic correction
    return tot_effError
