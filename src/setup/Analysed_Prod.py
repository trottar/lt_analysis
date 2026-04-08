#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-26 16:29:37 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

###################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
from pathlib import Path

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with -  RunNumber ParticleType ROOTPrefix\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number, particle type, and max number of events
runNum = sys.argv[1]
ParticleType = sys.argv[2]
ROOTPrefix = sys.argv[3]
MaxEvent = "-1"

##############################################################################################################################################
'''
Define and set up cuts
'''

cut_f = '/DB/CUTS/run_type/coin_prod.cuts'

# defining Cuts
if ParticleType == "kaon":
    cuts = ["coin_ek_cut_all_noRF","coin_ek_cut_prompt_noRF","coin_ek_cut_rand_noRF","coin_ek_cut_all_RF","coin_ek_cut_prompt_RF","coin_ek_cut_rand_RF"]
if ParticleType == "pion":
    cuts = ["coin_epi_cut_all_noRF","coin_epi_cut_prompt_noRF","coin_epi_cut_rand_noRF","coin_epi_cut_all_RF","coin_epi_cut_prompt_RF","coin_epi_cut_rand_RF"]
if ParticleType == "proton":
    cuts = ["coin_ep_cut_all_noRF","coin_ep_cut_prompt_noRF","coin_ep_cut_rand_noRF","coin_ep_cut_all_RF","coin_ep_cut_prompt_RF","coin_ep_cut_rand_RF"]    

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for progress bar
from ltsep import Misc
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))
from farm_env.ltsep_paths import create_ltsep_root
    
lt=create_ltsep_root(os.path.realpath(__file__),"Prod",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
SKIMPATH=lt.SKIMPATH
if str(SKIMPATH).endswith(f"{ANATYPE}LT"):
    SKIM_OUTPATH = str(SKIMPATH)
elif "None" in str(SKIMPATH):
    SKIM_OUTPATH = str(SKIMPATH).replace("None", f"{ANATYPE}LT")
else:
    SKIM_OUTPATH = os.path.join(str(SKIMPATH), f"{ANATYPE}LT")
os.makedirs(SKIM_OUTPATH, exist_ok=True)

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

#################################################################################################################################################################

def _build_structured_array_from_columns(column_data, headers, tree_name):
    if not column_data:
        print("WARNING: Skipping {} because the data is empty.".format(tree_name))
        return None

    if not headers or len(headers) != len(column_data):
        print("WARNING: Skipping {} due to issues with column names.".format(tree_name))
        return None

    numeric_columns = [np.asarray(column, dtype=np.float64) for column in column_data]
    if not numeric_columns or len(numeric_columns[0]) == 0:
        print("WARNING: Skipping {} because the data is empty.".format(tree_name))
        return None

    n_rows = len(numeric_columns[0])
    if any(len(column) != n_rows for column in numeric_columns):
        print("WARNING: Skipping {} due to mismatched column lengths.".format(tree_name))
        return None

    return np.core.records.fromarrays(numeric_columns, names=headers)


def _coerce_cut_indices(index_values, n_rows):
    raw_indices = np.asarray(index_values)
    if raw_indices.size == 0:
        return np.asarray([], dtype=np.int64)

    if np.issubdtype(raw_indices.dtype, np.integer):
        cut_indices = raw_indices.astype(np.int64, copy=False)
    elif np.issubdtype(raw_indices.dtype, np.floating):
        if not np.all(np.isfinite(raw_indices)):
            return None
        rounded_indices = np.rint(raw_indices)
        if not np.allclose(raw_indices, rounded_indices):
            return None
        cut_indices = rounded_indices.astype(np.int64, copy=False)
    else:
        return None

    if np.any(cut_indices < 0) or np.any(cut_indices >= n_rows):
        return None

    return cut_indices


def _apply_named_cuts(column_data, cut_names):
    numeric_columns = [np.asarray(column) for column in column_data]
    if not numeric_columns:
        return {cut_name: [] for cut_name in cut_names}

    n_rows = len(numeric_columns[0])
    if any(len(column) != n_rows for column in numeric_columns):
        raise ValueError("All columns must have the same number of rows before cuts are applied.")

    event_index = np.arange(n_rows, dtype=np.int64)
    reference_column = numeric_columns[0]
    cut_results = {}

    for i, cut_name in enumerate(cut_names):
        print("\nApplying {}...".format(cut_name), flush=True)
        fast_cut_indices = _coerce_cut_indices(c.add_cut(event_index, cut_name), n_rows)

        if fast_cut_indices is not None:
            reference_cut = np.asarray(c.add_cut(reference_column, cut_name))
            if np.array_equal(reference_column[fast_cut_indices], reference_cut):
                cut_results[cut_name] = [column[fast_cut_indices] for column in numeric_columns]
                Misc.progressBar(i, len(cut_names)-1, bar_length=25)
                sys.stdout.flush()
                continue

        cut_results[cut_name] = [c.add_cut(column, cut_name) for column in numeric_columns]
        Misc.progressBar(i, len(cut_names)-1, bar_length=25)
        sys.stdout.flush()

    return cut_results


#################################################################################################################################################################

def coin_kaon():

    # Define the array of arrays containing the relevant HMS and SHMS info                              

    NoCut_COIN_Kaon = [tree["H_gtr_yp"],tree["H_gtr_xp"],tree["H_dc_yp_fp"],tree["H_dc_xp_fp"],tree["H_dc_y_fp"],tree["H_dc_x_fp"],tree["P_gtr_yp"],tree["P_gtr_xp"],tree["P_dc_yp_fp"],tree["P_dc_xp_fp"],tree["P_dc_y_fp"],tree["P_dc_x_fp"],tree["P_dc_InsideDipoleExit"],tree["H_dc_InsideDipoleExit"],tree["H_gtr_beta"],  tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_eKCoinTime_ROC1"], tree["P_gtr_beta"],  tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMK"], tree["H_RF_Dist"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"],tree["ph_recoil"],tree["th_q"],tree["th_recoil"], tree["MandelT"], tree["emiss"],tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"],tree["Erecoil"],tree["emiss_nuc"],tree["Mrecoil"],tree["raster_x"],tree["raster_y"],tree["raster_z"]]

    Uncut_COIN_Kaon = list(NoCut_COIN_Kaon)

    # Create array of arrays of pions after cuts, all events, prompt and random          

    kaon_cut_results = _apply_named_cuts(NoCut_COIN_Kaon, [
        "coin_ek_cut_all_noRF",
        "coin_ek_cut_prompt_noRF",
        "coin_ek_cut_rand_noRF",
        "coin_ek_cut_all_RF",
        "coin_ek_cut_prompt_RF",
        "coin_ek_cut_rand_RF",
    ])

    Cut_Kaon_Events_all_noRF = kaon_cut_results["coin_ek_cut_all_noRF"]
    Cut_Kaon_Events_prompt_noRF = kaon_cut_results["coin_ek_cut_prompt_noRF"]
    Cut_Kaon_Events_rand_noRF = kaon_cut_results["coin_ek_cut_rand_noRF"]
    Cut_Kaon_Events_all_RF = kaon_cut_results["coin_ek_cut_all_RF"]
    Cut_Kaon_Events_prompt_RF = kaon_cut_results["coin_ek_cut_prompt_RF"]
    Cut_Kaon_Events_rand_RF = kaon_cut_results["coin_ek_cut_rand_RF"]
    
    COIN_Kaon = {
        "Uncut_Kaon_Events" : Uncut_COIN_Kaon,
        "Cut_Kaon_Events_all_noRF" : Cut_Kaon_Events_all_noRF,
        "Cut_Kaon_Events_prompt_noRF" : Cut_Kaon_Events_prompt_noRF,
        "Cut_Kaon_Events_rand_noRF" : Cut_Kaon_Events_rand_noRF,
        "Cut_Kaon_Events_all_RF" : Cut_Kaon_Events_all_RF,
        "Cut_Kaon_Events_prompt_RF" : Cut_Kaon_Events_prompt_RF,
        "Cut_Kaon_Events_rand_RF" : Cut_Kaon_Events_rand_RF,
        }

    return COIN_Kaon

##################################################################################################################################################################

def coin_pion():

    # Define the array of arrays containing the relevant HMS and SHMS info                              

    NoCut_COIN_Pion = [tree["H_gtr_yp"],tree["H_gtr_xp"],tree["H_dc_yp_fp"],tree["H_dc_xp_fp"],tree["H_dc_y_fp"],tree["H_dc_x_fp"],tree["P_gtr_yp"],tree["P_gtr_xp"],tree["P_dc_yp_fp"],tree["P_dc_xp_fp"],tree["P_dc_y_fp"],tree["P_dc_x_fp"],tree["P_dc_InsideDipoleExit"],tree["H_dc_InsideDipoleExit"],tree["H_gtr_beta"],  tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_ePiCoinTime_ROC1"], tree["P_gtr_beta"],  tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMK"], tree["H_RF_Dist"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"],tree["ph_recoil"],tree["th_q"],tree["th_recoil"], tree["MandelT"], tree["emiss"],tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"],tree["Erecoil"],tree["emiss_nuc"],tree["Mrecoil"],tree["raster_x"],tree["raster_y"],tree["raster_z"]]

    Uncut_COIN_Pion = list(NoCut_COIN_Pion)

    # Create array of arrays of pions after cuts, all events, prompt and random          

    pion_cut_results = _apply_named_cuts(NoCut_COIN_Pion, [
        "coin_epi_cut_all_noRF",
        "coin_epi_cut_prompt_noRF",
        "coin_epi_cut_rand_noRF",
        "coin_epi_cut_all_RF",
        "coin_epi_cut_prompt_RF",
        "coin_epi_cut_rand_RF",
    ])

    Cut_Pion_Events_all_noRF = pion_cut_results["coin_epi_cut_all_noRF"]
    Cut_Pion_Events_prompt_noRF = pion_cut_results["coin_epi_cut_prompt_noRF"]
    Cut_Pion_Events_rand_noRF = pion_cut_results["coin_epi_cut_rand_noRF"]
    Cut_Pion_Events_all_RF = pion_cut_results["coin_epi_cut_all_RF"]
    Cut_Pion_Events_prompt_RF = pion_cut_results["coin_epi_cut_prompt_RF"]
    Cut_Pion_Events_rand_RF = pion_cut_results["coin_epi_cut_rand_RF"]
    
    COIN_Pion = {
        "Uncut_Pion_Events" : Uncut_COIN_Pion,
        "Cut_Pion_Events_all_noRF" : Cut_Pion_Events_all_noRF,
        "Cut_Pion_Events_prompt_noRF" : Cut_Pion_Events_prompt_noRF,
        "Cut_Pion_Events_rand_noRF" : Cut_Pion_Events_rand_noRF,
        "Cut_Pion_Events_all_RF" : Cut_Pion_Events_all_RF,
        "Cut_Pion_Events_prompt_RF" : Cut_Pion_Events_prompt_RF,
        "Cut_Pion_Events_rand_RF" : Cut_Pion_Events_rand_RF,
        }

    return COIN_Pion

##################################################################################################################################################################

def coin_proton():

    # Define the array of arrays containing the relevant HMS and SHMS info                              

    NoCut_COIN_Proton = [tree["H_gtr_yp"],tree["H_gtr_xp"],tree["H_dc_yp_fp"],tree["H_dc_xp_fp"],tree["H_dc_y_fp"],tree["H_dc_x_fp"],tree["P_gtr_yp"],tree["P_gtr_xp"],tree["P_dc_yp_fp"],tree["P_dc_xp_fp"],tree["P_dc_y_fp"],tree["P_dc_x_fp"],tree["P_dc_InsideDipoleExit"],tree["H_dc_InsideDipoleExit"],tree["H_gtr_beta"],  tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"],  tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMK"], tree["H_RF_Dist"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"],tree["ph_recoil"],tree["th_q"],tree["th_recoil"], tree["MandelT"], tree["emiss"],tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"],tree["Erecoil"],tree["emiss_nuc"],tree["Mrecoil"],tree["raster_x"],tree["raster_y"],tree["raster_z"]]

    Uncut_COIN_Proton = list(NoCut_COIN_Proton)

    # Create array of arrays of protons after cuts, all events, prompt and random          

    proton_cut_results = _apply_named_cuts(NoCut_COIN_Proton, [
        "coin_ep_cut_all_noRF",
        "coin_ep_cut_prompt_noRF",
        "coin_ep_cut_rand_noRF",
        "coin_ep_cut_all_RF",
        "coin_ep_cut_prompt_RF",
        "coin_ep_cut_rand_RF",
    ])

    Cut_Proton_Events_all_noRF = proton_cut_results["coin_ep_cut_all_noRF"]
    Cut_Proton_Events_prompt_noRF = proton_cut_results["coin_ep_cut_prompt_noRF"]
    Cut_Proton_Events_rand_noRF = proton_cut_results["coin_ep_cut_rand_noRF"]
    Cut_Proton_Events_all_RF = proton_cut_results["coin_ep_cut_all_RF"]
    Cut_Proton_Events_prompt_RF = proton_cut_results["coin_ep_cut_prompt_RF"]
    Cut_Proton_Events_rand_RF = proton_cut_results["coin_ep_cut_rand_RF"]
    
    COIN_Proton = {
        "Uncut_Proton_Events" : Uncut_COIN_Proton,
        "Cut_Proton_Events_all_noRF" : Cut_Proton_Events_all_noRF,
        "Cut_Proton_Events_prompt_noRF" : Cut_Proton_Events_prompt_noRF,
        "Cut_Proton_Events_rand_noRF" : Cut_Proton_Events_rand_noRF,
        "Cut_Proton_Events_all_RF" : Cut_Proton_Events_all_RF,
        "Cut_Proton_Events_prompt_RF" : Cut_Proton_Events_prompt_RF,
        "Cut_Proton_Events_rand_RF" : Cut_Proton_Events_rand_RF,
        }

    return COIN_Proton

##################################################################################################################################################################

def main():
        
    print("Applying cuts for {}...".format(ParticleType))
    if ParticleType == "kaon":
        COIN_Data = coin_kaon()
    if ParticleType == "pion":
        COIN_Data = coin_pion()
    if ParticleType == "proton":
        COIN_Data = coin_proton()        

    print("\nCuts complete for {}. Starting ROOT output...".format(ParticleType), flush=True)

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)
    # RLT: I did what Stephen said not to so I can match naming scheme of simc root output

    COIN_Data_Header = ["hsyptar","hsxptar","hsypfp","hsxpfp","hsyfp","hsxfp","ssyptar","ssxptar","ssypfp","ssxpfp","ssyfp","ssxfp","P_dc_InsideDipoleExit","H_dc_InsideDipoleExit","H_gtr_beta", "hsdelta", "H_gtr_p", "H_hod_goodscinhit", "H_hod_goodstarttime", "H_cal_etotnorm", "H_cal_etottracknorm", "H_cer_npeSum", "CTime_ROC1", "P_gtr_beta", "P_gtr_p", "ssdelta", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_cal_etotnorm", "P_cal_etottracknorm", "P_aero_npeSum", "P_aero_xAtAero", "P_aero_yAtAero", "P_hgcer_npeSum", "P_hgcer_xAtCer", "P_hgcer_yAtCer", "MM", "H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon", "ph_q", "ph_recoil", "th_q", "th_recoil", "MandelT", "emiss", "pmiss", "pmx", "pmy", "pmz","Erecoil","emiss_nuc","Mrecoil","raster_x","raster_y","raster_z"]
    
    # Need to create a dict for all the branches we grab                                                
    data = {}
    data.update(COIN_Data)
    data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data

    term_search = ParticleType.capitalize()
    output_keys = [key for key in data_keys if term_search in key]
    
    out_f_file = "%s/%s_%s_%s_Raw_Data.root" % (SKIM_OUTPATH, ParticleType, runNum, MaxEvent)
        
    print("\n\nSaving data to new root files...", flush=True)
    output_index = 0
    for i in range (0, len(data_keys)):
        if(term_search in data_keys[i]):
            DFHeader=list(COIN_Data_Header)
            output_index += 1
        else:
            continue
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!                                                                      
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header
        # Old way changed 1/22/2023
        #if (i == 0):
        #    pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root(out_f_file, key ="%s" % data_keys[i])
        #elif (i != 0):
        #    pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root(out_f_file, key ="%s" % data_keys[i], mode ='a')
        
        print("\n[{}/{}] Preparing {}...".format(output_index, len(output_keys), data_keys[i]), flush=True)
        structured_array = _build_structured_array_from_columns(data.get(data_keys[i]), DFHeader, data_keys[i])
        if structured_array is None:
            Misc.progressBar(output_index-1, len(output_keys)-1, bar_length=25)
            sys.stdout.flush()
            continue

        # Save the structured array to ROOT file
        print("[{}/{}] Writing {} to {}...".format(output_index, len(output_keys), data_keys[i], out_f_file), flush=True)
        rnp.array2root(structured_array, out_f_file, mode='recreate' if i == 0 else 'update', treename=data_keys[i])
        
        Misc.progressBar(output_index-1, len(output_keys)-1, bar_length=25)
        sys.stdout.flush()

if __name__ == '__main__':
    main()
print ("Processing Complete")
