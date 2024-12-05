#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-05 05:46:39 trottar"
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
from ROOT import TH1D, TCutG, TFile
import sys, os, math

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
from utility import open_root_file

################################################################################################################################################

def particle_subtraction_cuts(subDict, inpDict, SubtractedParticle, hgcer_cutg=None):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    MM_offset_DUMMY = subDict["MM_offset_DUMMY"]
    MM_offset_RAND = subDict["MM_offset_RAND"]
    MM_offset_DUMMY_RAND = subDict["MM_offset_DUMMY_RAND"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################
    
    H_hsdelta_DATA = subDict["H_hsdelta_SUB_DATA"]
    H_hsxptar_DATA = subDict["H_hsxptar_SUB_DATA"]
    H_hsyptar_DATA = subDict["H_hsyptar_SUB_DATA"]
    H_ssxfp_DATA = subDict["H_ssxfp_SUB_DATA"]
    H_ssyfp_DATA = subDict["H_ssyfp_SUB_DATA"]
    H_ssxpfp_DATA = subDict["H_ssxpfp_SUB_DATA"]
    H_ssypfp_DATA = subDict["H_ssypfp_SUB_DATA"]
    H_hsxfp_DATA = subDict["H_hsxfp_SUB_DATA"]
    H_hsyfp_DATA = subDict["H_hsyfp_SUB_DATA"]
    H_hsxpfp_DATA = subDict["H_hsxpfp_SUB_DATA"]
    H_hsypfp_DATA = subDict["H_hsypfp_SUB_DATA"]
    H_ssdelta_DATA = subDict["H_ssdelta_SUB_DATA"]
    H_ssxptar_DATA = subDict["H_ssxptar_SUB_DATA"]
    H_ssyptar_DATA = subDict["H_ssyptar_SUB_DATA"]
    H_q_DATA = subDict["H_q_SUB_DATA"]
    H_Q2_DATA = subDict["H_Q2_SUB_DATA"]
    H_W_DATA = subDict["H_W_SUB_DATA"]
    H_t_DATA = subDict["H_t_SUB_DATA"]
    H_epsilon_DATA = subDict["H_epsilon_SUB_DATA"]
    H_MM_DATA = subDict["H_MM_SUB_DATA"]
    H_MM_nosub_DATA = subDict["H_MM_nosub_SUB_DATA"]
    H_th_DATA = subDict["H_th_SUB_DATA"]
    H_ph_DATA = subDict["H_ph_SUB_DATA"]
    H_ph_q_DATA = subDict["H_ph_q_SUB_DATA"]
    H_th_q_DATA = subDict["H_th_q_SUB_DATA"]
    H_ph_recoil_DATA = subDict["H_ph_recoil_SUB_DATA"]
    H_th_recoil_DATA = subDict["H_th_recoil_SUB_DATA"]
    H_pmiss_DATA = subDict["H_pmiss_SUB_DATA"]
    H_emiss_DATA = subDict["H_emiss_SUB_DATA"]
    H_pmx_DATA = subDict["H_pmx_SUB_DATA"]
    H_pmy_DATA = subDict["H_pmy_SUB_DATA"]
    H_pmz_DATA = subDict["H_pmz_SUB_DATA"]
    H_ct_DATA = subDict["H_ct_SUB_DATA"]
    H_cal_etottracknorm_DATA = subDict["H_cal_etottracknorm_SUB_DATA"]
    H_cer_npeSum_DATA = subDict["H_cer_npeSum_SUB_DATA"]
    P_cal_etottracknorm_DATA = subDict["P_cal_etottracknorm_SUB_DATA"]
    P_hgcer_npeSum_DATA = subDict["P_hgcer_npeSum_SUB_DATA"]
    P_aero_npeSum_DATA = subDict["P_aero_npeSum_SUB_DATA"]

    H_hsdelta_DUMMY = subDict["H_hsdelta_SUB_DUMMY"]
    H_hsxptar_DUMMY = subDict["H_hsxptar_SUB_DUMMY"]
    H_hsyptar_DUMMY = subDict["H_hsyptar_SUB_DUMMY"]
    H_ssxfp_DUMMY = subDict["H_ssxfp_SUB_DUMMY"]
    H_ssyfp_DUMMY = subDict["H_ssyfp_SUB_DUMMY"]
    H_ssxpfp_DUMMY = subDict["H_ssxpfp_SUB_DUMMY"]
    H_ssypfp_DUMMY = subDict["H_ssypfp_SUB_DUMMY"]
    H_hsxfp_DUMMY = subDict["H_hsxfp_SUB_DUMMY"]
    H_hsyfp_DUMMY = subDict["H_hsyfp_SUB_DUMMY"]
    H_hsxpfp_DUMMY = subDict["H_hsxpfp_SUB_DUMMY"]
    H_hsypfp_DUMMY = subDict["H_hsypfp_SUB_DUMMY"]
    H_ssdelta_DUMMY = subDict["H_ssdelta_SUB_DUMMY"]
    H_ssxptar_DUMMY = subDict["H_ssxptar_SUB_DUMMY"]
    H_ssyptar_DUMMY = subDict["H_ssyptar_SUB_DUMMY"]
    H_q_DUMMY = subDict["H_q_SUB_DUMMY"]
    H_Q2_DUMMY = subDict["H_Q2_SUB_DUMMY"]
    H_W_DUMMY = subDict["H_W_SUB_DUMMY"]
    H_t_DUMMY = subDict["H_t_SUB_DUMMY"]
    H_epsilon_DUMMY = subDict["H_epsilon_SUB_DUMMY"]
    H_MM_DUMMY = subDict["H_MM_SUB_DUMMY"]
    H_MM_nosub_DUMMY = subDict["H_MM_nosub_SUB_DUMMY"]    
    H_th_DUMMY = subDict["H_th_SUB_DUMMY"]
    H_ph_DUMMY = subDict["H_ph_SUB_DUMMY"]
    H_ph_q_DUMMY = subDict["H_ph_q_SUB_DUMMY"]
    H_th_q_DUMMY = subDict["H_th_q_SUB_DUMMY"]
    H_ph_recoil_DUMMY = subDict["H_ph_recoil_SUB_DUMMY"]
    H_th_recoil_DUMMY = subDict["H_th_recoil_SUB_DUMMY"]
    H_pmiss_DUMMY = subDict["H_pmiss_SUB_DUMMY"]
    H_emiss_DUMMY = subDict["H_emiss_SUB_DUMMY"]
    H_pmx_DUMMY = subDict["H_pmx_SUB_DUMMY"]
    H_pmy_DUMMY = subDict["H_pmy_SUB_DUMMY"]
    H_pmz_DUMMY = subDict["H_pmz_SUB_DUMMY"]
    H_ct_DUMMY = subDict["H_ct_SUB_DUMMY"]
    H_cal_etottracknorm_DUMMY = subDict["H_cal_etottracknorm_SUB_DUMMY"]
    H_cer_npeSum_DUMMY = subDict["H_cer_npeSum_SUB_DUMMY"]
    P_cal_etottracknorm_DUMMY = subDict["P_cal_etottracknorm_SUB_DUMMY"]
    P_hgcer_npeSum_DUMMY = subDict["P_hgcer_npeSum_SUB_DUMMY"]
    P_aero_npeSum_DUMMY = subDict["P_aero_npeSum_SUB_DUMMY"]

    H_hsdelta_RAND = subDict["H_hsdelta_SUB_RAND"]
    H_hsxptar_RAND = subDict["H_hsxptar_SUB_RAND"]
    H_hsyptar_RAND = subDict["H_hsyptar_SUB_RAND"]
    H_ssxfp_RAND = subDict["H_ssxfp_SUB_RAND"]
    H_ssyfp_RAND = subDict["H_ssyfp_SUB_RAND"]
    H_ssxpfp_RAND = subDict["H_ssxpfp_SUB_RAND"]
    H_ssypfp_RAND = subDict["H_ssypfp_SUB_RAND"]
    H_hsxfp_RAND = subDict["H_hsxfp_SUB_RAND"]
    H_hsyfp_RAND = subDict["H_hsyfp_SUB_RAND"]
    H_hsxpfp_RAND = subDict["H_hsxpfp_SUB_RAND"]
    H_hsypfp_RAND = subDict["H_hsypfp_SUB_RAND"]
    H_ssdelta_RAND = subDict["H_ssdelta_SUB_RAND"]
    H_ssxptar_RAND = subDict["H_ssxptar_SUB_RAND"]
    H_ssyptar_RAND = subDict["H_ssyptar_SUB_RAND"]
    H_q_RAND = subDict["H_q_SUB_RAND"]
    H_Q2_RAND = subDict["H_Q2_SUB_RAND"]
    H_W_RAND = subDict["H_W_SUB_RAND"]
    H_t_RAND = subDict["H_t_SUB_RAND"]
    H_epsilon_RAND = subDict["H_epsilon_SUB_RAND"]
    H_MM_RAND = subDict["H_MM_SUB_RAND"]
    H_MM_nosub_RAND = subDict["H_MM_nosub_SUB_RAND"]    
    H_th_RAND = subDict["H_th_SUB_RAND"]
    H_ph_RAND = subDict["H_ph_SUB_RAND"]
    H_ph_q_RAND = subDict["H_ph_q_SUB_RAND"]
    H_th_q_RAND = subDict["H_th_q_SUB_RAND"]
    H_ph_recoil_RAND = subDict["H_ph_recoil_SUB_RAND"]
    H_th_recoil_RAND = subDict["H_th_recoil_SUB_RAND"]
    H_pmiss_RAND = subDict["H_pmiss_SUB_RAND"]
    H_emiss_RAND = subDict["H_emiss_SUB_RAND"]
    H_pmx_RAND = subDict["H_pmx_SUB_RAND"]
    H_pmy_RAND = subDict["H_pmy_SUB_RAND"]
    H_pmz_RAND = subDict["H_pmz_SUB_RAND"]
    H_ct_RAND = subDict["H_ct_SUB_RAND"]
    H_cal_etottracknorm_RAND = subDict["H_cal_etottracknorm_SUB_RAND"]
    H_cer_npeSum_RAND = subDict["H_cer_npeSum_SUB_RAND"]
    P_cal_etottracknorm_RAND = subDict["P_cal_etottracknorm_SUB_RAND"]
    P_hgcer_npeSum_RAND = subDict["P_hgcer_npeSum_SUB_RAND"]
    P_aero_npeSum_RAND = subDict["P_aero_npeSum_SUB_RAND"]

    H_hsdelta_DUMMY_RAND = subDict["H_hsdelta_SUB_DUMMY_RAND"]
    H_hsxptar_DUMMY_RAND = subDict["H_hsxptar_SUB_DUMMY_RAND"]
    H_hsyptar_DUMMY_RAND = subDict["H_hsyptar_SUB_DUMMY_RAND"]
    H_ssxfp_DUMMY_RAND = subDict["H_ssxfp_SUB_DUMMY_RAND"]
    H_ssyfp_DUMMY_RAND = subDict["H_ssyfp_SUB_DUMMY_RAND"]
    H_ssxpfp_DUMMY_RAND = subDict["H_ssxpfp_SUB_DUMMY_RAND"]
    H_ssypfp_DUMMY_RAND = subDict["H_ssypfp_SUB_DUMMY_RAND"]
    H_hsxfp_DUMMY_RAND = subDict["H_hsxfp_SUB_DUMMY_RAND"]
    H_hsyfp_DUMMY_RAND = subDict["H_hsyfp_SUB_DUMMY_RAND"]
    H_hsxpfp_DUMMY_RAND = subDict["H_hsxpfp_SUB_DUMMY_RAND"]
    H_hsypfp_DUMMY_RAND = subDict["H_hsypfp_SUB_DUMMY_RAND"]
    H_ssdelta_DUMMY_RAND = subDict["H_ssdelta_SUB_DUMMY_RAND"]
    H_ssxptar_DUMMY_RAND = subDict["H_ssxptar_SUB_DUMMY_RAND"]
    H_ssyptar_DUMMY_RAND = subDict["H_ssyptar_SUB_DUMMY_RAND"]
    H_q_DUMMY_RAND = subDict["H_q_SUB_DUMMY_RAND"]
    H_Q2_DUMMY_RAND = subDict["H_Q2_SUB_DUMMY_RAND"]
    H_W_DUMMY_RAND = subDict["H_W_SUB_DUMMY_RAND"]
    H_t_DUMMY_RAND = subDict["H_t_SUB_DUMMY_RAND"]
    H_epsilon_DUMMY_RAND = subDict["H_epsilon_SUB_DUMMY_RAND"]
    H_MM_DUMMY_RAND = subDict["H_MM_SUB_DUMMY_RAND"]
    H_MM_nosub_DUMMY_RAND = subDict["H_MM_nosub_SUB_DUMMY_RAND"]    
    H_th_DUMMY_RAND = subDict["H_th_SUB_DUMMY_RAND"]
    H_ph_DUMMY_RAND = subDict["H_ph_SUB_DUMMY_RAND"]
    H_ph_q_DUMMY_RAND = subDict["H_ph_q_SUB_DUMMY_RAND"]
    H_th_q_DUMMY_RAND = subDict["H_th_q_SUB_DUMMY_RAND"]
    H_ph_recoil_DUMMY_RAND = subDict["H_ph_recoil_SUB_DUMMY_RAND"]
    H_th_recoil_DUMMY_RAND = subDict["H_th_recoil_SUB_DUMMY_RAND"]
    H_pmiss_DUMMY_RAND = subDict["H_pmiss_SUB_DUMMY_RAND"]
    H_emiss_DUMMY_RAND = subDict["H_emiss_SUB_DUMMY_RAND"]
    H_pmx_DUMMY_RAND = subDict["H_pmx_SUB_DUMMY_RAND"]
    H_pmy_DUMMY_RAND = subDict["H_pmy_SUB_DUMMY_RAND"]
    H_pmz_DUMMY_RAND = subDict["H_pmz_SUB_DUMMY_RAND"]
    H_ct_DUMMY_RAND = subDict["H_ct_SUB_DUMMY_RAND"]
    H_cal_etottracknorm_DUMMY_RAND = subDict["H_cal_etottracknorm_SUB_DUMMY_RAND"]
    H_cer_npeSum_DUMMY_RAND = subDict["H_cer_npeSum_SUB_DUMMY_RAND"]
    P_cal_etottracknorm_DUMMY_RAND = subDict["P_cal_etottracknorm_SUB_DUMMY_RAND"]
    P_hgcer_npeSum_DUMMY_RAND = subDict["P_hgcer_npeSum_SUB_DUMMY_RAND"]
    P_aero_npeSum_DUMMY_RAND = subDict["P_aero_npeSum_SUB_DUMMY_RAND"]

    MM_vs_CoinTime_DATA = subDict["MM_vs_CoinTime_SUB_DATA"]
    CoinTime_vs_beta_DATA = subDict["CoinTime_vs_beta_SUB_DATA"]
    MM_vs_beta_DATA = subDict["MM_vs_beta_SUB_DATA"]
    MM_vs_H_cer_DATA = subDict["MM_vs_H_cer_SUB_DATA"]
    MM_vs_H_cal_DATA = subDict["MM_vs_H_cal_SUB_DATA"]
    MM_vs_P_cal_DATA = subDict["MM_vs_P_cal_SUB_DATA"]
    MM_vs_P_hgcer_DATA = subDict["MM_vs_P_hgcer_SUB_DATA"]
    MM_vs_P_aero_DATA = subDict["MM_vs_P_aero_SUB_DATA"]
    phiq_vs_t_DATA = subDict["phiq_vs_t_SUB_DATA"]
    Q2_vs_W_DATA = subDict["Q2_vs_W_SUB_DATA"]
    Q2_vs_t_DATA = subDict["Q2_vs_t_SUB_DATA"]
    W_vs_t_DATA = subDict["W_vs_t_SUB_DATA"]
    EPS_vs_t_DATA = subDict["EPS_vs_t_SUB_DATA"]
    MM_vs_t_DATA = subDict["MM_vs_t_SUB_DATA"]
    P_hgcer_xAtCer_vs_yAtCer_DATA = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DATA"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DATA = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DATA"]
    P_hgcer_xAtCer_vs_MM_DATA = subDict["P_hgcer_xAtCer_vs_MM_SUB_DATA"]
    P_hgcer_nohole_xAtCer_vs_MM_DATA = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DATA"]
    P_hgcer_yAtCer_vs_MM_DATA = subDict["P_hgcer_yAtCer_vs_MM_SUB_DATA"]
    P_hgcer_nohole_yAtCer_vs_MM_DATA = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DATA"]

    MM_vs_CoinTime_DUMMY = subDict["MM_vs_CoinTime_SUB_DUMMY"]
    CoinTime_vs_beta_DUMMY = subDict["CoinTime_vs_beta_SUB_DUMMY"]
    MM_vs_beta_DUMMY = subDict["MM_vs_beta_SUB_DUMMY"]
    MM_vs_H_cer_DUMMY = subDict["MM_vs_H_cer_SUB_DUMMY"]
    MM_vs_H_cal_DUMMY = subDict["MM_vs_H_cal_SUB_DUMMY"]
    MM_vs_P_cal_DUMMY = subDict["MM_vs_P_cal_SUB_DUMMY"]
    MM_vs_P_hgcer_DUMMY = subDict["MM_vs_P_hgcer_SUB_DUMMY"]
    MM_vs_P_aero_DUMMY = subDict["MM_vs_P_aero_SUB_DUMMY"]
    phiq_vs_t_DUMMY = subDict["phiq_vs_t_SUB_DUMMY"]
    Q2_vs_W_DUMMY = subDict["Q2_vs_W_SUB_DUMMY"]
    Q2_vs_t_DUMMY = subDict["Q2_vs_t_SUB_DUMMY"]
    W_vs_t_DUMMY = subDict["W_vs_t_SUB_DUMMY"]
    EPS_vs_t_DUMMY = subDict["EPS_vs_t_SUB_DUMMY"]
    MM_vs_t_DUMMY = subDict["MM_vs_t_SUB_DUMMY"]
    P_hgcer_xAtCer_vs_yAtCer_DUMMY = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY"]
    P_hgcer_xAtCer_vs_MM_DUMMY = subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_nohole_xAtCer_vs_MM_DUMMY = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_yAtCer_vs_MM_DUMMY = subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY"]
    P_hgcer_nohole_yAtCer_vs_MM_DUMMY = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY"]

    MM_vs_CoinTime_RAND = subDict["MM_vs_CoinTime_SUB_RAND"]
    CoinTime_vs_beta_RAND = subDict["CoinTime_vs_beta_SUB_RAND"]
    MM_vs_beta_RAND = subDict["MM_vs_beta_SUB_RAND"]
    MM_vs_H_cer_RAND = subDict["MM_vs_H_cer_SUB_RAND"]
    MM_vs_H_cal_RAND = subDict["MM_vs_H_cal_SUB_RAND"]
    MM_vs_P_cal_RAND = subDict["MM_vs_P_cal_SUB_RAND"]
    MM_vs_P_hgcer_RAND = subDict["MM_vs_P_hgcer_SUB_RAND"]
    MM_vs_P_aero_RAND = subDict["MM_vs_P_aero_SUB_RAND"]
    phiq_vs_t_RAND = subDict["phiq_vs_t_SUB_RAND"]
    Q2_vs_W_RAND = subDict["Q2_vs_W_SUB_RAND"]
    Q2_vs_t_RAND = subDict["Q2_vs_t_SUB_RAND"]
    W_vs_t_RAND = subDict["W_vs_t_SUB_RAND"]
    EPS_vs_t_RAND = subDict["EPS_vs_t_SUB_RAND"]
    MM_vs_t_RAND = subDict["MM_vs_t_SUB_RAND"]
    P_hgcer_xAtCer_vs_yAtCer_RAND = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_RAND"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_RAND = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_RAND"]
    P_hgcer_xAtCer_vs_MM_RAND = subDict["P_hgcer_xAtCer_vs_MM_SUB_RAND"]
    P_hgcer_nohole_xAtCer_vs_MM_RAND = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_RAND"]
    P_hgcer_yAtCer_vs_MM_RAND = subDict["P_hgcer_yAtCer_vs_MM_SUB_RAND"]
    P_hgcer_nohole_yAtCer_vs_MM_RAND = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_RAND"]

    MM_vs_CoinTime_DUMMY_RAND = subDict["MM_vs_CoinTime_SUB_DUMMY_RAND"]
    CoinTime_vs_beta_DUMMY_RAND = subDict["CoinTime_vs_beta_SUB_DUMMY_RAND"]
    MM_vs_beta_DUMMY_RAND = subDict["MM_vs_beta_SUB_DUMMY_RAND"]
    MM_vs_H_cer_DUMMY_RAND = subDict["MM_vs_H_cer_SUB_DUMMY_RAND"]
    MM_vs_H_cal_DUMMY_RAND = subDict["MM_vs_H_cal_SUB_DUMMY_RAND"]
    MM_vs_P_cal_DUMMY_RAND = subDict["MM_vs_P_cal_SUB_DUMMY_RAND"]
    MM_vs_P_hgcer_DUMMY_RAND = subDict["MM_vs_P_hgcer_SUB_DUMMY_RAND"]
    MM_vs_P_aero_DUMMY_RAND = subDict["MM_vs_P_aero_SUB_DUMMY_RAND"]
    phiq_vs_t_DUMMY_RAND = subDict["phiq_vs_t_SUB_DUMMY_RAND"]
    Q2_vs_W_DUMMY_RAND = subDict["Q2_vs_W_SUB_DUMMY_RAND"]
    Q2_vs_t_DUMMY_RAND = subDict["Q2_vs_t_SUB_DUMMY_RAND"]
    W_vs_t_DUMMY_RAND = subDict["W_vs_t_SUB_DUMMY_RAND"]
    EPS_vs_t_DUMMY_RAND = subDict["EPS_vs_t_SUB_DUMMY_RAND"]
    MM_vs_t_DUMMY_RAND = subDict["MM_vs_t_SUB_DUMMY_RAND"]
    P_hgcer_xAtCer_vs_yAtCer_DUMMY_RAND = subDict["P_hgcer_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"]
    P_hgcer_nohole_xAtCer_vs_yAtCer_DUMMY_RAND = subDict["P_hgcer_nohole_xAtCer_vs_yAtCer_SUB_DUMMY_RAND"]
    P_hgcer_xAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_xAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_nohole_xAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_nohole_xAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_yAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_yAtCer_vs_MM_SUB_DUMMY_RAND"]
    P_hgcer_nohole_yAtCer_vs_MM_DUMMY_RAND = subDict["P_hgcer_nohole_yAtCer_vs_MM_SUB_DUMMY_RAND"]    
    
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
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DATA
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOHOLECUTS = apply_data_sub_cuts(evt)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
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

          # Phase shift to right setting
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
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOHOLECUTS = apply_data_sub_cuts(evt)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
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

          # Phase shift to right setting
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
          
          H_cal_etottracknorm_DUMMY.Fill(evt.H_cal_etottracknorm)
          H_cer_npeSum_DUMMY.Fill(evt.H_cer_npeSum)

          P_cal_etottracknorm_DUMMY.Fill(evt.P_cal_etottracknorm)
          P_hgcer_npeSum_DUMMY.Fill(evt.P_hgcer_npeSum)
          P_aero_npeSum_DUMMY.Fill(evt.P_aero_npeSum)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp
        
        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_RAND
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOHOLECUTS = apply_data_sub_cuts(evt)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
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

          # Phase shift to right setting
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
          #H_MM_RAND.Fill(pow(adj_MM, 2))  
          #H_MM_RAND.Fill(evt.Mrecoil)
          
          H_cal_etottracknorm_RAND.Fill(evt.H_cal_etottracknorm)
          H_cer_npeSum_RAND.Fill(evt.H_cer_npeSum)

          P_cal_etottracknorm_RAND.Fill(evt.P_cal_etottracknorm)
          P_hgcer_npeSum_RAND.Fill(evt.P_hgcer_npeSum)
          P_aero_npeSum_RAND.Fill(evt.P_aero_npeSum)
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY_RAND
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOHOLECUTS = apply_data_sub_cuts(evt)
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
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

          # Phase shift to right setting
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
          #H_MM_DUMMY_RAND.Fill(pow(adj_MM, 2))
          #H_MM_DUMMY_RAND.Fill(evt.Mrecoil)
          
          H_cal_etottracknorm_DUMMY_RAND.Fill(evt.H_cal_etottracknorm)
          H_cer_npeSum_DUMMY_RAND.Fill(evt.H_cer_npeSum)

          P_cal_etottracknorm_DUMMY_RAND.Fill(evt.P_cal_etottracknorm)
          P_hgcer_npeSum_DUMMY_RAND.Fill(evt.P_hgcer_npeSum)
          P_aero_npeSum_DUMMY_RAND.Fill(evt.P_aero_npeSum)

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
    H_t_DUMMY_RAND.Scale(1/nWindows)
    H_epsilon_DUMMY_RAND.Scale(1/nWindows)
    H_MM_DUMMY_RAND.Scale(1/nWindows)
    H_MM_nosub_DUMMY_RAND.Scale(1/nWindows)    
    H_pmiss_DUMMY_RAND.Scale(1/nWindows)
    H_emiss_DUMMY_RAND.Scale(1/nWindows)
    H_pmx_DUMMY_RAND.Scale(1/nWindows)
    H_pmy_DUMMY_RAND.Scale(1/nWindows)
    H_pmz_DUMMY_RAND.Scale(1/nWindows)
    H_W_DUMMY_RAND.Scale(1/nWindows)
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
    H_t_DUMMY.Add(H_t_DUMMY_RAND,-1)
    H_epsilon_DUMMY.Add(H_epsilon_DUMMY_RAND,-1)
    H_MM_DUMMY.Add(H_MM_DUMMY_RAND,-1)
    H_MM_nosub_DUMMY.Add(H_MM_nosub_DUMMY_RAND,-1)    
    H_pmiss_DUMMY.Add(H_pmiss_DUMMY_RAND,-1)
    H_emiss_DUMMY.Add(H_emiss_DUMMY_RAND,-1)
    H_pmx_DUMMY.Add(H_pmx_DUMMY_RAND,-1)
    H_pmy_DUMMY.Add(H_pmy_DUMMY_RAND,-1)
    H_pmz_DUMMY.Add(H_pmz_DUMMY_RAND,-1)
    H_W_DUMMY.Add(H_W_DUMMY_RAND,-1)
    H_ct_DUMMY.Add(H_ct_DUMMY_RAND,-1)

################################################################################################################################################

def particle_subtraction_ave(t_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg=None):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    MM_offset_DUMMY = subDict["MM_offset_DUMMY"]
    MM_offset_RAND = subDict["MM_offset_RAND"]
    MM_offset_DUMMY_RAND = subDict["MM_offset_DUMMY_RAND"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################
    # Defin of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################

    hist_dict = {}
    
    for j in range(len(t_bins)-1):
        
        hist_dict["H_Q2_DATA_{}".format(j)] = subDict["H_Q2_SUB_DATA_{}".format(j)]
        hist_dict["H_W_DATA_{}".format(j)] = subDict["H_W_SUB_DATA_{}".format(j)]
        hist_dict["H_t_DATA_{}".format(j)] = subDict["H_t_SUB_DATA_{}".format(j)]
        hist_dict["H_epsilon_DATA_{}".format(j)] = subDict["H_epsilon_SUB_DATA_{}".format(j)]
        hist_dict["H_MM_DATA_{}".format(j)] = subDict["H_MM_SUB_DATA_{}".format(j)]
        hist_dict["H_MM_nosub_DATA_{}".format(j)] = subDict["H_MM_nosub_SUB_DATA_{}".format(j)]        

        hist_dict["H_Q2_DUMMY_{}".format(j)] = subDict["H_Q2_SUB_DUMMY_{}".format(j)]
        hist_dict["H_W_DUMMY_{}".format(j)] = subDict["H_W_SUB_DUMMY_{}".format(j)]
        hist_dict["H_t_DUMMY_{}".format(j)] = subDict["H_t_SUB_DUMMY_{}".format(j)]
        hist_dict["H_epsilon_DUMMY_{}".format(j)] = subDict["H_epsilon_SUB_DUMMY_{}".format(j)]
        hist_dict["H_MM_DUMMY_{}".format(j)] = subDict["H_MM_SUB_DUMMY_{}".format(j)]
        hist_dict["H_MM_nosub_DUMMY_{}".format(j)] = subDict["H_MM_nosub_SUB_DUMMY_{}".format(j)]        

        hist_dict["H_Q2_RAND_{}".format(j)] = subDict["H_Q2_SUB_RAND_{}".format(j)]
        hist_dict["H_W_RAND_{}".format(j)] = subDict["H_W_SUB_RAND_{}".format(j)]
        hist_dict["H_t_RAND_{}".format(j)] = subDict["H_t_SUB_RAND_{}".format(j)]
        hist_dict["H_epsilon_RAND_{}".format(j)] = subDict["H_epsilon_SUB_RAND_{}".format(j)]
        hist_dict["H_MM_RAND_{}".format(j)] = subDict["H_MM_SUB_RAND_{}".format(j)]
        hist_dict["H_MM_nosub_RAND_{}".format(j)] = subDict["H_MM_nosub_SUB_RAND_{}".format(j)]        

        hist_dict["H_Q2_DUMMY_RAND_{}".format(j)] = subDict["H_Q2_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_W_DUMMY_RAND_{}".format(j)] = subDict["H_W_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_t_DUMMY_RAND_{}".format(j)] = subDict["H_t_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)] = subDict["H_epsilon_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_MM_DUMMY_RAND_{}".format(j)] = subDict["H_MM_SUB_DUMMY_RAND_{}".format(j)]
        hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)] = subDict["H_MM_nosub_SUB_DUMMY_RAND_{}".format(j)]        
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

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
            
    ##############
    ##############        
    ##############

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DATA
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_MM_nosub_DATA_{}".format(j)].Fill(adj_MM)            
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_Q2_DATA_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DATA_{}".format(j)].Fill(-evt.MandelT)
                    hist_dict["H_W_DATA_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DATA_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DATA_{}".format(j)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_MM_nosub_DUMMY_{}".format(j)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_Q2_DUMMY_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DUMMY_{}".format(j)].Fill(-evt.MandelT)
                    hist_dict["H_W_DUMMY_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DUMMY_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DUMMY_{}".format(j)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_RAND
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_MM_nosub_RAND_{}".format(j)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_Q2_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_RAND_{}".format(j)].Fill(-evt.MandelT)
                    hist_dict["H_W_RAND_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_RAND_{}".format(j)].Fill(adj_MM)
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY_RAND
        
        ##############
        ##############        
        ##############
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:                
                    hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Fill(adj_MM)                                
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                    hist_dict["H_Q2_DUMMY_RAND_{}".format(j)].Fill(evt.Q2)
                    hist_dict["H_t_DUMMY_RAND_{}".format(j)].Fill(-evt.MandelT)
                    hist_dict["H_W_DUMMY_RAND_{}".format(j)].Fill(evt.W)
                    hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Fill(evt.epsilon)
                    hist_dict["H_MM_DUMMY_RAND_{}".format(j)].Fill(adj_MM)
                  
    for j in range(len(t_bins)-1):

        # Data Random subtraction window
        hist_dict["H_Q2_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_W_RAND_{}".format(j)].Scale(1/nWindows)    
        hist_dict["H_t_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_epsilon_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_nosub_RAND_{}".format(j)].Scale(1/nWindows)        

        # Data Dummy_Random subtraction window
        hist_dict["H_Q2_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_W_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_t_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)
        hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)].Scale(1/nWindows)        

        ###
        # Data Random subtraction
        hist_dict["H_Q2_DATA_{}".format(j)].Add(hist_dict["H_Q2_RAND_{}".format(j)],-1)
        hist_dict["H_W_DATA_{}".format(j)].Add(hist_dict["H_W_RAND_{}".format(j)],-1)
        hist_dict["H_t_DATA_{}".format(j)].Add(hist_dict["H_t_RAND_{}".format(j)],-1)
        hist_dict["H_epsilon_DATA_{}".format(j)].Add(hist_dict["H_epsilon_RAND_{}".format(j)],-1)
        hist_dict["H_MM_DATA_{}".format(j)].Add(hist_dict["H_MM_RAND_{}".format(j)],-1)
        hist_dict["H_MM_nosub_DATA_{}".format(j)].Add(hist_dict["H_MM_nosub_RAND_{}".format(j)],-1)        

        ###
        # Dummy Random subtraction
        hist_dict["H_Q2_DUMMY_{}".format(j)].Add(hist_dict["H_Q2_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_W_DUMMY_{}".format(j)].Add(hist_dict["H_W_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_t_DUMMY_{}".format(j)].Add(hist_dict["H_t_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_epsilon_DUMMY_{}".format(j)].Add(hist_dict["H_epsilon_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_MM_DUMMY_{}".format(j)].Add(hist_dict["H_MM_DUMMY_RAND_{}".format(j)],-1)
        hist_dict["H_MM_nosub_DUMMY_{}".format(j)].Add(hist_dict["H_MM_nosub_DUMMY_RAND_{}".format(j)],-1)        
        
################################################################################################################################################

def particle_subtraction_yield(t_bins, phi_bins, subDict, inpDict, SubtractedParticle, hgcer_cutg=None):

    W = inpDict["W"] 
    Q2 = inpDict["Q2"]
    EPSSET = inpDict["EPSSET"]
    ParticleType = inpDict["ParticleType"]

    mm_min = inpDict["mm_min"] 
    mm_max = inpDict["mm_max"]
    
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 

    nWindows = subDict["nWindows"]
    phi_setting = subDict["phi_setting"]

    MM_offset_DATA = subDict["MM_offset_DATA"]
    MM_offset_DUMMY = subDict["MM_offset_DUMMY"]
    MM_offset_RAND = subDict["MM_offset_RAND"]
    MM_offset_DUMMY_RAND = subDict["MM_offset_DUMMY_RAND"]
    
    ################################################################################################################################################
    # Import function to define cut bools
    from apply_cuts import apply_data_cuts, apply_data_sub_cuts, set_val
    set_val(inpDict) # Set global variables for optimization
    
    ################################################################################################################################################
    # Define data root file trees of interest

    rootFileData = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDATAFilename}.root"
    if not os.path.isfile(rootFileData):
        print("\n\nERROR: No data file found called {}\n\n".format(rootFileData))
        sys.exit(2)

    InFile_DATA = open_root_file(rootFileData)

    TBRANCH_DATA  = InFile_DATA.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))

    TBRANCH_RAND  = InFile_DATA.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################
    # Define dummy root file trees of interest

    rootFileDummy = f"{OUTPATH}/{phi_setting}_{SubtractedParticle}_{InDUMMYFilename}.root"
    if not os.path.isfile(rootFileDummy):
        print("\n\nERROR: No dummy file found called {}\n\n".format(rootFileDummy))
        sys.exit(2)

    InFile_DUMMY = open_root_file(rootFileDummy)

    TBRANCH_DUMMY  = InFile_DUMMY.Get("Cut_{}_Events_prompt_noRF".format(SubtractedParticle.capitalize()))
    
    TBRANCH_DUMMY_RAND  = InFile_DUMMY.Get("Cut_{}_Events_rand_noRF".format(SubtractedParticle.capitalize()))

    ################################################################################################################################################

    hist_dict = {}
    
    # Loop through bins in t_data and identify events in specified bins
    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
    
            hist_dict["H_t_DATA_{}_{}".format(j, k)] = subDict["H_t_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_MM_DATA_{}_{}".format(j, k)] = subDict["H_MM_SUB_DATA_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DATA_{}_{}".format(j, k)]            

            hist_dict["H_t_DUMMY_{}_{}".format(j, k)] = subDict["H_t_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)] = subDict["H_MM_SUB_DUMMY_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DUMMY_{}_{}".format(j, k)]            

            hist_dict["H_t_RAND_{}_{}".format(j, k)] = subDict["H_t_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_RAND_{}_{}".format(j, k)] = subDict["H_MM_SUB_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_RAND_{}_{}".format(j, k)]            

            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_t_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_MM_SUB_DUMMY_RAND_{}_{}".format(j, k)]
            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)] = subDict["H_MM_nosub_SUB_DUMMY_RAND_{}_{}".format(j, k)]            
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1,0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    c0_dict = {}

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
            
    ##############
    ##############        
    ##############

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction data...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DATA):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DATA.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DATA
        
        ##############
        ##############        
        ##############
        
        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):            
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_t_DATA_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY

        ##############
        ##############        
        ##############

        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:                
                            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_t_DUMMY_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)].Fill(adj_MM)

    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_RAND

        ##############
        ##############        
        ##############

        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:                
                            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_t_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_dict["H_MM_RAND_{}_{}".format(j, k)].Fill(adj_MM)
          
    ################################################################################################################################################
    # Fill histograms for various trees called above

    print("\nGrabbing {} {} subtraction dummy random...".format(phi_setting,SubtractedParticle))
    for i,evt in enumerate(TBRANCH_DUMMY_RAND):

        # Progress bar
        Misc.progressBar(i, TBRANCH_DUMMY_RAND.GetEntries(),bar_length=25)        

        ##############
        # HARD CODED #
        ##############

        adj_hsdelta = evt.hsdelta + c0_dict["Q{}W{}_{}e".format(Q2,W,EPSSET)]*evt.hsxpfp

        # Check if variable shift branch exists
        try:
            adj_MM = evt.MM_shift
        except AttributeError:
            adj_MM = evt.MM + MM_offset_DUMMY_RAND

        ##############
        ##############        
        ##############


        # Phase shift to right setting
        #phi_shift = (evt.ph_q+math.pi)
        phi_shift = (evt.ph_q)        
        
        if ParticleType == "kaon":
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
            NOMMCUTS = apply_data_sub_cuts(evt) and not hgcer_cutg.IsInside(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer) and evt.P_hgcer_npeSum > 2.0
        else:
            ALLCUTS = apply_data_cuts(evt, mm_min, mm_max)
            NOMMCUTS = apply_data_sub_cuts(evt)

        if(NOMMCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)
            
        if(ALLCUTS):
            for j in range(len(t_bins)-1):
                for k in range(len(phi_bins)-1):
                    if t_bins[j] <= -evt.MandelT <= t_bins[j+1]:
                        if phi_bins[k] <= (phi_shift)*(180 / math.pi) <= phi_bins[k+1]:
                            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Fill(-evt.MandelT)
                            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Fill(adj_MM)

    for j in range(len(t_bins)-1):
        for k in range(len(phi_bins)-1):
            
            # Data Random subtraction window    
            hist_dict["H_t_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            

            # Data Dummy_Random subtraction window
            hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)
            hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)].Scale(1/nWindows)            

            ###
            # Data Random subtraction
            hist_dict["H_t_DATA_{}_{}".format(j, k)].Add(hist_dict["H_t_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_nosub_DATA_{}_{}".format(j, k)].Add(hist_dict["H_MM_nosub_RAND_{}_{}".format(j, k)],-1)            

            ###
            # Dummy Random subtraction
            hist_dict["H_t_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_t_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_MM_DUMMY_RAND_{}_{}".format(j, k)],-1)
            hist_dict["H_MM_nosub_DUMMY_{}_{}".format(j, k)].Add(hist_dict["H_MM_nosub_DUMMY_RAND_{}_{}".format(j, k)],-1)                        
