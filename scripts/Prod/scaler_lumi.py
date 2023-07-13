#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-07-13 10:48:29 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import uproot as up
import numpy as np
import math

################################################################################################################################################

def scaler(runNum, s_tree, thres_curr=3):

    '''
    SCALER TREE, TSP
    '''
    s_evts = s_tree.array("P.BCM4A.scaler")

    P_BCM4A_scalerCharge = s_tree.array("P.BCM4A.scalerCharge")
    P_BCM2_scalerCharge = s_tree.array("P.BCM2.scalerCharge")
    P_BCM4B_scalerCharge = s_tree.array("P.BCM4B.scalerCharge")
    P_BCM1_scalerCharge = s_tree.array("P.BCM1.scalerCharge")
    P_BCM4C_scalerCharge = s_tree.array("P.BCM4C.scalerCharge")

    P_BCM4A_scalerCurrent = s_tree.array("P.BCM4A.scalerCurrent")
    P_BCM2_scalerCurrent = s_tree.array("P.BCM2.scalerCurrent")
    P_BCM4B_scalerCurrent = s_tree.array("P.BCM4B.scalerCurrent")
    P_BCM1_scalerCurrent = s_tree.array("P.BCM1.scalerCurrent")
    P_BCM4C_scalerCurrent = s_tree.array("P.BCM4C.scalerCurrent")

    P_1Mhz_scalerTime = s_tree.array("P.1MHz.scalerTime")

    P_pTRIG1_scaler = s_tree.array("P.pTRIG1.scaler")
    P_pTRIG2_scaler = s_tree.array("P.pTRIG2.scaler")
    P_pTRIG3_scaler = s_tree.array("P.pTRIG3.scaler")
    P_pTRIG4_scaler = s_tree.array("P.pTRIG4.scaler")
    P_pTRIG5_scaler = s_tree.array("P.pTRIG5.scaler")
    P_pTRIG6_scaler = s_tree.array("P.pTRIG6.scaler")

    P_pL1ACCP_scaler = s_tree.array("P.pL1ACCP.scaler")
    P_pPRE40_scaler = s_tree.array("P.pPRE40.scaler")
    P_pPRE100_scaler = s_tree.array("P.pPRE100.scaler")
    P_pPRE150_scaler = s_tree.array("P.pPRE150.scaler")
    P_pPRE200_scaler = s_tree.array("P.pPRE200.scaler")
    P_pPRE40_scaler = s_tree.array("P.pPRE40.scaler")
    P_pPRE100_scaler = s_tree.array("P.pPRE100.scaler")
    P_pPRE150_scaler = s_tree.array("P.pPRE150.scaler")
    P_pPRE200_scaler = s_tree.array("P.pPRE200.scaler")

    P_pEL_LO_LO_scaler = s_tree.array("P.pEL_LO_LO.scaler")
    P_pEL_LO_scaler = s_tree.array("P.pEL_LO.scaler")
    P_pEL_HI_scaler = s_tree.array("P.pEL_HI.scaler")
    P_pEL_REAL_scaler = s_tree.array("P.pEL_REAL.scaler")
    P_pEL_CLEAN_scaler = s_tree.array("P.pEL_CLEAN.scaler")
    P_pSTOF_scaler = s_tree.array("P.pSTOF.scaler")

    P_pEL_LO_LO_scaler = s_tree.array("P.pEL_LO_LO.scaler")
    P_pEL_LO_scaler = s_tree.array("P.pEL_LO.scaler")
    P_pEL_HI_scaler = s_tree.array("P.pEL_HI.scaler")
    P_pEL_REAL_scaler = s_tree.array("P.pEL_REAL.scaler")
    P_pEL_CLEAN_scaler = s_tree.array("P.pEL_CLEAN.scaler")
    P_pSTOF_scaler = s_tree.array("P.pSTOF.scaler")
    P_pPRHI_scaler = s_tree.array("P.PRHI.scaler")
    P_pPRLO_scaler = s_tree.array("P.PRLO.scaler")

    P_EDTM_scaler = s_tree.array("P.EDTM.scaler")
    
    NBCM = 5
    NTRIG = 6
    NPRE = 4
    NRATE = 6
    SHMSNRATE = 8

    bcm_name = ["BCM1 ", "BCM2 ", "BCM4A", "BCM4B", "BCM4C"]

    trig_name = ["TRIG1", "TRIG2", "TRIG3", "TRIG4", "TRIG5", "TRIG6"]

    PRE_name = ["40", "100", "150", "200"]

    rate_name = ["EL_LO_LO", "EL_LO", "EL_HI", "EL_REAL", "EL_CLEAN", "STOF"]

    SHMS_rate_name = ["EL_LO_LO", "EL_LO", "EL_HI",
                      "EL_REAL", "EL_CLEAN", "STOF", "PR_HI", "PR_LO"]

    bcm_value = [P_BCM1_scalerCharge, P_BCM2_scalerCharge,
                 P_BCM4A_scalerCharge, P_BCM4B_scalerCharge, P_BCM4C_scalerCharge]

    time_value = P_1Mhz_scalerTime

    current = [P_BCM1_scalerCurrent, P_BCM2_scalerCurrent,
                 P_BCM4A_scalerCurrent, P_BCM4B_scalerCurrent, P_BCM4C_scalerCurrent]

    trig_value = [P_pTRIG1_scaler, P_pTRIG2_scaler, P_pTRIG3_scaler,P_pTRIG4_scaler, P_pTRIG5_scaler, P_pTRIG6_scaler]
    #trig_value = [P_pTRIG1_scaler, P_pTRIG2_scaler, P_pEL_CLEAN_scaler,P_pTRIG4_scaler, P_pTRIG5_scaler, P_pTRIG6_scaler]

    hms_el_clean = P_pEL_CLEAN_scaler

    acctrig_value = P_pL1ACCP_scaler

    PRE_value = [P_pPRE40_scaler, P_pPRE100_scaler,
                 P_pPRE150_scaler, P_pPRE200_scaler]

    SHMS_PRE_value = [P_pPRE40_scaler, P_pPRE100_scaler,
                      P_pPRE150_scaler, P_pPRE200_scaler]

    rate_value = [P_pEL_LO_LO_scaler, P_pEL_LO_scaler, P_pEL_HI_scaler,
                  P_pEL_REAL_scaler, P_pEL_CLEAN_scaler, P_pSTOF_scaler]

    SHMS_rate_value = [P_pEL_LO_LO_scaler, P_pEL_LO_scaler, P_pEL_HI_scaler,
                       P_pEL_REAL_scaler, P_pEL_CLEAN_scaler, P_pSTOF_scaler, P_pPRHI_scaler, P_pPRLO_scaler]

    EDTM_value = P_EDTM_scaler

    # Variables useful in Process
    # To find total charge
    name = [0]*NBCM
    charge_sum = [0]*NBCM
    time_sum = [0]*NBCM
    previous_charge = [0]*NBCM
    previous_time = 0
    previous_time = [0]*NBCM
    current_I = 0
    current_time = 0

    # To determine computer livetime
    name = [0]*NTRIG
    trig_sum = [0]*NTRIG
    hms_el_clean_sum = [0]
    previous_trig = [0]*NTRIG
    previous_hms_el_clean = [0]
    pretrigger = 0
    previous_pretrigger = 0
    acctrig_sum = 0
    previous_acctrig = 0

    # To determine HMS electronic livetime
    name = [0]*NPRE
    PRE_sum = [0]*NPRE
    previous_PRE = [0]*NPRE

    # To determine SHMS electronic livetime
    SHMS_PRE_sum = [0]*NPRE
    SHMS_previous_PRE = [0]*NPRE

    # To determine HMS trigger rates
    name = [0]*NRATE
    rate_sum = [0]*NRATE
    previous_rate = [0]*NRATE

    # To determine SHMS trigger rates
    rate_name = [0]*SHMSNRATE
    SHMS_rate_sum = [0]*SHMSNRATE
    SHMS_previous_rate = [0]*SHMSNRATE

    # To determine number of EDTM events
    EDTM_sum = 0
    EDTM_current = 0
    previous_EDTM = 0
    
    # Set bcm to use (0-bcm1, 1-bcm2, 2-bcm4A, 3-bcm4B, 4-bcm4C)
    bcm_ix = 0
    
    for ibcm in range(0, 5):
        previous_acctrig = (acctrig_value[0] - EDTM_current)
        previous_EDTM = EDTM_value[0]
        for itrig in range(0, NTRIG):
            previous_trig[itrig] = trig_value[itrig][0]
            previous_hms_el_clean = hms_el_clean[0]
        for iPRE in range(0, NPRE):
            previous_PRE[iPRE] = PRE_value[iPRE][0]
            SHMS_previous_PRE[iPRE] = SHMS_PRE_value[iPRE][0]
        for iRATE in range(0, NRATE):
            previous_rate[iRATE] = rate_value[iRATE][0]
        for iRATE in range(0, SHMSNRATE):
            SHMS_previous_rate[iRATE] = SHMS_rate_value[iRATE][0]
        previous_time[ibcm] = time_value[0]
        previous_charge[ibcm] = bcm_value[ibcm][0]
        # Iterate over all scaler events to get various scaler values
        for i, evt in enumerate(s_evts):
            if (time_value[i] != previous_time[ibcm]):
                # Current calculation using iterative charge and time values.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                #current_I = (bcm_value[ibcm][i] - previous_charge[ibcm])/(time_value[i] - previous_time[ibcm])
                current_I = current[ibcm][i]
            if (current[ibcm][i] > thres_curr ):
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                charge_sum[ibcm] += (bcm_value[ibcm][i] - previous_charge[ibcm])
                time_sum[ibcm] += (time_value[i] - previous_time[ibcm])
            # Current cuts and selection of BCM1
            if (ibcm == bcm_ix and current[ibcm][i] > thres_curr):
                # EDTM scaler iteration.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                EDTM_current = (EDTM_value[i] - previous_EDTM)
                EDTM_sum += EDTM_current
                # Accquired trigger sum calculation using iterative level 1 accepted values.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                # Changing acctrig to not subtract EDTM to get CPULT for all events
                # acctrig_sum += ((acctrig_value[i] - EDTM_current) - previous_acctrig)
                acctrig_sum += ((acctrig_value[i]) - previous_acctrig)
                hms_el_clean_sum += (hms_el_clean[i] - previous_hms_el_clean) # HMS EL_CLEAN
                for itrig in range(0, NTRIG):
                    # Trigger scaler iteration.
                    # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                    trig_sum[itrig] += (trig_value[itrig][i] - previous_trig[itrig])
                    #print("trig_value[%s] = " %(itrig),trig_value[itrig][i])
                    #print("previous_trig[%s] = " %(itrig),previous_trig[itrig])
                for iPRE in range(0, NPRE):
                    # Pre-trig scaler iteration. Used in electronic LT calculations.
                    # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                    PRE_sum[iPRE] += (PRE_value[iPRE][i] - previous_PRE[iPRE])
                    SHMS_PRE_sum[iPRE] += (SHMS_PRE_value[iPRE][i] - SHMS_previous_PRE[iPRE])
                for iRATE in range(0, NRATE):
                    rate_sum[iRATE] += (rate_value[iRATE][i] - previous_rate[iRATE])
                for iRATE in range(0, SHMSNRATE):
                    SHMS_rate_sum[iRATE] += (SHMS_rate_value[iRATE][i] - SHMS_previous_rate[iRATE])
            # Changing acctrig to not subtract EDTM to get CPULT for all events
            # previous_acctrig = (acctrig_value[i] - EDTM_current)
            previous_acctrig = (acctrig_value[i])
            previous_EDTM = EDTM_value[i]
            for itrig in range(0, NTRIG):
                previous_trig[itrig] = trig_value[itrig][i]
            for iPRE in range(0, NPRE):
                previous_PRE[iPRE] = PRE_value[iPRE][i]
                SHMS_previous_PRE[iPRE] = SHMS_PRE_value[iPRE][i]
            for iRATE in range(0, NRATE):
                previous_rate[iRATE] = rate_value[iRATE][i]
            for iRATE in range(0, SHMSNRATE):
                SHMS_previous_rate[iRATE] = SHMS_rate_value[iRATE][i]
            previous_time[ibcm] = time_value[i]
            previous_charge[ibcm] = bcm_value[ibcm][i]

        
    # Creates a dictionary for the calculated luminosity values 
    scalers = {
        "run number" : runNum,
        "time": time_sum[bcm_ix],
        "charge": charge_sum[bcm_ix],
        "current": charge_sum[bcm_ix]/time_sum[bcm_ix],
        "sent_edtm": EDTM_sum,
        "HMS_EL_CLEAN_scaler" : hms_el_clean_sum,
        "HMS_EL_CLEAN_scaler_accp" : hms_el_clean_sum,
        "rate_HMS_EL_CLEAN": hms_el_clean_sum/time_sum[bcm_ix],
        "yield_HMS_scaler": (hms_el_clean_sum-EDTM_sum)/charge_sum[bcm_ix],
        "uncern_yieldRel_HMS_scaler": np.sqrt(hms_el_clean_sum)/hms_el_clean_sum,
            
    }
    
    return scalers
