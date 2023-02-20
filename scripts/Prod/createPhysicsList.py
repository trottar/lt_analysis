#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-20 03:59:15 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os
import ROOT
import uproot as up

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=30:
    print("!!!!! ERROR !!!!!\n Expected 29 arguments\n Usage is with - Q2 POL EPSVAL TMIN TMAX NumtBins Kset runNumRight runNumLeft runNumCenter pThetaValRight pThetaValLeft pThetaValCenter EbeamValRight EbeamValLeft EbeamValCenter EffValRight EffValLeft EffValCenter EffErrRight EffErrLeft EffErrCenter ChargeValRight ChargeValLeft ChargeValCenter ChargeErrRight ChargeErrLeft ChargeErrCenter kinematics OutFullAnalysisFilename\n!!!!! ERROR !!!!!")
    sys.exit(1)

################################################################################################################################################
'''
User Inputs
'''

Q2 = sys.argv[1].replace("p",".")
POL = sys.argv[2]
EPSVAL = sys.argv[3]
TMIN = sys.argv[4]
TMAX = sys.argv[5]
NumtBins = sys.argv[6]
Kset = sys.argv[7]

runNumRight = list(sys.argv[8].split(" "))
runNumLeft = list(sys.argv[9].split(" "))
runNumCenter = list(sys.argv[10].split(" "))

pThetaValRight = list(sys.argv[11].split(" "))
pThetaValLeft = list(sys.argv[12].split(" "))
pThetaValCenter = list(sys.argv[13].split(" "))

EbeamValRight = list(sys.argv[14].split(" "))
EbeamValLeft = list(sys.argv[15].split(" "))
EbeamValCenter = list(sys.argv[16].split(" "))

EffValRight = list(sys.argv[17].split(" "))
EffValLeft = list(sys.argv[18].split(" "))
EffValCenter = list(sys.argv[19].split(" "))
EffErrRight = list(sys.argv[20].split(" "))
EffErrLeft = list(sys.argv[21].split(" "))
EffErrCenter = list(sys.argv[22].split(" "))

ChargeValRight = list(sys.argv[23].split(" "))
ChargeValLeft = list(sys.argv[24].split(" "))
ChargeValCenter = list(sys.argv[25].split(" "))
ChargeErrRight = list(sys.argv[26].split(" "))
ChargeErrLeft = list(sys.argv[27].split(" "))
ChargeErrCenter = list(sys.argv[28].split(" "))

kinematics = sys.argv[29].split("_")

OutFilename = sys.argv[30]

InSIMCFilenameRight = "Prod_Coin_{}.root".format(kinematics[0]+"right_"+kinematics[1])
InSIMCFilenameLeft = "Prod_Coin_{}.root".format(kinematics[0]+"left_"+kinematics[1])
InSIMCFilenameCenter = "Prod_Coin_{}.root".format(kinematics[0]+"center_"+kinematics[1])

particle = "kaon"

if particle == "kaon":
    PID = "k"
elif particle == "pion":
    PID = "pi"

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

InDATAFilename = OUTPATH+"/" + OutFilename + ".root"

InFile_DATA = up.open(InDATAFilename)

###############################################################################################################################################

if float(runNumRight[0]) != 0:
    yield_right_data = InFile_DATA["Right Data/H_yield_DATA"].values
    phibin_right_data = InFile_DATA["Right Data/H_phibin_DATA"].values
    tbin_right_data = InFile_DATA["Right Data/H_tbin_DATA"].values

if float(runNumLeft[0]) != 0:
    yield_left_data = InFile_DATA["Left Data/H_yield_DATA"].values
    phibin_left_data = InFile_DATA["Left Data/H_phibin_DATA"].values
    tbin_left_data = InFile_DATA["Left Data/H_tbin_DATA"].values

if float(runNumCenter[0]) != 0:
    yield_center_data = InFile_DATA["Center Data/H_yield_DATA"].values
    phibin_center_data = InFile_DATA["Center Data/H_phibin_DATA"].values
    tbin_center_data = InFile_DATA["Center Data/H_tbin_DATA"].values

InFile_DATA.Close()

print("\n\n~~~~~~~~~",yield_left_data)
print("~~~~~~~~~",tbin_left_data)
print("~~~~~~~~~",phibin_left_data)

print("\n\n~~~~~~~~~",len(yield_left_data))
print("~~~~~~~~~",len(tbin_left_data))
print("~~~~~~~~~",len(phibin_left_data))

print("\n\n~~~~~~~~~",yield_center_data)
print("~~~~~~~~~",tbin_center_data)
print("~~~~~~~~~",phibin_center_data)

################################################################################################################################################

def write_to_file(f_out,line):
    # Open a file in append mode
    with open(f_out, 'a') as f:
        # Write the value of the variable to the file
        f.write(line)

################################################################################################################################################

# Define thpq vector relative to middle setting
if float(runNumRight[0]) != 0:
    thpq_right = abs(float(pThetaValCenter[0])-float(pThetaValRight[0]))
if float(runNumLeft[0]) != 0:
    thpq_left = abs(float(pThetaValCenter[0])-float(pThetaValLeft[0]))
if float(runNumCenter[0]) != 0:
    thpq_center = 0.000

################################################################################################################################################
            
f_list_settings = '{}/src/list.settings'.format(LTANAPATH)
if not os.path.exists(f_list_settings):
    open(f_list_settings, "w").close()
# First check if line exists
with open(f_list_settings, 'r') as f:
    lines = f.readlines()
    if float(runNumRight[0]) != 0:
        check_line = "{} {} {} -{:.3f} {} {} {} {}\n".format(POL, Q2, EPSVAL, thpq_right, TMIN, TMAX, NumtBins, Kset)
        # Check if the line already exists
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
    if float(runNumLeft[0]) != 0:
        check_line = "{} {} {} +{:.3f} {} {} {} {}\n".format(POL, Q2, EPSVAL, thpq_left, TMIN, TMAX, NumtBins, Kset)
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
    if float(runNumCenter[0]) != 0:
        check_line = "{} {} {} {:.3f} {} {} {} {}\n".format(POL, Q2, EPSVAL, thpq_center, TMIN, TMAX, NumtBins, Kset)
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
            
################################################################################################################################################

if float(runNumRight[0]) != 0:
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_{:.3f}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100, thpq_right)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_right_data):
            if relyield == 0.0:
                # convert uC to C (10^-6C=1uC)
                check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_right_data[i], tbin_right_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)

if float(runNumLeft[0]) != 0:
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_{:.3f}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100, thpq_left)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_left_data):
            if relyield == 0.0:
                # convert uC to C (10^-6C=1uC)
                check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_left_data[i], tbin_left_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)

if float(runNumCenter[0]) != 0:
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_{:.3f}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100, thpq_center)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_center_data):
            if relyield == 0.0:
                # convert uC to C (10^-6C=1uC)
                check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_center_data[i], tbin_center_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)
                
################################################################################################################################################
