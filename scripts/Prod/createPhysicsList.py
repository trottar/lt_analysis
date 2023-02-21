#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-21 01:11:16 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os
import ROOT
import numpy as np

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=31:
    print("!!!!! ERROR !!!!!\n Expected 31 arguments\n Usage is with - Q2 POL EPSVAL TMIN TMAX NumtBins NumPhiBins Kset runNumRight runNumLeft runNumCenter pThetaValRight pThetaValLeft pThetaValCenter EbeamValRight EbeamValLeft EbeamValCenter EffValRight EffValLeft EffValCenter EffErrRight EffErrLeft EffErrCenter ChargeValRight ChargeValLeft ChargeValCenter ChargeErrRight ChargeErrLeft ChargeErrCenter kinematics OutFullAnalysisFilename\n!!!!! ERROR !!!!!")
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
NumtBins = int(sys.argv[6])
NumPhiBins = int(sys.argv[7])
Kset = sys.argv[8]

runNumRight = list(sys.argv[9].split(" "))
runNumLeft = list(sys.argv[10].split(" "))
runNumCenter = list(sys.argv[11].split(" "))

pThetaValRight = list(sys.argv[12].split(" "))
pThetaValLeft = list(sys.argv[13].split(" "))
pThetaValCenter = list(sys.argv[14].split(" "))
EbeamValRight = list(sys.argv[15].split(" "))
EbeamValLeft = list(sys.argv[16].split(" "))
EbeamValCenter = list(sys.argv[17].split(" "))

EffValRight = list(sys.argv[18].split(" "))
EffValLeft = list(sys.argv[19].split(" "))
EffValCenter = list(sys.argv[20].split(" "))
EffErrRight = list(sys.argv[21].split(" "))
EffErrLeft = list(sys.argv[22].split(" "))
EffErrCenter = list(sys.argv[23].split(" "))

ChargeValRight = list(sys.argv[24].split(" "))
ChargeValLeft = list(sys.argv[25].split(" "))
ChargeValCenter = list(sys.argv[26].split(" "))
ChargeErrRight = list(sys.argv[27].split(" "))
ChargeErrLeft = list(sys.argv[28].split(" "))
ChargeErrCenter = list(sys.argv[29].split(" "))

kinematics = sys.argv[30].split("_")

OutFilename = sys.argv[31]

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

InFile_DATA = ROOT.TFile.Open(InDATAFilename,"READ")

###############################################################################################################################################

if float(runNumRight[0]) != 0:
    yield_right_data = []
    yield_right_simc = []
    phibin_right_data = []
    tbin_right_data = []
    averQ2_right_data = []
    averW_right_data = []
    avert_right_data = []
    TBRANCH_RIGHT = InFile_DATA.Get("Right")
    for i, evt in enumerate(TBRANCH_RIGHT):
        if i <= NumtBins*NumPhiBins:
            yield_right_data.append(evt.yield_data)
            yield_right_simc.append(evt.yield_simc)
            phibin_right_data.append(evt.phibins)
            tbin_right_data.append(evt.tbins)
            averQ2_right_data.append(evt.aver_Q2)
            averW_right_data.append(evt.aver_W)
            avert_right_data.append(evt.aver_t)
    relyield_right = [0.0 if s == 0.0 else 0.0 if np.isnan(d) else 0.0 if np.isnan(s) else d/s for d,s in zip(yield_right_data,yield_right_simc)] 

if float(runNumLeft[0]) != 0:
    yield_left_data = []
    yield_left_simc = []
    phibin_left_data = []
    tbin_left_data = []
    averQ2_left_data = []
    averW_left_data = []
    avert_left_data = []
    TBRANCH_LEFT = InFile_DATA.Get("Left")
    for i, evt in enumerate(TBRANCH_LEFT):
        if i <= NumtBins*NumPhiBins:
            yield_left_data.append(evt.yield_data)
            yield_left_simc.append(evt.yield_simc)
            phibin_left_data.append(evt.phibins)
            tbin_left_data.append(evt.tbins)
            averQ2_left_data.append(evt.aver_Q2)
            averW_left_data.append(evt.aver_W)
            avert_left_data.append(evt.aver_t)
    relyield_left = [0.0 if s == 0.0 else 0.0 if np.isnan(d) else 0.0 if np.isnan(s) else d/s for d,s in zip(yield_left_data,yield_left_simc)] 

if float(runNumCenter[0]) != 0:
    yield_center_data = []
    yield_center_simc = []
    phibin_center_data = []
    tbin_center_data = []
    averQ2_center_data = []
    averW_center_data = []
    avert_center_data = []
    TBRANCH_CENTER = InFile_DATA.Get("Center")
    for i, evt in enumerate(TBRANCH_CENTER):
        if i <= NumtBins*NumPhiBins:
            yield_center_data.append(evt.yield_data)
            yield_center_simc.append(evt.yield_simc)
            phibin_center_data.append(evt.phibins)
            tbin_center_data.append(evt.tbins)
            averQ2_center_data.append(evt.aver_Q2)
            averW_center_data.append(evt.aver_W)
            avert_center_data.append(evt.aver_t)
    relyield_center = [0.0 if s == 0.0 else 0.0 if np.isnan(d) else 0.0 if np.isnan(s) else d/s for d,s in zip(yield_center_data,yield_center_simc)] 
    
print("\n\n~~~~~~~~~",relyield_left)
print("~~~~~~~~~",tbin_left_data)
print("~~~~~~~~~",phibin_left_data)
print("~~~~~~~~~",averQ2_left_data)
print("~~~~~~~~~",averW_left_data)
print("~~~~~~~~~",avert_left_data)

print("\n\n~~~~~~~~~",len(relyield_left))
print("~~~~~~~~~",len(tbin_left_data))
print("~~~~~~~~~",len(phibin_left_data))
print("~~~~~~~~~",len(averQ2_left_data))
print("~~~~~~~~~",len(averW_left_data))
print("~~~~~~~~~",len(avert_left_data))

print("\n\n~~~~~~~~~",relyield_center)
print("~~~~~~~~~",tbin_center_data)
print("~~~~~~~~~",phibin_center_data)
print("~~~~~~~~~",averQ2_center_data)
print("~~~~~~~~~",averW_center_data)
print("~~~~~~~~~",avert_center_data)

InFile_DATA.Close()

################################################################################################################################################

def write_to_file(f_out,line):
    # Open a file in append mode
    with open(f_out, 'a') as f:
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
            thpq_right = abs(float(pThetaValCenter[i])-float(pThetaValRight[i]))            
        else:
            continue
        
if float(runNumLeft[0]) != 0:
    runNums= runNumLeft
    for i, run in enumerate(runNumLeft):
        runNum = run
        pid_log = "%s/log/Analysed_Prod_%s.log" % (LTANAPATH,runNum)
        if os.path.exists(pid_log):
            thpq_left = abs(float(pThetaValCenter[i])-float(pThetaValLeft[i]))            
        else:
            continue

if float(runNumCenter[0]) != 0:
    thpq_center = 0.000

print(thpq_left)
    
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

    f_list = '{}/src/averages/aver.{}_{}_{}_{:.0f}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        if float(runNumRight[0]) != 0:
            for i, relyield in enumerate(relyield_right):
                if relyield != 0.0:
                    check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_right_data[i], tbin_right_data[i])
                    # Check if the line already exists
                    if check_line not in lines:
                        write_to_file(f_list,check_line)

        if float(runNumLeft[0]) != 0:
            for i, relyield in enumerate(relyield_left):
                if relyield != 0.0:
                    check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_left_data[i], tbin_left_data[i])
                    # Check if the line already exists
                    if check_line not in lines:
                        write_to_file(f_list,check_line)

        if float(runNumCenter[0]) != 0:
            for i, relyield in enumerate(relyield_center):
                if relyield != 0.0:
                    check_line = "{:.4f} {:.4f} {} {}\n".format(relyield, 1.0000, phibin_center_data[i], tbin_center_data[i])
                    # Check if the line already exists
                    if check_line not in lines:
                        write_to_file(f_list,check_line)
                
################################################################################################################################################

if float(runNumRight[0]) != 0:
    #f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_-{}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100, int(thpq_right*1000))
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_-3000.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_right_data):
            if relyield != 0.0:
                check_line = "{:.4f} {:.4f} {} {}\n".format(averQ2_right_data[i], 1.0, averW_right_data[i], 1.0, avert_right_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)

if float(runNumLeft[0]) != 0:
    #f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_+{}.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100, int(thpq_left*1000))
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_+3000.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_left_data):
            if relyield != 0.0:
                check_line = "{:.4f} {:.4f} {} {}\n".format(averQ2_left_data[i], 1.0, averW_left_data[i], 1.0, avert_left_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)

if float(runNumCenter[0]) != 0:
    f_list = '{}/src/kindata/kindata.{}_{}_{}_{:.0f}_+0000.dat'.format(LTANAPATH, PID, POL, Q2.replace(".",""), float(EPSVAL)*100)

    if not os.path.exists(f_list):
        open(f_list, "w").close()    
    # Open a file in read mode
    with open(f_list, 'r') as f:
        lines = f.readlines()
        for i, relyield in enumerate(relyield_center_data):
            if relyield != 0.0:
                check_line = "{:.4f} {:.4f} {} {}\n".format(averQ2_center_data[i], 1.0, averW_center_data[i], 1.0, avert_center_data[i])
                # Check if the line already exists
                if check_line not in lines:
                    write_to_file(f_list,check_line)

################################################################################################################################################

