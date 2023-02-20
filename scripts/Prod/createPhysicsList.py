#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-02-19 23:01:34 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os
import ROOT

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

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

###############################################################################################################################################


InDATAFilename = OUTPATH+"/" + OutFilename + ".root"

InFile_DATA = ROOT.TFile.Open(InDATAFilename, "OPEN")

TBRANCH_RIGHT_DATA  = InFile_DATA.Get("Right Data")
TBRANCH_LEFT_DATA  = InFile_DATA.Get("Left Data")
TBRANCH_CENTER_DATA  = InFile_DATA.Get("Center Data")

###############################################################################################################################################

relyield_right_data = []
print("\nGrabbing right data yield...")
for i,evt in enumerate(TBRANCH_RIGHT_DATA):
    # Progress bar
      Misc.progressBar(i, TBRANCH_RIGHT_DATA.GetEntries(),bar_length=25)      
      relyield_right_data.append(evt.H_relyield_DATA)
      
relyield_left_data = []
print("\nGrabbing left data yield...")
for i,evt in enumerate(TBRANCH_LEFT_DATA):
    # Progress bar
      Misc.progressBar(i, TBRANCH_LEFT_DATA.GetEntries(),bar_length=25)      
      relyield_left_data.append(evt.H_relyield_DATA)
      
relyield_center_data = []
print("\nGrabbing center data yield...")
for i,evt in enumerate(TBRANCH_CENTER_DATA):
    # Progress bar
      Misc.progressBar(i, TBRANCH_CENTER_DATA.GetEntries(),bar_length=25)      
      relyield_center_data.append(evt.H_relyield_DATA)
      
###############################################################################################################################################

# Grabs simc number of events and normalizaton factor
simc_right_hist = "%s/OUTPUT/Analysis/%sLT/%s" % (LTANAPATH,ANATYPE,InSIMCFilenameRight.replace('.root','.hist'))
try:
    f_simc_right = open(simc_right_hist)
    angle_flag = False
    for line in f_simc_right:
        #print(line)
        if "Ngen" in line:
            val = line.split("=")
            simc_right_nevents = int(val[1])
        if "normfac" in line:
            val = line.split("=")
            simc_right_normfactor = float(val[1])
    if 'simc_right_nevents' and 'simc_right_normfactor' in locals():
        print('\n\nsimc_right_nevents = ',simc_right_nevents,'\nsimc_right_normfactor = ',simc_right_normfactor,'\n\n')
    else:
        print("ERROR: Invalid simc right hist file %s" % simc_right_hist)
        #sys.exit(1)
    f_simc_right.close()
except FileNotFoundError:
    print("No right simc file found...")
    #sys.exit(1)
    
simc_left_hist = "%s/OUTPUT/Analysis/%sLT/%s" % (LTANAPATH,ANATYPE,InSIMCFilenameLeft.replace('.root','.hist'))
try:
    f_simc_left = open(simc_left_hist)
    angle_flag = False
    for line in f_simc_left:
        #print(line)
        if "Ngen" in line:
            val = line.split("=")
            simc_left_nevents = int(val[1])
        if "normfac" in line:
            val = line.split("=")
            simc_left_normfactor = float(val[1])
    if 'simc_left_nevents' and 'simc_left_normfactor' in locals():
        print('\n\nsimc_left_nevents = ',simc_left_nevents,'\nsimc_left_normfactor = ',simc_left_normfactor,'\n\n')
    else:
        print("ERROR: Invalid simc left hist file %s" % simc_left_hist)
        #sys.exit(1)
    f_simc_left.close()
except FileNotFoundError:
    print("No left simc file found...")
    sys.exit(1)
    
try:    
    simc_center_hist = "%s/OUTPUT/Analysis/%sLT/%s" % (LTANAPATH,ANATYPE,InSIMCFilenameCenter.replace('.root','.hist'))
    f_simc_center = open(simc_center_hist)
    angle_flag = False
    for line in f_simc_center:
        #print(line)
        if "Ngen" in line:
            val = line.split("=")
            simc_center_nevents = int(val[1])
        if "normfac" in line:
            val = line.split("=")
            simc_center_normfactor = float(val[1])
    if 'simc_center_nevents' and 'simc_center_normfactor' in locals():
        print('\n\nsimc_center_nevents = ',simc_center_nevents,'\nsimc_center_normfactor = ',simc_center_normfactor,'\n\n')
    else:
        print("ERROR: Invalid simc center hist file %s" % simc_center_hist)
        #sys.exit(1)
    f_simc_center.close()
except FileNotFoundError:
    print("No center simc file found...")
    sys.exit(1)
    
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
thpq_left = abs(float(pThetaValCenter[0])-float(pThetaValLeft[0]))
thpq_center = 0.000

################################################################################################################################################

f_list = '{}/src/root_ana/list.simc'.format(LTANAPATH)
if not os.path.exists(f_list):
    open(f_list, "w").close()    
# Open a file in read mode
with open(f_list, 'r') as f:
    lines = f.readlines()
    try:
        if float(runNumRight[0]) != 0:
            # Write the value of the variable to the file
            check_line = "{} {:.5f} {} -{:.3f} {} {}\n".format(Q2,float(EbeamValRight[0]),EPSVAL,thpq_right,simc_right_normfactor,simc_right_nevents)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list,check_line)
    except NameError:
        print("")                
    if float(runNumLeft[0]) != 0:
        check_line = "{} {:.5f} {} +{:.3f} {} {}\n".format(Q2,float(EbeamValLeft[0]),EPSVAL,thpq_left,simc_left_normfactor,simc_left_nevents)
        if check_line not in lines:
            write_to_file(f_list,check_line)
    if float(runNumCenter[0]) != 0:
        check_line = "{} {:.5f} {} {:.3f} {} {}\n".format(Q2,float(EbeamValCenter[0]),EPSVAL,thpq_center,simc_center_normfactor,simc_center_nevents)
        if check_line not in lines:
            write_to_file(f_list,check_line)
            
################################################################################################################################################
            
f_list_settings = '{}/src/root_ana/list.settings'.format(LTANAPATH)
if not os.path.exists(f_list_settings):
    open(f_list_settings, "w").close()
# First check if line exists
with open(f_list_settings, 'r') as f:
    lines = f.readlines()
    if float(runNumRight[0]) != 0:
        check_line = "{} {} {} -{:.3f} {} {} {} {}\n".format(POL,Q2,EPSVAL,thpq_right,TMIN,TMAX,NumtBins,Kset)
        # Check if the line already exists
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
    if float(runNumLeft[0]) != 0:
        check_line = "{} {} {} +{:.3f} {} {} {} {}\n".format(POL,Q2,EPSVAL,thpq_left,TMIN,TMAX,NumtBins,Kset)
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
    if float(runNumCenter[0]) != 0:
        check_line = "{} {} {} {:.3f} {} {} {} {}\n".format(POL,Q2,EPSVAL,thpq_center,TMIN,TMAX,NumtBins,Kset)
        if check_line not in lines:
            write_to_file(f_list_settings,check_line)
            
################################################################################################################################################
            
f_list = '{}/src/lists/list.{}_{:.0f}'.format(LTANAPATH,Q2.replace(".",""),float(EPSVAL)*100)

if not os.path.exists(f_list):
    open(f_list, "w").close()    
# Open a file in read mode
with open(f_list, 'r') as f:
    lines = f.readlines()
    if float(runNumRight[0]) != 0:
        # Write the value of the variable to the file
        for i,thpq in enumerate(EbeamValRight):
            # convert uC to C (10^-6C=1uC)
            check_line = "{} {} {} {} {} {:.5f} {} {} -{:.3f}\n" \
                    .format(runNumRight[i],Q2,EbeamValRight[i],float(ChargeValRight[i])/1000000,ChargeErrRight[i], \
                            float(EffValRight[i]),EffErrRight[i],EPSVAL,thpq_right)
            # Check if the line already exists
            if check_line not in lines:
                write_to_file(f_list,check_line)
    if float(runNumLeft[0]) != 0:
        for i,thpq in enumerate(EbeamValLeft):
            check_line = "{} {} {} {} {} {:.5f} {} {} +{:.3f}\n" \
                    .format(runNumLeft[i],Q2,EbeamValLeft[i],float(ChargeValLeft[i])/1000000,ChargeErrLeft[i], \
                            float(EffValLeft[i]),EffErrLeft[i],EPSVAL,thpq_left)
            if check_line not in lines:
                write_to_file(f_list,check_line)
    if float(runNumCenter[0]) != 0:
        for i,thpq in enumerate(EbeamValCenter):
            check_line = "{} {} {} {} {} {:.5f} {} {} {:.3f}\n" \
                    .format(runNumCenter[i],Q2,EbeamValCenter[i],float(ChargeValCenter[i])/1000000,ChargeErrCenter[i], \
                            float(EffValCenter[i]),EffErrCenter[i],EPSVAL,thpq_center)
            if check_line not in lines:
                write_to_file(f_list,check_line)
                
################################################################################################################################################
