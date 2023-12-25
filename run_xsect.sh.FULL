#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2023-10-08 13:08:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ ${HOSTNAME} = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f16`
LTANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f17`

# Flag definitions (flags: h)
while getopts 'h' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_xsect.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Get the average kinematics per t,phi bin."
	echo "             Then find the unsep xsect for given kinematic."
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
        echo "    -h, help"
	echo "     Q2=arg1, W=arg2"
	echo
	echo " Avaliable Kinematics..."	
	echo "                      Q2=5p5, W=3p02"
	echo "                      Q2=4p4, W=2p74"
	echo "                      Q2=3p0, W=3p14"
	echo "                      Q2=3p0, W=2p32"
	echo "                      Q2=2p1, W=2p95"
	echo "                      Q2=0p5, W=2p40"
        exit 0
        ;;
	a) a_flag='true' ;;
        t) t_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

Q2=$1
W=$2
echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5]"
echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40]"
if [[ -z "$1" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid Q2..."
    while true; do
	echo ""
	read -p "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5] - or press ctrl-c to exit : " Q2
	case $Q2 in
	    '');; # If blank, prompt again
	    '5p5'|'4p4'|'3p0'|'2p1'|'0p5') break;; # If a valid option, break the loop and continue
	esac
    done
fi
if [[ -z "$2" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid W..."
    while true; do
	echo ""
	read -p "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40] - or press ctrl-c to exit : " W
	case $W in
	    '');; # If blank, prompt again
	    '3p02'|'2p74'|'3p14'|'2p32'|'2p95'|'2p40') break;; # If a valid option, break the loop and continue
	esac
    done
fi

ParticleType=$3
POL=$4

NumtBins=$5
NumPhiBins=$6

##############
# HARD CODED #
##############

if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
    LOEPS=0.1838
    HIEPS=0.5291
fi

if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
    LOEPS=0.4805
    HIEPS=0.7148
fi

if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
    LOEPS=0.3935
    HIEPS=0.6668
fi

if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
    LOEPS=0.5736
    HIEPS=0.8791
fi

if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
    LOEPS=0.2477
    HIEPS=0.7864
fi

if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
    LOEPS=0.4515
    HIEPS=0.6979
fi

##############
##############
##############

echo
echo "---------------------------------------------------------"
echo
echo "Finding unesp xsect for Q2=${Q2}, W=${W}  setting..."
echo
echo "                        LOEPS = ${LOEPS}"
echo "                        HIEPS = ${HIEPS}"
echo
echo "---------------------------------------------------------"
echo

# Replace p with '.'
Q2=${Q2//p/.}

cd "${LTANAPATH}/src/"
echo
echo "Compiling average_kinematics.f..."
eval "gfortran -o average_kinematics average_kinematics.f"
# Check the exit status of the Fortran script
if [ $? -ne 0 ]; then
    echo
    echo
    echo "2 ERROR: Fortran script failed!"
    echo "       See error above..."
    exit 1
fi
echo
echo "Running average_kinematics..."
./average_kinematics.expect ${ParticleType} ${POL} ${Q2} ${LOEPS} ${HIEPS}

echo
echo "Compiling calc_xsect.f..."
eval "gfortran -o calc_xsect calc_xsect.f"
# Check the exit status of the Fortran script
if [ $? -ne 0 ]; then
    echo
    echo
    echo "2 ERROR: Fortran script failed!"
    echo "       See error above..."
    exit 1
fi
echo
echo "Running calc_xsect..."
./calc_xsect.expect ${ParticleType} ${POL} ${Q2} ${LOEPS} ${HIEPS}

# Replace p with '.'
Q2=${Q2//./p}

KIN="Q${Q2}W${W}"

# Define input and output file names
OutUnsepxsectsFilename="${ParticleType}_xsects_${KIN}"

cd "${LTANAPATH}/src/plotting/"
python3 plot_xsects.py ${ParticleType} ${POL} ${Q2} ${W} ${LOEPS} ${HIEPS} ${NumtBins} ${NumPhiBins} ${KIN} ${OutUnsepxsectsFilename}
# Check the exit status of the Python script
if [ $? -ne 0 ]; then
    echo
    echo
    echo "2 ERROR: Python script failed!"
    echo "       See error above..."
    exit 1
fi
