#! /bin/bash

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ ${HOSTNAME} = *"farm"* ]]; then
    PATHFILE_INFO=`python3 $replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
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

# Flag definitions (flags: h, c, a, s)
while getopts 'h' flag; do
    case "${flag}" in
        h) 
        echo "----------------------------------------------------------"
        echo "./set_ProdInput.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Runs simc and updates simc input files"
        echo "----------------------------------------------------------"
        echo
        echo "The following flags can be called for the production analysis..."
	echo "    EPSILON=arg1, Q2=arg2, W=arg3, ParticleType=arg4"
        echo "    -h, help"
	echo
	echo " Avaliable Kinematics..."	
	echo "                      Q2=5p5, W=3p02"
	echo "                      Q2=4p4, W=2p74"
	echo "                      Q2=3p0, W=3p14"
	echo "                      Q2=3p0, W=2p32"
	echo "                      Q2=2p1, W=2p95"
	echo "                      Q2=0p5, W=2p40"
	echo "                      Q2=0p4, W=2p20 (Q2 = 0.38, pion only)"		
        exit 0
        ;;
        *) print_usage
        exit 1 ;;
    esac
done

EPSILON=$1
Q2=$2
W=$3
ParticleType=$4
echo "Epsilon must be - high - low - Case Sensitive!"
echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4]"
echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20]"
if [[ -z "$1" || ! "$EPSILON" =~ high|low ]]; then # Check the 1st argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid epsilon..."
    while true; do
	echo ""
	read -p "Epsilon must be - high - low - Case Sensitive! - or press ctrl-c to exit : " EPSILON
	case $EPSILON in
	    '');; # If blank, prompt again
	    'high'|'low') break;; # If a valid option, break the loop and continue
	esac
    done
fi
if [[ -z "$2" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5|0p4 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid Q2..."
    while true; do
	echo ""
	read -p "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4] - or press ctrl-c to exit : " Q2
	case $Q2 in
	    '');; # If blank, prompt again
	    '5p5'|'4p4'|'3p0'|'2p1'|'0p5'|'0p4') break;; # If a valid option, break the loop and continue
	esac
    done
fi
if [[ -z "$3" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40|2p20 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid W..."
    while true; do
	echo ""
	read -p "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20] - or press ctrl-c to exit : " W
	case $W in
	    '');; # If blank, prompt again
	    '3p02'|'2p74'|'3p14'|'2p32'|'2p95'|'2p40'|'2p20') break;; # If a valid option, break the loop and continue
	esac
    done
fi
if [[ -z "$4" || ! "$ParticleType" =~ kaon|pion|proton ]]; then # Check the 3rd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid particle type..."
    while true; do
	echo ""
	read -p "Particle type must be one of - [kaon - pion - proton] - or press ctrl-c to exit : " ParticleType
	case $ParticleType in
	    '');; # If blank, prompt again
	    'kaon'|'pion'|'proton') break;; # If a valid option, break the loop and continue
	esac
    done
fi

InputSIMC_right="Prod_Coin_Q${Q2}W${W}right_${EPSILON}e"
InputSIMC_left="Prod_Coin_Q${Q2}W${W}left_${EPSILON}e"
InputSIMC_center="Prod_Coin_Q${Q2}W${W}center_${EPSILON}e"

cd ${LTANAPATH}/src/setup

./run_simc_tree "${InputSIMC_center}" "production"
