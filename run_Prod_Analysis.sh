#! /bin/bash

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

# Flag definitions (flags: h, a, o, s)
while getopts 'hdat' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_Prod_Analysis.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Plots data vs simc"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
        echo "    -h, help"
	echo "    -d, debug"	
        echo "    -a, combine data for each phi setting"	
        echo "    -t, set t-bin (!!!required for script!!!)"
	echo "        EPSILON=arg1, Q2=arg2, W=arg3, target=arg4"
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
	d) d_flag='true' ;;
	a) a_flag='true' ;;
        t) t_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# When any flag is used then the user input changes argument order
if [[ $t_flag = "true" || $a_flag = "true" ]]; then

    EPSILON=$(echo "$2" | tr '[:upper:]' '[:lower:]')
    Q2=$3
    W=$4
    echo "Epsilon must be - high - low - Case Sensitive!"
    echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5]"
    echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40]"
    if [[ -z "$2" || ! "$EPSILON" =~ high|low ]]; then # Check the 1st argument was provided and that it's one of the valid options
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
    if [[ -z "$3" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
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
    if [[ -z "$4" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
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
    
else
    
    EPSILON=$(echo "$1" | tr '[:upper:]' '[:lower:]')
    Q2=$2
    W=$3
    echo "Epsilon must be - high - low - Case Sensitive!"
    echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5]"
    echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40]"
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
    if [[ -z "$2" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
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
    if [[ -z "$3" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
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
fi

##############
# HARD CODED #
##############
NumtBins=5
NumPhiBins=16

# Define global variables for lt_analysis scripts
POL="+1" # All KaonLT is positive polarity
TMIN=0.010
TMAX=0.990
KSet=1 # Arbitrary value

# Efficiency csv file
#EffData="coin_production_Prod_efficiency_data_2022_12_05.csv"
#EffData="coin_production_Prod_efficiency_data_2022_12_30.csv"
EffData="coin_production_Prod_efficiency_data_2023_01_01.csv"

TreeNames="Uncut_Kaon_Events Cut_Kaon_Events_all_noRF Cut_Kaon_Events_prompt_noRF Cut_Kaon_Events_rand_noRF Cut_Kaon_Events_all_RF Cut_Kaon_Events_prompt_RF Cut_Kaon_Events_rand_RF"

# Function that calls python script to grab run numbers
grab_runs () {
    RunList=$1
    INPDIR="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/${RunList}"
    if [[ -e $INPDIR ]]; then
	cd "${LTANAPATH}/scripts"
	RunNumArr=$(python3 getRunNumbers.py $INPDIR)
	echo $RunNumArr
    else
	exit
    fi
}

echo
echo "---------------------------------------------------------"
echo
echo "Beginning analysis for Q2=${Q2}, W=${W}, ${EPSILON} setting..."
echo
echo "                       Number of t bins: ${NumtBins}"
echo "                       Range of t: ${TMIN} - ${TMAX}"
echo "                       Number of Phi bins: ${NumPhiBins}"
echo
echo "---------------------------------------------------------"
echo

data_right=()
data_left=()
data_center=()
dummy_right=()
dummy_left=()
dummy_center=()
# Get run numbers for left, right, and, center settings
declare -a PHI=("RIGHT" "LEFT" "CENTER")
for i in "${PHI[@]}"
do
    if [[ $d_flag = "true" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right="Prod_Test"
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs Prod_Test )"             # RIGHT, Q2=5p5, W=3p02
	    echo "Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left="Prod_Test"
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs Prod_Test )"		 # LEFT, Q2=5p5, W=3p02
	    echo "Run Numbers: [${data_left[@]}]"
	    echo	    
	elif [[ $i = "CENTER" ]]; then
	    file_center="Prod_Test"
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs Prod_Test )"		 # CENTER, Q2=5p5, W=3p02
	    echo "Run Numbers: [${data_center[@]}]"
	    echo
	fi
	KIN="Prod_Test"	
    fi

    if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    # Define run list based off kinematics selected
	    file_right_dummy="Q5p5W3p02right_${EPSILON}e_dummy"
	    file_right="Q5p5W3p02right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    # Converts python output to bash array
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=5p5, W=3p02
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"             # RIGHT, Q2=5p5, W=3p02
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q5p5W3p02left_${EPSILON}e_dummy"
	    file_left="Q5p5W3p02left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=5p5, W=3p02
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=5p5, W=3p02
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q5p5W3p02center_${EPSILON}e_dummy"
	    file_center="Q5p5W3p02center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=5p5, W=3p02
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=5p5, W=3p02
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.1838
	else
	    EPSVAL=0.5291
	fi
	KIN_dummy="Q5p5W3p02_${EPSILON}e_dummy"
	KIN="Q5p5W3p02_${EPSILON}e"
    fi
    
    if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right_dummy="Q4p4W2p74right_${EPSILON}e_dummy"
	    file_right="Q4p4W2p74right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=4p4, W=2p74
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=4p4, W=2p74
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q4p4W2p74left_${EPSILON}e_dummy"
	    file_left="Q4p4W2p74left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=4p4, W=2p74
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=4p4, W=2p74
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q4p4W2p74center_${EPSILON}e_dummy"
	    file_center="Q4p4W2p74center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=4p4, W=2p74
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=4p4, W=2p74
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.4805
	else
	    EPSVAL=0.7148
	fi
	KIN_dummy="Q4p4W2p74_${EPSILON}e_dummy"
	KIN="Q4p4W2p74_${EPSILON}e"
    fi
    
    if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right_dummy="Q3p0W3p14right_${EPSILON}e_dummy"
	    file_right="Q3p0W3p14right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=3p14
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=3p14
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q3p0W3p14left_${EPSILON}e_dummy"
	    file_left="Q3p0W3p14left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=3p14
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=3p14
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q3p0W3p14center_${EPSILON}e_dummy"
	    file_center="Q3p0W3p14center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=3p14
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=3p14
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.3935
	else
	    EPSVAL=0.6668
	fi
	KIN_dummy="Q3p0W3p14_${EPSILON}e_dummy"
	KIN="Q3p0W3p14_${EPSILON}e"
    fi
    
    if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right_dummy="Q3p0W2p32right_${EPSILON}e_dummy"
	    file_right="Q3p0W2p32right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=2p32
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=2p32
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q3p0W2p32left_${EPSILON}e_dummy"
	    file_left="Q3p0W2p32left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=2p32
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=2p32
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q3p0W2p32center_${EPSILON}e_dummy"
	    file_center="Q3p0W2p32center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=2p32
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=2p32
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.5736
	else
	    EPSVAL=0.8791
	fi
	KIN_dummy="Q3p0W2p32_${EPSILON}e_dummy"
	KIN="Q3p0W2p32_${EPSILON}e"
    fi
    
    if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right_dummy="Q2p1W2p95right_${EPSILON}e_dummy"
	    file_right="Q2p1W2p95right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=2p1, W=2p95
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=2p1, W=2p95
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q2p1W2p95left_${EPSILON}e_dummy"
	    file_left="Q2p1W2p95left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=2p1, W=2p95
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=2p1, W=2p95
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q2p1W2p95center_${EPSILON}e_dummy"
	    file_center="Q2p1W2p95center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=2p1, W=2p95
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=2p1, W=2p95
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.2477
	else
	    EPSVAL=0.7864
	fi
	KIN_dummy="Q2p1W2p95_${EPSILON}e_dummy"
	KIN="Q2p1W2p95_${EPSILON}e"
    fi
    
    if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
	if [[ $i = "RIGHT" ]]; then
	    file_right_dummy="Q0p5W2p40right_${EPSILON}e_dummy"
	    file_right="Q0p5W2p40right_${EPSILON}e"
	    echo "Reading in run numbers for right file ${file_right_dummy}..."
	    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=0p5, W=2p40
	    echo "Dummy Run Numbers: [${dummy_right[@]}]"	    
	    echo "Reading in run numbers for right file ${file_right}..."
	    IFS=', ' read -r -a data_right <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=0p5, W=2p40
	    echo "Data Run Numbers: [${data_right[@]}]"
	    echo
	elif [[ $i = "LEFT" ]]; then
	    file_left_dummy="Q0p5W2p40left_${EPSILON}e_dummy"
	    file_left="Q0p5W2p40left_${EPSILON}e"
	    echo "Reading in run numbers for left file ${file_left_dummy}..."
	    IFS=', ' read -r -a dummy_left <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=0p5, W=2p40
	    echo "Dummy Run Numbers: [${dummy_left[@]}]"	    
	    echo "Reading in run numbers for left file ${file_left}..."
	    IFS=', ' read -r -a data_left <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=0p5, W=2p40
	    echo "Data Run Numbers: [${data_left[@]}]"
	    echo
	elif [[ $i = "CENTER" ]]; then
	    file_center_dummy="Q0p5W2p40center_${EPSILON}e_dummy"
	    file_center="Q0p5W2p40center_${EPSILON}e"
	    echo "Reading in run numbers for center file ${file_center_dummy}..."
	    IFS=', ' read -r -a dummy_center <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=0p5, W=2p40
	    echo "Dummy Run Numbers: [${dummy_center[@]}]"	    
	    echo "Reading in run numbers for center file ${file_center}..."
	    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=0p5, W=2p40
	    echo "Data Run Numbers: [${data_center[@]}]"
	    echo
	fi
	if [[ ${EPSILON} == "low" ]]; then
	    EPSVAL=0.4515
	else
	    EPSVAL=0.6979
	fi
	KIN_dummy="Q0p5W2p40_${EPSILON}e_dummy"
	KIN="Q0p5W2p40_${EPSILON}e"
    fi
    
done

# Define input and output file names
InDATAFilename="Proc_Data_${KIN}.root"
OutDATAFilename="Analysed_Data_${KIN}"
InDUMMYFilename="Proc_Dummy_${KIN_dummy}.root"
OutDUMMYFilename="Analysed_Dummy_${KIN_dummy}"
OutFullAnalysisFilename="FullAnalysis_${KIN}"

# When analysis flag is used then the analysis script (Analysed_Prod.py)
# will create a new root file per run number which are combined using hadd
if [[ $a_flag = "true" ]]; then

    # Checks that array isn't empty
    if [ ${#dummy_right[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Right.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Right.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Right.root"
	fi
	echo
	echo "Combining right dummy..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Dummy" ${TreeNames} "${OutDUMMYFilename}_Right" "${dummy_right[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Right.root" ]; then
	    for i in "${dummy_right[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Dummy.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Dummy to ${i}_Proc_Dummy..."
		    #mv ${i}_-1_Raw_Dummy.root ${i}_-1_Proc_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Dummy.root ${i}_-1_Raw_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Dummy.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi

    # Checks that array isn't empty
    if [ ${#data_right[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Right.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Right.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Right.root"
	fi
	echo
	echo "Combining right data..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" ${TreeNames} "${OutDATAFilename}_Right" "${data_right[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Right.root" ]; then
	    for i in "${data_right[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Data.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Data to ${i}_Proc_Data..."
		    #mv ${i}_-1_Raw_Data.root ${i}_-1_Proc_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Data.root ${i}_-1_Raw_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Data.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi
    
    # Checks that array isn't empty
    if [ ${#dummy_left[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Left.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Left.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Left.root"
	fi
	echo
	echo "Combining left dummy..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Dummy" ${TreeNames} "${OutDUMMYFilename}_Left" "${dummy_left[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Left.root" ]; then
	    for i in "${dummy_left[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Dummy.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Dummy to ${i}_Proc_Dummy..."
		    #mv ${i}_-1_Raw_Dummy.root ${i}_-1_Proc_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Dummy.root ${i}_-1_Raw_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Dummy.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi

    # Checks that array isn't empty
    if [ ${#data_left[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root"
	fi
	echo
	echo "Combining left data..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" ${TreeNames} "${OutDATAFilename}_Left" "${data_left[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root" ]; then
	    for i in "${data_left[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Data.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Data to ${i}_Proc_Data..."
		    #mv ${i}_-1_Raw_Data.root ${i}_-1_Proc_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Data.root ${i}_-1_Raw_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Data.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi
    
    # Checks that array isn't empty
    if [ ${#dummy_center[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Center.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Center.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Center.root"
	fi
	echo
	echo "Combining center dummy..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Dummy" ${TreeNames} "${OutDUMMYFilename}_Center" "${dummy_center[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDUMMYFilename}_Center.root" ]; then
	    for i in "${dummy_center[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Dummy.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Dummy to ${i}_Proc_Dummy..."
		    #mv ${i}_-1_Raw_Dummy.root ${i}_-1_Proc_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Dummy.root ${i}_-1_Raw_Dummy.root # <runNum>_-1_Proc_Dummy.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Dummy.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi

    # Checks that array isn't empty
    if [ ${#data_center[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Center.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Center.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Center.root"
	fi
	echo
	echo "Combining center data..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" ${TreeNames} "${OutDATAFilename}_Center" "${data_center[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Center.root" ]; then
	    for i in "${data_center[@]}"
	    do       
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Data.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Data to ${i}_Proc_Data..."
		    #mv ${i}_-1_Raw_Data.root ${i}_-1_Proc_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Data.root ${i}_-1_Raw_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Data.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi
    
    # Checks that array isn't empty
    if [ ${#data_left[@]} -ne 0 ]; then
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root" ]; then
	    echo
	    echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root exists already, deleting..."
	    rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root"
	fi
	echo
	echo "Combining left data..."
	echo
	cd "${LTANAPATH}/scripts"
	python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" ${TreeNames} "${OutDATAFilename}_Left" "${data_left[*]}"
	echo
	if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Left.root" ]; then
	    for i in "${data_left[@]}"
	    do
		if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${i}_-1_Raw_Data.root" ]; then
		    cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    echo "Renaming ${i}_Raw_Data to ${i}_Proc_Data..."
		    #mv ${i}_-1_Raw_Data.root ${i}_-1_Proc_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		    #mv ${i}_-1_Proc_Data.root ${i}_-1_Raw_Data.root # <runNum>_-1_Proc_Data.root is used in later LT_analysis
		else
		    echo "WARNING: ${i}_Raw_Data.root does not exist!"
		fi
	    done	 
	fi
	echo
    fi    
    
fi

cd "${LTANAPATH}/scripts"

# Checks that array isn't empty
if [[ ${#data_right[@]} -ne 0 ]]; then
    echo
    echo "Calculating data total effective charge right..."
    PYRIGHTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_right[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYRIGHTSTRING"
    DataChargeValRight=("${arr1[@]}")
    DataChargeErrRight=("${arr2[@]}")
    DataEffValRight=("${arr3[@]}")
    DataEffErrRight=("${arr4[@]}")
    DatapThetaValRight=("${arr5[@]}")
    DataEbeamValRight=("${arr6[@]}")
    #echo ${DataChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DataChargeSumRight=$(IFS=+; echo "$((${DataChargeValRight[*]}))") # Only works for integers
    echo "Total Charge Right: ${DataChargeSumRight} uC"
fi

# Checks that array isn't empty
if [[ ${#dummy_right[@]} -ne 0 ]]; then
    echo
    echo "Calculating dummy total effective charge right..."
    PYRIGHTSTRING=$(python3 findEffectiveCharge.py ${EffDummy} "${dummy_right[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYRIGHTSTRING"
    DummyChargeValRight=("${arr1[@]}")
    DummyChargeErrRight=("${arr2[@]}")
    DummyEffValRight=("${arr3[@]}")
    DummyEffErrRight=("${arr4[@]}")
    DummypThetaValRight=("${arr5[@]}")
    DummyEbeamValRight=("${arr6[@]}")
    #echo ${DummyChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DummyChargeSumRight=$(IFS=+; echo "$((${DummyChargeValRight[*]}))") # Only works for integers
    echo "Total Dummy Charge Right: ${DummyChargeSumRight} uC"
fi

# Checks that array isn't empty
if [[ ${#data_left[@]} -ne 0 ]]; then
    echo
    echo "Calculating data total effective charge left..."
    PYLEFTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_left[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYLEFTSTRING"
    DataChargeValLeft=("${arr1[@]}")
    DataChargeErrLeft=("${arr2[@]}")
    DataEffValLeft=("${arr3[@]}")
    DataEffErrLeft=("${arr4[@]}")
    DatapThetaValLeft=("${arr5[@]}")
    DataEbeamValLeft=("${arr6[@]}")
    #echo ${DataChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DataChargeSumLeft=$(IFS=+; echo "$((${DataChargeValLeft[*]}))") # Only works for integers
    echo "Total Charge Left: ${DataChargeSumLeft} uC"
fi

# Checks that array isn't empty
if [[ ${#dummy_left[@]} -ne 0 ]]; then
    echo
    echo "Calculating dummy total effective charge left..."
    PYLEFTSTRING=$(python3 findEffectiveCharge.py ${EffDummy} "${dummy_left[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYLEFTSTRING"
    DummyChargeValLeft=("${arr1[@]}")
    DummyChargeErrLeft=("${arr2[@]}")
    DummyEffValLeft=("${arr3[@]}")
    DummyEffErrLeft=("${arr4[@]}")
    DummypThetaValLeft=("${arr5[@]}")
    DummyEbeamValLeft=("${arr6[@]}")
    #echo ${DummyChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DummyChargeSumLeft=$(IFS=+; echo "$((${DummyChargeValLeft[*]}))") # Only works for integers
    echo "Total Dummy Charge Left: ${DummyChargeSumLeft} uC"
fi

# Checks that array isn't empty
if [[ ${#data_center[@]} -ne 0 ]]; then
    echo
    echo "Calculating data total effective charge center..."
    PYCENTERSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_center[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYCENTERSTRING"
    DataChargeValCenter=("${arr1[@]}")
    DataChargeErrCenter=("${arr2[@]}")
    DataEffValCenter=("${arr3[@]}")
    DataEffErrCenter=("${arr4[@]}")
    DatapThetaValCenter=("${arr5[@]}")
    DataEbeamValCenter=("${arr6[@]}")
    #echo ${DataChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DataChargeSumCenter=$(IFS=+; echo "$((${DataChargeValCenter[*]}))") # Only works for integers
    echo "Total Charge Center: ${DataChargeSumCenter} uC"
fi

# Checks that array isn't empty
if [[ ${#dummy_center[@]} -ne 0 ]]; then
    echo
    echo "Calculating dummy total effective charge center..."
    PYCENTERSTRING=$(python3 findEffectiveCharge.py ${EffDummy} "${dummy_center[*]}")
    arr1=()
    arr2=()
    arr3=()
    arr4=()
    arr5=()
    arr6=()
    itt=0
    while read line; do
	itt=$((itt+1))
	# split the line into an array based on space
	IFS=' ' read -ra line_array <<< "$line"
	# store the elements in the corresponding array
	eval "arr$itt=(\"\${line_array[@]}\")"
    done <<< "$PYCENTERSTRING"
    DummyChargeValCenter=("${arr1[@]}")
    DummyChargeErrCenter=("${arr2[@]}")
    DummyEffValCenter=("${arr3[@]}")
    DummyEffErrCenter=("${arr4[@]}")
    DummypThetaValCenter=("${arr5[@]}")
    DummyEbeamValCenter=("${arr6[@]}")
    #echo ${DummyChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DummyChargeSumCenter=$(IFS=+; echo "$((${DummyChargeValCenter[*]}))") # Only works for integers
    echo "Total Dummy Charge Center: ${DummyChargeSumCenter} uC"
fi

# Run the plotting script if t-flag enabled
# Checks that array isn't empty
if [[ $t_flag = "true" || $d_flag = "true" ]]; then
    echo
    echo
    echo
    echo "Finding t-bins..."
    cd "${LTANAPATH}/scripts/Prod"    
    if [ ${#data_right[@]} -eq 0 ]; then
	python3 find_tBinRange.py ${KIN} ${W} ${Q2} ${EPSVAL} ${OutDATAFilename} ${OutFullAnalysisFilename} ${TMIN} ${TMAX} ${NumtBins} ${NumPhiBins} "0" "${data_left[*]}" "${data_center[*]}" "0" "${dummy_left[*]}" "${dummy_center[*]}" "0" ${DataChargeSumLeft} ${DataChargeSumCenter} "0" ${DummyChargeSumLeft} ${DummyChargeSumCenter} "0" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" ${EffData}
    else
	python3 find_tBinRange.py ${KIN} ${W} ${Q2} ${EPSVAL} ${OutDATAFilename} ${OutFullAnalysisFilename} ${TMIN} ${TMAX} ${NumtBins} ${NumPhiBins} "${data_right[*]}" "${data_left[*]}" "${data_center[*]}" "${dummy_right[*]}" "${dummy_left[*]}" "${dummy_center[*]}" ${DataChargeSumRight} ${DataChargeSumLeft} ${DataChargeSumCenter} ${DummyChargeSumRight} ${DummyChargeSumLeft} ${DummyChargeSumCenter} "${DataEffValRight[*]}" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" ${EffData}
    fi
fi

echo
echo
echo
echo "Creating analysis lists..."
cd "${LTANAPATH}/scripts/Prod/"

# Create input for lt_analysis code
if [ ${#data_right[@]} -eq 0 ]; then
    python3 createPhysicsList.py ${Q2} ${POL} ${EPSVAL} ${TMIN} ${TMAX} ${NumtBins} ${KSet} "0" "${data_left[*]}" "${data_center[*]}" "0" "${DatapThetaValLeft[*]}" "${DatapThetaValCenter[*]}" "0" "${DataEbeamValLeft[*]}" "${DataEbeamValCenter[*]}" "0" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" "0" "${DataEffErrLeft[*]}" "${DataEffErrCenter[*]}" "0" "${DataChargeValLeft[*]}" "${DataChargeValCenter[*]}" "0" "${DataChargeErrLeft[*]}" "${DataChargeErrCenter[*]}" ${KIN}
else
    python3 createPhysicsList.py ${Q2} ${POL} ${EPSVAL} ${TMIN} ${TMAX} ${NumtBins} ${KSet} "${data_right[*]}" "${data_left[*]}" "${data_center[*]}" "${DatapThetaValRight[*]}" "${DatapThetaValLeft[*]}" "${DatapThetaValCenter[*]}" "${DataEbeamValRight[*]}" "${DataEbeamValLeft[*]}" "${DataEbeamValCenter[*]}" "${DataEffValRight[*]}" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" "${DataEffErrRight[*]}" "${DataEffErrLeft[*]}" "${DataEffErrCenter[*]}" "${DataChargeValRight[*]}" "${DataChargeValLeft[*]}" "${DataChargeValCenter[*]}" "${DataChargeErrRight[*]}" "${DataChargeErrLeft[*]}" "${DataChargeErrCenter[*]}" ${KIN}
fi

if [[ $t_flag = "true" || $d_flag = "true" ]]; then
    cd "${LTANAPATH}"
    evince "OUTPUT/Analysis/${ANATYPE}LT/${OutFullAnalysisFilename}.pdf"
fi

echo
echo
echo
echo "Script Complete!"
