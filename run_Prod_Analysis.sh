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

# Flag definitions (flags: h, a, b, p)
while getopts 'hcip' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_Prod_Analysis.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Plots data vs simc"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
	echo "    Q2=arg1, W=arg2"
        echo "    -h, help"
        echo "    -c, combine all runs for each setting"
	echo "    -i, iterate SIMC to find proper weight"
	echo "    -p, specify particle type (kaon, pion, or proton). Otherwise runs for all."
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
	c) c_flag='true' ;;
	i) i_flag='true' ;;
	p) p_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# - `date` is a command that prints or sets the system date and time.
# - `+%H` extracts the hour in 24-hour format.
# - `+%M` extracts the minute.
# - `+%S` extracts the second.
# - `+%Y` extracts the year.
# - `+%B` extracts the full month name.
# - `%d` extracts the day of the month.
formatted_date=$(date +H%HM%MS%S_%Y%B%d)

##############
# HARD CODED #
##############

#DEBUG="False" # Flag for no plot splash
DEBUG="True" # Flag for plot splash

if [[ $p_flag != "true" ]]; then
    ParticleType="kaon"
    #ParticleType="pion"
else
    echo "Please select a particle type..."
    while true; do
	echo ""
	read -p "Particle type must be one of - [kaon - pion - proton] - or press ctrl-c to exit : " ParticleType
	case $ParticleType in
	    '');; # If blank, prompt again
	    'kaon'|'pion'|'proton') break;; # If a valid option, break the loop and continue
	esac
    done
fi

# Efficiency csv file
#EffData="coin_production_Prod_efficiency_data_2022_12_05.csv"
#EffData="coin_production_Prod_efficiency_data_2022_12_30.csv"
#EffData="coin_production_Prod_efficiency_data_2023_01_01.csv"
#EffData="coin_production_Prod_efficiency_data_2023_12_18.csv"
#EffData="coin_production_Prod_efficiency_data_2024_01_11.csv"
EffData="coin_production_Prod_efficiency_data_2024_01_14.csv"
##############
##############
##############

declare -a EPS=("low" "high")
for j in "${EPS[@]}"
do
    # When any flag is used then the user input changes argument order
    if [[ $i_flag = "true" || $p_flag = "true" || $c_flag = "true" ]]; then

	EPSILON=$j
	Q2=$2
	W=$3
	echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5]"
	echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40]"
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

    else

	EPSILON=$j
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
    fi

    # When analysis flag is used then the analysis script (Analysed_Prod.py)
    # will create a new root file per run number which are combined using hadd
    if [[ $c_flag = "true" ]]; then
	
	# Function that calls python script to grab run numbers
	grab_runs () {
	    RunList=$1
	    # Location of list of run lists
	    INPDIR="${UTILPATH}/run_list/${ANATYPE}LT/${RunList}"
	    if [[ -e $INPDIR ]]; then
		cd "${LTANAPATH}/src/setup"
		RunNumArr=$(python3 getRunNumbers.py $INPDIR)
		echo $RunNumArr
	    else
		exit
	    fi
	}
	
	echo
	echo "---------------------------------------------------------"
	echo
	echo "Combining data for Q2=${Q2}, W=${W}, ${EPSILON} setting..."
	echo
	echo "---------------------------------------------------------"
	echo
	
	data_right_tmp=()
	data_left_tmp=()
	data_center_tmp=()
	dummy_right_tmp=()
	dummy_left_tmp=()
	dummy_center_tmp=()
	# Get run numbers for left, right, and, center settings
	declare -a PHI=("RIGHT" "LEFT" "CENTER")
	for i in "${PHI[@]}"
	do
	    ##############
	    # HARD CODED #
	    ##############
	    if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    # Define run list based off kinematics selected
		    file_right_dummy="Q5p5W3p02right_${EPSILON}e_dummy"
		    file_right="Q5p5W3p02right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    # Converts python output to bash array
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"             # RIGHT, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q5p5W3p02left_${EPSILON}e_dummy"
		    file_left="Q5p5W3p02left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q5p5W3p02center_${EPSILON}e_dummy"
		    file_center="Q5p5W3p02center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.1838
		HIEPS=0.5291
		KIN="Q5p5W3p02_${EPSILON}e"
	    fi

	    if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q4p4W2p74right_${EPSILON}e_dummy"
		    file_right="Q4p4W2p74right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q4p4W2p74left_${EPSILON}e_dummy"
		    file_left="Q4p4W2p74left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q4p4W2p74center_${EPSILON}e_dummy"
		    file_center="Q4p4W2p74center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.4805
		HIEPS=0.7148
		KIN="Q4p4W2p74_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W3p14right_${EPSILON}e_dummy"
		    file_right="Q3p0W3p14right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q3p0W3p14left_${EPSILON}e_dummy"
		    file_left="Q3p0W3p14left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q3p0W3p14center_${EPSILON}e_dummy"
		    file_center="Q3p0W3p14center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.3935
		HIEPS=0.6668
		KIN="Q3p0W3p14_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W2p32right_${EPSILON}e_dummy"
		    file_right="Q3p0W2p32right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q3p0W2p32left_${EPSILON}e_dummy"
		    file_left="Q3p0W2p32left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q3p0W2p32center_${EPSILON}e_dummy"
		    file_center="Q3p0W2p32center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.5736
		HIEPS=0.8791
		KIN="Q3p0W2p32_${EPSILON}e"
	    fi

	    if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q2p1W2p95right_${EPSILON}e_dummy"
		    file_right="Q2p1W2p95right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q2p1W2p95left_${EPSILON}e_dummy"
		    file_left="Q2p1W2p95left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q2p1W2p95center_${EPSILON}e_dummy"
		    file_center="Q2p1W2p95center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.2477
		HIEPS=0.7864
		KIN="Q2p1W2p95_${EPSILON}e"
	    fi

	    if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q0p5W2p40right_${EPSILON}e_dummy"
		    file_right="Q0p5W2p40right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q0p5W2p40left_${EPSILON}e_dummy"
		    file_left="Q0p5W2p40left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q0p5W2p40center_${EPSILON}e_dummy"
		    file_center="Q0p5W2p40center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.4515
		HIEPS=0.6979
		KIN="Q0p5W2p40_${EPSILON}e"
	    fi
	    ##############
	    ##############
	    ##############
	done

	# Define input and output file names
	OutDATAFilename="Analysed_Data_${KIN}"
	OutDUMMYFilename="Analysed_Dummy_${KIN}"
	OutFullAnalysisFilename="FullAnalysis_${KIN}"
	
	if [[ $ParticleType = "kaon" ]]; then
	    TreeNames="Uncut_Kaon_Events Cut_Kaon_Events_all_noRF Cut_Kaon_Events_prompt_noRF Cut_Kaon_Events_rand_noRF Cut_Kaon_Events_all_RF Cut_Kaon_Events_prompt_RF Cut_Kaon_Events_rand_RF"
	fi
	if [[ $ParticleType = "pion" ]]; then
	    TreeNames="Uncut_Pion_Events Cut_Pion_Events_all_noRF Cut_Pion_Events_prompt_noRF Cut_Pion_Events_rand_noRF Cut_Pion_Events_all_RF Cut_Pion_Events_prompt_RF Cut_Pion_Events_rand_RF"
	fi
	if [[ $ParticleType = "proton" ]]; then
	    TreeNames="Uncut_Proton_Events Cut_Proton_Events_all_noRF Cut_Proton_Events_prompt_noRF Cut_Proton_Events_rand_noRF Cut_Proton_Events_all_RF Cut_Proton_Events_prompt_RF Cut_Proton_Events_rand_RF"
	fi

	# Checks that array isn't empty
	if [ ${#dummy_right_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Right.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Right.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Right.root"
	    fi
	    echo
	    echo "Combining right ${ParticleType} dummy..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDUMMYFilename}_Right" "${dummy_right_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Right.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Right.root" ]; then
		for i in "${dummy_right_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Right.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${#data_right_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Right.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Right.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Right.root"
	    fi
	    echo
	    echo "Combining right ${ParticleType} data..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDATAFilename}_Right" "${data_right_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Right.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Right.root" ]; then
		for i in "${data_right_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Right.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${#dummy_left_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Left.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Left.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Left.root"
	    fi
	    echo
	    echo "Combining left ${ParticleType} dummy..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDUMMYFilename}_Left" "${dummy_left_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Left.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Left.root" ]; then
		for i in "${dummy_left_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Left.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${#data_left_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Left.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Left.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Left.root"
	    fi
	    echo
	    echo "Combining left ${ParticleType} data..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDATAFilename}_Left" "${data_left_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Left.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Left.root" ]; then
		for i in "${data_left_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Left.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${#dummy_center_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Center.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Center.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Center.root"
	    fi
	    echo
	    echo "Combining center ${ParticleType} dummy..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDUMMYFilename}_Center" "${dummy_center_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Center.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDUMMYFilename}_Center.root" ]; then
		for i in "${dummy_center_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDUMMYFilename}_Center.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${#data_center_tmp[@]} -ne 0 ]; then
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Center.root" ]; then
		echo
		echo "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Center.root exists already, deleting..."
		rm -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${OutDATAFilename}_Center.root"
	    fi
	    echo
	    echo "Combining center ${ParticleType} data..."
	    echo
	    cd "${LTANAPATH}/src/setup"
	    python3 mergeRootFiles.py "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "_-1_Raw_Data" "${TreeNames}" "${ParticleType}_${OutDATAFilename}_Center" "${data_center_tmp[*]}" "${ParticleType}" "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Center.err"
	    echo
	    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${OutDATAFilename}_Center.root" ]; then
		for i in "${data_center_tmp[@]}"
		do       
		    if [ -f "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/${ParticleType}_${i}_-1_Raw_Data.root" ]; then
			cd "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT"
		    else
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >> "${LTANAPATH}/log/${ParticleType}_${OutDATAFilename}_Center.err"
			echo "WARNING: ${ParticleType}_${i}_Raw_Data.root does not exist!" >&2 # Redirect to stderr
		    fi
		done	 
	    fi
	    echo
	fi
    fi
done

##############
# HARD CODED #
##############

# tbins should not exceed 8 (major drop in statistics)
# TMIN should not equal zero (unless calc_xsect.f is adapted)
# Make sure 3 sig figs (no more)

if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
    # Q2=2p1, W=2p95
    NumtBins=3
    NumPhiBins=10
    TMIN=0.100
    TMAX=0.600
elif [[ $Q2 = "3p0" && $W = "2p32" ]]; then
    # Q2=3p0, W=2p32
    NumtBins=5
    NumPhiBins=10
    TMIN=0.100
    TMAX=0.990
elif [[ $Q2 = "3p0" && $W = "3p14" ]]; then
    # Q2=3p0, W=3p14
    NumtBins=5
    NumPhiBins=16
    TMIN=0.100
    TMAX=0.850
elif [[ $Q2 = "4p4" && $W = "2p74" ]]; then
    # Q2=4p4, W=2p74
    NumtBins=2
    NumPhiBins=8
    TMIN=0.001
    TMAX=1.300    
elif [[ $Q2 = "5p5" && $W = "3p02" ]]; then
    # Q2=5p5, W=3p02
    NumtBins=2
    NumPhiBins=8
    TMIN=0.001
    TMAX=1.300
else
    # For testing
    NumtBins=1
    NumPhiBins=1
    #NumtBins=2
    #NumPhiBins=8
    TMIN=0.001
    TMAX=0.990    
fi

# Define global variables for lt_analysis scripts
POL="+1" # All KaonLT is positive polarity

##############
##############
##############

if [[ $i_flag != "true" ]]; then
    # Need to rerun loop separately so that the combined files for high and low epsilon exists for diamond cut script
    declare -a EPS=("low" "high")
    for j in "${EPS[@]}"
    do

	# Redefine epsilon based on loop
	EPSILON=$j

	# Function that calls python script to grab run numbers
	grab_runs () {
	    RunList=$1
	    # Location of list of run lists
	    INPDIR="${UTILPATH}/run_list/${ANATYPE}LT/${RunList}"
	    if [[ -e $INPDIR ]]; then
		cd "${LTANAPATH}/src/setup"
		RunNumArr=$(python3 getRunNumbers.py $INPDIR)
		echo $RunNumArr
	    else
		exit
	    fi
	}

	echo
	echo "---------------------------------------------------------"
	echo
	echo "${formatted_date}"
	echo
	echo "Beginning analysis for Q2=${Q2}, W=${W}, ${EPSILON} setting..."
	echo
	echo "                       Number of t bins: ${NumtBins}"
	echo "                       Range of t: ${TMIN} - ${TMAX}"
	echo "                       Number of Phi bins: ${NumPhiBins}"
	echo
	echo "---------------------------------------------------------"
	echo

	data_right_tmp=()
	data_left_tmp=()
	data_center_tmp=()
	dummy_right_tmp=()
	dummy_left_tmp=()
	dummy_center_tmp=()
	# Get run numbers for left, right, and, center settings
	declare -a PHI=("RIGHT" "LEFT" "CENTER")
	for i in "${PHI[@]}"
	do
	    ##############
	    # HARD CODED #
	    ##############	    
	    if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    # Define run list based off kinematics selected
		    file_right_dummy="Q5p5W3p02right_${EPSILON}e_dummy"
		    file_right="Q5p5W3p02right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    # Converts python output to bash array
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"             # RIGHT, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q5p5W3p02left_${EPSILON}e_dummy"
		    file_left="Q5p5W3p02left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q5p5W3p02center_${EPSILON}e_dummy"
		    file_center="Q5p5W3p02center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.1838
		HIEPS=0.5291
		KIN="Q5p5W3p02_${EPSILON}e"
	    fi

	    if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q4p4W2p74right_${EPSILON}e_dummy"
		    file_right="Q4p4W2p74right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q4p4W2p74left_${EPSILON}e_dummy"
		    file_left="Q4p4W2p74left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q4p4W2p74center_${EPSILON}e_dummy"
		    file_center="Q4p4W2p74center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.4805
		HIEPS=0.7148
		KIN="Q4p4W2p74_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W3p14right_${EPSILON}e_dummy"
		    file_right="Q3p0W3p14right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q3p0W3p14left_${EPSILON}e_dummy"
		    file_left="Q3p0W3p14left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q3p0W3p14center_${EPSILON}e_dummy"
		    file_center="Q3p0W3p14center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.3935
		HIEPS=0.6668
		KIN="Q3p0W3p14_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W2p32right_${EPSILON}e_dummy"
		    file_right="Q3p0W2p32right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q3p0W2p32left_${EPSILON}e_dummy"
		    file_left="Q3p0W2p32left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q3p0W2p32center_${EPSILON}e_dummy"
		    file_center="Q3p0W2p32center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.5736
		HIEPS=0.8791
		KIN="Q3p0W2p32_${EPSILON}e"
	    fi

	    if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q2p1W2p95right_${EPSILON}e_dummy"
		    file_right="Q2p1W2p95right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q2p1W2p95left_${EPSILON}e_dummy"
		    file_left="Q2p1W2p95left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q2p1W2p95center_${EPSILON}e_dummy"
		    file_center="Q2p1W2p95center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.2477
		HIEPS=0.7864
		KIN="Q2p1W2p95_${EPSILON}e"
	    fi

	    if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q0p5W2p40right_${EPSILON}e_dummy"
		    file_right="Q0p5W2p40right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right_tmp <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_right_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for right file ${file_right}..."
		    IFS=', ' read -r -a data_right_tmp <<< "$( grab_runs ${file_right} )"		 # RIGHT, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_right_tmp[@]}]"
		    echo
		elif [[ $i = "LEFT" ]]; then
		    file_left_dummy="Q0p5W2p40left_${EPSILON}e_dummy"
		    file_left="Q0p5W2p40left_${EPSILON}e"
		    echo "Reading in run numbers for left file ${file_left_dummy}..."
		    IFS=', ' read -r -a dummy_left_tmp <<< "$( grab_runs ${file_left_dummy} )"             # LEFT, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_left_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for left file ${file_left}..."
		    IFS=', ' read -r -a data_left_tmp <<< "$( grab_runs ${file_left} )"		 # LEFT, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_left_tmp[@]}]"
		    echo
		elif [[ $i = "CENTER" ]]; then
		    file_center_dummy="Q0p5W2p40center_${EPSILON}e_dummy"
		    file_center="Q0p5W2p40center_${EPSILON}e"
		    echo "Reading in run numbers for center file ${file_center_dummy}..."
		    IFS=', ' read -r -a dummy_center_tmp <<< "$( grab_runs ${file_center_dummy} )"             # CENTER, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_center_tmp[@]}]"
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center_tmp <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_center_tmp[@]}]"
		    echo
		fi
		LOEPS=0.4515
		HIEPS=0.6979
		KIN="Q0p5W2p40_${EPSILON}e"
	    fi
	    ##############
	    ##############
	    ##############
	done

	# Define input and output file names
	OutDATAFilename="Analysed_Data_${KIN}"
	OutDUMMYFilename="Analysed_Dummy_${KIN}"
	OutFullAnalysisFilename="FullAnalysis_${KIN}"

	cd "${LTANAPATH}/src/setup"
	
	# Checks that array isn't empty
	if [[ ${#data_right_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating data total effective charge right..."
	    PYRIGHTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_right_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYRIGHTSTRING"
	    DataEffChargeValRight=("${arr1[@]}")
	    DataEffChargeErrRight=("${arr2[@]}")
	    DataEffValRight=("${arr3[@]}")
	    DataEffErrRight=("${arr4[@]}")
	    DatapThetaValRight=("${arr5[@]}")
	    DataEbeamValRight=("${arr6[@]}")
	    data_right=("${arr7[@]}")
	    if [ "${#data_right_tmp[@]}" -ne "${#data_right[@]}" ]; then
		echo "Removing bad right data runs..."
		for run in "${data_right_tmp[@]}"; do
		    if [[ ! " ${data_right[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DataEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYRIGHTTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DataEffChargeValRight[*]}" "${DataEffChargeErrRight[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYRIGHTTOTEFFCHARGE"
	    TotDataEffChargeValRight=("${arr1[@]}")
	    TotDataEffChargeErrRight=("${arr2[@]}")
	    echo "Total Effective Charge Right: ${TotDataEffChargeValRight} mC"
	    echo "Total Effective Charge Right Error: ${TotDataEffChargeErrRight}"
	    echo "Run numbers: [${data_right[@]}]"
	    #echo "Effective Charge per Run: [${DataEffChargeValRight[@]}]"
	    #echo "Effective Charge Error per Run: [${DataEffChargeErrRight[@]}]"
	    #echo "Efficiency per Run: [${DataEffValRight[@]}]"
	    #echo "Efficiency Error per Run: [${DataEffErrRight[@]}]"
	    #echo "Theta per Run: [${DatapThetaValRight[@]}]"
	    #echo "Beam Energy per Run: [${DataEbeamValRight[@]}]"
	fi

	# Checks that array isn't empty
	if [[ ${#dummy_right_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating dummy total effective charge right..."
	    PYRIGHTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${dummy_right_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYRIGHTSTRING"
	    DummyEffChargeValRight=("${arr1[@]}")
	    DummyEffChargeErrRight=("${arr2[@]}")
	    DummyEffValRight=("${arr3[@]}")
	    DummyEffErrRight=("${arr4[@]}")
	    DummypThetaValRight=("${arr5[@]}")
	    DummyEbeamValRight=("${arr6[@]}")
	    dummy_right=("${arr7[@]}")
	    if [ "${#dummy_right_tmp[@]}" -ne "${#dummy_right[@]}" ]; then
		echo "Removing bad right dummy runs..."
		for run in "${dummy_right_tmp[@]}"; do
		    if [[ ! " ${dummy_right[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DummyEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYRIGHTTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DummyEffChargeValRight[*]}" "${DummyEffChargeErrRight[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYRIGHTTOTEFFCHARGE"
	    TotDummyEffChargeValRight=("${arr1[@]}")
	    TotDummyEffChargeErrRight=("${arr2[@]}")
	    echo "Total Effective Charge Right: ${TotDummyEffChargeValRight} mC"
	    echo "Total Effective Charge Right Error: ${TotDummyEffChargeErrRight}"
	    echo "Run numbers: [${dummy_right[@]}]"
	    #echo "Effective Charge per Run: [${DummyEffChargeValRight[@]}]"
	    #echo "Effective Charge Error per Run: [${DummyEffChargeErrRight[@]}]"
	    #echo "Efficiency per Run: [${DummyEffValRight[@]}]"
	    #echo "Efficiency Error per Run: [${DummyEffErrRight[@]}]"
	    #echo "Theta per Run: [${DummypThetaValRight[@]}]"
	    #echo "Beam Energy per Run: [${DummyEbeamValRight[@]}]"
	fi

	# Checks that array isn't empty
	if [[ ${#data_left_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating data total effective charge left..."
	    PYLEFTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_left_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYLEFTSTRING"
	    DataEffChargeValLeft=("${arr1[@]}")
	    DataEffChargeErrLeft=("${arr2[@]}")
	    DataEffValLeft=("${arr3[@]}")
	    DataEffErrLeft=("${arr4[@]}")
	    DatapThetaValLeft=("${arr5[@]}")
	    DataEbeamValLeft=("${arr6[@]}")
	    data_left=("${arr7[@]}")
	    if [ "${#data_left_tmp[@]}" -ne "${#data_left[@]}" ]; then
		echo "Removing bad left data runs..."
		for run in "${data_left_tmp[@]}"; do
		    if [[ ! " ${data_left[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DataEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYLEFTTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DataEffChargeValLeft[*]}" "${DataEffChargeErrLeft[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYLEFTTOTEFFCHARGE"
	    TotDataEffChargeValLeft=("${arr1[@]}")
	    TotDataEffChargeErrLeft=("${arr2[@]}")
	    echo "Total Effective Charge Left: ${TotDataEffChargeValLeft} mC"
	    echo "Total Effective Charge Left Error: ${TotDataEffChargeErrLeft}"
	    echo "Run numbers: [${data_left[@]}]"
	    #echo "Effective Charge per Run: [${DataEffChargeValLeft[@]}]"
	    #echo "Effective Charge Error per Run: [${DataEffChargeErrLeft[@]}]"
	    #echo "Efficiency per Run: [${DataEffValLeft[@]}]"
	    #echo "Efficiency Error per Run: [${DataEffErrLeft[@]}]"
	    #echo "Theta per Run: [${DatapThetaValLeft[@]}]"
	    #echo "Beam Energy per Run: [${DataEbeamValLeft[@]}]"
	fi

	# Checks that array isn't empty
	if [[ ${#dummy_left_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating dummy total effective charge left..."
	    PYLEFTSTRING=$(python3 findEffectiveCharge.py ${EffData} "${dummy_left_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYLEFTSTRING"
	    DummyEffChargeValLeft=("${arr1[@]}")
	    DummyEffChargeErrLeft=("${arr2[@]}")
	    DummyEffValLeft=("${arr3[@]}")
	    DummyEffErrLeft=("${arr4[@]}")
	    DummypThetaValLeft=("${arr5[@]}")
	    DummyEbeamValLeft=("${arr6[@]}")
	    dummy_left=("${arr7[@]}")
	    if [ "${#dummy_left_tmp[@]}" -ne "${#dummy_left[@]}" ]; then
		echo "Removing bad left dummy runs..."
		for run in "${dummy_left_tmp[@]}"; do
		    if [[ ! " ${dummy_left[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DummyEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYLEFTTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DummyEffChargeValLeft[*]}" "${DummyEffChargeErrLeft[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYLEFTTOTEFFCHARGE"
	    TotDummyEffChargeValLeft=("${arr1[@]}")
	    TotDummyEffChargeErrLeft=("${arr2[@]}")
	    echo "Total Effective Charge Left: ${TotDummyEffChargeValLeft} mC"
	    echo "Total Effective Charge Left Error: ${TotDummyEffChargeErrLeft}"
	    echo "Run numbers: [${dummy_left[@]}]"
	    #echo "Effective Charge per Run: [${DummyEffChargeValLeft[@]}]"
	    #echo "Effective Charge Error per Run: [${DummyEffChargeErrLeft[@]}]"
	    #echo "Efficiency per Run: [${DummyEffValLeft[@]}]"
	    #echo "Efficiency Error per Run: [${DummyEffErrLeft[@]}]"
	    #echo "Theta per Run: [${DummypThetaValLeft[@]}]"
	    #echo "Beam Energy per Run: [${DummyEbeamValLeft[@]}]"
	fi

	# Checks that array isn't empty
	if [[ ${#data_center_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating data total effective charge center..."
	    PYCENTERSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data_center_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYCENTERSTRING"
	    DataEffChargeValCenter=("${arr1[@]}")
	    DataEffChargeErrCenter=("${arr2[@]}")
	    DataEffValCenter=("${arr3[@]}")
	    DataEffErrCenter=("${arr4[@]}")
	    DatapThetaValCenter=("${arr5[@]}")
	    DataEbeamValCenter=("${arr6[@]}")
	    data_center=("${arr7[@]}")
	    if [ "${#data_center_tmp[@]}" -ne "${#data_center[@]}" ]; then
		echo "Removing bad center data runs..."
		for run in "${data_center_tmp[@]}"; do
		    if [[ ! " ${data_center[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DataEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYCENTERTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DataEffChargeValCenter[*]}" "${DataEffChargeErrCenter[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYCENTERTOTEFFCHARGE"
	    TotDataEffChargeValCenter=("${arr1[@]}")
	    TotDataEffChargeErrCenter=("${arr2[@]}")
	    echo "Total Effective Charge Center: ${TotDataEffChargeValCenter} mC"
	    echo "Total Effective Charge Center Error: ${TotDataEffChargeErrCenter}"
	    echo "Run numbers: [${data_center[@]}]"
	    #echo "Effective Charge per Run: [${DataEffChargeValCenter[@]}]"
	    #echo "Effective Charge Error per Run: [${DataEffChargeErrCenter[@]}]"
	    #echo "Efficiency per Run: [${DataEffValCenter[@]}]"
	    #echo "Efficiency Error per Run: [${DataEffErrCenter[@]}]"
	    #echo "Theta per Run: [${DatapThetaValCenter[@]}]"
	    #echo "Beam Energy per Run: [${DataEbeamValCenter[@]}]"
	fi

	# Checks that array isn't empty
	if [[ ${#dummy_center_tmp[@]} -ne 0 ]]; then
	    echo
	    echo "Calculating dummy total effective charge center..."
	    PYCENTERSTRING=$(python3 findEffectiveCharge.py ${EffData} "${dummy_center_tmp[*]}" ${ParticleType})
	    arr1=()
	    arr2=()
	    arr3=()
	    arr4=()
	    arr5=()
	    arr6=()
	    arr7=()	    
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYCENTERSTRING"
	    DummyEffChargeValCenter=("${arr1[@]}")
	    DummyEffChargeErrCenter=("${arr2[@]}")
	    DummyEffValCenter=("${arr3[@]}")
	    DummyEffErrCenter=("${arr4[@]}")
	    DummypThetaValCenter=("${arr5[@]}")
	    DummyEbeamValCenter=("${arr6[@]}")
	    dummy_center=("${arr7[@]}")
	    if [ "${#dummy_center_tmp[@]}" -ne "${#dummy_center[@]}" ]; then
		echo "Removing bad center dummy runs..."
		for run in "${dummy_center_tmp[@]}"; do
		    if [[ ! " ${dummy_center[@]} " =~ " $run " ]]; then
			echo "        Removing run: $run"
		    fi
		done
	    fi
	    #echo ${DummyEffChargeVal[*]}
	    # Sums the array to get the total effective charge and effective charge error
	    PYCENTERTOTEFFCHARGE=$(python3 calcTotalEffectiveCharge.py "${DummyEffChargeValCenter[*]}" "${DummyEffChargeErrCenter[*]}")
	    arr1=()
	    arr2=()
	    itt=0
	    while read line; do
		itt=$((itt+1))
		# split the line into an array based on space
		IFS=' ' read -ra line_array <<< "$line"
		# store the elements in the corresponding array
		eval "arr$itt=(\"\${line_array[@]}\")"
	    done <<< "$PYCENTERTOTEFFCHARGE"
	    TotDummyEffChargeValCenter=("${arr1[@]}")
	    TotDummyEffChargeErrCenter=("${arr2[@]}")
	    echo "Total Effective Charge Center: ${TotDummyEffChargeValCenter} mC"
	    echo "Total Effective Charge Center Error: ${TotDummyEffChargeErrCenter}"
	    echo "Run numbers: [${dummy_center[@]}]"
	    #echo "Effective Charge per Run: [${DummyEffChargeValCenter[@]}]"
	    #echo "Effective Charge Error per Run: [${DummyEffChargeErrCenter[@]}]"
	    #echo "Efficiency per Run: [${DummyEffValCenter[@]}]"
	    #echo "Efficiency Error per Run: [${DummyEffErrCenter[@]}]"
	    #echo "Theta per Run: [${DummypThetaValCenter[@]}]"
	    #echo "Beam Energy per Run: [${DummyEbeamValCenter[@]}]"
	fi
	
	cd "${LTANAPATH}/src"

	if [ $j = "low" ]; then
	    echo
	    echo "Finding t/phi bins for low epsilon..."
	else
	    echo
	    echo "Using low epsilon t/phi bins for high epsilon..."
	fi
	
	if [ ${#data_right[@]} -eq 0 ]; then
	    python3 main.py ${KIN} ${W} ${Q2} ${LOEPS} ${HIEPS} ${OutDATAFilename} ${OutDUMMYFilename} ${OutFullAnalysisFilename} ${TMIN} ${TMAX} ${NumtBins} ${NumPhiBins} "0" "${data_left[*]}" "${data_center[*]}" "0" ${TotDataEffChargeValLeft} ${TotDataEffChargeValCenter} "0" ${TotDummyEffChargeValLeft} ${TotDummyEffChargeValCenter} "0" ${TotDataEffChargeErrLeft} ${TotDataEffChargeErrCenter} "0" ${TotDummyEffChargeErrLeft} ${TotDummyEffChargeErrCenter} "0" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" "0" "${DataEffErrLeft[*]}" "${DataEffErrCenter[*]}" ${EffData} ${ParticleType} $j "0" "${DatapThetaValLeft[*]}" "${DatapThetaValCenter[*]}" "0" "${DataEbeamValLeft[*]}" "${DataEbeamValCenter[*]}" ${POL} ${formatted_date} ${DEBUG}
	    # Check the exit status of the Python script
	    if [ $? -ne 0 ]; then
		echo
		echo
		echo "1 ERROR: Python script failed!"
		echo "       See error above..."
		exit 1
	    fi
	else
	    python3 main.py ${KIN} ${W} ${Q2} ${LOEPS} ${HIEPS} ${OutDATAFilename} ${OutDUMMYFilename} ${OutFullAnalysisFilename} ${TMIN} ${TMAX} ${NumtBins} ${NumPhiBins} "${data_right[*]}" "${data_left[*]}" "${data_center[*]}" ${TotDataEffChargeValRight} ${TotDataEffChargeValLeft} ${TotDataEffChargeValCenter} ${TotDummyEffChargeValRight} ${TotDummyEffChargeValLeft} ${TotDummyEffChargeValCenter} ${TotDataEffChargeErrRight} ${TotDataEffChargeErrLeft} ${TotDataEffChargeErrCenter} ${TotDummyEffChargeErrRight} ${TotDummyEffChargeErrLeft} ${TotDummyEffChargeErrCenter} "${DataEffValRight[*]}" "${DataEffValLeft[*]}" "${DataEffValCenter[*]}" "${DataEffErrRight[*]}" "${DataEffErrLeft[*]}" "${DataEffErrCenter[*]}" ${EffData} ${ParticleType} $j "${DatapThetaValRight[*]}" "${DatapThetaValLeft[*]}" "${DatapThetaValCenter[*]}" "${DataEbeamValRight[*]}" "${DataEbeamValLeft[*]}" "${DataEbeamValCenter[*]}" ${POL} ${formatted_date} ${DEBUG}
	    # Check the exit status of the Python script
	    if [ $? -ne 0 ]; then
		echo
		echo
		echo "1 ERROR: Python script failed!"
		echo "       See error above..."
		exit 1
	    fi
	fi	

	if [ $j = "low" ]; then
	    echo
	    echo
	    echo
	    echo "Low Epsilon Completed!"
	else
	    echo
	    echo
	    echo
	    echo "High Epsilon Completed!"	
	fi
    done
else
    # Need to rerun loop separately so that the combined files for high and low epsilon exists for diamond cut script
    declare -a EPS=("low" "high")
    for j in "${EPS[@]}"
    do

	# Redefine epsilon based on loop
	EPSILON=$j

	# Function that calls python script to grab run numbers
	grab_runs () {
	    RunList=$1
	    # Location of list of run lists
	    INPDIR="${UTILPATH}/run_list/${ANATYPE}LT/${RunList}"
	    if [[ -e $INPDIR ]]; then
		cd "${LTANAPATH}/src/setup"
		RunNumArr=$(python3 getRunNumbers.py $INPDIR)
		echo $RunNumArr
	    else
		exit
	    fi
	}

	echo
	echo "---------------------------------------------------------"
	echo
	echo "${formatted_date}"
	echo
	echo "Running weight iteration analysis for Q2=${Q2}, W=${W}, ${EPSILON} setting..."
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
	    ##############
	    # HARD CODED #
	    ##############
	    if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    # Define run list based off kinematics selected
		    file_right_dummy="Q5p5W3p02right_${EPSILON}e_dummy"
		    file_right="Q5p5W3p02right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    # Converts python output to bash array
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=5p5, W=3p02
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=5p5, W=3p02
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.1838
		HIEPS=0.5291
		KIN="Q5p5W3p02_${EPSILON}e"
	    fi

	    if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q4p4W2p74right_${EPSILON}e_dummy"
		    file_right="Q4p4W2p74right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=4p4, W=2p74
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=4p4, W=2p74
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.4805
		HIEPS=0.7148
		KIN="Q4p4W2p74_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W3p14right_${EPSILON}e_dummy"
		    file_right="Q3p0W3p14right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=3p14
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=3p14
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.3935
		HIEPS=0.6668
		KIN="Q3p0W3p14_${EPSILON}e"
	    fi

	    if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q3p0W2p32right_${EPSILON}e_dummy"
		    file_right="Q3p0W2p32right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=3p0, W=2p32
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=3p0, W=2p32
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.5736
		HIEPS=0.8791
		KIN="Q3p0W2p32_${EPSILON}e"
	    fi

	    if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q2p1W2p95right_${EPSILON}e_dummy"
		    file_right="Q2p1W2p95right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=2p1, W=2p95
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=2p1, W=2p95
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.2477
		HIEPS=0.7864
		KIN="Q2p1W2p95_${EPSILON}e"
	    fi

	    if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
		if [[ $i = "RIGHT" ]]; then
		    file_right_dummy="Q0p5W2p40right_${EPSILON}e_dummy"
		    file_right="Q0p5W2p40right_${EPSILON}e"
		    echo "Reading in run numbers for right file ${file_right_dummy}..."
		    IFS=', ' read -r -a dummy_right <<< "$( grab_runs ${file_right_dummy} )"             # RIGHT, Q2=0p5, W=2p40
		    echo "Dummy Run Numbers: [${dummy_right[@]}]"
		    echo
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
		    echo
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
		    echo
		    echo "Reading in run numbers for center file ${file_center}..."
		    IFS=', ' read -r -a data_center <<< "$( grab_runs ${file_center} )"		 # CENTER, Q2=0p5, W=2p40
		    echo "Data Run Numbers: [${data_center[@]}]"
		    echo
		fi
		LOEPS=0.4515
		HIEPS=0.6979
		KIN="Q0p5W2p40_${EPSILON}e"
	    fi
	    ##############
	    ##############
	    ##############
	done

	# Define input and output file names
	OutDATAFilename="Analysed_Data_${KIN}"
	OutDUMMYFilename="Analysed_Dummy_${KIN}"
	OutFullAnalysisFilename="FullAnalysis_${KIN}"
	
	cd "${LTANAPATH}/src"

	if [ $j = "low" ]; then
	    echo
	    echo "Finding new simc weight for low epsilon..."
	else
	    echo
	    echo "Finding new simc weight for for high epsilon..."
	fi

	python3 main_iter.py ${KIN} ${W} ${Q2} ${LOEPS} ${HIEPS} ${ParticleType} $j ${POL} ${OutFullAnalysisFilename} ${formatted_date} ${NumtBins} ${NumPhiBins} ${DEBUG}

	# Check the exit status of the Python script
	if [ $? -ne 0 ]; then
	    echo
	    echo
	    echo "1 ERROR: Python script failed!"
	    echo "       See error above..."
	    exit 1
	fi	
	
	if [ $j = "low" ]; then
	    echo
	    echo
	    echo
	    echo "Low Epsilon Completed!"
	else
	    echo
	    echo
	    echo
	    echo "High Epsilon Completed!"	
	fi	    
    done
fi
    
echo
echo
echo
echo "Script Complete!"
