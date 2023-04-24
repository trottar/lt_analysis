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

# Flag definitions (flags: h, p)
while getopts 'hp' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./applyCuts_Prod.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Creates root files with PID and CT cuts applied"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
        echo "    -h, help"
	echo "     EPSILON=arg1, PHIVAL=arg2, Q2=arg3, W=arg4, target=arg5, RUNNUM=arg6"
	echo "    -p, specify particle type (kaon, pion, or proton). Otherwise runs for all."
	echo "        EPSILON=arg1, PHIVAL=arg2, Q2=arg3, W=arg4, target=arg5, RUNNUM=arg6, ParticleType=arg7"
	echo
	echo " Avaliable Kinematics..."
	echo "                      EPSILON={high,low}"
	echo "                      PHIVAL={right,left,center}"
	echo "                      Q2=5p5, W=3p02"
	echo "                      Q2=4p4, W=2p74"
	echo "                      Q2=3p0, W=3p14"
	echo "                      Q2=3p0, W=2p32"
	echo "                      Q2=2p1, W=2p95"
	echo "                      Q2=0p5, W=2p40"
        exit 0
        ;;
	p) p_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ $p_flag = "true" ]]; then

    EPSILON=$(echo "$2" | tr '[:upper:]' '[:lower:]')
    PHIVAL=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    Q2=$4
    W=$5
    TargetType=$(echo "$6" | tr '[:upper:]' '[:lower:]')
    RUNNUM=$7
    ParticleType=$8
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
    if [[ -z "$3" || ! "$PHIVAL" =~ right|left|center ]]; then # Check the 1st argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid phi value..."
	while true; do
	    echo ""
	    read -p "Phi value must be - right - left - center - Case Sensitive! - or press ctrl-c to exit : " PHIVAL
	    case $PHIVAL in
		'');; # If blank, prompt again
		'right'|'left'|'center') break;; # If a valid option, break the loop and continue
	    esac
	done
    fi
    if [[ -z "$4" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
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
    if [[ -z "$5" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
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
    if [[ -z "$6" || ! "$TargetType" =~ lh2|dummy ]]; then # Check the 3rd argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid target type..."
	while true; do
	    echo ""
	    read -p "Target type must be one of - [lh2 - dummy] - or press ctrl-c to exit : " TargetType
	    case $TargetType in
		'');; # If blank, prompt again
		'lh2'|'dummy') break;; # If a valid option, break the loop and continue
	    esac
	done
    fi
    if [[ -z "$8" || ! "$ParticleType" =~ kaon|pion|proton ]]; then # Check the 3rd argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid target type..."
	while true; do
	    echo ""
	    read -p "Particle type must be one of - [kaon - pion - proton] - or press ctrl-c to exit : " ParticleType
	    case $ParticleType in
		'');; # If blank, prompt again
		'kaon'|'pion'|'proton') break;; # If a valid option, break the loop and continue
	    esac
	done
    fi
    
    echo
    echo "---------------------------------------------------------"
    echo
    echo "Creating analysis root file for Q2=${Q2}, W=${W}, EPSILON=${EPSILON}, PHIVAL=${PHIVAL} setting..."
    echo
    echo "Only running for ${ParticleType}"
    echo
    echo "---------------------------------------------------------"
    echo    
    
else
    
    EPSILON=$(echo "$1" | tr '[:upper:]' '[:lower:]')
    PHIVAL=$(echo "$2" | tr '[:upper:]' '[:lower:]')
    Q2=$3
    W=$4
    TargetType=$(echo "$5" | tr '[:upper:]' '[:lower:]')
    RUNNUM=$6
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
    if [[ -z "$2" || ! "$PHIVAL" =~ right|left|center ]]; then # Check the 1st argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid phi value..."
	while true; do
	    echo ""
	    read -p "Phi value must be - right - left - center - Case Sensitive! - or press ctrl-c to exit : " PHIVAL
	    case $PHIVAL in
		'');; # If blank, prompt again
		'right'|'left'|'center') break;; # If a valid option, break the loop and continue
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
    if [[ -z "$5" || ! "$TargetType" =~ lh2|dummy ]]; then # Check the 3rd argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid target type..."
	while true; do
	    echo ""
	    read -p "Target type must be one of - [lh2 - dummy] - or press ctrl-c to exit : " TargetType
	    case $TargetType in
		'');; # If blank, prompt again
		'lh2'|'dummy') break;; # If a valid option, break the loop and continue
	    esac
	done
    fi

    echo
    echo "---------------------------------------------------------"
    echo
    echo "Creating analysis root file for Q2=${Q2}, W=${W}, EPSILON=${EPSILON}, PHIVAL=${PHIVAL} setting..."
    echo 
    echo "Running for kaon, pion, proton"    
    echo
    echo "---------------------------------------------------------"
    echo
fi

if [[ $Q2 = "5p5" && $W = "3p02" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    # Define run list based off kinematics selected
	    file_right="Q5p5W3p02right_${EPSILON}e_dummy"
	else
	    file_right="Q5p5W3p02right_${EPSILON}e"		   
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q5p5W3p02left_${EPSILON}e_dummy"
	else
	    file_left="Q5p5W3p02left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q5p5W3p02center_${EPSILON}e_dummy"
	else
	    file_center="Q5p5W3p02center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.1838
    else
	EPSVAL=0.5291
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q5p5W3p02_${EPSILON}e_dummy"
    else
	KIN="Q5p5W3p02_${EPSILON}e"
    fi
fi
if [[ $Q2 = "4p4" && $W = "2p74" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_right="Q4p4W2p74right_${EPSILON}e_dummy"
	else
	    file_right="Q4p4W2p74right_${EPSILON}e"
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q4p4W2p74left_${EPSILON}e_dummy"
	else
	    file_left="Q4p4W2p74left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q4p4W2p74center_${EPSILON}e_dummy"
	else
	    file_center="Q4p4W2p74center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo	
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.4805
    else
	EPSVAL=0.7148
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q4p4W2p74_${EPSILON}e_dummy"
    else
	KIN="Q4p4W2p74_${EPSILON}e"
    fi
fi
if [[ $Q2 = "3p0" && $W = "3p14" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_right="Q3p0W3p14right_${EPSILON}e_dummy"
	else
	    file_right="Q3p0W3p14right_${EPSILON}e"
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q3p0W3p14left_${EPSILON}e_dummy"
	else
	    file_left="Q3p0W3p14left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo	
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q3p0W3p14center_${EPSILON}e_dummy"
	else
	    file_center="Q3p0W3p14center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.3935
    else
	EPSVAL=0.6668
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q3p0W3p14_${EPSILON}e_dummy"
    else
	KIN="Q3p0W3p14_${EPSILON}e"
    fi
fi
if [[ $Q2 = "3p0" && $W = "2p32" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_right="Q3p0W2p32right_${EPSILON}e_dummy"
	else
	    file_right="Q3p0W2p32right_${EPSILON}e"
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q3p0W2p32left_${EPSILON}e_dummy"
	else
	    file_left="Q3p0W2p32left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q3p0W2p32center_${EPSILON}e_dummy"
	else
	    file_center="Q3p0W2p32center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.5736
    else
	EPSVAL=0.8791
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q3p0W2p32_${EPSILON}e_dummy"
    else
	KIN="Q3p0W2p32_${EPSILON}e"
    fi
fi
if [[ $Q2 = "2p1" && $W = "2p95" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_right="Q2p1W2p95right_${EPSILON}e_dummy"
	else
	    file_right="Q2p1W2p95right_${EPSILON}e"
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q2p1W2p95left_${EPSILON}e_dummy"
	else
	    file_left="Q2p1W2p95left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q2p1W2p95center_${EPSILON}e_dummy"
	else
	    file_center="Q2p1W2p95center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.2477
    else
	EPSVAL=0.7864
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q2p1W2p95_${EPSILON}e_dummy"
    else
	KIN="Q2p1W2p95_${EPSILON}e"
    fi
fi        
if [[ $Q2 = "0p5" && $W = "2p40" ]]; then
    if [[ $PHIVAL = "right" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_right="Q0p5W2p40right_${EPSILON}e_dummy"
	else
	    file_right="Q0p5W2p40right_${EPSILON}e"
	fi
	echo "Reading in run numbers for right file ${file_right}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "left" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_left="Q0p5W2p40left_${EPSILON}e_dummy"
	else
	    file_left="Q0p5W2p40left_${EPSILON}e"
	fi
	echo "Reading in run numbers for left file ${file_left}..."
	echo "Run Number: ${RUNNUM}"
	echo
    elif [[ $PHIVAL = "center" ]]; then
	if [[ $TargetType = "dummy" ]]; then
	    file_center="Q0p5W2p40center_${EPSILON}e_dummy"
	else
	    file_center="Q0p5W2p40center_${EPSILON}e"
	fi
	echo "Reading in run numbers for center file ${file_center}..."
	echo "Run Number: ${RUNNUM}"
	echo
    fi
    if [[ ${EPSILON} == "low" ]]; then
	EPSVAL=0.4515
    else
	EPSVAL=0.6979
    fi
    if [[ $TargetType = "dummy" ]]; then
	KIN="Q0p5W2p40_${EPSILON}e_dummy"
    else
	KIN="Q0p5W2p40_${EPSILON}e"
    fi
fi    


if [[ $p_flag = "true" ]]; then

    # Define input and output file names
    InDATAFilename="Proc_Data_${ParticleType}_${KIN}.root"
    OutDATAFilename="Analysed_Data_${ParticleType}_${KIN}"
    OutFullAnalysisFilename="FullAnalysis_${ParticleType}_${KIN}"

    # The analysis script (Analysed_Prod.py) will create a new root file per run number
    if [ ${PHIVAL} = "right" ]; then
	echo
	echo "Analysing right data..."
	echo
	echo
	echo "--------------------------------"
	echo "Analysing right data run ${RUNNUM}..."
	echo "--------------------------------"
	echo
	cd "${LTANAPATH}/scripts/Prod"
	python3 Analysed_Prod.py "${RUNNUM}" "${ParticleType}" | tee ../../log/Right_Analysed_Prod_${ParticleType}_${RUNNUM}.log
	echo
    fi

    # Checks that array isn't empty
    if [ ${PHIVAL} = "left" ]; then
	echo
	echo "Analysing left data..."
	echo
	echo
	echo "--------------------------------"
	echo "Analysing left data run ${RUNNUM}..."
	echo "--------------------------------"
	echo
	cd "${LTANAPATH}/scripts/Prod"
	python3 Analysed_Prod.py "${RUNNUM}" "${ParticleType}" | tee ../../log/Left_Analysed_Prod_${ParticleType}_${RUNNUM}.log
    fi

    # Checks that array isn't empty
    if [ ${PHIVAL} = "center" ]; then
	echo
	echo "Analysing center data..."
	echo
	echo
	echo "--------------------------------"
	echo "Analysing center data run ${RUNNUM}..."
	echo "--------------------------------"
	echo
	cd "${LTANAPATH}/scripts/Prod"
	python3 Analysed_Prod.py "${RUNNUM}" "${ParticleType}" | tee ../../log/Center_Analysed_Prod_${ParticleType}_${RUNNUM}.log
    fi

else
    declare -a ParticleTypes=("kaon" "pion" "proton")
    for i in "${ParticleTypes[@]}"
    do

	# Define input and output file names
	InDATAFilename="Proc_Data_${i}_${KIN}.root"
	OutDATAFilename="Analysed_Data_${i}_${KIN}"
	OutFullAnalysisFilename="FullAnalysis_${i}_${KIN}"

	# The analysis script (Analysed_Prod.py) will create a new root file per run number
	if [ ${PHIVAL} = "right" ]; then
	    echo
	    echo "Analysing right data..."
	    echo
	    echo
	    echo "--------------------------------"
	    echo "Analysing right data run ${RUNNUM}..."
	    echo "--------------------------------"
	    echo
	    cd "${LTANAPATH}/scripts/Prod"
	    python3 Analysed_Prod.py "${RUNNUM}" "${i}" | tee ../../log/Right_Analysed_Prod_${i}_${RUNNUM}.log
	    echo
	fi

	# Checks that array isn't empty
	if [ ${PHIVAL} = "left" ]; then
	    echo
	    echo "Analysing left data..."
	    echo
	    echo
	    echo "--------------------------------"
	    echo "Analysing left data run ${RUNNUM}..."
	    echo "--------------------------------"
	    echo
	    cd "${LTANAPATH}/scripts/Prod"
	    python3 Analysed_Prod.py "${RUNNUM}" "${i}" | tee ../../log/Left_Analysed_Prod_${i}_${RUNNUM}.log
	fi

	# Checks that array isn't empty
	if [ ${PHIVAL} = "center" ]; then
	    echo
	    echo "Analysing center data..."
	    echo
	    echo
	    echo "--------------------------------"
	    echo "Analysing center data run ${RUNNUM}..."
	    echo "--------------------------------"
	    echo
	    cd "${LTANAPATH}/scripts/Prod"
	    python3 Analysed_Prod.py "${RUNNUM}" "${i}" | tee ../../log/Center_Analysed_Prod_${i}_${RUNNUM}.log
	    
	fi
    done
fi
	
echo
echo
echo
echo "Script Complete!"
