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
while getopts 'haos' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_HeeP_Analysis.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Plots data vs simc"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
        echo "    -h, help"
        echo "    -a, analyze"
	echo "        coin -> KIN=arg1"
	echo "        sing -> SPEC=arg1 KIN=arg2 (requires -s flag)"
	echo "    -s, single arm"
	echo "    -o, offset to replay applied "
	echo "        (this is a label, offsets must be applied explicitly in replay before this step)"
        exit 0
        ;;
        a) a_flag='true' ;;
        o) o_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# When any flag is used then the user input changes argument order
if [[ $a_flag = "true" || $o_flag = "true" || $s_flag = "true" ]]; then
    if [[ $s_flag = "true" ]]; then
	spec=$(echo "$2" | tr '[:upper:]' '[:lower:]')
	SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
	KIN=$3
	ROOTPREFIX=replay_${spec}_heep
    else
	KIN=$2
	ROOTPREFIX=replay_coin_heep
    fi
else
    KIN=$1
    ROOTPREFIX=replay_coin_heep
fi

# When just single flag used...
if [[ $s_flag = "true" ]]; then
    ##############
    # HARD CODED #
    ##############
    # Defines efficiency table to use
    if [[ $SPEC = "HMS" ]]; then
	EffData="hms_heep_HeePSing_efficiency_data_2022_07_28.csv"
    else
	EffData="shms_heep_HeePSing_efficiency_data_2022_07_28.csv"
    fi
    InDATAFilename="Raw_Data_${SPEC}_${KIN}.root"
    InDUMMYFilename="Raw_DummyData_${SPEC}_${KIN}.root"
    if [[ $o_flag = "true" ]]; then
	InSIMCFilename="Heep_${SPEC}_${KIN}_Offset.root"
    else
	InSIMCFilename="Heep_${SPEC}_${KIN}.root"
    fi
    OutDATAFilename="Analysed_Data_${SPEC}_${KIN}"
    OutDUMMYFilename="Analysed_DummyData_${SPEC}_${KIN}"
    if [[ $o_flag = "true" ]]; then
	OutFullAnalysisFilename="FullAnalysis_${SPEC}_${KIN}_Offset"
    else
	OutFullAnalysisFilename="FullAnalysis_${SPEC}_${KIN}"
    fi
else
    ##############
    # HARD CODED #
    ##############
    # Defines efficiency table to use
    #EffData="coin_heep_HeePCoin_efficiency_data_2022_09_09.csv"
    EffData="coin_heep_HeePCoin_efficiency_data_2022_12_02.csv"
    InDATAFilename="Raw_Data_${KIN}.root"
    InDUMMYFilename="Raw_DummyData_${KIN}.root"
    if [[ $o_flag = "true" ]]; then
	InSIMCFilename="Heep_Coin_${KIN}_Offset.root"
    else
	InSIMCFilename="Heep_Coin_${KIN}.root"
    fi    
    OutDATAFilename="Analysed_Data_${KIN}"
    OutDUMMYFilename="Analysed_DummyData_${KIN}"
    if [[ $o_flag = "true" ]]; then
	OutFullAnalysisFilename="FullAnalysis_${KIN}_Offset"
    else
	OutFullAnalysisFilename="FullAnalysis_${KIN}"
    fi
fi

# Function that calls python script to grab run numbers
grab_runs () {
    RunList=$1
    INPDIR="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/${RunList}"
    if [[ -e $INPDIR ]]; then
	cd "${LTANAPATH}/src/setup"
	RunNumArr=$(python3 getRunNumbers.py $INPDIR)
	echo $RunNumArr
    else
	exit
    fi
}

##############
# HARD CODED #
##############
# Define heep run numbers for a particular setting
if [[ $KIN = "10p6" && $s_flag != "true" ]]; then
    # Define run list based off kinematics selected
    file_dummy="HeePCoin_10p6_Autumn18_dummy"
    file="HeePCoin_10p6_Autumn18"
    echo "Reading in run numbers for file ${file_dummy}..."
    # Converts python output to bash array
    IFS=', ' read -r -a dummydata <<< "$( grab_runs ${file_dummy} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Dummy Run Numbers: [${dummydata[@]}]"
    echo
    echo "Reading in run numbers for file ${file}..."
    IFS=', ' read -r -a data <<< "$( grab_runs ${file} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Data Run Numbers: [${data[@]}]"
    echo
elif [[ $KIN = "8p2" && $s_flag != "true" ]]; then
    # Define run list based off kinematics selected
    file_dummy="HeePCoin_8p2_Spring19_dummy"
    file="HeePCoin_8p2_Spring19"
    echo "Reading in run numbers for file ${file_dummy}..."
    # Converts python output to bash array
    IFS=', ' read -r -a dummydata <<< "$( grab_runs ${file_dummy} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Dummy Run Numbers: [${dummydata[@]}]"
    echo
    echo "Reading in run numbers for file ${file}..."
    IFS=', ' read -r -a data <<< "$( grab_runs ${file} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Data Run Numbers: [${data[@]}]"
    echo
elif [[ $KIN = "6p2" && $s_flag != "true" ]]; then
    # Define run list based off kinematics selected
    file_dummy="HeePCoin_6p2_Spring19_dummy"
    file="HeePCoin_6p2_Spring19"
    echo "Reading in run numbers for file ${file_dummy}..."
    # Converts python output to bash array
    IFS=', ' read -r -a dummydata <<< "$( grab_runs ${file_dummy} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Dummy Run Numbers: [${dummydata[@]}]"
    echo
    echo "Reading in run numbers for file ${file}..."
    IFS=', ' read -r -a data <<< "$( grab_runs ${file} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Data Run Numbers: [${data[@]}]"
    echo
elif [[ $KIN = "4p9" && $s_flag != "true" ]]; then
    # Define run list based off kinematics selected
    file_dummy="HeePCoin_4p9_Autumn18_dummy"
    file="HeePCoin_4p9_Autumn18"
    echo "Reading in run numbers for file ${file_dummy}..."
    # Converts python output to bash array
    IFS=', ' read -r -a dummydata <<< "$( grab_runs ${file_dummy} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Dummy Run Numbers: [${dummydata[@]}]"
    echo
    echo "Reading in run numbers for file ${file}..."
    IFS=', ' read -r -a data <<< "$( grab_runs ${file} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Data Run Numbers: [${data[@]}]"
    echo
elif [[ $KIN = "3p8" && $s_flag != "true" ]]; then
    # Define run list based off kinematics selected
    file_dummy="HeePCoin_3p8_Autumn18_dummy"
    file="HeePCoin_3p8_Autumn18"
    echo "Reading in run numbers for file ${file_dummy}..."
    # Converts python output to bash array
    IFS=', ' read -r -a dummydata <<< "$( grab_runs ${file_dummy} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Dummy Run Numbers: [${dummydata[@]}]"
    echo
    echo "Reading in run numbers for file ${file}..."
    IFS=', ' read -r -a data <<< "$( grab_runs ${file} )"             # RIGHT, Q2=5p5, W=3p02
    echo "Data Run Numbers: [${data[@]}]"
    echo
elif [[ $KIN = "10p6" && $s_flag = "true" ]]; then
    declare -a data=(7974 7975 7976)
    #    declare -a data=(7974) # Just one test run
    declare -a dummydata=(7977)    
elif [[ $KIN = "8p2" && $s_flag != "true" ]]; then
    declare -a data=(111)
    declare -a dummydata=(111)
elif [[ $KIN = "4p9" && $s_flag = "true" ]]; then
    declare -a data=(111)
    declare -a dummydata=(111)
elif [[ $KIN = "3p8" && $s_flag = "true" ]]; then
    declare -a data=(111)
    declare -a dummydata=(111)
else
    echo "Invalid kinematic setting, ${KIN}"
    exit 128
fi

# When analysis flag is used then the analysis script (Analysed_{SPEC}.py)
# will create a new root file per run number which are combined using hadd
if [[ $a_flag = "true" ]]; then
    if [[ $s_flag = "true" ]]; then
	cd "${LTANAPATH}/src/HeeP/SING"
	echo
	echo "Analysing ${SPEC} data..."
	echo

	for i in "${data[@]}"
	do
	    echo
	    echo "-----------------------------"
	    echo "Analysing data run $i..."
	    echo "-----------------------------"
	    echo
	    python3 Analysed_SING.py "$i" ${SPEC} | tee ../../../log/Analysed_HeeP_SING_${SPEC}_$i.log
	    #root -l <<EOF 
	    #.x $LTANAPATH/Analysed_SING.C("$InDATAFilename","$OutDATAFilename")
	    #EOF
	done
	cd "${LTANAPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."  
	hadd -f ${OutDATAFilename}.root *_-1_${SPEC}_Raw_Data.root
	#rm -f *_-1_${SPEC}_Raw_Data.root

	cd "${LTANAPATH}/src/HeeP/SING"    
	echo
	echo "Analysing ${SPEC} dummy data..."
	echo

	for i in "${dummydata[@]}"
	do
	    echo
	    echo "-----------------------------------"
	    echo "Analysing dummy data run $i..."
	    echo "-----------------------------------"
	    echo
	    python3 Analysed_SING.py "$i" ${SPEC} | tee ../../../log/Analysed_HeeP_SING_${SPEC}_$i.log
	    #root -l <<EOF 
	    #.x $LTANAPATH/Analysed_SING.C("$InDUMMYFilename","$OutDUMMYFilename")
	    #EOF
	done
	cd "${LTANAPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."
	hadd -f ${OutDUMMYFilename}.root *_-1_${SPEC}_Raw_Data.root
	#rm -f *_-1_${SPEC}_Raw_Data.root	
    else
	cd "${LTANAPATH}/src/HeeP/COIN"
	echo
	echo "Analysing data..."
	echo

	for i in "${data[@]}"
	do
	    echo
	    echo "-----------------------------"
	    echo "Analysing data run $i..."
	    echo "-----------------------------"
	    echo
	    python3 Analysed_COIN.py "$i" | tee ../../../log/Analysed_HeeP_COIN_$i.log
	    #root -l <<EOF 
	    #.x $LTANAPATH/Analysed_COIN.C("$InDATAFilename","$OutDATAFilename")
	    #EOF
	done
	cd "${LTANAPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."  
	hadd -f ${OutDATAFilename}.root *_-1_Raw_Data.root
	for i in *_-1_Raw_Data.root; do mv -- "$i" "${i%_-1_Raw_Data.root}_-1_Raw_Target.root"; done

	cd "${LTANAPATH}/src/HeeP/COIN"    
	echo
	echo "Analysing dummy data..."
	echo

	for i in "${dummydata[@]}"
	do
	    echo
	    echo "-----------------------------------"
	    echo "Analysing dummy data run $i..."
	    echo "-----------------------------------"
	    echo
	    python3 Analysed_COIN.py "$i" | tee ../../../log/Analysed_HeeP_COIN_$i.log
	    #root -l <<EOF 
	    #.x $LTANAPATH/Analysed_COIN.C("$InDUMMYFilename","$OutDUMMYFilename")
	    #EOF
	done
	cd "${LTANAPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."
	hadd -f ${OutDUMMYFilename}.root *_-1_Raw_Data.root
	for i in *_-1_Raw_Data.root; do mv -- "$i" "${i%_-1_Raw_Data.root}_-1_Raw_Dummy.root"; done
    fi
fi

cd "${LTANAPATH}/src/setup"

# Checks that array isn't empty
if [[ ${#data[@]} -ne 0 ]]; then
    echo
    echo "Calculating data total effective charge ..."
    PYSTRING=$(python3 findEffectiveCharge.py ${EffData} "${data[*]}")
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
    done <<< "$PYSTRING"
    DataChargeVal=("${arr1[@]}")
    DataChargeErr=("${arr2[@]}")
    DataEffVal=("${arr3[@]}")
    DataEffErr=("${arr4[@]}")
    DatapThetaVal=("${arr5[@]}")
    DataEbeamVal=("${arr6[@]}")
    #echo ${DataChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DataChargeSum=$(IFS=+; echo "$((${DataChargeVal[*]}))") # Only works for integers
    echo "Total Charge : ${DataChargeSum} uC"
fi

# Checks that array isn't empty
if [[ ${#dummydata[@]} -ne 0 ]]; then
    echo
    echo "Calculating dummy data total effective charge ..."
    PYSTRING=$(python3 findEffectiveCharge.py ${EffData} "${dummydata[*]}")
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
    done <<< "$PYSTRING"
    DummyChargeVal=("${arr1[@]}")
    DummyChargeErr=("${arr2[@]}")
    DummyEffVal=("${arr3[@]}")
    DummyEffErr=("${arr4[@]}")
    DummypThetaVal=("${arr5[@]}")
    DummyEbeamVal=("${arr6[@]}")
    #echo ${DummyChargeVal[*]}
    # Sums the array to get the total effective charge
    # Note: this must be done as an array! This is why uC is used at this step
    #       and later converted to C
    DummyChargeSum=$(IFS=+; echo "$((${DummyChargeVal[*]}))") # Only works for integers
    echo "Total Charge : ${DummyChargeSum} uC"
fi

# Finally, run the plotting script
if [[ $s_flag = "true" ]]; then
    cd "${LTANAPATH}/src/HeeP/SING"
    python3 HeepSing.py ${KIN} "${OutDATAFilename}.root" $DataChargeSum "${DataEffVal[*]}" "${OutDUMMYFilename}.root" $DummyChargeSum "${DummyEffVal[*]}" ${InSIMCFilename} ${OutFullAnalysisFilename} ${EffData} ${SPEC}
else
    cd "${LTANAPATH}/src/HeeP/COIN"
    python3 HeepCoin.py ${KIN} "${OutDATAFilename}.root" $DataChargeSum "${DataEffVal[*]}" "${data[*]}" "${OutDUMMYFilename}.root" $DummyChargeSum "${DummyEffVal[*]}" "${dummydata[*]}" ${InSIMCFilename} ${OutFullAnalysisFilename} ${EffData}
fi

cd "${LTANAPATH}"
evince "OUTPUT/Analysis/HeeP/${OutFullAnalysisFilename}.pdf"

echo
echo
echo
echo "Script Complete!"
