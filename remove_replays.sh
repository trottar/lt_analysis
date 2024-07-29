#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-07-28 20:26:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
#!/bin/tcsh


# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ ${HOSTNAME} = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /group/c-kaonlt/USERS/${USER}/hallc_replay_lt/UTIL_KAONLT/bin/python/ltsep replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
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

file_name="Prod_Test"
numbers_to_match=()
IFS=', ' read -r -a numbers_to_match <<< "$( grab_runs ${file_name} )"
echo "Run Numbers: [${numbers_to_match[@]}]"

replay_root_path="${ROOTPATH}/${ANATYPE}LT"

while true; do
    # Prompt for confirmation before proceeding
    read -p "Are you sure you want to remove files with specified numbers in the filename? (yes/no): " answer

    case "$answer" in
        [Yy]* ) 
            # Loop through each number in the list
            for number in "${numbers_to_match[@]}"
            do
                # Remove files with the specified number in the filename
                rm "$replay_root_path"/*"$number"*
                echo "Removed files from $replay_root_path with number $number"
            done
            break ;;
        [Nn]* ) 
            echo "Operation aborted."
            exit ;;
        * ) 
            echo "Please answer yes or no." ;;
    esac
done
