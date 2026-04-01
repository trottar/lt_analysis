#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2025-06-21 15:31:51 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#


# Runs a repo-local ltsep wrapper so batch jobs do not depend on upstream
# getPathDict.py calling os.getlogin() on worker nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
if [[ ! -f "${LTSEP_FIELDS_SCRIPT}" ]]; then
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
fi
PATHFILE_INFO="$(python3 "${LTSEP_FIELDS_SCRIPT}" "$PWD")"
path_rc=$?
if [[ "${path_rc}" -ne 0 || -z "${PATHFILE_INFO}" ]]; then
    echo "ERROR: failed to resolve ltsep paths via ${LTSEP_FIELDS_SCRIPT}" >&2
    exit 1
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
SKIMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f16`
SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f17`
LTANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f18`

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

# Pion only
#Q2="0p4" #0p38
#W="2p20"

# TODO
#Q2="0p5"
#W="2p40"
# TODO
#Q2="2p1"
#W="2p95"
# ONGOING **
#Q2="3p0"
#W="2p32"
# DONE
#Q2="3p0"
#W="3p14"
# DONE
#Q2="4p4"
#W="2p74"
# TODO
Q2="5p5"
W="3p02"

KIN="Q${Q2}W${W}"

TARGET=("LH2" "dummy")
EPS=("high" "low")
PHISET=("center" "left" "right")

#TARGET=("dummy")
#EPS=("high" "low")
#PHISET=("center" "left" "right")

for t in "${TARGET[@]}"; do
    for e in "${EPS[@]}"; do
        for p in "${PHISET[@]}"; do

	    if [ $t = "dummy" ]; then
		file_name="${KIN}${p}_${e}e_dummy"
	    else
		file_name="${KIN}${p}_${e}e"
	    fi

	    numbers_to_match=()
	    IFS=', ' read -r -a numbers_to_match <<< "$( grab_runs ${file_name} )"
	    echo
	    echo "${file_name}"
	    echo "Run Numbers: [${numbers_to_match[@]}]"

            # Check if numbers_to_match is empty and break the loop
            if [ ${#numbers_to_match[@]} -eq 0 ]; then
                echo "No run numbers to process. Exiting loop."
                break
            fi
	    
	    inpFile="${UTILPATH}/run_list/${ANATYPE}LT/${file_name}"

	    replay_root_path="${ROOTPATH}/${ANATYPE}LT"
	    while true; do
		# Prompt for confirmation before proceeding
		#read -p "Are you sure you want to remove files with specified numbers in the filename? (yes/no): " answer
		answer="yes"

		case "$answer" in
		    [Yy]* )
			# Finds number of lines of inpFile
			numlines=$(eval "wc -l < ${inpFile}")
			# Loop through each number in the list
			for number in "${numbers_to_match[@]}"
			do
			    echo
			    echo "Running ${number}"
			    cd "${LTANAPATH}"
			    rootfile=/cache/hallc/kaonlt/Pass3_Dec_2023/ROOTfiles/Analysis/KaonLT/Kaon_coin_replay_production_${number}_-1.root
				jcache get ${rootfile}
			    ./applyCuts_Prod.sh -pm ${e} ${p} ${Q2} ${W} ${t} ${number} kaon
				./applyCuts_Prod.sh -p ${e} ${p} ${Q2} ${W} ${t} ${number} pion
			done
			break ;;
		    [Nn]* ) 
			echo "Operation aborted."
			break ;;
		    * ) 
			echo "Please answer yes or no." ;;
		esac
	    done
	done
    done
done
