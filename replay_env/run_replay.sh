#! /bin/bash
###########################################################################################################################
# Created - 20/July/21, Author - Muhammad Junaid (mjo147@uregina.ca), University of Regina, Canada (Copyright (c) junaid) #
# 28/11/21 - Version 2 - Utilises new ltsep package by Richard Trotta
# 19/01/22 - Version 3 - Target type is now an argument
# 20/03/26 - Version 4 - RLT, moved to lt_analysis and slimmed unnecessary parts
###########################################################################################################################
# To run this script, execute ./scriptname $RUNNUMBER$

#################################################################################################################################################

echo "Starting analysis of Kaon events"
echo "I take as arguments the run number and max number of events!"
# Input params - run number and max number of events
RUNNUMBER=$1
if [[ -z "$1" ]]; then
    echo "I need an input run number"
    echo "Please provide a run number as input"
fi
MAXEVENTS=-1

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

# #################################################################################################################################################

normalize_ltsep_dir() {
    local base_dir="$1"
    local leaf="${ANATYPE}LT"
    if [[ "${base_dir}" == *"/${leaf}" ]]; then
        printf '%s\n' "${base_dir}"
    elif [[ "${base_dir}" == *"None"* ]]; then
        printf '%s\n' "${base_dir/None/${leaf}}"
    else
        printf '%s\n' "${base_dir}/${leaf}"
    fi
}

REPLAY_OUTPUT_DIR="$(normalize_ltsep_dir "${ROOTPATH}")"
SCALER_OUTPUT_DIR="${UTILPATH}/ROOTfiles/Scalers"
SCALER_OUTPUT_FILE="${SCALER_OUTPUT_DIR}/coin_${ANATYPE}LT_replay_scalers_${RUNNUMBER}_${MAXEVENTS}.root"
BCM_PARAM_FILE="bcmcurrent_${RUNNUMBER}_.param"

# Source stuff depending upon hostname. Change or add more as needed  
if [[ "${HOST}" = *"farm"* ]]; then
    if [[ "${HOST}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
	source /apps/root/6.18.04/setroot_CUE.bash
    fi
    cd "$HCANAPATH"
    source "$HCANAPATH/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOST}" = *"qcd"* ]]; then
    source "$REPLAYPATH/setup.sh" 
fi

cd $REPLAYPATH
echo "Replay path is $REPLAYPATH"
# ###################################################################################################################################################
###################################################################################################################################################
mkdir -p "${SCALER_OUTPUT_DIR}"

if [ ! -f "${SCALER_OUTPUT_FILE}" ]; then
    if ! "${REPLAYPATH}/hcana" -l -q -b "${REPLAYPATH}/SCRIPTS/COIN/SCALERS/replay_${ANATYPE}LT_scalers.C(${RUNNUMBER},${MAXEVENTS})"; then
        echo "ERROR: scaler replay failed for run ${RUNNUMBER}"
        exit 1
    fi
    if [ ! -f "${SCALER_OUTPUT_FILE}" ]; then
        echo "ERROR: scaler replay did not create ${SCALER_OUTPUT_FILE}"
        exit 1
    fi
    cd "$REPLAYPATH"
    root -b -l<<EOF 
.L ${REPLAYPATH}/ScalerCalib.C
.x ${REPLAYPATH}/run.C("${SCALER_OUTPUT_FILE}")
.q  
EOF
    if [ ! -f "${BCM_PARAM_FILE}" ]; then
        echo "ERROR: bcm calibration did not create ${BCM_PARAM_FILE}"
        exit 1
    fi
    mv "${BCM_PARAM_FILE}" "$REPLAYPATH/PARAM/HMS/BCM/CALIB/bcmcurrent_$RUNNUMBER.param"
    cd $REPLAYPATH
else echo "Scaler replayfile already found for this run in ${SCALER_OUTPUT_DIR} - Skipping scaler replay step"
fi

sleep 3

if [ ! -f "${REPLAY_OUTPUT_DIR}/${ANATYPE}_coin_replay_production_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
			"${REPLAYPATH}/hcana" -l -q -b "${REPLAYPATH}/SCRIPTS/COIN/PRODUCTION/FullReplay_${ANATYPE}LT_Phys_Prod.C(${RUNNUMBER},${MAXEVENTS})" | tee $UTILPATH/REPORT_OUTPUT/Analysis/${ANATYPE}LT/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
			replay_rc=${PIPESTATUS[0]}
			if [ "${replay_rc}" -ne 0 ]; then
			    echo "ERROR: full replay failed for run ${RUNNUMBER}"
			    exit "${replay_rc}"
			fi
else echo "Replayfile already found for this run in ${REPLAY_OUTPUT_DIR}/ - Skipping replay step"
fi
