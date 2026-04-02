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
echo "I take as arguments the run number!"
# Input params - run number and max number of events
RUNNUMBER=$1
if [[ -z "$1" ]]; then
    echo "I need an input run number"
    echo "Please provide a run number as input"
fi
MAXEVENTS=-1
JOB_LAUNCH_DIR="$(pwd)"

# Runs a repo-local ltsep wrapper so batch jobs do not depend on upstream
# getPathDict.py calling os.getlogin() on worker nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"

resolve_ltsep_python() {
    local shell_user="${USER:-$(id -un 2>/dev/null)}"
    local -a candidates=()
    local import_probe='import importlib; m=importlib.import_module("ltsep"); getattr(m, "Root", None) or getattr(importlib.import_module("ltsep.ltsep"), "Root")'

    if [[ -n "${LTSEP_PYTHON:-}" ]]; then
        candidates+=("${LTSEP_PYTHON}")
    fi
    if [[ -n "${REPLAY_LT_ENV:-}" ]]; then
        candidates+=("${REPLAY_LT_ENV}/bin/python3" "${REPLAY_LT_ENV}/bin/python")
    fi
    if [[ -n "${replay_lt_env:-}" ]]; then
        candidates+=("${replay_lt_env}/bin/python3" "${replay_lt_env}/bin/python")
    fi
    if [[ -n "${shell_user}" ]]; then
        candidates+=(
            "/u/group/c-kaonlt/USERS/${shell_user}/replay_lt_env/bin/python3"
            "/u/group/c-kaonlt/USERS/${shell_user}/replay_lt_env/bin/python"
            "/group/c-kaonlt/USERS/${shell_user}/replay_lt_env/bin/python3"
            "/group/c-kaonlt/USERS/${shell_user}/replay_lt_env/bin/python"
        )
    fi
    candidates+=("python3" "python" "/usr/bin/python3" "/usr/bin/python")

    local candidate resolved
    for candidate in "${candidates[@]}"; do
        if [[ "${candidate}" == */* ]]; then
            [[ -x "${candidate}" ]] || continue
            if "${candidate}" -c "${import_probe}" >/dev/null 2>&1; then
                printf '%s\n' "${candidate}"
                return 0
            fi
        else
            resolved="$(command -v "${candidate}" 2>/dev/null)" || continue
            if [[ -n "${resolved}" ]] && "${resolved}" -c "${import_probe}" >/dev/null 2>&1; then
                printf '%s\n' "${resolved}"
                return 0
            fi
        fi
    done
    return 1
}

LTSEP_PYTHON="$(resolve_ltsep_python)"
if [[ -z "${LTSEP_PYTHON}" ]]; then
    echo "ERROR: could not find a Python interpreter that imports ltsep" >&2
    exit 1
fi

PATHFILE_INFO="$("${LTSEP_PYTHON}" "${LTSEP_FIELDS_SCRIPT}" "${REPO_ROOT}")"
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
SCALER_REPORT_DIR="${UTILPATH}/REPORT_OUTPUT/Scalers"
SCALER_REPORT_FILE="${SCALER_REPORT_DIR}/${ANATYPE}_output_coin_scalers_Summary_${RUNNUMBER}_${MAXEVENTS}.report"
SCALER_MACRO_FILE="${REPLAYPATH}/SCRIPTS/COIN/SCALERS/replay_${ANATYPE}LT_coin_scalers.C"
FULL_REPLAY_MACRO_FILE="${REPLAYPATH}/SCRIPTS/COIN/PRODUCTION/FullReplay_${ANATYPE}LT_Phys_Prod.C"
BCM_PARAM_FILE="bcmcurrent_${RUNNUMBER}_.param"
BCM_CALIB_DIR="${REPLAYPATH}/CALIBRATION/bcm_current_map"
FULL_REPLAY_REPORT_DIR="${UTILPATH}/REPORT_OUTPUT/Analysis/${ANATYPE}LT"
FULL_REPLAY_REPORT_FILE="${FULL_REPLAY_REPORT_DIR}/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report"
FULL_REPLAY_OUTPUT_FILE="${REPLAY_OUTPUT_DIR}/${ANATYPE}_coin_replay_production_${RUNNUMBER}_${MAXEVENTS}.root"
REPORT_TARBALL_BASENAME="${SWIF_REPORT_TARBALL_BASENAME:-FullReplay_run_${RUNNUMBER}_${ANATYPE}LT.tar}"
SWIF_OUTPUT_FILE="${JOB_LAUNCH_DIR}/$(basename "${FULL_REPLAY_OUTPUT_FILE}")"
SWIF_REPORT_OUTPUT_FILE="${JOB_LAUNCH_DIR}/$(basename "${FULL_REPLAY_REPORT_FILE}")"
SWIF_REPORT_TARBALL_FILE="${JOB_LAUNCH_DIR}/${REPORT_TARBALL_BASENAME}"

# Source farm environment when available.
if [[ -f /site/12gev_phys/softenv.sh ]]; then
    source /site/12gev_phys/softenv.sh 2.3
fi
if [[ -f /apps/root/6.18.04/setroot_CUE.bash ]]; then
    source /apps/root/6.18.04/setroot_CUE.bash
fi

cd "$REPLAYPATH" || exit 1
if [ ! -f "${REPLAYPATH}/setup.sh" ]; then
    echo "ERROR: replay setup script not found at ${REPLAYPATH}/setup.sh"
    exit 1
fi
source "${REPLAYPATH}/setup.sh"

if [ ! -x "${REPLAYPATH}/hcana" ]; then
    echo "ERROR: hcana not found or not executable at ${REPLAYPATH}/hcana"
    exit 1
fi
if [ ! -f "${SCALER_MACRO_FILE}" ]; then
    echo "ERROR: scaler macro not found at ${SCALER_MACRO_FILE}"
    exit 1
fi
if [ ! -f "${FULL_REPLAY_MACRO_FILE}" ]; then
    echo "ERROR: full replay macro not found at ${FULL_REPLAY_MACRO_FILE}"
    exit 1
fi
# ###################################################################################################################################################
###################################################################################################################################################
mkdir -p "${SCALER_OUTPUT_DIR}"
mkdir -p "${SCALER_REPORT_DIR}"
mkdir -p "${FULL_REPLAY_REPORT_DIR}"

if [ ! -f "${SCALER_OUTPUT_FILE}" ]; then
    "${REPLAYPATH}/hcana" -l -q -b "${SCALER_MACRO_FILE}(${RUNNUMBER},${MAXEVENTS})" |& tee "${SCALER_REPORT_FILE}"
    scaler_rc=${PIPESTATUS[0]}
    if [ "${scaler_rc}" -ne 0 ]; then
        echo "ERROR: scaler replay failed for run ${RUNNUMBER}"
        echo "ERROR: see ${SCALER_REPORT_FILE}"
        exit "${scaler_rc}"
    fi
    if [ ! -f "${SCALER_OUTPUT_FILE}" ]; then
        echo "ERROR: scaler replay did not create ${SCALER_OUTPUT_FILE}"
        echo "ERROR: see ${SCALER_REPORT_FILE}"
        exit 1
    fi
    if ! command -v root >/dev/null 2>&1; then
        echo "ERROR: root command not found after sourcing replay environment"
        exit 1
    fi
    if [ ! -d "${BCM_CALIB_DIR}" ]; then
        echo "ERROR: bcm calibration directory not found at ${BCM_CALIB_DIR}"
        exit 1
    fi
    if [ ! -f "${BCM_CALIB_DIR}/ScalerCalib.C" ]; then
        echo "ERROR: bcm calibration macro not found at ${BCM_CALIB_DIR}/ScalerCalib.C"
        exit 1
    fi
    if [ ! -f "${BCM_CALIB_DIR}/run.C" ]; then
        echo "ERROR: bcm runner macro not found at ${BCM_CALIB_DIR}/run.C"
        exit 1
    fi
    cd "${BCM_CALIB_DIR}" || exit 1
    root -b -l <<EOF
.L ScalerCalib.C
.x run.C("${SCALER_OUTPUT_FILE}")
.q
EOF
    if [ ! -f "${BCM_PARAM_FILE}" ]; then
        echo "ERROR: bcm calibration did not create ${BCM_PARAM_FILE}"
        exit 1
    fi
    mv "${BCM_PARAM_FILE}" "$REPLAYPATH/PARAM/HMS/BCM/CALIB/bcmcurrent_$RUNNUMBER.param"
    cd "$REPLAYPATH" || exit 1
else echo "Scaler replayfile already found for this run in ${SCALER_OUTPUT_DIR} - Skipping scaler replay step"
fi

sleep 3

if [ ! -f "${REPLAY_OUTPUT_DIR}/${ANATYPE}_coin_replay_production_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
				"${REPLAYPATH}/hcana" -l -q -b "${FULL_REPLAY_MACRO_FILE}(${RUNNUMBER},${MAXEVENTS})" |& tee "${FULL_REPLAY_REPORT_FILE}"
					replay_rc=${PIPESTATUS[0]}
					if [ "${replay_rc}" -ne 0 ]; then
					    echo "ERROR: full replay failed for run ${RUNNUMBER}"
					    exit "${replay_rc}"
					fi
else echo "Replayfile already found for this run in ${REPLAY_OUTPUT_DIR}/ - Skipping replay step"
fi

if [ ! -f "${FULL_REPLAY_OUTPUT_FILE}" ]; then
    echo "ERROR: replay output file not found at ${FULL_REPLAY_OUTPUT_FILE}"
    exit 1
fi

if [ ! -f "${FULL_REPLAY_REPORT_FILE}" ]; then
    echo "ERROR: replay report file not found at ${FULL_REPLAY_REPORT_FILE}"
    exit 1
fi

shopt -s nullglob
REPORT_ARTIFACTS=( "${FULL_REPLAY_REPORT_DIR}"/*"${RUNNUMBER}_${MAXEVENTS}"* )
shopt -u nullglob
if [ "${#REPORT_ARTIFACTS[@]}" -eq 0 ]; then
    echo "ERROR: no replay report artifacts found for run ${RUNNUMBER} in ${FULL_REPLAY_REPORT_DIR}"
    exit 1
fi

REPORT_ARTIFACT_NAMES=()
for report_artifact in "${REPORT_ARTIFACTS[@]}"; do
    if [ -f "${report_artifact}" ]; then
        REPORT_ARTIFACT_NAMES+=( "$(basename "${report_artifact}")" )
    fi
done

if [ "${#REPORT_ARTIFACT_NAMES[@]}" -eq 0 ]; then
    echo "ERROR: no replay report files found for run ${RUNNUMBER} in ${FULL_REPLAY_REPORT_DIR}"
    exit 1
fi

cp -f "${FULL_REPLAY_OUTPUT_FILE}" "${SWIF_OUTPUT_FILE}"
copy_rc=$?
if [ "${copy_rc}" -ne 0 ]; then
    echo "ERROR: failed to stage replay output into ${SWIF_OUTPUT_FILE}"
    exit "${copy_rc}"
fi

cp -f "${FULL_REPLAY_REPORT_FILE}" "${SWIF_REPORT_OUTPUT_FILE}"
report_copy_rc=$?
if [ "${report_copy_rc}" -ne 0 ]; then
    echo "ERROR: failed to stage replay report into ${SWIF_REPORT_OUTPUT_FILE}"
    exit "${report_copy_rc}"
fi

tar -cf "${SWIF_REPORT_TARBALL_FILE}" -C "${FULL_REPLAY_REPORT_DIR}" "${REPORT_ARTIFACT_NAMES[@]}"
tar_rc=$?
if [ "${tar_rc}" -ne 0 ]; then
    echo "ERROR: failed to create replay report tarball at ${SWIF_REPORT_TARBALL_FILE}"
    exit "${tar_rc}"
fi

echo "Staged SWIF output copy at ${SWIF_OUTPUT_FILE}"
echo "Staged SWIF report copy at ${SWIF_REPORT_OUTPUT_FILE}"
echo "Staged SWIF report tarball at ${SWIF_REPORT_TARBALL_FILE}"
