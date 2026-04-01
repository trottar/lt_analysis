#! /bin/bash

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

# Define global variables for lt_analysis scripts
Q2="4p4"
W="2p74"
TMIN=0.400
TMAX=0.750
# Combined KaonLT 4p4+3p0 (2p74+2p32)
#Q2="4p3"
#W="2p70"
#TMIN=0.400
#TMAX=0.750
ParticleType="kaon"
#
#Q2="1p6"
#W="2p22"
#TMIN=0.001
#TMAX=0.300
#
#Q2="2p4"
#W="2p22"
#TMIN=0.100
#TMAX=0.600
# Combined Fpi-2
#Q2="2p0"
#W="2p22"
#TMIN=0.001
#TMAX=0.600
#ParticleType="pion"

POL="+1" # Positive polarity

DEBUG="False" # Flag for no plot splash
#DEBUG="True" # Flag for plot splash

KIN="Q${Q2}W${W}"
OutFilename="Testing_${KIN}"
# - `date` is a command that prints or sets the system date and time.
# - `+%H` extracts the hour in 24-hour format.
# - `+%M` extracts the minute.
# - `+%S` extracts the second.
# - `+%Y` extracts the year.
# - `+%B` extracts the full month name.
# - `%d` extracts the day of the month.
formatted_date=$(date +%Y%B%d_H%HM%MS%S) 

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Running fitting algorithm for ${ParticleType} at Kinematics: ${KIN} with POL: ${POL}"
echo "Output will be saved to: ${OutFilename}"
echo "Current time and date: ${formatted_date}"
echo "DEBUG mode is: ${DEBUG}"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo

cd "${LTANAPATH}/testing"
python3 main.py ${KIN} ${W} ${Q2} ${TMIN} ${TMAX} ${ParticleType} ${POL} ${OutFilename} ${formatted_date} ${DEBUG}
