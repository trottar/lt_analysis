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
PATHFILE_INFO="$(python3 "${LTSEP_FIELDS_SCRIPT}" "${REPO_ROOT}")"
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

# Flag definitions (flags: h, c, s)
while getopts 'hcs' flag; do
    case "${flag}" in
        h) 
        echo "----------------------------------------------------------"
        echo "./check_Offsets.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: Check the effect of offsets to dW, dE_m, dp_m(par), dp_m(perp)"
        echo "----------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
	echo "    If no flags called arguments are..."
	echo "        coin -> KIN=arg1"
	echo "        sing -> SPEC=arg1 KIN=arg2 (requires -s flag)"
        echo "    -h, help"
        echo "    -c, compile fortran code"
	echo "        coin -> KIN=arg1"
	echo "        sing -> SPEC=arg1 KIN=arg2 (requires -s flag)"	
	echo "    -s, single arm"
        exit 0
        ;;
        c) c_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

HEEPFOR="heepcheck"

cd ${LTANAPATH}/src/HeeP
# When any flag is used then the user input changes argument order
if [[ $c_flag = "true" && $s_flag = "true" ]]; then
    echo "Compiling ${HEEPFOR}.f..."
    eval "gfortran -o  ${HEEPFOR} ${HEEPFOR}.f"
    spec=$2
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    KIN=$3
    InputSIMC="Heep_${SPEC}_${KIN}"
elif [[ $c_flag = "true" && $s_flag != "true" ]]; then
    echo "Compiling ${HEEPFOR}.f..."
    eval "gfortran -o  ${HEEPFOR} ${HEEPFOR}.f"
    KIN=$2
    InputSIMC="Heep_Coin_${KIN}"
elif [[ $s_flag = "true" ]]; then
    spec=$2
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    KIN=$3
    InputSIMC="Heep_${SPEC}_${KIN}"
else
    KIN=$1
    InputSIMC="Heep_Coin_${KIN}"
fi

cd ${LTANAPATH}/src/setup
# Python script that gets current values of simc input file
SIMCINP=`python3 getSetting.py ${InputSIMC}`

# From getSetting.py define variables for beam and theta to
# be used as inputs for fortran elastics script
BEAMINP=`echo ${SIMCINP} | cut -d ',' -f1`
THETAINP=`echo ${SIMCINP} | cut -d ',' -f2`

cd ${LTANAPATH}/src/HeeP
# Runs fortran code using 'expect' which takes the user input
# value then runs the heepcheck script to check offsets
# (Fotran script is run in background)
./${HEEPFOR}.expect $(echo "${BEAMINP}*1000"|bc) ${THETAINP} # piping bc allows float arithmetic
