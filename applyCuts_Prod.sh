#! /bin/bash

JOB_LAUNCH_DIR="${SWIF_JOB_WORK_DIR:-${SWIF_JOB_STAGE_DIR:-}}"
RUNNING_UNDER_SWIF=0
if [[ -n "${JOB_LAUNCH_DIR}" ]]; then
    RUNNING_UNDER_SWIF=1
fi

# Runs a repo-local ltsep wrapper so batch jobs do not depend on upstream
# getPathDict.py calling os.getlogin() on worker nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
if [[ ! -f "${LTSEP_FIELDS_SCRIPT}" ]]; then
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
fi

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
    echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4]"
    echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20]"
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
    if [[ -z "$4" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5|0p4 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid Q2..."
	while true; do
	    echo ""
	    read -p "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4] - or press ctrl-c to exit : " Q2
	    case $Q2 in
		'');; # If blank, prompt again
		'5p5'|'4p4'|'3p0'|'2p1'|'0p5'|'0p4') break;; # If a valid option, break the loop and continue
	    esac
	done
    fi
    if [[ -z "$5" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40|2p20 ]]; then # Check the 3rd argument was provided and that it's one of the valid options
	echo ""
	echo "I need a valid W..."
	while true; do
	    echo ""
	    read -p "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20] - or press ctrl-c to exit : " W
	    case $W in
		'');; # If blank, prompt again
		'3p02'|'2p74'|'3p14'|'2p32'|'2p95'|'2p40'|'2p20') break;; # If a valid option, break the loop and continue
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
    echo "Running for kaon and pion"
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

if [[ $TargetType = "dummy" ]]; then
    DataType="dummy"
else
    DataType="data"
fi

cd "${LTANAPATH}/src/setup"
SKIM_OUTPUT_DIR="$(normalize_ltsep_dir "${SKIMPATH}")"
mkdir -p "${SKIM_OUTPUT_DIR}"
PATH_DIAGNOSTICS_DONE=0

resolve_real_path() {
    local input_path="$1"
    local resolved_path=""
    if resolved_path="$(readlink -f -- "${input_path}" 2>/dev/null)" && [[ -n "${resolved_path}" ]]; then
        printf '%s\n' "${resolved_path}"
    else
        printf '%s\n' "${input_path}"
    fi
}

build_replay_input_dir() {
    local analysis_leaf="${ANATYPE}LT"
    local base_path="${1%/}"

    case "${base_path}" in
        */Analysis/"${analysis_leaf}")
            printf '%s\n' "${base_path}"
            ;;
        */Analysis)
            printf '%s/%s\n' "${base_path}" "${analysis_leaf}"
            ;;
        */ROOTfiles)
            printf '%s/Analysis/%s\n' "${base_path}" "${analysis_leaf}"
            ;;
        */None)
            printf '%s/Analysis/%s\n' "${base_path%/None}" "${analysis_leaf}"
            ;;
        *)
            printf '%s/Analysis/%s\n' "${base_path}" "${analysis_leaf}"
            ;;
    esac
}

build_replay_input_file() {
    local replay_basename="${ANATYPE}_coin_replay_production_${RUNNUM}_-1.root"
    local replay_input_dir
    replay_input_dir="$(build_replay_input_dir "$1")"
    printf '%s/%s\n' "${replay_input_dir}" "${replay_basename}"
}

cache_path_to_mss_path() {
    local cache_path="$1"
    local suffix=""
    case "${cache_path}" in
        /cache/*)
            suffix="${cache_path#/cache/}"
            ;;
        /lustre*/expphy/cache/*)
            suffix="${cache_path#*/expphy/cache/}"
            ;;
        *)
            return 1
            ;;
    esac
    printf '/mss/%s\n' "${suffix}"
}

print_applycuts_path_diagnostics() {
    local replay_input_file="$1"
    local replay_input_real="$2"
    local replay_input_dir
    local replay_input_real_dir
    replay_input_dir="$(dirname "${replay_input_file}")"
    replay_input_real_dir="$(dirname "${replay_input_real}")"
    echo "Replay ROOTPATH   : ${ROOTPATH}"
    echo "Replay root real  : $(resolve_real_path "${ROOTPATH}")"
    echo "Replay input dir  : ${replay_input_dir}"
    echo "Replay input real : ${replay_input_real_dir}"
    echo "Replay input file : ${replay_input_file}"
    echo "Replay file real  : ${replay_input_real}"
    echo "Skim base path    : ${SKIMPATH}"
    echo "Skim output dir   : ${SKIM_OUTPUT_DIR}"
    echo "Skim output real  : $(resolve_real_path "${SKIM_OUTPUT_DIR}")"
}

ensure_cache_backed_replay_input_ready() {
    local replay_input_real="$1"
    local replay_input_mss

    case "${replay_input_real}" in
        /cache/*|/lustre*/expphy/cache/*)
            if [[ -e "${replay_input_real}" && -s "${replay_input_real}" ]]; then
                echo "Cache-backed replay input is available at ${replay_input_real}"
                return 0
            fi
            if ! replay_input_mss="$(cache_path_to_mss_path "${replay_input_real}")"; then
                echo "ERROR: replay input resolves into cache, but MSS path could not be derived from ${replay_input_real}"
                return 2
            fi
            if [[ "${replay_input_mss}" != /mss/hallc/kaonlt/*/ROOTfiles/Analysis/*/*.root ]]; then
                echo "ERROR: derived MSS recovery path looks invalid: ${replay_input_mss}"
                return 2
            fi
            if ! command -v jcache >/dev/null 2>&1; then
                echo "ERROR: replay input resolves into cache, but jcache is not available to request ${replay_input_mss}"
                return 2
            fi
            echo "Replay input resolves into cache and is not currently available."
            echo "Resolved cache replay file : ${replay_input_real}"
            printf 'Derived MSS recovery file  : <%s>\n' "${replay_input_mss}"
            echo "Requesting recovery from tape with:"
            printf '  jcache get <%s>\n' "${replay_input_mss}"
            jcache get "${replay_input_mss}"
            local jcache_rc=$?
            if [[ "${jcache_rc}" -ne 0 ]]; then
                echo "ERROR: jcache get failed for ${replay_input_mss}"
                return "${jcache_rc}"
            fi
            echo "Replay cache recovery requested. Rerun applyCuts after staging completes."
            return 2
            ;;
        *)
            return 0
            ;;
    esac
}

stage_swif_output_copy() {
    local source_file="$1"
    local staged_file="${JOB_LAUNCH_DIR}/$(basename "${source_file}")"

    if [[ "${RUNNING_UNDER_SWIF}" -ne 1 ]]; then
        return 0
    fi

    if [[ ! -f "${source_file}" ]]; then
        echo "ERROR: expected skim output not found at ${source_file}"
        return 1
    fi

    cp -f "${source_file}" "${staged_file}"
    local copy_rc=$?
    if [[ "${copy_rc}" -ne 0 ]]; then
        echo "ERROR: failed to stage skim output into ${staged_file}"
        return "${copy_rc}"
    fi

    echo "Staged SWIF skim output copy at ${staged_file}"
    return 0
}

run_applycuts_particle() {
    local phi_label="$1"
    local particle="$2"
    local replay_input_file
    local replay_input_real
    replay_input_file="$(build_replay_input_file "${ROOTPATH}")"
    replay_input_real="$(build_replay_input_file "$(resolve_real_path "${ROOTPATH}")")"
    local out_f_file="${SKIM_OUTPUT_DIR}/${particle}_${RUNNUM}_-1_Raw_Data.root"
    local log_file="${LTANAPATH}/log/${phi_label}_${particle}_${RUNNUM}_${KIN}.log"

    if [[ "${PATH_DIAGNOSTICS_DONE}" -ne 1 ]]; then
        print_applycuts_path_diagnostics "${replay_input_file}" "${replay_input_real}"
        PATH_DIAGNOSTICS_DONE=1
    fi

    if [ -e "$out_f_file" ]; then
        echo "$out_f_file already exists. Skipping analysis for ${particle}."
        stage_swif_output_copy "$out_f_file"
        return $?
    fi

    ensure_cache_backed_replay_input_ready "${replay_input_real}"
    local cache_ready_rc=$?
    if [[ "${cache_ready_rc}" -ne 0 ]]; then
        return "${cache_ready_rc}"
    fi

    rm -f "$log_file"
    "${LTSEP_PYTHON}" Analysed_Prod.py "${RUNNUM}" "${particle}" "${ANATYPE}_coin_replay_production" |& tee -a "$log_file"
    local rc=${PIPESTATUS[0]}
    if [[ $rc -ne 0 ]]; then
        return $rc
    fi

    if [[ ! -f "$out_f_file" ]]; then
        echo "ERROR: expected skim output not found at ${out_f_file}"
        return 1
    fi

    stage_swif_output_copy "$out_f_file"
    return $?
}

if [[ $p_flag = "true" ]]; then

    # Define input and output file names
    InDATAFilename="Proc_Data_${ParticleType}_${KIN}.root"
    OutDATAFilename="Analysed_Data_${ParticleType}_${KIN}"
    OutFullAnalysisFilename="FullAnalysis_${ParticleType}_${KIN}"

    # The analysis script (Analysed_Prod.py) will create a new root file per run number
    if [ ${PHIVAL} = "right" ]; then
	echo
	echo "Analysing right ${DataType}..."
	echo
	echo
	echo "------------------------------------------------------"
	echo "Analysing right ${DataType} ${ParticleType} run ${RUNNUM}..."
	echo "------------------------------------------------------"
	echo
	run_applycuts_particle "Right" "${ParticleType}"
	rc=$?
	if [[ $rc -ne 0 ]]; then
	    exit $rc
	fi
	echo
    fi

    # Checks that array isn't empty
    if [ ${PHIVAL} = "left" ]; then
	echo
	echo "Analysing left ${DataType}..."
	echo
	echo
	echo "------------------------------------------------------"
	echo "Analysing left ${DataType} ${ParticleType} run ${RUNNUM}..."
	echo "------------------------------------------------------"
	echo
	run_applycuts_particle "Left" "${ParticleType}"
	rc=$?
	if [[ $rc -ne 0 ]]; then
	    exit $rc
	fi
    fi

    # Checks that array isn't empty
    if [ ${PHIVAL} = "center" ]; then
	echo
	echo "Analysing center ${DataType}..."
	echo
	echo
	echo "------------------------------------------------------"
	echo "Analysing center ${DataType} ${ParticleType} run ${RUNNUM}..."
	echo "------------------------------------------------------"
	echo
	run_applycuts_particle "Center" "${ParticleType}"
	rc=$?
	if [[ $rc -ne 0 ]]; then
	    exit $rc
	fi
    fi

else
    declare -a ParticleTypes=("kaon" "pion")
    for i in "${ParticleTypes[@]}"
    do 
	# Define input and output file names
	InDATAFilename="Proc_Data_${i}_${KIN}.root"
	OutDATAFilename="Analysed_Data_${i}_${KIN}"
	OutFullAnalysisFilename="FullAnalysis_${i}_${KIN}"

	# The analysis script (Analysed_Prod.py) will create a new root file per run number
	if [ ${PHIVAL} = "right" ]; then
	    echo
	    echo "Analysing right ${DataType}..."
	    echo
	    echo
	    echo "------------------------------------------------------"
	    echo "Analysing right ${DataType} ${i} run ${RUNNUM}..."
	    echo "------------------------------------------------------"
	    echo
	    run_applycuts_particle "Right" "${i}"
	    rc=$?
	    if [[ $rc -ne 0 ]]; then
		exit $rc
	    fi
	    echo
	fi

	# Checks that array isn't empty
	if [ ${PHIVAL} = "left" ]; then
	    echo
	    echo "Analysing left ${DataType}..."
	    echo
	    echo
	    echo "------------------------------------------------------"
	    echo "Analysing left ${DataType} ${i} run ${RUNNUM}..."
	    echo "------------------------------------------------------"
	    echo
	    run_applycuts_particle "Left" "${i}"
	    rc=$?
	    if [[ $rc -ne 0 ]]; then
		exit $rc
	    fi
	fi

	# Checks that array isn't empty
	if [ ${PHIVAL} = "center" ]; then
	    echo
	    echo "Analysing center ${DataType}..."
	    echo
	    echo
	    echo "------------------------------------------------------"
	    echo "Analysing center ${DataType} ${i} run ${RUNNUM}..."
	    echo "------------------------------------------------------"
	    echo
	    run_applycuts_particle "Center" "${i}"
	    rc=$?
	    if [[ $rc -ne 0 ]]; then
		exit $rc
	    fi
	fi
    done
fi

echo
echo
echo
echo "Script Complete!"
