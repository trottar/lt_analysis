#! /bin/bash

# Runs a repo-local ltsep wrapper so batch jobs do not depend on upstream
# getPathDict.py calling os.getlogin() on worker nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
if [[ ! -f "${LTSEP_FIELDS_SCRIPT}" ]]; then
    echo "ERROR: failed to locate ${LTSEP_FIELDS_SCRIPT}" >&2
    exit 1
fi

PATHFILE_INFO="$(python3 "${LTSEP_FIELDS_SCRIPT}" "${REPO_ROOT}")"
path_rc=$?
if [[ "${path_rc}" -ne 0 || -z "${PATHFILE_INFO}" ]]; then
    echo "ERROR: failed to resolve ltsep paths via ${LTSEP_FIELDS_SCRIPT}" >&2
    exit 1
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1`
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

sync_symlink() {
    local link_path="$1"
    local target_path="$2"
    local label="$3"
    local current_target=""

    if [ -L "${link_path}" ]; then
        current_target=$(readlink "${link_path}")
        if [ "${current_target}" = "${target_path}" ] && [ -e "${link_path}" ]; then
            return 0
        fi

        rm "${link_path}"
    elif [ -e "${link_path}" ]; then
        echo "ERROR: ${link_path} exists and is not a symlink" >&2
        exit 1
    fi

    ln -s "${target_path}" "${link_path}"
    echo "Updated ${label}: ${link_path}-->${target_path}"
}

resolve_physical_path() {
    local input_path="$1"

    if command -v readlink >/dev/null 2>&1; then
        readlink -f "${input_path}" 2>/dev/null && return 0
    fi

    cd "${input_path}" >/dev/null 2>&1 && pwd -P
}

normalize_toggle() {
    local value
    value=$(printf '%s' "${1}" | tr '[:upper:]' '[:lower:]')
    case "${value}" in
        true|1|yes|on)
            printf 'true\n'
            ;;
        false|0|no|off|'')
            printf 'false\n'
            ;;
        *)
            return 1
            ;;
    esac
}

append_background_if_enabled() {
    local background_name="$1"
    local toggle_value="$2"
    local normalized=""

    if ! normalized=$(normalize_toggle "${toggle_value}"); then
        echo "ERROR: invalid toggle value '${toggle_value}' for ${background_name}" >&2
        exit 1
    fi

    if [ "${normalized}" = "true" ]; then
        ACTIVE_BACKGROUNDS+=("${background_name}")
    fi
}

run_background_setting() {
    local background_name="$1"
    local input_dir="${LTANAPATH}/input/backgrounds/${background_name}"
    local output_dir="${BACKGROUND_OUTPUT_BASE}/${background_name}"
    local log_dir="${LTANAPATH}/log/background_samples/${background_name}"
    local phi_setting=""
    local input_name=""
    local input_file=""

    if [ ! -d "${input_dir}" ]; then
        echo "Skipping ${background_name}: input directory ${input_dir} does not exist"
        return 0
    fi

    mkdir -p "${output_dir}" "${log_dir}"
    sync_symlink "${BACKGROUND_SIMC_RUNTIME_PATH}/input" "${input_dir}/" "${background_name} input"
    sync_symlink "${BACKGROUND_SIMC_RUNTIME_PATH}/OUTPUTS" "${output_dir}/" "${background_name} output"
    sync_symlink "${BACKGROUND_SIMC_RUNTIME_PATH}/worksim" "${VOLATILEPATH}/worksim/" "${background_name} worksim"

    for phi_setting in right left center; do
        input_name="Prod_Coin_Q${Q2}W${W}${phi_setting}_${EPSILON}e"
        input_file="${input_dir}/${input_name}.inp"

        if [ ! -f "${input_file}" ]; then
            echo "Skipping missing input ${input_file}"
            continue
        fi

        echo
        echo "Running ${background_name} background for ${input_name}..."
        echo
        if SIMC_LOG_DIR="${log_dir}" "${LTANAPATH}/src/setup/run_simc_tree" "${input_name}" "${BACKGROUND_REACTION}"; then
            COMPLETED_RUNS+=("${background_name}:${input_name}")
        else
            echo "WARNING: ${background_name} background failed for ${input_name}; continuing"
            FAILED_RUNS+=("${background_name}:${input_name}")
        fi
    done
}

# Flag definitions
while getopts 'h' flag; do
    case "${flag}" in
        h)
        echo "---------------------------------------------------------------------"
        echo "./background_samples/set_BackgroundInput.sh {variable arguments, see help}"
        echo
        echo "Description: Runs the background SIMC samples for one epsilon setting"
        echo "---------------------------------------------------------------------"
        echo
        echo "The following arguments can be called for the production backgrounds..."
        echo "    EPSILON=arg1, Q2=arg2, W=arg3"
        echo "    -h, help"
        echo
        echo "Avaliable Kinematics..."
        echo "                      Q2=5p5, W=3p02"
        echo "                      Q2=4p4, W=2p74"
        echo "                      Q2=3p0, W=3p14"
        echo "                      Q2=3p0, W=2p32"
        echo "                      Q2=2p1, W=2p95"
        echo "                      Q2=0p5, W=2p40"
        echo "                      Q2=0p4, W=2p20 (Q2 = 0.38, pion only)"
        exit 0
        ;;
        *)
        exit 1
        ;;
    esac
done

EPSILON=$1
Q2=$2
W=$3

CONFIG_FILE="${SCRIPT_DIR}/background_samples.conf"
BACKGROUND_SAMPLE_DIR="${LTANAPATH}/background_samples"
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: missing configuration file ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

BACKGROUND_OUTPUT_BASE="${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/background_samples"
BACKGROUND_OUTPUT_LINK="${BACKGROUND_SAMPLE_DIR}/OUTPUTS"
BACKGROUND_SIMC_LINK="${LTANAPATH}/background_samples/simc_gfortran"
BACKGROUND_REACTION="${BACKGROUND_REACTION:-production}"

echo "Epsilon must be - high - low - Case Sensitive!"
echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4]"
echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20]"
if [[ -z "$1" || ! "$EPSILON" =~ high|low ]]; then
    echo ""
    echo "I need a valid epsilon..."
    while true; do
        echo ""
        read -p "Epsilon must be - high - low - Case Sensitive! - or press ctrl-c to exit : " EPSILON
        case $EPSILON in
            '');;
            'high'|'low') break;;
        esac
    done
fi
if [[ -z "$2" || ! "$Q2" =~ 5p5|4p4|3p0|2p1|0p5|0p4 ]]; then
    echo ""
    echo "I need a valid Q2..."
    while true; do
        echo ""
        read -p "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5 - 0p4] - or press ctrl-c to exit : " Q2
        case $Q2 in
            '');;
            '5p5'|'4p4'|'3p0'|'2p1'|'0p5'|'0p4') break;;
        esac
    done
fi
if [[ -z "$3" || ! "$W" =~ 3p02|2p74|3p14|2p32|2p95|2p40|2p20 ]]; then
    echo ""
    echo "I need a valid W..."
    while true; do
        echo ""
        read -p "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40 - 2p20] - or press ctrl-c to exit : " W
        case $W in
            '');;
            '3p02'|'2p74'|'3p14'|'2p32'|'2p95'|'2p40'|'2p20') break;;
        esac
    done
fi

mkdir -p "${BACKGROUND_OUTPUT_BASE}" "${LTANAPATH}/log/background_samples"
sync_symlink "${BACKGROUND_OUTPUT_LINK}" "${BACKGROUND_OUTPUT_BASE}/" "background output"

if [ -n "${BACKGROUND_SIMC_PATH}" ]; then
    sync_symlink "${BACKGROUND_SIMC_LINK}" "${BACKGROUND_SIMC_PATH}" "background simc"
fi

if [ ! -L "${BACKGROUND_SIMC_LINK}" ] || [ ! -e "${BACKGROUND_SIMC_LINK}" ]; then
    echo "ERROR: background SIMC symlink is missing or broken at ${BACKGROUND_SIMC_LINK}" >&2
    echo "       Set BACKGROUND_SIMC_PATH in ${CONFIG_FILE} and rerun ./set_SymLinks.sh" >&2
    exit 1
fi

BACKGROUND_SIMC_RUNTIME_PATH="$(resolve_physical_path "${BACKGROUND_SIMC_LINK}")"
if [ -z "${BACKGROUND_SIMC_RUNTIME_PATH}" ] || [ ! -d "${BACKGROUND_SIMC_RUNTIME_PATH}" ]; then
    echo "ERROR: failed to resolve a physical background SIMC path from ${BACKGROUND_SIMC_LINK}" >&2
    exit 1
fi

export SIMC_OVERRIDE_PATH="${BACKGROUND_SIMC_RUNTIME_PATH}"
# Background SIMC uses its own model definitions, so skip the production
# physics_iterate.f/par.pl update helpers used by set_ProdInput.sh.

declare -a ACTIVE_BACKGROUNDS=()
declare -a COMPLETED_RUNS=()
declare -a FAILED_RUNS=()
append_background_if_enabled "neutron" "${RUN_NEUTRON:-true}"
append_background_if_enabled "delta" "${RUN_DELTA:-true}"
append_background_if_enabled "sidis" "${RUN_SIDIS:-true}"

if [ "${#ACTIVE_BACKGROUNDS[@]}" -eq 0 ]; then
    echo "No backgrounds are enabled in ${CONFIG_FILE}"
    exit 0
fi

for background_name in "${ACTIVE_BACKGROUNDS[@]}"; do
    run_background_setting "${background_name}"
done

echo
echo "Completed ${#COMPLETED_RUNS[@]} background sample runs for Q2=${Q2}, W=${W}, ${EPSILON} epsilon"

if [ "${#FAILED_RUNS[@]}" -gt 0 ]; then
    echo
    echo "The following background sample runs failed:"
    for failed_run in "${FAILED_RUNS[@]}"; do
        echo "  ${failed_run}"
    done
    exit 1
fi

echo "Background sample generation complete"
