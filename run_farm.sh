#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2026-03-30"
# ================================================================
#
# Author:  OpenAI Codex
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPLAY_SUBMIT_SCRIPT="${SCRIPT_DIR}/farm_env/submit_replay.py"
APPLYCUTS_SUBMIT_SCRIPT="${SCRIPT_DIR}/farm_env/submit_applycuts.py"
REBALANCE_SCRIPT="${SCRIPT_DIR}/farm_env/rebalance_swif.py"
WORKFLOW_PREFIX="kaonlt"
DEFAULT_MANIFEST_DIR="${SCRIPT_DIR}/input/kaon"

# Flag definitions (flags: h, s, r, a, n, c, w, m)
while getopts 'hsrancw:m:' flag; do
    case "${flag}" in
        h)
        echo "--------------------------------------------------------------"
        echo "./run_farm.sh -{flags} {variable arguments, see help}"
        echo
        echo "Description: Submit or rebalance KaonLT replay or applyCuts"
        echo "             SWIF2 workflows for a Q2/W family."
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called..."
        echo "    Q2=arg1, W=arg2"
        echo "    -h, help"
        echo "    -s, actually submit jobs (default is dry-run)"
        echo "    -r, rebalance an existing workflow instead of submitting jobs"
        echo "    -a, with -r, actually apply modify-jobs commands (default is dry-run)"
        echo "    -n, do not call 'swif2 run' after submit/rebalance"
        echo "    -c, use applyCuts mode instead of replay mode"
        echo "    -w, override workflow name"
        echo "    -m, override manifest directory (default: input/kaon)"
        echo
        echo "Notes..."
        echo "    Replay mode scans all matching JSON variants in the selected manifest"
        echo "    directory and submits"
        echo "    one replay job per unique run."
        echo "    applyCuts mode scans the same selected manifest directory and"
        echo "    submits one job per"
        echo "    manifest variant + run, but only when replay output exists."
        echo
        echo "Examples..."
        echo "    ./run_farm.sh 3p0 3p14"
        echo "    ./run_farm.sh -s 3p0 3p14"
        echo "    ./run_farm.sh -c 3p0 3p14"
        echo "    ./run_farm.sh -c -s 3p0 3p14"
        echo "    ./run_farm.sh -m input/kaon_test -s 3p0 3p14"
        echo "    ./run_farm.sh -r 3p0 3p14"
        echo "    ./run_farm.sh -r -c -a -n 3p0 3p14"
        echo
        echo "Available Kinematics..."
        echo "                      Q2=5p5, W=3p02"
        echo "                      Q2=4p4, W=2p74"
        echo "                      Q2=3p0, W=3p14"
        echo "                      Q2=3p0, W=2p32"
        echo "                      Q2=2p1, W=2p95"
        echo "                      Q2=0p5, W=2p40"
        exit 0
        ;;
        s) submit_flag='true' ;;
        r) rebalance_flag='true' ;;
        a) apply_flag='true' ;;
        n) no_run_flag='true' ;;
        c) cuts_flag='true' ;;
        w) workflow_override="${OPTARG}" ;;
        m) manifest_dir="${OPTARG}" ;;
        *)
        exit 1
        ;;
    esac
done

shift $((OPTIND - 1))

if [[ ! -f "${REPLAY_SUBMIT_SCRIPT}" ]]; then
    echo "Replay submit helper not found: ${REPLAY_SUBMIT_SCRIPT}"
    exit 1
fi

if [[ ! -f "${APPLYCUTS_SUBMIT_SCRIPT}" ]]; then
    echo "applyCuts submit helper not found: ${APPLYCUTS_SUBMIT_SCRIPT}"
    exit 1
fi

if [[ ! -f "${REBALANCE_SCRIPT}" ]]; then
    echo "Rebalance helper not found: ${REBALANCE_SCRIPT}"
    exit 1
fi

if [[ "${rebalance_flag}" = "true" && "${submit_flag}" = "true" ]]; then
    echo "Please choose either submit mode or rebalance mode, not both."
    exit 1
fi

if [[ "${rebalance_flag}" != "true" && "${apply_flag}" = "true" ]]; then
    echo "The -a flag is only valid together with -r."
    exit 1
fi

Q2=$1
W=$2

echo "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5]"
echo "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40]"
if [[ -z "${Q2}" || ! "${Q2}" =~ ^(5p5|4p4|3p0|2p1|0p5)$ ]]; then
    echo ""
    echo "I need a valid Q2..."
    while true; do
        echo ""
        read -p "Q2 must be one of - [5p5 - 4p4 - 3p0 - 2p1 - 0p5] - or press ctrl-c to exit : " Q2
        case "${Q2}" in
            '');;
            '5p5'|'4p4'|'3p0'|'2p1'|'0p5') break;;
        esac
    done
fi
if [[ -z "${W}" || ! "${W}" =~ ^(3p02|2p74|3p14|2p32|2p95|2p40)$ ]]; then
    echo ""
    echo "I need a valid W..."
    while true; do
        echo ""
        read -p "W must be one of - [3p02 - 2p74 - 3p14 - 2p32 - 2p95 - 2p40] - or press ctrl-c to exit : " W
        case "${W}" in
            '');;
            '3p02'|'2p74'|'3p14'|'2p32'|'2p95'|'2p40') break;;
        esac
    done
fi

FAMILY="Q${Q2}W${W}"
USER_NAME="${USER:-user}"
if [[ "${cuts_flag}" = "true" ]]; then
    MODE="applyCuts"
    MODE_SCRIPT="${APPLYCUTS_SUBMIT_SCRIPT}"
    WORKFLOW_SUFFIX="_applycuts_${USER_NAME}"
else
    MODE="replay"
    MODE_SCRIPT="${REPLAY_SUBMIT_SCRIPT}"
    WORKFLOW_SUFFIX="_${USER_NAME}"
fi

if [[ -n "${workflow_override}" ]]; then
    WORKFLOW="${workflow_override}"
else
    WORKFLOW="${WORKFLOW_PREFIX}_${FAMILY}${WORKFLOW_SUFFIX}"
fi

if [[ -z "${manifest_dir}" ]]; then
    manifest_dir="${DEFAULT_MANIFEST_DIR}"
fi

if [[ "${rebalance_flag}" = "true" ]]; then
    echo
    echo "---------------------------------------------------------"
    echo
    echo "Inspecting ${MODE} SWIF2 workflow for Q2=${Q2}, W=${W}..."
    echo
    echo "                        WORKFLOW = ${WORKFLOW}"
    echo
    echo "---------------------------------------------------------"
    echo

    rebalance_cmd=(python3 "${REBALANCE_SCRIPT}" "${WORKFLOW}")
    if [[ "${apply_flag}" = "true" ]]; then
        rebalance_cmd+=(--apply)
    fi
    if [[ "${no_run_flag}" = "true" ]]; then
        rebalance_cmd+=(--no-run)
    fi

    echo "Running: ${rebalance_cmd[*]}"
    "${rebalance_cmd[@]}"
    exit $?
fi

echo
echo "---------------------------------------------------------"
echo
echo "Preparing ${MODE} submission for Q2=${Q2}, W=${W}..."
echo
echo "                        WORKFLOW = ${WORKFLOW}"
echo "                    MANIFEST DIR = ${manifest_dir}"
echo
echo "---------------------------------------------------------"
echo

submit_cmd=(python3 "${MODE_SCRIPT}" "${Q2}" "${W}" --workflow-name "${WORKFLOW}" --manifest-dir "${manifest_dir}")
if [[ "${submit_flag}" = "true" ]]; then
    submit_cmd+=(--submit)
fi
if [[ "${no_run_flag}" = "true" ]]; then
    submit_cmd+=(--no-run)
fi

echo "Running: ${submit_cmd[*]}"
"${submit_cmd[@]}"
