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
SUBMIT_SCRIPT="${SCRIPT_DIR}/farm_env/submit_replay.py"
REBALANCE_SCRIPT="${SCRIPT_DIR}/farm_env/rebalance_swif.py"
WORKFLOW_PREFIX="kaonlt"

# Flag definitions (flags: h, s, r, a, n, w)
while getopts 'hsranw:' flag; do
    case "${flag}" in
        h)
        echo "--------------------------------------------------------------"
        echo "./run_farm.sh -{flags} {variable arguments, see help}"
        echo
        echo "Description: Submit or rebalance KaonLT replay SWIF2 workflows"
        echo "             for a Q2/W family using the farm_env helpers."
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called..."
        echo "    Q2=arg1, W=arg2"
        echo "    -h, help"
        echo "    -s, actually submit jobs via submit_replay.py (default is dry-run)"
        echo "    -r, rebalance an existing workflow instead of submitting jobs"
        echo "    -a, with -r, actually apply modify-jobs commands (default is dry-run)"
        echo "    -n, do not call 'swif2 run' after submit/rebalance"
        echo "    -w, override workflow name"
        echo
        echo "Notes..."
        echo "    Submit mode automatically scans all matching JSON variants in input/kaon"
        echo "    for the requested Q2/W family, so all targets/settings are included."
        echo
        echo "Examples..."
        echo "    ./run_farm.sh 3p0 3p14"
        echo "    ./run_farm.sh -s 3p0 3p14"
        echo "    ./run_farm.sh -r 3p0 3p14"
        echo "    ./run_farm.sh -r -a -n 3p0 3p14"
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
        w) workflow_override="${OPTARG}" ;;
        *)
        exit 1
        ;;
    esac
done

shift $((OPTIND - 1))

if [[ ! -f "${SUBMIT_SCRIPT}" ]]; then
    echo "Submit helper not found: ${SUBMIT_SCRIPT}"
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
if [[ -n "${workflow_override}" ]]; then
    WORKFLOW="${workflow_override}"
else
    WORKFLOW="${WORKFLOW_PREFIX}_${FAMILY}_${USER_NAME}"
fi

if [[ "${rebalance_flag}" = "true" ]]; then
    echo
    echo "---------------------------------------------------------"
    echo
    echo "Inspecting SWIF2 workflow for Q2=${Q2}, W=${W}..."
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
echo "Preparing replay submission for Q2=${Q2}, W=${W}..."
echo
echo "                        WORKFLOW = ${WORKFLOW}"
echo
echo "---------------------------------------------------------"
echo

submit_cmd=(python3 "${SUBMIT_SCRIPT}" "${Q2}" "${W}" --workflow-name "${WORKFLOW}")
if [[ "${submit_flag}" = "true" ]]; then
    submit_cmd+=(--submit)
fi
if [[ "${no_run_flag}" = "true" ]]; then
    submit_cmd+=(--no-run)
fi

echo "Running: ${submit_cmd[*]}"
"${submit_cmd[@]}"
