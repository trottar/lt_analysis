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
JASMINE_SCRIPT="${SCRIPT_DIR}/farm_env/jasmine_put_from_manifest.py"
DIAGNOSE_SCRIPT="${SCRIPT_DIR}/farm_env/diagnose_swif_failures.py"
WORKFLOW_PREFIX="kaonlt"
DEFAULT_MANIFEST_DIR="${SCRIPT_DIR}/input/kaon"

# Flag definitions (flags: h, s, r, a, n, c, j, d, w, m, A, P)
while getopts 'hsrancjdw:m:A:P:' flag; do
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
        echo "    -j, use interactive Jasmine upload mode"
        echo "    -d, diagnose failed jobs in an existing workflow"
        echo "    -w, override workflow name"
        echo "    -m, override manifest directory (default: input/kaon)"
        echo "    -A, override SWIF/Slurm account (default from helper)"
        echo "    -P, override SWIF/Slurm partition (default from helper)"
        echo
        echo "Notes..."
        echo "    Replay mode scans all matching JSON variants in the selected manifest"
        echo "    directory and submits one replay job per unique run, first"
        echo "    requesting raw-file cache staging for missing coin_all data."
        echo "    applyCuts mode scans the same selected manifest directory and"
        echo "    submits one job per manifest variant + run, but only when the"
        echo "    replay ROOT file exists on MSS and is already present in cache"
        echo "    (requesting jcache staging first when needed)."
        echo "    Jasmine uploads are run separately from an interactive ifarm"
        echo "    session, but can be launched from this same wrapper with -j."
        echo "    In Jasmine mode, replay upload is the default and -c switches"
        echo "    the upload product kind to skim."
        echo
        echo "Examples..."
        echo "    ./run_farm.sh 3p0 3p14"
        echo "    ./run_farm.sh -s 3p0 3p14"
        echo "    ./run_farm.sh -c 3p0 3p14"
        echo "    ./run_farm.sh -c -s 3p0 3p14"
        echo "    ./run_farm.sh -m input/kaon_test -s 3p0 3p14"
        echo "    ./run_farm.sh -r 3p0 3p14"
        echo "    ./run_farm.sh -r -c -a -n 3p0 3p14"
        echo "    ./run_farm.sh -j -m input/kaon_test 3p0 3p14"
        echo "    ./run_farm.sh -j -c -m input/kaon_test -s 3p0 3p14"
        echo "    ./run_farm.sh -d 3p0 3p14"
        echo "    ./run_farm.sh -d -c 3p0 3p14"
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
        j) upload_flag='true' ;;
        d) diagnose_flag='true' ;;
        w) workflow_override="${OPTARG}" ;;
        m) manifest_dir="${OPTARG}" ;;
        A) account_override="${OPTARG}" ;;
        P) partition_override="${OPTARG}" ;;
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

if [[ ! -f "${JASMINE_SCRIPT}" ]]; then
    echo "Jasmine helper not found: ${JASMINE_SCRIPT}"
    exit 1
fi

if [[ ! -f "${DIAGNOSE_SCRIPT}" ]]; then
    echo "Diagnose helper not found: ${DIAGNOSE_SCRIPT}"
    exit 1
fi

if [[ "${rebalance_flag}" = "true" && "${submit_flag}" = "true" ]]; then
    echo "Please choose either submit mode or rebalance mode, not both."
    exit 1
fi

if [[ "${upload_flag}" = "true" && "${rebalance_flag}" = "true" ]]; then
    echo "Please choose either Jasmine upload mode or rebalance mode, not both."
    exit 1
fi

if [[ "${diagnose_flag}" = "true" && "${rebalance_flag}" = "true" ]]; then
    echo "Please choose either diagnose mode or rebalance mode, not both."
    exit 1
fi

if [[ "${diagnose_flag}" = "true" && "${upload_flag}" = "true" ]]; then
    echo "Please choose either diagnose mode or Jasmine upload mode, not both."
    exit 1
fi

if [[ "${rebalance_flag}" != "true" && "${apply_flag}" = "true" ]]; then
    echo "The -a flag is only valid together with -r."
    exit 1
fi

if [[ "${upload_flag}" = "true" && "${apply_flag}" = "true" ]]; then
    echo "The -a flag is not used in Jasmine upload mode."
    exit 1
fi

if [[ "${diagnose_flag}" = "true" && "${apply_flag}" = "true" ]]; then
    echo "The -a flag is not used in diagnose mode."
    exit 1
fi

if [[ "${upload_flag}" = "true" && "${no_run_flag}" = "true" ]]; then
    echo "The -n flag is not used in Jasmine upload mode."
    exit 1
fi

if [[ "${diagnose_flag}" = "true" && "${submit_flag}" = "true" ]]; then
    echo "The -s flag is not used in diagnose mode."
    exit 1
fi

if [[ "${diagnose_flag}" = "true" && "${no_run_flag}" = "true" ]]; then
    echo "The -n flag is not used in diagnose mode."
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
    JASMINE_PRODUCT_KIND="skim"
else
    MODE="replay"
    MODE_SCRIPT="${REPLAY_SUBMIT_SCRIPT}"
    WORKFLOW_SUFFIX="_${USER_NAME}"
    JASMINE_PRODUCT_KIND="replay"
fi

if [[ -n "${workflow_override}" ]]; then
    WORKFLOW="${workflow_override}"
else
    WORKFLOW="${WORKFLOW_PREFIX}_${FAMILY}${WORKFLOW_SUFFIX}"
fi

if [[ -z "${manifest_dir}" ]]; then
    manifest_dir="${DEFAULT_MANIFEST_DIR}"
elif [[ "${manifest_dir}" != /* ]]; then
    manifest_dir="${SCRIPT_DIR}/${manifest_dir}"
fi

if [[ "${diagnose_flag}" = "true" ]]; then
    echo
    echo "---------------------------------------------------------"
    echo
    echo "Diagnosing failed ${MODE} jobs for Q2=${Q2}, W=${W}..."
    echo
    echo "                        WORKFLOW = ${WORKFLOW}"
    echo
    echo "---------------------------------------------------------"
    echo

    diagnose_cmd=(python3 "${DIAGNOSE_SCRIPT}" "${WORKFLOW}")
    echo "Running: ${diagnose_cmd[*]}"
    "${diagnose_cmd[@]}"
    exit $?
fi

if [[ "${upload_flag}" = "true" ]]; then
    shopt -s nullglob
    matched_manifests=( "${manifest_dir}/${FAMILY}"*.json )
    shopt -u nullglob

    if [[ ${#matched_manifests[@]} -eq 0 ]]; then
        echo "No manifests found for ${FAMILY} in ${manifest_dir}"
        exit 1
    fi

    echo
    echo "---------------------------------------------------------"
    echo
    echo "Preparing interactive Jasmine ${JASMINE_PRODUCT_KIND} upload for Q2=${Q2}, W=${W}..."
    echo
    echo "                    MANIFEST DIR = ${manifest_dir}"
    echo "                   PRODUCT KIND = ${JASMINE_PRODUCT_KIND}"
    echo
    echo "---------------------------------------------------------"
    echo

    for manifest_path in "${matched_manifests[@]}"; do
        upload_cmd=(python3 "${JASMINE_SCRIPT}" --manifest-path "${manifest_path}" --product-kind "${JASMINE_PRODUCT_KIND}")
        if [[ "${submit_flag}" = "true" ]]; then
            upload_cmd+=(--submit)
        fi

        echo "Running: ${upload_cmd[*]}"
        "${upload_cmd[@]}" || exit $?
        echo
    done
    exit 0
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
if [[ -n "${account_override}" ]]; then
    submit_cmd+=(--account "${account_override}")
fi
if [[ -n "${partition_override}" ]]; then
    submit_cmd+=(--partition "${partition_override}")
fi

echo "Running: ${submit_cmd[*]}"
"${submit_cmd[@]}"
