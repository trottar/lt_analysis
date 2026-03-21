#! /bin/bash

set -euo pipefail

if [[ $# -ne 10 ]]; then
    echo "Usage: ./run_shift.sh <ltanapath> <particle_type> <simc_root> <data_root> <Q2> <W> <theta_cm_deg> <beam_energy_gev> <mm_min> <mm_max>" >&2
    exit 1
fi

LTANAPATH="$1"
PARTICLE_TYPE="$2"
SIMC_ROOT="$3"
DATA_ROOT="$4"
Q2="$5"
W="$6"
THETA_CM_DEG="$7"
BEAM_ENERGY_GEV="$8"
MM_MIN="$9"
MM_MAX="${10}"

if command -v cygpath >/dev/null 2>&1; then
    LTANAPATH="$(cygpath -u "${LTANAPATH}" 2>/dev/null || printf '%s' "${LTANAPATH}")"
fi

SRC_DIR="${LTANAPATH}/src"
TSHIFT_EXE="none"

if [[ "${PARTICLE_TYPE}" == "kaon" ]]; then
    TSHIFT_SRC="${SRC_DIR}/tshift_kaonff.f"
    TSHIFT_INC="${SRC_DIR}/linux_suppl.inc"
    TSHIFT_EXE="${SRC_DIR}/tshift_kaonff"

    echo "Compiling tshift_kaonff.f..." >&2
    (
        cd "${SRC_DIR}"
        gfortran -o tshift_kaonff tshift_kaonff.f
    ) >&2
fi

python3 "${SRC_DIR}/setup/calc_shift_values.py" \
    "${PARTICLE_TYPE}" \
    "${SIMC_ROOT}" \
    "${DATA_ROOT}" \
    "${Q2}" \
    "${W}" \
    "${THETA_CM_DEG}" \
    "${BEAM_ENERGY_GEV}" \
    "${MM_MIN}" \
    "${MM_MAX}" \
    "${TSHIFT_EXE}"
