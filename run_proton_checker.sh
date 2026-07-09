#!/usr/bin/env bash

set -euo pipefail

runner_version="2"

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <phi_setting> <Q2> <W> <eps_setting>" >&2
    echo "Example: $0 Center 3p0 3p14 high" >&2
    exit 1
fi

phi_setting="$1"
Q2="$2"
W="$3"
eps_setting="$4"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
macro="${script_dir}/check_slow_protons.C"

if [[ ! -f "$macro" ]]; then
    echo "Error:  ROOT macro not found: $macro" >&2
    exit 2
fi

echo "Running check_slow_protons.C (runner version ${runner_version}; RF timing first, CTime_ROC1 fallback)"
root -l -b -q \
    "${macro}(\"${phi_setting}\",\"${Q2}\",\"${W}\",\"${eps_setting}\")"
