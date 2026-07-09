#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <phi_setting> <Q2> <W> <eps_setting>" >&2
    echo "Example: $0 Center 3p0 3p14 low" >&2
    exit 1
fi

phi_setting="$1"
Q2="$2"
W="$3"
eps_setting="$4"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
macro="${script_dir}/check_slow_protons.C+"

if [[ ! -f "$macro" ]]; then
    echo "Error: ROOT macro not found: $macro" >&2
    exit 2
fi

echo "Running check_slow_protons.C"
root -l -b -q \
    ".L ${macro}(\"${phi_setting}\",\"${Q2}\",\"${W}\",\"${eps_setting}\")"
