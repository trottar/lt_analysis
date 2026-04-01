#!/bin/bash

set -euo pipefail

if [[ "$#" -lt 1 ]]; then
    echo "ERROR: run_with_softenv.sh needs a command to execute" >&2
    exit 2
fi

if [[ -f /site/12gev_phys/softenv.sh ]]; then
    source /site/12gev_phys/softenv.sh 2.3 >/dev/null 2>&1
fi

exec "$@"
