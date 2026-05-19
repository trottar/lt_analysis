#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"

cd "${repo_root}"

bash -n run_Prod_Analysis.sh

python -m py_compile src/utility/background_config.py
python -m py_compile src/utility/frozen_manifest.py
python -m py_compile src/utility/correction_ledger.py
python -m py_compile src/utility/epsilon_correction_compare.py
python -m py_compile src/utility/nonklambda_crosscheck.py
python -m py_compile src/utility/run_bg_systematics.py
python -m py_compile src/utility/write_final_analysis_summary.py
python -m py_compile src/utility/mm_background_subtraction.py
python -m py_compile src/utility/check_nominal_regression.py

python - <<'PY'
import os
import sys

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

import background_config
import frozen_manifest
import epsilon_correction_compare
import nonklambda_crosscheck
import run_bg_systematics
import write_final_analysis_summary
import mm_background_subtraction
import check_nominal_regression

print("analysis hardening pure-python imports OK")
PY
