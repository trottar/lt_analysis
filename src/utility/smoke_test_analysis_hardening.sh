#!/usr/bin/env bash

set -euo pipefail
export MPLBACKEND="${MPLBACKEND:-Agg}"

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

python - <<'PY'
import os
import sys
import tempfile

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

from epsilon_correction_compare import write_epsilon_empirical_compare
from nonklambda_crosscheck import write_nonklambda_crosscheck
from write_final_analysis_summary import write_final_analysis_summary

with tempfile.TemporaryDirectory() as tmpdir:
    write_epsilon_empirical_compare(
        {
            "particle_type": "kaon",
            "q2": "test",
            "w": "test",
            "active_profile": "fit1_only",
            "rows": [],
            "summary": {},
            "warnings": [],
        },
        tmpdir,
        "kaon",
        "test",
        "test",
        active_profile="fit1_only",
    )
    write_final_analysis_summary(
        {
            "particle_type": "kaon",
            "q2": "test",
            "w": "test",
            "active_profile": "fit1_only",
            "common_epsilon_scale_behavior": "independent",
            "low_epsilon_corrections": {},
            "high_epsilon_corrections": {},
            "epsilon_empirical_compare": {},
            "xsect_outputs": {},
        },
        tmpdir,
        "kaon",
        "test",
        "test",
        active_profile="fit1_only",
    )
    write_nonklambda_crosscheck(
        {
            "comparison_mode": "preprocessed_prediction",
            "active_profile": "fit1_only",
            "rows": [],
            "warnings": [],
        },
        tmpdir,
        "kaon",
        "test",
        "test",
        active_profile="fit1_only",
    )

print("analysis hardening profile-aware writer smoke OK")
PY

python - <<'PY'
import os
import sys
import tempfile

try:
    import ROOT  # noqa: F401
except Exception:
    print("ROOT unavailable; skipping correction-ledger runtime smoke")
    raise SystemExit(0)

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

from correction_ledger import write_correction_ledger

payload = {
    "particle_type": "kaon",
    "epsset": "low",
    "q2": "test",
    "w": "test",
    "outfilename": "test",
    "mm_cut_window": [1.0, 1.2],
    "active_profile": "fit1_only",
    "resolved_optimizer_settings": {},
    "settings": [],
    "combined_totals": {
        "combined_stage_yields": {},
        "all_settings_have_bin_level_ledger": False,
    },
}

with tempfile.TemporaryDirectory() as tmpdir:
    write_correction_ledger(payload, tmpdir, "kaon", "test", active_profile="fit1_only")

print("analysis hardening correction-ledger runtime smoke OK")
PY
