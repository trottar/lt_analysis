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

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

import frozen_manifest

strict_hash_paths, warn_only_hash_paths = frozen_manifest.get_iteration_manifest_hash_policy("src/main_iter.py")

manifest = {
    "q2": "test",
    "w": "test",
    "epsilon_values": {"low": "0.1"},
    "t_bin_edges": [0.0, 1.0],
    "phi_bin_edges": [0.0, 90.0],
    "key_file_hashes": {
        "src/main_auto.py": {"sha256": "expected-main-auto"},
        "src/main_iter.py": {"sha256": "expected-main-iter"},
    },
}
inp_dict = {
    "Q2": "test",
    "W": "test",
    "EPSSET": "low",
    "EPSVAL": "0.1",
    "NumtBins": 1,
    "NumPhiBins": 1,
}
original_allow = frozen_manifest.ALLOW_CONFIG_DRIFT
original_builder = frozen_manifest.build_key_file_hashes
try:
    frozen_manifest.ALLOW_CONFIG_DRIFT = False
    frozen_manifest.build_key_file_hashes = lambda repo_root: {
        "src/main_auto.py": {"sha256": "different-main-auto"},
        "src/main_iter.py": {"sha256": "different-main-iter"},
    }

    try:
        frozen_manifest.validate_iteration_inputs_against_manifest(
            manifest,
            dict(inp_dict),
            [0.0, 1.0],
            [0.0, 90.0],
            ".",
            strict_hash_paths=strict_hash_paths,
            warn_only_hash_paths=warn_only_hash_paths,
        )
    except ValueError as exc:
        if "src/main_iter.py" not in str(exc):
            raise
    else:
        raise AssertionError("Expected current-driver hash mismatch to fail validation")

    frozen_manifest.build_key_file_hashes = lambda repo_root: {
        "src/main_auto.py": {"sha256": "different-main-auto"},
        "src/main_iter.py": {"sha256": "expected-main-iter"},
    }
    validated = dict(inp_dict)
    frozen_manifest.validate_iteration_inputs_against_manifest(
        manifest,
        validated,
        [0.0, 1.0],
        [0.0, 90.0],
        ".",
        strict_hash_paths=strict_hash_paths,
        warn_only_hash_paths=warn_only_hash_paths,
    )
    warnings = validated.get("manifest_validation_warnings", [])
    if not any("src/main_auto.py" in warning for warning in warnings):
        raise AssertionError("Expected sibling-driver drift to warn, not fail")
finally:
    frozen_manifest.ALLOW_CONFIG_DRIFT = original_allow
    frozen_manifest.build_key_file_hashes = original_builder

print("analysis hardening manifest hash-policy smoke OK")
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
