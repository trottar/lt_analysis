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
python -m py_compile src/utility/regenerate_workflow_diagnostic_plots.py

python - <<'PY'
import contextlib
import io
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
import regenerate_workflow_diagnostic_plots

print("analysis hardening pure-python imports OK")
PY

python - <<'PY'
import math
import os
import sys

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

from background_config import (
    get_particle_subtraction_window_config,
    resolve_particle_subtraction_windows,
)

config = get_particle_subtraction_window_config("kaon", "pion")
if config is None:
    raise AssertionError("Expected kaon->pion subtraction-window config to exist")
if tuple(config["windows"]["pi_n"]) != (0.88, 0.94):
    raise AssertionError("Unexpected pi_n subtraction window: {}".format(config["windows"]["pi_n"]))
if tuple(config["windows"]["pi_delta"]) != (1.35, 1.45):
    raise AssertionError("Unexpected pi_delta subtraction window: {}".format(config["windows"]["pi_delta"]))

resolved_nominal = resolve_particle_subtraction_windows("kaon", "pion", 0.0)
resolved_shifted = resolve_particle_subtraction_windows("kaon", "pion", 0.02)
for key, expected in {"pi_n": (0.88, 0.94), "pi_delta": (1.35, 1.45)}.items():
    actual = resolved_nominal.get(key)
    if actual is None or any(not math.isclose(actual[idx], expected[idx], rel_tol=0.0, abs_tol=1e-12) for idx in range(2)):
        raise AssertionError("Unexpected nominal resolved window for {}: {}".format(key, actual))
    shifted = resolved_shifted.get(key)
    if shifted is None or any(
        not math.isclose(shifted[idx], expected[idx] + 0.02, rel_tol=0.0, abs_tol=1e-12)
        for idx in range(2)
    ):
        raise AssertionError("Unexpected shifted resolved window for {}: {}".format(key, shifted))

print("analysis hardening subtraction-window resolver smoke OK")
PY

python - <<'PY'
import contextlib
import io
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
            "particle_subtraction_windows": {
                "apply_mm_offset_data": True,
                "windows": {"pi_n": [0.88, 0.94], "pi_delta": [1.35, 1.45]},
            },
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

python - "$@" <<'PY'
import os
import sys

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

import frozen_manifest


def _is_manifest_filename(name):
    return name.endswith(".json") and "_frozen_manifest" in name


def _discover_default_manifest_paths(repo_root):
    analysis_root = os.path.join(repo_root, "OUTPUT", "Analysis")
    manifests = []
    if not os.path.isdir(analysis_root):
        return manifests
    for analysis_name in sorted(os.listdir(analysis_root)):
        analysis_dir = os.path.join(analysis_root, analysis_name)
        if not os.path.isdir(analysis_dir):
            continue
        for entry in sorted(os.listdir(analysis_dir)):
            if _is_manifest_filename(entry):
                manifests.append(os.path.join(analysis_dir, entry))
    return manifests


def _discover_manifest_paths(raw_paths):
    manifests = []
    for raw_path in raw_paths:
        abs_path = os.path.abspath(raw_path)
        if os.path.isfile(abs_path):
            if not _is_manifest_filename(os.path.basename(abs_path)):
                raise FileNotFoundError(
                    "Path is not a frozen manifest JSON: {}".format(abs_path)
                )
            manifests.append(abs_path)
            continue
        if not os.path.isdir(abs_path):
            raise FileNotFoundError("Manifest/cache path does not exist: {}".format(abs_path))

        discovered = []
        for candidate_root in (abs_path, os.path.join(abs_path, "json")):
            if not os.path.isdir(candidate_root):
                continue
            for entry in sorted(os.listdir(candidate_root)):
                if _is_manifest_filename(entry):
                    discovered.append(os.path.join(candidate_root, entry))
        if not discovered:
            for walk_root, _, files in os.walk(abs_path):
                for entry in sorted(files):
                    if _is_manifest_filename(entry):
                        discovered.append(os.path.join(walk_root, entry))
        if not discovered:
            raise FileNotFoundError(
                "No frozen manifest JSON files found under {}".format(abs_path)
            )
        manifests.extend(discovered)

    unique_manifests = []
    seen = set()
    for manifest_path in manifests:
        normalized = os.path.abspath(manifest_path)
        if normalized in seen:
            continue
        seen.add(normalized)
        unique_manifests.append(normalized)
    return unique_manifests


repo_root = os.getcwd()
requested_paths = sys.argv[1:]
if requested_paths:
    manifest_paths = _discover_manifest_paths(requested_paths)
else:
    manifest_paths = _discover_default_manifest_paths(repo_root)

if not manifest_paths:
    print(
        "analysis hardening cache-manifest preflight skipped "
        "(no staged top-level manifests found; pass a cache dir or manifest path to check a frozen run explicitly)"
    )
    raise SystemExit(0)

strict_hash_paths, warn_only_hash_paths = frozen_manifest.get_iteration_manifest_hash_policy()
validated = 0
warning_count = 0
for manifest_path in manifest_paths:
    manifest = frozen_manifest.load_frozen_manifest(manifest_path)
    try:
        result = frozen_manifest.validate_manifest_hashes_against_repo(
            manifest,
            repo_root,
            strict_hash_paths=strict_hash_paths,
            warn_only_hash_paths=warn_only_hash_paths,
            allow_config_drift=frozen_manifest.ALLOW_CONFIG_DRIFT,
            emit_warnings=False,
            manifest_path=manifest_path,
        )
    except (FileNotFoundError, ValueError) as exc:
        print("analysis hardening cache-manifest preflight FAILED")
        print("Manifest: {}".format(manifest_path))
        print("")
        print(exc)
        raise SystemExit(2)
    warnings = result.get("warnings", [])
    warning_count += len(warnings)
    validated += 1
    if warnings:
        print(
            "analysis hardening cache-manifest preflight WARN: {} ({} warning-only drift item(s))".format(
                manifest_path,
                len(warnings),
            )
        )

print(
    "analysis hardening cache-manifest preflight OK ({} manifest(s), {} warning-only drift item(s))".format(
        validated,
        warning_count,
    )
)
PY

python - <<'PY'
import contextlib
import io
import os
import sys
import tempfile
import json

utility_path = os.path.join(os.getcwd(), "src", "utility")
if utility_path not in sys.path:
    sys.path.append(utility_path)

import frozen_manifest
import regenerate_workflow_diagnostic_plots

strict_hash_paths, warn_only_hash_paths = frozen_manifest.get_iteration_manifest_hash_policy("src/main_iter.py")

manifest = {
    "q2": "test",
    "w": "test",
    "epsilon_values": {"low": "0.1"},
    "t_bin_edges": [0.0, 1.0],
    "phi_bin_edges": [0.0, 90.0],
    "key_file_hashes": {
        "src/utility/background_config.py": {"sha256": "expected-background-config"},
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
        "src/utility/background_config.py": {"sha256": "different-background-config"},
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
        if "src/utility/background_config.py" not in str(exc):
            raise
    else:
        raise AssertionError("Expected background-config hash mismatch to fail validation")

    frozen_manifest.build_key_file_hashes = lambda repo_root: {
        "src/utility/background_config.py": {"sha256": "expected-background-config"},
        "src/main_auto.py": {"sha256": "different-main-auto"},
        "src/main_iter.py": {"sha256": "different-main-iter"},
    }
    validated = dict(inp_dict)
    with contextlib.redirect_stdout(io.StringIO()):
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
    if not any("src/main_iter.py" in warning for warning in warnings):
        raise AssertionError("Expected current-driver code drift to warn, not fail")
finally:
    frozen_manifest.ALLOW_CONFIG_DRIFT = original_allow
    frozen_manifest.build_key_file_hashes = original_builder

print("analysis hardening manifest hash-policy smoke OK")

with tempfile.TemporaryDirectory() as tmpdir:
    plot_dir = os.path.join(tmpdir, "plots")
    json_dir = os.path.join(tmpdir, "json")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(json_dir, exist_ok=True)

    particle = "kaon"
    q2 = "test"
    w = "test"
    active_profile = "fit1_only"
    artifact_paths = frozen_manifest.get_analysis_artifact_paths(tmpdir, particle, q2, w, active_profile=active_profile)
    ledger_low_paths = frozen_manifest.get_correction_ledger_paths(tmpdir, particle, "test_low", active_profile=active_profile)
    ledger_high_paths = frozen_manifest.get_correction_ledger_paths(tmpdir, particle, "test_high", active_profile=active_profile)

    manifest_payload = {
        "particle_type": particle,
        "q2": q2,
        "w": w,
        "active_profile": active_profile,
        "epsilon_values": {"low": 0.1, "high": 0.2},
        "mm_cut_window": [1.10, 1.15],
        "particle_subtraction_windows": {
            "apply_mm_offset_data": True,
            "windows": {"pi_n": [0.88, 0.94], "pi_delta": [1.35, 1.45]},
        },
        "t_bin_edges": [0.1, 0.2, 0.3],
        "phi_bin_edges": [-180.0, -60.0, 60.0, 180.0],
        "selected_bg_scale1s": {"low": {"low:Center": 0.4}, "high": {"high:Center": 0.5}},
        "selected_bg_scale2s": {"low": {"low:Center": 0.2}, "high": {"high:Center": 0.3}},
        "resolved_optimizer_settings": {"selection_mode": "weighted"},
        "common_epsilon_scale_behavior": "independent",
        "created_at": "2026-05-20T00:00:00Z",
    }
    bundle_payload = {
        "particle_type": particle,
        "q2": q2,
        "w": w,
        "active_profile": active_profile,
        "low": {"argv": ["src/main.py"], "inp_dict": {"OutFilename": "test_low", "EPSSET": "low", "EPSVAL": 0.1, "NumtBins": 2, "NumPhiBins": 3, "ParticleType": particle, "POL": 1}},
        "high": {"argv": ["src/main.py"], "inp_dict": {"OutFilename": "test_high", "EPSSET": "high", "EPSVAL": 0.2, "NumtBins": 2, "NumPhiBins": 3, "ParticleType": particle, "POL": 1}},
    }

    def build_ledger(epsset, fit1_total, fit2_total):
        stage_yields = {
            "raw_prompt": 100.0,
            "after_random_subtraction": 90.0,
            "after_dummy_subtraction": 82.0,
            "after_pion_subtraction": 76.0,
            "after_empirical_fit1": 76.0 - fit1_total,
            "after_empirical_fit2": 76.0 - fit1_total - fit2_total,
            "final_lambda_window": 76.0 - fit1_total - fit2_total,
        }
        return {
            "particle_type": particle,
            "epsset": epsset,
            "q2": q2,
            "w": w,
            "outfilename": "test_{}".format(epsset),
            "particle_subtraction_windows": {
                "apply_mm_offset_data": True,
                "windows": {"pi_n": [0.88, 0.94], "pi_delta": [1.35, 1.45]},
            },
            "active_profile": active_profile,
            "resolved_optimizer_settings": {"selection_mode": "weighted"},
            "settings": [
                {
                    "phi_setting": "Center",
                    "bin_level_ledger_available": True,
                    "combined_stage_yields": dict(stage_yields),
                    "fit1_fractional_correction": fit1_total / 76.0,
                    "fit2_fractional_correction": fit2_total / (76.0 - fit1_total),
                    "total_empirical_fraction": (fit1_total + fit2_total) / 76.0,
                    "empirical_fit_uncertainty_frac": 0.02,
                    "failed_fit_count": 0,
                    "zero_background_fallback_count": 0,
                    "oversub_diagnostics": {
                        "fit1": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.0},
                        "fit2": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.0},
                    },
                    "t_phi_bins": [
                        {
                            "bin_key": "t_bin1phi_bin1",
                            "stage_yields": dict(stage_yields),
                            "fit1_fractional_correction": fit1_total / 76.0,
                            "fit2_fractional_correction": fit2_total / (76.0 - fit1_total),
                            "total_empirical_fraction": (fit1_total + fit2_total) / 76.0,
                            "empirical_fit_uncertainty_frac": 0.02,
                            "fit1_failed": False,
                            "fit2_failed": False,
                            "fit1_zero_background_fallback": False,
                            "fit2_zero_background_fallback": False,
                            "oversub_diagnostics": {
                                "fit1": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.0},
                                "fit2": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.0},
                            },
                        },
                        {
                            "bin_key": "t_bin2phi_bin3",
                            "stage_yields": {
                                "raw_prompt": 50.0,
                                "after_random_subtraction": 45.0,
                                "after_dummy_subtraction": 41.0,
                                "after_pion_subtraction": 38.0,
                                "after_empirical_fit1": 36.5,
                                "after_empirical_fit2": 35.5,
                                "final_lambda_window": 35.5,
                            },
                            "fit1_fractional_correction": 1.5 / 38.0,
                            "fit2_fractional_correction": 1.0 / 36.5,
                            "total_empirical_fraction": 2.5 / 38.0,
                            "empirical_fit_uncertainty_frac": 0.03,
                            "fit1_failed": False,
                            "fit2_failed": False,
                            "fit1_zero_background_fallback": False,
                            "fit2_zero_background_fallback": False,
                            "oversub_diagnostics": {
                                "fit1": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.01},
                                "fit2": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.02},
                            },
                        },
                    ],
                }
            ],
            "combined_totals": {
                "combined_stage_yields": dict(stage_yields),
                "all_settings_have_bin_level_ledger": True,
                "fit1_fractional_correction": fit1_total / 76.0,
                "fit2_fractional_correction": fit2_total / (76.0 - fit1_total),
                "total_empirical_fraction": (fit1_total + fit2_total) / 76.0,
                "empirical_fit_uncertainty_frac": 0.02,
                "failed_fit_count": 0,
                "zero_background_fallback_count": 0,
                "oversub_diagnostics": {
                    "fit1": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.01},
                    "fit2": {"oversub_bin_count": 0, "max_unclamped_ratio": 1.02},
                },
            },
        }

    low_ledger = build_ledger("low", 4.0, 2.0)
    high_ledger = build_ledger("high", 5.0, 1.0)
    epsilon_compare = {
        "particle_type": particle,
        "q2": q2,
        "w": w,
        "active_profile": active_profile,
        "rows": [
            {"phi_setting": "Center", "bin_key": "t_bin1phi_bin1", "f_emp_low": 0.08, "f_emp_high": 0.09, "delta_f_emp": 0.01},
            {"phi_setting": "Center", "bin_key": "t_bin2phi_bin3", "f_emp_low": 0.07, "f_emp_high": 0.05, "delta_f_emp": -0.02},
        ],
        "mean_f_emp_low": 0.075,
        "mean_f_emp_high": 0.07,
        "rms_delta_f_emp": 0.0158,
        "max_abs_delta_f_emp": 0.02,
        "warning_threshold": 0.03,
        "fail_threshold": None,
    }
    final_summary = {
        "particle_type": particle,
        "q2": q2,
        "w": w,
        "active_profile": active_profile,
        "common_epsilon_scale_behavior": "independent",
        "particle_subtraction_windows": manifest_payload["particle_subtraction_windows"],
        "resolved_optimizer_settings": {"selection_mode": "weighted"},
        "fit1_scales": manifest_payload["selected_bg_scale1s"],
        "fit2_scales": manifest_payload["selected_bg_scale2s"],
        "low_epsilon_corrections": low_ledger["combined_totals"],
        "high_epsilon_corrections": high_ledger["combined_totals"],
        "epsilon_empirical_compare": epsilon_compare,
        "xsect_outputs": {
            "separated_csv": [{"sigma_L": "0.11", "sigma_T": "0.22", "sigma_LT": "0.03", "sigma_TT": "0.04"}],
            "unseparated_csv": [{"xsec": "0.33"}],
        },
        "manifest_hash": "synthetic",
        "iteration_count": 2,
        "convergence_status": "synthetic",
    }
    systematics_summary = {
        "nominal_profile": "nominal_weighted",
        "generated_at": "2026-05-20T00:00:00Z",
        "rows": [
            {"profile_name": "nominal_weighted", "lambda_yield": 100.0, "percent_deviation_from_nominal": 0.0, "pass_fail_status": "passed", "sigma_L": 0.1, "sigma_T": 0.2, "sigma_LT": 0.03, "sigma_TT": 0.04, "error_message": None},
            {"profile_name": "fit1_only", "lambda_yield": 97.0, "percent_deviation_from_nominal": -3.0, "pass_fail_status": "passed", "sigma_L": 0.09, "sigma_T": 0.19, "sigma_LT": 0.025, "sigma_TT": 0.035, "error_message": None},
        ],
        "lambda_yield_rms": 2.121,
        "lambda_yield_envelope": 3.0,
        "preflight": {"passed": True},
    }
    nonklambda_summary = {
        "comparison_mode": "preprocessed_prediction",
        "active_profile": active_profile,
        "prediction_count": 2,
        "missing_empirical_bin_count": 0,
        "max_abs_lambda_window_correction_difference": 0.02,
        "warnings": [],
        "rows": [
            {"channel_name": "nonKL", "status": "ok", "lambda_window_correction_difference": 0.02},
            {"channel_name": "nonKL", "status": "ok", "lambda_window_correction_difference": -0.01},
        ],
    }

    for path, payload in [
        (os.path.join(json_dir, os.path.basename(artifact_paths["manifest_profile"])), manifest_payload),
        (os.path.join(json_dir, os.path.basename(artifact_paths["input_bundle_profile"])), bundle_payload),
        (os.path.join(json_dir, os.path.basename(artifact_paths["epsilon_compare_json_profile"])), epsilon_compare),
        (os.path.join(json_dir, os.path.basename(artifact_paths["final_summary_json_profile"])), final_summary),
        (os.path.join(json_dir, os.path.basename(artifact_paths["systematics_json"])), systematics_summary),
        (os.path.join(json_dir, os.path.basename(artifact_paths["nonklambda_json_profile"])), nonklambda_summary),
        (os.path.join(json_dir, os.path.basename(ledger_low_paths["json_profile"])), low_ledger),
        (os.path.join(json_dir, os.path.basename(ledger_high_paths["json_profile"])), high_ledger),
    ]:
        with open(path, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)

    status = regenerate_workflow_diagnostic_plots.main([tmpdir, "--profile", active_profile])
    if status != 0:
        raise AssertionError("Expected workflow diagnostic plot regeneration smoke to succeed")

    expected_pdf = os.path.join(plot_dir, "{}_Q{}W{}_workflow_diagnostics_{}.pdf".format(particle, q2, w, active_profile))
    if not os.path.exists(expected_pdf):
        raise AssertionError("Expected workflow diagnostic PDF was not created: {}".format(expected_pdf))

print("analysis hardening workflow-plot smoke OK")
PY

python - <<'PY'
from pathlib import Path

for rel_path in ("src/main_iter.py", "src/main_auto.py"):
    text = Path(rel_path).read_text()
    validate_index = text.find("validate_iteration_inputs_against_manifest(")
    xfit_index = text.find("x_fit_in_t(")
    if validate_index < 0 or xfit_index < 0:
        raise AssertionError("Could not find validation/x_fit markers in {}".format(rel_path))
    if validate_index > xfit_index:
        raise AssertionError(
            "Frozen-manifest validation occurs after x_fit_in_t in {}; fail-fast ordering regressed".format(rel_path)
        )

print("analysis hardening iteration-preflight ordering OK")
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
    "particle_subtraction_windows": {
        "apply_mm_offset_data": True,
        "windows": {"pi_n": [0.88, 0.94], "pi_delta": [1.35, 1.45]},
    },
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
