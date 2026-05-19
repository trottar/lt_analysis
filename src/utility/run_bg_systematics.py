#! /usr/bin/python

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from datetime import datetime

from background_config import BG_SYSTEMATIC_PROFILES
from frozen_manifest import (
    get_analysis_artifact_paths,
    get_correction_ledger_paths,
    load_zeroth_iteration_input_bundle,
)


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def _write_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _extract_scalar(rows, candidate_keys):
    if not rows:
        return None
    first_row = rows[0]
    for key in candidate_keys:
        if key in first_row and first_row[key] not in ("", None):
            try:
                return float(first_row[key])
            except Exception:
                return first_row[key]
    return None


def _copy_if_exists(src, dst_dir):
    if not src or not os.path.exists(src):
        return
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy(src, os.path.join(dst_dir, os.path.basename(src)))


def _run_profile(argv, profile_name, repo_root, snapshot_path):
    env = os.environ.copy()
    env["LT_BG_PROFILE_SNAPSHOT_JSON"] = snapshot_path
    subprocess.run(
        [sys.executable] + list(argv),
        cwd=repo_root,
        check=True,
        env=env,
    )


def _profile_output_dir(outpath, particle_type, q2, w, profile_name):
    return os.path.join(
        outpath,
        "bg_systematics",
        "{}_Q{}W{}".format(particle_type, q2, w),
        profile_name,
    )


def run_bg_systematics(bundle, repo_root, outpath, particle_type, q2, w):
    results = []
    nominal_row = None
    for profile_name in BG_SYSTEMATIC_PROFILES:
        profile_dir = _profile_output_dir(outpath, particle_type, q2, w, profile_name)
        os.makedirs(profile_dir, exist_ok=True)
        snapshot_path = os.path.join(profile_dir, "bg_profile_snapshot.json")
        _write_json(snapshot_path, {"BG_OPT_ACTIVE_PROFILE": profile_name})

        row = {
            "profile_name": profile_name,
            "selected_t_bins": None,
            "selected_phi_bins": None,
            "selected_fit1_scales": None,
            "selected_fit2_scales": None,
            "lambda_yield": None,
            "unseparated_cross_section": None,
            "sigma_L": None,
            "sigma_T": None,
            "sigma_LT": None,
            "sigma_TT": None,
            "percent_deviation_from_nominal": None,
            "pass_fail_status": "failed",
            "error_message": None,
        }
        try:
            low_entry = bundle.get("low", {})
            high_entry = bundle.get("high", {})
            if not low_entry or not high_entry:
                raise ValueError("0th iteration input bundle must contain both low and high entries")

            _run_profile(low_entry["argv"], profile_name, repo_root, snapshot_path)
            _run_profile(high_entry["argv"], profile_name, repo_root, snapshot_path)

            artifact_paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=profile_name)
            final_summary = _load_json(artifact_paths["final_summary_json_profile"])
            row["selected_t_bins"] = final_summary.get("t_bin_edges")
            row["selected_phi_bins"] = final_summary.get("phi_bin_edges")
            row["selected_fit1_scales"] = final_summary.get("fit1_scales")
            row["selected_fit2_scales"] = final_summary.get("fit2_scales")

            low_ledger_path = get_correction_ledger_paths(
                outpath,
                particle_type,
                bundle["low"]["inp_dict"]["OutFilename"],
                active_profile=profile_name,
            )["json_profile"]
            high_ledger_path = get_correction_ledger_paths(
                outpath,
                particle_type,
                bundle["high"]["inp_dict"]["OutFilename"],
                active_profile=profile_name,
            )["json_profile"]
            low_ledger = _load_json(low_ledger_path)
            high_ledger = _load_json(high_ledger_path)
            row["lambda_yield"] = float(low_ledger["combined_totals"]["combined_stage_yields"]["final_lambda_window"]) + float(
                high_ledger["combined_totals"]["combined_stage_yields"]["final_lambda_window"]
            )

            xsect_outputs = final_summary.get("xsect_outputs", {})
            separated_rows = xsect_outputs.get("separated_csv") or []
            unsep_rows = xsect_outputs.get("unseparated_csv") or []
            row["unseparated_cross_section"] = _extract_scalar(unsep_rows, ["xsec", "unsep_xsec", "cross_section"])
            row["sigma_L"] = _extract_scalar(separated_rows, ["sigma_L", "sigL", "L"])
            row["sigma_T"] = _extract_scalar(separated_rows, ["sigma_T", "sigT", "T"])
            row["sigma_LT"] = _extract_scalar(separated_rows, ["sigma_LT", "sigLT", "LT"])
            row["sigma_TT"] = _extract_scalar(separated_rows, ["sigma_TT", "sigTT", "TT"])
            row["pass_fail_status"] = "passed"

            for artifact in [
                artifact_paths["manifest"],
                artifact_paths["manifest_profile"],
                artifact_paths["input_bundle"],
                artifact_paths["input_bundle_profile"],
                artifact_paths["epsilon_compare_json"],
                artifact_paths["epsilon_compare_json_profile"],
                artifact_paths["epsilon_compare_csv"],
                artifact_paths["epsilon_compare_csv_profile"],
                artifact_paths["final_summary_json"],
                artifact_paths["final_summary_json_profile"],
                artifact_paths["final_summary_csv"],
                artifact_paths["final_summary_csv_profile"],
                artifact_paths["final_summary_md"],
                artifact_paths["final_summary_md_profile"],
                low_ledger_path,
                high_ledger_path,
            ]:
                _copy_if_exists(artifact, profile_dir)
        except Exception as exc:
            row["error_message"] = str(exc)
        results.append(row)
        if profile_name == "nominal_weighted" and row["pass_fail_status"] == "passed":
            nominal_row = deepcopy(row)

    nominal_value = None if nominal_row is None else nominal_row.get("lambda_yield")
    passed_values = []
    for row in results:
        if row["pass_fail_status"] != "passed":
            continue
        if nominal_value not in (None, 0):
            row["percent_deviation_from_nominal"] = 100.0 * (float(row["lambda_yield"]) - float(nominal_value)) / float(nominal_value)
        passed_values.append(float(row["lambda_yield"]))

    summary = {
        "particle_type": particle_type,
        "q2": q2,
        "w": w,
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "nominal_profile": "nominal_weighted",
        "rows": results,
        "lambda_yield_rms": None,
        "lambda_yield_envelope": None,
    }
    if nominal_value not in (None, 0) and passed_values:
        deltas = [value - float(nominal_value) for value in passed_values]
        summary["lambda_yield_rms"] = (sum(delta * delta for delta in deltas) / len(deltas)) ** 0.5
        summary["lambda_yield_envelope"] = max((abs(delta) for delta in deltas), default=0.0)
    return summary


def write_bg_systematics_summary(summary, outpath, particle_type, q2, w, active_profile=None):
    paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    _write_json(paths["systematics_json"], summary)
    fieldnames = [
        "profile_name",
        "selected_t_bins",
        "selected_phi_bins",
        "selected_fit1_scales",
        "selected_fit2_scales",
        "lambda_yield",
        "unseparated_cross_section",
        "sigma_L",
        "sigma_T",
        "sigma_LT",
        "sigma_TT",
        "percent_deviation_from_nominal",
        "pass_fail_status",
        "error_message",
    ]
    with open(paths["systematics_csv"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary.get("rows", []):
            writer.writerow({key: row.get(key) for key in fieldnames})
    return [paths["systematics_json"], paths["systematics_csv"]]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Replay configured background profiles from a 0th-iteration input bundle.")
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument("--outpath", required=True)
    parser.add_argument("--particle", required=True)
    parser.add_argument("--q2", required=True)
    parser.add_argument("--w", required=True)
    args = parser.parse_args()

    bundle = load_zeroth_iteration_input_bundle(args.outpath, args.particle, args.q2, args.w)
    if not bundle:
        raise SystemExit("Could not find 0th-iteration input bundle for {} Q{}W{}".format(args.particle, args.q2, args.w))
    summary = run_bg_systematics(bundle, args.repo_root, args.outpath, args.particle, args.q2, args.w)
    write_bg_systematics_summary(summary, args.outpath, args.particle, args.q2, args.w)
