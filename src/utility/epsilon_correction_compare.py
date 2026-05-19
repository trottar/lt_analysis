#! /usr/bin/python

from __future__ import annotations

import csv
import json
import math
import os

from background_config import EMP_EPS_DIFF_FAIL_THRESHOLD, EMP_EPS_DIFF_WARN_THRESHOLD
from frozen_manifest import get_analysis_artifact_paths, write_json_with_aliases


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def _safe_div(numerator, denominator):
    denominator = float(denominator)
    if abs(denominator) <= 1.0e-12:
        return 0.0
    return float(numerator) / denominator


def _iter_bin_entries(ledger_payload):
    for setting in ledger_payload.get("settings", []):
        phi_setting = setting.get("phi_setting")
        for entry in setting.get("t_phi_bins", []):
            stage_yields = entry.get("stage_yields", {})
            after_pion = float(stage_yields.get("after_pion_subtraction") or 0.0)
            final_yield = float(stage_yields.get("final_lambda_window") or 0.0)
            yield {
                "phi_setting": phi_setting,
                "bin_key": entry.get("bin_key"),
                "after_pion": after_pion,
                "final_yield": final_yield,
                "f_emp": _safe_div(after_pion - final_yield, after_pion),
            }


def build_epsilon_empirical_compare(low_ledger, high_ledger):
    low_bins = {
        (entry["phi_setting"], entry["bin_key"]): entry
        for entry in _iter_bin_entries(low_ledger)
    }
    high_bins = {
        (entry["phi_setting"], entry["bin_key"]): entry
        for entry in _iter_bin_entries(high_ledger)
    }
    shared_keys = sorted(set(low_bins.keys()) & set(high_bins.keys()))
    rows = []
    deltas = []
    warnings = []
    failures = []
    for phi_setting, bin_key in shared_keys:
        low_entry = low_bins[(phi_setting, bin_key)]
        high_entry = high_bins[(phi_setting, bin_key)]
        delta = float(high_entry["f_emp"] - low_entry["f_emp"])
        deltas.append(delta)
        row = {
            "phi_setting": phi_setting,
            "bin_key": bin_key,
            "f_emp_low": low_entry["f_emp"],
            "f_emp_high": high_entry["f_emp"],
            "delta_f_emp": delta,
        }
        rows.append(row)
        if abs(delta) > float(EMP_EPS_DIFF_WARN_THRESHOLD):
            warnings.append(row)
        if EMP_EPS_DIFF_FAIL_THRESHOLD is not None and abs(delta) > float(EMP_EPS_DIFF_FAIL_THRESHOLD):
            failures.append(row)

    mean_low = _safe_div(sum(row["f_emp_low"] for row in rows), len(rows)) if rows else 0.0
    mean_high = _safe_div(sum(row["f_emp_high"] for row in rows), len(rows)) if rows else 0.0
    rms_delta = math.sqrt(sum(delta * delta for delta in deltas) / len(deltas)) if deltas else 0.0
    max_abs_delta = max((abs(delta) for delta in deltas), default=0.0)
    return {
        "particle_type": low_ledger.get("particle_type"),
        "q2": low_ledger.get("q2"),
        "w": low_ledger.get("w"),
        "active_profile": high_ledger.get("active_profile") or low_ledger.get("active_profile"),
        "rows": rows,
        "mean_f_emp_low": float(mean_low),
        "mean_f_emp_high": float(mean_high),
        "rms_delta_f_emp": float(rms_delta),
        "max_abs_delta_f_emp": float(max_abs_delta),
        "warning_threshold": EMP_EPS_DIFF_WARN_THRESHOLD,
        "fail_threshold": EMP_EPS_DIFF_FAIL_THRESHOLD,
        "warning_bins": warnings,
        "failure_bins": failures,
    }


def write_epsilon_empirical_compare(payload, outpath, particle_type, q2, w, active_profile=None):
    paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    written = []
    written.extend(write_json_with_aliases(payload, paths["epsilon_compare_json"]))

    fieldnames = [
        "phi_setting",
        "bin_key",
        "f_emp_low",
        "f_emp_high",
        "delta_f_emp",
    ]
    with open(paths["epsilon_compare_csv"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in payload.get("rows", []):
            writer.writerow({key: row.get(key) for key in fieldnames})
    written.append(paths["epsilon_compare_csv"])
    return written


def load_ledgers_and_compare(low_ledger_path, high_ledger_path, outpath=None, particle_type=None, q2=None, w=None, active_profile=None):
    low_ledger = _load_json(low_ledger_path)
    high_ledger = _load_json(high_ledger_path)
    payload = build_epsilon_empirical_compare(low_ledger, high_ledger)
    if outpath and particle_type and q2 and w:
        write_epsilon_empirical_compare(payload, outpath, particle_type, q2, w, active_profile=active_profile)
    return payload


def validate_epsilon_compare(payload):
    if payload.get("failure_bins"):
        raise ValueError(
            "Empirical epsilon-difference threshold exceeded in {} bins".format(
                len(payload.get("failure_bins", []))
            )
        )
    for warning in payload.get("warning_bins", []):
        print(
            "WARNING: Large empirical epsilon difference for {} {} "
            "(delta_f_emp={:.4f})".format(
                warning.get("phi_setting"),
                warning.get("bin_key"),
                float(warning.get("delta_f_emp", 0.0) or 0.0),
            )
        )
    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare empirical low/high epsilon corrections from correction ledgers.")
    parser.add_argument("low_ledger")
    parser.add_argument("high_ledger")
    parser.add_argument("--outpath")
    parser.add_argument("--particle")
    parser.add_argument("--q2")
    parser.add_argument("--w")
    parser.add_argument("--profile")
    args = parser.parse_args()

    result = load_ledgers_and_compare(
        args.low_ledger,
        args.high_ledger,
        outpath=args.outpath,
        particle_type=args.particle,
        q2=args.q2,
        w=args.w,
        active_profile=args.profile,
    )
    validate_epsilon_compare(result)
