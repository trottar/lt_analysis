#! /usr/bin/python

from __future__ import annotations

import csv
import json
import math
import os

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

from frozen_manifest import get_analysis_artifact_paths, write_json_with_aliases


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def _safe_div(numerator, denominator):
    denominator = float(denominator)
    if abs(denominator) <= 1.0e-12:
        return 0.0
    return float(numerator) / denominator


def _shape_chi2_ndf(empirical_shape, predicted_shape):
    if not empirical_shape or not predicted_shape or len(empirical_shape) != len(predicted_shape):
        return None
    chi2 = 0.0
    n_used = 0
    for empirical, predicted in zip(empirical_shape, predicted_shape):
        emp_val = float(empirical.get("value", 0.0) or 0.0)
        pred_val = float(predicted.get("value", 0.0) or 0.0)
        variance = float(empirical.get("variance", abs(emp_val) + abs(pred_val)) or (abs(emp_val) + abs(pred_val)))
        if variance <= 0.0:
            continue
        chi2 += ((emp_val - pred_val) ** 2) / variance
        n_used += 1
    if n_used <= 1:
        return None
    return float(chi2 / (n_used - 1))


def _load_empirical_bin_map(low_ledger, high_ledger):
    result = {}
    for eps_name, ledger in (("low", low_ledger), ("high", high_ledger)):
        for setting in ledger.get("settings", []):
            for entry in setting.get("t_phi_bins", []):
                result[(eps_name, setting.get("phi_setting"), entry.get("bin_key"))] = {
                    "fit1_fraction": entry.get("fit1_fractional_correction"),
                    "fit2_fraction": entry.get("fit2_fractional_correction"),
                    "total_fraction": entry.get("total_empirical_fraction"),
                }
    return result


def build_nonklambda_crosscheck(config_payload):
    low_ledger = _load_json(config_payload["low_ledger"])
    high_ledger = _load_json(config_payload["high_ledger"])
    empirical_map = _load_empirical_bin_map(low_ledger, high_ledger)
    rows = []
    for channel in config_payload.get("channels", []):
        channel_name = channel.get("channel_name")
        label = channel.get("systematic_label")
        normalization_uncertainty = channel.get("normalization_uncertainty")
        for prediction in channel.get("predictions", []):
            epsset = str(prediction.get("epsilon_setting")).lower()
            phi_setting = prediction.get("phi_setting")
            bin_key = prediction.get("bin_key")
            empirical = empirical_map.get((epsset, phi_setting, bin_key), {})
            predicted_total = float(prediction.get("predicted_total_correction", 0.0) or 0.0)
            empirical_total = float(empirical.get("total_fraction", 0.0) or 0.0)
            empirical_fit1 = float(empirical.get("fit1_fraction", 0.0) or 0.0)
            empirical_fit2 = float(empirical.get("fit2_fraction", 0.0) or 0.0)
            row = {
                "channel_name": channel_name,
                "systematic_label": label,
                "epsilon_setting": epsset,
                "phi_setting": phi_setting,
                "bin_key": bin_key,
                "normalization_uncertainty": normalization_uncertainty,
                "predicted_total_correction": predicted_total,
                "empirical_total_correction": empirical_total,
                "lambda_window_correction_difference": empirical_total - predicted_total,
                "empirical_fit1_minus_prediction": empirical_fit1 - float(prediction.get("predicted_fit1_correction", 0.0) or 0.0),
                "empirical_fit2_minus_prediction": empirical_fit2 - float(prediction.get("predicted_fit2_correction", 0.0) or 0.0),
                "shape_chi2_ndf": _shape_chi2_ndf(
                    prediction.get("empirical_shape"),
                    prediction.get("predicted_shape"),
                ),
                "integral_difference": empirical_total - predicted_total,
            }
            rows.append(row)
    summary = {
        "analysis_manifest": config_payload.get("analysis_manifest"),
        "rows": rows,
        "channel_count": len(config_payload.get("channels", [])),
        "max_abs_lambda_window_correction_difference": max(
            (abs(float(row.get("lambda_window_correction_difference", 0.0) or 0.0)) for row in rows),
            default=0.0,
        ),
    }
    return summary


def write_nonklambda_crosscheck(summary, outpath, particle_type, q2, w, active_profile=None):
    paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    written = []
    written.extend(write_json_with_aliases(summary, paths["nonklambda_json"]))

    fieldnames = [
        "channel_name",
        "systematic_label",
        "epsilon_setting",
        "phi_setting",
        "bin_key",
        "normalization_uncertainty",
        "predicted_total_correction",
        "empirical_total_correction",
        "lambda_window_correction_difference",
        "empirical_fit1_minus_prediction",
        "empirical_fit2_minus_prediction",
        "shape_chi2_ndf",
        "integral_difference",
    ]
    with open(paths["nonklambda_csv"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary.get("rows", []):
            writer.writerow({key: row.get(key) for key in fieldnames})
    written.append(paths["nonklambda_csv"])

    if plt is not None:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_title("Non-KLambda cross-check")
        ax.set_xlabel("Row index")
        ax.set_ylabel("Lambda-window correction difference")
        yvals = [float(row.get("lambda_window_correction_difference", 0.0) or 0.0) for row in summary.get("rows", [])]
        ax.plot(range(len(yvals)), yvals, marker="o", linestyle="-")
        ax.axhline(0.0, color="black", linewidth=1.0)
        fig.tight_layout()
        fig.savefig(paths["nonklambda_pdf"])
        plt.close(fig)
        written.append(paths["nonklambda_pdf"])
    return written


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run the non-KLambda cross-check utility.")
    parser.add_argument("config_json")
    parser.add_argument("--outpath")
    parser.add_argument("--particle")
    parser.add_argument("--q2")
    parser.add_argument("--w")
    parser.add_argument("--profile")
    args = parser.parse_args()

    payload = _load_json(args.config_json)
    result = build_nonklambda_crosscheck(payload)
    if args.outpath and args.particle and args.q2 and args.w:
        write_nonklambda_crosscheck(result, args.outpath, args.particle, args.q2, args.w, active_profile=args.profile)
