#! /usr/bin/python

from __future__ import annotations

import csv
import json
import os
import shutil

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

from frozen_manifest import (
    get_analysis_artifact_paths,
    resolve_profile_output_paths,
    write_json_with_aliases,
)


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


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


def _prediction_has_supported_correction(prediction):
    return any(
        key in prediction
        for key in (
            "predicted_total_correction",
            "predicted_fit1_correction",
            "predicted_fit2_correction",
        )
    )


def _validate_json_file(path, label):
    if not path:
        raise ValueError("Missing required '{}' path in non-KLambda cross-check config".format(label))
    if not os.path.exists(path):
        raise FileNotFoundError("Non-KLambda cross-check {} file was not found: {}".format(label, path))
    try:
        _load_json(path)
    except Exception as exc:
        raise ValueError("Could not read non-KLambda cross-check {} JSON {}: {}".format(label, path, exc))


def validate_nonklambda_config(config_payload):
    if not isinstance(config_payload, dict):
        raise ValueError("Non-KLambda cross-check config must be a JSON object")

    _validate_json_file(config_payload.get("low_ledger"), "low_ledger")
    _validate_json_file(config_payload.get("high_ledger"), "high_ledger")

    analysis_manifest = config_payload.get("analysis_manifest")
    if analysis_manifest and not os.path.exists(analysis_manifest):
        raise FileNotFoundError("Non-KLambda analysis manifest was not found: {}".format(analysis_manifest))

    channels = config_payload.get("channels")
    if channels is None:
        raise ValueError("Non-KLambda cross-check config must define 'channels'")
    if not isinstance(channels, list):
        raise ValueError("Non-KLambda cross-check 'channels' must be a list")

    warnings = []
    if not channels:
        warnings.append("Non-KLambda cross-check channel list is empty")

    for channel_index, channel in enumerate(channels):
        if not isinstance(channel, dict):
            raise ValueError("Channel entry {} must be a JSON object".format(channel_index))
        if not channel.get("channel_name"):
            raise ValueError("Channel entry {} is missing 'channel_name'".format(channel_index))
        predictions = channel.get("predictions")
        if predictions is None:
            raise ValueError("Channel '{}' is missing 'predictions'".format(channel["channel_name"]))
        if not isinstance(predictions, list):
            raise ValueError("Channel '{}' predictions must be a list".format(channel["channel_name"]))
        for prediction_index, prediction in enumerate(predictions):
            if not isinstance(prediction, dict):
                raise ValueError(
                    "Prediction {} in channel '{}' must be a JSON object".format(
                        prediction_index,
                        channel["channel_name"],
                    )
                )
            for required_key in ("epsilon_setting", "phi_setting", "bin_key"):
                if prediction.get(required_key) in (None, ""):
                    raise ValueError(
                        "Prediction {} in channel '{}' is missing '{}'".format(
                            prediction_index,
                            channel["channel_name"],
                            required_key,
                        )
                    )
            if not _prediction_has_supported_correction(prediction):
                raise ValueError(
                    "Prediction {} in channel '{}' must include at least one of "
                    "'predicted_total_correction', 'predicted_fit1_correction', or "
                    "'predicted_fit2_correction'".format(
                        prediction_index,
                        channel["channel_name"],
                    )
                )
    return warnings


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


def _float_or_none(value):
    if value in (None, ""):
        return None
    return float(value)


def _predicted_total_correction(prediction):
    if prediction.get("predicted_total_correction") not in (None, ""):
        return float(prediction.get("predicted_total_correction"))
    fit1 = _float_or_none(prediction.get("predicted_fit1_correction"))
    fit2 = _float_or_none(prediction.get("predicted_fit2_correction"))
    if fit1 is None and fit2 is None:
        return None
    return float((fit1 or 0.0) + (fit2 or 0.0))


def build_nonklambda_crosscheck(config_payload, active_profile=None):
    warnings = validate_nonklambda_config(config_payload)
    low_ledger = _load_json(config_payload["low_ledger"])
    high_ledger = _load_json(config_payload["high_ledger"])
    empirical_map = _load_empirical_bin_map(low_ledger, high_ledger)
    rows = []
    prediction_count = 0
    missing_empirical_bin_count = 0

    for channel in config_payload.get("channels", []):
        channel_name = channel.get("channel_name")
        label = channel.get("systematic_label")
        normalization_uncertainty = channel.get("normalization_uncertainty")
        for prediction in channel.get("predictions", []):
            prediction_count += 1
            epsset = str(prediction.get("epsilon_setting")).lower()
            phi_setting = prediction.get("phi_setting")
            bin_key = prediction.get("bin_key")
            empirical = empirical_map.get((epsset, phi_setting, bin_key))
            predicted_total = _predicted_total_correction(prediction)
            predicted_fit1 = _float_or_none(prediction.get("predicted_fit1_correction"))
            predicted_fit2 = _float_or_none(prediction.get("predicted_fit2_correction"))
            row = {
                "channel_name": channel_name,
                "systematic_label": label,
                "epsilon_setting": epsset,
                "phi_setting": phi_setting,
                "bin_key": bin_key,
                "normalization_uncertainty": normalization_uncertainty,
                "predicted_total_correction": predicted_total,
                "predicted_fit1_correction": predicted_fit1,
                "predicted_fit2_correction": predicted_fit2,
                "shape_chi2_ndf": _shape_chi2_ndf(
                    prediction.get("empirical_shape"),
                    prediction.get("predicted_shape"),
                ),
                "comparison_mode": "preprocessed_prediction",
                "active_profile": active_profile or config_payload.get("active_profile") or high_ledger.get("active_profile") or low_ledger.get("active_profile"),
            }
            if empirical is None:
                missing_empirical_bin_count += 1
                row.update(
                    {
                        "status": "missing_empirical_bin",
                        "empirical_total_correction": None,
                        "empirical_fit1_correction": None,
                        "empirical_fit2_correction": None,
                        "lambda_window_correction_difference": None,
                        "empirical_fit1_minus_prediction": None,
                        "empirical_fit2_minus_prediction": None,
                        "integral_difference": None,
                    }
                )
            else:
                empirical_total = _float_or_none(empirical.get("total_fraction"))
                empirical_fit1 = _float_or_none(empirical.get("fit1_fraction"))
                empirical_fit2 = _float_or_none(empirical.get("fit2_fraction"))
                row.update(
                    {
                        "status": "ok",
                        "empirical_total_correction": empirical_total,
                        "empirical_fit1_correction": empirical_fit1,
                        "empirical_fit2_correction": empirical_fit2,
                        "lambda_window_correction_difference": None
                        if empirical_total is None or predicted_total is None
                        else empirical_total - predicted_total,
                        "empirical_fit1_minus_prediction": None
                        if empirical_fit1 is None or predicted_fit1 is None
                        else empirical_fit1 - predicted_fit1,
                        "empirical_fit2_minus_prediction": None
                        if empirical_fit2 is None or predicted_fit2 is None
                        else empirical_fit2 - predicted_fit2,
                        "integral_difference": None
                        if empirical_total is None or predicted_total is None
                        else empirical_total - predicted_total,
                    }
                )
            rows.append(row)

    if warnings:
        for warning in warnings:
            print("WARNING: {}".format(warning))

    active_profile = active_profile or config_payload.get("active_profile") or high_ledger.get("active_profile") or low_ledger.get("active_profile")
    summary = {
        "comparison_mode": "preprocessed_prediction",
        "active_profile": active_profile,
        "analysis_manifest": config_payload.get("analysis_manifest"),
        "low_ledger_active_profile": low_ledger.get("active_profile"),
        "high_ledger_active_profile": high_ledger.get("active_profile"),
        "channel_count": len(config_payload.get("channels", [])),
        "prediction_count": prediction_count,
        "missing_empirical_bin_count": missing_empirical_bin_count,
        "rows": rows,
        "warnings": warnings,
        "max_abs_lambda_window_correction_difference": max(
            (
                abs(float(row.get("lambda_window_correction_difference", 0.0) or 0.0))
                for row in rows
                if row.get("lambda_window_correction_difference") is not None
            ),
            default=0.0,
        ),
    }
    return summary


def write_nonklambda_crosscheck(summary, outpath, particle_type, q2, w, active_profile=None):
    paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    written = []
    written.extend(
        write_json_with_aliases(
            summary,
            paths["nonklambda_json"],
            paths["nonklambda_json_profile"],
            active_profile=active_profile,
        )
    )

    fieldnames = [
        "channel_name",
        "systematic_label",
        "epsilon_setting",
        "phi_setting",
        "bin_key",
        "status",
        "normalization_uncertainty",
        "predicted_total_correction",
        "predicted_fit1_correction",
        "predicted_fit2_correction",
        "empirical_total_correction",
        "empirical_fit1_correction",
        "empirical_fit2_correction",
        "lambda_window_correction_difference",
        "empirical_fit1_minus_prediction",
        "empirical_fit2_minus_prediction",
        "shape_chi2_ndf",
        "integral_difference",
        "comparison_mode",
        "active_profile",
    ]
    primary_csv_path, alias_csv_path = resolve_profile_output_paths(
        paths["nonklambda_csv"],
        profile_path=paths["nonklambda_csv_profile"],
        active_profile=active_profile,
    )
    with open(primary_csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary.get("rows", []):
            writer.writerow({key: row.get(key) for key in fieldnames})
    written.append(primary_csv_path)
    if alias_csv_path and os.path.abspath(alias_csv_path) != os.path.abspath(primary_csv_path):
        shutil.copy(primary_csv_path, alias_csv_path)
        written.append(alias_csv_path)

    if plt is not None:
        primary_pdf_path, alias_pdf_path = resolve_profile_output_paths(
            paths["nonklambda_pdf"],
            profile_path=paths["nonklambda_pdf_profile"],
            active_profile=active_profile,
        )
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_title("Non-KLambda cross-check")
        ax.set_xlabel("Row index")
        ax.set_ylabel("Lambda-window correction difference")
        yvals = [
            float(row.get("lambda_window_correction_difference", 0.0) or 0.0)
            for row in summary.get("rows", [])
            if row.get("lambda_window_correction_difference") is not None
        ]
        if yvals:
            ax.plot(range(len(yvals)), yvals, marker="o", linestyle="-")
        else:
            ax.text(0.5, 0.5, "No comparable empirical bins", ha="center", va="center", transform=ax.transAxes)
        ax.axhline(0.0, color="black", linewidth=1.0)
        fig.tight_layout()
        fig.savefig(primary_pdf_path)
        plt.close(fig)
        written.append(primary_pdf_path)
        if alias_pdf_path and os.path.abspath(alias_pdf_path) != os.path.abspath(primary_pdf_path):
            shutil.copy(primary_pdf_path, alias_pdf_path)
            written.append(alias_pdf_path)
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
    result = build_nonklambda_crosscheck(payload, active_profile=args.profile)
    if args.outpath and args.particle and args.q2 and args.w:
        write_nonklambda_crosscheck(result, args.outpath, args.particle, args.q2, args.w, active_profile=args.profile)
