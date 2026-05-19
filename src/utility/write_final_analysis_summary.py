#! /usr/bin/python

from __future__ import annotations

import csv
import hashlib
import json
import os
from copy import deepcopy

from frozen_manifest import get_analysis_artifact_paths, write_json_with_aliases


def _read_json_if_exists(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r") as handle:
        return json.load(handle)


def _read_csv_if_exists(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def _read_dat_if_exists(path):
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            rows.append(stripped.split())
    return rows


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_final_analysis_summary(
    manifest_payload,
    low_ledger_payload,
    high_ledger_payload,
    epsilon_compare_payload=None,
    systematics_payload=None,
    nonklambda_payload=None,
    xsect_payload=None,
    current_inp_dict=None,
    manifest_path=None,
):
    current_inp_dict = current_inp_dict or {}
    summary = {
        "particle_type": manifest_payload.get("particle_type"),
        "q2": manifest_payload.get("q2"),
        "w": manifest_payload.get("w"),
        "mm_cut_window": manifest_payload.get("mm_cut_window"),
        "t_bin_edges": manifest_payload.get("t_bin_edges"),
        "phi_bin_edges": manifest_payload.get("phi_bin_edges"),
        "active_profile": manifest_payload.get("active_profile"),
        "resolved_optimizer_settings": deepcopy(manifest_payload.get("resolved_optimizer_settings", {})),
        "fit1_scales": deepcopy(manifest_payload.get("selected_bg_scale1s", {})),
        "fit2_scales": deepcopy(manifest_payload.get("selected_bg_scale2s", {})),
        "low_epsilon_corrections": deepcopy((low_ledger_payload or {}).get("combined_totals", {})),
        "high_epsilon_corrections": deepcopy((high_ledger_payload or {}).get("combined_totals", {})),
        "epsilon_empirical_compare": deepcopy(epsilon_compare_payload or {}),
        "systematics_summary": deepcopy(systematics_payload or {}),
        "nonklambda_crosscheck": deepcopy(nonklambda_payload or {}),
        "xsect_outputs": deepcopy(xsect_payload or {}),
        "manifest_hash": _sha256_file(manifest_path) if manifest_path and os.path.exists(manifest_path) else None,
        "iteration_count": current_inp_dict.get("iter_num", 0),
        "convergence_status": current_inp_dict.get("convergence_status", "unknown"),
    }
    return summary


def _markdown_lines(summary):
    lines = [
        "# Final Analysis Summary",
        "",
        "| Field | Value |",
        "| --- | --- |",
        "| Particle | {} |".format(summary.get("particle_type", "")),
        "| Q2/W | Q{} W{} |".format(summary.get("q2", ""), summary.get("w", "")),
        "| Active profile | {} |".format(summary.get("active_profile", "")),
        "| MM cut | {} |".format(summary.get("mm_cut_window", "")),
        "| t bins | {} |".format(summary.get("t_bin_edges", "")),
        "| phi bins | {} |".format(summary.get("phi_bin_edges", "")),
        "| Manifest hash | {} |".format(summary.get("manifest_hash", "")),
        "| Iteration count | {} |".format(summary.get("iteration_count", "")),
        "| Convergence status | {} |".format(summary.get("convergence_status", "")),
        "",
        "## Empirical Corrections",
        "",
        "| Epsilon | Fit 1 frac | Fit 2 frac | Total empirical frac | Fit uncertainty frac |",
        "| --- | ---: | ---: | ---: | ---: |",
        "| Low | {fit1:.6f} | {fit2:.6f} | {total:.6f} | {unc:.6f} |".format(
            fit1=float(summary.get("low_epsilon_corrections", {}).get("fit1_fractional_correction", 0.0) or 0.0),
            fit2=float(summary.get("low_epsilon_corrections", {}).get("fit2_fractional_correction", 0.0) or 0.0),
            total=float(summary.get("low_epsilon_corrections", {}).get("total_empirical_fraction", 0.0) or 0.0),
            unc=float(summary.get("low_epsilon_corrections", {}).get("empirical_fit_uncertainty_frac", 0.0) or 0.0),
        ),
        "| High | {fit1:.6f} | {fit2:.6f} | {total:.6f} | {unc:.6f} |".format(
            fit1=float(summary.get("high_epsilon_corrections", {}).get("fit1_fractional_correction", 0.0) or 0.0),
            fit2=float(summary.get("high_epsilon_corrections", {}).get("fit2_fractional_correction", 0.0) or 0.0),
            total=float(summary.get("high_epsilon_corrections", {}).get("total_empirical_fraction", 0.0) or 0.0),
            unc=float(summary.get("high_epsilon_corrections", {}).get("empirical_fit_uncertainty_frac", 0.0) or 0.0),
        ),
        "",
    ]
    epsilon_compare = summary.get("epsilon_empirical_compare", {})
    if epsilon_compare:
        lines.extend(
            [
                "## High/Low Epsilon Compare",
                "",
                "| Mean low | Mean high | RMS delta | Max abs delta |",
                "| ---: | ---: | ---: | ---: |",
                "| {low:.6f} | {high:.6f} | {rms:.6f} | {maxabs:.6f} |".format(
                    low=float(epsilon_compare.get("mean_f_emp_low", 0.0) or 0.0),
                    high=float(epsilon_compare.get("mean_f_emp_high", 0.0) or 0.0),
                    rms=float(epsilon_compare.get("rms_delta_f_emp", 0.0) or 0.0),
                    maxabs=float(epsilon_compare.get("max_abs_delta_f_emp", 0.0) or 0.0),
                ),
                "",
            ]
        )
    if summary.get("xsect_outputs"):
        lines.extend(
            [
                "## Cross Section Files",
                "",
                "```json",
                json.dumps(summary.get("xsect_outputs"), indent=2, sort_keys=True),
                "```",
            ]
        )
    return "\n".join(lines) + "\n"


def write_final_analysis_summary(summary, outpath, particle_type, q2, w, active_profile=None):
    artifact_paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    written = []
    written.extend(write_json_with_aliases(summary, artifact_paths["final_summary_json"]))

    csv_rows = [
        {
            "particle_type": summary.get("particle_type"),
            "q2": summary.get("q2"),
            "w": summary.get("w"),
            "active_profile": summary.get("active_profile"),
            "manifest_hash": summary.get("manifest_hash"),
            "iteration_count": summary.get("iteration_count"),
            "convergence_status": summary.get("convergence_status"),
            "low_fit1_fractional_correction": summary.get("low_epsilon_corrections", {}).get("fit1_fractional_correction"),
            "low_fit2_fractional_correction": summary.get("low_epsilon_corrections", {}).get("fit2_fractional_correction"),
            "low_total_empirical_fraction": summary.get("low_epsilon_corrections", {}).get("total_empirical_fraction"),
            "high_fit1_fractional_correction": summary.get("high_epsilon_corrections", {}).get("fit1_fractional_correction"),
            "high_fit2_fractional_correction": summary.get("high_epsilon_corrections", {}).get("fit2_fractional_correction"),
            "high_total_empirical_fraction": summary.get("high_epsilon_corrections", {}).get("total_empirical_fraction"),
            "epsilon_compare_rms_delta_f_emp": summary.get("epsilon_empirical_compare", {}).get("rms_delta_f_emp"),
            "epsilon_compare_max_abs_delta_f_emp": summary.get("epsilon_empirical_compare", {}).get("max_abs_delta_f_emp"),
        }
    ]
    with open(artifact_paths["final_summary_csv"], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        for row in csv_rows:
            writer.writerow(row)
    written.append(artifact_paths["final_summary_csv"])

    with open(artifact_paths["final_summary_md"], "w") as handle:
        handle.write(_markdown_lines(summary))
    written.append(artifact_paths["final_summary_md"])
    return written


def collect_xsect_payload(outpath, ltanapath, particle_type, q2, w, anatype, pol):
    pol_str = "pl" if str(pol) == "1" else "mn"
    payload = {}
    payload["separated_csv"] = _read_csv_if_exists(
        os.path.join(outpath, "{}LT_Q{}W{}.csv".format(anatype, q2.replace(".", "p"), w.replace(".", "p")))
    )
    payload["unseparated_csv"] = _read_csv_if_exists(
        os.path.join(ltanapath, "src", particle_type, "xsects", "unsep_Q{}W{}.csv".format(q2.replace("p", ""), w.replace("p", "")))
    )
    payload["x_sep_dat"] = _read_dat_if_exists(
        os.path.join(ltanapath, "src", particle_type, "xsects", "x_sep.{}_Q{}W{}.dat".format(
            pol_str,
            q2.replace("p", ""),
            w.replace("p", ""),
        ))
    )
    return payload
