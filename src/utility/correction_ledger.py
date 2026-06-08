#! /usr/bin/python

from __future__ import annotations

import csv
import json
import math
import os
from copy import deepcopy

from utility import integrate_hist_range
from frozen_manifest import get_correction_ledger_paths, write_json_with_aliases
from background_config import get_particle_subtraction_window_config


STAGE_ORDER = [
    "raw_prompt",
    "after_random_subtraction",
    "after_dummy_subtraction",
    "after_pion_subtraction",
    "after_empirical_fit1",
    "after_empirical_fit2",
    "final_lambda_window",
]

TOP_LEVEL_STAGE_HISTS = {
    "raw_prompt": "H_MM_rand_dummy_DATA",
    "after_random_subtraction": "H_MM_dummy_DATA",
    "after_dummy_subtraction": "H_MM_nosub_DATA",
    "after_pion_subtraction": "H_MM_pisub_DATA",
    "after_empirical_fit1": "H_MM_fit1sub_DATA",
    "after_empirical_fit2": "H_MM_full_DATA",
    "final_lambda_window": "H_MM_full_DATA",
}


def _safe_div(numerator, denominator):
    denominator = float(denominator)
    if abs(denominator) <= 1.0e-12:
        return 0.0
    return float(numerator) / denominator


def _pct_change(current_value, previous_value):
    return 100.0 * _safe_div(current_value - previous_value, previous_value)


def _hist_window_yield(hist, mm_min, mm_max):
    if hist is None:
        return None
    return float(integrate_hist_range(hist, mm_min, mm_max))


def _top_level_stage_yields(hist, mm_min, mm_max):
    return {
        stage_name: _hist_window_yield(hist.get(hist_key), mm_min, mm_max)
        for stage_name, hist_key in TOP_LEVEL_STAGE_HISTS.items()
    }


def _percent_change_map(stage_yields):
    pct = {}
    previous_value = None
    for stage_name in STAGE_ORDER:
        current_value = stage_yields.get(stage_name)
        if current_value is None or previous_value is None:
            pct[stage_name] = None
        else:
            pct[stage_name] = _pct_change(current_value, previous_value)
        previous_value = current_value if current_value is not None else previous_value
    return pct


def _empirical_corrections(stage_yields):
    after_pion = float(stage_yields.get("after_pion_subtraction") or 0.0)
    after_fit1 = float(stage_yields.get("after_empirical_fit1") or 0.0)
    final_yield = float(stage_yields.get("final_lambda_window") or 0.0)
    fit1_correction = after_pion - after_fit1
    fit2_correction = after_fit1 - final_yield
    total_correction = after_pion - final_yield
    return {
        "fit1_correction": float(fit1_correction),
        "fit2_correction": float(fit2_correction),
        "total_empirical_correction": float(total_correction),
        "fit1_fractional_correction": _safe_div(fit1_correction, after_pion),
        "fit2_fractional_correction": _safe_div(fit2_correction, after_fit1),
        "total_empirical_fraction": _safe_div(total_correction, after_pion),
    }


def _aggregate_oversub_diagnostics(entries):
    summary = {
        "fit1": {
            "oversub_bin_count": 0,
            "max_unclamped_ratio": 0.0,
            "oversub_integral": 0.0,
            "affected_lambda_fraction": 0.0,
            "affected_mm_range": None,
        },
        "fit2": {
            "oversub_bin_count": 0,
            "max_unclamped_ratio": 0.0,
            "oversub_integral": 0.0,
            "affected_lambda_fraction": 0.0,
            "affected_mm_range": None,
        },
    }
    for entry in entries:
        diagnostics = entry.get("oversub_diagnostics", {})
        for fit_name in ("fit1", "fit2"):
            fit_diag = diagnostics.get(fit_name) or {}
            current = summary[fit_name]
            current["oversub_bin_count"] += int(fit_diag.get("oversub_bin_count", 0) or 0)
            current["max_unclamped_ratio"] = max(
                float(current.get("max_unclamped_ratio", 0.0) or 0.0),
                float(fit_diag.get("max_unclamped_ratio", 0.0) or 0.0),
            )
            current["oversub_integral"] += float(fit_diag.get("oversub_integral", 0.0) or 0.0)
            current["affected_lambda_fraction"] = max(
                float(current.get("affected_lambda_fraction", 0.0) or 0.0),
                float(fit_diag.get("affected_lambda_fraction", 0.0) or 0.0),
            )
            affected_range = fit_diag.get("affected_mm_range")
            if affected_range:
                if current["affected_mm_range"] is None:
                    current["affected_mm_range"] = list(affected_range)
                else:
                    current["affected_mm_range"][0] = min(current["affected_mm_range"][0], float(affected_range[0]))
                    current["affected_mm_range"][1] = max(current["affected_mm_range"][1], float(affected_range[1]))
    return summary


def _build_bin_entry(bin_key, entry):
    stage_yields = deepcopy(entry.get("stage_window_yields", {}))
    empirical = _empirical_corrections(stage_yields)
    fit1_frac_err = float(entry.get("bg_fit1_frac_err", 0.0) or 0.0)
    fit2_frac_err = float(entry.get("bg_fit2_frac_err", 0.0) or 0.0)
    empirical_unc_abs = math.sqrt(
        (float(empirical["fit1_correction"]) * fit1_frac_err) ** 2
        + (float(empirical["fit2_correction"]) * fit2_frac_err) ** 2
    )
    final_yield = float(stage_yields.get("final_lambda_window") or 0.0)
    fit1_zero_bg = fit1_frac_err >= 1.0 and abs(float(empirical["fit1_correction"])) <= 1.0e-12
    fit2_zero_bg = fit2_frac_err >= 1.0 and abs(float(empirical["fit2_correction"])) <= 1.0e-12
    return {
        "bin_key": bin_key,
        "stage_yields": stage_yields,
        "percent_change_by_stage": _percent_change_map(stage_yields),
        "scale_factor": float(entry.get("scale_factor", 0.0) or 0.0),
        "bg_fit1_frac_err": fit1_frac_err,
        "bg_fit2_frac_err": fit2_frac_err,
        "empirical_fit_uncertainty_abs": float(empirical_unc_abs),
        "empirical_fit_uncertainty_frac": _safe_div(empirical_unc_abs, final_yield),
        "fit1_failed": bool(fit1_frac_err >= 1.0),
        "fit2_failed": bool(fit2_frac_err >= 1.0),
        "fit1_zero_background_fallback": bool(fit1_zero_bg),
        "fit2_zero_background_fallback": bool(fit2_zero_bg),
        "oversub_diagnostics": deepcopy(entry.get("oversub_diagnostics", {})),
        **empirical,
    }


def _combine_stage_yields(bin_entries):
    combined = {stage_name: 0.0 for stage_name in STAGE_ORDER}
    for bin_entry in bin_entries:
        for stage_name in STAGE_ORDER:
            combined[stage_name] += float(bin_entry["stage_yields"].get(stage_name) or 0.0)
    return combined


def _combine_setting_stage_yields(setting_entries):
    combined = {stage_name: 0.0 for stage_name in STAGE_ORDER}
    for setting_entry in setting_entries:
        stage_yields = setting_entry.get("combined_stage_yields", {})
        for stage_name in STAGE_ORDER:
            combined[stage_name] += float(stage_yields.get(stage_name) or 0.0)
    return combined


def _summarize_setting(hist, processed_dict, mm_min, mm_max):
    bin_level_available = bool(processed_dict)
    bin_entries = []
    for bin_key in sorted(processed_dict.keys()):
        bin_entries.append(_build_bin_entry(bin_key, processed_dict[bin_key]))
    combined_stage_yields = _combine_stage_yields(bin_entries) if bin_entries else _top_level_stage_yields(hist, mm_min, mm_max)
    empirical = _empirical_corrections(combined_stage_yields)
    total_unc_abs = math.sqrt(
        sum(float(bin_entry.get("empirical_fit_uncertainty_abs", 0.0) or 0.0) ** 2 for bin_entry in bin_entries)
    )
    final_yield = float(combined_stage_yields.get("final_lambda_window") or 0.0)
    return {
        "phi_setting": hist.get("phi_setting"),
        "bin_level_ledger_available": bool(bin_level_available),
        "overall_stage_yields": _top_level_stage_yields(hist, mm_min, mm_max),
        "combined_stage_yields": combined_stage_yields,
        "percent_change_by_stage": _percent_change_map(combined_stage_yields),
        "fit1_fractional_correction": empirical["fit1_fractional_correction"],
        "fit2_fractional_correction": empirical["fit2_fractional_correction"],
        "total_empirical_fraction": empirical["total_empirical_fraction"],
        "fit1_correction": empirical["fit1_correction"],
        "fit2_correction": empirical["fit2_correction"],
        "total_empirical_correction": empirical["total_empirical_correction"],
        "empirical_fit_uncertainty_abs": float(total_unc_abs),
        "empirical_fit_uncertainty_frac": _safe_div(total_unc_abs, final_yield),
        "failed_fit_count": int(
            sum(1 for entry in bin_entries if entry.get("fit1_failed") or entry.get("fit2_failed"))
        ),
        "zero_background_fallback_count": int(
            sum(
                1
                for entry in bin_entries
                if entry.get("fit1_zero_background_fallback") or entry.get("fit2_zero_background_fallback")
            )
        ),
        "oversub_diagnostics": _aggregate_oversub_diagnostics(bin_entries),
        "t_phi_bins": bin_entries,
    }


def build_correction_ledger(histlist, inp_dict, bg_summary=None):
    mm_min = float(inp_dict["mm_min"])
    mm_max = float(inp_dict["mm_max"])
    setting_entries = []
    for hist in histlist:
        processed_dict = hist.get("_yield_data_processed_dict", {})
        setting_entries.append(_summarize_setting(hist, processed_dict, mm_min, mm_max))

    combined_bin_entries = []
    for setting_entry in setting_entries:
        combined_bin_entries.extend(setting_entry.get("t_phi_bins", []))
    combined_stage_yields = (
        _combine_stage_yields(combined_bin_entries)
        if combined_bin_entries
        else _combine_setting_stage_yields(setting_entries)
    )
    combined_empirical = _empirical_corrections(combined_stage_yields)
    combined_unc_sources = combined_bin_entries or setting_entries
    combined_unc_abs = math.sqrt(
        sum(float(entry.get("empirical_fit_uncertainty_abs", 0.0) or 0.0) ** 2 for entry in combined_unc_sources)
    )
    final_yield = float(combined_stage_yields.get("final_lambda_window") or 0.0)

    payload = {
        "particle_type": inp_dict.get("ParticleType"),
        "epsset": inp_dict.get("EPSSET"),
        "q2": inp_dict.get("Q2"),
        "w": inp_dict.get("W"),
        "outfilename": inp_dict.get("OutFilename"),
        "mm_cut_window": [mm_min, mm_max],
        "particle_subtraction_windows": get_particle_subtraction_window_config(
            inp_dict.get("ParticleType"),
            "pion",
        ),
        "active_profile": inp_dict.get("bg_active_profile"),
        "resolved_optimizer_settings": deepcopy(inp_dict.get("bg_resolved_profile", {})),
        "bg_optimization_summary": deepcopy(bg_summary or inp_dict.get("bg_optimization_summary", {})),
        "settings": setting_entries,
        "combined_totals": {
            "combined_stage_yields": combined_stage_yields,
            "all_settings_have_bin_level_ledger": bool(
                setting_entries and all(setting_entry.get("bin_level_ledger_available") for setting_entry in setting_entries)
            ),
            "percent_change_by_stage": _percent_change_map(combined_stage_yields),
            "fit1_fractional_correction": combined_empirical["fit1_fractional_correction"],
            "fit2_fractional_correction": combined_empirical["fit2_fractional_correction"],
            "total_empirical_fraction": combined_empirical["total_empirical_fraction"],
            "fit1_correction": combined_empirical["fit1_correction"],
            "fit2_correction": combined_empirical["fit2_correction"],
            "total_empirical_correction": combined_empirical["total_empirical_correction"],
            "empirical_fit_uncertainty_abs": float(combined_unc_abs),
            "empirical_fit_uncertainty_frac": _safe_div(combined_unc_abs, final_yield),
            "failed_fit_count": int(
                sum(int(setting_entry.get("failed_fit_count", 0) or 0) for setting_entry in setting_entries)
            ),
            "zero_background_fallback_count": int(
                sum(int(setting_entry.get("zero_background_fallback_count", 0) or 0) for setting_entry in setting_entries)
            ),
            "oversub_diagnostics": _aggregate_oversub_diagnostics(combined_bin_entries),
        },
    }
    return payload


def _ledger_rows(payload):
    rows = []
    resolved_profile_json = json.dumps(payload.get("resolved_optimizer_settings", {}), sort_keys=True)
    for setting_entry in payload.get("settings", []):
        row = {
            "row_kind": "setting_total",
            "epsset": payload.get("epsset"),
            "phi_setting": setting_entry.get("phi_setting"),
            "bin_key": "",
            "active_profile": payload.get("active_profile"),
            "resolved_profile_json": resolved_profile_json,
        }
        row.update(setting_entry.get("combined_stage_yields", {}))
        row.update({
            "bin_level_ledger_available": setting_entry.get("bin_level_ledger_available"),
            "fit1_fractional_correction": setting_entry.get("fit1_fractional_correction"),
            "fit2_fractional_correction": setting_entry.get("fit2_fractional_correction"),
            "total_empirical_fraction": setting_entry.get("total_empirical_fraction"),
            "empirical_fit_uncertainty_frac": setting_entry.get("empirical_fit_uncertainty_frac"),
            "failed_fit_count": setting_entry.get("failed_fit_count"),
            "zero_background_fallback_count": setting_entry.get("zero_background_fallback_count"),
            "fit1_oversub_bin_count": setting_entry.get("oversub_diagnostics", {}).get("fit1", {}).get("oversub_bin_count"),
            "fit2_oversub_bin_count": setting_entry.get("oversub_diagnostics", {}).get("fit2", {}).get("oversub_bin_count"),
            "fit1_max_unclamped_ratio": setting_entry.get("oversub_diagnostics", {}).get("fit1", {}).get("max_unclamped_ratio"),
            "fit2_max_unclamped_ratio": setting_entry.get("oversub_diagnostics", {}).get("fit2", {}).get("max_unclamped_ratio"),
        })
        rows.append(row)
        for bin_entry in setting_entry.get("t_phi_bins", []):
            row = {
                "row_kind": "t_phi_bin",
                "epsset": payload.get("epsset"),
                "phi_setting": setting_entry.get("phi_setting"),
                "bin_key": bin_entry.get("bin_key"),
                "active_profile": payload.get("active_profile"),
                "resolved_profile_json": resolved_profile_json,
            }
            row.update(bin_entry.get("stage_yields", {}))
            row.update({
                "fit1_fractional_correction": bin_entry.get("fit1_fractional_correction"),
                "fit2_fractional_correction": bin_entry.get("fit2_fractional_correction"),
                "total_empirical_fraction": bin_entry.get("total_empirical_fraction"),
                "empirical_fit_uncertainty_frac": bin_entry.get("empirical_fit_uncertainty_frac"),
                "failed_fit_count": int(bool(bin_entry.get("fit1_failed") or bin_entry.get("fit2_failed"))),
                "zero_background_fallback_count": int(bool(bin_entry.get("fit1_zero_background_fallback") or bin_entry.get("fit2_zero_background_fallback"))),
                "fit1_oversub_bin_count": bin_entry.get("oversub_diagnostics", {}).get("fit1", {}).get("oversub_bin_count"),
                "fit2_oversub_bin_count": bin_entry.get("oversub_diagnostics", {}).get("fit2", {}).get("oversub_bin_count"),
                "fit1_max_unclamped_ratio": bin_entry.get("oversub_diagnostics", {}).get("fit1", {}).get("max_unclamped_ratio"),
                "fit2_max_unclamped_ratio": bin_entry.get("oversub_diagnostics", {}).get("fit2", {}).get("max_unclamped_ratio"),
            })
            rows.append(row)

    combined = payload.get("combined_totals", {})
    row = {
        "row_kind": "combined_total",
        "epsset": payload.get("epsset"),
        "phi_setting": "Combined",
        "bin_key": "",
        "active_profile": payload.get("active_profile"),
        "resolved_profile_json": resolved_profile_json,
    }
    row.update(combined.get("combined_stage_yields", {}))
    row.update({
        "all_settings_have_bin_level_ledger": combined.get("all_settings_have_bin_level_ledger"),
        "fit1_fractional_correction": combined.get("fit1_fractional_correction"),
        "fit2_fractional_correction": combined.get("fit2_fractional_correction"),
        "total_empirical_fraction": combined.get("total_empirical_fraction"),
        "empirical_fit_uncertainty_frac": combined.get("empirical_fit_uncertainty_frac"),
        "failed_fit_count": combined.get("failed_fit_count"),
        "zero_background_fallback_count": combined.get("zero_background_fallback_count"),
        "fit1_oversub_bin_count": combined.get("oversub_diagnostics", {}).get("fit1", {}).get("oversub_bin_count"),
        "fit2_oversub_bin_count": combined.get("oversub_diagnostics", {}).get("fit2", {}).get("oversub_bin_count"),
        "fit1_max_unclamped_ratio": combined.get("oversub_diagnostics", {}).get("fit1", {}).get("max_unclamped_ratio"),
        "fit2_max_unclamped_ratio": combined.get("oversub_diagnostics", {}).get("fit2", {}).get("max_unclamped_ratio"),
    })
    rows.append(row)
    return rows


def write_correction_ledger(payload, outpath, particle_type, outfilename, active_profile=None):
    paths = get_correction_ledger_paths(outpath, particle_type, outfilename, active_profile=active_profile)
    written = []
    written.extend(
        write_json_with_aliases(
            payload,
            paths["json"],
            paths["json_profile"],
            active_profile=active_profile,
        )
    )

    rows = _ledger_rows(payload)
    fieldnames = [
        "row_kind",
        "epsset",
        "phi_setting",
        "bin_key",
        "active_profile",
        "resolved_profile_json",
        *STAGE_ORDER,
        "bin_level_ledger_available",
        "all_settings_have_bin_level_ledger",
        "fit1_fractional_correction",
        "fit2_fractional_correction",
        "total_empirical_fraction",
        "empirical_fit_uncertainty_frac",
        "failed_fit_count",
        "zero_background_fallback_count",
        "fit1_oversub_bin_count",
        "fit2_oversub_bin_count",
        "fit1_max_unclamped_ratio",
        "fit2_max_unclamped_ratio",
    ]
    primary_csv_path = paths["csv_profile"]
    with open(primary_csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})
    written.append(primary_csv_path)
    if active_profile in (None, "", "nominal_weighted") and os.path.abspath(paths["csv"]) != os.path.abspath(primary_csv_path):
        with open(paths["csv"], "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key) for key in fieldnames})
        written.append(paths["csv"])
    return written
