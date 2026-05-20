#! /usr/bin/python

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import shutil
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from frozen_manifest import (
    get_analysis_artifact_paths,
    get_correction_ledger_paths,
    should_write_generic_alias,
)


STAGE_ORDER = [
    "raw_prompt",
    "after_random_subtraction",
    "after_dummy_subtraction",
    "after_pion_subtraction",
    "after_empirical_fit1",
    "after_empirical_fit2",
    "final_lambda_window",
]

STAGE_LABELS = {
    "raw_prompt": "Raw",
    "after_random_subtraction": "After Rand",
    "after_dummy_subtraction": "After Dummy",
    "after_pion_subtraction": "After Pion",
    "after_empirical_fit1": "After Fit 1",
    "after_empirical_fit2": "After Fit 2",
    "final_lambda_window": "Final",
}

BIN_KEY_RE = re.compile(r"t_bin(\d+)phi_bin(\d+)")
PHI_SETTING_COLORS = {
    "Center": "#1f77b4",
    "Left": "#d62728",
    "Right": "#2ca02c",
}


def _infer_run_dir(path):
    path = Path(path).expanduser()
    if path.is_file():
        parent = path.parent
        if parent.name.lower() in {"json", "csv", "plots"}:
            return parent.parent
        return parent
    if path.is_dir() and path.name.lower() in {"json", "csv", "plots"}:
        return path.parent
    return path


def _infer_output_dir_for_manifest(manifest_path):
    return _infer_run_dir(manifest_path) / "plots"


def _resolve_search_root(path, pattern):
    path = Path(path).expanduser()
    if path.is_file():
        if "_frozen_manifest" in path.name and path.suffix.lower() == ".json":
            return [path]
        run_dir = _infer_run_dir(path)
        if (run_dir / "json").is_dir():
            return sorted((run_dir / "json").rglob(pattern))
        return sorted(run_dir.rglob(pattern))

    if not path.is_dir():
        return sorted(Path().glob(str(path)))

    lower_name = path.name.lower()
    if lower_name == "plots" and (path.parent / "json").is_dir():
        return sorted((path.parent / "json").rglob(pattern))
    if lower_name == "json":
        return sorted(path.rglob(pattern))
    if (path / "json").is_dir():
        return sorted((path / "json").rglob(pattern))
    return sorted(path.rglob(pattern))


def _collect_manifest_paths(inputs, pattern):
    manifest_paths = []
    seen = set()
    for raw_input in inputs:
        for candidate in _resolve_search_root(raw_input, pattern):
            if not candidate.is_file():
                continue
            if "_frozen_manifest" not in candidate.name or candidate.suffix.lower() != ".json":
                continue
            resolved = candidate.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            manifest_paths.append(candidate)
    return sorted(manifest_paths)


def _describe_input_behavior(raw_input):
    path = Path(raw_input).expanduser()
    if path.is_file():
        if "_frozen_manifest" in path.name:
            return "{} -> direct manifest".format(path)
        return "{} -> infer run directory and search for manifests".format(path)
    if path.is_dir():
        lower_name = path.name.lower()
        if lower_name == "plots" and (path.parent / "json").is_dir():
            return "{} -> uses sibling json/ and writes back to plots/".format(path)
        if lower_name == "json":
            return "{} -> recursive manifest search, writes to sibling plots/".format(path)
        if (path / "json").is_dir():
            return "{} -> uses embedded json/ and writes to embedded plots/".format(path)
        return "{} -> recursive search, writes to each run's plots/".format(path)
    return "{} -> glob expansion".format(raw_input)


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def _maybe_load_json(path):
    if not path or not os.path.exists(path):
        return None
    try:
        return _load_json(path)
    except Exception as exc:
        print("WARNING: Failed to load JSON {}: {}".format(path, exc))
        return None


def _profile_suffix_score(path, active_profile):
    if not active_profile:
        return 0
    safe_profile = str(active_profile).replace(" ", "_")
    return int(path.stem.endswith("_{}".format(safe_profile)))


def _discover_manifest_runs(manifest_paths, profile_filter=None):
    discovered = {}
    for manifest_path in manifest_paths:
        payload = _maybe_load_json(manifest_path)
        if not isinstance(payload, dict):
            continue
        particle = payload.get("particle_type")
        q2 = str(payload.get("q2", ""))
        w = str(payload.get("w", ""))
        active_profile = payload.get("active_profile") or "nominal_weighted"
        if profile_filter and str(active_profile) != str(profile_filter):
            continue
        run_dir = _infer_run_dir(manifest_path)
        key = (str(run_dir.resolve()), str(particle), q2, w, str(active_profile))
        score = (_profile_suffix_score(manifest_path, active_profile), int(manifest_path.parent.name.lower() == "json"))
        incumbent = discovered.get(key)
        if incumbent is None or score > incumbent["score"]:
            discovered[key] = {
                "run_dir": run_dir,
                "manifest_path": manifest_path,
                "manifest": payload,
                "particle_type": particle,
                "q2": q2,
                "w": w,
                "active_profile": active_profile,
                "score": score,
            }
    return [entry for entry in sorted(discovered.values(), key=lambda item: item["manifest_path"])]


def _find_existing_path(run_dir, paths, subdir=None):
    run_dir = Path(run_dir)
    candidates = []
    seen = set()
    for path in paths:
        if not path:
            continue
        path = Path(path)
        for candidate in [path, run_dir / path.name, run_dir / subdir / path.name if subdir else None]:
            if candidate is None:
                continue
            candidate = Path(candidate)
            candidate_key = str(candidate.resolve()) if candidate.exists() else str(candidate)
            if candidate_key in seen:
                continue
            seen.add(candidate_key)
            candidates.append(candidate)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _find_json_artifact(run_dir, primary_path=None, alias_path=None):
    return _find_existing_path(run_dir, [primary_path, alias_path], subdir="json")


def _find_bg_opt_csvs(run_dir):
    run_dir = Path(run_dir)
    csv_dir = run_dir / "csv"
    if csv_dir.is_dir():
        return sorted(csv_dir.glob("*_bg_opt_*.csv"))
    return sorted(run_dir.glob("*_bg_opt_*.csv"))


def _search_correction_ledgers(run_dir, particle_type, q2, w, active_profile):
    ledgers = {}
    for base_dir in [Path(run_dir) / "json", Path(run_dir)]:
        if not base_dir.is_dir():
            continue
        for candidate in sorted(base_dir.glob("*_correction_ledger*.json")):
            payload = _maybe_load_json(candidate)
            if not isinstance(payload, dict):
                continue
            if str(payload.get("particle_type")) != str(particle_type):
                continue
            if str(payload.get("q2")) != str(q2) or str(payload.get("w")) != str(w):
                continue
            profile_name = payload.get("active_profile") or "nominal_weighted"
            if str(profile_name) != str(active_profile):
                continue
            epsset = str(payload.get("epsset", "")).lower()
            if epsset not in {"low", "high"}:
                continue
            score = _profile_suffix_score(candidate, active_profile)
            incumbent = ledgers.get(epsset)
            if incumbent is None or score >= incumbent["score"]:
                ledgers[epsset] = {"path": candidate, "payload": payload, "score": score}
    return ledgers


def _load_run_artifacts(run_entry):
    run_dir = Path(run_entry["run_dir"])
    particle_type = run_entry["particle_type"]
    q2 = run_entry["q2"]
    w = run_entry["w"]
    active_profile = run_entry["active_profile"]
    artifact_paths = get_analysis_artifact_paths(str(run_dir), particle_type, q2, w, active_profile=active_profile)

    bundle_path = _find_json_artifact(
        run_dir,
        artifact_paths["input_bundle_profile"],
        artifact_paths["input_bundle"],
    )
    input_bundle = _maybe_load_json(bundle_path)

    ledgers = {}
    if isinstance(input_bundle, dict):
        for eps_name in ("low", "high"):
            outfilename = (
                input_bundle.get(eps_name, {})
                .get("inp_dict", {})
                .get("OutFilename")
            )
            if not outfilename:
                continue
            ledger_paths = get_correction_ledger_paths(str(run_dir), particle_type, outfilename, active_profile=active_profile)
            ledger_path = _find_json_artifact(run_dir, ledger_paths["json_profile"], ledger_paths["json"])
            payload = _maybe_load_json(ledger_path)
            if payload is not None:
                ledgers[eps_name] = {"path": ledger_path, "payload": payload}

    if "low" not in ledgers or "high" not in ledgers:
        search_ledgers = _search_correction_ledgers(run_dir, particle_type, q2, w, active_profile)
        for eps_name, entry in search_ledgers.items():
            ledgers.setdefault(eps_name, {"path": entry["path"], "payload": entry["payload"]})

    epsilon_compare_path = _find_json_artifact(
        run_dir,
        artifact_paths["epsilon_compare_json_profile"],
        artifact_paths["epsilon_compare_json"],
    )
    final_summary_path = _find_json_artifact(
        run_dir,
        artifact_paths["final_summary_json_profile"],
        artifact_paths["final_summary_json"],
    )
    nonklambda_path = _find_json_artifact(
        run_dir,
        artifact_paths["nonklambda_json_profile"],
        artifact_paths["nonklambda_json"],
    )
    systematics_path = _find_json_artifact(
        run_dir,
        artifact_paths["systematics_json"],
        artifact_paths["systematics_json"],
    )

    return {
        "manifest_path": run_entry["manifest_path"],
        "manifest": run_entry["manifest"],
        "bundle_path": bundle_path,
        "input_bundle": input_bundle,
        "low_ledger_path": ledgers.get("low", {}).get("path"),
        "low_ledger": ledgers.get("low", {}).get("payload"),
        "high_ledger_path": ledgers.get("high", {}).get("path"),
        "high_ledger": ledgers.get("high", {}).get("payload"),
        "epsilon_compare_path": epsilon_compare_path,
        "epsilon_compare": _maybe_load_json(epsilon_compare_path),
        "final_summary_path": final_summary_path,
        "final_summary": _maybe_load_json(final_summary_path),
        "nonklambda_path": nonklambda_path,
        "nonklambda": _maybe_load_json(nonklambda_path),
        "systematics_path": systematics_path,
        "systematics": _maybe_load_json(systematics_path),
        "bg_opt_csvs": _find_bg_opt_csvs(run_dir),
    }


def _plot_output_paths(run_dir, particle_type, q2, w, active_profile, output_dir=None, suffix=None):
    safe_profile = None if not active_profile else str(active_profile).replace(" ", "_")
    plots_dir = Path(output_dir).expanduser() if output_dir else _infer_output_dir_for_manifest(run_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    generic_name = "{}_Q{}W{}_workflow_diagnostics".format(particle_type, q2, w)
    profile_name = "{}_Q{}W{}_workflow_diagnostics{}".format(
        particle_type,
        q2,
        w,
        "" if not safe_profile else "_{}".format(safe_profile),
    )
    if suffix:
        generic_name = "{}_{}".format(generic_name, suffix)
        profile_name = "{}_{}".format(profile_name, suffix)
    generic_path = plots_dir / "{}.pdf".format(generic_name)
    profile_path = plots_dir / "{}.pdf".format(profile_name)
    if should_write_generic_alias(active_profile):
        return profile_path, generic_path
    return profile_path, None


def _float_or_nan(value):
    try:
        return float(value)
    except Exception:
        return float("nan")


def _short_profile_json(payload):
    if not payload:
        return "{}"
    return json.dumps(payload, indent=2, sort_keys=True)


def _safe_label(text):
    if text is None:
        return ""
    return str(text)


def _wrap_lines(lines, width=108):
    wrapped = []
    for line in lines:
        if not line:
            wrapped.append("")
            continue
        wrapped.extend(textwrap.wrap(str(line), width=width) or [""])
    return wrapped


def _text_page(pdf, title, lines, fontsize=9):
    fig, ax = plt.subplots(figsize=(11.0, 8.5))
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold")
    ax.text(
        0.01,
        0.98,
        "\n".join(_wrap_lines(lines)),
        transform=ax.transAxes,
        va="top",
        ha="left",
        family="monospace",
        fontsize=fontsize,
    )
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _format_value(value, precision=5):
    if value is None:
        return "None"
    try:
        numeric = float(value)
        if math.isnan(numeric):
            return "nan"
        return "{:.{prec}f}".format(numeric, prec=precision)
    except Exception:
        return str(value)


def _phi_color(phi_setting):
    return PHI_SETTING_COLORS.get(str(phi_setting), "#7f7f7f")


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


def _parse_bin_key(bin_key):
    match = BIN_KEY_RE.fullmatch(str(bin_key or ""))
    if not match:
        return None
    return int(match.group(1)) - 1, int(match.group(2)) - 1


def _matrix_from_entries(entries, value_getter):
    parsed = []
    max_t = -1
    max_phi = -1
    for entry in entries:
        indices = _parse_bin_key(entry.get("bin_key"))
        if indices is None:
            continue
        t_index, phi_index = indices
        max_t = max(max_t, t_index)
        max_phi = max(max_phi, phi_index)
        parsed.append((t_index, phi_index, entry))
    if not parsed:
        return None
    matrix = np.full((max_t + 1, max_phi + 1), np.nan)
    for t_index, phi_index, entry in parsed:
        matrix[t_index, phi_index] = _float_or_nan(value_getter(entry))
    return matrix


def _plot_heatmap(ax, matrix, title, colorbar_label, cmap="viridis", center_zero=False):
    if matrix is None:
        ax.axis("off")
        ax.set_title(title)
        ax.text(0.5, 0.5, "No bin-level data", ha="center", va="center", transform=ax.transAxes)
        return
    if center_zero:
        finite_vals = matrix[np.isfinite(matrix)]
        vmax = np.max(np.abs(finite_vals)) if finite_vals.size else 1.0
        vmin = -vmax
    else:
        vmin = None
        vmax = None
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Phi bin")
    ax.set_ylabel("t bin")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels([str(index + 1) for index in range(matrix.shape[1])], fontsize=8)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels([str(index + 1) for index in range(matrix.shape[0])], fontsize=8)
    if matrix.size <= 30:
        for t_index in range(matrix.shape[0]):
            for phi_index in range(matrix.shape[1]):
                value = matrix[t_index, phi_index]
                if np.isfinite(value):
                    ax.text(phi_index, t_index, "{:.3f}".format(float(value)), ha="center", va="center", fontsize=7, color="white")
    plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label=colorbar_label)


def _plot_manifest_pages(pdf, artifacts):
    manifest = artifacts["manifest"]
    profile_settings = manifest.get("resolved_optimizer_settings", {})
    lines = [
        "Manifest: {}".format(artifacts["manifest_path"]),
        "Particle: {}".format(manifest.get("particle_type")),
        "Q2/W: Q{} W{}".format(manifest.get("q2"), manifest.get("w")),
        "Active profile: {}".format(manifest.get("active_profile")),
        "Common epsilon scale behavior: {}".format(manifest.get("common_epsilon_scale_behavior", "independent")),
        "MM cut window: {}".format(manifest.get("mm_cut_window")),
        "Low epsilon: {}".format((manifest.get("epsilon_values") or {}).get("low")),
        "High epsilon: {}".format((manifest.get("epsilon_values") or {}).get("high")),
        "Git commit: {}".format(manifest.get("git_commit_hash")),
        "Created at: {}".format(manifest.get("created_at")),
        "",
        "Resolved optimizer settings:",
        _short_profile_json(profile_settings),
    ]
    _text_page(pdf, "0th Iteration Freeze Contract", lines, fontsize=8)

    t_edges = manifest.get("t_bin_edges") or []
    phi_edges = manifest.get("phi_bin_edges") or []
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5))
    fig.suptitle("Frozen Binning", fontsize=15, fontweight="bold")
    for axis, edges, label in [
        (axes[0, 0], t_edges, "t edges"),
        (axes[0, 1], phi_edges, "phi edges"),
    ]:
        if edges:
            axis.step(range(len(edges)), edges, where="mid", linewidth=2.0)
            axis.scatter(range(len(edges)), edges, s=30)
            axis.set_xticks(range(len(edges)))
        axis.set_title(label)
        axis.set_xlabel("Edge index")
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    for axis, edges, label in [
        (axes[1, 0], t_edges, "t bin widths"),
        (axes[1, 1], phi_edges, "phi bin widths"),
    ]:
        widths = np.diff(np.asarray(edges, dtype=float)) if len(edges) >= 2 else np.asarray([])
        axis.bar(range(1, len(widths) + 1), widths, color="#4c78a8")
        axis.set_title(label)
        axis.set_xlabel("Bin index")
        axis.set_ylabel("Width")
        axis.grid(alpha=0.25, axis="y")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.5), sharex=True)
    fig.suptitle("Selected Empirical Scale Maps", fontsize=15, fontweight="bold")
    for axis, fit_key, title in [
        (axes[0], "selected_bg_scale1s", "Fit 1 scales"),
        (axes[1], "selected_bg_scale2s", "Fit 2 scales"),
    ]:
        items = []
        for eps_name in ("low", "high"):
            scale_map = manifest.get(fit_key, {}).get(eps_name, {}) or {}
            for scale_key, scale_value in sorted(scale_map.items()):
                label = scale_key if ":" in str(scale_key) else "{}:{}".format(eps_name, scale_key)
                items.append((label, _float_or_nan(scale_value), eps_name))
        if items:
            xvals = np.arange(len(items))
            axis.bar(
                xvals,
                [item[1] for item in items],
                color=["#4c78a8" if item[2] == "low" else "#f58518" for item in items],
            )
            axis.set_xticks(xvals)
            axis.set_xticklabels([item[0].replace(":", "\n") for item in items], rotation=0, fontsize=8)
        axis.set_ylabel("Scale")
        axis.set_title(title)
        axis.grid(alpha=0.25, axis="y")
    axes[1].set_xlabel("Epsilon / phi setting")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _plot_input_bundle_page(pdf, artifacts):
    bundle = artifacts.get("input_bundle")
    if not isinstance(bundle, dict):
        return
    lines = [
        "Input bundle: {}".format(artifacts.get("bundle_path")),
        "Bundle active profile: {}".format(bundle.get("active_profile")),
        "",
    ]
    for eps_name in ("low", "high"):
        entry = bundle.get(eps_name) or {}
        inp_dict = entry.get("inp_dict") or {}
        lines.extend(
            [
                "{} epsilon".format(eps_name.capitalize()),
                "  argv: {}".format(" ".join(str(token) for token in entry.get("argv", []))),
                "  OutFilename: {}".format(inp_dict.get("OutFilename")),
                "  EPSSET / EPSVAL: {} / {}".format(inp_dict.get("EPSSET"), inp_dict.get("EPSVAL")),
                "  NumtBins / NumPhiBins: {} / {}".format(inp_dict.get("NumtBins"), inp_dict.get("NumPhiBins")),
                "  ParticleType / POL: {} / {}".format(inp_dict.get("ParticleType"), inp_dict.get("POL")),
                "",
            ]
        )
    _text_page(pdf, "0th Iteration Replay Inputs", lines, fontsize=8)


def _plot_correction_ledger_pages(pdf, ledger, title_prefix):
    if not isinstance(ledger, dict):
        return
    settings = ledger.get("settings", [])
    combined = ledger.get("combined_totals", {})
    stage_labels = [STAGE_LABELS[stage] for stage in STAGE_ORDER]
    combined_stage_yields = [combined.get("combined_stage_yields", {}).get(stage, 0.0) for stage in STAGE_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 8.5))
    fig.suptitle("{}: Stage Yields".format(title_prefix), fontsize=15, fontweight="bold")
    axes[0].plot(stage_labels, combined_stage_yields, marker="o", linewidth=2.0, color="#4c78a8")
    axes[0].set_title("Combined totals")
    axes[0].set_ylabel("Lambda-window yield")
    axes[0].tick_params(axis="x", rotation=35)
    axes[0].grid(alpha=0.25)

    for setting_entry in settings:
        phi_setting = _safe_label(setting_entry.get("phi_setting"))
        stage_yields = [setting_entry.get("combined_stage_yields", {}).get(stage, 0.0) for stage in STAGE_ORDER]
        axes[1].plot(
            stage_labels,
            stage_yields,
            marker="o",
            linewidth=1.8,
            label=phi_setting,
            color=_phi_color(phi_setting),
        )
    axes[1].set_title("Per-setting totals")
    axes[1].set_ylabel("Lambda-window yield")
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].grid(alpha=0.25)
    if settings:
        axes[1].legend(fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    labels = [_safe_label(setting_entry.get("phi_setting")) for setting_entry in settings] + ["Combined"]
    fit1_frac = [setting_entry.get("fit1_fractional_correction", 0.0) for setting_entry in settings] + [combined.get("fit1_fractional_correction", 0.0)]
    fit2_frac = [setting_entry.get("fit2_fractional_correction", 0.0) for setting_entry in settings] + [combined.get("fit2_fractional_correction", 0.0)]
    total_frac = [setting_entry.get("total_empirical_fraction", 0.0) for setting_entry in settings] + [combined.get("total_empirical_fraction", 0.0)]
    failed_counts = [setting_entry.get("failed_fit_count", 0) for setting_entry in settings] + [combined.get("failed_fit_count", 0)]
    zero_bg_counts = [setting_entry.get("zero_background_fallback_count", 0) for setting_entry in settings] + [combined.get("zero_background_fallback_count", 0)]
    fit1_max_ratio = [
        (setting_entry.get("oversub_diagnostics", {}).get("fit1", {}) or {}).get("max_unclamped_ratio", 0.0)
        for setting_entry in settings
    ] + [(combined.get("oversub_diagnostics", {}).get("fit1", {}) or {}).get("max_unclamped_ratio", 0.0)]
    fit2_max_ratio = [
        (setting_entry.get("oversub_diagnostics", {}).get("fit2", {}) or {}).get("max_unclamped_ratio", 0.0)
        for setting_entry in settings
    ] + [(combined.get("oversub_diagnostics", {}).get("fit2", {}) or {}).get("max_unclamped_ratio", 0.0)]
    uncertainty_frac = [setting_entry.get("empirical_fit_uncertainty_frac", 0.0) for setting_entry in settings] + [combined.get("empirical_fit_uncertainty_frac", 0.0)]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5))
    fig.suptitle("{}: Correction Diagnostics".format(title_prefix), fontsize=15, fontweight="bold")
    xvals = np.arange(len(labels))
    width = 0.24
    axes[0, 0].bar(xvals - width, fit1_frac, width=width, label="Fit 1")
    axes[0, 0].bar(xvals, fit2_frac, width=width, label="Fit 2")
    axes[0, 0].bar(xvals + width, total_frac, width=width, label="Total")
    axes[0, 0].set_title("Empirical correction fractions")
    axes[0, 0].set_xticks(xvals)
    axes[0, 0].set_xticklabels(labels, rotation=25)
    axes[0, 0].grid(alpha=0.25, axis="y")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].bar(xvals - width / 2.0, failed_counts, width=width, label="Failed fits", color="#e45756")
    axes[0, 1].bar(xvals + width / 2.0, zero_bg_counts, width=width, label="Zero-bg fallback", color="#72b7b2")
    axes[0, 1].set_title("Failure / fallback counts")
    axes[0, 1].set_xticks(xvals)
    axes[0, 1].set_xticklabels(labels, rotation=25)
    axes[0, 1].grid(alpha=0.25, axis="y")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].bar(xvals - width / 2.0, fit1_max_ratio, width=width, label="Fit 1", color="#4c78a8")
    axes[1, 0].bar(xvals + width / 2.0, fit2_max_ratio, width=width, label="Fit 2", color="#f58518")
    axes[1, 0].axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    axes[1, 0].set_title("Max unclamped MM background ratio")
    axes[1, 0].set_xticks(xvals)
    axes[1, 0].set_xticklabels(labels, rotation=25)
    axes[1, 0].grid(alpha=0.25, axis="y")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].bar(xvals, uncertainty_frac, color="#54a24b")
    axes[1, 1].set_title("Empirical fit uncertainty fraction")
    axes[1, 1].set_xticks(xvals)
    axes[1, 1].set_xticklabels(labels, rotation=25)
    axes[1, 1].grid(alpha=0.25, axis="y")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    n_settings = max(1, len(settings))
    fig, axes = plt.subplots(2, n_settings, figsize=(4.0 * n_settings, 8.5), squeeze=False)
    fig.suptitle("{}: Bin-Level Heatmaps".format(title_prefix), fontsize=15, fontweight="bold")
    for index, setting_entry in enumerate(settings):
        entries = setting_entry.get("t_phi_bins", [])
        total_matrix = _matrix_from_entries(entries, lambda entry: entry.get("total_empirical_fraction"))
        final_matrix = _matrix_from_entries(entries, lambda entry: (entry.get("stage_yields", {}) or {}).get("final_lambda_window"))
        phi_setting = _safe_label(setting_entry.get("phi_setting"))
        _plot_heatmap(
            axes[0, index],
            total_matrix,
            "{} total empirical fraction".format(phi_setting),
            "fraction",
            cmap="coolwarm",
            center_zero=True,
        )
        _plot_heatmap(
            axes[1, index],
            final_matrix,
            "{} final Lambda yield".format(phi_setting),
            "yield",
            cmap="viridis",
            center_zero=False,
        )
    if not settings:
        for axis in axes.flatten():
            axis.axis("off")
            axis.text(0.5, 0.5, "No setting-level correction ledger entries", ha="center", va="center", transform=axis.transAxes)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _plot_epsilon_compare_pages(pdf, compare_payload):
    if not isinstance(compare_payload, dict):
        return
    rows = compare_payload.get("rows", [])
    if not rows:
        _text_page(
            pdf,
            "High/Low Epsilon Compare",
            [
                "No shared low/high bin-level empirical correction rows were found.",
                "Source: {}".format(compare_payload),
            ],
        )
        return

    phi_settings = sorted({row.get("phi_setting") for row in rows})
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 8.5))
    fig.suptitle("High/Low Epsilon Empirical Correction Compare", fontsize=15, fontweight="bold")
    for phi_setting in phi_settings:
        subset = [row for row in rows if row.get("phi_setting") == phi_setting]
        axes[0].scatter(
            [float(row.get("f_emp_low", 0.0) or 0.0) for row in subset],
            [float(row.get("f_emp_high", 0.0) or 0.0) for row in subset],
            label=phi_setting,
            color=_phi_color(phi_setting),
            s=40,
        )
    low_vals = [float(row.get("f_emp_low", 0.0) or 0.0) for row in rows]
    high_vals = [float(row.get("f_emp_high", 0.0) or 0.0) for row in rows]
    lo = min(low_vals + high_vals)
    hi = max(low_vals + high_vals)
    axes[0].plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1.0)
    axes[0].set_title("f_emp low vs high")
    axes[0].set_xlabel("f_emp low")
    axes[0].set_ylabel("f_emp high")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    xvals = np.arange(len(rows))
    axes[1].bar(
        xvals,
        [float(row.get("delta_f_emp", 0.0) or 0.0) for row in rows],
        color=[_phi_color(row.get("phi_setting")) for row in rows],
    )
    warn_threshold = compare_payload.get("warning_threshold")
    fail_threshold = compare_payload.get("fail_threshold")
    if warn_threshold is not None:
        axes[1].axhline(float(warn_threshold), color="#e45756", linestyle="--", linewidth=1.0)
        axes[1].axhline(-float(warn_threshold), color="#e45756", linestyle="--", linewidth=1.0)
    if fail_threshold is not None:
        axes[1].axhline(float(fail_threshold), color="#7f0000", linestyle=":", linewidth=1.0)
        axes[1].axhline(-float(fail_threshold), color="#7f0000", linestyle=":", linewidth=1.0)
    axes[1].set_title("delta_f_emp by shared bin")
    axes[1].set_xlabel("Shared bin index")
    axes[1].set_ylabel("delta_f_emp")
    axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, max(1, len(phi_settings)), figsize=(4.0 * max(1, len(phi_settings)), 4.6), squeeze=False)
    fig.suptitle(
        "High - Low Empirical Fraction Difference Heatmaps\nmean_low={} mean_high={} rms={} max_abs={}".format(
            _format_value(compare_payload.get("mean_f_emp_low"), precision=4),
            _format_value(compare_payload.get("mean_f_emp_high"), precision=4),
            _format_value(compare_payload.get("rms_delta_f_emp"), precision=4),
            _format_value(compare_payload.get("max_abs_delta_f_emp"), precision=4),
        ),
        fontsize=14,
        fontweight="bold",
    )
    for index, phi_setting in enumerate(phi_settings):
        subset = [row for row in rows if row.get("phi_setting") == phi_setting]
        matrix = _matrix_from_entries(subset, lambda row: row.get("delta_f_emp"))
        _plot_heatmap(
            axes[0, index],
            matrix,
            "{} delta_f_emp".format(phi_setting),
            "delta_f_emp",
            cmap="coolwarm",
            center_zero=True,
        )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _plot_final_summary_pages(pdf, summary):
    if not isinstance(summary, dict):
        return
    low_corr = summary.get("low_epsilon_corrections", {}) or {}
    high_corr = summary.get("high_epsilon_corrections", {}) or {}
    xsect_outputs = summary.get("xsect_outputs", {}) or {}
    separated_rows = xsect_outputs.get("separated_csv") or []
    unsep_rows = xsect_outputs.get("unseparated_csv") or []

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 8.5))
    fig.suptitle("Final Analysis Summary", fontsize=15, fontweight="bold")
    metric_labels = ["Fit 1", "Fit 2", "Total"]
    xvals = np.arange(len(metric_labels))
    width = 0.35
    axes[0].bar(
        xvals - width / 2.0,
        [
            low_corr.get("fit1_fractional_correction", 0.0),
            low_corr.get("fit2_fractional_correction", 0.0),
            low_corr.get("total_empirical_fraction", 0.0),
        ],
        width=width,
        label="Low",
        color="#4c78a8",
    )
    axes[0].bar(
        xvals + width / 2.0,
        [
            high_corr.get("fit1_fractional_correction", 0.0),
            high_corr.get("fit2_fractional_correction", 0.0),
            high_corr.get("total_empirical_fraction", 0.0),
        ],
        width=width,
        label="High",
        color="#f58518",
    )
    axes[0].set_title("Low / high empirical correction fractions")
    axes[0].set_xticks(xvals)
    axes[0].set_xticklabels(metric_labels)
    axes[0].grid(alpha=0.25, axis="y")
    axes[0].legend()

    xsect_labels = ["Unsep", "sigma_L", "sigma_T", "sigma_LT", "sigma_TT"]
    xsect_values = [
        _extract_scalar(unsep_rows, ["xsec", "unsep_xsec", "cross_section"]),
        _extract_scalar(separated_rows, ["sigma_L", "sigL", "L"]),
        _extract_scalar(separated_rows, ["sigma_T", "sigT", "T"]),
        _extract_scalar(separated_rows, ["sigma_LT", "sigLT", "LT"]),
        _extract_scalar(separated_rows, ["sigma_TT", "sigTT", "TT"]),
    ]
    finite_xsect = [value for value in xsect_values if isinstance(value, (int, float)) and math.isfinite(float(value))]
    if finite_xsect:
        axes[1].bar(range(len(xsect_labels)), [0.0 if value is None else float(value) for value in xsect_values], color="#54a24b")
        axes[1].set_xticks(range(len(xsect_labels)))
        axes[1].set_xticklabels(xsect_labels, rotation=20)
        axes[1].set_title("Cross-section scalars from final summary")
        axes[1].grid(alpha=0.25, axis="y")
    else:
        axes[1].axis("off")
        axes[1].text(
            0.02,
            0.98,
            "\n".join(
                [
                    "No numeric cross-section scalars were available.",
                    "Active profile: {}".format(summary.get("active_profile")),
                    "Common epsilon scale behavior: {}".format(summary.get("common_epsilon_scale_behavior")),
                    "Manifest hash: {}".format(summary.get("manifest_hash")),
                    "Iteration count: {}".format(summary.get("iteration_count")),
                    "Convergence status: {}".format(summary.get("convergence_status")),
                ]
            ),
            transform=axes[1].transAxes,
            va="top",
            ha="left",
            family="monospace",
            fontsize=9,
        )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    _text_page(
        pdf,
        "Final Summary Metadata",
        [
            "Final summary: {}".format(summary),
            "Active profile: {}".format(summary.get("active_profile")),
            "Common epsilon scale behavior: {}".format(summary.get("common_epsilon_scale_behavior")),
            "Manifest hash: {}".format(summary.get("manifest_hash")),
            "Iteration count: {}".format(summary.get("iteration_count")),
            "Convergence status: {}".format(summary.get("convergence_status")),
            "",
            "Resolved optimizer settings:",
            _short_profile_json(summary.get("resolved_optimizer_settings")),
            "",
            "Fit 1 scales:",
            _short_profile_json(summary.get("fit1_scales")),
            "",
            "Fit 2 scales:",
            _short_profile_json(summary.get("fit2_scales")),
        ],
        fontsize=8,
    )


def _plot_systematics_pages(pdf, systematics):
    if not isinstance(systematics, dict):
        return
    rows = systematics.get("rows", [])
    if not rows:
        return

    labels = [row.get("profile_name") for row in rows]
    lambda_values = [None if row.get("lambda_yield") is None else float(row.get("lambda_yield")) for row in rows]
    deviations = [
        None if row.get("percent_deviation_from_nominal") is None else float(row.get("percent_deviation_from_nominal"))
        for row in rows
    ]
    colors = ["#4c78a8" if row.get("pass_fail_status") == "passed" else "#9d9d9d" for row in rows]

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.5))
    fig.suptitle("Background Systematic Replay Summary", fontsize=15, fontweight="bold")
    axes[0].bar(range(len(rows)), [0.0 if value is None else value for value in lambda_values], color=colors)
    axes[0].set_xticks(range(len(rows)))
    axes[0].set_xticklabels(labels, rotation=25, ha="right")
    axes[0].set_ylabel("Lambda yield")
    axes[0].set_title("Per-profile Lambda yield")
    axes[0].grid(alpha=0.25, axis="y")

    axes[1].bar(range(len(rows)), [0.0 if value is None else value for value in deviations], color=colors)
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[1].set_xticks(range(len(rows)))
    axes[1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1].set_ylabel("% deviation from nominal")
    axes[1].set_title(
        "Deviation envelope: RMS={} envelope={}".format(
            _format_value(systematics.get("lambda_yield_rms"), precision=4),
            _format_value(systematics.get("lambda_yield_envelope"), precision=4),
        )
    )
    axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    sigma_keys = [
        ("sigma_L", "sigma_L"),
        ("sigma_T", "sigma_T"),
        ("sigma_LT", "sigma_LT"),
        ("sigma_TT", "sigma_TT"),
    ]
    finite_sigma = any(
        row.get(key) not in (None, "")
        for row in rows
        for key, _ in sigma_keys
    )
    if finite_sigma:
        fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5))
        fig.suptitle("Systematic Replay Cross-Section Scalars", fontsize=15, fontweight="bold")
        for axis, (key, label) in zip(axes.flatten(), sigma_keys):
            values = [0.0 if row.get(key) in (None, "") else float(row.get(key)) for row in rows]
            axis.bar(range(len(rows)), values, color=colors)
            axis.set_title(label)
            axis.set_xticks(range(len(rows)))
            axis.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
            axis.grid(alpha=0.25, axis="y")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    lines = [
        "Systematics summary: {}".format(systematics.get("generated_at")),
        "Nominal profile: {}".format(systematics.get("nominal_profile")),
        "Preflight: {}".format(_short_profile_json(systematics.get("preflight"))),
        "",
    ]
    for row in rows:
        lines.append(
            "{} -> {} | deviation={} | error={}".format(
                row.get("profile_name"),
                row.get("pass_fail_status"),
                _format_value(row.get("percent_deviation_from_nominal"), precision=4),
                row.get("error_message"),
            )
        )
    _text_page(pdf, "Systematic Replay Status", lines, fontsize=8)


def _plot_nonklambda_pages(pdf, nonklambda):
    if not isinstance(nonklambda, dict):
        return
    rows = nonklambda.get("rows", [])
    if not rows:
        _text_page(
            pdf,
            "Non-KLambda Cross-Check",
            [
                "No non-KLambda rows were available.",
                "Warnings: {}".format(nonklambda.get("warnings")),
            ],
        )
        return

    statuses = {}
    for row in rows:
        statuses[row.get("status", "unknown")] = statuses.get(row.get("status", "unknown"), 0) + 1

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 8.5))
    fig.suptitle("Non-KLambda Cross-Check", fontsize=15, fontweight="bold")
    axes[0].bar(range(len(statuses)), list(statuses.values()), color="#4c78a8")
    axes[0].set_xticks(range(len(statuses)))
    axes[0].set_xticklabels(list(statuses.keys()), rotation=20)
    axes[0].set_title("Status counts")
    axes[0].grid(alpha=0.25, axis="y")

    comparable_rows = [row for row in rows if row.get("lambda_window_correction_difference") is not None]
    yvals = [float(row.get("lambda_window_correction_difference", 0.0) or 0.0) for row in comparable_rows]
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    if yvals:
        axes[1].plot(range(len(yvals)), yvals, marker="o", linestyle="-", color="#f58518")
    else:
        axes[1].text(0.5, 0.5, "No comparable empirical bins", ha="center", va="center", transform=axes[1].transAxes)
    axes[1].set_title("Lambda-window correction difference")
    axes[1].set_xlabel("Comparable row index")
    axes[1].set_ylabel("Empirical - prediction")
    axes[1].grid(alpha=0.25, axis="y")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    channel_max = {}
    for row in comparable_rows:
        channel = row.get("channel_name")
        channel_max[channel] = max(
            channel_max.get(channel, 0.0),
            abs(float(row.get("lambda_window_correction_difference", 0.0) or 0.0)),
        )
    lines = [
        "Comparison mode: {}".format(nonklambda.get("comparison_mode")),
        "Active profile: {}".format(nonklambda.get("active_profile")),
        "Prediction count: {}".format(nonklambda.get("prediction_count")),
        "Missing empirical bin count: {}".format(nonklambda.get("missing_empirical_bin_count")),
        "Max abs Lambda-window correction difference: {}".format(
            _format_value(nonklambda.get("max_abs_lambda_window_correction_difference"), precision=4)
        ),
        "",
        "Channel max abs differences:",
    ]
    for channel, value in sorted(channel_max.items()):
        lines.append("  {} -> {}".format(channel, _format_value(value, precision=4)))
    if nonklambda.get("warnings"):
        lines.extend(["", "Warnings:"] + ["  {}".format(warning) for warning in nonklambda.get("warnings", [])])
    _text_page(pdf, "Non-KLambda Cross-Check Summary", lines, fontsize=8)


def _plot_cover_page(pdf, run_entry, artifacts):
    manifest = artifacts["manifest"]
    run_dir = Path(run_entry["run_dir"])
    available = [
        ("Manifest", artifacts.get("manifest_path")),
        ("Input bundle", artifacts.get("bundle_path")),
        ("Low correction ledger", artifacts.get("low_ledger_path")),
        ("High correction ledger", artifacts.get("high_ledger_path")),
        ("Epsilon compare", artifacts.get("epsilon_compare_path")),
        ("Final summary", artifacts.get("final_summary_path")),
        ("Systematics summary", artifacts.get("systematics_path")),
        ("Non-KLambda cross-check", artifacts.get("nonklambda_path")),
    ]
    lines = [
        "Workflow diagnostics plot regeneration",
        "",
        "Run directory: {}".format(run_dir),
        "Particle: {}".format(manifest.get("particle_type")),
        "Q2/W: Q{} W{}".format(manifest.get("q2"), manifest.get("w")),
        "Active profile: {}".format(manifest.get("active_profile")),
        "Manifest source: {}".format(artifacts.get("manifest_path")),
        "",
        "Available workflow artifacts:",
    ]
    for label, path in available:
        lines.append("  {}: {}".format(label, path if path else "missing"))
    lines.extend(
        [
            "",
            "BG optimization CSVs found: {}".format(len(artifacts.get("bg_opt_csvs", []))),
            "Use regenerate_bg_opt_plots.py for detailed Step-4 candidate pages.",
        ]
    )
    _text_page(pdf, "Production Workflow Diagnostic Plots", lines, fontsize=8)


def _generate_workflow_pdf(run_entry, artifacts, output_dir=None, suffix=None):
    primary_path, alias_path = _plot_output_paths(
        run_entry["manifest_path"],
        run_entry["particle_type"],
        run_entry["q2"],
        run_entry["w"],
        run_entry["active_profile"],
        output_dir=output_dir,
        suffix=suffix,
    )
    with PdfPages(primary_path) as pdf:
        _plot_cover_page(pdf, run_entry, artifacts)
        _plot_input_bundle_page(pdf, artifacts)
        _plot_manifest_pages(pdf, artifacts)
        if artifacts.get("low_ledger"):
            _plot_correction_ledger_pages(pdf, artifacts["low_ledger"], "Low epsilon correction ledger")
        if artifacts.get("high_ledger"):
            _plot_correction_ledger_pages(pdf, artifacts["high_ledger"], "High epsilon correction ledger")
        if artifacts.get("epsilon_compare"):
            _plot_epsilon_compare_pages(pdf, artifacts["epsilon_compare"])
        if artifacts.get("systematics"):
            _plot_systematics_pages(pdf, artifacts["systematics"])
        if artifacts.get("nonklambda"):
            _plot_nonklambda_pages(pdf, artifacts["nonklambda"])
        if artifacts.get("final_summary"):
            _plot_final_summary_pages(pdf, artifacts["final_summary"])
    written = [primary_path]
    if alias_path and os.path.abspath(alias_path) != os.path.abspath(primary_path):
        shutil.copy(primary_path, alias_path)
        written.append(alias_path)
    return written


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate workflow-stage diagnostic PDFs from cached manifest/ledger/systematics JSON artifacts."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Path(s) to an archived run directory, a json/ directory, a plots/ directory, "
            "a single manifest JSON file, or a glob."
        ),
    )
    parser.add_argument(
        "--pattern",
        default="*_frozen_manifest*.json",
        help="Filename pattern to use when an input is a directory. Default: *_frozen_manifest*.json",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Optional active-profile filter when scanning directories.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional directory to place regenerated PDFs. Default: infer sibling plots/ from each run.",
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help="Optional suffix to append to regenerated PDF filenames before .pdf.",
    )
    args = parser.parse_args(argv)

    manifest_paths = _collect_manifest_paths(args.inputs, args.pattern)
    if not manifest_paths:
        print("No matching frozen manifest JSON files found.")
        return 1

    run_entries = _discover_manifest_runs(manifest_paths, profile_filter=args.profile)
    if not run_entries:
        print("No manifest groups matched the requested profile filter.")
        return 1

    print("Found {} manifest-backed workflow run(s).".format(len(run_entries)))
    for raw_input in args.inputs:
        print("  {}".format(_describe_input_behavior(raw_input)))

    failures = []
    for run_entry in run_entries:
        try:
            artifacts = _load_run_artifacts(run_entry)
            written_paths = _generate_workflow_pdf(run_entry, artifacts, output_dir=args.output_dir, suffix=args.suffix)
            for path in written_paths:
                print("Wrote {}".format(path))
        except Exception as exc:
            failures.append((run_entry, exc))
            print(
                "FAILED {} Q{}W{} {} -> {}".format(
                    run_entry["particle_type"],
                    run_entry["q2"],
                    run_entry["w"],
                    run_entry["active_profile"],
                    exc,
                )
            )

    if failures:
        print("\n{} workflow plot regeneration job(s) failed:".format(len(failures)))
        for run_entry, exc in failures:
            print(
                "  {} Q{}W{} {}: {}".format(
                    run_entry["particle_type"],
                    run_entry["q2"],
                    run_entry["w"],
                    run_entry["active_profile"],
                    exc,
                )
            )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
