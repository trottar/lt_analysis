#! /usr/bin/python

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from copy import deepcopy
from datetime import datetime, timezone

import numpy as np

from background_config import ALLOW_CONFIG_DRIFT, BG_MODELS


KEY_FILE_HASH_TARGETS = [
    "src/utility/background_config.py",
    "src/utility/bg_optimization.py",
    "src/main.py",
    "src/main_iter.py",
    "src/main_auto.py",
    "src/cuts/background_fit.py",
    "src/cuts/rand_sub.py",
    "src/binning/calculate_yield.py",
]


def get_iteration_manifest_hash_policy(current_script_relpath=None):
    """
    Later iterations should fail on live background-config drift and only warn
    on code-file drift.

    The frozen manifest already captures the 0th-iteration products. Once those
    products are materialized and validated, a later SIMC-weight iteration
    should not be blocked by bookkeeping or driver-code updates in
    ``main_iter.py`` / ``main_auto.py``. Keep the signature so existing callers
    do not need to change, but the only strict path is the shared analysis
    config.
    """
    strict_hash_paths = {"src/utility/background_config.py"}
    warn_only_hash_paths = set(KEY_FILE_HASH_TARGETS) - strict_hash_paths
    return strict_hash_paths, warn_only_hash_paths


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(val) for val in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_git_commit_hash(repo_root):
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=repo_root,
                stderr=subprocess.DEVNULL,
                text=True,
            )
            .strip()
        )
    except Exception:
        return None


def build_key_file_hashes(repo_root):
    hashes = {}
    for rel_path in KEY_FILE_HASH_TARGETS:
        abs_path = os.path.join(repo_root, rel_path)
        hashes[rel_path] = {
            "exists": os.path.exists(abs_path),
            "sha256": _sha256_file(abs_path) if os.path.exists(abs_path) else None,
        }
    return hashes


def get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=None):
    stem = os.path.join(outpath, "{}_Q{}W{}".format(particle_type, q2, w))
    safe_profile = None
    if active_profile:
        safe_profile = str(active_profile).replace(" ", "_")
    profile_suffix = "" if not safe_profile else "_{}".format(safe_profile)
    return {
        "manifest": stem + "_frozen_manifest.json",
        "manifest_profile": stem + "_frozen_manifest{}.json".format(profile_suffix),
        "input_bundle": stem + "_0th_iteration_input_bundle.json",
        "input_bundle_profile": stem + "_0th_iteration_input_bundle{}.json".format(profile_suffix),
        "epsilon_compare_json": stem + "_epsilon_empirical_compare.json",
        "epsilon_compare_csv": stem + "_epsilon_empirical_compare.csv",
        "epsilon_compare_json_profile": stem + "_epsilon_empirical_compare{}.json".format(profile_suffix),
        "epsilon_compare_csv_profile": stem + "_epsilon_empirical_compare{}.csv".format(profile_suffix),
        "final_summary_json": stem + "_final_analysis_summary.json",
        "final_summary_csv": stem + "_final_analysis_summary.csv",
        "final_summary_md": stem + "_final_analysis_summary.md",
        "final_summary_json_profile": stem + "_final_analysis_summary{}.json".format(profile_suffix),
        "final_summary_csv_profile": stem + "_final_analysis_summary{}.csv".format(profile_suffix),
        "final_summary_md_profile": stem + "_final_analysis_summary{}.md".format(profile_suffix),
        "systematics_json": stem + "_bg_systematics_summary.json",
        "systematics_csv": stem + "_bg_systematics_summary.csv",
        "nonklambda_json": stem + "_nonklambda_crosscheck.json",
        "nonklambda_csv": stem + "_nonklambda_crosscheck.csv",
        "nonklambda_pdf": stem + "_nonklambda_crosscheck.pdf",
        "nonklambda_json_profile": stem + "_nonklambda_crosscheck{}.json".format(profile_suffix),
        "nonklambda_csv_profile": stem + "_nonklambda_crosscheck{}.csv".format(profile_suffix),
        "nonklambda_pdf_profile": stem + "_nonklambda_crosscheck{}.pdf".format(profile_suffix),
    }


def get_correction_ledger_paths(outpath, particle_type, outfilename, active_profile=None):
    stem = os.path.join(outpath, "{}_{}".format(particle_type, outfilename))
    safe_profile = None
    if active_profile:
        safe_profile = str(active_profile).replace(" ", "_")
    profile_suffix = "" if not safe_profile else "_{}".format(safe_profile)
    return {
        "json": stem + "_correction_ledger.json",
        "csv": stem + "_correction_ledger.csv",
        "json_profile": stem + "_correction_ledger{}.json".format(profile_suffix),
        "csv_profile": stem + "_correction_ledger{}.csv".format(profile_suffix),
    }


def _read_json_if_exists(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r") as handle:
        return json.load(handle)


def should_write_generic_alias(active_profile):
    return active_profile in (None, "", "nominal_weighted")


def resolve_profile_output_paths(generic_path, profile_path=None, active_profile=None):
    primary_path = profile_path or generic_path
    alias_path = None
    if (
        generic_path
        and should_write_generic_alias(active_profile)
        and os.path.abspath(generic_path) != os.path.abspath(primary_path)
    ):
        alias_path = generic_path
    return primary_path, alias_path


def write_json_with_aliases(payload, generic_path, profile_path=None, active_profile=None):
    primary_path, alias_path = resolve_profile_output_paths(
        generic_path,
        profile_path=profile_path,
        active_profile=active_profile,
    )
    serializable = _json_ready(payload)
    with open(primary_path, "w") as handle:
        json.dump(serializable, handle, indent=2, sort_keys=True)
    written = [primary_path]
    if alias_path and os.path.abspath(alias_path) != os.path.abspath(primary_path):
        with open(alias_path, "w") as handle:
            json.dump(serializable, handle, indent=2, sort_keys=True)
        written.append(alias_path)
    return written


def write_text_with_aliases(text, generic_path, profile_path=None, active_profile=None):
    primary_path, alias_path = resolve_profile_output_paths(
        generic_path,
        profile_path=profile_path,
        active_profile=active_profile,
    )
    with open(primary_path, "w") as handle:
        handle.write(text)
    written = [primary_path]
    if alias_path and os.path.abspath(alias_path) != os.path.abspath(primary_path):
        with open(alias_path, "w") as handle:
            handle.write(text)
        written.append(alias_path)
    return written


def _extract_eps_entry(inp_dict):
    serializable = deepcopy(inp_dict)
    return _json_ready(serializable)


def write_zeroth_iteration_input_bundle(
    outpath,
    particle_type,
    q2,
    w,
    epsset,
    argv,
    inp_dict,
    active_profile=None,
):
    artifact_paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    bundle = (
        _read_json_if_exists(artifact_paths["input_bundle_profile"])
        or _read_json_if_exists(artifact_paths["input_bundle"])
        or {
        "particle_type": particle_type,
        "q2": q2,
        "w": w,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    )
    bundle["active_profile"] = active_profile
    bundle[str(epsset).lower()] = {
        "argv": list(argv),
        "inp_dict": _extract_eps_entry(inp_dict),
        "recorded_at": datetime.now(timezone.utc).isoformat(),
    }
    return write_json_with_aliases(
        bundle,
        artifact_paths["input_bundle"],
        artifact_paths["input_bundle_profile"],
        active_profile=active_profile,
    )


def load_zeroth_iteration_input_bundle(outpath, particle_type, q2, w, active_profile=None):
    artifact_paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    preferred_paths = [artifact_paths["input_bundle_profile"], artifact_paths["input_bundle"]]
    for candidate in preferred_paths:
        payload = _read_json_if_exists(candidate)
        if payload is not None:
            return payload
    return None


def _infer_fit_metadata(q2, w, phi_settings, eps_values):
    by_fit = {}
    q2w_key = "Q{}W{}".format(q2, w)
    q2w_models = BG_MODELS.get(q2w_key, {})
    for fit_name, prefix in (("Fit 1", "fixquad"), ("Fit 2", "cheb2")):
        fit_entries = {}
        for epsset, eps_token in eps_values.items():
            for phi_setting in phi_settings:
                model_key = "{}_{}_{}e".format(prefix, phi_setting, eps_token)
                model_config = (
                    q2w_models.get(fit_name, {}).get(model_key)
                    if isinstance(q2w_models.get(fit_name), dict)
                    else None
                )
                fit_entries["{}:{}".format(epsset, phi_setting)] = {
                    "model_key": model_key,
                    "func_expr": None if model_config is None else model_config.get("func_expr"),
                    "n_par": None if model_config is None else model_config.get("n_par"),
                    "sidebands": None if model_config is None else deepcopy(model_config.get("sidebands")),
                }
        by_fit[fit_name] = fit_entries
    return by_fit


def build_frozen_manifest_payload(
    particle_type,
    q2,
    w,
    low_eps,
    high_eps,
    histlist,
    low_bg_summary,
    high_bg_summary,
    low_bundle_entry,
    high_bundle_entry,
    current_inp_dict,
    repo_root,
):
    phi_settings = [hist.get("phi_setting", "") for hist in histlist]
    t_bins = list(np.asarray(histlist[0]["t_bins"], dtype=float)) if histlist else []
    phi_bins = list(np.asarray(histlist[0]["phi_bins"], dtype=float)) if histlist else []
    active_profile = (
        current_inp_dict.get("bg_active_profile")
        or (high_bg_summary or {}).get("active_profile")
        or (low_bg_summary or {}).get("active_profile")
    )
    resolved_profile = (
        current_inp_dict.get("bg_resolved_profile")
        or (high_bg_summary or {}).get("resolved_profile")
        or (low_bg_summary or {}).get("resolved_profile")
        or {}
    )
    eps_values = {
        "low": "low",
        "high": "high",
    }
    payload = {
        "particle_type": particle_type,
        "polarity": current_inp_dict.get("POL"),
        "q2": q2,
        "w": w,
        "epsilon_values": {
            "low": low_eps,
            "high": high_eps,
        },
        "mm_cut_window": [
            float(current_inp_dict.get("mm_min")),
            float(current_inp_dict.get("mm_max")),
        ],
        "t_bin_edges": t_bins,
        "phi_bin_edges": phi_bins,
        "background_models": _infer_fit_metadata(q2, w, phi_settings, eps_values),
        "selected_bg_scale1s": {
            "low": (low_bg_summary or {}).get("selected_bg_scale1s", {}),
            "high": (high_bg_summary or {}).get("selected_bg_scale1s", {}),
        },
        "selected_bg_scale2s": {
            "low": (low_bg_summary or {}).get("selected_bg_scale2s", {}),
            "high": (high_bg_summary or {}).get("selected_bg_scale2s", {}),
        },
        "active_profile": active_profile,
        "resolved_optimizer_settings": _json_ready(resolved_profile),
        "optimizer_selection_mode": (high_bg_summary or {}).get("selection_mode"),
        "optimizer_metric_weights": deepcopy(resolved_profile.get("metric_weights", {})),
        "use_common_epsilon_scales": bool(
            current_inp_dict.get("bg_use_common_epsilon_scales")
            or resolved_profile.get("use_common_epsilon_scales", False)
        ),
        "common_epsilon_scale_behavior": (
            current_inp_dict.get("bg_common_epsilon_scale_behavior")
            or resolved_profile.get("common_epsilon_scale_behavior")
            or "independent"
        ),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git_commit_hash": get_git_commit_hash(repo_root),
        "key_file_hashes": build_key_file_hashes(repo_root),
        "zeroth_iteration_inputs": {
            "low": low_bundle_entry,
            "high": high_bundle_entry,
        },
        "bg_optimization_summaries": {
            "low": _json_ready(low_bg_summary),
            "high": _json_ready(high_bg_summary),
        },
    }
    return payload


def write_frozen_manifest(outpath, particle_type, q2, w, payload, active_profile=None):
    artifact_paths = get_analysis_artifact_paths(outpath, particle_type, q2, w, active_profile=active_profile)
    return write_json_with_aliases(
        payload,
        artifact_paths["manifest"],
        artifact_paths["manifest_profile"],
        active_profile=active_profile,
    )


def read_interval_file(path):
    with open(path, "r") as handle:
        lines = handle.readlines()
    if len(lines) < 2:
        raise ValueError("Interval file {} does not contain a second line".format(path))
    raw_values = [value for value in lines[1].strip().split("\t") if value]
    if raw_values and raw_values[0].lower().startswith(("t", "phi")):
        raw_values = raw_values[1:]
    return np.asarray([float(value) for value in raw_values], dtype=float)


def _array_matches(lhs, rhs, atol=1.0e-9):
    try:
        return np.allclose(np.asarray(lhs, dtype=float), np.asarray(rhs, dtype=float), atol=atol, rtol=0.0)
    except Exception:
        return False


def find_frozen_manifest_path(base_dir, particle_type, q2, w):
    candidates = [
        os.path.join(base_dir, "{}_Q{}W{}_frozen_manifest.json".format(particle_type, q2, w)),
        os.path.join(base_dir, "json", "{}_Q{}W{}_frozen_manifest.json".format(particle_type, q2, w)),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        "Frozen manifest not found in {} for {} Q{}W{}".format(
            base_dir,
            particle_type,
            q2,
            w,
        )
    )


def load_frozen_manifest(manifest_path):
    with open(manifest_path, "r") as handle:
        return json.load(handle)


def validate_iteration_inputs_against_manifest(
    manifest,
    current_inp_dict,
    t_bins,
    phi_bins,
    repo_root,
    required_paths=None,
    strict_hash_paths=None,
    warn_only_hash_paths=None,
):
    required_paths = required_paths or []
    for path in required_paths:
        if not os.path.exists(path):
            raise FileNotFoundError("Required cached artifact is missing: {}".format(path))

    expected_q2 = str(manifest.get("q2"))
    expected_w = str(manifest.get("w"))
    if str(current_inp_dict.get("Q2")) != expected_q2 or str(current_inp_dict.get("W")) != expected_w:
        raise ValueError(
            "Manifest Q2/W mismatch: expected Q{}W{}, got Q{}W{}".format(
                expected_q2,
                expected_w,
                current_inp_dict.get("Q2"),
                current_inp_dict.get("W"),
            )
        )

    epsset = str(current_inp_dict.get("EPSSET", "")).lower()
    manifest_eps = manifest.get("epsilon_values", {}).get(epsset)
    current_eps = current_inp_dict.get("EPSVAL")
    if manifest_eps is not None and current_eps is not None and str(manifest_eps) != str(current_eps):
        raise ValueError(
            "Manifest epsilon mismatch for {}: expected {}, got {}".format(
                epsset,
                manifest_eps,
                current_eps,
            )
        )

    if not _array_matches(manifest.get("t_bin_edges", []), t_bins):
        raise ValueError("t-bin edges do not match frozen manifest")
    if not _array_matches(manifest.get("phi_bin_edges", []), phi_bins):
        raise ValueError("phi-bin edges do not match frozen manifest")
    if int(current_inp_dict.get("NumtBins", len(t_bins) - 1)) != len(t_bins) - 1:
        raise ValueError("Requested NumtBins does not match frozen manifest bin count")
    if int(current_inp_dict.get("NumPhiBins", len(phi_bins) - 1)) != len(phi_bins) - 1:
        raise ValueError("Requested NumPhiBins does not match frozen manifest bin count")

    manifest_hashes = manifest.get("key_file_hashes", {})
    current_hashes = build_key_file_hashes(repo_root)
    current_inp_dict["manifest_validation_warnings"] = []
    if strict_hash_paths is not None:
        strict_hash_paths = set(strict_hash_paths)
    if warn_only_hash_paths is not None:
        warn_only_hash_paths = set(warn_only_hash_paths)
    for rel_path, manifest_entry in manifest_hashes.items():
        manifest_hash = None if not isinstance(manifest_entry, dict) else manifest_entry.get("sha256")
        current_hash = current_hashes.get(rel_path, {}).get("sha256")
        if not manifest_hash:
            continue
        if manifest_hash == current_hash:
            continue
        warning = (
            "Config/code hash drift detected for {} "
            "(expected {}, got {})".format(rel_path, manifest_hash, current_hash or "missing")
        )
        should_fail = not ALLOW_CONFIG_DRIFT
        if strict_hash_paths is not None:
            should_fail = rel_path in strict_hash_paths and rel_path not in (warn_only_hash_paths or set())
        if should_fail:
            raise ValueError(warning)
        print("WARNING: {}".format(warning))
        current_inp_dict["manifest_validation_warnings"].append(warning)

    return {
        "manifest_hashes": current_hashes,
        "warnings": list(current_inp_dict.get("manifest_validation_warnings", [])),
    }
