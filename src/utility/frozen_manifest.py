#! /usr/bin/python

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from glob import glob
from copy import deepcopy
from datetime import datetime, timezone

import numpy as np

from background_config import (
    ALLOW_CONFIG_DRIFT,
    BG_MODELS,
    get_particle_subtraction_window_config,
)


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
        "particle_subtraction_windows": get_particle_subtraction_window_config(
            particle_type,
            "pion",
            inp_dict=current_inp_dict,
        ),
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
    wildcard_candidates = []
    for search_dir in (base_dir, os.path.join(base_dir, "json")):
        wildcard_candidates.extend(
            sorted(
                glob(
                    os.path.join(
                        search_dir,
                        "{}_Q{}W{}_frozen_manifest*.json".format(particle_type, q2, w),
                    )
                )
            )
        )
    for candidate in wildcard_candidates:
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


def _glob_first(patterns):
    for pattern in patterns:
        matches = sorted(glob(pattern))
        for match in matches:
            if os.path.exists(match):
                return match
    return None


def _find_cached_analysis_json(base_dir, particle_type, q2, w, eps_suffix):
    patterns = [
        os.path.join(base_dir, "json", "{}_*Q{}W{}_{}.json".format(particle_type, q2, w, eps_suffix)),
        os.path.join(base_dir, "{}_*Q{}W{}_{}.json".format(particle_type, q2, w, eps_suffix)),
    ]
    excluded_tokens = (
        "_bg_opt_",
        "_correction_ledger",
        "_final_analysis_summary",
        "_frozen_manifest",
        "_0th_iteration_input_bundle",
        "_epsilon_empirical_compare",
        "_bg_systematics_summary",
        "_nonklambda_crosscheck",
    )
    for pattern in patterns:
        for candidate in sorted(glob(pattern)):
            basename = os.path.basename(candidate)
            if any(token in basename for token in excluded_tokens):
                continue
            if os.path.exists(candidate):
                return candidate
    return None


def _infer_outfilename_from_json_path(json_path, particle_type):
    basename = os.path.basename(json_path)
    prefix = "{}_".format(particle_type)
    if basename.startswith(prefix):
        basename = basename[len(prefix):]
    if basename.endswith(".json"):
        basename = basename[:-5]
    return basename


def _read_interval_file_from_cache(base_dir, q2, w, prefix):
    rel_name = "{}_bin_interval_Q{}W{}".format(prefix, q2.replace("p", ""), w.replace("p", ""))
    candidates = [
        os.path.join(base_dir, rel_name),
        os.path.join(base_dir, "root", rel_name),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return read_interval_file(candidate)
    return None


def _load_cached_input_bundle(base_dir, particle_type, q2, w):
    bundle_path = _glob_first(
        [
            os.path.join(base_dir, "json", "{}_Q{}W{}_0th_iteration_input_bundle*.json".format(particle_type, q2, w)),
            os.path.join(base_dir, "{}_Q{}W{}_0th_iteration_input_bundle*.json".format(particle_type, q2, w)),
        ]
    )
    return (_read_json_if_exists(bundle_path), bundle_path) if bundle_path else (None, None)


def reconstruct_frozen_manifest_from_cache(base_dir, particle_type, q2, w):
    low_json = _find_cached_analysis_json(base_dir, particle_type, q2, w, "lowe")
    high_json = _find_cached_analysis_json(base_dir, particle_type, q2, w, "highe")
    if not low_json and not high_json:
        raise FileNotFoundError(
            "Could not reconstruct frozen manifest: no cached analysis JSON files were found in {}".format(base_dir)
        )

    low_payload = _read_json_if_exists(low_json) if low_json else None
    high_payload = _read_json_if_exists(high_json) if high_json else None
    primary_payload = high_payload or low_payload
    if primary_payload is None:
        raise FileNotFoundError(
            "Could not reconstruct frozen manifest: cached analysis JSON payloads were unreadable in {}".format(base_dir)
        )

    low_inp = (low_payload or {}).get("inpDict", {})
    high_inp = (high_payload or {}).get("inpDict", {})
    primary_inp = primary_payload.get("inpDict", {})
    primary_histlist = primary_payload.get("histlist", [])

    t_bins = None
    phi_bins = None
    if primary_histlist:
        t_bins = primary_histlist[0].get("t_bins")
        phi_bins = primary_histlist[0].get("phi_bins")
    if t_bins is None:
        t_bins = _read_interval_file_from_cache(base_dir, q2, w, "t")
    if phi_bins is None:
        phi_bins = _read_interval_file_from_cache(base_dir, q2, w, "phi")
    if t_bins is None or phi_bins is None:
        raise FileNotFoundError(
            "Could not reconstruct frozen manifest: cached t/phi interval files were not found in {}".format(base_dir)
        )

    input_bundle, input_bundle_path = _load_cached_input_bundle(base_dir, particle_type, q2, w)
    active_profile = (
        (input_bundle or {}).get("active_profile")
        or high_inp.get("bg_active_profile")
        or low_inp.get("bg_active_profile")
        or primary_inp.get("bg_active_profile")
    )

    synthesized_bundle = input_bundle or {
        "particle_type": particle_type,
        "q2": q2,
        "w": w,
        "active_profile": active_profile,
    }
    synthesized_bundle.setdefault(
        "low",
        {"inp_dict": {"OutFilename": _infer_outfilename_from_json_path(low_json, particle_type)}} if low_json else {},
    )
    synthesized_bundle.setdefault(
        "high",
        {"inp_dict": {"OutFilename": _infer_outfilename_from_json_path(high_json, particle_type)}} if high_json else {},
    )
    if low_json and "inp_dict" not in synthesized_bundle.get("low", {}):
        synthesized_bundle["low"]["inp_dict"] = {"OutFilename": _infer_outfilename_from_json_path(low_json, particle_type)}
    if high_json and "inp_dict" not in synthesized_bundle.get("high", {}):
        synthesized_bundle["high"]["inp_dict"] = {"OutFilename": _infer_outfilename_from_json_path(high_json, particle_type)}

    manifest_payload = {
        "particle_type": particle_type,
        "polarity": primary_inp.get("POL"),
        "q2": q2,
        "w": w,
        "epsilon_values": {
            "low": low_inp.get("EPSVAL", primary_inp.get("LOEPS")),
            "high": high_inp.get("EPSVAL", primary_inp.get("HIEPS")),
        },
        "mm_cut_window": [
            float(primary_inp.get("mm_min")) if primary_inp.get("mm_min") is not None else None,
            float(primary_inp.get("mm_max")) if primary_inp.get("mm_max") is not None else None,
        ],
        "particle_subtraction_windows": get_particle_subtraction_window_config(
            particle_type,
            "pion",
            inp_dict=primary_inp,
        ),
        "t_bin_edges": _json_ready(np.asarray(t_bins, dtype=float)),
        "phi_bin_edges": _json_ready(np.asarray(phi_bins, dtype=float)),
        "active_profile": active_profile,
        "zeroth_iteration_inputs": synthesized_bundle,
        "key_file_hashes": {},
        "created_at": datetime.now(timezone.utc).isoformat(),
        "reconstructed_from_cache": True,
        "reconstructed_from": {
            "base_dir": base_dir,
            "low_json": low_json,
            "high_json": high_json,
            "input_bundle": input_bundle_path,
        },
    }
    return manifest_payload


def load_or_reconstruct_frozen_manifest(base_dir, particle_type, q2, w):
    try:
        manifest_path = find_frozen_manifest_path(base_dir, particle_type, q2, w)
        return manifest_path, load_frozen_manifest(manifest_path), False
    except FileNotFoundError:
        manifest_payload = reconstruct_frozen_manifest_from_cache(base_dir, particle_type, q2, w)
        manifest_dir = os.path.join(base_dir, "json")
        if not os.path.isdir(manifest_dir):
            os.makedirs(manifest_dir, exist_ok=True)
        manifest_path = os.path.join(
            manifest_dir,
            "{}_Q{}W{}_frozen_manifest.json".format(particle_type, q2, w),
        )
        with open(manifest_path, "w") as handle:
            json.dump(_json_ready(manifest_payload), handle, indent=2, sort_keys=True)
        print("[FROZEN MANIFEST] Reconstructed missing manifest from cached artifacts: {}".format(manifest_path))
        return manifest_path, manifest_payload, True


def _build_hash_drift_message(
    rel_path,
    manifest_hash,
    current_hash,
    manifest=None,
    manifest_path=None,
    allow_config_drift=None,
):
    manifest = manifest or {}
    if allow_config_drift is None:
        allow_config_drift = ALLOW_CONFIG_DRIFT

    lines = [
        "Frozen-manifest validation failed: config/code hash drift detected.",
        "File: {}".format(rel_path),
        "Expected SHA256: {}".format(manifest_hash),
        "Current  SHA256: {}".format(current_hash or "missing"),
    ]
    if manifest_path:
        lines.append("Frozen manifest: {}".format(manifest_path))
    if manifest.get("particle_type") or manifest.get("q2") or manifest.get("w"):
        lines.append(
            "Frozen run: {} Q{}W{}".format(
                manifest.get("particle_type", "unknown"),
                manifest.get("q2", "unknown"),
                manifest.get("w", "unknown"),
            )
        )
    if manifest.get("active_profile"):
        lines.append("Frozen profile: {}".format(manifest.get("active_profile")))
    if manifest.get("created_at"):
        lines.append("Frozen at: {}".format(manifest.get("created_at")))
    lines.append("ALLOW_CONFIG_DRIFT = {}".format(bool(allow_config_drift)))

    if rel_path == "src/utility/background_config.py":
        lines.extend(
            [
                "",
                "This cached iteration chain was frozen with a different background_config.py.",
                "Later iterations stop here so they do not mix a new analysis config with old frozen data-side products.",
                "",
                "Next steps:",
                "1. Restore src/utility/background_config.py to the version used by this frozen run.",
                "2. Or rerun the 0th iteration to create a new frozen manifest/cache chain.",
                "3. Or, for exploratory work only, set ALLOW_CONFIG_DRIFT = True.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "This file differs from the version recorded in the frozen manifest.",
                "If you intended to continue the old frozen chain, restore the matching file or rerun the 0th iteration.",
            ]
        )

    return "\n".join(lines)


def validate_manifest_hashes_against_repo(
    manifest,
    repo_root,
    strict_hash_paths=None,
    warn_only_hash_paths=None,
    allow_config_drift=None,
    emit_warnings=True,
    manifest_path=None,
):
    manifest_hashes = manifest.get("key_file_hashes", {})
    current_hashes = build_key_file_hashes(repo_root)
    warnings = []
    if strict_hash_paths is not None:
        strict_hash_paths = set(strict_hash_paths)
    if warn_only_hash_paths is not None:
        warn_only_hash_paths = set(warn_only_hash_paths)
    if allow_config_drift is None:
        allow_config_drift = ALLOW_CONFIG_DRIFT

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
        detailed_message = _build_hash_drift_message(
            rel_path,
            manifest_hash,
            current_hash,
            manifest=manifest,
            manifest_path=manifest_path,
            allow_config_drift=allow_config_drift,
        )
        should_fail = not allow_config_drift
        if strict_hash_paths is not None:
            should_fail = rel_path in strict_hash_paths and rel_path not in (warn_only_hash_paths or set())
        if should_fail:
            raise ValueError(detailed_message)
        if emit_warnings:
            print("WARNING: {}".format(warning))
        warnings.append(warning)

    return {
        "manifest_hashes": current_hashes,
        "warnings": warnings,
    }


def validate_iteration_inputs_against_manifest(
    manifest,
    current_inp_dict,
    t_bins,
    phi_bins,
    repo_root,
    manifest_path=None,
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

    current_inp_dict["manifest_validation_warnings"] = []
    hash_validation = validate_manifest_hashes_against_repo(
        manifest,
        repo_root,
        strict_hash_paths=strict_hash_paths,
        warn_only_hash_paths=warn_only_hash_paths,
        allow_config_drift=ALLOW_CONFIG_DRIFT,
        emit_warnings=True,
        manifest_path=manifest_path,
    )
    current_inp_dict["manifest_validation_warnings"] = list(hash_validation.get("warnings", []))

    return {
        "manifest_hashes": hash_validation.get("manifest_hashes", {}),
        "warnings": list(current_inp_dict.get("manifest_validation_warnings", [])),
    }
