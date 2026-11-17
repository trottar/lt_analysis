#!/usr/bin/env python3
"""Bounded existing-input ROOT integration; no main/launcher or production writes.

prepare and check are ROOT-free. run is farm-only after independent review and
user publication. No automatic fallback, simulation generation or projection.
"""
from __future__ import annotations

import argparse
from array import array
import csv
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import subprocess
import sys
import traceback

sys.dont_write_bytecode = True
SCHEMA = "d2a_root_integration/v1"
SCOPE = {"setting": "Q4p4W2p74", "phi": "Center", "epsilon": "highe", "t_index": 2}
PIONS = ("pi_n", "pi_delta", "pi_sidis")
ROLES = (*PIONS, "k_lambda_signal", "k_sigma0_signal")
REL, ABS = 1e-10, 1e-12
REPO = Path(__file__).resolve().parents[1]
SEED = REPO / "testing/fixtures/d2a_center_highe_t3_integration_seed.json"


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")) + "\n").encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def load_json(path):
    def pairs(items):
        result = {}
        for k, v in items:
            require(k not in result, "duplicate_json_key:" + k)
            result[k] = v
        return result
    def invalid(value):
        raise ValueError("nonfinite_json:" + value)
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_constant=invalid)


def write_json(path, value):
    Path(path).write_bytes(encoded(value))


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def close(a, b, label):
    if isinstance(a, list) or isinstance(b, list):
        require(isinstance(a, list) and isinstance(b, list) and len(a) == len(b), label + ":length")
        for i, (x, y) in enumerate(zip(a, b)):
            close(x, y, label + ":" + str(i))
    else:
        require(finite(a) and finite(b) and math.isclose(a, b, rel_tol=REL, abs_tol=ABS), label)


def seed_checks(seed):
    require(seed["schema"] == SCHEMA and seed["scope"] == SCOPE, "wrong_seed_scope")
    require(len(seed["t_edges"]) == 4 and len(seed["phi_edges"]) == 10, "wrong_native_grid")
    edges = seed["mm_edges"]
    require(len(edges) == 101 and edges[0] == .7 and edges[-1] == 1.3, "wrong_MM_grid")
    cfg = seed["config"]
    require(cfg["Q2"] == "4p4" and cfg["W"] == "2p74" and cfg["EPSSET"] in ("high", "highe"), "wrong_config")
    require(cfg["particle_subtraction_mode"] == "simc_shape_components" and cfg["pion_subtraction_scope"] == "t_bin", "wrong_model_scope")
    require(seed["d1a_fixture"]["support"]["complete"] is True, "incomplete_original_support")
    for name in (*PIONS, "k_lambda_signal"):
        r = seed["inputs"][name]
        require(len(r["sha256"]) == 64 and r["bytes"] > 0 and r["path"], "unbound_input:" + name)
    for name, h in seed["experimental"].items():
        check_hist(h, edges, "experimental:" + name)
    close(seed["experimental"]["kaon"]["values"], seed["d1a_fixture"]["kaon"]["values"], "support_histogram_binding")


def plan_checks(plan, seed, seed_sha):
    seed_checks(seed)
    require(plan["schema"] == SCHEMA and plan["scope"] == SCOPE, "wrong_plan_scope")
    require(plan["seed_sha256"] == seed_sha and set(plan["inputs"]) == set(ROLES), "wrong_plan_binding")
    for name in ROLES:
        x, expected = plan["inputs"][name], seed["inputs"].get(name)
        if name == "k_sigma0_signal":
            if x is None:
                continue
            require(isinstance(x, dict) and x.get("explicit_existing_source") is True,
                    "sigma0_requires_explicit_existing_source")
            require(x.get("path") and len(x.get("sha256", "")) == 64 and x.get("bytes", 0) > 0,
                    "sigma0_identity_missing")
        else:
            require(isinstance(x, dict) and x["sha256"] == expected["sha256"] and x["bytes"] == expected["bytes"],
                    "input_identity_changed:" + name)
        require(PurePosixPath(x["path"]).is_absolute() or PureWindowsPath(x["path"]).is_absolute(), "input_path_not_absolute:" + name)


def safe_output(path, repo, inputs):
    out = Path(path).resolve()
    repo = Path(repo).resolve()
    require(not out.exists(), "output_already_exists")
    require(out != repo and repo not in out.parents, "output_inside_repository")
    for x in inputs:
        p = Path(x).resolve()
        require(out != p and out not in p.parents, "output_overlaps_input")
    require(out.parent.is_dir(), "output_parent_missing")
    return out


SELECTION_KEYS = ("Q2", "W", "EPSSET", "POL", "ParticleType", "tmin", "tmax",
                  "cut_poly", "cut_mode", "poly_points", "TMIN_RESOLUTION_THRESHOLD")


def selection_identity(cfg):
    return {k: deepcopy(cfg.get(k)) for k in SELECTION_KEYS}


def file_token(path):
    p = Path(path); st = p.stat()
    return [os.path.abspath(p), int(st.st_mtime_ns), int(st.st_size),
            int(st.st_ctime_ns), int(st.st_dev), int(st.st_ino)]


def resolve_input_bindings(shapes, plan, cfg):
    """Resolve every configured input before any scientific loader may read it."""
    bindings = {}
    for name in ROLES:
        item = plan["inputs"][name]
        requested = item["path"] if item else None
        if name == "k_lambda_signal":
            resolved, diagnostic = shapes.resolve_kaon_simc_signal_root_filename(requested, cfg, "Center")
        elif name == "k_sigma0_signal":
            resolved, diagnostic = shapes.resolve_kaon_simc_sigma0_root_filename(requested, cfg, "Center")
        else:
            resolved, diagnostic = requested, {}
        if item is None:
            require(resolved is None, "undeclared_Sigma0_source:prepare_an_explicit_identity_plan")
            continue
        require(resolved is not None and Path(resolved).resolve(strict=True) == Path(requested).resolve(strict=True),
                "resolved_input_substitution:" + name)
        token = list(shapes._get_root_file_cache_token(resolved))
        expected = file_token(Path(requested).resolve(strict=True))
        require(len(token) == 6 and token[1:] == expected[1:], "resolved_input_token:" + name)
        bindings[name] = dict(requested_path=requested, resolved_path=os.path.abspath(resolved),
                              canonical_path=str(Path(resolved).resolve(strict=True)), token=token,
                              resolution=deepcopy(diagnostic))
    require(shapes.resolve_simc_tree_name(cfg) == (str(cfg.get("simc_tree_name", "h10")).strip() or "h10"),
            "resolved_tree_changed")
    return bindings


def check_loader_binding(identity, binding, item, observed, cfg, name):
    """Bind the actual population to the preflight file, not just its plan label."""
    require(binding["requested_path"] == item["path"] and binding["canonical_path"] == observed["path"],
            "resolved_input_path:" + name)
    token = list(identity["source"])
    require(len(token) == 6 and token == binding["token"] and token[0] == binding["resolved_path"]
            and token[1:] == observed["stat_token"][1:] and token[2] == item["bytes"], "loader_source_token:" + name)
    require(identity["tree"] == (str(cfg.get("simc_tree_name", "h10")).strip() or "h10"), "loader_tree:" + name)
    require(identity["selection"] == selection_identity(cfg), "loader_selection:" + name)


def check_staged_fit(snap, seed, enabled):
    require(snap["fallback_used"] is False, "full_fit_fallback")
    for side, suffix, key in (("pion", "pion_control", "H_pion_fit_step_overlays"),
                              ("kaon", "kaon_nosub", "H_kaon_fit_step_overlays")):
        require(snap["fit_status_" + side] == "success", "failed_full_fit:" + side)
        diag = snap["diagnostics"][side]
        require(diag["success"] is True and diag["fallback_used"] is False, "failed_fit_diagnostics:" + side)
        if side == "kaon" and seed["resolved_subtraction_config"]["kaon_nosub"]["pi_delta_signal_protected_fit"]["enabled"]:
            protected = diag["pi_delta_signal_protected_fit"]
            require(protected["enabled"] is True and protected["status"] == "success"
                    and protected["protected_fit_succeeded"] is True, "failed_protected_Delta_policy")
            # A successful, configured Lambda-only variant remains valid; a
            # zero-pi-delta failure policy is preserved as evidence, not success.
            require(protected["fit_variant"] in ("lambda_sigma0_protected", "lambda_only_protected_fallback"),
                    "protected_Delta_variant")
        order = diag["fit_order"]
        require(isinstance(order, list) and len(set(order)) == len(order) and set(PIONS) <= set(order), "incomplete_fit_order:" + side)
        passes = diag["staged_pass_history"]
        require(type(diag["n_passes"]) is int and diag["n_passes"] > 0 and len(passes) == diag["n_passes"], "missing_stage_passes:" + side)
        steps = snap[key]
        require(len(steps) == len(passes) * len(order), "missing_component_stages:" + side)
        raw = snap["staged_amplitudes_raw"][suffix]
        require(all(finite(raw[c]) and raw[c] >= 0 for c in PIONS), "unavailable_raw_stage:" + side)
        for pi, summary in enumerate(passes, 1):
            for ci, component in enumerate(order):
                stage = summary[component]; step = steps[(pi - 1) * len(order) + ci]
                require(stage["success"] is True and finite(stage["amplitude"]) and stage["amplitude"] >= 0,
                        "failed_component_stage:" + side + ":" + component)
                require(step["component_name"] == component and step["pass_index"] == pi, "stage_identity:" + side)
                close(step["amplitude"], stage["amplitude"], "stage_amplitude")
                objective = step["fit_diagnostics"]
                bins = objective["fit_bin_indices"]
                require(isinstance(bins, list) and bins and len(set(bins)) == len(bins)
                        and all(type(i) is int and 1 <= i <= 100 for i in bins), "invalid_objective_bins")
                require(type(stage["n_fit_bins"]) is int and stage["n_fit_bins"] == len(bins)
                        and step["n_fit_bins"] == len(bins), "stage_bin_count")
                for field in ("fit_y", "fit_sigma", "fit_template"):
                    require(len(objective[field]) == len(bins) and all(finite(x) for x in objective[field]), "missing_objective:" + field)
                require(all(x > 0 for x in objective["fit_sigma"]), "invalid_objective_sigma")
                if objective["estimator"] == "least_squares":
                    num = math.fsum(t*y/(e*e) for t,y,e in zip(objective["fit_template"], objective["fit_y"], objective["fit_sigma"]))
                    den = math.fsum((t/e)**2 for t,e in zip(objective["fit_template"], objective["fit_sigma"]))
                    require(den > 0, "objective_zero_support")
                    close(objective["weighted_numerator"], num, "objective_numerator")
                    close(objective["weighted_denominator"], den, "objective_denominator")
                    close(stage["amplitude"], max(num/den, 0.), "objective_amplitude")
                else:
                    require(objective["estimator"] == "window_integral", "unknown_objective")
                    den, num = math.fsum(objective["fit_template"]), math.fsum(objective["fit_y"])
                    require(den > 0, "objective_zero_support")
                    close(objective["template_sum"], den, "objective_template_sum")
                    close(objective["target_sum"], num, "objective_target_sum")
                    close(stage["amplitude"], max(num/den, 0.), "objective_amplitude")
                marker = objective.get("d1a_original_support")
                require((marker is not None) is (enabled and side == "kaon" and component == "pi_sidis"), "misplaced_D1a_marker")
                if marker is not None:
                    require(marker["candidate"] == "D1.a" and marker["original_support_schema"] == "kaon_sidis_original_support/v1"
                            and marker["minimum_fit_bins"] == 2 and len(bins) >= 2, "invalid_D1a_marker")
                    residual = step["H_residual_input"];check_hist(residual, seed["mm_edges"], "D1a_residual")
                    windows = step["amplitude_windows"] or step["anchor_windows"]
                    candidates = [i for i in range(1,101) if diag["fit_min"] <= (seed["mm_edges"][i-1]+seed["mm_edges"][i])/2 <= diag["fit_max"]
                                  and (not windows or any(a <= (seed["mm_edges"][i-1]+seed["mm_edges"][i])/2 <= b for a,b in windows))
                                  and not any(a <= (seed["mm_edges"][i-1]+seed["mm_edges"][i])/2 <= b for a,b in step["excluded_windows"])]
                    sources = seed["d1a_fixture"]["support"]["sources"]
                    excluded = [i for i in candidates if all(v["entries"][i] == 0 and v["sumw2"][i] == 0 for v in sources.values())
                                and seed["experimental"]["kaon"]["values"][i] == 0 and seed["experimental"]["kaon"]["stored_sumw2"][i] == 0
                                and residual["stored_sumw2"][i] > 0]
                    require(marker["excluded_bin_indices"] == excluded and bins == [i for i in candidates if i not in excluded], "D1a_objective_decision")


def identify_inputs(plan):
    records = {}
    for name, item in plan["inputs"].items():
        if item is None:
            records[name] = {"status": "not_configured"}
            continue
        p = Path(item["path"]).resolve(strict=True)
        require(p.is_file() and p.stat().st_size == item["bytes"], "input_size:" + name)
        h = file_hash(p)
        require(h == item["sha256"], "input_hash:" + name)
        records[name] = {"status": "available", "path": str(p), "bytes": p.stat().st_size, "sha256": h, "stat_token": file_token(p)}
    return records


def git(repo, *args):
    return subprocess.check_output(["git", "--no-optional-locks", "-C", str(repo), *args], stderr=subprocess.PIPE)


def source_state(repo, seed, require_tracked=True):
    files = {}
    for name, expected in seed["scientific_source_LF_sha256"].items():
        b = (repo / name).read_bytes()
        h = hashlib.sha256(b.replace(b"\r\n", b"\n")).hexdigest()
        require(h == expected, "scientific_source_changed:" + name)
        files[name] = h
    if require_tracked:
        for name in ("testing/run_d2a_root_integration.py", "testing/fixtures/d2a_center_highe_t3_integration_seed.json"):
            git(repo, "ls-files", "--error-unmatch", name)
            require(git(repo, "show", "HEAD:" + name).replace(b"\r\n", b"\n") == (repo / name).read_bytes().replace(b"\r\n", b"\n"), "mechanism_unpublished_or_dirty:" + name)
    index_path = Path(git(repo, "rev-parse", "--git-path", "index").decode().strip())
    if not index_path.is_absolute():
        index_path = repo / index_path
    return {"HEAD": git(repo, "rev-parse", "HEAD").decode().strip(),
            "refs_sha256": hashlib.sha256(git(repo, "show-ref")).hexdigest(),
            "status_sha256": hashlib.sha256(git(repo, "status", "--porcelain=v2", "--untracked-files=all")).hexdigest(),
            "index_sha256": file_hash(index_path), "scientific_source_LF_sha256": files}


def check_hist(h, edges, label):
    require(isinstance(h, dict) and h.get("status") == "available", label + ":unavailable")
    close(h["edges"], edges, label + ":axis")
    require(len(h["values"]) == len(edges) + 1, label + ":bins")
    close([v * v for v in h["errors"]], h["stored_sumw2"], label + ":Sumw2")
    require(all(finite(x) for x in h["values"]) and all(x >= 0 for x in h["stored_sumw2"]), label + ":finite")


def moments(m, n, label):
    require(type(m["count"]) is int and m["count"] >= 0, label + ":count")
    require(len(m["values"]) == n and len(m["variances"]) == n, label + ":bins")
    require(all(finite(x) for x in m["values"] + m["variances"]), label + ":finite")
    require(all(x >= 0 for x in m["variances"]), label + ":negative_variance")
    close(math.fsum(m["values"]), m["signed"], label + ":signed")
    close(math.fsum(m["variances"]), m["sumw2"], label + ":sumw2")
    close(math.fsum(m["values"][1:-1]), m["normal_sumw"], label + ":normal")
    close(math.fsum(m["variances"][1:-1]), m["normal_sumw2"], label + ":normal_variance")
    require(m["absolute"] + ABS >= abs(m["signed"]), label + ":absolute")
    require(m["status"] == ("available" if m["count"] and m["normal_sumw"] > 0 else "unsupported"), label + ":availability")


def sum_moments(whole, parts, label):
    require(whole["count"] == sum(p["count"] for p in parts), label + ":count")
    for key in ("signed", "absolute", "sumw2"):
        close(whole[key], math.fsum(p[key] for p in parts), label + ":" + key)
    for key in ("values", "variances"):
        close(whole[key], [math.fsum(v) for v in zip(*(p[key] for p in parts))], label + ":" + key)


def normalize_check(raw, hist, edges, label):
    check_hist(hist, edges, label)
    f = raw["normal_sumw"] if raw["normal_sumw"] > 0 else 1.
    close(hist["values"], [v / f for v in raw["values"]], label + ":one_normalization")
    close(hist["stored_sumw2"], [v / (f * f) for v in raw["variances"]], label + ":variance_normalization")


def compatible(actual, expected, path="config"):
    """Compare retained fields only; explicit old unavailable fields stay unknown."""
    if isinstance(expected, dict) and expected.get("status") == "unavailable":
        return
    if isinstance(expected, dict):
        require(isinstance(actual, dict), path)
        for k, v in expected.items():
            # Availability is a result of the explicitly declared optional input.
            if k in ("sigma0_active", "sigma0_requested"):
                continue
            require(k in actual, path + ":missing:" + k)
            compatible(actual[k], v, path + ":" + k)
    elif isinstance(expected, list):
        require(isinstance(actual, (list, tuple)) and len(actual) == len(expected), path)
        for i, (a, e) in enumerate(zip(actual, expected)):
            compatible(a, e, path + ":" + str(i))
    elif isinstance(expected, bool):
        require(actual is expected, path)
    elif finite(expected):
        close(actual, expected, path)
    else:
        require(actual == expected, path)


def check_result(result, plan, seed, seed_sha):
    plan_checks(plan, seed, seed_sha)
    require(result["schema"] == SCHEMA and result["scope"] == SCOPE, "wrong_result_scope")
    require(result["runtime"]["kind"] == "PyROOT" and result["runtime"]["ROOT_version"], "ROOT_execution_not_reported")
    require(result["plan_digest"] == digest(plan) and result["seed_sha256"] == seed_sha, "result_binding")
    require(result["preservation"]["before"] == result["preservation"]["after"], "source_or_input_changed")
    require(result["preservation"]["before"]["source"]["scientific_source_LF_sha256"] == seed["scientific_source_LF_sha256"], "source_identity")
    for name in ROLES:
        expected = plan["inputs"][name]
        observed = result["preservation"]["before"]["inputs"][name]
        if expected is None:
            require(observed == {"status": "not_configured"}, "optional_input_hidden")
        else:
            require(observed["status"] == "available" and observed["sha256"] == expected["sha256"] and observed["bytes"] == expected["bytes"], "input_binding:" + name)
    require(set(result["input_bindings"]) == {n for n in ROLES if plan["inputs"][n] is not None}, "resolved_input_inventory")
    require(set(result["templates"]) == set(ROLES), "five_channel_inventory")
    for name, payload in result["templates"].items():
        if plan["inputs"][name] is None:
            require(payload["status"] == "not_configured" and payload["scope_usable"] is False, "optional_scope_substitution")
            continue
        identity = payload["identity"]
        check_loader_binding(identity, result["input_bindings"][name], plan["inputs"][name],
                             result["preservation"]["before"]["inputs"][name], seed["config"], name)
        require(identity["schema"] == 2 and identity["population_definition"] == 3 and identity["reuse_safe"] is True, "population_identity:" + name)
        require(identity["component"] == name and identity["particle_type"] == "kaon" and identity["phi_setting"] == "Center", "identity_semantics:" + name)
        require(identity["normalization"] == "normal_bins_once_after_event_accumulation" and identity["weight"] == "getattr(evt, iter_weight, 1.0)", "weight_normalization_semantics:" + name)
        close(identity["mm_selection_window"], [seed["config"]["mm_min"], seed["config"]["mm_max"]], "identity:MM_selection")
        close(identity["t_edges"], seed["t_edges"], "identity:t")
        close(identity["phi_edges"], seed["phi_edges"], "identity:phi")
        close(identity["hgcer_selection"]["points"], seed["hgcer_points"], "HGCer_geometry")
        require(identity["mm_axis"] == {"nbins": 100, "min": .7, "max": 1.3, "flows": True}, "identity:MM")
        support = payload["support"]
        for k in ("setting", "outside_t"):
            moments(support[k], 102, name + ":" + k)
        require(set(support["parents"]) == {"t_bin1", "t_bin2", "t_bin3"}, "three_parent_inventory")
        for ti in (1, 2, 3):
            key = "t_bin" + str(ti)
            parent, children, outside = support["parents"][key], support["children"][key], support["outside_phi"][key]
            moments(parent, 102, name + key);moments(outside, 102, name + key + "outside")
            require(set(children) == {"phi_bin" + str(i) for i in range(1, 10)}, "nine_child_inventory")
            for pi, child in children.items():
                moments(child, 102, name + key + pi)
                normalize_check(child, payload["children"][key][pi], seed["mm_edges"], name + key + pi)
            sum_moments(parent, [*children.values(), outside], name + key + ":closure")
            normalize_check(parent, payload["parents"][key], seed["mm_edges"], name + key)
        sum_moments(support["setting"], [*support["parents"].values(), support["outside_t"]], name + ":setting_closure")
        normalize_check(support["setting"], payload["setting"], seed["mm_edges"], name + ":setting")
        if name in PIONS:
            expected = seed["expected_pion_support"][name]
            for key in ("setting", "outside_t"):
                for field in ("count", "signed", "absolute", "sumw2", "values", "variances"):
                    close(support[key][field], expected[key][field], name + key + field)
            for key in support["parents"]:
                for field in ("count", "signed", "absolute", "sumw2", "values", "variances"):
                    close(support["parents"][key][field], expected["parents"][key][field], name + key + field)
                    close(support["outside_phi"][key][field], expected["outside_phi"][key][field], name + key + "outside" + field)
                for pi in support["children"][key]:
                    for field in ("count", "values", "variances"):
                        close(support["children"][key][pi][field], expected["children"][key][pi][field], name + key + pi + field)
        require(payload["scope_usable"] is True, "required_t3_scope_unsupported:" + name)
        cache = result["cache_checks"][name]
        require(cache["cold_entries"] > cache["before_entries"] and cache["warm_entries"] == cache["cold_entries"] and cache["probe_entries"] > cache["warm_entries"], "cache_key_invalidation:" + name)
        require(cache["warm_equal"] is True and cache["clone_independent"] is True and cache["probe_parent_equal"] is True and cache["probe_phi_edges"] != seed["phi_edges"], "cache_clone_values:" + name)
    for mode in ("disabled", "enabled"):
        fit = result["fits"][mode]
        check_staged_fit(fit["snapshot"], seed, mode == "enabled")
        require(fit["d1a_enabled"] is (mode == "enabled") and fit["input_unchanged"] is True, "D1a_input_mutation")
        for key, source in (("H_kaon_nosub_input", "kaon"), ("H_pion_control_input", "pion")):
            h = fit["snapshot"][key];check_hist(h, seed["mm_edges"], key)
            close(h["values"], seed["experimental"][source]["values"], key + ":original")
            close(h["stored_sumw2"], seed["experimental"][source]["stored_sumw2"], key + ":original_errors")
        for key in ("A_n", "A_delta", "A_sidis", "B_n", "B_delta", "B_sidis"):
            require(finite(fit["snapshot"][key]), "missing_actual_fit_amplitude:" + key)
        for side in ("kaon_nosub", "pion_control"):
            require(set(PIONS) <= set(fit["snapshot"]["staged_amplitudes_raw"][side]), "missing_staged_results")
        for key in ("H_kaon_fit_step_overlays", "H_pion_fit_step_overlays"):
            require(isinstance(fit["snapshot"].get(key), list) and fit["snapshot"][key], "missing_stage_history")
        require(fit["snapshot"]["pion_component_alignment"] == seed["alignment"], "frozen_alignment_changed")
        def support_markers(value):
            if isinstance(value, dict):
                return ([value["d1a_original_support"]] if "d1a_original_support" in value else []) + [x for v in value.values() for x in support_markers(v)]
            if isinstance(value, list):
                return [x for v in value for x in support_markers(v)]
            return []
        markers = support_markers(fit["snapshot"]["H_kaon_fit_step_overlays"])
        require(bool(markers) is (mode == "enabled"), "D1a_full_fit_objective_marker")
        compatible(fit["snapshot"]["resolved_subtraction_config"], seed["resolved_subtraction_config"])
    d1 = result["d1a_boundary"]
    require(d1["physical_histograms_unchanged"] is True, "D1a_physical_bins_changed")
    require(d1["disabled_bins"] == list(range(61, 69)) and d1["enabled_bins"] == [61, 63, 64, 65, 66, 67, 68] and d1["excluded_bins"] == [62], "D1a_empty_bin_objective")
    close(d1["disabled_sidis"], 0., "D1a_disabled_boundary")
    close(d1["enabled_sidis"], .218867157746124, "D1a_fixed_fixture_boundary")
    require(result["templates_unchanged_after_fits"] is True, "fit_mutated_source_templates")
    require(result["cut_state_restored"] is True and result["full_analysis_executed"] is False, "owner_scope")
    return {"status": "PASS", "scope": "detached_corrected_shared_entry_ROOT_integration", "physical_background_accepted": False}


def restore_hist(ROOT, saved, name):
    edges = saved["edges"]
    hist = ROOT.TH1D(name, name, len(edges) - 1, array("d", edges))
    hist.SetDirectory(0);hist.Sumw2()
    for i, (v, e) in enumerate(zip(saved["values"], saved["errors"])):
        hist.SetBinContent(i, v);hist.SetBinError(i, e)
    # Entry metadata is unavailable in saved spectra and unused by pinned fitter.
    hist.SetEntries(0)
    return hist


@contextmanager
def preserved_cuts(cuts):
    before = dict(vars(cuts))
    try:
        yield
    finally:
        vars(cuts).clear();vars(cuts).update(before)


def root_payload(payload, histogram, resolver):
    return {"status": "available", "identity": deepcopy(payload.get("event_population_identity")),
            "support": deepcopy(payload.get("event_population_support")),
            "diagnostics": deepcopy(payload.get("diagnostics")),
            "setting": histogram(payload["setting_shape_full"]),
            "parents": {k: histogram(v) for k, v in payload["t_shapes"].items()},
            "children": {k: {p: histogram(h) for p, h in v.items()} for k, v in payload["binned_shapes"].items()},
            "scope_usable": resolver(payload, "t_bin3", t_bin_index=2) is not None}


def compute_root(plan, seed):
    for name in ("src/cuts", "src/utility", "src/binning"):
        sys.path.insert(0, str(REPO / name))
    import ROOT
    import numpy as np
    import scipy
    import apply_cuts as cuts
    import pion_component_shapes as shapes
    import pion_component_fits as fits
    from pion_background_chain_audit import snapshot, histogram
    ROOT.gROOT.SetBatch(True)
    require(ROOT.gROOT.GetVersion() and getattr(ROOT, "__file__", None), "real_PyROOT_required")
    cfg = deepcopy(seed["config"])
    bindings = resolve_input_bindings(shapes, plan, cfg)
    before_cuts = dict(vars(cuts))
    templates, cache_checks, actual_payloads = {}, {}, {}
    with preserved_cuts(cuts):
        for name in ROLES:
            item = plan["inputs"][name]
            if item is None:
                missing = shapes.load_kaon_simc_sigma0_shape(None, cfg, "Center", t_bins=seed["t_edges"], phi_bins=seed["phi_edges"])
                require(fits.resolve_scope_single_shape(missing, "t_bin3", t_bin_index=2) is None, "optional_signal_substitution")
                templates[name] = {"status": "not_configured", "scope_usable": False, "diagnostics": snapshot(missing["diagnostics"])}
                actual_payloads[name] = missing
                continue
            def load(phi):
                kwargs = dict(t_bins=seed["t_edges"], phi_bins=phi, context="d2a_integration_" + name)
                if name in PIONS:
                    return shapes.load_pion_simc_component_shape(item["path"], cfg, "Center", "kaon", name, **kwargs)
                fn = shapes.load_kaon_simc_signal_shape if name == "k_lambda_signal" else shapes.load_kaon_simc_sigma0_shape
                return fn(item["path"], cfg, "Center", **kwargs)
            before_n = len(shapes._SIMC_COMPONENT_SHAPE_CACHE)
            cold = load(seed["phi_edges"]);cold_n = len(shapes._SIMC_COMPONENT_SHAPE_CACHE)
            warm = load(seed["phi_edges"]);warm_n = len(shapes._SIMC_COMPONENT_SHAPE_CACHE)
            cold_record = root_payload(cold, histogram, fits.resolve_scope_single_shape)
            check_loader_binding(cold_record["identity"], bindings[name], item,
                                 {"path": bindings[name]["canonical_path"], "stat_token": file_token(Path(item["path"]).resolve(strict=True))}, cfg, name)
            warm_record = root_payload(warm, histogram, fits.resolve_scope_single_shape)
            warm_equal = (digest(cold_record) == digest(warm_record))
            original = cold["t_shapes"]["t_bin3"].GetBinContent(1)
            warm["t_shapes"]["t_bin3"].SetBinContent(1, original + 123.)
            third = load(seed["phi_edges"])
            independent = (cold["t_shapes"]["t_bin3"].GetBinContent(1) == original and third["t_shapes"]["t_bin3"].GetBinContent(1) == original)
            probe_phi = [seed["phi_edges"][0], seed["phi_edges"][-1]]
            probe = load(probe_phi);probe_n = len(shapes._SIMC_COMPONENT_SHAPE_CACHE)
            probe_parent_equal = all(histogram(cold["t_shapes"][k]) == histogram(probe["t_shapes"][k]) for k in cold["t_shapes"])
            cache_checks[name] = dict(before_entries=before_n, cold_entries=cold_n, warm_entries=warm_n, probe_entries=probe_n, warm_equal=warm_equal, clone_independent=independent, probe_parent_equal=probe_parent_equal, probe_phi_edges=probe_phi)
            templates[name] = cold_record;actual_payloads[name] = cold
    restored = dict(vars(cuts)) == before_cuts
    pion_payload = {"components": {c: actual_payloads[c] for c in PIONS}}
    components = fits.resolve_scope_component_shapes(pion_payload, "t_bin3", t_bin_index=2)
    lam = fits.resolve_scope_single_shape(actual_payloads["k_lambda_signal"], "t_bin3", t_bin_index=2)
    sig = fits.resolve_scope_single_shape(actual_payloads["k_sigma0_signal"], "t_bin3", t_bin_index=2)
    fit_results = {}
    for mode in ("disabled", "enabled"):
        hc = restore_hist(ROOT, seed["experimental"]["pion"], "d2a_pion_" + mode)
        hk = restore_hist(ROOT, seed["experimental"]["kaon"], "d2a_kaon_" + mode)
        inputs_before = (histogram(hc), histogram(hk))
        inp = deepcopy(cfg);inp["d1a_kaon_sidis_support_enabled"] = mode == "enabled"
        fit = fits.build_particle_subtraction_component_result(hc, hk, components, inp, "t_bin3", kaon_signal_shape=lam, kaon_sigma0_shape=sig,
            kaon_sigma0_source_diagnostics=actual_payloads["k_sigma0_signal"]["diagnostics"], mm_offset_data=seed["template_mm_offset_data"], phi_setting="Center", context="d2a_integrate_" + mode,
            pion_component_alignment=deepcopy(seed["alignment"]), kaon_experimental_support=deepcopy(seed["d1a_fixture"]["support"]))
        fit_results[mode] = {"d1a_enabled": mode == "enabled", "input_unchanged": inputs_before == (histogram(hc), histogram(hk)), "snapshot": snapshot(fit)}
    f = seed["d1a_fixture"]
    k, n, sid = [restore_hist(ROOT, f[name], "d1a_original_" + name) for name in ("kaon", "neutron", "sidis")]
    unchanged = [histogram(x) for x in (k, n, sid)]
    windows = {x["component_name"]: x["anchor_windows"] for x in f["stages"]};windows["pi_delta"] = windows["pi_sidis"]
    passes = []
    for support in (None, deepcopy(f["support"])):
        amps, errors, fitted, steps = fits._run_staged_component_pass(k, {"pi_n": n, "pi_sidis": sid, "pi_delta": sid}, ["pi_n", "pi_sidis"], .7, 1.3, windows, {}, {}, [], {}, "A", 0, kaon_sidis_experimental_support=support)
        passes.append({"amps": snapshot(amps), "steps": snapshot(steps)})
    d1 = {"disabled_bins": passes[0]["steps"][1]["fit_diagnostics"]["fit_bin_indices"], "enabled_bins": passes[1]["steps"][1]["fit_diagnostics"]["fit_bin_indices"],
          "excluded_bins": passes[1]["steps"][1]["fit_diagnostics"]["d1a_original_support"]["excluded_bin_indices"], "disabled_sidis": passes[0]["amps"]["pi_sidis"], "enabled_sidis": passes[1]["amps"]["pi_sidis"],
          "physical_histograms_unchanged": unchanged == [histogram(x) for x in (k, n, sid)], "fixed_fixture_results": passes, "meaning": "Historical fixed-input ROOT objective parity, not the new D2.a fitted amplitude"}
    template_preservation = all(root_payload(actual_payloads[c], histogram, fits.resolve_scope_single_shape) == templates[c] for c in ROLES if plan["inputs"][c] is not None)
    imports = {}
    for name in ("ROOT", "numpy", "scipy", "ltsep", "pion_component_shapes", "pion_component_fits", "utility"):
        module = sys.modules.get(name)
        filename = getattr(module, "__file__", None)
        if filename and Path(filename).is_file():
            imports[name] = {"status": "available", "lexical_path": str(filename), "resolved_path": str(Path(filename).resolve()), "entry_file_sha256": file_hash(filename)}
        else:
            imports[name] = {"status": "unavailable", "reason": "module entry file not reported"}
    return dict(runtime={"kind": "PyROOT", "ROOT_version": ROOT.gROOT.GetVersion(), "ROOT_module": str(Path(ROOT.__file__).resolve()), "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "imports": imports}, input_bindings=bindings, templates=templates, cache_checks=cache_checks, fits=fit_results,
                d1a_boundary=d1, cut_state_restored=restored, templates_unchanged_after_fits=template_preservation, full_analysis_executed=False)


def render_outputs(result, seed, output):
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figures = output / "figures";figures.mkdir()
    def save(fig, name):
        fig.tight_layout()
        for suffix in ("png", "svg"):
            fig.savefig(figures / (name + "." + suffix), dpi=240)
        plt.close(fig)
    edges = seed["mm_edges"]
    for c in PIONS:
        p = result["templates"][c];raw = p["support"]["parents"]["t_bin3"]
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.stairs(p["parents"]["t_bin3"]["values"][1:-1], edges)
        ax.set(xlabel="Missing mass [GeV]", ylabel="Unit shape / native bin", title="Center/high t3 " + c + " — actual ROOT direct-event template")
        ax.grid(alpha=.2);save(fig, c + "_direct_t3")
        parts = [*p["support"]["children"]["t_bin3"].values(), p["support"]["outside_phi"]["t_bin3"]]
        diff = [a - math.fsum(x) for a, x in zip(raw["values"], zip(*(x["values"] for x in parts)))]
        fig, ax = plt.subplots(figsize=(9, 5));ax.stairs(diff[1:-1], edges)
        ax.set(xlabel="Missing mass [GeV]", ylabel="Raw parent − native children − outside", title=c + " ROOT closure; flows retained in JSON")
        ax.grid(alpha=.2);save(fig, c + "_closure")
    for side, key in (("kaon", "H_kaon_fit_total_refined"), ("pion", "H_pion_fit_total_refined")):
        fig, ax = plt.subplots(figsize=(9, 5));ax.stairs(seed["experimental"][side]["values"][1:-1], edges, color="black", label="Retained experimental input")
        for mode, style in (("disabled", "--"), ("enabled", "-")):
            h = result["fits"][mode]["snapshot"][key];ax.stairs(h["values"][1:-1], edges, linestyle=style, label="D1.a " + mode)
        ax.set(xlabel="Missing mass [GeV]", ylabel="Existing normalized yield / native bin", title="Center/high t3 " + side + " — conditional integration fit")
        ax.legend();ax.grid(alpha=.2);save(fig, side + "_D1a_modes")
    with (output / "fit_scalars.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f);keys=("A_n", "A_delta", "A_sidis", "B_n", "B_delta", "B_sidis")
        w.writerow(("D1a_mode", *keys))
        for mode in ("disabled", "enabled"):
            w.writerow((mode, *(result["fits"][mode]["snapshot"][k] for k in keys)))


def verify_directory(output, seed_path=SEED):
    out = Path(output).resolve();manifest = load_json(out / "SHA256SUMS.json")
    expected = {"owner.json", "result.json", "plan.json", "seed.json", "run.log", "checker.json", "fit_scalars.csv"}
    names = [c + suffix for c in PIONS for suffix in ("_direct_t3", "_closure")] + ["kaon_D1a_modes", "pion_D1a_modes"]
    expected |= {"figures/" + n + "." + ext for n in names for ext in ("png", "svg")}
    require(set(manifest) == expected, "artifact_inventory")
    for name, item in manifest.items():
        p = out / name
        require(out in p.resolve().parents and not p.is_symlink() and p.is_file() and p.stat().st_size == item["bytes"] and file_hash(p) == item["sha256"], "artifact_integrity:" + name)
    require(not (out / "receipt.json").is_symlink() and not (out / "SHA256SUMS.json").is_symlink(), "external_integrity_metadata")
    require({p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()} == expected | {"SHA256SUMS.json", "receipt.json"}, "unexpected_artifacts")
    owner = load_json(out / "owner.json")
    require(owner["status"] == "passed" and owner["schema"] == SCHEMA and owner["scope"] == SCOPE, "owner_not_passed")
    seed = load_json(seed_path);seed_sha = file_hash(seed_path)
    require(file_hash(out / "seed.json") == seed_sha, "seed_authority")
    checked = check_result(load_json(out / "result.json"), load_json(out / "plan.json"), seed, seed_sha)
    require(load_json(out / "checker.json") == checked, "checker_disagrees")
    receipt = load_json(out / "receipt.json")
    require(receipt["status"] == "passed" and receipt["manifest_sha256"] == file_hash(out / "SHA256SUMS.json"), "receipt_binding")
    for name in names:
        png = (out / "figures" / (name + ".png")).read_bytes()
        require(png.startswith(b"\x89PNG\r\n\x1a\n") and len(png) >= 24, "invalid_PNG:" + name)
        width, height = int.from_bytes(png[16:20], "big"), int.from_bytes(png[20:24], "big")
        require(width >= 800 and height >= 400, "figure_too_small:" + name)
        svg = (out / "figures" / (name + ".svg")).read_text(encoding="utf-8")
        require("<svg" in svg and "</svg>" in svg, "invalid_SVG:" + name)
    with (out / "fit_scalars.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 2 and [r["D1a_mode"] for r in rows] == ["disabled", "enabled"], "scalar_csv_scope")
    result = load_json(out / "result.json")
    for row in rows:
        for key in ("A_n", "A_delta", "A_sidis", "B_n", "B_delta", "B_sidis"):
            close(float(row[key]), result["fits"][row["D1a_mode"]]["snapshot"][key], "scalar_csv_binding")
    return checked


def run_owner(plan_path, plan_sha, output, seed_path=SEED, repo=REPO,
              backend=compute_root, renderer=render_outputs, state_reader=source_state):
    require(file_hash(plan_path) == plan_sha, "plan_hash")
    seed, plan = load_json(seed_path), load_json(plan_path);seed_sha = file_hash(seed_path)
    plan_checks(plan, seed, seed_sha)
    out = safe_output(output, repo, [plan_path, seed_path, *(x["path"] for x in plan["inputs"].values() if x)])
    out.mkdir()
    write_json(out / "owner.json", {"status": "running", "schema": SCHEMA, "scope": SCOPE})
    before = None
    try:
        before = {"source": state_reader(repo, seed), "inputs": identify_inputs(plan)}
        write_json(out / "plan.json", plan);(out / "seed.json").write_bytes(Path(seed_path).read_bytes())
        with (out / "run.log").open("w", encoding="utf-8") as log, redirect_stdout(log), redirect_stderr(log):
            result = backend(plan, seed)
        after = {"source": state_reader(repo, seed), "inputs": identify_inputs(plan)}
        require(before == after, "source_or_input_changed_during_run")
        result.update(schema=SCHEMA, scope=SCOPE, plan_digest=digest(plan), seed_sha256=seed_sha,
                      preservation={"before": before, "after": after})
        write_json(out / "result.json", result)  # Retain failed/policy fit evidence, unsealed.
        checked = check_result(result, plan, seed, seed_sha)
        renderer(result, seed, out)
        require(before == {"source": state_reader(repo, seed), "inputs": identify_inputs(plan)}, "source_or_input_changed_before_seal")
        write_json(out / "checker.json", checked)
        write_json(out / "owner.json", {"status": "passed", "schema": SCHEMA, "scope": SCOPE})
        members = {p.relative_to(out).as_posix(): {"bytes": p.stat().st_size, "sha256": file_hash(p)} for p in out.rglob("*") if p.is_file()}
        write_json(out / "SHA256SUMS.json", members)
        write_json(out / "receipt.json", {"status": "passed", "schema": SCHEMA, "manifest_sha256": file_hash(out / "SHA256SUMS.json"), "physical_background_accepted": False})
        verify_directory(out, seed_path)
        return out
    except BaseException as exc:
        preservation = {"before": before, "after": None}
        if before is not None:
            try:
                preservation["after"] = {"source": state_reader(repo, seed), "inputs": identify_inputs(plan)}
            except BaseException as guard_error:
                preservation["guard_error"] = str(guard_error)
        write_json(out / "owner.json", {"status": "failed", "schema": SCHEMA, "scope": SCOPE, "error_type": type(exc).__name__, "error": str(exc), "preservation": preservation})
        # A failed owner never leaves a success marker, even if final checking fails.
        receipt = out / "receipt.json"
        if receipt.exists():
            write_json(receipt, {"status": "failed", "error": str(exc)})
        raise


def prepare(args):
    seed = load_json(SEED);seed_checks(seed)
    inputs = deepcopy(seed["inputs"])
    if args.pion_root_dir:
        for c, subdir in (("pi_n", "neutron"), ("pi_delta", "delta"), ("pi_sidis", "sidis")):
            inputs[c]["path"] = str((args.pion_root_dir / subdir / Path(inputs[c]["path"]).name).resolve())
    if args.lambda_root:
        inputs["k_lambda_signal"]["path"] = str(args.lambda_root.resolve())
    if args.sigma0_root:
        require(args.sigma0_sha256 is not None, "sigma0_hash_required")
        p = args.sigma0_root.resolve(strict=True)
        require(file_hash(p) == args.sigma0_sha256, "sigma0_hash")
        inputs["k_sigma0_signal"] = {"path": str(p), "bytes": p.stat().st_size, "sha256": args.sigma0_sha256, "explicit_existing_source": True}
    else:
        require(args.sigma0_sha256 is None, "sigma0_path_required")
    plan = {"schema": SCHEMA, "scope": SCOPE, "seed_sha256": file_hash(SEED), "inputs": inputs}
    plan_checks(plan, seed, file_hash(SEED))
    path = safe_output(args.output, REPO, [SEED, *(x["path"] for x in inputs.values() if x)])
    with path.open("xb") as stream:
        stream.write(encoded(plan))
    print(json.dumps({"plan": str(path), "sha256": file_hash(path), "ROOT_executed": False, "input_availability": "Checked only by run preflight; original recorded locations may need explicit same-hash relocation"}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare");p.add_argument("--output", type=Path, required=True)
    p.add_argument("--pion-root-dir", type=Path);p.add_argument("--lambda-root", type=Path)
    p.add_argument("--sigma0-root", type=Path);p.add_argument("--sigma0-sha256")
    p = sub.add_parser("run");p.add_argument("--plan", type=Path, required=True);p.add_argument("--plan-sha256", required=True);p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("check");p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":prepare(args)
        elif args.command == "run":print(json.dumps({"output": str(run_owner(args.plan, args.plan_sha256, args.output)), "status": "passed"}))
        else:print(json.dumps(verify_directory(args.output)))
        return 0
    except Exception as exc:
        print(type(exc).__name__ + ": " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
