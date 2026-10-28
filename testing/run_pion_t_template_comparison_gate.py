#!/usr/bin/env python3
"""Own one isolated Center/Left lowe template comparison, then verify and collect.

Farm execution requires independent diff, publication, synchronization and
readiness review. This module imports operational helpers only, never physics.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
from contextlib import contextmanager, redirect_stdout
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import zipfile

try:
    from . import run_e8_2_left_lowe_scientific_audit_gate as debug
except ImportError:
    import run_e8_2_left_lowe_scientific_audit_gate as debug

isolation = debug.isolation
collector, require, sha256, git = debug.collector, debug.require, debug.sha256, debug.git
BASE_HEAD = "554ad84630ad77920ce17954d99f73df8ed2e6ac"
REVIEWED_HEAD = "1144fe56b68b0980f7f7ce59887987b10e0f3d9a"
CLEANUP_REPAIR_BASE = "f77b73786960f8a4136267372efa2b607f064706"
# Exact Git blob identities of the amended cleanup helper and its tests.
CLEANUP_REPAIR_BLOBS = {
    "testing/run_e8_4_fix5_canonical_five_plot_gate.py": "a05f041ceaf2a1ef5798002272a1e29a22083b83",
    "testing/test_run_e8_4_fix5_canonical_five_plot_gate.py": "48cc342829933da42cf925f453ad19fe352a97c7",
}
# Approved co-published companions at d6a923d; future edits remain frozen.
REVIEWED_PUBLICATION_BLOBS = {
    "testing/pion_local_t_event_validation_replay.py": "8addd03229a95fba8753c539ee966dc4ec74d01b",
    "testing/pion_local_t_statistical_validation.py": "93ddf7b329598485e4a6b5236f64abcc7e8aa063",
    "testing/pion_t_template_comparison_gate_profile.json": "5020d2e71e6ea36e2bce6a8e97b8200c9255132f",
    "testing/render_pion_t_template_comparison_pdf.py": "3039b71a0030e2de71ed15c20a66af1192d81970",
    "testing/test_pion_local_t_statistical_validation.py": "7b0dfda5a6c66b0064159121fd5b79ed122354c3",
    "testing/test_render_pion_t_template_comparison_pdf.py": "c4a64c090833565a575713c41b72ae11a3f26fa8",
}
PROFILE = "testing/pion_t_phi_template_comparison_gate_profile_v2.json"
PROFILE_ID = "pion_t_phi_template_comparison_center_left_lowe/v2"
ENABLE_ENV = "LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON"
CHILD_ENABLE_ENV = "LT_ANALYSIS_PION_T_PHI_TEMPLATE_COMPARISON"
COMMAND = ["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"]
SETTINGS = [dict(row) for row in isolation.CANONICAL_SETTINGS]
COMPONENTS = ("pi_n", "pi_delta", "pi_sidis")
SOURCE_FILES = (
    "src/cuts/pion_t_template_comparison.py", "src/utility/pion_component_shapes.py",
    "src/cuts/rand_sub.py", "src/cuts/pion_component_fits.py", "src/cuts/pion_t_bin_parents.py",
    "src/cuts/apply_cuts.py", "src/utility/background_config.py", "src/utility/canonical_binning.py",
    "src/cuts/hgcer_hole.py", "src/cuts/data_coordinates.py",
    "src/utility/root_histogram_ownership.py", "run_Prod_Analysis.sh")
ARTIFACTS = {"global": [
    {"key": setting.lower() + "_comparison_" + extension,
     "basename_template": f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.{extension}",
     "kind": "json" if extension == "json" else "file", "required": True}
    for setting in ("Center", "Left") for extension in ("json", "pdf")] + [
    {"key": setting.lower() + "_child_comparison_json",
     "basename_template": f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json",
     "kind": "json", "required": True} for setting in ("Center", "Left")], "settings": []}
# Exact consumed scientific/reporting dependencies and their regression tests.
# Unrelated pending statistics/replay implementation remains outside this scope.
INTEGRATION_FILES = {
    "src/cuts/pion_t_phi_template_comparison.py", "src/binning/calculate_yield.py",
    "src/cuts/pion_t_template_comparison.py",
    "testing/test_pion_t_phi_template_comparison.py",
    "testing/test_pion_t_template_comparison.py",
    "testing/test_pion_local_t_event_validation.py",
    "testing/render_pion_t_phi_template_comparison.py",
}
BOUNDARIES = {"diagnostic_only": True, "scientifically_accepted": False,
    "production_applied": False, "ordinary_production_outputs_may_be_written": False,
    "whole_output_root_immutable": False, "private_scientific_outputs": True,
    "pdf_legibility_verified": False,
    "method_a_detached": True, "method_b_numerically_excluded": True}

# Infrastructure entry points, never replacement environment configuration.
# /volatile may resolve through trusted aliases above the physical subtree.
VOLATILE_MOUNT = Path("/volatile")
GROUP_MOUNT = Path("/group")


def directory_identity(path, private=False, owned=True):
    path = Path(path)
    require(path.is_absolute() and path.resolve(strict=True) == path and
            all(not p.is_symlink() for p in (path, *path.parents)) and path.is_dir(),
            "runtime_directory_alias_or_invalid:" + str(path))
    info = path.stat()
    require((not owned or info.st_uid == os.geteuid()) and
            (not private or stat.S_IMODE(info.st_mode) in (0o700, 0o2700)),
            "runtime_directory_owner_or_mode:" + str(path))
    return {"path": str(path), "canonical": str(path.resolve()), "device": info.st_dev, "inode": info.st_ino,
            "uid": info.st_uid, "mode": stat.S_IMODE(info.st_mode)}


def volatile_alias_identity(path):
    """Snapshot infrastructure links, including links in their target paths.

    Only root-controlled parents without group/other write, or read-only
    synthetic infrastructure parents, may supply trusted aliases. User-owned
    writable directories never extend this trust into the configured base.
    """
    path = Path(path)
    canonical = path.resolve(strict=True)
    trusted_uid = VOLATILE_MOUNT.lstat().st_uid
    links, visited, pending = {}, set(), [path]
    while pending:
        candidate = pending.pop()
        if candidate in visited:
            continue
        visited.add(candidate)
        require(len(visited) <= 64, "volatile_alias_chain_invalid")
        for entry in reversed((candidate, *candidate.parents)):
            if not entry.is_symlink():
                continue
            info, parent = entry.lstat(), entry.parent.resolve(strict=True)
            parent_info = parent.stat()
            mode = stat.S_IMODE(parent_info.st_mode)
            require((parent_info.st_uid == 0 and mode & 0o022 == 0) or mode & 0o222 == 0,
                    "volatile_alias_parent_writable:" + str(entry))
            require(info.st_uid == parent_info.st_uid == trusted_uid,
                    "volatile_alias_owner_invalid:" + str(entry))
            target = os.readlink(entry)
            links[str(entry)] = {"target": target, "canonical": str(entry.resolve(strict=True)),
                "device": info.st_dev, "inode": info.st_ino, "uid": info.st_uid,
                "parent": directory_identity(parent, owned=False)}
            pending.append(entry.parent / target)
    require(path.resolve(strict=True) == canonical, "volatile_alias_changed")
    return canonical, links


def volatile_base(repo, baseline):
    raw = baseline.get("VOLATILEPATH", "")
    require(isinstance(raw, str) and raw and Path(raw).is_absolute(), "volatile_base_not_absolute")
    require(VOLATILE_MOUNT.is_absolute(), "volatile_anchor_invalid")
    anchor, aliases = volatile_alias_identity(VOLATILE_MOUNT)
    base, base_aliases = volatile_alias_identity(Path(raw))
    require(all(aliases.get(key) == value for key, value in base_aliases.items()),
            "volatile_untrusted_base_alias")
    anchor_identity = directory_identity(anchor, owned=False)
    identity = directory_identity(base)
    require(base != anchor and base.is_relative_to(anchor) and
            identity["device"] == anchor_identity["device"] and os.access(base, os.W_OK | os.X_OK),
            "volatile_base_invalid")
    for forbidden in (Path(repo).resolve(), GROUP_MOUNT.resolve()):
        require(not base.is_relative_to(forbidden) and not forbidden.is_relative_to(base),
                "volatile_base_overlap")
    if GROUP_MOUNT.exists():
        require(base.stat().st_dev != GROUP_MOUNT.stat().st_dev, "volatile_on_group_device")
    require(base.stat().st_dev != Path(repo).stat().st_dev, "volatile_on_repository_device")
    return {"resolver_VOLATILEPATH": raw, "base": identity,
            "anchor": anchor_identity, "aliases": aliases,
            "mount_device": anchor_identity["device"]}


def verify_runtime_parent(repo, record):
    current = volatile_base(repo, {"VOLATILEPATH": record["resolver_VOLATILEPATH"]})
    require(current == {k: record[k] for k in current}, "volatile_base_identity_changed")
    parent = Path(record["parent"]["path"])
    require(parent.parent == Path(record["base"]["path"]) and
            parent.name.startswith("kaonlt-ab-analysis-") and
            directory_identity(parent, private=True) == record["parent"],
            "runtime_parent_identity_changed")
    return parent


def runtime_directory(repo, record, path, private=False):
    parent = verify_runtime_parent(repo, record)
    path = Path(path)
    require(path != parent and path.is_relative_to(parent), "runtime_path_outside_volatile_parent")
    identity = directory_identity(path, private=private)
    require(identity["device"] == record["mount_device"], "runtime_device_changed")
    return identity


@contextmanager
def private_data_parent(repo, record):
    # No implicit recursive cleanup after a failed first identity check.
    require(shutil.rmtree.avoids_symlink_attacks, "runtime_cleanup_not_symlink_safe")
    require(volatile_base(repo, {"VOLATILEPATH": record["resolver_VOLATILEPATH"]}) == record,
            "volatile_base_identity_changed")
    parent = Path(tempfile.mkdtemp(prefix="kaonlt-ab-analysis-", dir=record["base"]["path"]))
    created = None
    try:
        created = parent.lstat()
        record["parent"] = directory_identity(parent, private=True)
        verify_runtime_parent(repo, record)
        yield parent
    finally:
        if "parent" not in record:
            current = volatile_base(repo, {"VOLATILEPATH": record["resolver_VOLATILEPATH"]})
            require(current == {k: record[k] for k in current}, "volatile_base_identity_changed")
            require(created is not None and stat.S_ISDIR(created.st_mode) and
                    parent.parent == Path(record["base"]["path"]) and
                    parent.name.startswith("kaonlt-ab-analysis-"),
                    "runtime_initial_parent_cleanup_unsafe:" + str(parent))
            identity = directory_identity(parent)
            require((identity["device"], identity["inode"], identity["uid"], identity["mode"]) ==
                    (created.st_dev, created.st_ino, created.st_uid, stat.S_IMODE(created.st_mode)) and
                    created.st_dev == record["mount_device"],
                    "runtime_initial_parent_cleanup_unsafe:" + str(parent))
            parent.rmdir()
        else:
            verify_runtime_parent(repo, record)
            if "code_parent" in record:
                require(not os.path.lexists(record["code_parent"]["path"]), "runtime_code_identity_changed_or_cleanup_incomplete")
            for location, identity in record.get("locations", {}).items():
                require(runtime_directory(repo, record, Path(location)) == identity,
                        "runtime_cleanup_location_changed:" + location)
            require(shutil.rmtree.avoids_symlink_attacks, "runtime_cleanup_not_symlink_safe")
            shutil.rmtree(parent)
        require(not os.path.lexists(parent), "runtime_parent_cleanup_incomplete")
        current = volatile_base(repo, {"VOLATILEPATH": record["resolver_VOLATILEPATH"]})
        require(current == {k: record[k] for k in current}, "volatile_base_identity_changed")


def verify_code_layout(repo, worktree, record):
    isolation.cleanup_path(repo, worktree.parent, worktree, "analysis")
    require(worktree.parent.parent == Path(repo).resolve().parent and
            directory_identity(worktree.parent, private=True) == record["code_parent"] and
            directory_identity(worktree) == record["analysis"] and
            worktree.stat().st_dev == Path(repo).stat().st_dev != record["mount_device"],
            "runtime_code_identity_changed")
    for location, identity in record.get("code_locations", {}).items():
        require(directory_identity(Path(location)) == identity,
                "runtime_code_location_changed:" + location)


@contextmanager
def analysis_worktree(repo, source_commit, record):
    with private_data_parent(repo, record):
        with isolation.owned_worktree(repo, source_commit, "analysis") as worktree:
            record["code_parent"] = directory_identity(worktree.parent, private=True)
            record["analysis"] = directory_identity(worktree)
            verify_code_layout(repo, worktree, record)
            try:
                yield worktree
            finally:
                verify_code_layout(repo, worktree, record)
                verify_runtime_parent(repo, record)
                if "locations" in record:
                    verify_private(worktree, Path(record["parent"]["path"]) / "ab-output/Analysis/KaonLT")
        require(not os.path.lexists(worktree.parent), "runtime_code_cleanup_incomplete")


def verify_runtime_layout(repo, worktree, record, output, env):
    verify_code_layout(repo, worktree, record)
    parent = verify_runtime_parent(repo, record)
    require(output == parent / "ab-output/Analysis/KaonLT", "private_output_path_changed")
    verify_private(worktree, output)
    code_paths = [worktree.parent / "runtime-python", worktree.parent / "runtime-python/ltsep"]
    code = {str(p): directory_identity(p) for p in code_paths}
    require(all(identity["device"] == record["analysis"]["device"] for identity in code.values()),
            "runtime_code_device_changed")
    if "code_locations" in record:
        require(code == record["code_locations"], "runtime_code_location_changed")
    record["code_locations"] = code
    locations = [parent / "ab-output", output.parent, output]
    locations += [Path(env[k]) for k in ("TMPDIR", "XDG_CACHE_HOME", "MPLCONFIGDIR", "PYTHONPYCACHEPREFIX")]
    identities = {str(p): runtime_directory(repo, record, p, private=p == locations[0] or p in locations[3:])
                  for p in locations}
    require(not output.is_relative_to(Path(record["shared_output"])) and
            not Path(record["shared_output"]).is_relative_to(output), "runtime_shared_output_overlap")
    if "locations" in record:
        require(identities == record["locations"], "runtime_location_identity_changed")
    record["locations"] = identities
    return identities


def prepare_runtime_overlay(repo, worktree, placement):
    """Use the established installed-package algorithm in its group geometry."""
    verify_code_layout(repo, worktree, placement)
    env, overlay = isolation.prepare_runtime_overlay(repo, worktree)
    require(overlay["baseline_paths"] == placement["baseline_paths"], "runtime_baseline_changed")
    verify_code_layout(repo, worktree, placement)
    return env, overlay


def strict_json(raw):
    def invalid(value):
        raise ValueError("nonfinite_json:" + value)
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate_json_key:" + key)
            result[key] = value
        return result
    result = json.loads(raw, parse_constant=invalid, object_pairs_hook=pairs)
    def finite(value):
        if isinstance(value, float):
            require(math.isfinite(value), "nonfinite_json_number")
        elif isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(result)
    return result


def resolved_profile(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None, "source_commit_requires_full_sha")
    profile = collector.load_validation_profile(Path(repo) / PROFILE)
    expected = {"schema_version": "pion_hgcer_validation_bundle_profile/v4",
        "validation_profile": PROFILE_ID, "collection_mode": "generic_artifacts",
        "settings": SETTINGS, "artifacts": ARTIFACTS,
        "source_identity": {"required_analysis_commit": REVIEWED_HEAD, "allowed_committed_files": [],
                            "allowed_non_analysis_path_prefixes": ["docs/memory/"]}}
    require(profile == expected, "profile_contract_mismatch")
    profile = deepcopy(profile)
    profile["source_identity"]["required_analysis_commit"] = source_commit
    return profile


def names():
    return {row["basename_template"] for row in ARTIFACTS["global"]}


def preflight(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None, "source_commit_requires_full_sha")
    before = isolation.snapshot(repo)
    require(before["branch"] == "test", "wrong_branch")
    require(before["head"] == source_commit, "wrong_head")
    require(before["origin_test"] == source_commit, "source_not_observed_pushed")
    require(source_commit != BASE_HEAD, "owner_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, source_commit)
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, REVIEWED_HEAD)
    git(repo, "merge-base", "--is-ancestor", REVIEWED_HEAD, source_commit)
    for line in before["porcelain"].splitlines():
        path = line[3:]
        require(path in isolation.FARM_OUTPUTS or path.startswith("OUTPUT/") or
                (line[:2] == "??" and re.fullmatch(r"kaonlt_review(?:\([^/]+\))?\.diff", path)),
                "dirty_gate_source:" + line)
    # Historical ancestry and cleanup pins stay strict; only exact reviewed paths.
    changed = git(repo, "diff", "--name-only", REVIEWED_HEAD, source_commit).splitlines()
    allowed = {PROFILE, "testing/run_pion_t_template_comparison_gate.py",
               "testing/test_run_pion_t_template_comparison_gate.py", *CLEANUP_REPAIR_BLOBS,
               *INTEGRATION_FILES, *REVIEWED_PUBLICATION_BLOBS}
    require(all(p in allowed or p.startswith("docs/memory/") for p in changed),
            "unexpected_committed_source")
    git(repo, "merge-base", "--is-ancestor", CLEANUP_REPAIR_BASE, source_commit)
    for path, blob in CLEANUP_REPAIR_BLOBS.items():
        require(git(repo, "rev-parse", source_commit + ":" + path) == blob,
                "cleanup_repair_blob_changed:" + path)
    for path, blob in REVIEWED_PUBLICATION_BLOBS.items():
        require(git(repo, "rev-parse", source_commit + ":" + path) == blob,
                "reviewed_publication_blob_changed:" + path)
    return before


def child_source_inventory(repo):
    """Read the observer's declared inventory without importing physics/ROOT."""
    base_tree = ast.parse((Path(repo) / "src/cuts/pion_t_template_comparison.py").read_text())
    base = [n.value for n in base_tree.body if isinstance(n, ast.Assign) and
            any(isinstance(t, ast.Name) and t.id == "REQUIRED_SOURCE_FILES" for t in n.targets)]
    require(len(base) == 1 and ast.literal_eval(base[0]) == SOURCE_FILES,
            "observer_base_source_inventory_changed")
    tree = ast.parse((Path(repo) / "src/cuts/pion_t_phi_template_comparison.py").read_text())
    declarations = [n.value for n in ast.walk(tree) if isinstance(n, ast.Assign) and
                    any(isinstance(t, ast.Name) and t.id == "files" for t in n.targets)]
    require(len(declarations) == 1 and isinstance(declarations[0], ast.Tuple),
            "observer_child_source_inventory_invalid")
    elements = declarations[0].elts
    require(isinstance(elements[0], ast.Starred) and isinstance(elements[0].value, ast.Name)
            and elements[0].value.id == "REQUIRED_SOURCE_FILES", "observer_child_source_inventory_invalid")
    paths = SOURCE_FILES + tuple(ast.literal_eval(n) for n in elements[1:])
    require(len(set(paths)) == len(paths) and all(isinstance(p, str) and
            not Path(p).is_absolute() and ".." not in Path(p).parts for p in paths),
            "observer_child_source_inventory_invalid")
    return paths


def committed_identities(repo, source_commit):
    records = {}
    for path in child_source_inventory(repo):
        result = subprocess.run(["git", "show", source_commit + ":" + path], cwd=repo,
                                capture_output=True)
        require(result.returncode == 0 and result.stdout, "committed_source_missing:" + path)
        records[path] = {"sha256": hashlib.sha256(result.stdout).hexdigest(), "bytes": len(result.stdout)}
    return records


def source_preflight(repo, source_commit, profile):
    with isolation.owned_worktree(repo, source_commit, "source") as source:
        require(resolved_profile(source, source_commit) == profile, "detached_profile_mismatch")
        with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
            _, checks = module.collect_source_checks(source,
                required_analysis_commit=source_commit, allowed_committed_files=[])
            require(all(c["returncode"] == 0 for c in checks), "collector_source_check_failed")
            ancestor, _, unexpected = module._committed_identity(checks, [], ["docs/memory/"])
            require(ancestor and not unexpected, "collector_source_identity_invalid")
        identities = committed_identities(source, source_commit)
    return checks, identities


def child_environment(env):
    result = dict(env)
    for key in ("LT_ANALYSIS_DEBUG_LEFT_LOW", "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE",
                "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING"):
        result.pop(key, None)
    flags = (ENABLE_ENV, CHILD_ENABLE_ENV)
    require(not any(k.startswith("LT_ANALYSIS_") and k not in flags and v
                    for k, v in result.items()), "untrusted_analysis_environment")
    for flag in flags:
        require(result.get(flag, "1") == "1", "conflicting_comparison_environment")
        result[flag] = "1"
    return result


def verify_completion(log_path, returncode):
    require(returncode == 0, "analysis_failed")
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    require("All available low-epsilon settings completed in debug mode." in text and
            "Full high-epsilon processing is intentionally skipped in -d debug mode." in text,
            "all_low_completion_missing")
    require("High Epsilon Completed!" not in text, "unexpected_high_epsilon_execution")


def pdf_pages(path):
    # Existing Poppler backend; do not install packages or use regex page counts.
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    require(result.returncode == 0, "pdf_structure_unreadable:" + result.stderr)
    pages = re.findall(r"(?m)^Pages:\s+(\d+)\s*$", result.stdout)
    require(len(pages) == 1, "pdf_page_count_unavailable")
    return int(pages[0])


def fingerprint(hist):
    require(isinstance(hist, dict) and hist.get("status") == "available" and
            isinstance(hist.get("stored_sumw2"), list) and hist["stored_sumw2"] and
            re.fullmatch(r"[0-9a-f]{64}", hist.get("content_error_fingerprint", "")) is not None,
            "frozen_histogram_identity_invalid")
    return hist["content_error_fingerprint"]


def verify_payload(payload, setting, source_commit, identities, pdf_record):
    require(payload.get("schema_version") == "pion_t_template_comparison/v1" and
            payload.get("diagnostic_only") is True and payload.get("production_applied") is False and
            payload.get("scientifically_accepted") is False, "comparison_schema_or_authority_invalid")
    require(payload.get("source_HEAD") == source_commit and
            payload.get("source_identity_rejection_reasons") == [], "comparison_source_identity_invalid")
    require(payload.get("phi_setting") == setting and payload.get("setting") in (
        {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "ParticleType": "kaon"},
        {"Q2": "4p4", "W": "2p74", "EPSSET": "lowe", "ParticleType": "kaon"}) and
        payload.get("enable_environment_variable") == ENABLE_ENV, "comparison_setting_invalid")
    source = payload.get("source_files") or {}
    require(set(source) == set(SOURCE_FILES), "comparison_source_inventory_invalid")
    for path, identity in source.items():
        require(identity.get("status") == "available" and
                all(identity.get(k) == v for k, v in identities[path].items()),
                "committed_source_hash_mismatch:" + path)
    pdf = payload.get("pdf") or {}
    require(pdf.get("status") == "available" and pdf.get("sha256") == pdf_record["sha256"] and
            pdf.get("bytes") == pdf_record["bytes"] and
            Path(pdf.get("path", "")).name == pdf_record["basename"], "comparison_pdf_identity_invalid")
    edges = payload.get("canonical_t_edges")
    require(isinstance(edges, list) and len(edges) == 4 and
            all(type(x) in (int, float) and math.isfinite(x) for x in edges) and
            all(a < b for a, b in zip(edges, edges[1:])), "canonical_t_edges_invalid")
    parents = payload.get("parents")
    require(isinstance(parents, list) and len(parents) == 3 and
            sorted(p.get("t_bin_index", -1) for p in parents) == [0, 1, 2], "comparison_t_inventory_invalid")
    for parent in parents:
        ti = parent["t_bin_index"]
        require(type(ti) is int and parent.get("t_edges") == edges[ti:ti+2], "parent_t_edges_invalid")
        require(parent.get("authoritative_parent_unchanged") is True and
                parent.get("production_applied") is False and parent.get("scientifically_accepted") is False and
                parent.get("alignment_mode") == "fixed_parent_alignment_template_shape_test" and
                re.fullmatch(r"[0-9a-f]{64}", parent.get("baseline_parent_fingerprint", "")) is not None and
                all(parent.get(k) for k in ("pion_parent_id", "canonical_interval_pair_id",
                                           "canonical_interval_pair_hash", "coordinate_fingerprint")),
                "parent_authority_or_identity_invalid")
        baseline, candidate = parent.get("baseline_A"), parent.get("candidate_B")
        require(isinstance(baseline, dict) and baseline, "baseline_A_missing")
        for side, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
            frozen = fingerprint((parent.get("frozen_inputs") or {}).get(side))
            require(frozen == fingerprint(baseline.get(key)), "baseline_input_identity_mismatch")
            if candidate is not None:
                require(frozen == fingerprint(candidate.get(key)), "candidate_input_identity_mismatch")
        reasons = parent.get("rejection_reasons")
        require(isinstance(reasons, list) and all(isinstance(r, str) and r for r in reasons),
                "candidate_rejection_reasons_invalid")
        require(set(parent.get("templates") or {}) == set(COMPONENTS), "component_inventory_invalid")
        for component, evidence in parent["templates"].items():
            require(isinstance(evidence, dict) and isinstance(evidence.get("provenance"), dict) and
                    evidence["provenance"], "component_provenance_missing:" + component)
            provenance = evidence["provenance"]
            if provenance.get("status") == "available":
                require(re.fullmatch(r"[0-9a-f]{64}", provenance.get("sha256", "")) is not None and
                        provenance.get("path") and provenance.get("tree_name") and
                        provenance.get("weight_expression") == "getattr(evt, 'iter_weight', 1.0)" and
                        provenance.get("coordinate") == "missmass; -t; phipq*180/pi",
                        "component_provenance_invalid:" + component)
            else:
                require(provenance.get("status") == "unavailable" and provenance.get("reason") and
                        parent.get("status") == "unavailable" and reasons,
                        "component_provenance_invalid:" + component)
            if parent.get("status") == "available":
                require(evidence.get("status") == "genuinely_local" and
                        evidence.get("normalization_count") == 1 and
                        math.isclose(evidence.get("normalized_integral", 0), 1., abs_tol=1e-9) and
                        all(evidence.get(k) is True for k in ("canonical_selection_closure",
                            "weighted_child_to_parent_closure", "production_global_selection_closure")) and
                        evidence.get("t_edges") == edges[ti:ti+2], "component_normalization_invalid")
        if parent.get("status") == "available":
            require(isinstance(candidate, dict) and candidate and not reasons, "available_candidate_invalid")
            require(candidate.get("pion_component_alignment") == baseline.get("pion_component_alignment") and
                    candidate.get("resolved_subtraction_config") == baseline.get("resolved_subtraction_config"),
                    "candidate_alignment_or_configuration_changed")
        else:
            require(parent.get("status") == "unavailable" and candidate is None and reasons,
                    "unexplained_unavailable_candidate")
    expected = "available" if all(p["status"] == "available" for p in parents) else "partial_unavailable"
    require(payload.get("status") == expected, "comparison_status_invalid")
    return {"status": expected, "rejection_reasons_by_t": {
        str(p["t_bin_index"]): p["rejection_reasons"] for p in parents}, "canonical_t_edges": edges}


def verify_child_applications(cell, pion_observations):
    """Validate captured transfers, not fit physics or fitted-parameter errors."""
    applications = cell.get("applications")
    require(isinstance(applications, dict) and
            set(applications) == {"Global", "Local_t", "Local_t_phi"}, "child_application_inventory")
    ti, pi = cell["t_index"], cell["phi_index"]
    edges = pion_observations["nommcuts"]["edges"]
    size = len(edges) + 1  # Includes ROOT underflow/overflow.

    def vector(value, nonnegative=False):
        return (isinstance(value, list) and len(value) == size and all(
            type(v) in (int, float) and math.isfinite(v) and (not nonnegative or v >= 0)
            for v in value))

    for method, application in applications.items():
        require(isinstance(application, dict) and application.get("status") in ("available", "unavailable"),
                "child_application_status:" + method)
        available = application["status"] == "available"
        if method == "Local_t_phi" and not cell["fit_accepted"]:
            require(not available, "child_rejected_weights_applied")
        if not available:
            require(isinstance(application.get("reason"), str) and application["reason"].strip() and
                    not any(k in application for k in ("weights", "cuts")),
                    "child_application_unavailable_evidence:" + method)
            continue
        weights = application.get("weights")
        require(vector(weights) and weights[0] == weights[-1] == 0,
                "child_application_weights:" + method)
        cuts = application.get("cuts")
        require(isinstance(cuts, dict) and set(cuts) == {"allcuts", "nommcuts"},
                "child_application_cuts:" + method)
        for cut, spectrum in cuts.items():
            require(pion_observations[cut]["edges"] == edges, "child_application_axis:" + method)
            require(isinstance(spectrum, dict) and vector(spectrum.get("values")) and
                    vector(spectrum.get("variances"), nonnegative=True),
                    "child_application_spectrum:" + method + ":" + cut)
            observed = pion_observations[cut]["children"][ti][pi]
            for field, expected in (
                ("values", [v*w for v,w in zip(observed["values"], weights)]),
                ("variances", [v*w*w for v,w in zip(observed["variances"], weights)]),
            ):
                require(vector(expected) and all(math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-12)
                        for a,b in zip(spectrum[field], expected)),
                        "child_application_reconstruction:" + method + ":" + cut + ":" + field)


def verify_child_payload(payload, setting, source_commit, identities, ab, runtime_inputs=None):
    """Structural/provenance closure, independent of ROOT and scientific adequacy."""
    require(payload.get("schema_version") == "pion_t_phi_template_comparison/v1" and
            payload.get("status") == "captured_pending_independent_validation" and
            payload.get("diagnostic_only") is True and payload.get("production_applied") is False and
            payload.get("scientifically_accepted") is False, "child_schema_or_authority_invalid")
    require(payload.get("source_HEAD") == source_commit and payload.get("source_identity_complete") is True,
            "child_source_identity_invalid")
    sources = payload.get("source_files") or {}
    require(set(sources) == set(identities), "child_source_inventory_invalid")
    for path, identity in sources.items():
        require(identity.get("status") == "available" and all(identity.get(k) == v
                for k, v in identities[path].items()), "child_committed_source_hash_mismatch:" + path)
    require(payload.get("phi_setting") == setting and payload.get("t_edges") == ab["canonical_t_edges"] and
            payload.get("phi_edges") == ab.get("canonical_phi_edges"), "child_setting_or_edges_invalid")
    for edges, size in ((payload["t_edges"], 4), (payload["phi_edges"], 10)):
        require(isinstance(edges, list) and len(edges) == size and
                all(type(v) in (int, float) and math.isfinite(v) for v in edges) and
                all(a < b for a, b in zip(edges, edges[1:])), "child_edges_invalid")
    require(payload.get("parent_input_closure") == "PASS" and
            payload.get("production_parents_unchanged") is True, "child_parent_closure_invalid")
    inputs = payload.get("runtime_inputs") or {}
    require(set(inputs) == {"data", "dummy"} and
            set(ab.get("original_analysis_input_identity") or {}) == set(inputs), "child_runtime_input_inventory")
    for role, identity in inputs.items():
        original = ab["original_analysis_input_identity"][role]
        require(identity.get("status") == original.get("status") == "available" and
                all(identity.get(k) == original.get(k) for k in ("path", "sha256", "bytes")) and
                type(identity.get("bytes")) is int and identity["bytes"] > 0 and
                re.fullmatch(r"[0-9a-f]{64}", identity.get("sha256", "")) is not None,
                "child_runtime_input_mismatch:" + role)
        if runtime_inputs is not None:
            name = f"{setting}_kaon_Analysed_{role.title()}_Q4p4W2p74_lowe.root"
            actual = runtime_inputs.get(name) or {}
            require(str(Path(identity["path"]).resolve()) == actual.get("resolved") and
                    all(identity.get(k) == actual.get(k) for k in ("sha256", "bytes")),
                    "child_private_container_mismatch:" + role)

    def vector_equal(a, b, reason):
        require(isinstance(a, list) and isinstance(b, list) and len(a) == len(b) and
                all(type(x) in (float, int) and type(y) in (float, int) and
                    math.isfinite(x) and math.isfinite(y) and math.isclose(x, y, rel_tol=1e-8, abs_tol=1e-12)
                    for x, y in zip(a, b)), reason)

    observations = payload.get("observations") or {}
    require(set(observations) == {"pion", "kaon"}, "child_population_inventory")
    for pop, cuts in observations.items():
        require(set(cuts) == {"allcuts", "nommcuts"}, "child_cut_inventory")
        for cut, projection in cuts.items():
            edges = projection.get("edges")
            require(isinstance(edges, list) and len(edges) >= 2 and
                    all(type(v) in (int, float) and math.isfinite(v) for v in edges) and
                    all(a < b for a,b in zip(edges,edges[1:])), "child_mm_edges_invalid")
            require(projection.get("closure") == "PASS" and all(len(projection.get(k, [])) == 3
                    for k in ("parents", "children", "sources", "outside_phi")), "child_projection_inventory")
            for ti in range(3):
                children, source = projection["children"][ti], projection["sources"][ti]
                require(len(children) == 9 and set(source) == {"prompt", "rand", "dummy", "dummy_rand"},
                        "child_projection_inventory")
                parent = projection["parents"][ti]
                for cell in [parent, projection["outside_phi"][ti], *children, *source.values()]:
                    require(type(cell.get("count")) is int and cell["count"] >= 0, "child_count_invalid")
                    for field in ("values", "variances"):
                        v = cell.get(field)
                        require(isinstance(v, list) and len(v) == len(edges)+1 and all(
                            type(x) in (int, float) and math.isfinite(x) and (field != "variances" or x >= 0)
                            for x in v), "child_signed_sumw2_invalid")
                for group in ([*children, projection["outside_phi"][ti]], list(source.values())):
                    require(parent["count"] == sum(c["count"] for c in group), "child_count_closure")
                    for field in ("values", "variances"):
                        vector_equal(parent[field], [math.fsum(c[field][i] for c in group)
                            for i in range(len(edges)+1)], "child_signed_sumw2_closure")
                if cut == "nommcuts":
                    frozen = ab["parents"][ti]["frozen_inputs"][pop]
                    require(edges == frozen.get("edges"), "child_parent_mm_edges_mismatch")
                    vector_equal(parent["values"], frozen.get("values"), "child_parent_content_mismatch")
                    vector_equal(parent["variances"], frozen.get("stored_sumw2"), "child_parent_sumw2_mismatch")
    cells = payload.get("cells")
    require(isinstance(cells, list) and len(cells) == 27, "child_cell_inventory")
    keys = [(c.get("t_index"), c.get("phi_index")) for c in cells]
    require(all(type(t) is int and type(p) is int for t,p in keys) and
            Counter(keys) == Counter((t,p) for t in range(3) for p in range(9)), "child_cell_coordinates")
    require(payload.get("cell_status_counts") == dict(Counter(c.get("status") for c in cells)),
            "child_status_accounting")
    scopes = {"local_t_phi", "local_t_aggregate_fallback", "setting_fallback", "unavailable"}
    for cell in cells:
        ti, pi = cell["t_index"], cell["phi_index"]
        require(cell.get("t_edges") == payload["t_edges"][ti:ti+2] and
                cell.get("phi_edges") == payload["phi_edges"][pi:pi+2], "child_cell_edges")
        for field, components in (("template_scopes", set(COMPONENTS)),
                                  ("signal_scopes", {"k_lambda", "k_sigma0"})):
            mapping = cell.get(field) or {}
            require(set(mapping) == components and all(v in scopes for v in mapping.values()), "child_shape_scope")
        for pop in ("pion", "kaon"):
            h = cell["inputs"][pop]
            observed = observations[pop]["nommcuts"]
            require(h.get("edges") == observed["edges"], "child_input_mm_edges")
            fingerprint(h)
            vector_equal(h.get("values"), observed["children"][ti][pi]["values"], "child_input_content")
            vector_equal(h.get("stored_sumw2"), observed["children"][ti][pi]["variances"], "child_input_sumw2")
        returned, accepted = cell.get("fit_returned"), cell.get("fit_accepted")
        require(type(returned) is bool and type(accepted) is bool, "child_acceptance_missing")
        fit, evaluation = cell.get("fit"), cell.get("fit_acceptance")
        if returned:
            require(isinstance(fit, dict) and isinstance(evaluation, dict) and
                    evaluation.get("accepted") is accepted and isinstance(evaluation.get("reason"), str) and
                    isinstance(evaluation.get("diagnostics"), dict), "child_fit_gate_evidence")
            require(fit.get("analysis_scope") == f"t_bin{ti+1}phi_bin{pi+1}" and
                    all(isinstance(fit.get(k), str) and fit[k] for k in ("fit_status_pion", "fit_status_kaon")),
                    "child_protected_fit_status_missing")
            protected = ((fit.get("diagnostics") or {}).get("kaon") or {}).get("pi_delta_signal_protected_fit")
            require(cell.get("protected_model_status") == (protected or
                    {"status":"unavailable","reason":"protected_diagnostics_not_retained"}) and
                    isinstance(cell["protected_model_status"].get("status"), str) and
                    cell["protected_model_status"]["status"], "child_protected_model_status")
            for pop, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
                require(fingerprint(fit.get(key)) == fingerprint(cell["inputs"][pop]), "child_fit_input_identity")
                require(fit[key].get("edges") == cell["inputs"][pop]["edges"], "child_fit_input_edges")
                for field in ("values", "stored_sumw2"):
                    vector_equal(fit[key].get(field), cell["inputs"][pop].get(field), "child_fit_input_content")
            if accepted:
                diagnostics = evaluation["diagnostics"]
                require(not evaluation["reason"] and not fit.get("fallback_used") and
                        all(fit[k].strip().lower() == "success" for k in ("fit_status_pion", "fit_status_kaon")) and
                        all(diagnostics.get("fit_validation_" + pop) is True and
                            (((fit.get("diagnostics") or {}).get(pop) or {}).get("validation") or {}).get("accepted") is True
                            for pop in ("pion", "kaon")), "child_false_fit_acceptance")
                for key in ("A_n", "A_delta", "A_sidis", "B_n", "B_delta", "B_sidis"):
                    value = fit.get(key)
                    require(type(value) in (int, float) and math.isfinite(value) and value >= 0 and
                            diagnostics.get("amplitudes", {}).get(key) == value, "child_gate_amplitudes")
                expected = "fitted_local_templates" if all(v == "local_t_phi" for v in
                    cell["template_scopes"].values()) else "fitted_with_template_fallback"
            else:
                require(evaluation["reason"] and cell.get("reason") == evaluation["reason"], "child_rejection_reason")
                expected = "fit_rejected"
        else:
            # snapshot(None) is persisted as this exact explicit missing record.
            # Presence is required: omitted/malformed evidence must still fail.
            missing = {"status": "unavailable", "reason": "not retained by producer"}
            require(not accepted and all(key in cell and
                    (cell[key] is None or cell[key] == missing)
                    for key in ("fit", "fit_acceptance")) and
                    isinstance(cell.get("reason"), str) and cell["reason"].strip(),
                    "child_failed_fit_evidence")
            require(cell.get("protected_model_status") == {"status":"unavailable","reason":"fit_not_returned"},
                    "child_protected_model_status")
            expected = "empty_child" if cell["reason"] == "empty_observed_cell" else "fit_failed"
            if expected == "empty_child":
                require(not any(v > 0 for pop in ("pion", "kaon") for v in
                        cell["inputs"][pop]["stored_sumw2"]), "child_false_empty_status")
        require(cell.get("status") == expected, "child_fit_status_invalid")
        verify_child_applications(cell, observations["pion"])
    return {"collection_status": "verified", "cell_status_counts": payload["cell_status_counts"],
            "accepted_fits": sum(c["fit_accepted"] for c in cells), "physics_adequacy": "NOT VERIFIED"}


def verify_artifacts(outdir, source_commit, identities, started_ns, before, runtime_inputs=None):
    records, outcomes = {}, {}
    child_names = {n for n in names() if "pion-t-phi-" in n}
    require({p.name for p in outdir.glob("*pion-t-phi-template-comparison*")} == child_names,
            "child_artifact_inventory_invalid")
    for name in sorted(names()):
        path = outdir / name
        require(path.is_file() and not path.is_symlink() and path.stat().st_size > 0,
                "comparison_artifact_missing_or_empty:" + name)
        info = path.stat()
        require(info.st_mtime_ns >= started_ns and before[name] is None,
                "comparison_artifact_stale:" + name)
        records[name] = {"basename": name, "bytes": info.st_size, "sha256": sha256(path),
                         "mtime_ns": info.st_mtime_ns}
        if path.suffix == ".pdf":
            require(path.read_bytes().startswith(b"%PDF-"), "comparison_pdf_signature_invalid")
            require(pdf_pages(path) == 6, "comparison_pdf_page_count_invalid")
    for setting in ("Center", "Left"):
        stem = f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison"
        payload = strict_json((outdir / (stem + ".json")).read_text(encoding="utf-8"))
        outcomes[setting] = verify_payload(payload, setting, source_commit, identities, records[stem + ".pdf"])
        child = strict_json((outdir / (stem.replace("pion-t-", "pion-t-phi-") + ".json")).read_text(encoding="utf-8"))
        outcomes[setting]["child_comparison"] = verify_child_payload(
            child, setting, source_commit, identities, payload, runtime_inputs)
    require(outcomes["Center"]["canonical_t_edges"] == outcomes["Left"]["canonical_t_edges"],
            "mixed_canonical_t_edges")
    return records, outcomes


def verify_zip(path, source_commit, records):
    expected = {"manifest.json", "source_state.txt", "source_checks.txt", "global/"}
    expected.update(s["phi"] + "_" + s["epsilon"] + "/" for s in SETTINGS)
    expected.update("global/" + name for name in names())
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_crc_failed")
        require(Counter(archive.namelist()) == Counter(expected), "zip_inventory_invalid")
        manifest = strict_json(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [] and
                manifest.get("git_head") == manifest.get("required_analysis_commit") == source_commit and
                manifest.get("validation_profile") == PROFILE_ID and
                manifest.get("requested_kinematic") == "Q4p4W2p74" and
                manifest.get("requested_settings") == SETTINGS, "zip_manifest_invalid")
        require(len(manifest.get("settings", [])) == 5 and
                all(s.get("artifacts") == {} for s in manifest["settings"]), "zip_setting_artifacts_invalid")
        require(set(manifest.get("global_artifacts", {})) == {r["key"] for r in ARTIFACTS["global"]},
                "zip_global_inventory_invalid")
        for declaration in ARTIFACTS["global"]:
            name = declaration["basename_template"]
            raw = archive.read("global/" + name)
            record = manifest["global_artifacts"][declaration["key"]]
            require(len(raw) == records[name]["bytes"] and
                    hashlib.sha256(raw).hexdigest() == records[name]["sha256"] and
                    record.get("sha256") == records[name]["sha256"] and
                    record.get("byte_size") == records[name]["bytes"], "zip_artifact_identity_mismatch:" + name)


def destinations(outdir, output):
    return {"log": outdir / (output.stem + ".log"),
            "summary": outdir / (output.stem + "-run-summary.json"),
            "receipt": isolation.accepted.gate_status_path(outdir, output),
            "zip": output,
            "delivered_log": output.with_name(output.stem + ".log"),
            "delivered_summary": output.with_name(output.stem + "-run-summary.json"),
            "delivered_receipt": output.with_name(output.stem + "-gate-status.json")}


def file_identity(path):
    path = Path(path)
    info = path.stat()
    require(path.is_file(), "input_not_regular:" + str(path))
    return {"sha256": sha256(path), "bytes": info.st_size,
            "device": info.st_dev, "inode": info.st_ino, "mtime_ns": info.st_mtime_ns,
            "resolved": str(path.resolve()), "link": os.readlink(path) if path.is_symlink() else None}


def accepted_identities(outdir):
    records = {}
    for name in sorted(names()):
        p = Path(outdir) / name
        if not os.path.lexists(p):
            records[name] = None
            continue
        link = os.readlink(p) if p.is_symlink() else None
        require(p.is_file(), "accepted_artifact_not_regular_or_dangling:" + name)
        records[name] = {"identity": file_identity(p), "link": link, "resolved": str(p.resolve())}
    return records


def verify_accepted(outdir, before):
    require(accepted_identities(outdir) == before, "accepted_artifact_changed")
    return {"passed": True, "artifacts": before}


def private_output(worktree, paths, shared, data_parent):
    """Exclusive owner-local topology; the frozen shared-link helper is unused."""
    worktree, shared = Path(worktree), Path(shared).resolve()
    external_output = (Path(paths["VOLATILEPATH"]) / "OUTPUT").resolve()
    require(shared == (external_output / "Analysis/KaonLT").resolve(), "private_shared_authority")
    parent = Path(data_parent)
    private, link = parent / "ab-output", worktree / "OUTPUT"
    require(worktree.is_absolute() and worktree.resolve() == worktree and worktree.name == "analysis"
            and not worktree.is_symlink() and worktree.is_dir(), "private_worktree_invalid")
    require(parent.resolve() == parent and parent.stat().st_uid == os.geteuid()
            and stat.S_IMODE(parent.stat().st_mode) & 0o077 == 0, "private_parent_permissions")
    require(paths["ANATYPE"] == "Kaon" and Path(paths["LTANAPATH"]) == worktree
            and Path(paths["OUTPATH"]) == link / "Analysis/KaonLT", "private_path_authority")
    require(all(not private.is_relative_to(p) and not p.is_relative_to(private)
                for p in (external_output, shared)), "private_shared_overlap")
    require(not os.path.lexists(private) and not os.path.lexists(link), "private_output_collision")
    private.mkdir(mode=0o700)
    target = private / "Analysis/KaonLT"
    target.mkdir(parents=True, mode=0o700)
    link.symlink_to(private, target_is_directory=True)
    require(os.readlink(link) == str(private) and link.resolve() == private
            and Path(paths["OUTPATH"]).resolve() == target, "private_link_identity")
    return target


def verify_private(worktree, output):
    private = Path(output).parents[1]
    link = worktree / "OUTPUT"
    require(link.is_symlink() and os.readlink(link) == str(private) and link.resolve() == private,
            "private_link_changed")
    require(output == private / "Analysis/KaonLT" and all(p.is_dir() and not p.is_symlink()
            for p in (private, private / "Analysis", output)), "private_output_path_changed")
    return {"passed": True}


def skim_directories(repo, paths):
    """Only the existing pure normalizers; no ROOT or setup-module imports."""
    tree = ast.parse((Path(repo) / "src/setup/findEffectiveCharge.py").read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "resolve_skim_outpath"]
    require(len(nodes) == 1, "skim_normalizer_changed")
    ns = {"os": os}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "<skim-normalizer>", "exec"), ns)
    effective = ns["resolve_skim_outpath"](paths["SKIMPATH"], paths["ANATYPE"])
    text = (Path(repo) / "run_Prod_Analysis.sh").read_text()
    require(text.count("normalize_ltsep_dir() {") == 1, "launcher_skim_normalizer_changed")
    body = "normalize_ltsep_dir() {" + text.split("normalize_ltsep_dir() {", 1)[1].split("\n}", 1)[0] + "\n}\n"
    result = subprocess.run(["bash", "-c", body + '\nANATYPE="$1"\nnormalize_ltsep_dir "$2"',
                             "ab-skim-path", paths["ANATYPE"], paths["SKIMPATH"]], capture_output=True, text=True)
    require(result.returncode == 0 and result.stdout.strip() and not result.stderr, "launcher_skim_path_failed")
    return result.stdout.strip(), effective


def external_setup(worktree, paths, shared):
    external = isolation.external_symlink_preflight(worktree, paths)
    if external["background_simc_path"]:
        require(Path(external["background_simc_path"]).is_dir(), "background_simc_directory_missing")
    require((Path(shared) / "background_samples").is_dir(), "external_background_mkdir_required")
    skim, effective = skim_directories(worktree, paths)
    require(Path(skim).is_dir() and Path(effective).is_dir(), "external_skim_directory_missing")
    external["skim_directories"] = [skim, effective]
    external["directories"] = [{"path": str(p), "resolved": str(p.resolve()),
                                "inode": p.stat().st_ino, "device": p.stat().st_dev}
                               for p in (Path(shared) / "background_samples", Path(skim), Path(effective))]
    return external


def verify_external(record):
    for item in record["links"]:
        p = Path(item["path"])
        require(p.is_symlink() and p.exists() and os.readlink(p) == item["literal_target"]
                and str(p.resolve()) == item["resolved_target"], "external_link_changed:" + str(p))
    for item in record.get("directories", []):
        p = Path(item["path"])
        require(p.is_dir() and str(p.resolve()) == item["resolved"] and
                p.stat().st_ino == item["inode"] and p.stat().st_dev == item["device"], "external_directory_changed")
    return {"passed": True}


def reached_callers(repo):
    """Bounded source-only import walk for this unchanged debug entry chain."""
    repo = Path(repo)
    modules = {}
    for p in (repo / "src").rglob("*.py"):
        modules.setdefault(p.stem, []).append(p)
    pending = [(repo / p, cwd) for p, cwd in (("src/main.py", "src"),
        ("src/setup/findEffectiveCharge.py", "src/setup"), ("src/setup/calcTotalEffectiveCharge.py", "src/setup"),
        ("src/setup/set_sig_fortran.py", "src/setup"))]
    seen, callers = set(), set()
    while pending:
        p, cwd = pending.pop()
        if (p, cwd) in seen:
            continue
        seen.add((p, cwd))
        for n in ast.walk(ast.parse(p.read_text())):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in {"Root", "create_ltsep_root"}:
                require(n.args and not n.keywords and all(isinstance(m, ast.Constant) and isinstance(m.value, str)
                        for m in n.args[1:]), "resolver_call_changed:" + str(p))
                modes = tuple(m.value for m in n.args[1:])
                require(modes in ((), ("Plot_LTSep",)), "resolver_mode_changed:" + str(p))
                callers.add((p.relative_to(repo).as_posix(), modes, cwd))
            imports = ([n.module] if isinstance(n, ast.ImportFrom) and n.module else
                       [a.name for a in n.names] if isinstance(n, ast.Import) else [])
            for name in imports:
                candidates = modules.get(name.split(".")[-1], [])
                require(len(candidates) <= 1, "ambiguous_source_import:" + name)
                if candidates:
                    pending.append((candidates[0], cwd))
    return sorted(callers)


def effective_csv_path(tree, fields):
    source = ast.parse((Path(tree) / "src/setup/findEffectiveCharge.py").read_text())
    assignments = [n for n in source.body if isinstance(n, ast.Assign) and
                   any(isinstance(t, ast.Name) and t.id == "OUTPATH" for t in n.targets)]
    expected = ast.parse('OUTPATH = str(Path(LTANAPATH) / "OUTPUT" / "Analysis" / f"{ANATYPE}LT")').body[0]
    require(len(assignments) == 2 and ast.dump(assignments[-1]) == ast.dump(expected), "effective_charge_writer_changed")
    return Path(eval(compile(ast.Expression(assignments[-1].value), "<actual-csv-path>", "eval"),
        {"__builtins__": {}, "Path": Path, "str": str}, fields))


def private_callers(repo, worktree, env, baseline, output, installed_package):
    records = []
    copied = worktree.parent / "runtime-python/ltsep"
    for caller, modes, cwd in reached_callers(repo):
        observations = []
        for tree, environ, package in ((repo, {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, installed_package),
                                      (worktree, env, copied)):
            code = ('import json,os,sys,ltsep,ltsep.pathing; from farm_env.ltsep_paths import create_ltsep_root; '
                    'r=create_ltsep_root(os.path.realpath(sys.argv[1]),*json.loads(sys.argv[2])); '
                    'print(json.dumps({"paths":{k:str(getattr(r,k,"")) for k in json.loads(sys.argv[3])},'
                    '"package":ltsep.__file__,"pathing":ltsep.pathing.__file__}))')
            fields = isolation.python_json(tree / cwd, code,
                (str(tree / caller), json.dumps(modes), json.dumps(isolation.PATH_FIELDS)),
                {**environ, "PYTHONPATH": str(tree) + os.pathsep + environ.get("PYTHONPATH", "")})
            require(Path(fields["package"]).resolve() == package / "__init__.py" and
                    Path(fields["pathing"]).resolve().is_relative_to(package), "caller_package_identity:" + caller)
            require(Path(fields["paths"]["LTANAPATH"]).resolve() == tree.resolve(), "caller_ltanapath:" + caller)
            observations.append(fields["paths"])
        ordinary, runtime = observations
        require(all(ordinary[k] == runtime[k] == baseline[k] for k in isolation.PATH_FIELDS
                    if k not in {"OUTPATH", "LTANAPATH"}), "caller_stable_roles:" + caller)
        if modes:
            require(Path(ordinary["OUTPATH"]).resolve() == Path(baseline["OUTPATH"]).resolve() and
                    Path(runtime["OUTPATH"]) == worktree / "OUTPUT/Analysis/KaonLT" and
                    Path(runtime["OUTPATH"]).resolve() == output, "caller_private_outpath:" + caller)
        elif caller == "src/setup/findEffectiveCharge.py":
            require(effective_csv_path(repo, ordinary).resolve() == Path(baseline["OUTPATH"]).resolve() and
                    effective_csv_path(worktree, runtime).resolve() == output, "effective_charge_private_outpath")
        else:
            require(caller in {"src/setup/getDataTable.py", "src/setup/getEfficiencyValue.py", "src/setup/calcTotalEffectiveCharge.py"}
                    and not any(isinstance(n, ast.Attribute) and n.attr == "OUTPATH" for n in
                                ast.walk(ast.parse((worktree / caller).read_text()))), "unapproved_default_writer:" + caller)
        records.append({"caller": caller, "modes": modes, "cwd": cwd, "ordinary": ordinary, "runtime": runtime})
    require(records, "caller_inventory_empty")
    return records


def diamond_choices(root):
    import glob
    root = Path(root)
    def choose(tokens, eps, fallback=False, main=False):
        best = None
        for token in tokens:
            for p in glob.glob(str(root / ("*" + token + "*kaon*Q4p4W2p74*.root"))):
                low = p.lower()
                if "simc" in low or "shiftprep" in Path(p).name.lower() or eps not in (p if main else low):
                    continue
                if fallback and ("left" in low or "right" in low):
                    continue
                if best is None or len(p) < len(best):
                    best = p
        return None if best is None else Path(best).name
    return {f"{kind}:{s}:{eps}": choose(tokens, eps, kind == "fallback", kind == "main")
            for kind, settings in (("main", [(s, [s]) for s in ("Center", "Left", "Right")]),
                ("overlay", [(s, [s, s.lower()]) for s in ("Center", "Left", "Right")]), ("fallback", [("Center", [""])]))
            for s, tokens in settings for eps in ("low", "mid", "high")}


def copy_plan(shared, output):
    shared, output = Path(shared), Path(output)
    needed = set()
    for eps in ("lowe", "highe"):
        for s in ("Center", "Left", "Right"):
            pair = [f"{s}_kaon_Analysed_{r}_Q4p4W2p74_{eps}.root" for r in ("Data", "Dummy")]
            present = [(shared / n).is_file() for n in pair]
            if not any(present):
                require(eps != "lowe" or s not in ("Center", "Left"), "required_lowe_setting_missing:" + s)
                continue
            archived = shared / f"Prod_Coin_Q4p4W2p74{s.lower()}_{eps}.root"
            if not archived.exists() and not (eps == "lowe" and s in ("Center", "Left")):
                # Frozen shift_prep/main skip this optional setting. Its raw
                # ROOTs still enter the complete Diamond namespace below.
                continue
            require(all(present), "incomplete_data_dummy_pair:" + s + ":" + eps)
            needed.update(pair)
            # Low epsilon reaches full SIMC preparation; high is raw prepass only.
            for suffix in ((".root", ".hist") if eps == "lowe" else (".root",)):
                name = f"Prod_Coin_Q4p4W2p74{s.lower()}_{eps}{suffix}"
                require((shared / name).is_file(), "required_simc_input_missing:" + name)
                needed.add(name)
            if eps == "lowe":
                # Full kaon rand_sub consumes the existing merged pion control
                # pair; -d never runs the -c merge that produces these inputs.
                for role in ("Data", "Dummy"):
                    name = f"{s}_pion_Analysed_{role}_Q4p4W2p74_{eps}.root"
                    require((shared / name).is_file(), "required_pion_input_missing:" + name)
                    needed.add(name)
    # Preserve the complete actual Diamond candidate namespace and first-tie order.
    import glob
    needed.update(Path(p).name for p in glob.glob(str(shared / "*kaon*Q4p4W2p74*.root"))
                  if "simc" not in p.lower() and "shiftprep" not in Path(p).name.lower())
    plan = []
    for source in shared.iterdir():
        if source.name not in needed:
            continue
        target = output / source.name
        require(source.is_file() and not os.path.lexists(target), "input_missing_or_collision:" + source.name)
        require(target.parent == output and target.resolve().is_relative_to(output.resolve()), "copy_path_escape")
        plan.append({"source": str(source), "target": str(target), "bytes": source.stat().st_size})
    require({Path(p["source"]).name for p in plan} == needed, "copy_inventory_changed")
    return plan


def stage_copies(shared, output, plan, records):
    require(plan == copy_plan(shared, output), "copy_plan_changed")
    for entry in plan:
        source, target = Path(entry["source"]), Path(entry["target"])
        before = file_identity(source)
        records[source.name] = {"source": str(source), "original": before}
        with source.open("rb") as src, target.open("xb") as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        after = file_identity(target)
        require(file_identity(source) == before and before["bytes"] == entry["bytes"], "original_changed_during_copy")
        require(all(before[k] == after[k] for k in ("sha256", "bytes")) and
                (before["device"], before["inode"]) != (after["device"], after["inode"]), "private_copy_identity")
        records[source.name]["runtime_container_before"] = after
    require(diamond_choices(shared) == diamond_choices(output), "diamond_selection_changed")
    return records


def verify_originals(records):
    for item in records.values():
        require(file_identity(item["source"]) == item["original"], "original_input_changed:" + item["source"])
    return {"passed": True}


def private_environment(env, data_parent):
    result = child_environment(env)
    for key, leaf in (("TMPDIR", "tmp"), ("XDG_CACHE_HOME", "cache"), ("MPLCONFIGDIR", "mpl"),
                      ("PYTHONPYCACHEPREFIX", "pycache")):
        p = Path(data_parent) / ("ab-" + leaf)
        require(not os.path.lexists(p), "private_cache_collision:" + key)
        p.mkdir(mode=0o700)
        result[key] = str(p)
    result["PYTHONDONTWRITEBYTECODE"] = "1"
    return result


def skim_inventory(repo, paths):
    """Use the real existence-only predicate; never hash or open skim ROOTs."""
    _, skim = skim_directories(repo, paths)
    result = {}
    for eps in ("lowe", "highe"):
        for setting in ("center", "left", "right"):
            for suffix in ("", "_dummy"):
                p = Path(paths["UTILPATH"]) / "run_list/KaonLT" / f"Q4p4W2p74{setting}_{eps}{suffix}"
                if not p.exists():
                    result[str(p)] = {"kind": "runlist", "exists": False}
                    continue
                raw = p.read_bytes()
                runs = raw.decode().split()
                require(all(re.fullmatch(r"\d+", run) for run in runs), "run_list_not_integer")
                result[str(p)] = {"kind": "runlist", "exists": True, "sha256": hashlib.sha256(raw).hexdigest()}
                for run in runs:
                    target = Path(skim) / f"kaon_{run}_-1_Raw_Data.root"
                    result[str(target)] = {"kind": "skim_existence", "exists": target.exists()}
    return result


class GateStatus(isolation.accepted.GateStatus):
    def __init__(self, outdir, output, source_commit):
        super().__init__(outdir, output, source_commit)
        self.record.pop("phi", None)
        self.record.pop("epsilon", None)
        self.update(schema_version="pion_t_template_comparison_gate_status/v1",
                    requested_settings=[{"phi": s, "epsilon": "lowe"} for s in ("Center", "Left")],
                    boundaries=BOUNDARIES)


def execute_gate(repo, outdir, output, source_commit):
    repo, outdir, output = (Path(p).absolute() for p in (repo, outdir, output))
    # Use the established farm-owner boundary; no new machine-specific values.
    require(repo.resolve() == isolation.accepted.REPO.resolve() and
            outdir.resolve() == isolation.accepted.ARTIFACTS.resolve() and
            output.parent.resolve() == isolation.accepted.GLOBUS.resolve() and
            output.suffix == ".zip" and outdir.is_dir() and output.parent.is_dir(), "owner_path_boundary_invalid")
    require(not output.is_symlink(), "output_symlink_invalid")
    paths = destinations(outdir, output)
    require(len(set(map(str, paths.values()))) == len(paths), "companion_paths_overlap")
    for path in paths.values():
        require(not os.path.lexists(path), "output_collision:" + str(path))
    status = GateStatus(outdir, output, source_commit)
    before = overlay = installed = external = accepted = skims = None
    copies = {}
    succeeded = False
    def preservation():
        checks = []
        if before is not None:
            checks.append(("ordinary_checkout", lambda: isolation.verify_preservation(repo, before)))
        if installed is not None:
            checks.append(("installed_ltsep", lambda: isolation.verify_ltsep_preservation(installed)))
        if accepted is not None:
            checks.append(("accepted_artifacts", lambda: verify_accepted(outdir, accepted)))
        if external is not None:
            checks.append(("external_links", lambda: verify_external(external)))
        if skims is not None:
            checks.append(("skim_existence", lambda: require(skim_inventory(repo, baseline) == skims, "skim_existence_changed")))
        checks.append(("original_inputs", lambda: verify_originals(copies)))
        results, errors = {}, []
        for name, check in checks:
            try:
                results[name] = check()
            except Exception as exc:
                results[name] = {"passed": False, "reason": str(exc)}
                errors.append(name + ":" + str(exc))
        status.update(preservation=results, staged_inputs=copies)
        require(not errors, "preservation_failed:" + ";".join(errors))
    try:
        accepted = accepted_identities(outdir)
        before = preflight(repo, source_commit)
        profile = resolved_profile(repo, source_commit)
        child_environment(os.environ)  # Reject injected opt-ins before any path work.
        status.update(accepted_comparison_targets=accepted, destinations={k: str(v) for k,v in paths.items()})
        checks, identities = source_preflight(repo, source_commit, profile)
        status.update("runtime_isolation", collector_source_preflight_completed=True,
                      collector_source_checks=checks, committed_source_identities=identities)
        # Snapshot the original package before the frozen overlay helper: its
        # own return record is unavailable if copy/import/path verification fails.
        identity = isolation.import_identity(repo)
        package = Path(identity["package_file"]).resolve().parent
        protected = [Path(identity["package_file"]), Path(identity["pathing_file"]),
                     *sorted((package / "PATH_TO_DIR").glob("*.path"))]
        installed = {"original_source_sha256": {str(p): sha256(p) for p in protected}}
        baseline = isolation.path_fields(repo, repo)
        placement = volatile_base(repo, baseline)
        status.update(runtime_placement=placement)
        with analysis_worktree(repo, source_commit, placement) as worktree:
            placement["baseline_paths"] = baseline
            placement["shared_output"] = str(outdir.resolve())
            status.update(runtime_placement=placement)
            env, overlay = prepare_runtime_overlay(repo, worktree, placement)
            data_parent = Path(placement["parent"]["path"])
            env = private_environment(env, data_parent)
            probes = isolation.probe_paths(worktree, env, overlay["baseline_paths"])
            baseline = overlay["baseline_paths"]
            require(Path(baseline["OUTPATH"]).resolve() == outdir.resolve() ==
                    (Path(baseline["VOLATILEPATH"]) / "OUTPUT/Analysis/KaonLT").resolve(), "ordinary_output_authority")
            external = external_setup(worktree, probes[0]["paths"], outdir)
            private = private_output(worktree, probes[0]["paths"], outdir, data_parent)
            callers = private_callers(repo, worktree, env, baseline, private, package)
            skims = skim_inventory(worktree, baseline)
            plan = copy_plan(outdir, private)
            status.update("input_staging", copy_plan=plan, copy_files=len(plan),
                          copy_bytes=sum(p["bytes"] for p in plan), skim_existence_inventory=skims)
            print(json.dumps({"copy_plan": plan, "copy_bytes": sum(p["bytes"] for p in plan)}), file=sys.stderr, flush=True)
            verify_runtime_layout(repo, worktree, placement, private, env)
            status.update(runtime_placement=placement, private_output=str(private))
            stage_copies(outdir, private, plan, copies)
            verify_runtime_layout(repo, worktree, placement, private, env)
            verify_private(worktree, private)
            verify_external(external)
            verify_originals(copies)
            prior = {name: None if not os.path.lexists(private / name) else {"path": str(private / name)} for name in names()}
            require(all(v is None for v in prior.values()), "comparison_artifact_preexists")
            status.update("analysis", analysis_started=True, analysis_worktree=str(worktree),
                          isolation_probes=probes, private_output=str(private), private_callers=callers,
                          staged_inputs=copies, external_symlink_preflight=external)
            verify_runtime_layout(repo, worktree, placement, private, env)
            started = time.time_ns()
            returncode = debug.run_analysis(worktree, paths["log"], env)
            status.update("completion_markers", analysis_returncode=returncode,
                          log_identity={"sha256": sha256(paths["log"]), "bytes": paths["log"].stat().st_size})
            verify_runtime_layout(repo, worktree, placement, private, env)
            verify_completion(paths["log"], returncode)
            status.update("verify_artifacts", analysis_completed=True)
            runtime = {name: file_identity(private / name) for name in copies}
            records, outcomes = verify_artifacts(private, source_commit, identities, started, prior, runtime)
            require(skim_inventory(worktree, baseline) == skims, "skim_existence_changed")
            status.update(artifact_verification_completed=True, artifacts=records, diagnostic_outcomes=outcomes)
            preservation()
            status.update("collection", runtime_container_identities=runtime)
            verify_private(worktree, private)
            effective = worktree.parent / "ab-profile.json"
            with effective.open("x", encoding="utf-8") as handle:
                json.dump(profile, handle)
            with isolation.owned_worktree(repo, source_commit, "source") as source:
                require(source.resolve() != worktree.resolve() and not os.path.lexists(source / "OUTPUT"),
                        "collector_source_not_independent")
                require(resolved_profile(source, source_commit) == profile, "detached_profile_mismatch")
                with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
                    result = module.collect_validation_bundle(outdir=private, kinematic="Q4p4W2p74",
                        output=output, profile_path=effective, repo_root=source)
                require(result["returncode"] == 0, "collection_failed")
                status.update("verify_zip", collection_completed=True)
                verify_zip(output, source_commit, records)
                status.update(zip_verification_completed=True, zip_identity={"sha256": sha256(output),
                              "bytes": output.stat().st_size})
            require(not os.path.lexists(source), "source_cleanup_incomplete")
        require(not os.path.lexists(worktree.parent) and not os.path.lexists(data_parent), "private_cleanup_incomplete")
        status.update(worktree_cleanup_completed=True)
        preservation()
        summary = {"schema_version": "pion_t_template_comparison_run_summary/v1",
            "source_commit": source_commit, "analysis_command": COMMAND, "analysis_returncode": returncode,
            "diagnostic_outcomes": outcomes, "artifacts": records, "boundaries": BOUNDARIES,
            "runtime_placement": placement, "copy_bytes": sum(p["bytes"] for p in plan),
            "copy_plan": plan, "original_inputs": copies, "runtime_container_identities": runtime,
            "log_identity": status.record["log_identity"], "zip_identity": status.record["zip_identity"]}
        with paths["summary"].open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        isolation.deliver_evidence(status, paths["log"], paths["summary"], output)
        succeeded = True
        return output
    except Exception as exc:
        status.update(status="failed", failure_reason=str(exc))
        raise
    finally:
        if not succeeded:
            try:
                preservation()
            except Exception as exc:
                status.update(status="failed", preservation_failure_reason=str(exc))
                raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = execute_gate(args.repo, args.outdir, args.output, args.source_commit)
    except Exception as exc:
        print("Pion template gate failed: " + str(exc), file=sys.stderr)
        print("Owner gate status: " + str(isolation.accepted.gate_status_path(args.outdir, args.output)), file=sys.stderr)
        return 1
    print(result.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
