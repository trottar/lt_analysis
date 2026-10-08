#!/usr/bin/env python3
"""Own isolated Q4p4W2p74 full analysis, five-setting verification and collection.

Farm execution requires independent review, push, synchronization and readiness
approval. --source-commit is the exact reviewed/pushed execution identity.
Success stdout is one ZIP path; the per-attempt status owns final ZIP identity.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

try:
    from . import run_e8_4_fix5_left_lowe_plot_gate as accepted
except ImportError:
    import run_e8_4_fix5_left_lowe_plot_gate as accepted

collector = accepted.collector
require = accepted.require
sha256 = accepted.sha256
git = accepted.git
collection_module = accepted.collection_module
KINEMATIC = accepted.KINEMATIC
BASE_HEAD = "ed378e0f30a357c6293da64f8f8f3bfbe187d1fe"
PROFILE = "testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json"
PROFILE_ID = "phase_e8_4_fix5_canonical_five_isolated_runtime/v1"
SUMMARY = "Q4p4W2p74_e8_4_fix5_canonical_five_run-summary.json"
MATERIALIZATION_HEAD = "463d2657f696ecee33113edc3393ac51083a8944"
MATERIALIZATION_PREFIX = KINEMATIC + "_kaon_pion-background_hgcer-method-a-"
MATERIALIZATION_NAMES = {
    "comparison": MATERIALIZATION_PREFIX + "current-baseline-authority-comparison-input.json",
    "manifest": MATERIALIZATION_PREFIX + "current-baseline-authority-materialization-manifest.json",
    "f2": MATERIALIZATION_PREFIX + "acceptance-representation-current-baseline-candidate.json",
    "f3": MATERIALIZATION_PREFIX + "acceptance-map-current-baseline-candidate.json",
    "f4": MATERIALIZATION_PREFIX + "parent-preserving-correction-current-baseline-candidate.json"}
MATERIALIZATION_SHA256 = {
    "comparison": "4c3fdab5d05965a1b3f1dc835c9d8c529896e6eaf2da8eb5c8ac8fc26e95257a",
    "manifest": "e7c37f55b24739e8be0fb778a6c3491e17df746263e6d9d71a4f31ec65cd7908",
    "f2": "2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e",
    "f3": "c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228",
    "f4": "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7"}
CANDIDATES = {MATERIALIZATION_NAMES[k]: MATERIALIZATION_SHA256[k] for k in ("f2", "f3", "f4")}
F1_SHA256 = {
    "Left-lowe": "eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07",
    "Left-highe": "544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e",
    "Center-lowe": "2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16",
    "Center-highe": "c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941",
    "Right-highe": "e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652"}
HISTORICAL_INPUT_SHA256 = {
    "f2": "87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da",
    "f3": "04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95",
    "f4": "adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188"}
CANDIDATE_FINGERPRINTS = {
    "f2": {"representation_fingerprint": "e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216",
           "artifact_fingerprint": "87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6"},
    "f3": {"map_fingerprint": "6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548",
           "algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
           "artifact_fingerprint": "e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302"},
    "f4": {"correction_fingerprint": "71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98",
           "artifact_fingerprint": "0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4"}}
F3_RECONSTRUCTION = {KINEMATIC: {"source_file_sha256": MATERIALIZATION_SHA256["f3"],
    **CANDIDATE_FINGERPRINTS["f3"], "farm_source_head": "0" * 40}}
SCIENTIFIC_GATE = {"f2_scientific_payload_match": True, "f3_scientific_payload_match": True,
                   "f4_scientific_payload_match": False, "first_changed_stage": "F4"}
FARM_OUTPUTS = accepted.FARM_OUTPUTS
GATE_PATHS = ("src", "testing", "tools", "farm_env", "background_samples",
              "run_Prod_Analysis.sh", "set_SymLinks.sh")
# The unchanged collector's schema fixes serialization order. The scientific
# inventory is exactly these five pairs, with no Right/lowe.
CANONICAL_SETTINGS = accepted.CANONICAL_SETTINGS
COMMAND = ["./run_Prod_Analysis.sh", "4p4", "2p74"]
PATH_FIELDS = ("VOLATILEPATH", "ANALYSISPATH", "HCANAPATH", "REPLAYPATH",
               "UTILPATH", "PACKAGEPATH", "OUTPATH", "ROOTPATH", "SKIMPATH",
               "REPORTPATH", "CUTPATH", "PARAMPATH", "SCRIPTPATH", "ANATYPE",
               "USER", "HOST", "SIMCPATH", "LTANAPATH")


def validate_settings(settings):
    require(settings == list(CANONICAL_SETTINGS), "canonical_five_inventory_invalid")


def verify_candidate_materialization(directory):
    directory = Path(directory).resolve()
    require(directory.is_dir(), "candidate_materialization_directory_missing")
    payloads, identities = {}, {}
    for key, name in MATERIALIZATION_NAMES.items():
        path = directory / name
        require(path.is_file() and sha256(path) == MATERIALIZATION_SHA256[key],
                "materialization_hash_mismatch:" + key)
        payloads[key] = collector._strict_json_payload(path)
        require(isinstance(payloads[key], dict), "materialization_payload_invalid:" + key)
        identities[key] = {"path": str(path), "sha256": MATERIALIZATION_SHA256[key], "bytes": path.stat().st_size}
    manifest, comparison = payloads["manifest"], payloads["comparison"]
    require(manifest.get("schema_version") == "method_a_current_baseline_authority_materialization/v1" and
            manifest.get("source_head") == MATERIALIZATION_HEAD and manifest.get("kinematic_token") == KINEMATIC,
            "materialization_source_or_schema_mismatch")
    for key, expected in {"complete": True, "non_authoritative": True, "accepted_authority_mutated": False,
                          "production_objects_mutated": False, "production_application_performed": False,
                          "method_a_promoted": False}.items():
        require(manifest.get(key) is expected, "materialization_flag_mismatch:" + key)
    require(manifest.get("errors") == [], "materialization_errors_present")
    require(comparison.get("schema_version") == "method_a_current_baseline_authority_comparison/v1" and
            comparison.get("non_authoritative") is True, "comparison_schema_invalid")
    require(comparison.get("accepted_file_sha256") == HISTORICAL_INPUT_SHA256,
            "materialization_historical_inputs_mismatch")
    accepted_inputs = manifest.get("accepted_inputs", {})
    require(set(accepted_inputs) == set(HISTORICAL_INPUT_SHA256) and
            all(accepted_inputs[k].get("raw_sha256") == v for k, v in HISTORICAL_INPUT_SHA256.items()),
            "materialization_historical_inputs_mismatch")
    row = manifest.get("comparison_input", {})
    require(row.get("raw_sha256") == MATERIALIZATION_SHA256["comparison"] and
            row.get("copied_output_raw_sha256") == MATERIALIZATION_SHA256["comparison"] and
            row.get("copied_output_basename") == MATERIALIZATION_NAMES["comparison"] and
            row.get("schema_version") == comparison["schema_version"], "materialization_comparison_mismatch")
    f1 = manifest.get("current_f1_inputs")
    require(isinstance(f1, list) and len(f1) == 5 and f1 == comparison.get("f1_inputs") and
            {r.get("setting_id"): r.get("source_file_sha256") for r in f1} == F1_SHA256,
            "materialization_f1_inventory_mismatch")
    require(all(r.get("alias") == r.get("setting_id") and
                r.get("setting", {}).get("kinematic_token") == KINEMATIC for r in f1),
            "materialization_f1_setting_mismatch")
    for key, value in SCIENTIFIC_GATE.items():
        for gate in (manifest.get("scientific_gate", {}), comparison.get("summary", {})):
            require((gate.get(key) is value) if isinstance(value, bool) else gate.get(key) == value,
                    "materialization_scientific_gate_mismatch")
    require(comparison.get("candidate_serialized_sha256") ==
            {k: MATERIALIZATION_SHA256[k] for k in ("f2", "f3")}, "comparison_candidate_identity_mismatch")
    require(comparison.get("diagnostic_f3_authority_override") == F3_RECONSTRUCTION and
            manifest.get("diagnostic_f3_authority_override") == {
                "accepted_authority": False, "purpose": "diagnostic_candidate_construction_only",
                "record": F3_RECONSTRUCTION}, "materialization_f3_reconstruction_mismatch")
    outputs = manifest.get("candidate_outputs", {})
    require(set(outputs) == {"f2", "f3", "f4"}, "materialization_candidate_inventory_mismatch")
    for stage, body_name, fingerprint in (("f2", "representation", "representation_fingerprint"),
            ("f3", "acceptance_map", "map_fingerprint"), ("f4", "correction", "correction_fingerprint")):
        artifact, row = payloads[stage], outputs[stage]
        body = artifact.get(body_name, {})
        expected = {"basename": MATERIALIZATION_NAMES[stage], "raw_sha256": MATERIALIZATION_SHA256[stage],
                    "non_authoritative": True, **CANDIDATE_FINGERPRINTS[stage]}
        require(all(row.get(k) == v for k, v in expected.items()) and artifact.get("non_authoritative") is True and
                artifact.get("artifact_fingerprint") == expected["artifact_fingerprint"] and
                body.get("fingerprint") == expected[fingerprint], "materialization_candidate_identity_mismatch:" + stage)
        if stage == "f3":
            require(body.get("algorithm_fingerprint") == expected["algorithm_fingerprint"],
                    "materialization_candidate_identity_mismatch:f3")
        if stage == "f4":
            inherited = {"f3_source_file_sha256": MATERIALIZATION_SHA256["f3"],
                         **{"f3_" + k: v for k, v in CANDIDATE_FINGERPRINTS["f3"].items()}}
            require(all(row.get(k) == v and body.get(k) == v for k, v in inherited.items()),
                    "materialization_f4_inherited_identity_mismatch")
    return {"passed": True, "source_head": MATERIALIZATION_HEAD, "files": identities,
            "f1_sha256": dict(F1_SHA256), "scientific_gate": dict(SCIENTIFIC_GATE)}


def stage_candidates(directory, outdir, verification, status=None):
    require(verification.get("passed") is True, "candidate_materialization_not_verified")
    require(not outdir.resolve().is_relative_to(Path(directory).resolve()), "candidate_destination_inside_materialization")
    plans = []
    for name, digest in CANDIDATES.items():
        source, target = Path(directory).resolve() / name, outdir / name
        require(source.resolve() != target.resolve() and not target.is_symlink(), "candidate_source_target_alias")
        before = sha256(target) if target.is_file() else None
        require(not target.exists() or target.is_file(), "candidate_target_not_file")
        require(before in (None, accepted.CANDIDATES.get(name), digest), "unknown_candidate_target:" + name)
        require(sha256(source) == digest, "materialization_source_changed:" + name)
        plans.append((source, target, before, digest))
    records = []
    for source, target, before, digest in plans:
        action = "no-op" if before == digest else "installed" if before is None else "replaced"
        if before != digest:
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(prefix=".candidate-stage-", dir=outdir, delete=False) as handle:
                    temporary = Path(handle.name)
                    with source.open("rb") as stream: shutil.copyfileobj(stream, handle)
                    handle.flush(); os.fsync(handle.fileno())
                require(sha256(temporary) == digest, "candidate_staging_hash_mismatch")
                require(not target.is_symlink() and (sha256(target) if target.is_file() else None) == before,
                        "candidate_target_changed_during_staging")
                os.replace(temporary, target)
            finally:
                if temporary is not None and temporary.exists(): temporary.unlink()
        require(sha256(target) == digest, "candidate_installation_hash_mismatch")
        records.append({"path": str(target), "before_sha256": before, "after_sha256": digest, "action": action})
        if status is not None: status.update(candidate_installation={"passed": False, "files": list(records)})
    return {"passed": True, "files": records}


def validate_candidate_lineage(outdir, module):
    """Read JSON and run the existing calculators; retain aggregate evidence only."""
    paths = module.accepted_f6_3_artifact_paths(outdir, KINEMATIC)
    f1, f2, f3, f4, hashes = module.load_accepted_f6_3_authority(paths, include_f2=True)
    expected_f1 = module.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256
    require(expected_f1 == F1_SHA256 and
            {k: hashes.get(k) for k in expected_f1} == expected_f1, "lineage_f1_identity_mismatch")
    require(hashes.get("f2") == CANDIDATES[Path(paths["f2"]).name] and
            hashes.get("f3") == CANDIDATES[Path(paths["f3"]).name] and
            hashes.get("f4") == CANDIDATES[Path(paths["f4"]).name], "lineage_f2_f3_f4_identity_mismatch")
    settings = []
    for setting_id in F1_SHA256:
        factors, provenance, rows = module.reconstruct_transient_factor_map(f1, f3, f4,
            f2_artifact=f2, f2_input_file_sha256=hashes["f2"],
            f1_input_file_hashes={k: hashes[k] for k in expected_f1}, f3_input_file_sha256=hashes["f3"],
            f4_input_file_sha256=hashes["f4"], setting_id=setting_id)
        require(factors and set(factors) == set(rows) and
                all(math.isfinite(float(v)) and float(v) > 0 for v in factors.values()),
                "lineage_factor_population_invalid:" + setting_id)
        require(provenance.get("selected_setting_id") == setting_id and
                provenance.get("schema_version") == module.F6_3_PARALLEL_SCHEMA_VERSION and
                provenance.get("candidate_validation_source_head") == MATERIALIZATION_HEAD and
                provenance.get("accepted_f1_source_file_sha256") == expected_f1 and
                provenance.get("candidate_f3_source_file_sha256") == hashes["f3"] and
                provenance.get("candidate_f2_source_file_sha256") == hashes["f2"] and
                provenance.get("scientific_equivalence", {}).get("all_stages_passed") is True and
                provenance.get("scientific_equivalence", {}).get("mode") == "exact-lineage" and
                provenance.get("candidate_f4_source_file_sha256") == hashes["f4"] and
                provenance.get("transient_factor_population_count") == len(factors) and
                provenance.get("current_baseline_candidate_lineage") is True and
                provenance.get("branch_role") == "parallel_nonproduction_method_a_full_analysis" and
                provenance.get("candidate_f3_reconstruction_role") == "candidate_construction_sentinel_only" and
                re.fullmatch(r"[0-9a-f]{64}", str(provenance.get("transient_factor_identity_fingerprint"))) is not None and
                provenance.get("accepted_f4_runtime_authority", {}).get("accepted_authority_match") is True and
                provenance.get("accepted_f3_source_file_sha256") == hashes["f3"] and
                provenance.get("accepted_f4_source_file_sha256") == hashes["f4"] and
                all(provenance.get("accepted_f4_" + k) ==
                    module.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC[KINEMATIC][k]
                    for k in ("correction_fingerprint", "artifact_fingerprint")) and
                all(provenance.get("accepted_f3_" + k) == v for k, v in
                    module.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC[KINEMATIC].items()
                    if k in {"map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"}) and
                provenance.get("production_promotion_performed") is False and
                provenance.get("event_correction_persisted") is False,
                "lineage_setting_provenance_invalid:" + setting_id)
        settings.append({"setting_id": setting_id, "factor_count": len(factors),
                         "f4_shared_reproduction_passed": True, "provenance": provenance})
    return {"passed": True, "observed_sha256": hashes, "settings": settings,
            "f4_shared_reproduction_passed": True}


def lineage_preflight(worktree, outdir, env):
    return python_json(worktree,
        'import json,sys; from pathlib import Path; '
        'sys.path[:0]=[str(Path.cwd()/"src/cuts"),str(Path.cwd()/"src/utility")]; '
        'import pion_hgcer_method_a_parallel_full_procedure as f63; '
        'from testing import run_e8_4_fix5_canonical_five_plot_gate as owner; '
        'print(json.dumps(owner.validate_candidate_lineage(Path(sys.argv[1]),f63),allow_nan=False))',
        (str(outdir),), env=env)


def snapshot(repo):
    return {
        "branch": git(repo, "branch", "--show-current"),
        "head": git(repo, "rev-parse", "HEAD"),
        "origin_test": git(repo, "rev-parse", "refs/remotes/origin/test"),
        "porcelain": git(repo, "status", "--porcelain=v1", "--untracked-files=all"),
        "farm_local_sha256": {name: sha256(repo / name) for name in sorted(FARM_OUTPUTS)
                              if (repo / name).is_file()},
    }


def preflight(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None,
            "source_commit_requires_full_sha")
    state = snapshot(repo)
    require(state["branch"] == "test", "wrong_branch")
    require(state["head"] == source_commit, "wrong_head")
    require(state["origin_test"] == source_commit, "source_not_observed_pushed")
    require(source_commit != BASE_HEAD, "owner_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, source_commit)
    for line in state["porcelain"].splitlines():
        path = line[3:]
        require(path in FARM_OUTPUTS or path.startswith("OUTPUT/") or
                (line[:2] == "??" and re.fullmatch(r"kaonlt_review(?:\([^/]+\))?\.diff", path)),
                "dirty_gate_source:" + line)
    git(repo, "diff", "--check", "HEAD", "--", *GATE_PATHS,
        *(":(exclude)" + name for name in sorted(FARM_OUTPUTS)))
    return state


def verify_preservation(repo, before):
    after = snapshot(repo)
    for key in before:
        require(after[key] == before[key], "ordinary_checkout_changed:" + key)
    return {"passed": True, "snapshot": after}


def resolved_profile(repo, source_commit):
    profile = collector.load_validation_profile(repo / PROFILE)
    validate_settings(profile["settings"])
    require(profile["validation_profile"] == PROFILE_ID and
            profile["collection_mode"] == "generic_artifacts", "profile_mode_invalid")
    require(profile["source_identity"] == {
        "required_analysis_commit": BASE_HEAD, "allowed_committed_files": [],
        "allowed_non_analysis_path_prefixes": ["docs/memory/"]}, "profile_source_identity_invalid")
    profile["source_identity"]["required_analysis_commit"] = source_commit
    return profile


def cleanup_path(repo, parent, worktree, leaf):
    parent = parent.resolve()
    root = repo.resolve()
    require(parent.parent == root.parent and parent != root and not parent.is_relative_to(root) and
            not parent.is_symlink() and worktree.resolve() == parent / leaf and
            worktree.resolve() != root and not worktree.resolve().is_relative_to(root),
            "worktree_cleanup_path_invalid")


@contextmanager
def owned_worktree(repo, source_commit, leaf):
    require(leaf in {"analysis", "source"}, "worktree_role_invalid")
    with tempfile.TemporaryDirectory(prefix="kaonlt-e8-4-canonical-five-" + leaf + "-",
                                     dir=repo.resolve().parent) as directory:
        parent = Path(directory).resolve()
        worktree = parent / leaf
        cleanup_path(repo, parent, worktree, leaf)
        git(repo, "worktree", "add", "--detach", str(worktree), source_commit)
        try:
            require(git(worktree, "rev-parse", "HEAD") == source_commit, "worktree_head_mismatch")
            require(git(worktree, "branch", "--show-current") == "", "worktree_not_detached")
            require(git(worktree, "status", "--porcelain=v1", "--untracked-files=all") == "",
                    "worktree_dirty")
            yield worktree
        finally:
            cleanup_path(repo, parent, worktree, leaf)
            git(repo, "worktree", "remove", "--force", str(worktree))


def collector_source_preflight(repo, source_commit, profile):
    with owned_worktree(repo, source_commit, "source") as worktree:
        with collection_module(worktree) as module, redirect_stdout(sys.stderr):
            require(resolved_profile(worktree, source_commit) == profile, "collector_profile_mismatch")
            _, checks = module.collect_source_checks(worktree,
                required_analysis_commit=source_commit, allowed_committed_files=[])
            for check in checks:
                require(check["returncode"] == 0, "collector_source_check_failed:" +
                        check["name"] + ":" + check["stderr"])
            ancestor, _, unexpected = module._committed_identity(checks, [], ["docs/memory/"])
            require(ancestor and not unexpected, "collector_source_identity_invalid")
    return checks


def python_json(cwd, code, arguments=(), env=None):
    result = subprocess.run(["python3", "-B", "-c", code, *arguments], cwd=cwd,
                            env=env, capture_output=True, text=True)
    require(result.returncode == 0, "ltsep_python_probe_failed:" + result.stderr)
    return json.loads(result.stdout)


def import_identity(cwd, env=None):
    return python_json(cwd, 'import json,ltsep,ltsep.pathing; '
        'print(json.dumps({"package_file":ltsep.__file__,"pathing_file":ltsep.pathing.__file__}))', env=env)


def path_fields(tree, caller, env=None):
    result = subprocess.run(["python3", "-B", str(tree / "farm_env/print_ltsep_path_fields.py"),
                             str(caller)], cwd=tree, env=env, capture_output=True, text=True)
    require(result.returncode == 0, "ltsep_probe_failed:" + str(caller) + ":" + result.stderr)
    fields = result.stdout.strip().split(",")
    require(len(fields) == len(PATH_FIELDS), "ltsep_probe_fields_invalid")
    return dict(zip(PATH_FIELDS, fields))


def analysis_environment(overlay_parent, ambient=None):
    env = dict(os.environ if ambient is None else ambient)
    for key in ("LT_ANALYSIS_DEBUG_LEFT_LOW", "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE",
                "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING"):
        env.pop(key, None)
    previous = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(overlay_parent) + (os.pathsep + previous if previous else "")
    return env


def read_path_config(raw):
    """Parse ltsep's plain KEY=value .path format without changing its bytes."""
    values = {}
    for line in raw.decode("utf-8").splitlines():
        if not line.strip():
            continue
        key, separator, value = line.partition("=")
        require(separator and key.strip() and key.strip() not in values, "ltsep_config_invalid")
        values[key.strip()] = value.strip()
    return values


def prepare_runtime_overlay(repo, worktree):
    identity = import_identity(repo)
    package_file = Path(identity["package_file"]).resolve()
    package = package_file.parent
    pathing = Path(identity["pathing_file"]).resolve()
    require(package_file.name == "__init__.py" and pathing.is_relative_to(package),
            "ltsep_package_layout_invalid")
    baseline = path_fields(repo, repo)
    require(baseline["LTANAPATH"] and Path(baseline["LTANAPATH"]).is_absolute() and
            Path(baseline["LTANAPATH"]).resolve() == repo.resolve(), "baseline_ltanapath_not_ordinary_repo")
    parent = worktree.parent.resolve()
    cleanup_path(repo, parent, worktree, "analysis")
    overlay = parent / "runtime-python"
    require(not overlay.exists(), "ltsep_overlay_already_exists")
    overlay.mkdir()
    copied = overlay / "ltsep"
    shutil.copytree(package, copied, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    configs = []
    for file in sorted((copied / "PATH_TO_DIR").glob("*.path")):
        values = read_path_config(file.read_bytes())
        if "LTANAPATH" not in values:
            continue
        # Use the installed package's own raw SetPath lookup to account for
        # Root's derived OUTPATH/CUTPATH, without guessing those transformations.
        observed = python_json(repo,
            'import json,sys; from ltsep.pathing import SetPath; '
            'p=SetPath(sys.argv[1]); '
            'print(json.dumps({k:p.getPath(k) for k in json.loads(sys.argv[2])}))',
            (str(repo), json.dumps(list(values))))
        expanded = {key: value.replace("${USER}", baseline["USER"]) for key, value in values.items()}
        if expanded == observed:
            configs.append(file)
    require(len(configs) == 1, "ltsep_matching_config_not_unique")
    selected = configs[0]
    original_config = package / selected.relative_to(copied)
    original_bytes = selected.read_bytes()
    require(original_config.read_bytes() == original_bytes, "ltsep_copy_identity_mismatch")
    # Replace the value span only; keep all other fields, whitespace and line
    # endings byte-for-byte. Never patch installed files or Python code.
    matches = list(re.finditer(rb"(?m)^(LTANAPATH[ \t]*=[ \t]*)([^\r\n]*?)([ \t]*)(\r?$)", original_bytes))
    require(len(matches) == 1, "ltsep_ltanapath_assignment_not_unique")
    match = matches[0]
    patched = original_bytes[:match.start(2)] + str(worktree).encode("utf-8") + original_bytes[match.end(2):]
    selected.write_bytes(patched)
    protected = {str(file): sha256(file) for file in (package_file, pathing, original_config)}
    env = analysis_environment(overlay)
    imported = import_identity(worktree, env)
    require(Path(imported["package_file"]).resolve() == (copied / "__init__.py").resolve() and
            Path(imported["pathing_file"]).resolve().is_relative_to(copied.resolve()),
            "ltsep_overlay_import_identity_invalid")
    record = {"original_package_identity": identity, "baseline_paths": baseline,
        "overlay_parent": str(overlay), "selected_copied_config": str(selected),
        "original_config_path": str(original_config), "original_copied_file_sha256": hashlib.sha256(original_bytes).hexdigest(),
        "patched_copied_file_sha256": sha256(selected), "original_source_sha256": protected,
        "child_import_identity": imported, "import_isolation_passed": True}
    verify_ltsep_preservation(record)
    return env, record


def verify_ltsep_preservation(record):
    for name, digest in record["original_source_sha256"].items():
        require(Path(name).is_file() and sha256(name) == digest, "installed_ltsep_changed:" + name)
    return {"passed": True, "original_source_sha256": record["original_source_sha256"]}


def probe_paths(worktree, env=None, baseline=None):
    probes = []
    for caller in (worktree, worktree / "src/setup/set_sig_fortran.py",
                   worktree / "set_SymLinks.sh"):
        paths = path_fields(worktree, caller, env)
        require(paths["LTANAPATH"] and Path(paths["LTANAPATH"]).is_absolute() and
                Path(paths["LTANAPATH"]).resolve() == worktree.resolve(),
                "ltanapath_not_isolated:" + str(caller))
        require(all(paths[key] for key in ("SIMCPATH", "VOLATILEPATH", "ANATYPE")),
                "ltsep_external_paths_missing")
        probes.append({"caller": str(caller), "paths": paths, "isolation_passed": True})
    require(all(item["paths"] == probes[0]["paths"] for item in probes), "ltsep_probe_paths_disagree")
    if baseline is not None:
        stable = set(PATH_FIELDS) - {"LTANAPATH", "OUTPATH"}
        require(all(item["paths"][key] == baseline[key] for item in probes for key in stable),
                "ltsep_stable_fields_changed")
        expected_out = worktree / "OUTPUT/Analysis" / (baseline["ANATYPE"] + "LT")
        require(all(Path(item["paths"]["OUTPATH"]).is_absolute() and
                    Path(item["paths"]["OUTPATH"]).resolve() == expected_out.resolve()
                    for item in probes), "ltsep_plot_ltsep_outpath_invalid")
    return probes


def external_symlink_preflight(worktree, paths):
    """Read external links only; never invoke the mutation-bearing setup script."""
    simc = Path(paths["SIMCPATH"])
    require(simc.is_absolute() and simc.is_dir(), "simc_path_invalid")
    records = []
    for leaf in ("OUTPUTS", "input", "worksim"):
        link = simc / leaf
        require(link.is_symlink() and link.exists(), "external_simc_link_missing_or_broken:" + leaf)
        records.append({"path": str(link), "literal_target": os.readlink(link),
                        "resolved_target": str(link.resolve()), "mutation_required": False})
    # Source the exact tracked Bash configuration, as set_SymLinks does. The
    # current configuration contains assignments only; do not invent path values.
    result = subprocess.run(["bash", "-c",
        'source "$1" || exit; printf "%s" "${BACKGROUND_SIMC_PATH:-}"',
        "kaonlt-background-path", str(worktree / "background_samples/background_samples.conf")],
        cwd=worktree, env={**os.environ, "USER": paths["USER"]}, capture_output=True, text=True)
    require(result.returncode == 0, "background_config_resolution_failed:" + result.stderr)
    configured = result.stdout
    background = Path(configured) if configured else None
    if background is not None:
        require(background.is_absolute(), "background_simc_path_not_absolute")
        if background.is_dir():
            link = background / "worksim"
            expected = paths["VOLATILEPATH"] + "/worksim/"
            require(link.is_symlink() and link.exists() and os.readlink(link) == expected,
                    "background_simc_worksim_would_mutate")
            records.append({"path": str(link), "literal_target": os.readlink(link),
                            "resolved_target": str(link.resolve()), "mutation_required": False})
    return {"passed": True, "background_simc_path": configured, "links": records}


def run_analysis(worktree, log_path, env=None):
    with log_path.open("xb") as log:
        process = subprocess.Popen(COMMAND, cwd=worktree, env=env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b""):
            log.write(line)
            log.flush()
            sys.stderr.write(line.decode("utf-8", errors="replace"))
            sys.stderr.flush()
        return process.wait()


def verify_completion(log_path, returncode):
    require(returncode == 0, "analysis_failed")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    require("Low Epsilon Completed!" in text and "High Epsilon Completed!" in text,
            "full_analysis_completion_missing")
    require("Left / lowe debug analysis completed." not in text, "debug_analysis_not_canonical_five")


def pdf_name(phi, epsilon):
    return f"{phi}_kaon_rand_sub_{KINEMATIC}_{epsilon}_full-background-subtraction.pdf"


def verify_pages(payload, phi, epsilon):
    require({"phi": phi, "epsilon": epsilon} in CANONICAL_SETTINGS, "setting_not_canonical")
    require(payload.get("schema_version") == "full_background_subtraction_page_manifest/v1",
            "page_manifest_schema_invalid")
    expected = {"kinematic_token": KINEMATIC, "phi_setting": phi,
                "epsilon_filename_token": epsilon, "particle_type": "kaon",
                "epsilon_setting": {"lowe": "low", "highe": "high"}[epsilon]}
    setting = payload.get("setting")
    require(isinstance(setting, dict) and all(setting.get(k) == v for k, v in expected.items()),
            "page_manifest_setting_invalid")
    require(payload.get("pdf_basename") == pdf_name(phi, epsilon), "page_manifest_pdf_invalid")
    require(payload.get("renderer_failures") == [], "renderer_failures_nonempty")
    pages = payload.get("pages")
    require(isinstance(pages, list) and all(isinstance(page, dict) for page in pages),
            "page_inventory_invalid")
    ids = Counter(page.get("page_id") for page in pages)
    scopes = Counter((page.get("page_id"), page.get("scope")) for page in pages)
    require(all(ids[key] == 1 for key in accepted.NEW_PAGE_IDS), "new_page_missing_or_duplicate")
    require(all(scopes[key] == 1 for key in accepted.CRITICAL_PAGES), "critical_page_missing_or_duplicate")
    for page in pages:
        if page.get("page_id") == accepted.NEW_PAGE_IDS[-1]:
            require(page.get("scope") == "setting", "parent_closure_scope_invalid")
        elif page.get("page_id") in accepted.NEW_PAGE_IDS:
            index = int(page["page_id"][-1]) - 1
            require(page.get("scope") == "t" + str(index + 1), "new_page_scope_invalid")
            require(page.get("t_index") == index, "new_page_t_identity_invalid")
            require([child.get("phi_index") for child in page.get("represented_phi_inventory", [])]
                    == list(range(9)), "new_page_child_inventory_invalid")
    # These pages are emitted only after the unchanged producer has validated
    # live parity, nonproduction/Method-B flags and candidate parent closure.
    return len(pages)


def artifact_names(profile):
    names = set()
    for scope in ("global", "settings"):
        for entry in profile["artifacts"][scope]:
            for setting in CANONICAL_SETTINGS:
                names.add(entry["basename_template"].format(kinematic=KINEMATIC, **setting))
    return names


def verify_artifacts(outdir, profile, started_ns, before):
    records = {}
    for name in sorted(artifact_names(profile) - {SUMMARY}):
        path = outdir / name
        require(path.is_file() and path.stat().st_size > 0, "artifact_missing:" + name)
        digest = sha256(path)
        if name in CANDIDATES:
            require(digest == CANDIDATES[name], "candidate_identity_mismatch:" + name)
        else:
            require(path.stat().st_mtime_ns >= started_ns and
                    (path.stat().st_mtime_ns, path.stat().st_size, digest) != before.get(name),
                    "artifact_stale:" + name)
        if path.suffix == ".json":
            collector._strict_json_payload(path)
        records[name] = {"sha256": digest, "bytes": path.stat().st_size}
    settings = []
    for setting in CANONICAL_SETTINGS:
        phi, epsilon = setting["phi"], setting["epsilon"]
        pdf = outdir / pdf_name(phi, epsilon)
        require(pdf.read_bytes().startswith(b"%PDF-"), "procedure_pdf_invalid:" + pdf.name)
        pages = collector._strict_json_payload(pdf.with_name(pdf.stem + "-manifest.json"))
        count = verify_pages(pages, phi, epsilon)
        full = collector._strict_json_payload(outdir / f"kaon_FullAnalysis_{KINEMATIC}_{epsilon}.json")
        inp = full.get("inpDict", {})
        epsset = {"lowe": "low", "highe": "high"}[epsilon]
        require(all(inp.get(key) == value for key, value in {
            "ParticleType": "kaon", "EPSSET": epsset, "Q2": "4p4", "W": "2p74",
            "OutFilename": f"FullAnalysis_{KINEMATIC}_{epsilon}"}.items()),
            "full_analysis_identity_invalid")
        histlist = full.get("histlist")
        require(isinstance(histlist, list) and
                Counter(hist.get("phi_setting") for hist in histlist) ==
                Counter(s["phi"] for s in CANONICAL_SETTINGS if s["epsilon"] == epsilon),
                "full_analysis_setting_inventory_invalid")
        ledger = collector._strict_json_payload(outdir /
            f"kaon_FullAnalysis_{KINEMATIC}_{epsilon}_correction_ledger_no_empirical_residual.json")
        require(ledger.get("active_profile") == "no_empirical_residual" and
                ledger.get("particle_type") == "kaon" and
                ledger.get("epsset") == epsset and ledger.get("q2") == "4p4" and
                ledger.get("w") == "2p74" and
                Counter(row.get("phi_setting") for row in ledger.get("settings", [])) ==
                Counter(s["phi"] for s in CANONICAL_SETTINGS if s["epsilon"] == epsilon),
                "ledger_setting_inventory_invalid")
        csv_path = outdir / f"kaon_FullAnalysis_{KINEMATIC}_{epsilon}_correction_ledger_no_empirical_residual.csv"
        with csv_path.open(encoding="utf-8", newline="") as handle:
            totals = [row["phi_setting"] for row in csv.DictReader(handle) if row.get("row_kind") == "setting_total"]
        require(Counter(totals) == Counter(s["phi"] for s in CANONICAL_SETTINGS if s["epsilon"] == epsilon),
                "ledger_csv_setting_inventory_invalid")
        settings.append({**setting, "page_count": count, "renderer_failures": [],
                         "producer_e8_4_authority_and_parent_gates_passed": True})
    return records, settings


def verify_zip(path, source_commit, records):
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_integrity_failed")
        manifest = json.loads(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [], "bundle_incomplete")
        require(manifest.get("git_head") == source_commit and
                manifest.get("required_analysis_commit") == source_commit, "bundle_source_identity_mismatch")
        validate_settings(manifest.get("requested_settings"))
        require(manifest.get("validation_profile") == PROFILE_ID, "bundle_profile_mismatch")
        settings = manifest.get("settings", [])
        validate_settings([{key: row.get(key) for key in ("phi", "epsilon")} for row in settings])
        require(set(manifest.get("global_artifacts", {})) == {"candidate_f2", "candidate_f3", "candidate_f4", "run_summary"},
                "bundle_global_inventory_invalid")
        entries = list(manifest["global_artifacts"].values())
        for row in settings:
            require(set(row.get("artifacts", {})) == {"procedure_pdf", "page_manifest", "full_analysis",
                    "correction_ledger_json", "correction_ledger_csv"}, "bundle_setting_inventory_invalid")
            entries.extend(row["artifacts"].values())
        seen = set()
        for entry in entries:
            name = Path(entry["archive_path"]).name
            require(name in records, "bundle_artifact_inventory_invalid")
            raw = archive.read(entry["archive_path"])
            require(hashlib.sha256(raw).hexdigest() == entry["sha256"] == records[name]["sha256"]
                    and len(raw) == entry["byte_size"] == records[name]["bytes"], "bundle_artifact_hash_mismatch:" + name)
            seen.add(name)  # Epsilon-wide files intentionally occur in each phi directory.
        require(seen == set(records), "bundle_artifact_missing")


class GateStatus(accepted.GateStatus):
    def __init__(self, outdir, output, source_commit):
        self.path = accepted.gate_status_path(outdir, output)
        now = datetime.now(timezone.utc).isoformat()
        self.record = {
            "schema_version": "e8_4_fix5_canonical_five_owner_gate_status/v1",
            "source_commit": source_commit, "kinematic": KINEMATIC,
            "settings": list(CANONICAL_SETTINGS), "expected_zip_path": output.as_posix(),
            "analysis_log_path": (outdir / (output.stem + ".log")).as_posix(),
            "started_at_utc": now, "updated_at_utc": now, "status": "running",
            "stage": "preflight", "failure_reason": None,
            "analysis_started": False, "analysis_completed": False,
            "artifact_verification_completed": False, "collector_source_preflight_completed": False,
            "collection_completed": False, "zip_verification_completed": False}
        self._persist(initial=True)


def copy_companion(source, destination):
    """Exclusive delivery; a partial copy belongs to this failed attempt only."""
    with source.open("rb") as src, destination.open("xb") as dst:
        shutil.copyfileobj(src, dst)
        dst.flush()
        os.fsync(dst.fileno())


def deliver_evidence(status, log_path, summary_path, output):
    destinations = {
        "log": output.with_name(output.stem + ".log"),
        "run_summary": output.with_name(output.stem + "-run-summary.json"),
        "gate_status": output.with_name(output.stem + "-gate-status.json")}
    for destination in destinations.values():
        require(not os.path.lexists(destination), "companion_destination_exists:" + str(destination))
    status.update("companion_delivery")
    identities = {}
    for name, source in (("log", log_path), ("run_summary", summary_path)):
        destination = destinations[name]
        identity = {"path": str(destination), "sha256": sha256(source),
                    "bytes": source.stat().st_size}
        copy_companion(source, destination)
        require(sha256(destination) == identity["sha256"] and
                destination.stat().st_size == identity["bytes"], "companion_copy_identity_mismatch:" + name)
        identities[name] = identity
    status.update("complete", status="success", companion_evidence=identities,
                  delivered_gate_status_path=str(destinations["gate_status"]))
    # Final acceptance operation: transfer the already-published success receipt.
    # Any copy failure propagates to execute_gate's failed source-status handler.
    copy_companion(status.path, destinations["gate_status"])


def execute_gate(repo, outdir, output, source_commit, candidate_materialization_dir,
                 lineage_preflight_only=False):
    repo, outdir, output = Path(repo).resolve(), Path(outdir).resolve(), Path(output).resolve()
    status = GateStatus(outdir, output, source_commit)
    try:
        if lineage_preflight_only:
            status.update(mode="lineage-preflight-only",
                          role="separate lineage gate; not canonical-five runtime validation",
                          canonical_five_runtime_validation=False, analysis_started=False,
                          expected_zip_path=None, analysis_log_path=None,
                          receipt_path=str(status.path))
        require(not output.exists(), "output_zip_already_exists")
        require(output.parent.is_dir(), "output_directory_missing")
        provenance = preflight(repo, source_commit)
        if lineage_preflight_only:
            status.update(ordinary_checkout_before=provenance)
        profile = resolved_profile(repo, source_commit)
        status.update("candidate_materialization_verification")
        materialization = verify_candidate_materialization(candidate_materialization_dir)
        status.update(candidate_materialization=materialization)
        if not lineage_preflight_only:
            before = {name: (path.stat().st_mtime_ns, path.stat().st_size, sha256(path))
                      for name in artifact_names(profile) if (path := outdir / name).is_file()}
            log_path = outdir / (output.stem + ".log")
            require(not log_path.exists(), "run_log_already_exists")
            status.update("collector_source_preflight")
            checks = collector_source_preflight(repo, source_commit, profile)
            status.update(collector_source_preflight_completed=True, collector_source_checks=checks)
        status.update("candidate_staging")
        installation = stage_candidates(candidate_materialization_dir, outdir, materialization, status=status)
        status.update(candidate_installation=installation)
        with owned_worktree(repo, source_commit, "analysis") as worktree:
            status.update("path_isolation")
            env, overlay_record = prepare_runtime_overlay(repo, worktree)
            status.update(ltsep_runtime_overlay=overlay_record)
            probes = probe_paths(worktree, env, overlay_record["baseline_paths"])
            require(outdir == (Path(probes[0]["paths"]["VOLATILEPATH"]) /
                              "OUTPUT/Analysis/KaonLT").resolve(), "analysis_artifact_root_mismatch")
            status.update("external_symlink_preflight")
            external = external_symlink_preflight(worktree, probes[0]["paths"])
            identity = {"path": str(worktree), "head": source_commit, "detached": True,
                        "clean_before_execution": True, "owner_created": True}
            status.update("f6_3_lineage_preflight")
            lineage = lineage_preflight(worktree, outdir, env)
            require(lineage.get("passed") is True and lineage.get("f4_shared_reproduction_passed") is True and
                    [row.get("setting_id") for row in lineage.get("settings", [])] == list(F1_SHA256),
                    "canonical_five_lineage_preflight_failed")
            status.update(f6_3_lineage_preflight=lineage)
            if lineage_preflight_only:
                status.update("worktree_cleanup", analysis_worktree=identity,
                              isolation_probes=probes, external_symlink_preflight=external)
            else:
                started_ns = time.time_ns()
                status.update("analysis", analysis_started=True, analysis_worktree=identity,
                              isolation_probes=probes, external_symlink_preflight=external)
                # Always check the primary checkout after the child, including an
                # exception/failure; never restore it or accept artifacts on mismatch.
                try:
                    returncode = run_analysis(worktree, log_path, env)
                finally:
                    status.update("ordinary_checkout_preservation")
                    try:
                        preservation = verify_preservation(repo, provenance)
                    finally:
                        ltsep_preservation = verify_ltsep_preservation(overlay_record)
                        status.update(installed_ltsep_preservation=ltsep_preservation)
                    status.update(ordinary_checkout_preservation=preservation)
                status.update("completion_markers", analysis_returncode=returncode)
                verify_completion(log_path, returncode)
                status.update(analysis_completed=True)
                status.update("verify_artifacts")
                records, setting_checks = verify_artifacts(outdir, profile, started_ns, before)
                status.update(artifact_verification_completed=True)
                summary = {"schema_version": "e8_4_fix5_canonical_five_run_summary/v1",
                    "source_commit": source_commit, "settings": list(CANONICAL_SETTINGS),
                    "analysis_command": COMMAND, "analysis_worktree": identity,
                    "isolation_probes": probes, "external_symlink_preflight": external,
                    "candidate_materialization": materialization, "candidate_installation": installation,
                    "f6_3_lineage_preflight": lineage,
                    "ltsep_runtime_overlay": overlay_record, "installed_ltsep_preservation": ltsep_preservation,
                    "ordinary_checkout_before": provenance, "ordinary_checkout_preservation": preservation,
                    "analysis_returncode": returncode, "log_sha256": sha256(log_path),
                    "setting_verification": setting_checks, "candidate_sha256": dict(CANDIDATES),
                    "collector_source_checks": checks, "artifacts": records,
                    # A packaged summary cannot contain its enclosing ZIP's own hash.
                    # The atomic success status supplies that final identity, without
                    # rewriting the immutable summary after collector hashing.
                    "zip_identity": {"path": str(output), "success_receipt": str(status.path)},
                    "boundaries": {"method_a": "detached/non-production", "method_b_numerically_excluded": True,
                        "production_promotion": False, "absolute_pion_misid_claim": False,
                        "absolute_simc_amplitude_claim": False}}
                summary_path = outdir / SUMMARY
                summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                records[SUMMARY] = {"sha256": sha256(summary_path), "bytes": summary_path.stat().st_size}
                status.update("collection")
                with tempfile.TemporaryDirectory(prefix="kaonlt-canonical-five-profile-") as directory:
                    effective = Path(directory) / "profile.json"
                    effective.write_text(json.dumps(profile), encoding="utf-8")
                    with owned_worktree(repo, source_commit, "source") as source:
                        with collection_module(source) as module, redirect_stdout(sys.stderr):
                            result = module.collect_validation_bundle(outdir=outdir, kinematic=KINEMATIC,
                                output=output, profile_path=effective, repo_root=source)
                        require(result["returncode"] == 0, "collection_failed")
                        status.update("source_recheck", collection_completed=True)
                        final_checks = collector_source_preflight(repo, source_commit, profile)
                        status.update(final_collector_source_checks=final_checks)
                        verify_preservation(repo, provenance)
                        verify_ltsep_preservation(overlay_record)
                        status.update("verify_zip")
                        verify_zip(output, source_commit, records)
                        status.update("worktree_cleanup", zip_verification_completed=True,
                                      zip_identity={"path": str(output), "sha256": sha256(output),
                                                    "bytes": output.stat().st_size})
        try:
            final_preservation = verify_preservation(repo, provenance)
        finally:
            final_ltsep_preservation = verify_ltsep_preservation(overlay_record)
        status.update(final_ordinary_checkout_preservation=final_preservation,
                      final_installed_ltsep_preservation=final_ltsep_preservation,
                      worktree_cleanup_completed=True)
        if lineage_preflight_only:
            # This atomic status is the directly returnable preflight receipt.
            # All cleanup and acceptance checks precede success publication.
            status.update("success", status="success", failure_reason=None)
            return status.path
        deliver_evidence(status, log_path, summary_path, output)
        return output
    except Exception as exc:
        status.update(status="failed", failure_reason=str(exc))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, required=True, help="active external analysis artifact root")
    parser.add_argument("--output", type=Path, required=True, help="fresh attempt output name; preflight-only returns its status receipt without creating a ZIP")
    parser.add_argument("--candidate-materialization-dir", type=Path, required=True,
                        help="reviewed detached fresh-F1 candidate materialization directory")
    parser.add_argument("--lineage-preflight-only", action="store_true",
                        help="run the separate five-setting lineage gate and return its receipt; no analysis")
    args = parser.parse_args(argv)
    try:
        result = execute_gate(args.repo, args.outdir, args.output, args.source_commit, args.candidate_materialization_dir,
                              lineage_preflight_only=args.lineage_preflight_only)
    except Exception as exc:
        print("Canonical-five gate failed: " + str(exc), file=sys.stderr)
        print("Owner gate status: " + str(accepted.gate_status_path(args.outdir, args.output)), file=sys.stderr)
        return 1
    print(result.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
