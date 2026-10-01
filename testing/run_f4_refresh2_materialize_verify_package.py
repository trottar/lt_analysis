"""Own the detached F.4.Refresh.2 materialize, verify, and package sequence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from testing import collect_pion_hgcer_validation_bundle as collector  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[1]
PROFILE_RELATIVE = "testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json"
MATERIALIZER = "testing/materialize_method_a_current_baseline_authority.py"
WRAPPER = "testing/package_pion_hgcer_validation_bundle.tcsh"
SOURCE_HEAD = "141a3d04f9e5d07be21dba14e0e63212c3990bf1"
COMPARISON_SHA256 = "c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5"
KINEMATIC = "Q4p4W2p74"
PROFILE_ID = "phase_f4_refresh2_current_baseline_candidate_materialization_farm_review/v1"
MATERIALIZATION_SCHEMA = "method_a_current_baseline_authority_materialization/v1"
ARTIFACT_ROOT = Path("/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT")
GLOBUS_ROOT = Path("/volatile/hallc/c-kaonlt/trottar/globus")
ALIASES = ("Left-lowe", "Left-highe", "Center-lowe", "Center-highe", "Right-highe")
SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
            ("Center", "highe"), ("Right", "highe"))
OUTPUT_SUFFIXES = {
    "f4_refresh1_comparison_input": "current-baseline-authority-comparison-input.json",
    "f4_refresh2_candidate_f2": "acceptance-representation-current-baseline-candidate.json",
    "f4_refresh2_candidate_f3": "acceptance-map-current-baseline-candidate.json",
    "f4_refresh2_candidate_f4": "parent-preserving-correction-current-baseline-candidate.json",
    "f4_refresh2_materialization_manifest": "current-baseline-authority-materialization-manifest.json",
}
STAGES = {
    "f2": "f4_refresh2_candidate_f2",
    "f3": "f4_refresh2_candidate_f3",
    "f4": "f4_refresh2_candidate_f4",
}
ALLOWED_COMMITTED = frozenset({
    PROFILE_RELATIVE,
    "testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py",
    "testing/run_f4_refresh2_materialize_verify_package.py",
    "testing/test_run_f4_refresh2_materialize_verify_package.py",
    "testing/test_memory_health.py",
    "tools/check_memory_health.py",
})


def _require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reject_constant(value: str) -> None:
    raise ValueError(f"nonstandard_json_constant:{value}")


def _unique_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate_json_key:{key}")
        result[key] = value
    return result


def strict_json(raw: bytes) -> dict:
    value = json.loads(raw.decode("utf-8"), parse_constant=_reject_constant,
                       object_pairs_hook=_unique_pairs)
    _require(isinstance(value, dict), "json_root_not_object")
    return value


def command_runner(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)


def _command(runner, command: list[str], cwd: Path, reason: str) -> str:
    result = runner(command, cwd)
    _require(result.returncode == 0,
             f"{reason}:exit={result.returncode}:{(result.stderr or result.stdout).strip()}")
    return result.stdout.strip()


def _f1_inputs(rows: list[str]) -> dict[str, Path]:
    _require(len(rows) == len(ALIASES), "canonical_five_f1_inputs_required")
    found: dict[str, Path] = {}
    for row in rows:
        alias, separator, name = row.partition("=")
        _require(separator == "=" and alias in ALIASES and alias not in found and bool(name),
                 f"invalid_or_duplicate_f1_alias:{alias}")
        found[alias] = Path(name)
    _require(set(found) == set(ALIASES), "canonical_five_f1_aliases_incomplete")
    return found


def _zip_path(output: str, globus_root: Path) -> Path:
    path = Path(output)
    if path.is_absolute():
        _require(path.parent == globus_root, "zip_must_be_in_canonical_globus_directory")
    else:
        _require(path.name == output, "zip_must_be_a_basename_or_canonical_absolute_path")
        path = globus_root / path
    _require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.zip", path.name) is not None,
             "unsafe_or_non_zip_output_name")
    return path


def _profile_outputs(profile: dict, artifact_dir: Path) -> dict[str, Path]:
    _require(profile.get("schema_version") == "pion_hgcer_validation_bundle_profile/v4"
             and profile.get("validation_profile") == PROFILE_ID
             and profile.get("collection_mode") == "generic_artifacts",
             "wrong_f4_refresh2_profile_identity")
    identity = profile.get("source_identity")
    _require(isinstance(identity, dict) and identity.get("required_analysis_commit") == SOURCE_HEAD,
             "wrong_required_analysis_commit")
    _require(set(identity.get("allowed_committed_files", [])) == ALLOWED_COMMITTED
             and identity.get("allowed_non_analysis_path_prefixes") == ["docs/memory/"],
             "profile_source_identity_not_narrow")
    _require(tuple((row.get("phi"), row.get("epsilon")) for row in profile.get("settings", []))
             == SETTINGS, "canonical_five_profile_settings_changed")
    _require(profile.get("artifacts", {}).get("settings") == [], "unexpected_setting_artifacts")
    rows = profile["artifacts"]["global"]
    _require(isinstance(rows, list) and len(rows) == 5, "five_global_artifacts_required")
    outputs: dict[str, Path] = {}
    for row in rows:
        _require(isinstance(row, dict), "global_artifact_declaration_invalid")
        key = row.get("key")
        _require(key in OUTPUT_SUFFIXES and key not in outputs and row.get("kind") == "json"
                 and row.get("required") is True, "global_artifact_declaration_invalid")
        expected = "{kinematic}_kaon_pion-background_hgcer-method-a-" + OUTPUT_SUFFIXES[key]
        _require(row.get("basename_template") == expected, "global_artifact_basename_changed")
        outputs[key] = artifact_dir / expected.format(kinematic=KINEMATIC)
    _require(set(outputs) == set(OUTPUT_SUFFIXES), "five_global_artifacts_required")
    return outputs


def preflight(args, *, repo_root: Path = REPO_ROOT, artifact_root: Path = ARTIFACT_ROOT,
              globus_root: Path = GLOBUS_ROOT, runner=command_runner,
              profile_path: Path | None = None) -> tuple[dict, dict[str, Path], dict[str, Path], Path]:
    _require(re.fullmatch(r"[0-9a-f]{40}", args.bundle_commit) is not None,
             "bundle_commit_must_be_full_lowercase_sha40")
    _require(args.kinematic == KINEMATIC, "unsupported_kinematic")
    _require(repo_root.is_dir() and (repo_root / ".git").exists(), "kaonlt_repository_missing")
    top = _command(runner, ["git", "rev-parse", "--show-toplevel"], repo_root, "git_top_level_failed")
    _require(Path(top).resolve() == repo_root.resolve(), "not_kaonlt_repository_root")
    head = _command(runner, ["git", "rev-parse", "HEAD"], repo_root, "git_head_failed")
    _require(head == args.bundle_commit, "bundle_commit_does_not_match_head")
    status = _command(runner, ["git", "status", "--porcelain=v1", "--untracked-files=all"],
                      repo_root, "git_status_failed")
    _require(not status, "dirty_repository_worktree")
    _command(runner, ["git", "merge-base", "--is-ancestor", SOURCE_HEAD, args.bundle_commit],
             repo_root, "materializer_source_not_ancestor")
    path = profile_path or repo_root / PROFILE_RELATIVE
    _require(path.is_file(), "f4_refresh2_profile_missing")
    profile = collector.load_validation_profile(path)
    _require(Path(args.artifact_dir).resolve() == artifact_root.resolve() and artifact_root.is_dir(),
             "artifact_dir_must_be_canonical_root")
    outputs = _profile_outputs(profile, artifact_root)
    changed = _command(runner, ["git", "diff", "--name-only", f"{SOURCE_HEAD}..{args.bundle_commit}"],
                       repo_root, "committed_range_diff_failed")
    unexpected = sorted({line for line in changed.splitlines() if line and line not in ALLOWED_COMMITTED
                         and not line.startswith("docs/memory/")})
    _require(not unexpected, f"unexpected_committed_paths:{','.join(unexpected)}")
    _require(Path(args.python).is_absolute() and Path(args.python).is_file()
             and os.access(args.python, os.X_OK), "python_unavailable_or_not_absolute")
    f1 = _f1_inputs(args.f1)
    supplied = [*f1.values(), args.accepted_f2, args.accepted_f3, args.accepted_f4, args.comparison]
    _require(all(Path(item).is_absolute() for item in supplied), "input_paths_must_be_absolute")
    _require(all(Path(item).is_file() for item in supplied), "required_input_missing_or_not_file")
    _require(len({Path(item).resolve() for item in supplied}) == len(supplied), "input_paths_not_distinct")
    _require(sha256_file(Path(args.comparison)) == COMPARISON_SHA256, "reviewed_comparison_sha_mismatch")
    zip_path = _zip_path(args.output, globus_root)
    _require(zip_path.parent.is_dir() and not (zip_path.exists() or zip_path.is_symlink()),
             "zip_destination_unavailable_or_exists")
    _require(all(not (item.exists() or item.is_symlink()) for item in outputs.values()),
             "stale_or_partial_materialization_output_exists")
    _require(all(item.resolve() not in {Path(source).resolve() for source in supplied}
                 for item in outputs.values()), "output_collides_with_input")
    return profile, outputs, f1, zip_path


def verify_materialization(outputs: dict[str, Path], comparison: Path) -> dict[str, str]:
    _require(set(outputs) == set(OUTPUT_SUFFIXES)
             and all(path.is_file() for path in outputs.values()), "materialization_output_missing")
    hashes = {key: sha256_file(path) for key, path in outputs.items()}
    _require(hashes["f4_refresh1_comparison_input"] == COMPARISON_SHA256,
             "comparison_copy_sha_mismatch")
    manifest = strict_json(outputs["f4_refresh2_materialization_manifest"].read_bytes())
    _require(manifest.get("schema_version") == MATERIALIZATION_SCHEMA, "materialization_schema_mismatch")
    for key, expected in {
        "non_authoritative": True, "accepted_authority_mutated": False,
        "production_objects_mutated": False, "production_application_performed": False,
        "method_a_promoted": False, "complete": True,
    }.items():
        _require(manifest.get(key) is expected, f"materialization_flag_mismatch:{key}")
    _require(manifest.get("errors") == [], "materialization_errors_present")
    _require(manifest.get("source_head") == SOURCE_HEAD, "materialization_source_head_mismatch")
    _require(manifest.get("kinematic_token") == KINEMATIC, "materialization_kinematic_mismatch")
    comparison_row = manifest.get("comparison_input")
    _require(isinstance(comparison_row, dict)
             and comparison_row.get("source_basename") == comparison.name
             and comparison_row.get("raw_sha256") == COMPARISON_SHA256
             and comparison_row.get("copied_output_basename") == outputs["f4_refresh1_comparison_input"].name
             and comparison_row.get("copied_output_raw_sha256") == COMPARISON_SHA256,
             "materialization_comparison_provenance_mismatch")
    gate = manifest.get("scientific_gate")
    _require(isinstance(gate, dict) and all(gate.get(key) is value for key, value in {
        "f2_scientific_payload_match": True, "f3_scientific_payload_match": True,
        "f4_scientific_payload_match": False,
    }.items()) and gate.get("first_changed_stage") == "F4", "materialization_scientific_gate_mismatch")
    candidates = manifest.get("candidate_outputs")
    _require(isinstance(candidates, dict) and set(candidates) == set(STAGES),
             "materialization_candidate_inventory_mismatch")
    for stage, key in STAGES.items():
        row = candidates[stage]
        _require(isinstance(row, dict) and row.get("basename") == outputs[key].name
                 and row.get("raw_sha256") == hashes[key]
                 and row.get("non_authoritative") is True,
                 f"materialization_candidate_hash_or_name_mismatch:{stage}")
    return hashes


def verify_zip(zip_path: Path, outputs: dict[str, Path],
               hashes: dict[str, str], bundle_commit: str) -> str:
    _require(zip_path.is_file(), "returned_zip_missing")
    with zipfile.ZipFile(zip_path) as archive:
        _require(archive.testzip() is None, "returned_zip_crc_failure")
        names = archive.namelist()
        _require(names.count("manifest.json") == 1, "returned_zip_manifest_missing_or_duplicate")
        manifest = strict_json(archive.read("manifest.json"))
        _require(manifest.get("complete") is True, "collector_bundle_incomplete")
        _require(manifest.get("validation_profile") == PROFILE_ID
                 and manifest.get("required_analysis_commit") == SOURCE_HEAD
                 and manifest.get("git_head") == bundle_commit
                 and manifest.get("requested_kinematic") == KINEMATIC,
                 "collector_bundle_identity_mismatch")
        _require(manifest.get("errors") == [], "collector_bundle_errors_present")
        records = manifest.get("global_artifacts")
        _require(isinstance(records, dict) and set(records) == set(outputs),
                 "collector_global_artifact_inventory_mismatch")
        for key, path in outputs.items():
            row = records[key]
            archive_path = "global/" + path.name
            _require(isinstance(row, dict) and row.get("status") == "exists"
                     and row.get("json_status") == "valid"
                     and row.get("archive_path") == archive_path
                     and row.get("sha256") == hashes[key]
                     and names.count(archive_path) == 1,
                     f"collector_global_artifact_mismatch:{key}")
            _require(hashlib.sha256(archive.read(archive_path)).hexdigest() == hashes[key],
                     f"packaged_artifact_hash_mismatch:{key}")
    return sha256_file(zip_path)


def run_operation(args, *, repo_root: Path = REPO_ROOT, artifact_root: Path = ARTIFACT_ROOT,
                  globus_root: Path = GLOBUS_ROOT, runner=command_runner,
                  profile_path: Path | None = None) -> tuple[Path, str]:
    _profile, outputs, f1, zip_path = preflight(
        args, repo_root=repo_root, artifact_root=artifact_root, globus_root=globus_root,
        runner=runner, profile_path=profile_path,
    )
    materializer = [args.python, MATERIALIZER]
    for alias in ALIASES:
        materializer.extend(["--f1", f"{alias}={f1[alias]}"])
    for stage in STAGES:
        materializer.extend([f"--accepted-{stage}", str(getattr(args, f"accepted_{stage}"))])
    materializer.extend(["--comparison", str(args.comparison),
                         "--expected-comparison-sha256", COMPARISON_SHA256,
                         "--source-head", SOURCE_HEAD, "--output-dir", str(artifact_root)])
    _command(runner, materializer, repo_root, "materializer_failed")
    hashes = verify_materialization(outputs, Path(args.comparison))
    package = ["tcsh", WRAPPER, "--python", args.python,
               "--bundle-commit", args.bundle_commit, "--profile", PROFILE_RELATIVE,
               "--artifact-dir", str(artifact_root), "--kinematic", KINEMATIC,
               "--output", args.output]
    for key in OUTPUT_SUFFIXES:
        package.extend(["--immutable", str(outputs[key]), hashes[key]])
    _command(runner, package, repo_root, "package_wrapper_failed")
    return zip_path, verify_zip(zip_path, outputs, hashes, args.bundle_commit)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True)
    parser.add_argument("--bundle-commit", required=True)
    parser.add_argument("--f1", action="append", required=True)
    for stage in STAGES:
        parser.add_argument(f"--accepted-{stage}", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--kinematic", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        zip_path, zip_sha = run_operation(args)
    except (OSError, ValueError, UnicodeError, json.JSONDecodeError, zipfile.BadZipFile,
            zipfile.LargeZipFile, KeyError, TypeError) as exc:
        parser.exit(1, f"STOP: F.4.Refresh.2 execution owner failed: {exc}\n")
    print(f"review ZIP: {zip_path}")
    print(f"review ZIP SHA-256: {zip_sha}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
