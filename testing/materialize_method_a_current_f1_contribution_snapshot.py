"""Own one detached current-F.1 F.2/F.3/F.4 candidate and verified ZIP.

Inputs are five explicit ROLE=PATH values or one directory with the five
canonical F.1 names. Existing external --outdir and --bundle-dir are required.
No historical comparison, canonical equivalence, farm authority or production
promotion is granted. Success prints only a separately labeled local locator
path; the locator is excluded from the transferable four-member ZIP.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import uuid
import zipfile

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT), str(REPO_ROOT / "testing"), str(REPO_ROOT / "src" / "cuts")]
from testing import compare_method_a_current_baseline_authority as shared
import analyze_pion_hgcer_method_a_current_snapshot_contribution_debug as downstream
import pion_hgcer_method_a_acceptance_contract as f1

f2, f3, f4 = shared.f2, shared.f3, shared.f4
debug = downstream.debug
require = debug.require
sha256 = downstream.sha256
json_bytes = downstream.json_bytes
KINEMATIC = debug.KINEMATIC
SETTINGS = debug.SETTINGS
# Explicitly NOT a farm source HEAD. The frozen builder requires 40 hex digits.
DETACHED_CONSTRUCTION_HEAD = "0" * 40
SCHEMA = "method_a_current_f1_contribution_snapshot_materialization/v1"
NAMES = {"f2": "candidate-f2.json", "f3": "candidate-f3.json", "f4": "candidate-f4.json"}
MEMBERS = (*NAMES.values(), "manifest.json")
WRITERS = {"f2": f2.write_pion_hgcer_method_a_acceptance_representation_json,
           "f3": f3.write_pion_hgcer_method_a_acceptance_map_json,
           "f4": f4.write_pion_hgcer_method_a_parent_preserving_correction_json}
BODY = {"f2": "representation", "f3": "acceptance_map", "f4": "correction"}


def strict_json(raw):
    value = downstream.parse_json(raw)
    # parse_constant alone does not reject overflowing JSON numbers (1e999).
    def finite(item):
        if isinstance(item, dict):
            for child in item.values():
                finite(child)
        elif isinstance(item, list):
            for child in item:
                finite(child)
        elif isinstance(item, float):
            require(math.isfinite(item), "nonfinite_json_number")
    finite(value)
    return value


def transferable(value):
    # Public candidate wrappers require an empty input_paths provenance field;
    # aggregate diagnostics prohibit that key even when empty. Preserve the
    # public schema and require the empty value; all event fields stay forbidden.
    def paths(item):
        if isinstance(item, dict):
            for key, child in item.items():
                if key == "input_paths":
                    require(child == {}, "private_input_paths_forbidden")
                else:
                    require(key.lower() not in debug.aggregate.FORBIDDEN,
                            "event_level_persistence_forbidden:" + key)
                paths(child)
        elif isinstance(item, list):
            for child in item:
                paths(child)
        elif isinstance(item, str):
            require(not item.startswith(("/", "\\")) and
                    re.match(r"^[A-Za-z]:[\\/]", item) is None,
                    "private_absolute_path_forbidden")
        else:
            debug.aggregate.guard(item)
    paths(value)


def canonical_name(role):
    phi, epsilon = role.split("-")
    return f1.pion_hgcer_method_a_acceptance_contract_filename(phi, "kaon", KINEMATIC, epsilon)


def input_paths(args):
    if args.f1_dir is not None:
        require(args.f1_dir.is_dir(), "f1_directory_missing")
        return {s: args.f1_dir / canonical_name(s) for s in SETTINGS}
    return shared._f1_inputs(args.f1)


def read_regular(path, role):
    """Check identity before/after reading; symlinks and special files fail."""
    require(not path.is_symlink(), "input_symlink:" + role)
    require(path.exists(), "input_missing:" + role)
    before = path.stat()
    require(stat.S_ISREG(before.st_mode), "input_not_regular:" + role)
    raw = path.read_bytes()
    after = path.stat()
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    require(identity(before) == identity(after), "input_changed_during_read:" + role)
    return raw, identity(after)


def load_f1(paths):
    require(set(paths) == set(SETTINGS), "f1_role_inventory_invalid")
    resolved, artifacts, hashes, identities = {}, [], {}, {}
    for role in SETTINGS:
        path = Path(paths[role])
        require(path.name == canonical_name(role), "f1_filename_invalid:" + role)
        raw, identity = read_regular(path, role)
        path = path.resolve(strict=True)
        require(all(path != old and not path.samefile(old) for old in resolved.values()),
                "f1_duplicate_file:" + role)
        artifact = strict_json(raw)
        setting = artifact.get("setting", {})
        require(setting.get("phi_setting", "") + "-" + setting.get("epsilon_filename_token", "") == role,
                "f1_role_assignment_invalid:" + role)
        require(setting.get("kinematic_token") == KINEMATIC, "f1_kinematic_invalid:" + role)
        resolved[role], hashes[role], identities[role] = path, sha256(raw), identity
        artifacts.append(artifact)
    f3._validate_f1_artifacts(artifacts)
    debug.f5._common_geometry(artifacts)
    return artifacts, hashes, resolved, identities


def recheck_f1(paths, hashes, identities):
    for role in SETTINGS:
        raw, identity = read_regular(paths[role], role)
        require(identity == identities[role] and sha256(raw) == hashes[role],
                "f1_changed_during_materialization:" + role)


def external_root(path):
    require(path.is_dir() and not path.is_symlink(), "external_root_missing_or_symlink")
    path = path.resolve(strict=True)
    ordinary = REPO_ROOT.resolve()
    require(path != ordinary and ordinary not in path.parents, "output_root_inside_checkout")
    return path


def construction_reference(artifact, raw_sha):
    body = artifact["acceptance_map"]
    return {"source_file_sha256": raw_sha, "map_fingerprint": body["fingerprint"],
            "algorithm_fingerprint": body["algorithm_fingerprint"],
            "artifact_fingerprint": artifact["artifact_fingerprint"],
            "farm_source_head": DETACHED_CONSTRUCTION_HEAD}


def available(artifact, stage):
    body = artifact[BODY[stage]]
    require(body.get("available") is True and body.get("status") == "available",
            stage + "_unavailable:" + str(body.get("invalid_parents", body.get("reason"))))
    transferable(artifact)
    require(artifact["provenance"]["input_paths"] == {}, stage + "_private_input_paths")


def serialize(stage, artifact, attempt):
    available(artifact, stage)
    path = attempt / NAMES[stage]
    require(not path.exists() and not path.is_symlink(), "candidate_output_exists:" + stage)
    WRITERS[stage](path, artifact)
    raw = path.read_bytes()
    require(raw == shared._writer_bytes(artifact), "public_writer_bytes_mismatch:" + stage)
    require(shared.first_mismatch(strict_json(raw), artifact) is None,
            "candidate_readback_mismatch:" + stage)
    return raw


def verify_stages(artifacts, hashes, raw):
    """Independent persisted-byte validators and exact frozen F.4 reproduction."""
    require(set(raw) == set(NAMES), "candidate_stage_inventory_invalid")
    candidates = {s: strict_json(data) for s, data in raw.items()}
    for stage, artifact in candidates.items():
        available(artifact, stage)
    parsed = f3._validate_f1_artifacts(artifacts)
    debug.f5._common_geometry(artifacts)
    representation, groups = f3._validate_f2_artifact(candidates["f2"], parsed,
                                                    f3._normalize_f1_hashes(parsed, hashes))
    map_value, models = f4._validate_f3_artifact(candidates["f3"], parsed, hashes, 1e-12)
    require(len(groups) == len(models) == 15, "f2_f3_parent_inventory_invalid")
    for key, expected in {"f2_source_file_sha256": sha256(raw["f2"]),
                          "f2_representation_fingerprint": representation["fingerprint"],
                          "f2_algorithm_fingerprint": representation["algorithm_fingerprint"]}.items():
        require(map_value[key] == expected, "f3_f2_binding_mismatch:" + key)
    reference = construction_reference(candidates["f3"], sha256(raw["f3"]))
    body = debug.validate_f4_wrapper(candidates["f4"])
    require(shared.first_mismatch(body["f3_runtime_authority"]["accepted"], reference) is None,
            "detached_f3_construction_reference_mismatch")
    # No canonical authority map is invoked. This verifies all 15 signed/source/
    # phi closures and positive finite transient factors using the frozen builder.
    snapshot, _discarded_aggregates = debug.reproduce_snapshot(
        artifacts, candidates["f3"], candidates["f4"],
        hashes={**hashes, "f3": sha256(raw["f3"]), "f4": sha256(raw["f4"])})
    return candidates, snapshot, reference


def source_observation():
    def git(*args):
        return subprocess.check_output(["git", "--no-optional-locks", *args], cwd=REPO_ROOT,
                                       text=True).rstrip("\r\n")
    head = git("rev-parse", "HEAD")
    require(re.fullmatch("[0-9a-f]{40}", head) is not None, "source_head_invalid")
    # Record status without transferable machine/user-specific filenames.
    status = git("status", "--porcelain=v1", "--untracked-files=all")
    counts = {}
    for row in status.splitlines():
        counts[row[:2]] = counts.get(row[:2], 0) + 1
    return {"git_head": head, "worktree_dirty": bool(status), "status_counts": counts}


def manifest_for(artifacts, hashes, raw, source):
    candidates, snapshot, reference = verify_stages(artifacts, hashes, raw)
    parsed = f3._validate_f1_artifacts(artifacts)
    rows = {p["setting_id"]: {"raw_sha256": hashes[p["setting_id"]],
             "stable_content_fingerprint": p["stable_content_fingerprint"],
             "contract_fingerprint": p["fingerprints"]["fingerprint"]} for p in parsed}
    outputs = {s: {"basename": NAMES[s], "raw_sha256": sha256(raw[s]),
                   "body_fingerprint": a[BODY[s]]["fingerprint"],
                   "artifact_fingerprint": a["artifact_fingerprint"]}
               for s, a in candidates.items()}
    for s in ("f2", "f3"):
        outputs[s]["algorithm_fingerprint"] = candidates[s][BODY[s]]["algorithm_fingerprint"]
    value = {"schema_version": SCHEMA, "kinematic_token": KINEMATIC,
             "role": "detached candidate materialization; no accepted scientific authority",
             "non_authoritative": True, "internal_snapshot_consistency": "PASS",
             "canonical_owner_equivalence": "BLOCKED", "canonical_owner_equivalence_evaluated": False,
             "method_a_promoted": False, "production_objects_mutated": False,
             "production_application_performed": False, "method_b_numerical_dependency": False,
             "source": source, "f1_inputs": rows, "candidate_outputs": outputs,
             "snapshot": snapshot, "detached_f3_construction_reference": {
                 "purpose": "internal_reproduction_only", "accepted_farm_authority": False,
                 "sentinel_head": DETACHED_CONSTRUCTION_HEAD, "record": reference},
             "checks": {"exact_f4_reproduction": True, "all_15_signed_parent_closures": True,
                        "cross_stage_bindings": True, "persisted_readback": True},
             "local_locator_in_transferable_zip": False}
    transferable(value)
    return debug.seal(value)


def verify_zip(path, artifacts, hashes, source, expected):
    with zipfile.ZipFile(path) as archive:
        require(archive.namelist() == list(MEMBERS) and archive.testzip() is None,
                "zip_inventory_or_crc_invalid")
        members = {name: archive.read(name) for name in MEMBERS}
    require(members == expected, "zip_member_bytes_mismatch")
    raw = {s: members[name] for s, name in NAMES.items()}
    persisted = strict_json(members["manifest.json"])
    require(shared.first_mismatch(persisted, manifest_for(artifacts, hashes, raw, source)) is None,
            "zip_manifest_or_fingerprint_mismatch")
    return sha256(path.read_bytes())


def create_zip(path, members):
    require(set(members) == set(MEMBERS), "zip_member_inventory_invalid")
    with zipfile.ZipFile(path, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in MEMBERS:
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, members[name])


def run(args):
    require(args.kinematic == KINEMATIC, "kinematic_unsupported")
    outroot, bundleroot = external_root(args.outdir), external_root(args.bundle_dir)
    artifacts, hashes, paths, identities = load_f1(input_paths(args))
    require(all(outroot != p and bundleroot != p for p in paths.values()), "output_input_alias")
    source = source_observation()
    token = uuid.uuid4().hex
    require(re.fullmatch(r"[0-9a-z]+", token) is not None, "attempt_token_invalid")
    attempt = outroot / ("kaonlt-current-f1-snapshot-" + token)
    final = bundleroot / (attempt.name + ".zip")
    require(not attempt.exists() and not attempt.is_symlink() and
            not final.exists() and not final.is_symlink(), "output_already_exists")
    # Exclusive ownership: collision failures never remove somebody else's files.
    attempt.mkdir()
    zip_created = False
    try:
        raw = {}
        candidate2 = f2.build_pion_hgcer_method_a_acceptance_representation_artifact(
            artifacts, input_file_hashes=hashes, input_paths={})
        raw["f2"] = serialize("f2", candidate2, attempt)
        candidate3 = f3.build_pion_hgcer_method_a_acceptance_map_artifact(
            artifacts, strict_json(raw["f2"]), f1_input_file_hashes=hashes,
            f2_input_file_sha256=sha256(raw["f2"]), input_paths={})
        raw["f3"] = serialize("f3", candidate3, attempt)
        reference = construction_reference(strict_json(raw["f3"]), sha256(raw["f3"]))
        candidate4 = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(
            artifacts, strict_json(raw["f3"]), f1_input_file_hashes=hashes,
            f3_input_file_sha256=sha256(raw["f3"]), input_paths={},
            accepted_f3_runtime_authority_by_kinematic={KINEMATIC: reference})
        raw["f4"] = serialize("f4", candidate4, attempt)
        manifest = manifest_for(artifacts, hashes, raw, source)
        members = {NAMES[s]: data for s, data in raw.items()}
        members["manifest.json"] = json_bytes(manifest)
        with (attempt / "manifest.json").open("xb") as handle:
            handle.write(members["manifest.json"])
        # Verify the persisted files, not just construction objects.
        require({n: (attempt/n).read_bytes() for n in MEMBERS} == members, "persisted_bytes_mismatch")
        staged_zip = attempt / "verification.zip"
        create_zip(staged_zip, members)
        expected_zip_sha = verify_zip(staged_zip, artifacts, hashes, source, members)
        recheck_f1(paths, hashes, identities)
        with final.open("xb") as target:
            zip_created = True
            with staged_zip.open("rb") as incoming:
                shutil.copyfileobj(incoming, target)
        require(final.resolve().parent == bundleroot and not final.is_symlink(), "final_zip_path_invalid")
        require(verify_zip(final, artifacts, hashes, source, members) == expected_zip_sha,
                "final_zip_sha256_mismatch")
        locator = attempt / "LOCAL_FARM_ONLY_inputs.json"
        value = {"f1": {s: str(p) for s, p in paths.items()},
                 "f3": str(attempt / NAMES["f3"]), "f4": str(attempt / NAMES["f4"])}
        with locator.open("xb") as handle:
            handle.write(json_bytes(value))
        require(strict_json(locator.read_bytes()) == value, "local_locator_readback_mismatch")
        require({n: (attempt/n).read_bytes() for n in MEMBERS} == members, "final_candidate_bytes_mismatch")
        staged_zip.unlink()
        # Last gate immediately before success publication; original inputs never copied.
        recheck_f1(paths, hashes, identities)
        return locator
    except Exception:
        if zip_created:
            final.unlink()
        # Only the exclusively created attempt is owned. No ordinary checkout cleanup.
        shutil.rmtree(attempt)
        raise


def parser():
    value = argparse.ArgumentParser(description=__doc__)
    inputs = value.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--f1-dir", type=Path)
    inputs.add_argument("--f1", action="append", help="ROLE=PATH, exactly five canonical roles")
    value.add_argument("--kinematic", required=True, choices=[KINEMATIC])
    value.add_argument("--outdir", type=Path, required=True)
    value.add_argument("--bundle-dir", type=Path, required=True)
    return value


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args)
    except (OSError, subprocess.SubprocessError, zipfile.BadZipFile) as exc:
        print("materialization_failed:" + type(exc).__name__, file=sys.stderr)
        return 1
    except (ValueError, KeyError, TypeError, IndexError, OverflowError) as exc:
        print("materialization_failed:" + str(exc), file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
