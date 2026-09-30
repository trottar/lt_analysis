"""Materialize a detached current-baseline F.2/F.3/F.4 candidate lineage."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from testing import compare_method_a_current_baseline_authority as comparison  # noqa: E402


SCHEMA = "method_a_current_baseline_authority_materialization/v1"
PREFIX = "_kaon_pion-background_hgcer-method-a-"
SUFFIXES = {
    "comparison": "current-baseline-authority-comparison-input.json",
    "f2": "acceptance-representation-current-baseline-candidate.json",
    "f3": "acceptance-map-current-baseline-candidate.json",
    "f4": "parent-preserving-correction-current-baseline-candidate.json",
    "manifest": "current-baseline-authority-materialization-manifest.json",
}
BUILDERS = (
    ("f2", "representation", comparison.f2, comparison.F2_PROVENANCE),
    ("f3", "acceptance_map", comparison.f3, comparison.F3_PROVENANCE),
    ("f4", "correction", comparison.f4, comparison.F4_PROVENANCE),
)
WRITER_NAMES = {
    "f2": "write_pion_hgcer_method_a_acceptance_representation_json",
    "f3": "write_pion_hgcer_method_a_acceptance_map_json",
    "f4": "write_pion_hgcer_method_a_parent_preserving_correction_json",
}


def _hex(value, length, label):
    if not isinstance(value, str) or re.fullmatch(rf"[0-9a-fA-F]{{{length}}}", value) is None:
        raise ValueError(f"{label} must be {length} hex characters")
    return value.lower()


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _mapping(value, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _comparison_section(payload, name):
    return _mapping(payload.get(name), f"comparison.{name}")


def _stage_fingerprints(stage, body):
    if stage == "f2":
        return {"representation_fingerprint": body["fingerprint"],
                "algorithm_fingerprint": body["algorithm_fingerprint"]}
    if stage == "f3":
        return {"map_fingerprint": body["fingerprint"],
                "algorithm_fingerprint": body["algorithm_fingerprint"],
                "accepted_basis": body["accepted_basis"]}
    return {"correction_fingerprint": body["fingerprint"]}


def _scientific_result(accepted, candidate, body_name, excluded):
    accepted_body = comparison._body(accepted, body_name)
    candidate_body = comparison._body(candidate, body_name)
    mismatch = comparison.first_mismatch(
        comparison.scientific_projection(accepted_body, excluded),
        comparison.scientific_projection(candidate_body, excluded),
    )
    return {
        "scientific_payload_match": mismatch is None,
        "first_mismatch_path": None if mismatch is None else mismatch["path"],
        "first_mismatch": mismatch,
    }


def _f1_inventory(paths, expected):
    expected_rows = expected.get("f1_inputs")
    _require(isinstance(expected_rows, list) and len(expected_rows) == len(comparison.ALIASES), "comparison F.1 inventory invalid")
    expected_by_id = {}
    for row in expected_rows:
        row = _mapping(row, "comparison F.1 row")
        sid = row.get("setting_id")
        _require(isinstance(sid, str) and sid in comparison.ALIASES and sid not in expected_by_id, "comparison F.1 setting inventory invalid")
        expected_by_id[sid] = row
    _require(set(expected_by_id) == set(comparison.ALIASES), "comparison F.1 setting inventory incomplete")
    artifacts, hashes, provisional = [], {}, []
    for alias in comparison.ALIASES:
        artifact, sha = comparison._read(paths[alias])
        setting = _mapping(artifact.get("setting"), "F.1 setting")
        _require(f"{setting.get('phi_setting')}-{setting.get('epsilon_filename_token')}" == alias, f"F.1 alias/artifact mismatch: {alias}")
        _require(_hex(expected_by_id[alias].get("source_file_sha256"), 64, "comparison F.1 SHA") == sha, f"current F.1 SHA mismatch: {alias}")
        artifacts.append(artifact)
        hashes[alias] = sha
        provisional.append({"alias": alias, "setting_id": alias, "source_file_sha256": sha})
    kinematics = {artifact["setting"].get("kinematic_token") for artifact in artifacts}
    _require(len(kinematics) == 1 and isinstance(next(iter(kinematics)), str), "current F.1 kinematics differ")
    kinematic = next(iter(kinematics))
    _require(kinematic and re.fullmatch(r"[A-Za-z0-9_]+", kinematic) is not None, "kinematic token unsafe for output filename")
    for row in expected_rows:
        _require(_mapping(row.get("setting"), "comparison F.1 setting").get("kinematic_token") == kinematic, "comparison/current F.1 kinematic mismatch")
    return artifacts, hashes, provisional, expected_by_id, kinematic


def _complete_f1_inventory(artifacts, provisional, expected_by_id, candidate_f2):
    f2_rows = comparison._body(candidate_f2, "representation").get("input_fingerprints")
    _require(isinstance(f2_rows, list) and len(f2_rows) == len(comparison.ALIASES), "candidate F.2 F.1 fingerprints invalid")
    result = []
    for index, (item, fp) in enumerate(zip(provisional, f2_rows)):
        fp = _mapping(fp, "candidate F.2 F.1 fingerprint")
        alias = item["setting_id"]
        _require(fp.get("setting_id") == alias and fp.get("source_file_sha256") == item["source_file_sha256"], "candidate F.2 F.1 fingerprint inventory mismatch")
        contract = _mapping(artifacts[index].get("contract"), "F.1 contract")
        training = contract.get("method_a_training_records")
        application = contract.get("application_records")
        _require(isinstance(training, list) and isinstance(application, list), "F.1 record populations invalid")
        row = {**item, "setting": fp["setting"],
               "stable_f1_content_fingerprint": fp["stable_f1_content_fingerprint"],
               "f1_contract_fingerprint": fp["fingerprint"],
               "training_record_count": len(training), "application_record_count": len(application),
               "training_population_fingerprint": fp["method_a_training_population_fingerprint"],
               "application_population_fingerprint": fp["application_population_fingerprint"]}
        _require(row == expected_by_id[alias], f"current F.1 fingerprint/provenance mismatch: {alias}")
        result.append(row)
    return result


def _accepted_inputs(args, expected):
    expected_shas = _comparison_section(expected, "accepted_file_sha256")
    artifacts, manifest_rows = {}, {}
    for stage in ("f2", "f3", "f4"):
        path = getattr(args, f"accepted_{stage}")
        artifact, sha = comparison._read(path)
        _require(sha == _hex(expected_shas.get(stage), 64, f"comparison accepted {stage} SHA"), f"accepted {stage} raw SHA mismatch")
        body_name = {"f2": "representation", "f3": "acceptance_map", "f4": "correction"}[stage]
        body = comparison._body(artifact, body_name)
        _require(artifact.get("non_authoritative") is True and isinstance(artifact.get("artifact_fingerprint"), str), f"accepted {stage} artifact identity invalid")
        artifacts[stage] = artifact
        manifest_rows[stage] = {"input_basename": path.name, "raw_sha256": sha,
                                "artifact_fingerprint": artifact["artifact_fingerprint"],
                                **_stage_fingerprints(stage, body)}
    return artifacts, manifest_rows


def _names(kinematic):
    names = {key: kinematic + PREFIX + suffix for key, suffix in SUFFIXES.items()}
    _require(len(set(names.values())) == 5, "candidate output names collide")
    return names


def _preflight(args):
    paths = comparison._f1_inputs(args.f1)
    inputs = [*paths.values(), args.accepted_f2, args.accepted_f3, args.accepted_f4, args.comparison]
    _require(all(path.is_file() for path in inputs), "every input must exist as a file")
    resolved = [path.resolve() for path in inputs]
    _require(len(set(resolved)) == len(resolved), "input paths must be distinct")
    _require(args.output_dir.is_dir(), "output directory must already exist")
    expected_sha = _hex(args.expected_comparison_sha256, 64, "expected comparison SHA")
    source_head = _hex(args.source_head, 40, "source head")
    expected, comparison_sha = comparison._read(args.comparison)
    _require(comparison_sha == expected_sha, "comparison raw SHA mismatch")
    _require(expected.get("schema_version") == "method_a_current_baseline_authority_comparison/v1" and expected.get("non_authoritative") is True, "comparison schema/non-authoritative status invalid")
    artifacts, hashes, provisional, expected_f1, kinematic = _f1_inventory(paths, expected)
    names = _names(kinematic)
    targets = {key: args.output_dir / name for key, name in names.items()}
    for key, target in targets.items():
        _require(target.name not in {path.name for path in (args.accepted_f2, args.accepted_f3, args.accepted_f4)}, f"{key} candidate basename collides with accepted artifact")
        _require(target.resolve() not in resolved, f"{key} output resolves to an input")
        _require(args.overwrite or not (target.exists() or target.is_symlink()), f"output already exists: {target}")
        _require(not target.exists() or target.is_file(), f"output is not a file: {target}")
    accepted, accepted_rows = _accepted_inputs(args, expected)
    return artifacts, hashes, provisional, expected_f1, kinematic, names, targets, expected, comparison_sha, source_head, accepted, accepted_rows


def _build(args, preflight):
    artifacts, hashes, provisional, expected_f1, kinematic, names, targets, expected, comparison_sha, source_head, accepted, accepted_rows = preflight
    expected_candidate_shas = _comparison_section(expected, "candidate_serialized_sha256")
    candidate_f2 = comparison.f2.build_pion_hgcer_method_a_acceptance_representation_artifact(artifacts, input_file_hashes=hashes, input_paths={})
    f2_bytes = comparison._writer_bytes(candidate_f2)
    f2_sha = hashlib.sha256(f2_bytes).hexdigest()
    _require(f2_sha == _hex(expected_candidate_shas.get("f2"), 64, "comparison candidate F.2 SHA"), "candidate F.2 writer SHA mismatch")
    f1_rows = _complete_f1_inventory(artifacts, provisional, expected_f1, candidate_f2)
    f2_result = _scientific_result(accepted["f2"], candidate_f2, "representation", comparison.F2_PROVENANCE)
    _require(f2_result == {"scientific_payload_match": True, "first_mismatch_path": None, "first_mismatch": None}, "independent F.2 scientific equality failed")
    _require(all(_comparison_section(expected, "f2").get(key) == value for key, value in f2_result.items()), "reviewed comparison F.2 gate invalid")

    candidate_f3 = comparison.f3.build_pion_hgcer_method_a_acceptance_map_artifact(artifacts, candidate_f2, f1_input_file_hashes=hashes, f2_input_file_sha256=f2_sha, input_paths={})
    f3_bytes = comparison._writer_bytes(candidate_f3)
    f3_sha = hashlib.sha256(f3_bytes).hexdigest()
    _require(f3_sha == _hex(expected_candidate_shas.get("f3"), 64, "comparison candidate F.3 SHA"), "candidate F.3 writer SHA mismatch")
    f3_result = _scientific_result(accepted["f3"], candidate_f3, "acceptance_map", comparison.F3_PROVENANCE)
    _require(f3_result == {"scientific_payload_match": True, "first_mismatch_path": None, "first_mismatch": None}, "independent F.3 scientific equality failed")
    _require(all(_comparison_section(expected, "f3").get(key) == value for key, value in f3_result.items()), "reviewed comparison F.3 gate invalid")

    f3_map = comparison._body(candidate_f3, "acceptance_map")
    override = {kinematic: {"source_file_sha256": f3_sha, "map_fingerprint": f3_map["fingerprint"],
                            "algorithm_fingerprint": f3_map["algorithm_fingerprint"],
                            "artifact_fingerprint": candidate_f3["artifact_fingerprint"],
                            "farm_source_head": "0" * 40}}
    _require(override == expected.get("diagnostic_f3_authority_override"), "comparison candidate F.3 diagnostic override mismatch")
    candidate_f4 = comparison.f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, candidate_f3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={}, accepted_f3_runtime_authority_by_kinematic=override)
    f4_result = _scientific_result(accepted["f4"], candidate_f4, "correction", comparison.F4_PROVENANCE)
    f4_details = comparison.compare_f4(comparison._body(accepted["f4"], "correction"), comparison._body(candidate_f4, "correction"))
    f4_result.update(f4_details)
    expected_f4 = _comparison_section(expected, "f4")
    for key in ("scientific_payload_match", "first_mismatch_path", "first_mismatch", "parents", "global_maxima", "setting_maxima"):
        _require(key in expected_f4 and f4_result[key] == expected_f4[key], f"F.4 comparison reproduction mismatch: {key}")
    summary = {f"{stage}_scientific_payload_match": result["scientific_payload_match"] for stage, result in (("f2", f2_result), ("f3", f3_result), ("f4", f4_result))}
    summary["first_changed_stage"] = comparison.first_changed_stage(*(summary[f"{stage}_scientific_payload_match"] for stage in ("f2", "f3", "f4")))
    _require(summary == {"f2_scientific_payload_match": True, "f3_scientific_payload_match": True, "f4_scientific_payload_match": False, "first_changed_stage": "F4"}, "materialization scientific gate not F4-only")
    _require(summary == expected.get("summary"), "comparison summary reproduction mismatch")
    candidates = {"f2": candidate_f2, "f3": candidate_f3, "f4": candidate_f4}
    bytes_by_stage = {"f2": f2_bytes, "f3": f3_bytes, "f4": comparison._writer_bytes(candidate_f4)}
    return candidates, bytes_by_stage, override, f1_rows, accepted_rows, summary, f4_result


def _publish(args, preflight, built):
    artifacts, hashes, provisional, expected_f1, kinematic, names, targets, expected, comparison_sha, source_head, accepted, accepted_rows = preflight
    candidates, expected_bytes, override, f1_rows, _built_accepted_rows, summary, f4_result = built
    output_rows = {}
    # All preflight and scientific gates have completed before a temporary file is created.
    with tempfile.TemporaryDirectory(prefix="method-a-current-baseline-", dir=args.output_dir) as temp_name:
        temp = Path(temp_name)
        copied = temp / names["comparison"]
        copied.write_bytes(args.comparison.read_bytes())
        _require(copied.read_bytes() == args.comparison.read_bytes() and hashlib.sha256(copied.read_bytes()).hexdigest() == comparison_sha, "comparison copy is not byte-identical")
        for stage, body_name, _module, _excluded in BUILDERS:
            artifact = candidates[stage]
            _require(artifact.get("non_authoritative") is True, f"candidate {stage} is authoritative")
            path = temp / names[stage]
            getattr(_module, WRITER_NAMES[stage])(path, artifact)
            raw = path.read_bytes()
            _require(raw == expected_bytes[stage], f"candidate {stage} public-writer bytes differ")
            sha = hashlib.sha256(raw).hexdigest()
            body = comparison._body(artifact, body_name)
            output_rows[stage] = {"basename": names[stage], "raw_sha256": sha,
                                  "artifact_fingerprint": artifact["artifact_fingerprint"],
                                  **_stage_fingerprints(stage, body),
                                  "non_authoritative": True}
        f4_body = comparison._body(candidates["f4"], "correction")
        output_rows["f4"].update({"f3_source_file_sha256": f4_body["f3_source_file_sha256"],
                                  "f3_map_fingerprint": f4_body["f3_map_fingerprint"],
                                  "f3_algorithm_fingerprint": f4_body["f3_algorithm_fingerprint"],
                                  "f3_artifact_fingerprint": f4_body["f3_artifact_fingerprint"]})
        _require(output_rows["f2"]["raw_sha256"] == expected["candidate_serialized_sha256"]["f2"] and output_rows["f3"]["raw_sha256"] == expected["candidate_serialized_sha256"]["f3"], "candidate written SHA drift")
        manifest = {"schema_version": SCHEMA, "non_authoritative": True,
                    "accepted_authority_mutated": False, "production_objects_mutated": False,
                    "production_application_performed": False, "method_a_promoted": False,
                    "source_head": source_head, "kinematic_token": kinematic,
                    "comparison_input": {"source_basename": args.comparison.name, "raw_sha256": comparison_sha,
                                         "copied_output_basename": names["comparison"],
                                         "copied_output_raw_sha256": comparison_sha,
                                         "schema_version": expected["schema_version"]},
                    "current_f1_inputs": f1_rows, "accepted_inputs": accepted_rows,
                    "candidate_outputs": output_rows,
                    "scientific_gate": {**summary, "f4_first_mismatch_path": f4_result["first_mismatch_path"],
                                        "f4_global_maxima": f4_result["global_maxima"]},
                    "diagnostic_f3_authority_override": {"accepted_authority": False,
                                                         "purpose": "diagnostic_candidate_construction_only",
                                                         "record": override},
                    "complete": True, "errors": []}
        manifest_bytes = comparison._writer_bytes(manifest)
        (temp / names["manifest"]).write_bytes(manifest_bytes)
        _require((temp / names["manifest"]).read_bytes() == manifest_bytes, "manifest temporary bytes drift")
        # Recheck exact targets before publishing. The manifest is the final marker.
        for key, target in targets.items():
            _require(target.resolve() not in {path.resolve() for path in [*comparison._f1_inputs(args.f1).values(), args.accepted_f2, args.accepted_f3, args.accepted_f4, args.comparison]}, "output/input alias appeared")
            _require(args.overwrite or not (target.exists() or target.is_symlink()), f"output appeared during materialization: {target}")
        previous_manifest = temp / "previous-manifest.json"
        had_manifest = args.overwrite and targets["manifest"].exists()
        if had_manifest:
            os.replace(targets["manifest"], previous_manifest)
        published_candidate = False
        try:
            for key in ("comparison", "f2", "f3", "f4"):
                os.replace(temp / names[key], targets[key])
                published_candidate = True
            _require(targets["comparison"].read_bytes() == args.comparison.read_bytes(), "published comparison copy drift")
            for stage in ("f2", "f3", "f4"):
                _require(targets[stage].read_bytes() == expected_bytes[stage], f"published candidate {stage} bytes drift")
            os.replace(temp / names["manifest"], targets["manifest"])
        except Exception:
            if published_candidate:
                targets["manifest"].unlink(missing_ok=True)
            elif had_manifest:
                os.replace(previous_manifest, targets["manifest"])
            raise
    return manifest


def run(args):
    preflight = _preflight(args)
    built = _build(args, preflight)
    return _publish(args, preflight, built)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--f1", action="append", required=True, help="Canonical ALIAS=PATH; repeat five times")
    for stage in ("f2", "f3", "f4"):
        parser.add_argument(f"--accepted-{stage}", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--expected-comparison-sha256", required=True)
    parser.add_argument("--source-head", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, KeyError, TypeError, OverflowError, comparison.f2.MethodAAcceptanceRepresentationError, comparison.f3.MethodAAcceptanceMapError, comparison.f4.MethodAParentPreservingCorrectionError) as exc:
        parser.exit(1, f"materialization failed: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
