"""Detached, read-only F.2/F.3/F.4 current-baseline authority comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "cuts"))
import pion_hgcer_method_a_acceptance_representation as f2  # noqa: E402
import pion_hgcer_method_a_acceptance_map as f3  # noqa: E402
import pion_hgcer_method_a_parent_preserving_correction as f4  # noqa: E402

ALIASES = ("Left-lowe", "Left-highe", "Center-lowe", "Center-highe", "Right-highe")
F2_PROVENANCE = frozenset(("input_fingerprints", "fingerprint_inputs", "fingerprint"))
F3_PROVENANCE = F2_PROVENANCE | frozenset(("f2_representation_fingerprint", "f2_source_file_sha256"))
F4_PROVENANCE = F2_PROVENANCE | frozenset(("f3_source_file_sha256", "f3_map_fingerprint", "f3_artifact_fingerprint", "f3_runtime_authority"))
REQUIRED = {
    "representation": ("candidate_definitions", "algorithm_config", "algorithm_fingerprint", "response_support", "groups", "candidate_summaries", "recommendation"),
    "acceptance_map": ("accepted_basis", "ordered_features", "algorithm_config", "algorithm_fingerprint", "models"),
    "correction": ("accepted_basis", "algorithm_config", "parents"),
}
PARENT_METRICS = (
    "application_event_count", "zero_baseline_weight_count", "in_support_count", "ood_count", "ood_fraction",
    "baseline_parent_sum", "absolute_baseline_parent_sum", "baseline_cancellation_ratio", "raw_shape_parent_sum",
    "parent_normalization", "adjusted_parent_sum", "closure_residual", "ood_final_correction",
    "raw_shape_factor_summary", "correction_factor_summary", "source_diagnostics", "canonical_phi_diagnostics",
)


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _bad_constant(value):
    raise ValueError(f"non-finite JSON constant: {value}")


def _read(path):
    raw = Path(path).read_bytes()
    obj = json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_bad_constant)
    if not isinstance(obj, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return obj, hashlib.sha256(raw).hexdigest()


def _writer_bytes(artifact):
    # Same encoding, key order, indentation, and trailing newline as the public writers.
    return (json.dumps(artifact, sort_keys=True, indent=2, allow_nan=False) + "\n").replace("\n", os.linesep).encode("utf-8")


def _body(artifact, name):
    if not isinstance(artifact, dict) or not isinstance(artifact.get(name), dict):
        raise ValueError(f"missing {name} payload")
    body = artifact[name]
    for key in REQUIRED[name]:
        if key not in body:
            raise ValueError(f"missing {name}.{key}")
    if name == "representation" and (not isinstance(body["groups"], list) or len(body["groups"]) != 15):
        raise ValueError("F.2 group inventory must contain 15 groups")
    if name == "acceptance_map" and (not isinstance(body["models"], list) or len(body["models"]) != 15):
        raise ValueError("F.3 model inventory must contain 15 models")
    return body


def scientific_projection(body, excluded):
    """Drop only the declared provenance fields; retain future scientific fields."""
    return {key: value for key, value in body.items() if key not in excluded}


def first_mismatch(left, right, path="$"):
    if type(left) is not type(right):
        return {"path": path, "accepted": left, "candidate": right}
    if isinstance(left, dict):
        for key in sorted(left.keys() | right.keys()):
            child = f"{path}.{key}"
            if key not in left or key not in right:
                return {"path": child, "accepted": left.get(key), "candidate": right.get(key), "missing_side": "accepted" if key not in left else "candidate"}
            found = first_mismatch(left[key], right[key], child)
            if found:
                return found
        return None
    if isinstance(left, list):
        if len(left) != len(right):
            return {"path": f"{path}.length", "accepted": len(left), "candidate": len(right)}
        for index, (a, b) in enumerate(zip(left, right)):
            found = first_mismatch(a, b, f"{path}[{index}]")
            if found:
                return found
        return None
    if left != right:
        return {"path": path, "accepted": left, "candidate": right}
    return None


def _delta(accepted, candidate):
    if isinstance(accepted, (int, float)) and not isinstance(accepted, bool) and isinstance(candidate, (int, float)) and not isinstance(candidate, bool):
        absolute = abs(candidate - accepted)
        relative = absolute / abs(accepted) if accepted != 0 else None
        return {"accepted": accepted, "candidate": candidate, "absolute_difference": absolute, "relative_difference": relative}
    if isinstance(accepted, dict) and isinstance(candidate, dict):
        return {key: _delta(accepted[key], candidate[key]) for key in sorted(accepted.keys() | candidate.keys()) if key in accepted and key in candidate} | {
            key: {"accepted": accepted.get(key), "candidate": candidate.get(key), "missing_side": "accepted" if key not in accepted else "candidate"}
            for key in sorted(accepted.keys() ^ candidate.keys())
        }
    if isinstance(accepted, list) and isinstance(candidate, list):
        return [
            _delta(accepted[index], candidate[index])
            if index < len(accepted) and index < len(candidate)
            else {
                "accepted": accepted[index] if index < len(accepted) else None,
                "candidate": candidate[index] if index < len(candidate) else None,
                "missing_side": "accepted" if index >= len(accepted) else "candidate",
            }
            for index in range(max(len(accepted), len(candidate)))
        ]
    return {"accepted": accepted, "candidate": candidate, "match": accepted == candidate}


def _numeric_deltas(tree):
    if isinstance(tree, dict):
        if "absolute_difference" in tree:
            yield tree
        else:
            for value in tree.values():
                yield from _numeric_deltas(value)
    elif isinstance(tree, list):
        for value in tree:
            yield from _numeric_deltas(value)


def _parent_index(body):
    parents = body["parents"]
    if not isinstance(parents, list) or len(parents) != 15:
        raise ValueError("F.4 parent inventory must contain 15 parents")
    result = {}
    for row in parents:
        if not isinstance(row, dict) or not isinstance(row.get("setting_id"), str) or type(row.get("canonical_t_index")) is not int:
            raise ValueError("malformed F.4 parent identity")
        key = (row["setting_id"], row["canonical_t_index"])
        if key in result:
            raise ValueError(f"duplicate F.4 parent: {key}")
        for metric in PARENT_METRICS:
            if metric not in row:
                raise ValueError(f"missing F.4 parent metric: {metric}")
        result[key] = row
    return result


def compare_f4(accepted, candidate):
    a, c = _parent_index(accepted), _parent_index(candidate)
    if a.keys() != c.keys():
        raise ValueError("accepted/candidate F.4 parent inventories differ")
    comparisons = []
    maxima = {name: {"absolute_difference": 0, "relative_difference": 0} for name in ("parent_normalization", "correction_factor_summary", "baseline_parent_sum")}
    by_setting = {}
    for key in sorted(a):
        metrics = {metric: _delta(a[key][metric], c[key][metric]) for metric in PARENT_METRICS}
        comparisons.append({"setting_id": key[0], "canonical_t_index": key[1], "metrics": metrics})
        setting_max = by_setting.setdefault(key[0], {name: {"absolute_difference": 0, "relative_difference": 0} for name in maxima})
        for name in maxima:
            for delta in _numeric_deltas(metrics[name]):
                for label in ("absolute_difference", "relative_difference"):
                    value = delta[label]
                    if value is not None:
                        maxima[name][label] = max(maxima[name][label], value)
                        setting_max[name][label] = max(setting_max[name][label], value)
    return {"parents": comparisons, "global_maxima": maxima, "setting_maxima": by_setting}


def first_changed_stage(f2_match, f3_match, f4_match):
    return next((stage for stage, matched in (("F2", f2_match), ("F3", f3_match), ("F4", f4_match)) if not matched), "none")


def _f1_inputs(entries):
    paths = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError("--f1 requires ALIAS=PATH")
        alias, path = entry.split("=", 1)
        if alias not in ALIASES or alias in paths or not path:
            raise ValueError(f"unexpected or duplicate F.1 alias: {alias}")
        paths[alias] = Path(path)
    if set(paths) != set(ALIASES):
        raise ValueError(f"canonical-five F.1 inventory required: {ALIASES}")
    return paths


def run(args):
    paths = _f1_inputs(args.f1)
    all_inputs = [*paths.values(), args.accepted_f2, args.accepted_f3, args.accepted_f4]
    resolved = [path.resolve() for path in all_inputs]
    if len(set(resolved)) != len(resolved) or args.output.resolve() in resolved:
        raise ValueError("input paths must be distinct and output must not be an input")
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"output exists: {args.output}")
    f1_artifacts, hashes, f1_provenance = [], {}, []
    for alias in ALIASES:
        artifact, sha = _read(paths[alias])
        setting = artifact.get("setting")
        if not isinstance(setting, dict) or f"{setting.get('phi_setting')}-{setting.get('epsilon_filename_token')}" != alias:
            raise ValueError(f"F.1 alias/artifact setting mismatch: {alias}")
        f1_artifacts.append(artifact)
        f1_provenance.append({"alias": alias, "source_file_sha256": sha})
    accepted = {}
    for name, path in (("f2", args.accepted_f2), ("f3", args.accepted_f3), ("f4", args.accepted_f4)):
        accepted[name], sha = _read(path)
        accepted[name + "_sha256"] = sha
    # The public builders enforce the complete current F.1 contract and downstream continuity.
    # Use explicit byte hashes indexed by the validated F.1 setting IDs.
    for artifact, item in zip(f1_artifacts, f1_provenance):
        sid = f"{artifact['setting']['phi_setting']}-{artifact['setting']['epsilon_filename_token']}"
        hashes[sid] = item["source_file_sha256"]
        item["setting_id"] = sid
    candidate_f2 = f2.build_pion_hgcer_method_a_acceptance_representation_artifact(f1_artifacts, input_file_hashes=hashes, input_paths={})
    f2_sha = hashlib.sha256(_writer_bytes(candidate_f2)).hexdigest()
    candidate_f3 = f3.build_pion_hgcer_method_a_acceptance_map_artifact(f1_artifacts, candidate_f2, f1_input_file_hashes=hashes, f2_input_file_sha256=f2_sha, input_paths={})
    f3_sha = hashlib.sha256(_writer_bytes(candidate_f3)).hexdigest()
    f3_map = _body(candidate_f3, "acceptance_map")
    kinematics = {item["setting"]["kinematic_token"] for item in f1_artifacts}
    if len(kinematics) != 1:
        raise ValueError("F.1 kinematics differ")
    kinematic = kinematics.pop()
    # F.4 requires a 40-hex farm head in its authority record. All zeros is
    # an explicit diagnostic sentinel: the candidate has no farm authority.
    override = {kinematic: {"source_file_sha256": f3_sha, "map_fingerprint": f3_map["fingerprint"], "algorithm_fingerprint": f3_map["algorithm_fingerprint"], "artifact_fingerprint": candidate_f3["artifact_fingerprint"], "farm_source_head": "0" * 40}}
    candidate_f4 = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(f1_artifacts, candidate_f3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={}, accepted_f3_runtime_authority_by_kinematic=override)
    for item, fp in zip(f1_provenance, _body(candidate_f2, "representation")["input_fingerprints"]):
        if item["setting_id"] != fp["setting_id"] or item["source_file_sha256"] != fp["source_file_sha256"]:
            raise ValueError("candidate F.2 F.1 fingerprint inventory mismatch")
        item.update({"setting": fp["setting"], "stable_f1_content_fingerprint": fp["stable_f1_content_fingerprint"], "f1_contract_fingerprint": fp["fingerprint"], "training_record_count": len(f1_artifacts[ALIASES.index(item["alias"])]["contract"]["method_a_training_records"]), "application_record_count": len(f1_artifacts[ALIASES.index(item["alias"])]["contract"]["application_records"]), "training_population_fingerprint": fp["method_a_training_population_fingerprint"], "application_population_fingerprint": fp["application_population_fingerprint"]})
    comparisons = {}
    for name, exclusion, stage in (("representation", F2_PROVENANCE, "f2"), ("acceptance_map", F3_PROVENANCE, "f3"), ("correction", F4_PROVENANCE, "f4")):
        a = scientific_projection(_body(accepted[stage], name), exclusion)
        c = scientific_projection(_body({"f2": candidate_f2, "f3": candidate_f3, "f4": candidate_f4}[stage], name), exclusion)
        mismatch = first_mismatch(a, c)
        comparisons[stage] = {"scientific_payload_match": mismatch is None, "first_mismatch_path": None if mismatch is None else mismatch["path"], "first_mismatch": mismatch}
    f4_details = compare_f4(_body(accepted["f4"], "correction"), _body(candidate_f4, "correction"))
    summary = {f"{stage}_scientific_payload_match": comparisons[stage]["scientific_payload_match"] for stage in ("f2", "f3", "f4")}
    summary["first_changed_stage"] = first_changed_stage(*(summary[f"{stage}_scientific_payload_match"] for stage in ("f2", "f3", "f4")))
    result = {"schema_version": "method_a_current_baseline_authority_comparison/v1", "non_authoritative": True, "f1_inputs": f1_provenance, "accepted_file_sha256": {stage: accepted[stage + "_sha256"] for stage in ("f2", "f3", "f4")}, "candidate_serialized_sha256": {"f2": f2_sha, "f3": f3_sha}, "diagnostic_f3_authority_override": override, "f2": comparisons["f2"], "f3": comparisons["f3"], "f4": {**comparisons["f4"], **f4_details}, "summary": summary}
    args.output.write_bytes(_writer_bytes(result))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--f1", action="append", required=True, help="Canonical ALIAS=PATH; repeat five times")
    for name in ("f2", "f3", "f4"):
        parser.add_argument(f"--accepted-{name}", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, KeyError, TypeError, OverflowError, f2.MethodAAcceptanceRepresentationError, f3.MethodAAcceptanceMapError, f4.MethodAParentPreservingCorrectionError) as exc:
        parser.exit(1, f"comparison failed: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
