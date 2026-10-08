"""Fail-closed transient candidate F.4 lineage for the private F.6.3 branch.

This module deliberately contains no ROOT import and no tree traversal.  It
reuses the accepted F.4 shared calculator, joins its transient factors to the
already-existing yield cache by the frozen F.1 identity, and returns aggregate
provenance separately from the short-lived identity-to-factor map.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections.abc import Mapping, Sequence

import numpy as np

F6_3_PARALLEL_SCHEMA_VERSION = "f6_3_parallel_method_a_authority/v1"
_CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))
_TOLERANCE = 1.0e-12

# Detached F.6.3 current-baseline candidate lineage only. These records do not
# replace the historical accepted F.3/F.4 authorities or promote Method A.
F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD = "463d2657f696ecee33113edc3393ac51083a8944"
F6_3_CANDIDATE_BASENAMES_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "f2": "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json",
        "f3": "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json",
        "f4": "Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json",
    },
}
F6_3_CANDIDATE_F2_AUTHORITY = {
    "source_file_sha256": "2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e",
    "representation_fingerprint": "e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216",
    "artifact_fingerprint": "87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6",
}

# Closed, top-level provenance exclusions audited against the public builders.
# Fingerprint inputs duplicate scientific body fields as well as lineage IDs.
# Unknown fields stay compared. This is shared with the detached comparator.
SCIENTIFIC_EQUIVALENCE_SCHEMA = "f6_3_scientific_equivalence/v1"
F2_PROVENANCE = frozenset(("input_fingerprints", "fingerprint_inputs", "fingerprint"))
F3_PROVENANCE = F2_PROVENANCE | frozenset(("f2_representation_fingerprint", "f2_source_file_sha256"))
F4_PROVENANCE = F2_PROVENANCE | frozenset(("f3_source_file_sha256", "f3_map_fingerprint", "f3_artifact_fingerprint", "f3_runtime_authority"))
SCIENTIFIC_PROVENANCE_EXCLUSIONS = {"f2": F2_PROVENANCE, "f3": F3_PROVENANCE, "f4": F4_PROVENANCE}


def scientific_projection(body, excluded):
    """Remove only named top-level provenance; retain every other field."""
    return {key: value for key, value in body.items() if key not in excluded}


def first_mismatch(left, right, path="$"):
    """Exact recursive equality, including types, inventory and unknown fields."""
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


def _writer_bytes(artifact):
    """Match public artifact writers; transient serialization is never persisted."""
    return (json.dumps(artifact, sort_keys=True, indent=2, allow_nan=False) + "\n").replace("\n", os.linesep).encode("utf-8")
F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256 = {
    "Left-lowe": "eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07",
    "Left-highe": "544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e",
    "Center-lowe": "2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16",
    "Center-highe": "c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941",
    "Right-highe": "e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652",
}
# The zero farm head is a candidate-construction sentinel only, required to
# reproduce the materialized candidate F.4 exactly. It is not farm F.3 authority.
F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "source_file_sha256": "c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228",
        "map_fingerprint": "6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548",
        "algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
        "artifact_fingerprint": "e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302",
        "farm_source_head": "0000000000000000000000000000000000000000",
    },
}
F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "source_file_sha256": "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7",
        "correction_fingerprint": "71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98",
        "artifact_fingerprint": "0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4",
        "farm_source_head": F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
        "f3_source_file_sha256": F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["source_file_sha256"],
        "f3_map_fingerprint": F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["map_fingerprint"],
        "f3_algorithm_fingerprint": F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["algorithm_fingerprint"],
        "f3_artifact_fingerprint": F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["artifact_fingerprint"],
        "f1_source_file_sha256": F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256,
    },
}


class MethodAParallelFullProcedureError(RuntimeError):
    """The non-production Method-A branch cannot safely be constructed."""


def _runtime_dependencies():
    """Import accepted scientific calculators only when F.6.3 is requested."""
    import pion_hgcer_method_a_parent_preserving_correction as f4
    import pion_hgcer_method_a_tphi_propagation as f5
    from pion_hgcer_method_a_acceptance_contract import pion_hgcer_method_a_acceptance_contract_filename
    from pion_hgcer_method_a_acceptance_map import pion_hgcer_method_a_acceptance_map_filename
    from pion_hgcer_method_a_parent_preserving_correction import pion_hgcer_method_a_parent_preserving_correction_filename
    return f4, f5, pion_hgcer_method_a_acceptance_contract_filename, pion_hgcer_method_a_acceptance_map_filename, pion_hgcer_method_a_parent_preserving_correction_filename


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_json(value: object) -> str:
    return _sha256_bytes(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8"))


def _f4_scaled_close(left: float, right: float) -> bool:
    """Use F.4's frozen signed-baseline comparison semantics exactly."""
    return abs(left - right) <= _TOLERANCE * max(1.0, abs(left), abs(right))


def _load_json(path: str) -> tuple[dict[str, object], str]:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def bad_constant(value):
        raise ValueError("nonfinite JSON constant:" + value)

    try:
        with open(path, "rb") as handle:
            raw = handle.read()
        payload = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs, parse_constant=bad_constant)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise MethodAParallelFullProcedureError("f6_3_accepted_artifact_load_failed:{}".format(os.path.basename(path))) from exc
    if not isinstance(payload, dict):
        raise MethodAParallelFullProcedureError("f6_3_accepted_artifact_payload_invalid")
    return payload, _sha256_bytes(raw)


def _setting_id(phi_setting: object, epsilon_token: object) -> str:
    setting = (str(phi_setting), str(epsilon_token).strip().lower())
    if setting not in _CANONICAL_SETTINGS:
        raise MethodAParallelFullProcedureError("f6_3_setting_unsupported")
    return "{}-{}".format(*setting)


def accepted_f6_3_artifact_paths(outpath: object, kinematic_token: object) -> dict[str, object]:
    """Return current F.1 and exact candidate F.2/F.3/F.4 paths; never fallback."""
    root = os.fspath(outpath)
    token = str(kinematic_token)
    candidate = F6_3_CANDIDATE_BASENAMES_BY_KINEMATIC.get(token)
    if candidate is None:
        raise MethodAParallelFullProcedureError("f6_3_kinematic_unsupported")
    _f4, _f5, f1_filename, _f3_filename, _f4_filename = _runtime_dependencies()
    f1 = {
        "{}-{}".format(phi, epsilon): os.path.join(
            root,
            f1_filename(phi, "kaon", token, epsilon),
        )
        for phi, epsilon in _CANONICAL_SETTINGS
    }
    return {
        "f1": f1,
        "f2": os.path.join(root, candidate["f2"]),
        "f3": os.path.join(root, candidate["f3"]),
        "f4": os.path.join(root, candidate["f4"]),
    }


def load_accepted_f6_3_authority(paths: Mapping[str, object], *, include_f2: bool = False):
    """Load exact paths; F.6.3 opts into F.2, older read-only diagnostics keep their API."""
    raw_f1 = paths.get("f1")
    if not isinstance(raw_f1, Mapping) or set(raw_f1) != {"{}-{}".format(*item) for item in _CANONICAL_SETTINGS}:
        raise MethodAParallelFullProcedureError("f6_3_f1_path_inventory_invalid")
    f1, f1_hashes = [], {}
    for phi, epsilon in _CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon)
        payload, digest = _load_json(os.fspath(raw_f1[setting_id]))
        f1.append(payload); f1_hashes[setting_id] = digest
    f3, f3_hash = _load_json(os.fspath(paths.get("f3", "")))
    f4, f4_hash = _load_json(os.fspath(paths.get("f4", "")))
    if include_f2:
        f2, f2_hash = _load_json(os.fspath(paths.get("f2", "")))
        return f1, f2, f3, f4, {**f1_hashes, "f2": f2_hash, "f3": f3_hash, "f4": f4_hash}
    return f1, f3, f4, {**f1_hashes, "f3": f3_hash, "f4": f4_hash}


def _validate_candidate_f2_f3(f2_artifact, f3_artifact, f2_sha, f3_sha, f4):
    """Authenticate reviewed wrappers independently of regenerated F.1 lineage."""
    f3 = f4._f3
    definitions = (
        ("f2", f2_artifact, f2_sha, F6_3_CANDIDATE_F2_AUTHORITY,
         "representation", "representation_fingerprint", f3._F2_ARTIFACT_SCHEMA,
         f3._F2_SCHEMA, f3._F2_FINGERPRINT_SCHEMA,
         {"probe_only": True, "probability_map_constructed": False, "basis_frozen": False}),
        ("f3", f3_artifact, f3_sha, F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"],
         "acceptance_map", "map_fingerprint", f4._F3_ARTIFACT_SCHEMA,
         f4._F3_SCHEMA, f4._F3_FINGERPRINT_SCHEMA,
         {"map_constructed": True, "map_applied": False, "absolute_probability_constructed": False,
          "parent_normalization_constructed": False, "correction_constructed": False, "basis_frozen": True}),
    )
    common = {"non_authoritative": True, "manual_review_required": True,
              "production_objects_mutated": False, "production_application_performed": False,
              "method_b_numerical_dependency": False, "event_probability_persisted": False,
              "event_application_performed": False}
    for stage, artifact, raw_sha, pin, body_name, fp_name, schema, body_schema, fp_schema, flags in definitions:
        def require(value, reason):
            if not value:
                raise MethodAParallelFullProcedureError(f"{stage}_runtime_authority_{reason}_mismatch")
        require(isinstance(artifact, Mapping), "wrapper")
        body = artifact.get(body_name)
        require(isinstance(body, Mapping), "body")
        require(raw_sha == pin["source_file_sha256"], "source_file_sha256")
        require(artifact.get("schema_version") == schema, "schema")
        require(body.get("schema_version") == body_schema and body.get("fingerprint_schema_version") == fp_schema, "body_schema")
        for name, expected in (common | flags).items():
            require(artifact.get(name) is expected and body.get(name) is expected, name)
        require(body.get("available") is True and body.get("status") == "available" and
                body.get("diagnostic_stage") == "complete", "available")
        require(body.get("fingerprint") == pin[fp_name], fp_name)
        require(artifact.get("artifact_fingerprint") == pin["artifact_fingerprint"], "artifact_fingerprint")
        inputs = body.get("fingerprint_inputs")
        require(isinstance(inputs, Mapping) and _sha256_json(inputs) == body["fingerprint"], "fingerprint_inputs")
        bound = ("candidate_definitions", "algorithm_config", "algorithm_fingerprint", "response_support", "groups", "candidate_summaries", "recommendation") if stage == "f2" else ("algorithm_config", "algorithm_fingerprint", "f2_representation_fingerprint", "f2_algorithm_fingerprint", "f2_source_file_sha256", "models")
        require(all(inputs.get(name) == body.get(name) for name in bound), "fingerprint_content")
        provenance = artifact.get("provenance")
        require(isinstance(provenance, Mapping), "provenance")
        require(_sha256_json({"schema_version": schema, fp_name: body["fingerprint"],
                             "input_paths": provenance.get("input_paths")}) == artifact["artifact_fingerprint"], "artifact_fingerprint")
        if stage == "f3":
            require(body.get("algorithm_fingerprint") == pin["algorithm_fingerprint"], "algorithm_fingerprint")
            require(body.get("f2_source_file_sha256") == f2_sha and
                    body.get("f2_representation_fingerprint") == f2_artifact["representation"]["fingerprint"] and
                    body.get("f2_algorithm_fingerprint") == f2_artifact["representation"]["algorithm_fingerprint"], "f2_identity")


def _compare_stage(stage, reviewed, current):
    exclusions = SCIENTIFIC_PROVENANCE_EXCLUSIONS[stage]
    mismatch = first_mismatch(scientific_projection(reviewed, exclusions),
                              scientific_projection(current, exclusions))
    if mismatch is not None:
        raise MethodAParallelFullProcedureError(f"f6_3_scientific_equivalence_mismatch:{stage}:{mismatch['path']}")
    return {"passed": True, "first_mismatch_path": None}


def _lineage_reconstruction(f1_artifacts, f2_artifact, f3_artifact, persisted, parsed,
                            hashes, f2_sha, f3_sha, f4):
    """Sequential exact science gates; factors remain private until all pass."""
    f3 = f4._f3
    import pion_hgcer_method_a_acceptance_representation as f2
    expected = F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["f1_source_file_sha256"]
    exact = hashes == expected
    decision = {"schema_version": SCIENTIFIC_EQUIVALENCE_SCHEMA,
                "mode": "exact-lineage" if exact else "equivalent-new-lineage",
                "provenance_exclusions": {stage: sorted(fields) for stage, fields in SCIENTIFIC_PROVENANCE_EXCLUSIONS.items()}}
    if exact:
        # Keep the ordinary persisted-lineage validators and original F.4 path.
        f3._validate_f2_artifact(f2_artifact, parsed,
                               [{"setting_id": item["setting_id"], "sha256": hashes[item["setting_id"]]} for item in parsed])
        current2, current3 = f2_artifact, f3_artifact
        decision["f2"] = _compare_stage("f2", f2_artifact["representation"], current2["representation"])
        decision["f3"] = _compare_stage("f3", f3_artifact["acceptance_map"], current3["acceptance_map"])
        correction, review = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            f1_artifacts, f3_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha,
            accepted_f3_runtime_authority_by_kinematic=F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC)
        if correction != persisted or correction.get("fingerprint") != persisted.get("fingerprint"):
            raise MethodAParallelFullProcedureError("f6_3_f4_shared_reproduction_mismatch")
        current4 = f4._artifact_wrapper(correction, input_paths={}, generated_at_utc=None, git_head=None, git_status_short=None)
    else:
        current2 = f2.build_pion_hgcer_method_a_acceptance_representation_artifact(
            f1_artifacts, input_file_hashes=hashes, input_paths={})
        f2_sha = _sha256_bytes(_writer_bytes(current2))
        decision["f2"] = _compare_stage("f2", f2_artifact["representation"], current2["representation"])
        current3 = f3.build_pion_hgcer_method_a_acceptance_map_artifact(
            f1_artifacts, current2, f1_input_file_hashes=hashes, f2_input_file_sha256=f2_sha, input_paths={})
        f3_sha = _sha256_bytes(_writer_bytes(current3))
        decision["f3"] = _compare_stage("f3", f3_artifact["acceptance_map"], current3["acceptance_map"])
        body3 = current3["acceptance_map"]
        # Explicit construction sentinel, not accepted farm F.3 authority.
        transient_authority = {"Q4p4W2p74": {"source_file_sha256": f3_sha,
            "map_fingerprint": body3["fingerprint"], "algorithm_fingerprint": body3["algorithm_fingerprint"],
            "artifact_fingerprint": current3["artifact_fingerprint"], "farm_source_head": "0" * 40}}
        current4, review = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data(
            f1_artifacts, current3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha,
            accepted_f3_runtime_authority_by_kinematic=transient_authority, input_paths={})
    correction = current4["correction"]
    decision["f4"] = _compare_stage("f4", persisted, correction)
    expected_parents = {(sid, index) for sid in hashes for index in range(3)}
    parents = correction.get("parents", [])
    if len(parents) != 15 or {(p.get("setting_id"), p.get("canonical_t_index")) for p in parents} != expected_parents:
        raise MethodAParallelFullProcedureError("f6_3_current_parent_inventory_invalid")
    if len(review) != 15 or {(p.get("setting_id"), p.get("canonical_t_index")) for p in review} != expected_parents:
        raise MethodAParallelFullProcedureError("f6_3_current_review_inventory_invalid")
    decision["all_stages_passed"] = True
    current = {"reconstruction_role": "transient_current_lineage_no_farm_authority",
               "f1_source_file_sha256": dict(sorted(hashes.items())),
               "f1_stable_content_fingerprints": {p["setting_id"]: p["stable_content_fingerprint"] for p in parsed}}
    for stage, artifact, body_name, sha in (("f2", current2, "representation", f2_sha),
            ("f3", current3, "acceptance_map", f3_sha),
            ("f4", current4, "correction", _sha256_bytes(_writer_bytes(current4)))):
        current[stage] = {"source_file_sha256": sha, "scientific_fingerprint": artifact[body_name]["fingerprint"],
                          "artifact_fingerprint": artifact["artifact_fingerprint"]}
    return correction, review, current, decision


def reconstruct_transient_factor_map(
    f1_artifacts: Sequence[Mapping[str, object]],
    f3_artifact: Mapping[str, object],
    f4_artifact: Mapping[str, object],
    *,
    f2_artifact: Mapping[str, object],
    f2_input_file_sha256: str,
    f1_input_file_hashes: Mapping[str, object],
    f3_input_file_sha256: str,
    f4_input_file_sha256: str,
    setting_id: str,
) -> tuple[
    dict[tuple[str, int], float],
    dict[str, object],
    dict[tuple[str, int], dict[str, object]],
]:
    """Reproduce F.4 and expose factors only for immediate branch filling."""
    if setting_id not in {"{}-{}".format(*item) for item in _CANONICAL_SETTINGS}:
        raise MethodAParallelFullProcedureError("f6_3_setting_unsupported")
    _f4, _f5, _f1_filename, _f3_filename, _f4_filename = _runtime_dependencies()
    try:
        _artifact, persisted, authority = _f5._validate_f4_artifact(
            f4_artifact, str(f4_input_file_sha256),
            F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC,
        )
        parsed = _f4._f3._validate_f1_artifacts(f1_artifacts)
        hashes = {str(key): str(value) for key, value in f1_input_file_hashes.items()}
        if set(hashes) != {str(item["setting_id"]) for item in parsed}:
            raise MethodAParallelFullProcedureError("f6_3_f1_hash_inventory_invalid")
        if set(hashes) != {"{}-{}".format(*item) for item in _CANONICAL_SETTINGS} or any(p["setting"]["kinematic_token"] != "Q4p4W2p74" for p in parsed):
            raise MethodAParallelFullProcedureError("f6_3_f1_setting_inventory_invalid")
        _validate_candidate_f2_f3(f2_artifact, f3_artifact, str(f2_input_file_sha256), str(f3_input_file_sha256), _f4)
        recomputed, review_data, current_lineage, equivalence = _lineage_reconstruction(
            f1_artifacts, f2_artifact, f3_artifact, persisted, parsed, hashes,
            str(f2_input_file_sha256), str(f3_input_file_sha256), _f4)
        rows_by_parent = _f4._raw_application_rows(f1_artifacts, parsed, _TOLERANCE)
        # F.4's sanitized calculation rows deliberately omit these baseline
        # coordinates. Join them from the same fully audited raw F.1 population
        # for F.6.3 parity only; do not change the shared F.4 calculator.
        baseline_coordinates = {
            (str(row["source_label"]), int(row["entry_index"])): {
                name: row[name] for name in ("analysis_MM", "analysis_t")
            }
            for artifact in f1_artifacts
            if "{}-{}".format(
                artifact["setting"]["phi_setting"],
                artifact["setting"]["epsilon_filename_token"],
            ) == setting_id
            for row in artifact["contract"]["application_records"]
        }
    except MethodAParallelFullProcedureError:
        raise
    except Exception as exc:
        raise MethodAParallelFullProcedureError("f6_3_f4_shared_reproduction_failed:{}".format(exc)) from exc
    factors: dict[tuple[str, int], float] = {}
    accepted_rows: dict[tuple[str, int], dict[str, object]] = {}
    for review in review_data:
        if str(review.get("setting_id")) != setting_id:
            continue
        key = (setting_id, int(review.get("canonical_t_index")))
        rows = rows_by_parent.get(key)
        values = np.asarray(review.get("correction_factors"), dtype=float)
        if rows is None or values.ndim != 1 or values.size != len(rows):
            raise MethodAParallelFullProcedureError("f6_3_factor_row_alignment_invalid")
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise MethodAParallelFullProcedureError("f6_3_factor_nonfinite_or_nonpositive")
        for row, value in zip(rows, values):
            identity = (str(row.get("source_label")), int(row.get("entry_index")))
            if identity in factors:
                raise MethodAParallelFullProcedureError("f6_3_factor_identity_duplicate")
            factors[identity] = float(value)
            accepted_rows[identity] = {
                name: row[name]
                for name in (
                    "t_index", "phi_index",
                    "signed_source_coefficient", "baseline_pion_weight_w0",
                    "signed_baseline_event_contribution",
                )
            }
            accepted_rows[identity].update(baseline_coordinates[identity])
    if not factors:
        raise MethodAParallelFullProcedureError("f6_3_selected_factor_population_missing")
    identities = sorted((source, entry) for source, entry in factors)
    provenance = {
        "schema_version": F6_3_PARALLEL_SCHEMA_VERSION,
        "branch_role": "parallel_nonproduction_method_a_full_analysis",
        "current_baseline_candidate_lineage": True,
        "candidate_validation_source_head": F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
        "candidate_f3_reconstruction_role": "candidate_construction_sentinel_only",
        "candidate_f3_source_file_sha256": str(f3_input_file_sha256),
        "candidate_f2_source_file_sha256": str(f2_input_file_sha256),
        "candidate_f4_source_file_sha256": str(f4_input_file_sha256),
        "production_promotion_performed": False,
        "selected_setting_id": setting_id,
        "accepted_f1_source_file_sha256": dict(sorted(hashes.items())),
        "accepted_f3_source_file_sha256": str(f3_input_file_sha256),
        "accepted_f3_map_fingerprint": str(persisted.get("f3_map_fingerprint")),
        "accepted_f3_algorithm_fingerprint": str(persisted.get("f3_algorithm_fingerprint")),
        "accepted_f3_artifact_fingerprint": str(persisted.get("f3_artifact_fingerprint")),
        "accepted_f4_source_file_sha256": str(f4_input_file_sha256),
        "accepted_f4_correction_fingerprint": str(persisted.get("fingerprint")),
        "accepted_f4_artifact_fingerprint": str(f4_artifact.get("artifact_fingerprint")),
        "accepted_f4_runtime_authority": authority,
        "transient_factor_population_count": len(factors),
        "transient_factor_identity_fingerprint": _sha256_json(identities),
        "event_correction_persisted": False,
        "reviewed_candidate": {
            "validation_materialization_source_head": F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
            "f1_source_file_sha256": dict(F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]["f1_source_file_sha256"]),
            "f2": dict(F6_3_CANDIDATE_F2_AUTHORITY),
            "f3": dict(F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]),
            "f4": dict(F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC["Q4p4W2p74"]),
        },
        "current_runtime_lineage": current_lineage,
        "scientific_equivalence": equivalence,
    }
    return factors, provenance, accepted_rows


def validate_live_cache_parity(
    transient_factors: Mapping[tuple[str, int], float],
    live_rows: Sequence[Mapping[str, object]],
    accepted_rows: Mapping[tuple[str, int], Mapping[str, object]],
) -> dict[str, object]:
    """Require exact F.1 identity/baseline parity before any multiplier is used."""
    seen: dict[tuple[str, int], Mapping[str, object]] = {}
    required = (
        "source_label", "entry_index", "t_index", "phi_index", "analysis_MM",
        "analysis_t", "signed_source_coefficient", "baseline_pion_weight_w0",
        "signed_baseline_event_contribution",
    )
    for row in live_rows:
        if any(name not in row for name in required):
            raise MethodAParallelFullProcedureError("f6_3_live_cache_field_missing")
        identity = (str(row["source_label"]), int(row["entry_index"]))
        if identity in seen:
            raise MethodAParallelFullProcedureError("f6_3_live_cache_identity_duplicate")
        values = (row["analysis_MM"], row["analysis_t"], row["signed_source_coefficient"], row["baseline_pion_weight_w0"], row["signed_baseline_event_contribution"])
        if not all(math.isfinite(float(value)) for value in values):
            raise MethodAParallelFullProcedureError("f6_3_live_cache_nonfinite")
        if float(row["baseline_pion_weight_w0"]) < 0.0:
            raise MethodAParallelFullProcedureError("f6_3_live_cache_w0_negative")
        if not _f4_scaled_close(
            float(row["signed_baseline_event_contribution"]),
            float(row["signed_source_coefficient"])
            * float(row["baseline_pion_weight_w0"]),
        ):
            raise MethodAParallelFullProcedureError("f6_3_live_cache_baseline_identity_mismatch")
        seen[identity] = row
    expected = set(transient_factors)
    observed = set(seen)
    if expected != observed:
        raise MethodAParallelFullProcedureError("f6_3_live_cache_identity_inventory_mismatch")
    if set(accepted_rows) != expected:
        raise MethodAParallelFullProcedureError("f6_3_accepted_identity_inventory_mismatch")
    for identity, expected_row in accepted_rows.items():
        observed_row = seen[identity]
        for name in ("t_index", "phi_index"):
            if int(observed_row[name]) != int(expected_row[name]):
                raise MethodAParallelFullProcedureError("f6_3_live_cache_f1_{}_mismatch".format(name))
        for name in (
            "analysis_MM", "analysis_t", "signed_source_coefficient",
            "baseline_pion_weight_w0", "signed_baseline_event_contribution",
        ):
            if not _f4_scaled_close(
                float(observed_row[name]), float(expected_row[name]),
            ):
                raise MethodAParallelFullProcedureError("f6_3_live_cache_f1_{}_mismatch".format(name))
    for identity, factor in transient_factors.items():
        if not math.isfinite(float(factor)) or float(factor) <= 0.0:
            raise MethodAParallelFullProcedureError("f6_3_factor_nonfinite_or_nonpositive")
    return {
        "live_cache_identity_count": len(seen),
        "live_cache_identity_fingerprint": _sha256_json(sorted((source, entry) for source, entry in seen)),
        "live_cache_parity_passed": True,
    }


def unavailable_parallel_source(reason: object, *, setting_id: object = None) -> dict[str, object]:
    """Return the safe consumer boundary when the optional branch cannot run."""
    return {
        "schema_version": "f6_3_parallel_method_a_source/v1",
        "available": False,
        "reason": str(reason),
        "selected_setting_id": None if setting_id is None else str(setting_id),
        "branch_role": "parallel_nonproduction_method_a_full_analysis",
        "baseline_production_mutated": False,
        "production_promotion_performed": False,
        "method_b_numerical_dependency": False,
        "empirical_residual_used": False,
        "event_correction_persisted": False,
        "canonical_child_renormalization_performed": False,
        "baseline_public_output_unchanged": True,
    }
