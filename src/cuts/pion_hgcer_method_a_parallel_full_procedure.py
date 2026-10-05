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
        "f3": "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json",
        "f4": "Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json",
    },
}
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
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
        payload = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
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
    """Return current F.1 and exact candidate F.3/F.4 paths; never fallback."""
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
        "f3": os.path.join(root, candidate["f3"]),
        "f4": os.path.join(root, candidate["f4"]),
    }


def load_accepted_f6_3_authority(paths: Mapping[str, object]) -> tuple[list[dict[str, object]], dict[str, object], dict[str, object], dict[str, str]]:
    """Load the exact accepted authority set with its observed file hashes."""
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
    return f1, f3, f4, {**f1_hashes, "f3": f3_hash, "f4": f4_hash}


def reconstruct_transient_factor_map(
    f1_artifacts: Sequence[Mapping[str, object]],
    f3_artifact: Mapping[str, object],
    f4_artifact: Mapping[str, object],
    *,
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
        recomputed, review_data = _f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            f1_artifacts,
            f3_artifact,
            f1_input_file_hashes=hashes,
            f3_input_file_sha256=str(f3_input_file_sha256),
            accepted_f3_runtime_authority_by_kinematic=F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC,
        )
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
    if recomputed != persisted or recomputed.get("fingerprint") != persisted.get("fingerprint"):
        raise MethodAParallelFullProcedureError("f6_3_f4_shared_reproduction_mismatch")

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
