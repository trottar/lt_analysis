"""Fail-closed transient F.4 authority for the private F.6.3 branch.

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
    """Return deterministic accepted F.1/F.3/F.4 paths; never search/fallback."""
    root = os.fspath(outpath)
    token = str(kinematic_token)
    _f4, _f5, f1_filename, f3_filename, f4_filename = _runtime_dependencies()
    f1 = {
        "{}-{}".format(phi, epsilon): os.path.join(
            root,
            f1_filename(phi, "kaon", token, epsilon),
        )
        for phi, epsilon in _CANONICAL_SETTINGS
    }
    return {
        "f1": f1,
        "f3": os.path.join(root, f3_filename(token)),
        "f4": os.path.join(root, f4_filename(token)),
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
            f4_artifact, str(f4_input_file_sha256), None,
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
        )
        rows_by_parent = _f4._raw_application_rows(f1_artifacts, parsed, _TOLERANCE)
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
                    "t_index", "phi_index", "analysis_MM", "analysis_t",
                    "signed_source_coefficient", "baseline_pion_weight_w0",
                    "signed_baseline_event_contribution",
                )
            }
    if not factors:
        raise MethodAParallelFullProcedureError("f6_3_selected_factor_population_missing")
    identities = sorted((source, entry) for source, entry in factors)
    provenance = {
        "schema_version": F6_3_PARALLEL_SCHEMA_VERSION,
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
