"""Detached Phase F.5 signed ``(t, phi)`` propagation of accepted F.4 output.

This module is deliberately downstream of the frozen F.4 public calculator.
It constructs review-only aggregate templates and never exposes a correction to
the normal KaonLT runtime.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os

import numpy as np

import pion_hgcer_method_a_parent_preserving_correction as _f4


METHOD_A_TPHI_PROPAGATION_SCHEMA_VERSION = "pion_hgcer_method_a_tphi_propagation/v1"
METHOD_A_TPHI_PROPAGATION_FINGERPRINT_SCHEMA_VERSION = "pion_hgcer_method_a_tphi_propagation_fingerprint/v1"
METHOD_A_TPHI_PROPAGATION_ARTIFACT_SCHEMA_VERSION = "pion_hgcer_method_a_tphi_propagation_artifact/v1"
_F4_ARTIFACT_SCHEMA = "pion_hgcer_method_a_parent_preserving_correction_artifact/v1"
_F4_SCHEMA = "pion_hgcer_method_a_parent_preserving_correction/v1"
_F4_FINGERPRINT_SCHEMA = "pion_hgcer_method_a_parent_preserving_correction_fingerprint/v1"
_CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))
_T_EDGES = (0.4, 0.5666666666666667, 0.7333333333333334, 0.9)
_PHI_EDGES = (-180.0, -140.0, -100.0, -60.0, -20.0, 20.0, 60.0, 100.0, 140.0, 180.0)
_TOLERANCE = 1.0e-12

# Source-owned acceptance boundary.  This is intentionally not discovered from
# docs, the network, an archive, or a command-line override.
ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "source_file_sha256": "adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188",
        "correction_fingerprint": "362241005c02f2149e260c391b5c3d35793287573128b42cf5ed693419d9d2f3",
        "artifact_fingerprint": "c4b9f513d5918ca77179bbe5fab28c5d9d14e63a440338545e73961fa50b9d67",
        "farm_source_head": "67e0298c51759c7a5ba693464d2c2655bf39250d",
        "f3_source_file_sha256": "04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95",
        "f3_map_fingerprint": "81b2a1e89ef9689b24c6dd145f53b26f7ac8fee9d5666cbc2ae8caa9da2830e6",
        "f3_algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
        "f3_artifact_fingerprint": "f8a12313bbd81aca48402a8c7eb4773c2dbe3ad0c7c6f64a2f212e98202d6ee2",
        "f1_source_file_sha256": {
            "Left-lowe": "f16c5e89f62e848ec221335bd4916265757040af9970b005d8f67cc800e0077d",
            "Left-highe": "203f4c76f1a251e3e8f231fa3a5e50c9a4fa337420efffffd803205c6d7ea218",
            "Center-lowe": "1593e22b55382b4a9e831d3a1114584e2e3057fcbc4aeeea4aa74c948edf39f8",
            "Center-highe": "5de64b850735ebe70040a703bac999bdd1ac84dc821aba8e1a21ec86ae4b9db3",
            "Right-highe": "77d98006ff81e772c466bf8d10ec83088509a0b8a4d430bded1172bc118bc15f",
        },
    },
}


class MethodATPhiPropagationError(ValueError):
    """The frozen F.1/F.3/F.4 authority chain or propagation is invalid."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodATPhiPropagationError("{}_not_json_safe".format(label)) from exc


def _map(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodATPhiPropagationError("{}_invalid".format(label))
    return value


def _seq(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodATPhiPropagationError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodATPhiPropagationError("{}_nonfinite".format(label))
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodATPhiPropagationError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(result):
        raise MethodATPhiPropagationError("{}_nonfinite".format(label))
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodATPhiPropagationError("{}_invalid".format(label))
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodATPhiPropagationError("{}_invalid".format(label)) from exc
    if result != value:
        raise MethodATPhiPropagationError("{}_invalid".format(label))
    return result


def _hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value.lower()):
        raise MethodATPhiPropagationError("{}_invalid".format(label))
    return value.lower()


def _close(left: float, right: float, tolerance: float = _TOLERANCE) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def _require(value: Mapping[str, object], name: str, expected: object, label: str) -> None:
    if value.get(name) != expected:
        raise MethodATPhiPropagationError("{}_{}".format(label, name))


def _authority(kinematic: str, observed_file_sha256: str, correction: Mapping[str, object], artifact_fingerprint: str, authority_by_kinematic: Mapping[str, object] | None) -> dict[str, object]:
    records = ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC if authority_by_kinematic is None else _map(authority_by_kinematic, "f4_runtime_authority_records")
    raw = records.get(kinematic)
    if raw is None:
        raise MethodATPhiPropagationError("f4_runtime_authority_kinematic_unsupported")
    accepted_raw = _map(raw, "f4_runtime_authority")
    accepted = {
        "source_file_sha256": _hash(accepted_raw.get("source_file_sha256"), "f4_runtime_authority_source_file_sha256"),
        "correction_fingerprint": _hash(accepted_raw.get("correction_fingerprint"), "f4_runtime_authority_correction_fingerprint"),
        "artifact_fingerprint": _hash(accepted_raw.get("artifact_fingerprint"), "f4_runtime_authority_artifact_fingerprint"),
        "farm_source_head": str(accepted_raw.get("farm_source_head")),
    }
    if len(accepted["farm_source_head"]) != 40 or any(char not in "0123456789abcdef" for char in accepted["farm_source_head"].lower()):
        raise MethodATPhiPropagationError("f4_runtime_authority_farm_source_head_invalid")
    observed = {
        "source_file_sha256": _hash(observed_file_sha256, "f4_observed_source_file_sha256"),
        "correction_fingerprint": _hash(correction.get("fingerprint"), "f4_observed_correction_fingerprint"),
        "artifact_fingerprint": _hash(artifact_fingerprint, "f4_observed_artifact_fingerprint"),
    }
    for name in observed:
        if accepted[name] != observed[name]:
            raise MethodATPhiPropagationError("f4_runtime_authority_{}_mismatch".format(name))
    # The default boundary pins inherited F.3/F.1 evidence too.  Test-only
    # authorities may omit these records at this builder boundary.
    for name in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint"):
        if name in accepted_raw and accepted_raw[name] != correction.get(name):
            raise MethodATPhiPropagationError("f4_runtime_authority_{}_mismatch".format(name))
    if "f1_source_file_sha256" in accepted_raw:
        expected_f1 = _map(accepted_raw["f1_source_file_sha256"], "f4_runtime_authority_f1_source_file_sha256")
        observed_f1 = {str(item["setting_id"]): item["source_file_sha256"] for item in _seq(correction.get("input_fingerprints"), "f4_input_fingerprints") if isinstance(item, Mapping)}
        if dict(expected_f1) != observed_f1:
            raise MethodATPhiPropagationError("f4_runtime_authority_f1_source_file_sha256_mismatch")
    return {"kinematic_token": kinematic, "accepted": accepted, "observed": observed, "accepted_authority_match": True}


def _validate_f4_artifact(value: object, source_file_sha256: str, authority_by_kinematic: Mapping[str, object] | None) -> tuple[Mapping[str, object], Mapping[str, object], dict[str, object]]:
    artifact = _map(value, "f4_artifact")
    _require(artifact, "schema_version", _F4_ARTIFACT_SCHEMA, "f4_artifact")
    for name, expected in (("non_authoritative", True), ("relative_map_consumed", True), ("relative_map_modified", False), ("parent_normalization_constructed", True), ("correction_constructed", True), ("correction_applied_to_production", False), ("event_correction_evaluated_for_detached_diagnostic", True), ("event_correction_persisted", False), ("child_renormalization_performed", False), ("downstream_template_application_performed", False), ("absolute_probability_constructed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False), ("basis_frozen", True), ("manual_review_required", True)):
        _require(artifact, name, expected, "f4_artifact")
    correction = _map(artifact.get("correction"), "f4_correction")
    _require(correction, "schema_version", _F4_SCHEMA, "f4_correction")
    _require(correction, "fingerprint_schema_version", _F4_FINGERPRINT_SCHEMA, "f4_correction")
    for name, expected in (("status", "available"), ("available", True), ("reason", None), ("diagnostic_stage", "complete"), ("non_authoritative", True), ("accepted_basis", "hgcer3"), ("basis_frozen", True), ("relative_map_consumed", True), ("relative_map_modified", False), ("parent_normalization_constructed", True), ("correction_constructed", True), ("correction_applied_to_production", False), ("event_correction_evaluated_for_detached_diagnostic", True), ("event_correction_persisted", False), ("child_renormalization_performed", False), ("downstream_template_application_performed", False), ("absolute_probability_constructed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False)):
        _require(correction, name, expected, "f4_correction")
    inputs = _map(correction.get("fingerprint_inputs"), "f4_fingerprint_inputs")
    if correction.get("fingerprint") != _sha256(inputs):
        raise MethodATPhiPropagationError("f4_correction_fingerprint_mismatch")
    provenance = _map(artifact.get("provenance"), "f4_provenance")
    artifact_fingerprint = _hash(artifact.get("artifact_fingerprint"), "f4_artifact_fingerprint")
    expected_wrapper = _sha256({"schema_version": _F4_ARTIFACT_SCHEMA, "correction_fingerprint": correction.get("fingerprint"), "input_paths": provenance.get("input_paths")})
    if artifact_fingerprint != expected_wrapper:
        raise MethodATPhiPropagationError("f4_artifact_fingerprint_mismatch")
    parents = _seq(correction.get("parents"), "f4_parents")
    if len(parents) != 15:
        raise MethodATPhiPropagationError("f4_parent_inventory_invalid")
    seen: set[tuple[str, int]] = set()
    for raw in parents:
        parent = _map(raw, "f4_parent"); setting = _map(parent.get("setting"), "f4_parent_setting")
        key = ("{}-{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token")), _integer(parent.get("canonical_t_index"), "f4_parent_t_index"))
        if key in seen or key[0] not in {"{}-{}".format(*item) for item in _CANONICAL_SETTINGS} or key[1] not in (0, 1, 2):
            raise MethodATPhiPropagationError("f4_parent_identity_invalid")
        seen.add(key)
        if parent.get("closure_passed") is not True:
            raise MethodATPhiPropagationError("f4_parent_closure_failed")
        for name in ("baseline_parent_sum", "raw_shape_parent_sum", "parent_normalization", "adjusted_parent_sum"):
            if _finite(parent.get(name), "f4_parent_{}".format(name)) <= 0.0:
                raise MethodATPhiPropagationError("f4_parent_{}_invalid".format(name))
        support = _map(parent.get("f3_support"), "f4_parent_f3_support")
        if support.get("statistically_sparse") is not False or support.get("support_gate_passed") is not True:
            raise MethodATPhiPropagationError("f4_parent_support_invalid")
        if not isinstance(parent.get("canonical_phi_diagnostics"), Sequence):
            raise MethodATPhiPropagationError("f4_parent_phi_diagnostics_invalid")
    expected = {("{}-{}".format(phi, epsilon), index) for phi, epsilon in _CANONICAL_SETTINGS for index in range(3)}
    if seen != expected:
        raise MethodATPhiPropagationError("f4_parent_inventory_invalid")
    kinematic = str(_map(_seq(correction.get("input_fingerprints"), "f4_input_fingerprints")[0], "f4_input_fingerprint").get("setting", {}).get("kinematic_token"))
    authority = _authority(kinematic, source_file_sha256, correction, artifact_fingerprint, authority_by_kinematic)
    return artifact, correction, authority


def _raw_by_setting(f1_artifacts: Sequence[Mapping[str, object]]) -> dict[str, Mapping[str, object]]:
    result: dict[str, Mapping[str, object]] = {}
    for raw in f1_artifacts:
        setting = _map(_map(raw, "f1_artifact").get("setting"), "f1_setting")
        key = "{}-{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"))
        if key in result:
            raise MethodATPhiPropagationError("f1_setting_duplicate")
        result[key] = raw
    return result


def _common_geometry(f1_artifacts: Sequence[Mapping[str, object]]) -> None:
    for raw in f1_artifacts:
        contract = _map(_map(raw, "f1_artifact").get("contract"), "f1_contract")
        t_edges = tuple(_finite(value, "f1_t_edge") for value in _seq(contract.get("t_edges"), "f1_t_edges"))
        phi_edges = tuple(_finite(value, "f1_phi_edge") for value in _seq(contract.get("phi_edges"), "f1_phi_edges"))
        if t_edges != _T_EDGES or phi_edges != _PHI_EDGES:
            raise MethodATPhiPropagationError("f1_common_tphi_geometry_invalid")


def _matrix(values: Sequence[Sequence[Sequence[float]]]) -> list[list[float]]:
    return [[math.fsum(cell) for cell in row] for row in values]


def _parent_phi_diagnostics(parent: Mapping[str, object]) -> dict[int, Mapping[str, object]]:
    result: dict[int, Mapping[str, object]] = {}
    for raw in _seq(parent.get("canonical_phi_diagnostics"), "f4_parent_phi_diagnostics"):
        child = _map(raw, "f4_phi_diagnostic"); index = _integer(child.get("phi_index"), "f4_phi_index")
        if index in result or index < 0 or index >= 9:
            raise MethodATPhiPropagationError("f4_phi_diagnostic_invalid")
        if _finite(child.get("phi_low"), "f4_phi_low") != _PHI_EDGES[index] or _finite(child.get("phi_high"), "f4_phi_high") != _PHI_EDGES[index + 1]:
            raise MethodATPhiPropagationError("f4_phi_diagnostic_geometry_invalid")
        result[index] = child
    return result


def _build_propagation(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None) -> dict[str, object]:
    f1_hashes = {str(key): _hash(value, "f1_input_file_sha256") for key, value in _map(f1_input_file_hashes, "f1_input_file_hashes").items()}
    _common_geometry(f1_artifacts)
    artifact, persisted, authority = _validate_f4_artifact(f4_artifact, _hash(f4_input_file_sha256, "f4_input_file_sha256"), accepted_f4_runtime_authority_by_kinematic)
    try:
        recomputed, review_data = _f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_hashes, f3_input_file_sha256=_hash(f3_input_file_sha256, "f3_input_file_sha256"), accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
        parsed = _f4._f3._validate_f1_artifacts(f1_artifacts)
        rows_by_parent = _f4._raw_application_rows(f1_artifacts, parsed, _TOLERANCE)
    except (ValueError, _f4.MethodAParentPreservingCorrectionError) as exc:
        raise MethodATPhiPropagationError("f4_shared_reproduction_failed:{}".format(exc)) from exc
    if recomputed != persisted:
        raise MethodATPhiPropagationError("f4_shared_reproduction_mismatch")
    if recomputed.get("fingerprint") != persisted.get("fingerprint"):
        raise MethodATPhiPropagationError("f4_shared_fingerprint_mismatch")
    persisted_parents = {(str(row["setting_id"]), int(row["canonical_t_index"])): _map(row, "f4_parent") for row in _seq(persisted.get("parents"), "f4_parents")}
    review_by_parent = {(str(row.get("setting_id")), _integer(row.get("canonical_t_index"), "f4_review_t_index")): _map(row, "f4_review") for row in review_data}
    if len(review_by_parent) != 15 or set(review_by_parent) != set(persisted_parents):
        raise MethodATPhiPropagationError("f4_review_inventory_invalid")
    raw_settings = _raw_by_setting(f1_artifacts)
    setting_payloads: list[dict[str, object]] = []
    all_parent_summaries: list[dict[str, object]] = []
    for phi, epsilon in _CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon)
        raw = raw_settings.get(setting_id)
        if raw is None:
            raise MethodATPhiPropagationError("f1_setting_missing")
        setting = _map(raw.get("setting"), "f1_setting")
        count = [[0 for _ in range(9)] for _ in range(3)]
        baseline_values: list[list[list[float]]] = [[[ ] for _ in range(9)] for _ in range(3)]
        adjusted_values: list[list[list[float]]] = [[[ ] for _ in range(9)] for _ in range(3)]
        parent_summaries: list[dict[str, object]] = []
        for t_index in range(3):
            key = (setting_id, t_index); parent = persisted_parents.get(key); review = review_by_parent.get(key); rows = rows_by_parent.get(key)
            if parent is None or review is None or rows is None:
                raise MethodATPhiPropagationError("f4_parent_or_f1_rows_missing")
            factors = np.asarray(review.get("correction_factors"), dtype=float)
            if factors.ndim != 1 or factors.size != len(rows) or not np.all(np.isfinite(factors)) or np.any(factors <= 0.0):
                raise MethodATPhiPropagationError("f4_review_factor_alignment_invalid")
            source_values: dict[str, list[tuple[float, float]]] = {}
            # Exact row/factor length equality is enforced immediately above;
            # plain zip preserves that validated pairing on farm Python versions
            # without zip(strict=...).
            for row, factor in zip(rows, factors):
                phi_index = _integer(row.get("phi_index"), "f1_phi_index")
                if phi_index < 0 or phi_index >= 9:
                    raise MethodATPhiPropagationError("f1_phi_index_invalid")
                baseline = _finite(row.get("signed_baseline_event_contribution"), "f1_signed_baseline")
                adjusted_value = baseline * float(factor)
                count[t_index][phi_index] += 1; baseline_values[t_index][phi_index].append(baseline); adjusted_values[t_index][phi_index].append(adjusted_value)
                source_values.setdefault(str(row.get("source_label")), []).append((baseline, adjusted_value))
            baseline = _matrix(baseline_values)[t_index]; adjusted = _matrix(adjusted_values)[t_index]
            delta = [adjusted[index] - baseline[index] for index in range(9)]
            parent_baseline = math.fsum(baseline); parent_adjusted = math.fsum(adjusted)
            f4_baseline = _finite(parent.get("baseline_parent_sum"), "f4_parent_baseline")
            f4_adjusted = _finite(parent.get("adjusted_parent_sum"), "f4_parent_adjusted")
            if not _close(parent_baseline, f4_baseline) or not _close(parent_adjusted, f4_adjusted) or not _close(math.fsum(delta), 0.0):
                raise MethodATPhiPropagationError("f4_parent_template_closure_mismatch")
            f4_sources: dict[str, Mapping[str, object]] = {}
            for raw_source in _seq(parent.get("source_diagnostics"), "f4_source_diagnostics"):
                source = _map(raw_source, "f4_source_diagnostic"); label = str(source.get("source_label"))
                if label in f4_sources:
                    raise MethodATPhiPropagationError("f4_source_diagnostic_duplicate")
                f4_sources[label] = source
            if set(f4_sources) != set(source_values):
                raise MethodATPhiPropagationError("f4_source_diagnostic_inventory_mismatch")
            for label, values in source_values.items():
                source = f4_sources[label]; source_baseline = math.fsum(item[0] for item in values); source_adjusted = math.fsum(item[1] for item in values)
                if _integer(source.get("event_count"), "f4_source_event_count") != len(values) or not _close(_finite(source.get("baseline_signed_sum"), "f4_source_baseline"), source_baseline) or not _close(_finite(source.get("adjusted_signed_sum"), "f4_source_adjusted"), source_adjusted) or not _close(_finite(source.get("signed_delta"), "f4_source_delta"), source_adjusted - source_baseline):
                    raise MethodATPhiPropagationError("f4_source_diagnostic_mismatch")
            diagnostics = _parent_phi_diagnostics(parent)
            for phi_index in range(9):
                child = diagnostics.get(phi_index)
                if count[t_index][phi_index] == 0:
                    if child is not None:
                        raise MethodATPhiPropagationError("f4_empty_phi_diagnostic_present")
                    continue
                if child is None:
                    raise MethodATPhiPropagationError("f4_occupied_phi_diagnostic_missing")
                if _integer(child.get("event_count"), "f4_phi_event_count") != count[t_index][phi_index] or not _close(_finite(child.get("baseline_signed_sum"), "f4_phi_baseline"), baseline[phi_index]) or not _close(_finite(child.get("adjusted_signed_sum"), "f4_phi_adjusted"), adjusted[phi_index]) or not _close(_finite(child.get("signed_delta"), "f4_phi_delta"), delta[phi_index]):
                    raise MethodATPhiPropagationError("f4_phi_diagnostic_mismatch")
            redistribution = [item / f4_baseline for item in delta]
            parent_summaries.append({"canonical_t_index": t_index, "baseline_parent_sum": parent_baseline, "adjusted_parent_sum": parent_adjusted, "closure_residual": parent_adjusted - parent_baseline, "max_abs_child_delta": max(abs(item) for item in delta), "sum_abs_child_delta": math.fsum(abs(item) for item in delta), "max_abs_redistribution": max(abs(item) for item in redistribution), "sum_abs_redistribution": math.fsum(abs(item) for item in redistribution), "max_abs_redistribution_phi_index": max(range(9), key=lambda index: abs(redistribution[index]))})
        baseline_matrix = _matrix(baseline_values); adjusted_matrix = _matrix(adjusted_values)
        delta_matrix = [[adjusted_matrix[t][p] - baseline_matrix[t][p] for p in range(9)] for t in range(3)]
        base_share = [[baseline_matrix[t][p] / parent_summaries[t]["baseline_parent_sum"] for p in range(9)] for t in range(3)]
        adjusted_share = [[adjusted_matrix[t][p] / parent_summaries[t]["adjusted_parent_sum"] for p in range(9)] for t in range(3)]
        redistribution = [[delta_matrix[t][p] / parent_summaries[t]["baseline_parent_sum"] for p in range(9)] for t in range(3)]
        setting_payloads.append({"setting": _copy(setting, "setting"), "setting_id": setting_id, "t_edges": list(_T_EDGES), "phi_edges": list(_PHI_EDGES), "event_counts": count, "baseline_signed_contents": baseline_matrix, "adjusted_signed_contents": adjusted_matrix, "signed_delta_contents": delta_matrix, "baseline_share_of_parent": base_share, "adjusted_share_of_parent": adjusted_share, "redistribution_fraction_of_parent": redistribution, "parents": parent_summaries, "setting_baseline_sum": math.fsum(item["baseline_parent_sum"] for item in parent_summaries), "setting_adjusted_sum": math.fsum(item["adjusted_parent_sum"] for item in parent_summaries), "setting_closure_residual": math.fsum(item["closure_residual"] for item in parent_summaries), "setting_max_abs_redistribution": max(item["max_abs_redistribution"] for item in parent_summaries)})
        all_parent_summaries.extend({"setting_id": setting_id, **item} for item in parent_summaries)
    input_fingerprints = _copy(persisted.get("input_fingerprints"), "f4_input_fingerprints")
    core: dict[str, object] = {"schema_version": METHOD_A_TPHI_PROPAGATION_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_TPHI_PROPAGATION_FINGERPRINT_SCHEMA_VERSION, "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete", "non_authoritative": True, "accepted_basis": "hgcer3", "basis_frozen": True, "f4_correction_consumed": True, "f4_correction_modified": False, "event_level_propagation_performed": True, "event_correction_persisted": False, "baseline_template_constructed": True, "adjusted_template_constructed": True, "parallel_templates_constructed": True, "canonical_child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "manual_review_required": True, "f4_source_file_sha256": _hash(f4_input_file_sha256, "f4_input_file_sha256"), "f4_correction_fingerprint": persisted["fingerprint"], "f4_artifact_fingerprint": artifact["artifact_fingerprint"], "f4_runtime_authority": authority, "f3_inherited_authority": {name: persisted[name] for name in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint", "f3_runtime_authority")}, "input_fingerprints": input_fingerprints, "t_edges": list(_T_EDGES), "phi_edges": list(_PHI_EDGES), "setting_templates": setting_payloads, "parent_closure": all_parent_summaries}
    fingerprint_inputs = {"schema_version": core["schema_version"], "fingerprint_schema_version": core["fingerprint_schema_version"], "f4_source_file_sha256": core["f4_source_file_sha256"], "f4_correction_fingerprint": core["f4_correction_fingerprint"], "f4_artifact_fingerprint": core["f4_artifact_fingerprint"], "f4_runtime_authority": core["f4_runtime_authority"], "f3_inherited_authority": core["f3_inherited_authority"], "input_content": [{"setting_id": row["setting_id"], "stable_f1_content_fingerprint": row["stable_f1_content_fingerprint"], "source_file_sha256": row["source_file_sha256"]} for row in input_fingerprints], "t_edges": core["t_edges"], "phi_edges": core["phi_edges"], "setting_templates": core["setting_templates"], "parent_closure": core["parent_closure"], "policy": {"formula": "sum_b_and_sum_b_times_f4_C_by_canonical_tphi_cell/v1", "child_normalization": "forbidden/v1", "f4_reuse": "public_shared_calculator_required/v1", "closure_relative_tolerance": _TOLERANCE}}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    return _copy(core, "tphi_propagation")  # type: ignore[return-value]


def build_pion_hgcer_method_a_tphi_propagation(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> dict[str, object]:
    """Build F.5 aggregate templates by invoking F.4's public calculator once."""
    return _build_propagation(f1_artifacts, f3_artifact, f4_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, f4_input_file_sha256=f4_input_file_sha256, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)


def build_pion_hgcer_method_a_tphi_propagation_artifact(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, input_paths: Mapping[str, object] | None = None, generated_at_utc: str | None = None, git_head: str | None = None, git_status_short: str | None = None, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> dict[str, object]:
    propagation = build_pion_hgcer_method_a_tphi_propagation(f1_artifacts, f3_artifact, f4_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, f4_input_file_sha256=f4_input_file_sha256, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    provenance = {"generated_at_utc": generated_at_utc, "git_head": git_head, "git_status_short": git_status_short, "input_paths": _copy({} if input_paths is None else input_paths, "input_paths")}
    return {"schema_version": METHOD_A_TPHI_PROPAGATION_ARTIFACT_SCHEMA_VERSION, "propagation": propagation, "non_authoritative": True, "f4_correction_consumed": True, "f4_correction_modified": False, "event_level_propagation_performed": True, "event_correction_persisted": False, "baseline_template_constructed": True, "adjusted_template_constructed": True, "parallel_templates_constructed": True, "canonical_child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "basis_frozen": True, "manual_review_required": True, "provenance": provenance, "artifact_fingerprint": _sha256({"schema_version": METHOD_A_TPHI_PROPAGATION_ARTIFACT_SCHEMA_VERSION, "propagation_fingerprint": propagation["fingerprint"], "input_paths": provenance["input_paths"]})}


def pion_hgcer_method_a_tphi_propagation_filename(kinematic: object) -> str:
    token = str(kinematic)
    if not token or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise MethodATPhiPropagationError("tphi_propagation_filename_token_invalid")
    return "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.json".format(token)


def write_pion_hgcer_method_a_tphi_propagation_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    with open(os.fspath(path), "w", encoding="utf-8") as handle:
        json.dump(_copy(payload, "tphi_propagation_artifact"), handle, sort_keys=True, indent=2, allow_nan=False); handle.write("\n")
    return os.fspath(path)
