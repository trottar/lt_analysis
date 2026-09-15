"""Detached Phase F.4 parent-preserving Method-A-only correction candidate.

Consumes frozen F.1 application contracts and an accepted F.3 relative map.
All event-level values are transient; JSON persists only parent/source/child
aggregates.  Nothing in this module is a production correction.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os

import numpy as np
from scipy.spatial import cKDTree

import pion_hgcer_method_a_acceptance_map as _f3


METHOD_A_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION = "pion_hgcer_method_a_parent_preserving_correction/v1"
METHOD_A_PARENT_PRESERVING_CORRECTION_FINGERPRINT_SCHEMA_VERSION = "pion_hgcer_method_a_parent_preserving_correction_fingerprint/v1"
METHOD_A_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION = "pion_hgcer_method_a_parent_preserving_correction_artifact/v1"
_F3_ARTIFACT_SCHEMA = "pion_hgcer_method_a_acceptance_map_artifact/v1"
_F3_SCHEMA = "pion_hgcer_method_a_acceptance_map/v1"
_F3_FINGERPRINT_SCHEMA = "pion_hgcer_method_a_acceptance_map_fingerprint/v1"
_CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))
_FEATURES = ("SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer")
DEFAULT_ALGORITHM_CONFIG = {
    "algorithm_version": "method_a_parent_preserving_correction/v1",
    "correction_formula_version": "signed_parent_preserving_relative_shape/v1",
    "ood_raw_shape_policy": "neutral_one_outside_f3_p99_support/v1",
    "parent_normalization_policy": "full_signed_canonical_t_application_population/v1",
    "closure_relative_tolerance": 1.0e-12,
    "f3_support_absolute_tolerance": 1.0e-12,
    "signed_baseline_identity_relative_tolerance": 1.0e-12,
}

# This is a source-owned acceptance boundary, not a dynamic discovery rule.
# F.4 may consume only the F.3 artifact accepted by the detached farm gate.
ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "source_file_sha256": "04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95",
        "map_fingerprint": "81b2a1e89ef9689b24c6dd145f53b26f7ac8fee9d5666cbc2ae8caa9da2830e6",
        "algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
        "artifact_fingerprint": "f8a12313bbd81aca48402a8c7eb4773c2dbe3ad0c7c6f64a2f212e98202d6ee2",
        "farm_source_head": "5382cfc1994b078c620b32c043938134c33ffa39",
    },
}


class MethodAParentPreservingCorrectionError(ValueError):
    """Frozen F.1/F.3 authority or correction mathematics is invalid."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodAParentPreservingCorrectionError("{}_not_json_safe".format(label)) from exc


def _map(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    return value


def _seq(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodAParentPreservingCorrectionError("{}_nonfinite".format(label))
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodAParentPreservingCorrectionError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(result):
        raise MethodAParentPreservingCorrectionError("{}_nonfinite".format(label))
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label)) from exc
    if result != value:
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    return result


def _require(value: Mapping[str, object], name: str, expected: object, label: str) -> None:
    if value.get(name) != expected:
        raise MethodAParentPreservingCorrectionError("{}_{}".format(label, name))


def _percentiles(values: Sequence[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise MethodAParentPreservingCorrectionError("correction_summary_invalid")
    return {"min": float(np.min(array)), "p01": float(np.percentile(array, 1.0)), "p50": float(np.percentile(array, 50.0)), "p99": float(np.percentile(array, 99.0)), "max": float(np.max(array))}


def _close(left: float, right: float, tolerance: float) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def _resolve_config(value: Mapping[str, object] | None) -> dict[str, object]:
    config = dict(DEFAULT_ALGORITHM_CONFIG)
    if value is not None:
        supplied = _map(value, "algorithm_config")
        unexpected = sorted(set(supplied) - set(config))
        if unexpected:
            raise MethodAParentPreservingCorrectionError("algorithm_config_unknown:{}".format(",".join(unexpected)))
        config.update(supplied)
    for name in ("closure_relative_tolerance", "f3_support_absolute_tolerance", "signed_baseline_identity_relative_tolerance"):
        if _finite(config[name], "algorithm_config_{}".format(name)) <= 0.0:
            raise MethodAParentPreservingCorrectionError("algorithm_config_{}_invalid".format(name))
    for name, expected in (("correction_formula_version", "signed_parent_preserving_relative_shape/v1"), ("ood_raw_shape_policy", "neutral_one_outside_f3_p99_support/v1"), ("parent_normalization_policy", "full_signed_canonical_t_application_population/v1")):
        if config.get(name) != expected:
            raise MethodAParentPreservingCorrectionError("algorithm_config_{}_invalid".format(name))
    return _copy(config, "algorithm_config")  # type: ignore[return-value]


def _validate_hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value.lower()):
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    return value.lower()


def _strict_edges(value: object, label: str) -> tuple[float, ...]:
    edges = tuple(_finite(item, "{}_edge".format(label)) for item in _seq(value, label))
    if len(edges) < 2 or any(left >= right for left, right in zip(edges, edges[1:])):
        raise MethodAParentPreservingCorrectionError("{}_invalid".format(label))
    return edges


def _f3_runtime_authority(kinematic: str, observed_source_file_sha256: str, map_value: Mapping[str, object], artifact_fingerprint: str, authority_by_kinematic: Mapping[str, object] | None) -> dict[str, object]:
    authority_records = ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC if authority_by_kinematic is None else _map(authority_by_kinematic, "f3_runtime_authority_records")
    raw_authority = authority_records.get(kinematic)
    if raw_authority is None:
        raise MethodAParentPreservingCorrectionError("f3_runtime_authority_kinematic_unsupported")
    authority = _map(raw_authority, "f3_runtime_authority")
    accepted = {
        "source_file_sha256": _validate_hash(authority.get("source_file_sha256"), "f3_runtime_authority_source_file_sha256"),
        "map_fingerprint": _validate_hash(authority.get("map_fingerprint"), "f3_runtime_authority_map_fingerprint"),
        "algorithm_fingerprint": _validate_hash(authority.get("algorithm_fingerprint"), "f3_runtime_authority_algorithm_fingerprint"),
        "artifact_fingerprint": _validate_hash(authority.get("artifact_fingerprint"), "f3_runtime_authority_artifact_fingerprint"),
        "farm_source_head": str(authority.get("farm_source_head")),
    }
    if len(accepted["farm_source_head"]) != 40 or any(char not in "0123456789abcdef" for char in accepted["farm_source_head"].lower()):
        raise MethodAParentPreservingCorrectionError("f3_runtime_authority_farm_source_head_invalid")
    observed = {
        "source_file_sha256": observed_source_file_sha256,
        "map_fingerprint": _validate_hash(map_value.get("fingerprint"), "f3_observed_map_fingerprint"),
        "algorithm_fingerprint": _validate_hash(map_value.get("algorithm_fingerprint"), "f3_observed_algorithm_fingerprint"),
        "artifact_fingerprint": artifact_fingerprint,
    }
    for name in ("source_file_sha256", "map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"):
        if observed[name] != accepted[name]:
            raise MethodAParentPreservingCorrectionError("f3_runtime_authority_{}_mismatch".format(name))
    return {"kinematic_token": kinematic, "accepted": accepted, "observed": observed, "accepted_authority_match": True}


def _preflight_raw_phi_degrees(f1_artifacts: Sequence[Mapping[str, object]]) -> None:
    """Fail with an F.4 semantic error before canonical JSON rejects NaN/Inf."""
    for artifact_index, raw_artifact in enumerate(f1_artifacts):
        artifact = _map(raw_artifact, "f1_artifact_{}".format(artifact_index))
        contract = _map(artifact.get("contract"), "f1_contract_{}".format(artifact_index))
        for row_index, raw_row in enumerate(_seq(contract.get("application_records"), "f1_application_records_{}".format(artifact_index))):
            row = _map(raw_row, "f1_application_{}_{}".format(artifact_index, row_index))
            _finite(row.get("phi_degrees"), "f1_application_phi_degrees")


def _raw_application_rows(f1_artifacts: Sequence[Mapping[str, object]], parsed: Sequence[Mapping[str, object]], tolerance: float) -> dict[tuple[str, int], list[dict[str, object]]]:
    """Add F.4-only authoritative signed/phi fields after F.3's full F.1 audit."""
    raw_by_setting: dict[str, Mapping[str, object]] = {}
    for raw in f1_artifacts:
        setting = _map(raw.get("setting"), "f1_setting")
        raw_by_setting["{}-{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"))] = raw
    result: dict[tuple[str, int], list[dict[str, object]]] = {}
    for item in parsed:
        setting_id = str(item["setting_id"])
        raw_artifact = raw_by_setting.get(setting_id)
        if raw_artifact is None:
            raise MethodAParentPreservingCorrectionError("f1_raw_setting_missing")
        contract = _map(raw_artifact.get("contract"), "f1_contract")
        phi_edges = _strict_edges(contract.get("phi_edges"), "f1_phi_edges")
        raw_rows = _seq(contract.get("application_records"), "f1_application_records")
        sanitized = _seq(item["application"], "f1_application")
        if len(raw_rows) != len(sanitized):
            raise MethodAParentPreservingCorrectionError("f1_application_length_mismatch")
        by_identity: dict[tuple[str, int], Mapping[str, object]] = {}
        for index, raw_row in enumerate(raw_rows):
            row = _map(raw_row, "f1_application_{}".format(index))
            source = row.get("source_label"); entry = row.get("entry_index")
            if not isinstance(source, str):
                raise MethodAParentPreservingCorrectionError("f1_application_identity_invalid")
            identity = (source, _integer(entry, "f1_application_entry_index"))
            if identity in by_identity:
                raise MethodAParentPreservingCorrectionError("f1_application_identity_duplicate")
            by_identity[identity] = row
        for clean in sanitized:
            row = _map(clean, "f1_application")
            identity = (str(row["source_label"]), int(row["entry_index"]))
            raw = by_identity.get(identity)
            if raw is None:
                raise MethodAParentPreservingCorrectionError("f1_application_identity_missing")
            if raw.get("phi_status") != "inside_phi":
                raise MethodAParentPreservingCorrectionError("f1_application_phi_status_invalid")
            phi_index = _integer(raw.get("phi_index"), "f1_application_phi_index")
            phi_low = _finite(raw.get("phi_low"), "f1_application_phi_low")
            phi_high = _finite(raw.get("phi_high"), "f1_application_phi_high")
            phi_degrees = _finite(raw.get("phi_degrees"), "f1_application_phi_degrees")
            if phi_index < 0 or phi_index >= len(phi_edges) - 1 or phi_low != phi_edges[phi_index] or phi_high != phi_edges[phi_index + 1]:
                raise MethodAParentPreservingCorrectionError("f1_application_phi_geometry_invalid")
            if not phi_low <= phi_degrees < phi_high and not (phi_index == len(phi_edges) - 2 and phi_degrees == phi_high):
                raise MethodAParentPreservingCorrectionError("f1_application_phi_assignment_invalid")
            if str(raw.get("source_label")) != str(row["source_label"]) or _integer(raw.get("entry_index"), "f1_application_raw_entry_index") != int(row["entry_index"]) or _integer(raw.get("t_index"), "f1_application_raw_t_index") != int(row["t_index"]) or _finite(raw.get("t_low"), "f1_application_raw_t_low") != float(row["t_low"]) or _finite(raw.get("t_high"), "f1_application_raw_t_high") != float(row["t_high"]):
                raise MethodAParentPreservingCorrectionError("f1_application_identity_assignment_mismatch")
            w0 = _finite(raw.get("baseline_pion_weight_w0"), "f1_application_w0")
            coefficient = _finite(raw.get("signed_source_coefficient"), "f1_application_source_coefficient")
            baseline = _finite(raw.get("signed_baseline_event_contribution"), "f1_application_baseline_contribution")
            if w0 < 0.0:
                raise MethodAParentPreservingCorrectionError("f1_application_w0_negative")
            if not _close(baseline, coefficient * w0, tolerance):
                raise MethodAParentPreservingCorrectionError("f1_application_signed_baseline_identity_mismatch")
            combined = {**row, "phi_index": phi_index, "phi_low": phi_low, "phi_high": phi_high, "phi_degrees": phi_degrees, "phi_status": "inside_phi", "baseline_pion_weight_w0": w0, "signed_source_coefficient": coefficient, "signed_baseline_event_contribution": baseline}
            result.setdefault((setting_id, int(row["t_index"])), []).append(combined)
    for rows in result.values():
        rows.sort(key=lambda row: (str(row["source_label"]), int(row["entry_index"])))
    return result


def _validate_f3_artifact(value: object, parsed_f1: Sequence[Mapping[str, object]], f1_hashes: Mapping[str, str], tolerance: float) -> tuple[Mapping[str, object], dict[tuple[str, int], Mapping[str, object]]]:
    artifact = _map(value, "f3_artifact")
    _require(artifact, "schema_version", _F3_ARTIFACT_SCHEMA, "f3_artifact")
    for name, expected in (("non_authoritative", True), ("map_constructed", True), ("map_applied", False), ("absolute_probability_constructed", False), ("parent_normalization_constructed", False), ("correction_constructed", False), ("event_probability_persisted", False), ("event_application_performed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False), ("basis_frozen", True), ("manual_review_required", True)):
        _require(artifact, name, expected, "f3_artifact")
    provenance = _map(artifact.get("provenance"), "f3_provenance")
    map_value = _map(artifact.get("acceptance_map"), "f3_map")
    if artifact.get("artifact_fingerprint") != _sha256({"schema_version": _F3_ARTIFACT_SCHEMA, "map_fingerprint": map_value.get("fingerprint"), "input_paths": provenance.get("input_paths")}):
        raise MethodAParentPreservingCorrectionError("f3_artifact_fingerprint_mismatch")
    _require(map_value, "schema_version", _F3_SCHEMA, "f3_map"); _require(map_value, "fingerprint_schema_version", _F3_FINGERPRINT_SCHEMA, "f3_map")
    for name, expected in (("status", "available"), ("available", True), ("reason", None), ("diagnostic_stage", "complete"), ("non_authoritative", True), ("map_constructed", True), ("map_applied", False), ("absolute_probability_constructed", False), ("parent_normalization_constructed", False), ("correction_constructed", False), ("event_probability_persisted", False), ("event_application_performed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False), ("accepted_basis", "hgcer3"), ("ordered_features", list(_FEATURES)), ("basis_frozen", True), ("manual_review_required", True)):
        _require(map_value, name, expected, "f3_map")
    inputs = _map(map_value.get("fingerprint_inputs"), "f3_fingerprint_inputs")
    if map_value.get("fingerprint") != _sha256(inputs):
        raise MethodAParentPreservingCorrectionError("f3_map_fingerprint_mismatch")
    if map_value.get("algorithm_config") != _f3.DEFAULT_ALGORITHM_CONFIG:
        raise MethodAParentPreservingCorrectionError("f3_algorithm_config_invalid")
    for name in ("algorithm_config", "algorithm_fingerprint", "f2_representation_fingerprint", "f2_algorithm_fingerprint", "f2_source_file_sha256", "models"):
        if inputs.get(name) != map_value.get(name):
            raise MethodAParentPreservingCorrectionError("f3_fingerprint_content_mismatch:{}".format(name))
    expected_algorithm = _sha256({"map_schema_version": _F3_SCHEMA, "fingerprint_schema_version": _F3_FINGERPRINT_SCHEMA, "accepted_basis": "hgcer3", "ordered_features": list(_FEATURES), "algorithm_config": map_value.get("algorithm_config")})
    if map_value.get("algorithm_fingerprint") != expected_algorithm:
        raise MethodAParentPreservingCorrectionError("f3_algorithm_fingerprint_mismatch")
    f1_by_setting = {str(item["setting_id"]): item for item in parsed_f1}
    f3_inputs = _seq(map_value.get("input_fingerprints"), "f3_input_fingerprints")
    if len(f3_inputs) != len(f1_by_setting):
        raise MethodAParentPreservingCorrectionError("f3_input_inventory_invalid")
    expected_input_content = [{"setting_id": str(item["setting_id"]), "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "source_file_sha256": f1_hashes[str(item["setting_id"])]} for item in parsed_f1]
    if inputs.get("input_content") != expected_input_content:
        raise MethodAParentPreservingCorrectionError("f3_fingerprint_input_content_mismatch")
    observed_input_settings: set[str] = set()
    for raw_input in f3_inputs:
        item = _map(raw_input, "f3_input_fingerprint"); setting_id = str(item.get("setting_id")); source = f1_by_setting.get(setting_id)
        if setting_id in observed_input_settings or source is None or item.get("setting") != source["setting"] or item.get("source_file_sha256") != f1_hashes.get(setting_id) or item.get("stable_f1_content_fingerprint") != source["stable_content_fingerprint"] or item.get("f1_contract_fingerprint") != _map(source["fingerprints"], "f1_fingerprints").get("fingerprint"):
            raise MethodAParentPreservingCorrectionError("f3_f1_authority_mismatch")
        observed_input_settings.add(setting_id)
    if observed_input_settings != set(f1_by_setting):
        raise MethodAParentPreservingCorrectionError("f3_input_inventory_invalid")
    models = _seq(map_value.get("models"), "f3_models")
    if len(models) != 15:
        raise MethodAParentPreservingCorrectionError("f3_model_count_invalid")
    indexed: dict[tuple[str, int], Mapping[str, object]] = {}
    for raw_model in models:
        model = _map(raw_model, "f3_model"); setting_id = str(model.get("setting_id")); t_index = _integer(model.get("canonical_t_index"), "f3_t_index")
        source = f1_by_setting.get(setting_id)
        if source is None or (setting_id, t_index) in indexed or model.get("setting") != source["setting"] or model.get("ordered_features") != list(_FEATURES):
            raise MethodAParentPreservingCorrectionError("f3_model_identity_invalid")
        geometry = _seq(source["canonical_t_geometry"], "f1_geometry")
        if t_index < 0 or t_index >= len(geometry) or model.get("canonical_t_low") != geometry[t_index][1] or model.get("canonical_t_high") != geometry[t_index][2]:
            raise MethodAParentPreservingCorrectionError("f3_model_geometry_invalid")
        fit = _map(model.get("logistic_fit"), "f3_logistic_fit"); scaler = _map(model.get("scaler"), "f3_scaler"); support = _map(model.get("application_support"), "f3_support"); continuity = _map(model.get("f2_support_continuity"), "f3_f2_continuity")
        if fit.get("converged") is not True or fit.get("optimizer") != _f3.DEFAULT_ALGORITHM_CONFIG["optimizer"] or fit.get("regularization_lambda") != _f3.DEFAULT_ALGORITHM_CONFIG["regularization_lambda"] or fit.get("optimizer_maxiter") != _f3.DEFAULT_ALGORITHM_CONFIG["optimizer_maxiter"] or fit.get("optimizer_gtol") != _f3.DEFAULT_ALGORITHM_CONFIG["optimizer_gtol"] or scaler.get("method") != _f3.DEFAULT_ALGORITHM_CONFIG["scaling"] or model.get("relative_response_definition") != "exp(coefficients_dot_robust_scaled_hgcer3_without_intercept)" or model.get("relative_response_at_training_median") != 1.0 or support.get("valid") is not True or support.get("statistically_sparse") is not False or support.get("support_gate_passed") is not True or continuity.get("passed") is not True or continuity.get("absolute_tolerance") != 1.0e-12 or continuity.get("mismatched_fields") != []:
            raise MethodAParentPreservingCorrectionError("f3_model_gate_invalid")
        for name in ("median", "divisor"):
            values = _seq(scaler.get(name), "f3_scaler_{}".format(name))
            if len(values) != 3 or any(not math.isfinite(float(item)) for item in values):
                raise MethodAParentPreservingCorrectionError("f3_scaler_invalid")
        if any(float(item) <= 0.0 for item in _seq(scaler.get("divisor"), "f3_divisor")) or len(_seq(fit.get("coefficients"), "f3_coefficients")) != 3 or any(not math.isfinite(float(item)) for item in _seq(fit.get("coefficients"), "f3_coefficients")) or not math.isfinite(_finite(fit.get("intercept_provenance_only"), "f3_intercept")) or not math.isfinite(_finite(fit.get("objective"), "f3_objective")):
            raise MethodAParentPreservingCorrectionError("f3_model_parameters_invalid")
        for name in ("training_nn_p50", "training_nn_p95", "training_nn_p99", "application_nn_p50", "application_nn_p95", "application_nn_p99", "application_nn_max", "application_ood_fraction", "support_distance_threshold"):
            _finite(support.get(name), "f3_support_{}".format(name))
        indexed[(setting_id, t_index)] = model
    expected_keys = {(str(item["setting_id"]), int(index)) for item in parsed_f1 for index, _, _ in _seq(item["canonical_t_geometry"], "geometry")}
    if set(indexed) != expected_keys:
        raise MethodAParentPreservingCorrectionError("f3_model_inventory_invalid")
    return map_value, indexed


def _support_continuity(training: Sequence[Mapping[str, object]], application: Sequence[Mapping[str, object]], model: Mapping[str, object], tolerance: float) -> tuple[cKDTree, np.ndarray, np.ndarray, float, dict[str, object]]:
    scaler = _map(model.get("scaler"), "f3_scaler"); support = _map(model.get("application_support"), "f3_support")
    values = _f3._feature_array(training, _FEATURES); recomputed = _f3._robust_scale(values)
    if recomputed is None:
        raise MethodAParentPreservingCorrectionError("f3_support_scaling_invalid")
    median, divisor = recomputed
    stored_median = np.asarray(scaler.get("median"), dtype=float); stored_divisor = np.asarray(scaler.get("divisor"), dtype=float)
    if stored_median.shape != (3,) or stored_divisor.shape != (3,) or not np.allclose(median, stored_median, rtol=0.0, atol=tolerance) or not np.allclose(divisor, stored_divisor, rtol=0.0, atol=tolerance):
        raise MethodAParentPreservingCorrectionError("f3_scaler_continuity_failed")
    tree = cKDTree((values - median) / divisor)
    distances = np.asarray(tree.query((values - median) / divisor, k=2)[0][:, 1], dtype=float)
    threshold = float(np.percentile(distances, 99.0))
    nonprompt = [row for row in application if row["source_label"] != "prompt"]
    app = np.asarray(tree.query((_f3._feature_array(nonprompt, _FEATURES) - median) / divisor, k=1)[0], dtype=float) if nonprompt else np.asarray([], dtype=float)
    summary = _f3._percentiles(distances); application_summary = _f3._percentiles(app)
    observed = {"nonprompt_application_count": int(app.size), "training_nn_p50": summary["p50"], "training_nn_p95": summary["p95"], "training_nn_p99": summary["p99"], "application_nn_p50": application_summary["p50"], "application_nn_p95": application_summary["p95"], "application_nn_p99": application_summary["p99"], "application_nn_max": application_summary["max"], "application_ood_count": int(np.sum(app > threshold)), "application_ood_fraction": float(np.sum(app > threshold) / app.size) if app.size else 0.0, "statistically_sparse": int(app.size) < 20, "support_gate_passed": bool((float(np.sum(app > threshold) / app.size) if app.size else 0.0) <= 0.10), "support_distance_threshold": threshold}
    for name, actual in observed.items():
        stored = support.get(name)
        if isinstance(actual, float):
            if stored is None or abs(actual - float(stored)) > tolerance:
                raise MethodAParentPreservingCorrectionError("f3_support_continuity_failed:{}".format(name))
        elif stored != actual:
            raise MethodAParentPreservingCorrectionError("f3_support_continuity_failed:{}".format(name))
    return tree, median, divisor, threshold, observed


def _unavailable(parsed: Sequence[Mapping[str, object]], map_value: Mapping[str, object], f1_hashes: Mapping[str, str], f3_sha256: str, f3_artifact_fingerprint: str, f3_runtime_authority: Mapping[str, object], config: Mapping[str, object], invalid: Sequence[Mapping[str, object]]) -> dict[str, object]:
    inputs = [{"setting": _copy(item["setting"], "setting"), "setting_id": item["setting_id"], "source_file_sha256": f1_hashes[item["setting_id"]], "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "f1_contract_fingerprint": _map(item["fingerprints"], "fingerprints")["fingerprint"]} for item in parsed]
    core = {"schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_FINGERPRINT_SCHEMA_VERSION, "status": "unavailable", "available": False, "reason": "one_or_more_parent_corrections_invalid", "diagnostic_stage": "incomplete", "non_authoritative": True, "accepted_basis": "hgcer3", "basis_frozen": True, "relative_map_consumed": True, "relative_map_modified": False, "parent_normalization_constructed": False, "correction_constructed": False, "correction_applied_to_production": False, "event_correction_evaluated_for_detached_diagnostic": False, "event_correction_persisted": False, "child_renormalization_performed": False, "downstream_template_application_performed": False, "absolute_probability_constructed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "algorithm_config": _copy(config, "config"), "f3_source_file_sha256": f3_sha256, "f3_map_fingerprint": map_value["fingerprint"], "f3_algorithm_fingerprint": map_value["algorithm_fingerprint"], "f3_artifact_fingerprint": f3_artifact_fingerprint, "f3_runtime_authority": _copy(f3_runtime_authority, "f3_runtime_authority"), "input_fingerprints": inputs, "invalid_parents": _copy(list(invalid), "invalid_parents"), "parents": []}
    fp = {"schema_version": core["schema_version"], "fingerprint_schema_version": core["fingerprint_schema_version"], "accepted_basis": "hgcer3", "algorithm_config": core["algorithm_config"], "f3_source_file_sha256": f3_sha256, "f3_map_fingerprint": core["f3_map_fingerprint"], "f3_algorithm_fingerprint": core["f3_algorithm_fingerprint"], "f3_artifact_fingerprint": core["f3_artifact_fingerprint"], "f3_runtime_authority": core["f3_runtime_authority"], "input_content": [{"setting_id": item["setting_id"], "stable_f1_content_fingerprint": item["stable_f1_content_fingerprint"], "source_file_sha256": item["source_file_sha256"]} for item in inputs], "invalid_parents": core["invalid_parents"], "parents": []}
    core["fingerprint_inputs"] = fp; core["fingerprint"] = _sha256(fp)
    return _copy(core, "unavailable_correction")  # type: ignore[return-value]


def _build_pion_hgcer_method_a_parent_preserving_correction(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, algorithm_config: Mapping[str, object] | None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None, include_review_data: bool) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Shared F.4 calculation; review arrays remain transient at this boundary."""
    config = _resolve_config(algorithm_config); _preflight_raw_phi_degrees(f1_artifacts); parsed = _f3._validate_f1_artifacts(f1_artifacts)
    if set(f1_input_file_hashes) != {item["setting_id"] for item in parsed}:
        raise MethodAParentPreservingCorrectionError("f1_input_file_hashes_settings_invalid")
    hashes = {str(key): _validate_hash(value, "f1_input_file_hash") for key, value in f1_input_file_hashes.items()}
    f3_sha = _validate_hash(f3_input_file_sha256, "f3_input_file_hash")
    raw_rows = _raw_application_rows(f1_artifacts, parsed, float(config["signed_baseline_identity_relative_tolerance"]))
    map_value, models = _validate_f3_artifact(f3_artifact, parsed, hashes, float(config["f3_support_absolute_tolerance"]))
    f3_artifact_fingerprint = _validate_hash(f3_artifact.get("artifact_fingerprint"), "f3_artifact_fingerprint")
    kinematics = {str(_map(item["setting"], "f1_setting").get("kinematic_token")) for item in parsed}
    if len(kinematics) != 1:
        raise MethodAParentPreservingCorrectionError("f1_kinematic_identity_invalid")
    f3_runtime_authority = _f3_runtime_authority(next(iter(kinematics)), f3_sha, map_value, f3_artifact_fingerprint, accepted_f3_runtime_authority_by_kinematic)
    parents: list[dict[str, object]] = []; invalid: list[dict[str, object]] = []; review_data: list[dict[str, object]] = []
    for item in parsed:
        setting_id = str(item["setting_id"]); training_by_t: dict[int, list[Mapping[str, object]]] = {}
        for row in _seq(item["training"], "training"):
            clean = _map(row, "training"); training_by_t.setdefault(int(clean["t_index"]), []).append(clean)
        for t_index, t_low, t_high in _seq(item["canonical_t_geometry"], "geometry"):
            key = (setting_id, int(t_index)); model = models[key]; application = raw_rows.get(key, []); training = training_by_t.get(int(t_index), [])
            identity = {"setting": _copy(item["setting"], "setting"), "setting_id": setting_id, "canonical_t_index": int(t_index), "canonical_t_low": float(t_low), "canonical_t_high": float(t_high)}
            try:
                tree, median, divisor, threshold, support = _support_continuity(training, application, model, float(config["f3_support_absolute_tolerance"]))
                fit = _map(model.get("logistic_fit"), "f3_fit"); coefficients = np.asarray(fit.get("coefficients"), dtype=float)
                if coefficients.shape != (3,): raise MethodAParentPreservingCorrectionError("f3_coefficients_invalid")
                factors: list[float] = []; corrections: list[float] = []; baseline: list[float] = []; shaped: list[float] = []; adjusted: list[float] = []; source: dict[str, list[tuple[float, float]]] = {}; child: dict[tuple[int, float, float], list[tuple[float, float]]] = {}; in_support = 0; ood = 0; zero_w0 = 0
                feature_values = np.asarray([[float(row[name]) for name in _FEATURES] for row in application], dtype=float)
                scaled_values = (feature_values - median) / divisor
                distances = np.asarray(tree.query(scaled_values, k=1)[0], dtype=float)
                inside_mask = distances <= threshold
                raw_factors = np.ones(len(application), dtype=float)
                raw_factors[inside_mask] = np.exp((scaled_values @ coefficients)[inside_mask])
                if not np.all(np.isfinite(raw_factors)) or np.any(raw_factors <= 0.0):
                    raise MethodAParentPreservingCorrectionError("relative_response_invalid")
                in_support = int(np.sum(inside_mask)); ood = int(len(application) - in_support)
                for row, raw_factor in zip(application, raw_factors):
                    b = float(row["signed_baseline_event_contribution"]); factors.append(raw_factor); baseline.append(b); shaped.append(b * raw_factor)
                    if float(row["baseline_pion_weight_w0"]) == 0.0: zero_w0 += 1
                    source.setdefault(str(row["source_label"]), []).append((b, raw_factor))
                    child.setdefault((int(row["phi_index"]), float(row["phi_low"]), float(row["phi_high"])), []).append((b, raw_factor))
                B = math.fsum(baseline); U = math.fsum(shaped)
                if not math.isfinite(B) or B <= 0.0: raise MethodAParentPreservingCorrectionError("baseline_parent_sum_invalid")
                if not math.isfinite(U) or U <= 0.0: raise MethodAParentPreservingCorrectionError("raw_shape_parent_sum_invalid")
                normalization = U / B
                if not math.isfinite(normalization) or normalization <= 0.0: raise MethodAParentPreservingCorrectionError("parent_normalization_invalid")
                for b, factor in zip(baseline, factors):
                    correction = factor / normalization
                    if not math.isfinite(correction) or correction <= 0.0: raise MethodAParentPreservingCorrectionError("correction_invalid")
                    corrections.append(correction); adjusted.append(b * correction)
                adjusted_sum = math.fsum(adjusted); residual = adjusted_sum - B; tolerance = float(config["closure_relative_tolerance"]) * max(1.0, abs(B))
                if abs(residual) > tolerance: raise MethodAParentPreservingCorrectionError("parent_closure_failed")
                source_summary = [{"source_label": label, "event_count": len(rows), "baseline_signed_sum": math.fsum(value[0] for value in rows), "adjusted_signed_sum": math.fsum(value[0] * value[1] / normalization for value in rows), "signed_delta": math.fsum(value[0] * value[1] / normalization for value in rows) - math.fsum(value[0] for value in rows)} for label, rows in sorted(source.items())]
                child_summary = [{"phi_index": child_key[0], "phi_low": child_key[1], "phi_high": child_key[2], "event_count": len(rows), "baseline_signed_sum": math.fsum(value[0] for value in rows), "adjusted_signed_sum": math.fsum(value[0] * value[1] / normalization for value in rows), "signed_delta": math.fsum(value[0] * value[1] / normalization for value in rows) - math.fsum(value[0] for value in rows)} for child_key, rows in sorted(child.items())]
                if not _close(math.fsum(float(row["baseline_signed_sum"]) for row in child_summary), B, float(config["closure_relative_tolerance"])) or not _close(math.fsum(float(row["adjusted_signed_sum"]) for row in child_summary), adjusted_sum, float(config["closure_relative_tolerance"])) or not _close(math.fsum(float(row["signed_delta"]) for row in child_summary), 0.0, float(config["closure_relative_tolerance"])):
                    raise MethodAParentPreservingCorrectionError("child_closure_failed")
                parents.append({**identity, "application_event_count": len(application), "zero_baseline_weight_count": zero_w0, "in_support_count": in_support, "ood_count": ood, "ood_fraction": float(ood / len(application)) if application else 0.0, "baseline_parent_sum": B, "absolute_baseline_parent_sum": math.fsum(abs(value) for value in baseline), "baseline_cancellation_ratio": abs(B) / math.fsum(abs(value) for value in baseline), "raw_shape_parent_sum": U, "parent_normalization": normalization, "adjusted_parent_sum": adjusted_sum, "closure_residual": residual, "closure_relative_scale": tolerance, "closure_passed": True, "raw_shape_factor_summary": _percentiles(factors), "correction_factor_summary": _percentiles(corrections), "ood_final_correction": 1.0 / normalization, "f3_support": support, "source_diagnostics": source_summary, "canonical_phi_diagnostics": child_summary})
                if include_review_data:
                    review_data.append({"setting_id": setting_id, "canonical_t_index": int(t_index), "raw_shape_factors": np.asarray(raw_factors, dtype=float).copy(), "correction_factors": np.asarray(corrections, dtype=float).copy(), "in_support_mask": np.asarray(inside_mask, dtype=bool).copy()})
            except MethodAParentPreservingCorrectionError as exc:
                invalid.append({**identity, "invalid_reason": str(exc)})
    if invalid:
        return _unavailable(parsed, map_value, hashes, f3_sha, f3_artifact_fingerprint, f3_runtime_authority, config, invalid), []
    parents.sort(key=lambda row: (_CANONICAL_SETTINGS.index((str(_map(row["setting"], "setting")["phi_setting"]), str(_map(row["setting"], "setting")["epsilon_filename_token"]))), int(row["canonical_t_index"])))
    inputs = [{"setting": _copy(item["setting"], "setting"), "setting_id": item["setting_id"], "source_file_sha256": hashes[item["setting_id"]], "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "f1_contract_fingerprint": _map(item["fingerprints"], "fingerprints")["fingerprint"]} for item in parsed]
    fingerprint_inputs = {"schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_FINGERPRINT_SCHEMA_VERSION, "accepted_basis": "hgcer3", "algorithm_config": config, "f3_source_file_sha256": f3_sha, "f3_map_fingerprint": map_value["fingerprint"], "f3_algorithm_fingerprint": map_value["algorithm_fingerprint"], "f3_artifact_fingerprint": f3_artifact_fingerprint, "f3_runtime_authority": f3_runtime_authority, "input_content": [{"setting_id": item["setting_id"], "stable_f1_content_fingerprint": item["stable_f1_content_fingerprint"], "source_file_sha256": item["source_file_sha256"]} for item in inputs], "parents": parents}
    core = {"schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_FINGERPRINT_SCHEMA_VERSION, "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete", "non_authoritative": True, "accepted_basis": "hgcer3", "basis_frozen": True, "relative_map_consumed": True, "relative_map_modified": False, "parent_normalization_constructed": True, "correction_constructed": True, "correction_applied_to_production": False, "event_correction_evaluated_for_detached_diagnostic": True, "event_correction_persisted": False, "child_renormalization_performed": False, "downstream_template_application_performed": False, "absolute_probability_constructed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "algorithm_config": config, "f3_source_file_sha256": f3_sha, "f3_map_fingerprint": map_value["fingerprint"], "f3_algorithm_fingerprint": map_value["algorithm_fingerprint"], "f3_artifact_fingerprint": f3_artifact_fingerprint, "f3_runtime_authority": f3_runtime_authority, "input_fingerprints": inputs, "parents": parents}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    review_data.sort(key=lambda value: (_CANONICAL_SETTINGS.index(tuple(str(part) for part in str(value["setting_id"]).split("-", 1))), int(value["canonical_t_index"])))
    return _copy(core, "parent_preserving_correction"), review_data  # type: ignore[return-value]


def build_pion_hgcer_method_a_parent_preserving_correction(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, algorithm_config: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> dict[str, object]:
    """Build aggregate-only F.4 output and discard all transient event factors."""
    correction, _review_data = _build_pion_hgcer_method_a_parent_preserving_correction(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, algorithm_config=algorithm_config, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic, include_review_data=False)
    return correction


def build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, algorithm_config: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Build F.4 once, retaining correction values only for immediate review rendering."""
    return _build_pion_hgcer_method_a_parent_preserving_correction(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, algorithm_config=algorithm_config, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic, include_review_data=True)


def _artifact_wrapper(correction: Mapping[str, object], *, input_paths: Mapping[str, object] | None, generated_at_utc: str | None, git_head: str | None, git_status_short: str | None) -> dict[str, object]:
    provenance = {"generated_at_utc": generated_at_utc, "git_head": git_head, "git_status_short": git_status_short, "input_paths": _copy({} if input_paths is None else input_paths, "input_paths")}
    return {"schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION, "correction": correction, "non_authoritative": True, "relative_map_consumed": True, "relative_map_modified": False, "parent_normalization_constructed": bool(correction["parent_normalization_constructed"]), "correction_constructed": bool(correction["correction_constructed"]), "correction_applied_to_production": False, "event_correction_evaluated_for_detached_diagnostic": bool(correction["event_correction_evaluated_for_detached_diagnostic"]), "event_correction_persisted": False, "child_renormalization_performed": False, "downstream_template_application_performed": False, "absolute_probability_constructed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "basis_frozen": True, "manual_review_required": True, "provenance": provenance, "artifact_fingerprint": _sha256({"schema_version": METHOD_A_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION, "correction_fingerprint": correction["fingerprint"], "input_paths": provenance["input_paths"]})}


def build_pion_hgcer_method_a_parent_preserving_correction_artifact(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, input_paths: Mapping[str, object] | None = None, algorithm_config: Mapping[str, object] | None = None, generated_at_utc: str | None = None, git_head: str | None = None, git_status_short: str | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> dict[str, object]:
    correction = build_pion_hgcer_method_a_parent_preserving_correction(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, algorithm_config=algorithm_config, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    return _artifact_wrapper(correction, input_paths=input_paths, generated_at_utc=generated_at_utc, git_head=git_head, git_status_short=git_status_short)


def build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data(f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, input_paths: Mapping[str, object] | None = None, algorithm_config: Mapping[str, object] | None = None, generated_at_utc: str | None = None, git_head: str | None = None, git_status_short: str | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> tuple[dict[str, object], list[dict[str, object]]]:
    correction, review_data = build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, algorithm_config=algorithm_config, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    return _artifact_wrapper(correction, input_paths=input_paths, generated_at_utc=generated_at_utc, git_head=git_head, git_status_short=git_status_short), review_data


def pion_hgcer_method_a_parent_preserving_correction_filename(kinematic: object) -> str:
    token = str(kinematic)
    if not token or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise MethodAParentPreservingCorrectionError("correction_filename_token_invalid")
    return "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json".format(token)


def write_pion_hgcer_method_a_parent_preserving_correction_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    with open(os.fspath(path), "w", encoding="utf-8") as handle:
        json.dump(_copy(payload, "correction_artifact"), handle, sort_keys=True, indent=2, allow_nan=False); handle.write("\n")
    return os.fspath(path)
