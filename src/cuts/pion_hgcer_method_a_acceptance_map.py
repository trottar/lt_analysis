"""Detached Phase F.3 Method-A ``hgcer3`` relative response map.

This module consumes only serialized, authority-checked F.1 and F.2 artifacts.
It constructs fifteen support-aware *relative* response definitions for review;
it never constructs an absolute probability, correction, weight, normalization,
or production-side object.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.special import expit


METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION = "pion_hgcer_method_a_acceptance_map/v1"
METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_map_fingerprint/v1"
)
METHOD_A_ACCEPTANCE_MAP_ARTIFACT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_map_artifact/v1"
)
_F1_ARTIFACT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract_artifact/v2"
_F1_CONTRACT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract/v2"
_F1_FINGERPRINT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2"
_F2_SCHEMA = "pion_hgcer_method_a_acceptance_representation/v1"
_F2_FINGERPRINT_SCHEMA = "pion_hgcer_method_a_acceptance_representation_fingerprint/v1"
_F2_ARTIFACT_SCHEMA = "pion_hgcer_method_a_acceptance_representation_artifact/v1"
_PRIMARY_FEATURES = (
    "SHMS_delta", "SHMS_xptar", "SHMS_yptar", "P_hgcer_xAtCer", "P_hgcer_yAtCer",
)
ACCEPTED_BASIS = "hgcer3"
ACCEPTED_FEATURES = ("SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer")
_CANONICAL_SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)
_PROMPT_SOURCE = "prompt"
_LOW_NPE_UPPER_BOUND = 2.0
_F1_FEATURE_METADATA = {
    "primary_acceptance_features": list(_PRIMARY_FEATURES),
    "training_population": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
    "training_low_definition": "0_lt_P_hgcer_npeSum_le_2",
    "training_control_definition": "P_hgcer_npeSum_gt_2",
    "application_population": "authoritative_physical_pion_control_P_hgcer_npeSum_gt_2",
    "parent_coordinate": "canonical_t",
    "downstream_yield_coordinates": ["canonical_t", "canonical_phi"],
    "phi_is_training_feature": False,
    "method_b_numerical_dependency": False,
    "probability_map_constructed": False,
    "weight_adjustment_constructed": False,
    "future_normalization_policy": "future_parent_t_only_no_tphi_child_renormalization",
    "absolute_leakage_probability_claimed": False,
}
_F1_PROVENANCE_FINGERPRINT_NAMES = (
    "phase_a_contract_fingerprint", "phase_a_pion_event_population_fingerprint",
    "method_a_fingerprint", "method_a_event_population_fingerprint",
    "part1_config_fingerprint", "coordinate_fingerprint",
)
_F1_CHILD_ASSIGNMENT_FIELDS = (
    "source_label", "entry_index", "t_index", "phi_index", "phi_low", "phi_high", "phi_status",
)

DEFAULT_ALGORITHM_CONFIG = {
    "algorithm_version": "method_a_acceptance_map/v1",
    "minimum_low_count": 25,
    "minimum_control_count": 100,
    "scaling": "full_parent_training_median_iqr_else_std/v1",
    "optimizer": "scipy_lbfgsb_analytic_gradient/v1",
    "regularization_lambda": 1.0e-3,
    "optimizer_maxiter": 1000,
    "optimizer_gtol": 1.0e-8,
    "support_nn_percentile": 99.0,
    "support_ood_fraction_max": 0.10,
    "support_sparse_nonprompt_count": 20,
    "f2_support_continuity_absolute_tolerance": 1.0e-12,
    "review_grid_side": 101,
}
_F2_FROZEN_CONFIG = {
    "algorithm_version": "method_a_acceptance_representation_probe/v1",
    "minimum_low_count": 25,
    "minimum_control_count": 100,
    "fold_count": 5,
    "fold_assignment": "stable_identity_sorted_round_robin_per_class/v1",
    "scaling": "fold_training_median_iqr_else_std/v1",
    "optimizer": "scipy_lbfgsb_analytic_gradient/v1",
    "regularization_lambda": 1.0e-3,
    "optimizer_maxiter": 1000,
    "optimizer_gtol": 1.0e-8,
    "information_median_auc_loss_max": 0.02,
    "information_max_auc_loss_max": 0.05,
    "information_median_balanced_log_loss_penalty_max": 0.01,
    "information_max_balanced_log_loss_penalty_max": 0.02,
    "support_nn_percentile": 99.0,
    "support_ood_fraction_max": 0.10,
    "support_sparse_nonprompt_count": 20,
    "recommendation_logic_version": "minimum_dimension_without_metric_tiebreak/v1",
}


class MethodAAcceptanceMapError(ValueError):
    """Raised when frozen F.1/F.2 authority cannot support Phase F.3."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _json_copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceMapError("{}_not_json_safe".format(label)) from exc


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodAAcceptanceMapError("{}_invalid".format(label))
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodAAcceptanceMapError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodAAcceptanceMapError("{}_nonfinite".format(label))
    try:
        resolved = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceMapError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(resolved):
        raise MethodAAcceptanceMapError("{}_nonfinite".format(label))
    return resolved


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodAAcceptanceMapError("{}_invalid".format(label))
    try:
        resolved = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceMapError("{}_invalid".format(label)) from exc
    if resolved != value:
        raise MethodAAcceptanceMapError("{}_invalid".format(label))
    return resolved


def _string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise MethodAAcceptanceMapError("{}_invalid".format(label))
    return value


def _require(value: Mapping[str, object], name: str, expected: object, label: str) -> None:
    if name not in value or value.get(name) != expected:
        raise MethodAAcceptanceMapError("{}_{}".format(label, name))


def _identity(record: Mapping[str, object], label: str) -> tuple[str, int]:
    source = _string(record.get("source_label"), "{}_source_label".format(label))
    entry = _integer(record.get("entry_index"), "{}_entry_index".format(label))
    if entry < 0:
        raise MethodAAcceptanceMapError("{}_entry_index_invalid".format(label))
    return source, entry


def _setting_id(phi: str, epsilon: str) -> str:
    return "{}-{}".format(phi, epsilon)


def _validate_setting(value: object) -> dict[str, object]:
    setting = _mapping(value, "f1_setting")
    phi = _string(setting.get("phi_setting"), "f1_setting_phi")
    epsilon = _string(setting.get("epsilon_filename_token"), "f1_setting_epsilon").lower()
    if (phi, epsilon) not in _CANONICAL_SETTINGS:
        raise MethodAAcceptanceMapError("f1_setting_not_canonical")
    if _string(setting.get("particle_type"), "f1_setting_particle").lower() != "kaon":
        raise MethodAAcceptanceMapError("f1_setting_particle_invalid")
    kinematic = _string(setting.get("kinematic_token"), "f1_setting_kinematic")
    return {
        "phi_setting": phi, "epsilon_filename_token": epsilon,
        "epsilon_setting": _string(setting.get("epsilon_setting"), "f1_setting_epsilon_setting").lower(),
        "particle_type": "kaon", "kinematic_token": kinematic,
        "Q2": _json_copy(setting.get("Q2"), "f1_setting_Q2"),
        "W": _json_copy(setting.get("W"), "f1_setting_W"),
    }


def _validate_metadata(value: object) -> None:
    if dict(_mapping(value, "f1_feature_metadata")) != _F1_FEATURE_METADATA:
        raise MethodAAcceptanceMapError("f1_feature_metadata_invalid")


def _edge_array(value: object, label: str) -> tuple[float, ...]:
    raw_edges = _sequence(value, label)
    if len(raw_edges) < 2:
        raise MethodAAcceptanceMapError("{}_length_invalid".format(label))
    edges = tuple(_finite(edge, "{}_{}".format(label, index)) for index, edge in enumerate(raw_edges))
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise MethodAAcceptanceMapError("{}_not_strictly_increasing".format(label))
    return edges


def _validate_t_assignment(record: object, t_edges: Sequence[float], label: str) -> None:
    value = _mapping(record, label)
    index = _integer(value.get("t_index"), "{}_t_index".format(label))
    if index < 0 or index >= len(t_edges) - 1:
        raise MethodAAcceptanceMapError("{}_t_index_invalid".format(label))
    if (_finite(value.get("t_low"), "{}_t_low".format(label)), _finite(value.get("t_high"), "{}_t_high".format(label))) != (t_edges[index], t_edges[index + 1]):
        raise MethodAAcceptanceMapError("{}_t_geometry_mismatch".format(label))


def _reconstruct_f1_fingerprints(contract: Mapping[str, object], raw_training: Sequence[object], raw_application: Sequence[object], feature_metadata: Mapping[str, object], label: str) -> dict[str, str]:
    provenance = {name: _string(contract.get(name), "{}_{}".format(label, name)) for name in _F1_PROVENANCE_FINGERPRINT_NAMES}
    host_state = _string(contract.get("host_state"), "{}_host_state".format(label))
    _require(contract, "source_target_state", "post_proton_noRF", label)
    try:
        closure = _mapping(contract.get("method_a_training_summary"), "{}_training_summary".format(label))["by_t_delta"]
    except KeyError as exc:
        raise MethodAAcceptanceMapError("{}_method_a_closure_invalid".format(label)) from exc
    child_projection: list[dict[str, object]] = []
    for index, raw in enumerate(raw_application):
        row = _mapping(raw, "{}_application_{}".format(label, index))
        missing = [name for name in _F1_CHILD_ASSIGNMENT_FIELDS if name not in row]
        if missing:
            raise MethodAAcceptanceMapError("{}_application_child_projection_missing:{}".format(label, ",".join(missing)))
        child_projection.append({name: row[name] for name in _F1_CHILD_ASSIGNMENT_FIELDS})
    values = {
        "method_a_training_population_fingerprint": _sha256(raw_training),
        "application_population_fingerprint": _sha256(raw_application),
        "acceptance_feature_metadata_fingerprint": _sha256(feature_metadata),
        "application_child_assignment_projection_fingerprint": _sha256(child_projection),
    }
    for name, expected in values.items():
        if contract.get(name) != expected:
            raise MethodAAcceptanceMapError("{}_{}_mismatch".format(label, name))
    fingerprint_inputs = {
        "schema_version": _F1_CONTRACT_SCHEMA, "fingerprint_schema_version": _F1_FINGERPRINT_SCHEMA,
        **provenance, "host_state": host_state, "source_target_state": "post_proton_noRF",
        "t_edges": contract.get("t_edges"), "delta_edges": contract.get("delta_edges"), "phi_edges": contract.get("phi_edges"),
        **values, "method_a_closure": closure, "feature_metadata": feature_metadata,
    }
    if contract.get("fingerprint_inputs") != fingerprint_inputs:
        raise MethodAAcceptanceMapError("{}_fingerprint_inputs_mismatch".format(label))
    if contract.get("fingerprint") != _sha256(fingerprint_inputs):
        raise MethodAAcceptanceMapError("{}_fingerprint_mismatch".format(label))
    return {**provenance, **values, "fingerprint": str(contract["fingerprint"]), "host_state": host_state, "source_target_state": "post_proton_noRF"}


def _stable_f1_content_fingerprint(setting: Mapping[str, object], contract: Mapping[str, object], raw_training: Sequence[object], raw_application: Sequence[object]) -> str:
    names = (
        "schema_version", "fingerprint_schema_version", "status", "available", "reason", "diagnostic_stage", "non_authoritative", "production_objects_mutated", "refinement_applied", "production_application_performed", "event_application_performed", "method_b_numerical_dependency", "future_weight_adjustment_constructed", "phase_a_contract_fingerprint", "phase_a_pion_event_population_fingerprint", "method_a_fingerprint", "method_a_event_population_fingerprint", "part1_config_fingerprint", "coordinate_fingerprint", "host_state", "source_target_state", "t_edges", "delta_edges", "phi_edges", "method_a_training_summary", "feature_metadata",
    )
    return _sha256({"setting": setting, "contract_authority": {name: contract.get(name) for name in names}, "method_a_training_records": sorted((_json_copy(row, "f1_training_record") for row in raw_training), key=_canonical_json), "application_records": sorted((_json_copy(row, "f1_application_record") for row in raw_application), key=_canonical_json)})


def _sanitize_training_record(value: object, label: str) -> dict[str, object]:
    row = _mapping(value, label)
    source, entry = _identity(row, label)
    if source != _PROMPT_SOURCE or row.get("nommcuts") is not True:
        raise MethodAAcceptanceMapError("{}_selection_invalid".format(label))
    npe = _finite(row.get("P_hgcer_npeSum"), "{}_npe".format(label))
    response = _string(row.get("response_class"), "{}_response_class".format(label))
    expected = "low" if 0.0 < npe <= _LOW_NPE_UPPER_BOUND else "control" if npe > _LOW_NPE_UPPER_BOUND else None
    if response != expected:
        raise MethodAAcceptanceMapError("{}_response_classification_invalid".format(label))
    result: dict[str, object] = {"source_label": source, "entry_index": entry, "response_class": response, "P_hgcer_npeSum": npe, "t_index": _integer(row.get("t_index"), "{}_t_index".format(label)), "t_low": _finite(row.get("t_low"), "{}_t_low".format(label)), "t_high": _finite(row.get("t_high"), "{}_t_high".format(label))}
    if result["t_low"] >= result["t_high"]:
        raise MethodAAcceptanceMapError("{}_t_geometry_invalid".format(label))
    for feature in _PRIMARY_FEATURES:
        result[feature] = _finite(row.get(feature), "{}_{}".format(label, feature))
    return result


def _sanitize_application_record(value: object, label: str) -> dict[str, object]:
    row = _mapping(value, label)
    source, entry = _identity(row, label)
    npe = _finite(row.get("P_hgcer_npeSum"), "{}_npe".format(label))
    if npe <= _LOW_NPE_UPPER_BOUND:
        raise MethodAAcceptanceMapError("{}_npe_not_physical_control".format(label))
    result: dict[str, object] = {"source_label": source, "entry_index": entry, "P_hgcer_npeSum": npe, "t_index": _integer(row.get("t_index"), "{}_t_index".format(label)), "t_low": _finite(row.get("t_low"), "{}_t_low".format(label)), "t_high": _finite(row.get("t_high"), "{}_t_high".format(label))}
    if result["t_low"] >= result["t_high"]:
        raise MethodAAcceptanceMapError("{}_t_geometry_invalid".format(label))
    for feature in _PRIMARY_FEATURES:
        result[feature] = _finite(row.get(feature), "{}_{}".format(label, feature))
    return result


def _validate_f1_artifacts(f1_artifacts: object) -> list[dict[str, object]]:
    raw_artifacts = _sequence(f1_artifacts, "f1_artifacts")
    if len(raw_artifacts) != len(_CANONICAL_SETTINGS):
        raise MethodAAcceptanceMapError("f1_setting_count_invalid")
    parsed: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    kinematics: set[str] = set()
    geometry_reference: tuple[tuple[int, float, float], ...] | None = None
    for artifact_index, raw_artifact in enumerate(raw_artifacts):
        artifact = _mapping(raw_artifact, "f1_artifact_{}".format(artifact_index))
        _require(artifact, "schema_version", _F1_ARTIFACT_SCHEMA, "f1_artifact")
        for name in ("non_authoritative", "production_objects_mutated", "refinement_applied", "production_application_performed", "event_application_performed"):
            _require(artifact, name, True if name == "non_authoritative" else False, "f1_artifact")
        setting = _validate_setting(artifact.get("setting"))
        key = (str(setting["phi_setting"]), str(setting["epsilon_filename_token"]))
        if key in seen:
            raise MethodAAcceptanceMapError("f1_setting_duplicate")
        seen.add(key); kinematics.add(str(setting["kinematic_token"]))
        contract = _mapping(artifact.get("contract"), "f1_contract")
        _require(contract, "schema_version", _F1_CONTRACT_SCHEMA, "f1_contract")
        _require(contract, "fingerprint_schema_version", _F1_FINGERPRINT_SCHEMA, "f1_contract")
        for name, expected in (("status", "available"), ("available", True), ("diagnostic_stage", "complete"), ("reason", None)):
            _require(contract, name, expected, "f1_contract")
        for name in ("non_authoritative", "method_b_numerical_dependency", "event_application_performed", "production_application_performed", "production_objects_mutated", "future_weight_adjustment_constructed", "refinement_applied"):
            _require(contract, name, True if name == "non_authoritative" else False, "f1_contract")
        metadata = _mapping(contract.get("feature_metadata"), "f1_feature_metadata")
        _validate_metadata(metadata)
        t_edges = _edge_array(contract.get("t_edges"), "f1_t_edges")
        _edge_array(contract.get("delta_edges"), "f1_delta_edges"); _edge_array(contract.get("phi_edges"), "f1_phi_edges")
        raw_training = _sequence(contract.get("method_a_training_records"), "f1_training_records")
        raw_application = _sequence(contract.get("application_records"), "f1_application_records")
        if not raw_training or not raw_application:
            raise MethodAAcceptanceMapError("f1_population_empty")
        for index, row in enumerate(raw_training): _validate_t_assignment(row, t_edges, "f1_training_{}_{}".format(artifact_index, index))
        for index, row in enumerate(raw_application): _validate_t_assignment(row, t_edges, "f1_application_{}_{}".format(artifact_index, index))
        fingerprints = _reconstruct_f1_fingerprints(contract, raw_training, raw_application, metadata, "f1_contract_{}".format(artifact_index))
        training = [_sanitize_training_record(row, "f1_training_{}_{}".format(artifact_index, index)) for index, row in enumerate(raw_training)]
        application = [_sanitize_application_record(row, "f1_application_{}_{}".format(artifact_index, index)) for index, row in enumerate(raw_application)]
        training.sort(key=lambda row: (str(row["source_label"]), int(row["entry_index"])))
        application.sort(key=lambda row: (str(row["source_label"]), int(row["entry_index"])))
        training_ids = [(str(row["source_label"]), int(row["entry_index"])) for row in training]
        application_ids = [(str(row["source_label"]), int(row["entry_index"])) for row in application]
        if len(training_ids) != len(set(training_ids)) or len(application_ids) != len(set(application_ids)):
            raise MethodAAcceptanceMapError("f1_record_identity_duplicate")
        controls = {(str(row["source_label"]), int(row["entry_index"])) for row in training if row["response_class"] == "control"}
        if any((str(row["source_label"]), int(row["entry_index"])) not in controls for row in application if row["source_label"] == _PROMPT_SOURCE):
            raise MethodAAcceptanceMapError("f1_prompt_application_identity_missing_training_control")
        geometry = tuple((index, t_edges[index], t_edges[index + 1]) for index in range(len(t_edges) - 1))
        if geometry_reference is None: geometry_reference = geometry
        elif geometry_reference != geometry: raise MethodAAcceptanceMapError("f1_canonical_t_geometry_mixed")
        parsed.append({"setting": setting, "setting_id": _setting_id(*key), "training": training, "application": application, "fingerprints": fingerprints, "stable_content_fingerprint": _stable_f1_content_fingerprint(setting, contract, raw_training, raw_application), "canonical_t_geometry": geometry})
    if seen != set(_CANONICAL_SETTINGS): raise MethodAAcceptanceMapError("f1_canonical_settings_missing_or_unexpected")
    if len(kinematics) != 1: raise MethodAAcceptanceMapError("f1_kinematic_mixed")
    parsed.sort(key=lambda item: _CANONICAL_SETTINGS.index((str(_mapping(item["setting"], "setting")["phi_setting"]), str(_mapping(item["setting"], "setting")["epsilon_filename_token"]))))
    return parsed


def _resolve_config(algorithm_config: Mapping[str, object] | None) -> dict[str, object]:
    config = dict(DEFAULT_ALGORITHM_CONFIG)
    if algorithm_config is not None:
        supplied = _mapping(algorithm_config, "algorithm_config")
        unexpected = sorted(set(supplied) - set(config))
        if unexpected:
            raise MethodAAcceptanceMapError("algorithm_config_unknown:{}".format(",".join(unexpected)))
        config.update(supplied)
    for key in ("minimum_low_count", "minimum_control_count", "optimizer_maxiter", "support_sparse_nonprompt_count", "review_grid_side"):
        if _integer(config[key], "algorithm_config_{}".format(key)) <= 0:
            raise MethodAAcceptanceMapError("algorithm_config_{}_invalid".format(key))
    for key in ("regularization_lambda", "optimizer_gtol", "support_nn_percentile", "support_ood_fraction_max", "f2_support_continuity_absolute_tolerance"):
        if _finite(config[key], "algorithm_config_{}".format(key)) < 0.0:
            raise MethodAAcceptanceMapError("algorithm_config_{}_invalid".format(key))
    if config["optimizer"] != "scipy_lbfgsb_analytic_gradient/v1":
        raise MethodAAcceptanceMapError("algorithm_config_optimizer_invalid")
    if config["scaling"] != "full_parent_training_median_iqr_else_std/v1":
        raise MethodAAcceptanceMapError("algorithm_config_scaling_invalid")
    if not 0.0 < float(config["support_nn_percentile"]) <= 100.0:
        raise MethodAAcceptanceMapError("algorithm_config_support_nn_percentile_invalid")
    return _json_copy(config, "algorithm_config")  # type: ignore[return-value]


def _feature_array(records: Sequence[Mapping[str, object]], features: Sequence[str]) -> np.ndarray:
    return np.asarray([[float(row[name]) for name in features] for row in records], dtype=float)


def _robust_scale(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    median = np.asarray(np.median(values, axis=0), dtype=float)
    scale = np.asarray(np.percentile(values, 75.0, axis=0) - np.percentile(values, 25.0, axis=0), dtype=float)
    invalid = ~np.isfinite(scale) | (scale <= 0.0)
    if np.any(invalid):
        standard_deviation = np.asarray(np.std(values, axis=0), dtype=float)
        scale[invalid] = standard_deviation[invalid]
    if not np.all(np.isfinite(median)) or not np.all(np.isfinite(scale)) or np.any(scale <= 0.0):
        return None
    return median, scale


def _logistic_objective(parameters: np.ndarray, values: np.ndarray, labels: np.ndarray, weights: np.ndarray, regularization_lambda: float) -> tuple[float, np.ndarray]:
    intercept, coefficients = parameters[0], parameters[1:]
    logits = intercept + values @ coefficients
    losses = np.logaddexp(0.0, logits) - labels * logits
    objective = float(np.mean(weights * losses) + 0.5 * regularization_lambda * np.dot(coefficients, coefficients))
    residual = weights * (expit(logits) - labels)
    gradient = np.empty_like(parameters)
    gradient[0] = float(np.mean(residual))
    gradient[1:] = values.T @ residual / values.shape[0] + regularization_lambda * coefficients
    return objective, gradient


def _percentiles(values: np.ndarray) -> dict[str, float | None]:
    if values.size == 0:
        return {"p50": None, "p95": None, "p99": None, "max": None}
    return {"p50": float(np.percentile(values, 50.0)), "p95": float(np.percentile(values, 95.0)), "p99": float(np.percentile(values, 99.0)), "max": float(np.max(values))}


def _support_summary(training: Sequence[Mapping[str, object]], application: Sequence[Mapping[str, object]], config: Mapping[str, object]) -> tuple[dict[str, object], cKDTree | None, np.ndarray | None, np.ndarray | None]:
    nonprompt = [row for row in application if row["source_label"] != _PROMPT_SOURCE]
    sparse = len(nonprompt) < int(config["support_sparse_nonprompt_count"])
    values = _feature_array(training, ACCEPTED_FEATURES)
    scale = _robust_scale(values)
    if scale is None:
        return ({"valid": False, "invalid_reason": "training_scaling_invalid", "nonprompt_application_count": len(nonprompt), "training_nn_p50": None, "training_nn_p95": None, "training_nn_p99": None, "application_nn_p50": None, "application_nn_p95": None, "application_nn_p99": None, "application_nn_max": None, "application_ood_count": None, "application_ood_fraction": None, "statistically_sparse": sparse, "support_gate_passed": False}, None, None, None)
    median, divisor = scale
    scaled_training = (values - median) / divisor
    if scaled_training.shape[0] < 2:
        return ({"valid": False, "invalid_reason": "training_neighbor_reference_invalid", "nonprompt_application_count": len(nonprompt), "training_nn_p50": None, "training_nn_p95": None, "training_nn_p99": None, "application_nn_p50": None, "application_nn_p95": None, "application_nn_p99": None, "application_nn_max": None, "application_ood_count": None, "application_ood_fraction": None, "statistically_sparse": sparse, "support_gate_passed": False}, None, None, None)
    tree = cKDTree(scaled_training)
    training_distances = np.asarray(tree.query(scaled_training, k=2)[0][:, 1], dtype=float)
    threshold = float(np.percentile(training_distances, float(config["support_nn_percentile"])))
    app_distances = np.asarray(tree.query((_feature_array(nonprompt, ACCEPTED_FEATURES) - median) / divisor, k=1)[0], dtype=float) if nonprompt else np.asarray([], dtype=float)
    train = _percentiles(training_distances)
    applied = _percentiles(app_distances)
    ood_count = int(np.sum(app_distances > threshold))
    ood_fraction = float(ood_count / app_distances.size) if app_distances.size else 0.0
    return ({"valid": True, "invalid_reason": None, "nonprompt_application_count": int(app_distances.size), "training_nn_p50": train["p50"], "training_nn_p95": train["p95"], "training_nn_p99": train["p99"], "application_nn_p50": applied["p50"], "application_nn_p95": applied["p95"], "application_nn_p99": applied["p99"], "application_nn_max": applied["max"], "application_ood_count": ood_count, "application_ood_fraction": ood_fraction, "statistically_sparse": sparse, "support_gate_passed": bool(ood_fraction <= float(config["support_ood_fraction_max"])), "support_distance_threshold": threshold}, tree, median, divisor)


def _normalize_f1_hashes(parsed: Sequence[Mapping[str, object]], values: Mapping[str, object] | None) -> list[dict[str, object]]:
    supplied = {} if values is None else _mapping(values, "f1_input_file_hashes")
    expected = {str(item["setting_id"]) for item in parsed}
    if set(supplied) != expected:
        raise MethodAAcceptanceMapError("f1_input_file_hashes_settings_invalid")
    result = []
    for item in parsed:
        sha = supplied.get(str(item["setting_id"]))
        if not isinstance(sha, str) or len(sha) != 64 or any(char not in "0123456789abcdef" for char in sha.lower()):
            raise MethodAAcceptanceMapError("f1_input_file_hash_invalid")
        result.append({"setting_id": str(item["setting_id"]), "sha256": sha.lower()})
    return result


def _f2_candidate_definition() -> dict[str, object]:
    return {"candidate_id": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "dimension": 3, "automatic_recommendation_eligible": True, "role": "reduced_detector_coordinate"}


def _validate_f2_artifact(value: object, parsed_f1: Sequence[Mapping[str, object]], f1_hashes: Sequence[Mapping[str, object]]) -> tuple[Mapping[str, object], dict[tuple[str, int], Mapping[str, object]]]:
    artifact = _mapping(value, "f2_artifact")
    _require(artifact, "schema_version", _F2_ARTIFACT_SCHEMA, "f2_artifact")
    for name, expected in (("non_authoritative", True), ("probe_only", True), ("probability_map_constructed", False), ("event_probability_persisted", False), ("event_application_performed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False), ("basis_frozen", False), ("manual_review_required", True)):
        _require(artifact, name, expected, "f2_artifact")
    provenance = _mapping(artifact.get("provenance"), "f2_artifact_provenance")
    expected_artifact_fingerprint = _sha256({"schema_version": _F2_ARTIFACT_SCHEMA, "representation_fingerprint": artifact.get("representation", {}).get("fingerprint") if isinstance(artifact.get("representation"), Mapping) else None, "input_paths": provenance.get("input_paths")})
    if artifact.get("artifact_fingerprint") != expected_artifact_fingerprint:
        raise MethodAAcceptanceMapError("f2_artifact_fingerprint_mismatch")
    representation = _mapping(artifact.get("representation"), "f2_representation")
    _require(representation, "schema_version", _F2_SCHEMA, "f2_representation")
    _require(representation, "fingerprint_schema_version", _F2_FINGERPRINT_SCHEMA, "f2_representation")
    for name, expected in (("status", "available"), ("available", True), ("diagnostic_stage", "complete"), ("non_authoritative", True), ("probe_only", True), ("probability_map_constructed", False), ("event_probability_persisted", False), ("event_application_performed", False), ("production_application_performed", False), ("production_objects_mutated", False), ("method_b_numerical_dependency", False), ("basis_frozen", False), ("manual_review_required", True)):
        _require(representation, name, expected, "f2_representation")
    if representation.get("primary_acceptance_features") != list(_PRIMARY_FEATURES):
        raise MethodAAcceptanceMapError("f2_primary_features_invalid")
    if dict(_mapping(representation.get("algorithm_config"), "f2_algorithm_config")) != _F2_FROZEN_CONFIG:
        raise MethodAAcceptanceMapError("f2_algorithm_config_invalid")
    fingerprint_inputs = _mapping(representation.get("fingerprint_inputs"), "f2_fingerprint_inputs")
    if representation.get("fingerprint") != _sha256(fingerprint_inputs):
        raise MethodAAcceptanceMapError("f2_representation_fingerprint_mismatch")
    for name in ("candidate_definitions", "algorithm_config", "algorithm_fingerprint", "response_support", "groups", "candidate_summaries", "recommendation"):
        if fingerprint_inputs.get(name) != representation.get(name):
            raise MethodAAcceptanceMapError("f2_fingerprint_content_mismatch:{}".format(name))
    algorithm_fingerprint = _sha256({"representation_schema_version": _F2_SCHEMA, "fingerprint_schema_version": _F2_FINGERPRINT_SCHEMA, "candidate_definitions": representation.get("candidate_definitions"), "algorithm_config": representation.get("algorithm_config")})
    if representation.get("algorithm_fingerprint") != algorithm_fingerprint:
        raise MethodAAcceptanceMapError("f2_algorithm_fingerprint_mismatch")
    definitions = _sequence(representation.get("candidate_definitions"), "f2_candidate_definitions")
    if _f2_candidate_definition() not in definitions:
        raise MethodAAcceptanceMapError("f2_hgcer3_definition_invalid")
    recommendation = _mapping(representation.get("recommendation"), "f2_recommendation")
    for name, expected in (("recommendation_status", "unique_supported_reduced_basis"), ("recommended_basis", ACCEPTED_BASIS), ("basis_frozen", False), ("manual_review_required", True)):
        _require(recommendation, name, expected, "f2_recommendation")
    summaries = _sequence(representation.get("candidate_summaries"), "f2_candidate_summaries")
    hgcer_summary = next((item for item in summaries if isinstance(item, Mapping) and item.get("candidate_id") == ACCEPTED_BASIS), None)
    if not isinstance(hgcer_summary, Mapping) or hgcer_summary.get("ordered_features") != list(ACCEPTED_FEATURES) or hgcer_summary.get("information_gate_passed") is not True or hgcer_summary.get("application_support_gate_passed") is not True or hgcer_summary.get("overall_candidate_passed") is not True:
        raise MethodAAcceptanceMapError("f2_hgcer3_summary_invalid")
    expected_f1 = {str(item["setting_id"]): item for item in parsed_f1}
    hashes = {str(item["setting_id"]): str(item["sha256"]) for item in f1_hashes}
    inputs = _sequence(representation.get("input_fingerprints"), "f2_input_fingerprints")
    if len(inputs) != len(expected_f1):
        raise MethodAAcceptanceMapError("f2_input_fingerprint_count_invalid")
    for input_item in inputs:
        item = _mapping(input_item, "f2_input_fingerprint")
        setting_id = _string(item.get("setting_id"), "f2_input_setting_id")
        source = expected_f1.get(setting_id)
        if source is None or item.get("setting") != source["setting"] or item.get("source_file_sha256") != hashes[setting_id] or item.get("stable_f1_content_fingerprint") != source["stable_content_fingerprint"]:
            raise MethodAAcceptanceMapError("f2_input_fingerprint_stale")
        fingerprints = _mapping(source["fingerprints"], "f1_fingerprints")
        for name, expected in fingerprints.items():
            if item.get(name) != expected:
                raise MethodAAcceptanceMapError("f2_input_contract_fingerprint_stale")
    if fingerprint_inputs.get("input_content") != [{"setting_id": str(item["setting_id"]), "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "source_file_sha256": hashes[str(item["setting_id"])]} for item in parsed_f1]:
        raise MethodAAcceptanceMapError("f2_input_content_invalid")
    groups = _sequence(representation.get("groups"), "f2_groups")
    if len(groups) != 15:
        raise MethodAAcceptanceMapError("f2_group_count_invalid")
    indexed: dict[tuple[str, int], Mapping[str, object]] = {}
    for raw_group in groups:
        group = _mapping(raw_group, "f2_group")
        setting_id = _string(group.get("setting_id"), "f2_group_setting_id")
        index = _integer(group.get("canonical_t_index"), "f2_group_t_index")
        if (setting_id, index) in indexed:
            raise MethodAAcceptanceMapError("f2_group_duplicate")
        source = expected_f1.get(setting_id)
        if source is None or index < 0 or index >= len(_sequence(source["canonical_t_geometry"], "f1_geometry")):
            raise MethodAAcceptanceMapError("f2_group_identity_invalid")
        geometry = _sequence(source["canonical_t_geometry"], "f1_geometry")[index]
        training = [row for row in _sequence(source["training"], "f1_training") if _mapping(row, "f1_training")["t_index"] == index]
        application = [row for row in _sequence(source["application"], "f1_application") if _mapping(row, "f1_application")["t_index"] == index]
        if group.get("setting") != source["setting"] or group.get("canonical_t_low") != geometry[1] or group.get("canonical_t_high") != geometry[2] or group.get("training_positive_count") != len(training) or group.get("low_count") != sum(row["response_class"] == "low" for row in training) or group.get("control_count") != sum(row["response_class"] == "control" for row in training) or group.get("nonprompt_application_count") != sum(row["source_label"] != _PROMPT_SOURCE for row in application):
            raise MethodAAcceptanceMapError("f2_group_geometry_invalid")
        metrics = _mapping(group.get("candidate_metrics"), "f2_candidate_metrics")
        candidate = _mapping(metrics.get(ACCEPTED_BASIS), "f2_hgcer3_metrics")
        support = _mapping(candidate.get("application_support"), "f2_hgcer3_support")
        if candidate.get("low_count") != group.get("low_count") or candidate.get("control_count") != group.get("control_count") or candidate.get("valid") is not True or support.get("valid") is not True or support.get("statistically_sparse") is not False or support.get("support_gate_passed") is not True:
            raise MethodAAcceptanceMapError("f2_hgcer3_group_invalid")
        for name in ("nonprompt_application_count", "application_ood_count"):
            _integer(support.get(name), "f2_support_{}".format(name))
        for name in ("training_nn_p50", "training_nn_p95", "training_nn_p99", "application_nn_p50", "application_nn_p95", "application_nn_p99", "application_nn_max", "application_ood_fraction"):
            _finite(support.get(name), "f2_support_{}".format(name))
        indexed[(setting_id, index)] = group
    if set(indexed) != {(str(item["setting_id"]), index) for item in parsed_f1 for index, _, _ in _sequence(item["canonical_t_geometry"], "f1_geometry")}:
        raise MethodAAcceptanceMapError("f2_group_inventory_invalid")
    return representation, indexed


def _support_continuity(f2_group: Mapping[str, object], support: Mapping[str, object], tolerance: float) -> dict[str, object]:
    metrics = _mapping(f2_group.get("candidate_metrics"), "f2_candidate_metrics")
    f2_support = _mapping(_mapping(metrics.get(ACCEPTED_BASIS), "f2_hgcer3_metrics").get("application_support"), "f2_hgcer3_support")
    fields = ("nonprompt_application_count", "training_nn_p50", "training_nn_p95", "training_nn_p99", "application_nn_p50", "application_nn_p95", "application_nn_p99", "application_nn_max", "application_ood_count", "application_ood_fraction", "statistically_sparse", "support_gate_passed")
    mismatches: list[str] = []
    for name in fields:
        left, right = support.get(name), f2_support.get(name)
        if isinstance(left, float) or isinstance(right, float):
            if left is None or right is None or abs(float(left) - float(right)) > tolerance:
                mismatches.append(name)
        elif left != right:
            mismatches.append(name)
    return {"passed": not mismatches, "absolute_tolerance": tolerance, "compared_fields": list(fields), "mismatched_fields": mismatches, "f2_training_nn_p99": f2_support.get("training_nn_p99")}


def _unavailable_result(parsed: Sequence[Mapping[str, object]], config: Mapping[str, object], f2_representation: Mapping[str, object], f1_hashes: Sequence[Mapping[str, object]], f2_file_sha256: str, reasons: Sequence[Mapping[str, object]]) -> dict[str, object]:
    inputs = [{"setting": _json_copy(item["setting"], "input_setting"), "setting_id": item["setting_id"], "source_file_sha256": next(entry["sha256"] for entry in f1_hashes if entry["setting_id"] == item["setting_id"]), "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "f1_contract_fingerprint": _mapping(item["fingerprints"], "f1_fingerprints")["fingerprint"]} for item in parsed]
    core = {"schema_version": METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION, "status": "unavailable", "available": False, "reason": "one_or_more_parent_models_invalid", "diagnostic_stage": "incomplete", "non_authoritative": True, "map_constructed": False, "map_applied": False, "absolute_probability_constructed": False, "parent_normalization_constructed": False, "correction_constructed": False, "event_probability_persisted": False, "event_application_performed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "accepted_basis": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "basis_frozen": True, "manual_review_required": True, "algorithm_config": _json_copy(config, "algorithm_config"), "f2_representation_fingerprint": f2_representation["fingerprint"], "f2_algorithm_fingerprint": f2_representation["algorithm_fingerprint"], "f2_source_file_sha256": f2_file_sha256, "input_fingerprints": inputs, "invalid_parents": _json_copy(list(reasons), "invalid_parents"), "models": []}
    fingerprint_inputs = {"map_schema_version": METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION, "accepted_basis": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "algorithm_config": core["algorithm_config"], "f2_representation_fingerprint": core["f2_representation_fingerprint"], "f2_algorithm_fingerprint": core["f2_algorithm_fingerprint"], "f2_source_file_sha256": f2_file_sha256, "input_content": [{"setting_id": item["setting_id"], "stable_f1_content_fingerprint": item["stable_f1_content_fingerprint"], "source_file_sha256": item["source_file_sha256"]} for item in inputs], "invalid_parents": core["invalid_parents"], "models": []}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    return _json_copy(core, "unavailable_map")  # type: ignore[return-value]


def build_pion_hgcer_method_a_acceptance_map(f1_artifacts: Sequence[Mapping[str, object]], f2_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f2_input_file_sha256: str, algorithm_config: Mapping[str, object] | None = None) -> dict[str, object]:
    """Build the aggregate-only, support-aware Phase F.3 relative map."""
    parsed = _validate_f1_artifacts(f1_artifacts)
    config = _resolve_config(algorithm_config)
    f1_hashes = _normalize_f1_hashes(parsed, f1_input_file_hashes)
    if not isinstance(f2_input_file_sha256, str) or len(f2_input_file_sha256) != 64 or any(char not in "0123456789abcdef" for char in f2_input_file_sha256.lower()):
        raise MethodAAcceptanceMapError("f2_input_file_hash_invalid")
    representation, f2_groups = _validate_f2_artifact(f2_artifact, parsed, f1_hashes)
    models: list[dict[str, object]] = []
    invalid: list[dict[str, object]] = []
    for item in parsed:
        by_t_training: dict[int, list[Mapping[str, object]]] = {}
        by_t_application: dict[int, list[Mapping[str, object]]] = {}
        for row in _sequence(item["training"], "training"):
            parsed_row = _mapping(row, "training"); by_t_training.setdefault(int(parsed_row["t_index"]), []).append(parsed_row)
        for row in _sequence(item["application"], "application"):
            parsed_row = _mapping(row, "application"); by_t_application.setdefault(int(parsed_row["t_index"]), []).append(parsed_row)
        for t_index, t_low, t_high in _sequence(item["canonical_t_geometry"], "geometry"):
            training, application = by_t_training.get(int(t_index), []), by_t_application.get(int(t_index), [])
            low_count = sum(row["response_class"] == "low" for row in training)
            control_count = sum(row["response_class"] == "control" for row in training)
            identity = {"setting_id": item["setting_id"], "canonical_t_index": int(t_index), "canonical_t_low": float(t_low), "canonical_t_high": float(t_high), "low_count": int(low_count), "control_count": int(control_count)}
            if low_count < int(config["minimum_low_count"]) or control_count < int(config["minimum_control_count"]):
                invalid.append({**identity, "invalid_reason": "insufficient_response_support"}); continue
            support, _tree, median, divisor = _support_summary(training, application, config)
            if not bool(support["valid"]):
                invalid.append({**identity, "invalid_reason": support["invalid_reason"]}); continue
            continuity = _support_continuity(f2_groups[(str(item["setting_id"]), int(t_index))], support, float(config["f2_support_continuity_absolute_tolerance"]))
            if not continuity["passed"]:
                invalid.append({**identity, "invalid_reason": "f2_application_support_continuity_failed", "support_continuity": continuity}); continue
            if bool(support["statistically_sparse"]) or not bool(support["support_gate_passed"]):
                invalid.append({**identity, "invalid_reason": "application_support_gate_failed", "support": support}); continue
            assert median is not None and divisor is not None
            values = _feature_array(training, ACCEPTED_FEATURES)
            scaled = (values - median) / divisor
            labels = np.asarray([1.0 if row["response_class"] == "low" else 0.0 for row in training], dtype=float)
            weights = np.where(labels == 1.0, labels.size / (2.0 * low_count), labels.size / (2.0 * control_count))
            result = minimize(_logistic_objective, np.zeros(scaled.shape[1] + 1, dtype=float), args=(scaled, labels, weights, float(config["regularization_lambda"])), method="L-BFGS-B", jac=True, options={"maxiter": int(config["optimizer_maxiter"]), "gtol": float(config["optimizer_gtol"])})
            if not result.success or not np.all(np.isfinite(result.x)) or not math.isfinite(float(result.fun)):
                invalid.append({**identity, "invalid_reason": "optimizer_not_converged"}); continue
            models.append({"setting": _json_copy(item["setting"], "model_setting"), **identity, "training_positive_count": len(training), "nonprompt_application_count": support["nonprompt_application_count"], "ordered_features": list(ACCEPTED_FEATURES), "scaler": {"method": str(config["scaling"]), "median": [float(value) for value in median], "divisor": [float(value) for value in divisor]}, "logistic_fit": {"optimizer": str(config["optimizer"]), "regularization_lambda": float(config["regularization_lambda"]), "optimizer_maxiter": int(config["optimizer_maxiter"]), "optimizer_gtol": float(config["optimizer_gtol"]), "converged": True, "objective": float(result.fun), "intercept_provenance_only": float(result.x[0]), "coefficients": [float(value) for value in result.x[1:]]}, "relative_response_definition": "exp(coefficients_dot_robust_scaled_hgcer3_without_intercept)", "relative_response_at_training_median": 1.0, "training_feature_minimum": [float(value) for value in np.min(values, axis=0)], "training_feature_maximum": [float(value) for value in np.max(values, axis=0)], "application_support": support, "f2_support_continuity": continuity})
    if invalid:
        return _unavailable_result(parsed, config, representation, f1_hashes, f2_input_file_sha256.lower(), invalid)
    models.sort(key=lambda model: (_CANONICAL_SETTINGS.index((str(_mapping(model["setting"], "model_setting")["phi_setting"]), str(_mapping(model["setting"], "model_setting")["epsilon_filename_token"]))), int(model["canonical_t_index"])))
    inputs = [{"setting": _json_copy(item["setting"], "input_setting"), "setting_id": item["setting_id"], "source_file_sha256": next(entry["sha256"] for entry in f1_hashes if entry["setting_id"] == item["setting_id"]), "stable_f1_content_fingerprint": item["stable_content_fingerprint"], "f1_contract_fingerprint": _mapping(item["fingerprints"], "f1_fingerprints")["fingerprint"]} for item in parsed]
    algorithm_fingerprint = _sha256({"map_schema_version": METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION, "accepted_basis": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "algorithm_config": config})
    fingerprint_inputs = {"map_schema_version": METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION, "accepted_basis": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "algorithm_config": config, "algorithm_fingerprint": algorithm_fingerprint, "f2_representation_fingerprint": representation["fingerprint"], "f2_algorithm_fingerprint": representation["algorithm_fingerprint"], "f2_source_file_sha256": f2_input_file_sha256.lower(), "input_content": [{"setting_id": item["setting_id"], "stable_f1_content_fingerprint": item["stable_f1_content_fingerprint"], "source_file_sha256": item["source_file_sha256"]} for item in inputs], "models": models}
    core = {"schema_version": METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION, "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete", "non_authoritative": True, "map_constructed": True, "map_applied": False, "absolute_probability_constructed": False, "parent_normalization_constructed": False, "correction_constructed": False, "event_probability_persisted": False, "event_application_performed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "accepted_basis": ACCEPTED_BASIS, "ordered_features": list(ACCEPTED_FEATURES), "basis_frozen": True, "manual_review_required": True, "algorithm_config": config, "algorithm_fingerprint": algorithm_fingerprint, "f2_representation_fingerprint": representation["fingerprint"], "f2_algorithm_fingerprint": representation["algorithm_fingerprint"], "f2_source_file_sha256": f2_input_file_sha256.lower(), "input_fingerprints": inputs, "models": models}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    return _json_copy(core, "acceptance_map")  # type: ignore[return-value]


def build_pion_hgcer_method_a_acceptance_map_artifact(f1_artifacts: Sequence[Mapping[str, object]], f2_artifact: Mapping[str, object], *, f1_input_file_hashes: Mapping[str, object], f2_input_file_sha256: str, input_paths: Mapping[str, object] | None = None, algorithm_config: Mapping[str, object] | None = None, generated_at_utc: str | None = None, git_head: str | None = None, git_status_short: str | None = None) -> dict[str, object]:
    """Wrap the detached map with non-scientific command-line provenance."""
    acceptance_map = build_pion_hgcer_method_a_acceptance_map(f1_artifacts, f2_artifact, f1_input_file_hashes=f1_input_file_hashes, f2_input_file_sha256=f2_input_file_sha256, algorithm_config=algorithm_config)
    provenance = {"generated_at_utc": generated_at_utc, "git_head": git_head, "git_status_short": git_status_short, "input_paths": _json_copy({} if input_paths is None else input_paths, "input_paths")}
    return {"schema_version": METHOD_A_ACCEPTANCE_MAP_ARTIFACT_SCHEMA_VERSION, "acceptance_map": acceptance_map, "non_authoritative": True, "map_constructed": bool(acceptance_map["map_constructed"]), "map_applied": False, "absolute_probability_constructed": False, "parent_normalization_constructed": False, "correction_constructed": False, "event_probability_persisted": False, "event_application_performed": False, "production_application_performed": False, "production_objects_mutated": False, "method_b_numerical_dependency": False, "basis_frozen": True, "manual_review_required": True, "provenance": provenance, "artifact_fingerprint": _sha256({"schema_version": METHOD_A_ACCEPTANCE_MAP_ARTIFACT_SCHEMA_VERSION, "map_fingerprint": acceptance_map["fingerprint"], "input_paths": provenance["input_paths"]})}


def _filename_token(value: object) -> str:
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise MethodAAcceptanceMapError("method_a_acceptance_map_filename_token_invalid")
    token = str(value)
    if not token or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise MethodAAcceptanceMapError("method_a_acceptance_map_filename_token_invalid")
    return token


def pion_hgcer_method_a_acceptance_map_filename(kinematic_token: object) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(_filename_token(kinematic_token))


def write_pion_hgcer_method_a_acceptance_map_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    target = os.fspath(path)
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(_json_copy(payload, "acceptance_map_artifact"), handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return target


def evaluate_relative_response_grid(model: Mapping[str, object], raw_values: np.ndarray, training_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate a transient plotting grid and mask points outside model support."""
    values = np.asarray(raw_values, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(ACCEPTED_FEATURES) or not np.all(np.isfinite(values)):
        raise MethodAAcceptanceMapError("map_grid_values_invalid")
    scaler = _mapping(model.get("scaler"), "model_scaler")
    fit = _mapping(model.get("logistic_fit"), "model_logistic_fit")
    support = _mapping(model.get("application_support"), "model_support")
    median = np.asarray(scaler.get("median"), dtype=float); divisor = np.asarray(scaler.get("divisor"), dtype=float)
    coefficients = np.asarray(fit.get("coefficients"), dtype=float)
    if median.shape != (3,) or divisor.shape != (3,) or coefficients.shape != (3,) or np.any(divisor <= 0.0):
        raise MethodAAcceptanceMapError("map_model_parameters_invalid")
    scaled = (values - median) / divisor
    training = np.asarray(training_values, dtype=float)
    if training.ndim != 2 or training.shape[1] != len(ACCEPTED_FEATURES) or training.shape[0] < 2 or not np.all(np.isfinite(training)):
        raise MethodAAcceptanceMapError("map_model_training_values_invalid")
    threshold = _finite(support.get("support_distance_threshold"), "map_support_distance_threshold")
    tree = cKDTree((training - median) / divisor)
    in_support = np.asarray(tree.query(scaled, k=1)[0], dtype=float) <= threshold
    response = np.full(values.shape[0], np.nan, dtype=float)
    response[in_support] = np.exp((scaled @ coefficients)[in_support])
    if not np.all(np.isfinite(response[in_support])):
        raise MethodAAcceptanceMapError("map_grid_response_nonfinite")
    return response, in_support


__all__ = (
    "ACCEPTED_BASIS", "ACCEPTED_FEATURES", "DEFAULT_ALGORITHM_CONFIG",
    "METHOD_A_ACCEPTANCE_MAP_ARTIFACT_SCHEMA_VERSION", "METHOD_A_ACCEPTANCE_MAP_FINGERPRINT_SCHEMA_VERSION", "METHOD_A_ACCEPTANCE_MAP_SCHEMA_VERSION", "MethodAAcceptanceMapError", "build_pion_hgcer_method_a_acceptance_map", "build_pion_hgcer_method_a_acceptance_map_artifact", "evaluate_relative_response_grid", "pion_hgcer_method_a_acceptance_map_filename", "write_pion_hgcer_method_a_acceptance_map_json",
)
