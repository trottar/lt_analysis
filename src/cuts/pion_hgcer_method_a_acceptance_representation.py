"""Detached global Method-A acceptance-representation audit.

Phase F.2 consumes only serialized F.1 v2 artifacts.  It compares a fixed
set of acceptance-coordinate bases for a future response map, but constructs
no map, probability, correction, weight, or production-side object.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import hashlib
import json
import math
import os

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.special import expit


METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_representation/v1"
)
METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_representation_fingerprint/v1"
)
METHOD_A_ACCEPTANCE_REPRESENTATION_ARTIFACT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_representation_artifact/v1"
)

_F1_ARTIFACT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract_artifact/v2"
_F1_CONTRACT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract/v2"
_F1_FINGERPRINT_SCHEMA = "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2"
_PRIMARY_FEATURES = (
    "SHMS_delta",
    "SHMS_xptar",
    "SHMS_yptar",
    "P_hgcer_xAtCer",
    "P_hgcer_yAtCer",
)
_CANONICAL_SETTINGS = (
    ("Left", "lowe"),
    ("Left", "highe"),
    ("Center", "lowe"),
    ("Center", "highe"),
    ("Right", "highe"),
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
    "phase_a_contract_fingerprint",
    "phase_a_pion_event_population_fingerprint",
    "method_a_fingerprint",
    "method_a_event_population_fingerprint",
    "part1_config_fingerprint",
    "coordinate_fingerprint",
)
_F1_CHILD_ASSIGNMENT_FIELDS = (
    "source_label", "entry_index", "t_index", "phi_index", "phi_low",
    "phi_high", "phi_status",
)

CANDIDATE_REPRESENTATIONS = (
    {
        "candidate_id": "delta_only",
        "ordered_features": ("SHMS_delta",),
        "dimension": 1,
        "automatic_recommendation_eligible": True,
        "role": "reduced_reference",
    },
    {
        "candidate_id": "track3",
        "ordered_features": ("SHMS_delta", "SHMS_xptar", "SHMS_yptar"),
        "dimension": 3,
        "automatic_recommendation_eligible": True,
        "role": "reduced_track_coordinate",
    },
    {
        "candidate_id": "hgcer3",
        "ordered_features": (
            "SHMS_delta",
            "P_hgcer_xAtCer",
            "P_hgcer_yAtCer",
        ),
        "dimension": 3,
        "automatic_recommendation_eligible": True,
        "role": "reduced_detector_coordinate",
    },
    {
        "candidate_id": "full5_reference",
        "ordered_features": _PRIMARY_FEATURES,
        "dimension": 5,
        "automatic_recommendation_eligible": False,
        "role": "diagnostic_information_reference_only",
    },
)

DEFAULT_ALGORITHM_CONFIG = {
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


class MethodAAcceptanceRepresentationError(ValueError):
    """Raised when frozen F.1 input cannot support an F.2 audit."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _json_copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRepresentationError("{}_not_json_safe".format(label)) from exc


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label))
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodAAcceptanceRepresentationError("{}_nonfinite".format(label))
    try:
        resolved = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRepresentationError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(resolved):
        raise MethodAAcceptanceRepresentationError("{}_nonfinite".format(label))
    return resolved


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label))
    try:
        resolved = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label)) from exc
    if resolved != value:
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label))
    return resolved


def _string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise MethodAAcceptanceRepresentationError("{}_invalid".format(label))
    return value


def _require(value: Mapping[str, object], name: str, expected: object, label: str) -> None:
    if name not in value or value.get(name) != expected:
        raise MethodAAcceptanceRepresentationError("{}_{}".format(label, name))


def _identity(record: Mapping[str, object], label: str) -> tuple[str, int]:
    source = _string(record.get("source_label"), "{}_source_label".format(label))
    entry = _integer(record.get("entry_index"), "{}_entry_index".format(label))
    if entry < 0:
        raise MethodAAcceptanceRepresentationError("{}_entry_index_invalid".format(label))
    return source, entry


def _setting_id(phi: str, epsilon: str) -> str:
    return "{}-{}".format(phi, epsilon)


def _candidate_definitions() -> list[dict[str, object]]:
    return [
        {
            "candidate_id": candidate["candidate_id"],
            "ordered_features": list(candidate["ordered_features"]),
            "dimension": candidate["dimension"],
            "automatic_recommendation_eligible": candidate["automatic_recommendation_eligible"],
            "role": candidate["role"],
        }
        for candidate in CANDIDATE_REPRESENTATIONS
    ]


def _resolve_config(algorithm_config: Mapping[str, object] | None) -> dict[str, object]:
    config = dict(DEFAULT_ALGORITHM_CONFIG)
    if algorithm_config is not None:
        supplied = _mapping(algorithm_config, "algorithm_config")
        unexpected = sorted(set(supplied) - set(config))
        if unexpected:
            raise MethodAAcceptanceRepresentationError(
                "algorithm_config_unknown:{}".format(",".join(unexpected))
            )
        config.update(supplied)
    if config["fold_count"] != 5:
        raise MethodAAcceptanceRepresentationError("algorithm_config_fold_count_invalid")
    if config["optimizer"] != "scipy_lbfgsb_analytic_gradient/v1":
        raise MethodAAcceptanceRepresentationError("algorithm_config_optimizer_invalid")
    for key in (
        "minimum_low_count", "minimum_control_count", "optimizer_maxiter",
        "support_sparse_nonprompt_count",
    ):
        if _integer(config[key], "algorithm_config_{}".format(key)) <= 0:
            raise MethodAAcceptanceRepresentationError("algorithm_config_{}_invalid".format(key))
    for key in (
        "regularization_lambda", "optimizer_gtol",
        "information_median_auc_loss_max", "information_max_auc_loss_max",
        "information_median_balanced_log_loss_penalty_max",
        "information_max_balanced_log_loss_penalty_max", "support_nn_percentile",
        "support_ood_fraction_max",
    ):
        if _finite(config[key], "algorithm_config_{}".format(key)) < 0.0:
            raise MethodAAcceptanceRepresentationError("algorithm_config_{}_invalid".format(key))
    if not 0.0 < float(config["support_nn_percentile"]) <= 100.0:
        raise MethodAAcceptanceRepresentationError("algorithm_config_support_nn_percentile_invalid")
    return _json_copy(config, "algorithm_config")  # type: ignore[return-value]


def _validate_setting(value: object) -> dict[str, object]:
    setting = _mapping(value, "f1_setting")
    phi = _string(setting.get("phi_setting"), "f1_setting_phi")
    epsilon = _string(setting.get("epsilon_filename_token"), "f1_setting_epsilon").lower()
    if (phi, epsilon) not in _CANONICAL_SETTINGS:
        raise MethodAAcceptanceRepresentationError("f1_setting_not_canonical")
    if _string(setting.get("particle_type"), "f1_setting_particle").lower() != "kaon":
        raise MethodAAcceptanceRepresentationError("f1_setting_particle_invalid")
    if _string(setting.get("kinematic_token"), "f1_setting_kinematic") != setting.get("kinematic_token"):
        raise MethodAAcceptanceRepresentationError("f1_setting_kinematic_invalid")
    return {
        "phi_setting": phi,
        "epsilon_filename_token": epsilon,
        "epsilon_setting": _string(setting.get("epsilon_setting"), "f1_setting_epsilon_setting").lower(),
        "particle_type": "kaon",
        "kinematic_token": setting["kinematic_token"],
        "Q2": _json_copy(setting.get("Q2"), "f1_setting_Q2"),
        "W": _json_copy(setting.get("W"), "f1_setting_W"),
    }


def _validate_metadata(value: object) -> None:
    metadata = _mapping(value, "f1_feature_metadata")
    if dict(metadata) != _F1_FEATURE_METADATA:
        raise MethodAAcceptanceRepresentationError("f1_feature_metadata_invalid")


def _edge_array(value: object, label: str) -> tuple[float, ...]:
    raw_edges = _sequence(value, label)
    if len(raw_edges) < 2:
        raise MethodAAcceptanceRepresentationError("{}_length_invalid".format(label))
    edges = tuple(_finite(edge, "{}_{}".format(label, index)) for index, edge in enumerate(raw_edges))
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise MethodAAcceptanceRepresentationError("{}_not_strictly_increasing".format(label))
    return edges


def _validate_t_assignment(record: object, t_edges: Sequence[float], label: str) -> None:
    value = _mapping(record, label)
    t_index = _integer(value.get("t_index"), "{}_t_index".format(label))
    if t_index < 0 or t_index >= len(t_edges) - 1:
        raise MethodAAcceptanceRepresentationError("{}_t_index_invalid".format(label))
    t_low = _finite(value.get("t_low"), "{}_t_low".format(label))
    t_high = _finite(value.get("t_high"), "{}_t_high".format(label))
    if t_low != t_edges[t_index] or t_high != t_edges[t_index + 1]:
        raise MethodAAcceptanceRepresentationError("{}_t_geometry_mismatch".format(label))


def _reconstruct_f1_fingerprints(
    contract: Mapping[str, object], raw_training: Sequence[object], raw_application: Sequence[object],
    feature_metadata: Mapping[str, object], label: str,
) -> dict[str, str]:
    """Verify the v2 producer's exact serialized F.1 fingerprint inputs."""
    provenance = {
        name: _string(contract.get(name), "{}_{}".format(label, name))
        for name in _F1_PROVENANCE_FINGERPRINT_NAMES
    }
    host_state = _string(contract.get("host_state"), "{}_host_state".format(label))
    _require(contract, "source_target_state", "post_proton_noRF", label)
    try:
        method_a_closure = _mapping(contract.get("method_a_training_summary"), "{}_training_summary".format(label))["by_t_delta"]
    except KeyError as exc:
        raise MethodAAcceptanceRepresentationError("{}_method_a_closure_invalid".format(label)) from exc
    child_projection: list[dict[str, object]] = []
    for record_index, raw_record in enumerate(raw_application):
        record = _mapping(raw_record, "{}_application_{}".format(label, record_index))
        missing = [name for name in _F1_CHILD_ASSIGNMENT_FIELDS if name not in record]
        if missing:
            raise MethodAAcceptanceRepresentationError(
                "{}_application_child_projection_missing:{}".format(label, ",".join(missing))
            )
        child_projection.append({name: record[name] for name in _F1_CHILD_ASSIGNMENT_FIELDS})
    values = {
        "method_a_training_population_fingerprint": _sha256(raw_training),
        "application_population_fingerprint": _sha256(raw_application),
        "acceptance_feature_metadata_fingerprint": _sha256(feature_metadata),
        "application_child_assignment_projection_fingerprint": _sha256(child_projection),
    }
    for name, expected in values.items():
        if contract.get(name) != expected:
            raise MethodAAcceptanceRepresentationError("{}_{}_mismatch".format(label, name))
    fingerprint_inputs = {
        "schema_version": _F1_CONTRACT_SCHEMA,
        "fingerprint_schema_version": _F1_FINGERPRINT_SCHEMA,
        **provenance,
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
        "t_edges": contract.get("t_edges"),
        "delta_edges": contract.get("delta_edges"),
        "phi_edges": contract.get("phi_edges"),
        **values,
        "method_a_closure": method_a_closure,
        "feature_metadata": feature_metadata,
    }
    if contract.get("fingerprint_inputs") != fingerprint_inputs:
        raise MethodAAcceptanceRepresentationError("{}_fingerprint_inputs_mismatch".format(label))
    if contract.get("fingerprint") != _sha256(fingerprint_inputs):
        raise MethodAAcceptanceRepresentationError("{}_fingerprint_mismatch".format(label))
    return {
        **provenance,
        **values,
        "coordinate_fingerprint": provenance["coordinate_fingerprint"],
        "fingerprint": str(contract["fingerprint"]),
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
    }


def _stable_f1_content_fingerprint(
    setting: Mapping[str, object], contract: Mapping[str, object],
    raw_training: Sequence[object], raw_application: Sequence[object],
) -> str:
    """Bind F.2's scientific fingerprint to F.1 content, not record order.

    The frozen F.1 producer fingerprints its serialized population order.  F.2
    verifies those producer fingerprints above, but its own deterministic audit
    fingerprint must remain invariant if the same valid records are reordered.
    """
    authority_names = (
        "schema_version", "fingerprint_schema_version", "status", "available", "reason",
        "diagnostic_stage", "non_authoritative", "production_objects_mutated",
        "refinement_applied", "production_application_performed",
        "event_application_performed", "method_b_numerical_dependency",
        "future_weight_adjustment_constructed", "phase_a_contract_fingerprint",
        "phase_a_pion_event_population_fingerprint", "method_a_fingerprint",
        "method_a_event_population_fingerprint", "part1_config_fingerprint",
        "coordinate_fingerprint", "host_state", "source_target_state", "t_edges",
        "delta_edges", "phi_edges", "method_a_training_summary", "feature_metadata",
    )
    return _sha256({
        "setting": setting,
        "contract_authority": {name: contract.get(name) for name in authority_names},
        "method_a_training_records": sorted(
            (_json_copy(row, "f1_training_record") for row in raw_training), key=_canonical_json,
        ),
        "application_records": sorted(
            (_json_copy(row, "f1_application_record") for row in raw_application), key=_canonical_json,
        ),
    })


def _sanitize_training_record(value: object, label: str) -> dict[str, object]:
    record = _mapping(value, label)
    source, entry = _identity(record, label)
    if source != _PROMPT_SOURCE or record.get("nommcuts") is not True:
        raise MethodAAcceptanceRepresentationError("{}_selection_invalid".format(label))
    npe = _finite(record.get("P_hgcer_npeSum"), "{}_npe".format(label))
    response_class = _string(record.get("response_class"), "{}_response_class".format(label))
    expected = "low" if 0.0 < npe <= _LOW_NPE_UPPER_BOUND else "control" if npe > _LOW_NPE_UPPER_BOUND else None
    if expected is None or response_class != expected:
        raise MethodAAcceptanceRepresentationError("{}_response_classification_invalid".format(label))
    result: dict[str, object] = {
        "source_label": source,
        "entry_index": entry,
        "response_class": response_class,
        "P_hgcer_npeSum": npe,
        "t_index": _integer(record.get("t_index"), "{}_t_index".format(label)),
        "t_low": _finite(record.get("t_low"), "{}_t_low".format(label)),
        "t_high": _finite(record.get("t_high"), "{}_t_high".format(label)),
    }
    if result["t_low"] >= result["t_high"]:
        raise MethodAAcceptanceRepresentationError("{}_t_geometry_invalid".format(label))
    for feature in _PRIMARY_FEATURES:
        result[feature] = _finite(record.get(feature), "{}_{}".format(label, feature))
    return result


def _sanitize_application_record(value: object, label: str) -> dict[str, object]:
    record = _mapping(value, label)
    source, entry = _identity(record, label)
    npe = _finite(record.get("P_hgcer_npeSum"), "{}_npe".format(label))
    if npe <= _LOW_NPE_UPPER_BOUND:
        raise MethodAAcceptanceRepresentationError("{}_npe_not_physical_control".format(label))
    result: dict[str, object] = {
        "source_label": source,
        "entry_index": entry,
        "P_hgcer_npeSum": npe,
        "t_index": _integer(record.get("t_index"), "{}_t_index".format(label)),
        "t_low": _finite(record.get("t_low"), "{}_t_low".format(label)),
        "t_high": _finite(record.get("t_high"), "{}_t_high".format(label)),
    }
    if result["t_low"] >= result["t_high"]:
        raise MethodAAcceptanceRepresentationError("{}_t_geometry_invalid".format(label))
    for feature in _PRIMARY_FEATURES:
        result[feature] = _finite(record.get(feature), "{}_{}".format(label, feature))
    return result


def _validate_f1_artifacts(f1_artifacts: object) -> list[dict[str, object]]:
    raw_artifacts = _sequence(f1_artifacts, "f1_artifacts")
    if len(raw_artifacts) != len(_CANONICAL_SETTINGS):
        raise MethodAAcceptanceRepresentationError("f1_setting_count_invalid")
    parsed: list[dict[str, object]] = []
    seen_settings: set[tuple[str, str]] = set()
    kinematics: set[str] = set()
    geometry_reference: tuple[tuple[int, float, float], ...] | None = None
    for artifact_index, raw_artifact in enumerate(raw_artifacts):
        artifact = _mapping(raw_artifact, "f1_artifact_{}".format(artifact_index))
        _require(artifact, "schema_version", _F1_ARTIFACT_SCHEMA, "f1_artifact")
        for name in (
            "non_authoritative", "production_objects_mutated", "refinement_applied",
            "production_application_performed", "event_application_performed",
        ):
            expected = True if name == "non_authoritative" else False
            _require(artifact, name, expected, "f1_artifact")
        setting = _validate_setting(artifact.get("setting"))
        setting_key = (str(setting["phi_setting"]), str(setting["epsilon_filename_token"]))
        if setting_key in seen_settings:
            raise MethodAAcceptanceRepresentationError("f1_setting_duplicate")
        seen_settings.add(setting_key)
        kinematics.add(str(setting["kinematic_token"]))
        contract = _mapping(artifact.get("contract"), "f1_contract")
        _require(contract, "schema_version", _F1_CONTRACT_SCHEMA, "f1_contract")
        _require(contract, "fingerprint_schema_version", _F1_FINGERPRINT_SCHEMA, "f1_contract")
        _require(contract, "status", "available", "f1_contract")
        _require(contract, "available", True, "f1_contract")
        _require(contract, "diagnostic_stage", "complete", "f1_contract")
        _require(contract, "reason", None, "f1_contract")
        for name in (
            "non_authoritative", "method_b_numerical_dependency",
            "event_application_performed", "production_application_performed",
            "production_objects_mutated", "future_weight_adjustment_constructed",
            "refinement_applied",
        ):
            expected = True if name == "non_authoritative" else False
            _require(contract, name, expected, "f1_contract")
        feature_metadata = _mapping(contract.get("feature_metadata"), "f1_feature_metadata")
        _validate_metadata(feature_metadata)
        t_edges = _edge_array(contract.get("t_edges"), "f1_t_edges")
        _edge_array(contract.get("delta_edges"), "f1_delta_edges")
        _edge_array(contract.get("phi_edges"), "f1_phi_edges")
        raw_training = _sequence(contract.get("method_a_training_records"), "f1_training_records")
        raw_application = _sequence(contract.get("application_records"), "f1_application_records")
        if not raw_training or not raw_application:
            raise MethodAAcceptanceRepresentationError("f1_population_empty")
        for row_index, row in enumerate(raw_training):
            _validate_t_assignment(row, t_edges, "f1_training_{}_{}".format(artifact_index, row_index))
        for row_index, row in enumerate(raw_application):
            _validate_t_assignment(row, t_edges, "f1_application_{}_{}".format(artifact_index, row_index))
        fingerprints = _reconstruct_f1_fingerprints(
            contract, raw_training, raw_application, feature_metadata,
            "f1_contract_{}".format(artifact_index),
        )
        stable_content_fingerprint = _stable_f1_content_fingerprint(
            setting, contract, raw_training, raw_application,
        )
        training = [
            _sanitize_training_record(row, "f1_training_{}_{}".format(artifact_index, row_index))
            for row_index, row in enumerate(raw_training)
        ]
        application = [
            _sanitize_application_record(row, "f1_application_{}_{}".format(artifact_index, row_index))
            for row_index, row in enumerate(raw_application)
        ]
        training.sort(key=lambda row: (str(row["source_label"]), int(row["entry_index"])))
        application.sort(key=lambda row: (str(row["source_label"]), int(row["entry_index"])))
        training_ids = [(str(row["source_label"]), int(row["entry_index"])) for row in training]
        application_ids = [(str(row["source_label"]), int(row["entry_index"])) for row in application]
        if len(training_ids) != len(set(training_ids)) or len(application_ids) != len(set(application_ids)):
            raise MethodAAcceptanceRepresentationError("f1_record_identity_duplicate")
        control_ids = {
            (str(row["source_label"]), int(row["entry_index"]))
            for row in training if row["response_class"] == "control"
        }
        if any(
            (str(row["source_label"]), int(row["entry_index"])) not in control_ids
            for row in application if row["source_label"] == _PROMPT_SOURCE
        ):
            raise MethodAAcceptanceRepresentationError("f1_prompt_application_identity_missing_training_control")
        ordered_geometry = tuple(
            (index, t_edges[index], t_edges[index + 1])
            for index in range(len(t_edges) - 1)
        )
        if geometry_reference is None:
            geometry_reference = ordered_geometry
        elif geometry_reference != ordered_geometry:
            raise MethodAAcceptanceRepresentationError("f1_canonical_t_geometry_mixed")
        parsed.append({
            "setting": setting,
            "setting_id": _setting_id(*setting_key),
            "training": training,
            "application": application,
            "fingerprints": fingerprints,
            "stable_content_fingerprint": stable_content_fingerprint,
            "canonical_t_geometry": ordered_geometry,
        })
    if seen_settings != set(_CANONICAL_SETTINGS):
        raise MethodAAcceptanceRepresentationError("f1_canonical_settings_missing_or_unexpected")
    if len(kinematics) != 1:
        raise MethodAAcceptanceRepresentationError("f1_kinematic_mixed")
    parsed.sort(key=lambda item: _CANONICAL_SETTINGS.index((
        str(_mapping(item["setting"], "setting")["phi_setting"]),
        str(_mapping(item["setting"], "setting")["epsilon_filename_token"]),
    )))
    return parsed


def _feature_array(records: Sequence[Mapping[str, object]], features: Sequence[str]) -> np.ndarray:
    return np.asarray([[float(record[name]) for name in features] for record in records], dtype=float)


def _robust_scale(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    median = np.median(values, axis=0)
    iqr = np.percentile(values, 75.0, axis=0) - np.percentile(values, 25.0, axis=0)
    scale = np.asarray(iqr, dtype=float)
    zero = ~np.isfinite(scale) | (scale <= 0.0)
    if np.any(zero):
        std = np.std(values, axis=0)
        scale[zero] = std[zero]
    if not np.all(np.isfinite(median)) or not np.all(np.isfinite(scale)) or np.any(scale <= 0.0):
        return None
    return np.asarray(median, dtype=float), scale


def _stable_folds(records: Sequence[Mapping[str, object]], fold_count: int) -> np.ndarray:
    folds = np.empty(len(records), dtype=int)
    for response_class in ("low", "control"):
        indexes = [
            index for index, record in enumerate(records)
            if record["response_class"] == response_class
        ]
        indexes.sort(key=lambda index: (
            str(records[index]["source_label"]), int(records[index]["entry_index"])
        ))
        for rank, index in enumerate(indexes):
            folds[index] = rank % fold_count
    return folds


def _auc(labels: np.ndarray, scores: np.ndarray) -> float:
    positive = labels == 1.0
    n_positive = int(np.sum(positive))
    n_negative = int(labels.size - n_positive)
    if n_positive == 0 or n_negative == 0:
        raise MethodAAcceptanceRepresentationError("probe_auc_class_missing")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, dtype=float)
    start = 0
    while start < scores.size:
        end = start + 1
        while end < scores.size and scores[order[end]] == scores[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * ((start + 1) + end)
        start = end
    return float((np.sum(ranks[positive]) - n_positive * (n_positive + 1) / 2.0) / (n_positive * n_negative))


def _logistic_objective(
    parameters: np.ndarray, values: np.ndarray, labels: np.ndarray, weights: np.ndarray,
    regularization_lambda: float,
) -> tuple[float, np.ndarray]:
    intercept = parameters[0]
    coefficients = parameters[1:]
    logits = intercept + values @ coefficients
    losses = np.logaddexp(0.0, logits) - labels * logits
    objective = float(np.mean(weights * losses) + 0.5 * regularization_lambda * np.dot(coefficients, coefficients))
    residual = weights * (expit(logits) - labels)
    gradient = np.empty_like(parameters)
    gradient[0] = float(np.mean(residual))
    gradient[1:] = values.T @ residual / values.shape[0] + regularization_lambda * coefficients
    return objective, gradient


def _probe_candidate(
    records: Sequence[Mapping[str, object]], candidate: Mapping[str, object], config: Mapping[str, object],
) -> dict[str, object]:
    features = tuple(candidate["ordered_features"])
    fold_count = int(config["fold_count"])
    labels = np.asarray([1.0 if row["response_class"] == "low" else 0.0 for row in records], dtype=float)
    folds = _stable_folds(records, fold_count)
    scores = np.empty(len(records), dtype=float)
    heldout_logits = np.empty(len(records), dtype=float)
    for fold in range(fold_count):
        training_index = folds != fold
        heldout_index = folds == fold
        values_train = _feature_array([records[index] for index in np.flatnonzero(training_index)], features)
        scale = _robust_scale(values_train)
        if scale is None:
            return {
                "low_count": int(np.sum(labels == 1.0)),
                "control_count": int(np.sum(labels == 0.0)),
                "roc_auc": None,
                "balanced_log_loss": None,
                "fold_count": fold_count,
                "fit_converged_all_folds": False,
                "valid": False,
                "invalid_reason": "fold_scaling_invalid",
            }
        median, divisor = scale
        scaled_train = (values_train - median) / divisor
        train_labels = labels[training_index]
        low_count = int(np.sum(train_labels == 1.0))
        control_count = int(np.sum(train_labels == 0.0))
        if low_count == 0 or control_count == 0:
            return {
                "low_count": int(np.sum(labels == 1.0)),
                "control_count": int(np.sum(labels == 0.0)),
                "roc_auc": None,
                "balanced_log_loss": None,
                "fold_count": fold_count,
                "fit_converged_all_folds": False,
                "valid": False,
                "invalid_reason": "fold_class_missing",
            }
        weights = np.where(
            train_labels == 1.0,
            train_labels.size / (2.0 * low_count),
            train_labels.size / (2.0 * control_count),
        )
        result = minimize(
            _logistic_objective,
            np.zeros(scaled_train.shape[1] + 1, dtype=float),
            args=(scaled_train, train_labels, weights, float(config["regularization_lambda"])),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": int(config["optimizer_maxiter"]), "gtol": float(config["optimizer_gtol"])},
        )
        if not result.success or not np.all(np.isfinite(result.x)):
            return {
                "low_count": int(np.sum(labels == 1.0)),
                "control_count": int(np.sum(labels == 0.0)),
                "roc_auc": None,
                "balanced_log_loss": None,
                "fold_count": fold_count,
                "fit_converged_all_folds": False,
                "valid": False,
                "invalid_reason": "optimizer_not_converged",
            }
        values_heldout = _feature_array([records[index] for index in np.flatnonzero(heldout_index)], features)
        logits = result.x[0] + ((values_heldout - median) / divisor) @ result.x[1:]
        heldout_logits[heldout_index] = logits
        scores[heldout_index] = expit(logits)
    losses = np.logaddexp(0.0, heldout_logits) - labels * heldout_logits
    return {
        "low_count": int(np.sum(labels == 1.0)),
        "control_count": int(np.sum(labels == 0.0)),
        "roc_auc": _auc(labels, scores),
        "balanced_log_loss": float(0.5 * np.mean(losses[labels == 1.0]) + 0.5 * np.mean(losses[labels == 0.0])),
        "fold_count": fold_count,
        "fit_converged_all_folds": True,
        "valid": True,
        "invalid_reason": None,
    }


def _percentiles(values: np.ndarray) -> dict[str, float | None]:
    if values.size == 0:
        return {"p50": None, "p95": None, "p99": None, "max": None}
    return {
        "p50": float(np.percentile(values, 50.0)),
        "p95": float(np.percentile(values, 95.0)),
        "p99": float(np.percentile(values, 99.0)),
        "max": float(np.max(values)),
    }


def _support_candidate(
    training: Sequence[Mapping[str, object]], application: Sequence[Mapping[str, object]],
    candidate: Mapping[str, object], config: Mapping[str, object],
) -> dict[str, object]:
    features = tuple(candidate["ordered_features"])
    nonprompt_count = sum(row["source_label"] != _PROMPT_SOURCE for row in application)
    statistically_sparse = nonprompt_count < int(config["support_sparse_nonprompt_count"])
    training_values = _feature_array(training, features)
    scale = _robust_scale(training_values)
    if scale is None:
        return {
            "valid": False, "invalid_reason": "training_scaling_invalid",
            "nonprompt_application_count": nonprompt_count,
            "training_nn_p50": None, "training_nn_p95": None, "training_nn_p99": None,
            "application_nn_p50": None, "application_nn_p95": None, "application_nn_p99": None,
            "application_nn_max": None, "application_ood_count": None,
            "application_ood_fraction": None, "statistically_sparse": statistically_sparse,
            "support_gate_passed": False,
        }
    median, divisor = scale
    scaled_training = (training_values - median) / divisor
    if scaled_training.shape[0] < 2:
        return {
            "valid": False, "invalid_reason": "training_neighbor_reference_invalid",
            "nonprompt_application_count": nonprompt_count,
            "training_nn_p50": None, "training_nn_p95": None, "training_nn_p99": None,
            "application_nn_p50": None, "application_nn_p95": None, "application_nn_p99": None,
            "application_nn_max": None, "application_ood_count": None,
            "application_ood_fraction": None, "statistically_sparse": statistically_sparse,
            "support_gate_passed": False,
        }
    tree = cKDTree(scaled_training)
    training_distances = np.asarray(tree.query(scaled_training, k=2)[0][:, 1], dtype=float)
    threshold = float(np.percentile(training_distances, float(config["support_nn_percentile"])))
    nonprompt = [row for row in application if row["source_label"] != _PROMPT_SOURCE]
    app_distances = (
        np.asarray(tree.query((_feature_array(nonprompt, features) - median) / divisor, k=1)[0], dtype=float)
        if nonprompt else np.asarray([], dtype=float)
    )
    training_summary = _percentiles(training_distances)
    application_summary = _percentiles(app_distances)
    ood_count = int(np.sum(app_distances > threshold))
    ood_fraction = float(ood_count / app_distances.size) if app_distances.size else 0.0
    return {
        "valid": True, "invalid_reason": None,
        "nonprompt_application_count": int(app_distances.size),
        "training_nn_p50": training_summary["p50"],
        "training_nn_p95": training_summary["p95"],
        "training_nn_p99": training_summary["p99"],
        "application_nn_p50": application_summary["p50"],
        "application_nn_p95": application_summary["p95"],
        "application_nn_p99": application_summary["p99"],
        "application_nn_max": application_summary["max"],
        "application_ood_count": ood_count,
        "application_ood_fraction": ood_fraction,
        "statistically_sparse": statistically_sparse,
        "support_gate_passed": bool(ood_fraction <= float(config["support_ood_fraction_max"])),
    }


def _normalize_hashes(
    parsed: Sequence[Mapping[str, object]], input_file_hashes: Mapping[str, object] | None,
) -> list[dict[str, object]]:
    supplied = {} if input_file_hashes is None else _mapping(input_file_hashes, "input_file_hashes")
    expected = {str(item["setting_id"]) for item in parsed}
    if supplied and set(supplied) != expected:
        raise MethodAAcceptanceRepresentationError("input_file_hashes_settings_invalid")
    result = []
    for item in parsed:
        setting_id = str(item["setting_id"])
        value = supplied.get(setting_id)
        if value is not None and (not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value.lower())):
            raise MethodAAcceptanceRepresentationError("input_file_hash_invalid")
        result.append({"setting_id": setting_id, "sha256": value})
    return result


def _group_inventory(parsed: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    groups: list[dict[str, object]] = []
    for item in parsed:
        training_by_t: dict[int, list[Mapping[str, object]]] = {}
        application_by_t: dict[int, list[Mapping[str, object]]] = {}
        for record in _sequence(item["training"], "training"):
            training_by_t.setdefault(int(_mapping(record, "training")["t_index"]), []).append(_mapping(record, "training"))
        for record in _sequence(item["application"], "application"):
            application_by_t.setdefault(int(_mapping(record, "application")["t_index"]), []).append(_mapping(record, "application"))
        for t_index, t_low, t_high in _sequence(item["canonical_t_geometry"], "geometry"):
            training = training_by_t.get(int(t_index), [])
            application = application_by_t.get(int(t_index), [])
            low_count = sum(record["response_class"] == "low" for record in training)
            control_count = sum(record["response_class"] == "control" for record in training)
            groups.append({
                "setting": _json_copy(item["setting"], "setting"),
                "setting_id": item["setting_id"],
                "canonical_t_index": int(t_index), "canonical_t_low": float(t_low), "canonical_t_high": float(t_high),
                "training_positive_count": len(training), "low_count": int(low_count),
                "control_count": int(control_count),
                "nonprompt_application_count": sum(record["source_label"] != _PROMPT_SOURCE for record in application),
                "_training": training, "_application": application,
                "candidate_metrics": {},
            })
    return groups


def _candidate_summary(
    candidate: Mapping[str, object], groups: Sequence[Mapping[str, object]], config: Mapping[str, object],
) -> dict[str, object]:
    candidate_id = str(candidate["candidate_id"])
    reference_id = "full5_reference"
    losses, penalties = [], []
    invalid_information = 0
    invalid_support = 0
    sparse_groups = 0
    support_group_passes = []
    for group in groups:
        metrics = _mapping(_mapping(group["candidate_metrics"], "candidate_metrics")[candidate_id], "candidate_metrics")
        reference = _mapping(_mapping(group["candidate_metrics"], "candidate_metrics")[reference_id], "reference_metrics")
        if candidate_id == reference_id:
            metrics["auc_loss_relative_to_full5_reference"] = 0.0 if metrics["valid"] else None
            metrics["balanced_log_loss_penalty_relative_to_full5_reference"] = 0.0 if metrics["valid"] else None
        elif metrics["valid"] and reference["valid"]:
            metrics["auc_loss_relative_to_full5_reference"] = float(reference["roc_auc"]) - float(metrics["roc_auc"])
            metrics["balanced_log_loss_penalty_relative_to_full5_reference"] = float(metrics["balanced_log_loss"]) - float(reference["balanced_log_loss"])
        else:
            metrics["auc_loss_relative_to_full5_reference"] = None
            metrics["balanced_log_loss_penalty_relative_to_full5_reference"] = None
        if metrics["auc_loss_relative_to_full5_reference"] is None:
            invalid_information += 1
        else:
            losses.append(float(metrics["auc_loss_relative_to_full5_reference"]))
            penalties.append(float(metrics["balanced_log_loss_penalty_relative_to_full5_reference"]))
        support = _mapping(metrics["application_support"], "application_support")
        if not support["valid"]:
            invalid_support += 1
        else:
            support_group_passes.append(bool(support["support_gate_passed"]))
            sparse_groups += int(bool(support["statistically_sparse"]))
    information_passed = (
        invalid_information == 0
        and bool(losses)
        and float(np.median(losses)) <= float(config["information_median_auc_loss_max"])
        and float(np.max(losses)) <= float(config["information_max_auc_loss_max"])
        and float(np.median(penalties)) <= float(config["information_median_balanced_log_loss_penalty_max"])
        and float(np.max(penalties)) <= float(config["information_max_balanced_log_loss_penalty_max"])
    )
    support_passed = invalid_support == 0 and all(support_group_passes)
    eligible = bool(candidate["automatic_recommendation_eligible"])
    return {
        "candidate_id": candidate_id,
        "ordered_features": list(candidate["ordered_features"]),
        "dimension": int(candidate["dimension"]),
        "automatic_recommendation_eligible": eligible,
        "information_retention": {
            "median_auc_loss": float(np.median(losses)) if losses else None,
            "max_auc_loss": float(np.max(losses)) if losses else None,
            "median_balanced_log_loss_penalty": float(np.median(penalties)) if penalties else None,
            "max_balanced_log_loss_penalty": float(np.max(penalties)) if penalties else None,
            "valid_group_count": len(losses), "invalid_group_count": invalid_information,
            "information_gate_passed": information_passed,
        },
        "application_support": {
            "valid_group_count": len(groups) - invalid_support,
            "invalid_group_count": invalid_support,
            "statistically_sparse_group_count": sparse_groups,
            "application_support_gate_passed": support_passed,
        },
        "information_gate_passed": information_passed,
        "application_support_gate_passed": support_passed,
        "overall_candidate_passed": bool(eligible and information_passed and support_passed),
    }


def _recommendation(groups: Sequence[Mapping[str, object]], summaries: Sequence[Mapping[str, object]], support_ok: bool) -> dict[str, object]:
    if not support_ok:
        return {
            "recommendation_status": "no_supported_reduced_basis",
            "recommended_basis": None, "basis_frozen": False, "manual_review_required": True,
        }
    passing = [
        summary for summary in summaries
        if summary.get("automatic_recommendation_eligible") is True
        and summary["overall_candidate_passed"]
    ]
    if not passing:
        return {
            "recommendation_status": "no_supported_reduced_basis",
            "recommended_basis": None, "basis_frozen": False, "manual_review_required": True,
        }
    if len(passing) == 1:
        return {
            "recommendation_status": "unique_supported_reduced_basis",
            "recommended_basis": passing[0]["candidate_id"], "basis_frozen": False,
            "manual_review_required": True,
        }
    minimum_dimension = min(int(summary["dimension"]) for summary in passing)
    minimum = [summary for summary in passing if int(summary["dimension"]) == minimum_dimension]
    if len(minimum) == 1:
        return {
            "recommendation_status": "unique_minimum_dimension_supported_basis",
            "recommended_basis": minimum[0]["candidate_id"], "basis_frozen": False,
            "manual_review_required": True,
        }
    return {
        "recommendation_status": "multiple_supported_reduced_bases",
        "recommended_basis": None, "basis_frozen": False, "manual_review_required": True,
    }


def build_pion_hgcer_method_a_acceptance_representation(
    f1_artifacts: Sequence[Mapping[str, object]], *, algorithm_config: Mapping[str, object] | None = None,
    input_file_hashes: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Build the detached, deterministic global F.2 representation audit.

    Expected frozen-input defects raise :class:`MethodAAcceptanceRepresentationError`.
    Insufficient but otherwise valid response support produces a non-recommending
    F.2 result instead of pooling canonical-t parents.
    """
    parsed = _validate_f1_artifacts(f1_artifacts)
    config = _resolve_config(algorithm_config)
    hashes = _normalize_hashes(parsed, input_file_hashes)
    groups = _group_inventory(parsed)
    support_ok = all(
        int(group["low_count"]) >= int(config["minimum_low_count"])
        and int(group["control_count"]) >= int(config["minimum_control_count"])
        for group in groups
    )
    insufficient_response_groups = [
        {
            "setting_id": group["setting_id"],
            "canonical_t_index": group["canonical_t_index"],
            "low_count": group["low_count"],
            "control_count": group["control_count"],
        }
        for group in groups
        if int(group["low_count"]) < int(config["minimum_low_count"])
        or int(group["control_count"]) < int(config["minimum_control_count"])
    ]
    if support_ok:
        for group in groups:
            metrics = _mapping(group["candidate_metrics"], "candidate_metrics")
            for candidate in CANDIDATE_REPRESENTATIONS:
                candidate_metrics = _probe_candidate(
                    _sequence(group["_training"], "group_training"), candidate, config
                )
                candidate_metrics["application_support"] = _support_candidate(
                    _sequence(group["_training"], "group_training"),
                    _sequence(group["_application"], "group_application"), candidate, config,
                )
                metrics[str(candidate["candidate_id"])] = candidate_metrics
    else:
        for group in groups:
            metrics = _mapping(group["candidate_metrics"], "candidate_metrics")
            for candidate in CANDIDATE_REPRESENTATIONS:
                metrics[str(candidate["candidate_id"])] = {
                    "low_count": group["low_count"], "control_count": group["control_count"],
                    "roc_auc": None, "balanced_log_loss": None,
                    "fold_count": int(config["fold_count"]), "fit_converged_all_folds": False,
                    "valid": False, "invalid_reason": "insufficient_response_support",
                    "application_support": {
                        "valid": False, "invalid_reason": "not_evaluated_insufficient_response_support",
                        "nonprompt_application_count": group["nonprompt_application_count"],
                        "training_nn_p50": None, "training_nn_p95": None, "training_nn_p99": None,
                        "application_nn_p50": None, "application_nn_p95": None, "application_nn_p99": None,
                        "application_nn_max": None, "application_ood_count": None,
                        "application_ood_fraction": None,
                        "statistically_sparse": int(group["nonprompt_application_count"])
                        < int(config["support_sparse_nonprompt_count"]),
                        "support_gate_passed": False,
                    },
                }
    summaries = [_candidate_summary(candidate, groups, config) for candidate in CANDIDATE_REPRESENTATIONS]
    recommendation = _recommendation(groups, summaries, support_ok)
    public_groups = []
    for group in groups:
        public = {key: value for key, value in group.items() if not key.startswith("_")}
        public_groups.append(public)
    input_fingerprints = [
        {
            "setting": _json_copy(item["setting"], "input_setting"),
            "setting_id": item["setting_id"],
            "stable_f1_content_fingerprint": item["stable_content_fingerprint"],
            **_json_copy(item["fingerprints"], "input_fingerprints"),
            "source_file_sha256": next(entry["sha256"] for entry in hashes if entry["setting_id"] == item["setting_id"]),
        }
        for item in parsed
    ]
    candidate_definitions = _candidate_definitions()
    algorithm_fingerprint = _sha256({
        "representation_schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION,
        "fingerprint_schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION,
        "candidate_definitions": candidate_definitions,
        "algorithm_config": config,
    })
    fingerprint_inputs = {
        "representation_schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION,
        "fingerprint_schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION,
        "candidate_definitions": candidate_definitions,
        "algorithm_config": config,
        "algorithm_fingerprint": algorithm_fingerprint,
        "input_content": [
            {
                "setting_id": item["setting_id"],
                "stable_f1_content_fingerprint": item["stable_content_fingerprint"],
                "source_file_sha256": next(
                    entry["sha256"] for entry in hashes if entry["setting_id"] == item["setting_id"]
                ),
            }
            for item in parsed
        ],
        "response_support": {
            "all_groups_satisfy_minimum": support_ok,
            "minimum_low_count": int(config["minimum_low_count"]),
            "minimum_control_count": int(config["minimum_control_count"]),
            "insufficient_group_count": len(insufficient_response_groups),
            "insufficient_groups": insufficient_response_groups,
        },
        "groups": public_groups,
        "candidate_summaries": summaries,
        "recommendation": recommendation,
    }
    result_core = {
        "schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION,
        "fingerprint_schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION,
        "status": "available", "available": True, "diagnostic_stage": "complete",
        "non_authoritative": True, "probe_only": True,
        "probability_map_constructed": False, "event_probability_persisted": False,
        "event_application_performed": False, "production_application_performed": False,
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "basis_frozen": False, "manual_review_required": True,
        "primary_acceptance_features": list(_PRIMARY_FEATURES),
        "candidate_definitions": candidate_definitions,
        "algorithm_config": config,
        "algorithm_fingerprint": algorithm_fingerprint,
        "input_fingerprints": input_fingerprints,
        "response_support": fingerprint_inputs["response_support"],
        "groups": public_groups,
        "candidate_summaries": summaries,
        "recommendation": recommendation,
    }
    result = _json_copy(result_core, "representation_result")
    result["fingerprint_inputs"] = _json_copy(fingerprint_inputs, "representation_fingerprint_inputs")
    result["fingerprint"] = _sha256(fingerprint_inputs)
    return result  # type: ignore[return-value]


def build_pion_hgcer_method_a_acceptance_representation_artifact(
    f1_artifacts: Sequence[Mapping[str, object]], *, algorithm_config: Mapping[str, object] | None = None,
    input_file_hashes: Mapping[str, object] | None = None, input_paths: Mapping[str, object] | None = None,
    generated_at_utc: str | None = None, git_head: str | None = None,
    git_status_short: str | None = None,
) -> dict[str, object]:
    """Wrap one F.2 result with non-scientific CLI provenance."""
    representation = build_pion_hgcer_method_a_acceptance_representation(
        f1_artifacts, algorithm_config=algorithm_config, input_file_hashes=input_file_hashes,
    )
    provenance = {
        "generated_at_utc": generated_at_utc,
        "git_head": git_head,
        "git_status_short": git_status_short,
        "input_paths": _json_copy({} if input_paths is None else input_paths, "input_paths"),
    }
    return {
        "schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_ARTIFACT_SCHEMA_VERSION,
        "representation": representation,
        "non_authoritative": True, "probe_only": True,
        "probability_map_constructed": False, "event_probability_persisted": False,
        "event_application_performed": False, "production_application_performed": False,
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "basis_frozen": False, "manual_review_required": True,
        "provenance": provenance,
        "artifact_fingerprint": _sha256({
            "schema_version": METHOD_A_ACCEPTANCE_REPRESENTATION_ARTIFACT_SCHEMA_VERSION,
            "representation_fingerprint": representation["fingerprint"],
            "input_paths": provenance["input_paths"],
        }),
    }


def _filename_token(value: object) -> str:
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise ValueError("method_a_acceptance_representation_filename_token_invalid")
    token = str(value)
    if not token or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise ValueError("method_a_acceptance_representation_filename_token_invalid")
    return token


def pion_hgcer_method_a_acceptance_representation_filename(kinematic_token: object) -> str:
    """Return the deterministic global F.2 JSON basename."""
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-representation.json".format(
        _filename_token(kinematic_token)
    )


def write_pion_hgcer_method_a_acceptance_representation_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    """Write a JSON-safe deterministic F.2 artifact with one trailing newline."""
    target = os.fspath(path)
    serializable = _json_copy(payload, "representation_artifact")
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(serializable, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return target


__all__ = (
    "CANDIDATE_REPRESENTATIONS",
    "DEFAULT_ALGORITHM_CONFIG",
    "METHOD_A_ACCEPTANCE_REPRESENTATION_ARTIFACT_SCHEMA_VERSION",
    "METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION",
    "METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION",
    "MethodAAcceptanceRepresentationError",
    "build_pion_hgcer_method_a_acceptance_representation",
    "build_pion_hgcer_method_a_acceptance_representation_artifact",
    "pion_hgcer_method_a_acceptance_representation_filename",
    "write_pion_hgcer_method_a_acceptance_representation_json",
)
