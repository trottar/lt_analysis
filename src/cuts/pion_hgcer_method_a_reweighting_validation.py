"""Detached Phase F.6.1 validation of accepted Method-A event reweighting.

This module consumes frozen F.1/F.3/F.4/F.5 artifacts.  It evaluates F.4
event factors only while constructing aggregate validation material; factors,
event identities, and feature rows never enter its persisted payload.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os

import numpy as np

import pion_hgcer_method_a_parent_preserving_correction as _f4
import pion_hgcer_method_a_tphi_propagation as _f5


METHOD_A_REWEIGHTING_VALIDATION_SCHEMA_VERSION = "pion_hgcer_method_a_reweighting_validation/v1"
METHOD_A_REWEIGHTING_VALIDATION_FINGERPRINT_SCHEMA_VERSION = "pion_hgcer_method_a_reweighting_validation_fingerprint/v1"
METHOD_A_REWEIGHTING_VALIDATION_ARTIFACT_SCHEMA_VERSION = "pion_hgcer_method_a_reweighting_validation_artifact/v1"
_F5_ARTIFACT_SCHEMA = "pion_hgcer_method_a_tphi_propagation_artifact/v1"
_F5_SCHEMA = "pion_hgcer_method_a_tphi_propagation/v1"
_CANONICAL_SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)
_MODEL_VARIABLES = ("SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer")
_INDEPENDENT_VARIABLES = ("SHMS_xptar", "SHMS_yptar", "phi", "analysis_MM", "Q2", "W")
_REQUIRED_TRAINING_VARIABLES = _MODEL_VARIABLES + _INDEPENDENT_VARIABLES
_TOLERANCE = 1.0e-12
DISPLAY_POLICY = {
    "one_dimensional_nonphi_bin_count": 40,
    "hgcer_xy_bin_count_per_axis": 30,
    "phi_edges": "frozen_canonical_phi_edges/v1",
    "degenerate_finite_range_half_width": "max(1,abs(center))*1e-6/v1",
    "nonphi_range": "common_finite_low_and_prompt_physical_control/v1",
    "signed_missing_mass_range": "complete_physical_application_population/v1",
    "percentile_trimming": False,
    "clipping": False,
    "winsorization": False,
}


# Source-owned authority: these values must not be discovered dynamically by
# an analyzer, archive, or command-line argument.
ACCEPTED_RUNTIME_AUTHORITY_BY_KINEMATIC = {
    "Q4p4W2p74": {
        "f1_source_file_sha256": {
            "Left-lowe": "f16c5e89f62e848ec221335bd4916265757040af9970b005d8f67cc800e0077d",
            "Left-highe": "203f4c76f1a251e3e8f231fa3a5e50c9a4fa337420efffffd803205c6d7ea218",
            "Center-lowe": "1593e22b55382b4a9e831d3a1114584e2e3057fcbc4aeeea4aa74c948edf39f8",
            "Center-highe": "5de64b850735ebe70040a703bac999bdd1ac84dc821aba8e1a21ec86ae4b9db3",
            "Right-highe": "77d98006ff81e772c466bf8d10ec83088509a0b8a4d430bded1172bc118bc15f",
        },
        "f3_source_file_sha256": "04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95",
        "f3_map_fingerprint": "81b2a1e89ef9689b24c6dd145f53b26f7ac8fee9d5666cbc2ae8caa9da2830e6",
        "f3_algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
        "f3_artifact_fingerprint": "f8a12313bbd81aca48402a8c7eb4773c2dbe3ad0c7c6f64a2f212e98202d6ee2",
        "f4_source_file_sha256": "adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188",
        "f4_correction_fingerprint": "362241005c02f2149e260c391b5c3d35793287573128b42cf5ed693419d9d2f3",
        "f4_artifact_fingerprint": "c4b9f513d5918ca77179bbe5fab28c5d9d14e63a440338545e73961fa50b9d67",
        "f5_source_file_sha256": "143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be",
        "f5_propagation_fingerprint": "d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa",
        "f5_artifact_fingerprint": "261968ee7d9590d7d0afe0cd15155ef745a4169f95f52ffffd392aadd320de63",
        "farm_source_head": "6634e9cb470cf35f21f5d475ec6ce33b524cd233",
    },
}


class MethodAReweightingValidationError(ValueError):
    """Raised when detached F.6.1 authority or closure is unavailable."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodAReweightingValidationError("{}_not_json_safe".format(label)) from exc


def _map(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodAReweightingValidationError("{}_invalid".format(label))
    return value


def _seq(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodAReweightingValidationError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodAReweightingValidationError("{}_nonfinite".format(label))
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodAReweightingValidationError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(result):
        raise MethodAReweightingValidationError("{}_nonfinite".format(label))
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodAReweightingValidationError("{}_invalid".format(label))
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodAReweightingValidationError("{}_invalid".format(label)) from exc
    if result != value:
        raise MethodAReweightingValidationError("{}_invalid".format(label))
    return result


def _hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value.lower()):
        raise MethodAReweightingValidationError("{}_invalid".format(label))
    return value.lower()


def _close(left: float, right: float) -> bool:
    return abs(left - right) <= _TOLERANCE * max(1.0, abs(left), abs(right))


def _identity(row: Mapping[str, object], label: str) -> tuple[str, int]:
    source = row.get("source_label")
    if not isinstance(source, str) or not source:
        raise MethodAReweightingValidationError("{}_source_invalid".format(label))
    entry = _integer(row.get("entry_index"), "{}_entry_index".format(label))
    if entry < 0:
        raise MethodAReweightingValidationError("{}_entry_invalid".format(label))
    return source, entry


def _setting_id(setting: Mapping[str, object]) -> str:
    return "{}-{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"))


def _runtime_authority(
    kinematic: str,
    f1_hashes: Mapping[str, str],
    f3_artifact: Mapping[str, object],
    f3_file_sha256: str,
    f4_artifact: Mapping[str, object],
    f4_file_sha256: str,
    f5_artifact: Mapping[str, object],
    f5_file_sha256: str,
    authority_by_kinematic: Mapping[str, object] | None,
) -> dict[str, object]:
    records = ACCEPTED_RUNTIME_AUTHORITY_BY_KINEMATIC if authority_by_kinematic is None else _map(authority_by_kinematic, "f6_1_runtime_authority_records")
    accepted_raw = _map(records.get(kinematic), "f6_1_runtime_authority")
    f3_map = _map(f3_artifact.get("acceptance_map"), "f3_acceptance_map")
    f4_correction = _map(f4_artifact.get("correction"), "f4_correction")
    f5_propagation = _map(f5_artifact.get("propagation"), "f5_propagation")
    observed = {
        "f1_source_file_sha256": dict(f1_hashes),
        "f3_source_file_sha256": _hash(f3_file_sha256, "f3_source_file_sha256"),
        "f3_map_fingerprint": _hash(f3_map.get("fingerprint"), "f3_map_fingerprint"),
        "f3_algorithm_fingerprint": _hash(f3_map.get("algorithm_fingerprint"), "f3_algorithm_fingerprint"),
        "f3_artifact_fingerprint": _hash(f3_artifact.get("artifact_fingerprint"), "f3_artifact_fingerprint"),
        "f4_source_file_sha256": _hash(f4_file_sha256, "f4_source_file_sha256"),
        "f4_correction_fingerprint": _hash(f4_correction.get("fingerprint"), "f4_correction_fingerprint"),
        "f4_artifact_fingerprint": _hash(f4_artifact.get("artifact_fingerprint"), "f4_artifact_fingerprint"),
        "f5_source_file_sha256": _hash(f5_file_sha256, "f5_source_file_sha256"),
        "f5_propagation_fingerprint": _hash(f5_propagation.get("fingerprint"), "f5_propagation_fingerprint"),
        "f5_artifact_fingerprint": _hash(f5_artifact.get("artifact_fingerprint"), "f5_artifact_fingerprint"),
    }
    accepted: dict[str, object] = {}
    for name, value in observed.items():
        if name == "f1_source_file_sha256":
            supplied = _map(accepted_raw.get(name), "f6_1_runtime_authority_{}".format(name))
            normalized = {str(key): _hash(item, "f6_1_runtime_authority_f1_hash") for key, item in supplied.items()}
            if normalized != value:
                raise MethodAReweightingValidationError("f6_1_runtime_authority_{}_mismatch".format(name))
            accepted[name] = normalized
        else:
            normalized = _hash(accepted_raw.get(name), "f6_1_runtime_authority_{}".format(name))
            if normalized != value:
                raise MethodAReweightingValidationError("f6_1_runtime_authority_{}_mismatch".format(name))
            accepted[name] = normalized
    head = accepted_raw.get("farm_source_head")
    if not isinstance(head, str) or len(head) != 40 or any(char not in "0123456789abcdef" for char in head.lower()):
        raise MethodAReweightingValidationError("f6_1_runtime_authority_farm_source_head_invalid")
    accepted["farm_source_head"] = head.lower()
    return {"kinematic_token": kinematic, "accepted": accepted, "observed": observed, "accepted_authority_match": True}


def _validate_f5_artifact(value: object) -> tuple[Mapping[str, object], Mapping[str, object]]:
    artifact = _map(value, "f5_artifact")
    if artifact.get("schema_version") != _F5_ARTIFACT_SCHEMA:
        raise MethodAReweightingValidationError("f5_artifact_schema_version")
    for name, expected in (
        ("non_authoritative", True), ("f4_correction_consumed", True),
        ("f4_correction_modified", False), ("event_level_propagation_performed", True),
        ("event_correction_persisted", False), ("canonical_child_renormalization_performed", False),
        ("smoothing_or_interpolation_performed", False), ("yield_constructed", False),
        ("cross_section_constructed", False), ("root_object_constructed", False),
        ("production_application_performed", False), ("production_objects_mutated", False),
        ("method_b_numerical_dependency", False),
    ):
        if artifact.get(name) != expected:
            raise MethodAReweightingValidationError("f5_artifact_{}".format(name))
    propagation = _map(artifact.get("propagation"), "f5_propagation")
    if propagation.get("schema_version") != _F5_SCHEMA or propagation.get("available") is not True:
        raise MethodAReweightingValidationError("f5_propagation_invalid")
    fingerprint_inputs = _map(propagation.get("fingerprint_inputs"), "f5_fingerprint_inputs")
    if propagation.get("fingerprint") != _f5._sha256(fingerprint_inputs):
        raise MethodAReweightingValidationError("f5_propagation_fingerprint_mismatch")
    provenance = _map(artifact.get("provenance"), "f5_provenance")
    expected_artifact = _f5._sha256({"schema_version": _F5_ARTIFACT_SCHEMA, "propagation_fingerprint": propagation.get("fingerprint"), "input_paths": provenance.get("input_paths")})
    if artifact.get("artifact_fingerprint") != expected_artifact:
        raise MethodAReweightingValidationError("f5_artifact_fingerprint_mismatch")
    templates = _seq(propagation.get("setting_templates"), "f5_setting_templates")
    parents = _seq(propagation.get("parent_closure"), "f5_parent_closure")
    if len(templates) != 5 or len(parents) != 15:
        raise MethodAReweightingValidationError("f5_inventory_invalid")
    return artifact, propagation


def _raw_records(f1_artifacts: Sequence[Mapping[str, object]]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for raw_artifact in f1_artifacts:
        artifact = _map(raw_artifact, "f1_artifact")
        setting = _map(artifact.get("setting"), "f1_setting")
        setting_id = _setting_id(setting)
        if setting_id in result:
            raise MethodAReweightingValidationError("f1_setting_duplicate")
        contract = _map(artifact.get("contract"), "f1_contract")
        training: dict[tuple[str, int], Mapping[str, object]] = {}
        application: dict[tuple[str, int], Mapping[str, object]] = {}
        for index, raw in enumerate(_seq(contract.get("method_a_training_records"), "f1_training_records")):
            row = _map(raw, "f1_training_{}".format(index)); key = _identity(row, "f1_training_{}".format(index))
            if key in training:
                raise MethodAReweightingValidationError("f6_1_training_identity_duplicate")
            training[key] = row
        for index, raw in enumerate(_seq(contract.get("application_records"), "f1_application_records")):
            row = _map(raw, "f1_application_{}".format(index)); key = _identity(row, "f1_application_{}".format(index))
            if key in application:
                raise MethodAReweightingValidationError("f6_1_application_identity_duplicate")
            application[key] = row
        result[setting_id] = {"setting": setting, "contract": contract, "training": training, "application": application}
    expected = {"{}-{}".format(phi, epsilon) for phi, epsilon in _CANONICAL_SETTINGS}
    if set(result) != expected:
        raise MethodAReweightingValidationError("f1_setting_inventory_invalid")
    return result


def _phi_assignment(value: float, edges: np.ndarray, index: int) -> bool:
    return edges[index] <= value < edges[index + 1] or (index == len(edges) - 2 and value == edges[index + 1])


def _check_parity(training: Mapping[str, object], application: Mapping[str, object]) -> None:
    for name in ("t_index", "SHMS_delta", "P_hgcer_npeSum", "P_hgcer_xAtCer", "P_hgcer_yAtCer", "SHMS_xptar", "SHMS_yptar", "analysis_t", "analysis_MM", "Q2", "W", "epsilon"):
        if not _close(_finite(training.get(name), "f6_1_training_{}".format(name)), _finite(application.get(name), "f6_1_application_{}".format(name))):
            raise MethodAReweightingValidationError("f6_1_prompt_parity_mismatch:{}".format(name))
    for name in ("allcuts", "nommcuts"):
        if training.get(name) != application.get(name):
            raise MethodAReweightingValidationError("f6_1_prompt_parity_mismatch:{}".format(name))


def _required_training_value(row: Mapping[str, object], name: str) -> float:
    return _finite(row.get(name), "f6_1_required_training_variable_missing:{}".format(name))


def _display_edges(values: Sequence[float], bins: int, label: str) -> list[float]:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise MethodAReweightingValidationError("{}_values_invalid".format(label))
    low, high = float(np.min(array)), float(np.max(array))
    if low == high:
        half = max(1.0, abs(low)) * 1.0e-6
        low, high = low - half, high + half
    return [float(value) for value in np.linspace(low, high, bins + 1)]


def _normalized_histogram(values: Sequence[float], weights: Sequence[float], edges: Sequence[float], label: str) -> list[float]:
    histogram = np.asarray(np.histogram(np.asarray(values, dtype=float), bins=np.asarray(edges, dtype=float), weights=np.asarray(weights, dtype=float))[0], dtype=float)
    total = math.fsum(float(value) for value in histogram)
    if not math.isfinite(total) or total <= 0.0:
        raise MethodAReweightingValidationError("{}_normalization_invalid".format(label))
    result = histogram / total
    if not np.all(np.isfinite(result)):
        raise MethodAReweightingValidationError("{}_normalization_invalid".format(label))
    return [float(value) for value in result]


def _hellinger(left: Sequence[float], right: Sequence[float]) -> float:
    result = float(math.sqrt(0.5 * math.fsum((math.sqrt(max(0.0, a)) - math.sqrt(max(0.0, b))) ** 2 for a, b in zip(left, right))))
    if not math.isfinite(result):
        raise MethodAReweightingValidationError("f6_1_hellinger_invalid")
    return result


def _value(row: Mapping[str, object], variable: str, *, application: bool) -> float:
    name = "phi_degrees" if application and variable == "phi" else variable
    return _finite(row.get(name), "f6_1_required_{}_variable_missing:{}".format("application" if application else "training", variable))


def _shape_comparisons(low_rows: Sequence[Mapping[str, object]], prompt_pairs: Sequence[Mapping[str, object]], phi_edges: Sequence[float]) -> tuple[dict[str, object], dict[str, object]]:
    comparisons: dict[str, object] = {}
    for variable in _MODEL_VARIABLES + _INDEPENDENT_VARIABLES:
        low_values = [_value(row, variable, application=False) for row in low_rows]
        control_values = [_value(_map(pair["row"], "f6_1_pair_row"), variable, application=True) for pair in prompt_pairs]
        edges = list(phi_edges) if variable == "phi" else _display_edges(low_values + control_values, int(DISPLAY_POLICY["one_dimensional_nonphi_bin_count"]), "f6_1_{}".format(variable))
        baseline_weights = [_finite(_map(pair["row"], "f6_1_pair_row").get("baseline_pion_weight_w0"), "f6_1_w0") for pair in prompt_pairs]
        adjusted_weights = [weight * _finite(pair.get("factor"), "f6_1_correction") for weight, pair in zip(baseline_weights, prompt_pairs)]
        low = _normalized_histogram(low_values, [1.0] * len(low_values), edges, "f6_1_low_{}".format(variable))
        baseline = _normalized_histogram(control_values, baseline_weights, edges, "f6_1_baseline_{}".format(variable))
        method_a = _normalized_histogram(control_values, adjusted_weights, edges, "f6_1_method_a_{}".format(variable))
        baseline_h = _hellinger(low, baseline); method_h = _hellinger(low, method_a)
        comparisons[variable] = {
            "edges": edges, "low_response_unit_area": low, "baseline_prompt_control_unit_area": baseline,
            "method_a_prompt_control_unit_area": method_a, "hellinger_baseline": baseline_h,
            "hellinger_method_a": method_h, "delta_hellinger": method_h - baseline_h,
        }
    x_low = [_value(row, "P_hgcer_xAtCer", application=False) for row in low_rows]
    y_low = [_value(row, "P_hgcer_yAtCer", application=False) for row in low_rows]
    x_control = [_value(_map(pair["row"], "f6_1_pair_row"), "P_hgcer_xAtCer", application=True) for pair in prompt_pairs]
    y_control = [_value(_map(pair["row"], "f6_1_pair_row"), "P_hgcer_yAtCer", application=True) for pair in prompt_pairs]
    x_edges = _display_edges(x_low + x_control, int(DISPLAY_POLICY["hgcer_xy_bin_count_per_axis"]), "f6_1_hgcer_x")
    y_edges = _display_edges(y_low + y_control, int(DISPLAY_POLICY["hgcer_xy_bin_count_per_axis"]), "f6_1_hgcer_y")
    weights = [_finite(_map(pair["row"], "f6_1_pair_row").get("baseline_pion_weight_w0"), "f6_1_w0") for pair in prompt_pairs]
    adjusted = [weight * _finite(pair.get("factor"), "f6_1_correction") for weight, pair in zip(weights, prompt_pairs)]
    matrices: dict[str, list[list[float]]] = {}
    for name, xs, ys, selected_weights in (("low_response_unit_area", x_low, y_low, [1.0] * len(x_low)), ("baseline_prompt_control_unit_area", x_control, y_control, weights), ("method_a_prompt_control_unit_area", x_control, y_control, adjusted)):
        matrix = np.asarray(np.histogram2d(xs, ys, bins=(x_edges, y_edges), weights=selected_weights)[0], dtype=float)
        total = math.fsum(float(value) for value in matrix.flat)
        if not math.isfinite(total) or total <= 0.0:
            raise MethodAReweightingValidationError("f6_1_hgcer_xy_normalization_invalid")
        matrices[name] = [[float(value) for value in row] for row in matrix / total]
    return comparisons, {"x_edges": x_edges, "y_edges": y_edges, **matrices}


def _signed_histogram(values: Sequence[float], baseline: Sequence[float], adjusted: Sequence[float], edges: Sequence[float]) -> dict[str, list[float]]:
    base = np.asarray(np.histogram(np.asarray(values, dtype=float), bins=np.asarray(edges, dtype=float), weights=np.asarray(baseline, dtype=float))[0], dtype=float)
    shifted = np.asarray(np.histogram(np.asarray(values, dtype=float), bins=np.asarray(edges, dtype=float), weights=np.asarray(adjusted, dtype=float))[0], dtype=float)
    return {"edges": [float(value) for value in edges], "baseline_signed_contents": [float(value) for value in base], "method_a_signed_contents": [float(value) for value in shifted], "signed_delta_contents": [float(value) for value in shifted - base]}


def _same(left: object, right: object) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left == right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return _close(float(left), float(right))
    if isinstance(left, Sequence) and not isinstance(left, (str, bytes, bytearray)) and isinstance(right, Sequence) and not isinstance(right, (str, bytes, bytearray)):
        return len(left) == len(right) and all(_same(a, b) for a, b in zip(left, right))
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return set(left) == set(right) and all(_same(left[key], right[key]) for key in left)
    return left == right


def _independent_tphi(
    pairs_by_parent: Mapping[tuple[str, int], Sequence[Mapping[str, object]]],
    persisted_f4: Mapping[str, object],
    raw: Mapping[str, Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    f4_parents = {(str(item["setting_id"]), int(item["canonical_t_index"])): _map(item, "f4_parent") for item in _seq(persisted_f4.get("parents"), "f4_parents")}
    settings: list[dict[str, object]] = []; closures: list[dict[str, object]] = []
    for phi, epsilon in _CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); state = raw[setting_id]; contract = _map(state["contract"], "f1_contract")
        phi_edges = [_finite(value, "f6_1_phi_edge") for value in _seq(contract.get("phi_edges"), "f6_1_phi_edges")]
        t_edges = [_finite(value, "f6_1_t_edge") for value in _seq(contract.get("t_edges"), "f6_1_t_edges")]
        counts = [[0 for _ in range(9)] for _ in range(3)]
        base_values: list[list[list[float]]] = [[[ ] for _ in range(9)] for _ in range(3)]
        adjusted_values: list[list[list[float]]] = [[[ ] for _ in range(9)] for _ in range(3)]
        parent_rows: list[dict[str, object]] = []
        for t_index in range(3):
            parent = f4_parents.get((setting_id, t_index)); pairs = pairs_by_parent.get((setting_id, t_index))
            if parent is None or pairs is None:
                raise MethodAReweightingValidationError("f6_1_tphi_parent_missing")
            for pair in pairs:
                row = _map(pair["row"], "f6_1_pair_row"); phi_index = _integer(row.get("phi_index"), "f6_1_phi_index")
                if phi_index < 0 or phi_index >= 9:
                    raise MethodAReweightingValidationError("f6_1_phi_index_invalid")
                baseline = _finite(row.get("signed_baseline_event_contribution"), "f6_1_signed_baseline")
                value = baseline * _finite(pair.get("factor"), "f6_1_correction")
                counts[t_index][phi_index] += 1; base_values[t_index][phi_index].append(baseline); adjusted_values[t_index][phi_index].append(value)
            baseline_row = [math.fsum(values) for values in base_values[t_index]]
            adjusted_row = [math.fsum(values) for values in adjusted_values[t_index]]
            delta_row = [adjusted_row[index] - baseline_row[index] for index in range(9)]
            parent_baseline = math.fsum(baseline_row); parent_adjusted = math.fsum(adjusted_row)
            if not _close(parent_baseline, _finite(parent.get("baseline_parent_sum"), "f4_parent_baseline")) or not _close(parent_adjusted, _finite(parent.get("adjusted_parent_sum"), "f4_parent_adjusted")):
                raise MethodAReweightingValidationError("f6_1_f4_parent_reproduction_mismatch")
            redistribution = [value / parent_baseline for value in delta_row]
            parent_rows.append({"canonical_t_index": t_index, "baseline_parent_sum": parent_baseline, "adjusted_parent_sum": parent_adjusted, "closure_residual": parent_adjusted - parent_baseline, "max_abs_child_delta": max(abs(value) for value in delta_row), "sum_abs_child_delta": math.fsum(abs(value) for value in delta_row), "max_abs_redistribution": max(abs(value) for value in redistribution), "sum_abs_redistribution": math.fsum(abs(value) for value in redistribution), "max_abs_redistribution_phi_index": max(range(9), key=lambda index: abs(redistribution[index]))})
        baseline = [[math.fsum(values) for values in row] for row in base_values]
        adjusted = [[math.fsum(values) for values in row] for row in adjusted_values]
        delta = [[adjusted[t][p] - baseline[t][p] for p in range(9)] for t in range(3)]
        setting = {"setting": _copy(state["setting"], "setting"), "setting_id": setting_id, "t_edges": t_edges, "phi_edges": phi_edges, "event_counts": counts, "baseline_signed_contents": baseline, "adjusted_signed_contents": adjusted, "signed_delta_contents": delta, "baseline_share_of_parent": [[baseline[t][p] / parent_rows[t]["baseline_parent_sum"] for p in range(9)] for t in range(3)], "adjusted_share_of_parent": [[adjusted[t][p] / parent_rows[t]["adjusted_parent_sum"] for p in range(9)] for t in range(3)], "redistribution_fraction_of_parent": [[delta[t][p] / parent_rows[t]["baseline_parent_sum"] for p in range(9)] for t in range(3)], "parents": parent_rows, "setting_baseline_sum": math.fsum(row["baseline_parent_sum"] for row in parent_rows), "setting_adjusted_sum": math.fsum(row["adjusted_parent_sum"] for row in parent_rows), "setting_closure_residual": math.fsum(row["closure_residual"] for row in parent_rows), "setting_max_abs_redistribution": max(row["max_abs_redistribution"] for row in parent_rows)}
        settings.append(setting); closures.extend({"setting_id": setting_id, **row} for row in parent_rows)
    return settings, closures


def _verify_f5_continuity(
    pairs_by_parent: Mapping[tuple[str, int], Sequence[Mapping[str, object]]],
    persisted_f4: Mapping[str, object],
    raw: Mapping[str, Mapping[str, object]],
    persisted_f5: Mapping[str, object],
) -> dict[str, object]:
    settings, closures = _independent_tphi(pairs_by_parent, persisted_f4, raw)
    observed = {str(item["setting_id"]): item for item in settings}
    expected = {str(item["setting_id"]): _map(item, "f5_setting") for item in _seq(persisted_f5.get("setting_templates"), "f5_setting_templates")}
    if set(observed) != set(expected):
        raise MethodAReweightingValidationError("f6_1_f5_setting_inventory_mismatch")
    names = ("event_counts", "baseline_signed_contents", "adjusted_signed_contents", "signed_delta_contents", "baseline_share_of_parent", "adjusted_share_of_parent", "redistribution_fraction_of_parent", "setting_baseline_sum", "setting_adjusted_sum", "setting_closure_residual")
    for setting_id, item in observed.items():
        for name in names:
            if not _same(item[name], expected[setting_id].get(name)):
                raise MethodAReweightingValidationError("f6_1_f5_aggregate_mismatch:{}".format(name))
    expected_closure = {(str(item["setting_id"]), int(item["canonical_t_index"])): _map(item, "f5_parent") for item in _seq(persisted_f5.get("parent_closure"), "f5_parent_closure")}
    for row in closures:
        expected_row = expected_closure.get((str(row["setting_id"]), int(row["canonical_t_index"])))
        if expected_row is None or not _same(row, expected_row):
            raise MethodAReweightingValidationError("f6_1_f5_parent_closure_mismatch")
    return {"checked": True, "tolerance": _TOLERANCE, "setting_count": 5, "parent_count": 15, "cell_count": 135, "settings": settings, "parent_closure": closures}


def _validate_and_pair_rows(
    f1_artifacts: Sequence[Mapping[str, object]],
    f3_artifact: Mapping[str, object],
    f1_hashes: Mapping[str, str],
    f3_file_sha256: str,
    persisted_f4: Mapping[str, object],
    accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None,
) -> tuple[dict[tuple[str, int], list[dict[str, object]]], dict[tuple[str, int], dict[str, object]], Mapping[str, object], Mapping[str, object]]:
    try:
        recomputed_f4, review = _f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(f1_artifacts, f3_artifact, f1_input_file_hashes=f1_hashes, f3_input_file_sha256=f3_file_sha256, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
        parsed = _f4._f3._validate_f1_artifacts(f1_artifacts)
        rows_by_parent = _f4._raw_application_rows(f1_artifacts, parsed, _TOLERANCE)
    except (ValueError, _f4.MethodAParentPreservingCorrectionError) as exc:
        raise MethodAReweightingValidationError("f6_1_f4_shared_reproduction_failed:{}".format(exc)) from exc
    if recomputed_f4 != persisted_f4 or recomputed_f4.get("fingerprint") != persisted_f4.get("fingerprint"):
        raise MethodAReweightingValidationError("f6_1_f4_shared_reproduction_mismatch")
    review_by_parent = {(str(item.get("setting_id")), _integer(item.get("canonical_t_index"), "f6_1_review_t_index")): _map(item, "f6_1_review") for item in review}
    persisted_parents = {(str(item["setting_id"]), int(item["canonical_t_index"])): _map(item, "f4_parent") for item in _seq(persisted_f4.get("parents"), "f4_parents")}
    if len(review_by_parent) != 15 or set(review_by_parent) != set(persisted_parents) or set(rows_by_parent) != set(persisted_parents):
        raise MethodAReweightingValidationError("f6_1_f4_review_inventory_invalid")
    raw = _raw_records(f1_artifacts); pairs_by_parent: dict[tuple[str, int], list[dict[str, object]]] = {}
    for key, rows in rows_by_parent.items():
        review_row = review_by_parent[key]; factors = np.asarray(review_row.get("correction_factors"), dtype=float); support = np.asarray(review_row.get("in_support_mask"), dtype=bool)
        if factors.ndim != 1 or support.ndim != 1 or factors.size != len(rows) or support.size != len(rows) or not np.all(np.isfinite(factors)) or np.any(factors <= 0.0):
            raise MethodAReweightingValidationError("f6_1_f4_review_factor_alignment_invalid")
        application = _map(raw[key[0]]["application"], "f6_1_raw_application")
        paired: list[dict[str, object]] = []
        for row, factor, mask in zip(rows, factors, support):
            identity = _identity(row, "f6_1_application")
            raw_row = application.get(identity)
            if raw_row is None:
                raise MethodAReweightingValidationError("f6_1_application_raw_identity_missing")
            combined = {**row, **_map(raw_row, "f6_1_raw_application_row")}
            paired.append({"row": combined, "factor": float(factor), "in_support": bool(mask)})
        paired.sort(key=lambda item: _identity(_map(item["row"], "f6_1_pair_row"), "f6_1_pair"))
        identities = [_identity(_map(item["row"], "f6_1_pair_row"), "f6_1_pair") for item in paired]
        if len(identities) != len(set(identities)):
            raise MethodAReweightingValidationError("f6_1_application_identity_duplicate")
        pairs_by_parent[key] = paired
    return pairs_by_parent, {key: dict(value) for key, value in persisted_parents.items()}, recomputed_f4, raw


def _parent_payloads(
    pairs_by_parent: Mapping[tuple[str, Sequence[Mapping[str, object]]], Sequence[Mapping[str, object]]],
    raw: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    payloads: list[dict[str, object]] = []
    for phi, epsilon in _CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); state = raw[setting_id]; training = _map(state["training"], "f6_1_training")
        phi_edges = [_finite(value, "f6_1_phi_edge") for value in _seq(_map(state["contract"], "f6_1_contract").get("phi_edges"), "f6_1_phi_edges")]
        if len(phi_edges) != 10:
            raise MethodAReweightingValidationError("f6_1_canonical_phi_geometry_invalid")
        for t_index in range(3):
            pairs = pairs_by_parent.get((setting_id, t_index))
            if pairs is None:
                raise MethodAReweightingValidationError("f6_1_parent_pairs_missing")
            controls: list[Mapping[str, object]] = []; low: list[Mapping[str, object]] = []
            for key, row in training.items():
                if _integer(row.get("t_index"), "f6_1_training_t_index") != t_index or key[0] != "prompt":
                    continue
                response = row.get("response_class"); npe = _finite(row.get("P_hgcer_npeSum"), "f6_1_training_npe")
                if response == "low":
                    if row.get("nommcuts") is not True or not (0.0 < npe <= 2.0):
                        raise MethodAReweightingValidationError("f6_1_low_population_invalid")
                    for variable in _REQUIRED_TRAINING_VARIABLES:
                        _required_training_value(row, variable)
                    low.append(row)
            prompt_pairs: list[Mapping[str, object]] = []
            for pair in pairs:
                row = _map(pair["row"], "f6_1_pair_row")
                if str(row.get("source_label")) != "prompt":
                    continue
                identity = _identity(row, "f6_1_prompt_application")
                matched = training.get(identity)
                if matched is None or matched.get("response_class") != "control":
                    raise MethodAReweightingValidationError("f6_1_prompt_application_identity_missing_training_control")
                _check_parity(matched, row)
                training_phi = _finite(matched.get("phi"), "f6_1_training_phi")
                application_phi = _finite(row.get("phi_degrees"), "f6_1_application_phi_degrees")
                index = _integer(row.get("phi_index"), "f6_1_application_phi_index")
                if row.get("phi_status") != "inside_phi" or index < 0 or index >= 9 or not _close(training_phi, application_phi) or not _phi_assignment(application_phi, np.asarray(phi_edges, dtype=float), index):
                    raise MethodAReweightingValidationError("f6_1_phi_semantic_parity_mismatch")
                for variable in _REQUIRED_TRAINING_VARIABLES:
                    _required_training_value(matched, variable); _value(row, variable, application=True)
                controls.append(matched); prompt_pairs.append(pair)
            for row in low:
                value = _finite(row.get("phi"), "f6_1_low_phi")
                if not (phi_edges[0] <= value <= phi_edges[-1]):
                    raise MethodAReweightingValidationError("f6_1_low_phi_outside_canonical_domain")
            if not low or not prompt_pairs:
                raise MethodAReweightingValidationError("f6_1_prompt_population_empty")
            comparisons, detector_xy = _shape_comparisons(low, prompt_pairs, phi_edges)
            full_values = [_finite(_map(pair["row"], "f6_1_pair_row").get("analysis_MM"), "f6_1_full_analysis_MM") for pair in pairs]
            baseline = [_finite(_map(pair["row"], "f6_1_pair_row").get("signed_baseline_event_contribution"), "f6_1_full_baseline") for pair in pairs]
            adjusted = [value * _finite(pair.get("factor"), "f6_1_correction") for value, pair in zip(baseline, pairs)]
            phi_values = [_finite(_map(pair["row"], "f6_1_pair_row").get("phi_degrees"), "f6_1_full_phi") for pair in pairs]
            source_counts = {name: 0 for name in ("prompt", "rand", "dummy", "dummy_rand")}
            for pair in pairs:
                label = str(_map(pair["row"], "f6_1_pair_row").get("source_label"))
                if label not in source_counts:
                    raise MethodAReweightingValidationError("f6_1_application_source_invalid")
                source_counts[label] += 1
            payloads.append({"setting": _copy(state["setting"], "setting"), "setting_id": setting_id, "canonical_t_index": t_index, "population_counts": {"low_response_prompt_training": len(low), "prompt_physical_control": len(prompt_pairs), "full_physical_application": len(pairs), "training_control_not_in_physical_prompt_application": len(controls) - len(prompt_pairs), "full_application_by_source": source_counts}, "identity_and_phi_parity": {"prompt_control_count": len(prompt_pairs), "all_prompt_controls_matched": True, "all_common_fields_matched": True, "phi_semantics_closed_without_conversion": True}, "support": {"in_support_count": sum(bool(pair["in_support"]) for pair in pairs), "ood_count": sum(not bool(pair["in_support"]) for pair in pairs)}, "prompt_shape_comparisons": comparisons, "hgcer_xy": detector_xy, "signed_background": {"analysis_MM": _signed_histogram(full_values, baseline, adjusted, _display_edges(full_values, int(DISPLAY_POLICY["one_dimensional_nonphi_bin_count"]), "f6_1_signed_mm")), "phi": _signed_histogram(phi_values, baseline, adjusted, phi_edges)}})
    return payloads


def build_pion_hgcer_method_a_reweighting_validation(
    f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], f5_artifact: Mapping[str, object], *,
    f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, f5_input_file_sha256: str,
    accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
) -> dict[str, object]:
    try:
        parsed = _f4._f3._validate_f1_artifacts(f1_artifacts)
    except ValueError as exc:
        raise MethodAReweightingValidationError("f6_1_f1_validation_failed:{}".format(exc)) from exc
    f1_hashes = {str(key): _hash(value, "f1_input_file_sha256") for key, value in _map(f1_input_file_hashes, "f1_input_file_hashes").items()}
    if set(f1_hashes) != {str(item["setting_id"]) for item in parsed}:
        raise MethodAReweightingValidationError("f6_1_f1_input_hash_inventory_invalid")
    artifact_f4, persisted_f4, _f4_authority = _f5._validate_f4_artifact(f4_artifact, _hash(f4_input_file_sha256, "f4_input_file_sha256"), accepted_f4_runtime_authority_by_kinematic)
    artifact_f5, persisted_f5 = _validate_f5_artifact(f5_artifact)
    setting = _map(parsed[0]["setting"], "f1_setting"); kinematic = str(setting.get("kinematic_token"))
    authority = _runtime_authority(kinematic, f1_hashes, f3_artifact, _hash(f3_input_file_sha256, "f3_input_file_sha256"), artifact_f4, _hash(f4_input_file_sha256, "f4_input_file_sha256"), artifact_f5, _hash(f5_input_file_sha256, "f5_input_file_sha256"), accepted_runtime_authority_by_kinematic)
    pairs, _parents, recomputed_f4, raw = _validate_and_pair_rows(f1_artifacts, f3_artifact, f1_hashes, _hash(f3_input_file_sha256, "f3_input_file_sha256"), persisted_f4, accepted_f3_runtime_authority_by_kinematic)
    try:
        recomputed_f5 = _f5.build_pion_hgcer_method_a_tphi_propagation(f1_artifacts, f3_artifact, artifact_f4, f1_input_file_hashes=f1_hashes, f3_input_file_sha256=_hash(f3_input_file_sha256, "f3_input_file_sha256"), f4_input_file_sha256=_hash(f4_input_file_sha256, "f4_input_file_sha256"), accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    except (ValueError, _f5.MethodATPhiPropagationError) as exc:
        raise MethodAReweightingValidationError("f6_1_f5_shared_reproduction_failed:{}".format(exc)) from exc
    if recomputed_f5 != persisted_f5:
        raise MethodAReweightingValidationError("f6_1_f5_shared_reproduction_mismatch")
    if recomputed_f5.get("fingerprint") != authority["accepted"]["f5_propagation_fingerprint"]:
        raise MethodAReweightingValidationError("f6_1_f5_fingerprint_mismatch")
    parent_payloads = _parent_payloads(pairs, raw)
    continuity = _verify_f5_continuity(pairs, persisted_f4, raw, persisted_f5)
    inputs = [{"setting_id": str(item["setting_id"]), "source_file_sha256": f1_hashes[str(item["setting_id"])], "stable_f1_content_fingerprint": item["stable_content_fingerprint"]} for item in parsed]
    core: dict[str, object] = {"schema_version": METHOD_A_REWEIGHTING_VALIDATION_SCHEMA_VERSION, "fingerprint_schema_version": METHOD_A_REWEIGHTING_VALIDATION_FINGERPRINT_SCHEMA_VERSION, "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete", "non_authoritative": True, "validation_only": True, "f4_correction_consumed": True, "f4_correction_modified": False, "f5_continuity_checked": True, "event_correction_evaluated_for_detached_validation": True, "event_correction_persisted": False, "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "method_b_numerical_dependency": False, "shape_improvement_gate_applied": False, "manual_review_required": True, "display_policy": _copy(DISPLAY_POLICY, "display_policy"), "runtime_authority": authority, "input_fingerprints": inputs, "f4_reproduction": {"exact_payload_match": True, "correction_fingerprint": recomputed_f4["fingerprint"]}, "f5_reproduction": {"exact_payload_match": True, "propagation_fingerprint": recomputed_f5["fingerprint"]}, "parents": parent_payloads, "f5_continuity": continuity}
    fingerprint_inputs = {"schema_version": core["schema_version"], "fingerprint_schema_version": core["fingerprint_schema_version"], "display_policy": core["display_policy"], "runtime_authority": core["runtime_authority"], "input_fingerprints": core["input_fingerprints"], "f4_reproduction": core["f4_reproduction"], "f5_reproduction": core["f5_reproduction"], "parents": core["parents"], "f5_continuity": core["f5_continuity"]}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    serialized = _canonical_json(core)
    for forbidden in ("correction_factors", "raw_shape_factors", "in_support_mask", "event_corrections", "entry_index", "application_records", "method_a_training_records"):
        if forbidden in serialized:
            raise MethodAReweightingValidationError("f6_1_forbidden_event_persistence:{}".format(forbidden))
    return _copy(core, "reweighting_validation")  # type: ignore[return-value]


def build_pion_hgcer_method_a_reweighting_validation_artifact(
    f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], f5_artifact: Mapping[str, object], *,
    f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, f5_input_file_sha256: str,
    input_paths: Mapping[str, object] | None = None, generated_at_utc: str | None = None, git_head: str | None = None, git_status_short: str | None = None,
    accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
) -> dict[str, object]:
    validation = build_pion_hgcer_method_a_reweighting_validation(f1_artifacts, f3_artifact, f4_artifact, f5_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, f4_input_file_sha256=f4_input_file_sha256, f5_input_file_sha256=f5_input_file_sha256, accepted_runtime_authority_by_kinematic=accepted_runtime_authority_by_kinematic, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    provenance = {"generated_at_utc": generated_at_utc, "git_head": git_head, "git_status_short": git_status_short, "input_paths": _copy({} if input_paths is None else input_paths, "input_paths")}
    artifact = {"schema_version": METHOD_A_REWEIGHTING_VALIDATION_ARTIFACT_SCHEMA_VERSION, "validation": validation, "non_authoritative": True, "validation_only": True, "f4_correction_consumed": True, "f4_correction_modified": False, "f5_continuity_checked": True, "event_correction_evaluated_for_detached_validation": True, "event_correction_persisted": False, "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "method_b_numerical_dependency": False, "shape_improvement_gate_applied": False, "manual_review_required": True, "provenance": provenance}
    artifact["artifact_fingerprint"] = _sha256({"schema_version": artifact["schema_version"], "validation_fingerprint": validation["fingerprint"], "input_paths": provenance["input_paths"]})
    return _copy(artifact, "reweighting_validation_artifact")  # type: ignore[return-value]


def pion_hgcer_method_a_reweighting_validation_filename(kinematic: object) -> str:
    token = str(kinematic)
    if not token or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise MethodAReweightingValidationError("reweighting_validation_filename_token_invalid")
    return "{}_kaon_pion-background_hgcer_method-a-reweighting-validation.json".format(token)


def write_pion_hgcer_method_a_reweighting_validation_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    with open(os.fspath(path), "w", encoding="utf-8") as handle:
        json.dump(_copy(payload, "reweighting_validation_artifact"), handle, sort_keys=True, indent=2, allow_nan=False); handle.write("\n")
    return os.fspath(path)
