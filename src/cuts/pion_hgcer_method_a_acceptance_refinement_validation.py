"""Detached F.6.2 acceptance-correlated Method-A refinement validation.

The module consumes the frozen F.1--F.6.1 artifacts, evaluates F.4 factors
only transiently, and persists aggregate distributions and diagnostics only.
It intentionally has no normal-analysis, ROOT, yield, or Method-B ownership.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

import pion_hgcer_method_a_parent_preserving_correction as _f4
import pion_hgcer_method_a_reweighting_validation as _f61


_UTILITY_DIRECTORY = Path(__file__).resolve().parents[1] / "utility"
if str(_UTILITY_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(_UTILITY_DIRECTORY))
from background_config import resolve_analysis_runtime_config  # noqa: E402


METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_refinement_validation/v1"
)
METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_FINGERPRINT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_refinement_validation_fingerprint/v1"
)
METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_ARTIFACT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_refinement_validation_artifact/v1"
)
_F6_1_ARTIFACT_SCHEMA = "pion_hgcer_method_a_reweighting_validation_artifact/v1"
_F6_1_SCHEMA = "pion_hgcer_method_a_reweighting_validation/v1"
_KINEMATIC = "Q4p4W2p74"
_CANONICAL_SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)
_ONE_DIMENSIONAL_VARIABLES = (
    "SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer",
    "SHMS_xptar", "SHMS_yptar", "analysis_MM",
)
_MODEL_VARIABLES = _ONE_DIMENSIONAL_VARIABLES[:3]
_INDEPENDENT_VARIABLES = _ONE_DIMENSIONAL_VARIABLES[3:5]
_JOINT_VARIABLES = (
    ("analysis_MM", "SHMS_xptar"), ("analysis_MM", "SHMS_yptar"),
    ("SHMS_delta", "SHMS_xptar"), ("SHMS_delta", "SHMS_yptar"),
)
_SOURCE_LABELS = ("prompt", "rand", "dummy", "dummy_rand")
_TOLERANCE = 1.0e-12
_BOOTSTRAP_POLICY = {
    "policy_version": "F6.2_BOOTSTRAP/v1",
    "replicas": 2000,
    "global_seed": 20260916,
    "confidence_interval": {"method": "linear_percentile/v1", "low_percent": 2.5, "high_percent": 97.5},
    "control_resampling": "paired_baseline_and_method_a/v1",
    "low_response_resampling": "independent_of_control/v1",
    "signed_resampling": "fixed_count_per_source_stratum_paired_baseline_and_method_a/v1",
    "minimum_valid_fraction_for_kappa_and_rho": 0.90,
}


# Source-owned authority.  Analyzer callers use this exact record; focused
# tests may inject a self-consistent synthetic record into the builder.
ACCEPTED_F6_1_ARTIFACT_AUTHORITY_BY_KINEMATIC = {
    _KINEMATIC: {
        "source_file_sha256": "62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6",
        "validation_fingerprint": "d992d789b3434897d76691df27c190a0f51479f67c0d5b496e84835e517c77f8",
        "artifact_fingerprint": "377872a218a780481347402e4410c81bd568a4fb7682cd49f0f3352b499bad41",
    }
}


class MethodAAcceptanceRefinementValidationError(ValueError):
    """The frozen validation chain or detached diagnostic is invalid."""


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _copy(value: object, label: str) -> object:
    try:
        return json.loads(_canonical_json(value))
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRefinementValidationError("{}_not_json_safe".format(label)) from exc


def _map(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    return value


def _seq(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise MethodAAcceptanceRefinementValidationError("{}_nonfinite".format(label))
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRefinementValidationError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(result):
        raise MethodAAcceptanceRefinementValidationError("{}_nonfinite".format(label))
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label)) from exc
    if result != value:
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    return result


def _hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value.lower()):
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    return value.lower()


def _identity(row: Mapping[str, object], label: str) -> tuple[str, int]:
    source = row.get("source_label")
    if not isinstance(source, str) or source not in _SOURCE_LABELS:
        raise MethodAAcceptanceRefinementValidationError("{}_source_label_invalid".format(label))
    entry = _integer(row.get("entry_index"), "{}_entry_index".format(label))
    if entry < 0:
        raise MethodAAcceptanceRefinementValidationError("{}_entry_index_invalid".format(label))
    return source, entry


def _strict_edges(value: object, label: str, *, expected_length: int | None = None) -> list[float]:
    edges = [_finite(item, "{}_edge".format(label)) for item in _seq(value, label)]
    if expected_length is not None and len(edges) != expected_length:
        raise MethodAAcceptanceRefinementValidationError("{}_length_invalid".format(label))
    if len(edges) < 2 or any(right <= left for left, right in zip(edges, edges[1:])):
        raise MethodAAcceptanceRefinementValidationError("{}_not_strictly_increasing".format(label))
    return edges


def _metric(value: float | None = None, reason: str | None = None) -> dict[str, object]:
    if value is None:
        return {"available": False, "value": None, "reason": reason or "unavailable"}
    if not math.isfinite(value):
        return {"available": False, "value": None, "reason": "nonfinite_input"}
    return {"available": True, "value": float(value), "reason": None}


def _population(contents: np.ndarray, reason: str | None = None) -> dict[str, object]:
    total = float(np.sum(contents))
    if not np.all(np.isfinite(contents)):
        return {"available": False, "unit_area": [0.0] * int(contents.size), "reason": "nonfinite_input"}
    if reason is not None:
        return {"available": False, "unit_area": [0.0] * int(contents.size), "reason": reason}
    if not math.isfinite(total) or total <= 0.0:
        return {"available": False, "unit_area": [0.0] * int(contents.size), "reason": "normalization_invalid"}
    return {"available": True, "unit_area": [float(item) for item in (contents / total)], "reason": None}


def _population_array(population: Mapping[str, object]) -> np.ndarray | None:
    if population.get("available") is not True:
        return None
    values = np.asarray(population.get("unit_area"), dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or not math.isclose(float(np.sum(values)), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        return None
    return values


def _population_matrix(population: Mapping[str, object]) -> np.ndarray | None:
    if population.get("available") is not True:
        return None
    values = np.asarray(population.get("unit_area"), dtype=float)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or not math.isclose(float(np.sum(values)), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        return None
    return values


def _unavailable_reason(low: Mapping[str, object], baseline: Mapping[str, object], method_a: Mapping[str, object]) -> str:
    for candidate in (low, baseline, method_a):
        reason = candidate.get("reason")
        if isinstance(reason, str) and reason:
            return reason
    return "normalization_invalid"


def _hellinger(left: np.ndarray, right: np.ndarray) -> float:
    return float(math.sqrt(0.5 * float(np.sum((np.sqrt(left) - np.sqrt(right)) ** 2))))


def _total_variation(left: np.ndarray, right: np.ndarray) -> float:
    return float(0.5 * float(np.sum(np.abs(left - right))))


def _alignment_metrics(low: np.ndarray | None, baseline: np.ndarray | None, method_a: np.ndarray | None) -> dict[str, object]:
    if low is None or baseline is None or method_a is None:
        reason = _unavailable_reason(
            {"reason": "low_population_empty" if low is None else None},
            {"reason": "control_population_empty" if baseline is None else None},
            {"reason": "control_population_empty" if method_a is None else None},
        )
        return {"kappa": _metric(reason=reason), "rho": _metric(reason=reason)}
    residual, correction = low - baseline, method_a - baseline
    residual_norm = float(np.linalg.norm(residual.ravel()))
    correction_norm = float(np.linalg.norm(correction.ravel()))
    if not math.isfinite(residual_norm) or not math.isfinite(correction_norm):
        return {"kappa": _metric(reason="nonfinite_input"), "rho": _metric(reason="nonfinite_input")}
    if residual_norm == 0.0:
        return {"kappa": _metric(reason="residual_norm_zero"), "rho": _metric(reason="residual_norm_zero")}
    rho = _metric(correction_norm / residual_norm)
    if correction_norm == 0.0:
        return {"kappa": _metric(reason="correction_norm_zero"), "rho": rho}
    return {
        "kappa": _metric(float(np.dot(residual.ravel(), correction.ravel()) / (residual_norm * correction_norm))),
        "rho": rho,
    }


def _comparison_metrics(low: Mapping[str, object], baseline: Mapping[str, object], method_a: Mapping[str, object]) -> dict[str, object]:
    l_array, b_array, a_array = _population_array(low), _population_array(baseline), _population_array(method_a)
    common_reason = _unavailable_reason(low, baseline, method_a)
    hb = _metric(_hellinger(l_array, b_array)) if l_array is not None and b_array is not None else _metric(reason=common_reason)
    ha = _metric(_hellinger(l_array, a_array)) if l_array is not None and a_array is not None else _metric(reason=common_reason)
    delta = _metric(float(ha["value"]) - float(hb["value"])) if hb["available"] and ha["available"] else _metric(reason=common_reason)
    tv_lb = _metric(_total_variation(l_array, b_array)) if l_array is not None and b_array is not None else _metric(reason=common_reason)
    tv_la = _metric(_total_variation(l_array, a_array)) if l_array is not None and a_array is not None else _metric(reason=common_reason)
    refinement = _metric(_total_variation(a_array, b_array)) if a_array is not None and b_array is not None else _metric(reason=common_reason)
    return {"H_B": hb, "H_A": ha, "DeltaH": delta, "TV_LB": tv_lb, "TV_LA": tv_la, "R": refinement, **_alignment_metrics(l_array, b_array, a_array)}


def _joint_metrics(low: Mapping[str, object], baseline: Mapping[str, object], method_a: Mapping[str, object]) -> dict[str, object]:
    l_array, b_array, a_array = _population_matrix(low), _population_matrix(baseline), _population_matrix(method_a)
    common_reason = _unavailable_reason(low, baseline, method_a)
    tv_lb = _metric(_total_variation(l_array, b_array)) if l_array is not None and b_array is not None else _metric(reason=common_reason)
    tv_la = _metric(_total_variation(l_array, a_array)) if l_array is not None and a_array is not None else _metric(reason=common_reason)
    refinement = _metric(_total_variation(a_array, b_array)) if a_array is not None and b_array is not None else _metric(reason=common_reason)
    return {"TV_LB": tv_lb, "TV_LA": tv_la, "R": refinement, "kappa": _alignment_metrics(l_array, b_array, a_array)["kappa"]}


def _require_within_edges(values: np.ndarray, edges: Sequence[float], label: str) -> None:
    if values.size and (not np.all(np.isfinite(values)) or np.any(values < edges[0]) or np.any(values > edges[-1])):
        raise MethodAAcceptanceRefinementValidationError("{}_outside_frozen_edges".format(label))


def _one_dimensional_payload(low_values: Sequence[float], control_values: Sequence[float], baseline_weights: Sequence[float], method_a_weights: Sequence[float], edges: Sequence[float], label: str) -> dict[str, object]:
    low = np.asarray(low_values, dtype=float); control = np.asarray(control_values, dtype=float)
    baseline = np.asarray(baseline_weights, dtype=float); adjusted = np.asarray(method_a_weights, dtype=float)
    if control.size != baseline.size or control.size != adjusted.size:
        raise MethodAAcceptanceRefinementValidationError("{}_control_weight_alignment_invalid".format(label))
    _require_within_edges(low, edges, "{}_low".format(label)); _require_within_edges(control, edges, "{}_control".format(label))
    if baseline.size and (not np.all(np.isfinite(baseline)) or np.any(baseline <= 0.0)):
        raise MethodAAcceptanceRefinementValidationError("{}_baseline_weight_invalid".format(label))
    if adjusted.size and (not np.all(np.isfinite(adjusted)) or np.any(adjusted <= 0.0)):
        raise MethodAAcceptanceRefinementValidationError("{}_method_a_weight_invalid".format(label))
    l_contents = np.histogram(low, bins=np.asarray(edges, dtype=float))[0].astype(float)
    b_contents = np.histogram(control, bins=np.asarray(edges, dtype=float), weights=baseline)[0].astype(float)
    a_contents = np.histogram(control, bins=np.asarray(edges, dtype=float), weights=adjusted)[0].astype(float)
    l_pop = _population(l_contents, "low_population_empty" if low.size == 0 else None)
    b_pop = _population(b_contents, "control_population_empty" if control.size == 0 else None)
    a_pop = _population(a_contents, "control_population_empty" if control.size == 0 else None)
    return {"edges": [float(item) for item in edges], "L": l_pop, "B": b_pop, "A": a_pop, "metrics": _comparison_metrics(l_pop, b_pop, a_pop)}


def _two_dimensional_payload(low_x: Sequence[float], low_y: Sequence[float], control_x: Sequence[float], control_y: Sequence[float], baseline_weights: Sequence[float], method_a_weights: Sequence[float], x_edges: Sequence[float], y_edges: Sequence[float], label: str) -> dict[str, object]:
    lx, ly = np.asarray(low_x, dtype=float), np.asarray(low_y, dtype=float)
    cx, cy = np.asarray(control_x, dtype=float), np.asarray(control_y, dtype=float)
    baseline, adjusted = np.asarray(baseline_weights, dtype=float), np.asarray(method_a_weights, dtype=float)
    if lx.size != ly.size or cx.size != cy.size or cx.size != baseline.size or cx.size != adjusted.size:
        raise MethodAAcceptanceRefinementValidationError("{}_shape_alignment_invalid".format(label))
    _require_within_edges(lx, x_edges, "{}_low_x".format(label)); _require_within_edges(ly, y_edges, "{}_low_y".format(label))
    _require_within_edges(cx, x_edges, "{}_control_x".format(label)); _require_within_edges(cy, y_edges, "{}_control_y".format(label))
    if baseline.size and (not np.all(np.isfinite(baseline)) or np.any(baseline <= 0.0)):
        raise MethodAAcceptanceRefinementValidationError("{}_baseline_weight_invalid".format(label))
    if adjusted.size and (not np.all(np.isfinite(adjusted)) or np.any(adjusted <= 0.0)):
        raise MethodAAcceptanceRefinementValidationError("{}_method_a_weight_invalid".format(label))
    x_bins, y_bins = np.asarray(x_edges, dtype=float), np.asarray(y_edges, dtype=float)
    l_contents = np.histogram2d(lx, ly, bins=(x_bins, y_bins))[0].astype(float)
    b_contents = np.histogram2d(cx, cy, bins=(x_bins, y_bins), weights=baseline)[0].astype(float)
    a_contents = np.histogram2d(cx, cy, bins=(x_bins, y_bins), weights=adjusted)[0].astype(float)
    l_pop = _population(l_contents.ravel(), "low_population_empty" if lx.size == 0 else None)
    b_pop = _population(b_contents.ravel(), "control_population_empty" if cx.size == 0 else None)
    a_pop = _population(a_contents.ravel(), "control_population_empty" if cx.size == 0 else None)
    for population, contents in ((l_pop, l_contents), (b_pop, b_contents), (a_pop, a_contents)):
        if population["available"]:
            population["unit_area"] = [[float(item) for item in row] for row in (contents / float(np.sum(contents)))]
    return {"x_edges": [float(item) for item in x_edges], "y_edges": [float(item) for item in y_edges], "L": l_pop, "B": b_pop, "A": a_pop, "metrics": _joint_metrics(l_pop, b_pop, a_pop)}


def _coarsen_edges(edges: Sequence[float], label: str) -> list[float]:
    values = _strict_edges(edges, label, expected_length=41)
    selected = [values[index] for index in range(0, 41, 4)]
    if len(selected) != 11:
        raise MethodAAcceptanceRefinementValidationError("{}_coarsening_invalid".format(label))
    return selected


def _phi_index(value: float, edges: Sequence[float]) -> int | None:
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        if low <= value < high or (index == len(edges) - 2 and value == high):
            return index
    return None


def _validate_f6_1_artifact(artifact_value: object, f6_1_file_sha256: str, authority_by_kinematic: Mapping[str, object] | None) -> tuple[Mapping[str, object], Mapping[str, object], dict[str, object]]:
    artifact = _map(artifact_value, "f6_1_artifact")
    if artifact.get("schema_version") != _F6_1_ARTIFACT_SCHEMA:
        raise MethodAAcceptanceRefinementValidationError("f6_1_artifact_schema_invalid")
    validation = _map(artifact.get("validation"), "f6_1_validation")
    if validation.get("schema_version") != _F6_1_SCHEMA or validation.get("available") is not True:
        raise MethodAAcceptanceRefinementValidationError("f6_1_validation_invalid")
    if validation.get("fingerprint") != _f61._sha256(_map(validation.get("fingerprint_inputs"), "f6_1_fingerprint_inputs")):
        raise MethodAAcceptanceRefinementValidationError("f6_1_validation_fingerprint_invalid")
    provenance = _map(artifact.get("provenance"), "f6_1_provenance")
    expected_artifact = _f61._sha256({"schema_version": _F6_1_ARTIFACT_SCHEMA, "validation_fingerprint": validation.get("fingerprint"), "input_paths": provenance.get("input_paths")})
    if artifact.get("artifact_fingerprint") != expected_artifact:
        raise MethodAAcceptanceRefinementValidationError("f6_1_artifact_fingerprint_invalid")
    records = ACCEPTED_F6_1_ARTIFACT_AUTHORITY_BY_KINEMATIC if authority_by_kinematic is None else _map(authority_by_kinematic, "f6_2_f6_1_authority_records")
    raw = _map(records.get(_KINEMATIC), "f6_2_f6_1_authority")
    observed = {"source_file_sha256": _hash(f6_1_file_sha256, "f6_1_source_file_sha256"), "validation_fingerprint": _hash(validation.get("fingerprint"), "f6_1_validation_fingerprint"), "artifact_fingerprint": _hash(artifact.get("artifact_fingerprint"), "f6_1_artifact_fingerprint")}
    accepted: dict[str, str] = {}
    for name, value in observed.items():
        expected = _hash(raw.get(name), "f6_2_f6_1_authority_{}".format(name))
        if expected != value:
            raise MethodAAcceptanceRefinementValidationError("f6_2_f6_1_authority_{}_mismatch".format(name))
        accepted[name] = expected
    return artifact, validation, {"accepted": accepted, "observed": observed, "accepted_authority_match": True}


def _bootstrap_policy(config: Mapping[str, object] | None) -> dict[str, object]:
    policy = dict(_BOOTSTRAP_POLICY)
    policy["confidence_interval"] = dict(_BOOTSTRAP_POLICY["confidence_interval"])
    if config is None:
        return policy
    supplied = _map(config, "bootstrap_test_config")
    if set(supplied) != {"replicas"}:
        raise MethodAAcceptanceRefinementValidationError("bootstrap_test_config_invalid")
    replicas = _integer(supplied.get("replicas"), "bootstrap_test_replicas")
    if replicas <= 0 or replicas >= int(_BOOTSTRAP_POLICY["replicas"]):
        raise MethodAAcceptanceRefinementValidationError("bootstrap_test_replicas_invalid")
    policy["replicas"] = replicas
    policy["test_only_override"] = True
    return policy


def _child_seed(setting_id: str, t_index: int, phi_index: int, global_seed: int) -> int:
    material = "F6.2_BOOTSTRAP/v1|{}|{}|{}|{}".format(global_seed, setting_id, t_index, phi_index).encode("utf-8")
    return int.from_bytes(hashlib.sha256(material).digest()[:8], byteorder="big", signed=False)


def _linear_percentile(values: Sequence[float], percent: float) -> float:
    clean = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not clean:
        raise MethodAAcceptanceRefinementValidationError("linear_percentile_empty")
    if not (0.0 <= percent <= 100.0):
        raise MethodAAcceptanceRefinementValidationError("linear_percentile_percent_invalid")
    location = (len(clean) - 1) * (percent / 100.0)
    lower, upper = int(math.floor(location)), int(math.ceil(location))
    if lower == upper:
        return float(clean[lower])
    fraction = location - lower
    return float(clean[lower] + fraction * (clean[upper] - clean[lower]))


def _bootstrap_summary(values: Sequence[float], requested: int, name: str, policy: Mapping[str, object]) -> dict[str, object]:
    valid = [float(value) for value in values if math.isfinite(float(value))]
    invalid = requested - len(valid)
    needs_ninety_percent = name in {"kappa", "rho"}
    if not valid:
        return {"requested_replicas": requested, "valid_replicas": 0, "invalid_replicas": invalid, "available": False, "reason": "no_valid_bootstrap_replicas", "ci_low": None, "ci_high": None}
    if needs_ninety_percent and len(valid) < math.ceil(float(policy["minimum_valid_fraction_for_kappa_and_rho"]) * requested):
        return {"requested_replicas": requested, "valid_replicas": len(valid), "invalid_replicas": invalid, "available": False, "reason": "insufficient_valid_bootstrap_replicas", "ci_low": None, "ci_high": None}
    interval = _map(policy["confidence_interval"], "bootstrap_confidence_interval")
    return {"requested_replicas": requested, "valid_replicas": len(valid), "invalid_replicas": invalid, "available": True, "reason": None, "ci_low": _linear_percentile(valid, _finite(interval.get("low_percent"), "bootstrap_ci_low")), "ci_high": _linear_percentile(valid, _finite(interval.get("high_percent"), "bootstrap_ci_high"))}


def _metric_value(shape: Mapping[str, object], name: str) -> float | None:
    metrics = _map(shape.get("metrics"), "bootstrap_metrics")
    value = _map(metrics.get(name), "bootstrap_metric_{}".format(name))
    return float(value["value"]) if value.get("available") is True else None


def _effective_sample_size(weights: Sequence[float], label: str) -> dict[str, object]:
    values = np.asarray(weights, dtype=float)
    if values.size == 0:
        return _metric(reason="control_population_empty")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise MethodAAcceptanceRefinementValidationError("{}_invalid".format(label))
    denominator = float(np.sum(values * values)); numerator = float(np.sum(values)) ** 2
    return _metric(numerator / denominator) if denominator > 0.0 else _metric(reason="normalization_invalid")


def _runtime_kaon_window() -> dict[str, object]:
    try:
        config = resolve_analysis_runtime_config("4p4", "2p74")
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceRefinementValidationError("kaon_window_runtime_config_unavailable") from exc
    low, high = _finite(config.get("mm_min"), "kaon_window_mm_min"), _finite(config.get("mm_max"), "kaon_window_mm_max")
    if low != 1.10 or high != 1.16 or high <= low:
        raise MethodAAcceptanceRefinementValidationError("kaon_window_runtime_config_mismatch")
    return {"mm_min": low, "mm_max": high, "membership": "[low, high)", "runtime_config_source": str(config.get("source")), "runtime_config_hash": _hash(config.get("config_hash"), "kaon_window_config_hash")}


def _value(row: Mapping[str, object], name: str, label: str) -> float:
    return _finite(row.get(name), "{}_{}".format(label, name))


def _parent_display_edges(f6_1_validation: Mapping[str, object]) -> dict[tuple[str, int], dict[str, list[float]]]:
    result: dict[tuple[str, int], dict[str, list[float]]] = {}
    parents = _seq(f6_1_validation.get("parents"), "f6_1_parents")
    if len(parents) != 15:
        raise MethodAAcceptanceRefinementValidationError("f6_1_parent_inventory_invalid")
    for item in parents:
        parent = _map(item, "f6_1_parent")
        key = (str(parent.get("setting_id")), _integer(parent.get("canonical_t_index"), "f6_1_parent_t_index"))
        comparisons = _map(parent.get("prompt_shape_comparisons"), "f6_1_prompt_shape_comparisons")
        variables: dict[str, list[float]] = {}
        for variable in _ONE_DIMENSIONAL_VARIABLES:
            comparison = _map(comparisons.get(variable), "f6_1_comparison_{}".format(variable))
            variables[variable] = _strict_edges(comparison.get("edges"), "f6_1_edges_{}".format(variable), expected_length=41)
        signed = _map(parent.get("signed_background"), "f6_1_signed_background")
        signed_mm = _map(signed.get("analysis_MM"), "f6_1_signed_analysis_MM")
        variables["signed_analysis_MM"] = _strict_edges(signed_mm.get("edges"), "f6_1_signed_analysis_MM_edges")
        if key in result:
            raise MethodAAcceptanceRefinementValidationError("f6_1_parent_duplicate")
        result[key] = variables
    return result


def _raw_state(f1_artifacts: Sequence[Mapping[str, object]]) -> dict[str, dict[str, object]]:
    try:
        return _f61._raw_records(f1_artifacts)
    except (ValueError, _f61.MethodAReweightingValidationError) as exc:
        raise MethodAAcceptanceRefinementValidationError("f6_2_raw_f1_records_invalid:{}".format(exc)) from exc


def _reproduce_and_pair(
    f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f1_hashes: Mapping[str, str], f3_sha: str,
    persisted_f4: Mapping[str, object], accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None,
) -> tuple[dict[tuple[str, int], list[dict[str, object]]], dict[str, dict[str, object]], Mapping[str, object]]:
    try:
        pairs, _parents, recomputed_f4, raw = _f61._validate_and_pair_rows(f1_artifacts, f3_artifact, f1_hashes, f3_sha, persisted_f4, accepted_f3_runtime_authority_by_kinematic)
    except (ValueError, _f61.MethodAReweightingValidationError) as exc:
        raise MethodAAcceptanceRefinementValidationError("f6_2_f4_transient_review_reproduction_failed:{}".format(exc)) from exc
    return pairs, raw, recomputed_f4


def _geometry(contract: Mapping[str, object], setting_id: str) -> tuple[dict[int, tuple[float, float]], list[float], list[float]]:
    phi_edges = _strict_edges(contract.get("phi_edges"), "{}_phi_edges".format(setting_id), expected_length=10)
    delta_edges = _strict_edges(contract.get("delta_edges"), "{}_delta_edges".format(setting_id))
    t_geometry: dict[int, tuple[float, float]] = {}
    t_edges = _strict_edges(contract.get("t_edges"), "{}_t_edges".format(setting_id), expected_length=4)
    for index, (low, high) in enumerate(zip(t_edges, t_edges[1:])):
        t_geometry[index] = (low, high)
    if set(t_geometry) != {0, 1, 2}:
        raise MethodAAcceptanceRefinementValidationError("{}_canonical_t_inventory_invalid".format(setting_id))
    return t_geometry, phi_edges, delta_edges


def _child_rows(
    state: Mapping[str, object], pairs: Sequence[Mapping[str, object]], t_index: int, phi_index: int, phi_edges: Sequence[float], label: str,
) -> tuple[list[Mapping[str, object]], list[Mapping[str, object]], list[Mapping[str, object]]]:
    training = _map(state.get("training"), "{}_training".format(label))
    low: list[Mapping[str, object]] = []
    controls: list[Mapping[str, object]] = []
    full: list[Mapping[str, object]] = []
    for row in training.values():
        if _integer(row.get("t_index"), "{}_training_t_index".format(label)) != t_index or row.get("source_label") != "prompt":
            continue
        if row.get("response_class") != "low":
            continue
        npe = _value(row, "P_hgcer_npeSum", "{}_low".format(label))
        if row.get("nommcuts") is not True or not (0.0 < npe <= 2.0):
            raise MethodAAcceptanceRefinementValidationError("{}_low_population_invalid".format(label))
        degrees = math.degrees(_value(row, "phi", "{}_training".format(label)))
        assignment = _phi_index(degrees, phi_edges)
        if assignment is None:
            raise MethodAAcceptanceRefinementValidationError("{}_low_phi_outside_canonical_geometry".format(label))
        if assignment == phi_index:
            low.append(row)
    for item in pairs:
        pair = _map(item, "{}_pair".format(label)); row = _map(pair.get("row"), "{}_pair_row".format(label))
        source, entry = _identity(row, "{}_pair".format(label))
        row_t, row_phi = _integer(row.get("t_index"), "{}_application_t_index".format(label)), _integer(row.get("phi_index"), "{}_application_phi_index".format(label))
        degrees = _value(row, "phi_degrees", "{}_application".format(label))
        if row.get("phi_status") != "inside_phi" or _phi_index(degrees, phi_edges) != row_phi:
            raise MethodAAcceptanceRefinementValidationError("{}_application_phi_geometry_invalid".format(label))
        factor = _finite(pair.get("factor"), "{}_correction_factor".format(label))
        if factor <= 0.0 or not isinstance(pair.get("in_support"), bool):
            raise MethodAAcceptanceRefinementValidationError("{}_review_factor_or_support_invalid".format(label))
        if row_t != t_index or row_phi != phi_index:
            continue
        if source not in _SOURCE_LABELS:
            raise MethodAAcceptanceRefinementValidationError("{}_source_label_invalid".format(label))
        _value(row, "analysis_MM", "{}_full".format(label)); _value(row, "signed_baseline_event_contribution", "{}_full".format(label))
        full.append(pair)
        if source == "prompt":
            matched = training.get((source, entry))
            if matched is None or matched.get("response_class") != "control":
                raise MethodAAcceptanceRefinementValidationError("{}_prompt_identity_missing_training_control".format(label))
            training_degrees = math.degrees(_value(matched, "phi", "{}_training_control".format(label)))
            if not math.isclose(training_degrees, degrees, rel_tol=0.0, abs_tol=_TOLERANCE):
                raise MethodAAcceptanceRefinementValidationError("{}_phi_semantic_parity_mismatch".format(label))
            for name in _ONE_DIMENSIONAL_VARIABLES:
                if not math.isclose(_value(matched, name, "{}_training_control".format(label)), _value(row, name, "{}_application_control".format(label)), rel_tol=0.0, abs_tol=_TOLERANCE):
                    raise MethodAAcceptanceRefinementValidationError("{}_prompt_parity_mismatch:{}".format(label, name))
            w0 = _value(row, "baseline_pion_weight_w0", "{}_control".format(label))
            if w0 <= 0.0:
                raise MethodAAcceptanceRefinementValidationError("{}_baseline_weight_invalid".format(label))
            controls.append(pair)
    return low, controls, full


def _values(rows: Sequence[Mapping[str, object]], variable: str, *, training: bool, label: str) -> list[float]:
    values: list[float] = []
    for item in rows:
        row = item if training else _map(item.get("row"), "{}_row".format(label))
        values.append(_value(row, variable, label))
    return values


def _weights(control: Sequence[Mapping[str, object]], label: str) -> tuple[list[float], list[float]]:
    baseline, adjusted = [], []
    for item in control:
        row = _map(item.get("row"), "{}_row".format(label)); w0 = _value(row, "baseline_pion_weight_w0", label); factor = _finite(item.get("factor"), "{}_factor".format(label))
        if w0 <= 0.0 or factor <= 0.0:
            raise MethodAAcceptanceRefinementValidationError("{}_weight_invalid".format(label))
        baseline.append(w0); adjusted.append(w0 * factor)
    return baseline, adjusted


def _support_summary(control: Sequence[Mapping[str, object]], full: Sequence[Mapping[str, object]]) -> dict[str, object]:
    def summarize(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
        in_support = sum(bool(row.get("in_support")) for row in rows); total = len(rows)
        return {"in_support_count": in_support, "ood_count": total - in_support, "ood_fraction": float((total - in_support) / total) if total else None}
    return {"prompt_control": summarize(control), "full_physical_application": summarize(full)}


def _signed_payload(full: Sequence[Mapping[str, object]], mm_edges: Sequence[float], window: Mapping[str, object], label: str) -> tuple[dict[str, object], list[dict[str, float | str]]]:
    mm, baseline, adjusted = [], [], []
    by_source: dict[str, dict[str, float | int]] = {source: {"event_count": 0, "baseline_signed_sum": 0.0, "method_a_signed_sum": 0.0} for source in _SOURCE_LABELS}
    transient: list[dict[str, float | str]] = []
    for item in full:
        row = _map(item.get("row"), "{}_row".format(label)); source = str(row.get("source_label")); value = _value(row, "analysis_MM", label); b = _value(row, "signed_baseline_event_contribution", label); factor = _finite(item.get("factor"), "{}_factor".format(label))
        if source not in by_source or factor <= 0.0:
            raise MethodAAcceptanceRefinementValidationError("{}_signed_source_or_factor_invalid".format(label))
        a = b * factor
        if not math.isfinite(a):
            raise MethodAAcceptanceRefinementValidationError("{}_signed_method_a_nonfinite".format(label))
        mm.append(value); baseline.append(b); adjusted.append(a); transient.append({"source": source, "mm": value, "baseline": b, "adjusted": a})
        by_source[source]["event_count"] = int(by_source[source]["event_count"]) + 1
        by_source[source]["baseline_signed_sum"] = float(by_source[source]["baseline_signed_sum"]) + b
        by_source[source]["method_a_signed_sum"] = float(by_source[source]["method_a_signed_sum"]) + a
    mm_array, b_array, a_array = np.asarray(mm, dtype=float), np.asarray(baseline, dtype=float), np.asarray(adjusted, dtype=float)
    _require_within_edges(mm_array, mm_edges, "{}_signed_mm".format(label))
    b_contents = np.histogram(mm_array, bins=np.asarray(mm_edges, dtype=float), weights=b_array)[0].astype(float)
    a_contents = np.histogram(mm_array, bins=np.asarray(mm_edges, dtype=float), weights=a_array)[0].astype(float)
    low, high = _finite(window.get("mm_min"), "window_low"), _finite(window.get("mm_max"), "window_high")
    selected = (mm_array >= low) & (mm_array < high)
    pb, pa = float(np.sum(b_array[selected])), float(np.sum(a_array[selected]))
    fraction = _metric((pa - pb) / pb) if math.isfinite(pb) and pb != 0.0 else _metric(reason="baseline_window_sum_zero_or_nonfinite")
    vb, va = float(np.sum(b_array * b_array)), float(np.sum(a_array * a_array))
    variance = {"V_B": _metric(vb), "V_A": _metric(va), "R_V": _metric(va / vb) if vb != 0.0 else _metric(reason="baseline_variance_zero")}
    return ({"source_aggregates": [{"source_label": source, **by_source[source]} for source in _SOURCE_LABELS], "analysis_MM": {"edges": [float(item) for item in mm_edges], "baseline_signed_contents": [float(item) for item in b_contents], "method_a_signed_contents": [float(item) for item in a_contents]}, "kaon_window": {"P_B_K": _metric(pb), "P_A_K": _metric(pa), "DeltaP_K": _metric(pa - pb), "f_refine_K": fraction}, "variance_proxy": variance}, transient)


def _bootstrap_child(low: Sequence[Mapping[str, object]], control: Sequence[Mapping[str, object]], full: Sequence[Mapping[str, object]], variable_edges: Mapping[str, Sequence[float]], delta_edges: Sequence[float], policy: Mapping[str, object], setting_id: str, t_index: int, phi_index: int, window: Mapping[str, object]) -> dict[str, object]:
    requested = _integer(policy.get("replicas"), "bootstrap_replicas")
    seed = _child_seed(setting_id, t_index, phi_index, _integer(policy.get("global_seed"), "bootstrap_global_seed")); rng = np.random.default_rng(seed)
    low_values = {name: np.asarray(_values(low, name, training=True, label="bootstrap_low"), dtype=float) for name in _ONE_DIMENSIONAL_VARIABLES}
    control_values = {name: np.asarray(_values(control, name, training=False, label="bootstrap_control"), dtype=float) for name in _ONE_DIMENSIONAL_VARIABLES}
    baseline, adjusted = _weights(control, "bootstrap_control")
    b_array, a_array = np.asarray(baseline, dtype=float), np.asarray(adjusted, dtype=float)
    signed: dict[str, list[dict[str, float | str]]] = {source: [] for source in _SOURCE_LABELS}
    for item in full:
        row = _map(item.get("row"), "bootstrap_full_row"); source = str(row.get("source_label")); signed[source].append({"source": source, "mm": _value(row, "analysis_MM", "bootstrap_full"), "baseline": _value(row, "signed_baseline_event_contribution", "bootstrap_full"), "adjusted": _value(row, "signed_baseline_event_contribution", "bootstrap_full") * _finite(item.get("factor"), "bootstrap_full_factor")})
    values: dict[str, list[float]] = {"{}:{}".format(variable, metric): [] for variable in _ONE_DIMENSIONAL_VARIABLES for metric in ("DeltaH", "kappa", "rho")}
    joint_values: dict[str, list[float]] = {"{}__{}".format(x_name, y_name): [] for x_name, y_name in _JOINT_VARIABLES[:2]}
    window_values: list[float] = []
    for _ in range(requested):
        low_index = rng.integers(0, len(low), size=len(low)) if low else np.asarray([], dtype=int)
        control_index = rng.integers(0, len(control), size=len(control)) if control else np.asarray([], dtype=int)
        for variable in _ONE_DIMENSIONAL_VARIABLES:
            shape = _one_dimensional_payload(low_values[variable][low_index], control_values[variable][control_index], b_array[control_index], a_array[control_index], variable_edges[variable], "bootstrap_{}".format(variable))
            for metric in ("DeltaH", "kappa", "rho"):
                value = _metric_value(shape, metric)
                if value is not None:
                    values["{}:{}".format(variable, metric)].append(value)
        for x_name, y_name in _JOINT_VARIABLES[:2]:
            shape = _two_dimensional_payload(low_values[x_name][low_index], low_values[y_name][low_index], control_values[x_name][control_index], control_values[y_name][control_index], b_array[control_index], a_array[control_index], _coarsen_edges(variable_edges[x_name], "bootstrap_{}".format(x_name)), _coarsen_edges(variable_edges[y_name], "bootstrap_{}".format(y_name)), "bootstrap_{}__{}".format(x_name, y_name))
            value = _metric_value(shape, "kappa")
            if value is not None:
                joint_values["{}__{}".format(x_name, y_name)].append(value)
        selected_b, selected_a = [], []
        window_low, window_high = _finite(window.get("mm_min"), "bootstrap_window_low"), _finite(window.get("mm_max"), "bootstrap_window_high")
        for source in _SOURCE_LABELS:
            rows = signed[source]
            if rows:
                indices = rng.integers(0, len(rows), size=len(rows))
                for index in indices:
                    row = rows[int(index)]
                    if window_low <= float(row["mm"]) < window_high:
                        selected_b.append(float(row["baseline"])); selected_a.append(float(row["adjusted"]))
        if full:
            window_values.append(float(np.sum(selected_a) - np.sum(selected_b)))
    one_dimensional = {variable: {metric: _bootstrap_summary(values["{}:{}".format(variable, metric)], requested, metric, policy) for metric in ("DeltaH", "kappa", "rho")} for variable in _ONE_DIMENSIONAL_VARIABLES}
    joints = {name: {"kappa": _bootstrap_summary(raw, requested, "kappa", policy)} for name, raw in joint_values.items()}
    return {"child_seed": seed, "policy": {"replicas": requested, "global_seed": int(policy["global_seed"]), "confidence_interval": _copy(policy["confidence_interval"], "bootstrap_ci")}, "one_dimensional": one_dimensional, "joint_missing_mass_acceptance": joints, "kaon_window": {"DeltaP_K": _bootstrap_summary(window_values, requested, "DeltaP_K", policy)}}


def _child_payload(state: Mapping[str, object], pairs: Sequence[Mapping[str, object]], parent_edges: Mapping[str, Sequence[float]], t_index: int, phi_index: int, t_bounds: tuple[float, float], phi_edges: Sequence[float], delta_edges: Sequence[float], policy: Mapping[str, object], window: Mapping[str, object]) -> dict[str, object]:
    setting_id = str(state.get("setting_id")); label = "{}_t{}_phi{}".format(setting_id, t_index, phi_index)
    low, controls, full = _child_rows(state, pairs, t_index, phi_index, phi_edges, label)
    baseline, adjusted = _weights(controls, label)
    one_dimensional = {variable: _one_dimensional_payload(_values(low, variable, training=True, label="{}_low".format(label)), _values(controls, variable, training=False, label="{}_control".format(label)), baseline, adjusted, parent_edges[variable], "{}_{}".format(label, variable)) for variable in _ONE_DIMENSIONAL_VARIABLES}
    joint: dict[str, object] = {}
    for x_name, y_name in _JOINT_VARIABLES:
        x_edges = delta_edges if x_name == "SHMS_delta" else _coarsen_edges(parent_edges[x_name], "{}_{}".format(label, x_name))
        y_edges = _coarsen_edges(parent_edges[y_name], "{}_{}".format(label, y_name))
        payload = _two_dimensional_payload(_values(low, x_name, training=True, label="{}_low".format(label)), _values(low, y_name, training=True, label="{}_low".format(label)), _values(controls, x_name, training=False, label="{}_control".format(label)), _values(controls, y_name, training=False, label="{}_control".format(label)), baseline, adjusted, x_edges, y_edges, "{}_{}__{}".format(label, x_name, y_name))
        payload["x_variable"] = x_name; payload["y_variable"] = y_name
        joint["{}__{}".format(x_name, y_name)] = payload
    signed, _transient = _signed_payload(full, parent_edges["signed_analysis_MM"], window, label)
    return {"setting_id": setting_id, "canonical_t_index": t_index, "canonical_t_low": t_bounds[0], "canonical_t_high": t_bounds[1], "phi_index": phi_index, "phi_low": phi_edges[phi_index], "phi_high": phi_edges[phi_index + 1], "population_counts": {"N_low": len(low), "N_control": len(controls), "N_full": len(full), "full_by_source": {source: sum(1 for item in full if _map(item.get("row"), "full_row").get("source_label") == source) for source in _SOURCE_LABELS}}, "support": _support_summary(controls, full), "effective_sample_size": {"baseline_w0": _effective_sample_size(baseline, "{}_baseline_neff".format(label)), "method_a_w0_times_C": _effective_sample_size(adjusted, "{}_method_a_neff".format(label))}, "one_dimensional": one_dimensional, "joint_distributions": joint, "signed_background": signed, "bootstrap": _bootstrap_child(low, controls, full, parent_edges, delta_edges, policy, setting_id, t_index, phi_index, window), "automatic_case_assignment": False, "numerical_case_thresholds_applied": False, "manual_physics_review_required": True}


def _assert_aggregate_only_persistence(value: object) -> None:
    forbidden = {"entry_index", "event_id", "event_identity", "application_records", "method_a_training_records", "correction_factors", "raw_shape_factors", "in_support_mask", "event_corrections", "lookup_table", "raw_rows", "review_data"}
    def visit(item: object) -> None:
        if isinstance(item, Mapping):
            for key, child in item.items():
                if str(key) in forbidden:
                    raise MethodAAcceptanceRefinementValidationError("forbidden_event_persistence:{}".format(key))
                visit(child)
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray)):
            for child in item:
                visit(child)
    visit(value)


def build_pion_hgcer_method_a_acceptance_refinement_validation(
    f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], f5_artifact: Mapping[str, object], f6_1_artifact: Mapping[str, object], *,
    f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, f5_input_file_sha256: str, f6_1_input_file_sha256: str,
    accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f6_1_artifact_authority_by_kinematic: Mapping[str, object] | None = None,
    bootstrap_test_config: Mapping[str, object] | None = None,
) -> dict[str, object]:
    policy = _bootstrap_policy(bootstrap_test_config)
    f1_hashes = {str(key): _hash(value, "f1_input_file_sha256") for key, value in _map(f1_input_file_hashes, "f1_input_file_hashes").items()}
    f6_artifact, persisted_f6, f6_authority = _validate_f6_1_artifact(f6_1_artifact, f6_1_input_file_sha256, accepted_f6_1_artifact_authority_by_kinematic)
    try:
        recomputed_f6 = _f61.build_pion_hgcer_method_a_reweighting_validation(f1_artifacts, f3_artifact, f4_artifact, f5_artifact, f1_input_file_hashes=f1_hashes, f3_input_file_sha256=_hash(f3_input_file_sha256, "f3_input_file_sha256"), f4_input_file_sha256=_hash(f4_input_file_sha256, "f4_input_file_sha256"), f5_input_file_sha256=_hash(f5_input_file_sha256, "f5_input_file_sha256"), accepted_runtime_authority_by_kinematic=accepted_runtime_authority_by_kinematic, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
    except (ValueError, _f61.MethodAReweightingValidationError) as exc:
        raise MethodAAcceptanceRefinementValidationError("f6_2_f6_1_reproduction_failed:{}".format(exc)) from exc
    if recomputed_f6 != persisted_f6:
        raise MethodAAcceptanceRefinementValidationError("f6_2_f6_1_exact_reproduction_mismatch")
    if recomputed_f6.get("fingerprint") != f6_authority["accepted"]["validation_fingerprint"]:
        raise MethodAAcceptanceRefinementValidationError("f6_2_f6_1_validation_fingerprint_mismatch")
    persisted_f4 = _map(f4_artifact.get("correction"), "f4_correction")
    pairs_by_parent, raw, recomputed_f4 = _reproduce_and_pair(f1_artifacts, f3_artifact, f1_hashes, _hash(f3_input_file_sha256, "f3_input_file_sha256"), persisted_f4, accepted_f3_runtime_authority_by_kinematic)
    if recomputed_f4 != persisted_f4:
        raise MethodAAcceptanceRefinementValidationError("f6_2_f4_exact_reproduction_mismatch")
    display_edges = _parent_display_edges(recomputed_f6); window = _runtime_kaon_window()
    parents: list[dict[str, object]] = []
    for phi, epsilon in _CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); state = _map(raw.get(setting_id), "{}_raw_state".format(setting_id)); state = {**state, "setting_id": setting_id}
        t_geometry, phi_edges, delta_edges = _geometry(_map(state.get("contract"), "{}_contract".format(setting_id)), setting_id)
        for t_index in range(3):
            pairs = pairs_by_parent.get((setting_id, t_index))
            parent_edges = display_edges.get((setting_id, t_index))
            if pairs is None or parent_edges is None:
                raise MethodAAcceptanceRefinementValidationError("f6_2_parent_inventory_invalid")
            children = [_child_payload(state, pairs, parent_edges, t_index, phi_index, t_geometry[t_index], phi_edges, delta_edges, policy, window) for phi_index in range(9)]
            parents.append({"setting": _copy(state.get("setting"), "setting"), "setting_id": setting_id, "canonical_t_index": t_index, "canonical_t_low": t_geometry[t_index][0], "canonical_t_high": t_geometry[t_index][1], "phi_edges": phi_edges, "children": children})
    if len(parents) != 15 or sum(len(_seq(parent["children"], "children")) for parent in parents) != 135:
        raise MethodAAcceptanceRefinementValidationError("f6_2_child_inventory_invalid")
    core: dict[str, object] = {
        "schema_version": METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_SCHEMA_VERSION,
        "fingerprint_schema_version": METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_FINGERPRINT_SCHEMA_VERSION,
        "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete",
        "non_authoritative": True, "validation_only": True,
        "f6_1_consumed": True, "f6_1_modified": False, "f4_correction_consumed": True, "f4_correction_modified": False,
        "event_correction_evaluated_for_detached_validation": True, "event_correction_persisted": False,
        "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False,
        "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False,
        "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False,
        "method_b_numerical_dependency": False, "automatic_case_assignment": False,
        "numerical_case_thresholds_applied": False, "final_yield_uncertainty_claimed": False,
        "manual_review_required": True, "baseline_disagreement_is_not_automatic_failure": True,
        "kaon_region_leakage_is_permitted": True, "low_response_is_reference_shape_not_absolute_probability": True,
        "bootstrap_policy": policy, "kaon_window": window,
        "runtime_authority": _copy(recomputed_f6.get("runtime_authority"), "runtime_authority"),
        "f6_1_authority": f6_authority,
        "f6_1_reproduction": {"exact_payload_match": True, "validation_fingerprint": recomputed_f6.get("fingerprint"), "artifact_fingerprint": f6_artifact.get("artifact_fingerprint")},
        "f4_reproduction": {"exact_payload_match": True, "correction_fingerprint": recomputed_f4.get("fingerprint")},
        "f5_reproduction": _copy(recomputed_f6.get("f5_reproduction"), "f5_reproduction"),
        "f5_continuity": _copy(recomputed_f6.get("f5_continuity"), "f5_continuity"),
        "parents": parents,
    }
    fingerprint_inputs = {key: core[key] for key in ("schema_version", "fingerprint_schema_version", "bootstrap_policy", "kaon_window", "runtime_authority", "f6_1_authority", "f6_1_reproduction", "f4_reproduction", "f5_reproduction", "f5_continuity", "parents")}
    core["fingerprint_inputs"] = fingerprint_inputs; core["fingerprint"] = _sha256(fingerprint_inputs)
    _assert_aggregate_only_persistence(core)
    return _copy(core, "acceptance_refinement_validation")  # type: ignore[return-value]


def build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(
    f1_artifacts: Sequence[Mapping[str, object]], f3_artifact: Mapping[str, object], f4_artifact: Mapping[str, object], f5_artifact: Mapping[str, object], f6_1_artifact: Mapping[str, object], *,
    f1_input_file_hashes: Mapping[str, object], f3_input_file_sha256: str, f4_input_file_sha256: str, f5_input_file_sha256: str, f6_1_input_file_sha256: str,
    input_paths: Mapping[str, object] | None = None, generated_at_utc: str | None = None,
    accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None,
    accepted_f6_1_artifact_authority_by_kinematic: Mapping[str, object] | None = None,
    bootstrap_test_config: Mapping[str, object] | None = None,
) -> dict[str, object]:
    validation = build_pion_hgcer_method_a_acceptance_refinement_validation(f1_artifacts, f3_artifact, f4_artifact, f5_artifact, f6_1_artifact, f1_input_file_hashes=f1_input_file_hashes, f3_input_file_sha256=f3_input_file_sha256, f4_input_file_sha256=f4_input_file_sha256, f5_input_file_sha256=f5_input_file_sha256, f6_1_input_file_sha256=f6_1_input_file_sha256, accepted_runtime_authority_by_kinematic=accepted_runtime_authority_by_kinematic, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic, accepted_f6_1_artifact_authority_by_kinematic=accepted_f6_1_artifact_authority_by_kinematic, bootstrap_test_config=bootstrap_test_config)
    provenance = {"generated_at_utc": generated_at_utc, "input_paths": _copy({} if input_paths is None else input_paths, "input_paths")}
    artifact: dict[str, object] = {"schema_version": METHOD_A_ACCEPTANCE_REFINEMENT_VALIDATION_ARTIFACT_SCHEMA_VERSION, "validation": validation, "non_authoritative": True, "validation_only": True, "f6_1_consumed": True, "f6_1_modified": False, "event_correction_persisted": False, "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "method_b_numerical_dependency": False, "automatic_case_assignment": False, "numerical_case_thresholds_applied": False, "final_yield_uncertainty_claimed": False, "provenance": provenance}
    artifact["artifact_fingerprint"] = _sha256({"schema_version": artifact["schema_version"], "validation_fingerprint": validation["fingerprint"], "input_paths": provenance["input_paths"]})
    return _copy(artifact, "acceptance_refinement_validation_artifact")  # type: ignore[return-value]


def pion_hgcer_method_a_acceptance_refinement_validation_filename(kinematic: object) -> str:
    token = str(kinematic)
    if token != _KINEMATIC or token != token.strip() or any(char.isspace() for char in token) or any(char in token for char in "\\\\/:" ) or ".." in token:
        raise MethodAAcceptanceRefinementValidationError("acceptance_refinement_validation_filename_token_invalid")
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json".format(token)


def write_pion_hgcer_method_a_acceptance_refinement_validation_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    with open(os.fspath(path), "w", encoding="utf-8") as handle:
        json.dump(_copy(payload, "acceptance_refinement_validation_artifact"), handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return os.fspath(path)
