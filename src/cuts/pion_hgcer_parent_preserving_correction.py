"""Detached same-t parent-preserving map from frozen E.7 and Phase-A data.

The map is scientific-review evidence only.  It neither changes the existing
pion subtraction nor touches ROOT, events, histograms, yields, or weights.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os
from types import MappingProxyType


PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION = (
    "pion_hgcer_parent_preserving_correction/v1"
)
PION_HGCER_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION = (
    "pion_hgcer_parent_preserving_correction_artifact/v1"
)
PARENT_NORMALIZATION_EPSILON = 1.0e-12
_BASELINE_ALGEBRA_RELATIVE_TOLERANCE = 1.0e-12
_CLOSURE_TOLERANCE_DEFINITION = (
    "max(1e-12,1e-10*max(1,abs(canonical_baseline_before),"
    "abs(canonical_baseline_after)))"
)
_PROTOTYPE_SCHEMA = "pion_hgcer_ab_combination_prototype/v1"
_PHASE_A_SCHEMA = "pion_hgcer_event_contract/v1"
_PHASE_A_FINGERPRINT_SCHEMA = "pion_hgcer_event_contract_fingerprint/v2"
_PHASE_A_FINGERPRINT_EPHEMERAL_PROVENANCE_FIELDS = frozenset((
    "canonical_interval_pair_id",
))
_PHASE_D_SCHEMA = "pion_hgcer_phase_d_checkpoint/v1"
_SOURCE_TARGET_STATE = "post_proton_noRF"
_AVAILABILITY_STATES = {
    "both_comparable",
    "both_present_not_comparable",
    "a_only",
    "b_only",
    "neither_available",
}
_UNAVAILABLE_E7_REASONS = {
    "a_only": "method_b_unavailable",
    "b_only": "method_a_unavailable",
    "both_present_not_comparable": "frozen_ab_not_comparable",
    "neither_available": "both_methods_unavailable",
}


def _mapping(value):
    return value if isinstance(value, Mapping) else {}


def _sequence(value):
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _nonempty_string(value):
    return isinstance(value, str) and bool(value)


def _integer(value):
    if isinstance(value, bool):
        return None
    try:
        result = int(value)
    except (TypeError, ValueError):
        return None
    return result if value == result else None


def _finite(value):
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _strict_edges(value):
    if not _sequence(value) or len(value) < 2:
        return None
    edges = []
    for raw in value:
        scalar = _finite(raw)
        if scalar is None:
            return None
        edges.append(scalar)
    if any(right <= left for left, right in zip(edges, edges[1:])):
        return None
    return tuple(edges)


def _canonical_json(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _fingerprint(value):
    try:
        encoded = _canonical_json(value).encode("ascii")
    except (TypeError, ValueError, OverflowError):
        return None
    return hashlib.sha256(encoded).hexdigest()


def _serialized_equal(left, right):
    try:
        return _canonical_json(left) == _canonical_json(right)
    except (TypeError, ValueError, OverflowError):
        return False


def _json_copy(value, context="payload"):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("{}_contains_nonfinite_float".format(context))
        return float(value)
    if isinstance(value, MappingProxyType):
        value = dict(value)
    if isinstance(value, Mapping):
        return {
            str(key): _json_copy(child, "{}.{}".format(context, key))
            for key, child in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [
            _json_copy(child, "{}[{}]".format(context, index))
            for index, child in enumerate(value)
        ]
    try:
        if hasattr(value, "item"):
            return _json_copy(value.item(), context)
        if hasattr(value, "tolist"):
            return _json_copy(value.tolist(), context)
    except Exception as exc:
        raise ValueError("{}_is_not_json_safe".format(context)) from exc
    raise ValueError("{}_is_not_json_safe".format(context))


def _phase_a_json_ready(value):
    """Mirror the frozen Phase-A JSON projection for record fingerprints."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, dict):
        return {
            str(key): _phase_a_json_ready(child)
            for key, child in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_phase_a_json_ready(child) for child in value]
    try:
        if hasattr(value, "item"):
            return _phase_a_json_ready(value.item())
        if hasattr(value, "tolist"):
            return _phase_a_json_ready(value.tolist())
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return float(numeric) if math.isfinite(numeric) else None


def _phase_a_pion_event_population_fingerprint(records):
    """Reproduce Phase-A's complete ordered pion-record fingerprint exactly."""
    projected = []
    for record in records or ():
        fingerprint_record = _phase_a_json_ready(record)
        if isinstance(fingerprint_record, dict):
            fingerprint_record = {
                key: value
                for key, value in fingerprint_record.items()
                if key not in _PHASE_A_FINGERPRINT_EPHEMERAL_PROVENANCE_FIELDS
            }
        projected.append(fingerprint_record)
    return _fingerprint(projected)


def _unavailable(reason):
    """Return a stable, non-authoritative unavailable E.7.1 payload."""
    return {
        "schema_version": PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "method": "same_t_signed_baseline_parent_preserving_rescale",
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "event_application_performed": False,
        "source_ab_combination_prototype_fingerprint": None,
        "source_phase_d_checkpoint_schema": None,
        "source_checkpoint_payload_fingerprint": None,
        "phase_a_contract_fingerprint": None,
        "phase_a_event_population_fingerprint": None,
        "coordinate_fingerprint": None,
        "host_state": None,
        "source_target_state": None,
        "t_edges": [],
        "delta_edges": [],
        "parent_normalization_epsilon": PARENT_NORMALIZATION_EPSILON,
        "closure_tolerance_definition": _CLOSURE_TOLERANCE_DEFINITION,
        "parents": [],
        "cells": [],
        "uncertainty_model": "not_defined",
        "fingerprint_inputs": {},
        "fingerprint": None,
    }


def _prototype_contract(value):
    prototype = _mapping(value)
    required = (
        "schema_version", "status", "available", "reason", "method",
        "log_weight_method_a", "log_weight_method_b", "non_authoritative",
        "production_objects_mutated", "refinement_applied",
        "production_application_performed",
        "prototype_statistical_optimality_claimed", "prototype_uncertainty_model",
        "source_phase_d_checkpoint_schema", "source_checkpoint_payload_fingerprint",
        "phase_a_contract_fingerprint", "coordinate_fingerprint",
        "method_a_comparison_fingerprint", "method_b_comparison_fingerprint",
        "source_method_a_comparison_payload_fingerprint",
        "source_method_b_comparison_payload_fingerprint", "ab_comparison_fingerprint",
        "host_state", "source_target_state", "t_edges", "delta_edges", "cells",
        "fingerprint_inputs", "fingerprint",
    )
    if not prototype or any(key not in prototype for key in required):
        return None, None, None, "e7_prototype_contract_invalid"
    t_edges = _strict_edges(prototype["t_edges"])
    delta_edges = _strict_edges(prototype["delta_edges"])
    if (
        prototype["schema_version"] != _PROTOTYPE_SCHEMA
        or prototype["status"] != "available"
        or prototype["available"] is not True
        or prototype["reason"] is not None
        or prototype["method"] != "equal_weight_log_geometric_mean"
        or prototype["log_weight_method_a"] != 0.5
        or prototype["log_weight_method_b"] != 0.5
        or prototype["non_authoritative"] is not True
        or prototype["production_objects_mutated"] is not False
        or prototype["refinement_applied"] is not False
        or prototype["production_application_performed"] is not False
        or prototype["prototype_statistical_optimality_claimed"] is not False
        or prototype["prototype_uncertainty_model"] != "not_defined"
        or prototype["source_phase_d_checkpoint_schema"] != _PHASE_D_SCHEMA
        or not all(
            _nonempty_string(prototype[key])
            for key in (
                "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
                "coordinate_fingerprint", "method_a_comparison_fingerprint",
                "method_b_comparison_fingerprint",
                "source_method_a_comparison_payload_fingerprint",
                "source_method_b_comparison_payload_fingerprint",
                "ab_comparison_fingerprint", "fingerprint",
            )
        )
        or prototype["host_state"] not in {
            "proton_cleaned", "identity_no_proton_cleaning",
        }
        or prototype["source_target_state"] != _SOURCE_TARGET_STATE
        or t_edges is None
        or delta_edges is None
        or not _sequence(prototype["cells"])
        or not isinstance(prototype["fingerprint_inputs"], Mapping)
    ):
        return None, None, None, "e7_prototype_contract_invalid"
    expected_inputs = {
        "schema_version": _PROTOTYPE_SCHEMA,
        "method": "equal_weight_log_geometric_mean",
        "log_weight_method_a": 0.5,
        "log_weight_method_b": 0.5,
        "source_phase_d_checkpoint_schema": prototype["source_phase_d_checkpoint_schema"],
        "source_checkpoint_payload_fingerprint": prototype[
            "source_checkpoint_payload_fingerprint"
        ],
        "phase_a_contract_fingerprint": prototype["phase_a_contract_fingerprint"],
        "coordinate_fingerprint": prototype["coordinate_fingerprint"],
        "method_a_comparison_fingerprint": prototype[
            "method_a_comparison_fingerprint"
        ],
        "method_b_comparison_fingerprint": prototype[
            "method_b_comparison_fingerprint"
        ],
        "source_method_a_comparison_payload_fingerprint": prototype[
            "source_method_a_comparison_payload_fingerprint"
        ],
        "source_method_b_comparison_payload_fingerprint": prototype[
            "source_method_b_comparison_payload_fingerprint"
        ],
        "ab_comparison_fingerprint": prototype["ab_comparison_fingerprint"],
        "host_state": prototype["host_state"],
        "source_target_state": prototype["source_target_state"],
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "cells": _json_copy(prototype["cells"], "prototype.cells"),
    }
    if (
        not _serialized_equal(prototype["fingerprint_inputs"], expected_inputs)
        or prototype["fingerprint"] != _fingerprint(expected_inputs)
    ):
        return None, None, None, "e7_prototype_contract_invalid"
    return prototype, t_edges, delta_edges, None


def _phase_a_contract(value):
    phase_a = _mapping(value)
    required = (
        "schema_version", "fingerprint_schema_version", "status", "available",
        "immutable_record_contract", "production_objects_mutated",
        "refinement_applied", "rf_restoration_applied", "pion_closure",
        "host_closure", "contract_fingerprint", "fingerprint_inputs",
        "coordinate_fingerprint", "pion_event_population_fingerprint",
        "canonical_t_edges", "delta_edges", "pion_records", "host_state",
        "source_target_state",
    )
    if not phase_a or any(key not in phase_a for key in required):
        return None, None, None, "phase_a_contract_invalid"
    t_edges = _strict_edges(phase_a["canonical_t_edges"])
    delta_edges = _strict_edges(phase_a["delta_edges"])
    pion_closure = _mapping(phase_a["pion_closure"])
    host_closure = _mapping(phase_a["host_closure"])
    if (
        phase_a["schema_version"] != _PHASE_A_SCHEMA
        or phase_a["fingerprint_schema_version"] != _PHASE_A_FINGERPRINT_SCHEMA
        or phase_a["status"] != "available"
        or phase_a["available"] is not True
        or phase_a["immutable_record_contract"] is not True
        or phase_a["production_objects_mutated"] is not False
        or phase_a["refinement_applied"] is not False
        or phase_a["rf_restoration_applied"] is not False
        or pion_closure.get("passed") is not True
        or host_closure.get("passed") is not True
        or not all(
            _nonempty_string(phase_a[key])
            for key in (
                "contract_fingerprint", "coordinate_fingerprint",
                "pion_event_population_fingerprint",
            )
        )
        or not isinstance(phase_a["fingerprint_inputs"], Mapping)
        or phase_a["contract_fingerprint"] != _fingerprint(
            phase_a["fingerprint_inputs"]
        )
        or phase_a["host_state"] not in {
            "proton_cleaned", "identity_no_proton_cleaning",
        }
        or phase_a["source_target_state"] != _SOURCE_TARGET_STATE
        or t_edges is None
        or delta_edges is None
        or not _sequence(phase_a["pion_records"])
    ):
        return None, None, None, "phase_a_contract_invalid"
    return phase_a, t_edges, delta_edges, None


def _prototype_cell(value, t_edges, delta_edges):
    cell = _mapping(value)
    required = (
        "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
        "method_a", "method_b", "comparison", "prototype_status",
        "prototype_reason", "prototype_log_scale", "prototype_relative_scale",
    )
    if not cell or any(key not in cell for key in required):
        return None, "e7_prototype_cell_contract_invalid"
    t_index = _integer(cell["t_index"])
    delta_index = _integer(cell["delta_index"])
    if (
        t_index is None
        or delta_index is None
        or not 0 <= t_index < len(t_edges) - 1
        or not 0 <= delta_index < len(delta_edges) - 1
        or not _serialized_equal(cell["t_low"], t_edges[t_index])
        or not _serialized_equal(cell["t_high"], t_edges[t_index + 1])
        or not _serialized_equal(cell["delta_low"], delta_edges[delta_index])
        or not _serialized_equal(cell["delta_high"], delta_edges[delta_index + 1])
    ):
        return None, "e7_prototype_cell_geometry_invalid"
    method_a = _mapping(cell["method_a"])
    method_b = _mapping(cell["method_b"])
    comparison = _mapping(cell["comparison"])
    if (
        any(key not in method_a for key in ("present", "candidate", "low", "high", "status"))
        or any(key not in method_b for key in ("present", "candidate", "uncertainty", "status"))
        or any(
            key not in comparison
            for key in (
                "availability", "ratio_B_over_A", "log_ratio_B_over_A",
                "diagnostic_interval_relation",
            )
        )
        or not isinstance(method_a["present"], bool)
        or not isinstance(method_b["present"], bool)
        or comparison["availability"] not in _AVAILABILITY_STATES
    ):
        return None, "e7_prototype_cell_contract_invalid"
    availability = comparison["availability"]
    a_present = method_a["present"]
    b_present = method_b["present"]
    a_candidate = _finite(method_a["candidate"])
    a_low = _finite(method_a["low"])
    a_high = _finite(method_a["high"])
    b_candidate = _finite(method_b["candidate"])
    b_uncertainty = _finite(method_b["uncertainty"])
    if a_present:
        if (
            method_a["status"] not in {"available", "marginal"}
            or None in (a_candidate, a_low, a_high)
            or a_candidate < 0.0
            or a_low < 0.0
            or a_low > a_candidate
            or a_candidate > a_high
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif (
        method_a["status"] != "unavailable"
        or any(method_a[key] is not None for key in ("candidate", "low", "high"))
    ):
        return None, "e7_prototype_cell_contract_invalid"
    if b_present:
        if (
            method_b["status"] != "available_multi_region"
            or b_candidate is None
            or b_candidate <= 0.0
            or b_uncertainty is None
            or b_uncertainty <= 0.0
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif (
        method_b["status"] not in {
            "single_region_only", "unavailable", "region_marginal",
            "region_inconsistent", "shape_poor_veto",
        }
        or any(method_b[key] is not None for key in ("candidate", "uncertainty"))
    ):
        return None, "e7_prototype_cell_contract_invalid"
    relation = comparison["diagnostic_interval_relation"]
    ratio = comparison["ratio_B_over_A"]
    log_ratio = comparison["log_ratio_B_over_A"]
    if relation not in {"overlap", "disjoint", "not_evaluable"}:
        return None, "e7_prototype_cell_contract_invalid"
    if availability == "both_comparable":
        if (
            not a_present
            or not b_present
            or a_candidate <= 0.0
            or _finite(ratio) is None
            or _finite(log_ratio) is None
            or relation not in {"overlap", "disjoint"}
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif availability == "both_present_not_comparable":
        if (
            not a_present
            or not b_present
            or a_candidate != 0.0
            or ratio is not None
            or log_ratio is not None
            or relation != "not_evaluable"
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif availability == "a_only":
        if (
            not a_present
            or b_present
            or ratio is not None
            or log_ratio is not None
            or relation != "not_evaluable"
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif availability == "b_only":
        if (
            a_present
            or not b_present
            or ratio is not None
            or log_ratio is not None
            or relation != "not_evaluable"
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif (
        a_present
        or b_present
        or ratio is not None
        or log_ratio is not None
        or relation != "not_evaluable"
    ):
        return None, "e7_prototype_cell_contract_invalid"
    scale = _finite(cell["prototype_relative_scale"])
    log_scale = _finite(cell["prototype_log_scale"])
    if availability == "both_comparable":
        if (
            cell["prototype_status"] != "available"
            or cell["prototype_reason"] is not None
            or scale is None
            or scale <= 0.0
            or log_scale is None
            or not a_present
            or not b_present
        ):
            return None, "e7_prototype_cell_contract_invalid"
    elif (
        cell["prototype_status"] != "unavailable"
        or cell["prototype_reason"] != _UNAVAILABLE_E7_REASONS[availability]
        or cell["prototype_relative_scale"] is not None
        or cell["prototype_log_scale"] is not None
    ):
        return None, "e7_prototype_cell_contract_invalid"
    return {
        "t_index": int(t_index),
        "t_low": float(t_edges[t_index]),
        "t_high": float(t_edges[t_index + 1]),
        "delta_index": int(delta_index),
        "delta_low": float(delta_edges[delta_index]),
        "delta_high": float(delta_edges[delta_index + 1]),
        "e7_prototype_status": str(cell["prototype_status"]),
        "e7_prototype_reason": cell["prototype_reason"],
        "e7_raw_relative_scale": scale if scale is not None else None,
        "e7_comparison_availability": str(availability),
        "refinable": bool(
            cell["prototype_status"] == "available"
            and availability == "both_comparable"
            and scale is not None
            and scale > 0.0
        ),
    }, None


def _new_baseline_metrics():
    return {
        "baseline_record_count": 0,
        "baseline_signed_sum": 0.0,
        "baseline_absolute_support": 0.0,
        "baseline_sumw2": 0.0,
    }


def _finalize_baseline_metrics(metrics):
    result = dict(metrics)
    sumw2 = float(result["baseline_sumw2"])
    support = float(result["baseline_absolute_support"])
    result["baseline_neff"] = support * support / sumw2 if sumw2 > 0.0 else 0.0
    return result


def _record_baseline_metrics(metrics, contribution):
    metrics["baseline_record_count"] += 1
    metrics["baseline_signed_sum"] += contribution
    metrics["baseline_absolute_support"] += abs(contribution)
    metrics["baseline_sumw2"] += contribution * contribution


def _algebra_matches(contribution, coefficient, weight):
    reconstructed = coefficient * weight
    tolerance = max(
        1.0e-12,
        _BASELINE_ALGEBRA_RELATIVE_TOLERANCE
        * max(1.0, abs(contribution), abs(reconstructed)),
    )
    return abs(contribution - reconstructed) <= tolerance


def _collect_baseline_records(records, t_edges, delta_edges, cells):
    cell_metrics = {
        coordinate: _new_baseline_metrics() for coordinate in cells
    }
    outside_metrics = [
        _new_baseline_metrics() for _unused in range(len(t_edges) - 1)
    ]
    for record_number, raw_record in enumerate(records):
        record = _mapping(raw_record)
        if not record or "nommcuts" not in record or not isinstance(record["nommcuts"], bool):
            return None, None, "phase_a_pion_record_contract_invalid:{}".format(
                record_number
            )
        if record["nommcuts"] is not True:
            continue
        required = (
            "canonical_t_index", "canonical_t_lower_edge", "canonical_t_upper_edge",
            "delta_index", "delta_lower_edge", "delta_upper_edge",
            "signed_baseline_event_contribution", "signed_source_coefficient",
            "baseline_pion_weight_w0", "noRF_provenance",
        )
        if any(key not in record for key in required):
            return None, None, "phase_a_pion_record_contract_invalid:{}".format(
                record_number
            )
        t_index = _integer(record["canonical_t_index"])
        contribution = _finite(record["signed_baseline_event_contribution"])
        coefficient = _finite(record["signed_source_coefficient"])
        weight = _finite(record["baseline_pion_weight_w0"])
        if (
            t_index is None
            or not 0 <= t_index < len(t_edges) - 1
            or contribution is None
            or coefficient is None
            or weight is None
            or record["noRF_provenance"] != "noRF"
            or not _serialized_equal(record["canonical_t_lower_edge"], t_edges[t_index])
            or not _serialized_equal(record["canonical_t_upper_edge"], t_edges[t_index + 1])
            or not _algebra_matches(contribution, coefficient, weight)
        ):
            return None, None, "phase_a_pion_record_contract_invalid:{}".format(
                record_number
            )
        if record["delta_index"] is None:
            if record["delta_lower_edge"] is not None or record["delta_upper_edge"] is not None:
                return None, None, "phase_a_pion_record_contract_invalid:{}".format(
                    record_number
                )
            _record_baseline_metrics(outside_metrics[t_index], contribution)
            continue
        delta_index = _integer(record["delta_index"])
        if (
            delta_index is None
            or not 0 <= delta_index < len(delta_edges) - 1
            or not _serialized_equal(record["delta_lower_edge"], delta_edges[delta_index])
            or not _serialized_equal(record["delta_upper_edge"], delta_edges[delta_index + 1])
        ):
            return None, None, "phase_a_pion_record_contract_invalid:{}".format(
                record_number
            )
        _record_baseline_metrics(cell_metrics[(t_index, delta_index)], contribution)
    return cell_metrics, outside_metrics, None


def _closure_tolerance(before, after):
    return max(1.0e-12, 1.0e-10 * max(1.0, abs(before), abs(after)))


def _apply_identity(cells, refinable_indices, status, reason):
    refinable = set(refinable_indices)
    for cell in cells:
        if cell["delta_index"] in refinable:
            cell["C_final"] = 1.0
            cell["C_final_status"] = status
            cell["C_final_reason"] = reason
        else:
            cell["C_final"] = 1.0
            cell["C_final_status"] = "identity_unmodified"
            cell["C_final_reason"] = "e7_cell_not_refinable"
        cell["parent_normalization_divisor"] = None
        cell["uncertainty"] = None


def _parent_payload(t_index, t_edges, cells, outside_metrics):
    """Construct one parent while preserving only its signed canonical baseline."""
    refinable = [cell for cell in cells if cell["refinable"]]
    refinable_indices = [cell["delta_index"] for cell in refinable]
    canonical_before = sum(cell["baseline_signed_sum"] for cell in cells)
    outside = float(outside_metrics["baseline_signed_sum"])
    refinable_before = sum(cell["baseline_signed_sum"] for cell in refinable)
    refinable_scaled = sum(
        cell["baseline_signed_sum"] * cell["e7_raw_relative_scale"]
        for cell in refinable
    )
    attempted = {
        "refinable_baseline_signed_sum": refinable_before,
        "refinable_raw_scaled_signed_sum": refinable_scaled,
        "normalization_divisor": None,
        "canonical_baseline_after": None,
        "closure_difference": None,
        "closure_tolerance": None,
        "closure_passed": None,
    }
    parent_status = "available_parent_preserved"
    parent_reason = None
    if not refinable:
        parent_status = "identity_no_refinable_cells"
        parent_reason = "no_e7_refinable_cells"
        _apply_identity(cells, (), "identity_unmodified", "e7_cell_not_refinable")
    elif len(refinable) == 1:
        parent_status = "identity_single_refinable_cell"
        parent_reason = "parent_preservation_requires_multiple_refinable_cells"
        _apply_identity(
            cells,
            refinable_indices,
            "identity_single_refinable_cell",
            parent_reason,
        )
    else:
        normalization = None
        if refinable_before <= PARENT_NORMALIZATION_EPSILON:
            parent_status = "identity_parent_normalization_unavailable"
            parent_reason = "refinable_baseline_signed_sum_nonpositive"
        elif refinable_scaled <= PARENT_NORMALIZATION_EPSILON:
            parent_status = "identity_parent_normalization_unavailable"
            parent_reason = "refinable_raw_scaled_signed_sum_nonpositive"
        else:
            normalization = refinable_scaled / refinable_before
            attempted["normalization_divisor"] = normalization
            if not math.isfinite(normalization) or normalization <= 0.0:
                parent_status = "identity_parent_normalization_unavailable"
                parent_reason = "parent_normalization_divisor_invalid"
        if parent_reason is not None:
            _apply_identity(
                cells,
                refinable_indices,
                "identity_parent_normalization_unavailable",
                parent_reason,
            )
        else:
            for cell in cells:
                cell["uncertainty"] = None
                if cell["refinable"]:
                    cell["parent_normalization_divisor"] = normalization
                    cell["C_final"] = cell["e7_raw_relative_scale"] / normalization
                    cell["C_final_status"] = "available_parent_preserved"
                    cell["C_final_reason"] = None
                else:
                    cell["parent_normalization_divisor"] = normalization
                    cell["C_final"] = 1.0
                    cell["C_final_status"] = "identity_unmodified"
                    cell["C_final_reason"] = "e7_cell_not_refinable"
    canonical_after = sum(
        cell["baseline_signed_sum"] * cell["C_final"] for cell in cells
    )
    tolerance = _closure_tolerance(canonical_before, canonical_after)
    difference = canonical_after - canonical_before
    closure_passed = abs(difference) <= tolerance
    attempted.update({
        "canonical_baseline_after": canonical_after,
        "closure_difference": difference,
        "closure_tolerance": tolerance,
        "closure_passed": closure_passed,
    })
    if parent_status == "available_parent_preserved" and not closure_passed:
        parent_status = "identity_parent_normalization_unavailable"
        parent_reason = "parent_preservation_closure_failed"
        _apply_identity(
            cells,
            refinable_indices,
            "identity_parent_normalization_unavailable",
            parent_reason,
        )
        canonical_after = sum(
            cell["baseline_signed_sum"] * cell["C_final"] for cell in cells
        )
        tolerance = _closure_tolerance(canonical_before, canonical_after)
        difference = canonical_after - canonical_before
        closure_passed = abs(difference) <= tolerance
    final_normalization = (
        attempted["normalization_divisor"]
        if parent_status == "available_parent_preserved"
        else None
    )
    return {
        "t_index": int(t_index),
        "t_low": float(t_edges[t_index]),
        "t_high": float(t_edges[t_index + 1]),
        "refinable_delta_indices": list(refinable_indices),
        "refinable_cell_count": len(refinable_indices),
        "identity_delta_indices": [
            cell["delta_index"] for cell in cells if not cell["refinable"]
        ],
        "refinable_baseline_signed_sum": refinable_before,
        "refinable_raw_scaled_signed_sum": refinable_scaled,
        "normalization_divisor": final_normalization,
        "canonical_baseline_before": canonical_before,
        "canonical_baseline_after": canonical_after,
        "outside_delta_baseline": outside,
        "full_parent_baseline_before": canonical_before + outside,
        "full_parent_baseline_after": canonical_after + outside,
        "closure_difference": difference,
        "closure_tolerance": tolerance,
        "closure_passed": closure_passed,
        "negative_baseline_cell_count": sum(
            cell["baseline_signed_sum"] < 0.0 for cell in cells
        ),
        "zero_baseline_cell_count": sum(
            cell["baseline_signed_sum"] == 0.0 for cell in cells
        ),
        "parent_status": parent_status,
        "parent_reason": parent_reason,
        "attempted_normalization": attempted,
    }


def build_pion_hgcer_parent_preserving_correction(
    ab_combination_prototype,
    phase_a_event_contract,
):
    """Build the detached parent-preserving E.7.1 map from frozen parents."""
    try:
        prototype, t_edges, delta_edges, reason = _prototype_contract(
            ab_combination_prototype
        )
        if reason is not None:
            return _unavailable(reason)
        phase_a, phase_t_edges, phase_delta_edges, reason = _phase_a_contract(
            phase_a_event_contract
        )
        if reason is not None:
            return _unavailable(reason)
        if (
            prototype["phase_a_contract_fingerprint"] != phase_a["contract_fingerprint"]
            or prototype["coordinate_fingerprint"] != phase_a["coordinate_fingerprint"]
            or not _serialized_equal(t_edges, phase_t_edges)
            or not _serialized_equal(delta_edges, phase_delta_edges)
            or prototype["host_state"] != phase_a["host_state"]
            or prototype["source_target_state"] != phase_a["source_target_state"]
        ):
            return _unavailable("e7_phase_a_provenance_mismatch")
        expected_count = (len(t_edges) - 1) * (len(delta_edges) - 1)
        if len(prototype["cells"]) != expected_count:
            return _unavailable("e7_prototype_cell_grid_invalid")
        cells_by_coordinate = {}
        for source in prototype["cells"]:
            cell, reason = _prototype_cell(source, t_edges, delta_edges)
            if reason is not None:
                return _unavailable(reason)
            coordinate = (cell["t_index"], cell["delta_index"])
            if coordinate in cells_by_coordinate:
                return _unavailable("e7_prototype_cell_grid_invalid")
            cells_by_coordinate[coordinate] = cell
        if len(cells_by_coordinate) != expected_count:
            return _unavailable("e7_prototype_cell_grid_invalid")
        phase_a_population_fingerprint = phase_a[
            "pion_event_population_fingerprint"
        ]
        phase_a_input_population_fingerprint = phase_a["fingerprint_inputs"].get(
            "pion_event_population_fingerprint"
        )
        actual_population_fingerprint = _phase_a_pion_event_population_fingerprint(
            phase_a["pion_records"]
        )
        if (
            not all(
                _nonempty_string(value)
                for value in (
                    phase_a_population_fingerprint,
                    phase_a_input_population_fingerprint,
                    actual_population_fingerprint,
                )
            )
            or phase_a_population_fingerprint != phase_a_input_population_fingerprint
            or phase_a_population_fingerprint != actual_population_fingerprint
        ):
            return _unavailable(
                "phase_a_pion_event_population_fingerprint_mismatch"
            )
        baseline_metrics, outside_metrics, reason = _collect_baseline_records(
            phase_a["pion_records"], t_edges, delta_edges, cells_by_coordinate
        )
        if reason is not None:
            return _unavailable(reason)
        grouped = [[] for _unused in range(len(t_edges) - 1)]
        for coordinate in sorted(cells_by_coordinate):
            cell = dict(cells_by_coordinate[coordinate])
            cell.update(_finalize_baseline_metrics(baseline_metrics[coordinate]))
            cell["parent_normalization_divisor"] = None
            cell["C_final"] = None
            cell["C_final_status"] = None
            cell["C_final_reason"] = None
            cell["uncertainty"] = None
            grouped[cell["t_index"]].append(cell)
        parents = []
        cells = []
        for t_index, group in enumerate(grouped):
            group.sort(key=lambda row: row["delta_index"])
            if len(group) != len(delta_edges) - 1:
                return _unavailable("e7_prototype_cell_grid_invalid")
            parents.append(
                _parent_payload(
                    t_index,
                    t_edges,
                    group,
                    _finalize_baseline_metrics(outside_metrics[t_index]),
                )
            )
            cells.extend(group)
        fingerprint_inputs = {
            "schema_version": PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION,
            "method": "same_t_signed_baseline_parent_preserving_rescale",
            "source_ab_combination_prototype_fingerprint": prototype["fingerprint"],
            "source_checkpoint_payload_fingerprint": prototype[
                "source_checkpoint_payload_fingerprint"
            ],
            "phase_a_contract_fingerprint": phase_a["contract_fingerprint"],
            "phase_a_event_population_fingerprint": phase_a[
                "pion_event_population_fingerprint"
            ],
            "coordinate_fingerprint": phase_a["coordinate_fingerprint"],
            "host_state": phase_a["host_state"],
            "source_target_state": phase_a["source_target_state"],
            "t_edges": list(t_edges),
            "delta_edges": list(delta_edges),
            "parent_normalization_epsilon": PARENT_NORMALIZATION_EPSILON,
            "closure_tolerance_definition": _CLOSURE_TOLERANCE_DEFINITION,
            "parents": _json_copy(parents, "parents"),
            "cells": _json_copy(cells, "cells"),
        }
        fingerprint = _fingerprint(fingerprint_inputs)
        if fingerprint is None:
            return _unavailable("parent_preserving_fingerprint_inputs_invalid")
        return {
            "schema_version": PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION,
            "status": "available",
            "available": True,
            "reason": None,
            "method": "same_t_signed_baseline_parent_preserving_rescale",
            "non_authoritative": True,
            "production_objects_mutated": False,
            "refinement_applied": False,
            "production_application_performed": False,
            "event_application_performed": False,
            "source_ab_combination_prototype_fingerprint": prototype["fingerprint"],
            "source_phase_d_checkpoint_schema": prototype[
                "source_phase_d_checkpoint_schema"
            ],
            "source_checkpoint_payload_fingerprint": prototype[
                "source_checkpoint_payload_fingerprint"
            ],
            "phase_a_contract_fingerprint": phase_a["contract_fingerprint"],
            "phase_a_event_population_fingerprint": phase_a[
                "pion_event_population_fingerprint"
            ],
            "coordinate_fingerprint": phase_a["coordinate_fingerprint"],
            "host_state": phase_a["host_state"],
            "source_target_state": phase_a["source_target_state"],
            "t_edges": list(t_edges),
            "delta_edges": list(delta_edges),
            "parent_normalization_epsilon": PARENT_NORMALIZATION_EPSILON,
            "closure_tolerance_definition": _CLOSURE_TOLERANCE_DEFINITION,
            "parents": _json_copy(parents, "parents"),
            "cells": _json_copy(cells, "cells"),
            "uncertainty_model": "not_defined",
            "fingerprint_inputs": _json_copy(fingerprint_inputs, "fingerprint_inputs"),
            "fingerprint": fingerprint,
        }
    except Exception:
        return _unavailable("unexpected_parent_preserving_build_failure")


def _filename_token(value):
    if isinstance(value, (bool, list, tuple, dict, set, MappingProxyType)) or isinstance(
        value, Mapping
    ):
        raise ValueError("parent_preserving_correction_filename_token_invalid")
    if not isinstance(value, (str, int, float)):
        raise ValueError("parent_preserving_correction_filename_token_invalid")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("parent_preserving_correction_filename_token_invalid")
    raw = str(value)
    token = raw.strip()
    if (
        not token
        or token != raw
        or any(character.isspace() for character in token)
        or any(character in token for character in "\\\\/:")
        or ".." in token
    ):
        raise ValueError("parent_preserving_correction_filename_token_invalid")
    return token


def _setting(value):
    setting = _mapping(value)
    required = (
        "kinematic_token", "Q2", "W", "epsilon_setting",
        "epsilon_filename_token", "phi_setting", "particle_type",
    )
    if not setting or any(key not in setting for key in required):
        raise ValueError("parent_preserving_correction_setting_invalid")
    result = {
        "kinematic_token": _filename_token(setting["kinematic_token"]),
        "Q2": _json_copy(setting["Q2"], "setting.Q2"),
        "W": _json_copy(setting["W"], "setting.W"),
        "epsilon_setting": _filename_token(setting["epsilon_setting"]).lower(),
        "epsilon_filename_token": _filename_token(
            setting["epsilon_filename_token"]
        ).lower(),
        "phi_setting": _filename_token(setting["phi_setting"]),
        "particle_type": _filename_token(setting["particle_type"]).lower(),
    }
    if result["particle_type"] != "kaon":
        raise ValueError("parent_preserving_correction_particle_type_invalid")
    if result["epsilon_setting"] not in {"high", "low"}:
        raise ValueError("parent_preserving_correction_epsilon_setting_invalid")
    if result["epsilon_filename_token"] != "{}e".format(
        result["epsilon_setting"]
    ):
        raise ValueError("parent_preserving_correction_epsilon_filename_token_invalid")
    return result


def pion_hgcer_parent_preserving_correction_filename(
    phi_setting,
    particle_type,
    kinematic_token,
    epsilon_filename_token,
):
    """Return the deterministic E.7.1 downstream artifact basename."""
    phi = _filename_token(phi_setting)
    particle = _filename_token(particle_type).lower()
    kinematic = _filename_token(kinematic_token)
    epsilon = _filename_token(epsilon_filename_token).lower()
    if particle != "kaon":
        raise ValueError("parent_preserving_correction_particle_type_invalid")
    if epsilon not in {"highe", "lowe"}:
        raise ValueError("parent_preserving_correction_epsilon_filename_token_invalid")
    return "{}_{}_pion-background_hgcer_parent-preserving-correction_{}_{}.json".format(
        phi, particle, kinematic, epsilon
    )


def build_pion_hgcer_parent_preserving_correction_artifact(*, setting, correction):
    """Wrap one detached map with setting identity for JSON persistence only."""
    payload = _json_copy(correction, "correction")
    if not isinstance(payload, Mapping) or payload.get("schema_version") != (
        PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION
    ):
        raise ValueError("parent_preserving_correction_artifact_contract_invalid")
    if (
        payload.get("non_authoritative") is not True
        or payload.get("production_objects_mutated") is not False
        or payload.get("refinement_applied") is not False
        or payload.get("production_application_performed") is not False
        or payload.get("event_application_performed") is not False
    ):
        raise ValueError("parent_preserving_correction_artifact_authority_invalid")
    return {
        "schema_version": PION_HGCER_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION,
        "setting": _setting(setting),
        "correction": payload,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "event_application_performed": False,
    }


def write_pion_hgcer_parent_preserving_correction_json(path, payload):
    """Write one deterministic JSON-only E.7.1 artifact with a trailing newline."""
    target = os.fspath(path)
    serialized = _json_copy(payload, "artifact")
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(serialized, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return target


__all__ = (
    "PARENT_NORMALIZATION_EPSILON",
    "PION_HGCER_PARENT_PRESERVING_CORRECTION_ARTIFACT_SCHEMA_VERSION",
    "PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION",
    "build_pion_hgcer_parent_preserving_correction",
    "build_pion_hgcer_parent_preserving_correction_artifact",
    "pion_hgcer_parent_preserving_correction_filename",
    "write_pion_hgcer_parent_preserving_correction_json",
)
