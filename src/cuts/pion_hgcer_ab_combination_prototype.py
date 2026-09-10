"""Detached equal-log Method-A/legacy-Method-B scale-map prototype.

This module consumes only the frozen Phase-D checkpoint.  It is deliberately
non-authoritative: it creates no application weights and mutates no analysis
or ROOT object.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math


PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION = (
    "pion_hgcer_ab_combination_prototype/v1"
)

_PHASE_D_SCHEMA = "pion_hgcer_phase_d_checkpoint/v1"
_METHOD_A_COMPARISON_SCHEMA = "pion_hgcer_method_a_comparison/v1"
_METHOD_B_COMPARISON_SCHEMA = "pion_hgcer_method_b_comparison/v1"
_AB_COMPARISON_SCHEMA = "pion_hgcer_ab_comparison/v1"
_METHOD_A_COMPARISON_METHOD = "method_a_same_t_comparison_representation"
_METHOD_B_COMPARISON_METHOD = "method_b_comparison_representation"
_SOURCE_TARGET_STATE = "post_proton_noRF"
_AVAILABILITY_STATES = {
    "both_comparable",
    "both_present_not_comparable",
    "a_only",
    "b_only",
    "neither_available",
}
_UNAVAILABLE_REASONS = {
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


def _unavailable(reason):
    """Return a stable detached unavailable prototype payload."""
    return {
        "schema_version": PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "method": "equal_weight_log_geometric_mean",
        "log_weight_method_a": 0.5,
        "log_weight_method_b": 0.5,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "prototype_statistical_optimality_claimed": False,
        "prototype_uncertainty_model": "not_defined",
        "source_phase_d_checkpoint_schema": None,
        "source_checkpoint_payload_fingerprint": None,
        "phase_a_contract_fingerprint": None,
        "coordinate_fingerprint": None,
        "method_a_comparison_fingerprint": None,
        "method_b_comparison_fingerprint": None,
        "source_method_a_comparison_payload_fingerprint": None,
        "source_method_b_comparison_payload_fingerprint": None,
        "ab_comparison_fingerprint": None,
        "host_state": None,
        "source_target_state": None,
        "t_edges": [],
        "delta_edges": [],
        "cells": [],
        "fingerprint_inputs": {},
        "fingerprint": None,
    }


def _checkpoint_contract(value):
    checkpoint = _mapping(value)
    required = (
        "schema_version", "status", "available", "source_checkpoint_payload_fingerprint",
        "non_authoritative", "comparison_performed", "classification_performed",
        "classification_scope", "decision_performed", "statistical_compatibility_claimed",
        "production_objects_mutated", "refinement_applied", "method_a_comparison",
        "method_b_comparison", "ab_comparison",
    )
    if not checkpoint or any(key not in checkpoint for key in required):
        return None, "phase_d_checkpoint_contract_invalid"
    if checkpoint["schema_version"] != _PHASE_D_SCHEMA:
        return None, "phase_d_checkpoint_schema_invalid"
    if checkpoint["available"] is not True or checkpoint["status"] != "available":
        return None, "phase_d_checkpoint_unavailable"
    if (
        checkpoint["non_authoritative"] is not True
        or checkpoint["comparison_performed"] is not True
        or checkpoint["classification_performed"] is not True
        or checkpoint["classification_scope"] != "availability_only_non_prescriptive"
        or checkpoint["decision_performed"] is not False
        or checkpoint["statistical_compatibility_claimed"] is not False
        or checkpoint["production_objects_mutated"] is not False
        or checkpoint["refinement_applied"] is not False
    ):
        return None, "phase_d_checkpoint_authority_invalid"
    if not _nonempty_string(checkpoint["source_checkpoint_payload_fingerprint"]):
        return None, "phase_d_checkpoint_provenance_invalid"
    return checkpoint, None


def _representation_contract(value, schema, dependency_name, source_payload_name, reason):
    representation = _mapping(value)
    required = (
        "schema_version", "method", "status", "available", "source_checkpoint_payload_fingerprint",
        "phase_a_contract_fingerprint", "coordinate_fingerprint", source_payload_name,
        "canonical_t_edges", "delta_edges", "fingerprint_inputs", "fingerprint",
        "non_authoritative", dependency_name, "comparison_performed",
        "classification_performed", "production_objects_mutated", "refinement_applied",
    )
    if not representation or any(key not in representation for key in required):
        return None, reason
    expected_method = (
        _METHOD_A_COMPARISON_METHOD
        if schema == _METHOD_A_COMPARISON_SCHEMA else _METHOD_B_COMPARISON_METHOD
    )
    if (
        representation["schema_version"] != schema
        or representation["method"] != expected_method
        or representation["status"] != "available"
        or representation["available"] is not True
        or not all(
            _nonempty_string(representation[key])
            for key in (
                "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
                "coordinate_fingerprint", source_payload_name, "fingerprint",
            )
        )
        or _strict_edges(representation["canonical_t_edges"]) is None
        or _strict_edges(representation["delta_edges"]) is None
        or not isinstance(representation["fingerprint_inputs"], Mapping)
        or representation["fingerprint"] != _fingerprint(representation["fingerprint_inputs"])
        or representation["non_authoritative"] is not True
        or representation[dependency_name] is not False
        or representation["comparison_performed"] is not False
        or representation["classification_performed"] is not False
        or representation["production_objects_mutated"] is not False
        or representation["refinement_applied"] is not False
    ):
        return None, reason
    return representation, None


def _ab_contract(value):
    comparison = _mapping(value)
    required = (
        "schema_version", "method", "status", "available",
        "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
        "coordinate_fingerprint", "method_a_comparison_fingerprint",
        "method_b_comparison_fingerprint", "source_method_a_comparison_payload_fingerprint",
        "source_method_b_comparison_payload_fingerprint", "canonical_t_edges", "delta_edges",
        "host_state", "source_target_state", "cells", "fingerprint_inputs", "fingerprint",
        "non_authoritative", "comparison_performed", "classification_performed",
        "classification_scope", "decision_performed", "statistical_compatibility_claimed",
        "production_objects_mutated", "refinement_applied",
    )
    if not comparison or any(key not in comparison for key in required):
        return None, "ab_comparison_contract_invalid"
    if (
        comparison["schema_version"] != _AB_COMPARISON_SCHEMA
        or comparison["method"] != "non_authoritative_ab_comparison"
        or comparison["status"] != "available"
        or comparison["available"] is not True
        or not all(
            _nonempty_string(comparison[key])
            for key in (
                "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
                "coordinate_fingerprint", "method_a_comparison_fingerprint",
                "method_b_comparison_fingerprint",
                "source_method_a_comparison_payload_fingerprint",
                "source_method_b_comparison_payload_fingerprint", "fingerprint",
            )
        )
        or _strict_edges(comparison["canonical_t_edges"]) is None
        or _strict_edges(comparison["delta_edges"]) is None
        or comparison["host_state"] not in {"proton_cleaned", "identity_no_proton_cleaning"}
        or comparison["source_target_state"] != _SOURCE_TARGET_STATE
        or not _sequence(comparison["cells"])
        or not isinstance(comparison["fingerprint_inputs"], Mapping)
        or comparison["fingerprint"] != _fingerprint(comparison["fingerprint_inputs"])
        or comparison["non_authoritative"] is not True
        or comparison["comparison_performed"] is not True
        or comparison["classification_performed"] is not True
        or comparison["classification_scope"] != "availability_only_non_prescriptive"
        or comparison["decision_performed"] is not False
        or comparison["statistical_compatibility_claimed"] is not False
        or comparison["production_objects_mutated"] is not False
        or comparison["refinement_applied"] is not False
    ):
        return None, "ab_comparison_contract_invalid"
    return comparison, None


def _source_linkage(checkpoint, method_a, method_b, comparison):
    if (
        method_a["fingerprint"] != comparison["method_a_comparison_fingerprint"]
        or method_b["fingerprint"] != comparison["method_b_comparison_fingerprint"]
    ):
        return "ab_comparison_representation_fingerprint_mismatch"
    if (
        _fingerprint(method_a) != comparison["source_method_a_comparison_payload_fingerprint"]
        or _fingerprint(method_b) != comparison["source_method_b_comparison_payload_fingerprint"]
    ):
        return "ab_comparison_representation_payload_fingerprint_mismatch"
    source_fingerprint = checkpoint["source_checkpoint_payload_fingerprint"]
    if any(
        source.get("source_checkpoint_payload_fingerprint") != source_fingerprint
        for source in (method_a, method_b, comparison)
    ):
        return "ab_comparison_source_checkpoint_fingerprint_mismatch"
    if any(
        source["phase_a_contract_fingerprint"] != comparison["phase_a_contract_fingerprint"]
        or source["coordinate_fingerprint"] != comparison["coordinate_fingerprint"]
        for source in (method_a, method_b)
    ):
        return "ab_comparison_provenance_mismatch"
    if any(
        not _serialized_equal(source[key], comparison[key])
        for source in (method_a, method_b)
        for key in ("canonical_t_edges", "delta_edges")
    ):
        return "ab_comparison_geometry_mismatch"
    if (
        method_b.get("host_state") != comparison["host_state"]
        or method_b.get("source_target_state") != comparison["source_target_state"]
    ):
        return "ab_comparison_host_state_mismatch"
    return None


def _cell(source, t_edges, delta_edges):
    cell = _mapping(source)
    required = (
        "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
        "method_a", "method_b", "comparison",
    )
    if not cell or any(key not in cell for key in required):
        return None, "ab_comparison_cell_contract_invalid"
    t_index = _integer(cell["t_index"])
    delta_index = _integer(cell["delta_index"])
    if (
        t_index is None or delta_index is None
        or not 0 <= t_index < len(t_edges) - 1
        or not 0 <= delta_index < len(delta_edges) - 1
        or not _serialized_equal(cell["t_low"], t_edges[t_index])
        or not _serialized_equal(cell["t_high"], t_edges[t_index + 1])
        or not _serialized_equal(cell["delta_low"], delta_edges[delta_index])
        or not _serialized_equal(cell["delta_high"], delta_edges[delta_index + 1])
    ):
        return None, "ab_comparison_cell_geometry_invalid"
    method_a = _mapping(cell["method_a"])
    method_b = _mapping(cell["method_b"])
    relation = _mapping(cell["comparison"])
    required_a = (
        "present", "comparison_candidate", "comparison_candidate_low",
        "comparison_candidate_high", "comparison_candidate_status",
    )
    required_b = (
        "present", "comparison_candidate", "comparison_candidate_uncertainty",
        "comparison_candidate_status",
    )
    required_relation = (
        "availability", "ratio_B_over_A", "log_ratio_B_over_A",
        "diagnostic_interval_relation",
    )
    if (
        any(key not in method_a for key in required_a)
        or any(key not in method_b for key in required_b)
        or any(key not in relation for key in required_relation)
        or not isinstance(method_a["present"], bool)
        or not isinstance(method_b["present"], bool)
    ):
        return None, "ab_comparison_cell_contract_invalid"
    a_present = method_a["present"]
    b_present = method_b["present"]
    a_candidate = _finite(method_a["comparison_candidate"])
    a_low = _finite(method_a["comparison_candidate_low"])
    a_high = _finite(method_a["comparison_candidate_high"])
    if a_present:
        if (
            method_a["comparison_candidate_status"] not in {"available", "marginal"}
            or None in (a_candidate, a_low, a_high)
            or a_candidate < 0.0 or a_low < 0.0
            or a_low > a_candidate or a_candidate > a_high
        ):
            return None, "ab_comparison_method_a_cell_invalid"
    elif (
        method_a["comparison_candidate_status"] != "unavailable"
        or any(value is not None for value in (
            method_a["comparison_candidate"], method_a["comparison_candidate_low"],
            method_a["comparison_candidate_high"],
        ))
    ):
        return None, "ab_comparison_method_a_cell_invalid"
    b_candidate = _finite(method_b["comparison_candidate"])
    b_uncertainty = _finite(method_b["comparison_candidate_uncertainty"])
    if b_present:
        if (
            method_b["comparison_candidate_status"] != "available_multi_region"
            or b_candidate is None or b_candidate <= 0.0
            or b_uncertainty is None or b_uncertainty <= 0.0
        ):
            return None, "ab_comparison_method_b_cell_invalid"
    elif (
        method_b["comparison_candidate_status"] not in {
            "single_region_only", "unavailable", "region_marginal", "region_inconsistent",
            "shape_poor_veto",
        }
        or any(value is not None for value in (
            method_b["comparison_candidate"], method_b["comparison_candidate_uncertainty"],
        ))
    ):
        return None, "ab_comparison_method_b_cell_invalid"
    availability = relation["availability"]
    ratio = relation["ratio_B_over_A"]
    log_ratio = relation["log_ratio_B_over_A"]
    interval_relation = relation["diagnostic_interval_relation"]
    if availability not in _AVAILABILITY_STATES or interval_relation not in {
        "overlap", "disjoint", "not_evaluable",
    }:
        return None, "ab_comparison_availability_invalid"
    if availability == "both_comparable":
        if (
            not a_present or not b_present or a_candidate <= 0.0
            or _finite(ratio) is None or _finite(log_ratio) is None
            or interval_relation not in {"overlap", "disjoint"}
        ):
            return None, "ab_comparison_availability_invalid"
    elif availability == "both_present_not_comparable":
        if (
            not a_present or not b_present or a_candidate != 0.0
            or ratio is not None or log_ratio is not None or interval_relation != "not_evaluable"
        ):
            return None, "ab_comparison_availability_invalid"
    elif availability == "a_only":
        if (
            not a_present or b_present or ratio is not None or log_ratio is not None
            or interval_relation != "not_evaluable"
        ):
            return None, "ab_comparison_availability_invalid"
    elif availability == "b_only":
        if (
            a_present or not b_present or ratio is not None or log_ratio is not None
            or interval_relation != "not_evaluable"
        ):
            return None, "ab_comparison_availability_invalid"
    elif (
        a_present or b_present or ratio is not None or log_ratio is not None
        or interval_relation != "not_evaluable"
    ):
        return None, "ab_comparison_availability_invalid"
    return {
        "t_index": t_index,
        "t_low": float(t_edges[t_index]),
        "t_high": float(t_edges[t_index + 1]),
        "delta_index": delta_index,
        "delta_low": float(delta_edges[delta_index]),
        "delta_high": float(delta_edges[delta_index + 1]),
        "method_a": {
            "present": a_present,
            "candidate": a_candidate if a_present else None,
            "low": a_low if a_present else None,
            "high": a_high if a_present else None,
            "status": str(method_a["comparison_candidate_status"]),
        },
        "method_b": {
            "present": b_present,
            "candidate": b_candidate if b_present else None,
            "uncertainty": b_uncertainty if b_present else None,
            "status": str(method_b["comparison_candidate_status"]),
        },
        "comparison": {
            "availability": str(availability),
            "ratio_B_over_A": _finite(ratio) if ratio is not None else None,
            "log_ratio_B_over_A": _finite(log_ratio) if log_ratio is not None else None,
            "diagnostic_interval_relation": str(interval_relation),
        },
    }, None


def _prototype_cell(cell):
    """Create the fixed equal-log central value or an explicit unavailable cell."""
    availability = cell["comparison"]["availability"]
    prototype_status = "unavailable"
    prototype_reason = _UNAVAILABLE_REASONS.get(availability)
    log_scale = None
    relative_scale = None
    if availability == "both_comparable":
        method_a = cell["method_a"]["candidate"]
        method_b = cell["method_b"]["candidate"]
        log_scale = 0.5 * (math.log(method_a) + math.log(method_b))
        relative_scale = math.exp(log_scale)
        prototype_status = "available"
        prototype_reason = None
    return {
        **cell,
        "method_a": dict(cell["method_a"]),
        "method_b": dict(cell["method_b"]),
        "comparison": dict(cell["comparison"]),
        "prototype_status": prototype_status,
        "prototype_reason": prototype_reason,
        "prototype_log_scale": log_scale,
        "prototype_relative_scale": relative_scale,
    }


def build_pion_hgcer_ab_combination_prototype(phase_d_checkpoint):
    """Build a detached equal-log central-scale map from frozen Phase-D only."""
    try:
        checkpoint, reason = _checkpoint_contract(phase_d_checkpoint)
        if reason is not None:
            return _unavailable(reason)
        method_a, reason = _representation_contract(
            checkpoint["method_a_comparison"], _METHOD_A_COMPARISON_SCHEMA,
            "method_b_numerical_dependency", "source_method_a_payload_fingerprint",
            "method_a_comparison_contract_invalid",
        )
        if reason is not None:
            return _unavailable(reason)
        method_b, reason = _representation_contract(
            checkpoint["method_b_comparison"], _METHOD_B_COMPARISON_SCHEMA,
            "method_a_numerical_dependency", "source_method_b_payload_fingerprint",
            "method_b_comparison_contract_invalid",
        )
        if reason is not None:
            return _unavailable(reason)
        comparison, reason = _ab_contract(checkpoint["ab_comparison"])
        if reason is not None:
            return _unavailable(reason)
        reason = _source_linkage(checkpoint, method_a, method_b, comparison)
        if reason is not None:
            return _unavailable(reason)
        t_edges = _strict_edges(comparison["canonical_t_edges"])
        delta_edges = _strict_edges(comparison["delta_edges"])
        expected_count = (len(t_edges) - 1) * (len(delta_edges) - 1)
        if len(comparison["cells"]) != expected_count:
            return _unavailable("ab_comparison_cell_grid_invalid")
        cells = []
        seen = set()
        for source in comparison["cells"]:
            cell, reason = _cell(source, t_edges, delta_edges)
            if reason is not None:
                return _unavailable(reason)
            coordinate = (cell["t_index"], cell["delta_index"])
            if coordinate in seen:
                return _unavailable("ab_comparison_cell_grid_invalid")
            seen.add(coordinate)
            cells.append(_prototype_cell(cell))
        if len(seen) != expected_count:
            return _unavailable("ab_comparison_cell_grid_invalid")
        cells.sort(key=lambda row: (row["t_index"], row["delta_index"]))
        fingerprint_inputs = {
            "schema_version": PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION,
            "method": "equal_weight_log_geometric_mean",
            "log_weight_method_a": 0.5,
            "log_weight_method_b": 0.5,
            "source_phase_d_checkpoint_schema": checkpoint["schema_version"],
            "source_checkpoint_payload_fingerprint": checkpoint[
                "source_checkpoint_payload_fingerprint"
            ],
            "phase_a_contract_fingerprint": comparison["phase_a_contract_fingerprint"],
            "coordinate_fingerprint": comparison["coordinate_fingerprint"],
            "method_a_comparison_fingerprint": comparison[
                "method_a_comparison_fingerprint"
            ],
            "method_b_comparison_fingerprint": comparison[
                "method_b_comparison_fingerprint"
            ],
            "source_method_a_comparison_payload_fingerprint": comparison[
                "source_method_a_comparison_payload_fingerprint"
            ],
            "source_method_b_comparison_payload_fingerprint": comparison[
                "source_method_b_comparison_payload_fingerprint"
            ],
            "ab_comparison_fingerprint": comparison["fingerprint"],
            "host_state": comparison["host_state"],
            "source_target_state": comparison["source_target_state"],
            "t_edges": list(t_edges),
            "delta_edges": list(delta_edges),
            "cells": json.loads(_canonical_json(cells)),
        }
        fingerprint = _fingerprint(fingerprint_inputs)
        if fingerprint is None:
            return _unavailable("prototype_fingerprint_inputs_invalid")
        return {
            "schema_version": PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION,
            "status": "available",
            "available": True,
            "reason": None,
            "method": "equal_weight_log_geometric_mean",
            "log_weight_method_a": 0.5,
            "log_weight_method_b": 0.5,
            "non_authoritative": True,
            "production_objects_mutated": False,
            "refinement_applied": False,
            "production_application_performed": False,
            "prototype_statistical_optimality_claimed": False,
            "prototype_uncertainty_model": "not_defined",
            "source_phase_d_checkpoint_schema": checkpoint["schema_version"],
            "source_checkpoint_payload_fingerprint": checkpoint[
                "source_checkpoint_payload_fingerprint"
            ],
            "phase_a_contract_fingerprint": comparison["phase_a_contract_fingerprint"],
            "coordinate_fingerprint": comparison["coordinate_fingerprint"],
            "method_a_comparison_fingerprint": comparison[
                "method_a_comparison_fingerprint"
            ],
            "method_b_comparison_fingerprint": comparison[
                "method_b_comparison_fingerprint"
            ],
            "source_method_a_comparison_payload_fingerprint": comparison[
                "source_method_a_comparison_payload_fingerprint"
            ],
            "source_method_b_comparison_payload_fingerprint": comparison[
                "source_method_b_comparison_payload_fingerprint"
            ],
            "ab_comparison_fingerprint": comparison["fingerprint"],
            "host_state": comparison["host_state"],
            "source_target_state": comparison["source_target_state"],
            "t_edges": list(t_edges),
            "delta_edges": list(delta_edges),
            "cells": json.loads(_canonical_json(cells)),
            "fingerprint_inputs": json.loads(_canonical_json(fingerprint_inputs)),
            "fingerprint": fingerprint,
        }
    except Exception:
        return _unavailable("unexpected_prototype_build_failure")


__all__ = (
    "PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION",
    "build_pion_hgcer_ab_combination_prototype",
)
