"""Read-only Phase-D comparison-input contract for frozen HGCer diagnostics.

This module deliberately validates only the detached Phase-C checkpoint.  It
does not choose a Method-A observable, transform Method-B, compare numerical
results, or own any analysis object.
"""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import json
import math


COMPARISON_INPUT_SCHEMA_VERSION = "pion_hgcer_ab_comparison_input/v1"

_CHECKPOINT_SCHEMA_VERSION = "pion_hgcer_refinement_checkpoint/v1"
_PHASE_A_SCHEMA_VERSION = "pion_hgcer_event_contract/v1"
_PHASE_A_FINGERPRINT_SCHEMA_VERSION = "pion_hgcer_event_contract_fingerprint/v2"
_METHOD_A_SCHEMA_VERSION = "pion_hgcer_method_a/v1"
_METHOD_B_SCHEMA_VERSION = "pion_hgcer_method_b/v1"
_HOST_STATES = {"proton_cleaned", "identity_no_proton_cleaning"}
_SOURCE_TARGET_STATE = "post_proton_noRF"


class _ComparisonInputUnavailable(Exception):
    """Expected validation failure carrying the public reason and stage."""

    def __init__(self, reason, stage):
        super().__init__(reason)
        self.reason = str(reason)
        self.stage = str(stage)


def _canonical_json(value):
    """Return the one JSON representation used for snapshots and identity."""
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _unavailable(reason, stage, exception=None):
    """Return the detached, non-authoritative failure shape."""
    result = {
        "schema_version": COMPARISON_INPUT_SCHEMA_VERSION,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "non_authoritative": True,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    if exception is not None:
        result["exception_type"] = type(exception).__name__
        result["exception_message"] = str(exception)
    return result


def _fail(reason, stage):
    raise _ComparisonInputUnavailable(reason, stage)


def _is_finite_number(value):
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
    )


def _is_integer(value):
    return not isinstance(value, bool) and isinstance(value, int)


def _nonempty_string(value):
    return isinstance(value, str) and bool(value)


def _serialized_equal(left, right):
    """Compare values exactly as serialized checkpoint data, not numerically."""
    return _canonical_json(left) == _canonical_json(right)


def _mapping(value, reason, stage):
    if not isinstance(value, Mapping):
        _fail(reason, stage)
    return value


def _require_keys(mapping, keys, reason, stage):
    if any(key not in mapping for key in keys):
        _fail(reason, stage)


def _validate_checkpoint(snapshot):
    stage = "checkpoint_validation"
    _mapping(snapshot, "checkpoint_contract_invalid", stage)
    _require_keys(
        snapshot,
        (
            "schema_version",
            "setting",
            "phase_a",
            "method_a",
            "method_b",
            "canonical_t_edges",
            "delta_edges",
            "host_state_summary",
            "non_authoritative",
            "production_objects_mutated",
            "refinement_applied",
        ),
        "checkpoint_contract_invalid",
        stage,
    )
    if snapshot["schema_version"] != _CHECKPOINT_SCHEMA_VERSION:
        _fail("checkpoint_schema_mismatch", stage)
    if (
        snapshot["non_authoritative"] is not True
        or snapshot["production_objects_mutated"] is not False
        or snapshot["refinement_applied"] is not False
    ):
        _fail("checkpoint_authority_contract_invalid", stage)


def _validate_setting(setting):
    stage = "setting_validation"
    setting = _mapping(setting, "setting_contract_invalid", stage)
    _require_keys(
        setting,
        (
            "kinematic_token",
            "Q2",
            "W",
            "epsilon_setting",
            "epsilon_filename_token",
            "phi_setting",
            "particle_type",
        ),
        "setting_contract_invalid",
        stage,
    )
    if not _nonempty_string(setting["kinematic_token"]):
        _fail("setting_contract_invalid", stage)
    if not _is_finite_number(setting["Q2"]) or not _is_finite_number(setting["W"]):
        _fail("setting_contract_invalid", stage)
    if not _nonempty_string(setting["phi_setting"]):
        _fail("setting_contract_invalid", stage)
    if setting["particle_type"] != "kaon":
        _fail("setting_contract_invalid", stage)
    epsilon = setting["epsilon_setting"]
    filename_token = setting["epsilon_filename_token"]
    if not isinstance(epsilon, str) or epsilon not in {"low", "high"}:
        _fail("setting_contract_invalid", stage)
    if not isinstance(filename_token, str) or filename_token != "{}e".format(epsilon):
        _fail("setting_contract_invalid", stage)


def _strict_edges(value, invalid_reason, stage):
    if not isinstance(value, list) or len(value) < 2:
        _fail(invalid_reason, stage)
    previous = None
    for edge in value:
        if not _is_finite_number(edge):
            _fail(invalid_reason, stage)
        if previous is not None and edge <= previous:
            _fail(invalid_reason, stage)
        previous = edge
    return value


def _reported_summary_edges(summary, key):
    """Return geometry only when an unavailable method actually reports it."""
    if key not in summary or summary[key] is None or summary[key] == []:
        return None
    return summary[key]


def _validate_phase_a(phase):
    stage = "phase_a_provenance"
    phase = _mapping(phase, "phase_a_unavailable", stage)
    summary = _mapping(phase.get("summary"), "phase_a_unavailable", stage)
    if summary.get("status") != "available" or summary.get("available") is not True:
        _fail("phase_a_unavailable", stage)
    contract_fingerprint = phase.get("contract_fingerprint")
    if not _nonempty_string(contract_fingerprint):
        _fail("phase_a_contract_fingerprint_missing", stage)
    if summary.get("contract_fingerprint") != contract_fingerprint:
        _fail("phase_a_contract_fingerprint_mismatch", stage)
    coordinate_fingerprint = phase.get("coordinate_fingerprint")
    if not _nonempty_string(coordinate_fingerprint):
        _fail("coordinate_fingerprint_missing", stage)
    if summary.get("coordinate_fingerprint") != coordinate_fingerprint:
        _fail("coordinate_fingerprint_mismatch", stage)
    if summary.get("schema_version") != _PHASE_A_SCHEMA_VERSION:
        _fail("phase_a_schema_mismatch", stage)
    if (
        summary.get("fingerprint_schema_version")
        != _PHASE_A_FINGERPRINT_SCHEMA_VERSION
    ):
        _fail("phase_a_fingerprint_schema_mismatch", stage)
    if (
        summary.get("pion_closure_passed") is not True
        or summary.get("host_closure_passed") is not True
    ):
        _fail("phase_a_closure_unavailable", stage)
    if (
        summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_restoration_applied") is not False
    ):
        _fail("phase_a_authority_contract_invalid", stage)
    return phase, summary, contract_fingerprint, coordinate_fingerprint


def _validate_geometry(snapshot, phase_summary, method_a, method_b):
    stage = "canonical_geometry"
    t_edges = _strict_edges(
        snapshot["canonical_t_edges"], "canonical_t_edges_invalid", stage
    )
    delta_edges = _strict_edges(snapshot["delta_edges"], "delta_edges_invalid", stage)
    if not _serialized_equal(phase_summary.get("canonical_t_edges"), t_edges):
        _fail("canonical_t_edges_mismatch", stage)
    if not _serialized_equal(phase_summary.get("delta_edges"), delta_edges):
        _fail("delta_edges_mismatch", stage)
    for method in (method_a, method_b):
        summary = method.get("summary") if isinstance(method, Mapping) else None
        if not isinstance(summary, Mapping):
            continue
        reported_t_edges = _reported_summary_edges(summary, "t_edges")
        if reported_t_edges is not None and not _serialized_equal(reported_t_edges, t_edges):
            _fail("canonical_t_edges_mismatch", stage)
        reported_delta_edges = _reported_summary_edges(summary, "delta_edges")
        if (
            reported_delta_edges is not None
            and not _serialized_equal(reported_delta_edges, delta_edges)
        ):
            _fail("delta_edges_mismatch", stage)
    return t_edges, delta_edges


def _validate_host_state(snapshot, phase, phase_summary):
    stage = "host_provenance"
    host_summary = _mapping(
        snapshot.get("host_state_summary"), "host_state_invalid", stage
    )
    phase_host_state = phase.get("host_state")
    serialized_states = (
        phase_host_state,
        phase_summary.get("host_state"),
        host_summary.get("phase_a_host_state"),
    )
    if any(
        not isinstance(state, str) or state not in _HOST_STATES
        for state in serialized_states
    ):
        _fail("host_state_invalid", stage)
    if any(state != phase_host_state for state in serialized_states):
        _fail("host_state_mismatch", stage)
    source_states = (
        phase.get("source_target_state"),
        phase_summary.get("source_target_state"),
        host_summary.get("source_target_state"),
    )
    if any(state != _SOURCE_TARGET_STATE for state in source_states):
        _fail("source_target_state_mismatch", stage)
    method_b_host_state = host_summary.get("method_b_host_state")
    if method_b_host_state is not None and method_b_host_state != phase_host_state:
        _fail("host_state_mismatch", stage)
    return phase_host_state


def _validate_method_structure(method, name):
    stage = "{}_provenance".format(name)
    reason = "{}_status_contract_invalid".format(name)
    method = _mapping(method, reason, stage)
    required = ("status", "available", "reason", "fingerprint", "summary", "cells")
    _require_keys(method, required, reason, stage)
    if name == "method_b":
        _require_keys(method, ("parent_region_references",), reason, stage)
    status = method["status"]
    available = method["available"]
    if (
        not isinstance(status, str)
        or status not in {"available", "unavailable"}
        or not isinstance(available, bool)
        or (status == "available") != available
    ):
        _fail(reason, stage)
    summary = _mapping(method["summary"], reason, stage)
    if summary.get("status") != status or summary.get("available") is not available:
        _fail(reason, stage)
    if not isinstance(method["cells"], list):
        _fail("{}_cell_grid_mismatch".format(name), "{}_geometry".format(name))
    if name == "method_b" and not isinstance(method["parent_region_references"], list):
        _fail("method_b_parent_reference_geometry_mismatch", "method_b_geometry")
    return method, summary, available


def _validate_method_a_provenance(method, summary, available, phase_fingerprint, coordinate):
    if not available:
        return
    stage = "method_a_provenance"
    fingerprint = method.get("fingerprint")
    if not _nonempty_string(fingerprint):
        _fail("method_a_fingerprint_missing", stage)
    if summary.get("fingerprint") != fingerprint:
        _fail("method_a_fingerprint_mismatch", stage)
    if summary.get("schema_version") != _METHOD_A_SCHEMA_VERSION:
        _fail("method_a_schema_mismatch", stage)
    if summary.get("phase_a_contract_fingerprint") != phase_fingerprint:
        _fail("method_a_phase_a_fingerprint_mismatch", stage)
    if summary.get("coordinate_fingerprint") != coordinate:
        _fail("method_a_coordinate_fingerprint_mismatch", stage)
    if (
        summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_ct_required") is not False
    ):
        _fail("method_a_authority_contract_invalid", stage)


def _validate_method_b_provenance(
    method, summary, available, phase_fingerprint, coordinate, host_state
):
    if not available:
        return
    stage = "method_b_provenance"
    fingerprint = method.get("fingerprint")
    if not _nonempty_string(fingerprint):
        _fail("method_b_fingerprint_missing", stage)
    if summary.get("fingerprint") != fingerprint:
        _fail("method_b_fingerprint_mismatch", stage)
    if summary.get("schema_version") != _METHOD_B_SCHEMA_VERSION:
        _fail("method_b_schema_mismatch", stage)
    if summary.get("phase_a_contract_fingerprint") != phase_fingerprint:
        _fail("method_b_phase_a_fingerprint_mismatch", stage)
    if summary.get("coordinate_fingerprint") != coordinate:
        _fail("method_b_coordinate_fingerprint_mismatch", stage)
    if summary.get("host_state") != host_state:
        _fail("method_b_host_state_mismatch", stage)
    if (
        summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_ct_required") is not False
        or summary.get("interpolation_used") is not False
    ):
        _fail("method_b_authority_contract_invalid", stage)


def _cell_indexes(cell, method_name, t_count, delta_count, required):
    grid_reason = "{}_cell_grid_mismatch".format(method_name)
    stage = "{}_geometry".format(method_name)
    if not isinstance(cell, Mapping):
        _fail(grid_reason, stage)
    if required and ("t_index" not in cell or "delta_index" not in cell):
        _fail(grid_reason, stage)
    t_index = cell.get("t_index")
    delta_index = cell.get("delta_index")
    if t_index is not None:
        if not _is_integer(t_index) or not 0 <= t_index < t_count:
            _fail(grid_reason, stage)
    if delta_index is not None:
        if not _is_integer(delta_index) or not 0 <= delta_index < delta_count:
            _fail(grid_reason, stage)
    return t_index, delta_index


def _validate_cell_geometry(cell, method_name, t_index, delta_index, t_edges, delta_edges, required):
    reason = "{}_cell_geometry_mismatch".format(method_name)
    stage = "{}_geometry".format(method_name)
    geometry = (
        ("t_low", t_index, t_edges, lambda index: t_edges[index]),
        ("t_high", t_index, t_edges, lambda index: t_edges[index + 1]),
        ("delta_low", delta_index, delta_edges, lambda index: delta_edges[index]),
        ("delta_high", delta_index, delta_edges, lambda index: delta_edges[index + 1]),
    )
    for key, index, _edges, expected in geometry:
        if key not in cell:
            if required:
                _fail(reason, stage)
            continue
        if index is None or not _serialized_equal(cell[key], expected(index)):
            _fail(reason, stage)


def _validate_method_cells(method, name, available, t_edges, delta_edges, host_state):
    cells = method["cells"]
    t_count = len(t_edges) - 1
    delta_count = len(delta_edges) - 1
    if available and len(cells) != t_count * delta_count:
        _fail("{}_cell_grid_mismatch".format(name), "{}_geometry".format(name))
    seen = set()
    for cell in cells:
        t_index, delta_index = _cell_indexes(
            cell, name, t_count, delta_count, required=available
        )
        _validate_cell_geometry(
            cell,
            name,
            t_index,
            delta_index,
            t_edges,
            delta_edges,
            required=available,
        )
        if name == "method_b" and "host_state" in cell and cell["host_state"] != host_state:
            _fail("method_b_cell_host_state_mismatch", "method_b_geometry")
        if available:
            coordinate = (t_index, delta_index)
            if coordinate in seen:
                _fail("{}_cell_grid_mismatch".format(name), "{}_geometry".format(name))
            seen.add(coordinate)


def _validate_parent_region_references(references, t_edges):
    stage = "method_b_geometry"
    t_count = len(t_edges) - 1
    for reference in references:
        if not isinstance(reference, Mapping) or "t_index" not in reference:
            continue
        t_index = reference["t_index"]
        if not _is_integer(t_index) or not 0 <= t_index < t_count:
            _fail("method_b_parent_reference_geometry_mismatch", stage)


def _available(snapshot, source_fingerprint):
    """Build the final detached D.1 payload from the validated snapshot."""
    return {
        "schema_version": COMPARISON_INPUT_SCHEMA_VERSION,
        "status": "available",
        "available": True,
        "reason": None,
        "diagnostic_stage": "complete",
        "source_checkpoint_schema_version": snapshot["schema_version"],
        "source_checkpoint_payload_fingerprint": source_fingerprint,
        "setting": snapshot["setting"],
        "canonical_t_edges": snapshot["canonical_t_edges"],
        "delta_edges": snapshot["delta_edges"],
        "phase_a": snapshot["phase_a"],
        "method_a": snapshot["method_a"],
        "method_b": snapshot["method_b"],
        "host_state_summary": snapshot["host_state_summary"],
        "non_authoritative": True,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }


def build_pion_hgcer_comparison_input_contract(checkpoint_payload):
    """Validate and detach the frozen Phase-C A/B comparison input.

    A successful result establishes only that the existing checkpoint can be
    consumed later.  It intentionally does not compare, rank, or modify either
    diagnostic method.
    """
    if not isinstance(checkpoint_payload, Mapping):
        return _unavailable("checkpoint_contract_invalid", "checkpoint_validation")
    try:
        encoded = _canonical_json(checkpoint_payload)
        snapshot = json.loads(encoded)
    except (TypeError, ValueError, OverflowError):
        return _unavailable("checkpoint_snapshot_invalid", "snapshot_validation")
    source_fingerprint = hashlib.sha256(encoded.encode("ascii")).hexdigest()
    try:
        _validate_checkpoint(snapshot)
        _validate_setting(snapshot["setting"])
        phase, phase_summary, phase_fingerprint, coordinate_fingerprint = _validate_phase_a(
            snapshot["phase_a"]
        )
        method_a, method_a_summary, method_a_available = _validate_method_structure(
            snapshot["method_a"], "method_a"
        )
        method_b, method_b_summary, method_b_available = _validate_method_structure(
            snapshot["method_b"], "method_b"
        )
        t_edges, delta_edges = _validate_geometry(
            snapshot, phase_summary, method_a, method_b
        )
        host_state = _validate_host_state(snapshot, phase, phase_summary)
        _validate_method_a_provenance(
            method_a,
            method_a_summary,
            method_a_available,
            phase_fingerprint,
            coordinate_fingerprint,
        )
        _validate_method_a_cells(
            method_a, method_a_available, t_edges, delta_edges, host_state
        )
        _validate_method_b_provenance(
            method_b,
            method_b_summary,
            method_b_available,
            phase_fingerprint,
            coordinate_fingerprint,
            host_state,
        )
        _validate_method_b_cells(
            method_b, method_b_available, t_edges, delta_edges, host_state
        )
        _validate_parent_region_references(method_b["parent_region_references"], t_edges)
        return _available(snapshot, source_fingerprint)
    except _ComparisonInputUnavailable as exc:
        return _unavailable(exc.reason, exc.stage)
    except Exception as exc:
        return _unavailable(
            "unexpected_comparison_input_build_failure",
            "unexpected_exception",
            exception=exc,
        )


def _validate_method_a_cells(method, available, t_edges, delta_edges, host_state):
    _validate_method_cells(
        method, "method_a", available, t_edges, delta_edges, host_state
    )


def _validate_method_b_cells(method, available, t_edges, delta_edges, host_state):
    _validate_method_cells(
        method, "method_b", available, t_edges, delta_edges, host_state
    )


__all__ = (
    "COMPARISON_INPUT_SCHEMA_VERSION",
    "build_pion_hgcer_comparison_input_contract",
)
