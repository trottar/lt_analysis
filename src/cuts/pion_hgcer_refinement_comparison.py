"""Read-only Phase-D comparison contracts for frozen HGCer diagnostics.

This module validates the detached Phase-C A/B comparison input and builds
non-authoritative Method-A and Method-B comparison representations. It does
not perform a Method-A/Method-B numerical comparison, construct a correction,
or own or mutate analysis objects.
"""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import json
import math


COMPARISON_INPUT_SCHEMA_VERSION = "pion_hgcer_ab_comparison_input/v1"
METHOD_A_COMPARISON_SCHEMA_VERSION = "pion_hgcer_method_a_comparison/v1"
METHOD_B_COMPARISON_SCHEMA_VERSION = "pion_hgcer_method_b_comparison/v1"

_CHECKPOINT_SCHEMA_VERSION = "pion_hgcer_refinement_checkpoint/v1"
_PHASE_A_SCHEMA_VERSION = "pion_hgcer_event_contract/v1"
_PHASE_A_FINGERPRINT_SCHEMA_VERSION = "pion_hgcer_event_contract_fingerprint/v2"
_METHOD_A_SCHEMA_VERSION = "pion_hgcer_method_a/v1"
_METHOD_B_SCHEMA_VERSION = "pion_hgcer_method_b/v1"
_HOST_STATES = {"proton_cleaned", "identity_no_proton_cleaning"}
_SOURCE_TARGET_STATE = "post_proton_noRF"
_METHOD_A_COMPARISON_METHOD = "method_a_same_t_comparison_representation"
_METHOD_B_COMPARISON_METHOD = "method_b_comparison_representation"
_METHOD_A_WILSON_Z_95 = 1.959963984540054
_METHOD_A_PARENT_MINIMUM_DELTA_CELLS = 2
_METHOD_A_SUPPORT_THRESHOLD_KEYS = (
    "supported_positive_count",
    "supported_low_count",
    "supported_control_count",
    "marginal_positive_count",
    "minimum_control_count_for_ratio",
    "minimum_low_count_for_ratio",
)
_METHOD_A_D2_CELL_FIELDS = (
    "t_index",
    "t_low",
    "t_high",
    "delta_index",
    "delta_low",
    "delta_high",
    "prompt_positive_count",
    "prompt_low_count",
    "prompt_control_count",
    "partition_closure_passed",
    "f_low",
    "f_low_low",
    "f_low_high",
    "R_low_control",
    "R_low_control_low",
    "R_low_control_high",
    "signed_R_low_control",
    "signed_R_low_control_sigma",
    "prompt_vs_signed_status",
    "nommcuts_vs_allcuts_status",
    "support_class",
    "method_A_status",
    "method_A_reason",
)
_METHOD_B_D3_CELL_FIELDS = (
    "t_index",
    "t_low",
    "t_high",
    "delta_index",
    "delta_low",
    "delta_high",
    "host_state",
    "host_record_count",
    "baseline_record_count",
    "region_consistency_status",
    "region_consistency_reason",
    "shape_chi2",
    "shape_ndf",
    "shape_chi2_ndf",
    "shape_max_abs_pull",
    "shape_usable_bin_count",
    "shape_status",
    "shape_reason",
    "candidate_L_B",
    "candidate_L_B_uncertainty",
    "candidate_L_B_status",
    "method_B_status",
    "method_B_reason",
)
_METHOD_B_D3_REGION_FIELDS = (
    "region_name",
    "available",
    "reason",
    "region_role",
    "window_source",
    "mm_low",
    "mm_high",
    "protected_signal_overlap",
    "support_status",
    "support_reason",
    "host_record_count",
    "host_yield",
    "host_sumw2",
    "host_neff",
    "baseline_record_count",
    "baseline_pion_yield",
    "baseline_pion_sumw2",
    "baseline_pion_neff",
    "baseline_pion_significance",
    "residual",
    "residual_sigma",
    "fractional_residual",
    "raw_ratio",
    "raw_ratio_sigma",
    "parent_reference_ratio",
    "parent_reference_sigma",
    "parent_relative_ratio",
    "parent_relative_sigma",
    "parent_relative_status",
    "parent_relative_reason",
)
_METHOD_B_D3_PARENT_REFERENCE_FIELDS = (
    "t_index",
    "t_low",
    "t_high",
    "region_name",
    "usable_delta_cell_count",
    "contributing_delta_indices",
    "combined_host_abs_support",
    "combined_host_sumw2",
    "combined_host_neff",
    "combined_baseline_abs_support",
    "combined_baseline_sumw2",
    "combined_baseline_neff",
    "parent_reference_ratio",
    "parent_reference_uncertainty",
    "parent_reference_status",
    "parent_reference_reason",
    "weighting",
)
_METHOD_B_CELL_STATUSES = {
    "available",
    "marginal",
    "internally_inconsistent",
    "shape_inconsistent",
    "unavailable",
}
_METHOD_B_CANDIDATE_STATUSES = {
    "available_multi_region",
    "single_region_only",
    "unavailable",
    "region_marginal",
    "region_inconsistent",
    "shape_poor_veto",
}
_METHOD_B_CONSISTENCY_STATUSES = {
    "region_consistent",
    "region_marginal",
    "region_inconsistent",
    "insufficient_regions",
}
_METHOD_B_SHAPE_STATUSES = {"good", "marginal", "poor", "unavailable"}


class _ComparisonInputUnavailable(Exception):
    """Expected validation failure carrying the public reason and stage."""

    def __init__(self, reason, stage):
        super().__init__(reason)
        self.reason = str(reason)
        self.stage = str(stage)


class _MethodAComparisonUnavailable(Exception):
    """Expected D.2 validation failure carrying the public reason and stage."""

    def __init__(self, reason, stage):
        super().__init__(reason)
        self.reason = str(reason)
        self.stage = str(stage)


class _MethodBComparisonUnavailable(Exception):
    """Expected D.3 validation failure carrying the public reason and stage."""

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


def _is_safe_scalar_metadata(value):
    """Accept only exact finite numeric or safe string setting scalars."""
    if _is_finite_number(value):
        return True
    if not isinstance(value, str) or not value or value != value.strip():
        return False
    if ".." in value or any(
        character.isspace()
        or ord(character) < 32
        or ord(character) == 127
        or character in "/\\:"
        for character in value
    ):
        return False
    return True


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
    if not _is_safe_scalar_metadata(setting["Q2"]) or not _is_safe_scalar_metadata(
        setting["W"]
    ):
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


def _prevalidate_source_scalar_metadata(checkpoint_payload):
    """Report malformed Q2/W metadata before JSON rejects non-finite values."""
    setting = checkpoint_payload.get("setting")
    if not isinstance(setting, Mapping):
        return
    for key in ("Q2", "W"):
        if key in setting and not _is_safe_scalar_metadata(setting[key]):
            _fail("setting_contract_invalid", "setting_validation")


def build_pion_hgcer_comparison_input_contract(checkpoint_payload):
    """Validate and detach the frozen Phase-C A/B comparison input.

    A successful result establishes only that the existing checkpoint can be
    consumed later.  It intentionally does not compare, rank, or modify either
    diagnostic method.
    """
    if not isinstance(checkpoint_payload, Mapping):
        return _unavailable("checkpoint_contract_invalid", "checkpoint_validation")
    try:
        _prevalidate_source_scalar_metadata(checkpoint_payload)
    except _ComparisonInputUnavailable as exc:
        return _unavailable(exc.reason, exc.stage)
    except Exception as exc:
        return _unavailable(
            "unexpected_comparison_input_build_failure",
            "unexpected_exception",
            exception=exc,
        )
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


def _method_a_comparison_unavailable(reason, stage, exception=None):
    """Return the detached, Method-A-only D.2 failure shape."""
    result = {
        "schema_version": METHOD_A_COMPARISON_SCHEMA_VERSION,
        "method": _METHOD_A_COMPARISON_METHOD,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "non_authoritative": True,
        "method_b_numerical_dependency": False,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    if exception is not None:
        result["exception_type"] = type(exception).__name__
        result["exception_message"] = str(exception)
    return result


def _method_a_comparison_fail(reason, stage):
    raise _MethodAComparisonUnavailable(reason, stage)


def _d2_mapping(value, reason, stage):
    if not isinstance(value, Mapping):
        _method_a_comparison_fail(reason, stage)
    return value


def _d2_require_keys(mapping, keys, reason, stage):
    if any(key not in mapping for key in keys):
        _method_a_comparison_fail(reason, stage)


def _d2_strict_edges(value, invalid_reason, stage):
    if not isinstance(value, list) or len(value) < 2:
        _method_a_comparison_fail(invalid_reason, stage)
    previous = None
    for edge in value:
        if not _is_finite_number(edge):
            _method_a_comparison_fail(invalid_reason, stage)
        if previous is not None and edge <= previous:
            _method_a_comparison_fail(invalid_reason, stage)
        previous = edge
    return value


def _method_a_wilson_interval(successes, total):
    """Return the exact Wilson-95 construction frozen by Method A."""
    if total <= 0 or successes < 0 or successes > total:
        return None, None
    fraction = float(successes) / float(total)
    z2 = _METHOD_A_WILSON_Z_95 * _METHOD_A_WILSON_Z_95
    denominator = 1.0 + z2 / total
    center = (fraction + z2 / (2.0 * total)) / denominator
    spread = (
        _METHOD_A_WILSON_Z_95
        * math.sqrt(
            fraction * (1.0 - fraction) / total + z2 / (4.0 * total * total)
        )
        / denominator
    )
    lower = 0.0 if successes == 0 else max(0.0, center - spread)
    upper = 1.0 if successes == total else min(1.0, center + spread)
    return lower, upper


def _method_a_ratio_from_fraction(value):
    if value is None or value >= 1.0:
        return None
    return float(value) / (1.0 - float(value))


def _method_a_support_class(positive, low, control, thresholds):
    if (
        positive >= thresholds["supported_positive_count"]
        and low >= thresholds["supported_low_count"]
        and control >= thresholds["supported_control_count"]
    ):
        return "supported"
    if (
        positive >= thresholds["marginal_positive_count"]
        and control >= thresholds["minimum_control_count_for_ratio"]
        and low >= thresholds["minimum_low_count_for_ratio"]
    ):
        return "marginal"
    return "unsupported"


def _validate_method_a_comparison_input(contract):
    """Validate and JSON-detach only the D.2 Method-A input projection."""
    stage = "comparison_input_validation"
    _d2_mapping(contract, "comparison_input_contract_invalid", stage)
    _d2_require_keys(
        contract,
        (
            "schema_version",
            "status",
            "available",
            "source_checkpoint_payload_fingerprint",
            "phase_a",
            "method_a",
            "canonical_t_edges",
            "delta_edges",
            "non_authoritative",
            "comparison_performed",
            "classification_performed",
            "production_objects_mutated",
            "refinement_applied",
        ),
        "comparison_input_contract_invalid",
        stage,
    )
    if contract["schema_version"] != COMPARISON_INPUT_SCHEMA_VERSION:
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)
    status = contract["status"]
    available = contract["available"]
    if not isinstance(status, str) or not isinstance(available, bool):
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)
    if status == "unavailable" and available is False:
        _method_a_comparison_fail("comparison_input_unavailable", stage)
    if status != "available" or available is not True:
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)
    if (
        contract["non_authoritative"] is not True
        or contract["comparison_performed"] is not False
        or contract["classification_performed"] is not False
        or contract["production_objects_mutated"] is not False
        or contract["refinement_applied"] is not False
    ):
        _method_a_comparison_fail(
            "comparison_input_authority_contract_invalid", "comparison_input_authority"
        )
    if not _nonempty_string(contract["source_checkpoint_payload_fingerprint"]):
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)

    # This projection deliberately omits the D.1 Method-B subtree.  The
    # serialized projection also provides all returned D.2 objects with their
    # own JSON-safe, non-aliased storage.
    projection = {
        "source_checkpoint_payload_fingerprint": contract[
            "source_checkpoint_payload_fingerprint"
        ],
        "phase_a": contract["phase_a"],
        "method_a": contract["method_a"],
        "canonical_t_edges": contract["canonical_t_edges"],
        "delta_edges": contract["delta_edges"],
    }
    try:
        return json.loads(_canonical_json(projection))
    except (TypeError, ValueError, OverflowError):
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)


def _validate_method_a_comparison_phase_a(phase):
    stage = "phase_a_provenance"
    phase = _d2_mapping(phase, "comparison_input_contract_invalid", stage)
    contract_fingerprint = phase.get("contract_fingerprint")
    if not _nonempty_string(contract_fingerprint):
        _method_a_comparison_fail("phase_a_fingerprint_missing", stage)
    coordinate_fingerprint = phase.get("coordinate_fingerprint")
    if not _nonempty_string(coordinate_fingerprint):
        _method_a_comparison_fail("coordinate_fingerprint_missing", stage)
    summary = phase.get("summary")
    if not isinstance(summary, Mapping):
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)
    if (
        summary.get("status") != "available"
        or summary.get("available") is not True
        or summary.get("contract_fingerprint") != contract_fingerprint
        or summary.get("coordinate_fingerprint") != coordinate_fingerprint
        or summary.get("schema_version") != _PHASE_A_SCHEMA_VERSION
        or summary.get("fingerprint_schema_version")
        != _PHASE_A_FINGERPRINT_SCHEMA_VERSION
        or summary.get("pion_closure_passed") is not True
        or summary.get("host_closure_passed") is not True
        or summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_restoration_applied") is not False
    ):
        _method_a_comparison_fail("comparison_input_contract_invalid", stage)
    return contract_fingerprint, coordinate_fingerprint


def _validate_method_a_comparison_thresholds(summary):
    stage = "method_a_summary"
    thresholds = summary.get("support_thresholds")
    if not isinstance(thresholds, Mapping) or set(thresholds) != set(
        _METHOD_A_SUPPORT_THRESHOLD_KEYS
    ):
        _method_a_comparison_fail("method_a_support_thresholds_invalid", stage)
    for key in _METHOD_A_SUPPORT_THRESHOLD_KEYS:
        value = thresholds[key]
        if not _is_integer(value) or value < 0:
            _method_a_comparison_fail("method_a_support_thresholds_invalid", stage)
    return dict(thresholds)


def _validate_method_a_comparison_summary(
    method_a, phase_fingerprint, coordinate_fingerprint, t_edges, delta_edges
):
    stage = "method_a_provenance"
    method_a = _d2_mapping(method_a, "method_a_summary_invalid", stage)
    if method_a.get("status") == "unavailable" and method_a.get("available") is False:
        _method_a_comparison_fail("method_a_unavailable", stage)
    if method_a.get("status") != "available" or method_a.get("available") is not True:
        _method_a_comparison_fail("method_a_summary_invalid", stage)
    fingerprint = method_a.get("fingerprint")
    if not _nonempty_string(fingerprint):
        _method_a_comparison_fail("method_a_fingerprint_missing", stage)
    summary = method_a.get("summary")
    if not isinstance(summary, Mapping):
        _method_a_comparison_fail("method_a_summary_invalid", "method_a_summary")
    if (
        summary.get("schema_version") != _METHOD_A_SCHEMA_VERSION
        or summary.get("status") != "available"
        or summary.get("available") is not True
        or summary.get("fingerprint") != fingerprint
        or summary.get("phase_a_contract_fingerprint") != phase_fingerprint
        or summary.get("coordinate_fingerprint") != coordinate_fingerprint
        or summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_ct_required") is not False
    ):
        _method_a_comparison_fail("method_a_summary_invalid", "method_a_summary")
    if summary.get("uncertainty_method") != "wilson_95_percent":
        _method_a_comparison_fail(
            "method_a_uncertainty_contract_invalid", "method_a_summary"
        )
    if not _serialized_equal(summary.get("t_edges"), t_edges) or not _serialized_equal(
        summary.get("delta_edges"), delta_edges
    ):
        _method_a_comparison_fail("canonical_geometry_mismatch", "canonical_geometry")
    return fingerprint, _validate_method_a_comparison_thresholds(summary)


def _validate_method_a_comparison_geometry(cell, t_edges, delta_edges):
    t_index = cell["t_index"]
    delta_index = cell["delta_index"]
    if (
        not _is_integer(t_index)
        or not _is_integer(delta_index)
        or not 0 <= t_index < len(t_edges) - 1
        or not 0 <= delta_index < len(delta_edges) - 1
        or not _serialized_equal(cell["t_low"], t_edges[t_index])
        or not _serialized_equal(cell["t_high"], t_edges[t_index + 1])
        or not _serialized_equal(cell["delta_low"], delta_edges[delta_index])
        or not _serialized_equal(cell["delta_high"], delta_edges[delta_index + 1])
    ):
        _method_a_comparison_fail("canonical_geometry_mismatch", "canonical_geometry")


def _validate_method_a_comparison_counts(cell):
    for key in (
        "prompt_positive_count",
        "prompt_low_count",
        "prompt_control_count",
    ):
        if not _is_integer(cell[key]) or cell[key] < 0:
            _method_a_comparison_fail(
                "method_a_cell_count_contract_invalid", "method_a_cells"
            )
    if (
        cell["partition_closure_passed"] is not True
        or cell["prompt_positive_count"]
        != cell["prompt_low_count"] + cell["prompt_control_count"]
    ):
        _method_a_comparison_fail(
            "method_a_cell_count_contract_invalid", "method_a_cells"
        )


def _validate_method_a_comparison_available_ratio(cell):
    positive = cell["prompt_positive_count"]
    low = cell["prompt_low_count"]
    control = cell["prompt_control_count"]
    if positive <= 0 or control <= 0:
        _method_a_comparison_fail("method_a_cell_ratio_contract_invalid", "method_a_cells")
    expected_f_low = float(low) / float(positive)
    expected_f_low_low, expected_f_low_high = _method_a_wilson_interval(low, positive)
    expected_ratio = float(low) / float(control)
    expected_ratio_low = _method_a_ratio_from_fraction(expected_f_low_low)
    expected_ratio_high = _method_a_ratio_from_fraction(expected_f_low_high)
    values = (
        cell["f_low"],
        cell["f_low_low"],
        cell["f_low_high"],
        cell["R_low_control"],
        cell["R_low_control_low"],
        cell["R_low_control_high"],
    )
    if (
        any(not _is_finite_number(value) for value in values)
        or not 0.0 <= cell["f_low_low"] <= cell["f_low"] <= cell["f_low_high"] <= 1.0
        or not 0.0
        <= cell["R_low_control_low"]
        <= cell["R_low_control"]
        <= cell["R_low_control_high"]
        or not _serialized_equal(cell["f_low"], expected_f_low)
        or not _serialized_equal(cell["f_low_low"], expected_f_low_low)
        or not _serialized_equal(cell["f_low_high"], expected_f_low_high)
        or not _serialized_equal(cell["R_low_control"], expected_ratio)
        or not _serialized_equal(cell["R_low_control_low"], expected_ratio_low)
        or not _serialized_equal(cell["R_low_control_high"], expected_ratio_high)
    ):
        _method_a_comparison_fail("method_a_cell_ratio_contract_invalid", "method_a_cells")


def _validate_method_a_comparison_cells(method_a, thresholds, t_edges, delta_edges):
    stage = "method_a_cells"
    cells = method_a.get("cells")
    if not isinstance(cells, list):
        _method_a_comparison_fail("method_a_cells_invalid", stage)
    expected_count = (len(t_edges) - 1) * (len(delta_edges) - 1)
    if len(cells) != expected_count:
        _method_a_comparison_fail("method_a_cell_count_contract_invalid", stage)
    selected_cells = []
    seen = set()
    for cell in cells:
        if not isinstance(cell, Mapping) or any(
            key not in cell for key in _METHOD_A_D2_CELL_FIELDS
        ):
            _method_a_comparison_fail("method_a_cells_invalid", stage)
        _validate_method_a_comparison_geometry(cell, t_edges, delta_edges)
        coordinate = (cell["t_index"], cell["delta_index"])
        if coordinate in seen:
            _method_a_comparison_fail("method_a_cell_count_contract_invalid", stage)
        seen.add(coordinate)
        _validate_method_a_comparison_counts(cell)
        support_class = _method_a_support_class(
            cell["prompt_positive_count"],
            cell["prompt_low_count"],
            cell["prompt_control_count"],
            thresholds,
        )
        if cell["support_class"] != support_class or cell["method_A_status"] not in {
            "available",
            "unavailable",
        }:
            _method_a_comparison_fail(
                "method_a_cell_status_contract_invalid", stage
            )
        expected_status = "unavailable" if support_class == "unsupported" else "available"
        if cell["method_A_status"] != expected_status:
            _method_a_comparison_fail(
                "method_a_cell_status_contract_invalid", stage
            )
        if cell["method_A_status"] == "available":
            _validate_method_a_comparison_available_ratio(cell)
        selected_cells.append({key: cell[key] for key in _METHOD_A_D2_CELL_FIELDS})
    return sorted(selected_cells, key=lambda cell: (cell["t_index"], cell["delta_index"]))


def _method_a_parent_reference(t_index, t_edges, cells, thresholds):
    """Build one same-t aggregate prompt-count parent reference."""
    contributors = [
        cell
        for cell in cells
        if cell["method_A_status"] == "available"
        and cell["support_class"] in {"supported", "marginal"}
    ]
    positive = sum(cell["prompt_positive_count"] for cell in contributors)
    low = sum(cell["prompt_low_count"] for cell in contributors)
    control = sum(cell["prompt_control_count"] for cell in contributors)
    closure = positive == low + control
    support_class = _method_a_support_class(positive, low, control, thresholds)
    f_low = float(low) / float(positive) if positive > 0 else None
    f_low_low, f_low_high = _method_a_wilson_interval(low, positive)
    ratio = float(low) / float(control) if control > 0 else None
    ratio_low = _method_a_ratio_from_fraction(f_low_low)
    ratio_high = _method_a_ratio_from_fraction(f_low_high)
    parent = {
        "t_index": t_index,
        "t_low": t_edges[t_index],
        "t_high": t_edges[t_index + 1],
        "contributing_delta_indices": [
            cell["delta_index"] for cell in contributors
        ],
        "contributing_delta_cell_count": len(contributors),
        "supported_contributing_cell_count": sum(
            cell["support_class"] == "supported" for cell in contributors
        ),
        "marginal_contributing_cell_count": sum(
            cell["support_class"] == "marginal" for cell in contributors
        ),
        "prompt_positive_count": positive,
        "prompt_low_count": low,
        "prompt_control_count": control,
        "partition_closure_passed": closure,
        "f_low": f_low,
        "f_low_low": f_low_low,
        "f_low_high": f_low_high,
        "R_low_control": ratio,
        "R_low_control_low": ratio_low,
        "R_low_control_high": ratio_high,
        "support_class": support_class,
        "parent_reference_status": "unavailable",
        "parent_reference_reason": None,
        "uncertainty_method": "wilson_95_percent_aggregated_prompt_counts",
        "minimum_contributing_delta_cells": _METHOD_A_PARENT_MINIMUM_DELTA_CELLS,
    }
    if len(contributors) < _METHOD_A_PARENT_MINIMUM_DELTA_CELLS:
        parent["parent_reference_reason"] = "insufficient_contributing_delta_cells"
    elif support_class == "unsupported":
        parent["parent_reference_reason"] = "parent_support_insufficient"
    elif not closure or low <= 0 or control <= 0 or ratio is None or ratio <= 0.0:
        parent["parent_reference_reason"] = "parent_reference_nonpositive"
    elif (
        not _is_finite_number(ratio_low)
        or not _is_finite_number(ratio_high)
        or ratio_low <= 0.0
        or ratio_low > ratio
        or ratio > ratio_high
    ):
        parent["parent_reference_reason"] = "parent_interval_invalid"
    elif support_class == "supported":
        parent["parent_reference_status"] = "available"
    else:
        parent["parent_reference_status"] = "marginal"
    return parent


def _method_a_comparison_cells(cells, parents):
    parent_by_t_index = {parent["t_index"]: parent for parent in parents}
    result = []
    for source_cell in cells:
        parent = parent_by_t_index[source_cell["t_index"]]
        cell = dict(source_cell)
        cell.update(
            {
                "parent_reference_R_low_control": parent["R_low_control"],
                "parent_reference_R_low_control_low": parent["R_low_control_low"],
                "parent_reference_R_low_control_high": parent["R_low_control_high"],
                "parent_reference_status": parent["parent_reference_status"],
                "method_a_comparison_candidate": None,
                "method_a_comparison_candidate_low": None,
                "method_a_comparison_candidate_high": None,
                "method_a_comparison_candidate_status": "unavailable",
                "method_a_comparison_candidate_reason": None,
                "candidate_interval_method": None,
                "candidate_covariance_treatment": None,
            }
        )
        if source_cell["method_A_status"] != "available":
            cell["method_a_comparison_candidate_reason"] = "cell_support_insufficient"
        elif parent["parent_reference_status"] == "unavailable":
            cell["method_a_comparison_candidate_reason"] = (
                parent["parent_reference_reason"] or "parent_reference_unavailable"
            )
        elif (
            parent["R_low_control"] is None
            or parent["R_low_control"] <= 0.0
        ):
            cell["method_a_comparison_candidate_reason"] = "parent_reference_nonpositive"
        elif (
            parent["R_low_control_low"] is None
            or parent["R_low_control_low"] <= 0.0
            or parent["R_low_control_high"] is None
        ):
            cell["method_a_comparison_candidate_reason"] = "parent_interval_invalid"
        else:
            cell["method_a_comparison_candidate"] = (
                source_cell["R_low_control"] / parent["R_low_control"]
            )
            cell["method_a_comparison_candidate_low"] = (
                source_cell["R_low_control_low"] / parent["R_low_control_high"]
            )
            cell["method_a_comparison_candidate_high"] = (
                source_cell["R_low_control_high"] / parent["R_low_control_low"]
            )
            cell["method_a_comparison_candidate_status"] = (
                "available"
                if (
                    source_cell["support_class"] == "supported"
                    and parent["parent_reference_status"] == "available"
                )
                else "marginal"
            )
            cell["candidate_interval_method"] = "ratio_envelope_from_wilson_bounds"
            cell["candidate_covariance_treatment"] = (
                "shared_parent_covariance_not_modeled"
            )
        result.append(cell)
    return result


def build_pion_hgcer_method_a_comparison(comparison_input_contract):
    """Build a detached same-t Method-A comparison representation.

    This consumes an already successful D.1 snapshot and uses no Method-B
    value.  It creates a descriptive parent-relative view of the frozen
    Method-A low/control observable without changing any upstream object.
    """
    if not isinstance(comparison_input_contract, Mapping):
        return _method_a_comparison_unavailable(
            "comparison_input_contract_invalid", "comparison_input_validation"
        )
    try:
        snapshot = _validate_method_a_comparison_input(comparison_input_contract)
        phase_fingerprint, coordinate_fingerprint = _validate_method_a_comparison_phase_a(
            snapshot["phase_a"]
        )
        t_edges = _d2_strict_edges(
            snapshot["canonical_t_edges"],
            "canonical_geometry_invalid",
            "canonical_geometry",
        )
        delta_edges = _d2_strict_edges(
            snapshot["delta_edges"], "canonical_geometry_invalid", "canonical_geometry"
        )
        method_a_fingerprint, thresholds = _validate_method_a_comparison_summary(
            snapshot["method_a"],
            phase_fingerprint,
            coordinate_fingerprint,
            t_edges,
            delta_edges,
        )
        cells = _validate_method_a_comparison_cells(
            snapshot["method_a"], thresholds, t_edges, delta_edges
        )
        parents = [
            _method_a_parent_reference(
                t_index,
                t_edges,
                [cell for cell in cells if cell["t_index"] == t_index],
                thresholds,
            )
            for t_index in range(len(t_edges) - 1)
        ]
        comparison_cells = _method_a_comparison_cells(cells, parents)
        source_method_a_payload_fingerprint = hashlib.sha256(
            _canonical_json(snapshot["method_a"]).encode("ascii")
        ).hexdigest()
        parent_definition = "same_t_aggregate_prompt_low_control_counts"
        uncertainty_definition = "wilson_95_percent_aggregated_prompt_counts"
        candidate_definition = "same_t_parent_ratio_with_ratio_envelope_bounds"
        fingerprint_inputs = {
            "schema_version": METHOD_A_COMPARISON_SCHEMA_VERSION,
            "method": _METHOD_A_COMPARISON_METHOD,
            "phase_a_contract_fingerprint": phase_fingerprint,
            "coordinate_fingerprint": coordinate_fingerprint,
            "method_a_fingerprint": method_a_fingerprint,
            "source_method_a_payload_fingerprint": source_method_a_payload_fingerprint,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "minimum_parent_delta_cells": _METHOD_A_PARENT_MINIMUM_DELTA_CELLS,
            "support_thresholds": thresholds,
            "parent_definition": parent_definition,
            "uncertainty_definition": uncertainty_definition,
            "candidate_definition": candidate_definition,
        }
        fingerprint = hashlib.sha256(
            _canonical_json(fingerprint_inputs).encode("ascii")
        ).hexdigest()
        return {
            "schema_version": METHOD_A_COMPARISON_SCHEMA_VERSION,
            "method": _METHOD_A_COMPARISON_METHOD,
            "status": "available",
            "available": True,
            "reason": None,
            "diagnostic_stage": "complete",
            "source_checkpoint_payload_fingerprint": snapshot[
                "source_checkpoint_payload_fingerprint"
            ],
            "source_method_a_payload_fingerprint": source_method_a_payload_fingerprint,
            "phase_a_contract_fingerprint": phase_fingerprint,
            "coordinate_fingerprint": coordinate_fingerprint,
            "method_a_fingerprint": method_a_fingerprint,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "support_thresholds": thresholds,
            "parent_definition": parent_definition,
            "uncertainty_definition": uncertainty_definition,
            "candidate_definition": candidate_definition,
            "parent_references": parents,
            "cells": comparison_cells,
            "fingerprint_inputs": fingerprint_inputs,
            "fingerprint": fingerprint,
            "non_authoritative": True,
            "method_b_numerical_dependency": False,
            "comparison_performed": False,
            "classification_performed": False,
            "production_objects_mutated": False,
            "refinement_applied": False,
        }
    except _MethodAComparisonUnavailable as exc:
        return _method_a_comparison_unavailable(exc.reason, exc.stage)
    except Exception as exc:
        return _method_a_comparison_unavailable(
            "unexpected_method_a_comparison_build_failure",
            "unexpected_exception",
            exception=exc,
        )


def _method_b_comparison_unavailable(reason, stage, exception=None):
    """Return the detached, Method-B-only D.3 failure shape."""
    result = {
        "schema_version": METHOD_B_COMPARISON_SCHEMA_VERSION,
        "method": _METHOD_B_COMPARISON_METHOD,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "non_authoritative": True,
        "method_a_numerical_dependency": False,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    if exception is not None:
        result["exception_type"] = type(exception).__name__
        result["exception_message"] = str(exception)
    return result


def _method_b_comparison_fail(reason, stage):
    raise _MethodBComparisonUnavailable(reason, stage)


def _d3_mapping(value, reason, stage):
    if not isinstance(value, Mapping):
        _method_b_comparison_fail(reason, stage)
    return value


def _d3_require_keys(mapping, keys, reason, stage):
    if any(key not in mapping for key in keys):
        _method_b_comparison_fail(reason, stage)


def _d3_strict_edges(value):
    if not isinstance(value, list) or len(value) < 2:
        _method_b_comparison_fail("canonical_geometry_invalid", "canonical_geometry")
    previous = None
    for edge in value:
        if not _is_finite_number(edge) or (
            previous is not None and edge <= previous
        ):
            _method_b_comparison_fail("canonical_geometry_invalid", "canonical_geometry")
        previous = edge
    return value


def _validate_method_b_comparison_input(contract):
    """Validate and detach only the D.3 Method-B input projection."""
    stage = "comparison_input_validation"
    _d3_mapping(contract, "comparison_input_contract_invalid", stage)
    _d3_require_keys(
        contract,
        (
            "schema_version",
            "status",
            "available",
            "source_checkpoint_payload_fingerprint",
            "phase_a",
            "method_b",
            "canonical_t_edges",
            "delta_edges",
            "host_state_summary",
            "non_authoritative",
            "comparison_performed",
            "classification_performed",
            "production_objects_mutated",
            "refinement_applied",
        ),
        "comparison_input_contract_invalid",
        stage,
    )
    if contract["schema_version"] != COMPARISON_INPUT_SCHEMA_VERSION:
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)
    if not isinstance(contract["status"], str) or not isinstance(
        contract["available"], bool
    ):
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)
    if contract["status"] == "unavailable" and contract["available"] is False:
        _method_b_comparison_fail("comparison_input_unavailable", stage)
    if contract["status"] != "available" or contract["available"] is not True:
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)
    if (
        contract["non_authoritative"] is not True
        or contract["comparison_performed"] is not False
        or contract["classification_performed"] is not False
        or contract["production_objects_mutated"] is not False
        or contract["refinement_applied"] is not False
    ):
        _method_b_comparison_fail(
            "comparison_input_authority_contract_invalid", "comparison_input_authority"
        )
    if not _nonempty_string(contract["source_checkpoint_payload_fingerprint"]):
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)
    projection = {
        "source_checkpoint_payload_fingerprint": contract[
            "source_checkpoint_payload_fingerprint"
        ],
        "phase_a": contract["phase_a"],
        "method_b": contract["method_b"],
        "canonical_t_edges": contract["canonical_t_edges"],
        "delta_edges": contract["delta_edges"],
        "host_state_summary": contract["host_state_summary"],
    }
    try:
        return json.loads(_canonical_json(projection))
    except (TypeError, ValueError, OverflowError):
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)


def _validate_method_b_comparison_phase_host(phase, host_summary):
    stage = "phase_a_provenance"
    phase = _d3_mapping(phase, "comparison_input_contract_invalid", stage)
    summary = _d3_mapping(phase.get("summary"), "comparison_input_contract_invalid", stage)
    host_summary = _d3_mapping(
        host_summary, "host_state_invalid", "host_provenance"
    )
    phase_fingerprint = phase.get("contract_fingerprint")
    if not _nonempty_string(phase_fingerprint):
        _method_b_comparison_fail("phase_a_fingerprint_missing", stage)
    coordinate_fingerprint = phase.get("coordinate_fingerprint")
    if not _nonempty_string(coordinate_fingerprint):
        _method_b_comparison_fail("coordinate_fingerprint_missing", stage)
    if (
        summary.get("status") != "available"
        or summary.get("available") is not True
        or summary.get("schema_version") != _PHASE_A_SCHEMA_VERSION
        or summary.get("fingerprint_schema_version")
        != _PHASE_A_FINGERPRINT_SCHEMA_VERSION
        or summary.get("contract_fingerprint") != phase_fingerprint
        or summary.get("coordinate_fingerprint") != coordinate_fingerprint
        or summary.get("pion_closure_passed") is not True
        or summary.get("host_closure_passed") is not True
        or summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_restoration_applied") is not False
    ):
        _method_b_comparison_fail("comparison_input_contract_invalid", stage)
    host_state = phase.get("host_state")
    if (
        host_state not in _HOST_STATES
        or summary.get("host_state") != host_state
        or host_summary.get("phase_a_host_state") != host_state
        or host_summary.get("method_b_host_state") != host_state
    ):
        _method_b_comparison_fail("host_state_invalid", "host_provenance")
    source_target_state = phase.get("source_target_state")
    if (
        source_target_state != _SOURCE_TARGET_STATE
        or summary.get("source_target_state") != source_target_state
        or host_summary.get("source_target_state") != source_target_state
    ):
        _method_b_comparison_fail(
            "source_target_state_mismatch", "host_provenance"
        )
    return phase_fingerprint, coordinate_fingerprint, host_state, source_target_state, summary


def _validate_method_b_comparison_summary(
    method_b, phase_fingerprint, coordinate_fingerprint, host_state, t_edges, delta_edges
):
    stage = "method_b_provenance"
    method_b = _d3_mapping(method_b, "method_b_summary_invalid", stage)
    if method_b.get("status") == "unavailable" and method_b.get("available") is False:
        _method_b_comparison_fail("method_b_unavailable", stage)
    if method_b.get("status") != "available" or method_b.get("available") is not True:
        _method_b_comparison_fail("method_b_summary_invalid", stage)
    fingerprint = method_b.get("fingerprint")
    if not _nonempty_string(fingerprint):
        _method_b_comparison_fail("method_b_fingerprint_missing", stage)
    summary = _d3_mapping(method_b.get("summary"), "method_b_summary_invalid", stage)
    if (
        summary.get("schema_version") != _METHOD_B_SCHEMA_VERSION
        or summary.get("status") != "available"
        or summary.get("available") is not True
        or summary.get("fingerprint") != fingerprint
        or summary.get("phase_a_contract_fingerprint") != phase_fingerprint
        or summary.get("coordinate_fingerprint") != coordinate_fingerprint
        or summary.get("host_state") != host_state
        or summary.get("production_objects_mutated") is not False
        or summary.get("refinement_applied") is not False
        or summary.get("rf_ct_required") is not False
        or summary.get("interpolation_used") is not False
    ):
        _method_b_comparison_fail("method_b_summary_invalid", stage)
    if not _serialized_equal(summary.get("t_edges"), t_edges) or not _serialized_equal(
        summary.get("delta_edges"), delta_edges
    ):
        _method_b_comparison_fail("canonical_geometry_mismatch", "canonical_geometry")
    return method_b, fingerprint


def _validate_method_b_region(row):
    stage = "method_b_regions"
    if not isinstance(row, Mapping) or any(
        key not in row for key in _METHOD_B_D3_REGION_FIELDS
    ):
        _method_b_comparison_fail("method_b_region_contract_invalid", stage)
    if (
        not _nonempty_string(row["region_name"])
        or not isinstance(row["available"], bool)
        or not _is_finite_number(row["mm_low"])
        or not _is_finite_number(row["mm_high"])
        or row["mm_low"] >= row["mm_high"]
        or not isinstance(row["protected_signal_overlap"], bool)
        or row["support_status"] not in {"usable", "unavailable"}
        or row["parent_relative_status"] not in {"available", "unavailable"}
    ):
        _method_b_comparison_fail("method_b_region_status_invalid", stage)
    for key in ("host_record_count", "baseline_record_count"):
        if not _is_integer(row[key]) or row[key] < 0:
            _method_b_comparison_fail("method_b_region_contract_invalid", stage)
    if row["support_status"] == "usable" and (
        not _is_finite_number(row["raw_ratio"])
        or not _is_finite_number(row["raw_ratio_sigma"])
        or row["raw_ratio_sigma"] <= 0.0
    ):
        _method_b_comparison_fail("method_b_region_contract_invalid", stage)
    if row["parent_relative_status"] == "available" and (
        not _is_finite_number(row["parent_relative_ratio"])
        or row["parent_relative_ratio"] <= 0.0
        or not _is_finite_number(row["parent_relative_sigma"])
        or row["parent_relative_sigma"] <= 0.0
        or not _is_finite_number(row["parent_reference_ratio"])
        or row["parent_reference_ratio"] <= 0.0
        or not _is_finite_number(row["parent_reference_sigma"])
        or row["parent_reference_sigma"] < 0.0
    ):
        _method_b_comparison_fail(
            "method_b_parent_relative_contract_invalid", stage
        )
    return {key: row[key] for key in _METHOD_B_D3_REGION_FIELDS}


def _validate_method_b_comparison_cells(method_b, t_edges, delta_edges, host_state):
    stage = "method_b_cells"
    cells = method_b.get("cells")
    if not isinstance(cells, list):
        _method_b_comparison_fail("method_b_cells_invalid", stage)
    expected_count = (len(t_edges) - 1) * (len(delta_edges) - 1)
    if len(cells) != expected_count:
        _method_b_comparison_fail("method_b_cell_grid_invalid", stage)
    selected_cells = []
    seen = set()
    region_names = None
    for cell in cells:
        if not isinstance(cell, Mapping) or any(
            key not in cell for key in _METHOD_B_D3_CELL_FIELDS
        ) or "regions" not in cell:
            _method_b_comparison_fail("method_b_cells_invalid", stage)
        t_index = cell["t_index"]
        delta_index = cell["delta_index"]
        if (
            not _is_integer(t_index)
            or not _is_integer(delta_index)
            or not 0 <= t_index < len(t_edges) - 1
            or not 0 <= delta_index < len(delta_edges) - 1
            or not _serialized_equal(cell["t_low"], t_edges[t_index])
            or not _serialized_equal(cell["t_high"], t_edges[t_index + 1])
            or not _serialized_equal(cell["delta_low"], delta_edges[delta_index])
            or not _serialized_equal(cell["delta_high"], delta_edges[delta_index + 1])
        ):
            _method_b_comparison_fail(
                "method_b_cell_geometry_invalid", "canonical_geometry"
            )
        coordinate = (t_index, delta_index)
        if coordinate in seen:
            _method_b_comparison_fail("method_b_cell_grid_invalid", stage)
        seen.add(coordinate)
        if cell["host_state"] != host_state:
            _method_b_comparison_fail(
                "method_b_cell_host_state_mismatch", stage
            )
        for key in ("host_record_count", "baseline_record_count"):
            if not _is_integer(cell[key]) or cell[key] < 0:
                _method_b_comparison_fail("method_b_cells_invalid", stage)
        if (
            cell["method_B_status"] not in _METHOD_B_CELL_STATUSES
            or cell["candidate_L_B_status"] not in _METHOD_B_CANDIDATE_STATUSES
            or cell["region_consistency_status"] not in _METHOD_B_CONSISTENCY_STATUSES
            or cell["shape_status"] not in _METHOD_B_SHAPE_STATUSES
        ):
            _method_b_comparison_fail("method_b_status_contract_invalid", stage)
        if (
            not _is_integer(cell["shape_ndf"])
            or cell["shape_ndf"] < 0
            or not _is_integer(cell["shape_usable_bin_count"])
            or cell["shape_usable_bin_count"] < 0
            or cell["shape_ndf"] != cell["shape_usable_bin_count"]
        ):
            _method_b_comparison_fail("method_b_shape_contract_invalid", stage)
        if cell["shape_status"] in {"good", "marginal", "poor"} and any(
            not _is_finite_number(cell[key])
            for key in ("shape_chi2", "shape_chi2_ndf", "shape_max_abs_pull")
        ):
            _method_b_comparison_fail("method_b_shape_contract_invalid", stage)
        candidate_status = cell["candidate_L_B_status"]
        if candidate_status == "available_multi_region":
            if (
                not _is_finite_number(cell["candidate_L_B"])
                or cell["candidate_L_B"] <= 0.0
                or not _is_finite_number(cell["candidate_L_B_uncertainty"])
                or cell["candidate_L_B_uncertainty"] <= 0.0
                or cell["method_B_status"] != "available"
                or cell["method_B_reason"] is not None
                or cell["region_consistency_status"] != "region_consistent"
                or cell["shape_status"] == "poor"
            ):
                _method_b_comparison_fail(
                    "method_b_candidate_contract_invalid", stage
                )
        elif cell["candidate_L_B"] is not None or cell["candidate_L_B_uncertainty"] is not None:
            _method_b_comparison_fail("method_b_candidate_contract_invalid", stage)
        if (
            cell["region_consistency_status"] == "region_inconsistent"
            and cell["method_B_status"] != "internally_inconsistent"
        ) or (
            cell["region_consistency_status"] == "region_marginal"
            and cell["method_B_status"] != "marginal"
        ) or (
            candidate_status == "shape_poor_veto"
            and (
                cell["method_B_status"] != "shape_inconsistent"
                or cell["shape_status"] != "poor"
            )
        ) or (
            candidate_status in {"single_region_only", "unavailable"}
            and cell["method_B_status"] != "unavailable"
        ):
            _method_b_comparison_fail("method_b_status_contract_invalid", stage)
        regions = cell["regions"]
        if not isinstance(regions, list) or not regions:
            _method_b_comparison_fail("method_b_region_contract_invalid", stage)
        selected_regions = [_validate_method_b_region(region) for region in regions]
        names = tuple(region["region_name"] for region in selected_regions)
        if len(names) != len(set(names)):
            _method_b_comparison_fail("method_b_region_contract_invalid", stage)
        if region_names is None:
            region_names = names
        elif names != region_names:
            _method_b_comparison_fail("method_b_region_contract_invalid", stage)
        selected = {key: cell[key] for key in _METHOD_B_D3_CELL_FIELDS}
        selected["regions"] = selected_regions
        selected_cells.append(selected)
    return selected_cells, region_names


def _validate_method_b_parent_references(references, region_names, t_edges, delta_edges):
    stage = "method_b_parent_references"
    if not isinstance(references, list):
        _method_b_comparison_fail("method_b_parent_references_invalid", stage)
    expected_count = (len(t_edges) - 1) * len(region_names)
    if len(references) != expected_count:
        _method_b_comparison_fail("method_b_parent_references_invalid", stage)
    selected_references = []
    lookup = {}
    for reference in references:
        if not isinstance(reference, Mapping) or any(
            key not in reference for key in _METHOD_B_D3_PARENT_REFERENCE_FIELDS
        ):
            _method_b_comparison_fail("method_b_parent_references_invalid", stage)
        t_index = reference["t_index"]
        region_name = reference["region_name"]
        if (
            not _is_integer(t_index)
            or not 0 <= t_index < len(t_edges) - 1
            or region_name not in region_names
            or not _serialized_equal(reference["t_low"], t_edges[t_index])
            or not _serialized_equal(reference["t_high"], t_edges[t_index + 1])
        ):
            _method_b_comparison_fail(
                "method_b_parent_reference_geometry_invalid", stage
            )
        deltas = reference["contributing_delta_indices"]
        if (
            not isinstance(deltas, list)
            or not _is_integer(reference["usable_delta_cell_count"])
            or reference["usable_delta_cell_count"] < 0
            or reference["usable_delta_cell_count"] != len(deltas)
            or len(deltas) != len(set(deltas))
            or any(
                not _is_integer(delta) or not 0 <= delta < len(delta_edges) - 1
                for delta in deltas
            )
        ):
            _method_b_comparison_fail(
                "method_b_parent_reference_geometry_invalid", stage
            )
        key = (t_index, region_name)
        if key in lookup:
            _method_b_comparison_fail(
                "method_b_parent_references_invalid", stage
            )
        if reference["weighting"] != "inverse_variance":
            _method_b_comparison_fail(
                "method_b_parent_reference_status_invalid", stage
            )
        if reference["parent_reference_status"] == "available":
            if (
                not _is_finite_number(reference["parent_reference_ratio"])
                or reference["parent_reference_ratio"] <= 0.0
                or not _is_finite_number(reference["parent_reference_uncertainty"])
                or reference["parent_reference_uncertainty"] < 0.0
                or reference["parent_reference_reason"] is not None
            ):
                _method_b_comparison_fail(
                    "method_b_parent_reference_status_invalid", stage
                )
        elif reference["parent_reference_status"] == "unavailable":
            if not _nonempty_string(reference["parent_reference_reason"]):
                _method_b_comparison_fail(
                    "method_b_parent_reference_status_invalid", stage
                )
        else:
            _method_b_comparison_fail(
                "method_b_parent_reference_status_invalid", stage
            )
        selected = {
            key: reference[key] for key in _METHOD_B_D3_PARENT_REFERENCE_FIELDS
        }
        lookup[key] = selected
        selected_references.append(selected)
    if set(lookup) != {
        (t_index, region_name)
        for t_index in range(len(t_edges) - 1)
        for region_name in region_names
    }:
        _method_b_comparison_fail("method_b_parent_references_invalid", stage)
    return selected_references, lookup


def _validate_method_b_parent_reference_linkage(cells, references):
    for cell in cells:
        for row in cell["regions"]:
            reference = references[(cell["t_index"], row["region_name"])]
            if not _serialized_equal(
                row["parent_reference_ratio"], reference["parent_reference_ratio"]
            ) or not _serialized_equal(
                row["parent_reference_sigma"],
                reference["parent_reference_uncertainty"],
            ):
                _method_b_comparison_fail(
                    "method_b_parent_reference_linkage_mismatch",
                    "method_b_parent_references",
                )


def _method_b_comparison_cells(cells):
    result = []
    for source_cell in cells:
        cell = dict(source_cell)
        if source_cell["candidate_L_B_status"] == "available_multi_region":
            cell.update(
                {
                    "method_b_comparison_candidate": source_cell["candidate_L_B"],
                    "method_b_comparison_candidate_uncertainty": source_cell[
                        "candidate_L_B_uncertainty"
                    ],
                    "method_b_comparison_candidate_status": "available_multi_region",
                    "method_b_comparison_candidate_reason": None,
                }
            )
        else:
            cell.update(
                {
                    "method_b_comparison_candidate": None,
                    "method_b_comparison_candidate_uncertainty": None,
                    "method_b_comparison_candidate_status": source_cell[
                        "candidate_L_B_status"
                    ],
                    "method_b_comparison_candidate_reason": source_cell[
                        "method_B_reason"
                    ],
                }
            )
        result.append(cell)
    return result


def build_pion_hgcer_method_b_comparison(comparison_input_contract):
    """Build a detached projection of the frozen Method-B diagnostics.

    D.3 validates representation integrity and copies Method-B-owned values;
    it neither derives a second Method-B result nor accesses Method A.
    """
    if not isinstance(comparison_input_contract, Mapping):
        return _method_b_comparison_unavailable(
            "comparison_input_contract_invalid", "comparison_input_validation"
        )
    try:
        snapshot = _validate_method_b_comparison_input(comparison_input_contract)
        (
            phase_fingerprint,
            coordinate_fingerprint,
            host_state,
            source_target_state,
            phase_summary,
        ) = _validate_method_b_comparison_phase_host(
            snapshot["phase_a"], snapshot["host_state_summary"]
        )
        t_edges = _d3_strict_edges(snapshot["canonical_t_edges"])
        delta_edges = _d3_strict_edges(snapshot["delta_edges"])
        if not _serialized_equal(phase_summary.get("canonical_t_edges"), t_edges) or not _serialized_equal(
            phase_summary.get("delta_edges"), delta_edges
        ):
            _method_b_comparison_fail(
                "canonical_geometry_mismatch", "canonical_geometry"
            )
        method_b, method_b_fingerprint = _validate_method_b_comparison_summary(
            snapshot["method_b"],
            phase_fingerprint,
            coordinate_fingerprint,
            host_state,
            t_edges,
            delta_edges,
        )
        cells, region_names = _validate_method_b_comparison_cells(
            method_b, t_edges, delta_edges, host_state
        )
        parents, parent_lookup = _validate_method_b_parent_references(
            method_b.get("parent_region_references"),
            region_names,
            t_edges,
            delta_edges,
        )
        _validate_method_b_parent_reference_linkage(cells, parent_lookup)
        comparison_cells = _method_b_comparison_cells(cells)
        source_method_b_payload_fingerprint = hashlib.sha256(
            _canonical_json(snapshot["method_b"]).encode("ascii")
        ).hexdigest()
        fingerprint_inputs = {
            "schema_version": METHOD_B_COMPARISON_SCHEMA_VERSION,
            "method": _METHOD_B_COMPARISON_METHOD,
            "phase_a_contract_fingerprint": phase_fingerprint,
            "coordinate_fingerprint": coordinate_fingerprint,
            "method_b_fingerprint": method_b_fingerprint,
            "source_method_b_payload_fingerprint": source_method_b_payload_fingerprint,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "host_state": host_state,
            "source_target_state": source_target_state,
            "candidate_projection_definition": "exact_native_candidate_L_B_copy",
            "regional_projection_definition": "compact_frozen_regional_evidence",
            "parent_reference_projection_definition": "compact_frozen_same_t_parent_references",
            "shape_projection_definition": "compact_frozen_shape_summary",
        }
        fingerprint = hashlib.sha256(
            _canonical_json(fingerprint_inputs).encode("ascii")
        ).hexdigest()
        return {
            "schema_version": METHOD_B_COMPARISON_SCHEMA_VERSION,
            "method": _METHOD_B_COMPARISON_METHOD,
            "status": "available",
            "available": True,
            "reason": None,
            "diagnostic_stage": "complete",
            "source_checkpoint_payload_fingerprint": snapshot[
                "source_checkpoint_payload_fingerprint"
            ],
            "source_method_b_payload_fingerprint": source_method_b_payload_fingerprint,
            "phase_a_contract_fingerprint": phase_fingerprint,
            "coordinate_fingerprint": coordinate_fingerprint,
            "method_b_fingerprint": method_b_fingerprint,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "host_state": host_state,
            "source_target_state": source_target_state,
            "parent_region_references": parents,
            "cells": comparison_cells,
            "fingerprint_inputs": fingerprint_inputs,
            "fingerprint": fingerprint,
            "non_authoritative": True,
            "method_a_numerical_dependency": False,
            "comparison_performed": False,
            "classification_performed": False,
            "production_objects_mutated": False,
            "refinement_applied": False,
        }
    except _MethodBComparisonUnavailable as exc:
        return _method_b_comparison_unavailable(exc.reason, exc.stage)
    except Exception as exc:
        return _method_b_comparison_unavailable(
            "unexpected_method_b_comparison_build_failure",
            "unexpected_exception",
            exception=exc,
        )


__all__ = (
    "COMPARISON_INPUT_SCHEMA_VERSION",
    "METHOD_A_COMPARISON_SCHEMA_VERSION",
    "METHOD_B_COMPARISON_SCHEMA_VERSION",
    "build_pion_hgcer_comparison_input_contract",
    "build_pion_hgcer_method_a_comparison",
    "build_pion_hgcer_method_b_comparison",
)
