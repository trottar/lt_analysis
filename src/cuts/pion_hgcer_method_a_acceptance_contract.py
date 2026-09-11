"""Detached dual-population Method-A acceptance event contract.

This JSON-only sidecar separates the uncensored positive-HGCer response
population used to train a future Method-A acceptance study from the frozen
physical pion-control population used downstream.  It constructs neither a
correction nor a production weight.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
import os


METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_event_contract/v2"
)
METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2"
)
METHOD_A_ACCEPTANCE_EVENT_CONTRACT_ARTIFACT_SCHEMA_VERSION = (
    "pion_hgcer_method_a_acceptance_event_contract_artifact/v2"
)

_PHASE_A_SCHEMA = "pion_hgcer_event_contract/v1"
_PHASE_A_FINGERPRINT_SCHEMA = "pion_hgcer_event_contract_fingerprint/v2"
_SOURCE_TARGET_STATE = "post_proton_noRF"
_METHOD_A_SCHEMA = "pion_hgcer_method_a/v1"
_METHOD_A_METHOD = "observed_positive_hgcer_response"
_PRIMARY_FEATURES = (
    "SHMS_delta",
    "SHMS_xptar",
    "SHMS_yptar",
    "P_hgcer_xAtCer",
    "P_hgcer_yAtCer",
)
_PROMPT_SOURCE = "prompt"
_LOW_RESPONSE_UPPER_BOUND = 2.0
_PARENT_CACHE_REQUIRED_FIELDS = (
    "source_label", "entry_index", "coordinate_fingerprint", "rf_state",
    "source_tree_name", "t_index", "adj_t", "adj_MM", "coefficient",
    "allcuts", "nommcuts", "ssdelta", "delta_index", "P_hgcer_npeSum",
    "P_hgcer_xAtCer", "P_hgcer_yAtCer", "ssxptar", "ssyptar",
    "hsxptar", "hsyptar", "phi_degrees", "phi_index",
)


class MethodAAcceptanceContractUnavailable(RuntimeError):
    """Expected frozen-input validation failure for the F.1 sidecar."""


def _mapping(value):
    return value if isinstance(value, Mapping) else {}


def _sequence(value):
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _materialize_edge_values(value, label, *, allow_none=False):
    """Detach an ordinary or indexable one-dimensional edge collection."""
    if value is None:
        if allow_none:
            return []
        raise MethodAAcceptanceContractUnavailable("{}_edges_invalid".format(label))
    if isinstance(value, (str, bytes, bytearray, bool)):
        raise MethodAAcceptanceContractUnavailable("{}_edges_invalid".format(label))
    if isinstance(value, (list, tuple)):
        return list(value)
    try:
        length = len(value)
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceContractUnavailable(
            "{}_edges_invalid".format(label)
        ) from exc
    if isinstance(length, bool) or not isinstance(length, int) or length < 0:
        raise MethodAAcceptanceContractUnavailable("{}_edges_invalid".format(label))
    try:
        return [value[index] for index in range(length)]
    except (TypeError, KeyError, IndexError, AttributeError, ValueError) as exc:
        raise MethodAAcceptanceContractUnavailable(
            "{}_edges_invalid".format(label)
        ) from exc


def _finite(value):
    if isinstance(value, bool):
        return None
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    return scalar if math.isfinite(scalar) else None


def _integer(value):
    if isinstance(value, bool):
        return None
    try:
        integer = int(value)
    except (TypeError, ValueError):
        return None
    return integer if value == integer else None


def _optional_finite(value):
    if value is None:
        return None
    return _finite(value)


def _nonempty_string(value):
    return value if isinstance(value, str) and value else None


def _json_copy(value, label):
    try:
        serialized = json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise MethodAAcceptanceContractUnavailable(
            "{}_not_json_safe".format(label)
        ) from exc
    return json.loads(serialized)


def _hash(value):
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _strict_edges(value, label):
    raw_edges = _materialize_edge_values(value, label)
    if len(raw_edges) < 2:
        raise MethodAAcceptanceContractUnavailable("{}_edges_invalid".format(label))
    edges = []
    for item in raw_edges:
        scalar = _finite(item)
        if scalar is None:
            raise MethodAAcceptanceContractUnavailable(
                "{}_edges_nonfinite".format(label)
            )
        edges.append(scalar)
    if any(edges[index] >= edges[index + 1] for index in range(len(edges) - 1)):
        raise MethodAAcceptanceContractUnavailable(
            "{}_edges_not_strictly_increasing".format(label)
        )
    return edges


def _edge_pair(edges, index, label):
    resolved = _integer(index)
    if resolved is None or not 0 <= resolved < len(edges) - 1:
        raise MethodAAcceptanceContractUnavailable("{}_index_invalid".format(label))
    return resolved, edges[resolved], edges[resolved + 1]


def _equal(left, right):
    if isinstance(left, bool) or isinstance(right, bool):
        return left is right
    if left is None or right is None:
        return left is right
    left_scalar, right_scalar = _finite(left), _finite(right)
    if left_scalar is not None and right_scalar is not None:
        return left_scalar == right_scalar
    return left == right


def _identity(record, label):
    source = _nonempty_string(record.get("source_label"))
    entry = _integer(record.get("entry_index"))
    if source is None or entry is None or entry < 0:
        raise MethodAAcceptanceContractUnavailable("{}_identity_invalid".format(label))
    return source, entry


def _phase_a_contract(value):
    phase = _mapping(value)
    required = (
        "schema_version", "fingerprint_schema_version", "status", "available",
        "contract_fingerprint", "pion_event_population_fingerprint",
        "coordinate_fingerprint", "canonical_t_edges", "delta_edges",
        "pion_records", "host_state", "source_target_state",
        "rf_restoration_applied", "immutable_record_contract",
        "production_objects_mutated", "refinement_applied",
    )
    if any(field not in phase for field in required):
        raise MethodAAcceptanceContractUnavailable("phase_a_contract_invalid")
    if (
        phase.get("schema_version") != _PHASE_A_SCHEMA
        or phase.get("fingerprint_schema_version") != _PHASE_A_FINGERPRINT_SCHEMA
        or phase.get("status") != "available"
        or phase.get("available") is not True
    ):
        raise MethodAAcceptanceContractUnavailable("phase_a_contract_unavailable")
    if (
        phase.get("rf_restoration_applied") is not False
        or phase.get("immutable_record_contract") is not True
        or phase.get("production_objects_mutated") is not False
        or phase.get("refinement_applied") is not False
        or phase.get("source_target_state") != _SOURCE_TARGET_STATE
    ):
        raise MethodAAcceptanceContractUnavailable("phase_a_authority_invalid")
    for field in (
        "contract_fingerprint", "pion_event_population_fingerprint",
        "coordinate_fingerprint", "host_state",
    ):
        if _nonempty_string(phase.get(field)) is None:
            raise MethodAAcceptanceContractUnavailable(
                "phase_a_{}_missing".format(field)
            )
    if not _sequence(phase.get("pion_records")):
        raise MethodAAcceptanceContractUnavailable("phase_a_pion_records_invalid")
    return phase


def _cache_parent_rows(cache):
    rows = cache.get("records")
    if not _sequence(rows):
        raise MethodAAcceptanceContractUnavailable("pion_cache_parent_records_invalid")
    index = {}
    for source in rows:
        row = _mapping(source)
        if any(field not in row for field in _PARENT_CACHE_REQUIRED_FIELDS):
            raise MethodAAcceptanceContractUnavailable(
                "pion_cache_parent_fields_missing"
            )
        if _nonempty_string(row.get("source_tree_name")) is None:
            raise MethodAAcceptanceContractUnavailable(
                "pion_cache_parent_provenance_invalid"
            )
        key = _identity(row, "pion_cache_parent")
        if key in index:
            raise MethodAAcceptanceContractUnavailable("pion_cache_parent_identity_duplicate")
        index[key] = row
    return index


def _column_values(value):
    try:
        return [_detach_child_scalar(value[index]) for index in range(len(value))]
    except MethodAAcceptanceContractUnavailable:
        raise
    except (TypeError, KeyError, IndexError):
        raise MethodAAcceptanceContractUnavailable("pion_cache_child_columns_invalid")


def _detach_child_scalar(value):
    """Detach one child-cache scalar from an optional NumPy-like wrapper."""
    return _detach_scalar(value, "pion_cache_child_scalar_invalid")


def _detach_scalar(value, reason):
    """Return one ordinary scalar without truth-testing an array-like value."""
    if value is None or type(value) in (str, bool, int, float):
        return value
    item = getattr(value, "item", None)
    if not callable(item):
        raise MethodAAcceptanceContractUnavailable(reason)
    try:
        detached = item()
    except Exception as exc:
        raise MethodAAcceptanceContractUnavailable(reason) from exc
    if detached is None or type(detached) in (str, bool, int, float):
        return detached
    raise MethodAAcceptanceContractUnavailable(reason)


def _cache_child_rows(cache):
    source_sections = _mapping(cache.get("child_event_cache"))
    if not source_sections:
        raise MethodAAcceptanceContractUnavailable("pion_cache_child_records_invalid")
    required = (
        "source_label", "entry_index", "coefficient", "coordinate_fingerprint",
        "adj_t", "adj_MM", "ssdelta", "P_hgcer_npeSum",
        "P_hgcer_xAtCer", "P_hgcer_yAtCer", "ssxptar", "ssyptar",
        "hsxptar", "hsyptar", "allcuts", "nommcuts", "t_index",
        "phi_index", "phi_degrees", "Q2", "W", "epsilon", "theta_cm_deg",
    )
    index = {}
    for section_name, section_value in source_sections.items():
        section = _mapping(section_value)
        if any(field not in section for field in required):
            raise MethodAAcceptanceContractUnavailable("pion_cache_child_fields_missing")
        columns = {field: _column_values(section[field]) for field in required}
        lengths = {len(column) for column in columns.values()}
        if len(lengths) != 1:
            raise MethodAAcceptanceContractUnavailable("pion_cache_child_columns_misaligned")
        for position in range(lengths.pop()):
            row = {field: columns[field][position] for field in required}
            source = _nonempty_string(row.get("source_label"))
            if source != str(section_name):
                raise MethodAAcceptanceContractUnavailable("pion_cache_child_source_mismatch")
            key = _identity(row, "pion_cache_child")
            if key in index:
                raise MethodAAcceptanceContractUnavailable("pion_cache_child_identity_duplicate")
            index[key] = row
    return index


def _record_geometry(record, t_edges, delta_edges):
    t_index, t_low, t_high = _edge_pair(
        t_edges, record.get("canonical_t_index"), "canonical_t"
    )
    if (
        not _equal(record.get("canonical_t_lower_edge"), t_low)
        or not _equal(record.get("canonical_t_upper_edge"), t_high)
    ):
        raise MethodAAcceptanceContractUnavailable("phase_a_t_geometry_mismatch")
    delta_index = record.get("delta_index")
    if delta_index is None:
        if (
            record.get("delta_lower_edge") is not None
            or record.get("delta_upper_edge") is not None
        ):
            raise MethodAAcceptanceContractUnavailable("phase_a_outside_delta_edges_invalid")
        return t_index, t_low, t_high, None, None, None
    resolved_delta, delta_low, delta_high = _edge_pair(
        delta_edges, delta_index, "delta"
    )
    if (
        not _equal(record.get("delta_lower_edge"), delta_low)
        or not _equal(record.get("delta_upper_edge"), delta_high)
    ):
        raise MethodAAcceptanceContractUnavailable("phase_a_delta_geometry_mismatch")
    return t_index, t_low, t_high, resolved_delta, delta_low, delta_high


def _required_feature(row, name):
    value = _finite(row.get(name))
    if value is None:
        raise MethodAAcceptanceContractUnavailable(
            "acceptance_feature_nonfinite:{}".format(name)
        )
    return value


def _parent_parity(phase, parent, t_index, delta_index):
    paired = (
        ("coordinate_fingerprint", "coordinate_fingerprint"),
        ("canonical_t_index", "t_index"),
        ("analysis_abs_t", "adj_t"),
        ("analysis_MM", "adj_MM"),
        ("signed_source_coefficient", "coefficient"),
        ("allcuts", "allcuts"),
        ("nommcuts", "nommcuts"),
        ("SHMS_delta", "ssdelta"),
        ("P_hgcer_npeSum", "P_hgcer_npeSum"),
        ("P_hgcer_xAtCer", "P_hgcer_xAtCer"),
        ("P_hgcer_yAtCer", "P_hgcer_yAtCer"),
    )
    for phase_name, parent_name in paired:
        if not _equal(phase.get(phase_name), parent.get(parent_name)):
            raise MethodAAcceptanceContractUnavailable(
                "phase_a_parent_parity_mismatch:{}".format(phase_name)
            )
    cache_delta = parent.get("delta_index")
    if (None if cache_delta is None else _integer(cache_delta)) != delta_index:
        raise MethodAAcceptanceContractUnavailable("phase_a_parent_delta_index_mismatch")
    if _integer(parent.get("t_index")) != t_index:
        raise MethodAAcceptanceContractUnavailable("phase_a_parent_t_index_mismatch")


def _parent_no_rf_provenance(phase, parent):
    """Close the frozen Phase-A to authoritative-cache no-RF link exactly."""
    phase_no_rf = phase.get("noRF_provenance")
    parent_no_rf = parent.get("rf_state")
    if phase_no_rf != "noRF":
        raise MethodAAcceptanceContractUnavailable(
            "phase_a_noRF_provenance_invalid"
        )
    if parent_no_rf != "noRF":
        raise MethodAAcceptanceContractUnavailable(
            "pion_cache_parent_noRF_provenance_invalid"
        )
    if phase_no_rf != parent_no_rf:
        raise MethodAAcceptanceContractUnavailable(
            "phase_a_parent_noRF_provenance_mismatch"
        )


def _child_parity(parent, child, phi_index):
    paired = (
        ("coordinate_fingerprint", "coordinate_fingerprint"),
        ("t_index", "t_index"),
        ("coefficient", "coefficient"),
        ("adj_t", "adj_t"),
        ("adj_MM", "adj_MM"),
        ("ssdelta", "ssdelta"),
        ("P_hgcer_npeSum", "P_hgcer_npeSum"),
        ("P_hgcer_xAtCer", "P_hgcer_xAtCer"),
        ("P_hgcer_yAtCer", "P_hgcer_yAtCer"),
        ("ssxptar", "ssxptar"),
        ("ssyptar", "ssyptar"),
        ("hsxptar", "hsxptar"),
        ("hsyptar", "hsyptar"),
        ("allcuts", "allcuts"),
        ("nommcuts", "nommcuts"),
        ("phi_degrees", "phi_degrees"),
    )
    for parent_name, child_name in paired:
        if not _equal(parent.get(parent_name), child.get(child_name)):
            raise MethodAAcceptanceContractUnavailable(
                "parent_child_parity_mismatch:{}".format(parent_name)
            )
    if _integer(child.get("phi_index")) != phi_index:
        raise MethodAAcceptanceContractUnavailable("parent_child_phi_index_mismatch")


def _phi_metadata(parent, phi_edges, child_index):
    phi_degrees = _required_feature(parent, "phi_degrees")
    raw_index = parent.get("phi_index")
    if raw_index is None:
        if child_index is not None:
            raise MethodAAcceptanceContractUnavailable("outside_phi_child_present")
        return phi_degrees, None, None, None, "outside_phi"
    phi_index, phi_low, phi_high = _edge_pair(phi_edges, raw_index, "phi")
    if child_index != phi_index:
        raise MethodAAcceptanceContractUnavailable("inside_phi_child_missing")
    if not phi_low <= phi_degrees < phi_high and not (
        phi_index == len(phi_edges) - 2 and phi_degrees == phi_high
    ):
        raise MethodAAcceptanceContractUnavailable("stored_phi_assignment_mismatch")
    return phi_degrees, phi_index, phi_low, phi_high, "inside_phi"


def _boolean(value, label):
    if type(value) is not bool:
        raise MethodAAcceptanceContractUnavailable(label)
    return value


def _normalized_part1_row(value, label):
    row = _mapping(value)
    if not row:
        raise MethodAAcceptanceContractUnavailable("{}_invalid".format(label))
    return {
        str(name): _detach_scalar(item, "{}_scalar_invalid".format(label))
        for name, item in row.items()
    }


def _part1_index(rows, label):
    if not _sequence(rows) or not rows:
        raise MethodAAcceptanceContractUnavailable("{}_invalid".format(label))
    index = {}
    for value in rows:
        row = _normalized_part1_row(value, label)
        key = _identity(row, label)
        if key in index:
            raise MethodAAcceptanceContractUnavailable("{}_identity_duplicate".format(label))
        index[key] = row
    return index


def _part1_diagnostic(value, phase, t_edges, delta_edges):
    diagnostic = _mapping(value)
    if (
        diagnostic.get("status") != "available"
        or diagnostic.get("rf_restoration_applied") is not False
        or _nonempty_string(diagnostic.get("coordinate_fingerprint")) is None
    ):
        raise MethodAAcceptanceContractUnavailable("part1_diagnostic_unavailable")
    if diagnostic.get("coordinate_fingerprint") != phase["coordinate_fingerprint"]:
        raise MethodAAcceptanceContractUnavailable("part1_coordinate_fingerprint_mismatch")
    if _strict_edges(diagnostic.get("t_edges"), "part1_t") != t_edges:
        raise MethodAAcceptanceContractUnavailable("part1_t_geometry_mismatch")
    if _strict_edges(diagnostic.get("delta_edges"), "part1_delta") != delta_edges:
        raise MethodAAcceptanceContractUnavailable("part1_delta_geometry_mismatch")
    if _nonempty_string(diagnostic.get("config_fingerprint")) is None:
        raise MethodAAcceptanceContractUnavailable("part1_config_fingerprint_missing")
    provenance = _mapping(_mapping(diagnostic.get("source_provenance")).get("pion"))
    if not provenance or _mapping(provenance.get(_PROMPT_SOURCE)).get("rf_state") != "noRF":
        raise MethodAAcceptanceContractUnavailable("part1_pion_provenance_invalid")
    if any(_mapping(item).get("rf_state") != "noRF" for item in provenance.values()):
        raise MethodAAcceptanceContractUnavailable("part1_pion_provenance_invalid")
    records = _mapping(diagnostic.get("records"))
    sidecar = _mapping(diagnostic.get("phase_e_acceptance_records"))
    response = _part1_index(records.get("pion"), "part1_response_records")
    acceptance = _part1_index(sidecar.get("pion"), "part1_acceptance_records")
    if set(response) != set(acceptance):
        if set(response) - set(acceptance):
            raise MethodAAcceptanceContractUnavailable("part1_acceptance_match_missing")
        raise MethodAAcceptanceContractUnavailable("part1_acceptance_match_extra")
    return diagnostic, response, acceptance


def _part1_geometry(row, t_edges, delta_edges):
    t_index, t_low, t_high = _edge_pair(
        t_edges, row.get("canonical_t_index"), "part1_canonical_t"
    )
    delta_index, delta_low, delta_high = _edge_pair(
        delta_edges, row.get("delta_index"), "part1_delta"
    )
    return t_index, t_low, t_high, delta_index, delta_low, delta_high


def _part1_parity(response, acceptance):
    paired = (
        "coordinate_fingerprint", "analysis_t", "canonical_t_index", "ssdelta",
        "delta_index", "P_hgcer_npeSum", "allcuts", "nommcuts",
        "diagnostic_weight", "rf_applied_to_diagnostic",
    )
    for name in paired:
        if name not in response or name not in acceptance:
            raise MethodAAcceptanceContractUnavailable(
                "part1_response_acceptance_fields_missing:{}".format(name)
            )
        if not _equal(response[name], acceptance[name]):
            raise MethodAAcceptanceContractUnavailable(
                "part1_response_acceptance_parity_mismatch:{}".format(name)
            )
    if (
        response.get("coordinate_fingerprint") is None
        or response.get("rf_applied_to_diagnostic") is not False
        or acceptance.get("rf_applied_to_diagnostic") is not False
    ):
        raise MethodAAcceptanceContractUnavailable("part1_response_provenance_invalid")


def _method_a_cells(value, phase, t_edges, delta_edges):
    method = _mapping(value)
    required = (
        "schema_version", "method", "status", "available", "non_authoritative",
        "production_objects_mutated", "refinement_applied", "fingerprint",
        "event_population_fingerprint", "coordinate_fingerprint",
        "phase_a_contract_fingerprint", "t_edges", "delta_edges", "cells", "summary",
        "response_population_definition", "physical_control_definition",
        "low_response_definition",
    )
    if any(name not in method for name in required):
        raise MethodAAcceptanceContractUnavailable("method_a_contract_invalid")
    if (
        method.get("schema_version") != _METHOD_A_SCHEMA
        or method.get("method") != _METHOD_A_METHOD
        or method.get("status") != "available"
        or method.get("available") is not True
    ):
        raise MethodAAcceptanceContractUnavailable("method_a_contract_unavailable")
    if (
        method.get("non_authoritative") is not True
        or method.get("production_objects_mutated") is not False
        or method.get("refinement_applied") is not False
    ):
        raise MethodAAcceptanceContractUnavailable("method_a_authority_invalid")
    for name in ("fingerprint", "event_population_fingerprint", "coordinate_fingerprint", "phase_a_contract_fingerprint"):
        if _nonempty_string(method.get(name)) is None:
            raise MethodAAcceptanceContractUnavailable("method_a_{}_missing".format(name))
    if (
        method["coordinate_fingerprint"] != phase["coordinate_fingerprint"]
        or method["phase_a_contract_fingerprint"] != phase["contract_fingerprint"]
    ):
        raise MethodAAcceptanceContractUnavailable("method_a_provenance_mismatch")
    if _strict_edges(method["t_edges"], "method_a_t") != t_edges:
        raise MethodAAcceptanceContractUnavailable("method_a_t_geometry_mismatch")
    if _strict_edges(method["delta_edges"], "method_a_delta") != delta_edges:
        raise MethodAAcceptanceContractUnavailable("method_a_delta_geometry_mismatch")
    if (
        method["response_population_definition"] != "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0"
        or method["physical_control_definition"] != "P_hgcer_npeSum_gt_2"
        or method["low_response_definition"] != "0_lt_P_hgcer_npeSum_le_2"
    ):
        raise MethodAAcceptanceContractUnavailable("method_a_response_definition_invalid")
    if not _sequence(method["cells"]):
        raise MethodAAcceptanceContractUnavailable("method_a_cells_invalid")
    cells = {}
    for value in method["cells"]:
        row = _mapping(value)
        fields = (
            "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
            "prompt_positive_count", "prompt_low_count", "prompt_control_count",
            "partition_closure_passed",
        )
        if any(name not in row for name in fields):
            raise MethodAAcceptanceContractUnavailable("method_a_cell_fields_missing")
        t_index, t_low, t_high = _edge_pair(t_edges, row["t_index"], "method_a_cell_t")
        delta_index, delta_low, delta_high = _edge_pair(delta_edges, row["delta_index"], "method_a_cell_delta")
        if (
            not _equal(row["t_low"], t_low) or not _equal(row["t_high"], t_high)
            or not _equal(row["delta_low"], delta_low) or not _equal(row["delta_high"], delta_high)
        ):
            raise MethodAAcceptanceContractUnavailable("method_a_cell_geometry_mismatch")
        key = (t_index, delta_index)
        if key in cells:
            raise MethodAAcceptanceContractUnavailable("method_a_cell_identity_duplicate")
        cells[key] = row
    expected = {
        (t_index, delta_index)
        for t_index in range(len(t_edges) - 1)
        for delta_index in range(len(delta_edges) - 1)
    }
    if set(cells) != expected:
        raise MethodAAcceptanceContractUnavailable("method_a_cell_lattice_invalid")
    return method, cells


def _training_population(response_index, acceptance_index, method_cells, t_edges, delta_edges):
    rows = []
    audit = {
        "observed_nonpositive_response_record_count": 0,
        "observed_prompt_nommcuts_nonpositive_response_count": 0,
        "zero_or_nonpositive_included_in_training": False,
        "absolute_leakage_probability_claimed": False,
    }
    for key in sorted(response_index):
        response, acceptance = response_index[key], acceptance_index[key]
        _part1_parity(response, acceptance)
        t_index, t_low, t_high, delta_index, delta_low, delta_high = _part1_geometry(
            response, t_edges, delta_edges
        )
        acceptance_geometry = _part1_geometry(acceptance, t_edges, delta_edges)
        if acceptance_geometry != (t_index, t_low, t_high, delta_index, delta_low, delta_high):
            raise MethodAAcceptanceContractUnavailable("part1_response_acceptance_geometry_mismatch")
        npe = _required_feature(response, "P_hgcer_npeSum")
        source_label = _nonempty_string(response.get("source_label"))
        if source_label is None:
            raise MethodAAcceptanceContractUnavailable("part1_response_source_invalid")
        nommcuts = _boolean(response.get("nommcuts"), "part1_response_nommcuts_invalid")
        allcuts = _boolean(response.get("allcuts"), "part1_response_allcuts_invalid")
        if npe <= 0.0:
            audit["observed_nonpositive_response_record_count"] += 1
            if source_label == _PROMPT_SOURCE and nommcuts:
                audit["observed_prompt_nommcuts_nonpositive_response_count"] += 1
            continue
        if source_label != _PROMPT_SOURCE or nommcuts is not True:
            continue
        response_class = "low" if npe <= _LOW_RESPONSE_UPPER_BOUND else "control"
        rows.append({
            "source_label": key[0], "entry_index": key[1],
            "coordinate_fingerprint": response["coordinate_fingerprint"],
            "t_index": t_index, "t_low": t_low, "t_high": t_high,
            "SHMS_delta": _required_feature(response, "ssdelta"),
            "delta_index": delta_index, "delta_low": delta_low, "delta_high": delta_high,
            "P_hgcer_npeSum": npe, "response_class": response_class,
            "P_hgcer_xAtCer": _required_feature(response, "P_hgcer_xAtCer"),
            "P_hgcer_yAtCer": _required_feature(response, "P_hgcer_yAtCer"),
            "SHMS_xptar": _required_feature(acceptance, "ssxptar"),
            "SHMS_yptar": _required_feature(acceptance, "ssyptar"),
            "allcuts": allcuts, "nommcuts": nommcuts,
            "analysis_t": _required_feature(response, "analysis_t"),
            "analysis_MM": _required_feature(response, "analysis_MM"),
            "diagnostic_weight": _required_feature(response, "diagnostic_weight"),
            "Q2": _optional_finite(response.get("Q2")),
            "W": _optional_finite(response.get("W")),
            "epsilon": _optional_finite(response.get("epsilon")),
            "phi": _optional_finite(response.get("phi")),
        })
    rows.sort(key=lambda row: (row["source_label"], row["entry_index"]))
    by_cell = {}
    closure_cells = []
    for t_index in range(len(t_edges) - 1):
        for delta_index in range(len(delta_edges) - 1):
            local = [row for row in rows if (row["t_index"], row["delta_index"]) == (t_index, delta_index)]
            counts = {
                "positive_count": len(local),
                "low_count": sum(row["response_class"] == "low" for row in local),
                "control_count": sum(row["response_class"] == "control" for row in local),
            }
            if counts["positive_count"] != counts["low_count"] + counts["control_count"]:
                raise MethodAAcceptanceContractUnavailable("training_partition_closure_invalid")
            method_cell = method_cells[(t_index, delta_index)]
            if method_cell.get("partition_closure_passed") is not True:
                raise MethodAAcceptanceContractUnavailable("method_a_partition_closure_invalid")
            comparisons = (
                ("positive_count", "prompt_positive_count"),
                ("low_count", "prompt_low_count"),
                ("control_count", "prompt_control_count"),
            )
            for training_name, method_name in comparisons:
                if _integer(method_cell.get(method_name)) != counts[training_name]:
                    raise MethodAAcceptanceContractUnavailable(
                        "method_a_training_cell_closure_mismatch:{}:{}:{}".format(
                            t_index, delta_index, training_name
                        )
                    )
            by_cell[(t_index, delta_index)] = counts
            closure_cells.append({
                "t_index": t_index, "t_low": t_edges[t_index], "t_high": t_edges[t_index + 1],
                "delta_index": delta_index, "delta_low": delta_edges[delta_index],
                "delta_high": delta_edges[delta_index + 1], **counts,
                "method_a_partition_closure_passed": True,
                "method_a_closure_passed": True,
            })
    return rows, audit, closure_cells


def _training_summary(rows, audit, closure_cells, method):
    positive = sum(row["positive_count"] for row in closure_cells)
    low = sum(row["low_count"] for row in closure_cells)
    control = sum(row["control_count"] for row in closure_cells)
    if _integer(_mapping(method.get("summary")).get("prompt_positive_nommcuts_records")) != positive:
        raise MethodAAcceptanceContractUnavailable("method_a_training_setting_closure_mismatch")
    if positive != low + control:
        raise MethodAAcceptanceContractUnavailable("training_setting_partition_closure_invalid")
    return {
        "training_record_count": len(rows),
        "prompt_positive_count": positive,
        "prompt_low_count": low,
        "prompt_control_count": control,
        "partition_closure_passed": True,
        "method_a_closure_passed": True,
        "by_t_delta": closure_cells,
        "response_threshold_audit": audit,
    }


def _application_summary(records, t_count, phi_count, unmatched_parent, unmatched_child):
    per_t = [{
        "t_index": index, "record_count": 0, "inside_phi_count": 0,
        "outside_phi_count": 0, "prompt_record_count": 0,
        "non_prompt_record_count": 0,
    } for index in range(t_count)]
    per_t_phi = [{
        "t_index": t_index, "phi_index": phi_index, "record_count": 0,
        "prompt_record_count": 0, "non_prompt_record_count": 0,
    } for t_index in range(t_count) for phi_index in range(phi_count)]
    groups = {(row["t_index"], row["phi_index"]): row for row in per_t_phi}
    for record in records:
        group = per_t[record["t_index"]]
        group["record_count"] += 1
        group["{}_count".format(record["phi_status"])] += 1
        source_kind = "prompt" if record["source_label"] == _PROMPT_SOURCE else "non_prompt"
        group["{}_record_count".format(source_kind)] += 1
        if record["phi_index"] is not None:
            phi_group = groups[(record["t_index"], record["phi_index"])]
            phi_group["record_count"] += 1
            phi_group["{}_record_count".format(source_kind)] += 1
    return {
        "phase_a_record_count": len(records), "matched_parent_record_count": len(records),
        "unmatched_parent_cache_record_count": len(unmatched_parent),
        "unmatched_child_cache_record_count": len(unmatched_child),
        "inside_phi_count": sum(row["inside_phi_count"] for row in per_t),
        "outside_phi_count": sum(row["outside_phi_count"] for row in per_t),
        "prompt_record_count": sum(row["prompt_record_count"] for row in per_t),
        "non_prompt_record_count": sum(row["non_prompt_record_count"] for row in per_t),
        "by_t": per_t, "by_t_phi": per_t_phi,
    }


def _application_population(phase, cache, t_edges, delta_edges, phi_edges):
    if str(cache.get("coordinate_fingerprint") or "") != phase["coordinate_fingerprint"]:
        raise MethodAAcceptanceContractUnavailable("pion_cache_coordinate_fingerprint_mismatch")
    if _strict_edges(cache.get("delta_edges"), "pion_cache_delta") != delta_edges:
        raise MethodAAcceptanceContractUnavailable("pion_cache_delta_geometry_mismatch")
    parent_index = _cache_parent_rows(cache)
    child_index = _cache_child_rows(cache)
    records, phase_keys, child_keys = [], set(), set()
    for source in phase["pion_records"]:
        phase_record = _mapping(source)
        key = _identity(phase_record, "phase_a_pion")
        if key in phase_keys:
            raise MethodAAcceptanceContractUnavailable("phase_a_pion_identity_duplicate")
        phase_keys.add(key)
        parent = parent_index.get(key)
        if parent is None:
            raise MethodAAcceptanceContractUnavailable("phase_a_parent_match_missing")
        t_index, t_low, t_high, delta_index, delta_low, delta_high = _record_geometry(
            phase_record, t_edges, delta_edges
        )
        _parent_no_rf_provenance(phase_record, parent)
        _parent_parity(phase_record, parent, t_index, delta_index)
        npe = _required_feature(phase_record, "P_hgcer_npeSum")
        if npe <= _LOW_RESPONSE_UPPER_BOUND:
            raise MethodAAcceptanceContractUnavailable("application_population_npe_not_physical_control")
        for feature_name in ("P_hgcer_xAtCer", "P_hgcer_yAtCer", "ssxptar", "ssyptar"):
            _required_feature(parent, feature_name)
        child = child_index.get(key)
        phi_degrees, phi_index, phi_low, phi_high, phi_status = _phi_metadata(
            parent, phi_edges, None if child is None else _integer(child.get("phi_index"))
        )
        if phi_index is not None:
            _child_parity(parent, child, phi_index)
            child_keys.add(key)
        elif child is not None:
            raise MethodAAcceptanceContractUnavailable("outside_phi_child_present")
        records.append({
            "source_label": key[0], "entry_index": key[1],
            "coordinate_fingerprint": phase["coordinate_fingerprint"],
            "t_index": t_index, "t_low": t_low, "t_high": t_high,
            "phi_degrees": phi_degrees, "phi_index": phi_index,
            "phi_low": phi_low, "phi_high": phi_high, "phi_status": phi_status,
            "SHMS_delta": _required_feature(phase_record, "SHMS_delta"),
            "delta_index": delta_index, "delta_low": delta_low, "delta_high": delta_high,
            "P_hgcer_npeSum": npe,
            "P_hgcer_xAtCer": _required_feature(parent, "P_hgcer_xAtCer"),
            "P_hgcer_yAtCer": _required_feature(parent, "P_hgcer_yAtCer"),
            "SHMS_xptar": _required_feature(parent, "ssxptar"),
            "SHMS_yptar": _required_feature(parent, "ssyptar"),
            "HMS_xptar": _optional_finite(parent.get("hsxptar")),
            "HMS_yptar": _optional_finite(parent.get("hsyptar")),
            "signed_source_coefficient": _required_feature(phase_record, "signed_source_coefficient"),
            "baseline_pion_weight_w0": _required_feature(phase_record, "baseline_pion_weight_w0"),
            "signed_baseline_event_contribution": _required_feature(phase_record, "signed_baseline_event_contribution"),
            "allcuts": _boolean(phase_record.get("allcuts"), "phase_a_allcuts_invalid"),
            "nommcuts": _boolean(phase_record.get("nommcuts"), "phase_a_nommcuts_invalid"),
            "analysis_t": _required_feature(phase_record, "analysis_abs_t"),
            "analysis_MM": _required_feature(phase_record, "analysis_MM"),
            "Q2": _optional_finite(None if child is None else child.get("Q2")),
            "W": _optional_finite(None if child is None else child.get("W")),
            "epsilon": _optional_finite(None if child is None else child.get("epsilon")),
            "theta_cm_deg": _optional_finite(None if child is None else child.get("theta_cm_deg")),
        })
    unmatched_parent = sorted(set(parent_index) - phase_keys)
    unmatched_child = sorted(set(child_index) - child_keys)
    return records, _application_summary(
        records, len(t_edges) - 1, len(phi_edges) - 1, unmatched_parent, unmatched_child
    ), unmatched_parent, unmatched_child


def _unavailable(reason, *, stage="validation", phase=None, phi_edges=None, diagnostic=None, method_a=None):
    source = _mapping(phase)
    part1 = _mapping(diagnostic)
    method = _mapping(method_a)
    def unavailable_edges(value, label):
        try:
            return [float(item) for item in _materialize_edge_values(value, label, allow_none=True)]
        except (MethodAAcceptanceContractUnavailable, TypeError, ValueError):
            return []
    return {
        "schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION,
        "fingerprint_schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
        "status": "unavailable", "available": False, "reason": str(reason),
        "diagnostic_stage": str(stage), "non_authoritative": True,
        "production_objects_mutated": False, "refinement_applied": False,
        "production_application_performed": False, "event_application_performed": False,
        "method_b_numerical_dependency": False, "future_weight_adjustment_constructed": False,
        "phase_a_contract_fingerprint": source.get("contract_fingerprint"),
        "phase_a_pion_event_population_fingerprint": source.get("pion_event_population_fingerprint"),
        "method_a_fingerprint": method.get("fingerprint"),
        "method_a_event_population_fingerprint": method.get("event_population_fingerprint"),
        "part1_config_fingerprint": part1.get("config_fingerprint"),
        "coordinate_fingerprint": source.get("coordinate_fingerprint"),
        "host_state": source.get("host_state"), "source_target_state": source.get("source_target_state"),
        "t_edges": unavailable_edges(source.get("canonical_t_edges"), "canonical_t"),
        "delta_edges": unavailable_edges(source.get("delta_edges"), "delta"),
        "phi_edges": unavailable_edges(phi_edges, "phi"),
        "method_a_training_records": [], "application_records": [],
        "method_a_training_summary": {"by_t_delta": [], "response_threshold_audit": {}},
        "application_summary": {"by_t": [], "by_t_phi": []},
        "method_a_training_population_fingerprint": None,
        "application_population_fingerprint": None,
        "acceptance_feature_metadata_fingerprint": None,
        "application_child_assignment_projection_fingerprint": None,
        "fingerprint_inputs": {}, "fingerprint": None,
    }


def build_pion_hgcer_method_a_acceptance_event_contract(
    pion_hgcer_tdelta_diagnostic, pion_hgcer_method_a, phase_a_event_contract,
    pion_control_cache, *, phi_edges,
):
    """Build the detached F.1 v2 Method-A-training/application sidecar."""
    phase = _mapping(phase_a_event_contract)
    diagnostic = _mapping(pion_hgcer_tdelta_diagnostic)
    method = _mapping(pion_hgcer_method_a)
    try:
        phase = _phase_a_contract(phase)
        t_edges = _strict_edges(phase["canonical_t_edges"], "canonical_t")
        delta_edges = _strict_edges(phase["delta_edges"], "delta")
        resolved_phi_edges = _strict_edges(phi_edges, "phi")
        diagnostic, response_index, acceptance_index = _part1_diagnostic(
            diagnostic, phase, t_edges, delta_edges
        )
        method, method_cells = _method_a_cells(method, phase, t_edges, delta_edges)
        training_records, response_audit, closure_cells = _training_population(
            response_index, acceptance_index, method_cells, t_edges, delta_edges
        )
        training_summary = _training_summary(
            training_records, response_audit, closure_cells, method
        )
        cache = _mapping(pion_control_cache)
        if not cache:
            raise MethodAAcceptanceContractUnavailable("pion_cache_invalid")
        application_records, application_summary, unmatched_parent, unmatched_child = _application_population(
            phase, cache, t_edges, delta_edges, resolved_phi_edges
        )
        feature_metadata = {
            "primary_acceptance_features": list(_PRIMARY_FEATURES),
            "training_population": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
            "training_low_definition": "0_lt_P_hgcer_npeSum_le_2",
            "training_control_definition": "P_hgcer_npeSum_gt_2",
            "application_population": "authoritative_physical_pion_control_P_hgcer_npeSum_gt_2",
            "parent_coordinate": "canonical_t",
            "downstream_yield_coordinates": ["canonical_t", "canonical_phi"],
            "phi_is_training_feature": False, "method_b_numerical_dependency": False,
            "probability_map_constructed": False, "weight_adjustment_constructed": False,
            "future_normalization_policy": "future_parent_t_only_no_tphi_child_renormalization",
            "absolute_leakage_probability_claimed": False,
        }
        training_fingerprint = _hash(training_records)
        application_fingerprint = _hash(application_records)
        feature_fingerprint = _hash(feature_metadata)
        child_projection = [{
            "source_label": row["source_label"], "entry_index": row["entry_index"],
            "t_index": row["t_index"], "phi_index": row["phi_index"],
            "phi_low": row["phi_low"], "phi_high": row["phi_high"],
            "phi_status": row["phi_status"],
        } for row in application_records]
        child_fingerprint = _hash(child_projection)
        fingerprint_inputs = {
            "schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION,
            "fingerprint_schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
            "phase_a_contract_fingerprint": phase["contract_fingerprint"],
            "phase_a_pion_event_population_fingerprint": phase["pion_event_population_fingerprint"],
            "method_a_fingerprint": method["fingerprint"],
            "method_a_event_population_fingerprint": method["event_population_fingerprint"],
            "part1_config_fingerprint": diagnostic["config_fingerprint"],
            "coordinate_fingerprint": phase["coordinate_fingerprint"],
            "host_state": phase["host_state"], "source_target_state": phase["source_target_state"],
            "t_edges": t_edges, "delta_edges": delta_edges, "phi_edges": resolved_phi_edges,
            "method_a_training_population_fingerprint": training_fingerprint,
            "application_population_fingerprint": application_fingerprint,
            "acceptance_feature_metadata_fingerprint": feature_fingerprint,
            "application_child_assignment_projection_fingerprint": child_fingerprint,
            "method_a_closure": closure_cells, "feature_metadata": feature_metadata,
        }
        return {
            "schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION,
            "fingerprint_schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
            "status": "available", "available": True, "reason": None,
            "diagnostic_stage": "complete", "non_authoritative": True,
            "production_objects_mutated": False, "refinement_applied": False,
            "production_application_performed": False, "event_application_performed": False,
            "method_b_numerical_dependency": False, "future_weight_adjustment_constructed": False,
            "phase_a_contract_fingerprint": phase["contract_fingerprint"],
            "phase_a_pion_event_population_fingerprint": phase["pion_event_population_fingerprint"],
            "method_a_fingerprint": method["fingerprint"],
            "method_a_event_population_fingerprint": method["event_population_fingerprint"],
            "part1_config_fingerprint": diagnostic["config_fingerprint"],
            "coordinate_fingerprint": phase["coordinate_fingerprint"],
            "host_state": phase["host_state"], "source_target_state": phase["source_target_state"],
            "t_edges": t_edges, "delta_edges": delta_edges, "phi_edges": resolved_phi_edges,
            "feature_metadata": feature_metadata,
            "method_a_training_records": training_records,
            "application_records": application_records,
            "method_a_training_summary": training_summary,
            "application_summary": application_summary,
            "unmatched_parent_cache_identities": [list(key) for key in unmatched_parent],
            "unmatched_child_cache_identities": [list(key) for key in unmatched_child],
            "method_a_training_population_fingerprint": training_fingerprint,
            "application_population_fingerprint": application_fingerprint,
            "acceptance_feature_metadata_fingerprint": feature_fingerprint,
            "application_child_assignment_projection_fingerprint": child_fingerprint,
            "fingerprint_inputs": fingerprint_inputs, "fingerprint": _hash(fingerprint_inputs),
        }
    except MethodAAcceptanceContractUnavailable as exc:
        return _unavailable(
            str(exc), phase=phase, phi_edges=phi_edges, diagnostic=diagnostic, method_a=method
        )
    except Exception:
        return _unavailable(
            "unexpected_method_a_acceptance_event_contract_build_failure",
            stage="unexpected_exception", phase=phase, phi_edges=phi_edges,
            diagnostic=diagnostic, method_a=method,
        )


def _filename_token(value):
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise ValueError("method_a_acceptance_contract_filename_token_invalid")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("method_a_acceptance_contract_filename_token_invalid")
    token = str(value)
    if (
        not token or token != token.strip() or any(char.isspace() for char in token)
        or any(char in token for char in "\\\\/:" ) or ".." in token
    ):
        raise ValueError("method_a_acceptance_contract_filename_token_invalid")
    return token


def pion_hgcer_method_a_acceptance_contract_filename(
    phi_setting, particle_type, kinematic_token, epsilon_filename_token
):
    """Return the deterministic F.1 contract basename."""
    phi = _filename_token(phi_setting)
    particle = _filename_token(particle_type).lower()
    kinematic = _filename_token(kinematic_token)
    epsilon = _filename_token(epsilon_filename_token).lower()
    if particle != "kaon" or epsilon not in {"lowe", "highe"}:
        raise ValueError("method_a_acceptance_contract_filename_setting_invalid")
    return "{}_{}_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(
        phi, particle, kinematic, epsilon
    )


def _artifact_setting(value):
    setting = _mapping(value)
    required = (
        "kinematic_token", "Q2", "W", "epsilon_setting",
        "epsilon_filename_token", "phi_setting", "particle_type",
    )
    if any(key not in setting for key in required):
        raise ValueError("method_a_acceptance_contract_artifact_setting_invalid")
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
    if result["particle_type"] != "kaon" or result["epsilon_filename_token"] not in {"lowe", "highe"}:
        raise ValueError("method_a_acceptance_contract_artifact_setting_invalid")
    return result


def build_pion_hgcer_method_a_acceptance_contract_artifact(*, setting, contract):
    """Wrap one detached F.1 contract with its resolved setting identity."""
    payload = _json_copy(contract, "contract")
    if payload.get("schema_version") != METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION:
        raise ValueError("method_a_acceptance_contract_artifact_contract_invalid")
    if (
        payload.get("non_authoritative") is not True
        or payload.get("production_objects_mutated") is not False
        or payload.get("refinement_applied") is not False
        or payload.get("production_application_performed") is not False
        or payload.get("event_application_performed") is not False
    ):
        raise ValueError("method_a_acceptance_contract_artifact_authority_invalid")
    return {
        "schema_version": METHOD_A_ACCEPTANCE_EVENT_CONTRACT_ARTIFACT_SCHEMA_VERSION,
        "setting": _artifact_setting(setting),
        "contract": payload,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "event_application_performed": False,
    }


def write_pion_hgcer_method_a_acceptance_contract_json(path, payload):
    """Write a deterministic, JSON-safe F.1 artifact with one trailing newline."""
    target = os.fspath(path)
    serializable = _json_copy(payload, "artifact")
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(serializable, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return target


__all__ = (
    "METHOD_A_ACCEPTANCE_EVENT_CONTRACT_ARTIFACT_SCHEMA_VERSION",
    "METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION",
    "METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION",
    "build_pion_hgcer_method_a_acceptance_event_contract",
    "build_pion_hgcer_method_a_acceptance_contract_artifact",
    "pion_hgcer_method_a_acceptance_contract_filename",
    "write_pion_hgcer_method_a_acceptance_contract_json",
)
