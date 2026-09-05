"""Detached D.6 through D.11 and E.2 procedure pages for the procedure PDF.

This module is presentation-only.  It receives already-built proton-cleaning
objects, clones only what it draws, and never rebuilds a fit, event lookup, or
weight.  Later procedure phases can append pages between the explicit PDF
open/close helpers without inheriting the technical diagnostic-PDF lifecycle.
"""

from __future__ import annotations

from array import array
import hashlib
import json
import math
import os
from collections.abc import Mapping, Sequence

from canonical_binning import find_canonical_bin
from pion_component_subtraction import simc_shape_pion_weight_from_value


D6_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d6/v1"
D7_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d7/v1"
D8_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d8/v1"
D9_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d9/v1"
D10_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d10/v1"
D11_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d11/v1"
E2_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_e2/v1"
FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX = "_full-background-subtraction"

_TIMING_T_METHOD = "timing_t_event_weight"
_CTIME_AERO_METHOD = "ctime_aero_event_weight"
_D9_SIDES = ("kaon", "pion")
_D9_METHOD_A_COMPARISON_SCHEMA = "pion_hgcer_method_a_comparison/v1"
_D10_PHASE_A_SCHEMA = "pion_hgcer_event_contract/v1"
_D10_METHOD_B_SCHEMA = "pion_hgcer_method_b/v1"
_D10_METHOD_B_COMPARISON_SCHEMA = "pion_hgcer_method_b_comparison/v1"
_D10_SOURCE_TARGET_STATE = "post_proton_noRF"
_D11_PHASE_D_CHECKPOINT_SCHEMA = "pion_hgcer_phase_d_checkpoint/v1"
_D11_METHOD_A_COMPARISON_SCHEMA = "pion_hgcer_method_a_comparison/v1"
_D11_METHOD_B_COMPARISON_SCHEMA = "pion_hgcer_method_b_comparison/v1"
_D11_AB_COMPARISON_SCHEMA = "pion_hgcer_ab_comparison/v1"
_D11_SOURCE_TARGET_STATE = "post_proton_noRF"
_D11_AVAILABILITY_LABELS = {
    "both_comparable": "Both methods available",
    "both_present_not_comparable": "Both present; ratio undefined",
    "a_only": "Method A only",
    "b_only": "Method B only",
    "neither_available": "Neither method available",
}


def _import_root():
    try:
        import ROOT
    except Exception:
        return None
    return ROOT


def _mapping(value):
    return value if isinstance(value, Mapping) else {}


def _strict_edges(value):
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return None
    try:
        edges = [float(edge) for edge in value]
    except (TypeError, ValueError):
        return None
    if len(edges) < 2 or any(not math.isfinite(edge) for edge in edges):
        return None
    if any(right <= left for left, right in zip(edges, edges[1:])):
        return None
    return edges


def _unavailable(reason):
    return {
        "schema_version": D6_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "method": None,
        "method_supported": False,
        "t_edges": [],
        "delta_edges": [],
        "per_t": [],
    }


def _d7_unavailable(reason):
    return {
        "schema_version": D7_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "method": None,
        "t_edges": [],
        "delta_edges": [],
        "per_t": [],
    }


def _d8_unavailable(reason):
    return {
        "schema_version": D8_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "t_edges": [],
        "delta_edges": [],
        "per_t": [],
    }


def _d9_unavailable(reason):
    return {
        "schema_version": D9_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "t_edges": [],
        "delta_edges": [],
        "coordinate_fingerprint": None,
        "method_a_thresholds": {"available": False},
        "per_t": [],
    }


def _d10_unavailable(reason):
    """Return the detached D.10 unavailable payload without source aliases."""
    return {
        "schema_version": D10_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "phase_a_contract_fingerprint": None,
        "coordinate_fingerprint": None,
        "method_b_fingerprint": None,
        "host_state": None,
        "host_label": None,
        "source_target_state": None,
        "t_edges": [],
        "delta_edges": [],
        "mm_edges": [],
        "mm_regions": [],
        "protected_regions": [],
        "d3_available": False,
        "d3_reason": None,
        "per_t": [],
    }


def _d11_unavailable(reason):
    """Return a detached unavailable D.11 payload without source aliases."""
    return {
        "schema_version": D11_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "source_checkpoint_payload_fingerprint": None,
        "phase_a_contract_fingerprint": None,
        "coordinate_fingerprint": None,
        "method_a_comparison_fingerprint": None,
        "method_b_comparison_fingerprint": None,
        "ab_comparison_fingerprint": None,
        "host_state": None,
        "source_target_state": None,
        "t_edges": [],
        "delta_edges": [],
        "per_t": [],
    }


def _timing_axis_label(branch):
    text = str(branch or "").strip()
    lowered = text.lower()
    if "ctime" in lowered or "coincidence" in lowered:
        return "Coincidence time [ns]"
    if "rf" in lowered:
        return "RF timing [ns]"
    return "Selected timing observable [ns]"


def _product_t_geometry(products):
    t_edges = []
    normalized = []
    for expected_index, product in enumerate(products):
        product = _mapping(product)
        try:
            t_index = int(product.get("t_index"))
        except (TypeError, ValueError):
            return None, None, "canonical_t_index_invalid"
        if t_index != expected_index:
            return None, None, "canonical_t_index_invalid"
        edges = _strict_edges(product.get("t_edges"))
        if edges is None or len(edges) != 2:
            return None, None, "canonical_t_edges_invalid"
        if t_edges and edges[0] != t_edges[-1]:
            return None, None, "canonical_t_edges_noncontiguous"
        if not t_edges:
            t_edges.append(edges[0])
        t_edges.append(edges[1])
        normalized.append((t_index, edges, product))
    if not normalized:
        return None, None, "canonical_t_products_missing"
    return t_edges, normalized, None


def _timing_t_pid_entry(result, delta_count, t_index):
    cells_by_delta = result.get("H_delta_timing_t_cells")
    if not isinstance(cells_by_delta, Sequence) or isinstance(cells_by_delta, (str, bytes)):
        return {"available": False, "reason": "timing_cell_histograms_missing"}
    if len(cells_by_delta) != delta_count:
        return {"available": False, "reason": "timing_cell_geometry_mismatch"}
    cells = []
    for delta_cells in cells_by_delta:
        if not isinstance(delta_cells, Sequence) or isinstance(delta_cells, (str, bytes)):
            return {"available": False, "reason": "timing_cell_histograms_invalid"}
        if t_index >= len(delta_cells) or delta_cells[t_index] is None:
            return {"available": False, "reason": "timing_cell_histogram_missing"}
        cells.append(delta_cells[t_index])
    return {
        "available": True,
        "kind": "timing_t_cells",
        "timing_axis_label": _timing_axis_label(result.get("selected_timing_branch")),
        "cell_histograms": tuple(cells),
        "aerogel_validation_available": bool(
            _mapping(result.get("diagnostics")).get("aerogel_vs_t_validation")
            and _mapping(result.get("_aerogel_vs_t_root_payload"))
        ),
    }


def _legacy_pid_entry(result):
    histogram = result.get("H_global_pid")
    if histogram is None:
        return {"available": False, "reason": "ctime_aerogel_pid_histogram_missing"}
    return {
        "available": True,
        "kind": "ctime_aero",
        "timing_axis_label": "Coincidence time [ns]",
        "histogram": histogram,
        "shared_across_t": True,
    }


def build_full_background_subtraction_d6_payload(cleaning_result, cleaning_application):
    """Select D.6 display inputs without recalculating proton cleaning.

    ROOT objects remain references because the renderer immediately clones them.
    Geometry and scalar metadata are copied into fresh lists/dictionaries so
    callers cannot mutate the resolved cleaning mappings through this result.
    """
    result = _mapping(cleaning_result)
    application = _mapping(cleaning_application)
    products = application.get("canonical_t_products")
    if not isinstance(products, Sequence) or isinstance(products, (str, bytes)):
        return _unavailable("canonical_t_products_missing")
    t_edges, normalized_products, reason = _product_t_geometry(products)
    if reason is not None:
        return _unavailable(reason)
    delta_edges = _strict_edges(result.get("delta_edges"))
    delta_geometry_available = delta_edges is not None
    if not delta_geometry_available:
        # The raw-MM page is independent of proton-delta geometry.  Keep it
        # available and make only delta-dependent pages unavailable.
        delta_edges = []

    method = str(result.get("method") or "").strip()
    method_supported = method in {_TIMING_T_METHOD, _CTIME_AERO_METHOD}
    if method == _TIMING_T_METHOD:
        weight_source = application.get("H_proton_weight_vs_delta_t")
    elif method == _CTIME_AERO_METHOD:
        weight_source = application.get("H_proton_weight_vs_delta_aero")
    else:
        weight_source = None

    per_t = []
    delta_count = len(delta_edges) - 1
    for t_index, product_edges, product in normalized_products:
        raw_source = _mapping(product.get("raw_targets")).get("h_mm_nosub")
        raw = {
            "available": raw_source is not None,
            "reason": None if raw_source is not None else "raw_mm_source_missing",
            "histogram": raw_source,
        }
        if method == _TIMING_T_METHOD:
            pid = (
                _timing_t_pid_entry(result, delta_count, t_index)
                if delta_geometry_available
                else {"available": False, "reason": "proton_delta_edges_invalid"}
            )
        elif method == _CTIME_AERO_METHOD:
            pid = _legacy_pid_entry(result)
        else:
            pid = {"available": False, "reason": "proton_cleaning_method_unsupported"}
        weight = {
            "available": bool(
                method_supported and delta_geometry_available and weight_source is not None
            ),
            "reason": (
                None
                if method_supported and delta_geometry_available and weight_source is not None
                else (
                    "proton_delta_edges_invalid"
                    if method_supported and not delta_geometry_available
                    else "proton_weight_map_missing"
                    if method_supported
                    else "proton_cleaning_method_unsupported"
                )
            ),
            "histogram": weight_source,
            "kind": method if method_supported else None,
        }
        per_t.append({
            "t_index": int(t_index),
            "t_low": float(product_edges[0]),
            "t_high": float(product_edges[1]),
            "raw_mm": raw,
            "pid": pid,
            "weight": weight,
        })

    return {
        "schema_version": D6_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "method": method,
        "method_supported": bool(method_supported),
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "delta_geometry_available": bool(delta_geometry_available),
        "per_t": per_t,
    }


def _d7_mm_source(product, target_key, reason):
    histogram = _mapping(product.get(target_key)).get("h_mm_nosub")
    return {
        "available": histogram is not None,
        "reason": None if histogram is not None else reason,
        "histogram": histogram,
    }


def _d7_empty_projection(reason, delta_count=0):
    return {
        "available": False,
        "reason": str(reason),
        "rows_by_delta": tuple(tuple() for _ in range(max(0, int(delta_count)))),
        "exclusions": {},
    }


def _d7_projection_groups(normalized_products, delta_edges):
    return {
        t_index: [list() for _ in range(len(delta_edges) - 1)]
        for t_index, _edges, _product in normalized_products
    }


def _d7_legacy_t_index(entry, t_edges):
    """Delegate legacy CTime/aerogel membership to the shared contract."""
    try:
        return int(find_canonical_bin(entry.get("adj_t"), t_edges))
    except Exception:
        return -1


def _build_d7_delta_projection(
    cleaning_result, prepared_source_bundle, normalized_products, t_edges, delta_edges
):
    prepared_sources = _mapping(prepared_source_bundle).get("prepared_sources")
    if not isinstance(prepared_sources, Mapping):
        return _d7_empty_projection("prepared_source_rows_missing", len(delta_edges) - 1)
    lookup = _mapping(cleaning_result).get("_prepared_event_weight_lookup")
    if not isinstance(lookup, Mapping):
        return _d7_empty_projection("frozen_applied_lookup_missing", len(delta_edges) - 1)

    method = str(_mapping(cleaning_result).get("method") or "").strip()
    groups = _d7_projection_groups(normalized_products, delta_edges)
    exclusions = {
        "prepared_rows_seen": 0,
        "nommcuts_rows": 0,
        "invalid_source_entry_index": 0,
        "missing_frozen_lookup_row": 0,
        "missing_mm_or_weight": 0,
        "missing_applied_factor": 0,
        "invalid_frozen_delta_index": 0,
        "invalid_canonical_t_index": 0,
        "legacy_t_membership_rows": 0,
        "selected_rows": 0,
    }
    for source_label, source_spec in prepared_sources.items():
        source_spec = _mapping(source_spec)
        coefficient = source_spec.get("coefficient")
        try:
            physical_weight = float(coefficient)
        except (TypeError, ValueError):
            physical_weight = float("nan")
        entries = source_spec.get("entries")
        if not isinstance(entries, Mapping):
            continue
        for entry_index, entry in entries.items():
            exclusions["prepared_rows_seen"] += 1
            entry = _mapping(entry)
            if not bool(entry.get("nommcuts")):
                continue
            exclusions["nommcuts_rows"] += 1
            try:
                numeric_entry_index = int(entry_index)
            except (TypeError, ValueError):
                exclusions["invalid_source_entry_index"] += 1
                continue
            frozen = _mapping(lookup.get("{}:{}".format(source_label, numeric_entry_index)))
            if not frozen:
                exclusions["missing_frozen_lookup_row"] += 1
                continue
            try:
                missing_mass = float(entry.get("adj_mm"))
            except (TypeError, ValueError):
                missing_mass = float("nan")
            if not (math.isfinite(missing_mass) and math.isfinite(physical_weight)):
                exclusions["missing_mm_or_weight"] += 1
                continue
            try:
                proton_probability = float(frozen.get("proton_weight"))
                cleaned_factor = float(frozen.get("cleaned_factor"))
            except (TypeError, ValueError):
                proton_probability = cleaned_factor = float("nan")
            if not (math.isfinite(proton_probability) and math.isfinite(cleaned_factor)):
                exclusions["missing_applied_factor"] += 1
                continue
            delta_index = frozen.get("delta_index")
            if isinstance(delta_index, bool) or not isinstance(delta_index, int):
                exclusions["invalid_frozen_delta_index"] += 1
                continue
            if not 0 <= delta_index < len(delta_edges) - 1:
                exclusions["invalid_frozen_delta_index"] += 1
                continue
            t_index = frozen.get("t_index")
            if isinstance(t_index, bool) or not isinstance(t_index, int):
                if method != _CTIME_AERO_METHOD:
                    exclusions["invalid_canonical_t_index"] += 1
                    continue
                t_index = _d7_legacy_t_index(entry, t_edges)
                exclusions["legacy_t_membership_rows"] += 1
            if t_index not in groups:
                exclusions["invalid_canonical_t_index"] += 1
                continue
            groups[t_index][delta_index].append({
                "missing_mass": float(missing_mass),
                "raw_weight": float(physical_weight),
                "proton_contribution": float(physical_weight * proton_probability),
                "cleaned_contribution": float(physical_weight * cleaned_factor),
            })
            exclusions["selected_rows"] += 1
    return {
        "available": True,
        "reason": None,
        "rows_by_t_delta": {
            t_index: tuple(tuple(row) for row in rows_by_delta)
            for t_index, rows_by_delta in groups.items()
        },
        "exclusions": exclusions,
    }


def build_full_background_subtraction_d7_payload(
    cleaning_result, cleaning_application, prepared_source_bundle
):
    """Select frozen D.7 MM inputs without rerunning proton cleaning.

    The canonical MM products remain references until the renderer clones
    them.  Row values are copied scalar records from prepared entries and the
    committed applied lookup; this builder never accesses source trees.
    """
    result = _mapping(cleaning_result)
    application = _mapping(cleaning_application)
    products = application.get("canonical_t_products")
    if not isinstance(products, Sequence) or isinstance(products, (str, bytes)):
        return _d7_unavailable("canonical_t_products_missing")
    t_edges, normalized_products, reason = _product_t_geometry(products)
    if reason is not None:
        return _d7_unavailable(reason)
    delta_edges = _strict_edges(result.get("delta_edges"))
    if delta_edges is None:
        delta_edges = []
        projection = _d7_empty_projection("proton_delta_edges_invalid")
    else:
        projection = _build_d7_delta_projection(
            result, prepared_source_bundle, normalized_products, t_edges, delta_edges
        )

    per_t = []
    for t_index, product_edges, product in normalized_products:
        raw_mm = _d7_mm_source(product, "raw_targets", "raw_mm_source_missing")
        proton_mm = _d7_mm_source(
            product, "proton_targets", "proton_mm_source_missing"
        )
        cleaned_mm = _d7_mm_source(
            product,
            "cleaned_targets_pre_rf",
            "cleaned_pre_rf_mm_source_missing",
        )
        rows_by_delta = tuple(
            (projection.get("rows_by_t_delta") or {}).get(
                t_index, tuple(tuple() for _ in range(len(delta_edges) - 1))
            )
        )
        delta_available = bool(
            projection.get("available") and raw_mm.get("available")
        )
        per_t.append({
            "t_index": int(t_index),
            "t_low": float(product_edges[0]),
            "t_high": float(product_edges[1]),
            "raw_mm": raw_mm,
            "proton_mm": proton_mm,
            "cleaned_mm": cleaned_mm,
            "delta_projection": {
                "available": delta_available,
                "reason": (
                    None
                    if delta_available
                    else projection.get("reason") or raw_mm.get("reason")
                ),
                "rows_by_delta": rows_by_delta,
                "exclusions": dict(projection.get("exclusions") or {}),
            },
        })
    return {
        "schema_version": D7_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "method": str(result.get("method") or "").strip(),
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "projection_exclusions": dict(projection.get("exclusions") or {}),
        "per_t": per_t,
    }


def _d8_histogram_source(application, key, reason):
    histogram = _mapping(application).get(key)
    return {
        "available": histogram is not None,
        "reason": None if histogram is not None else str(reason),
        "histogram": histogram,
    }


def _d8_final_application_reason(parent):
    status = _mapping(_mapping(parent).get("final_diagnostic_application_status"))
    return str(
        status.get("final_reason")
        or status.get("reason")
        or status.get("detail")
        or "final_parent_application_unavailable"
    )


def _d8_parent_cache_geometry(pion_parents, pion_control_cache):
    if not isinstance(pion_parents, Sequence) or isinstance(pion_parents, (str, bytes)):
        return None, None, "pion_parent_collection_missing"
    cache = _mapping(pion_control_cache)
    by_t = cache.get("by_t")
    if not isinstance(by_t, Sequence) or isinstance(by_t, (str, bytes)):
        return None, None, "pion_control_cache_missing"
    if not pion_parents or len(pion_parents) != len(by_t):
        return None, None, "pion_parent_cache_count_mismatch"

    normalized = []
    t_edges = []
    for expected_index, (parent, cache_entry) in enumerate(zip(pion_parents, by_t)):
        parent = _mapping(parent)
        cache_entry = _mapping(cache_entry)
        parent_index = parent.get("t_bin_index")
        cache_index = cache_entry.get("t_index")
        if (
            isinstance(parent_index, bool)
            or isinstance(cache_index, bool)
            or not isinstance(parent_index, int)
            or not isinstance(cache_index, int)
            or parent_index != expected_index
            or cache_index != expected_index
        ):
            return None, None, "pion_parent_cache_t_index_mismatch"
        parent_edges = _strict_edges(parent.get("t_edges"))
        cache_edges = _strict_edges(cache_entry.get("t_edges"))
        if (
            parent_edges is None
            or cache_edges is None
            or len(parent_edges) != 2
            or len(cache_edges) != 2
            or parent_edges != cache_edges
        ):
            return None, None, "pion_parent_cache_t_geometry_mismatch"
        if t_edges and parent_edges[0] != t_edges[-1]:
            return None, None, "pion_parent_cache_t_geometry_mismatch"
        if not t_edges:
            t_edges.append(parent_edges[0])
        t_edges.append(parent_edges[1])
        normalized.append((expected_index, parent_edges, parent, cache_entry))
    return t_edges, normalized, None


def _d8_empty_delta_projection(reason, delta_count=0):
    return {
        "available": False,
        "reason": str(reason),
        "rows_by_delta": tuple(tuple() for _ in range(max(0, int(delta_count)))),
        "exclusions": {},
        "closure": {
            "status": "unavailable",
            "content_passed": False,
            "error_passed": None,
            "coverage_complete": False,
        },
    }


def _d8_histogram_bin_count(histogram):
    try:
        count = int(histogram.GetNbinsX())
    except Exception:
        return None
    return count if count >= 1 else None


def _d8_compare_projection_bins(template, contents, variances):
    bin_count = _d8_histogram_bin_count(template)
    if bin_count is None or len(contents) != bin_count + 2:
        return {
            "available": False,
            "content_passed": False,
            "error_passed": None,
            "content_max_abs_difference": None,
            "error_max_abs_difference": None,
        }
    content_max_abs_difference = 0.0
    content_passed = True
    error_available = hasattr(template, "GetBinError")
    error_max_abs_difference = 0.0 if error_available else None
    error_passed = True if error_available else None
    for bin_index in range(bin_count + 2):
        try:
            reference = float(template.GetBinContent(bin_index))
        except Exception:
            return {
                "available": False,
                "content_passed": False,
                "error_passed": None,
                "content_max_abs_difference": None,
                "error_max_abs_difference": None,
            }
        comparison = float(contents[bin_index])
        difference = abs(reference - comparison)
        content_max_abs_difference = max(content_max_abs_difference, difference)
        if difference > 1.0e-12 * max(1.0, abs(reference), abs(comparison)):
            content_passed = False
        if error_available:
            try:
                reference_variance = float(template.GetBinError(bin_index)) ** 2
            except Exception:
                error_available = False
                error_max_abs_difference = None
                error_passed = None
                continue
            comparison_variance = float(variances[bin_index])
            variance_difference = abs(reference_variance - comparison_variance)
            error_max_abs_difference = max(error_max_abs_difference, variance_difference)
            if variance_difference > 1.0e-12 * max(
                1.0, abs(reference_variance), abs(comparison_variance)
            ):
                error_passed = False
    return {
        "available": True,
        "content_passed": bool(content_passed),
        "error_passed": error_passed,
        "content_max_abs_difference": float(content_max_abs_difference),
        "error_max_abs_difference": error_max_abs_difference,
    }


def _build_d8_delta_projection(final_application, cache_entry, delta_edges):
    application = _mapping(final_application)
    cache_t = _mapping(cache_entry)
    template = application.get("H_pion_subtraction_template_MM_nosub")
    reference = application.get("H_pion_control_model")
    weights = application.get("weights")
    records = cache_t.get("records")
    if template is None or reference is None or weights is None:
        return _d8_empty_delta_projection(
            "final_pion_template_or_weight_reference_missing", len(delta_edges) - 1
        )
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        return _d8_empty_delta_projection("pion_control_records_missing", len(delta_edges) - 1)
    bin_count = _d8_histogram_bin_count(template)
    if bin_count is None or not hasattr(template, "GetXaxis"):
        return _d8_empty_delta_projection("baseline_template_geometry_invalid", len(delta_edges) - 1)
    try:
        template_axis = template.GetXaxis()
    except Exception:
        return _d8_empty_delta_projection("baseline_template_geometry_invalid", len(delta_edges) - 1)

    rows_by_delta = [list() for _ in range(len(delta_edges) - 1)]
    contents = [0.0] * (bin_count + 2)
    variances = [0.0] * (bin_count + 2)
    exclusions = {
        "records_seen": 0,
        "non_nommcuts_records": 0,
        "invalid_frozen_delta_index": 0,
        "missing_record_fields": 0,
        "invalid_template_bin": 0,
        "assigned_records": 0,
        "nonzero_contributions": 0,
    }
    for record in records:
        exclusions["records_seen"] += 1
        record = _mapping(record)
        if not bool(record.get("nommcuts")):
            exclusions["non_nommcuts_records"] += 1
            continue
        delta_index = record.get("delta_index")
        if isinstance(delta_index, bool) or not isinstance(delta_index, int):
            exclusions["invalid_frozen_delta_index"] += 1
            continue
        if not 0 <= delta_index < len(rows_by_delta):
            exclusions["invalid_frozen_delta_index"] += 1
            continue
        try:
            missing_mass = float(record.get("adj_MM"))
            coefficient = float(record.get("coefficient"))
        except (TypeError, ValueError):
            missing_mass = coefficient = float("nan")
        if not (math.isfinite(missing_mass) and math.isfinite(coefficient)):
            exclusions["missing_record_fields"] += 1
            continue
        contribution = coefficient * simc_shape_pion_weight_from_value(
            missing_mass, reference, weights
        )
        if not math.isfinite(contribution):
            exclusions["missing_record_fields"] += 1
            continue
        try:
            template_bin = int(template_axis.FindBin(missing_mass))
        except Exception:
            template_bin = -1
        if not 0 <= template_bin <= bin_count + 1:
            exclusions["invalid_template_bin"] += 1
            continue
        rows_by_delta[delta_index].append({
            "missing_mass": float(missing_mass),
            "baseline_contribution": float(contribution),
        })
        contents[template_bin] += contribution
        variances[template_bin] += contribution * contribution
        exclusions["assigned_records"] += 1
        exclusions["nonzero_contributions"] += int(contribution != 0.0)

    comparison = _d8_compare_projection_bins(template, contents, variances)
    coverage_complete = exclusions["invalid_frozen_delta_index"] == 0
    numerical_match = bool(
        comparison.get("available")
        and comparison.get("content_passed")
        and comparison.get("error_passed") is not False
    )
    if not comparison.get("available"):
        available = False
        reason = "baseline_template_closure_unavailable"
        status = "unavailable"
    elif not coverage_complete:
        available = True
        reason = None
        status = "incomplete_frozen_delta_coverage"
    elif numerical_match:
        available = True
        reason = None
        status = "closed"
    else:
        available = False
        reason = "baseline_template_closure_mismatch"
        status = "mismatch"
    return {
        "available": available,
        "reason": reason,
        "rows_by_delta": tuple(tuple(rows) for rows in rows_by_delta),
        "exclusions": exclusions,
        "closure": {
            "status": status,
            "coverage_complete": coverage_complete,
            **comparison,
        },
    }


def build_full_background_subtraction_d8_payload(pion_parents, pion_control_cache):
    """Select final authoritative-pion display products without changing them."""
    t_edges, normalized, reason = _d8_parent_cache_geometry(
        pion_parents, pion_control_cache
    )
    if reason is not None:
        return _d8_unavailable(reason)
    cache = _mapping(pion_control_cache)
    delta_edges = _strict_edges(cache.get("delta_edges"))
    delta_geometry_available = delta_edges is not None
    if not delta_geometry_available:
        delta_edges = []

    per_t = []
    for t_index, parent_edges, parent, cache_entry in normalized:
        final_application = parent.get("final_diagnostic_application_result")
        final_reason = _d8_final_application_reason(parent)
        application_available = isinstance(final_application, Mapping)
        before = _d8_histogram_source(
            final_application,
            "H_MM_nosub_before_pion_subtraction",
            final_reason if not application_available else "before_pion_mm_source_missing",
        )
        baseline = _d8_histogram_source(
            final_application,
            "H_pion_subtraction_template_MM_nosub",
            final_reason if not application_available else "baseline_pion_mm_source_missing",
        )
        after = _d8_histogram_source(
            final_application,
            "H_MM_nosub_after_pion_subtraction",
            final_reason if not application_available else "after_pion_mm_source_missing",
        )
        if not application_available:
            projection = _d8_empty_delta_projection(final_reason)
        elif not delta_geometry_available:
            projection = _d8_empty_delta_projection("pion_control_delta_edges_invalid")
        else:
            projection = _build_d8_delta_projection(
                final_application, cache_entry, delta_edges
            )
        per_t.append({
            "t_index": int(t_index),
            "t_low": float(parent_edges[0]),
            "t_high": float(parent_edges[1]),
            "before_pion_mm": before,
            "baseline_pion_mm": baseline,
            "after_pion_mm": after,
            "delta_projection": {
                "available": bool(projection.get("available") and baseline.get("available")),
                "reason": (
                    projection.get("reason")
                    if not projection.get("available")
                    else baseline.get("reason")
                ),
                "rows_by_delta": tuple(projection.get("rows_by_delta") or ()),
                "exclusions": dict(projection.get("exclusions") or {}),
                "closure": dict(projection.get("closure") or {}),
            },
        })
    return {
        "schema_version": D8_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "delta_geometry_available": bool(delta_geometry_available),
        "per_t": per_t,
    }


def _d9_integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _d9_part1_rows(records, side, t_count, delta_count):
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        return None, "part1_{}_records_missing".format(side)
    per_t = [[] for _ in range(t_count)]
    for source_record in records:
        record = _mapping(source_record)
        if not record or record.get("side") != side:
            return None, "part1_{}_record_contract_invalid".format(side)
        t_index = record.get("canonical_t_index")
        delta_index = record.get("delta_index")
        if (
            not _d9_integer(t_index)
            or not _d9_integer(delta_index)
            or not 0 <= t_index < t_count
            or not 0 <= delta_index < delta_count
        ):
            return None, "part1_{}_record_membership_invalid".format(side)
        try:
            delta = float(record.get("ssdelta"))
            npe = float(record.get("P_hgcer_npeSum"))
            diagnostic_weight = float(record.get("diagnostic_weight"))
        except (TypeError, ValueError):
            return None, "part1_{}_record_scalar_invalid".format(side)
        if not all(math.isfinite(value) for value in (delta, npe, diagnostic_weight)):
            return None, "part1_{}_record_scalar_invalid".format(side)
        per_t[t_index].append({
            "t_index": int(t_index),
            "delta_index": int(delta_index),
            "side": str(side),
            "delta": delta,
            "npe": npe,
            "diagnostic_weight": diagnostic_weight,
        })
    return tuple(tuple(rows) for rows in per_t), None


def _d9_thresholds(method_a):
    method = _mapping(method_a)
    if method.get("status") != "available" or method.get("available") is not True:
        return {"available": False}, "method_a_unavailable"
    configuration = _mapping(method.get("configuration"))
    try:
        positive = float(configuration.get("positive_response_threshold"))
        low_upper = float(configuration.get("low_response_upper_threshold"))
    except (TypeError, ValueError):
        return {"available": False}, "method_a_thresholds_invalid"
    if not (math.isfinite(positive) and math.isfinite(low_upper)) or positive != 0.0 or low_upper != 2.0:
        return {"available": False}, "method_a_thresholds_invalid"
    return {
        "available": True,
        "positive_response_threshold": positive,
        "low_response_upper_threshold": low_upper,
    }, None


def _d9_method_a_relative_cells(method_a, comparison, t_edges, delta_edges, coordinate_fingerprint):
    """Validate and copy D.2 cells without rebuilding a same-t reference."""
    method = _mapping(method_a)
    projection = _mapping(comparison)
    if (
        method.get("status") != "available"
        or method.get("available") is not True
        or method.get("non_authoritative") is not True
        or method.get("production_objects_mutated") is not False
        or method.get("refinement_applied") is not False
    ):
        return None, "method_a_unavailable"
    if (
        not isinstance(method.get("fingerprint"), str)
        or not method.get("fingerprint")
        or method.get("coordinate_fingerprint") != coordinate_fingerprint
        or _strict_edges(method.get("t_edges")) != t_edges
        or _strict_edges(method.get("delta_edges")) != delta_edges
    ):
        return None, "method_a_provenance_mismatch"
    if (
        projection.get("schema_version") != _D9_METHOD_A_COMPARISON_SCHEMA
        or projection.get("status") != "available"
        or projection.get("available") is not True
        or projection.get("non_authoritative") is not True
        or projection.get("method_b_numerical_dependency") is not False
        or projection.get("comparison_performed") is not False
        or projection.get("classification_performed") is not False
        or projection.get("production_objects_mutated") is not False
        or projection.get("refinement_applied") is not False
        or projection.get("method_a_fingerprint") != method.get("fingerprint")
        or projection.get("coordinate_fingerprint") != coordinate_fingerprint
        or _strict_edges(projection.get("canonical_t_edges")) != t_edges
        or _strict_edges(projection.get("delta_edges")) != delta_edges
    ):
        return None, "method_a_comparison_provenance_mismatch"
    source_cells = projection.get("cells")
    expected_count = (len(t_edges) - 1) * (len(delta_edges) - 1)
    if not isinstance(source_cells, list) or len(source_cells) != expected_count:
        return None, "method_a_comparison_cell_grid_invalid"
    copied = [[] for _ in range(len(t_edges) - 1)]
    seen = set()
    for source_cell in source_cells:
        cell = _mapping(source_cell)
        required = (
            "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
            "method_a_comparison_candidate", "method_a_comparison_candidate_low",
            "method_a_comparison_candidate_high", "method_a_comparison_candidate_status",
            "support_class", "method_A_status", "method_A_reason",
            "method_a_comparison_candidate_reason",
        )
        if not cell or any(key not in cell for key in required):
            return None, "method_a_comparison_cell_contract_invalid"
        t_index = cell["t_index"]
        delta_index = cell["delta_index"]
        if (
            not _d9_integer(t_index)
            or not _d9_integer(delta_index)
            or not 0 <= t_index < len(t_edges) - 1
            or not 0 <= delta_index < len(delta_edges) - 1
            or cell["t_low"] != t_edges[t_index]
            or cell["t_high"] != t_edges[t_index + 1]
            or cell["delta_low"] != delta_edges[delta_index]
            or cell["delta_high"] != delta_edges[delta_index + 1]
            or (t_index, delta_index) in seen
        ):
            return None, "method_a_comparison_cell_geometry_invalid"
        seen.add((t_index, delta_index))
        status = cell["method_a_comparison_candidate_status"]
        if status not in {"available", "marginal", "unavailable"}:
            return None, "method_a_comparison_cell_status_invalid"
        if (
            cell["support_class"] not in {"supported", "marginal", "unsupported"}
            or cell["method_A_status"] not in {"available", "unavailable"}
        ):
            return None, "method_a_comparison_cell_status_invalid"
        candidate = cell["method_a_comparison_candidate"]
        low = cell["method_a_comparison_candidate_low"]
        high = cell["method_a_comparison_candidate_high"]
        if status in {"available", "marginal"}:
            try:
                candidate, low, high = float(candidate), float(low), float(high)
            except (TypeError, ValueError):
                return None, "method_a_comparison_candidate_invalid"
            if (
                not all(math.isfinite(value) for value in (candidate, low, high))
                or candidate < 0.0
                or low < 0.0
                or low > candidate
                or candidate > high
                or cell["method_A_status"] != "available"
                or cell["support_class"] not in {"supported", "marginal"}
            ):
                return None, "method_a_comparison_candidate_invalid"
        else:
            if candidate is not None or low is not None or high is not None:
                return None, "method_a_comparison_candidate_invalid"
        copied[t_index].append({
            "t_index": int(t_index),
            "t_low": float(cell["t_low"]),
            "t_high": float(cell["t_high"]),
            "delta_index": int(delta_index),
            "delta_low": float(cell["delta_low"]),
            "delta_high": float(cell["delta_high"]),
            "method_a_comparison_candidate": candidate,
            "method_a_comparison_candidate_low": low,
            "method_a_comparison_candidate_high": high,
            "method_a_comparison_candidate_status": str(status),
            "support_class": cell["support_class"],
            "method_A_status": cell["method_A_status"],
            "method_A_reason": cell["method_A_reason"],
            "method_a_comparison_candidate_reason": cell[
                "method_a_comparison_candidate_reason"
            ],
        })
    if len(seen) != expected_count:
        return None, "method_a_comparison_cell_grid_invalid"
    return tuple(
        tuple(sorted(rows, key=lambda cell: cell["delta_index"])) for rows in copied
    ), None


def build_full_background_subtraction_d9_payload(
    pion_hgcer_tdelta_diagnostic,
    pion_hgcer_method_a,
    pion_hgcer_method_a_comparison,
):
    """Select frozen Part-1/D.2 presentation inputs without recalculation."""
    response = _mapping(pion_hgcer_tdelta_diagnostic)
    if (
        response.get("status") not in {"available", "unavailable"}
        or response.get("non_authoritative") is not True
        or response.get("production_side_effect_free") is not True
        or response.get("production_hgcer_pid_unchanged") is not True
        or response.get("rf_restoration_applied") is not False
    ):
        return _d9_unavailable("part1_diagnostic_unavailable")
    t_edges = _strict_edges(response.get("t_edges"))
    delta_edges = _strict_edges(response.get("delta_edges"))
    coordinate_fingerprint = response.get("coordinate_fingerprint")
    if t_edges is None or delta_edges is None:
        return _d9_unavailable("part1_geometry_invalid")
    if not isinstance(coordinate_fingerprint, str) or not coordinate_fingerprint:
        return _d9_unavailable("part1_coordinate_fingerprint_missing")
    part1_available = response.get("status") == "available"
    records = _mapping(response.get("records")) if part1_available else {}
    histograms = _mapping(response.get("histograms")) if part1_available else {}
    row_groups = {}
    row_failures = {}
    for side in _D9_SIDES:
        if part1_available:
            row_groups[side], row_failures[side] = _d9_part1_rows(
                records.get(side), side, len(t_edges) - 1, len(delta_edges) - 1
            )
        else:
            row_groups[side] = tuple(
                tuple() for _unused in range(len(t_edges) - 1)
            )
            row_failures[side] = "part1_diagnostic_unavailable"
    template_keys = {
        "response": {
            "kaon": "H_hgcer_kaon_weighted",
            "pion": "H_hgcer_pion_weighted",
        },
        "delta": {
            "kaon": "H_hgcer_vs_delta_kaon_weighted",
            "pion": "H_hgcer_vs_delta_pion_weighted",
        },
    }
    thresholds, threshold_reason = _d9_thresholds(pion_hgcer_method_a)
    relative_groups, relative_reason = _d9_method_a_relative_cells(
        pion_hgcer_method_a,
        pion_hgcer_method_a_comparison,
        t_edges,
        delta_edges,
        coordinate_fingerprint,
    )
    per_t = []
    for t_index in range(len(t_edges) - 1):
        response_templates = {
            side: histograms.get(template_keys["response"][side]) for side in _D9_SIDES
        }
        delta_templates = {
            side: histograms.get(template_keys["delta"][side]) for side in _D9_SIDES
        }
        response_reason = next(
            (row_failures[side] for side in _D9_SIDES if row_failures[side] is not None),
            None,
        )
        if response_reason is None and any(
            response_templates[side] is None for side in _D9_SIDES
        ):
            response_reason = "part1_weighted_response_template_missing"
        delta_reason = next(
            (row_failures[side] for side in _D9_SIDES if row_failures[side] is not None),
            None,
        )
        if delta_reason is None and any(
            delta_templates[side] is None for side in _D9_SIDES
        ):
            delta_reason = "part1_weighted_delta_template_missing"
        per_t.append({
            "t_index": int(t_index),
            "t_low": float(t_edges[t_index]),
            "t_high": float(t_edges[t_index + 1]),
            "hgcer_response": {
                "available": response_reason is None,
                "reason": response_reason,
                "kaon_template": response_templates["kaon"],
                "pion_template": response_templates["pion"],
                "kaon_rows": tuple(row_groups.get("kaon") or ())[t_index]
                if row_groups.get("kaon") is not None else tuple(),
                "pion_rows": tuple(row_groups.get("pion") or ())[t_index]
                if row_groups.get("pion") is not None else tuple(),
            },
            "hgcer_delta_response": {
                "available": delta_reason is None,
                "reason": delta_reason,
                "kaon_template": delta_templates["kaon"],
                "pion_template": delta_templates["pion"],
                "kaon_rows": tuple(row_groups.get("kaon") or ())[t_index]
                if row_groups.get("kaon") is not None else tuple(),
                "pion_rows": tuple(row_groups.get("pion") or ())[t_index]
                if row_groups.get("pion") is not None else tuple(),
            },
            "method_a_relative": {
                "available": relative_reason is None,
                "reason": relative_reason,
                "cells": relative_groups[t_index]
                if relative_groups is not None else tuple(),
            },
        })
    return {
        "schema_version": D9_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "coordinate_fingerprint": str(coordinate_fingerprint),
        "method_a_thresholds": dict(thresholds),
        "method_a_threshold_reason": threshold_reason,
        "per_t": per_t,
    }


def full_background_subtraction_pdf_path(main_pdf):
    """Return the detached cumulative procedure-PDF destination."""
    root, extension = os.path.splitext(os.fspath(main_pdf))
    return "{}{}{}".format(
        root,
        FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX,
        extension or ".pdf",
    )


def _clone_display_histogram(histogram, name):
    if histogram is None or not hasattr(histogram, "Clone"):
        return None
    try:
        clone = histogram.Clone(str(name))
        if hasattr(clone, "SetDirectory"):
            clone.SetDirectory(0)
        return clone
    except Exception:
        return None


def _set_histogram_title(histogram, title):
    if histogram is not None and hasattr(histogram, "SetTitle"):
        histogram.SetTitle(str(title))


def _style_histogram(histogram, color=None):
    if histogram is None:
        return
    if color is not None and hasattr(histogram, "SetLineColor"):
        histogram.SetLineColor(color)
    if hasattr(histogram, "SetLineWidth"):
        histogram.SetLineWidth(2)
    if hasattr(histogram, "SetStats"):
        histogram.SetStats(0)


def _t_context(group):
    return "|t| = [{:.4f}, {:.4f}] GeV^2".format(
        float(group["t_low"]), float(group["t_high"])
    )


def _draw_small_note(ROOT, text):
    note = ROOT.TPaveText(0.14, 0.84, 0.86, 0.92, "NDC")
    note.SetFillStyle(0)
    note.SetBorderSize(0)
    note.SetTextAlign(22)
    note.SetTextSize(0.030)
    note.AddText(str(text))
    note.Draw()
    return note


def _draw_page_header(ROOT, canvas, heading, group):
    """Draw a retained, page-visible canonical-t header on a multipanel canvas."""
    canvas.cd()
    header = ROOT.TPaveText(0.16, 0.925, 0.84, 0.995, "NDC")
    header.SetFillStyle(0)
    header.SetBorderSize(0)
    header.SetTextAlign(22)
    header.SetTextSize(0.028)
    header.AddText(str(heading))
    header.AddText(_t_context(group))
    header.Draw()
    return header


def _histogram_visible_extrema(histogram):
    """Return finite ordinary-bin extrema without changing a histogram."""
    if histogram is None:
        return None
    try:
        bin_count = int(histogram.GetNbinsX())
    except Exception:
        return None
    minimum = maximum = None
    for bin_index in range(1, bin_count + 1):
        try:
            value = float(histogram.GetBinContent(bin_index))
        except Exception:
            continue
        if not math.isfinite(value):
            continue
        minimum = value if minimum is None else min(minimum, value)
        maximum = value if maximum is None else max(maximum, value)
    if minimum is None:
        return None
    return minimum, maximum


def _combined_histogram_y_range(histograms):
    """Return a padded signed-y range spanning visible bins and zero."""
    minimum = maximum = 0.0
    found_content = False
    for histogram in tuple(histograms or ()):
        extrema = _histogram_visible_extrema(histogram)
        if extrema is None:
            continue
        found_content = True
        minimum = min(minimum, extrema[0])
        maximum = max(maximum, extrema[1])
    span = maximum - minimum
    if not found_content or span <= 0.0:
        return -0.05, 0.05
    padding = 0.05 * span
    return minimum - padding, maximum + padding


def _apply_display_y_range(histogram, y_range):
    """Apply display-only limits to an axis-owning detached histogram clone."""
    if histogram is None or not isinstance(y_range, tuple) or len(y_range) != 2:
        return
    try:
        histogram.SetMinimum(float(y_range[0]))
        histogram.SetMaximum(float(y_range[1]))
    except Exception:
        return


def _render_raw_mm_page(ROOT, pdf_name, group):
    source = _mapping(group.get("raw_mm")).get("histogram")
    histogram = _clone_display_histogram(
        source, "H_full_background_d6_raw_t{}".format(group["t_index"] + 1)
    )
    if histogram is None:
        return False
    title = "Kaon-selected missing mass - {}".format(_t_context(group))
    _set_histogram_title(
        histogram, "{};Missing mass [GeV];Signed normalized yield".format(title)
    )
    _style_histogram(histogram, getattr(ROOT, "kBlack", 1))
    canvas = ROOT.TCanvas(
        "C_full_background_d6_raw_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    try:
        histogram.Draw("hist e")
        legend = ROOT.TLegend(0.66, 0.76, 0.88, 0.87)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(histogram, "Kaon-selected data", "l")
        legend.Draw()
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_timing_t_pid_page(ROOT, pdf_name, group, delta_edges):
    pid = _mapping(group.get("pid"))
    histograms = tuple(pid.get("cell_histograms") or ())
    if len(histograms) != len(delta_edges) - 1:
        return False
    panel_count = len(histograms)
    columns = min(3, max(1, panel_count))
    rows = int(math.ceil(float(panel_count) / float(columns)))
    title = "Proton-identification timing - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d6_pid_t{}".format(group["t_index"] + 1), title, 1300, 850
    )
    canvas.Divide(columns, rows)
    draw_objects = []
    try:
        for delta_index, source in enumerate(histograms):
            canvas.cd(delta_index + 1)
            histogram = _clone_display_histogram(
                source,
                "H_full_background_d6_pid_t{}_d{}".format(
                    group["t_index"] + 1, delta_index + 1
                ),
            )
            if histogram is None:
                return False
            panel_title = "delta = [{:.3f}, {:.3f}]%".format(
                delta_edges[delta_index], delta_edges[delta_index + 1]
            )
            _set_histogram_title(
                histogram,
                "{};{};Signed weighted yield".format(
                    panel_title, pid["timing_axis_label"]
                ),
            )
            _style_histogram(histogram, getattr(ROOT, "kBlue", 4))
            histogram.Draw("hist e")
            draw_objects.append(histogram)
        canvas.cd(1)
        draw_objects.append(_draw_small_note(ROOT, "Proton-identification timing"))
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Proton-identification timing", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_ctime_aero_pid_page(ROOT, pdf_name, group):
    source = _mapping(group.get("pid")).get("histogram")
    histogram = _clone_display_histogram(
        source, "H_full_background_d6_pid_ctime_t{}".format(group["t_index"] + 1)
    )
    if histogram is None:
        return False
    title = "Proton-identification timing vs aerogel - {}".format(_t_context(group))
    _set_histogram_title(
        histogram,
        "{};Aerogel NPE;Coincidence time [ns];Signed weighted yield".format(title),
    )
    if hasattr(histogram, "SetStats"):
        histogram.SetStats(0)
    canvas = ROOT.TCanvas(
        "C_full_background_d6_pid_ctime_t{}".format(group["t_index"] + 1),
        title,
        1200,
        800,
    )
    try:
        histogram.Draw("colz")
        _draw_small_note(ROOT, "Same timing/aerogel weighting geometry applies to this |t| bin")
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_timing_t_weight_page(ROOT, pdf_name, group):
    source = _mapping(group.get("weight")).get("histogram")
    display_map = _clone_display_histogram(
        source, "H_full_background_d6_weight_map_t{}".format(group["t_index"] + 1)
    )
    if display_map is None or not hasattr(display_map, "ProjectionX"):
        return False
    try:
        histogram = display_map.ProjectionX(
            "H_full_background_d6_weight_t{}".format(group["t_index"] + 1),
            int(group["t_index"]) + 1,
            int(group["t_index"]) + 1,
            "e",
        )
    except Exception:
        return False
    if hasattr(histogram, "SetDirectory"):
        histogram.SetDirectory(0)
    title = "Proton contamination weight - {}".format(_t_context(group))
    _set_histogram_title(histogram, "{};delta [%];Proton contamination weight".format(title))
    _style_histogram(histogram, getattr(ROOT, "kViolet", 6))
    canvas = ROOT.TCanvas(
        "C_full_background_d6_weight_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    try:
        histogram.Draw("hist e")
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_ctime_aero_weight_page(ROOT, pdf_name, group):
    source = _mapping(group.get("weight")).get("histogram")
    histogram = _clone_display_histogram(
        source, "H_full_background_d6_weight_ctime_t{}".format(group["t_index"] + 1)
    )
    if histogram is None:
        return False
    title = "Proton contamination weight - {}".format(_t_context(group))
    _set_histogram_title(
        histogram,
        "{};delta [%];Aerogel NPE;Proton contamination weight".format(title),
    )
    if hasattr(histogram, "SetStats"):
        histogram.SetStats(0)
    canvas = ROOT.TCanvas(
        "C_full_background_d6_weight_ctime_t{}".format(group["t_index"] + 1),
        title,
        1200,
        800,
    )
    try:
        histogram.Draw("colz")
        _draw_small_note(ROOT, "Same timing/aerogel weighting geometry applies to this |t| bin")
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def open_full_background_subtraction_pdf(pdf_name):
    """Open a multipage procedure PDF without emitting a page."""
    ROOT = _import_root()
    if ROOT is None:
        return False
    canvas = ROOT.TCanvas("C_full_background_subtraction_open", "procedure PDF", 1, 1)
    try:
        canvas.Print("{}[".format(pdf_name))
    finally:
        canvas.Close()
    return True


def close_full_background_subtraction_pdf(pdf_name):
    """Close a multipage procedure PDF without emitting a terminal page."""
    ROOT = _import_root()
    if ROOT is None:
        return False
    canvas = ROOT.TCanvas("C_full_background_subtraction_close", "procedure PDF", 1, 1)
    try:
        canvas.Print("{}]".format(pdf_name))
    finally:
        canvas.Close()
    return True


def _render_d6_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render the existing D.6 page group for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    if _mapping(group.get("raw_mm")).get("available"):
        if _render_raw_mm_page(ROOT, pdf_name, group):
            manifest.append({"page_id": "full_background.d6.raw_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.6 raw MM page unavailable for t{}".format(t_number))
    else:
        failures.append("D.6 raw MM input unavailable for t{}".format(t_number))

    method = presentation.get("method")
    if _mapping(group.get("pid")).get("available"):
        rendered_pid = (
            _render_timing_t_pid_page(ROOT, pdf_name, group, presentation.get("delta_edges") or ())
            if method == _TIMING_T_METHOD
            else _render_ctime_aero_pid_page(ROOT, pdf_name, group)
            if method == _CTIME_AERO_METHOD
            else False
        )
        if rendered_pid:
            manifest.append({"page_id": "full_background.d6.proton_pid", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.6 proton-identification page unavailable for t{}".format(t_number))
    else:
        failures.append("D.6 proton-identification input unavailable for t{}".format(t_number))

    if _mapping(group.get("weight")).get("available"):
        rendered_weight = (
            _render_timing_t_weight_page(ROOT, pdf_name, group)
            if method == _TIMING_T_METHOD
            else _render_ctime_aero_weight_page(ROOT, pdf_name, group)
            if method == _CTIME_AERO_METHOD
            else False
        )
        if rendered_weight:
            manifest.append({"page_id": "full_background.d6.proton_weight", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.6 proton-weight page unavailable for t{}".format(t_number))
    else:
        failures.append("D.6 proton-weight input unavailable for t{}".format(t_number))


def _render_d7_mm_overlay_page(
    ROOT, pdf_name, group, other_key, page_title, other_label, other_color, page_id
):
    raw = _mapping(group.get("raw_mm"))
    other = _mapping(group.get(other_key))
    raw_histogram = _clone_display_histogram(
        raw.get("histogram"), "H_full_background_d7_raw_{}_t{}".format(page_id.rsplit(".", 1)[-1], group["t_index"] + 1)
    )
    other_histogram = _clone_display_histogram(
        other.get("histogram"), "H_full_background_d7_other_{}_t{}".format(page_id.rsplit(".", 1)[-1], group["t_index"] + 1)
    )
    if raw_histogram is None or other_histogram is None:
        return False
    title = "{} - {}".format(page_title, _t_context(group))
    _set_histogram_title(
        raw_histogram, "{};Missing mass [GeV];Signed normalized yield".format(title)
    )
    _apply_display_y_range(
        raw_histogram,
        _combined_histogram_y_range((raw_histogram, other_histogram)),
    )
    _style_histogram(raw_histogram, getattr(ROOT, "kBlack", 1))
    _style_histogram(other_histogram, getattr(ROOT, other_color, 2))
    canvas = ROOT.TCanvas(
        "C_full_background_d7_{}_t{}".format(page_id.rsplit(".", 1)[-1], group["t_index"] + 1),
        title,
        1200,
        800,
    )
    try:
        raw_histogram.Draw("hist e")
        other_histogram.Draw("hist e same")
        legend = ROOT.TLegend(0.64, 0.72, 0.89, 0.87)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(raw_histogram, "Kaon-selected data", "l")
        legend.AddEntry(other_histogram, other_label, "l")
        legend.Draw()
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _new_d7_delta_histogram(source, name):
    """Create one detached empty local-MM display histogram from its source."""
    histogram = _clone_display_histogram(source, name)
    if histogram is None or not hasattr(histogram, "Reset"):
        return None
    try:
        histogram.Reset("ICES")
    except Exception:
        return None
    return histogram


def _render_d7_delta_page(ROOT, pdf_name, group, delta_edges):
    projection = _mapping(group.get("delta_projection"))
    rows_by_delta = tuple(projection.get("rows_by_delta") or ())
    if len(rows_by_delta) != len(delta_edges) - 1:
        return False
    panel_count = len(rows_by_delta)
    columns = min(3, max(1, panel_count))
    rows = int(math.ceil(float(panel_count) / float(columns)))
    title = "Proton subtraction across delta - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d7_delta_t{}".format(group["t_index"] + 1),
        title,
        1300,
        850,
    )
    canvas.Divide(columns, rows)
    draw_objects = []
    source = _mapping(group.get("raw_mm")).get("histogram")
    try:
        panel_histograms = []
        for delta_index, display_rows in enumerate(rows_by_delta):
            raw_histogram = _new_d7_delta_histogram(
                source,
                "H_full_background_d7_delta_raw_t{}_d{}".format(
                    group["t_index"] + 1, delta_index + 1
                ),
            )
            proton_histogram = _new_d7_delta_histogram(
                source,
                "H_full_background_d7_delta_proton_t{}_d{}".format(
                    group["t_index"] + 1, delta_index + 1
                ),
            )
            cleaned_histogram = _new_d7_delta_histogram(
                source,
                "H_full_background_d7_delta_cleaned_t{}_d{}".format(
                    group["t_index"] + 1, delta_index + 1
                ),
            )
            if raw_histogram is None or proton_histogram is None or cleaned_histogram is None:
                return False
            for row in display_rows:
                row = _mapping(row)
                raw_histogram.Fill(row["missing_mass"], row["raw_weight"])
                proton_histogram.Fill(row["missing_mass"], row["proton_contribution"])
                cleaned_histogram.Fill(row["missing_mass"], row["cleaned_contribution"])
            panel_histograms.append((raw_histogram, proton_histogram, cleaned_histogram))

        common_y_range = _combined_histogram_y_range(
            histogram
            for panel in panel_histograms
            for histogram in panel
        )
        for delta_index, histograms in enumerate(panel_histograms):
            canvas.cd(delta_index + 1)
            raw_histogram, proton_histogram, cleaned_histogram = histograms
            panel_title = "delta = [{:.3f}, {:.3f}] %".format(
                delta_edges[delta_index], delta_edges[delta_index + 1]
            )
            _set_histogram_title(
                raw_histogram,
                "{};Missing mass [GeV];Signed normalized yield".format(panel_title),
            )
            _apply_display_y_range(raw_histogram, common_y_range)
            _style_histogram(raw_histogram, getattr(ROOT, "kBlack", 1))
            _style_histogram(proton_histogram, getattr(ROOT, "kRed", 2))
            _style_histogram(cleaned_histogram, getattr(ROOT, "kBlue", 4))
            raw_histogram.Draw("hist e")
            proton_histogram.Draw("hist e same")
            cleaned_histogram.Draw("hist e same")
            if delta_index == 0:
                legend = ROOT.TLegend(0.53, 0.66, 0.89, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(raw_histogram, "Kaon-selected data", "l")
                legend.AddEntry(proton_histogram, "Proton contamination", "l")
                legend.AddEntry(cleaned_histogram, "Proton-cleaned kaon data", "l")
                legend.Draw()
                draw_objects.append(legend)
            draw_objects.extend((raw_histogram, proton_histogram, cleaned_histogram))
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Proton subtraction across delta", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d7_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render the D.7 result group for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    raw = _mapping(group.get("raw_mm"))
    proton = _mapping(group.get("proton_mm"))
    cleaned = _mapping(group.get("cleaned_mm"))
    if raw.get("available") and proton.get("available"):
        if _render_d7_mm_overlay_page(
            ROOT,
            pdf_name,
            group,
            "proton_mm",
            "Proton contamination in missing mass",
            "Proton contamination",
            "kRed",
            "full_background.d7.proton_mm",
        ):
            manifest.append({"page_id": "full_background.d7.proton_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.7 proton-contamination MM page unavailable for t{}".format(t_number))
    else:
        failures.append("D.7 proton-contamination MM input unavailable for t{}".format(t_number))
    if raw.get("available") and cleaned.get("available"):
        if _render_d7_mm_overlay_page(
            ROOT,
            pdf_name,
            group,
            "cleaned_mm",
            "Before and after proton subtraction",
            "Proton-cleaned kaon data",
            "kBlue",
            "full_background.d7.proton_cleaned_mm",
        ):
            manifest.append({"page_id": "full_background.d7.proton_cleaned_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.7 proton-cleaned MM page unavailable for t{}".format(t_number))
    else:
        failures.append("D.7 proton-cleaned MM input unavailable for t{}".format(t_number))
    projection = _mapping(group.get("delta_projection"))
    if projection.get("available"):
        if _render_d7_delta_page(ROOT, pdf_name, group, presentation.get("delta_edges") or ()):
            manifest.append({"page_id": "full_background.d7.proton_delta_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.7 proton-subtraction delta page unavailable for t{}".format(t_number))
    else:
        failures.append("D.7 proton-subtraction delta input unavailable for t{}".format(t_number))


def _render_d8_mm_overlay_page(
    ROOT, pdf_name, group, other_key, page_title, other_label, other_color, page_id
):
    before = _mapping(group.get("before_pion_mm"))
    other = _mapping(group.get(other_key))
    before_histogram = _clone_display_histogram(
        before.get("histogram"),
        "H_full_background_d8_before_{}_t{}".format(
            page_id.rsplit(".", 1)[-1], group["t_index"] + 1
        ),
    )
    other_histogram = _clone_display_histogram(
        other.get("histogram"),
        "H_full_background_d8_other_{}_t{}".format(
            page_id.rsplit(".", 1)[-1], group["t_index"] + 1
        ),
    )
    if before_histogram is None or other_histogram is None:
        return False
    title = "{} - {}".format(page_title, _t_context(group))
    _set_histogram_title(
        before_histogram, "{};Missing mass [GeV];Signed normalized yield".format(title)
    )
    _apply_display_y_range(
        before_histogram,
        _combined_histogram_y_range((before_histogram, other_histogram)),
    )
    _style_histogram(before_histogram, getattr(ROOT, "kBlack", 1))
    _style_histogram(other_histogram, getattr(ROOT, other_color, 2))
    canvas = ROOT.TCanvas(
        "C_full_background_d8_{}_t{}".format(
            page_id.rsplit(".", 1)[-1], group["t_index"] + 1
        ),
        title,
        1200,
        800,
    )
    try:
        before_histogram.Draw("hist e")
        other_histogram.Draw("hist e same")
        legend = ROOT.TLegend(0.60, 0.72, 0.89, 0.87)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(before_histogram, "Proton-cleaned kaon data", "l")
        legend.AddEntry(other_histogram, other_label, "l")
        legend.Draw()
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d8_delta_page(ROOT, pdf_name, group, delta_edges):
    projection = _mapping(group.get("delta_projection"))
    rows_by_delta = tuple(projection.get("rows_by_delta") or ())
    if len(rows_by_delta) != len(delta_edges) - 1:
        return False
    panel_count = len(rows_by_delta)
    columns = min(3, max(1, panel_count))
    rows = int(math.ceil(float(panel_count) / float(columns)))
    title = "Baseline pion background across delta - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d8_delta_t{}".format(group["t_index"] + 1),
        title,
        1300,
        850,
    )
    canvas.Divide(columns, rows)
    draw_objects = []
    source = _mapping(group.get("baseline_pion_mm")).get("histogram")
    try:
        panel_histograms = []
        for delta_index, display_rows in enumerate(rows_by_delta):
            histogram = _new_d7_delta_histogram(
                source,
                "H_full_background_d8_delta_baseline_t{}_d{}".format(
                    group["t_index"] + 1, delta_index + 1
                ),
            )
            if histogram is None:
                return False
            for row in display_rows:
                row = _mapping(row)
                histogram.Fill(row["missing_mass"], row["baseline_contribution"])
            panel_histograms.append(histogram)

        common_y_range = _combined_histogram_y_range(panel_histograms)
        for delta_index, histogram in enumerate(panel_histograms):
            canvas.cd(delta_index + 1)
            panel_title = "delta = [{:.3f}, {:.3f}] %".format(
                delta_edges[delta_index], delta_edges[delta_index + 1]
            )
            _set_histogram_title(
                histogram,
                "{};Missing mass [GeV];Signed normalized yield".format(panel_title),
            )
            _apply_display_y_range(histogram, common_y_range)
            _style_histogram(histogram, getattr(ROOT, "kOrange", 800) + 7)
            histogram.Draw("hist e")
            if delta_index == 0:
                legend = ROOT.TLegend(0.57, 0.76, 0.89, 0.87)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(histogram, "Baseline pion background", "l")
                legend.Draw()
                draw_objects.append(legend)
            draw_objects.append(histogram)
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Baseline pion background across delta", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _append_d8_projection_failure(group, failures):
    projection = _mapping(group.get("delta_projection"))
    exclusions = _mapping(projection.get("exclusions"))
    t_number = int(group.get("t_index", -1)) + 1
    invalid_delta = int(exclusions.get("invalid_frozen_delta_index", 0) or 0)
    if invalid_delta:
        failures.append(
            "D.8 frozen delta rows excluded for t{}: invalid_frozen_delta_index={}".format(
                t_number, invalid_delta
            )
        )
    closure = _mapping(projection.get("closure"))
    if closure.get("status") == "incomplete_frozen_delta_coverage":
        failures.append("D.8 delta projection has incomplete frozen coverage for t{}".format(t_number))


def _render_d8_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render the D.8 baseline-pion group for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    before = _mapping(group.get("before_pion_mm"))
    baseline = _mapping(group.get("baseline_pion_mm"))
    after = _mapping(group.get("after_pion_mm"))
    if before.get("available") and baseline.get("available"):
        if _render_d8_mm_overlay_page(
            ROOT,
            pdf_name,
            group,
            "baseline_pion_mm",
            "Baseline pion background in missing mass",
            "Baseline pion background",
            "kOrange",
            "full_background.d8.pion_background_mm",
        ):
            manifest.append({"page_id": "full_background.d8.pion_background_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.8 baseline pion-background MM page unavailable for t{}".format(t_number))
    else:
        failures.append("D.8 baseline pion-background MM input unavailable for t{}".format(t_number))
    if before.get("available") and after.get("available"):
        if _render_d8_mm_overlay_page(
            ROOT,
            pdf_name,
            group,
            "after_pion_mm",
            "Before and after baseline pion subtraction",
            "Baseline pion-subtracted kaon data",
            "kGreen",
            "full_background.d8.pion_subtracted_mm",
        ):
            manifest.append({"page_id": "full_background.d8.pion_subtracted_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.8 pion-subtracted MM page unavailable for t{}".format(t_number))
    else:
        failures.append("D.8 pion-subtracted MM input unavailable for t{}".format(t_number))
    projection = _mapping(group.get("delta_projection"))
    _append_d8_projection_failure(group, failures)
    if projection.get("available"):
        if _render_d8_delta_page(ROOT, pdf_name, group, presentation.get("delta_edges") or ()):
            manifest.append({"page_id": "full_background.d8.pion_delta_mm", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.8 baseline pion delta page unavailable for t{}".format(t_number))
    else:
        failures.append(
            "D.8 baseline pion delta input unavailable for t{}: {}".format(
                t_number, projection.get("reason")
            )
        )


def _d9_fresh_weighted_histogram(source, name, rows, *, versus_delta):
    """Fill a detached Part-1 weighted display clone from frozen scalar rows."""
    histogram = _clone_display_histogram(source, name)
    if histogram is None or not hasattr(histogram, "Reset") or not hasattr(histogram, "Fill"):
        return None
    try:
        histogram.Reset("ICES")
        for row in rows:
            row = _mapping(row)
            if versus_delta:
                histogram.Fill(row["delta"], row["npe"], row["diagnostic_weight"])
            else:
                histogram.Fill(row["npe"], row["diagnostic_weight"])
    except Exception:
        return None
    return histogram


def _histogram_visible_extrema_2d(histogram):
    """Return finite ordinary-cell extrema without changing a 2D histogram."""
    if histogram is None:
        return None
    try:
        x_count = int(histogram.GetNbinsX())
        y_count = int(histogram.GetNbinsY())
    except Exception:
        return None
    minimum = maximum = None
    for x_index in range(1, x_count + 1):
        for y_index in range(1, y_count + 1):
            try:
                value = float(histogram.GetBinContent(x_index, y_index))
            except Exception:
                continue
            if not math.isfinite(value):
                continue
            minimum = value if minimum is None else min(minimum, value)
            maximum = value if maximum is None else max(maximum, value)
    if minimum is None:
        return None
    return minimum, maximum


def _signed_z_range(histograms):
    """Return one signed display range that always contains zero."""
    minimum = maximum = 0.0
    found_content = False
    for histogram in histograms:
        extrema = _histogram_visible_extrema_2d(histogram)
        if extrema is None:
            continue
        found_content = True
        minimum = min(minimum, extrema[0])
        maximum = max(maximum, extrema[1])
    span = maximum - minimum
    if not found_content or span <= 0.0:
        return -0.05, 0.05
    padding = 0.05 * span
    return minimum - padding, maximum + padding


def _apply_display_z_range(histogram, z_range):
    if histogram is None or not isinstance(z_range, tuple) or len(z_range) != 2:
        return
    try:
        histogram.SetMinimum(float(z_range[0]))
        histogram.SetMaximum(float(z_range[1]))
    except Exception:
        return


def _d9_draw_threshold(ROOT, x_value, y_range, x_low, x_high, draw_objects):
    """Draw the frozen Method-A low-response boundary only on display clones."""
    if not hasattr(ROOT, "TLine"):
        return
    try:
        line = ROOT.TLine(float(x_value), float(y_range[0]), float(x_value), float(y_range[1]))
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        draw_objects.append(line)
    except Exception:
        return


def _d9_draw_delta_threshold(ROOT, y_value, x_low, x_high, draw_objects):
    if not hasattr(ROOT, "TLine"):
        return
    try:
        line = ROOT.TLine(float(x_low), float(y_value), float(x_high), float(y_value))
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        draw_objects.append(line)
    except Exception:
        return


def _d9_method_a_threshold_text(thresholds):
    """Return the display text for a valid frozen Method-A threshold mapping."""
    thresholds = _mapping(thresholds)
    if thresholds.get("available") is not True:
        return None
    try:
        positive = float(thresholds["positive_response_threshold"])
        low_upper = float(thresholds["low_response_upper_threshold"])
    except (KeyError, TypeError, ValueError):
        return None
    if not math.isfinite(positive) or not math.isfinite(low_upper):
        return None
    return "Method A low response: {:g} < HGCer NPE <= {:g}".format(
        positive, low_upper
    )


def _d9_draw_notice(ROOT, text, y_low, y_high):
    """Draw one retained D.9 notice at an explicitly separate NDC position."""
    notice = ROOT.TPaveText(0.14, float(y_low), 0.86, float(y_high), "NDC")
    notice.SetFillStyle(0)
    notice.SetBorderSize(0)
    notice.SetTextAlign(22)
    notice.SetTextSize(0.030)
    notice.AddText(str(text))
    notice.Draw()
    return notice


def _render_d9_hgcer_response_page(ROOT, pdf_name, group, thresholds):
    presentation = _mapping(group.get("hgcer_response"))
    kaon = _d9_fresh_weighted_histogram(
        presentation.get("kaon_template"),
        "H_full_background_d9_hgcer_kaon_t{}".format(group["t_index"] + 1),
        presentation.get("kaon_rows") or (),
        versus_delta=False,
    )
    pion = _d9_fresh_weighted_histogram(
        presentation.get("pion_template"),
        "H_full_background_d9_hgcer_pion_t{}".format(group["t_index"] + 1),
        presentation.get("pion_rows") or (),
        versus_delta=False,
    )
    if kaon is None or pion is None:
        return False
    title = "HGCer response - {}".format(_t_context(group))
    _set_histogram_title(kaon, "{};HGCer NPE;Signed weighted yield".format(title))
    y_range = _combined_histogram_y_range((kaon, pion))
    _apply_display_y_range(kaon, y_range)
    _style_histogram(kaon, getattr(ROOT, "kBlack", 1))
    _style_histogram(pion, getattr(ROOT, "kBlue", 4))
    canvas = ROOT.TCanvas(
        "C_full_background_d9_hgcer_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    draw_objects = []
    try:
        kaon.Draw("hist e")
        pion.Draw("hist e same")
        legend = ROOT.TLegend(0.62, 0.58, 0.89, 0.72)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(kaon, "Proton-cleaned kaon data", "l")
        legend.AddEntry(pion, "Pion-control sample", "l")
        legend.Draw()
        draw_objects.extend((kaon, pion, legend))
        threshold_text = _d9_method_a_threshold_text(thresholds)
        if threshold_text is not None:
            _d9_draw_threshold(
                ROOT,
                thresholds["low_response_upper_threshold"],
                y_range,
                0.0,
                thresholds["low_response_upper_threshold"],
                draw_objects,
            )
            draw_objects.append(
                _d9_draw_notice(ROOT, threshold_text, 0.75, 0.83)
            )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d9_hgcer_delta_page(ROOT, pdf_name, group, thresholds, delta_edges):
    presentation = _mapping(group.get("hgcer_delta_response"))
    kaon = _d9_fresh_weighted_histogram(
        presentation.get("kaon_template"),
        "H_full_background_d9_hgcer_delta_kaon_t{}".format(group["t_index"] + 1),
        presentation.get("kaon_rows") or (),
        versus_delta=True,
    )
    pion = _d9_fresh_weighted_histogram(
        presentation.get("pion_template"),
        "H_full_background_d9_hgcer_delta_pion_t{}".format(group["t_index"] + 1),
        presentation.get("pion_rows") or (),
        versus_delta=True,
    )
    if kaon is None or pion is None:
        return False
    z_range = _signed_z_range((kaon, pion))
    _apply_display_z_range(kaon, z_range)
    _apply_display_z_range(pion, z_range)
    title = "HGCer response across delta - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d9_hgcer_delta_t{}".format(group["t_index"] + 1),
        title,
        1300,
        800,
    )
    canvas.Divide(2, 1)
    draw_objects = []
    threshold_text = _d9_method_a_threshold_text(thresholds)
    try:
        for pad_index, (histogram, label) in enumerate(
            ((kaon, "Proton-cleaned kaon data"), (pion, "Pion-control sample")), 1
        ):
            canvas.cd(pad_index)
            _set_histogram_title(
                histogram, "{};delta [%];HGCer NPE;Signed weighted yield".format(label)
            )
            if hasattr(histogram, "SetStats"):
                histogram.SetStats(0)
            histogram.Draw("colz")
            if threshold_text is not None:
                _d9_draw_delta_threshold(
                    ROOT,
                    thresholds["low_response_upper_threshold"],
                    delta_edges[0],
                    delta_edges[-1],
                    draw_objects,
                )
            draw_objects.append(histogram)
        canvas.cd(1)
        if threshold_text is not None:
            draw_objects.append(_d9_draw_notice(ROOT, threshold_text, 0.75, 0.83))
        draw_objects.append(_draw_page_header(ROOT, canvas, "HGCer response across delta", group))
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _d9_relative_y_range(cells):
    values = [1.0]
    for cell in cells:
        cell = _mapping(cell)
        if cell.get("method_a_comparison_candidate_status") not in {"available", "marginal"}:
            continue
        for key in (
            "method_a_comparison_candidate_low",
            "method_a_comparison_candidate",
            "method_a_comparison_candidate_high",
        ):
            try:
                value = float(cell.get(key))
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                values.append(value)
    low, high = min(values), max(values)
    if high <= low:
        return low - 0.05, high + 0.05
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def _render_d9_method_a_relative_page(ROOT, pdf_name, group, delta_edges):
    presentation = _mapping(group.get("method_a_relative"))
    cells = tuple(presentation.get("cells") or ())
    if not hasattr(ROOT, "TH1D"):
        return False
    title = "Method A - HGCer response diagnostic - {}".format(_t_context(group))
    try:
        frame = ROOT.TH1D(
            "H_full_background_d9_method_a_relative_t{}".format(group["t_index"] + 1),
            "{};delta [%];Relative pion-background diagnostic".format(title),
            len(delta_edges) - 1,
            array("d", delta_edges),
        )
        frame.SetDirectory(0)
        if hasattr(frame, "SetStats"):
            frame.SetStats(0)
        y_range = _d9_relative_y_range(cells)
        frame.SetMinimum(y_range[0])
        frame.SetMaximum(y_range[1])
    except Exception:
        return False
    canvas = ROOT.TCanvas(
        "C_full_background_d9_method_a_relative_t{}".format(group["t_index"] + 1),
        title,
        1200,
        800,
    )
    draw_objects = [frame]
    try:
        frame.Draw("axis")
        graphs = {}
        for status, color, label in (
            ("available", getattr(ROOT, "kBlack", 1), "Supported"),
            ("marginal", getattr(ROOT, "kBlue", 4), "Marginal"),
        ):
            selected = [
                _mapping(cell) for cell in cells
                if _mapping(cell).get("method_a_comparison_candidate_status") == status
            ]
            if not selected or not hasattr(ROOT, "TGraphAsymmErrors"):
                continue
            graph = ROOT.TGraphAsymmErrors(len(selected))
            for point_index, cell in enumerate(selected):
                x_value = 0.5 * (float(cell["delta_low"]) + float(cell["delta_high"]))
                value = float(cell["method_a_comparison_candidate"])
                low = float(cell["method_a_comparison_candidate_low"])
                high = float(cell["method_a_comparison_candidate_high"])
                graph.SetPoint(point_index, x_value, value)
                graph.SetPointError(point_index, 0.0, 0.0, value - low, high - value)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetMarkerStyle(20 if status == "available" else 24)
            graph.SetMarkerSize(1.1)
            graph.Draw("P same")
            graphs[status] = (graph, label)
            draw_objects.append(graph)
        if hasattr(ROOT, "TLine"):
            unity = ROOT.TLine(float(delta_edges[0]), 1.0, float(delta_edges[-1]), 1.0)
            unity.SetLineStyle(2)
            unity.SetLineWidth(2)
            unity.Draw()
            draw_objects.append(unity)
        if graphs:
            legend = ROOT.TLegend(0.64, 0.72, 0.89, 0.87)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            for graph, label in graphs.values():
                legend.AddEntry(graph, label, "p")
            legend.Draw()
            draw_objects.append(legend)
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Method A - HGCer response diagnostic", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d9_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render D.9 HGCer/Method-A pages for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    thresholds = _mapping(presentation.get("method_a_thresholds"))
    response = _mapping(group.get("hgcer_response"))
    if response.get("available"):
        if _render_d9_hgcer_response_page(ROOT, pdf_name, group, thresholds):
            manifest.append({"page_id": "full_background.d9.hgcer_response", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.9 HGCer response page unavailable for t{}".format(t_number))
    else:
        failures.append("D.9 HGCer response input unavailable for t{}: {}".format(t_number, response.get("reason")))
    delta_response = _mapping(group.get("hgcer_delta_response"))
    if delta_response.get("available"):
        if _render_d9_hgcer_delta_page(
            ROOT, pdf_name, group, thresholds, presentation.get("delta_edges") or ()
        ):
            manifest.append({"page_id": "full_background.d9.hgcer_delta_response", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.9 HGCer delta-response page unavailable for t{}".format(t_number))
    else:
        failures.append("D.9 HGCer delta-response input unavailable for t{}: {}".format(t_number, delta_response.get("reason")))
    relative = _mapping(group.get("method_a_relative"))
    if relative.get("available"):
        if _render_d9_method_a_relative_page(
            ROOT, pdf_name, group, presentation.get("delta_edges") or ()
        ):
            manifest.append({"page_id": "full_background.d9.method_a_relative", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.9 Method-A relative page unavailable for t{}".format(t_number))
    else:
        failures.append("D.9 Method-A relative input unavailable for t{}: {}".format(t_number, relative.get("reason")))


def _append_d7_exclusion_failure(presentation, failures):
    exclusions = _mapping(presentation.get("projection_exclusions"))
    relevant = (
        "invalid_source_entry_index",
        "missing_frozen_lookup_row",
        "missing_mm_or_weight",
        "missing_applied_factor",
        "invalid_frozen_delta_index",
        "invalid_canonical_t_index",
    )
    parts = []
    for key in relevant:
        try:
            count = int(exclusions.get(key, 0) or 0)
        except (TypeError, ValueError):
            count = 0
        if count > 0:
            parts.append("{}={}".format(key, count))
    if parts:
        failures.append("D.7 frozen display rows excluded: {}".format(", ".join(parts)))


def _d10_integer(value):
    if isinstance(value, bool):
        return None
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    integer = int(scalar)
    if not math.isfinite(scalar) or scalar != integer:
        return None
    return integer


def _d10_finite(value):
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    return scalar if math.isfinite(scalar) else None


def _d10_nonempty_string(value):
    return isinstance(value, str) and bool(value.strip())


def _d10_sequence(value):
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _d10_region_rows(value, *, protected):
    """Copy only frozen Method-B region display metadata."""
    if not _d10_sequence(value) or not value:
        return None, "method_b_region_contract_invalid"
    result = []
    names = set()
    expected_role = "protected_signal" if protected else "pion_sensitive"
    for source in value:
        row = _mapping(source)
        required = (
            "region_name", "mm_low", "mm_high", "region_role", "window_source",
            "mm_offset_applied",
        )
        if not row or any(key not in row for key in required):
            return None, "method_b_region_contract_invalid"
        name = row["region_name"]
        low = _d10_finite(row["mm_low"])
        high = _d10_finite(row["mm_high"])
        offset = _d10_finite(row["mm_offset_applied"])
        if (
            not _d10_nonempty_string(name)
            or name in names
            or low is None
            or high is None
            or low >= high
            or offset is None
            or row["region_role"] != expected_role
            or not _d10_nonempty_string(row["window_source"])
        ):
            return None, "method_b_region_contract_invalid"
        copied = {
            "region_name": str(name),
            "mm_low": low,
            "mm_high": high,
            "region_role": expected_role,
            "window_source": str(row["window_source"]),
            "mm_offset_applied": offset,
        }
        if not protected:
            if (
                not isinstance(row.get("available"), bool)
                or not isinstance(row.get("protected_signal_overlap"), bool)
            ):
                return None, "method_b_region_contract_invalid"
            copied.update({
                "available": row["available"],
                "protected_signal_overlap": row["protected_signal_overlap"],
                "reason": row.get("reason"),
            })
        names.add(name)
        result.append(copied)
    if protected:
        if len(result) != 1:
            return None, "method_b_protected_region_contract_invalid"
        protected_row = result[0]
        if (
            protected_row["region_name"] != "KLambdaSigma0"
            or protected_row["mm_low"] != 1.10
            or protected_row["mm_high"] != 1.23
            or protected_row["mm_offset_applied"] != 0.0
        ):
            return None, "method_b_protected_region_contract_invalid"
    return tuple(result), None


def _d10_phase_a_contract(phase_a):
    phase = _mapping(phase_a)
    required = (
        "schema_version", "status", "available", "contract_fingerprint",
        "coordinate_fingerprint", "canonical_t_edges", "delta_edges", "host_state",
        "source_target_state", "pion_records", "kaon_host_records", "pion_closure",
        "host_closure", "production_objects_mutated", "refinement_applied",
        "rf_restoration_applied",
    )
    if not phase or any(key not in phase for key in required):
        return None, "phase_a_contract_invalid"
    t_edges = _strict_edges(phase["canonical_t_edges"])
    delta_edges = _strict_edges(phase["delta_edges"])
    pion_closure = _mapping(phase["pion_closure"])
    host_closure = _mapping(phase["host_closure"])
    if (
        phase["schema_version"] != _D10_PHASE_A_SCHEMA
        or phase["status"] != "available"
        or phase["available"] is not True
        or not _d10_nonempty_string(phase["contract_fingerprint"])
        or not _d10_nonempty_string(phase["coordinate_fingerprint"])
        or t_edges is None
        or delta_edges is None
        or phase["host_state"] not in {"proton_cleaned", "identity_no_proton_cleaning"}
        or phase["source_target_state"] != _D10_SOURCE_TARGET_STATE
        or pion_closure.get("passed") is not True
        or host_closure.get("passed") is not True
        or phase["production_objects_mutated"] is not False
        or phase["refinement_applied"] is not False
        or phase["rf_restoration_applied"] is not False
        or not _d10_sequence(phase["pion_records"])
        or not _d10_sequence(phase["kaon_host_records"])
    ):
        return None, "phase_a_contract_invalid"
    return {
        "contract_fingerprint": str(phase["contract_fingerprint"]),
        "coordinate_fingerprint": str(phase["coordinate_fingerprint"]),
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "host_state": phase["host_state"],
        "source_target_state": phase["source_target_state"],
        "pion_records": phase["pion_records"],
        "kaon_host_records": phase["kaon_host_records"],
    }, None


def _d10_method_b_contract(method_b, phase):
    method = _mapping(method_b)
    required = (
        "schema_version", "status", "available", "non_authoritative",
        "production_objects_mutated", "refinement_applied", "rf_ct_required",
        "interpolation_used", "phase_a_records_only", "method_a_numerical_dependency",
        "phase_a_contract_fingerprint", "coordinate_fingerprint", "host_state",
        "source_target_state", "t_edges", "delta_edges", "mm_binning", "mm_regions",
        "protected_regions", "fingerprint", "cells",
    )
    if not method or any(key not in method for key in required):
        return None, "method_b_contract_invalid"
    t_edges = _strict_edges(method["t_edges"])
    delta_edges = _strict_edges(method["delta_edges"])
    mm_edges = _strict_edges(method["mm_binning"])
    if (
        method["schema_version"] != _D10_METHOD_B_SCHEMA
        or method["status"] != "available"
        or method["available"] is not True
        or method["non_authoritative"] is not True
        or method["production_objects_mutated"] is not False
        or method["refinement_applied"] is not False
        or method["rf_ct_required"] is not False
        or method["interpolation_used"] is not False
        or method["phase_a_records_only"] is not True
        or method["method_a_numerical_dependency"] is not False
        or not _d10_nonempty_string(method["fingerprint"])
        or method["phase_a_contract_fingerprint"] != phase["contract_fingerprint"]
        or method["coordinate_fingerprint"] != phase["coordinate_fingerprint"]
        or list(t_edges or ()) != phase["t_edges"]
        or list(delta_edges or ()) != phase["delta_edges"]
        or method["host_state"] != phase["host_state"]
        or method["source_target_state"] != phase["source_target_state"]
        or method["source_target_state"] != _D10_SOURCE_TARGET_STATE
        or mm_edges is None
        or not _d10_sequence(method["cells"])
    ):
        return None, "method_b_contract_invalid"
    regions, reason = _d10_region_rows(method["mm_regions"], protected=False)
    if reason is not None:
        return None, reason
    protected, reason = _d10_region_rows(method["protected_regions"], protected=True)
    if reason is not None:
        return None, reason
    return {
        "fingerprint": str(method["fingerprint"]),
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "mm_edges": list(mm_edges),
        "mm_regions": regions,
        "protected_regions": protected,
        "cells": method["cells"],
    }, None


def _d10_event_rows(records, *, contribution_key, t_count, delta_count, host_state=None):
    copied = [[] for _unused in range(t_count)]
    failures = [None for _unused in range(t_count)]

    def fail(t_index, reason):
        if t_index is None:
            for index in range(t_count):
                failures[index] = failures[index] or reason
        else:
            failures[t_index] = failures[t_index] or reason

    for source in records:
        row = _mapping(source)
        if row.get("nommcuts") is not True:
            continue
        raw_t_index = row.get("canonical_t_index")
        raw_delta_index = row.get("delta_index")
        if raw_t_index is None:
            continue
        t_index = _d10_integer(raw_t_index)
        if t_index is None or not 0 <= t_index < t_count:
            fail(None, "phase_a_record_t_membership_invalid")
            continue
        if raw_delta_index is None:
            continue
        delta_index = _d10_integer(raw_delta_index)
        if delta_index is None or not 0 <= delta_index < delta_count:
            fail(t_index, "phase_a_record_delta_membership_invalid")
            continue
        missing_mass = _d10_finite(row.get("analysis_MM"))
        contribution = _d10_finite(row.get(contribution_key))
        if missing_mass is None or contribution is None:
            fail(t_index, "phase_a_record_scalar_invalid")
            continue
        if host_state is not None and row.get("host_state") != host_state:
            fail(t_index, "phase_a_host_record_state_invalid")
            continue
        copied[t_index].append({
            "t_index": t_index,
            "delta_index": delta_index,
            "missing_mass": missing_mass,
            "signed_contribution": contribution,
        })
    return tuple(tuple(rows) for rows in copied), tuple(failures)


def _d10_copy_cell_bins(value, mm_edges):
    if not _d10_sequence(value) or len(value) != len(mm_edges) + 1:
        return None, "method_b_cell_bins_invalid"
    copied = []
    regular_count = len(mm_edges) - 1
    for expected_index, source in enumerate(value):
        row = _mapping(source)
        required = (
            "index", "mm_low", "mm_high", "underflow", "overflow",
            "in_allowed_pion_sensitive_domain", "usable_for_shape",
            "host_record_count", "host_yield", "host_sumw2",
            "baseline_record_count", "baseline_yield", "baseline_sumw2",
            "residual", "pull",
        )
        if not row or any(key not in row for key in required):
            return None, "method_b_cell_bins_invalid"
        index = _d10_integer(row["index"])
        if index != expected_index:
            return None, "method_b_cell_bins_invalid"
        regular = 1 <= index <= regular_count
        if (
            row["underflow"] is not (index == 0)
            or row["overflow"] is not (index == regular_count + 1)
            or not isinstance(row["in_allowed_pion_sensitive_domain"], bool)
            or not isinstance(row["usable_for_shape"], bool)
        ):
            return None, "method_b_cell_bins_invalid"
        if regular:
            if row["mm_low"] != mm_edges[index - 1] or row["mm_high"] != mm_edges[index]:
                return None, "method_b_cell_mm_geometry_invalid"
        elif row["mm_low"] is not None or row["mm_high"] is not None:
            return None, "method_b_cell_mm_geometry_invalid"
        host_count = _d10_integer(row["host_record_count"])
        baseline_count = _d10_integer(row["baseline_record_count"])
        host_yield = _d10_finite(row["host_yield"])
        host_sumw2 = _d10_finite(row["host_sumw2"])
        baseline_yield = _d10_finite(row["baseline_yield"])
        baseline_sumw2 = _d10_finite(row["baseline_sumw2"])
        if (
            host_count is None or host_count < 0 or baseline_count is None or baseline_count < 0
            or host_yield is None or host_sumw2 is None or host_sumw2 < 0.0
            or baseline_yield is None or baseline_sumw2 is None or baseline_sumw2 < 0.0
        ):
            return None, "method_b_cell_bins_invalid"
        residual = row["residual"]
        pull = row["pull"]
        if (residual is not None and _d10_finite(residual) is None) or (
            pull is not None and _d10_finite(pull) is None
        ):
            return None, "method_b_cell_bins_invalid"
        copied.append({
            "index": index,
            "mm_low": row["mm_low"] if not regular else float(row["mm_low"]),
            "mm_high": row["mm_high"] if not regular else float(row["mm_high"]),
            "underflow": row["underflow"],
            "overflow": row["overflow"],
            "in_allowed_pion_sensitive_domain": row["in_allowed_pion_sensitive_domain"],
            "usable_for_shape": row["usable_for_shape"],
            "host_record_count": host_count,
            "host_yield": host_yield,
            "host_sumw2": host_sumw2,
            "baseline_record_count": baseline_count,
            "baseline_yield": baseline_yield,
            "baseline_sumw2": baseline_sumw2,
            "residual": None if residual is None else float(residual),
            "pull": None if pull is None else float(pull),
        })
    return tuple(copied), None


def _d10_local_closure_cells(cells, method, phase):
    t_edges = method["t_edges"]
    delta_edges = method["delta_edges"]
    t_count = len(t_edges) - 1
    delta_count = len(delta_edges) - 1
    if len(cells) != t_count * delta_count:
        return None, None, "method_b_cell_grid_invalid"
    grouped = [[] for _unused in range(t_count)]
    failures = [None for _unused in range(t_count)]
    seen = set()
    for source in cells:
        cell = _mapping(source)
        t_index = _d10_integer(cell.get("t_index"))
        delta_index = _d10_integer(cell.get("delta_index"))
        if (
            t_index is None or delta_index is None
            or not 0 <= t_index < t_count or not 0 <= delta_index < delta_count
            or (t_index, delta_index) in seen
            or cell.get("t_low") != t_edges[t_index]
            or cell.get("t_high") != t_edges[t_index + 1]
            or cell.get("delta_low") != delta_edges[delta_index]
            or cell.get("delta_high") != delta_edges[delta_index + 1]
            or cell.get("host_state") != phase["host_state"]
        ):
            return None, None, "method_b_cell_grid_invalid"
        seen.add((t_index, delta_index))
        if cell.get("mm_edges") != method["mm_edges"]:
            failures[t_index] = failures[t_index] or "method_b_cell_mm_geometry_invalid"
            continue
        bins, reason = _d10_copy_cell_bins(cell.get("bins"), method["mm_edges"])
        if reason is not None:
            failures[t_index] = failures[t_index] or reason
            continue
        shape_status = cell.get("shape_status")
        if shape_status not in {"good", "marginal", "poor", "unavailable"}:
            failures[t_index] = failures[t_index] or "method_b_cell_shape_invalid"
            continue
        grouped[t_index].append({
            "t_index": t_index,
            "delta_index": delta_index,
            "delta_low": float(delta_edges[delta_index]),
            "delta_high": float(delta_edges[delta_index + 1]),
            "mm_edges": list(method["mm_edges"]),
            "bins": bins,
            "shape_status": shape_status,
            "shape_reason": cell.get("shape_reason"),
        })
    if len(seen) != t_count * delta_count:
        return None, None, "method_b_cell_grid_invalid"
    for t_index, rows in enumerate(grouped):
        if len(rows) != delta_count:
            failures[t_index] = failures[t_index] or "method_b_cell_grid_invalid"
        rows.sort(key=lambda row: row["delta_index"])
    return tuple(tuple(rows) for rows in grouped), tuple(failures), None


def _d10_d3_relative_cells(comparison, phase, method):
    d3 = _mapping(comparison)
    if not d3 or d3.get("status") != "available" or d3.get("available") is not True:
        return None, "method_b_comparison_unavailable"
    required = (
        "schema_version", "non_authoritative", "method_a_numerical_dependency",
        "comparison_performed", "classification_performed", "production_objects_mutated",
        "refinement_applied", "phase_a_contract_fingerprint", "coordinate_fingerprint",
        "method_b_fingerprint", "canonical_t_edges", "delta_edges", "host_state",
        "source_target_state", "cells",
    )
    if any(key not in d3 for key in required):
        return None, "method_b_comparison_contract_invalid"
    if (
        d3["schema_version"] != _D10_METHOD_B_COMPARISON_SCHEMA
        or d3["non_authoritative"] is not True
        or d3["method_a_numerical_dependency"] is not False
        or d3["comparison_performed"] is not False
        or d3["classification_performed"] is not False
        or d3["production_objects_mutated"] is not False
        or d3["refinement_applied"] is not False
        or d3["phase_a_contract_fingerprint"] != phase["contract_fingerprint"]
        or d3["coordinate_fingerprint"] != phase["coordinate_fingerprint"]
        or d3["method_b_fingerprint"] != method["fingerprint"]
        or d3["canonical_t_edges"] != method["t_edges"]
        or d3["delta_edges"] != method["delta_edges"]
        or d3["host_state"] != phase["host_state"]
        or d3["source_target_state"] != phase["source_target_state"]
        or not _d10_sequence(d3["cells"])
    ):
        return None, "method_b_comparison_provenance_invalid"
    t_count = len(method["t_edges"]) - 1
    delta_count = len(method["delta_edges"]) - 1
    if len(d3["cells"]) != t_count * delta_count:
        return None, "method_b_comparison_cell_grid_invalid"
    grouped = [[] for _unused in range(t_count)]
    seen = set()
    statuses = {
        "available_multi_region", "single_region_only", "unavailable", "region_marginal",
        "region_inconsistent", "shape_poor_veto",
    }
    for source in d3["cells"]:
        cell = _mapping(source)
        required_cell = (
            "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
            "method_b_comparison_candidate", "method_b_comparison_candidate_uncertainty",
            "method_b_comparison_candidate_status", "method_b_comparison_candidate_reason",
            "method_B_status", "method_B_reason", "region_consistency_status",
            "region_consistency_reason", "shape_status", "shape_reason",
        )
        if not cell or any(key not in cell for key in required_cell):
            return None, "method_b_comparison_cell_contract_invalid"
        t_index = _d10_integer(cell["t_index"])
        delta_index = _d10_integer(cell["delta_index"])
        if (
            t_index is None or delta_index is None
            or not 0 <= t_index < t_count or not 0 <= delta_index < delta_count
            or (t_index, delta_index) in seen
            or cell["t_low"] != method["t_edges"][t_index]
            or cell["t_high"] != method["t_edges"][t_index + 1]
            or cell["delta_low"] != method["delta_edges"][delta_index]
            or cell["delta_high"] != method["delta_edges"][delta_index + 1]
            or cell["method_b_comparison_candidate_status"] not in statuses
            or cell["method_B_status"] not in {
                "available", "marginal", "internally_inconsistent", "shape_inconsistent", "unavailable",
            }
            or cell["region_consistency_status"] not in {
                "region_consistent", "region_marginal", "region_inconsistent", "insufficient_regions",
            }
            or cell["shape_status"] not in {"good", "marginal", "poor", "unavailable"}
        ):
            return None, "method_b_comparison_cell_contract_invalid"
        status = cell["method_b_comparison_candidate_status"]
        candidate = cell["method_b_comparison_candidate"]
        uncertainty = cell["method_b_comparison_candidate_uncertainty"]
        if status == "available_multi_region":
            candidate = _d10_finite(candidate)
            uncertainty = _d10_finite(uncertainty)
            if candidate is None or candidate <= 0.0 or uncertainty is None or uncertainty <= 0.0:
                return None, "method_b_comparison_candidate_invalid"
        elif candidate is not None or uncertainty is not None:
            return None, "method_b_comparison_candidate_invalid"
        seen.add((t_index, delta_index))
        grouped[t_index].append({
            "t_index": t_index,
            "t_low": float(cell["t_low"]),
            "t_high": float(cell["t_high"]),
            "delta_index": delta_index,
            "delta_low": float(cell["delta_low"]),
            "delta_high": float(cell["delta_high"]),
            "method_b_comparison_candidate": candidate,
            "method_b_comparison_candidate_uncertainty": uncertainty,
            "method_b_comparison_candidate_status": status,
            "method_b_comparison_candidate_reason": cell[
                "method_b_comparison_candidate_reason"
            ],
            "method_B_status": cell["method_B_status"],
            "method_B_reason": cell["method_B_reason"],
            "region_consistency_status": cell["region_consistency_status"],
            "region_consistency_reason": cell["region_consistency_reason"],
            "shape_status": cell["shape_status"],
            "shape_reason": cell["shape_reason"],
        })
    if len(seen) != t_count * delta_count:
        return None, "method_b_comparison_cell_grid_invalid"
    return tuple(tuple(sorted(rows, key=lambda row: row["delta_index"])) for rows in grouped), None


def build_full_background_subtraction_d10_payload(
    pion_hgcer_event_contract,
    pion_hgcer_method_b,
    pion_hgcer_method_b_comparison,
):
    """Select frozen Method-B explanation inputs without Method-B recalculation."""
    phase, reason = _d10_phase_a_contract(pion_hgcer_event_contract)
    if reason is not None:
        return _d10_unavailable(reason)
    method, reason = _d10_method_b_contract(pion_hgcer_method_b, phase)
    if reason is not None:
        return _d10_unavailable(reason)
    local_cells, local_failures, reason = _d10_local_closure_cells(
        method["cells"], method, phase
    )
    if reason is not None:
        return _d10_unavailable(reason)
    t_count = len(phase["t_edges"]) - 1
    delta_count = len(phase["delta_edges"]) - 1
    host_rows, host_failures = _d10_event_rows(
        phase["kaon_host_records"],
        contribution_key="signed_host_event_contribution",
        t_count=t_count,
        delta_count=delta_count,
        host_state=phase["host_state"],
    )
    pion_rows, pion_failures = _d10_event_rows(
        phase["pion_records"],
        contribution_key="signed_baseline_event_contribution",
        t_count=t_count,
        delta_count=delta_count,
    )
    relative_cells, d3_reason = _d10_d3_relative_cells(
        pion_hgcer_method_b_comparison, phase, method
    )
    host_label = (
        "Proton-cleaned kaon data"
        if phase["host_state"] == "proton_cleaned" else "Kaon-selected data"
    )
    per_t = []
    for t_index in range(t_count):
        inputs_reason = host_failures[t_index] or pion_failures[t_index]
        per_t.append({
            "t_index": t_index,
            "t_low": float(phase["t_edges"][t_index]),
            "t_high": float(phase["t_edges"][t_index + 1]),
            "mm_inputs": {
                "available": inputs_reason is None,
                "reason": inputs_reason,
                "host_rows": host_rows[t_index],
                "baseline_rows": pion_rows[t_index],
            },
            "local_closure": {
                "available": local_failures[t_index] is None,
                "reason": local_failures[t_index],
                "cells_by_delta": local_cells[t_index],
            },
            "method_b_relative": {
                "available": d3_reason is None,
                "reason": d3_reason,
                "cells": relative_cells[t_index] if relative_cells is not None else tuple(),
            },
        })
    return {
        "schema_version": D10_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "phase_a_contract_fingerprint": phase["contract_fingerprint"],
        "coordinate_fingerprint": phase["coordinate_fingerprint"],
        "method_b_fingerprint": method["fingerprint"],
        "host_state": phase["host_state"],
        "host_label": host_label,
        "source_target_state": phase["source_target_state"],
        "t_edges": list(phase["t_edges"]),
        "delta_edges": list(phase["delta_edges"]),
        "mm_edges": list(method["mm_edges"]),
        "mm_regions": [dict(row) for row in method["mm_regions"]],
        "protected_regions": [dict(row) for row in method["protected_regions"]],
        "d3_available": d3_reason is None,
        "d3_reason": d3_reason,
        "per_t": per_t,
    }


def _d10_new_histogram(ROOT, name, edges):
    if not hasattr(ROOT, "TH1D"):
        return None
    try:
        histogram = ROOT.TH1D(str(name), str(name), len(edges) - 1, array("d", edges))
        histogram.SetDirectory(0)
        if hasattr(histogram, "Sumw2"):
            histogram.Sumw2()
        return histogram
    except Exception:
        return None


def _d10_event_histogram(ROOT, name, edges, rows):
    histogram = _d10_new_histogram(ROOT, name, edges)
    if histogram is None:
        return None
    try:
        bin_count = len(edges) - 1
        contents = [0.0] * (bin_count + 2)
        sumw2 = [0.0] * (bin_count + 2)
        for row in rows:
            row = _mapping(row)
            value = float(row["missing_mass"])
            contribution = float(row["signed_contribution"])
            if value < edges[0]:
                bin_index = 0
            elif value > edges[-1]:
                bin_index = bin_count + 1
            elif value == edges[-1]:
                bin_index = bin_count
            else:
                bin_index = next(
                    index
                    for index in range(1, bin_count + 1)
                    if value < edges[index]
                )
            contents[bin_index] += contribution
            sumw2[bin_index] += contribution * contribution
        for bin_index, (content, variance) in enumerate(zip(contents, sumw2)):
            histogram.SetBinContent(bin_index, content)
            histogram.SetBinError(bin_index, math.sqrt(variance))
    except Exception:
        return None
    return histogram


def _d10_local_histograms(ROOT, name, cell):
    cell = _mapping(cell)
    edges = cell.get("mm_edges") or ()
    host = _d10_new_histogram(ROOT, "{}_host".format(name), edges)
    baseline = _d10_new_histogram(ROOT, "{}_baseline".format(name), edges)
    if host is None or baseline is None:
        return None, None
    try:
        for row in cell.get("bins") or ():
            row = _mapping(row)
            index = int(row["index"])
            if not 1 <= index <= len(edges) - 1:
                continue
            host.SetBinContent(index, float(row["host_yield"]))
            host.SetBinError(index, math.sqrt(float(row["host_sumw2"])))
            baseline.SetBinContent(index, float(row["baseline_yield"]))
            baseline.SetBinError(index, math.sqrt(float(row["baseline_sumw2"])))
    except Exception:
        return None, None
    return host, baseline


def _d10_draw_region_boundaries(ROOT, presentation, y_range, draw_objects):
    """Draw only frozen Method-B MM boundaries, never recomputed windows."""
    if not hasattr(ROOT, "TLine"):
        return {}
    lines = {}
    for region in tuple(presentation.get("protected_regions") or ()):
        region = _mapping(region)
        if region.get("region_name") != "KLambdaSigma0":
            continue
        for boundary in (region.get("mm_low"), region.get("mm_high")):
            try:
                line = ROOT.TLine(
                    float(boundary), float(y_range[0]), float(boundary), float(y_range[1])
                )
                line.SetLineColor(getattr(ROOT, "kRed", 2))
                line.SetLineStyle(3)
                line.SetLineWidth(1)
                line.Draw()
                draw_objects.append(line)
                lines.setdefault("protected", line)
            except Exception:
                continue
    return lines


def _d10_add_region_legend_entries(legend, lines):
    if "protected" in lines:
        legend.AddEntry(
            lines["protected"], "Pion-sensitive region (outside protected window)", "l"
        )
        legend.AddEntry(
            lines["protected"], "Lambda/Sigma protected region (1.10-1.23 GeV)", "l"
        )


def _render_d10_mm_inputs_page(ROOT, pdf_name, presentation, group):
    inputs = _mapping(group.get("mm_inputs"))
    edges = presentation.get("mm_edges") or ()
    host = _d10_event_histogram(
        ROOT,
        "H_full_background_d10_inputs_host_t{}".format(group["t_index"] + 1),
        edges,
        inputs.get("host_rows") or (),
    )
    baseline = _d10_event_histogram(
        ROOT,
        "H_full_background_d10_inputs_baseline_t{}".format(group["t_index"] + 1),
        edges,
        inputs.get("baseline_rows") or (),
    )
    if host is None or baseline is None:
        return False
    title = "Method B - Missing-mass closure inputs - {}".format(_t_context(group))
    _set_histogram_title(host, ";Missing mass [GeV];Signed normalized yield")
    y_range = _combined_histogram_y_range((host, baseline))
    _apply_display_y_range(host, y_range)
    _style_histogram(host, getattr(ROOT, "kBlack", 1))
    _style_histogram(baseline, getattr(ROOT, "kBlue", 4))
    canvas = ROOT.TCanvas(
        "C_full_background_d10_inputs_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    draw_objects = [host, baseline]
    try:
        host.Draw("hist e")
        baseline.Draw("hist e same")
        lines = _d10_draw_region_boundaries(ROOT, presentation, y_range, draw_objects)
        legend = ROOT.TLegend(0.55, 0.56, 0.89, 0.82)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(host, presentation.get("host_label", "Kaon-selected data"), "l")
        legend.AddEntry(baseline, "Baseline pion background", "l")
        _d10_add_region_legend_entries(legend, lines)
        legend.Draw()
        draw_objects.append(legend)
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Method B - Missing-mass closure inputs", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _d10_local_closure_y_range(histograms):
    """Return a signed display range from one local closure panel only."""
    values = [0.0]
    for histogram in tuple(histograms or ()):
        try:
            bin_count = int(histogram.GetNbinsX())
        except Exception:
            continue
        for bin_index in range(1, bin_count + 1):
            try:
                content = float(histogram.GetBinContent(bin_index))
                error = abs(float(histogram.GetBinError(bin_index)))
            except Exception:
                continue
            if math.isfinite(content) and math.isfinite(error):
                values.extend((content - error, content + error))
    low, high = min(values), max(values)
    if high <= low:
        return -0.05, 0.05
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def _render_d10_local_closure_page(ROOT, pdf_name, presentation, group):
    closure = _mapping(group.get("local_closure"))
    cells = tuple(closure.get("cells_by_delta") or ())
    if not cells:
        return False
    columns = min(3, len(cells))
    rows = int(math.ceil(float(len(cells)) / float(columns)))
    title = "Method B - Local missing-mass closure - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d10_local_t{}".format(group["t_index"] + 1), title, 1400, 900
    )
    canvas.Divide(columns, rows)
    draw_objects = []
    pairs = []
    try:
        for cell in cells:
            host, baseline = _d10_local_histograms(
                ROOT,
                "H_full_background_d10_local_t{}_d{}".format(
                    group["t_index"] + 1, int(cell["delta_index"]) + 1
                ),
                cell,
            )
            if host is None or baseline is None:
                return False
            pairs.append((cell, host, baseline))
        for panel_index, (cell, host, baseline) in enumerate(pairs, 1):
            canvas.cd(panel_index)
            panel_title = "delta = [{:.3f}, {:.3f}] %".format(
                float(cell["delta_low"]), float(cell["delta_high"])
            )
            _set_histogram_title(
                host, "{};Missing mass [GeV];Signed normalized yield".format(panel_title)
            )
            y_range = _d10_local_closure_y_range((host, baseline))
            _apply_display_y_range(host, y_range)
            _style_histogram(host, getattr(ROOT, "kBlack", 1))
            _style_histogram(baseline, getattr(ROOT, "kBlue", 4))
            host.Draw("hist e")
            baseline.Draw("hist e same")
            lines = _d10_draw_region_boundaries(ROOT, presentation, y_range, draw_objects)
            if panel_index == 1:
                legend = ROOT.TLegend(0.48, 0.54, 0.89, 0.84)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(host, presentation.get("host_label", "Kaon-selected data"), "l")
                legend.AddEntry(baseline, "Baseline pion background", "l")
                _d10_add_region_legend_entries(legend, lines)
                legend.Draw()
                draw_objects.append(legend)
            draw_objects.extend((host, baseline))
        canvas.cd(1)
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Method B - Local missing-mass closure", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _d10_relative_y_range(cells):
    values = [1.0]
    for cell in cells:
        cell = _mapping(cell)
        if cell.get("method_b_comparison_candidate_status") != "available_multi_region":
            continue
        candidate = _d10_finite(cell.get("method_b_comparison_candidate"))
        uncertainty = _d10_finite(cell.get("method_b_comparison_candidate_uncertainty"))
        if candidate is not None and uncertainty is not None:
            values.extend((candidate - uncertainty, candidate, candidate + uncertainty))
    low, high = min(values), max(values)
    if high <= low:
        return low - 0.05, high + 0.05
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def _render_d10_method_b_relative_page(ROOT, pdf_name, group, delta_edges):
    presentation = _mapping(group.get("method_b_relative"))
    cells = tuple(presentation.get("cells") or ())
    if not hasattr(ROOT, "TH1D"):
        return False
    title = "Method B - Missing-mass closure diagnostic - {}".format(_t_context(group))
    try:
        frame = ROOT.TH1D(
            "H_full_background_d10_relative_t{}".format(group["t_index"] + 1),
            "{};delta [%];Relative pion-background diagnostic".format(title),
            len(delta_edges) - 1,
            array("d", delta_edges),
        )
        frame.SetDirectory(0)
        if hasattr(frame, "SetStats"):
            frame.SetStats(0)
        y_range = _d10_relative_y_range(cells)
        frame.SetMinimum(y_range[0])
        frame.SetMaximum(y_range[1])
    except Exception:
        return False
    canvas = ROOT.TCanvas(
        "C_full_background_d10_relative_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    draw_objects = [frame]
    try:
        frame.Draw("axis")
        selected = [
            _mapping(cell) for cell in cells
            if _mapping(cell).get("method_b_comparison_candidate_status") == "available_multi_region"
        ]
        if selected and hasattr(ROOT, "TGraphErrors"):
            graph = ROOT.TGraphErrors(len(selected))
            for point_index, cell in enumerate(selected):
                x_value = 0.5 * (float(cell["delta_low"]) + float(cell["delta_high"]))
                graph.SetPoint(point_index, x_value, float(cell["method_b_comparison_candidate"]))
                graph.SetPointError(
                    point_index, 0.0, float(cell["method_b_comparison_candidate_uncertainty"])
                )
            graph.SetMarkerColor(getattr(ROOT, "kBlack", 1))
            graph.SetLineColor(getattr(ROOT, "kBlack", 1))
            graph.SetMarkerStyle(20)
            graph.SetMarkerSize(1.1)
            graph.Draw("P same")
            draw_objects.append(graph)
            legend = ROOT.TLegend(0.64, 0.74, 0.89, 0.86)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(graph, "Available local closure", "p")
            legend.Draw()
            draw_objects.append(legend)
        if hasattr(ROOT, "TLine"):
            unity = ROOT.TLine(float(delta_edges[0]), 1.0, float(delta_edges[-1]), 1.0)
            unity.SetLineStyle(2)
            unity.SetLineWidth(2)
            unity.Draw()
            draw_objects.append(unity)
        draw_objects.append(
            _draw_page_header(ROOT, canvas, "Method B - Missing-mass closure diagnostic", group)
        )
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d10_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render the detached D.10 Method-B explanation for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    inputs = _mapping(group.get("mm_inputs"))
    if inputs.get("available"):
        if _render_d10_mm_inputs_page(ROOT, pdf_name, presentation, group):
            manifest.append({"page_id": "full_background.d10.method_b_mm_inputs", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.10 Method-B MM-input page unavailable for t{}".format(t_number))
    else:
        failures.append("D.10 Method-B MM-input unavailable for t{}: {}".format(t_number, inputs.get("reason")))
    closure = _mapping(group.get("local_closure"))
    if closure.get("available"):
        if _render_d10_local_closure_page(ROOT, pdf_name, presentation, group):
            manifest.append({"page_id": "full_background.d10.method_b_local_closure", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.10 Method-B local-closure page unavailable for t{}".format(t_number))
    else:
        failures.append("D.10 Method-B local-closure input unavailable for t{}: {}".format(t_number, closure.get("reason")))
    relative = _mapping(group.get("method_b_relative"))
    if relative.get("available"):
        if _render_d10_method_b_relative_page(
            ROOT, pdf_name, group, presentation.get("delta_edges") or ()
        ):
            manifest.append({"page_id": "full_background.d10.method_b_relative", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.10 Method-B relative page unavailable for t{}".format(t_number))
    else:
        failures.append("D.10 Method-B relative input unavailable for t{}: {}".format(t_number, relative.get("reason")))


def render_full_background_subtraction_d6_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.6 pages in raw-MM, PID, weight order for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.6 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.6 procedure rendering unavailable: PyROOT not available")
        return result
    for group in tuple(presentation.get("per_t") or ()):
        _render_d6_t_pages(ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"])
    return result


def render_full_background_subtraction_d7_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.7 proton-subtraction pages for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.7 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.7 procedure rendering unavailable: PyROOT not available")
        return result
    _append_d7_exclusion_failure(presentation, result["failures"])
    for group in tuple(presentation.get("per_t") or ()):
        _render_d7_t_pages(ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"])
    return result


def render_full_background_subtraction_d8_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.8 baseline-pion pages for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.8 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.8 procedure rendering unavailable: PyROOT not available")
        return result
    for group in tuple(presentation.get("per_t") or ()):
        _render_d8_t_pages(
            ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"]
        )
    return result


def render_full_background_subtraction_d9_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.9 HGCer/Method-A pages for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.9 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.9 procedure rendering unavailable: PyROOT not available")
        return result
    for group in tuple(presentation.get("per_t") or ()):
        _render_d9_t_pages(
            ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"]
        )
    return result


def render_full_background_subtraction_d10_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.10 Method-B explanation pages for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.10 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.10 procedure rendering unavailable: PyROOT not available")
        return result
    for group in tuple(presentation.get("per_t") or ()):
        _render_d10_t_pages(
            ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"]
        )
    return result


def _d11_nonempty_string(value):
    return isinstance(value, str) and bool(value.strip())


def _d11_integer(value):
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if value == number else None


def _d11_finite(value):
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _d11_canonical_fingerprint(value):
    """Return the frozen comparison hash, or ``None`` for non-JSON data."""
    try:
        serialized = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError, OverflowError):
        return None
    return hashlib.sha256(serialized.encode("ascii")).hexdigest()


def _d11_serialized_equal(left, right):
    """Require the exact JSON geometry retained by the frozen checkpoint."""
    try:
        return json.dumps(
            left, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False,
        ) == json.dumps(
            right, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False,
        )
    except (TypeError, ValueError, OverflowError):
        return False


def _d11_checkpoint_contract(value):
    checkpoint = _mapping(value)
    if not checkpoint:
        return None, "phase_d_checkpoint_contract_invalid"
    if checkpoint.get("status") != "available" or checkpoint.get("available") is not True:
        return None, "phase_d_checkpoint_unavailable"
    required = (
        "schema_version", "source_checkpoint_payload_fingerprint", "non_authoritative",
        "comparison_performed", "classification_performed", "classification_scope",
        "decision_performed", "statistical_compatibility_claimed",
        "production_objects_mutated", "refinement_applied", "method_a_comparison",
        "method_b_comparison", "ab_comparison",
    )
    if any(key not in checkpoint for key in required):
        return None, "phase_d_checkpoint_contract_invalid"
    if (
        checkpoint["schema_version"] != _D11_PHASE_D_CHECKPOINT_SCHEMA
        or not _d11_nonempty_string(checkpoint["source_checkpoint_payload_fingerprint"])
        or checkpoint["non_authoritative"] is not True
        or checkpoint["comparison_performed"] is not True
        or checkpoint["classification_performed"] is not True
        or checkpoint["classification_scope"] != "availability_only_non_prescriptive"
        or checkpoint["decision_performed"] is not False
        or checkpoint["statistical_compatibility_claimed"] is not False
        or checkpoint["production_objects_mutated"] is not False
        or checkpoint["refinement_applied"] is not False
    ):
        return None, "phase_d_checkpoint_authority_invalid"
    return checkpoint, None


def _d11_representation_contract(
    value, schema, dependency_key, source_payload_key, reason
):
    representation = _mapping(value)
    required = (
        "schema_version", "status", "available", "source_checkpoint_payload_fingerprint",
        "phase_a_contract_fingerprint", "coordinate_fingerprint", source_payload_key,
        "canonical_t_edges", "delta_edges", "fingerprint_inputs", "fingerprint",
        "non_authoritative", dependency_key,
        "comparison_performed", "classification_performed", "production_objects_mutated",
        "refinement_applied",
    )
    if not representation or any(key not in representation for key in required):
        return None, reason
    t_edges = _strict_edges(representation["canonical_t_edges"])
    delta_edges = _strict_edges(representation["delta_edges"])
    if (
        representation["schema_version"] != schema
        or representation["status"] != "available"
        or representation["available"] is not True
        or not all(
            _d11_nonempty_string(representation[key])
            for key in (
                "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
                "coordinate_fingerprint", source_payload_key, "fingerprint",
            )
        )
        or t_edges is None
        or delta_edges is None
        or not isinstance(representation["fingerprint_inputs"], Mapping)
        or representation["fingerprint"] != _d11_canonical_fingerprint(
            representation["fingerprint_inputs"]
        )
        or representation["non_authoritative"] is not True
        or representation[dependency_key] is not False
        or representation["comparison_performed"] is not False
        or representation["classification_performed"] is not False
        or representation["production_objects_mutated"] is not False
        or representation["refinement_applied"] is not False
    ):
        return None, reason
    return representation, None


def _d11_ab_contract(value):
    comparison = _mapping(value)
    required = (
        "schema_version", "method", "status", "available",
        "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
        "coordinate_fingerprint", "method_a_comparison_fingerprint",
        "method_b_comparison_fingerprint", "source_method_a_comparison_payload_fingerprint",
        "source_method_b_comparison_payload_fingerprint", "canonical_t_edges", "delta_edges",
        "host_state", "source_target_state", "cells", "fingerprint_inputs", "fingerprint",
        "non_authoritative",
        "comparison_performed", "classification_performed", "classification_scope",
        "decision_performed", "statistical_compatibility_claimed",
        "production_objects_mutated", "refinement_applied",
    )
    if not comparison or any(key not in comparison for key in required):
        return None, "ab_comparison_contract_invalid"
    t_edges = _strict_edges(comparison["canonical_t_edges"])
    delta_edges = _strict_edges(comparison["delta_edges"])
    if (
        comparison["schema_version"] != _D11_AB_COMPARISON_SCHEMA
        or comparison["method"] != "non_authoritative_ab_comparison"
        or comparison["status"] != "available"
        or comparison["available"] is not True
        or not all(
            _d11_nonempty_string(comparison[key])
            for key in (
                "source_checkpoint_payload_fingerprint", "phase_a_contract_fingerprint",
                "coordinate_fingerprint", "method_a_comparison_fingerprint",
                "method_b_comparison_fingerprint",
                "source_method_a_comparison_payload_fingerprint",
                "source_method_b_comparison_payload_fingerprint", "fingerprint",
            )
        )
        or t_edges is None
        or delta_edges is None
        or comparison["host_state"] not in {"proton_cleaned", "identity_no_proton_cleaning"}
        or comparison["source_target_state"] != _D11_SOURCE_TARGET_STATE
        or not isinstance(comparison["cells"], Sequence)
        or isinstance(comparison["cells"], (str, bytes))
        or not isinstance(comparison["fingerprint_inputs"], Mapping)
        or comparison["fingerprint"] != _d11_canonical_fingerprint(
            comparison["fingerprint_inputs"]
        )
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


def _d11_source_linkage(checkpoint, method_a, method_b, comparison):
    if (
        method_a["fingerprint"] != comparison["method_a_comparison_fingerprint"]
        or method_b["fingerprint"] != comparison["method_b_comparison_fingerprint"]
    ):
        return "ab_comparison_representation_fingerprint_mismatch"
    if (
        _d11_canonical_fingerprint(method_a)
        != comparison["source_method_a_comparison_payload_fingerprint"]
        or _d11_canonical_fingerprint(method_b)
        != comparison["source_method_b_comparison_payload_fingerprint"]
    ):
        return "ab_comparison_representation_payload_fingerprint_mismatch"
    source_fingerprint = checkpoint["source_checkpoint_payload_fingerprint"]
    if any(
        source.get("source_checkpoint_payload_fingerprint") != source_fingerprint
        for source in (method_a, method_b, comparison)
    ):
        return "ab_comparison_source_checkpoint_fingerprint_mismatch"
    if any(
        source.get("phase_a_contract_fingerprint") != comparison["phase_a_contract_fingerprint"]
        or source.get("coordinate_fingerprint") != comparison["coordinate_fingerprint"]
        for source in (method_a, method_b)
    ):
        return "ab_comparison_provenance_mismatch"
    if any(not _d11_serialized_equal(source[key], comparison[key]) for source in (
        method_a, method_b
    ) for key in ("canonical_t_edges", "delta_edges")):
        return "ab_comparison_geometry_mismatch"
    if (
        method_b.get("host_state") != comparison["host_state"]
        or method_b.get("source_target_state") != comparison["source_target_state"]
    ):
        return "ab_comparison_host_state_mismatch"
    return None


def _d11_scalar_cell(source, t_edges, delta_edges):
    cell = _mapping(source)
    required = (
        "t_index", "t_low", "t_high", "delta_index", "delta_low", "delta_high",
        "method_a", "method_b", "comparison",
    )
    if not cell or any(key not in cell for key in required):
        return None, "ab_comparison_cell_contract_invalid"
    t_index = _d11_integer(cell["t_index"])
    delta_index = _d11_integer(cell["delta_index"])
    if (
        t_index is None or delta_index is None
        or not 0 <= t_index < len(t_edges) - 1
        or not 0 <= delta_index < len(delta_edges) - 1
        or not _d11_serialized_equal(cell["t_low"], t_edges[t_index])
        or not _d11_serialized_equal(cell["t_high"], t_edges[t_index + 1])
        or not _d11_serialized_equal(cell["delta_low"], delta_edges[delta_index])
        or not _d11_serialized_equal(cell["delta_high"], delta_edges[delta_index + 1])
    ):
        return None, "ab_comparison_cell_geometry_invalid"
    method_a = _mapping(cell["method_a"])
    method_b = _mapping(cell["method_b"])
    relation = _mapping(cell["comparison"])
    if any(
        key not in method_a
        for key in (
            "present", "comparison_candidate", "comparison_candidate_low",
            "comparison_candidate_high", "comparison_candidate_status",
        )
    ) or any(
        key not in method_b
        for key in (
            "present", "comparison_candidate", "comparison_candidate_uncertainty",
            "comparison_candidate_status",
        )
    ) or any(
        key not in relation
        for key in ("availability", "availability_reason", "ratio_B_over_A", "log_ratio_B_over_A")
    ):
        return None, "ab_comparison_cell_contract_invalid"
    a_present = method_a["present"]
    b_present = method_b["present"]
    if not isinstance(a_present, bool) or not isinstance(b_present, bool):
        return None, "ab_comparison_cell_contract_invalid"
    a_candidate = _d11_finite(method_a["comparison_candidate"])
    a_low = _d11_finite(method_a["comparison_candidate_low"])
    a_high = _d11_finite(method_a["comparison_candidate_high"])
    if a_present:
        if (
            method_a["comparison_candidate_status"] not in {"available", "marginal"}
            or None in (a_candidate, a_low, a_high)
            or a_candidate < 0.0 or a_low < 0.0 or a_low > a_candidate or a_candidate > a_high
        ):
            return None, "ab_comparison_method_a_cell_invalid"
    elif (
        method_a["comparison_candidate_status"] != "unavailable"
        or any(value is not None for value in (method_a["comparison_candidate"], method_a["comparison_candidate_low"], method_a["comparison_candidate_high"]))
    ):
        return None, "ab_comparison_method_a_cell_invalid"
    b_candidate = _d11_finite(method_b["comparison_candidate"])
    b_uncertainty = _d11_finite(method_b["comparison_candidate_uncertainty"])
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
        or any(value is not None for value in (method_b["comparison_candidate"], method_b["comparison_candidate_uncertainty"]))
    ):
        return None, "ab_comparison_method_b_cell_invalid"
    availability = relation["availability"]
    ratio = relation["ratio_B_over_A"]
    log_ratio = relation["log_ratio_B_over_A"]
    if availability not in _D11_AVAILABILITY_LABELS:
        return None, "ab_comparison_availability_invalid"
    if availability == "both_comparable":
        if not a_present or not b_present or a_candidate <= 0.0 or (
            _d11_finite(ratio) is None or _d11_finite(log_ratio) is None
        ):
            return None, "ab_comparison_availability_invalid"
    elif availability == "both_present_not_comparable":
        if not a_present or not b_present or a_candidate != 0.0 or ratio is not None or log_ratio is not None:
            return None, "ab_comparison_availability_invalid"
    elif availability == "a_only":
        if not a_present or b_present or ratio is not None or log_ratio is not None:
            return None, "ab_comparison_availability_invalid"
    elif availability == "b_only":
        if a_present or not b_present or ratio is not None or log_ratio is not None:
            return None, "ab_comparison_availability_invalid"
    elif a_present or b_present or ratio is not None or log_ratio is not None:
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
            "status": method_a["comparison_candidate_status"],
        },
        "method_b": {
            "present": b_present,
            "candidate": b_candidate if b_present else None,
            "uncertainty": b_uncertainty if b_present else None,
            "status": method_b["comparison_candidate_status"],
        },
        "comparison": {
            "availability": availability,
            "ratio_B_over_A": _d11_finite(ratio) if ratio is not None else None,
            "log_ratio_B_over_A": _d11_finite(log_ratio) if log_ratio is not None else None,
        },
    }, None


def _d11_page_cell(cell):
    """Return an independent scalar copy for one D.11 page group."""
    return {
        "t_index": cell["t_index"], "t_low": cell["t_low"], "t_high": cell["t_high"],
        "delta_index": cell["delta_index"], "delta_low": cell["delta_low"],
        "delta_high": cell["delta_high"],
        "method_a": dict(cell["method_a"]),
        "method_b": dict(cell["method_b"]),
        "comparison": dict(cell["comparison"]),
    }


def build_full_background_subtraction_d11_payload(phase_d_checkpoint):
    """Select only frozen Phase-D/D.4 display scalars without recomputation."""
    checkpoint, reason = _d11_checkpoint_contract(phase_d_checkpoint)
    if reason is not None:
        return _d11_unavailable(reason)
    method_a, reason = _d11_representation_contract(
        checkpoint["method_a_comparison"], _D11_METHOD_A_COMPARISON_SCHEMA,
        "method_b_numerical_dependency", "source_method_a_payload_fingerprint",
        "method_a_comparison_contract_invalid",
    )
    if reason is not None:
        return _d11_unavailable(reason)
    method_b, reason = _d11_representation_contract(
        checkpoint["method_b_comparison"], _D11_METHOD_B_COMPARISON_SCHEMA,
        "method_a_numerical_dependency", "source_method_b_payload_fingerprint",
        "method_b_comparison_contract_invalid",
    )
    if reason is not None:
        return _d11_unavailable(reason)
    comparison, reason = _d11_ab_contract(checkpoint["ab_comparison"])
    if reason is not None:
        return _d11_unavailable(reason)
    reason = _d11_source_linkage(checkpoint, method_a, method_b, comparison)
    if reason is not None:
        return _d11_unavailable(reason)
    t_edges = list(comparison["canonical_t_edges"])
    delta_edges = list(comparison["delta_edges"])
    expected = (len(t_edges) - 1) * (len(delta_edges) - 1)
    if len(comparison["cells"]) != expected:
        return _d11_unavailable("ab_comparison_cell_grid_invalid")
    grouped = [[] for _unused in range(len(t_edges) - 1)]
    seen = set()
    for source in comparison["cells"]:
        cell, reason = _d11_scalar_cell(source, t_edges, delta_edges)
        if reason is not None:
            return _d11_unavailable(reason)
        coordinate = (cell["t_index"], cell["delta_index"])
        if coordinate in seen:
            return _d11_unavailable("ab_comparison_cell_grid_invalid")
        seen.add(coordinate)
        grouped[cell["t_index"]].append(cell)
    if len(seen) != expected:
        return _d11_unavailable("ab_comparison_cell_grid_invalid")
    per_t = []
    for t_index, cells in enumerate(grouped):
        cells.sort(key=lambda row: row["delta_index"])
        if len(cells) != len(delta_edges) - 1:
            return _d11_unavailable("ab_comparison_cell_grid_invalid")
        overlay_cells = [_d11_page_cell(cell) for cell in cells if (
            cell["method_a"]["present"] or cell["method_b"]["present"]
        )]
        ratio_cells = [_d11_page_cell(cell) for cell in cells if (
            cell["comparison"]["availability"] == "both_comparable"
        )]
        availability_cells = []
        for cell in cells:
            copied = _d11_page_cell(cell)
            copied["availability_label"] = _D11_AVAILABILITY_LABELS[
                copied["comparison"]["availability"]
            ]
            availability_cells.append(copied)
        per_t.append({
            "t_index": t_index,
            "t_low": float(t_edges[t_index]),
            "t_high": float(t_edges[t_index + 1]),
            "ab_overlay": {
                "available": bool(overlay_cells),
                "reason": None if overlay_cells else "no_stored_method_candidates",
                "cells": overlay_cells,
            },
            "ratio_log": {
                "available": bool(ratio_cells),
                "reason": None if ratio_cells else "no_stored_both_comparable_cells",
                "cells": ratio_cells,
            },
            "availability": {"available": True, "reason": None, "cells": availability_cells},
        })
    return {
        "schema_version": D11_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "source_checkpoint_payload_fingerprint": checkpoint["source_checkpoint_payload_fingerprint"],
        "phase_a_contract_fingerprint": comparison["phase_a_contract_fingerprint"],
        "coordinate_fingerprint": comparison["coordinate_fingerprint"],
        "method_a_comparison_fingerprint": comparison["method_a_comparison_fingerprint"],
        "method_b_comparison_fingerprint": comparison["method_b_comparison_fingerprint"],
        "ab_comparison_fingerprint": comparison["fingerprint"],
        "host_state": comparison["host_state"],
        "source_target_state": comparison["source_target_state"],
        "t_edges": [float(edge) for edge in t_edges],
        "delta_edges": [float(edge) for edge in delta_edges],
        "per_t": per_t,
    }


def _d11_display_bounds(values, reference):
    finite = [float(reference)]
    finite.extend(float(value) for value in values if _d11_finite(value) is not None)
    low, high = min(finite), max(finite)
    if high <= low:
        return low - 0.05, high + 0.05
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def _d11_new_frame(ROOT, name, title, delta_edges, y_range):
    if not hasattr(ROOT, "TH1D"):
        return None
    try:
        frame = ROOT.TH1D(str(name), str(title), len(delta_edges) - 1, array("d", delta_edges))
        frame.SetDirectory(0)
        if hasattr(frame, "SetStats"):
            frame.SetStats(0)
        frame.SetMinimum(float(y_range[0]))
        frame.SetMaximum(float(y_range[1]))
        return frame
    except Exception:
        return None


def _render_d11_ab_overlay_page(ROOT, pdf_name, presentation, group):
    overlay = _mapping(group.get("ab_overlay"))
    cells = tuple(overlay.get("cells") or ())
    if not cells or not hasattr(ROOT, "TGraphAsymmErrors") or not hasattr(ROOT, "TGraphErrors"):
        return False
    values = [1.0]
    for cell in cells:
        cell = _mapping(cell)
        method_a = _mapping(cell.get("method_a"))
        method_b = _mapping(cell.get("method_b"))
        if method_a.get("present") is True:
            values.extend((method_a["low"], method_a["candidate"], method_a["high"]))
        if method_b.get("present") is True:
            values.extend((
                method_b["candidate"] - method_b["uncertainty"], method_b["candidate"],
                method_b["candidate"] + method_b["uncertainty"],
            ))
    y_range = _d11_display_bounds(values, 1.0)
    title = "Method A / Method B diagnostic comparison - {}".format(_t_context(group))
    frame = _d11_new_frame(
        ROOT, "H_full_background_d11_overlay_t{}".format(group["t_index"] + 1),
        "{};delta [%];Relative pion-background diagnostic".format(title),
        presentation["delta_edges"], y_range,
    )
    if frame is None:
        return False
    canvas = ROOT.TCanvas(
        "C_full_background_d11_overlay_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    draw_objects = [frame]
    try:
        frame.Draw("axis")
        a_cells = [cell for cell in cells if _mapping(cell).get("method_a", {}).get("present") is True]
        b_cells = [cell for cell in cells if _mapping(cell).get("method_b", {}).get("present") is True]
        legend = ROOT.TLegend(0.53, 0.61, 0.89, 0.79)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        if a_cells:
            graph_a = ROOT.TGraphAsymmErrors(len(a_cells))
            for index, cell in enumerate(a_cells):
                method_a = cell["method_a"]
                x_value = 0.5 * (cell["delta_low"] + cell["delta_high"])
                graph_a.SetPoint(index, x_value, method_a["candidate"])
                graph_a.SetPointError(
                    index, 0.0, 0.0, method_a["candidate"] - method_a["low"],
                    method_a["high"] - method_a["candidate"],
                )
            graph_a.SetMarkerColor(getattr(ROOT, "kBlue", 4))
            graph_a.SetLineColor(getattr(ROOT, "kBlue", 4))
            graph_a.SetMarkerStyle(20)
            graph_a.Draw("P same")
            legend.AddEntry(graph_a, "Method A - HGCer response diagnostic", "lep")
            draw_objects.append(graph_a)
        if b_cells:
            graph_b = ROOT.TGraphErrors(len(b_cells))
            for index, cell in enumerate(b_cells):
                method_b = cell["method_b"]
                x_value = 0.5 * (cell["delta_low"] + cell["delta_high"])
                graph_b.SetPoint(index, x_value, method_b["candidate"])
                graph_b.SetPointError(index, 0.0, method_b["uncertainty"])
            graph_b.SetMarkerColor(getattr(ROOT, "kRed", 2))
            graph_b.SetLineColor(getattr(ROOT, "kRed", 2))
            graph_b.SetMarkerStyle(24)
            graph_b.Draw("P same")
            legend.AddEntry(graph_b, "Method B - Missing-mass closure diagnostic", "lep")
            draw_objects.append(graph_b)
        legend.Draw()
        draw_objects.append(legend)
        unity = ROOT.TLine(
            float(presentation["delta_edges"][0]), 1.0,
            float(presentation["delta_edges"][-1]), 1.0,
        )
        unity.SetLineStyle(2)
        unity.SetLineWidth(2)
        unity.Draw()
        draw_objects.append(unity)
        draw_objects.append(_draw_page_header(ROOT, canvas, "Method A / Method B diagnostic comparison", group))
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d11_ratio_log_page(ROOT, pdf_name, presentation, group):
    ratio_log = _mapping(group.get("ratio_log"))
    cells = tuple(ratio_log.get("cells") or ())
    if not cells or not hasattr(ROOT, "TGraph"):
        return False
    title = "Method A / Method B relative comparison - {}".format(_t_context(group))
    canvas = ROOT.TCanvas(
        "C_full_background_d11_ratio_log_t{}".format(group["t_index"] + 1), title, 1400, 800
    )
    draw_objects = []
    try:
        canvas.Divide(2, 1)
        for panel, key, label, reference in (
            (1, "ratio_B_over_A", "B/A", 1.0),
            (2, "log_ratio_B_over_A", "ln(B/A)", 0.0),
        ):
            canvas.cd(panel)
            frame = _d11_new_frame(
                ROOT,
                "H_full_background_d11_{}_t{}".format(key, group["t_index"] + 1),
                "{};delta [%];{}".format(title, label), presentation["delta_edges"],
                _d11_display_bounds((cell["comparison"][key] for cell in cells), reference),
            )
            if frame is None:
                return False
            frame.Draw("axis")
            graph = ROOT.TGraph(len(cells))
            for index, cell in enumerate(cells):
                graph.SetPoint(
                    index, 0.5 * (cell["delta_low"] + cell["delta_high"]),
                    cell["comparison"][key],
                )
            graph.SetMarkerColor(getattr(ROOT, "kBlack", 1))
            graph.SetMarkerStyle(20)
            graph.Draw("P same")
            line = ROOT.TLine(
                float(presentation["delta_edges"][0]), float(reference),
                float(presentation["delta_edges"][-1]), float(reference),
            )
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            draw_objects.extend((frame, graph, line))
        canvas.cd(1)
        draw_objects.append(_draw_page_header(ROOT, canvas, "Method A / Method B relative comparison", group))
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _d11_availability_style(ROOT, availability):
    colors = {
        "both_comparable": getattr(ROOT, "kGreen", 3),
        "both_present_not_comparable": getattr(ROOT, "kOrange", 800),
        "a_only": getattr(ROOT, "kBlue", 4),
        "b_only": getattr(ROOT, "kRed", 2),
        "neither_available": getattr(ROOT, "kGray", 920),
    }
    return colors[availability]


def _render_d11_availability_page(ROOT, pdf_name, presentation, group):
    availability = _mapping(group.get("availability"))
    cells = tuple(availability.get("cells") or ())
    if not cells or not hasattr(ROOT, "TBox"):
        return False
    title = "Method availability across delta - {}".format(_t_context(group))
    frame = _d11_new_frame(
        ROOT, "H_full_background_d11_availability_t{}".format(group["t_index"] + 1),
        "{};delta [%];".format(title), presentation["delta_edges"], (0.0, 1.0),
    )
    if frame is None:
        return False
    canvas = ROOT.TCanvas(
        "C_full_background_d11_availability_t{}".format(group["t_index"] + 1), title, 1200, 800
    )
    draw_objects = [frame]
    try:
        frame.Draw("axis")
        for cell in cells:
            state = cell["comparison"]["availability"]
            box = ROOT.TBox(cell["delta_low"], 0.0, cell["delta_high"], 1.0)
            box.SetFillColor(_d11_availability_style(ROOT, state))
            box.SetLineColor(getattr(ROOT, "kBlack", 1))
            box.Draw()
            draw_objects.append(box)
        legend = ROOT.TLegend(0.55, 0.44, 0.89, 0.76)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        for state, label in _D11_AVAILABILITY_LABELS.items():
            exemplar = ROOT.TBox(0.0, 0.0, 0.0, 0.0)
            exemplar.SetFillColor(_d11_availability_style(ROOT, state))
            exemplar.SetLineColor(getattr(ROOT, "kBlack", 1))
            legend.AddEntry(exemplar, label, "f")
            draw_objects.append(exemplar)
        legend.Draw()
        draw_objects.append(legend)
        draw_objects.append(_draw_page_header(ROOT, canvas, "Method availability across delta", group))
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_d11_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Render the frozen D.11 explanation for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    overlay = _mapping(group.get("ab_overlay"))
    if overlay.get("available"):
        if _render_d11_ab_overlay_page(ROOT, pdf_name, presentation, group):
            manifest.append({"page_id": "full_background.d11.ab_overlay", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.11 A/B-overlay page unavailable for t{}".format(t_number))
    else:
        failures.append("D.11 A/B-overlay unavailable for t{}: {}".format(t_number, overlay.get("reason")))
    ratio_log = _mapping(group.get("ratio_log"))
    if ratio_log.get("available"):
        if _render_d11_ratio_log_page(ROOT, pdf_name, presentation, group):
            manifest.append({"page_id": "full_background.d11.ab_ratio_log", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.11 B/A-and-log page unavailable for t{}".format(t_number))
    else:
        failures.append("D.11 B/A-and-log unavailable for t{}: {}".format(t_number, ratio_log.get("reason")))
    availability = _mapping(group.get("availability"))
    if availability.get("available"):
        if _render_d11_availability_page(ROOT, pdf_name, presentation, group):
            manifest.append({"page_id": "full_background.d11.method_availability", "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("D.11 availability page unavailable for t{}".format(t_number))
    else:
        failures.append("D.11 availability unavailable for t{}: {}".format(t_number, availability.get("reason")))


def render_full_background_subtraction_d11_pages(pdf_name, payload, *, page_manifest=None):
    """Append D.11 stored A/B explanation pages for every canonical t bin."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    presentation = _mapping(payload)
    if not bool(presentation.get("available")):
        result["failures"].append(
            "D.11 procedure input unavailable: {}".format(presentation.get("reason"))
        )
        return result
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("D.11 procedure rendering unavailable: PyROOT not available")
        return result
    for group in tuple(presentation.get("per_t") or ()):
        _render_d11_t_pages(
            ROOT, pdf_name, presentation, _mapping(group), manifest, result["failures"]
        )
    return result


def render_full_background_subtraction_procedure_pages(
    pdf_name, d6_payload, d7_payload, d8_payload=None, d9_payload=None, d10_payload=None,
    d11_payload=None,
    *, e2_payload=None, page_manifest=None
):
    """Append available D.6 through D.11 and E.2 groups in canonical-t order."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    d6 = _mapping(d6_payload)
    d7 = _mapping(d7_payload)
    d8 = _mapping(d8_payload)
    d9 = _mapping(d9_payload)
    d10 = _mapping(d10_payload)
    d11 = _mapping(d11_payload)
    e2 = _mapping(e2_payload)
    d6_available = bool(d6.get("available"))
    d7_available = bool(d7.get("available"))
    d8_requested = d8_payload is not None
    d8_available = bool(d8.get("available"))
    d9_requested = d9_payload is not None
    d9_available = bool(d9.get("available"))
    d10_requested = d10_payload is not None
    d10_available = bool(d10.get("available"))
    d11_requested = d11_payload is not None
    d11_available = bool(d11.get("available"))
    e2_requested = e2_payload is not None
    e2_available = bool(e2.get("available"))
    if not d6_available and not d7_available and not d8_available and not d9_available and not d10_available and not d11_available and not e2_available:
        if d6_payload is not None:
            result["failures"].append(
                "D.6 procedure input unavailable: {}".format(d6.get("reason"))
            )
        if d7_payload is not None:
            result["failures"].append(
                "D.7 procedure input unavailable: {}".format(d7.get("reason"))
            )
        if d8_requested:
            result["failures"].append(
                "D.8 procedure input unavailable: {}".format(d8.get("reason"))
            )
        if d9_requested:
            result["failures"].append(
                "D.9 procedure input unavailable: {}".format(d9.get("reason"))
            )
        if d10_requested:
            result["failures"].append(
                "D.10 procedure input unavailable: {}".format(d10.get("reason"))
            )
        if d11_requested:
            result["failures"].append(
                "D.11 procedure input unavailable: {}".format(d11.get("reason"))
            )
        if e2_requested:
            result["failures"].append(
                "E.2 procedure input unavailable: {}".format(e2.get("reason"))
            )
        return result
    if not d6_available:
        result["failures"].append(
            "D.6 procedure input unavailable: {}".format(d6.get("reason"))
        )
    if not d7_available:
        result["failures"].append(
            "D.7 procedure input unavailable: {}".format(d7.get("reason"))
        )
    if d8_requested and not d8_available:
        result["failures"].append(
            "D.8 procedure input unavailable: {}".format(d8.get("reason"))
        )
    if d9_requested and not d9_available:
        result["failures"].append(
            "D.9 procedure input unavailable: {}".format(d9.get("reason"))
        )
    if d10_requested and not d10_available:
        result["failures"].append(
            "D.10 procedure input unavailable: {}".format(d10.get("reason"))
        )
    if d11_requested and not d11_available:
        result["failures"].append(
            "D.11 procedure input unavailable: {}".format(d11.get("reason"))
        )
    if e2_requested and not e2_available:
        result["failures"].append(
            "E.2 procedure input unavailable: {}".format(e2.get("reason"))
        )
    if d6_available and d7_available and list(d6.get("t_edges") or ()) != list(d7.get("t_edges") or ()):
        result["failures"].append("D.6/D.7 canonical t geometry mismatch")
        d7_available = False
    geometry_owner = (
        d6 if d6_available else d7 if d7_available else d8 if d8_available
        else d9 if d9_available else d10 if d10_available else d11 if d11_available
        else None
    )
    if (
        d8_available
        and geometry_owner is not None
        and list(geometry_owner.get("t_edges") or ()) != list(d8.get("t_edges") or ())
    ):
        result["failures"].append("D.8 canonical t geometry mismatch")
        d8_available = False
    if (
        d9_available
        and geometry_owner is not d9
        and list(geometry_owner.get("t_edges") or ()) != list(d9.get("t_edges") or ())
    ):
        result["failures"].append("D.9 canonical t geometry mismatch")
        d9_available = False
    if (
        d10_available
        and geometry_owner is not d10
        and list(geometry_owner.get("t_edges") or ()) != list(d10.get("t_edges") or ())
    ):
        result["failures"].append("D.10 canonical t geometry mismatch")
        d10_available = False
    if (
        d11_available
        and geometry_owner is not d11
        and list(geometry_owner.get("t_edges") or ()) != list(d11.get("t_edges") or ())
    ):
        result["failures"].append("D.11 canonical t geometry mismatch")
        d11_available = False
    for comparison_name, comparison_payload, comparison_available in (
        ("D.6", d6, d6_available),
        ("D.7", d7, d7_available),
    ):
        if (
            d8_available
            and comparison_available
            and list(comparison_payload.get("delta_edges") or ())
            and list(comparison_payload.get("delta_edges") or ()) != list(d8.get("delta_edges") or ())
        ):
            result["failures"].append(
                "D.8 cache delta geometry differs from {}".format(comparison_name)
            )
    for comparison_name, comparison_payload, comparison_available in (
        ("D.6", d6, d6_available),
        ("D.7", d7, d7_available),
        ("D.8", d8, d8_available),
    ):
        if (
            d9_available
            and comparison_available
            and list(comparison_payload.get("delta_edges") or ())
            and list(comparison_payload.get("delta_edges") or ()) != list(d9.get("delta_edges") or ())
        ):
            result["failures"].append(
                "D.9 frozen delta geometry differs from {}".format(comparison_name)
            )
    for comparison_name, comparison_payload, comparison_available in (
        ("D.6", d6, d6_available),
        ("D.7", d7, d7_available),
        ("D.8", d8, d8_available),
        ("D.9", d9, d9_available),
    ):
        if (
            d10_available
            and comparison_available
            and list(comparison_payload.get("delta_edges") or ())
            and list(comparison_payload.get("delta_edges") or ()) != list(d10.get("delta_edges") or ())
        ):
            result["failures"].append(
                "D.10 frozen delta geometry differs from {}".format(comparison_name)
            )
    for comparison_name, comparison_payload, comparison_available, hard_gate in (
        ("D.6", d6, d6_available, False),
        ("D.7", d7, d7_available, False),
        ("D.8", d8, d8_available, False),
        ("D.9", d9, d9_available, True),
        ("D.10", d10, d10_available, True),
    ):
        if (
            d11_available
            and comparison_available
            and list(comparison_payload.get("delta_edges") or ())
            and list(comparison_payload.get("delta_edges") or ()) != list(d11.get("delta_edges") or ())
        ):
            result["failures"].append(
                "D.11 frozen delta geometry differs from {}".format(comparison_name)
            )
            if hard_gate:
                d11_available = False
    if e2_available and geometry_owner is not None and (
        list(geometry_owner.get("t_edges") or ()) != list(e2.get("t_edges") or ())
        or list(geometry_owner.get("delta_edges") or ()) != list(e2.get("delta_edges") or ())
    ):
        result["failures"].append("E.2 frozen procedure geometry mismatch")
        e2_available = False
    ROOT = _import_root()
    if ROOT is None:
        result["failures"].append("full background-subtraction rendering unavailable: PyROOT not available")
        return result
    if d7_available:
        _append_d7_exclusion_failure(d7, result["failures"])
    d7_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(d7.get("per_t") or ())
        if isinstance(group, Mapping)
    }
    d8_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(d8.get("per_t") or ())
        if isinstance(group, Mapping)
    }
    d9_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(d9.get("per_t") or ())
        if isinstance(group, Mapping)
    }
    d10_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(d10.get("per_t") or ())
        if isinstance(group, Mapping)
    }
    d11_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(d11.get("per_t") or ())
        if isinstance(group, Mapping)
    }
    e2_by_index = {
        group.get("t_index"): _mapping(group)
        for group in tuple(e2.get("per_t") or ())
        if isinstance(group, Mapping)
    }

    def render_d8_group(t_index):
        if not d8_available:
            return
        d8_group = d8_by_index.get(t_index)
        if d8_group is None:
            result["failures"].append(
                "D.8 input missing canonical t{}".format(int(t_index) + 1)
            )
            return
        _render_d8_t_pages(ROOT, pdf_name, d8, d8_group, manifest, result["failures"])

    def render_d9_group(t_index):
        if not d9_available:
            return
        d9_group = d9_by_index.get(t_index)
        if d9_group is None:
            result["failures"].append(
                "D.9 input missing canonical t{}".format(int(t_index) + 1)
            )
            return
        _render_d9_t_pages(ROOT, pdf_name, d9, d9_group, manifest, result["failures"])

    def render_d10_group(t_index):
        if not d10_available:
            return
        d10_group = d10_by_index.get(t_index)
        if d10_group is None:
            result["failures"].append(
                "D.10 input missing canonical t{}".format(int(t_index) + 1)
            )
            return
        _render_d10_t_pages(ROOT, pdf_name, d10, d10_group, manifest, result["failures"])

    def render_d11_group(t_index):
        if not d11_available:
            return
        d11_group = d11_by_index.get(t_index)
        if d11_group is None:
            result["failures"].append(
                "D.11 input missing canonical t{}".format(int(t_index) + 1)
            )
            return
        _render_d11_t_pages(ROOT, pdf_name, d11, d11_group, manifest, result["failures"])

    def render_e2_group(t_index):
        if not e2_available:
            return
        e2_group = e2_by_index.get(t_index)
        if e2_group is None:
            result["failures"].append(
                "E.2 input missing canonical t{}".format(int(t_index) + 1)
            )
            return
        _render_e2_t_pages(ROOT, pdf_name, e2, e2_group, manifest, result["failures"])

    if d6_available:
        for group in tuple(d6.get("per_t") or ()):
            group = _mapping(group)
            _render_d6_t_pages(ROOT, pdf_name, d6, group, manifest, result["failures"])
            if d7_available:
                d7_group = d7_by_index.get(group.get("t_index"))
                if d7_group is None:
                    result["failures"].append(
                        "D.7 input missing canonical t{}".format(int(group.get("t_index", -1)) + 1)
                    )
                else:
                    _render_d7_t_pages(ROOT, pdf_name, d7, d7_group, manifest, result["failures"])
            render_d8_group(group.get("t_index"))
            render_d9_group(group.get("t_index"))
            render_d10_group(group.get("t_index"))
            render_d11_group(group.get("t_index"))
            render_e2_group(group.get("t_index"))
    elif d7_available:
        for group in tuple(d7.get("per_t") or ()):
            group = _mapping(group)
            _render_d7_t_pages(ROOT, pdf_name, d7, group, manifest, result["failures"])
            render_d8_group(group.get("t_index"))
            render_d9_group(group.get("t_index"))
            render_d10_group(group.get("t_index"))
            render_d11_group(group.get("t_index"))
            render_e2_group(group.get("t_index"))
    elif d8_available:
        for group in tuple(d8.get("per_t") or ()):
            group = _mapping(group)
            _render_d8_t_pages(ROOT, pdf_name, d8, group, manifest, result["failures"])
            render_d9_group(group.get("t_index"))
            render_d10_group(group.get("t_index"))
            render_d11_group(group.get("t_index"))
            render_e2_group(group.get("t_index"))
    elif d9_available:
        for group in tuple(d9.get("per_t") or ()):
            _render_d9_t_pages(ROOT, pdf_name, d9, _mapping(group), manifest, result["failures"])
            render_d10_group(_mapping(group).get("t_index"))
            render_d11_group(_mapping(group).get("t_index"))
            render_e2_group(_mapping(group).get("t_index"))
    elif d10_available:
        for group in tuple(d10.get("per_t") or ()):
            _render_d10_t_pages(ROOT, pdf_name, d10, _mapping(group), manifest, result["failures"])
            render_d11_group(_mapping(group).get("t_index"))
            render_e2_group(_mapping(group).get("t_index"))
    elif d11_available:
        for group in tuple(d11.get("per_t") or ()):
            _render_d11_t_pages(ROOT, pdf_name, d11, _mapping(group), manifest, result["failures"])
            render_e2_group(_mapping(group).get("t_index"))
    elif e2_available:
        for group in tuple(e2.get("per_t") or ()):
            _render_e2_t_pages(ROOT, pdf_name, e2, _mapping(group), manifest, result["failures"])
    return result


def _e2_unavailable(reason):
    """Return an E.2-local unavailable payload without source aliases."""
    return {
        "schema_version": E2_PRESENTATION_SCHEMA_VERSION,
        "available": False,
        "reason": str(reason),
        "non_authoritative": True,
        "production_objects_mutated": False,
        "phase_a_contract_fingerprint": None,
        "coordinate_fingerprint": None,
        "host_state": None,
        "host_label": None,
        "source_target_state": None,
        "physical_pion_control_threshold": None,
        "t_edges": [],
        "delta_edges": [],
        "per_t": [],
    }


def _e2_finite(value):
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    return scalar if math.isfinite(scalar) else None


def _e2_integer(value):
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return int(value)


def _e2_optional_coordinate(value):
    if value is None:
        return None, True
    scalar = _e2_finite(value)
    return scalar, scalar is not None


def _e2_sidecar_row(row, side, phase):
    """Validate one frozen acceptance row without assigning a new cell."""
    row = _mapping(row)
    required = (
        "side", "source_label", "entry_index", "coordinate_fingerprint",
        "analysis_t", "canonical_t_index", "ssdelta", "delta_index",
        "P_hgcer_npeSum", "ssxptar", "ssyptar", "diagnostic_weight",
        "allcuts", "nommcuts", "rf_applied_to_diagnostic",
    )
    if not row or any(key not in row for key in required):
        return None, "e2_sidecar_row_invalid"
    t_index = _e2_integer(row["canonical_t_index"])
    delta_index = _e2_integer(row["delta_index"])
    entry_index = _e2_integer(row["entry_index"])
    analysis_t = _e2_finite(row["analysis_t"])
    delta = _e2_finite(row["ssdelta"])
    npe = _e2_finite(row["P_hgcer_npeSum"])
    weight = _e2_finite(row["diagnostic_weight"])
    ssxptar, x_valid = _e2_optional_coordinate(row["ssxptar"])
    ssyptar, y_valid = _e2_optional_coordinate(row["ssyptar"])
    if (
        row["side"] != side
        or not _d10_nonempty_string(row["source_label"])
        or entry_index is None
        or row["coordinate_fingerprint"] != phase["coordinate_fingerprint"]
        or analysis_t is None
        or delta is None
        or npe is None
        or weight is None
        or not x_valid
        or not y_valid
        or type(row["allcuts"]) is not bool
        or type(row["nommcuts"]) is not bool
        or row["rf_applied_to_diagnostic"] is not False
        or t_index is None
        or delta_index is None
        or not 0 <= t_index < len(phase["t_edges"]) - 1
        or not 0 <= delta_index < len(phase["delta_edges"]) - 1
    ):
        return None, "e2_sidecar_row_invalid"
    if (
        find_canonical_bin(analysis_t, tuple(phase["t_edges"])) != t_index
        or find_canonical_bin(delta, tuple(phase["delta_edges"])) != delta_index
    ):
        return None, "e2_sidecar_cell_assignment_mismatch"
    return {
        "t_index": t_index,
        "delta_index": delta_index,
        "npe": npe,
        "ssxptar": ssxptar,
        "ssyptar": ssyptar,
        "absolute_support_weight": abs(weight),
        "nommcuts": row["nommcuts"],
    }, None


def build_full_background_subtraction_e2_payload(
    pion_hgcer_tdelta_diagnostic, pion_hgcer_event_contract,
):
    """Build detached SHMS-angle support shapes from frozen Part-1 rows."""
    phase, reason = _d10_phase_a_contract(pion_hgcer_event_contract)
    if reason is not None:
        return _e2_unavailable("e2_phase_a_contract_invalid")
    diagnostic = _mapping(pion_hgcer_tdelta_diagnostic)
    sidecar = _mapping(diagnostic.get("phase_e_acceptance_records"))
    t_edges = _strict_edges(diagnostic.get("t_edges"))
    delta_edges = _strict_edges(diagnostic.get("delta_edges"))
    configuration = _mapping(diagnostic.get("config"))
    threshold = _e2_finite(configuration.get("production_hgcer_threshold"))
    if (
        not diagnostic
        or diagnostic.get("status") != "available"
        or diagnostic.get("non_authoritative") is not True
        or diagnostic.get("production_side_effect_free") is not True
        or diagnostic.get("production_hgcer_pid_unchanged") is not True
        or diagnostic.get("rf_restoration_applied") is not False
        or not _d10_nonempty_string(diagnostic.get("coordinate_fingerprint"))
        or t_edges is None
        or delta_edges is None
        or list(t_edges) != phase["t_edges"]
        or list(delta_edges) != phase["delta_edges"]
        or diagnostic.get("coordinate_fingerprint") != phase["coordinate_fingerprint"]
        or threshold is None
        or abs(threshold - 2.0) > 1.0e-12
    ):
        return _e2_unavailable("e2_part1_contract_invalid")
    if set(sidecar) != {"kaon", "pion"}:
        return _e2_unavailable("e2_acceptance_sidecar_invalid")
    if any(
        not _d10_sequence(sidecar[side])
        for side in ("kaon", "pion")
    ):
        return _e2_unavailable("e2_acceptance_sidecar_invalid")

    cells_by_index = {
        (t_index, delta_index): {
            "delta_index": delta_index,
            "delta_low": float(delta_edges[delta_index]),
            "delta_high": float(delta_edges[delta_index + 1]),
            "kaon_rows": [],
            "pion_rows": [],
        }
        for t_index in range(len(t_edges) - 1)
        for delta_index in range(len(delta_edges) - 1)
    }
    for side in ("kaon", "pion"):
        for source_row in sidecar[side]:
            row, row_reason = _e2_sidecar_row(source_row, side, phase)
            if row_reason is not None:
                return _e2_unavailable(row_reason)
            if not row["nommcuts"]:
                continue
            if side == "pion" and not row["npe"] > threshold:
                continue
            cells_by_index[(row["t_index"], row["delta_index"])][
                "{}_rows".format(side)
            ].append({
                "ssxptar": row["ssxptar"],
                "ssyptar": row["ssyptar"],
                "absolute_support_weight": row["absolute_support_weight"],
            })

    per_t = []
    for t_index in range(len(t_edges) - 1):
        cells = []
        for delta_index in range(len(delta_edges) - 1):
            source = cells_by_index[(t_index, delta_index)]
            cells.append({
                "delta_index": delta_index,
                "delta_low": source["delta_low"],
                "delta_high": source["delta_high"],
                "kaon_rows": tuple(dict(row) for row in source["kaon_rows"]),
                "pion_rows": tuple(dict(row) for row in source["pion_rows"]),
            })
        per_t.append({
            "t_index": t_index,
            "t_low": float(t_edges[t_index]),
            "t_high": float(t_edges[t_index + 1]),
            "cells": tuple(cells),
        })
    host_state = phase["host_state"]
    return {
        "schema_version": E2_PRESENTATION_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "phase_a_contract_fingerprint": phase["contract_fingerprint"],
        "coordinate_fingerprint": phase["coordinate_fingerprint"],
        "host_state": host_state,
        "host_label": (
            "Proton-cleaned kaon sample"
            if host_state == "proton_cleaned" else "Kaon-selected sample"
        ),
        "source_target_state": phase["source_target_state"],
        "physical_pion_control_threshold": 2.0,
        "selection_definition": {
            "common": "noRF_nommcuts_canonical_t_delta",
            "kaon": "kaon_side_noRF_nommcuts",
            "pion": "pion_control_noRF_nommcuts_HGCer_gt_2",
            "weighting": "absolute_diagnostic_support",
            "normalization": "unit_area_independently_per_side_per_cell",
        },
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "per_t": tuple(per_t),
    }


def _e2_display_bounds(coordinate):
    """Return the fixed, presentation-only E.2 SHMS coordinate bounds."""
    if coordinate == "ssxptar":
        return -0.10, 0.10
    return -0.06, 0.06


def _e2_histogram(ROOT, name, title, coordinate, rows):
    """Fill and independently normalize one detached E.2 display histogram."""
    lower, upper = _e2_display_bounds(coordinate)
    try:
        histogram = ROOT.TH1D(
            str(name), str(title), 50, lower, upper,
        )
        histogram.SetDirectory(0)
        if hasattr(histogram, "SetStats"):
            histogram.SetStats(0)
    except Exception:
        return None
    for row in tuple(rows or ()):
        row = _mapping(row)
        value = _e2_finite(row.get(coordinate))
        weight = _e2_finite(row.get("absolute_support_weight"))
        if value is not None and weight is not None and weight > 0.0:
            histogram.Fill(value, weight)
    try:
        integral = float(histogram.Integral())
    except Exception:
        return None
    if not math.isfinite(integral) or integral <= 0.0:
        return None
    histogram.Scale(1.0 / integral)
    return histogram


def _e2_empty_frame(ROOT, name, title, coordinate):
    """Create an unfilled, detached frame for an empty E.2 canonical cell."""
    lower, upper = _e2_display_bounds(coordinate)
    try:
        frame = ROOT.TH1D(str(name), str(title), 50, lower, upper)
        frame.SetDirectory(0)
        if hasattr(frame, "SetStats"):
            frame.SetStats(0)
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.0)
        return frame
    except Exception:
        return None


def _e2_panel_notice(ROOT, text):
    try:
        note = ROOT.TPaveText(0.12, 0.38, 0.88, 0.62, "NDC")
        note.SetFillStyle(0)
        note.SetBorderSize(0)
        note.SetTextAlign(22)
        note.SetTextSize(0.045)
        note.AddText(str(text))
        note.Draw()
        return note
    except Exception:
        return None


def _render_e2_acceptance_page(ROOT, pdf_name, presentation, group, coordinate):
    cells = tuple(group.get("cells") or ())
    if not cells or not hasattr(ROOT, "TH1D"):
        return False
    x_label = "SHMS x'_{tar}" if coordinate == "ssxptar" else "SHMS y'_{tar}"
    heading = "{} acceptance across delta".format(x_label)
    title = "{} - {}".format(heading, _t_context(group))
    columns = min(3, len(cells))
    rows = int(math.ceil(float(len(cells)) / float(columns)))
    canvas = ROOT.TCanvas(
        "C_full_background_e2_{}_t{}".format(coordinate, group["t_index"] + 1),
        title,
        1400,
        900,
    )
    canvas.Divide(columns, rows)
    draw_objects = []
    legend_drawn = False
    try:
        for panel_index, cell in enumerate(cells, 1):
            cell = _mapping(cell)
            canvas.cd(panel_index)
            panel_title = "delta = [{:.3f}, {:.3f}] %;{};Normalized absolute support".format(
                float(cell["delta_low"]), float(cell["delta_high"]), x_label,
            )
            kaon = _e2_histogram(
                ROOT,
                "H_full_background_e2_{}_kaon_t{}_d{}".format(
                    coordinate, group["t_index"] + 1, int(cell["delta_index"]) + 1,
                ),
                panel_title,
                coordinate,
                cell.get("kaon_rows"),
            )
            pion = _e2_histogram(
                ROOT,
                "H_full_background_e2_{}_pion_t{}_d{}".format(
                    coordinate, group["t_index"] + 1, int(cell["delta_index"]) + 1,
                ),
                panel_title,
                coordinate,
                cell.get("pion_rows"),
            )
            if kaon is None and pion is None:
                frame = _e2_empty_frame(
                    ROOT,
                    "H_full_background_e2_{}_frame_t{}_d{}".format(
                        coordinate, group["t_index"] + 1, int(cell["delta_index"]) + 1,
                    ),
                    panel_title,
                    coordinate,
                )
                if frame is None:
                    return False
                frame.Draw("hist")
                draw_objects.append(frame)
                notice = _e2_panel_notice(ROOT, "No finite acceptance support")
                if notice is not None:
                    draw_objects.append(notice)
                continue
            frame = kaon if kaon is not None else pion
            maximum = max(
                float(histogram.GetMaximum())
                for histogram in (kaon, pion) if histogram is not None
            )
            frame.SetMinimum(0.0)
            frame.SetMaximum(max(1.0e-12, 1.15 * maximum))
            _style_histogram(kaon, getattr(ROOT, "kBlack", 1))
            _style_histogram(pion, getattr(ROOT, "kBlue", 4))
            frame.Draw("hist")
            other = pion if frame is kaon else kaon
            if other is not None:
                other.Draw("hist same")
            if not legend_drawn and hasattr(ROOT, "TLegend"):
                legend = ROOT.TLegend(0.48, 0.66, 0.89, 0.86)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                if kaon is not None:
                    legend.AddEntry(kaon, presentation["host_label"], "l")
                if pion is not None:
                    legend.AddEntry(pion, "Pion-control sample (HGCer NPE > 2)", "l")
                legend.Draw()
                draw_objects.append(legend)
                legend_drawn = True
            if kaon is None or pion is None:
                missing = presentation["host_label"] if kaon is None else "Pion-control sample"
                notice = _e2_panel_notice(ROOT, "{} unavailable".format(missing))
                if notice is not None:
                    draw_objects.append(notice)
            draw_objects.extend(histogram for histogram in (kaon, pion) if histogram is not None)
        canvas.cd(1)
        draw_objects.append(_draw_small_note(
            ROOT, "Shapes normalized independently; noRF, no MM cut; diagnostic only."
        ))
        draw_objects.append(_draw_page_header(ROOT, canvas, heading, group))
        canvas._full_background_e2_draw_objects = tuple(draw_objects)
        canvas.Print(pdf_name)
    finally:
        canvas.Close()
    return True


def _render_e2_t_pages(ROOT, pdf_name, presentation, group, manifest, failures):
    """Append detached E.2 SHMS acceptance pages for one canonical t bin."""
    t_number = int(group.get("t_index", -1)) + 1
    pages = (
        ("ssxptar", "full_background.e2.shms_xptar_acceptance"),
        ("ssyptar", "full_background.e2.shms_yptar_acceptance"),
    )
    for coordinate, page_id in pages:
        if _render_e2_acceptance_page(ROOT, pdf_name, presentation, group, coordinate):
            manifest.append({"page_id": page_id, "scope": "t{}".format(t_number), "authoritative": False})
        else:
            failures.append("E.2 {} page unavailable for t{}".format(coordinate, t_number))


__all__ = (
    "D6_PRESENTATION_SCHEMA_VERSION",
    "D7_PRESENTATION_SCHEMA_VERSION",
    "D8_PRESENTATION_SCHEMA_VERSION",
    "D9_PRESENTATION_SCHEMA_VERSION",
    "D10_PRESENTATION_SCHEMA_VERSION",
    "D11_PRESENTATION_SCHEMA_VERSION",
    "E2_PRESENTATION_SCHEMA_VERSION",
    "FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX",
    "build_full_background_subtraction_d6_payload",
    "build_full_background_subtraction_d7_payload",
    "build_full_background_subtraction_d8_payload",
    "build_full_background_subtraction_d9_payload",
    "build_full_background_subtraction_d10_payload",
    "build_full_background_subtraction_d11_payload",
    "build_full_background_subtraction_e2_payload",
    "close_full_background_subtraction_pdf",
    "full_background_subtraction_pdf_path",
    "open_full_background_subtraction_pdf",
    "render_full_background_subtraction_d6_pages",
    "render_full_background_subtraction_d7_pages",
    "render_full_background_subtraction_d8_pages",
    "render_full_background_subtraction_d9_pages",
    "render_full_background_subtraction_d10_pages",
    "render_full_background_subtraction_d11_pages",
    "render_full_background_subtraction_procedure_pages",
)
