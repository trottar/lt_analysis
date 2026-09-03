"""Detached D.6/D.7/D.8/D.9 procedure pages for the procedure PDF.

This module is presentation-only.  It receives already-built proton-cleaning
objects, clones only what it draws, and never rebuilds a fit, event lookup, or
weight.  Later procedure phases can append pages between the explicit PDF
open/close helpers without inheriting the technical diagnostic-PDF lifecycle.
"""

from __future__ import annotations

from array import array
import math
import os
from collections.abc import Mapping, Sequence

from canonical_binning import find_canonical_bin
from pion_component_subtraction import simc_shape_pion_weight_from_value


D6_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d6/v1"
D7_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d7/v1"
D8_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d8/v1"
D9_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d9/v1"
FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX = "_full-background-subtraction"

_TIMING_T_METHOD = "timing_t_event_weight"
_CTIME_AERO_METHOD = "ctime_aero_event_weight"
_D9_SIDES = ("kaon", "pion")
_D9_METHOD_A_COMPARISON_SCHEMA = "pion_hgcer_method_a_comparison/v1"


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
        legend = ROOT.TLegend(0.58, 0.67, 0.89, 0.87)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(kaon, "Proton-cleaned kaon data", "l")
        legend.AddEntry(pion, "Pion-control sample", "l")
        legend.Draw()
        draw_objects.extend((kaon, pion, legend))
        if thresholds.get("available"):
            _d9_draw_threshold(
                ROOT,
                thresholds["low_response_upper_threshold"],
                y_range,
                0.0,
                thresholds["low_response_upper_threshold"],
                draw_objects,
            )
            draw_objects.append(
                _draw_small_note(
                    ROOT, "Diagnostic only - Method A low response: 0 < HGCer NPE <= 2"
                )
            )
        else:
            draw_objects.append(
                _draw_small_note(ROOT, "Diagnostic only - Method A threshold unavailable")
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
            if thresholds.get("available"):
                _d9_draw_delta_threshold(
                    ROOT,
                    thresholds["low_response_upper_threshold"],
                    delta_edges[0],
                    delta_edges[-1],
                    draw_objects,
                )
            draw_objects.append(histogram)
        canvas.cd(1)
        notice = (
            "Diagnostic only - Method A low response: 0 < HGCer NPE <= 2"
            if thresholds.get("available")
            else "Diagnostic only - Method A threshold unavailable"
        )
        draw_objects.append(_draw_small_note(ROOT, notice))
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
    title = "Method A relative HGCer diagnostic - {}".format(_t_context(group))
    try:
        frame = ROOT.TH1D(
            "H_full_background_d9_method_a_relative_t{}".format(group["t_index"] + 1),
            "{};delta [%];Method A cell / same-|t| parent".format(title),
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
            _draw_small_note(ROOT, "Stored Method A values only - no refinement applied")
        )
        draw_objects.append(_draw_page_header(ROOT, canvas, "Method A relative diagnostic", group))
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


def render_full_background_subtraction_procedure_pages(
    pdf_name, d6_payload, d7_payload, d8_payload=None, d9_payload=None, *, page_manifest=None
):
    """Append available D.6 through D.9 groups in complete canonical-t order."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    d6 = _mapping(d6_payload)
    d7 = _mapping(d7_payload)
    d8 = _mapping(d8_payload)
    d9 = _mapping(d9_payload)
    d6_available = bool(d6.get("available"))
    d7_available = bool(d7.get("available"))
    d8_requested = d8_payload is not None
    d8_available = bool(d8.get("available"))
    d9_requested = d9_payload is not None
    d9_available = bool(d9.get("available"))
    if not d6_available and not d7_available and not d8_available and not d9_available:
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
    if d6_available and d7_available and list(d6.get("t_edges") or ()) != list(d7.get("t_edges") or ()):
        result["failures"].append("D.6/D.7 canonical t geometry mismatch")
        d7_available = False
    geometry_owner = (
        d6 if d6_available else d7 if d7_available else d8 if d8_available else d9
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
    elif d7_available:
        for group in tuple(d7.get("per_t") or ()):
            group = _mapping(group)
            _render_d7_t_pages(ROOT, pdf_name, d7, group, manifest, result["failures"])
            render_d8_group(group.get("t_index"))
            render_d9_group(group.get("t_index"))
    elif d8_available:
        for group in tuple(d8.get("per_t") or ()):
            group = _mapping(group)
            _render_d8_t_pages(ROOT, pdf_name, d8, group, manifest, result["failures"])
            render_d9_group(group.get("t_index"))
    elif d9_available:
        for group in tuple(d9.get("per_t") or ()):
            _render_d9_t_pages(ROOT, pdf_name, d9, _mapping(group), manifest, result["failures"])
    return result


__all__ = (
    "D6_PRESENTATION_SCHEMA_VERSION",
    "D7_PRESENTATION_SCHEMA_VERSION",
    "D8_PRESENTATION_SCHEMA_VERSION",
    "D9_PRESENTATION_SCHEMA_VERSION",
    "FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX",
    "build_full_background_subtraction_d6_payload",
    "build_full_background_subtraction_d7_payload",
    "build_full_background_subtraction_d8_payload",
    "build_full_background_subtraction_d9_payload",
    "close_full_background_subtraction_pdf",
    "full_background_subtraction_pdf_path",
    "open_full_background_subtraction_pdf",
    "render_full_background_subtraction_d6_pages",
    "render_full_background_subtraction_d7_pages",
    "render_full_background_subtraction_d8_pages",
    "render_full_background_subtraction_d9_pages",
    "render_full_background_subtraction_procedure_pages",
)
