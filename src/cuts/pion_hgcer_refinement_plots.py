"""Read-only PDF presentation for detached HGCer Phase A / Method A / Method B.

This module intentionally has no imports from the analysis, fitting, weighting, or
event-selection layers.  It consumes already-produced payloads and, when ROOT is
available, turns scalar records and display-only histogram clones into PDF pages.
It never evaluates a weight, changes a fit, or changes a production object.
"""

from __future__ import annotations

import math
import os
import textwrap


def build_pdf_destinations(main_pdf):
    """Return the deterministic main and automatic diagnostic PDF paths."""
    root, extension = os.path.splitext(os.fspath(main_pdf))
    extension = extension or ".pdf"
    return {
        "main": "{}{}".format(root, extension),
        "proton_debug": "{}_proton-debug{}".format(root, extension),
        "pion_fit_debug": "{}_pion-fit-debug{}".format(root, extension),
        "hgcer_debug": "{}_hgcer-debug{}".format(root, extension),
    }


def build_pdf_route_manifest(main_pdf):
    """Declare C.4 page ownership without using a mutable render-time policy."""
    destinations = build_pdf_destinations(main_pdf)
    routes = {
        "main": (
            "basic",
            "proton.summary.provenance_closure",
            "proton.summary.committed_mm",
            "proton.summary.commitment",
            "pion.canonical_t",
            "hgcer.phase_a.summary",
            "hgcer.method_a.support",
            "hgcer.method_a.f_low",
            "hgcer.method_b.status",
            "hgcer.method_b.candidate",
            "hgcer.method_b.regional_closure",
            "hgcer.method_b.shape_quality",
            "qa.setting_warnings",
        ),
        "proton_debug": ("proton.detail",),
        "pion_fit_debug": (
            "pion.setting_wide.detail",
            "pion.coordinate.detail",
            "pion.refinement.detail",
            "pion.dummy_random.detail",
        ),
        "hgcer_debug": ("hgcer.pid", "hgcer.part1", "hgcer.part1_5", "hgcer.part2"),
    }
    return {
        "destinations": destinations,
        "routes": {key: list(value) for key, value in routes.items()},
    }


def _mapping(value):
    return value if isinstance(value, dict) else {}


def _finite(value):
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    return scalar if math.isfinite(scalar) else None


def _short_fingerprint(value):
    text = str(value or "")
    return text[-16:] if text else "not recorded"


def _cell_key(cell):
    payload = _mapping(cell)
    try:
        return int(payload.get("t_index")), int(payload.get("delta_index"))
    except (TypeError, ValueError):
        return (10 ** 6, 10 ** 6)


def _edges(payload, name):
    values = list(_mapping(payload).get(name) or ())
    result = []
    for value in values:
        scalar = _finite(value)
        if scalar is None:
            return []
        result.append(scalar)
    return result if len(result) >= 2 else []


def _cell_status(cell, name, default="unavailable"):
    value = _mapping(cell).get(name)
    return str(value if value is not None else default).strip().lower() or default


_METHOD_A_SUPPORT_LABELS = {
    "unsupported": "unsupported",
    "marginal": "marginal",
    "supported": "supported",
}

_METHOD_B_STATUS_LABELS = {
    "unavailable": "unavailable",
    "marginal": "marginal",
    "shape_inconsistent": "shape inconsistent",
    "available": "available",
}

_METHOD_B_EXTRA_STATUS_LABELS = {
    "internally_inconsistent": "other stored inconsistency",
}

_SHAPE_STATUS_LABELS = {
    "unavailable": "unavailable",
    "marginal": "marginal",
    "good": "good",
    "poor": "poor",
}


def refinement_annotation_lines(candidate=False):
    """Return the stable visible diagnostic-only annotation for refinement pages."""
    lines = ["NON-AUTHORITATIVE DIAGNOSTIC / No refinement applied"]
    if candidate:
        lines.append("Diagnostic candidate only; no event correction")
    return tuple(lines)


def method_a_f_low_style(support_class):
    """Return distinct marker styles for the two displayed stored support classes."""
    return {"supported": 20, "marginal": 24}.get(str(support_class or "").lower())


def _stored_category_counts(cells, field, categories, summary_counts=None):
    """Count only serialized categories, preferring recorded cells when present."""
    rows = list(cells or ())
    if rows:
        counts = {category: 0 for category in categories}
        other = {}
        for cell in rows:
            status = _cell_status(cell, field)
            if status in counts:
                counts[status] += 1
            else:
                other[status] = other.get(status, 0) + 1
        return counts, other
    summary = _mapping(summary_counts)
    counts = {}
    other = {}
    for name, value in summary.items():
        try:
            count = int(value)
        except (TypeError, ValueError):
            count = 0
        if name in categories:
            counts[name] = count
        else:
            other[name] = count
    for category in categories:
        counts.setdefault(category, 0)
    return counts, other


def _display_value(value):
    """Render a recorded scalar or a compact recorded mapping without decisions."""
    if isinstance(value, dict):
        fields = ("status", "passed", "available", "reason", "count", "total")
        fragments = ["{}={}".format(name, value[name]) for name in fields if name in value]
        if not fragments:
            fragments = [
                "{}={}".format(name, item)
                for name, item in sorted(value.items())
                if not isinstance(item, (dict, list, tuple))
            ][:6]
        return ", ".join(fragments) if fragments else "recorded mapping"
    if value is None or value == "":
        return "not recorded"
    return str(value)


def phase_a_summary_payload(checkpoint, phase_a, display_context=None):
    """Extract only stored Phase-A contract fields for a compact overview."""
    checkpoint = _mapping(checkpoint)
    phase = _mapping(phase_a)
    setting = _mapping(checkpoint.get("setting"))
    phase_checkpoint = _mapping(checkpoint.get("phase_a"))
    summary = _mapping(phase.get("summary"))
    context = _mapping(display_context)
    if not summary:
        summary = _mapping(phase_checkpoint.get("summary"))
    return {
        "setting": dict(setting),
        "status": phase.get("status", summary.get("status", "unavailable")),
        "available": bool(phase.get("available", summary.get("available", False))),
        "reason": phase.get("reason", summary.get("reason")),
        "host_state": context.get("host_state", phase.get("host_state", phase_checkpoint.get("host_state"))),
        "source_target_state": phase.get("source_target_state", phase_checkpoint.get("source_target_state")),
        "pion_closure_passed": bool(
            _mapping(phase.get("pion_closure")).get("passed", summary.get("pion_closure_passed", False))
        ),
        "host_closure_passed": bool(
            _mapping(phase.get("host_closure")).get("passed", summary.get("host_closure_passed", False))
        ),
        "canonical_t_edges": _edges(phase, "canonical_t_edges") or list(checkpoint.get("canonical_t_edges") or ()),
        "delta_edges": _edges(phase, "delta_edges") or list(checkpoint.get("delta_edges") or ()),
        "lambda_gate_status": context.get("lambda_gate_status", phase.get("lambda_gate_status", summary.get("lambda_gate_status"))),
        "lambda_gate_production_action": context.get("production_action", phase.get(
            "lambda_gate_production_action", summary.get("lambda_gate_production_action")
        )),
        "proton_cleaning_committed": context.get("proton_cleaning_committed", "not recorded"),
        "contract_fingerprint": _short_fingerprint(
            phase.get("contract_fingerprint", phase_checkpoint.get("contract_fingerprint"))
        ),
        "coordinate_fingerprint": _short_fingerprint(
            phase.get("coordinate_fingerprint", phase_checkpoint.get("coordinate_fingerprint"))
        ),
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }


def method_a_plot_payload(method_a, checkpoint=None):
    """Return stored Method-A support and f_low cells without inference."""
    payload = _mapping(method_a)
    checkpoint_method = _mapping(_mapping(checkpoint).get("method_a"))
    summary = _mapping(payload.get("summary")) or _mapping(checkpoint_method.get("summary"))
    cells = [dict(_mapping(cell)) for cell in (payload.get("cells") or checkpoint_method.get("cells") or ())]
    cells.sort(key=_cell_key)
    support_counts, other_support_counts = _stored_category_counts(
        cells, "support_class", tuple(_METHOD_A_SUPPORT_LABELS), summary.get("support_counts")
    )
    return {
        "status": payload.get("status", checkpoint_method.get("status", "unavailable")),
        "available": bool(payload.get("available", checkpoint_method.get("available", False))),
        "reason": payload.get("reason", checkpoint_method.get("reason")),
        "t_edges": _edges(payload, "t_edges") or list(_mapping(checkpoint).get("canonical_t_edges") or ()),
        "delta_edges": _edges(payload, "delta_edges") or list(_mapping(checkpoint).get("delta_edges") or ()),
        "support_counts": support_counts,
        "other_support_counts": other_support_counts,
        "cells": cells,
        "non_authoritative": True,
        "no_interpolation": True,
        "no_correction_applied": True,
    }


def method_a_f_low_points(payload):
    """Select stored supported/marginal f_low observations without filling gaps."""
    points = []
    for cell in _mapping(payload).get("cells") or ():
        entry = _mapping(cell)
        if _cell_status(entry, "support_class") not in {"supported", "marginal"}:
            continue
        value = _finite(entry.get("f_low"))
        low = _finite(entry.get("f_low_low"))
        high = _finite(entry.get("f_low_high"))
        delta_low = _finite(entry.get("delta_low"))
        delta_high = _finite(entry.get("delta_high"))
        if None in (value, low, high, delta_low, delta_high):
            continue
        points.append({
            "t_index": entry.get("t_index"),
            "delta_index": entry.get("delta_index"),
            "delta_center": 0.5 * (delta_low + delta_high),
            "f_low": value,
            "f_low_low": low,
            "f_low_high": high,
            "support_class": _cell_status(entry, "support_class"),
        })
    return points


def method_b_plot_payload(method_b, checkpoint=None):
    """Return stored Method-B status, candidate, regional, and shape fields."""
    payload = _mapping(method_b)
    checkpoint_method = _mapping(_mapping(checkpoint).get("method_b"))
    summary = _mapping(payload.get("summary")) or _mapping(checkpoint_method.get("summary"))
    cells = [dict(_mapping(cell)) for cell in (payload.get("cells") or checkpoint_method.get("cells") or ())]
    cells.sort(key=_cell_key)
    method_status_counts, other_method_status_counts = _stored_category_counts(
        cells, "method_B_status", tuple(_METHOD_B_STATUS_LABELS), summary.get("method_B_status_counts")
    )
    return {
        "status": payload.get("status", checkpoint_method.get("status", "unavailable")),
        "available": bool(payload.get("available", checkpoint_method.get("available", False))),
        "reason": payload.get("reason", checkpoint_method.get("reason")),
        "t_edges": _edges(payload, "t_edges") or list(_mapping(checkpoint).get("canonical_t_edges") or ()),
        "delta_edges": _edges(payload, "delta_edges") or list(_mapping(checkpoint).get("delta_edges") or ()),
        "method_status_counts": method_status_counts,
        "other_method_status_counts": other_method_status_counts,
        "shape_status_counts": dict(summary.get("shape_status_counts") or {}),
        "candidate_status_counts": dict(summary.get("candidate_status_counts") or {}),
        "cells": cells,
        "non_authoritative": True,
        "frozen_pion_baseline": True,
        "no_refinement": True,
        "no_interpolation": True,
    }


def method_b_candidate_points(payload):
    """Select finite available multi-region candidate values; never substitute unity."""
    points = []
    for cell in _mapping(payload).get("cells") or ():
        entry = _mapping(cell)
        if _cell_status(entry, "candidate_L_B_status") != "available_multi_region":
            continue
        value = _finite(entry.get("candidate_L_B"))
        sigma = _finite(entry.get("candidate_L_B_uncertainty"))
        delta_low = _finite(entry.get("delta_low"))
        delta_high = _finite(entry.get("delta_high"))
        if None in (value, sigma, delta_low, delta_high):
            continue
        points.append({
            "t_index": entry.get("t_index"),
            "delta_index": entry.get("delta_index"),
            "delta_center": 0.5 * (delta_low + delta_high),
            "candidate_L_B": value,
            "candidate_L_B_uncertainty": sigma,
            "candidate_L_B_status": _cell_status(entry, "candidate_L_B_status"),
        })
    return points


def method_b_regional_rows(payload):
    """Return only the recorded C.4 regional closure rows for display."""
    selected = {"pi_n", "pi_sidis", "pi_delta_high"}
    rows = []
    for cell in _mapping(payload).get("cells") or ():
        entry = _mapping(cell)
        for region in entry.get("regions") or ():
            value = _mapping(region)
            if str(value.get("region_name") or "") not in selected:
                continue
            rows.append({
                "t_index": entry.get("t_index"),
                "delta_index": entry.get("delta_index"),
                "delta_low": entry.get("delta_low"),
                "delta_high": entry.get("delta_high"),
                "region_name": value.get("region_name"),
                "parent_relative_ratio": value.get("parent_relative_ratio"),
                "parent_relative_sigma": value.get("parent_relative_sigma"),
                "parent_relative_status": value.get("parent_relative_status"),
            })
    return rows


def method_b_regional_panels(payload):
    """Group recorded regional Qtilde ratios by their canonical-t parent."""
    regions = ("pi_n", "pi_sidis", "pi_delta_high")
    panels = {}
    for cell in _mapping(payload).get("cells") or ():
        try:
            t_index = int(_mapping(cell).get("t_index"))
        except (TypeError, ValueError):
            continue
        panels.setdefault(t_index, {"t_index": t_index, "series": {name: [] for name in regions}})
    for row in method_b_regional_rows(payload):
        if _cell_status(row, "parent_relative_status") != "available":
            continue
        low = _finite(row.get("delta_low"))
        high = _finite(row.get("delta_high"))
        ratio = _finite(row.get("parent_relative_ratio"))
        sigma = _finite(row.get("parent_relative_sigma"))
        if None in (low, high, ratio, sigma):
            continue
        try:
            t_index = int(row.get("t_index"))
            delta_index = int(row.get("delta_index"))
        except (TypeError, ValueError):
            continue
        region_name = str(row.get("region_name") or "")
        if region_name not in regions:
            continue
        panel = panels.setdefault(t_index, {"t_index": t_index, "series": {name: [] for name in regions}})
        panel["series"][region_name].append({
            "delta_index": delta_index,
            "delta_center": 0.5 * (low + high),
            "Qtilde": ratio,
            "Qtilde_uncertainty": sigma,
        })
    result = []
    for t_index in sorted(panels):
        panel = panels[t_index]
        for points in panel["series"].values():
            points.sort(key=lambda entry: entry["delta_index"])
        result.append(panel)
    return result


def unity_line_limits(payload):
    """Return the recorded full delta range for every reference-unity line."""
    return canonical_delta_frame_limits(payload)


def canonical_delta_frame_limits(payload):
    """Return the full recorded canonical delta domain for a plot frame."""
    edges = _edges(payload, "delta_edges")
    return (edges[0], edges[-1]) if len(edges) >= 2 else (None, None)


def display_y_range(points, value_key, uncertainty_key=None, *, include_values=()):
    """Return a finite display range from stored values/errors without selection changes."""
    limits = []
    for point in points or ():
        entry = _mapping(point)
        value = _finite(entry.get(value_key))
        if value is None:
            continue
        uncertainty = _finite(entry.get(uncertainty_key)) if uncertainty_key else None
        uncertainty = abs(uncertainty) if uncertainty is not None else 0.0
        limits.extend((value - uncertainty, value + uncertainty))
    for value in include_values or ():
        scalar = _finite(value)
        if scalar is not None:
            limits.append(scalar)
    if not limits:
        return (None, None)
    lower, upper = min(limits), max(limits)
    if upper == lower:
        margin = max(abs(lower) * 0.05, 0.05)
    else:
        margin = 0.10 * (upper - lower)
    return (lower - margin, upper + margin)


def method_a_f_low_frame_limits(points):
    """Use stored Method-A asymmetric intervals to make one parent frame."""
    limits = []
    for point in points or ():
        value = _finite(_mapping(point).get("f_low"))
        low = _finite(_mapping(point).get("f_low_low"))
        high = _finite(_mapping(point).get("f_low_high"))
        if None not in (value, low, high):
            limits.append({"value": value, "uncertainty": max(abs(value - low), abs(high - value))})
    return display_y_range(limits, "value", "uncertainty")


def method_b_candidate_frame_limits(points):
    """Use stored candidate intervals and retain the visible unity reference."""
    return display_y_range(points, "candidate_L_B", "candidate_L_B_uncertainty", include_values=(1.0,))


def method_b_regional_frame_limits(panel):
    """Use every displayed stored regional interval for one parent-local frame."""
    rows = []
    for points in _mapping(panel).get("series", {}).values():
        rows.extend(points or ())
    return display_y_range(rows, "Qtilde", "Qtilde_uncertainty", include_values=(1.0,))


def method_b_shape_rows(payload):
    """Extract stored shape metrics and their categorical state without thresholds."""
    return [
        {
            "t_index": _mapping(cell).get("t_index"),
            "delta_index": _mapping(cell).get("delta_index"),
            "shape_chi2_ndf": _mapping(cell).get("shape_chi2_ndf"),
            "shape_max_abs_pull": _mapping(cell).get("shape_max_abs_pull"),
            "shape_status": _mapping(cell).get("shape_status"),
            "candidate_L_B_status": _mapping(cell).get("candidate_L_B_status"),
            "shape_poor_veto": _cell_status(cell, "candidate_L_B_status") == "shape_poor_veto",
        }
        for cell in _mapping(payload).get("cells") or ()
    ]


def canonical_parent_lambda_summary(parent_rows):
    """Summarize only pre-existing canonical-parent K-Lambda display records."""
    entries = []
    for row in parent_rows or ():
        entry = _mapping(row)
        try:
            t_index = int(entry.get("t_index"))
        except (TypeError, ValueError):
            continue
        status = str(entry.get("status") or "not recorded").strip().lower()
        entries.append({
            "t_index": t_index,
            "status": status or "not recorded",
            "reason": entry.get("reason") or "not recorded",
        })
    entries.sort(key=lambda entry: entry["t_index"])
    available = sum(entry["status"] in ("available", "ok", "pass", "passed") for entry in entries)
    return {
        "entries": entries,
        "available": available,
        "total": len(entries),
        "unavailable": len(entries) - available,
    }


def _pass_fail_label(value):
    if value is True:
        return "PASS"
    if value is False:
        return "FAIL"
    return "not recorded"


def setting_qa_summary_payload(phase_a, method_a, method_b, part2=None, display_context=None, runtime_qa_context=None):
    """Collect existing setting-level coverage and QA state for terminal display."""
    phase = _mapping(phase_a)
    method_a_payload = method_a_plot_payload(method_a)
    method_b_payload = method_b_plot_payload(method_b)
    context = _mapping(display_context)
    runtime = _mapping(runtime_qa_context)
    lambda_summary = canonical_parent_lambda_summary(runtime.get("canonical_parent_k_lambda"))
    return {
        "method_a_coverage": dict(method_a_payload.get("support_counts") or {}),
        "method_b_coverage": dict(method_b_payload.get("method_status_counts") or {}),
        "method_b_other_coverage": dict(method_b_payload.get("other_method_status_counts") or {}),
        "phase_a_coordinate_status": phase.get("status") or "unavailable",
        "phase_a_pion_closure": _mapping(phase.get("pion_closure")).get("passed", "not recorded"),
        "phase_a_host_closure": _mapping(phase.get("host_closure")).get("passed", "not recorded"),
        "phase_host_state": context.get("host_state", phase.get("host_state", "not recorded")),
        "lambda_gate_status": context.get("lambda_gate_status", phase.get("lambda_gate_status", "not recorded")),
        "proton_action": context.get("production_action", phase.get("lambda_gate_production_action", "not recorded")),
        "proton_cleaning_committed": context.get("proton_cleaning_committed", "not recorded"),
        "aerogel_warnings": runtime.get("aerogel_warnings", "not available"),
        "proton_warnings": runtime.get("proton_warnings", "not available"),
        "canonical_parent_k_lambda": lambda_summary,
        "k_lambda_comparison": runtime.get("k_lambda_comparison", "not available"),
        "k_sigma0_protected_region": runtime.get("k_sigma0_protected_region", "not recorded"),
        "k_sigma0_availability": runtime.get("k_sigma0_availability", "not available"),
        "hgcer_diagnostic_availability": runtime.get("hgcer_diagnostic_availability", "not available"),
        "renderer_failures": runtime.get("renderer_failures", "not available"),
        "part2_status": _mapping(part2).get("status") or "not available",
    }


def _has_recorded_warning(value):
    if value is None or value in ("", "not available", "not recorded"):
        return False
    if isinstance(value, dict):
        return any(_has_recorded_warning(item) for item in value.values())
    return bool(value)


def _status_category(value, *, required=False, optional=False):
    token = str(value or "not recorded").strip().lower()
    if token in ("available", "ok", "pass", "passed", "true"):
        return "OK"
    if token in ("not available", "not recorded", "", "none"):
        return "INFORMATIONAL" if optional else ("WARNING" if required else "NOT_RECORDED")
    if token in ("warning", "marginal", "insufficient_support", "bypass"):
        return "WARNING"
    if token in ("unavailable", "failed", "fail", "failure", "error"):
        return "FAILURE" if required else "WARNING"
    return "WARNING" if required else "NOT_RECORDED"


def setting_qa_warning_states(qa, detached_warnings=()):
    """Classify existing presentation states without adding scientific criteria."""
    summary = _mapping(qa)
    states = []

    def add(label, category, detail):
        if category != "OK":
            states.append({"label": label, "category": category, "detail": detail})

    add("Phase A", _status_category(summary.get("phase_a_coordinate_status"), required=True), summary.get("phase_a_coordinate_status"))
    for label, value in (
        ("Phase-A pion closure", summary.get("phase_a_pion_closure")),
        ("Phase-A host closure", summary.get("phase_a_host_closure")),
    ):
        add(label, "OK" if value is True else ("FAILURE" if value is False else "WARNING"), _pass_fail_label(value))
    if _has_recorded_warning(summary.get("aerogel_warnings")):
        add("Aerogel diagnostics", "WARNING", _display_value(summary.get("aerogel_warnings")))
    if _has_recorded_warning(summary.get("proton_warnings")):
        add("Proton diagnostics", "WARNING", _display_value(summary.get("proton_warnings")))
    add(
        "HGCer diagnostics",
        _status_category(summary.get("hgcer_diagnostic_availability"), required=True),
        _display_value(summary.get("hgcer_diagnostic_availability")),
    )
    for parent in _mapping(summary.get("canonical_parent_k_lambda")).get("entries") or ():
        category = _status_category(parent.get("status"), required=True)
        add("Canonical-parent K-Lambda t{}".format(int(parent["t_index"]) + 1), category, parent.get("reason"))
    if not _mapping(summary.get("canonical_parent_k_lambda")).get("entries"):
        add("Canonical-parent K-Lambda", "WARNING", "not recorded")
    if _has_recorded_warning(summary.get("renderer_failures")):
        add("Required main-PDF renderer", "FAILURE", _display_value(summary.get("renderer_failures")))
    k_sigma0_detail = "protected region={}; explicit template={}".format(
        summary.get("k_sigma0_protected_region"), summary.get("k_sigma0_availability")
    )
    add("K-Sigma0 explicit template", "INFORMATIONAL", k_sigma0_detail)
    for entry in detached_warnings or ():
        item = _mapping(entry)
        scope = str(item.get("scope") or "detached diagnostic")
        status = str(item.get("status") or "not recorded").lower()
        if scope == "Part 2":
            category = "INFORMATIONAL"
        elif scope in ("Method A", "Method B") and status == "marginal":
            category = "INFORMATIONAL"
        else:
            category = _status_category(status, required=scope == "Phase A")
        add(scope, category, item.get("reason") or "not recorded")
    return states


def setting_qa_summary_lines(qa, detached_warnings=()):
    """Build the terminal QA text from stored display state and its classifier."""
    summary = _mapping(qa)
    lambda_summary = _mapping(summary.get("canonical_parent_k_lambda"))
    lines = ["Setting-level coverage + QA summary"]
    lines.extend((
        "Host/proton state: host={}; Lambda gate={}; action={}; committed={}".format(
            summary.get("phase_host_state"), summary.get("lambda_gate_status"),
            summary.get("proton_action"), summary.get("proton_cleaning_committed"),
        ),
        "Phase A: status={}; pion closure={}; host closure={}".format(
            summary.get("phase_a_coordinate_status"), _pass_fail_label(summary.get("phase_a_pion_closure")),
            _pass_fail_label(summary.get("phase_a_host_closure")),
        ),
        "Method A coverage: {}".format(_display_value(summary.get("method_a_coverage"))),
        "Method B coverage: {}".format(_display_value(summary.get("method_b_coverage"))),
        "Canonical-parent K-Lambda QA: {} / {} available".format(
            lambda_summary.get("available", 0), lambda_summary.get("total", 0),
        ),
    ))
    for parent in lambda_summary.get("entries") or ():
        lines.append("  t{}: {} — {}".format(int(parent["t_index"]) + 1, parent["status"], parent["reason"]))
    lines.extend((
        "K-Sigma0 protected region: {}; setting-wide explicit template: {}".format(
            summary.get("k_sigma0_protected_region"), summary.get("k_sigma0_availability"),
        ),
        "Aerogel warnings: {}; proton warnings: {}; HGCer diagnostics: {}".format(
            _display_value(summary.get("aerogel_warnings")), _display_value(summary.get("proton_warnings")),
            _display_value(summary.get("hgcer_diagnostic_availability")),
        ),
        "Renderer failures: {}".format(_display_value(summary.get("renderer_failures"))),
    ))
    if summary.get("method_b_other_coverage"):
        lines.append("Method B other recorded states: {}".format(_display_value(summary.get("method_b_other_coverage"))))
    states = setting_qa_warning_states(summary, detached_warnings)
    actionable = [state for state in states if state["category"] in ("WARNING", "FAILURE")]
    if actionable:
        lines.append("Outstanding setting-level QA states:")
        lines.extend("{}: {} — {}".format(state["category"], state["label"], state["detail"]) for state in actionable)
    else:
        lines.append("No outstanding setting-level QA warnings.")
    informational = [state for state in states if state["category"] in ("INFORMATIONAL", "NOT_RECORDED")]
    if informational:
        lines.append("Information:")
        lines.extend("{}: {} — {}".format(state["category"], state["label"], state["detail"]) for state in informational)
    return lines


def warning_payload(phase_a, method_a, method_b, part2=None, extra=None):
    """Preserve non-OK stored states for the final main-PDF warning page."""
    warnings = []
    for scope, payload in (
        ("Phase A", phase_a),
        ("Method A", method_a),
        ("Method B", method_b),
        ("Part 2", part2),
    ):
        entry = _mapping(payload)
        if not entry:
            continue
        status = str(entry.get("status") or "unavailable").lower()
        if status not in ("available", "ok", "pass", "passed"):
            warnings.append({
                "scope": scope,
                "status": status,
                "reason": entry.get("reason") or entry.get("diagnostic_stage") or "not recorded",
                "production_impact": "none; detached diagnostic only",
            })
    for entry in extra or ():
        item = _mapping(entry)
        if item:
            warnings.append(dict(item))
    return warnings


def _title(setting, section, diagnostic):
    setting = _mapping(setting)
    phi = setting.get("phi_setting") or "unknown phi"
    kinematic = setting.get("kinematic_token") or "unknown kinematic"
    epsilon = setting.get("epsilon_filename_token") or setting.get("epsilon_setting") or "unknown epsilon"
    return "{} | {} | {} — {}: {}".format(phi, kinematic, epsilon, section, diagnostic)


def _import_root():
    try:
        import ROOT
    except Exception:
        return None
    return ROOT


def _text_canvas(ROOT, name, title, lines, *, width=1200, height=800):
    canvas = ROOT.TCanvas(str(name), str(title), int(width), int(height))
    text = ROOT.TPaveText(0.04, 0.05, 0.96, 0.95, "NDC")
    text.SetFillStyle(0)
    text.SetBorderSize(0)
    text.SetTextAlign(12)
    text.SetTextSize(0.024)
    text.AddText(str(title))
    for line in lines:
        for wrapped in textwrap.wrap(str(line), width=112) or [""]:
            text.AddText(wrapped)
    text.Draw()
    return canvas, text


def _print_text_page(pdf_name, page_id, title, lines, manifest):
    ROOT = _import_root()
    if ROOT is None:
        return False
    canvas = None
    try:
        canvas, _text = _text_canvas(ROOT, "hgcer_plot_{}".format(page_id.replace(".", "_")), title, lines)
        canvas.Print(pdf_name)
        manifest.append({"page_id": page_id, "scope": "setting", "authoritative": False})
        return True
    finally:
        if canvas is not None:
            canvas.Close()


def open_diagnostic_pdf(pdf_name, checkpoint, *, role, main_pdf):
    """Open a supplement with a real provenance page, never an empty canvas."""
    setting = _mapping(_mapping(checkpoint).get("setting"))
    manifest = []
    _print_text_page(
        "{}(".format(pdf_name),
        "{}.provenance".format(role),
        _title(setting, "Diagnostic supplement", "provenance"),
        (
            "role: {}".format(role),
            "main PDF: {}".format(os.path.basename(os.fspath(main_pdf))),
            "authority: production subtraction remains the main PDF and frozen analysis payloads",
            "non_authoritative=true; production_objects_mutated=false; refinement_applied=false",
        ),
        manifest,
    )
    return manifest


def close_diagnostic_pdf(pdf_name, checkpoint, *, role, manifest):
    """Close a supplement on a real terminal provenance page."""
    setting = _mapping(_mapping(checkpoint).get("setting"))
    return _print_text_page(
        "{}".format(pdf_name) + ")",
        "{}.terminal".format(role),
        _title(setting, "Diagnostic supplement", "terminal status"),
        (
            "role: {}".format(role),
            "terminal presentation only; no analysis object was changed",
        ),
        manifest,
    )


def _clone_display(histogram, name):
    if histogram is None:
        return None
    try:
        clone = histogram.Clone(str(name))
        if hasattr(clone, "SetDirectory"):
            clone.SetDirectory(0)
        return clone
    except Exception:
        return None


def _display_frame(ROOT, name, title, x_limits, y_limits, y_label):
    """Create a detached ROOT frame with the prescribed visible ranges."""
    xmin, xmax = x_limits
    ymin, ymax = y_limits
    if None in (xmin, xmax, ymin, ymax):
        return None
    frame = ROOT.TH1D(str(name), "{};SHMS delta;{}".format(title, y_label), 1, float(xmin), float(xmax))
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(float(ymin))
    frame.SetMaximum(float(ymax))
    return frame


def _diagnostic_label(ROOT, *, candidate=False):
    """Create the mandatory visible detached-diagnostic label for a plot page."""
    label = ROOT.TPaveText(0.10, 0.90, 0.90, 0.985, "NDC")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.SetTextAlign(12)
    label.SetTextSize(0.022)
    for line in refinement_annotation_lines(candidate):
        label.AddText(line)
    return label


def _category_legend(ROOT, labels, category_codes):
    legend = ROOT.TPaveText(0.73, 0.10, 0.96, 0.30, "NDC")
    legend.SetFillColor(0)
    legend.SetBorderSize(1)
    legend.SetTextAlign(12)
    legend.SetTextSize(0.021)
    legend.AddText("categorical map key")
    for category, code in category_codes.items():
        legend.AddText("{}: {}".format(int(code), labels.get(category, category)))
    return legend


def _draw_map_page(ROOT, pdf_name, page_id, title, payload, *, value_key, category_codes=None, category_labels=None, manifest=None):
    t_edges = list(payload.get("t_edges") or ())
    delta_edges = list(payload.get("delta_edges") or ())
    cells = list(payload.get("cells") or ())
    if len(t_edges) < 2 or len(delta_edges) < 2:
        return _print_text_page(pdf_name, page_id, title, ("unavailable: canonical t/delta edges not recorded",) + refinement_annotation_lines(), manifest)
    from array import array
    histogram = ROOT.TH2D(
        "H_{}_{}".format(page_id.replace(".", "_"), abs(id(cells))),
        "{};SHMS delta;|t| [GeV^2]".format(title),
        len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges),
    )
    histogram.SetDirectory(0)
    for cell in cells:
        try:
            xbin = int(cell.get("delta_index")) + 1
            ybin = int(cell.get("t_index")) + 1
        except (TypeError, ValueError):
            continue
        value = cell.get(value_key)
        if category_codes is not None:
            value = category_codes.get(str(value or "unavailable").lower(), 0.0)
        value = _finite(value)
        if value is not None:
            histogram.SetBinContent(xbin, ybin, value)
    canvas = ROOT.TCanvas("C_{}".format(page_id.replace(".", "_")), title, 1100, 800)
    try:
        draw_objects = [histogram]
        histogram.Draw("colz text")
        label = _diagnostic_label(ROOT)
        draw_objects.append(label)
        label.Draw()
        if category_codes is not None:
            legend = _category_legend(ROOT, category_labels or {}, category_codes)
            draw_objects.append(legend)
            legend.Draw()
        canvas.Print(pdf_name)
        manifest.append({"page_id": page_id, "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def _draw_shape_quality_page(ROOT, pdf_name, title, payload, manifest):
    """Render stored shape values and recorded shape-veto state side by side."""
    t_edges = list(payload.get("t_edges") or ())
    delta_edges = list(payload.get("delta_edges") or ())
    if len(t_edges) < 2 or len(delta_edges) < 2:
        return _print_text_page(pdf_name, "hgcer.method_b.shape_quality", title, ("unavailable: canonical t/delta edges not recorded",) + refinement_annotation_lines(), manifest)
    from array import array
    suffix = abs(id(payload.get("cells")))
    maps = {
        "chi2": ROOT.TH2D("H_hgcer_method_b_shape_chi2_{}".format(suffix), "stored shape chi2/ndf;SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
        "pull": ROOT.TH2D("H_hgcer_method_b_shape_pull_{}".format(suffix), "stored maximum |pull|;SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
        "status": ROOT.TH2D("H_hgcer_method_b_shape_status_{}".format(suffix), "shape status (0 unavailable, 1 marginal, 2 good, 3 poor);SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
        "veto": ROOT.TH2D("H_hgcer_method_b_shape_veto_{}".format(suffix), "stored candidate shape-poor veto (0 no, 1 veto);SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
    }
    for histogram in maps.values():
        histogram.SetDirectory(0)
    codes = {"unavailable": 0.0, "marginal": 1.0, "good": 2.0, "poor": 3.0}
    for row in method_b_shape_rows(payload):
        try:
            xbin, ybin = int(row.get("delta_index")) + 1, int(row.get("t_index")) + 1
        except (TypeError, ValueError):
            continue
        chi2 = _finite(row.get("shape_chi2_ndf"))
        pull = _finite(row.get("shape_max_abs_pull"))
        if chi2 is not None:
            maps["chi2"].SetBinContent(xbin, ybin, chi2)
        if pull is not None:
            maps["pull"].SetBinContent(xbin, ybin, pull)
        maps["status"].SetBinContent(xbin, ybin, codes.get(str(row.get("shape_status") or "unavailable").lower(), 0.0))
        maps["veto"].SetBinContent(xbin, ybin, 1.0 if row.get("shape_poor_veto") else 0.0)
    canvas = ROOT.TCanvas("C_hgcer_method_b_shape_quality", title, 1900, 650)
    try:
        draw_objects = list(maps.values())
        canvas.Divide(4, 1)
        for pad, key in enumerate(("chi2", "pull", "status", "veto"), start=1):
            canvas.cd(pad)
            maps[key].Draw("colz text")
            label = _diagnostic_label(ROOT)
            draw_objects.append(label)
            label.Draw()
            if key == "status":
                legend = _category_legend(ROOT, _SHAPE_STATUS_LABELS, codes)
                draw_objects.append(legend)
                legend.Draw()
            elif key == "veto":
                legend = _category_legend(ROOT, {"not_vetoed": "not vetoed", "shape_poor_veto": "shape-poor veto"}, {"not_vetoed": 0.0, "shape_poor_veto": 1.0})
                draw_objects.append(legend)
                legend.Draw()
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.shape_quality", "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def _render_phase_a_page(pdf_name, checkpoint, phase_a, phase_a_display_context, manifest):
    payload = phase_a_summary_payload(checkpoint, phase_a, phase_a_display_context)
    setting = payload["setting"]
    lines = (
        "Phase A status: {} (available={})".format(payload["status"], payload["available"]),
        "host state: {}; source target: {}".format(payload["host_state"] or "not recorded", payload["source_target_state"] or "not recorded"),
        "pion closure: {}; host closure: {}".format(payload["pion_closure_passed"], payload["host_closure_passed"]),
        "canonical |t| edges: {}".format(payload["canonical_t_edges"] or "not recorded"),
        "delta edges: {}".format(payload["delta_edges"] or "not recorded"),
        "proton/Lambda gate: {} / {}; committed={}".format(payload["lambda_gate_status"] or "not recorded", payload["lambda_gate_production_action"] or "not recorded", payload["proton_cleaning_committed"]),
        "contract fingerprint: {}; coordinate fingerprint: {}".format(payload["contract_fingerprint"], payload["coordinate_fingerprint"]),
        "NON-AUTHORITATIVE DIAGNOSTIC / No refinement applied",
        "production_objects_mutated=false; refinement_applied=false",
        "reason: {}".format(payload["reason"] or "none"),
    )
    return _print_text_page(pdf_name, "hgcer.phase_a.summary", _title(setting, "Phase A", "HGCer Event/Population Contract"), lines, manifest)


def _render_method_a_pages(pdf_name, checkpoint, method_a, manifest):
    payload = method_a_plot_payload(method_a, checkpoint)
    setting = _mapping(_mapping(checkpoint).get("setting"))
    category_codes = {"unsupported": 0.0, "marginal": 1.0, "supported": 2.0}
    ROOT = _import_root()
    if ROOT is None:
        return False
    _draw_map_page(ROOT=ROOT, pdf_name=pdf_name, page_id="hgcer.method_a.support", title=_title(setting, "Method A", "support map"), payload=payload, value_key="support_class", category_codes=category_codes, category_labels=_METHOD_A_SUPPORT_LABELS, manifest=manifest)
    canvas = ROOT.TCanvas("C_hgcer_method_a_f_low", _title(setting, "Method A", "f_low by canonical t"), 1200, 800)
    try:
        draw_objects = []
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        points = method_a_f_low_points(payload)
        delta_limits = canonical_delta_frame_limits(payload)
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            parent_points = [point for point in points if int(point.get("t_index", -1)) == t_index]
            y_limits = method_a_f_low_frame_limits(parent_points)
            if not parent_points or None in delta_limits or None in y_limits:
                message = ROOT.TLatex()
                message.DrawLatexNDC(0.12, 0.55, "No supported or marginal f_low cells")
                label = _diagnostic_label(ROOT)
                label.Draw()
                draw_objects.extend((message, label))
                continue
            frame = _display_frame(
                ROOT,
                "H_hgcer_method_a_f_low_frame_{}".format(t_index),
                "{} t{}".format(_title(setting, "Method A", "f_low"), t_index + 1),
                delta_limits,
                y_limits,
                "f_low",
            )
            draw_objects.append(frame)
            frame.Draw("AXIS")
            legend = ROOT.TLegend(0.64, 0.72, 0.90, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            draw_objects.append(legend)
            has_series = False
            for support_class in ("supported", "marginal"):
                graph = ROOT.TGraphErrors()
                draw_objects.append(graph)
                point = 0
                for point_data in parent_points:
                    if int(point_data.get("t_index", -1)) != t_index or point_data.get("support_class") != support_class:
                        continue
                    value = point_data["f_low"]
                    graph.SetPoint(point, point_data["delta_center"], value)
                    graph.SetPointError(point, 0.0, max(value - point_data["f_low_low"], point_data["f_low_high"] - value))
                    point += 1
                if point:
                    graph.SetMarkerStyle(method_a_f_low_style(support_class))
                    graph.Draw("P same")
                    legend.AddEntry(graph, support_class, "p")
                    has_series = True
                    plotted = True
            if has_series:
                legend.Draw()
                label = _diagnostic_label(ROOT)
                label.Draw()
                draw_objects.append(label)
        if not plotted:
            canvas.cd(1)
            message = ROOT.TLatex()
            message.DrawLatexNDC(0.12, 0.55, "No supported or marginal f_low cells")
            label = _diagnostic_label(ROOT)
            label.Draw()
            draw_objects.extend((message, label))
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_a.f_low", "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def _render_method_b_pages(pdf_name, checkpoint, method_b, manifest):
    payload = method_b_plot_payload(method_b, checkpoint)
    setting = _mapping(_mapping(checkpoint).get("setting"))
    ROOT = _import_root()
    if ROOT is None:
        return False
    status_codes = {"unavailable": 0.0, "marginal": 1.0, "shape_inconsistent": 2.0, "available": 3.0, "internally_inconsistent": 4.0}
    status_labels = dict(_METHOD_B_STATUS_LABELS)
    status_labels.update(_METHOD_B_EXTRA_STATUS_LABELS)
    _draw_map_page(ROOT, pdf_name, "hgcer.method_b.status", _title(setting, "Method B", "cell status"), payload, value_key="method_B_status", category_codes=status_codes, category_labels=status_labels, manifest=manifest)

    canvas = ROOT.TCanvas("C_hgcer_method_b_candidate", _title(setting, "Method B", "candidate L_B"), 1200, 800)
    try:
        draw_objects = []
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        candidates = method_b_candidate_points(payload)
        delta_limits = canonical_delta_frame_limits(payload)
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            parent_points = [point for point in candidates if int(point.get("t_index", -1)) == t_index]
            y_limits = method_b_candidate_frame_limits(parent_points)
            if not parent_points or None in delta_limits or None in y_limits:
                message = ROOT.TLatex()
                message.DrawLatexNDC(0.12, 0.55, "No available multi-region candidate cells")
                label = _diagnostic_label(ROOT, candidate=True)
                label.Draw()
                draw_objects.extend((message, label))
                continue
            frame = _display_frame(
                ROOT,
                "H_hgcer_method_b_candidate_frame_{}".format(t_index),
                "{} t{}".format(_title(setting, "Method B", "candidate L_B"), t_index + 1),
                delta_limits,
                y_limits,
                "L_B candidate",
            )
            draw_objects.append(frame)
            frame.Draw("AXIS")
            unity = ROOT.TLine(delta_limits[0], 1.0, delta_limits[1], 1.0)
            draw_objects.append(unity)
            graph = ROOT.TGraphErrors()
            draw_objects.append(graph)
            point = 0
            for point_data in parent_points:
                if int(point_data.get("t_index", -1)) != t_index:
                    continue
                value, sigma = point_data["candidate_L_B"], point_data["candidate_L_B_uncertainty"]
                center = point_data["delta_center"]
                graph.SetPoint(point, center, value)
                graph.SetPointError(point, 0.0, sigma)
                point += 1
            if point:
                graph.SetMarkerStyle(20)
                graph.Draw("P same")
                unity.SetLineStyle(2)
                unity.Draw("same")
                label = _diagnostic_label(ROOT, candidate=True)
                label.Draw()
                draw_objects.append(label)
                plotted = True
        if not plotted:
            canvas.cd(1)
            message = ROOT.TLatex()
            message.DrawLatexNDC(0.12, 0.55, "No available multi-region candidate cells")
            label = _diagnostic_label(ROOT, candidate=True)
            label.Draw()
            draw_objects.extend((message, label))
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.candidate", "scope": "setting", "authoritative": False})
    finally:
        canvas.Close()

    regional_panels = method_b_regional_panels(payload)
    regional_canvas = ROOT.TCanvas("C_hgcer_method_b_regional", _title(setting, "Method B", "regional closure"), 1500, 800)
    try:
        draw_objects = []
        regions = ("pi_n", "pi_sidis", "pi_delta_high")
        regional_canvas.Divide(max(1, min(len(regional_panels), 4)), 1)
        plotted = False
        for pad, panel in enumerate(regional_panels[:4], start=1):
            regional_canvas.cd(pad)
            delta_limits = canonical_delta_frame_limits(payload)
            y_limits = method_b_regional_frame_limits(panel)
            panel_points = [point for points in panel["series"].values() for point in points]
            if not panel_points or None in delta_limits or None in y_limits:
                message = ROOT.TLatex()
                message.DrawLatexNDC(0.12, 0.55, "No stored available regional rows for this t parent")
                label = _diagnostic_label(ROOT)
                label.Draw()
                draw_objects.extend((message, label))
                continue
            frame = _display_frame(
                ROOT,
                "H_hgcer_method_b_regional_frame_{}".format(panel["t_index"]),
                "{} t{}".format(_title(setting, "Method B", "regional closure"), panel["t_index"] + 1),
                delta_limits,
                y_limits,
                "Qtilde stored ratio",
            )
            draw_objects.append(frame)
            frame.Draw("AXIS")
            legend = ROOT.TLegend(0.57, 0.68, 0.90, 0.89)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            draw_objects.append(legend)
            has_series = False
            for region_name, marker_style in zip(regions, (20, 24, 25)):
                graph = ROOT.TGraphErrors()
                draw_objects.append(graph)
                point = 0
                for row in panel["series"][region_name]:
                    graph.SetPoint(point, row["delta_center"], row["Qtilde"])
                    graph.SetPointError(point, 0.0, row["Qtilde_uncertainty"])
                    point += 1
                if point:
                    graph.SetMarkerStyle(marker_style)
                    graph.Draw("P same")
                    legend.AddEntry(graph, region_name, "p")
                    has_series = True
                    plotted = True
            if has_series:
                line = ROOT.TLine(delta_limits[0], 1.0, delta_limits[1], 1.0)
                line.SetLineStyle(2)
                draw_objects.append(line)
                line.Draw("same")
                legend.Draw()
                label = _diagnostic_label(ROOT)
                label.Draw()
                draw_objects.append(label)
        if not plotted:
            regional_canvas.cd(1)
            message = ROOT.TLatex()
            message.DrawLatexNDC(0.12, 0.55, "No stored available regional rows")
            label = _diagnostic_label(ROOT)
            label.Draw()
            draw_objects.extend((message, label))
        regional_canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.regional_closure", "scope": "setting", "authoritative": False})
    finally:
        regional_canvas.Close()

    _draw_shape_quality_page(ROOT, pdf_name, _title(setting, "Method B", "shape quality"), payload, manifest)
    return True


def render_pion_hgcer_refinement_pages(pdf_name, checkpoint, *, phase_a=None, method_a=None, method_b=None, phase_a_display_context=None, page_manifest=None):
    """Append all new detached Phase-A/Method-A/Method-B pages to an open PDF."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    if _import_root() is None:
        return manifest
    _render_phase_a_page(pdf_name, checkpoint, phase_a, phase_a_display_context, manifest)
    _render_method_a_pages(pdf_name, checkpoint, method_a, manifest)
    _render_method_b_pages(pdf_name, checkpoint, method_b, manifest)
    return manifest


def render_setting_warning_page(pdf_name, checkpoint, *, phase_a=None, method_a=None, method_b=None, part2=None, phase_a_display_context=None, runtime_qa_context=None, page_manifest=None, close_pdf=False):
    """Append the terminal setting-level QA summary to an already-open main PDF."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    warnings = warning_payload(phase_a, method_a, method_b, part2)
    qa = setting_qa_summary_payload(
        phase_a, method_a, method_b, part2,
        display_context=phase_a_display_context,
        runtime_qa_context=runtime_qa_context,
    )
    setting = _mapping(_mapping(checkpoint).get("setting"))
    lines = setting_qa_summary_lines(qa, warnings)
    target_pdf = "{}".format(pdf_name) + ")" if close_pdf else pdf_name
    _print_text_page(target_pdf, "qa.setting_warnings", _title(setting, "Warnings", "detached diagnostic status"), lines, manifest)
    return manifest


def proton_main_qa_payload(cleaning_result, cleaning_application, phase_a_display_context=None):
    """Select recorded proton QA and committed-spectrum objects for the main PDF."""
    result = _mapping(cleaning_result)
    diagnostics = _mapping(result.get("diagnostics"))
    application = _mapping(cleaning_application)
    application_diagnostics = _mapping(application.get("diagnostics"))
    gate = _mapping(diagnostics.get("lambda_preservation_gate"))
    context = _mapping(phase_a_display_context)
    host_state = context.get(
        "host_state",
        application.get("host_state", result.get("source_target_state", "not recorded")),
    )
    identity_host = host_state == "identity_no_proton_cleaning"
    if identity_host:
        spectra = {
            "raw": _mapping(application.get("raw_targets")).get("h_mm_nosub"),
            "estimated": _mapping(application.get("proton_targets")).get("h_mm_nosub"),
            "cleaned": _mapping(application.get("cleaned_targets_pre_rf")).get("h_mm_nosub"),
            "committed": _mapping(application.get("final_targets")).get("h_mm_nosub"),
        }
        numerical_closure = "not recorded"
        global_closure = "not recorded"
        identity_host_closure = application_diagnostics.get(
            "identity_host_closure", "not recorded"
        )
    else:
        spectra = {
            "raw": application.get("H_MM_before_proton_cleaning"),
            "estimated": application.get("H_MM_estimated_proton"),
            "cleaned": application.get("H_MM_after_proton_cleaning"),
            "committed": application.get("H_MM_after_proton_cleaning_final_rf"),
        }
        numerical_closure = {
            "proposed_pre_rf_closure_passed": gate.get("proposed_pre_rf_closure_passed"),
            "proposed_pre_rf_closure_difference": gate.get("proposed_pre_rf_closure_difference"),
            "final_applied_closure_passed": gate.get("final_applied_closure_passed"),
            "final_applied_pre_rf_closure_difference": gate.get("final_applied_pre_rf_closure_difference"),
        }
        global_closure = application_diagnostics.get(
            "canonical_t_global_closure", "not recorded"
        )
        identity_host_closure = "not recorded"
    return {
        "status": result.get("status") or diagnostics.get("status") or "not recorded",
        "method": result.get("method") or "not recorded",
        "canonical_t_binning": diagnostics.get("canonical_t_binning", "not recorded"),
        "shifted_t_consistency": diagnostics.get("cross_stage_t_consistency_summary", diagnostics.get("cross_stage_t_consistency", "not recorded")),
        "closure_mode": "identity_no_proton_cleaning" if identity_host else "proton_cleaned",
        "numerical_closure": numerical_closure,
        "global_closure": global_closure,
        "identity_host_closure": identity_host_closure,
        "lambda_gate_status": context.get("lambda_gate_status", gate.get("status", "not recorded")),
        "production_action": context.get("production_action", gate.get("production_action", application.get("production_action", "not recorded"))),
        "proton_cleaning_committed": context.get("proton_cleaning_committed", gate.get("proton_cleaning_committed", application.get("accepted", "not recorded"))),
        "host_state": host_state,
        "spectra": spectra,
    }


def proton_closure_summary_lines(qa):
    """Format the two already-recorded proton closure types without merging them."""
    payload = _mapping(qa)
    if payload.get("closure_mode") == "identity_no_proton_cleaning":
        identity = _mapping(payload.get("identity_host_closure"))
        transform = _mapping(identity.get("identity_transform_closure"))
        upstream = _mapping(identity.get("upstream_noRF_closure"))
        lines = [
            "Identity-host closure",
            "overall: {}".format(_pass_fail_label(identity.get("passed"))),
            "identity-transform: {}".format(_pass_fail_label(transform.get("passed"))),
            "upstream noRF: {}".format(_pass_fail_label(upstream.get("passed"))),
        ]
        for key, label in (("global_full", "global full"), ("global_cut", "global cut")):
            closure = _mapping(transform.get(key))
            if closure:
                lines.append("{}: {}".format(label, _pass_fail_label(closure.get("passed"))))
        return lines
    numerical = _mapping(payload.get("numerical_closure"))
    global_closure = _mapping(payload.get("global_closure"))
    lines = [
        "Numerical proton closure",
        "proposed: {} | raw - proton - cleaned = {}".format(
            _pass_fail_label(numerical.get("proposed_pre_rf_closure_passed")),
            _display_value(numerical.get("proposed_pre_rf_closure_difference")),
        ),
        "committed/applied: {} | raw - proton - cleaned = {}".format(
            _pass_fail_label(numerical.get("final_applied_closure_passed")),
            _display_value(numerical.get("final_applied_pre_rf_closure_difference")),
        ),
        "Canonical-|t| global closure",
    ]
    labels = {
        "raw": "raw",
        "proton_estimate": "proton estimate",
        "proton_cleaned_pre_rf": "cleaned pre-RF",
        "final_post_rf": "final",
    }
    if global_closure:
        for key, label in labels.items():
            closure = _mapping(global_closure.get(key))
            lines.append("{}: {}".format(label, _pass_fail_label(closure.get("passed"))))
    else:
        lines.append("not recorded")
    return lines


def proton_main_summary_lines(qa):
    """Build a compact host-state-aware main-PDF proton summary."""
    payload = _mapping(qa)
    lines = [
        "status: {}; method: {}".format(payload.get("status"), payload.get("method")),
        "host state: {}".format(payload.get("host_state")),
        "Lambda gate: {}".format(payload.get("lambda_gate_status")),
        "production action: {}; committed={}".format(
            payload.get("production_action"), payload.get("proton_cleaning_committed")
        ),
        "canonical-t provenance/consistency: {}".format(
            _display_value(payload.get("canonical_t_binning"))
        ),
        "shifted-t consistency: {}".format(
            _display_value(payload.get("shifted_t_consistency"))
        ),
    ]
    lines.extend(proton_closure_summary_lines(payload))
    if payload.get("closure_mode") == "identity_no_proton_cleaning":
        lines.append("No proton subtraction was applied.")
    lines.append("details are in the proton-debug supplement.")
    return tuple(lines)


def _render_proton_committed_mm_page(pdf_name, setting, qa, manifest):
    ROOT = _import_root()
    spectra = _mapping(qa.get("spectra"))
    ordered = (
        ("raw", "raw", 1),
        ("estimated", "estimated proton", 2),
        ("cleaned", "cleaned", 4),
        ("committed", "committed shifted MM", 6),
    )
    display = [
        (label, _clone_display(spectra.get(key), "H_proton_main_{}_{}".format(key, abs(id(qa)))), color)
        for key, label, color in ordered
    ]
    display = [(label, hist, color) for label, hist, color in display if hist is not None]
    if ROOT is None or not display:
        return _print_text_page(pdf_name, "proton.summary.committed_mm", _title(setting, "Proton", "committed shifted MM"), ("Final committed shifted MM overlay: not available", "Details remain in the proton-debug supplement."), manifest)
    identity_host = qa.get("closure_mode") == "identity_no_proton_cleaning"
    page_name = "identity-host committed MM" if identity_host else "committed shifted MM"
    canvas = ROOT.TCanvas("C_proton_main_committed_mm", _title(setting, "Proton", page_name), 1200, 800)
    try:
        draw_objects = [histogram for _label, histogram, _color in display]
        legend = ROOT.TLegend(0.64, 0.68, 0.90, 0.89)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        draw_objects.append(legend)
        host_annotation = ROOT.TPaveText(0.12, 0.82, 0.58, 0.91, "NDC")
        host_annotation.SetFillStyle(0)
        host_annotation.SetBorderSize(0)
        host_annotation.SetTextAlign(12)
        host_annotation.SetTextSize(0.022)
        host_annotation.AddText("host state: {}".format(qa.get("host_state")))
        host_annotation.AddText("proton-cleaning action: {}".format(qa.get("production_action")))
        draw_objects.append(host_annotation)
        for index, (label, histogram, color) in enumerate(display):
            histogram.SetLineColor(color)
            histogram.SetLineWidth(2 if label == "committed shifted MM" else 1)
            histogram.SetTitle("{};shifted missing mass;counts".format(_title(setting, "Proton", "committed shifted MM")))
            histogram.Draw("hist" if index == 0 else "hist same")
            legend.AddEntry(histogram, label, "l")
        legend.Draw()
        host_annotation.Draw()
        missing = [label for key, label, _color in ordered if spectra.get(key) is None]
        if missing:
            annotation = ROOT.TPaveText(0.12, 0.78, 0.58, 0.89, "NDC")
            annotation.SetFillStyle(0)
            annotation.SetBorderSize(0)
            annotation.SetTextAlign(12)
            annotation.SetTextSize(0.022)
            annotation.AddText("missing display states: {}".format(", ".join(missing)))
            draw_objects.append(annotation)
            annotation.Draw()
        canvas.Print(pdf_name)
        manifest.append({"page_id": "proton.summary.committed_mm", "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def render_proton_main_summary_pages(pdf_name, checkpoint, cleaning_result, cleaning_application, *, phase_a_display_context=None, page_manifest=None):
    """Emit concise main-PDF proton summaries from existing result/application objects."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    setting = _mapping(_mapping(checkpoint).get("setting"))
    qa = proton_main_qa_payload(cleaning_result, cleaning_application, phase_a_display_context)
    closure_lines = proton_main_summary_lines(qa)
    _print_text_page(
        pdf_name,
        "proton.summary.provenance_closure",
        _title(setting, "Proton", "provenance and closure"),
        closure_lines,
        manifest,
    )
    _render_proton_committed_mm_page(pdf_name, setting, qa, manifest)
    _print_text_page(pdf_name, "proton.summary.commitment", _title(setting, "Proton", "commitment summary"), (
        "Lambda gate: {}".format(qa["lambda_gate_status"]),
        "production action: {}; committed={}".format(qa["production_action"], qa["proton_cleaning_committed"]),
        "resulting host state: {}".format(qa["host_state"]),
        "presentation uses existing committed application objects only.",
    ), manifest)
    return manifest


__all__ = (
    "build_pdf_destinations",
    "build_pdf_route_manifest",
    "canonical_delta_frame_limits",
    "canonical_parent_lambda_summary",
    "close_diagnostic_pdf",
    "method_a_f_low_points",
    "method_a_f_low_frame_limits",
    "method_a_f_low_style",
    "method_a_plot_payload",
    "method_b_candidate_points",
    "method_b_candidate_frame_limits",
    "method_b_plot_payload",
    "method_b_regional_panels",
    "method_b_regional_frame_limits",
    "method_b_regional_rows",
    "method_b_shape_rows",
    "open_diagnostic_pdf",
    "phase_a_summary_payload",
    "proton_main_summary_lines",
    "proton_main_qa_payload",
    "proton_closure_summary_lines",
    "refinement_annotation_lines",
    "render_pion_hgcer_refinement_pages",
    "render_proton_main_summary_pages",
    "render_setting_warning_page",
    "setting_qa_summary_payload",
    "setting_qa_summary_lines",
    "setting_qa_warning_states",
    "unity_line_limits",
    "warning_payload",
)
