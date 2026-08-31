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
    edges = list(_mapping(payload).get("delta_edges") or ())
    return (edges[0], edges[-1]) if len(edges) >= 2 else (None, None)


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


def setting_qa_summary_payload(phase_a, method_a, method_b, part2=None, display_context=None, runtime_qa_context=None):
    """Collect existing setting-level coverage and QA state for terminal display."""
    phase = _mapping(phase_a)
    method_a_payload = method_a_plot_payload(method_a)
    method_b_payload = method_b_plot_payload(method_b)
    context = _mapping(display_context)
    runtime = _mapping(runtime_qa_context)
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
        "k_lambda_comparison": runtime.get("k_lambda_comparison", "not available"),
        "k_sigma0_availability": runtime.get("k_sigma0_availability", "not available"),
        "hgcer_diagnostic_availability": runtime.get("hgcer_diagnostic_availability", "not available"),
        "part2_status": _mapping(part2).get("status") or "not available",
    }


def _has_recorded_warning(value):
    if value is None or value in ("", "not available", "not recorded"):
        return False
    if isinstance(value, dict):
        return any(_has_recorded_warning(item) for item in value.values())
    return bool(value)


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
        histogram.Draw("colz text")
        label = _diagnostic_label(ROOT)
        label.Draw()
        if category_codes is not None:
            _category_legend(ROOT, category_labels or {}, category_codes).Draw()
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
        "veto": ROOT.TH2D("H_hgcer_method_b_shape_veto_{}".format(suffix), "stored candidate shape-pool veto (0 no, 1 veto);SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
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
        canvas.Divide(4, 1)
        for pad, key in enumerate(("chi2", "pull", "status", "veto"), start=1):
            canvas.cd(pad)
            maps[key].Draw("colz text")
            _diagnostic_label(ROOT).Draw()
            if key == "status":
                _category_legend(ROOT, _SHAPE_STATUS_LABELS, codes).Draw()
            elif key == "veto":
                _category_legend(ROOT, {"not_vetoed": "not vetoed", "shape_poor_veto": "shape-pool veto"}, {"not_vetoed": 0.0, "shape_poor_veto": 1.0}).Draw()
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
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            legend = ROOT.TLegend(0.64, 0.72, 0.90, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            first = True
            for support_class in ("supported", "marginal"):
                graph = ROOT.TGraphErrors()
                graph.SetTitle("{} t{};SHMS delta;f_low".format(_title(setting, "Method A", "f_low"), t_index + 1))
                point = 0
                for point_data in method_a_f_low_points(payload):
                    if int(point_data.get("t_index", -1)) != t_index or point_data.get("support_class") != support_class:
                        continue
                    value = point_data["f_low"]
                    graph.SetPoint(point, point_data["delta_center"], value)
                    graph.SetPointError(point, 0.0, max(value - point_data["f_low_low"], point_data["f_low_high"] - value))
                    point += 1
                if point:
                    graph.SetMarkerStyle(method_a_f_low_style(support_class))
                    graph.Draw("AP" if first else "P same")
                    legend.AddEntry(graph, support_class, "p")
                    first = False
                    plotted = True
            if not first:
                legend.Draw()
                _diagnostic_label(ROOT).Draw()
        if not plotted:
            canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No supported or marginal f_low cells")
            _diagnostic_label(ROOT).Draw()
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
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            delta_min, delta_max = unity_line_limits(payload)
            unity = ROOT.TLine(delta_min, 1.0, delta_max, 1.0) if delta_min is not None else None
            graph = ROOT.TGraphErrors()
            graph.SetTitle("{} t{};SHMS delta;L_B candidate".format(_title(setting, "Method B", "candidate"), t_index + 1))
            point = 0
            for point_data in method_b_candidate_points(payload):
                if int(point_data.get("t_index", -1)) != t_index:
                    continue
                value, sigma = point_data["candidate_L_B"], point_data["candidate_L_B_uncertainty"]
                center = point_data["delta_center"]
                graph.SetPoint(point, center, value)
                graph.SetPointError(point, 0.0, sigma)
                point += 1
            if point:
                graph.SetMarkerStyle(20)
                graph.Draw("AP")
                if unity is not None:
                    unity.SetLineStyle(2)
                    unity.Draw("same")
                _diagnostic_label(ROOT, candidate=True).Draw()
                plotted = True
        if not plotted:
            canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No available multi-region candidate cells")
            _diagnostic_label(ROOT, candidate=True).Draw()
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.candidate", "scope": "setting", "authoritative": False})
    finally:
        canvas.Close()

    regional_panels = method_b_regional_panels(payload)
    regional_canvas = ROOT.TCanvas("C_hgcer_method_b_regional", _title(setting, "Method B", "regional closure"), 1500, 800)
    try:
        regions = ("pi_n", "pi_sidis", "pi_delta_high")
        regional_canvas.Divide(max(1, min(len(regional_panels), 4)), 1)
        plotted = False
        for pad, panel in enumerate(regional_panels[:4], start=1):
            regional_canvas.cd(pad)
            legend = ROOT.TLegend(0.57, 0.68, 0.90, 0.89)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            first = True
            for region_name, marker_style in zip(regions, (20, 24, 25)):
                graph = ROOT.TGraphErrors()
                graph.SetTitle("{} t{};SHMS delta;Qtilde stored ratio".format(_title(setting, "Method B", "regional closure"), panel["t_index"] + 1))
                point = 0
                for row in panel["series"][region_name]:
                    graph.SetPoint(point, row["delta_center"], row["Qtilde"])
                    graph.SetPointError(point, 0.0, row["Qtilde_uncertainty"])
                    point += 1
                if point:
                    graph.SetMarkerStyle(marker_style)
                    graph.Draw("AP" if first else "P same")
                    legend.AddEntry(graph, region_name, "p")
                    first = False
                    plotted = True
            if not first:
                delta_min, delta_max = unity_line_limits(payload)
                if delta_min is not None:
                    line = ROOT.TLine(delta_min, 1.0, delta_max, 1.0)
                    line.SetLineStyle(2)
                    line.Draw("same")
                legend.Draw()
                _diagnostic_label(ROOT).Draw()
            else:
                ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No stored available regional rows for this t parent")
                _diagnostic_label(ROOT).Draw()
        if not plotted:
            regional_canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No stored available regional rows")
            _diagnostic_label(ROOT).Draw()
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
    """Append the final non-OK state summary to an already-open main PDF."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    warnings = warning_payload(phase_a, method_a, method_b, part2)
    qa = setting_qa_summary_payload(
        phase_a, method_a, method_b, part2,
        display_context=phase_a_display_context,
        runtime_qa_context=runtime_qa_context,
    )
    setting = _mapping(_mapping(checkpoint).get("setting"))
    lines = ["Setting-level coverage + QA warning summary"]
    lines.append("Method A coverage: {}".format(_display_value(qa["method_a_coverage"])))
    lines.append("Method B coverage: {}".format(_display_value(qa["method_b_coverage"])))
    if qa["method_b_other_coverage"]:
        lines.append("Method B other recorded states: {}".format(_display_value(qa["method_b_other_coverage"])))
    lines.append("Phase A coordinate/closure: status={}; pion={}; host={}".format(qa["phase_a_coordinate_status"], qa["phase_a_pion_closure"], qa["phase_a_host_closure"]))
    lines.append("Phase host: {}; Lambda gate: {}; proton action: {}; committed={}".format(qa["phase_host_state"], qa["lambda_gate_status"], qa["proton_action"], qa["proton_cleaning_committed"]))
    lines.append("Aerogel warnings: {}; proton warnings: {}".format(_display_value(qa["aerogel_warnings"]), _display_value(qa["proton_warnings"])))
    lines.append("K-Lambda comparison/fallback: {}".format(_display_value(qa["k_lambda_comparison"])))
    lines.append("K-Sigma0 availability: {}; HGCer diagnostic availability: {}; Part 2: {}".format(_display_value(qa["k_sigma0_availability"]), _display_value(qa["hgcer_diagnostic_availability"]), qa["part2_status"]))
    if not warnings and not _has_recorded_warning(qa["aerogel_warnings"]) and not _has_recorded_warning(qa["proton_warnings"]):
        lines.append("No outstanding setting-level QA warnings.")
    for entry in warnings:
        lines.append("{} | status={} | reason={} | production impact={}".format(entry.get("scope"), entry.get("status"), entry.get("reason"), entry.get("production_impact")))
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
    return {
        "status": result.get("status") or diagnostics.get("status") or "not recorded",
        "method": result.get("method") or "not recorded",
        "canonical_t_binning": diagnostics.get("canonical_t_binning", "not recorded"),
        "shifted_t_consistency": diagnostics.get("cross_stage_t_consistency_summary", diagnostics.get("cross_stage_t_consistency", "not recorded")),
        "numerical_closure": application_diagnostics.get("canonical_t_global_closure", diagnostics.get("timing_t_mm_diagnostics", "not recorded")),
        "global_closure": application_diagnostics.get("canonical_t_global_closure", "not recorded"),
        "lambda_gate_status": context.get("lambda_gate_status", gate.get("status", "not recorded")),
        "production_action": context.get("production_action", gate.get("production_action", application.get("production_action", "not recorded"))),
        "proton_cleaning_committed": context.get("proton_cleaning_committed", gate.get("proton_cleaning_committed", application.get("accepted", "not recorded"))),
        "host_state": context.get("host_state", application.get("host_state", result.get("source_target_state", "not recorded"))),
        "spectra": {
            "raw": application.get("H_MM_before_proton_cleaning"),
            "estimated": application.get("H_MM_estimated_proton"),
            "cleaned": application.get("H_MM_after_proton_cleaning"),
            "committed": application.get("H_MM_after_proton_cleaning_final_rf"),
        },
    }


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
    canvas = ROOT.TCanvas("C_proton_main_committed_mm", _title(setting, "Proton", "committed shifted MM"), 1200, 800)
    try:
        legend = ROOT.TLegend(0.64, 0.68, 0.90, 0.89)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        for index, (label, histogram, color) in enumerate(display):
            histogram.SetLineColor(color)
            histogram.SetLineWidth(2 if label == "committed shifted MM" else 1)
            histogram.SetTitle("{};shifted missing mass;counts".format(_title(setting, "Proton", "committed shifted MM")))
            histogram.Draw("hist" if index == 0 else "hist same")
            legend.AddEntry(histogram, label, "l")
        legend.Draw()
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
    _print_text_page(pdf_name, "proton.summary.provenance_closure", _title(setting, "Proton", "provenance and closure"), (
        "status: {}; method: {}".format(qa["status"], qa["method"]),
        "canonical-t provenance/consistency: {}".format(_display_value(qa["canonical_t_binning"])),
        "shifted-t consistency: {}".format(_display_value(qa["shifted_t_consistency"])),
        "numerical closure: {}".format(_display_value(qa["numerical_closure"])),
        "global closure: {}".format(_display_value(qa["global_closure"])),
        "details are in the proton-debug supplement.",
    ), manifest)
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
    "close_diagnostic_pdf",
    "method_a_f_low_points",
    "method_a_f_low_style",
    "method_a_plot_payload",
    "method_b_candidate_points",
    "method_b_plot_payload",
    "method_b_regional_panels",
    "method_b_regional_rows",
    "method_b_shape_rows",
    "open_diagnostic_pdf",
    "phase_a_summary_payload",
    "proton_main_qa_payload",
    "refinement_annotation_lines",
    "render_pion_hgcer_refinement_pages",
    "render_proton_main_summary_pages",
    "render_setting_warning_page",
    "setting_qa_summary_payload",
    "unity_line_limits",
    "warning_payload",
)
