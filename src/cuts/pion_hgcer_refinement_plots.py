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
            "proton.summary",
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


def phase_a_summary_payload(checkpoint, phase_a):
    """Extract only stored Phase-A contract fields for a compact overview."""
    checkpoint = _mapping(checkpoint)
    phase = _mapping(phase_a)
    setting = _mapping(checkpoint.get("setting"))
    phase_checkpoint = _mapping(checkpoint.get("phase_a"))
    summary = _mapping(phase.get("summary"))
    if not summary:
        summary = _mapping(phase_checkpoint.get("summary"))
    return {
        "setting": dict(setting),
        "status": phase.get("status", summary.get("status", "unavailable")),
        "available": bool(phase.get("available", summary.get("available", False))),
        "reason": phase.get("reason", summary.get("reason")),
        "host_state": phase.get("host_state", phase_checkpoint.get("host_state")),
        "source_target_state": phase.get("source_target_state", phase_checkpoint.get("source_target_state")),
        "pion_closure_passed": bool(
            _mapping(phase.get("pion_closure")).get("passed", summary.get("pion_closure_passed", False))
        ),
        "host_closure_passed": bool(
            _mapping(phase.get("host_closure")).get("passed", summary.get("host_closure_passed", False))
        ),
        "canonical_t_edges": _edges(phase, "canonical_t_edges") or list(checkpoint.get("canonical_t_edges") or ()),
        "delta_edges": _edges(phase, "delta_edges") or list(checkpoint.get("delta_edges") or ()),
        "lambda_gate_status": phase.get("lambda_gate_status", summary.get("lambda_gate_status")),
        "lambda_gate_production_action": phase.get(
            "lambda_gate_production_action", summary.get("lambda_gate_production_action")
        ),
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
    """Return faithful Method-A support and f_low cells without re-counting."""
    payload = _mapping(method_a)
    checkpoint_method = _mapping(_mapping(checkpoint).get("method_a"))
    summary = _mapping(payload.get("summary")) or _mapping(checkpoint_method.get("summary"))
    cells = [dict(_mapping(cell)) for cell in (payload.get("cells") or checkpoint_method.get("cells") or ())]
    cells.sort(key=_cell_key)
    return {
        "status": payload.get("status", checkpoint_method.get("status", "unavailable")),
        "available": bool(payload.get("available", checkpoint_method.get("available", False))),
        "reason": payload.get("reason", checkpoint_method.get("reason")),
        "t_edges": _edges(payload, "t_edges") or list(_mapping(checkpoint).get("canonical_t_edges") or ()),
        "delta_edges": _edges(payload, "delta_edges") or list(_mapping(checkpoint).get("delta_edges") or ()),
        "support_counts": dict(summary.get("support_counts") or {}),
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
    return {
        "status": payload.get("status", checkpoint_method.get("status", "unavailable")),
        "available": bool(payload.get("available", checkpoint_method.get("available", False))),
        "reason": payload.get("reason", checkpoint_method.get("reason")),
        "t_edges": _edges(payload, "t_edges") or list(_mapping(checkpoint).get("canonical_t_edges") or ()),
        "delta_edges": _edges(payload, "delta_edges") or list(_mapping(checkpoint).get("delta_edges") or ()),
        "method_status_counts": dict(summary.get("method_B_status_counts") or {}),
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
                "region_name": value.get("region_name"),
                "parent_relative_ratio": value.get("parent_relative_ratio"),
                "parent_relative_sigma": value.get("parent_relative_sigma"),
                "parent_relative_status": value.get("parent_relative_status"),
            })
    return rows


def method_b_shape_rows(payload):
    """Extract stored shape metrics and their categorical state without thresholds."""
    return [
        {
            "t_index": _mapping(cell).get("t_index"),
            "delta_index": _mapping(cell).get("delta_index"),
            "shape_chi2_ndf": _mapping(cell).get("shape_chi2_ndf"),
            "shape_max_abs_pull": _mapping(cell).get("shape_max_abs_pull"),
            "shape_status": _mapping(cell).get("shape_status"),
        }
        for cell in _mapping(payload).get("cells") or ()
    ]


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


def _draw_map_page(ROOT, pdf_name, page_id, title, payload, *, value_key, category_codes=None, manifest=None):
    t_edges = list(payload.get("t_edges") or ())
    delta_edges = list(payload.get("delta_edges") or ())
    cells = list(payload.get("cells") or ())
    if len(t_edges) < 2 or len(delta_edges) < 2:
        return _print_text_page(pdf_name, page_id, title, ("unavailable: canonical t/delta edges not recorded",), manifest)
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
        canvas.Print(pdf_name)
        manifest.append({"page_id": page_id, "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def _draw_shape_quality_page(ROOT, pdf_name, title, payload, manifest):
    """Render stored chi2/ndf, maximum pull, and shape category side by side."""
    t_edges = list(payload.get("t_edges") or ())
    delta_edges = list(payload.get("delta_edges") or ())
    if len(t_edges) < 2 or len(delta_edges) < 2:
        return _print_text_page(pdf_name, "hgcer.method_b.shape_quality", title, ("unavailable: canonical t/delta edges not recorded",), manifest)
    from array import array
    suffix = abs(id(payload.get("cells")))
    maps = {
        "chi2": ROOT.TH2D("H_hgcer_method_b_shape_chi2_{}".format(suffix), "stored shape chi2/ndf;SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
        "pull": ROOT.TH2D("H_hgcer_method_b_shape_pull_{}".format(suffix), "stored maximum |pull|;SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
        "status": ROOT.TH2D("H_hgcer_method_b_shape_status_{}".format(suffix), "shape status (0 unavailable, 1 marginal, 2 good, 3 poor);SHMS delta;|t| [GeV^2]", len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges)),
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
    canvas = ROOT.TCanvas("C_hgcer_method_b_shape_quality", title, 1700, 650)
    try:
        canvas.Divide(3, 1)
        for pad, key in enumerate(("chi2", "pull", "status"), start=1):
            canvas.cd(pad)
            maps[key].Draw("colz text")
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.shape_quality", "scope": "setting", "authoritative": False})
        return True
    finally:
        canvas.Close()


def _render_phase_a_page(pdf_name, checkpoint, phase_a, manifest):
    payload = phase_a_summary_payload(checkpoint, phase_a)
    setting = payload["setting"]
    lines = (
        "Phase A status: {} (available={})".format(payload["status"], payload["available"]),
        "host state: {}; source target: {}".format(payload["host_state"] or "not recorded", payload["source_target_state"] or "not recorded"),
        "pion closure: {}; host closure: {}".format(payload["pion_closure_passed"], payload["host_closure_passed"]),
        "canonical |t| edges: {}".format(payload["canonical_t_edges"] or "not recorded"),
        "delta edges: {}".format(payload["delta_edges"] or "not recorded"),
        "proton/Lambda gate: {} / {}".format(payload["lambda_gate_status"] or "not recorded", payload["lambda_gate_production_action"] or "not recorded"),
        "contract fingerprint: {}; coordinate fingerprint: {}".format(payload["contract_fingerprint"], payload["coordinate_fingerprint"]),
        "non_authoritative=true; production_objects_mutated=false; refinement_applied=false",
        "reason: {}".format(payload["reason"] or "none"),
    )
    return _print_text_page(pdf_name, "hgcer.phase_a.summary", _title(setting, "Phase A", "HGCer Event/Population Contract"), lines, manifest)


def _render_method_a_pages(pdf_name, checkpoint, method_a, manifest):
    payload = method_a_plot_payload(method_a, checkpoint)
    setting = _mapping(_mapping(checkpoint).get("setting"))
    category_codes = {"unsupported": 0.0, "marginal": 1.0, "supported": 2.0}
    _draw_map_page(ROOT=_import_root(), pdf_name=pdf_name, page_id="hgcer.method_a.support", title=_title(setting, "Method A", "support map"), payload=payload, value_key="support_class", category_codes=category_codes, manifest=manifest) if _import_root() is not None else None
    ROOT = _import_root()
    if ROOT is None:
        return False
    canvas = ROOT.TCanvas("C_hgcer_method_a_f_low", _title(setting, "Method A", "f_low by canonical t"), 1200, 800)
    try:
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            graph = ROOT.TGraphErrors()
            graph.SetTitle("{} t{};SHMS delta;f_low".format(_title(setting, "Method A", "f_low"), t_index + 1))
            point = 0
            for point_data in method_a_f_low_points(payload):
                if int(point_data.get("t_index", -1)) != t_index:
                    continue
                value = point_data["f_low"]
                graph.SetPoint(point, point_data["delta_center"], value)
                graph.SetPointError(point, 0.0, max(value - point_data["f_low_low"], point_data["f_low_high"] - value))
                point += 1
            if point:
                graph.SetMarkerStyle(24)
                graph.Draw("AP")
                plotted = True
        if not plotted:
            canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No supported or marginal f_low cells")
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
    status_codes = {"unavailable": 0.0, "marginal": 1.0, "shape_inconsistent": 2.0, "internally_inconsistent": 3.0, "available": 4.0}
    _draw_map_page(ROOT, pdf_name, "hgcer.method_b.status", _title(setting, "Method B", "cell status"), payload, value_key="method_B_status", category_codes=status_codes, manifest=manifest)

    canvas = ROOT.TCanvas("C_hgcer_method_b_candidate", _title(setting, "Method B", "candidate L_B"), 1200, 800)
    try:
        t_indices = sorted({int(cell.get("t_index")) for cell in payload["cells"] if _finite(cell.get("t_index")) is not None})
        canvas.Divide(max(1, min(len(t_indices), 4)), 1)
        plotted = False
        for pad, t_index in enumerate(t_indices[:4], start=1):
            canvas.cd(pad)
            unity = ROOT.TLine(0.0, 1.0, 1.0, 1.0)
            graph = ROOT.TGraphErrors()
            graph.SetTitle("{} t{};SHMS delta;L_B candidate".format(_title(setting, "Method B", "candidate"), t_index + 1))
            point = 0
            maximum_delta = 1.0
            for point_data in method_b_candidate_points(payload):
                if int(point_data.get("t_index", -1)) != t_index:
                    continue
                value, sigma = point_data["candidate_L_B"], point_data["candidate_L_B_uncertainty"]
                center = point_data["delta_center"]
                maximum_delta = max(maximum_delta, center)
                graph.SetPoint(point, center, value)
                graph.SetPointError(point, 0.0, sigma)
                point += 1
            if point:
                graph.SetMarkerStyle(20)
                graph.Draw("AP")
                unity.SetX2(maximum_delta)
                unity.SetLineStyle(2)
                unity.Draw("same")
                plotted = True
        if not plotted:
            canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No available multi-region candidate cells")
        canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.candidate", "scope": "setting", "authoritative": False})
    finally:
        canvas.Close()

    regional_rows = method_b_regional_rows(payload)
    regional_canvas = ROOT.TCanvas("C_hgcer_method_b_regional", _title(setting, "Method B", "regional closure"), 1500, 800)
    try:
        regions = ("pi_n", "pi_sidis", "pi_delta_high")
        regional_canvas.Divide(3, 1)
        plotted = False
        for pad, region_name in enumerate(regions, start=1):
            regional_canvas.cd(pad)
            graph = ROOT.TGraphErrors()
            graph.SetTitle("{} {};SHMS delta;Qtilde".format(_title(setting, "Method B", "regional closure"), region_name))
            point = 0
            for row in regional_rows:
                if row.get("region_name") != region_name or str(row.get("parent_relative_status") or "").lower() != "available":
                    continue
                cell = next((entry for entry in payload["cells"] if entry.get("t_index") == row.get("t_index") and entry.get("delta_index") == row.get("delta_index")), {})
                low, high = _finite(cell.get("delta_low")), _finite(cell.get("delta_high"))
                value, sigma = _finite(row.get("parent_relative_ratio")), _finite(row.get("parent_relative_sigma"))
                if None in (low, high, value, sigma):
                    continue
                graph.SetPoint(point, 0.5 * (low + high), value)
                graph.SetPointError(point, 0.0, sigma)
                point += 1
            if point:
                graph.SetMarkerStyle(20)
                graph.Draw("AP")
                plotted = True
            line = ROOT.TLine(0.0, 1.0, 1.0, 1.0)
            line.SetLineStyle(2)
            line.Draw("same")
        if not plotted:
            regional_canvas.cd(1)
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "No stored available regional rows")
        regional_canvas.Print(pdf_name)
        manifest.append({"page_id": "hgcer.method_b.regional_closure", "scope": "setting", "authoritative": False})
    finally:
        regional_canvas.Close()

    _draw_shape_quality_page(ROOT, pdf_name, _title(setting, "Method B", "shape quality"), payload, manifest)
    return True


def render_pion_hgcer_refinement_pages(pdf_name, checkpoint, *, phase_a=None, method_a=None, method_b=None, page_manifest=None):
    """Append all new detached Phase-A/Method-A/Method-B pages to an open PDF."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    if _import_root() is None:
        return manifest
    _render_phase_a_page(pdf_name, checkpoint, phase_a, manifest)
    _render_method_a_pages(pdf_name, checkpoint, method_a, manifest)
    _render_method_b_pages(pdf_name, checkpoint, method_b, manifest)
    return manifest


def render_setting_warning_page(pdf_name, checkpoint, *, phase_a=None, method_a=None, method_b=None, part2=None, page_manifest=None, close_pdf=False):
    """Append the final non-OK state summary to an already-open main PDF."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    warnings = warning_payload(phase_a, method_a, method_b, part2)
    setting = _mapping(_mapping(checkpoint).get("setting"))
    lines = ["No outstanding non-OK detached diagnostic states."] if not warnings else []
    for entry in warnings:
        lines.append("{} | status={} | reason={} | production impact={}".format(entry.get("scope"), entry.get("status"), entry.get("reason"), entry.get("production_impact")))
    target_pdf = "{}".format(pdf_name) + ")" if close_pdf else pdf_name
    _print_text_page(target_pdf, "qa.setting_warnings", _title(setting, "Warnings", "detached diagnostic status"), lines, manifest)
    return manifest


def render_proton_main_summary_pages(pdf_name, checkpoint, cleaning_result, cleaning_application, *, page_manifest=None):
    """Emit concise main-PDF proton summaries from existing result/application objects."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = _mapping(cleaning_result)
    diagnostics = _mapping(result.get("diagnostics"))
    setting = _mapping(_mapping(checkpoint).get("setting"))
    status = result.get("status") or diagnostics.get("status") or "not recorded"
    support = diagnostics.get("supported_delta_bins", "not recorded")
    marginal = diagnostics.get("marginal_delta_bins", "not recorded")
    _print_text_page(pdf_name, "proton.summary.timing_fit_support", _title(setting, "Proton", "timing and fit support"), (
        "status: {}; method: {}".format(status, result.get("method") or "not recorded"),
        "supported delta bins: {}; marginal delta bins: {}".format(support, marginal),
        "details are in the proton-debug supplement.",
    ), manifest)
    application = _mapping(cleaning_application)
    _print_text_page(pdf_name, "proton.summary.application", _title(setting, "Proton", "cleaning application"), (
        "accepted: {}; production action: {}".format(application.get("accepted", "not recorded"), application.get("production_action", "not recorded")),
        "raw / estimated proton / cleaned spectra are retained in the proton-debug supplement.",
        "presentation uses existing committed application objects only.",
    ), manifest)
    gate = _mapping(diagnostics.get("lambda_preservation_gate"))
    _print_text_page(pdf_name, "proton.summary.commitment", _title(setting, "Proton", "commitment summary"), (
        "Lambda gate: {}".format(gate.get("status") or "not recorded"),
        "production action: {}".format(gate.get("production_action") or "not recorded"),
        "host/closure state: {}".format(result.get("source_target_state") or diagnostics.get("source_target_state") or "not recorded"),
    ), manifest)
    return manifest


__all__ = (
    "build_pdf_destinations",
    "build_pdf_route_manifest",
    "close_diagnostic_pdf",
    "method_a_f_low_points",
    "method_a_plot_payload",
    "method_b_candidate_points",
    "method_b_plot_payload",
    "method_b_regional_rows",
    "method_b_shape_rows",
    "open_diagnostic_pdf",
    "phase_a_summary_payload",
    "render_pion_hgcer_refinement_pages",
    "render_proton_main_summary_pages",
    "render_setting_warning_page",
    "warning_payload",
)
