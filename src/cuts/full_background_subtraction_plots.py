"""Detached D.6 procedure pages for the full background-subtraction PDF.

This module is presentation-only.  It receives already-built proton-cleaning
objects, clones only what it draws, and never rebuilds a fit, event lookup, or
weight.  Later procedure phases can append pages between the explicit PDF
open/close helpers without inheriting the technical diagnostic-PDF lifecycle.
"""

from __future__ import annotations

import math
import os
from collections.abc import Mapping, Sequence


D6_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d6/v1"
FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX = "_full-background-subtraction"

_TIMING_T_METHOD = "timing_t_event_weight"
_CTIME_AERO_METHOD = "ctime_aero_event_weight"


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
        histogram, "{};Aerogel NPE;Coincidence time [ns]".format(title)
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
        histogram, "{};delta [%];Aerogel NPE".format(title)
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
    delta_edges = list(presentation.get("delta_edges") or ())
    method = presentation.get("method")
    for group in tuple(presentation.get("per_t") or ()):
        group = _mapping(group)
        t_number = int(group.get("t_index", -1)) + 1
        if _mapping(group.get("raw_mm")).get("available"):
            if _render_raw_mm_page(ROOT, pdf_name, group):
                manifest.append({"page_id": "full_background.d6.raw_mm", "scope": "t{}".format(t_number), "authoritative": False})
            else:
                result["failures"].append("D.6 raw MM page unavailable for t{}".format(t_number))
        else:
            result["failures"].append("D.6 raw MM input unavailable for t{}".format(t_number))

        if _mapping(group.get("pid")).get("available"):
            rendered_pid = (
                _render_timing_t_pid_page(ROOT, pdf_name, group, delta_edges)
                if method == _TIMING_T_METHOD
                else _render_ctime_aero_pid_page(ROOT, pdf_name, group)
                if method == _CTIME_AERO_METHOD
                else False
            )
            if rendered_pid:
                manifest.append({"page_id": "full_background.d6.proton_pid", "scope": "t{}".format(t_number), "authoritative": False})
            else:
                result["failures"].append("D.6 proton-identification page unavailable for t{}".format(t_number))
        else:
            result["failures"].append("D.6 proton-identification input unavailable for t{}".format(t_number))

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
                result["failures"].append("D.6 proton-weight page unavailable for t{}".format(t_number))
        else:
            result["failures"].append("D.6 proton-weight input unavailable for t{}".format(t_number))
    return result


__all__ = (
    "D6_PRESENTATION_SCHEMA_VERSION",
    "FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX",
    "build_full_background_subtraction_d6_payload",
    "close_full_background_subtraction_pdf",
    "full_background_subtraction_pdf_path",
    "open_full_background_subtraction_pdf",
    "render_full_background_subtraction_d6_pages",
)
