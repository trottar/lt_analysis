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

from canonical_binning import find_canonical_bin


D6_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d6/v1"
D7_PRESENTATION_SCHEMA_VERSION = "full_background_subtraction_d7/v1"
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
        for delta_index, display_rows in enumerate(rows_by_delta):
            canvas.cd(delta_index + 1)
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
            panel_title = "delta = [{:.3f}, {:.3f}] %".format(
                delta_edges[delta_index], delta_edges[delta_index + 1]
            )
            _set_histogram_title(
                raw_histogram,
                "{};Missing mass [GeV];Signed normalized yield".format(panel_title),
            )
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


def render_full_background_subtraction_procedure_pages(
    pdf_name, d6_payload, d7_payload, *, page_manifest=None
):
    """Append D.6 and D.7 page groups in complete canonical-t order."""
    manifest = page_manifest if isinstance(page_manifest, list) else []
    result = {"manifest": manifest, "failures": []}
    d6 = _mapping(d6_payload)
    d7 = _mapping(d7_payload)
    d6_available = bool(d6.get("available"))
    d7_available = bool(d7.get("available"))
    if not d6_available and not d7_available:
        result["failures"].extend((
            "D.6 procedure input unavailable: {}".format(d6.get("reason")),
            "D.7 procedure input unavailable: {}".format(d7.get("reason")),
        ))
        return result
    if not d6_available:
        result["failures"].append(
            "D.6 procedure input unavailable: {}".format(d6.get("reason"))
        )
    if not d7_available:
        result["failures"].append(
            "D.7 procedure input unavailable: {}".format(d7.get("reason"))
        )
    if d6_available and d7_available and list(d6.get("t_edges") or ()) != list(d7.get("t_edges") or ()):
        result["failures"].append("D.6/D.7 canonical t geometry mismatch")
        d7_available = False
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
    elif d7_available:
        for group in tuple(d7.get("per_t") or ()):
            _render_d7_t_pages(ROOT, pdf_name, d7, _mapping(group), manifest, result["failures"])
    return result


__all__ = (
    "D6_PRESENTATION_SCHEMA_VERSION",
    "D7_PRESENTATION_SCHEMA_VERSION",
    "FULL_BACKGROUND_SUBTRACTION_PDF_SUFFIX",
    "build_full_background_subtraction_d6_payload",
    "build_full_background_subtraction_d7_payload",
    "close_full_background_subtraction_pdf",
    "full_background_subtraction_pdf_path",
    "open_full_background_subtraction_pdf",
    "render_full_background_subtraction_d6_pages",
    "render_full_background_subtraction_d7_pages",
    "render_full_background_subtraction_procedure_pages",
)
