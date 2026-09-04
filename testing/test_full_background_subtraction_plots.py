"""Focused D.6 through D.11 procedure-PDF contract tests."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest
from copy import deepcopy
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import full_background_subtraction_plots as plots


EXPECTED_FULL_BACKGROUND_PAGE_IDS = (
    "full_background.d6.raw_mm",
    "full_background.d6.proton_pid",
    "full_background.d6.proton_weight",
    "full_background.d7.proton_mm",
    "full_background.d7.proton_cleaned_mm",
    "full_background.d7.proton_delta_mm",
    "full_background.d8.pion_background_mm",
    "full_background.d8.pion_subtracted_mm",
    "full_background.d8.pion_delta_mm",
    "full_background.d9.hgcer_response",
    "full_background.d9.hgcer_delta_response",
    "full_background.d9.method_a_relative",
    "full_background.d10.method_b_mm_inputs",
    "full_background.d10.method_b_local_closure",
    "full_background.d10.method_b_relative",
    "full_background.d11.ab_overlay",
    "full_background.d11.ab_ratio_log",
    "full_background.d11.method_availability",
)


def _assert_full_background_manifest_contract(test_case, manifest):
    """D.12: every retained page is a canonical non-authoritative t page."""
    seen_scopes = set()
    current_scope = None
    page_positions_by_scope = {}
    for page in manifest:
        page_id = page.get("page_id")
        scope = page.get("scope")
        test_case.assertIn(page_id, EXPECTED_FULL_BACKGROUND_PAGE_IDS)
        test_case.assertRegex(scope, r"^t[1-9][0-9]*$")
        test_case.assertIs(page.get("authoritative"), False)
        for forbidden in ("global", "status", "provenance", "warning", "qa"):
            test_case.assertNotIn(forbidden, page_id.lower())
        if scope != current_scope:
            test_case.assertNotIn(scope, seen_scopes)
            seen_scopes.add(scope)
            current_scope = scope
        page_positions_by_scope.setdefault(scope, []).append(
            EXPECTED_FULL_BACKGROUND_PAGE_IDS.index(page_id)
        )
    for positions in page_positions_by_scope.values():
        test_case.assertEqual(positions, sorted(positions))
        test_case.assertEqual(len(positions), len(set(positions)))


def _snapshot_root_th1(histogram):
    """Capture all regular-bin and axis state of a detached D.6/D.7/D.9 TH1."""
    axis = histogram.GetXaxis()
    bin_count = histogram.GetNbinsX()
    return {
        "nbins_x": bin_count,
        "contents": tuple(
            histogram.GetBinContent(index)
            for index in range(1, bin_count + 1)
        ),
        "errors": tuple(
            histogram.GetBinError(index)
            for index in range(1, bin_count + 1)
        ),
        "x_min": axis.GetXmin(),
        "x_max": axis.GetXmax(),
        "x_edges": tuple(
            axis.GetBinLowEdge(index)
            for index in range(1, bin_count + 2)
        ),
    }


def _snapshot_root_th2(histogram):
    """Capture all regular-bin and axis state of a detached D.6/D.9 TH2."""
    x_axis = histogram.GetXaxis()
    y_axis = histogram.GetYaxis()
    x_bins = histogram.GetNbinsX()
    y_bins = histogram.GetNbinsY()
    return {
        "nbins_x": x_bins,
        "nbins_y": y_bins,
        "contents": tuple(
            histogram.GetBinContent(x_index, y_index)
            for x_index in range(1, x_bins + 1)
            for y_index in range(1, y_bins + 1)
        ),
        "errors": tuple(
            histogram.GetBinError(x_index, y_index)
            for x_index in range(1, x_bins + 1)
            for y_index in range(1, y_bins + 1)
        ),
        "x_min": x_axis.GetXmin(),
        "x_max": x_axis.GetXmax(),
        "y_min": y_axis.GetXmin(),
        "y_max": y_axis.GetXmax(),
        "x_edges": tuple(
            x_axis.GetBinLowEdge(index)
            for index in range(1, x_bins + 2)
        ),
        "y_edges": tuple(
            y_axis.GetBinLowEdge(index)
            for index in range(1, y_bins + 2)
        ),
    }


class _Histogram:
    def __init__(self, label):
        self.label = label
        self.directory = object()
        self.title = "original"

    def Clone(self, name):
        clone = type(self)("{}:{}".format(self.label, name))
        clone.title = self.title
        return clone

    def SetDirectory(self, directory):
        self.directory = directory

    def SetTitle(self, title):
        self.title = title


class _BinnedHistogram(_Histogram):
    created = []

    def __init__(self, label, contents=()):
        super().__init__(label)
        self.contents = list(contents)
        self.display_minimum = None
        self.display_maximum = None

    def Clone(self, name):
        clone = type(self)("{}:{}".format(self.label, name), self.contents)
        type(self).created.append(clone)
        return clone

    def GetNbinsX(self):
        return len(self.contents)

    def GetBinContent(self, bin_index):
        return self.contents[int(bin_index) - 1]

    def Reset(self, _option):
        self.contents = [0.0] * len(self.contents)

    def Fill(self, value, weight):
        index = min(max(int(float(value)), 0), len(self.contents) - 1)
        self.contents[index] += float(weight)

    def SetMinimum(self, value):
        self.display_minimum = float(value)

    def SetMaximum(self, value):
        self.display_maximum = float(value)

    def SetLineColor(self, _color):
        return None

    def SetLineWidth(self, _width):
        return None

    def SetStats(self, _value):
        return None

    def Draw(self, _option):
        return None


class _BinnedHistogram2D(_Histogram):
    created = []

    def __init__(self, label, x_bins=2, y_bins=2):
        super().__init__(label)
        self.x_bins = int(x_bins)
        self.y_bins = int(y_bins)
        self.contents = [[0.0] * self.y_bins for _unused in range(self.x_bins)]
        self.display_minimum = None
        self.display_maximum = None

    def Clone(self, name):
        clone = type(self)("{}:{}".format(self.label, name), self.x_bins, self.y_bins)
        clone.contents = [list(row) for row in self.contents]
        type(self).created.append(clone)
        return clone

    def GetNbinsX(self):
        return self.x_bins

    def GetNbinsY(self):
        return self.y_bins

    def GetBinContent(self, x_index, y_index):
        return self.contents[int(x_index) - 1][int(y_index) - 1]

    def Reset(self, _option):
        self.contents = [[0.0] * self.y_bins for _unused in range(self.x_bins)]

    def Fill(self, x_value, y_value, weight):
        x_index = min(max(int(float(x_value)), 0), self.x_bins - 1)
        y_index = min(max(int(float(y_value)), 0), self.y_bins - 1)
        self.contents[x_index][y_index] += float(weight)

    def SetMinimum(self, value):
        self.display_minimum = float(value)

    def SetMaximum(self, value):
        self.display_maximum = float(value)

    def SetStats(self, _value):
        return None

    def Draw(self, _option):
        return None


class _ProjectionAxis:
    def __init__(self, bin_count):
        self._bin_count = int(bin_count)

    def FindBin(self, value):
        value = float(value)
        if value < 0.0:
            return 0
        if value >= self._bin_count:
            return self._bin_count + 1
        return int(value) + 1


class _ProjectionHistogram(_BinnedHistogram):
    def __init__(self, label, contents=(), errors=None):
        super().__init__(label, contents)
        self._underflow = 0.0
        self._overflow = 0.0
        self._variances = [float(value) ** 2 for value in (errors or [0.0] * len(self.contents))]
        self._underflow_variance = 0.0
        self._overflow_variance = 0.0
        self._axis = _ProjectionAxis(len(self.contents))

    def Clone(self, name):
        clone = type(self)(
            "{}:{}".format(self.label, name),
            self.contents,
            [variance ** 0.5 for variance in self._variances],
        )
        clone._underflow = self._underflow
        clone._overflow = self._overflow
        clone._underflow_variance = self._underflow_variance
        clone._overflow_variance = self._overflow_variance
        type(self).created.append(clone)
        return clone

    def GetXaxis(self):
        return self._axis

    def GetBinContent(self, bin_index):
        bin_index = int(bin_index)
        if bin_index == 0:
            return self._underflow
        if bin_index == len(self.contents) + 1:
            return self._overflow
        return super().GetBinContent(bin_index)

    def GetBinError(self, bin_index):
        bin_index = int(bin_index)
        if bin_index == 0:
            return self._underflow_variance ** 0.5
        if bin_index == len(self.contents) + 1:
            return self._overflow_variance ** 0.5
        return self._variances[bin_index - 1] ** 0.5

    def Reset(self, _option):
        super().Reset(_option)
        self._underflow = self._overflow = 0.0
        self._variances = [0.0] * len(self.contents)
        self._underflow_variance = self._overflow_variance = 0.0

    def Fill(self, value, weight):
        bin_index = self._axis.FindBin(value)
        weight = float(weight)
        if bin_index == 0:
            self._underflow += weight
            self._underflow_variance += weight * weight
        elif bin_index == len(self.contents) + 1:
            self._overflow += weight
            self._overflow_variance += weight * weight
        else:
            self.contents[bin_index - 1] += weight
            self._variances[bin_index - 1] += weight * weight


class _FakeCanvas:
    def __init__(self, *_args):
        return None

    def Divide(self, _columns, _rows):
        return None

    def cd(self, *_args):
        return self

    def Print(self, _name):
        return None

    def Close(self):
        return None


class _FakeLegend:
    def __init__(self, *coordinates):
        self.coordinates = tuple(float(value) for value in coordinates[:4])
        self.entries = []

    def SetBorderSize(self, _size):
        return None

    def SetFillStyle(self, _style):
        return None

    def AddEntry(self, _object, label, option):
        self.entries.append((str(label), str(option)))

    def Draw(self):
        return None


class _FakePaveText:
    def __init__(self, root, *coordinates):
        self._root = root
        self.coordinates = tuple(float(value) for value in coordinates[:4])
        self._text = []

    def SetFillStyle(self, _style):
        return None

    def SetBorderSize(self, _size):
        return None

    def SetTextAlign(self, _alignment):
        return None

    def SetTextSize(self, _size):
        return None

    def AddText(self, text):
        self._text.append(str(text))

    def Draw(self):
        self._root.drawn_text.append(tuple(self._text))


class _FakeROOT:
    kBlack = 1
    kBlue = 4
    kRed = 2

    def __init__(self):
        self.drawn_text = []
        self.legends = []
        self.pave_texts = []

    def TCanvas(self, *_args):
        return _FakeCanvas(*_args)

    def TLegend(self, *_args):
        legend = _FakeLegend(*_args)
        self.legends.append(legend)
        return legend

    def TPaveText(self, *_args):
        pave_text = _FakePaveText(self, *_args)
        self.pave_texts.append(pave_text)
        return pave_text


class _D10Axis:
    def __init__(self, edges):
        self.edges = tuple(float(edge) for edge in edges)

    def GetBinLowEdge(self, index):
        index = int(index)
        if index <= 0:
            return self.edges[0]
        if index >= len(self.edges):
            return self.edges[-1]
        return self.edges[index - 1]


class _D10DisplayHistogram:
    def __init__(self, name, title, bin_count, edges):
        self.name = str(name)
        self.title = str(title)
        self.edges = tuple(float(edge) for edge in edges)
        self.contents = [0.0] * (int(bin_count) + 2)
        self.errors = [0.0] * (int(bin_count) + 2)
        self.directory = object()
        self.display_minimum = None
        self.display_maximum = None
        self.root = None

    def SetDirectory(self, directory):
        self.directory = directory

    def Sumw2(self):
        return None

    def SetTitle(self, title):
        self.title = str(title)

    def SetMinimum(self, value):
        self.display_minimum = float(value)

    def SetMaximum(self, value):
        self.display_maximum = float(value)

    def SetStats(self, _value):
        return None

    def SetLineColor(self, _color):
        return None

    def SetLineWidth(self, _width):
        return None

    def Draw(self, option):
        if self.root is not None:
            self.root.drawn_histograms.append((self, str(option)))

    def SetBinContent(self, index, value):
        self.contents[int(index)] = float(value)

    def SetBinError(self, index, value):
        self.errors[int(index)] = float(value)

    def GetBinContent(self, index):
        return self.contents[int(index)]

    def GetBinError(self, index):
        return self.errors[int(index)]

    def GetNbinsX(self):
        return len(self.edges) - 1

    def GetXaxis(self):
        return _D10Axis(self.edges)


class _D10ROOT:
    def __init__(self):
        self.histograms = []

    def TH1D(self, *args):
        histogram = _D10DisplayHistogram(*args)
        self.histograms.append(histogram)
        return histogram


class _D10Line:
    def __init__(self, root, x1, y1, x2, y2):
        self.root = root
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)

    def SetLineColor(self, _color):
        return None

    def SetLineStyle(self, _style):
        return None

    def SetLineWidth(self, _width):
        return None

    def Draw(self):
        self.root.drawn_lines.append(self)


class _D10GraphErrors:
    def __init__(self, root, point_count):
        self.root = root
        self.points = [(None, None)] * int(point_count)
        self.errors = [(None, None)] * int(point_count)

    def SetPoint(self, index, x_value, y_value):
        self.points[int(index)] = (float(x_value), float(y_value))

    def SetPointError(self, index, x_error, y_error):
        self.errors[int(index)] = (float(x_error), float(y_error))

    def SetMarkerColor(self, _color):
        return None

    def SetLineColor(self, _color):
        return None

    def SetMarkerStyle(self, _style):
        return None

    def SetMarkerSize(self, _size):
        return None

    def Draw(self, _option):
        self.root.drawn_graphs.append(self)


class _D11GraphAsymmErrors:
    def __init__(self, root, point_count):
        self.root = root
        self.points = [(None, None)] * int(point_count)
        self.errors = [(None, None, None, None)] * int(point_count)

    def SetPoint(self, index, x_value, y_value):
        self.points[int(index)] = (float(x_value), float(y_value))

    def SetPointError(self, index, ex_low, ex_high, ey_low, ey_high):
        self.errors[int(index)] = (
            float(ex_low), float(ex_high), float(ey_low), float(ey_high)
        )

    def SetMarkerColor(self, _color):
        return None

    def SetLineColor(self, _color):
        return None

    def SetMarkerStyle(self, _style):
        return None

    def Draw(self, _option):
        self.root.drawn_asymm_graphs.append(self)


class _D11Graph:
    def __init__(self, root, point_count):
        self.root = root
        self.points = [(None, None)] * int(point_count)

    def SetPoint(self, index, x_value, y_value):
        self.points[int(index)] = (float(x_value), float(y_value))

    def SetMarkerColor(self, _color):
        return None

    def SetMarkerStyle(self, _style):
        return None

    def Draw(self, _option):
        self.root.drawn_plain_graphs.append(self)


class _D11Box:
    def __init__(self, root, x_low, y_low, x_high, y_high):
        self.root = root
        self.x_low = float(x_low)
        self.y_low = float(y_low)
        self.x_high = float(x_high)
        self.y_high = float(y_high)

    def SetFillColor(self, _color):
        return None

    def SetLineColor(self, _color):
        return None

    def Draw(self):
        self.root.drawn_boxes.append(self)


class _D10RenderROOT(_FakeROOT):
    """Test-only capture ROOT for the detached D.10 page renderers."""

    kGray = 920

    def __init__(self):
        super().__init__()
        self.histograms = []
        self.drawn_histograms = []
        self.drawn_lines = []
        self.drawn_graphs = []
        self.drawn_asymm_graphs = []
        self.drawn_plain_graphs = []
        self.drawn_boxes = []

    def TH1D(self, *args):
        histogram = _D10DisplayHistogram(*args)
        histogram.root = self
        self.histograms.append(histogram)
        return histogram

    def TLine(self, *args):
        return _D10Line(self, *args)

    def TGraphErrors(self, point_count):
        return _D10GraphErrors(self, point_count)

    def TGraphAsymmErrors(self, point_count):
        return _D11GraphAsymmErrors(self, point_count)

    def TGraph(self, point_count):
        return _D11Graph(self, point_count)

    def TBox(self, *args):
        return _D11Box(self, *args)


class _D11ROOTCapture:
    """Proxy real PyROOT constructors while retaining D.11 draw objects."""

    def __init__(self, root):
        self._root = root
        self.asymmetric_graphs = []
        self.symmetric_graphs = []
        self.central_graphs = []
        self.lines = []
        self.boxes = []

    def __getattr__(self, name):
        return getattr(self._root, name)

    def TGraphAsymmErrors(self, point_count):
        graph = self._root.TGraphAsymmErrors(int(point_count))
        self.asymmetric_graphs.append(graph)
        return graph

    def TGraphErrors(self, point_count):
        graph = self._root.TGraphErrors(int(point_count))
        self.symmetric_graphs.append(graph)
        return graph

    def TGraph(self, point_count):
        graph = self._root.TGraph(int(point_count))
        self.central_graphs.append(graph)
        return graph

    def TLine(self, *args):
        line = self._root.TLine(*args)
        self.lines.append(line)
        return line

    def TBox(self, *args):
        box = self._root.TBox(*args)
        self.boxes.append(box)
        return box


def _d10_visible_text(root):
    pieces = [histogram.title for histogram in root.histograms]
    pieces.extend(
        label for legend in root.legends for label, _option in legend.entries
    )
    pieces.extend(text for block in root.drawn_text for text in block)
    return "\n".join(pieces)


def _d10_histogram(root, name):
    return next(histogram for histogram in root.histograms if histogram.name == name)


def _rectangles_overlap(left, right):
    """Return whether two NDC rectangles have a positive-area overlap."""
    left_x1, left_y1, left_x2, left_y2 = (float(value) for value in left)
    right_x1, right_y1, right_x2, right_y2 = (float(value) for value in right)
    return (
        max(left_x1, right_x1) < min(left_x2, right_x2)
        and max(left_y1, right_y1) < min(left_y2, right_y2)
    )


class _ForbiddenTree:
    def __iter__(self):
        raise AssertionError("D.7 must not traverse source trees")


def _products(raw_sources, proton_sources=None, cleaned_sources=None):
    if proton_sources is None:
        proton_sources = tuple(
            _Histogram("proton-t{}".format(index)) for index, _raw in enumerate(raw_sources)
        )
    if cleaned_sources is None:
        cleaned_sources = tuple(
            _Histogram("cleaned-t{}".format(index)) for index, _raw in enumerate(raw_sources)
        )
    return tuple(
        {
            "t_index": index,
            "t_edges": [float(index), float(index + 1)],
            "raw_targets": {"h_mm_nosub": raw},
            "proton_targets": {"h_mm_nosub": proton_sources[index]},
            "cleaned_targets_pre_rf": {"h_mm_nosub": cleaned_sources[index]},
            "final_targets": {"h_mm_nosub": _Histogram("forbidden-final-{}".format(index))},
        }
        for index, raw in enumerate(raw_sources)
    )


def _timing_result():
    return {
        "method": "timing_t_event_weight",
        "selected_timing_branch": "P_RF_Dist_Track",
        "delta_edges": [-10.0, 0.0, 10.0],
        "H_delta_timing_t_cells": [
            (_Histogram("timing-d0-t0"), _Histogram("timing-d0-t1")),
            (_Histogram("timing-d1-t0"), _Histogram("timing-d1-t1")),
        ],
        "diagnostics": {"aerogel_vs_t_validation": {"enabled": True}},
        "_aerogel_vs_t_root_payload": {"global": {"aero": object()}},
    }


def _timing_application(raw_sources):
    return {
        "canonical_t_products": _products(raw_sources),
        "H_proton_weight_vs_delta_t": _Histogram("weight-delta-t"),
    }


def _legacy_result():
    return {
        "method": "ctime_aero_event_weight",
        "selected_timing_branch": "CTime_ROC1",
        "delta_edges": [-10.0, 0.0, 10.0],
        "H_global_pid": _Histogram("ctime-aero"),
    }


def _legacy_application(raw_sources):
    return {
        "canonical_t_products": _products(raw_sources),
        "H_proton_weight_vs_delta_aero": _Histogram("weight-delta-aero"),
    }


def _d7_prepared_bundle():
    return {
        "prepared_sources": {
            "prompt": {
                "coefficient": 2.0,
                "entries": {
                    0: {"nommcuts": True, "adj_mm": 0.2, "adj_t": 0.2},
                    1: {"nommcuts": True, "adj_mm": 1.2, "adj_t": 1.2},
                    2: {"nommcuts": False, "adj_mm": 0.4, "adj_t": 0.4},
                    3: {"nommcuts": True, "adj_mm": 0.6, "adj_t": 0.6},
                },
            },
            "dummy": {
                "coefficient": -1.0,
                "entries": {0: {"nommcuts": True, "adj_mm": 0.3, "adj_t": 0.3}},
            },
        }
    }


def _d7_timing_result():
    result = _timing_result()
    result["_prepared_event_weight_lookup"] = {
        "prompt:0": {
            "delta_index": 0,
            "t_index": 0,
            "proton_weight": 0.25,
            "cleaned_factor": 0.75,
        },
        "prompt:1": {
            "delta_index": 1,
            "t_index": 1,
            "proton_weight": 0.50,
            "cleaned_factor": 0.50,
        },
        "prompt:3": {
            "delta_index": 9,
            "t_index": 0,
            "proton_weight": 0.25,
            "cleaned_factor": 0.75,
        },
        "dummy:0": {
            "delta_index": 0,
            "t_index": 0,
            "proton_weight": 0.50,
            "cleaned_factor": 0.50,
        },
    }
    return result


def _d8_final_application(label, *, weights=(0.0, 2.0, 3.0, 0.0), zero=False):
    before = _ProjectionHistogram("{}-before".format(label), (8.0, 9.0))
    template = _ProjectionHistogram(
        "{}-template".format(label),
        (0.0, 3.0) if not zero else (0.0, 0.0),
        (8.0 ** 0.5, 3.0) if not zero else (0.0, 0.0),
    )
    after = _ProjectionHistogram("{}-after".format(label), (8.0, 6.0) if not zero else (8.0, 9.0))
    return {
        "H_MM_nosub_before_pion_subtraction": before,
        "H_pion_subtraction_template_MM_nosub": template,
        "H_MM_nosub_after_pion_subtraction": after,
        "H_pion_control_model": _ProjectionHistogram("{}-model".format(label), (1.0, 1.0)),
        "weights": list(weights),
        "H_MM_before_pion_subtraction": _ProjectionHistogram("{}-forbidden-cut-before".format(label), (99.0, 99.0)),
        "H_pion_subtraction_template_MM": _ProjectionHistogram("{}-forbidden-cut-template".format(label), (99.0, 99.0)),
        "H_MM_after_pion_subtraction": _ProjectionHistogram("{}-forbidden-cut-after".format(label), (99.0, 99.0)),
    }


def _d8_fixture():
    final_t0 = _d8_final_application("final-t0")
    final_t1 = _d8_final_application("final-t1", weights=(0.0, 0.0, 0.0, 0.0), zero=True)
    proposed_t0 = _d8_final_application("proposed-t0", weights=(0.0, 9.0, 9.0, 0.0))
    parents = (
        {
            "t_bin_index": 0,
            "t_edges": [0.0, 1.0],
            "proposed_diagnostic_application_result": proposed_t0,
            "final_diagnostic_application_result": final_t0,
            "final_diagnostic_application_status": {"final_status": "applied_component"},
        },
        {
            "t_bin_index": 1,
            "t_edges": [1.0, 2.0],
            "proposed_diagnostic_application_result": _d8_final_application("proposed-t1"),
            "final_diagnostic_application_result": final_t1,
            "final_diagnostic_application_status": {"final_status": "zero"},
        },
    )
    cache = {
        "delta_edges": [-10.0, 0.0, 10.0],
        "by_t": (
            {
                "t_index": 0,
                "t_edges": [0.0, 1.0],
                "records": (
                    {"nommcuts": True, "delta_index": 0, "adj_MM": 0.2, "coefficient": 1.0, "ssdelta": 5.0},
                    {"nommcuts": True, "delta_index": 1, "adj_MM": 0.2, "coefficient": -1.0, "ssdelta": -5.0},
                    {"nommcuts": True, "delta_index": 0, "adj_MM": 1.2, "coefficient": 1.0, "ssdelta": -5.0},
                    {"nommcuts": False, "delta_index": 0, "adj_MM": 0.2, "coefficient": 100.0},
                ),
            },
            {
                "t_index": 1,
                "t_edges": [1.0, 2.0],
                "records": (
                    {"nommcuts": True, "delta_index": 0, "adj_MM": 0.2, "coefficient": 1.0},
                ),
            },
        ),
    }
    return parents, cache


def _d8_single_scale_final_fallback_fixture():
    """Provide one frozen final single-scale parent without any resolver."""
    scale = 2.5
    final = {
        "H_MM_nosub_before_pion_subtraction": _ProjectionHistogram(
            "single-scale-final-before", (12.0, 8.0)
        ),
        "H_pion_subtraction_template_MM_nosub": _ProjectionHistogram(
            "single-scale-final-template", (5.0, -2.5), (5.0, 2.5)
        ),
        "H_MM_nosub_after_pion_subtraction": _ProjectionHistogram(
            "single-scale-final-after", (7.0, 10.5)
        ),
        "H_pion_control_model": _ProjectionHistogram(
            "single-scale-final-model", (1.0, 1.0)
        ),
        "weights": [0.0, scale, scale, 0.0],
    }
    proposed = {
        "H_MM_nosub_before_pion_subtraction": _ProjectionHistogram(
            "proposal-before-sentinel", (99.0, 98.0)
        ),
        "H_pion_subtraction_template_MM_nosub": _ProjectionHistogram(
            "proposal-template-sentinel", (97.0, 96.0)
        ),
        "H_MM_nosub_after_pion_subtraction": _ProjectionHistogram(
            "proposal-after-sentinel", (95.0, 94.0)
        ),
        "H_pion_control_model": _ProjectionHistogram(
            "proposal-model-sentinel", (1.0, 1.0)
        ),
        "weights": [0.0, 9.0, 9.0, 0.0],
    }
    parents = (
        {
            "t_bin_index": 0,
            "t_edges": [0.0, 1.0],
            "proposed_diagnostic_application_result": proposed,
            "final_diagnostic_application_result": final,
            "final_diagnostic_application_status": {
                "final_status": "applied_fallback",
                "fallback_mode": "single_scale",
            },
        },
    )
    cache = {
        "delta_edges": [-10.0, 0.0, 10.0],
        "by_t": (
            {
                "t_index": 0,
                "t_edges": [0.0, 1.0],
                "records": (
                    {
                        "nommcuts": True,
                        "delta_index": 0,
                        "adj_MM": 0.2,
                        "coefficient": 2.0,
                    },
                    {
                        "nommcuts": True,
                        "delta_index": 1,
                        "adj_MM": 1.2,
                        "coefficient": -1.0,
                    },
                ),
            },
        ),
    }
    return parents, cache


def _d9_fixture():
    """Return realistic detached Part-1, Method-A, and D.2 display inputs."""
    t_edges = [0.0, 1.0, 2.0]
    delta_edges = [-10.0, 0.0, 10.0]
    coordinate_fingerprint = "d9-coordinate-fingerprint"
    method_a_fingerprint = "d9-method-a-fingerprint"
    records = {}
    for side, sign in (("kaon", 1.0), ("pion", -1.0)):
        side_records = []
        for t_index in range(2):
            for delta_index in range(2):
                side_records.append({
                    "side": side,
                    "canonical_t_index": t_index,
                    "delta_index": delta_index,
                    # These sentinels must not replace the frozen indices above.
                    "ssdelta": 7.5 if delta_index == 0 else -7.5,
                    "P_hgcer_npeSum": 0.5 + float(t_index + delta_index),
                    "diagnostic_weight": sign * (3.25 + t_index + delta_index),
                    "adj_t": 99.0,
                    "reconstructed_weight": 12345.0,
                })
        records[side] = tuple(side_records)
    histograms = {
        "H_hgcer_kaon_weighted": _Histogram("d9-kaon-weighted"),
        "H_hgcer_pion_weighted": _Histogram("d9-pion-weighted"),
        "H_hgcer_vs_delta_kaon_weighted": _Histogram("d9-kaon-delta-weighted"),
        "H_hgcer_vs_delta_pion_weighted": _Histogram("d9-pion-delta-weighted"),
        "H_hgcer_kaon_absolute": _Histogram("forbidden-absolute"),
        "H_hgcer_pion_absolute": _Histogram("forbidden-absolute"),
    }
    diagnostic = {
        "status": "available",
        "non_authoritative": True,
        "production_side_effect_free": True,
        "production_hgcer_pid_unchanged": True,
        "rf_restoration_applied": False,
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "coordinate_fingerprint": coordinate_fingerprint,
        "records": records,
        "histograms": histograms,
    }
    method_a = {
        "status": "available",
        "available": True,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "fingerprint": method_a_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "configuration": {
            "positive_response_threshold": 0.0,
            "low_response_upper_threshold": 2.0,
        },
    }
    cells = []
    candidate_specs = {
        (0, 0): ("available", 0.0, 0.0, 0.25),
        (0, 1): ("marginal", 0.8, 0.5, 1.3),
        (1, 0): ("unavailable", None, None, None),
        (1, 1): ("available", 1.2, 0.9, 1.7),
    }
    for t_index in range(2):
        for delta_index in range(2):
            status, candidate, low, high = candidate_specs[(t_index, delta_index)]
            cells.append({
                "t_index": t_index,
                "t_low": t_edges[t_index],
                "t_high": t_edges[t_index + 1],
                "delta_index": delta_index,
                "delta_low": delta_edges[delta_index],
                "delta_high": delta_edges[delta_index + 1],
                "method_a_comparison_candidate": candidate,
                "method_a_comparison_candidate_low": low,
                "method_a_comparison_candidate_high": high,
                "method_a_comparison_candidate_status": status,
                "method_a_comparison_candidate_reason": None if candidate is not None else "parent_unavailable",
                "support_class": "supported" if status == "available" else "marginal",
                "method_A_status": "available" if candidate is not None else "unavailable",
                "method_A_reason": None if candidate is not None else "support_insufficient",
                "forbidden_recalculation_sentinel": 919.0,
            })
    comparison = {
        "schema_version": "pion_hgcer_method_a_comparison/v1",
        "status": "available",
        "available": True,
        "non_authoritative": True,
        "method_b_numerical_dependency": False,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "method_a_fingerprint": method_a_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "canonical_t_edges": t_edges,
        "delta_edges": delta_edges,
        "cells": cells,
    }
    return diagnostic, method_a, comparison


def _d10_fixture():
    """Return frozen Phase-A, Method-B, and D.3 inputs with source sentinels."""
    t_edges = [0.0, 1.0, 2.0]
    delta_edges = [-10.0, 0.0, 10.0]
    mm_edges = [0.80, 0.91, 1.03, 1.17, 1.35]
    phase_fingerprint = "d10-phase-a-fingerprint"
    coordinate_fingerprint = "d10-coordinate-fingerprint"
    method_b_fingerprint = "d10-method-b-fingerprint"
    host_rows = []
    pion_rows = []
    for t_index in range(2):
        for delta_index in range(2):
            host_rows.append({
                "canonical_t_index": t_index,
                "delta_index": delta_index,
                "analysis_MM": 0.85 + 0.10 * float(t_index + delta_index),
                # These raw coordinates and factors must never be used to reassign/rebuild.
                "analysis_abs_t": 99.0,
                "SHMS_delta": 99.0,
                "proton_cleaning_factor": 99.0,
                "nommcuts": True,
                "host_state": "proton_cleaned",
                "signed_host_event_contribution": 3.0 if delta_index == 0 else -1.5,
            })
            pion_rows.append({
                "canonical_t_index": t_index,
                "delta_index": delta_index,
                "analysis_MM": 0.87 + 0.10 * float(t_index + delta_index),
                "analysis_abs_t": -99.0,
                "SHMS_delta": -99.0,
                "signed_source_coefficient": 99.0,
                "baseline_pion_weight_w0": 99.0,
                "nommcuts": True,
                "signed_baseline_event_contribution": -2.5 if delta_index == 0 else 1.25,
            })

    def bins_for(t_index, delta_index):
        rows = []
        for index in range(len(mm_edges) + 1):
            regular = 1 <= index < len(mm_edges)
            rows.append({
                "index": index,
                "mm_low": mm_edges[index - 1] if regular else None,
                "mm_high": mm_edges[index] if regular else None,
                "underflow": index == 0,
                "overflow": index == len(mm_edges),
                "in_allowed_pion_sensitive_domain": regular,
                "usable_for_shape": regular,
                "host_record_count": 1 if regular else 0,
                "host_yield": 7.0 + t_index + delta_index if regular else 0.0,
                "host_sumw2": 9.0 if regular else 0.0,
                "baseline_record_count": 1 if regular else 0,
                "baseline_yield": -2.0 - t_index if regular else 0.0,
                "baseline_sumw2": 4.0 if regular else 0.0,
                "residual": 9.0 + delta_index if regular else 0.0,
                "pull": 2.0 if regular else None,
            })
        return rows

    method_b_cells = []
    d3_cells = []
    statuses = (
        "available_multi_region", "single_region_only", "region_marginal", "unavailable",
    )
    for t_index in range(2):
        for delta_index in range(2):
            status = statuses[t_index * 2 + delta_index]
            method_b_cells.append({
                "t_index": t_index,
                "t_low": t_edges[t_index],
                "t_high": t_edges[t_index + 1],
                "delta_index": delta_index,
                "delta_low": delta_edges[delta_index],
                "delta_high": delta_edges[delta_index + 1],
                "host_state": "proton_cleaned",
                "mm_edges": mm_edges,
                "bins": bins_for(t_index, delta_index),
                "shape_status": "good",
                "shape_reason": None,
                # Page 15 must not use these raw values.
                "candidate_L_B": 9.0,
                "candidate_L_B_uncertainty": 8.0,
            })
            candidate = 1.25 if status == "available_multi_region" else None
            uncertainty = 0.15 if status == "available_multi_region" else None
            d3_cells.append({
                "t_index": t_index,
                "t_low": t_edges[t_index],
                "t_high": t_edges[t_index + 1],
                "delta_index": delta_index,
                "delta_low": delta_edges[delta_index],
                "delta_high": delta_edges[delta_index + 1],
                "method_b_comparison_candidate": candidate,
                "method_b_comparison_candidate_uncertainty": uncertainty,
                "method_b_comparison_candidate_status": status,
                "method_b_comparison_candidate_reason": None if candidate is not None else status,
                "method_B_status": "available" if candidate is not None else "unavailable",
                "method_B_reason": None if candidate is not None else status,
                "region_consistency_status": "region_consistent",
                "region_consistency_reason": None,
                "shape_status": "good",
                "shape_reason": None,
            })

    phase_a = {
        "schema_version": "pion_hgcer_event_contract/v1",
        "status": "available",
        "available": True,
        "contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "canonical_t_edges": t_edges,
        "delta_edges": delta_edges,
        "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF",
        "pion_records": pion_rows,
        "kaon_host_records": host_rows,
        "pion_closure": {"passed": True},
        "host_closure": {"passed": True},
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_restoration_applied": False,
    }
    method_b = {
        "schema_version": "pion_hgcer_method_b/v1",
        "status": "available",
        "available": True,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "interpolation_used": False,
        "phase_a_records_only": True,
        "method_a_numerical_dependency": False,
        "phase_a_contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF",
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "mm_binning": mm_edges,
        "mm_regions": [
            {
                "region_name": "pion_sensitive_low",
                "mm_low": 0.80,
                "mm_high": 1.10,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
                "available": True,
                "protected_signal_overlap": False,
                "reason": None,
            },
            {
                "region_name": "pion_sensitive_high",
                "mm_low": 1.23,
                "mm_high": 1.35,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
                "available": True,
                "protected_signal_overlap": False,
                "reason": None,
            },
        ],
        "protected_regions": [{
            "region_name": "KLambdaSigma0",
            "mm_low": 1.10,
            "mm_high": 1.23,
            "region_role": "protected_signal",
            "window_source": "method_b.fixed_lambda_sigma_protected",
            "mm_offset_applied": 0.0,
        }],
        "fingerprint": method_b_fingerprint,
        "cells": method_b_cells,
    }
    comparison = {
        "schema_version": "pion_hgcer_method_b_comparison/v1",
        "status": "available",
        "available": True,
        "non_authoritative": True,
        "method_a_numerical_dependency": False,
        "comparison_performed": False,
        "classification_performed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "phase_a_contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "method_b_fingerprint": method_b_fingerprint,
        "canonical_t_edges": t_edges,
        "delta_edges": delta_edges,
        "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF",
        "cells": d3_cells,
    }
    return phase_a, method_b, comparison


def _d11_fixture(t_edges=(0.0, 1.0, 2.0), delta_edges=(-12.0, -5.0, 1.0, 6.0, 14.0, 23.0)):
    """Return a detached Phase-D checkpoint with stored D.4 sentinels."""
    def canonical_fingerprint(value):
        return hashlib.sha256(json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False,
        ).encode("ascii")).hexdigest()

    t_edges = list(t_edges)
    delta_edges = list(delta_edges)
    source_fingerprint = "d11-phase-c-source-fingerprint"
    phase_fingerprint = "d11-phase-a-fingerprint"
    coordinate_fingerprint = "d11-coordinate-fingerprint"
    method_a_fingerprint_inputs = {"representation": "d11-test-method-a"}
    method_b_fingerprint_inputs = {"representation": "d11-test-method-b"}
    method_a_fingerprint = canonical_fingerprint(method_a_fingerprint_inputs)
    method_b_fingerprint = canonical_fingerprint(method_b_fingerprint_inputs)
    definitions = (
        {
            "availability": "both_comparable", "availability_reason": None,
            "a": (True, 2.0, 1.0, 3.0, "available"),
            "b": (True, 4.0, 0.5, "available_multi_region"),
            # Deliberately not derivable from the A/B sentinels.
            "ratio": 7.25, "log": -3.5,
        },
        {
            "availability": "both_present_not_comparable",
            "availability_reason": "method_a_comparison_candidate_nonpositive",
            "a": (True, 0.0, 0.0, 0.25, "marginal"),
            "b": (True, 1.1, 0.2, "available_multi_region"),
            "ratio": None, "log": None,
        },
        {
            "availability": "a_only", "availability_reason": None,
            "a": (True, 0.8, 0.4, 1.2, "available"),
            "b": (False, None, None, "unavailable"), "ratio": None, "log": None,
        },
        {
            "availability": "b_only", "availability_reason": None,
            "a": (False, None, None, None, "unavailable"),
            "b": (True, 0.9, 0.1, "available_multi_region"), "ratio": None, "log": None,
        },
        {
            "availability": "neither_available", "availability_reason": None,
            "a": (False, None, None, None, "unavailable"),
            "b": (False, None, None, "unavailable"), "ratio": None, "log": None,
        },
    )
    cells = []
    for t_index in range(len(t_edges) - 1):
        for delta_index in range(len(delta_edges) - 1):
            definition = definitions[delta_index % len(definitions)]
            a_present, a_candidate, a_low, a_high, a_status = definition["a"]
            b_present, b_candidate, b_uncertainty, b_status = definition["b"]
            cells.append({
                "t_index": t_index,
                "t_low": t_edges[t_index],
                "t_high": t_edges[t_index + 1],
                "delta_index": delta_index,
                "delta_low": delta_edges[delta_index],
                "delta_high": delta_edges[delta_index + 1],
                "method_a": {
                    "present": a_present,
                    "comparison_candidate": a_candidate,
                    "comparison_candidate_low": a_low,
                    "comparison_candidate_high": a_high,
                    "comparison_candidate_status": a_status,
                },
                "method_b": {
                    "present": b_present,
                    "comparison_candidate": b_candidate,
                    "comparison_candidate_uncertainty": b_uncertainty,
                    "comparison_candidate_status": b_status,
                },
                "comparison": {
                    "availability": definition["availability"],
                    "availability_reason": definition["availability_reason"],
                    "ratio_B_over_A": definition["ratio"],
                    "log_ratio_B_over_A": definition["log"],
                    "diagnostic_interval_relation": "disjoint",
                },
            })
    method_a = {
        "schema_version": "pion_hgcer_method_a_comparison/v1",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source_fingerprint,
        "phase_a_contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "source_method_a_payload_fingerprint": "d11-source-method-a-payload",
        "fingerprint_inputs": method_a_fingerprint_inputs,
        "fingerprint": method_a_fingerprint,
        "non_authoritative": True, "method_b_numerical_dependency": False,
        "comparison_performed": False, "classification_performed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    method_b = {
        "schema_version": "pion_hgcer_method_b_comparison/v1",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source_fingerprint,
        "phase_a_contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "source_method_b_payload_fingerprint": "d11-source-method-b-payload",
        "fingerprint_inputs": method_b_fingerprint_inputs,
        "fingerprint": method_b_fingerprint,
        "host_state": "proton_cleaned", "source_target_state": "post_proton_noRF",
        "non_authoritative": True, "method_a_numerical_dependency": False,
        "comparison_performed": False, "classification_performed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    comparison_fingerprint_inputs = {"representation": "d11-test-ab"}
    ab_fingerprint = canonical_fingerprint(comparison_fingerprint_inputs)
    comparison = {
        "schema_version": "pion_hgcer_ab_comparison/v1",
        "method": "non_authoritative_ab_comparison",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source_fingerprint,
        "phase_a_contract_fingerprint": phase_fingerprint,
        "coordinate_fingerprint": coordinate_fingerprint,
        "method_a_comparison_fingerprint": method_a_fingerprint,
        "method_b_comparison_fingerprint": method_b_fingerprint,
        "source_method_a_comparison_payload_fingerprint": canonical_fingerprint(method_a),
        "source_method_b_comparison_payload_fingerprint": canonical_fingerprint(method_b),
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "host_state": "proton_cleaned", "source_target_state": "post_proton_noRF",
        "cells": cells, "fingerprint_inputs": comparison_fingerprint_inputs,
        "fingerprint": ab_fingerprint,
        "non_authoritative": True, "comparison_performed": True,
        "classification_performed": True,
        "classification_scope": "availability_only_non_prescriptive",
        "decision_performed": False, "statistical_compatibility_claimed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    return {
        "schema_version": "pion_hgcer_phase_d_checkpoint/v1",
        "source_checkpoint_payload_fingerprint": source_fingerprint,
        "method_a_comparison": method_a,
        "method_b_comparison": method_b,
        "ab_comparison": comparison,
        "status": "available", "available": True,
        "non_authoritative": True, "comparison_performed": True,
        "classification_performed": True,
        "classification_scope": "availability_only_non_prescriptive",
        "decision_performed": False, "statistical_compatibility_claimed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }


def _d11_cumulative_payload(label, t_edges=(0.0, 1.0), delta_edges=(-10.0, 0.0, 10.0)):
    """Minimal already-built page group for cumulative D.11 gate tests."""
    return {
        "label": str(label),
        "available": True,
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "per_t": [{"t_index": 0}],
    }


def _d12_cumulative_payload(
    label, t_edges=(0.0, 1.0, 2.0), delta_edges=(-10.0, 0.0, 10.0)
):
    """Detached two-bin presentation groups for D.12 cumulative-order tests."""
    return {
        "label": str(label),
        "available": True,
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "per_t": [
            {"t_index": index}
            for index in range(len(t_edges) - 1)
        ],
    }


def _d12_record_phase_pages(page_ids, omissions):
    """Record detached page availability without invoking ROOT page builders."""
    def recorder(_root, _pdf_name, _presentation, group, manifest, failures):
        scope = "t{}".format(int(group["t_index"]) + 1)
        for page_id in page_ids:
            if (scope, page_id) in omissions:
                failures.append("{} unavailable for {}".format(page_id, scope))
                continue
            manifest.append({
                "page_id": page_id,
                "scope": scope,
                "authoritative": False,
            })
    return recorder


class FullBackgroundSubtractionD6Tests(unittest.TestCase):
    def test_pdf_path_is_deterministic_and_detached(self):
        cases = (
            (
                r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe.pdf",
                r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction.pdf",
            ),
            (
                "/analysis/Left_kaon_rand_sub_Q4p4W2p74_highe.pdf",
                "/analysis/Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction.pdf",
            ),
            (
                "/analysis/Left_kaon_rand_sub_Q4p4W2p74_highe",
                "/analysis/Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction.pdf",
            ),
        )
        for main, expected in cases:
            with self.subTest(main=main):
                path = plots.full_background_subtraction_pdf_path(main)
                self.assertEqual(path, expected)
                self.assertNotEqual(path, main)
                self.assertNotIn("hgcer-ab-comparison", path)
                self.assertNotEqual(
                    path,
                    main.rsplit(".pdf", 1)[0] + "_hgcer-ab-comparison.pdf",
                )

    def test_timing_t_payload_uses_exact_products_and_native_inputs(self):
        raw = (_Histogram("raw-t0"), _Histogram("raw-t1"))
        result = _timing_result()
        application = _timing_application(raw)
        payload = plots.build_full_background_subtraction_d6_payload(result, application)

        self.assertTrue(payload["available"])
        self.assertEqual(payload["t_edges"], [0.0, 1.0, 2.0])
        self.assertEqual(payload["delta_edges"], [-10.0, 0.0, 10.0])
        self.assertEqual(payload["per_t"][0]["raw_mm"]["histogram"], raw[0])
        self.assertEqual(payload["per_t"][1]["raw_mm"]["histogram"], raw[1])
        self.assertEqual(
            payload["per_t"][0]["pid"]["timing_axis_label"], "RF timing [ns]"
        )
        self.assertEqual(
            payload["per_t"][0]["pid"]["cell_histograms"][1].label,
            "timing-d1-t0",
        )
        self.assertTrue(payload["per_t"][0]["pid"]["aerogel_validation_available"])
        self.assertEqual(
            payload["per_t"][1]["weight"]["histogram"],
            application["H_proton_weight_vs_delta_t"],
        )

    def test_legacy_payload_keeps_method_wide_inputs_explicit(self):
        raw = (_Histogram("raw-t0"), _Histogram("raw-t1"))
        result = _legacy_result()
        application = _legacy_application(raw)
        payload = plots.build_full_background_subtraction_d6_payload(result, application)

        self.assertTrue(payload["available"])
        for group in payload["per_t"]:
            self.assertTrue(group["pid"]["available"])
            self.assertEqual(group["pid"]["kind"], "ctime_aero")
            self.assertTrue(group["pid"]["shared_across_t"])
            self.assertEqual(group["pid"]["timing_axis_label"], "Coincidence time [ns]")
            self.assertEqual(group["weight"]["kind"], "ctime_aero_event_weight")

    def test_geometry_is_copied_but_source_objects_are_not_substituted(self):
        raw = (_Histogram("raw-t0"), _Histogram("raw-t1"))
        result = _timing_result()
        application = _timing_application(raw)
        payload = plots.build_full_background_subtraction_d6_payload(result, application)
        payload["delta_edges"][0] = -99.0
        payload["per_t"][0]["pid"]["timing_axis_label"] = "changed"

        self.assertEqual(result["delta_edges"][0], -10.0)
        self.assertEqual(result["H_delta_timing_t_cells"][0][0].title, "original")
        self.assertEqual(application["canonical_t_products"][0]["t_edges"], [0.0, 1.0])
        self.assertEqual(application["canonical_t_products"][0]["raw_targets"]["h_mm_nosub"], raw[0])

    def test_missing_raw_is_local_and_unsupported_method_has_no_pid_weight(self):
        application = _timing_application((_Histogram("raw-t0"), None))
        payload = plots.build_full_background_subtraction_d6_payload(
            {"method": "unknown", "delta_edges": [-1.0, 1.0]}, application
        )

        self.assertTrue(payload["available"])
        self.assertFalse(payload["method_supported"])
        self.assertTrue(payload["per_t"][0]["raw_mm"]["available"])
        self.assertFalse(payload["per_t"][1]["raw_mm"]["available"])
        for group in payload["per_t"]:
            self.assertFalse(group["pid"]["available"])
            self.assertFalse(group["weight"]["available"])

    def test_invalid_delta_is_local_but_invalid_t_geometry_is_unavailable(self):
        result = _timing_result()
        result["delta_edges"] = [-10.0, -10.0]
        invalid_delta = plots.build_full_background_subtraction_d6_payload(
            result, _timing_application((_Histogram("raw"), _Histogram("raw")))
        )
        self.assertTrue(invalid_delta["available"])
        self.assertFalse(invalid_delta["delta_geometry_available"])
        self.assertEqual(invalid_delta["delta_edges"], [])
        self.assertTrue(invalid_delta["per_t"][0]["raw_mm"]["available"])
        self.assertEqual(
            invalid_delta["per_t"][0]["pid"]["reason"],
            "proton_delta_edges_invalid",
        )
        self.assertEqual(
            invalid_delta["per_t"][0]["weight"]["reason"],
            "proton_delta_edges_invalid",
        )

        invalid_t = plots.build_full_background_subtraction_d6_payload(
            _timing_result(),
            {
                "canonical_t_products": (
                    {
                        "t_index": 0,
                        "t_edges": [0.0, 1.0],
                        "raw_targets": {"h_mm_nosub": _Histogram("raw")},
                    },
                    {
                        "t_index": 1,
                        "t_edges": [1.1, 2.0],
                        "raw_targets": {"h_mm_nosub": _Histogram("raw")},
                    },
                ),
                "H_proton_weight_vs_delta_t": _Histogram("weight"),
            },
        )
        self.assertFalse(invalid_t["available"])
        self.assertEqual(invalid_t["reason"], "canonical_t_edges_noncontiguous")

    def test_d7_uses_exact_pre_rf_sources_and_committed_signed_rows(self):
        raw = (_Histogram("raw-t0"), _Histogram("raw-t1"))
        proton = (_Histogram("proton-t0"), _Histogram("proton-t1"))
        cleaned = (_Histogram("cleaned-t0"), _Histogram("cleaned-t1"))
        application = {
            "canonical_t_products": _products(raw, proton, cleaned),
            "H_proton_weight_vs_delta_t": _Histogram("weight"),
        }
        result = _d7_timing_result()
        prepared = _d7_prepared_bundle()
        prepared["prepared_sources"]["prompt"]["tree"] = _ForbiddenTree()
        payload = plots.build_full_background_subtraction_d7_payload(
            result, application, prepared
        )

        self.assertTrue(payload["available"])
        self.assertEqual(payload["t_edges"], [0.0, 1.0, 2.0])
        self.assertEqual(payload["delta_edges"], [-10.0, 0.0, 10.0])
        self.assertIs(payload["per_t"][0]["raw_mm"]["histogram"], raw[0])
        self.assertIs(payload["per_t"][0]["proton_mm"]["histogram"], proton[0])
        self.assertIs(payload["per_t"][0]["cleaned_mm"]["histogram"], cleaned[0])
        self.assertIsNot(
            payload["per_t"][0]["raw_mm"]["histogram"],
            application["canonical_t_products"][0]["final_targets"]["h_mm_nosub"],
        )
        delta_zero = payload["per_t"][0]["delta_projection"]["rows_by_delta"][0]
        self.assertEqual(len(delta_zero), 2)
        self.assertEqual(delta_zero[0]["raw_weight"], 2.0)
        self.assertEqual(delta_zero[0]["proton_contribution"], 0.5)
        self.assertEqual(delta_zero[0]["cleaned_contribution"], 1.5)
        self.assertEqual(delta_zero[1]["raw_weight"], -1.0)
        self.assertEqual(delta_zero[1]["proton_contribution"], -0.5)
        self.assertEqual(delta_zero[1]["cleaned_contribution"], -0.5)
        exclusions = payload["per_t"][0]["delta_projection"]["exclusions"]
        self.assertEqual(exclusions["invalid_frozen_delta_index"], 1)
        self.assertEqual(exclusions["selected_rows"], 3)
        self.assertEqual(payload["projection_exclusions"]["invalid_frozen_delta_index"], 1)
        self.assertEqual(raw[0].title, "original")
        self.assertEqual(
            prepared["prepared_sources"]["prompt"]["entries"][0]["adj_mm"], 0.2
        )

    def test_d7_legacy_rows_use_shared_final_edge_membership(self):
        result = _legacy_result()
        result["_prepared_event_weight_lookup"] = {
            "prompt:0": {
                "delta_index": 1,
                "proton_weight": 0.25,
                "cleaned_factor": 0.75,
            }
        }
        prepared = {
            "prepared_sources": {
                "prompt": {
                    "coefficient": 1.0,
                    "entries": {0: {"nommcuts": True, "adj_mm": 1.5, "adj_t": 2.0}},
                }
            }
        }
        payload = plots.build_full_background_subtraction_d7_payload(
            result,
            _legacy_application((_Histogram("raw-t0"), _Histogram("raw-t1"))),
            prepared,
        )

        self.assertEqual(
            len(payload["per_t"][0]["delta_projection"]["rows_by_delta"][1]), 0
        )
        self.assertEqual(
            len(payload["per_t"][1]["delta_projection"]["rows_by_delta"][1]), 1
        )
        self.assertEqual(
            payload["per_t"][1]["delta_projection"]["exclusions"][
                "legacy_t_membership_rows"
            ],
            1,
        )

    def test_d7_missing_lookup_is_local_to_delta_page_and_has_no_aliases(self):
        raw = (_Histogram("raw-t0"), _Histogram("raw-t1"))
        application = _timing_application(raw)
        prepared = _d7_prepared_bundle()
        payload = plots.build_full_background_subtraction_d7_payload(
            _timing_result(), application, prepared
        )

        self.assertTrue(payload["available"])
        self.assertTrue(payload["per_t"][0]["raw_mm"]["available"])
        self.assertTrue(payload["per_t"][0]["proton_mm"]["available"])
        self.assertFalse(payload["per_t"][0]["delta_projection"]["available"])
        self.assertEqual(
            payload["per_t"][0]["delta_projection"]["reason"],
            "frozen_applied_lookup_missing",
        )
        payload["per_t"][0]["delta_projection"]["exclusions"]["changed"] = True
        self.assertNotIn("changed", prepared["prepared_sources"]["prompt"]["entries"][0])

    def test_d8_uses_only_final_parent_application_and_preserves_zero_fallback(self):
        parents, cache = _d8_fixture()
        payload = plots.build_full_background_subtraction_d8_payload(parents, cache)

        self.assertTrue(payload["available"])
        self.assertEqual(payload["t_edges"], [0.0, 1.0, 2.0])
        self.assertEqual(payload["delta_edges"], [-10.0, 0.0, 10.0])
        first = payload["per_t"][0]
        final = parents[0]["final_diagnostic_application_result"]
        proposed = parents[0]["proposed_diagnostic_application_result"]
        self.assertIs(first["before_pion_mm"]["histogram"], final["H_MM_nosub_before_pion_subtraction"])
        self.assertIs(first["baseline_pion_mm"]["histogram"], final["H_pion_subtraction_template_MM_nosub"])
        self.assertIs(first["after_pion_mm"]["histogram"], final["H_MM_nosub_after_pion_subtraction"])
        self.assertIsNot(first["baseline_pion_mm"]["histogram"], proposed["H_pion_subtraction_template_MM_nosub"])
        self.assertEqual(first["delta_projection"]["closure"]["status"], "closed")
        self.assertTrue(payload["per_t"][1]["delta_projection"]["available"])
        self.assertEqual(payload["per_t"][1]["delta_projection"]["closure"]["status"], "closed")
        payload["t_edges"][0] = -99.0
        first["delta_projection"]["rows_by_delta"][0][0]["baseline_contribution"] = 99.0
        self.assertEqual(parents[0]["t_edges"][0], 0.0)
        self.assertEqual(cache["by_t"][0]["records"][0]["coefficient"], 1.0)
        self.assertEqual(final["H_pion_subtraction_template_MM_nosub"].contents, [0.0, 3.0])

    def test_d8_uses_frozen_delta_index_and_signed_final_weights(self):
        parents, cache = _d8_fixture()
        payload = plots.build_full_background_subtraction_d8_payload(parents, cache)
        rows = payload["per_t"][0]["delta_projection"]["rows_by_delta"]

        self.assertEqual(len(rows[0]), 2)
        self.assertEqual(len(rows[1]), 1)
        self.assertEqual(rows[0][0]["baseline_contribution"], 2.0)
        self.assertEqual(rows[1][0]["baseline_contribution"], -2.0)
        self.assertEqual(rows[0][1]["baseline_contribution"], 3.0)
        self.assertEqual(payload["per_t"][0]["delta_projection"]["exclusions"]["non_nommcuts_records"], 1)

    def test_d8_single_scale_final_fallback_uses_frozen_final_application(self):
        parents, cache = _d8_single_scale_final_fallback_fixture()
        final = parents[0]["final_diagnostic_application_result"]
        proposed = parents[0]["proposed_diagnostic_application_result"]
        payload = plots.build_full_background_subtraction_d8_payload(parents, cache)

        self.assertTrue(payload["available"])
        self.assertEqual(
            parents[0]["final_diagnostic_application_status"],
            {"final_status": "applied_fallback", "fallback_mode": "single_scale"},
        )
        group = payload["per_t"][0]
        self.assertIs(
            group["before_pion_mm"]["histogram"],
            final["H_MM_nosub_before_pion_subtraction"],
        )
        self.assertIs(
            group["baseline_pion_mm"]["histogram"],
            final["H_pion_subtraction_template_MM_nosub"],
        )
        self.assertIs(
            group["after_pion_mm"]["histogram"],
            final["H_MM_nosub_after_pion_subtraction"],
        )
        for payload_key, final_key in (
            ("before_pion_mm", "H_MM_nosub_before_pion_subtraction"),
            ("baseline_pion_mm", "H_pion_subtraction_template_MM_nosub"),
            ("after_pion_mm", "H_MM_nosub_after_pion_subtraction"),
        ):
            self.assertIsNot(group[payload_key]["histogram"], proposed[final_key])

        projection = group["delta_projection"]
        self.assertTrue(projection["available"])
        self.assertEqual(projection["closure"]["status"], "closed")
        self.assertTrue(projection["closure"]["coverage_complete"])
        positive = projection["rows_by_delta"][0][0]["baseline_contribution"]
        negative = projection["rows_by_delta"][1][0]["baseline_contribution"]
        self.assertEqual(
            positive,
            2.0 * plots.simc_shape_pion_weight_from_value(
                0.2, final["H_pion_control_model"], final["weights"]
            ),
        )
        self.assertEqual(
            negative,
            -1.0 * plots.simc_shape_pion_weight_from_value(
                1.2, final["H_pion_control_model"], final["weights"]
            ),
        )
        self.assertEqual((positive, negative), (5.0, -2.5))
        self.assertGreater(positive, 0.0)
        self.assertLess(negative, 0.0)

    def test_d8_unavailable_final_never_substitutes_the_proposal(self):
        parents, cache = _d8_fixture()
        parents = list(parents)
        parents[0] = dict(
            parents[0],
            final_diagnostic_application_result=None,
            final_diagnostic_application_status={"final_reason": "skip_bin"},
        )
        payload = plots.build_full_background_subtraction_d8_payload(tuple(parents), cache)

        self.assertTrue(payload["available"])
        first = payload["per_t"][0]
        self.assertFalse(first["before_pion_mm"]["available"])
        self.assertEqual(first["before_pion_mm"]["reason"], "skip_bin")
        self.assertFalse(first["delta_projection"]["available"])
        self.assertEqual(first["delta_projection"]["reason"], "skip_bin")

    def test_d8_invalid_delta_coverage_and_numerical_mismatch_are_local(self):
        parents, cache = _d8_fixture()
        by_t = list(cache["by_t"])
        records = list(by_t[0]["records"])
        records[1] = dict(records[1], delta_index=9)
        by_t[0] = dict(by_t[0], records=tuple(records))
        incomplete = plots.build_full_background_subtraction_d8_payload(
            parents, dict(cache, by_t=tuple(by_t))
        )
        first = incomplete["per_t"][0]["delta_projection"]
        self.assertTrue(first["available"])
        self.assertEqual(first["closure"]["status"], "incomplete_frozen_delta_coverage")
        self.assertEqual(first["exclusions"]["invalid_frozen_delta_index"], 1)

        parents, cache = _d8_fixture()
        parents[0]["final_diagnostic_application_result"]["H_pion_subtraction_template_MM_nosub"].contents[1] = 99.0
        mismatch = plots.build_full_background_subtraction_d8_payload(parents, cache)
        first = mismatch["per_t"][0]
        self.assertTrue(first["before_pion_mm"]["available"])
        self.assertTrue(first["baseline_pion_mm"]["available"])
        self.assertFalse(first["delta_projection"]["available"])
        self.assertEqual(first["delta_projection"]["reason"], "baseline_template_closure_mismatch")

    def test_d8_parent_cache_geometry_mismatch_is_unavailable(self):
        parents, cache = _d8_fixture()
        bad_cache = dict(cache, by_t=(dict(cache["by_t"][0], t_edges=[0.0, 1.1]), cache["by_t"][1]))
        payload = plots.build_full_background_subtraction_d8_payload(parents, bad_cache)

        self.assertFalse(payload["available"])
        self.assertEqual(payload["reason"], "pion_parent_cache_t_geometry_mismatch")

    def test_d9_uses_only_frozen_part1_rows_weighted_templates_and_d2_values(self):
        diagnostic, method_a, comparison = _d9_fixture()
        payload = plots.build_full_background_subtraction_d9_payload(
            diagnostic, method_a, comparison
        )

        self.assertTrue(payload["available"])
        self.assertEqual(payload["schema_version"], plots.D9_PRESENTATION_SCHEMA_VERSION)
        self.assertEqual(payload["t_edges"], [0.0, 1.0, 2.0])
        self.assertEqual(payload["delta_edges"], [-10.0, 0.0, 10.0])
        self.assertEqual(payload["method_a_thresholds"], {
            "available": True,
            "positive_response_threshold": 0.0,
            "low_response_upper_threshold": 2.0,
        })
        first = payload["per_t"][0]
        self.assertIs(
            first["hgcer_response"]["kaon_template"],
            diagnostic["histograms"]["H_hgcer_kaon_weighted"],
        )
        self.assertIs(
            first["hgcer_delta_response"]["pion_template"],
            diagnostic["histograms"]["H_hgcer_vs_delta_pion_weighted"],
        )
        self.assertIsNot(
            first["hgcer_response"]["kaon_template"],
            diagnostic["histograms"]["H_hgcer_kaon_absolute"],
        )
        first_row = first["hgcer_response"]["kaon_rows"][0]
        self.assertEqual(first_row["t_index"], 0)
        self.assertEqual(first_row["delta_index"], 0)
        self.assertEqual(first_row["delta"], 7.5)
        self.assertEqual(first_row["diagnostic_weight"], 3.25)
        self.assertNotIn("reconstructed_weight", first_row)
        cells = first["method_a_relative"]["cells"]
        self.assertEqual(cells[0]["method_a_comparison_candidate"], 0.0)
        self.assertEqual(cells[0]["method_a_comparison_candidate_high"], 0.25)
        self.assertEqual(cells[1]["method_a_comparison_candidate_status"], "marginal")
        self.assertNotIn("forbidden_recalculation_sentinel", cells[0])

        payload["t_edges"][0] = -99.0
        first_row["diagnostic_weight"] = 99.0
        cells[0]["method_a_comparison_candidate"] = 99.0
        self.assertEqual(diagnostic["t_edges"][0], 0.0)
        self.assertEqual(diagnostic["records"]["kaon"][0]["diagnostic_weight"], 3.25)
        self.assertEqual(comparison["cells"][0]["method_a_comparison_candidate"], 0.0)

    def test_d9_threshold_failure_is_local_to_annotation_and_d2_failures_are_local(self):
        diagnostic, method_a, comparison = _d9_fixture()
        method_a["configuration"] = dict(
            method_a["configuration"], low_response_upper_threshold=2.1
        )
        threshold_payload = plots.build_full_background_subtraction_d9_payload(
            diagnostic, method_a, comparison
        )
        self.assertTrue(threshold_payload["available"])
        self.assertFalse(threshold_payload["method_a_thresholds"]["available"])
        self.assertEqual(
            threshold_payload["method_a_threshold_reason"], "method_a_thresholds_invalid"
        )
        self.assertTrue(threshold_payload["per_t"][0]["hgcer_response"]["available"])
        self.assertTrue(threshold_payload["per_t"][0]["method_a_relative"]["available"])

        diagnostic, method_a, comparison = _d9_fixture()
        comparison["method_a_fingerprint"] = "wrong"
        payload = plots.build_full_background_subtraction_d9_payload(
            diagnostic, method_a, comparison
        )
        self.assertTrue(payload["available"])
        self.assertTrue(payload["per_t"][0]["hgcer_response"]["available"])
        self.assertFalse(payload["per_t"][0]["method_a_relative"]["available"])
        self.assertEqual(
            payload["per_t"][0]["method_a_relative"]["reason"],
            "method_a_comparison_provenance_mismatch",
        )

    def test_d9_rejects_part1_coordinate_geometry_and_d2_grid_mismatches(self):
        for mutation, expected in (
            (lambda diagnostic, _method, _comparison: diagnostic.update(coordinate_fingerprint=""), "part1_coordinate_fingerprint_missing"),
            (lambda diagnostic, _method, _comparison: diagnostic.update(delta_edges=[-10.0, -10.0]), "part1_geometry_invalid"),
            (lambda _diagnostic, _method, comparison: comparison["cells"].pop(), None),
        ):
            diagnostic, method_a, comparison = _d9_fixture()
            mutation(diagnostic, method_a, comparison)
            payload = plots.build_full_background_subtraction_d9_payload(
                diagnostic, method_a, comparison
            )
            if expected is not None:
                self.assertFalse(payload["available"])
                self.assertEqual(payload["reason"], expected)
            else:
                self.assertTrue(payload["available"])
                self.assertEqual(
                    payload["per_t"][0]["method_a_relative"]["reason"],
                    "method_a_comparison_cell_grid_invalid",
                )

    def test_d9_part1_unavailability_requires_retained_geometry_for_page12(self):
        diagnostic, method_a, comparison = _d9_fixture()
        diagnostic["status"] = "unavailable"
        payload = plots.build_full_background_subtraction_d9_payload(
            diagnostic, method_a, comparison
        )
        self.assertTrue(payload["available"])
        self.assertFalse(payload["per_t"][0]["hgcer_response"]["available"])
        self.assertFalse(payload["per_t"][0]["hgcer_delta_response"]["available"])
        self.assertTrue(payload["per_t"][0]["method_a_relative"]["available"])

        diagnostic, method_a, comparison = _d9_fixture()
        diagnostic.update(status="unavailable", delta_edges=[])
        payload = plots.build_full_background_subtraction_d9_payload(
            diagnostic, method_a, comparison
        )
        self.assertFalse(payload["available"])
        self.assertEqual(payload["reason"], "part1_geometry_invalid")

    def test_d10_uses_frozen_phase_a_method_b_and_d3_sources_without_aliases(self):
        phase_a, method_b, comparison = _d10_fixture()
        phase_before = deepcopy(phase_a)
        method_before = deepcopy(method_b)
        comparison_before = deepcopy(comparison)
        payload = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )

        self.assertTrue(payload["available"])
        self.assertEqual(payload["schema_version"], plots.D10_PRESENTATION_SCHEMA_VERSION)
        self.assertEqual(payload["mm_edges"], [0.80, 0.91, 1.03, 1.17, 1.35])
        self.assertEqual(payload["host_label"], "Proton-cleaned kaon data")
        self.assertEqual(
            [
                (row["region_name"], row["mm_low"], row["mm_high"], row["mm_offset_applied"])
                for row in payload["mm_regions"]
            ],
            [
                ("pion_sensitive_low", 0.80, 1.10, 0.0),
                ("pion_sensitive_high", 1.23, 1.35, 0.0),
            ],
        )
        self.assertEqual(
            set(payload["protected_regions"][0]),
            {
                "region_name", "mm_low", "mm_high", "region_role", "window_source",
                "mm_offset_applied",
            },
        )
        self.assertEqual(
            payload["per_t"][0]["mm_inputs"]["host_rows"][0]["signed_contribution"], 3.0
        )
        self.assertEqual(
            payload["per_t"][0]["mm_inputs"]["baseline_rows"][0]["signed_contribution"], -2.5
        )
        first_bin = payload["per_t"][0]["local_closure"]["cells_by_delta"][0]["bins"][1]
        self.assertEqual(first_bin["host_yield"], 7.0)
        self.assertEqual(first_bin["host_sumw2"], 9.0)
        self.assertEqual(first_bin["baseline_yield"], -2.0)
        self.assertEqual(first_bin["baseline_sumw2"], 4.0)
        first_relative = payload["per_t"][0]["method_b_relative"]["cells"][0]
        self.assertEqual(first_relative["method_b_comparison_candidate"], 1.25)
        self.assertEqual(first_relative["method_b_comparison_candidate_uncertainty"], 0.15)
        self.assertNotEqual(first_relative["method_b_comparison_candidate"], 9.0)
        self.assertIsNone(
            payload["per_t"][0]["method_b_relative"]["cells"][1][
                "method_b_comparison_candidate"
            ]
        )

        payload["mm_regions"][0]["mm_low"] = -123.0
        payload["per_t"][0]["mm_inputs"]["host_rows"][0]["signed_contribution"] = 44.0
        payload["per_t"][0]["local_closure"]["cells_by_delta"][0]["bins"][1]["host_yield"] = 55.0
        self.assertEqual(phase_a, phase_before)
        self.assertEqual(method_b, method_before)
        self.assertEqual(comparison, comparison_before)

    def test_d10_none_membership_is_outside_geometry_but_non_none_invalid_is_local(self):
        phase_a, method_b, comparison = _d10_fixture()
        phase_a["kaon_host_records"].extend((
            {
                "canonical_t_index": None,
                "delta_index": None,
                "analysis_MM": 99.0,
                "analysis_abs_t": -999.0,
                "SHMS_delta": -999.0,
                "nommcuts": True,
                "host_state": "proton_cleaned",
                "signed_host_event_contribution": 44.0,
            },
            {
                "canonical_t_index": 0,
                "delta_index": None,
                "analysis_MM": 98.0,
                "analysis_abs_t": 999.0,
                "SHMS_delta": 999.0,
                "nommcuts": True,
                "host_state": "proton_cleaned",
                "signed_host_event_contribution": 55.0,
            },
        ))
        payload = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        self.assertTrue(payload["available"])
        self.assertTrue(payload["per_t"][0]["mm_inputs"]["available"])
        self.assertEqual(len(payload["per_t"][0]["mm_inputs"]["host_rows"]), 2)
        self.assertTrue(payload["per_t"][0]["local_closure"]["available"])
        self.assertTrue(payload["per_t"][0]["method_b_relative"]["available"])

        for key, invalid in (("canonical_t_index", 99), ("delta_index", 99), ("delta_index", 0.5)):
            with self.subTest(key=key, invalid=invalid):
                phase_a, method_b, comparison = _d10_fixture()
                phase_a["kaon_host_records"][0][key] = invalid
                payload = plots.build_full_background_subtraction_d10_payload(
                    phase_a, method_b, comparison
                )
                self.assertTrue(payload["available"])
                self.assertFalse(payload["per_t"][0]["mm_inputs"]["available"])
                self.assertTrue(payload["per_t"][0]["local_closure"]["available"])
                self.assertTrue(payload["per_t"][0]["method_b_relative"]["available"])

    def test_d10_page13_histogram_uses_native_upper_edge_and_signed_sumw2(self):
        root = _D10ROOT()
        edges = [0.80, 0.91, 1.03, 1.17, 1.35]
        histogram = plots._d10_event_histogram(
            root,
            "d10-upper-edge",
            edges,
            (
                {"missing_mass": 1.35, "signed_contribution": 6.25},
                {"missing_mass": 1.36, "signed_contribution": -2.0},
                {"missing_mass": 0.80, "signed_contribution": -3.0},
            ),
        )
        self.assertIsNotNone(histogram)
        self.assertEqual(histogram.GetNbinsX(), 4)
        self.assertEqual(
            tuple(histogram.GetXaxis().GetBinLowEdge(index) for index in range(1, 6)),
            tuple(edges),
        )
        self.assertEqual(histogram.GetBinContent(4), 6.25)
        self.assertEqual(histogram.GetBinError(4), 6.25)
        self.assertEqual(histogram.GetBinContent(5), -2.0)
        self.assertEqual(histogram.GetBinError(5), 2.0)
        self.assertEqual(histogram.GetBinContent(1), -3.0)
        self.assertEqual(histogram.GetBinError(1), 3.0)

    def test_d10_page13_renderer_locks_signed_range_boundaries_and_text(self):
        phase_a, method_b, comparison = _d10_fixture()
        phase_a["kaon_host_records"][0].update(
            analysis_MM=0.85, signed_host_event_contribution=-7.0
        )
        phase_a["kaon_host_records"][1].update(
            analysis_MM=0.95, signed_host_event_contribution=3.0
        )
        phase_a["pion_records"][0].update(
            analysis_MM=0.85, signed_baseline_event_contribution=-2.0
        )
        phase_a["pion_records"][1].update(
            analysis_MM=0.95, signed_baseline_event_contribution=11.0
        )
        presentation = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        root = _D10RenderROOT()
        self.assertTrue(
            plots._render_d10_mm_inputs_page(root, "ignored.pdf", presentation, presentation["per_t"][0])
        )

        host = _d10_histogram(root, "H_full_background_d10_inputs_host_t1")
        self.assertEqual(host.title, ";Missing mass [GeV];Signed normalized yield")
        self.assertNotIn("Method B - Missing-mass closure inputs", host.title)
        self.assertLessEqual(host.display_minimum, -7.0)
        self.assertGreaterEqual(host.display_maximum, 11.0)
        self.assertLessEqual(host.display_minimum, 0.0)
        self.assertGreaterEqual(host.display_maximum, 0.0)
        self.assertEqual(
            {line.x1 for line in root.drawn_lines},
            {1.10, 1.23},
        )
        self.assertTrue(all(line.x1 == line.x2 for line in root.drawn_lines))
        visible = _d10_visible_text(root)
        for expected in (
            "Method B - Missing-mass closure inputs",
            "Proton-cleaned kaon data",
            "Baseline pion background",
            "Pion-sensitive region (outside protected window)",
            "Lambda/Sigma protected region (1.10-1.23 GeV)",
            "|t| = [0.0000, 1.0000] GeV^2",
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected, visible)
        self.assertEqual(visible.count("Method B - Missing-mass closure inputs"), 1)
        for internal in ("host_state", "phase_a", "window_source", "candidate_L_B"):
            with self.subTest(internal=internal):
                self.assertNotIn(internal, visible)

    def test_d10_page13_renderer_uses_identity_host_physics_label(self):
        phase_a, method_b, comparison = _d10_fixture()
        for mapping in (phase_a, method_b, comparison):
            mapping["host_state"] = "identity_no_proton_cleaning"
        for record in phase_a["kaon_host_records"]:
            record["host_state"] = "identity_no_proton_cleaning"
        for cell in method_b["cells"]:
            cell["host_state"] = "identity_no_proton_cleaning"
        presentation = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        self.assertTrue(presentation["available"])
        root = _D10RenderROOT()
        self.assertTrue(
            plots._render_d10_mm_inputs_page(root, "ignored.pdf", presentation, presentation["per_t"][0])
        )
        visible = _d10_visible_text(root)
        self.assertIn("Kaon-selected data", visible)
        self.assertNotIn("Proton-cleaned kaon data", visible)
        self.assertNotIn("identity_no_proton_cleaning", visible)

    def test_d10_page14_renderer_uses_independent_ranges_delta_order_and_boundaries(self):
        phase_a, method_b, comparison = _d10_fixture()
        for cell in method_b["cells"]:
            if cell["t_index"] != 0:
                continue
            for row in cell["bins"][1:-1]:
                row.update(host_yield=0.0, host_sumw2=0.0, baseline_yield=0.0, baseline_sumw2=0.0)
        method_b["cells"][0]["bins"][1].update(
            host_yield=4.0, host_sumw2=9.0, baseline_yield=-3.0, baseline_sumw2=4.0
        )
        method_b["cells"][1]["bins"][1].update(
            host_yield=-9.0, host_sumw2=16.0, baseline_yield=12.0, baseline_sumw2=25.0
        )
        presentation = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        root = _D10RenderROOT()
        self.assertTrue(
            plots._render_d10_local_closure_page(
                root, "ignored.pdf", presentation, presentation["per_t"][0]
            )
        )

        hosts = [
            _d10_histogram(root, "H_full_background_d10_local_t1_d{}_host".format(index))
            for index in (1, 2)
        ]
        ranges = [(histogram.display_minimum, histogram.display_maximum) for histogram in hosts]
        self.assertNotEqual(ranges[0], ranges[1])
        self.assertLessEqual(ranges[0][0], -5.0)
        self.assertGreaterEqual(ranges[0][1], 7.0)
        self.assertLessEqual(ranges[1][0], -13.0)
        self.assertGreaterEqual(ranges[1][1], 17.0)
        for low, high in ranges:
            self.assertLessEqual(low, 0.0)
            self.assertGreaterEqual(high, 0.0)
        self.assertEqual(
            [histogram.title.split(";")[0] for histogram in hosts],
            ["delta = [-10.000, 0.000] %", "delta = [0.000, 10.000] %"],
        )
        self.assertTrue(all(histogram.edges == (0.80, 0.91, 1.03, 1.17, 1.35) for histogram in hosts))
        self.assertEqual(hosts[0].GetBinContent(1), 4.0)
        self.assertEqual(hosts[0].GetBinError(1), 3.0)
        baseline = _d10_histogram(root, "H_full_background_d10_local_t1_d1_baseline")
        self.assertEqual(baseline.GetBinContent(1), -3.0)
        self.assertEqual(baseline.GetBinError(1), 2.0)
        self.assertEqual(len(root.drawn_lines), 4)
        self.assertEqual(
            [{line.x1 for line in root.drawn_lines[offset:offset + 2]} for offset in (0, 2)],
            [{1.10, 1.23}] * 2,
        )
        visible = _d10_visible_text(root)
        self.assertIn("Method B - Local missing-mass closure", visible)
        self.assertNotIn("Diagnostic only - no refinement applied", visible)
        self.assertIn("|t| = [0.0000, 1.0000] GeV^2", visible)

    def test_d10_page15_renderer_locks_d3_points_status_unity_range_and_text(self):
        phase_a, method_b, comparison = _d10_fixture()
        presentation = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        root = _D10RenderROOT()
        self.assertTrue(
            plots._render_d10_method_b_relative_page(
                root, "ignored.pdf", presentation["per_t"][0], presentation["delta_edges"]
            )
        )
        self.assertEqual(len(root.drawn_graphs), 1)
        graph = root.drawn_graphs[0]
        self.assertEqual(graph.points, [(-5.0, 1.25)])
        self.assertEqual(graph.errors, [(0.0, 0.15)])
        self.assertNotIn((None, 9.0), graph.points)
        self.assertNotIn((0.0, 8.0), graph.errors)
        self.assertEqual([(line.x1, line.y1, line.x2, line.y2) for line in root.drawn_lines], [(-10.0, 1.0, 10.0, 1.0)])
        frame = _d10_histogram(root, "H_full_background_d10_relative_t1")
        self.assertLessEqual(frame.display_minimum, 1.10)
        self.assertGreaterEqual(frame.display_maximum, 1.40)
        visible = _d10_visible_text(root)
        for expected in (
            "Method B - Missing-mass closure diagnostic",
            "Relative pion-background diagnostic",
            "|t| = [0.0000, 1.0000] GeV^2",
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected, visible)

        isolated_cells = []
        statuses = (
            "available_multi_region", "single_region_only", "region_marginal",
            "region_inconsistent", "shape_poor_veto", "unavailable",
        )
        for index, status in enumerate(statuses):
            isolated_cells.append({
                "delta_low": -30.0 + 10.0 * index,
                "delta_high": -20.0 + 10.0 * index,
                "method_b_comparison_candidate_status": status,
                "method_b_comparison_candidate": 0.20 if status == "available_multi_region" else None,
                "method_b_comparison_candidate_uncertainty": 0.35 if status == "available_multi_region" else None,
            })
        isolated_group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "method_b_relative": {"available": True, "cells": isolated_cells},
        }
        isolated_root = _D10RenderROOT()
        self.assertTrue(
            plots._render_d10_method_b_relative_page(
                isolated_root, "ignored.pdf", isolated_group,
                [-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0],
            )
        )
        self.assertEqual(len(isolated_root.drawn_graphs), 1)
        self.assertEqual(isolated_root.drawn_graphs[0].points, [(-25.0, 0.20)])
        isolated_frame = _d10_histogram(isolated_root, "H_full_background_d10_relative_t1")
        self.assertLessEqual(isolated_frame.display_minimum, -0.15)
        self.assertGreaterEqual(isolated_frame.display_maximum, 1.0)

    def test_d11_uses_only_detached_stored_d4_scalars_and_availability(self):
        checkpoint = _d11_fixture()
        before = deepcopy(checkpoint)
        payload = plots.build_full_background_subtraction_d11_payload(checkpoint)

        self.assertTrue(payload["available"])
        self.assertEqual(payload["schema_version"], plots.D11_PRESENTATION_SCHEMA_VERSION)
        self.assertEqual(payload["t_edges"], [0.0, 1.0, 2.0])
        self.assertEqual(payload["delta_edges"], [-12.0, -5.0, 1.0, 6.0, 14.0, 23.0])
        first = payload["per_t"][0]
        self.assertEqual(first["ratio_log"]["cells"][0]["comparison"], {
            "availability": "both_comparable", "ratio_B_over_A": 7.25,
            "log_ratio_B_over_A": -3.5,
        })
        zero = first["ab_overlay"]["cells"][1]
        self.assertEqual(zero["method_a"]["candidate"], 0.0)
        self.assertEqual(zero["method_a"]["low"], 0.0)
        self.assertEqual(zero["method_a"]["high"], 0.25)
        self.assertEqual(zero["method_b"]["candidate"], 1.1)
        self.assertNotIn(zero, first["ratio_log"]["cells"])
        self.assertEqual(
            [cell["availability_label"] for cell in first["availability"]["cells"]],
            [
                "Both methods available", "Both present; ratio undefined", "Method A only",
                "Method B only", "Neither method available",
            ],
        )
        self.assertNotIn("diagnostic_interval_relation", first["availability"]["cells"][0]["comparison"])

        checkpoint["ab_comparison"]["cells"][0]["comparison"]["diagnostic_interval_relation"] = "overlap"
        self.assertEqual(
            plots.build_full_background_subtraction_d11_payload(checkpoint), payload
        )
        payload["per_t"][0]["ab_overlay"]["cells"][0]["method_a"]["candidate"] = 99.0
        self.assertEqual(before["ab_comparison"]["cells"][0]["method_a"]["comparison_candidate"], 2.0)
        self.assertEqual(checkpoint["method_a_comparison"], before["method_a_comparison"])

    def test_d11_rejects_checkpoint_provenance_geometry_and_grid_mismatches(self):
        mutations = (
            lambda checkpoint: checkpoint.update(source_checkpoint_payload_fingerprint="wrong"),
            lambda checkpoint: checkpoint["method_a_comparison"].update(fingerprint="wrong"),
            lambda checkpoint: checkpoint["method_a_comparison"]["fingerprint_inputs"].update(
                tampered=True
            ),
            lambda checkpoint: checkpoint["method_b_comparison"].update(fingerprint="wrong"),
            lambda checkpoint: checkpoint["ab_comparison"]["fingerprint_inputs"].update(
                tampered=True
            ),
            lambda checkpoint: checkpoint["method_a_comparison"].update(
                preserved_d2_value="tampered"
            ),
            lambda checkpoint: checkpoint["ab_comparison"].update(
                source_method_b_comparison_payload_fingerprint="wrong"
            ),
            lambda checkpoint: checkpoint["method_b_comparison"].update(phase_a_contract_fingerprint="wrong"),
            lambda checkpoint: checkpoint["method_b_comparison"].update(coordinate_fingerprint="wrong"),
            lambda checkpoint: checkpoint["method_b_comparison"].update(canonical_t_edges=[0.0, 1.0, 3.0]),
            lambda checkpoint: checkpoint["method_b_comparison"].update(delta_edges=[-12.0, 1.0, 23.0]),
            lambda checkpoint: checkpoint["method_b_comparison"].update(host_state="identity_no_proton_cleaning"),
            lambda checkpoint: checkpoint["method_b_comparison"].update(source_target_state="post_proton_RF"),
            lambda checkpoint: checkpoint["ab_comparison"]["cells"].pop(),
            lambda checkpoint: checkpoint["ab_comparison"]["cells"].append(
                deepcopy(checkpoint["ab_comparison"]["cells"][0])
            ),
            lambda checkpoint: checkpoint["ab_comparison"]["cells"][0].update(t_low=-1.0),
            lambda checkpoint: checkpoint["ab_comparison"]["cells"][0].update(delta_index=99),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                checkpoint = _d11_fixture()
                mutation(checkpoint)
                payload = plots.build_full_background_subtraction_d11_payload(checkpoint)
                self.assertFalse(payload["available"])
                self.assertEqual(payload["per_t"], [])

    def test_d11_page_local_availability_and_rendered_stored_values(self):
        checkpoint = _d11_fixture()
        for cell in checkpoint["ab_comparison"]["cells"]:
            if cell["t_index"] != 1:
                continue
            cell["method_a"].update(
                present=False,
                comparison_candidate=None,
                comparison_candidate_low=None,
                comparison_candidate_high=None,
                comparison_candidate_status="unavailable",
            )
            cell["method_b"].update(
                present=False,
                comparison_candidate=None,
                comparison_candidate_uncertainty=None,
                comparison_candidate_status="unavailable",
            )
            cell["comparison"].update(
                availability="neither_available", availability_reason=None,
                ratio_B_over_A=None, log_ratio_B_over_A=None,
            )
        payload = plots.build_full_background_subtraction_d11_payload(checkpoint)
        self.assertTrue(payload["available"])
        self.assertFalse(payload["per_t"][1]["ab_overlay"]["available"])
        self.assertFalse(payload["per_t"][1]["ratio_log"]["available"])
        self.assertTrue(payload["per_t"][1]["availability"]["available"])

        root = _D10RenderROOT()
        first = payload["per_t"][0]
        self.assertTrue(plots._render_d11_ab_overlay_page(root, "ignored.pdf", payload, first))
        self.assertEqual(len(root.drawn_asymm_graphs), 1)
        self.assertIn((-2.0, 0.0), root.drawn_asymm_graphs[0].points)
        zero_index = root.drawn_asymm_graphs[0].points.index((-2.0, 0.0))
        self.assertEqual(root.drawn_asymm_graphs[0].errors[zero_index], (0.0, 0.0, 0.0, 0.25))
        self.assertEqual(len(root.drawn_graphs), 1)
        self.assertIn((-8.5, 4.0), root.drawn_graphs[0].points)
        overlay_frame = _d10_histogram(root, "H_full_background_d11_overlay_t1")
        self.assertLessEqual(overlay_frame.display_minimum, 0.0)
        self.assertGreaterEqual(overlay_frame.display_maximum, 4.5)
        self.assertIn((-12.0, 1.0, 23.0, 1.0), [
            (line.x1, line.y1, line.x2, line.y2) for line in root.drawn_lines
        ])
        visible = _d10_visible_text(root)
        self.assertIn("Method A - HGCer response diagnostic", visible)
        self.assertIn("Method B - Missing-mass closure diagnostic", visible)
        self.assertNotIn("NON-AUTHORITATIVE DIAGNOSTIC", visible)
        self.assertNotIn("No refinement, correction, or method selection", visible)

        ratio_root = _D10RenderROOT()
        self.assertTrue(plots._render_d11_ratio_log_page(ratio_root, "ignored.pdf", payload, first))
        self.assertEqual(ratio_root.drawn_plain_graphs[0].points, [(-8.5, 7.25)])
        self.assertEqual(ratio_root.drawn_plain_graphs[1].points, [(-8.5, -3.5)])
        self.assertEqual(
            [(line.x1, line.y1, line.x2, line.y2) for line in ratio_root.drawn_lines],
            [(-12.0, 1.0, 23.0, 1.0), (-12.0, 0.0, 23.0, 0.0)],
        )

        availability_root = _D10RenderROOT()
        self.assertTrue(plots._render_d11_availability_page(availability_root, "ignored.pdf", payload, first))
        self.assertEqual(
            [(box.x_low, box.x_high) for box in availability_root.drawn_boxes],
            [(-12.0, -5.0), (-5.0, 1.0), (1.0, 6.0), (6.0, 14.0), (14.0, 23.0)],
        )
        labels = [label for legend in availability_root.legends for label, _option in legend.entries]
        self.assertEqual(labels, [
            "Both methods available", "Both present; ratio undefined", "Method A only",
            "Method B only", "Neither method available",
        ])
        legend = availability_root.legends[-1]
        self.assertEqual(legend.coordinates, (0.55, 0.44, 0.89, 0.76))
        self.assertFalse(any(
            "NON-AUTHORITATIVE DIAGNOSTIC" in pave._text
            or "No refinement, correction, or method selection" in pave._text
            for pave in availability_root.pave_texts
        ))

    def test_d12_cumulative_omissions_remain_local_and_t_ordered(self):
        """D.12: omitted pages never introduce placeholders or cross-t interleaving."""
        phase_page_ids = (
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[0:3],
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[3:6],
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[6:9],
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[9:12],
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[12:15],
            EXPECTED_FULL_BACKGROUND_PAGE_IDS[15:18],
        )
        omissions = {
            ("t1", "full_background.d6.proton_pid"),
            ("t1", "full_background.d8.pion_delta_mm"),
            ("t1", "full_background.d10.method_b_relative"),
            ("t2", "full_background.d9.method_a_relative"),
            ("t2", "full_background.d11.ab_ratio_log"),
        }
        payloads = {
            label: _d12_cumulative_payload(label)
            for label in ("D.6", "D.7", "D.8", "D.9", "D.10", "D.11")
        }
        with patch.object(plots, "_import_root", return_value=object()), patch.object(
            plots,
            "_render_d6_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[0], omissions),
        ), patch.object(
            plots,
            "_render_d7_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[1], omissions),
        ), patch.object(
            plots,
            "_render_d8_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[2], omissions),
        ), patch.object(
            plots,
            "_render_d9_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[3], omissions),
        ), patch.object(
            plots,
            "_render_d10_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[4], omissions),
        ), patch.object(
            plots,
            "_render_d11_t_pages",
            side_effect=_d12_record_phase_pages(phase_page_ids[5], omissions),
        ):
            rendered = plots.render_full_background_subtraction_procedure_pages(
                "ignored.pdf",
                payloads["D.6"],
                payloads["D.7"],
                payloads["D.8"],
                payloads["D.9"],
                payloads["D.10"],
                payloads["D.11"],
            )

        expected_manifest = [
            (scope, page_id)
            for scope in ("t1", "t2")
            for page_id in EXPECTED_FULL_BACKGROUND_PAGE_IDS
            if (scope, page_id) not in omissions
        ]
        self.assertEqual(
            [(page["scope"], page["page_id"]) for page in rendered["manifest"]],
            expected_manifest,
        )
        self.assertEqual(len(rendered["manifest"]), 36 - len(omissions))
        _assert_full_background_manifest_contract(self, rendered["manifest"])
        for scope, page_id in omissions:
            self.assertNotIn(
                (scope, page_id),
                [(page["scope"], page["page_id"]) for page in rendered["manifest"]],
            )
            self.assertIn(
                "{} unavailable for {}".format(page_id, scope),
                rendered["failures"],
            )

    def test_d11_early_display_delta_mismatches_do_not_suppress_d11(self):
        def record_d11(_root, _pdf_name, _presentation, group, manifest, _failures):
            for page_id in (
                "full_background.d11.ab_overlay",
                "full_background.d11.ab_ratio_log",
                "full_background.d11.method_availability",
            ):
                manifest.append({
                    "page_id": page_id,
                    "scope": "t{}".format(int(group["t_index"]) + 1),
                    "authoritative": False,
                })

        def noop(*_args):
            return None

        for mismatched_phase in ("D.6", "D.7", "D.8"):
            with self.subTest(mismatched_phase=mismatched_phase):
                payloads = {
                    phase: _d11_cumulative_payload(phase)
                    for phase in ("D.6", "D.7", "D.8", "D.9", "D.10", "D.11")
                }
                payloads[mismatched_phase]["delta_edges"] = [-12.0, -2.0, 8.0]
                with patch.object(plots, "_import_root", return_value=object()), patch.object(
                    plots, "_render_d6_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d7_t_pages", side_effect=noop), patch.object(
                    plots, "_render_d8_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d9_t_pages", side_effect=noop), patch.object(
                    plots, "_render_d10_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d11_t_pages", side_effect=record_d11):
                    rendered = plots.render_full_background_subtraction_procedure_pages(
                        "ignored.pdf", payloads["D.6"], payloads["D.7"], payloads["D.8"],
                        payloads["D.9"], payloads["D.10"], payloads["D.11"],
                    )
                self.assertIn(
                    "D.11 frozen delta geometry differs from {}".format(mismatched_phase),
                    rendered["failures"],
                )
                self.assertEqual(
                    [page["page_id"] for page in rendered["manifest"]],
                    [
                        "full_background.d11.ab_overlay",
                        "full_background.d11.ab_ratio_log",
                        "full_background.d11.method_availability",
                    ],
                )

    def test_d11_frozen_delta_and_t_mismatches_suppress_d11(self):
        def record_d11(_root, _pdf_name, _presentation, _group, manifest, _failures):
            manifest.append({"page_id": "full_background.d11.method_availability"})

        def noop(*_args):
            return None

        cases = (
            ("D.9", "delta_edges", [-10.0, 5.0, 10.0],
             "D.11 frozen delta geometry differs from D.9"),
            ("D.10", "delta_edges", [-10.0, 5.0, 10.0],
             "D.11 frozen delta geometry differs from D.10"),
            ("D.11", "t_edges", [0.0, 2.0], "D.11 canonical t geometry mismatch"),
        )
        for phase, key, value, expected_failure in cases:
            with self.subTest(phase=phase, key=key):
                payloads = {
                    label: _d11_cumulative_payload(label)
                    for label in ("D.6", "D.7", "D.8", "D.9", "D.10", "D.11")
                }
                payloads[phase][key] = value
                with patch.object(plots, "_import_root", return_value=object()), patch.object(
                    plots, "_render_d6_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d7_t_pages", side_effect=noop), patch.object(
                    plots, "_render_d8_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d9_t_pages", side_effect=noop), patch.object(
                    plots, "_render_d10_t_pages", side_effect=noop
                ), patch.object(plots, "_render_d11_t_pages", side_effect=record_d11):
                    rendered = plots.render_full_background_subtraction_procedure_pages(
                        "ignored.pdf", payloads["D.6"], payloads["D.7"], payloads["D.8"],
                        payloads["D.9"], payloads["D.10"], payloads["D.11"],
                    )
                self.assertIn(expected_failure, rendered["failures"])
                self.assertFalse(any(
                    page["page_id"].startswith("full_background.d11.")
                    for page in rendered["manifest"]
                ))

    def test_d10_d3_failures_are_local_to_the_relative_page(self):
        mutations = (
            ("method_b_fingerprint", "wrong-method-b"),
            ("coordinate_fingerprint", "wrong-coordinate"),
            ("canonical_t_edges", [0.0, 1.0, 3.0]),
            ("delta_edges", [-10.0, 5.0, 10.0]),
            ("host_state", "identity_no_proton_cleaning"),
        )
        for key, value in mutations:
            with self.subTest(key=key):
                phase_a, method_b, comparison = _d10_fixture()
                comparison[key] = value
                payload = plots.build_full_background_subtraction_d10_payload(
                    phase_a, method_b, comparison
                )
                self.assertTrue(payload["available"])
                self.assertTrue(payload["per_t"][0]["mm_inputs"]["available"])
                self.assertTrue(payload["per_t"][0]["local_closure"]["available"])
                self.assertFalse(payload["per_t"][0]["method_b_relative"]["available"])
        phase_a, method_b, comparison = _d10_fixture()
        comparison["cells"].pop()
        payload = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        self.assertTrue(payload["available"])
        self.assertFalse(payload["per_t"][0]["method_b_relative"]["available"])
        self.assertEqual(
            payload["per_t"][0]["method_b_relative"]["reason"],
            "method_b_comparison_cell_grid_invalid",
        )

    def test_d10_phase_method_b_failures_are_global_and_page_failures_are_local(self):
        for field, value in (
            ("phase_a_contract_fingerprint", "wrong-phase"),
            ("coordinate_fingerprint", "wrong-coordinate"),
            ("t_edges", [0.0, 1.0, 3.0]),
            ("delta_edges", [-10.0, 5.0, 10.0]),
            ("host_state", "identity_no_proton_cleaning"),
            ("source_target_state", "post_proton_RF"),
        ):
            with self.subTest(field=field):
                phase_a, method_b, comparison = _d10_fixture()
                method_b[field] = value
                payload = plots.build_full_background_subtraction_d10_payload(
                    phase_a, method_b, comparison
                )
                self.assertFalse(payload["available"])
                self.assertEqual(payload["reason"], "method_b_contract_invalid")

        for field, value in (
            ("region_name", "KLambda"),
            ("mm_low", 1.09),
            ("mm_high", 1.24),
            ("mm_offset_applied", 0.01),
        ):
            with self.subTest(protected_field=field):
                phase_a, method_b, comparison = _d10_fixture()
                method_b["protected_regions"][0][field] = value
                payload = plots.build_full_background_subtraction_d10_payload(
                    phase_a, method_b, comparison
                )
                self.assertFalse(payload["available"])
                self.assertEqual(payload["reason"], "method_b_protected_region_contract_invalid")

        phase_a, method_b, comparison = _d10_fixture()
        phase_a["kaon_host_records"][0]["delta_index"] = 99
        payload = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        self.assertTrue(payload["available"])
        self.assertFalse(payload["per_t"][0]["mm_inputs"]["available"])
        self.assertTrue(payload["per_t"][0]["local_closure"]["available"])
        self.assertTrue(payload["per_t"][0]["method_b_relative"]["available"])

        phase_a, method_b, comparison = _d10_fixture()
        method_b["cells"][0]["bins"][1]["host_sumw2"] = -1.0
        payload = plots.build_full_background_subtraction_d10_payload(
            phase_a, method_b, comparison
        )
        self.assertTrue(payload["available"])
        self.assertTrue(payload["per_t"][0]["mm_inputs"]["available"])
        self.assertFalse(payload["per_t"][0]["local_closure"]["available"])
        self.assertTrue(payload["per_t"][0]["method_b_relative"]["available"])

    def test_d10_accepts_each_native_d3_candidate_status_without_fabricating_points(self):
        statuses = (
            "available_multi_region", "single_region_only", "unavailable", "region_marginal",
            "region_inconsistent", "shape_poor_veto",
        )
        for status in statuses:
            with self.subTest(status=status):
                phase_a, method_b, comparison = _d10_fixture()
                cell = comparison["cells"][0]
                cell["method_b_comparison_candidate_status"] = status
                cell["method_b_comparison_candidate_reason"] = (
                    None if status == "available_multi_region" else status
                )
                cell["method_b_comparison_candidate"] = (
                    1.75 if status == "available_multi_region" else None
                )
                cell["method_b_comparison_candidate_uncertainty"] = (
                    0.20 if status == "available_multi_region" else None
                )
                payload = plots.build_full_background_subtraction_d10_payload(
                    phase_a, method_b, comparison
                )
                self.assertTrue(payload["available"])
                copied = payload["per_t"][0]["method_b_relative"]["cells"][0]
                self.assertEqual(copied["method_b_comparison_candidate_status"], status)
                if status == "available_multi_region":
                    self.assertEqual(copied["method_b_comparison_candidate"], 1.75)
                else:
                    self.assertIsNone(copied["method_b_comparison_candidate"])

    def test_d9_threshold_text_is_generated_from_frozen_payload_values(self):
        self.assertEqual(
            plots._d9_method_a_threshold_text({
                "available": True,
                "positive_response_threshold": 0.25,
                "low_response_upper_threshold": 1.75,
            }),
            "Method A low response: 0.25 < HGCer NPE <= 1.75",
        )
        self.assertIsNone(plots._d9_method_a_threshold_text({"available": False}))

    def test_d9_response_page_uses_detached_signed_weighted_clones_and_threshold_only(self):
        _BinnedHistogram.created = []
        root = _FakeROOT()
        kaon_source = _BinnedHistogram("d9-kaon", (11.0,))
        pion_source = _BinnedHistogram("d9-pion", (13.0,))
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "hgcer_response": {
                "kaon_template": kaon_source,
                "pion_template": pion_source,
                "kaon_rows": ({"npe": 0.5, "diagnostic_weight": -4.0},),
                "pion_rows": ({"npe": 0.5, "diagnostic_weight": 6.0},),
            },
        }
        self.assertTrue(
            plots._render_d9_hgcer_response_page(
                root,
                "unused.pdf",
                group,
                {
                    "available": True,
                    "positive_response_threshold": 0.0,
                    "low_response_upper_threshold": 2.0,
                },
            )
        )
        displays = [
            histogram for histogram in _BinnedHistogram.created
            if "H_full_background_d9_hgcer_" in histogram.label
        ]
        self.assertEqual(len(displays), 2)
        self.assertEqual(
            (displays[0].display_minimum, displays[0].display_maximum), (-4.5, 6.5)
        )
        self.assertEqual((displays[1].display_minimum, displays[1].display_maximum), (None, None))
        self.assertEqual(kaon_source.contents, [11.0])
        self.assertEqual(pion_source.contents, [13.0])
        visible_text = tuple(line for text in root.drawn_text for line in text)
        self.assertNotIn("Diagnostic only - no refinement applied", visible_text)
        self.assertIn("Method A low response: 0 < HGCer NPE <= 2", visible_text)
        self.assertEqual(len(root.legends), 1)
        self.assertEqual(len(root.pave_texts), 1)
        legend_rectangle = root.legends[0].coordinates
        threshold_rectangle = root.pave_texts[0].coordinates
        self.assertEqual(legend_rectangle, (0.62, 0.58, 0.89, 0.72))
        self.assertFalse(_rectangles_overlap(legend_rectangle, threshold_rectangle))

        unavailable_root = _FakeROOT()
        self.assertTrue(
            plots._render_d9_hgcer_response_page(
                unavailable_root, "unused.pdf", group, {"available": False}
            )
        )
        unavailable_text = tuple(
            line for text in unavailable_root.drawn_text for line in text
        )
        self.assertNotIn("Diagnostic only - no refinement applied", unavailable_text)
        self.assertNotIn("Method A low response: 0 < HGCer NPE <= 2", unavailable_text)
        self.assertFalse(any("threshold unavailable" in line for line in unavailable_text))
        self.assertEqual(len(unavailable_root.legends), 1)
        self.assertFalse(unavailable_root.pave_texts)

    def test_d9_delta_page_keeps_only_the_threshold_notice(self):
        _BinnedHistogram2D.created = []
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "hgcer_delta_response": {
                "kaon_template": _BinnedHistogram2D("d9-kaon-delta"),
                "pion_template": _BinnedHistogram2D("d9-pion-delta"),
                "kaon_rows": ({"delta": -5.0, "npe": 0.5, "diagnostic_weight": -4.0},),
                "pion_rows": ({"delta": 5.0, "npe": 1.5, "diagnostic_weight": 6.0},),
            },
        }
        root = _FakeROOT()
        thresholds = {
            "available": True,
            "positive_response_threshold": 0.0,
            "low_response_upper_threshold": 2.0,
        }
        self.assertTrue(
            plots._render_d9_hgcer_delta_page(
                root, "unused.pdf", group, thresholds, (-10.0, 0.0, 10.0)
            )
        )
        visible_text = tuple(line for text in root.drawn_text for line in text)
        self.assertNotIn("Diagnostic only - no refinement applied", visible_text)
        self.assertIn("Method A low response: 0 < HGCer NPE <= 2", visible_text)

        unavailable_root = _FakeROOT()
        self.assertTrue(
            plots._render_d9_hgcer_delta_page(
                unavailable_root, "unused.pdf", group, {"available": False},
                (-10.0, 0.0, 10.0),
            )
        )
        unavailable_text = tuple(
            line for text in unavailable_root.drawn_text for line in text
        )
        self.assertNotIn("Diagnostic only - no refinement applied", unavailable_text)
        self.assertFalse(any("Method A low response" in line for line in unavailable_text))
        self.assertFalse(any("threshold unavailable" in line for line in unavailable_text))

    def test_d8_delta_renderer_uses_common_range_and_visible_t_header(self):
        _ProjectionHistogram.created = []
        root = _FakeROOT()
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "baseline_pion_mm": {"histogram": _ProjectionHistogram("template", (0.0, 0.0))},
            "delta_projection": {
                "rows_by_delta": (
                    ({"missing_mass": 0.2, "baseline_contribution": -5.0},),
                    ({"missing_mass": 1.2, "baseline_contribution": 9.0},),
                )
            },
        }

        self.assertTrue(
            plots._render_d8_delta_page(root, "unused.pdf", group, (-10.0, 0.0, 10.0))
        )
        displays = [
            histogram
            for histogram in _ProjectionHistogram.created
            if "_d8_delta_baseline_" in histogram.label
        ]
        self.assertEqual(len(displays), 2)
        self.assertEqual(
            {(histogram.display_minimum, histogram.display_maximum) for histogram in displays},
            {(-5.7, 9.7)},
        )
        self.assertTrue(
            any(
                text == ("Baseline pion background across delta", "|t| = [0.0000, 1.0000] GeV^2")
                for text in root.drawn_text
            )
        )

    def test_signed_display_y_range_covers_every_histogram_and_zero(self):
        raw = _BinnedHistogram("raw", (1.0, 2.0))
        proton = _BinnedHistogram("proton", (-8.0, -6.0))
        cleaned = _BinnedHistogram("cleaned", (10.0, 12.0))

        display_range = plots._combined_histogram_y_range((raw, proton, cleaned))
        self.assertLessEqual(display_range[0], -8.0)
        self.assertGreaterEqual(display_range[1], 12.0)
        self.assertLessEqual(display_range[0], 0.0)
        self.assertGreaterEqual(display_range[1], 0.0)
        self.assertEqual(
            display_range,
            plots._combined_histogram_y_range((raw, proton, cleaned)),
        )
        self.assertEqual(plots._combined_histogram_y_range(()), (-0.05, 0.05))
        self.assertEqual(
            plots._combined_histogram_y_range((_BinnedHistogram("zero", (0.0, 0.0)),)),
            (-0.05, 0.05),
        )

    def test_d7_overlays_apply_shared_signed_range_to_the_axis_clone(self):
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "raw_mm": {"histogram": _BinnedHistogram("raw", (1.0, 2.0))},
            "proton_mm": {"histogram": _BinnedHistogram("proton", (-8.0, -6.0))},
            "cleaned_mm": {"histogram": _BinnedHistogram("cleaned", (10.0, 12.0))},
        }
        root = _FakeROOT()
        for other_key, other_label, color, page_id, minimum, maximum in (
            ("proton_mm", "Proton contamination", "kRed", "full_background.d7.proton_mm", -8.0, 2.0),
            ("cleaned_mm", "Proton-cleaned kaon data", "kBlue", "full_background.d7.proton_cleaned_mm", 1.0, 12.0),
        ):
            _BinnedHistogram.created = []
            self.assertTrue(
                plots._render_d7_mm_overlay_page(
                    root,
                    "unused.pdf",
                    group,
                    other_key,
                    "overlay",
                    other_label,
                    color,
                    page_id,
                )
            )
            raw_display = next(
                histogram
                for histogram in _BinnedHistogram.created
                if "_raw_" in histogram.label
            )
            self.assertLessEqual(raw_display.display_minimum, minimum)
            self.assertGreaterEqual(raw_display.display_maximum, maximum)
            self.assertLessEqual(raw_display.display_minimum, 0.0)
            self.assertGreaterEqual(raw_display.display_maximum, 0.0)

    def test_d7_delta_page_uses_one_common_signed_range_and_visible_t_header(self):
        _BinnedHistogram.created = []
        root = _FakeROOT()
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "raw_mm": {"histogram": _BinnedHistogram("template", (0.0, 0.0))},
            "delta_projection": {
                "rows_by_delta": (
                    (
                        {
                            "missing_mass": 0.2,
                            "raw_weight": 1.0,
                            "proton_contribution": -7.0,
                            "cleaned_contribution": 5.0,
                        },
                    ),
                    (
                        {
                            "missing_mass": 1.2,
                            "raw_weight": 2.0,
                            "proton_contribution": -3.0,
                            "cleaned_contribution": 11.0,
                        },
                    ),
                )
            },
        }

        self.assertTrue(
            plots._render_d7_delta_page(root, "unused.pdf", group, (-10.0, 0.0, 10.0))
        )
        raw_displays = [
            histogram
            for histogram in _BinnedHistogram.created
            if "_delta_raw_" in histogram.label
        ]
        self.assertEqual(len(raw_displays), 2)
        self.assertEqual(
            {(histogram.display_minimum, histogram.display_maximum) for histogram in raw_displays},
            {(-7.9, 11.9)},
        )
        self.assertTrue(
            any(
                text == ("Proton subtraction across delta", "|t| = [0.0000, 1.0000] GeV^2")
                for text in root.drawn_text
            )
        )

    def test_d6_timing_page_draws_visible_t_header(self):
        root = _FakeROOT()
        group = {
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "pid": {
                "timing_axis_label": "RF timing [ns]",
                "cell_histograms": (
                    _BinnedHistogram("timing-d0", (1.0, 2.0)),
                    _BinnedHistogram("timing-d1", (2.0, 3.0)),
                ),
            },
        }

        self.assertTrue(
            plots._render_timing_t_pid_page(root, "unused.pdf", group, (-10.0, 0.0, 10.0))
        )
        self.assertTrue(
            any(
                text == ("Proton-identification timing", "|t| = [0.0000, 1.0000] GeV^2")
                for text in root.drawn_text
            )
        )

    def test_clone_is_detached_without_source_mutation(self):
        source = _Histogram("source")
        clone = plots._clone_display_histogram(source, "display")
        self.assertIsNot(clone, source)
        self.assertEqual(source.directory != 0, True)
        self.assertEqual(clone.directory, 0)
        self.assertEqual(source.title, "original")

    def test_ctime_aerogel_colz_pages_define_their_z_axis_quantities(self):
        source = (REPO_ROOT / "src" / "cuts" / "full_background_subtraction_plots.py").read_text(encoding="utf-8")
        pid_start = source.index("def _render_ctime_aero_pid_page")
        weight_start = source.index("def _render_ctime_aero_weight_page")
        pid_source = source[pid_start:source.index("def _render_timing_t_weight_page", pid_start)]
        weight_source = source[weight_start:source.index("def open_full_background_subtraction_pdf", weight_start)]

        self.assertIn('histogram.Draw("colz")', pid_source)
        self.assertIn("Aerogel NPE;Coincidence time [ns];Signed weighted yield", pid_source)
        self.assertIn('histogram.Draw("colz")', weight_source)
        self.assertIn("delta [%];Aerogel NPE;Proton contamination weight", weight_source)

    def test_static_presentation_and_runtime_contracts(self):
        source = (REPO_ROOT / "src" / "cuts" / "full_background_subtraction_plots.py").read_text(encoding="utf-8")
        for heading in (
            "Kaon-selected missing mass",
            "Proton-identification timing",
            "Proton contamination weight",
            "Proton contamination in missing mass",
            "Before and after proton subtraction",
            "Proton subtraction across delta",
            "Baseline pion background in missing mass",
            "Before and after baseline pion subtraction",
            "Baseline pion background across delta",
            "HGCer response",
            "HGCer response across delta",
            "Method A - HGCer response diagnostic",
            "Method B - Missing-mass closure inputs",
            "Method B - Local missing-mass closure",
            "Method B - Missing-mass closure diagnostic",
            "Method A / Method B diagnostic comparison",
            "Method A / Method B relative comparison",
            "Method availability across delta",
        ):
            with self.subTest(heading=heading):
                self.assertIn(heading, source)
        for notice in (
            "Diagnostic only - no refinement applied",
            "Diagnostic only - no correction or method selection",
            "NON-AUTHORITATIVE DIAGNOSTIC",
            "No refinement, correction, or method selection",
        ):
            with self.subTest(notice=notice):
                self.assertNotIn(notice, source)
        technical_source = (
            REPO_ROOT / "src" / "cuts" / "pion_hgcer_refinement_plots.py"
        ).read_text(encoding="utf-8")
        for technical_notice in (
            "NON-AUTHORITATIVE DIAGNOSTIC / No refinement applied",
            "NON-AUTHORITATIVE DIAGNOSTIC",
        ):
            with self.subTest(technical_notice=technical_notice):
                self.assertIn(technical_notice, technical_source)
        for forbidden_phase_object in (
            "D12_PRESENTATION_SCHEMA_VERSION",
            "build_full_background_subtraction_d12_payload",
            "full_background.d12.",
            "D12_1_PRESENTATION_SCHEMA_VERSION",
            "build_full_background_subtraction_d12_1_payload",
            "full_background.d12_1.",
        ):
            with self.subTest(forbidden_phase_object=forbidden_phase_object):
                self.assertNotIn(forbidden_phase_object, source)
        open_helper_source = source[
            source.index("def open_full_background_subtraction_pdf"):
            source.index("def close_full_background_subtraction_pdf")
        ]
        close_helper_source = source[
            source.index("def close_full_background_subtraction_pdf"):
            source.index("def _render_d6_t_pages")
        ]
        self.assertIn('canvas.Print("{}[".format(pdf_name))', open_helper_source)
        self.assertIn('canvas.Print("{}]".format(pdf_name))', close_helper_source)
        self.assertNotIn("canvas.Print(pdf_name)", open_helper_source)
        self.assertNotIn("canvas.Print(pdf_name)", close_helper_source)
        pre_d9_source = source[:source.index("def _d9_integer")]
        for forbidden in (
            "import proton_contamination_weights",
            "build_kaon_proton_cleaning_result",
            "apply_kaon_proton_cleaning_to_targets",
            "_evaluate_event_proton_probability",
            "_build_timing_constraint_for_cell",
            "_fit_delta_timing_slice",
            "_resolve_post_proton_rf_application",
            "apply_low_epsilon_rf_after_proton_cleaning",
            "Rebin(",
            "Scale(",
            "Add(",
            "Multiply(",
            "SetBinContent(",
            "SetBinError(",
            "final_cleaned_factor",
            'get("sources")',
            "build_setting_t_bin_pion_parents",
            "build_particle_subtraction_component_result",
            "resolve_frozen_parent_application_policy",
            "evaluate_particle_subtraction_component_fit_result",
            "resolve_pion_hgcer_delta_edges",
            "build_pion_hgcer_method_a",
            "build_pion_hgcer_method_b",
            "single_scale",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, pre_d9_source)
        fresh_delta_histogram = source[
            source.index("def _new_d7_delta_histogram"):
            source.index("def _render_d7_delta_page")
        ]
        self.assertIn('histogram.Reset("ICES")', fresh_delta_histogram)
        for visible_internal in ("raw_targets", "prepared_event_weight_lookup", "support code"):
            self.assertNotIn(visible_internal, source[source.index("def _render_raw_mm_page"):])
        self.assertIn("from canonical_binning import find_canonical_bin", source)
        legacy_membership = source[
            source.index("def _d7_legacy_t_index"):
            source.index("def _build_d7_delta_projection")
        ]
        self.assertIn("find_canonical_bin", legacy_membership)
        d8_builder = source[
            source.index("def _build_d8_delta_projection"):
            source.index("def _d9_integer")
        ]
        self.assertIn("final_diagnostic_application_result", d8_builder)
        self.assertNotIn("proposed_diagnostic_application_result", d8_builder)
        self.assertIn("simc_shape_pion_weight_from_value", d8_builder)
        self.assertIn("invalid_frozen_delta_index", d8_builder)
        self.assertIn("H_pion_subtraction_template_MM_nosub", d8_builder)
        d9_source = source[
            source.index("def _d9_integer"):
            source.index("def _append_d7_exclusion_failure")
        ]
        for required in (
            "H_hgcer_kaon_weighted",
            "H_hgcer_pion_weighted",
            "H_hgcer_vs_delta_kaon_weighted",
            "H_hgcer_vs_delta_pion_weighted",
            "diagnostic_weight",
            "method_a_comparison_candidate",
            "TGraphAsymmErrors",
            "Signed weighted yield",
        ):
            self.assertIn(required, d9_source)
        d9_relative_renderer = d9_source[
            d9_source.index("def _render_d9_method_a_relative_page"):
            d9_source.index("def _render_d9_t_pages")
        ]
        self.assertIn("Method A - HGCer response diagnostic", d9_relative_renderer)
        self.assertIn("Relative pion-background diagnostic", d9_relative_renderer)
        self.assertNotIn(
            "Diagnostic only - no correction or method selection", d9_relative_renderer
        )
        self.assertNotIn("Method A cell / same-|t| parent", d9_relative_renderer)
        self.assertNotIn(
            "Stored Method A values only - no refinement applied", d9_relative_renderer
        )
        self.assertNotIn("Method A relative diagnostic", d9_relative_renderer)
        for forbidden in (
            "Method B",
            "ln(B/A)",
            "build_pion_hgcer_tdelta_diagnostic(",
            "resolve_pion_hgcer_delta_edges(",
            "build_pion_hgcer_method_a(",
            "build_pion_hgcer_method_a_comparison(",
            "build_pion_hgcer_method_b(",
            "build_pion_hgcer_method_b_comparison(",
            "build_pion_hgcer_ab_comparison(",
            "Rebin(",
            "Scale(",
            "Add(",
            "Multiply(",
            "SetBinContent(",
            "SetBinError(",
            "Normalized pion-control yield",
            "C_A",
            "C_B",
            "C_final",
            "use_A",
            "use_B",
            "combine_AB",
            "preferred_method",
            "selected_method",
        ):
            with self.subTest(d9_forbidden=forbidden):
                self.assertNotIn(forbidden, d9_source)

        d10_source = source[
            source.index("def _d10_integer"):
            source.index("def render_full_background_subtraction_d6_pages")
        ]
        for required in (
            "signed_host_event_contribution",
            "signed_baseline_event_contribution",
            "method_b_comparison_candidate",
            "Missing-mass closure inputs",
            "Local missing-mass closure",
            "Relative pion-background diagnostic",
            "KLambdaSigma0",
            "Pion-sensitive region (outside protected window)",
            "Lambda/Sigma protected region (1.10-1.23 GeV)",
        ):
            self.assertIn(required, d10_source)
        for forbidden in (
            "find_canonical_bin(",
            "build_pion_hgcer_event_contract(",
            "build_pion_hgcer_method_b(",
            "build_pion_hgcer_method_b_comparison(",
            "resolve_pion_hgcer_method_b_config(",
            "method_b_display_payload(",
            "phase_d_ab_plot_payload(",
            "candidate_L_B",
            "raw_ratio",
            "parent_reference_ratio",
            "build_pion_hgcer_method_a(",
            "build_pion_hgcer_ab_comparison(",
            "C_A",
            "C_B",
            "C_final",
            "use_A",
            "use_B",
            "combine_AB",
            "preferred_method",
        ):
            with self.subTest(d10_forbidden=forbidden):
                self.assertNotIn(forbidden, d10_source)
        d10_event_rows_source = source[
            source.index("def _d10_event_rows"):
            source.index("def _d10_copy_cell_bins")
        ]
        self.assertIn("raw_t_index is None", d10_event_rows_source)
        self.assertIn("raw_delta_index is None", d10_event_rows_source)
        self.assertNotIn("analysis_abs_t", d10_event_rows_source)
        self.assertNotIn("SHMS_delta", d10_event_rows_source)

        d11_labels_source = source[
            source.index("_D11_AVAILABILITY_LABELS"):
            source.index("def _import_root")
        ]
        for label in (
            "Both methods available",
            "Both present; ratio undefined",
            "Method A only",
            "Method B only",
            "Neither method available",
        ):
            with self.subTest(d11_availability_label=label):
                self.assertIn(label, d11_labels_source)
        d11_source = source[
            source.index("def _d11_nonempty_string"):
            source.index("def render_full_background_subtraction_procedure_pages")
        ]
        for required in (
            "D11_PRESENTATION_SCHEMA_VERSION",
            "ratio_B_over_A",
            "log_ratio_B_over_A",
            "TGraphAsymmErrors",
            "TGraphErrors",
            "TGraph",
            "TBox",
            "full_background.d11.ab_overlay",
            "full_background.d11.ab_ratio_log",
            "full_background.d11.method_availability",
        ):
            with self.subTest(d11_required=required):
                self.assertIn(required, d11_source)
        for forbidden in (
            "diagnostic_interval_relation",
            "build_pion_hgcer_comparison_input_contract(",
            "build_pion_hgcer_method_a_comparison(",
            "build_pion_hgcer_method_b_comparison(",
            "build_pion_hgcer_ab_comparison(",
            "phase_d_ab_plot_payload(",
            "find_canonical_bin(",
            "math.log(",
            "C_A",
            "C_B",
            "C_final",
            "use_A",
            "use_B",
            "combine_AB",
            "preferred_method",
            "selected_method",
        ):
            with self.subTest(d11_forbidden=forbidden):
                self.assertNotIn(forbidden, d11_source)
        d11_renderer_source = d11_source[
            d11_source.index("def _render_d11_ab_overlay_page"):
            d11_source.index("def _render_d11_t_pages")
        ]
        for notice in (
            "NON-AUTHORITATIVE DIAGNOSTIC",
            "No refinement, correction, or method selection",
        ):
            with self.subTest(d11_banner_absent=notice):
                self.assertNotIn(notice, d11_renderer_source)
        for forbidden in ("tension", "compatibility", "significance"):
            with self.subTest(d11_renderer_forbidden=forbidden):
                self.assertNotIn(forbidden, d11_renderer_source)
        cumulative_renderer_source = source[
            source.index("def render_full_background_subtraction_procedure_pages"):
            source.index("__all__", source.index("def render_full_background_subtraction_procedure_pages"))
        ]
        for forbidden in (
            "C_A", "C_B", "C_final", "use_A", "use_B", "combine_AB",
            "preferred_method", "selected_method", "chosen_method",
            "refined_pion_weight", "applied_refinement_weight",
            "tension", "compatibility", "significance", "chi2", "p_value",
        ):
            with self.subTest(cumulative_renderer_forbidden=forbidden):
                self.assertNotIn(forbidden, cumulative_renderer_source)

        runtime = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        start = runtime.index("# Phases D.6 through D.11 are terminal presentation only.")
        end = runtime.index("for supplement_key, role in (", start)
        block = runtime[start:end]
        for name in (
            "build_full_background_subtraction_d6_payload",
            "build_full_background_subtraction_d7_payload",
            "build_full_background_subtraction_d8_payload",
            "build_full_background_subtraction_d9_payload",
            "build_full_background_subtraction_d10_payload",
            "build_full_background_subtraction_d11_payload",
            "full_background_subtraction_pdf_path",
            "open_full_background_subtraction_pdf",
            "render_full_background_subtraction_procedure_pages",
            "close_full_background_subtraction_pdf",
        ):
            self.assertIn(name, block)
        for call in (
            "full_background_subtraction_pdf_path",
            "open_full_background_subtraction_pdf",
            "render_full_background_subtraction_procedure_pages",
            "close_full_background_subtraction_pdf",
        ):
            with self.subTest(lifecycle_call=call):
                self.assertEqual(block.count(call + "("), 1)
        self.assertNotIn("open_diagnostic_pdf", block)
        for forbidden in (
            "C_A", "C_B", "C_final", "use_A", "use_B", "combine_AB",
            "preferred_method", "selected_method", "chosen_method",
            "refined_pion_weight", "applied_refinement_weight",
            "tension", "compatibility", "significance", "chi2", "p_value",
        ):
            self.assertNotIn(forbidden, block)
        for forbidden_phase_object in (
            "D12_PRESENTATION_SCHEMA_VERSION",
            "build_full_background_subtraction_d12_payload",
            "full_background.d12.",
            "D12_1_PRESENTATION_SCHEMA_VERSION",
            "build_full_background_subtraction_d12_1_payload",
            "full_background.d12_1.",
        ):
            with self.subTest(runtime_forbidden_phase_object=forbidden_phase_object):
                self.assertNotIn(forbidden_phase_object, runtime)
        self.assertLess(
            runtime.index("build_full_background_subtraction_d6_payload(", start),
            runtime.index("build_full_background_subtraction_d7_payload(", start),
        )
        self.assertLess(
            runtime.index("build_full_background_subtraction_d7_payload(", start),
            runtime.index("build_full_background_subtraction_d8_payload(", start),
        )
        self.assertLess(
            runtime.index("build_full_background_subtraction_d8_payload(", start),
            runtime.index("build_full_background_subtraction_d9_payload(", start),
        )
        self.assertLess(
            runtime.index("build_full_background_subtraction_d9_payload(", start),
            runtime.index("build_full_background_subtraction_d10_payload(", start),
        )
        self.assertLess(
            runtime.index("build_full_background_subtraction_d10_payload(", start),
            runtime.index("build_full_background_subtraction_d11_payload(", start),
        )
        self.assertLess(
            runtime.index("build_full_background_subtraction_d11_payload(", start),
            runtime.index("render_full_background_subtraction_procedure_pages(", start),
        )
        self.assertLess(
            runtime.index("full_background_subtraction_pdf_path(", start),
            runtime.index("open_full_background_subtraction_pdf(", start),
        )
        self.assertLess(
            runtime.index("open_full_background_subtraction_pdf(", start),
            runtime.index("render_full_background_subtraction_procedure_pages(", start),
        )
        self.assertLess(
            runtime.index("render_full_background_subtraction_procedure_pages(", start),
            runtime.index("close_full_background_subtraction_pdf(", start),
        )
        self.assertLess(
            runtime.index("proton_cleaning_application = apply_kaon_proton_cleaning_to_targets("),
            start,
        )
        self.assertNotIn("full_background_subtraction_d6_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d7_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d8_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d9_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d10_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d11_payload\"]", block)
        self.assertNotIn("full_background_subtraction_d12_payload\"]", block)
        self.assertEqual(block.count('histDict["full_background_subtraction_'), 2)
        self.assertIn('histDict["full_background_subtraction_page_manifest"]', block)
        self.assertIn('histDict["full_background_subtraction_renderer_failures"]', block)

    @unittest.skipUnless(plots._import_root() is not None, "PyROOT not available")
    def test_d10_root_histograms_preserve_frozen_boundaries_and_stored_bins(self):
        """D.10.1: match frozen Method-B regular-bin ownership exactly."""
        ROOT = plots._import_root()
        ROOT.gROOT.SetBatch(True)
        edges = [0.80, 0.91, 1.03, 1.17, 1.35]
        event_histogram = plots._d10_event_histogram(
            ROOT,
            "H_d10_native_boundary",
            edges,
            (
                {"missing_mass": 1.35, "signed_contribution": 6.25},
                {"missing_mass": 1.36, "signed_contribution": -2.0},
                {"missing_mass": 0.80, "signed_contribution": -3.0},
            ),
        )
        self.assertIsNotNone(event_histogram)
        self.assertEqual(event_histogram.GetNbinsX(), 4)
        self.assertEqual(
            [event_histogram.GetXaxis().GetBinLowEdge(index) for index in range(1, 6)],
            edges,
        )
        self.assertAlmostEqual(event_histogram.GetBinContent(1), -3.0)
        self.assertAlmostEqual(event_histogram.GetBinError(1), 3.0)
        self.assertAlmostEqual(event_histogram.GetBinContent(4), 6.25)
        self.assertAlmostEqual(event_histogram.GetBinError(4), 6.25)
        self.assertAlmostEqual(event_histogram.GetBinContent(5), -2.0)
        self.assertAlmostEqual(event_histogram.GetBinError(5), 2.0)

        _phase_a, method_b, _comparison = _d10_fixture()
        host, baseline = plots._d10_local_histograms(
            ROOT, "H_d10_stored_local", method_b["cells"][0]
        )
        self.assertIsNotNone(host)
        self.assertIsNotNone(baseline)
        self.assertEqual(host.GetNbinsX(), 4)
        self.assertAlmostEqual(host.GetBinContent(1), 7.0)
        self.assertAlmostEqual(host.GetBinError(1), 3.0)
        self.assertAlmostEqual(baseline.GetBinContent(1), -2.0)
        self.assertAlmostEqual(baseline.GetBinError(1), 2.0)

    @unittest.skipUnless(plots._import_root() is not None, "PyROOT not available")
    def test_d11_root_graphs_boxes_and_stored_values_are_preserved(self):
        """D.11.1: real ROOT receives the frozen A/B scalar presentation unchanged."""
        ROOT = plots._import_root()
        ROOT.gROOT.SetBatch(True)
        checkpoint = _d11_fixture()
        checkpoint_before = deepcopy(checkpoint)
        payload = plots.build_full_background_subtraction_d11_payload(checkpoint)
        payload_before = deepcopy(payload)
        self.assertTrue(payload["available"])
        group = payload["per_t"][0]
        capture = _D11ROOTCapture(ROOT)
        with tempfile.TemporaryDirectory() as temporary:
            pdf_name = str(Path(temporary) / "d11-root-regression.pdf")
            self.assertTrue(
                plots._render_d11_ab_overlay_page(capture, pdf_name, payload, group)
            )
            self.assertTrue(
                plots._render_d11_ratio_log_page(capture, pdf_name, payload, group)
            )
            self.assertTrue(
                plots._render_d11_availability_page(capture, pdf_name, payload, group)
            )

        self.assertEqual(len(capture.asymmetric_graphs), 1)
        self.assertEqual(len(capture.symmetric_graphs), 1)
        self.assertEqual(len(capture.central_graphs), 2)
        asymmetric = capture.asymmetric_graphs[0]
        zero_index = next(
            index for index in range(asymmetric.GetN())
            if asymmetric.GetPointY(index) == 0.0
        )
        self.assertEqual(asymmetric.GetPointY(zero_index), 0.0)
        self.assertEqual(asymmetric.GetErrorYlow(zero_index), 0.0)
        self.assertEqual(asymmetric.GetErrorYhigh(zero_index), 0.25)
        symmetric = capture.symmetric_graphs[0]
        method_b_index = next(
            index for index in range(symmetric.GetN())
            if symmetric.GetPointY(index) == 4.0
        )
        self.assertEqual(symmetric.GetErrorY(method_b_index), 0.5)
        self.assertEqual(
            [graph.GetPointY(0) for graph in capture.central_graphs], [7.25, -3.5]
        )
        self.assertIn(
            (-12.0, 1.0, 23.0, 1.0),
            [(line.GetX1(), line.GetY1(), line.GetX2(), line.GetY2()) for line in capture.lines],
        )
        self.assertIn(
            (-12.0, 0.0, 23.0, 0.0),
            [(line.GetX1(), line.GetY1(), line.GetX2(), line.GetY2()) for line in capture.lines],
        )
        main_boxes = [
            box for box in capture.boxes
            if box.GetY1() == 0.0 and box.GetY2() == 1.0
        ]
        self.assertEqual(len(main_boxes), len(group["availability"]["cells"]))
        self.assertEqual(
            [(box.GetX1(), box.GetX2()) for box in main_boxes],
            [(-12.0, -5.0), (-5.0, 1.0), (1.0, 6.0), (6.0, 14.0), (14.0, 23.0)],
        )
        self.assertEqual(checkpoint, checkpoint_before)
        self.assertEqual(payload, payload_before)

    @unittest.skipUnless(plots._import_root() is not None, "PyROOT not available")
    def test_root_rendering_preserves_source_and_page_order(self):
        ROOT = plots._import_root()
        ROOT.gROOT.SetBatch(True)
        raw = ROOT.TH1D("H_d6_raw", "source", 3, 0.0, 3.0)
        raw.SetDirectory(0)
        raw.SetBinContent(1, 4.0)
        raw.SetBinError(1, 0.4)
        proton = ROOT.TH1D("H_d7_proton", "source", 3, 0.0, 3.0)
        proton.SetDirectory(0)
        proton.SetBinContent(1, 1.0)
        proton.SetBinError(1, 0.2)
        cleaned = ROOT.TH1D("H_d7_cleaned", "source", 3, 0.0, 3.0)
        cleaned.SetDirectory(0)
        cleaned.SetBinContent(1, 3.0)
        cleaned.SetBinError(1, 0.3)
        timing_cells = []
        for delta_index in range(2):
            row = []
            for t_index in range(2):
                histogram = ROOT.TH1D(
                    "H_d6_timing_{}_{}".format(delta_index, t_index), "source", 4, -2.0, 2.0
                )
                histogram.SetDirectory(0)
                histogram.SetBinContent(1, 1.0)
                histogram.SetBinError(1, 0.1 * (1 + delta_index + t_index))
                row.append(histogram)
            timing_cells.append(tuple(row))
        weight = ROOT.TH2D("H_d6_weight", "source", 2, -10.0, 10.0, 2, 0.0, 2.0)
        weight.SetDirectory(0)
        weight.SetBinContent(1, 1, 0.25)
        weight.SetBinError(1, 1, 0.025)
        result = _d7_timing_result()
        result["H_delta_timing_t_cells"] = timing_cells
        application = {
            "canonical_t_products": _products(
                (raw, raw),
                (proton, proton),
                (cleaned, cleaned),
            ),
            "H_proton_weight_vs_delta_t": weight,
        }
        d6_payload = plots.build_full_background_subtraction_d6_payload(result, application)
        d7_payload = plots.build_full_background_subtraction_d7_payload(
            result, application, _d7_prepared_bundle()
        )
        def d8_histogram(name, contents, errors=()):
            histogram = ROOT.TH1D(name, "source", 3, 0.0, 3.0)
            histogram.SetDirectory(0)
            histogram.Sumw2()
            for bin_index, content in enumerate(contents, start=1):
                histogram.SetBinContent(bin_index, float(content))
                if errors:
                    histogram.SetBinError(bin_index, float(errors[bin_index - 1]))
            return histogram

        d8_parents = []
        d8_cache_rows = []
        for t_index in range(2):
            template = d8_histogram(
                "H_d8_template_t{}".format(t_index),
                (2.0, 3.0, 0.0) if t_index == 0 else (0.0, 0.0, 0.0),
                (2.0, 3.0, 0.0) if t_index == 0 else (0.0, 0.0, 0.0),
            )
            before = d8_histogram(
                "H_d8_before_t{}".format(t_index), (4.0, 5.0, 1.0), (0.5, 0.75, 0.25)
            )
            after = d8_histogram(
                "H_d8_after_t{}".format(t_index), (4.0, 2.0, 1.0), (0.4, 0.6, 0.2)
            )
            model = d8_histogram("H_d8_model_t{}".format(t_index), (0.0, 0.0, 0.0))
            d8_parents.append({
                "t_bin_index": t_index,
                "t_edges": [float(t_index), float(t_index + 1)],
                "final_diagnostic_application_result": {
                    "H_MM_nosub_before_pion_subtraction": before,
                    "H_pion_subtraction_template_MM_nosub": template,
                    "H_MM_nosub_after_pion_subtraction": after,
                    "H_pion_control_model": model,
                    "weights": [0.0, 2.0, 3.0, 0.0],
                },
            })
            d8_cache_rows.append({
                "t_index": t_index,
                "t_edges": [float(t_index), float(t_index + 1)],
                "records": [
                    {
                        "nommcuts": True,
                        "delta_index": 0,
                        "adj_MM": 0.2,
                        "coefficient": 1.0 if t_index == 0 else 0.0,
                    },
                    {
                        "nommcuts": True,
                        "delta_index": 1,
                        "adj_MM": 1.2,
                        "coefficient": 1.0 if t_index == 0 else 0.0,
                    },
                ],
            })
        d8_payload = plots.build_full_background_subtraction_d8_payload(
            tuple(d8_parents), {"delta_edges": [-10.0, 0.0, 10.0], "by_t": tuple(d8_cache_rows)}
        )
        self.assertTrue(d8_payload["available"])
        self.assertEqual(
            [group["delta_projection"]["closure"]["status"] for group in d8_payload["per_t"]],
            ["closed", "closed"],
        )
        d9_diagnostic, d9_method_a, d9_comparison = _d9_fixture()
        d9_source_snapshots = []
        for side in ("kaon", "pion"):
            response_histogram = ROOT.TH1D(
                "H_d9_{}_response".format(side), "source", 4, 0.0, 4.0
            )
            response_histogram.SetDirectory(0)
            response_histogram.SetBinContent(1, 7.0 if side == "kaon" else -5.0)
            response_histogram.SetBinError(1, 0.7 if side == "kaon" else 0.5)
            delta_histogram = ROOT.TH2D(
                "H_d9_{}_delta".format(side), "source", 2, -10.0, 10.0, 4, 0.0, 4.0
            )
            delta_histogram.SetDirectory(0)
            delta_histogram.SetBinContent(1, 1, 3.0 if side == "kaon" else -2.0)
            delta_histogram.SetBinError(1, 1, 0.3 if side == "kaon" else 0.2)
            d9_diagnostic["histograms"]["H_hgcer_{}_weighted".format(side)] = response_histogram
            d9_diagnostic["histograms"]["H_hgcer_vs_delta_{}_weighted".format(side)] = delta_histogram
            d9_source_snapshots.append((
                response_histogram,
                _snapshot_root_th1(response_histogram),
            ))
            d9_source_snapshots.append((
                delta_histogram,
                _snapshot_root_th2(delta_histogram),
            ))
        d9_payload = plots.build_full_background_subtraction_d9_payload(
            d9_diagnostic, d9_method_a, d9_comparison
        )
        self.assertTrue(d9_payload["available"])
        d10_phase_a, d10_method_b, d10_comparison = _d10_fixture()
        d10_sources_before = deepcopy((d10_phase_a, d10_method_b, d10_comparison))
        d10_payload = plots.build_full_background_subtraction_d10_payload(
            d10_phase_a, d10_method_b, d10_comparison
        )
        self.assertTrue(d10_payload["available"])
        d11_checkpoint = _d11_fixture(
            t_edges=(0.0, 1.0, 2.0),
            delta_edges=(-10.0, 0.0, 10.0),
        )
        d11_checkpoint_before = deepcopy(d11_checkpoint)
        d11_payload = plots.build_full_background_subtraction_d11_payload(
            d11_checkpoint
        )
        self.assertTrue(d11_payload["available"])
        d6_d7_source_snapshots = tuple(
            (histogram, _snapshot_root_th1(histogram))
            for histogram in (raw, proton, cleaned) + tuple(
                histogram for row in timing_cells for histogram in row
            )
        )
        d6_weight_snapshot = _snapshot_root_th2(weight)
        d8_source_snapshots = []
        for parent in d8_parents:
            application = parent["final_diagnostic_application_result"]
            for key in (
                "H_MM_nosub_before_pion_subtraction",
                "H_pion_subtraction_template_MM_nosub",
                "H_MM_nosub_after_pion_subtraction",
            ):
                histogram = application[key]
                axis = histogram.GetXaxis()
                d8_source_snapshots.append((
                    histogram,
                    histogram.GetNbinsX(),
                    tuple(
                        (histogram.GetBinContent(index), histogram.GetBinError(index))
                        for index in range(1, histogram.GetNbinsX() + 1)
                    ),
                    axis.GetXmin(),
                    axis.GetXmax(),
                ))
        with tempfile.TemporaryDirectory() as temporary:
            pdf = str(Path(temporary) / "procedure.pdf")
            self.assertTrue(plots.open_full_background_subtraction_pdf(pdf))
            try:
                rendered = plots.render_full_background_subtraction_procedure_pages(
                    pdf, d6_payload, d7_payload, d8_payload, d9_payload, d10_payload,
                    d11_payload,
                )
            finally:
                plots.close_full_background_subtraction_pdf(pdf)
            self.assertTrue(Path(pdf).exists())
            self.assertGreater(Path(pdf).stat().st_size, 0)
            self.assertEqual(
                [page["page_id"] for page in rendered["manifest"]],
                list(EXPECTED_FULL_BACKGROUND_PAGE_IDS) * 2,
            )
            self.assertEqual(len(rendered["manifest"]), 36)
            self.assertEqual(
                [page["scope"] for page in rendered["manifest"]],
                ["t1"] * 18 + ["t2"] * 18,
            )
            self.assertTrue(all(page["authoritative"] is False for page in rendered["manifest"]))
            _assert_full_background_manifest_contract(self, rendered["manifest"])
            pdfinfo = shutil.which("pdfinfo")
            if pdfinfo is not None:
                pdfinfo_result = subprocess.run(
                    [pdfinfo, pdf],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(pdfinfo_result.returncode, 0, pdfinfo_result.stderr)
                page_count_line = next(
                    (
                        line for line in pdfinfo_result.stdout.splitlines()
                        if line.startswith("Pages:")
                    ),
                    None,
                )
                self.assertIsNotNone(page_count_line, pdfinfo_result.stdout)
                self.assertEqual(
                    int(page_count_line.split(":", 1)[1].strip()),
                    len(rendered["manifest"]),
                )
        for histogram, snapshot in d6_d7_source_snapshots:
            self.assertEqual(_snapshot_root_th1(histogram), snapshot)
        self.assertEqual(_snapshot_root_th2(weight), d6_weight_snapshot)
        for histogram, snapshot in d9_source_snapshots:
            if histogram.InheritsFrom("TH2"):
                self.assertEqual(_snapshot_root_th2(histogram), snapshot)
            else:
                self.assertEqual(_snapshot_root_th1(histogram), snapshot)
        for histogram, bin_count, bins, x_minimum, x_maximum in d8_source_snapshots:
            self.assertEqual(histogram.GetNbinsX(), bin_count)
            self.assertEqual(
                tuple(
                    (histogram.GetBinContent(index), histogram.GetBinError(index))
                    for index in range(1, histogram.GetNbinsX() + 1)
                ),
                bins,
            )
            self.assertEqual(histogram.GetXaxis().GetXmin(), x_minimum)
            self.assertEqual(histogram.GetXaxis().GetXmax(), x_maximum)
        self.assertEqual((d10_phase_a, d10_method_b, d10_comparison), d10_sources_before)
        self.assertEqual(d11_checkpoint, d11_checkpoint_before)


if __name__ == "__main__":
    unittest.main()
