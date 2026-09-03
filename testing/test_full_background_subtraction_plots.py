"""Focused D.6 through D.9 procedure-PDF contract tests."""

from __future__ import annotations

import tempfile
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import full_background_subtraction_plots as plots


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

    def SetBorderSize(self, _size):
        return None

    def SetFillStyle(self, _style):
        return None

    def AddEntry(self, *_args):
        return None

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


class FullBackgroundSubtractionD6Tests(unittest.TestCase):
    def test_pdf_path_is_deterministic_and_detached(self):
        main = r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe.pdf"
        path = plots.full_background_subtraction_pdf_path(main)
        self.assertEqual(
            path,
            r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction.pdf",
        )
        self.assertNotEqual(path, main)
        self.assertNotIn("hgcer-ab-comparison", path)

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

    def test_d9_response_page_uses_detached_signed_weighted_clones_and_notices(self):
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
        self.assertIn("Diagnostic only - no refinement applied", visible_text)
        self.assertIn("Method A low response: 0 < HGCer NPE <= 2", visible_text)
        self.assertEqual(len(root.legends), 1)
        self.assertEqual(len(root.pave_texts), 2)
        legend_rectangle = root.legends[0].coordinates
        authority_rectangle = root.pave_texts[0].coordinates
        threshold_rectangle = root.pave_texts[1].coordinates
        self.assertEqual(legend_rectangle, (0.62, 0.58, 0.89, 0.72))
        self.assertFalse(_rectangles_overlap(legend_rectangle, authority_rectangle))
        self.assertFalse(_rectangles_overlap(legend_rectangle, threshold_rectangle))
        self.assertFalse(_rectangles_overlap(authority_rectangle, threshold_rectangle))

        unavailable_root = _FakeROOT()
        self.assertTrue(
            plots._render_d9_hgcer_response_page(
                unavailable_root, "unused.pdf", group, {"available": False}
            )
        )
        unavailable_text = tuple(
            line for text in unavailable_root.drawn_text for line in text
        )
        self.assertIn("Diagnostic only - no refinement applied", unavailable_text)
        self.assertNotIn("Method A low response: 0 < HGCer NPE <= 2", unavailable_text)
        self.assertFalse(any("threshold unavailable" in line for line in unavailable_text))
        self.assertEqual(len(unavailable_root.legends), 1)
        self.assertEqual(len(unavailable_root.pave_texts), 1)
        self.assertFalse(
            _rectangles_overlap(
                unavailable_root.legends[0].coordinates,
                unavailable_root.pave_texts[0].coordinates,
            )
        )

    def test_d9_delta_page_keeps_authority_and_threshold_notices_separate(self):
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
        self.assertIn("Diagnostic only - no refinement applied", visible_text)
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
        self.assertIn("Diagnostic only - no refinement applied", unavailable_text)
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
        self.assertIn(
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
        ):
            with self.subTest(d9_forbidden=forbidden):
                self.assertNotIn(forbidden, d9_source)

        runtime = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        start = runtime.index("# Phases D.6/D.7/D.8/D.9 are terminal presentation only.")
        end = runtime.index("for supplement_key, role in (", start)
        block = runtime[start:end]
        for name in (
            "build_full_background_subtraction_d6_payload",
            "build_full_background_subtraction_d7_payload",
            "build_full_background_subtraction_d8_payload",
            "build_full_background_subtraction_d9_payload",
            "full_background_subtraction_pdf_path",
            "open_full_background_subtraction_pdf",
            "render_full_background_subtraction_procedure_pages",
            "close_full_background_subtraction_pdf",
        ):
            self.assertIn(name, block)
        self.assertNotIn("open_diagnostic_pdf", block)
        for forbidden in ("C_A", "C_B", "C_final", "use_A", "use_B", "combine_AB"):
            self.assertNotIn(forbidden, block)
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

    @unittest.skipUnless(plots._import_root() is not None, "PyROOT not available")
    def test_root_rendering_preserves_source_and_page_order(self):
        ROOT = plots._import_root()
        ROOT.gROOT.SetBatch(True)
        raw = ROOT.TH1D("H_d6_raw", "source", 3, 0.0, 3.0)
        raw.SetDirectory(0)
        raw.SetBinContent(1, 4.0)
        proton = ROOT.TH1D("H_d7_proton", "source", 3, 0.0, 3.0)
        proton.SetDirectory(0)
        proton.SetBinContent(1, 1.0)
        cleaned = ROOT.TH1D("H_d7_cleaned", "source", 3, 0.0, 3.0)
        cleaned.SetDirectory(0)
        cleaned.SetBinContent(1, 3.0)
        timing_cells = []
        for delta_index in range(2):
            row = []
            for t_index in range(2):
                histogram = ROOT.TH1D(
                    "H_d6_timing_{}_{}".format(delta_index, t_index), "source", 4, -2.0, 2.0
                )
                histogram.SetDirectory(0)
                histogram.SetBinContent(1, 1.0)
                row.append(histogram)
            timing_cells.append(tuple(row))
        weight = ROOT.TH2D("H_d6_weight", "source", 2, -10.0, 10.0, 2, 0.0, 2.0)
        weight.SetDirectory(0)
        weight.SetBinContent(1, 1, 0.25)
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
            delta_histogram = ROOT.TH2D(
                "H_d9_{}_delta".format(side), "source", 2, -10.0, 10.0, 4, 0.0, 4.0
            )
            delta_histogram.SetDirectory(0)
            delta_histogram.SetBinContent(1, 1, 3.0 if side == "kaon" else -2.0)
            d9_diagnostic["histograms"]["H_hgcer_{}_weighted".format(side)] = response_histogram
            d9_diagnostic["histograms"]["H_hgcer_vs_delta_{}_weighted".format(side)] = delta_histogram
            d9_source_snapshots.append((
                response_histogram,
                tuple(response_histogram.GetBinContent(index) for index in range(1, 5)),
            ))
            d9_source_snapshots.append((
                delta_histogram,
                tuple(
                    delta_histogram.GetBinContent(x_index, y_index)
                    for x_index in range(1, 3)
                    for y_index in range(1, 5)
                ),
            ))
        d9_payload = plots.build_full_background_subtraction_d9_payload(
            d9_diagnostic, d9_method_a, d9_comparison
        )
        self.assertTrue(d9_payload["available"])
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
                    pdf, d6_payload, d7_payload, d8_payload, d9_payload
                )
            finally:
                plots.close_full_background_subtraction_pdf(pdf)
            self.assertEqual(
                [page["page_id"] for page in rendered["manifest"]],
                [
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
                ],
            )
            self.assertEqual(
                [page["scope"] for page in rendered["manifest"]],
                ["t1"] * 12 + ["t2"] * 12,
            )
            self.assertTrue(all(page["authoritative"] is False for page in rendered["manifest"]))
            self.assertTrue(Path(pdf).exists())
        self.assertEqual(raw.GetBinContent(1), 4.0)
        self.assertEqual(raw.GetNbinsX(), 3)
        self.assertEqual(proton.GetBinContent(1), 1.0)
        self.assertEqual(cleaned.GetBinContent(1), 3.0)
        self.assertEqual(timing_cells[0][0].GetNbinsX(), 4)
        self.assertEqual(weight.GetNbinsX(), 2)
        self.assertEqual(weight.GetNbinsY(), 2)
        for histogram, contents in d9_source_snapshots:
            if histogram.InheritsFrom("TH2"):
                self.assertEqual(histogram.GetNbinsX(), 2)
                self.assertEqual(histogram.GetNbinsY(), 4)
                self.assertEqual(
                    tuple(
                        histogram.GetBinContent(x_index, y_index)
                        for x_index in range(1, 3)
                        for y_index in range(1, 5)
                    ),
                    contents,
                )
            else:
                self.assertEqual(histogram.GetNbinsX(), 4)
                self.assertEqual(
                    tuple(histogram.GetBinContent(index) for index in range(1, 5)),
                    contents,
                )
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


if __name__ == "__main__":
    unittest.main()
