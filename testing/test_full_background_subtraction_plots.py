"""Focused D.6/D.7/D.8 procedure-PDF contract tests."""

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
    def SetBorderSize(self, _size):
        return None

    def SetFillStyle(self, _style):
        return None

    def AddEntry(self, *_args):
        return None

    def Draw(self):
        return None


class _FakePaveText:
    def __init__(self, root):
        self._root = root
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

    def TCanvas(self, *_args):
        return _FakeCanvas(*_args)

    def TLegend(self, *_args):
        return _FakeLegend()

    def TPaveText(self, *_args):
        return _FakePaveText(self)


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
            "HGCer response",
            "Method A",
            "Method B",
            "ln(B/A)",
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
                self.assertNotIn(forbidden, source)
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
            source.index("def full_background_subtraction_pdf_path")
        ]
        self.assertIn("final_diagnostic_application_result", d8_builder)
        self.assertNotIn("proposed_diagnostic_application_result", d8_builder)
        self.assertIn("simc_shape_pion_weight_from_value", d8_builder)
        self.assertIn("invalid_frozen_delta_index", d8_builder)
        self.assertIn("H_pion_subtraction_template_MM_nosub", d8_builder)

        runtime = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        start = runtime.index("# Phases D.6/D.7/D.8 are terminal presentation only.")
        end = runtime.index("for supplement_key, role in (", start)
        block = runtime[start:end]
        for name in (
            "build_full_background_subtraction_d6_payload",
            "build_full_background_subtraction_d7_payload",
            "build_full_background_subtraction_d8_payload",
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
                    pdf, d6_payload, d7_payload, d8_payload
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
                    "full_background.d6.raw_mm",
                    "full_background.d6.proton_pid",
                    "full_background.d6.proton_weight",
                    "full_background.d7.proton_mm",
                    "full_background.d7.proton_cleaned_mm",
                    "full_background.d7.proton_delta_mm",
                    "full_background.d8.pion_background_mm",
                    "full_background.d8.pion_subtracted_mm",
                    "full_background.d8.pion_delta_mm",
                ],
            )
            self.assertEqual(
                [page["scope"] for page in rendered["manifest"]],
                ["t1"] * 9 + ["t2"] * 9,
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
