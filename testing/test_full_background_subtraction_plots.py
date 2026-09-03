"""Focused D.6 procedure-PDF contract tests."""

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
            "Baseline pion background",
            "HGCer response",
            "Method A",
            "Method B",
            "ln(B/A)",
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

        runtime = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        start = runtime.index("# Phases D.6/D.7 are terminal presentation only.")
        end = runtime.index("for supplement_key, role in (", start)
        block = runtime[start:end]
        for name in (
            "build_full_background_subtraction_d6_payload",
            "build_full_background_subtraction_d7_payload",
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
        with tempfile.TemporaryDirectory() as temporary:
            pdf = str(Path(temporary) / "procedure.pdf")
            self.assertTrue(plots.open_full_background_subtraction_pdf(pdf))
            try:
                rendered = plots.render_full_background_subtraction_procedure_pages(
                    pdf, d6_payload, d7_payload
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
                    "full_background.d6.raw_mm",
                    "full_background.d6.proton_pid",
                    "full_background.d6.proton_weight",
                    "full_background.d7.proton_mm",
                    "full_background.d7.proton_cleaned_mm",
                    "full_background.d7.proton_delta_mm",
                ],
            )
            self.assertEqual(
                [page["scope"] for page in rendered["manifest"]],
                ["t1"] * 6 + ["t2"] * 6,
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


if __name__ == "__main__":
    unittest.main()
