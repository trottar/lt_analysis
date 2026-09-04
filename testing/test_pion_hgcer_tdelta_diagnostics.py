"""Contracts for the non-authoritative pre-HGCer Part-1 diagnostic."""

from __future__ import annotations

import gc
from copy import deepcopy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_diagnostics as diagnostics
from data_coordinates import build_kaon_data_coordinate_contract
from background_config import get_pion_hgcer_diagnostic_config


class _Event:
    def __init__(self, *, mm, mandel_t, delta, npe):
        self.MM = float(mm)
        self.MandelT = float(mandel_t)
        self.ssdelta = float(delta)
        self.P_hgcer_npeSum = float(npe)
        self.P_hgcer_xAtCer = 0.0
        self.P_hgcer_yAtCer = 0.0
        self.Q2 = 4.4
        self.W = 2.74
        self.epsilon = 0.4
        self.ph_q = 0.0


class PionHGCerTDeltaPurePythonTests(unittest.TestCase):
    def test_histogram_title_catalog_is_brace_safe_for_all_sides_and_weightings(self):
        expected_keys = {
            "hgcer",
            "delta",
            "t",
            "mm",
            "hgcer_vs_delta",
            "hgcer_vs_t",
            "mm_vs_delta",
            "mm_vs_hgcer",
            "support_absolute",
            "support_effective",
            "boundary_effective",
            "support_class",
            "boundary_readiness",
        }
        side_titles = {
            "kaon": "proton-cleaned kaon",
            "pion": "pion control",
        }
        weighting_titles = {
            "weighted": "signed weighted yield",
            "absolute": "absolute support",
        }

        for side, side_title in side_titles.items():
            for weighting, weighting_title in weighting_titles.items():
                with self.subTest(side=side, weighting=weighting):
                    titles = diagnostics._histogram_title_catalog(side, weighting)
                    self.assertEqual(set(titles), expected_keys)
                    for title in titles.values():
                        self.assertTrue(title)
                        self.assertIn(";", title)
                        self.assertNotIn("{", title)
                        self.assertNotIn("}", title)
                    self.assertIn(side_title, titles["hgcer"])
                    self.assertIn(weighting_title, titles["hgcer"])
                    self.assertIn("GeV^2", titles["hgcer_vs_t"])
                    self.assertIn("GeV^2", titles["support_absolute"])

    def test_boundary_bands_use_adjacent_half_open_intervals(self):
        boundary = {"threshold": 2.0, "boundary_band_width": 0.5}
        for value in (1.50, 1.99):
            self.assertTrue(diagnostics.pion_hgcer_boundary_contains("kaon", value, boundary))
            self.assertFalse(diagnostics.pion_hgcer_boundary_contains("pion", value, boundary))
        for value in (2.00, 2.49):
            self.assertFalse(diagnostics.pion_hgcer_boundary_contains("kaon", value, boundary))
            self.assertTrue(diagnostics.pion_hgcer_boundary_contains("pion", value, boundary))
        self.assertFalse(diagnostics.pion_hgcer_boundary_contains("kaon", 1.49, boundary))
        self.assertFalse(diagnostics.pion_hgcer_boundary_contains("pion", 2.50, boundary))

    def test_boundary_readiness_requires_both_pid_sides(self):
        thresholds = {
            "ready_absolute_weight": 0.1,
            "marginal_absolute_weight": 0.01,
            "ready_effective_entries": 5,
            "marginal_effective_entries": 1,
        }
        ready = {"absolute_weight_support": 0.2, "sum_weight_squared": 0.008}
        marginal = {"absolute_weight_support": 0.02, "sum_weight_squared": 0.0004}
        absent = {"absolute_weight_support": 0.0, "sum_weight_squared": 0.0}
        self.assertEqual(
            diagnostics.classify_pion_hgcer_boundary_readiness(ready, ready, thresholds),
            "ready",
        )
        self.assertEqual(
            diagnostics.classify_pion_hgcer_boundary_readiness(marginal, marginal, thresholds),
            "marginal",
        )
        self.assertEqual(
            diagnostics.classify_pion_hgcer_boundary_readiness(ready, absent, thresholds),
            "insufficient",
        )

    def test_root_facing_display_strings_are_ascii_and_not_mojibake(self):
        strings = [
            diagnostics.pion_hgcer_display_text(kind)
            for kind in (
                "part1",
                "part1p5",
                "part1p5_provenance",
                "part1p5_boundary",
                "part1p5_population",
                "part1p5_boundary_support",
                "part1p5_readiness",
                "part1p5_summary",
                "compact_cell_status",
                "legend_kaon",
                "legend_pion",
                "projection_npe_weighted",
                "projection_npe_absolute",
                "projection_mm_weighted",
                "projection_mm_absolute",
                "npe_vs_delta",
                "part1_t_note",
                "cell_no_response",
                "cell_gate_note",
                "part1p5_norf_policy",
                "part1p5_summary_note",
                "kaon_pid_mm",
                "pion_pid_mm",
            )
        ]
        strings.extend((
            diagnostics.pion_hgcer_display_text("part1p5_t_boundary", t_index=1),
            diagnostics.pion_hgcer_display_text(
                "part1p5_cell_boundary", t_index=1, delta_index=2
            ),
        ))
        for text in strings:
            with self.subTest(text=text):
                text.encode("ascii")
                for forbidden in ("â", "Â", "Ã", "–", "—"):
                    self.assertNotIn(forbidden, text)

    def test_no_rf_provenance_is_explicit_and_rejects_rf_sources(self):
        kaon_bundle = {
            "sources": {
                "prompt": {
                    "tree": (object(),),
                    "tree_name": "Cut_Kaon_Events_prompt_noRF",
                    "coefficient": 1.0,
                }
            }
        }
        pion_bundle = {
            "prompt_tree": (object(),),
            "prompt_tree_name": "Cut_Pion_Events_prompt_noRF",
        }
        provenance = diagnostics.resolve_pion_hgcer_source_provenance(
            kaon_bundle, pion_bundle
        )
        self.assertEqual(provenance["kaon"]["prompt"]["pid_role"], "kaon_pid")
        self.assertEqual(provenance["pion"]["prompt"]["rf_state"], "noRF")
        pion_bundle["prompt_tree_name"] = "Cut_Pion_Events_prompt_RF"
        with self.assertRaisesRegex(ValueError, "pion_hgcer_norf_provenance_failed"):
            diagnostics.resolve_pion_hgcer_source_provenance(kaon_bundle, pion_bundle)

    def test_final_edges_use_the_shared_canonical_membership_rule(self):
        self.assertEqual(diagnostics.canonical_t_delta_index(0.9, (0.4, 0.6, 0.9)), 1)
        self.assertEqual(diagnostics.canonical_t_delta_index(20.0, (-10.0, 5.0, 20.0)), 1)
        self.assertIsNone(diagnostics.canonical_t_delta_index(20.00001, (-10.0, 5.0, 20.0)))

    def test_delta_edges_prefer_the_existing_proton_contract(self):
        edges, source = diagnostics.resolve_pion_hgcer_delta_edges(
            {"reuse_proton_delta_edges": True, "delta_range": (-10.0, 20.0), "delta_bins": 10},
            {"delta_edges": (-10.0, -2.0, 3.0, 20.0)},
        )
        self.assertEqual(edges, [-10.0, -2.0, 3.0, 20.0])
        self.assertEqual(source, "proton_cleaning_result.delta_edges")

    def test_numpy_backed_canonical_and_delta_edges_are_accepted(self):
        t_edges = np.array((0.4, 0.6, 0.9), dtype=float)
        delta_edges = np.array((-10.0, -2.0, 20.0), dtype=float)
        self.assertEqual(diagnostics.canonical_t_delta_index(0.9, t_edges), 1)
        self.assertEqual(diagnostics.canonical_t_delta_index(20.0, delta_edges), 1)
        self.assertEqual(diagnostics._float_edges(t_edges), [0.4, 0.6, 0.9])
        resolved, source = diagnostics.resolve_pion_hgcer_delta_edges(
            {"reuse_proton_delta_edges": True},
            {"delta_edges": delta_edges},
        )
        self.assertEqual(resolved, [-10.0, -2.0, 20.0])
        self.assertEqual(source, "proton_cleaning_result.delta_edges")

    def test_numpy_backed_edges_serialize_without_truth_value_coercion(self):
        serialized = diagnostics.serialize_pion_hgcer_tdelta_diagnostic(
            {
                "status": "unavailable",
                "t_edges": np.array((0.4, 0.6, 0.9), dtype=float),
                "delta_edges": np.array((-10.0, 20.0), dtype=float),
                "cells": np.array((), dtype=object),
            }
        )
        self.assertEqual(serialized["t_edges"], [0.4, 0.6, 0.9])
        self.assertEqual(serialized["delta_edges"], [-10.0, 20.0])
        json.dumps(serialized, allow_nan=False)

    def test_support_requires_non_cancelling_support_on_both_sides(self):
        thresholds = {
            "supported_absolute_weight": 0.1,
            "marginal_absolute_weight": 0.01,
            "supported_effective_entries": 5,
            "marginal_effective_entries": 1,
        }
        supported = {"absolute_weight_support": 0.2, "sum_weight_squared": 0.008}
        marginal = {"absolute_weight_support": 0.02, "sum_weight_squared": 0.0004}
        absent = {"absolute_weight_support": 0.0, "sum_weight_squared": 0.0}
        self.assertEqual(diagnostics.classify_pion_hgcer_support(supported, supported, thresholds), "supported")
        self.assertEqual(diagnostics.classify_pion_hgcer_support(marginal, marginal, thresholds), "marginal")
        self.assertEqual(diagnostics.classify_pion_hgcer_support(supported, absent, thresholds), "unsupported")

    def test_runtime_configuration_is_diagnostic_only_and_validated(self):
        config = get_pion_hgcer_diagnostic_config(
            {"Q2": "4p4", "W": "2p74", "pion_hgcer_diagnostic_config": {"delta_bins": 4}}
        )
        self.assertTrue(config["enabled"])
        self.assertEqual(config["delta_bins"], 4)
        self.assertEqual(config["emit_cell_pages"], "supported_marginal")
        self.assertEqual(config["hgcer_boundary"]["threshold"], 2.0)
        self.assertEqual(config["hgcer_boundary"]["boundary_band_width"], 0.5)
        self.assertEqual(
            config["boundary_readiness_thresholds"]["ready_effective_entries"],
            10.0,
        )
        with self.assertRaisesRegex(ValueError, "must match production_hgcer_threshold"):
            get_pion_hgcer_diagnostic_config({
                "pion_hgcer_diagnostic_config": {
                    "hgcer_boundary": {"threshold": 1.5}
                }
            })

    def test_serializer_never_exposes_root_histogram_payload(self):
        serialized = diagnostics.serialize_pion_hgcer_tdelta_diagnostic(
            {
                "status": "unavailable",
                "reason": "synthetic missing tree",
                "records": {"kaon": (), "pion": ()},
                "histograms": {"must_not": object()},
            }
        )
        self.assertEqual(serialized["reason"], "synthetic missing tree")
        self.assertNotIn("histograms", serialized)

    def test_serializer_excludes_phase_e_acceptance_record_sidecar(self):
        serialized = diagnostics.serialize_pion_hgcer_tdelta_diagnostic(
            {
                "status": "unavailable",
                "coordinate_fingerprint": "sidecar-freeze-coordinate",
                "t_edges": (0.0, 0.4, 1.0),
                "delta_edges": (-10.0, 0.0, 20.0),
                "records": {"kaon": (), "pion": ()},
                "phase_e_acceptance_records": {
                    "kaon": (
                        {
                            "canonical_t_index": 0,
                            "delta_index": 0,
                            "ssxptar": 0.01,
                            "ssyptar": -0.01,
                        },
                    ),
                    "pion": (),
                },
            }
        )
        self.assertNotIn("phase_e_acceptance_records", serialized)
        self.assertEqual(
            serialized["coordinate_fingerprint"], "sidecar-freeze-coordinate"
        )
        self.assertEqual(serialized["t_edges"], [0.0, 0.4, 1.0])
        self.assertEqual(serialized["delta_edges"], [-10.0, 0.0, 20.0])
        json.dumps(serialized, allow_nan=False)

    def test_unavailable_payload_preserves_concise_staged_error_metadata(self):
        original = IndexError("Replacement index 2 out of range for positional args tuple")
        failure = diagnostics._PionHGCerDiagnosticBuildFailure(
            "histogram_construction", original
        )
        self.assertEqual(failure.diagnostic_stage, "histogram_construction")
        self.assertIs(failure.original_exception, original)
        serialized = diagnostics.serialize_pion_hgcer_tdelta_diagnostic(
            {
                "status": "unavailable",
                "reason": str(failure.original_exception),
                "exception_type": type(failure.original_exception).__name__,
                "exception_message": str(failure.original_exception),
                "diagnostic_stage": failure.diagnostic_stage,
            }
        )
        self.assertEqual(serialized["exception_type"], "IndexError")
        self.assertEqual(serialized["exception_message"], str(original))
        self.assertEqual(serialized["diagnostic_stage"], "histogram_construction")

    def test_renderer_manifest_tracks_actual_canonical_t_and_cell_status_slots(self):
        for n_t in (2, 3, 5):
            payload = {
                "status": "available",
                "t_edges": tuple(float(index) for index in range(n_t + 1)),
                "config": {"emit_cell_pages": "supported_marginal", "emit_status_pages": True},
                "cells": tuple(
                    {
                        "t_index": t_index,
                        "delta_index": delta_index,
                        "support_class": ("supported", "marginal", "unsupported")[delta_index],
                    }
                    for t_index in range(n_t)
                    for delta_index in range(3)
                ),
            }
            manifest = diagnostics.expected_pion_hgcer_page_manifest(payload)
            page_ids = [entry["page_id"] for entry in manifest]
            self.assertEqual(len(page_ids), 3 + n_t + (3 * n_t) + 6 + n_t)
            self.assertIn("pion_hgcer.part1.t{}.delta3.status".format(n_t), page_ids)
            self.assertIn("pion_hgcer.part1p5.boundary", page_ids)
            self.assertIn("pion_hgcer.part1p5.summary", page_ids)
            self.assertEqual(len(page_ids), len(set(page_ids)))


@unittest.skipUnless(diagnostics.ROOT is not None, "PyROOT is unavailable")
class PionHGCerTDeltaRootTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._root_was_batch = bool(diagnostics.ROOT.gROOT.IsBatch())
        diagnostics.ROOT.gROOT.SetBatch(True)

    @classmethod
    def tearDownClass(cls):
        diagnostics.ROOT.gROOT.SetBatch(cls._root_was_batch)

    def _config(self):
        return {
            "enabled": True,
            "delta_range": (-10.0, 20.0),
            "delta_bins": 2,
            "reuse_proton_delta_edges": True,
            "hgcer_npe_range": (0.0, 10.0),
            "hgcer_npe_bins": 20,
            "mm_range": (0.7, 1.5),
            "mm_bins": 20,
            "support_thresholds": {
                "supported_absolute_weight": 0.0,
                "marginal_absolute_weight": 0.0,
                "supported_effective_entries": 0.0,
                "marginal_effective_entries": 0.0,
            },
            "emit_cell_pages": "none",
            "production_hgcer_threshold": 2.0,
            "hgcer_boundary": {
                "threshold": 2.0,
                "zoom_min": 0.0,
                "zoom_max": 4.0,
                "boundary_band_width": 0.5,
            },
            "boundary_readiness_thresholds": {
                "ready_absolute_weight": 0.0,
                "marginal_absolute_weight": 0.0,
                "ready_effective_entries": 0.0,
                "marginal_effective_entries": 0.0,
            },
            "emit_boundary_cell_pages": "none",
        }

    def test_histogram_construction_is_detached_sumw2_and_fully_labeled(self):
        t_edges = [0.4, 0.5666667, 0.7333333, 0.9]
        delta_edges = [-10.0, -7.0, -4.0, -1.0, 2.0, 5.0, 8.0, 11.0, 14.0, 17.0, 20.0]
        histograms = diagnostics._make_histograms(
            t_edges, delta_edges, self._config()
        )
        expected_keys = {"H_support_class", "H_boundary_readiness"}
        for side in ("kaon", "pion"):
            expected_keys.update(
                {
                    "H_support_absolute_{}".format(side),
                    "H_support_effective_{}".format(side),
                    "H_boundary_effective_{}".format(side),
                }
            )
            for weighting in ("weighted", "absolute"):
                for quantity in (
                    "H_hgcer",
                    "H_delta",
                    "H_t",
                    "H_mm",
                    "H_hgcer_vs_delta",
                    "H_hgcer_vs_t",
                    "H_mm_vs_delta",
                    "H_mm_vs_hgcer",
                ):
                    expected_keys.add("{}_{}_{}".format(quantity, side, weighting))
        self.assertEqual(set(histograms), expected_keys)
        for histogram in histograms.values():
            self.assertFalse(bool(histogram.GetDirectory()))
            self.assertGreater(histogram.GetSumw2N(), 0)
            self.assertTrue(histogram.GetTitle())
            self.assertTrue(histogram.GetXaxis().GetTitle())
        title = histograms["H_hgcer_kaon_weighted"].GetTitle()
        self.assertIn("P_hgcer_npeSum", title)
        self.assertIn("signed weighted yield", title)

    def _build_renderable_synthetic_result(self):
        contract = build_kaon_data_coordinate_contract(
            "Center", {"shift": 0.0}, {"shift": 0.0}, require_t_shift=True
        )
        events = [
            _Event(mm=1.08, mandel_t=-0.4, delta=-10.0, npe=0.5),
            _Event(mm=1.12, mandel_t=-0.9, delta=20.0, npe=3.0),
        ]
        proton_cleaning_result = {
            "accepted": True,
            "coordinate_fingerprint": contract["coordinate_fingerprint"],
            "delta_edges": (-10.0, 5.0, 20.0),
            "_prepared_event_weight_lookup": {
                "prompt:0": {"cleaned_factor": 0.5},
                "prompt:1": {"cleaned_factor": 0.5},
            },
        }
        pion_control_cache = {
            "coordinate_fingerprint": contract["coordinate_fingerprint"]
        }
        proton_cleaning_before = deepcopy(proton_cleaning_result)
        pion_control_before = deepcopy(pion_control_cache)
        kaon_bundle = {
            "sources": {
                "prompt": {
                    "tree": events,
                    "tree_name": "Cut_Kaon_Events_prompt_noRF",
                    "coefficient": 1.0,
                }
            },
            "prepared_sources": {
                "prompt": {
                    "entries": {
                        0: {"allcuts": True, "nommcuts": True},
                        1: {"allcuts": True, "nommcuts": True},
                    }
                }
            },
        }
        pion_bundle = {
            "prompt_tree": events,
            "rand_tree": (),
            "dummy_prompt_tree": (),
            "dummy_rand_tree": (),
            "prompt_tree_name": "Cut_Pion_Events_prompt_noRF",
            "rand_tree_name": "Cut_Pion_Events_rand_noRF",
        }
        config = self._config()
        config.update(
            {
                "emit_cell_pages": "supported_marginal",
                "emit_status_pages": True,
                "support_thresholds": {
                    "supported_absolute_weight": 0.75,
                    "marginal_absolute_weight": 0.1,
                    "supported_effective_entries": 1.0,
                    "marginal_effective_entries": 1.0,
                },
            }
        )
        return (
            diagnostics.build_pion_hgcer_tdelta_diagnostic(
                kaon_source_bundle=kaon_bundle,
                pion_tree_bundle=pion_bundle,
                proton_cleaning_result=proton_cleaning_result,
                pion_control_cache=pion_control_cache,
                coordinate_contract=contract,
                t_edges=(0.4, 0.6, 0.9),
                config=config,
                hole_contains=lambda *_args: False,
                evaluate_pion_event=lambda *_args, **_kwargs: (True, True, 0.0),
                mm_min=1.1,
                mm_max=1.16,
            ),
            proton_cleaning_result,
            pion_control_cache,
            proton_cleaning_before,
            pion_control_before,
        )

    def test_synthetic_build_serialize_and_render_cover_cells_and_final_edges(self):
        (
            result,
            proton_cleaning_result,
            pion_control_cache,
            proton_cleaning_before,
            pion_control_before,
        ) = (
            self._build_renderable_synthetic_result()
        )
        self.assertEqual(result["status"], "available")
        self.assertEqual(proton_cleaning_result, proton_cleaning_before)
        self.assertEqual(pion_control_cache, pion_control_before)
        self.assertEqual(
            [record["canonical_t_index"] for record in result["records"]["kaon"]],
            [0, 1],
        )
        self.assertEqual(
            [record["delta_index"] for record in result["records"]["pion"]],
            [0, 1],
        )
        self.assertEqual(
            result["source_audit"]["kaon"]["prompt"]["records_below_reference_boundary"],
            1,
        )
        self.assertEqual(
            result["source_audit"]["pion"]["prompt"]["records_at_or_above_reference_boundary"],
            1,
        )
        self.assertEqual(
            [record["proton_cleaning_factor"] for record in result["records"]["kaon"]],
            [0.5, 0.5],
        )
        self.assertTrue(
            all(
                record["proton_cleaning_factor"] is None
                for record in result["records"]["pion"]
            )
        )
        self.assertEqual(
            {cell["support_class"] for cell in result["cells"]},
            {"marginal", "unsupported"},
        )
        serialized = diagnostics.serialize_pion_hgcer_tdelta_diagnostic(result)
        json.dumps(serialized, allow_nan=False)
        self.assertNotIn("histograms", serialized)
        expected_manifest = diagnostics.expected_pion_hgcer_page_manifest(result)
        with tempfile.TemporaryDirectory() as temporary_directory:
            pdf_path = Path(temporary_directory) / "pion_hgcer_part1.pdf"
            opener = diagnostics.ROOT.TCanvas(
                diagnostics.unique_root_object_name(
                    "C_pion_hgcer_test_open", scope="pion_hgcer", role="test"
                ),
                "Part-1 test PDF opener",
                100,
                100,
            )
            opener.Print(str(pdf_path) + "(")
            emitted = diagnostics.render_pion_hgcer_tdelta_pages(
                str(pdf_path), result, close_pdf=True
            )
            self.assertEqual(emitted, expected_manifest)
            self.assertTrue(
                any(entry["page_id"].endswith(".status") for entry in emitted)
            )
            self.assertTrue(pdf_path.is_file())
            self.assertGreater(pdf_path.stat().st_size, 0)
        full_range_axis = result["histograms"]["H_hgcer_kaon_weighted"].GetXaxis()
        self.assertEqual(full_range_axis.GetFirst(), 1)
        self.assertEqual(
            full_range_axis.GetLast(),
            result["histograms"]["H_hgcer_kaon_weighted"].GetNbinsX(),
        )
        histogram = result["histograms"]["H_hgcer_kaon_absolute"]
        gc.collect()
        self.assertGreater(histogram.Integral(), 0.0)

    def test_pre_hgcer_populations_share_coordinates_and_preserve_factor_scope(self):
        contract = build_kaon_data_coordinate_contract(
            "Center", {"shift": 0.01}, {"shift": -0.02}, require_t_shift=True
        )
        kaon_events = [_Event(mm=1.08, mandel_t=-0.42, delta=-5.0, npe=0.5)]
        pion_events = [_Event(mm=1.08, mandel_t=-0.42, delta=-5.0, npe=0.75)]
        factor_lookup = {
            "prompt:0": {
                "cleaned_factor": 0.4,
                # A diagnostic must not use this RF-restored value.
                "final_cleaned_factor": 0.0,
            }
        }
        kaon_bundle = {
            "sources": {"prompt": {"tree": kaon_events, "tree_name": "Cut_Kaon_Events_prompt_noRF", "coefficient": 2.0}},
            "prepared_sources": {
                "prompt": {
                    "entries": {0: {"allcuts": True, "nommcuts": True}},
                    "coefficient": 2.0,
                }
            },
        }
        pion_bundle = {
            "prompt_tree": pion_events,
            "rand_tree": (),
            "dummy_prompt_tree": (),
            "dummy_rand_tree": (),
            "prompt_tree_name": "Cut_Pion_Events_prompt_noRF",
            "rand_tree_name": "Cut_Pion_Events_rand_noRF",
        }
        result = diagnostics.build_pion_hgcer_tdelta_diagnostic(
            kaon_source_bundle=kaon_bundle,
            pion_tree_bundle=pion_bundle,
            proton_cleaning_result={
                "accepted": True,
                "coordinate_fingerprint": contract["coordinate_fingerprint"],
                "delta_edges": (-10.0, 5.0, 20.0),
                "_prepared_event_weight_lookup": factor_lookup,
            },
            pion_control_cache={"coordinate_fingerprint": contract["coordinate_fingerprint"]},
            coordinate_contract=contract,
            t_edges=(0.4, 0.9),
            config=self._config(),
            hole_contains=lambda *_args: False,
            evaluate_pion_event=lambda *_args, **_kwargs: (True, True, 0.0),
            mm_min=1.1,
            mm_max=1.16,
        )
        self.assertEqual(result["status"], "available")
        self.assertEqual(len(result["records"]["kaon"]), 1)
        self.assertEqual(len(result["records"]["pion"]), 1)
        kaon_record = result["records"]["kaon"][0]
        pion_record = result["records"]["pion"][0]
        self.assertAlmostEqual(kaon_record["diagnostic_weight"], 0.8)
        self.assertIsNone(pion_record["proton_cleaning_factor"])
        self.assertAlmostEqual(pion_record["diagnostic_weight"], 2.0)
        self.assertAlmostEqual(kaon_record["analysis_MM"], 1.09)
        self.assertAlmostEqual(pion_record["analysis_t"], 0.40)
        self.assertEqual(result["source_audit"]["pion"]["prompt"]["records_below_reference_boundary"], 1)
        histogram = result["histograms"]["H_hgcer_kaon_absolute"]
        gc.collect()
        self.assertGreater(histogram.Integral(), 0.0)

    def test_boundary_metrics_keep_pid_sides_separate_on_the_common_cell_grid(self):
        contract = build_kaon_data_coordinate_contract(
            "Center", {"shift": 0.0}, {"shift": 0.0}, require_t_shift=True
        )
        kaon_event = _Event(mm=1.10, mandel_t=-0.50, delta=0.0, npe=1.50)
        pion_event = _Event(mm=1.10, mandel_t=-0.50, delta=0.0, npe=2.00)
        result = diagnostics.build_pion_hgcer_tdelta_diagnostic(
            kaon_source_bundle={
                "sources": {
                    "prompt": {
                        "tree": [kaon_event],
                        "tree_name": "Cut_Kaon_Events_prompt_noRF",
                        "coefficient": 1.0,
                    }
                },
                "prepared_sources": {
                    "prompt": {"entries": {0: {"allcuts": True, "nommcuts": True}}}
                },
            },
            pion_tree_bundle={
                "prompt_tree": [pion_event],
                "rand_tree": (),
                "dummy_prompt_tree": (),
                "dummy_rand_tree": (),
                "prompt_tree_name": "Cut_Pion_Events_prompt_noRF",
                "rand_tree_name": "Cut_Pion_Events_rand_noRF",
            },
            proton_cleaning_result={
                "accepted": True,
                "coordinate_fingerprint": contract["coordinate_fingerprint"],
                "delta_edges": (-10.0, 5.0, 20.0),
                "_prepared_event_weight_lookup": {"prompt:0": {"cleaned_factor": 1.0}},
            },
            pion_control_cache={"coordinate_fingerprint": contract["coordinate_fingerprint"]},
            coordinate_contract=contract,
            t_edges=(0.4, 0.9),
            config=self._config(),
            hole_contains=lambda *_args: False,
            evaluate_pion_event=lambda *_args, **_kwargs: (True, True, 0.0),
            mm_min=1.1,
            mm_max=1.16,
        )
        kaon_record = result["records"]["kaon"][0]
        pion_record = result["records"]["pion"][0]
        self.assertEqual(
            (kaon_record["canonical_t_index"], kaon_record["delta_index"]),
            (pion_record["canonical_t_index"], pion_record["delta_index"]),
        )
        cell = result["cells"][0]
        self.assertEqual(cell["boundary_support"]["kaon"]["record_count"], 1)
        self.assertEqual(cell["boundary_support"]["pion"]["record_count"], 1)
        self.assertEqual(cell["boundary_readiness"], "ready")
        histograms = result["histograms"]
        self.assertGreater(histograms["H_t_kaon_absolute"].Integral(), 0.0)
        self.assertGreater(histograms["H_t_pion_absolute"].Integral(), 0.0)
        self.assertEqual(
            histograms["H_hgcer_pion_absolute"].GetBinContent(
                histograms["H_hgcer_pion_absolute"].FindBin(1.50)
            ),
            0.0,
        )
        self.assertEqual(
            histograms["H_hgcer_kaon_absolute"].GetBinContent(
                histograms["H_hgcer_kaon_absolute"].FindBin(2.00)
            ),
            0.0,
        )

    def test_no_below_threshold_events_are_reported_as_a_source_limitation(self):
        contract = build_kaon_data_coordinate_contract(
            "Center", {"shift": 0.0}, {"shift": 0.0}, require_t_shift=True
        )
        event = _Event(mm=1.10, mandel_t=-0.50, delta=0.0, npe=3.0)
        kaon_bundle = {
            "sources": {"prompt": {"tree": [event], "tree_name": "Cut_Kaon_Events_prompt_noRF", "coefficient": 1.0}},
            "prepared_sources": {"prompt": {"entries": {0: {"allcuts": True, "nommcuts": True}}}},
        }
        pion_bundle = {
            "prompt_tree": [event],
            "rand_tree": (),
            "dummy_prompt_tree": (),
            "dummy_rand_tree": (),
            "prompt_tree_name": "Cut_Pion_Events_prompt_noRF",
            "rand_tree_name": "Cut_Pion_Events_rand_noRF",
        }
        result = diagnostics.build_pion_hgcer_tdelta_diagnostic(
            kaon_source_bundle=kaon_bundle,
            pion_tree_bundle=pion_bundle,
            proton_cleaning_result={
                "accepted": True,
                "coordinate_fingerprint": contract["coordinate_fingerprint"],
                "_prepared_event_weight_lookup": {"prompt:0": {"cleaned_factor": 1.0}},
            },
            pion_control_cache={"coordinate_fingerprint": contract["coordinate_fingerprint"]},
            coordinate_contract=contract,
            t_edges=(0.4, 0.9),
            config=self._config(),
            hole_contains=lambda *_args: False,
            evaluate_pion_event=lambda *_args, **_kwargs: (True, True, 0.0),
            mm_min=1.1,
            mm_max=1.16,
        )
        self.assertIn(
            "tree_may_be_truncated",
            result["source_audit"]["pion"]["prompt"]["source_tree_boundary_coverage"],
        )


class PionHGCerTDeltaIntegrationContracts(unittest.TestCase):
    def test_rand_sub_calls_the_diagnostic_after_authoritative_cache_and_keeps_it_out_of_parents(self):
        source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        call = source.index("build_pion_hgcer_tdelta_diagnostic(")
        cache = source.index("_build_authoritative_pion_control_source_cache(")
        parents = source.index("build_setting_t_bin_pion_parents(")
        self.assertGreater(call, cache)
        self.assertLess(call, parents)
        diagnostic_region = source[call:parents]
        self.assertIn("proton_cleaning_result=proton_cleaning_result", diagnostic_region)
        self.assertNotIn("_pion_t_bin_parent_results", diagnostic_region)
        self.assertIn("cleaned_factor", (REPO_ROOT / "src" / "cuts" / "pion_hgcer_diagnostics.py").read_text(encoding="utf-8"))
        self.assertNotIn("final_cleaned_factor", (REPO_ROOT / "src" / "cuts" / "pion_hgcer_diagnostics.py").read_text(encoding="utf-8"))

    def test_particle_background_tree_resolution_is_explicitly_no_rf(self):
        source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        self.assertIn('subtracted_particle, epsset, rf_state="noRF"', source)
        self.assertIn('ParticleType, EPSSET, rf_state="noRF"', source)


if __name__ == "__main__":
    unittest.main()
