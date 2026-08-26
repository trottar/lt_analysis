"""Contracts for the non-authoritative pre-HGCer Part-1 diagnostic."""

from __future__ import annotations

import gc
import importlib.util
import sys
import unittest
from pathlib import Path


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
            self.assertEqual(len(page_ids), 3 + n_t + (3 * n_t))
            self.assertIn("pion_hgcer.part1.t{}.delta3.status".format(n_t), page_ids)
            self.assertEqual(len(page_ids), len(set(page_ids)))


@unittest.skipUnless(diagnostics.ROOT is not None, "PyROOT is unavailable")
class PionHGCerTDeltaRootTests(unittest.TestCase):
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
        }

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
            "prompt_tree_name": "Cut_Pion_Events_prompt_RF",
            "rand_tree_name": "Cut_Pion_Events_rand_RF",
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
        self.assertEqual(result["source_audit"]["pion"]["prompt"]["hgc_below_threshold_records"], 1)
        histogram = result["histograms"]["H_hgcer_kaon_absolute"]
        gc.collect()
        self.assertGreater(histogram.Integral(), 0.0)

    def test_no_below_threshold_events_are_reported_as_a_source_limitation(self):
        contract = build_kaon_data_coordinate_contract(
            "Center", {"shift": 0.0}, {"shift": 0.0}, require_t_shift=True
        )
        event = _Event(mm=1.10, mandel_t=-0.50, delta=0.0, npe=3.0)
        kaon_bundle = {
            "sources": {"prompt": {"tree": [event], "coefficient": 1.0}},
            "prepared_sources": {"prompt": {"entries": {0: {"allcuts": True, "nommcuts": True}}}},
        }
        pion_bundle = {"prompt_tree": [event], "rand_tree": (), "dummy_prompt_tree": (), "dummy_rand_tree": ()}
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
            result["source_audit"]["pion"]["prompt"]["tree_hgcer_truncation_evidence"],
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


if __name__ == "__main__":
    unittest.main()
