"""Pure-Python contracts for C.4 detached HGCer PDF presentation."""

from __future__ import annotations

import ast
from pathlib import Path
import sys
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_refinement_plots as plots


def _checkpoint():
    return {
        "setting": {
            "phi_setting": "Left",
            "kinematic_token": "Q4p4W2p74",
            "epsilon_filename_token": "highe",
        },
        "canonical_t_edges": [0.1, 0.2, 0.3],
        "delta_edges": [-10.0, 0.0, 10.0],
        "phase_a": {
            "contract_fingerprint": "0123456789abcdef0123456789abcdef",
            "coordinate_fingerprint": "abcdef0123456789abcdef0123456789",
            "host_state": "proton_cleaned",
            "source_target_state": "post_proton_noRF",
        },
        "method_a": {"summary": {"support_counts": {"supported": 1, "marginal": 1, "unsupported": 2}}},
        "method_b": {
            "summary": {
                "method_B_status_counts": {"available": 1, "marginal": 1, "unavailable": 2},
                "shape_status_counts": {"good": 1, "poor": 1},
            }
        },
    }


def _method_a():
    return {
        "status": "available",
        "available": True,
        "t_edges": [0.1, 0.2, 0.3],
        "delta_edges": [-10.0, 0.0, 10.0],
        "summary": {"support_counts": {"supported": 1, "marginal": 1, "unsupported": 2}},
        "cells": [
            {
                "t_index": 0, "delta_index": 0, "delta_low": -10.0, "delta_high": 0.0,
                "support_class": "supported", "f_low": 0.25, "f_low_low": 0.20, "f_low_high": 0.30,
            },
            {
                "t_index": 0, "delta_index": 1, "delta_low": 0.0, "delta_high": 10.0,
                "support_class": "marginal", "f_low": 0.50, "f_low_low": 0.40, "f_low_high": 0.60,
            },
            {
                "t_index": 1, "delta_index": 0, "delta_low": -10.0, "delta_high": 0.0,
                "support_class": "unsupported", "f_low": None, "f_low_low": None, "f_low_high": None,
            },
        ],
    }


def _method_b():
    return {
        "status": "available",
        "available": True,
        "t_edges": [0.1, 0.2, 0.3],
        "delta_edges": [-10.0, 0.0, 10.0],
        "summary": {"method_B_status_counts": {"available": 1, "marginal": 1, "unavailable": 2}},
        "cells": [
            {
                "t_index": 0, "delta_index": 0, "delta_low": -10.0, "delta_high": 0.0,
                "method_B_status": "available", "candidate_L_B": 1.05,
                "candidate_L_B_uncertainty": 0.10, "candidate_L_B_status": "available_multi_region",
                "shape_chi2_ndf": 1.2, "shape_max_abs_pull": 1.5, "shape_status": "good",
                "regions": [
                    {"region_name": "pi_n", "parent_relative_ratio": 1.02, "parent_relative_sigma": 0.1, "parent_relative_status": "available"},
                    {"region_name": "pi_sidis", "parent_relative_ratio": 1.08, "parent_relative_sigma": 0.1, "parent_relative_status": "available"},
                    {"region_name": "pi_delta_high", "parent_relative_ratio": 1.03, "parent_relative_sigma": 0.1, "parent_relative_status": "available"},
                    {"region_name": "pi_delta_low", "parent_relative_ratio": 99.0, "parent_relative_sigma": 1.0, "parent_relative_status": "available"},
                ],
            },
            {
                "t_index": 0, "delta_index": 1, "delta_low": 0.0, "delta_high": 10.0,
                "method_B_status": "shape_inconsistent", "candidate_L_B": None,
                "candidate_L_B_uncertainty": None, "candidate_L_B_status": "shape_poor_veto",
                "shape_chi2_ndf": None, "shape_max_abs_pull": None, "shape_status": "poor",
                "regions": [],
            },
        ],
    }


class PionHGCerRefinementPlotTests(unittest.TestCase):
    def test_supplement_names_and_routes_are_deterministic(self):
        main = r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe.pdf"
        destinations = plots.build_pdf_destinations(main)
        self.assertEqual(destinations["main"], main)
        self.assertEqual(destinations["proton_debug"], r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe_proton-debug.pdf")
        self.assertEqual(destinations["pion_fit_debug"], r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe_pion-fit-debug.pdf")
        self.assertEqual(destinations["hgcer_debug"], r"C:\analysis\Left_kaon_rand_sub_Q4p4W2p74_highe_hgcer-debug.pdf")
        manifest = plots.build_pdf_route_manifest(main)
        self.assertIn("hgcer.phase_a.summary", manifest["routes"]["main"])
        self.assertEqual(
            [page_id for page_id in manifest["routes"]["main"] if page_id.startswith("proton.summary")],
            [
                "proton.summary.provenance_closure",
                "proton.summary.committed_mm",
                "proton.summary.commitment",
            ],
        )
        self.assertIn("pion.coordinate.detail", manifest["routes"]["pion_fit_debug"])
        self.assertIn("hgcer.part2", manifest["routes"]["hgcer_debug"])
        self.assertIn("proton.detail", manifest["routes"]["proton_debug"])

    def test_phase_a_summary_uses_stored_contract_and_checkpoint_fields(self):
        phase = {
            "status": "available", "available": True,
            "pion_closure": {"passed": True}, "host_closure": {"passed": True},
            "lambda_gate_status": "PASS", "lambda_gate_production_action": "committed",
        }
        payload = plots.phase_a_summary_payload(_checkpoint(), phase)
        self.assertEqual(payload["host_state"], "proton_cleaned")
        self.assertEqual(payload["source_target_state"], "post_proton_noRF")
        self.assertEqual(payload["canonical_t_edges"], [0.1, 0.2, 0.3])
        self.assertEqual(payload["delta_edges"], [-10.0, 0.0, 10.0])
        self.assertEqual(payload["lambda_gate_status"], "PASS")
        self.assertTrue(payload["non_authoritative"])
        self.assertFalse(payload["production_objects_mutated"])
        self.assertFalse(payload["refinement_applied"])

    def test_phase_a_runtime_display_context_is_explicit_and_read_only(self):
        payload = plots.phase_a_summary_payload(
            _checkpoint(),
            {"status": "available", "available": True},
            {
                "lambda_gate_status": "pass",
                "production_action": "apply",
                "proton_cleaning_committed": True,
                "host_state": "post_proton_noRF",
            },
        )
        self.assertEqual(payload["lambda_gate_status"], "pass")
        self.assertEqual(payload["lambda_gate_production_action"], "apply")
        self.assertTrue(payload["proton_cleaning_committed"])
        self.assertEqual(payload["host_state"], "post_proton_noRF")

    def test_method_a_keeps_fidelity_fields_nulls_and_support_counts(self):
        payload = plots.method_a_plot_payload(_method_a(), _checkpoint())
        self.assertEqual(payload["support_counts"], {"supported": 1, "marginal": 1, "unsupported": 1})
        self.assertEqual(payload["cells"][2]["f_low"], None)
        points = plots.method_a_f_low_points(payload)
        self.assertEqual(len(points), 2)
        self.assertEqual([point["support_class"] for point in points], ["supported", "marginal"])
        self.assertAlmostEqual(points[0]["f_low_low"], 0.20)
        self.assertTrue(payload["no_interpolation"])
        self.assertTrue(payload["no_correction_applied"])

    def test_method_b_keeps_status_candidate_regions_and_shape_without_substitution(self):
        payload = plots.method_b_plot_payload(_method_b(), _checkpoint())
        self.assertEqual(payload["method_status_counts"], {"available": 1, "marginal": 0, "shape_inconsistent": 1, "unavailable": 0})
        candidates = plots.method_b_candidate_points(payload)
        self.assertEqual(candidates, [{
            "t_index": 0, "delta_index": 0, "delta_center": -5.0,
            "candidate_L_B": 1.05, "candidate_L_B_uncertainty": 0.10,
            "candidate_L_B_status": "available_multi_region",
        }])
        regions = plots.method_b_regional_rows(payload)
        self.assertEqual({row["region_name"] for row in regions}, {"pi_n", "pi_sidis", "pi_delta_high"})
        self.assertNotIn("pi_delta_low", {row["region_name"] for row in regions})
        shapes = plots.method_b_shape_rows(payload)
        self.assertEqual(shapes[1]["shape_chi2_ndf"], None)
        self.assertEqual(shapes[1]["shape_status"], "poor")
        self.assertEqual(shapes[1]["candidate_L_B_status"], "shape_poor_veto")
        self.assertTrue(shapes[1]["shape_poor_veto"])
        self.assertTrue(payload["frozen_pion_baseline"])
        self.assertTrue(payload["no_refinement"])

    def test_method_b_regional_panels_are_parent_local_and_unity_uses_full_delta_edges(self):
        payload = plots.method_b_plot_payload(_method_b(), _checkpoint())
        panels = plots.method_b_regional_panels(payload)
        self.assertEqual([panel["t_index"] for panel in panels], [0])
        self.assertEqual(set(panels[0]["series"]), {"pi_n", "pi_sidis", "pi_delta_high"})
        self.assertTrue(all(len(points) == 1 for points in panels[0]["series"].values()))
        self.assertEqual(plots.unity_line_limits(payload), (-10.0, 10.0))

    def test_method_b_regional_panels_do_not_merge_different_t_parents(self):
        method_b = {
            "t_edges": [0.1, 0.2, 0.3],
            "delta_edges": [-10.0, 10.0],
            "cells": [
                {
                    "t_index": 0, "delta_index": 0, "delta_low": -10.0, "delta_high": 10.0,
                    "method_B_status": "available",
                    "regions": [{"region_name": "pi_n", "parent_relative_ratio": 1.0, "parent_relative_sigma": 0.1, "parent_relative_status": "available"}],
                },
                {
                    "t_index": 1, "delta_index": 0, "delta_low": -10.0, "delta_high": 10.0,
                    "method_B_status": "available",
                    "regions": [{"region_name": "pi_n", "parent_relative_ratio": 1.2, "parent_relative_sigma": 0.1, "parent_relative_status": "available"}],
                },
            ],
        }
        panels = plots.method_b_regional_panels(plots.method_b_plot_payload(method_b, _checkpoint()))
        self.assertEqual([panel["t_index"] for panel in panels], [0, 1])
        self.assertEqual(panels[0]["series"]["pi_n"][0]["Qtilde"], 1.0)
        self.assertEqual(panels[1]["series"]["pi_n"][0]["Qtilde"], 1.2)

    def test_single_region_candidate_does_not_create_method_b_availability(self):
        method_b = {
            "status": "unavailable",
            "cells": [{
                "t_index": 0, "delta_index": 0, "delta_low": -10.0, "delta_high": 10.0,
                "method_B_status": "unavailable", "candidate_L_B_status": "single_region_only",
                "candidate_L_B": 1.0, "candidate_L_B_uncertainty": 0.1,
            }],
        }
        payload = plots.method_b_plot_payload(method_b, _checkpoint())
        self.assertEqual(payload["method_status_counts"], {"available": 0, "marginal": 0, "shape_inconsistent": 0, "unavailable": 1})
        self.assertEqual(plots.method_b_candidate_points(payload), [])

    def test_shape_poor_veto_uses_the_stored_candidate_status_not_shape_thresholds(self):
        payload = plots.method_b_plot_payload(
            {
                "cells": [{
                    "t_index": 0, "delta_index": 0,
                    "method_B_status": "shape_inconsistent",
                    "candidate_L_B_status": "shape_poor_veto",
                    "shape_chi2_ndf": 0.2,
                    "shape_max_abs_pull": 0.1,
                    "shape_status": "good",
                }],
            },
            _checkpoint(),
        )
        row = plots.method_b_shape_rows(payload)[0]
        self.assertEqual(row["candidate_L_B_status"], "shape_poor_veto")
        self.assertTrue(row["shape_poor_veto"])
        self.assertEqual(row["shape_chi2_ndf"], 0.2)
        self.assertEqual(row["shape_status"], "good")

    def test_setting_qa_summary_reports_coverage_and_unavailable_optional_inputs(self):
        summary = plots.setting_qa_summary_payload(
            {"status": "available", "host_state": "post_proton_noRF", "pion_closure": {"passed": True}, "host_closure": {"passed": True}},
            _method_a(),
            _method_b(),
            {"status": "unavailable"},
            display_context={"lambda_gate_status": "pass", "production_action": "apply", "proton_cleaning_committed": True},
        )
        self.assertEqual(summary["method_a_coverage"], {"supported": 1, "marginal": 1, "unsupported": 1})
        self.assertEqual(summary["method_b_coverage"]["shape_inconsistent"], 1)
        self.assertEqual(summary["phase_host_state"], "post_proton_noRF")
        self.assertEqual(summary["lambda_gate_status"], "pass")
        self.assertEqual(summary["aerogel_warnings"], "not available")
        self.assertEqual(summary["k_sigma0_availability"], "not available")

    def test_proton_main_qa_keeps_existing_closure_and_commitment_records(self):
        qa = plots.proton_main_qa_payload(
            {
                "status": "available",
                "method": "timing_t_event_weight",
                "diagnostics": {
                    "canonical_t_binning": {"status": "available"},
                    "cross_stage_t_consistency_summary": {"passed": True},
                    "lambda_preservation_gate": {"status": "pass", "production_action": "apply"},
                },
            },
            {"accepted": True, "host_state": "post_proton_noRF", "diagnostics": {"canonical_t_global_closure": {"passed": True}}},
            {"proton_cleaning_committed": True},
        )
        self.assertEqual(qa["canonical_t_binning"], {"status": "available"})
        self.assertEqual(qa["shifted_t_consistency"], {"passed": True})
        self.assertEqual(qa["global_closure"], {"passed": True})
        self.assertTrue(qa["proton_cleaning_committed"])

    def test_category_labels_f_low_styles_and_annotation_are_explicit(self):
        self.assertEqual(set(plots._METHOD_A_SUPPORT_LABELS), {"supported", "marginal", "unsupported"})
        self.assertEqual(set(plots._METHOD_B_STATUS_LABELS), {"available", "marginal", "shape_inconsistent", "unavailable"})
        self.assertNotEqual(plots.method_a_f_low_style("supported"), plots.method_a_f_low_style("marginal"))
        self.assertIsNone(plots.method_a_f_low_style("unsupported"))
        self.assertEqual(
            plots.refinement_annotation_lines(candidate=True),
            (
                "NON-AUTHORITATIVE DIAGNOSTIC / No refinement applied",
                "Diagnostic candidate only; no event correction",
            ),
        )

    def test_warning_page_preserves_only_non_ok_state_and_production_impact(self):
        warnings = plots.warning_payload(
            {"status": "available"},
            {"status": "unavailable", "reason": "support_insufficient"},
            {"status": "marginal", "reason": "regional_overlap"},
            {"status": "unavailable", "reason": "current_source_missing"},
        )
        self.assertEqual([entry["scope"] for entry in warnings], ["Method A", "Method B", "Part 2"])
        self.assertTrue(all(entry["production_impact"] == "none; detached diagnostic only" for entry in warnings))

    def test_synthetic_renderer_smoke_is_a_noop_without_root(self):
        manifest = []
        with mock.patch.object(plots, "_import_root", return_value=None):
            result = plots.render_pion_hgcer_refinement_pages(
                "synthetic.pdf", _checkpoint(), phase_a={}, method_a=_method_a(), method_b=_method_b(), page_manifest=manifest
            )
        self.assertIs(result, manifest)
        self.assertEqual(manifest, [])

    def test_renderer_has_no_analysis_or_physics_imports_or_calls(self):
        source_path = REPO_ROOT / "src" / "cuts" / "pion_hgcer_refinement_plots.py"
        source = source_path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        self.assertTrue(set(imports).issubset({"__future__", "ROOT", "array", "math", "os", "textwrap"}))
        forbidden = {"method_comparison", "method_agreement", "C_B", "C_final", "refined_pion_weight", "applied_refinement_weight"}
        self.assertTrue(forbidden.isdisjoint(source))
        for marker in ("build_pion", "apply_pion", "particle_subtraction", "pion_component_fits", "pion_t_bin_parents"):
            self.assertNotIn(marker, source)
        self.assertIn("NON-AUTHORITATIVE DIAGNOSTIC / No refinement applied", source)
        self.assertIn("Diagnostic candidate only; no event correction", source)
        self.assertIn("unity_line_limits(payload)", source)
        self.assertNotIn("ROOT.TLine(0.0, 1.0", source)


if __name__ == "__main__":
    unittest.main()
