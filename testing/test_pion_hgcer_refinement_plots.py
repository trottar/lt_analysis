"""Pure-Python contracts for C.4 detached HGCer PDF presentation."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import sys
import types
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


def _left_low_method_b():
    """Reduced deterministic copy of the persisted Left-low Method-B fields."""
    delta_edges = [-10.0, -7.0, -4.0, -1.0, 2.0, 5.0, 8.0, 11.0, 14.0, 17.0, 20.0]
    candidates = {
        (0, 4): (0.9547608416143295, 0.12140029583276699),
        (1, 2): (1.0713191848003094, 0.20951468382894808),
    }
    remaining_statuses = (
        ["region_marginal"] * 3
        + ["single_region_only"] * 7
        + ["unavailable"] * 18
    )
    cells = []
    remaining_index = 0
    for t_index in range(3):
        for delta_index in range(10):
            low, high = delta_edges[delta_index:delta_index + 2]
            candidate = candidates.get((t_index, delta_index))
            if candidate is not None:
                candidate_status = "available_multi_region"
                method_status = "available"
                candidate_value, candidate_uncertainty = candidate
                shape_status = "good"
            else:
                candidate_status = remaining_statuses[remaining_index]
                remaining_index += 1
                method_status = "marginal" if candidate_status == "region_marginal" else "unavailable"
                candidate_value = None
                candidate_uncertainty = None
                shape_status = "marginal" if method_status == "marginal" else "unavailable"
            cells.append({
                "t_index": t_index,
                "delta_index": delta_index,
                "delta_low": low,
                "delta_high": high,
                "method_B_status": method_status,
                "candidate_L_B": candidate_value,
                "candidate_L_B_uncertainty": candidate_uncertainty,
                "candidate_L_B_status": candidate_status,
                "shape_chi2_ndf": None,
                "shape_max_abs_pull": None,
                "shape_status": shape_status,
                "regions": [],
            })
    return {
        "status": "available",
        "available": True,
        "t_edges": [0.1, 0.2, 0.3, 0.4],
        "delta_edges": delta_edges,
        "cells": cells,
        "summary": {
            "candidate_status_counts": {
                "available_multi_region": 2,
                "region_marginal": 3,
                "single_region_only": 7,
                "unavailable": 18,
            },
            "method_B_status_counts": {
                "available": 2,
                "marginal": 3,
                "shape_inconsistent": 0,
                "unavailable": 25,
            },
        },
    }


def _left_low_checkpoint():
    return {
        "setting": {"phi_setting": "Left", "kinematic_token": "Q4p4W2p74", "epsilon_filename_token": "low"},
        "canonical_t_edges": [0.1, 0.2, 0.3, 0.4],
        "delta_edges": [-10.0, -7.0, -4.0, -1.0, 2.0, 5.0, 8.0, 11.0, 14.0, 17.0, 20.0],
        "method_b": _left_low_method_b(),
    }


def _load_parent_lambda_renderer(comparison_renderer):
    """Load the parent presentation wrapper with inert non-rendering imports."""
    background_config = types.ModuleType("background_config")
    background_config.get_particle_subtraction_setting_key = lambda *_args, **_kwargs: "test"
    background_config.resolve_particle_subtraction_mode = lambda *_args, **_kwargs: "simc_shape_components"
    background_config.resolve_pion_subtraction_scope = lambda *_args, **_kwargs: "t_bin"

    component_fits = types.ModuleType("pion_component_fits")
    component_fits.print_particle_subtraction_kaon_lambda_comparison_page = comparison_renderer
    component_fits.record_particle_subtraction_page = lambda *_args, **_kwargs: None
    for name in (
        "build_particle_subtraction_component_result",
        "load_or_resolve_pion_component_alignment",
        "print_particle_subtraction_component_application_pages",
        "print_particle_subtraction_component_fit_pages",
        "resolve_scope_component_shapes",
        "resolve_scope_single_shape",
    ):
        setattr(component_fits, name, lambda *_args, **_kwargs: None)

    component_subtraction = types.ModuleType("pion_component_subtraction")
    for name in (
        "build_simc_shape_pion_control_weights",
        "build_t_bin_pion_parent_identity",
        "evaluate_particle_subtraction_component_fit_result",
        "fingerprint_histogram_content_error",
        "resolve_frozen_parent_application_policy",
        "validate_authoritative_t_bin_pion_parent",
        "validate_frozen_t_bin_pion_parent_collection",
    ):
        setattr(component_subtraction, name, lambda *_args, **_kwargs: None)

    root_ownership = types.ModuleType("root_histogram_ownership")
    root_ownership.clone_root_histogram = lambda histogram, **_kwargs: histogram
    coordinates = types.ModuleType("data_coordinates")
    coordinates.validate_kaon_data_coordinate_contract = lambda *_args, **_kwargs: {}
    modules = {
        "background_config": background_config,
        "pion_component_fits": component_fits,
        "pion_component_subtraction": component_subtraction,
        "root_histogram_ownership": root_ownership,
        "data_coordinates": coordinates,
    }
    module_path = REPO_ROOT / "src" / "cuts" / "pion_t_bin_parents.py"
    with mock.patch.dict(sys.modules, modules):
        spec = importlib.util.spec_from_file_location("_c43_parent_lambda_renderer", module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


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

    def test_left_low_checkpoint_fixture_retains_both_recorded_candidates(self):
        payload = plots.method_b_display_payload(None, _left_low_checkpoint())
        candidates = plots.method_b_candidate_points(payload)
        self.assertEqual(payload["source"], "checkpoint_method_b")
        self.assertEqual(len(candidates), 2)
        self.assertEqual(candidates, [
            {
                "t_index": 0, "delta_index": 4, "delta_center": 3.5,
                "candidate_L_B": 0.9547608416143295,
                "candidate_L_B_uncertainty": 0.12140029583276699,
                "candidate_L_B_status": "available_multi_region",
            },
            {
                "t_index": 1, "delta_index": 2, "delta_center": -2.5,
                "candidate_L_B": 1.0713191848003094,
                "candidate_L_B_uncertainty": 0.20951468382894808,
                "candidate_L_B_status": "available_multi_region",
            },
        ])
        self.assertEqual(
            plots.method_b_candidate_count_parity(payload),
            {"expected_count": 2, "selected_count": 2, "passed": True},
        )

    def test_checkpoint_wins_over_stale_runtime_without_field_merging(self):
        runtime = types.MappingProxyType({
            "status": "available",
            "available": True,
            "summary": {"method_B_status_counts": {"available": 2}},
            "cells": [types.MappingProxyType({
                "t_index": 0, "delta_index": 0,
                "method_B_status": "available",
                "candidate_L_B_status": "unavailable",
            })],
        })
        payload = plots.method_b_display_payload(runtime, _left_low_checkpoint())
        self.assertEqual(payload["source"], "checkpoint_method_b")
        self.assertTrue(payload["source_complete"])
        self.assertEqual(len(payload["cells"]), 30)
        self.assertEqual(len(plots.method_b_candidate_points(payload)), 2)
        self.assertEqual(payload["method_status_counts"], {
            "available": 2, "marginal": 3, "shape_inconsistent": 0, "unavailable": 25,
        })
        self.assertEqual(payload["summary"]["candidate_status_counts"]["available_multi_region"], 2)

    def test_runtime_method_b_is_used_only_when_checkpoint_method_b_is_absent(self):
        runtime = types.MappingProxyType(_left_low_method_b())
        payload = plots.method_b_display_payload(runtime, {
            "canonical_t_edges": [0.1, 0.2, 0.3, 0.4],
            "delta_edges": [-10.0, -7.0, -4.0, -1.0, 2.0, 5.0, 8.0, 11.0, 14.0, 17.0, 20.0],
        })
        self.assertEqual(payload["source"], "runtime_method_b_fallback")
        self.assertEqual(len(plots.method_b_candidate_points(payload)), 2)
        self.assertEqual(payload["method_status_counts"]["unavailable"], 25)

    def test_runtime_checkpoint_presentation_mismatch_keeps_checkpoint_display(self):
        runtime = _left_low_method_b()
        runtime["cells"][4]["candidate_L_B"] = 9.99
        checkpoint = _left_low_checkpoint()
        payload = plots.method_b_display_payload(runtime, checkpoint)
        parity = plots.method_b_display_source_parity(runtime, checkpoint)
        self.assertEqual(payload["source"], "checkpoint_method_b")
        self.assertEqual(plots.method_b_candidate_points(payload)[0]["candidate_L_B"], 0.9547608416143295)
        self.assertTrue(parity["checked"])
        self.assertFalse(parity["passed"])
        self.assertIn("candidate_tuples", parity["differences"])
        self.assertEqual(checkpoint["method_b"]["cells"][4]["candidate_L_B"], 0.9547608416143295)

    def test_method_b_coverage_retains_zero_categories_and_all_cells(self):
        payload = plots.method_b_display_payload(None, _left_low_checkpoint())
        self.assertEqual(payload["method_status_counts"], {
            "available": 2, "marginal": 3, "shape_inconsistent": 0, "unavailable": 25,
        })
        self.assertEqual(payload["coverage_parity"], {
            "checked": True, "passed": True, "cell_count": 30, "coverage_total": 30,
            "summary_checked": True, "summary_passed": True, "differences": (),
        })
        qa = plots.setting_qa_summary_payload(
            {"status": "available"}, {}, {}, method_b_display=payload,
        )
        self.assertEqual(qa["method_b_coverage"], payload["method_status_counts"])
        self.assertIn(
            "Method B coverage: available=2, marginal=3, shape_inconsistent=0, unavailable=25",
            plots.setting_qa_summary_lines(qa),
        )
        mismatched_checkpoint = _left_low_checkpoint()
        mismatched_checkpoint["method_b"]["summary"]["method_B_status_counts"]["available"] = 1
        mismatch = plots.method_b_display_payload(None, mismatched_checkpoint)["coverage_parity"]
        self.assertFalse(mismatch["passed"])
        self.assertIn("method_B_status_counts", mismatch["differences"])

    def test_method_b_coverage_keeps_explicit_other_statuses(self):
        method_b = _left_low_method_b()
        method_b["cells"][0]["method_B_status"] = "internally_inconsistent"
        method_b["summary"]["method_B_status_counts"] = {
            "available": 2, "marginal": 2, "shape_inconsistent": 0,
            "unavailable": 25, "internally_inconsistent": 1,
        }
        payload = plots.method_b_display_payload(method_b, {})
        self.assertEqual(payload["other_method_status_counts"], {"internally_inconsistent": 1})
        self.assertEqual(sum(payload["method_status_counts"].values()) + sum(payload["other_method_status_counts"].values()), 30)
        self.assertTrue(payload["coverage_parity"]["passed"])

    def test_candidate_page_state_uses_local_empty_panels_without_global_empty_message(self):
        payload = plots.method_b_display_payload(None, _left_low_checkpoint())
        state = plots.method_b_candidate_page_state(payload)
        self.assertEqual(state["candidate_count"], 2)
        self.assertEqual(state["parent_candidate_counts"], {0: 1, 1: 1, 2: 0})
        self.assertFalse(state["show_setting_empty"])

        local_empty = _left_low_method_b()
        local_empty["cells"] = [
            cell for cell in local_empty["cells"]
            if not (cell["t_index"] == 1 and cell["delta_index"] == 2)
        ]
        local_empty["cells"].append({
            "t_index": 1, "delta_index": 2, "delta_low": -4.0, "delta_high": -1.0,
            "method_B_status": "unavailable", "candidate_L_B": None,
            "candidate_L_B_uncertainty": None, "candidate_L_B_status": "unavailable",
            "shape_status": "unavailable", "regions": [],
        })
        local_state = plots.method_b_candidate_page_state(
            plots.method_b_display_payload(local_empty, {})
        )
        self.assertEqual(local_state["parent_candidate_counts"], {0: 1, 1: 0, 2: 0})
        self.assertFalse(local_state["show_setting_empty"])

    def test_one_display_payload_drives_method_b_pages_and_final_qa(self):
        payload = plots.method_b_display_payload(None, _left_low_checkpoint())
        qa = plots.setting_qa_summary_payload(
            {"status": "available"}, {}, {}, method_b_display=payload,
        )
        self.assertEqual(len(plots.method_b_shape_rows(payload)), 30)
        self.assertEqual(plots.method_b_candidate_page_state(payload)["candidate_count"], 2)
        self.assertEqual(qa["method_b_coverage"], payload["method_status_counts"])
        self.assertEqual(plots.unity_line_limits(payload), (-10.0, 20.0))

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

    def test_proton_cleaned_main_qa_retains_existing_compact_spectra_and_closures(self):
        raw, estimated, cleaned, committed = object(), object(), object(), object()
        qa = plots.proton_main_qa_payload(
            {
                "diagnostics": {
                    "lambda_preservation_gate": {
                        "proposed_pre_rf_closure_passed": True,
                        "final_applied_closure_passed": True,
                    },
                },
            },
            {
                "host_state": "proton_cleaned",
                "H_MM_before_proton_cleaning": raw,
                "H_MM_estimated_proton": estimated,
                "H_MM_after_proton_cleaning": cleaned,
                "H_MM_after_proton_cleaning_final_rf": committed,
                "diagnostics": {"canonical_t_global_closure": {"raw": {"passed": True}}},
            },
        )
        self.assertEqual(qa["closure_mode"], "proton_cleaned")
        self.assertEqual(qa["spectra"], {
            "raw": raw, "estimated": estimated, "cleaned": cleaned, "committed": committed,
        })
        self.assertEqual(qa["numerical_closure"]["proposed_pre_rf_closure_passed"], True)
        self.assertEqual(qa["global_closure"], {"raw": {"passed": True}})
        self.assertEqual(qa["identity_host_closure"], "not recorded")

    def test_identity_host_main_qa_uses_existing_target_maps_and_identity_closure(self):
        raw, proton_zero, cleaned, final = object(), object(), object(), object()
        identity_closure = {
            "passed": True,
            "identity_transform_closure": {
                "passed": True, "global_full": {"passed": True}, "global_cut": {"passed": True},
            },
            "upstream_noRF_closure": {"passed": True},
        }
        qa = plots.proton_main_qa_payload(
            {"diagnostics": {"lambda_preservation_gate": {"status": "bypass"}}},
            {
                "host_state": "identity_no_proton_cleaning",
                "raw_targets": {"h_mm_nosub": raw},
                "proton_targets": {"h_mm_nosub": proton_zero},
                "cleaned_targets_pre_rf": {"h_mm_nosub": cleaned},
                "final_targets": {"h_mm_nosub": final},
                "diagnostics": {"identity_host_closure": identity_closure},
            },
            {
                "host_state": "identity_no_proton_cleaning",
                "production_action": "bypass",
                "proton_cleaning_committed": False,
            },
        )
        self.assertEqual(qa["closure_mode"], "identity_no_proton_cleaning")
        self.assertEqual(qa["spectra"], {
            "raw": raw, "estimated": proton_zero, "cleaned": cleaned, "committed": final,
        })
        self.assertEqual(qa["numerical_closure"], "not recorded")
        self.assertEqual(qa["global_closure"], "not recorded")
        self.assertIs(qa["identity_host_closure"], identity_closure)
        lines = plots.proton_main_summary_lines(qa)
        self.assertIn("host state: identity_no_proton_cleaning", lines)
        self.assertIn("production action: bypass; committed=False", lines)
        self.assertIn("Identity-host closure", lines)
        self.assertIn("identity-transform: PASS", lines)
        self.assertIn("No proton subtraction was applied.", lines)

    def test_proton_numerical_and_canonical_global_closures_remain_distinct(self):
        qa = plots.proton_main_qa_payload(
            {
                "diagnostics": {
                    "lambda_preservation_gate": {
                        "proposed_pre_rf_closure_passed": True,
                        "proposed_pre_rf_closure_difference": 0.02,
                        "final_applied_closure_passed": False,
                        "final_applied_pre_rf_closure_difference": -0.03,
                    },
                },
            },
            {
                "diagnostics": {
                    "canonical_t_global_closure": {
                        "raw": {"passed": True},
                        "proton_estimate": {"passed": True},
                        "proton_cleaned_pre_rf": {"passed": True},
                        "final_post_rf": {"passed": False},
                    },
                },
            },
        )
        self.assertEqual(qa["numerical_closure"]["proposed_pre_rf_closure_difference"], 0.02)
        self.assertEqual(qa["numerical_closure"]["final_applied_pre_rf_closure_difference"], -0.03)
        self.assertNotEqual(qa["numerical_closure"], qa["global_closure"])
        lines = plots.proton_closure_summary_lines(qa)
        self.assertIn("Numerical proton closure", lines)
        self.assertIn("Canonical-|t| global closure", lines)
        self.assertIn("final: FAIL", lines)

    def test_warning_classifier_blocks_a_clean_statement_for_unavailable_hgcer(self):
        summary = plots.setting_qa_summary_payload(
            {"status": "available", "host_state": "post_proton_noRF", "pion_closure": {"passed": True}, "host_closure": {"passed": True}},
            _method_a(),
            _method_b(),
            {"status": "available"},
            display_context={"lambda_gate_status": "pass", "production_action": "apply", "proton_cleaning_committed": True},
            runtime_qa_context={
                "aerogel_warnings": [], "proton_warnings": [],
                "canonical_parent_k_lambda": [{"t_index": 0, "status": "available"}],
                "k_sigma0_protected_region": "active",
                "k_sigma0_availability": "not available",
                "hgcer_diagnostic_availability": "unavailable",
                "renderer_failures": [],
            },
        )
        states = plots.setting_qa_warning_states(summary)
        self.assertIn(
            {"label": "HGCer diagnostics", "category": "FAILURE", "detail": "unavailable"},
            states,
        )
        lines = plots.setting_qa_summary_lines(summary)
        self.assertIn("Outstanding setting-level QA states:", lines)
        self.assertNotIn("No outstanding setting-level QA warnings.", lines)
        self.assertTrue(any(
            line.startswith("INFORMATIONAL: K-Sigma0 explicit template")
            and "protected region=active; explicit template=not available" in line
            for line in lines
        ))

    def test_warning_classifier_blocks_a_clean_statement_for_proton_renderer_failure(self):
        summary = plots.setting_qa_summary_payload(
            {"status": "available", "pion_closure": {"passed": True}, "host_closure": {"passed": True}},
            _method_a(),
            _method_b(),
            {"status": "available"},
            runtime_qa_context={
                "aerogel_warnings": [],
                "proton_warnings": [],
                "canonical_parent_k_lambda": [{"t_index": 0, "status": "available"}],
                "hgcer_diagnostic_availability": "available",
                "renderer_failures": ["proton main summary: RuntimeError: test"],
            },
        )
        lines = plots.setting_qa_summary_lines(summary)
        self.assertIn("Outstanding setting-level QA states:", lines)
        self.assertNotIn("No outstanding setting-level QA warnings.", lines)
        self.assertTrue(any(
            line.startswith("FAILURE: Required main-PDF renderer") and "proton main summary" in line
            for line in lines
        ))

    def test_canonical_parent_lambda_summary_keeps_parent_identity_and_reason(self):
        summary = plots.canonical_parent_lambda_summary([
            {"t_index": 0, "status": "available"},
            {"t_index": 1, "status": "available"},
            {"t_index": 2, "status": "unavailable", "reason": "reference missing"},
        ])
        self.assertEqual((summary["available"], summary["total"], summary["unavailable"]), (2, 3, 1))
        self.assertEqual(summary["entries"][2], {"t_index": 2, "status": "unavailable", "reason": "reference missing"})
        lines = plots.setting_qa_summary_lines({
            "method_a_coverage": {}, "method_b_coverage": {}, "method_b_other_coverage": {},
            "phase_host_state": "not recorded", "lambda_gate_status": "not recorded",
            "proton_action": "not recorded", "proton_cleaning_committed": "not recorded",
            "phase_a_coordinate_status": "available", "phase_a_pion_closure": True,
            "phase_a_host_closure": True, "canonical_parent_k_lambda": summary,
            "k_sigma0_protected_region": "not recorded", "k_sigma0_availability": "not available",
            "aerogel_warnings": [], "proton_warnings": [],
            "hgcer_diagnostic_availability": "available", "renderer_failures": [],
        })
        self.assertIn("  t1: available", lines)
        self.assertIn("  t3: unavailable — reference missing", lines)

    def test_parent_lambda_render_outcome_uses_actual_comparison_result_not_provenance(self):
        module = _load_parent_lambda_renderer(lambda *_args, **_kwargs: True)
        parent = {
            "t_bin_index": 1,
            "analysis_scope": "canonical_t_1",
            "fit_result": {
                "diagnostics": {
                    "pi_delta_signal_protected_fit": {
                        "k_lambda_scope_template_availability": {"status": "unavailable"},
                    },
                },
            },
        }
        record = module._print_parent_lambda_page(
            "synthetic.pdf", parent, {"mm_min": 1.0, "mm_max": 1.3}, "test", [], "pion.t_bin.1"
        )
        self.assertEqual(record, {
            "t_index": 1, "status": "available", "reason": None, "page_recorded": True,
        })
        with mock.patch("builtins.print"):
            report = module._parent_plot_contract(
                [], [parent], setting_wide_enabled=False, canonical_parent_k_lambda_render=[record]
            )
        self.assertEqual(report["canonical_parent_k_lambda_render"], [record])

    def test_parent_lambda_unavailable_status_page_is_not_a_successful_comparison(self):
        def unavailable(*_args, **_kwargs):
            raise RuntimeError("reference unavailable")

        module = _load_parent_lambda_renderer(unavailable)
        parent = {"t_bin_index": 2, "analysis_scope": "canonical_t_2"}
        with mock.patch.object(module, "_print_parent_status_page", return_value=True):
            record = module._print_parent_lambda_page(
                "synthetic.pdf", parent, {"mm_min": 1.0, "mm_max": 1.3}, "test", [], "pion.t_bin.2"
            )
        self.assertEqual(record, {
            "t_index": 2,
            "status": "unavailable",
            "reason": "reference unavailable",
            "page_recorded": True,
        })

    def test_canonical_delta_frame_uses_edges_not_candidate_extent(self):
        payload = {"delta_edges": [-10.0, -7.0, -4.0, 5.0, 20.0]}
        points = [{"delta_center": -4.0}, {"delta_center": 5.0}]
        self.assertEqual(plots.canonical_delta_frame_limits(payload), (-10.0, 20.0))
        self.assertEqual(plots.unity_line_limits(payload), (-10.0, 20.0))
        self.assertEqual([point["delta_center"] for point in points], [-4.0, 5.0])

    def test_common_y_ranges_include_all_error_bars_and_unity(self):
        regional = {
            "series": {
                "pi_n": [{"Qtilde": 0.8, "Qtilde_uncertainty": 0.1}],
                "pi_sidis": [{"Qtilde": 1.0, "Qtilde_uncertainty": 0.2}],
                "pi_delta_high": [{"Qtilde": 1.7, "Qtilde_uncertainty": 0.3}],
            },
        }
        ymin, ymax = plots.method_b_regional_frame_limits(regional)
        self.assertLessEqual(ymin, 0.7)
        self.assertGreaterEqual(ymax, 2.0)
        self.assertLessEqual(ymin, 1.0)
        self.assertGreaterEqual(ymax, 1.0)
        a_min, a_max = plots.method_a_f_low_frame_limits([
            {"f_low": 0.3, "f_low_low": 0.1, "f_low_high": 0.4},
            {"f_low": 0.7, "f_low_low": 0.5, "f_low_high": 0.95},
        ])
        self.assertLessEqual(a_min, 0.1)
        self.assertGreaterEqual(a_max, 0.95)

    def test_renderer_retains_draw_objects_and_shape_poor_spelling(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_refinement_plots.py").read_text(encoding="utf-8")
        self.assertGreaterEqual(source.count("draw_objects"), 5)
        self.assertIn("shape-poor veto", source)
        self.assertNotIn("shape-pool veto", source)

    def test_rand_sub_contains_proton_renderer_failure_and_actual_parent_outcome_handoff(self):
        source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        proton_call = source.index("render_proton_main_summary_pages(")
        proton_failure = source.index("proton main summary: {}: {}")
        refinement_call = source.index("render_pion_hgcer_refinement_pages(")
        final_qa = source.index('"renderer_failures": setting_renderer_failures')
        warning_page = source.index("render_setting_warning_page(")
        method_b_display = source.index("method_b_display = method_b_display_payload(")
        self.assertLess(proton_call, proton_failure)
        self.assertLess(proton_failure, final_qa)
        self.assertLess(proton_failure, refinement_call)
        self.assertLess(refinement_call, warning_page)
        self.assertLess(method_b_display, final_qa)
        self.assertIn("Method-B display-source parity mismatch", source)
        self.assertIn("method_b_display=method_b_display", source)
        self.assertIn("canonical_parent_k_lambda_render", source)
        self.assertNotIn("k_lambda_scope_template_availability", source)
        self.assertNotIn("k_lambda_source_availability", source)

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
        self.assertTrue(set(imports).issubset({"__future__", "ROOT", "array", "collections.abc", "math", "os", "textwrap"}))
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
