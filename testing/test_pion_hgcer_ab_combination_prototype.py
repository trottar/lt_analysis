"""Focused contracts for the detached equal-log A/B prototype."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import inspect
import json
import math
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import pion_hgcer_ab_combination_prototype as prototype


def _fingerprint(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False,
    ).encode("ascii")).hexdigest()


def _method_a(present, candidate, low, high, status):
    return {
        "present": present,
        "comparison_candidate": candidate,
        "comparison_candidate_low": low,
        "comparison_candidate_high": high,
        "comparison_candidate_status": status,
    }


def _method_b(present, candidate, uncertainty, status):
    return {
        "present": present,
        "comparison_candidate": candidate,
        "comparison_candidate_uncertainty": uncertainty,
        "comparison_candidate_status": status,
    }


def _checkpoint():
    """Return one complete Phase-D checkpoint with all frozen D.4 states."""
    t_edges = [0.0, 1.0]
    delta_edges = [-5.0, -3.0, -1.0, 1.0, 3.0, 5.0]
    source = "phase-c-source"
    phase = "phase-a-contract"
    coordinate = "coordinate-contract"
    definitions = (
        (
            "both_comparable",
            _method_a(True, 4.0, 3.0, 5.0, "marginal"),
            _method_b(True, 9.0, 1.0, "available_multi_region"),
            2.25, math.log(2.25), "disjoint",
        ),
        (
            "both_present_not_comparable",
            _method_a(True, 0.0, 0.0, 0.25, "marginal"),
            _method_b(True, 1.2, 0.2, "available_multi_region"),
            None, None, "not_evaluable",
        ),
        (
            "a_only",
            _method_a(True, 1.1, 0.8, 1.4, "available"),
            _method_b(False, None, None, "unavailable"),
            None, None, "not_evaluable",
        ),
        (
            "b_only",
            _method_a(False, None, None, None, "unavailable"),
            _method_b(True, 0.8, 0.1, "available_multi_region"),
            None, None, "not_evaluable",
        ),
        (
            "neither_available",
            _method_a(False, None, None, None, "unavailable"),
            _method_b(False, None, None, "shape_poor_veto"),
            None, None, "not_evaluable",
        ),
    )
    cells = []
    for delta_index, (availability, method_a, method_b, ratio, log_ratio, relation) in enumerate(definitions):
        cells.append({
            "t_index": 0,
            "t_low": t_edges[0],
            "t_high": t_edges[1],
            "delta_index": delta_index,
            "delta_low": delta_edges[delta_index],
            "delta_high": delta_edges[delta_index + 1],
            "method_a": method_a,
            "method_b": method_b,
            "comparison": {
                "availability": availability,
                "availability_reason": None,
                "ratio_B_over_A": ratio,
                "log_ratio_B_over_A": log_ratio,
                "diagnostic_interval_relation": relation,
            },
        })
    method_a = {
        "schema_version": "pion_hgcer_method_a_comparison/v1",
        "method": "method_a_same_t_comparison_representation",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source,
        "phase_a_contract_fingerprint": phase,
        "coordinate_fingerprint": coordinate,
        "source_method_a_payload_fingerprint": "method-a-source",
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "fingerprint_inputs": {"representation": "method-a"},
        "non_authoritative": True, "method_b_numerical_dependency": False,
        "comparison_performed": False, "classification_performed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    method_a["fingerprint"] = _fingerprint(method_a["fingerprint_inputs"])
    method_b = {
        "schema_version": "pion_hgcer_method_b_comparison/v1",
        "method": "method_b_comparison_representation",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source,
        "phase_a_contract_fingerprint": phase,
        "coordinate_fingerprint": coordinate,
        "source_method_b_payload_fingerprint": "method-b-source",
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "host_state": "proton_cleaned", "source_target_state": "post_proton_noRF",
        "fingerprint_inputs": {"representation": "method-b"},
        "non_authoritative": True, "method_a_numerical_dependency": False,
        "comparison_performed": False, "classification_performed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    method_b["fingerprint"] = _fingerprint(method_b["fingerprint_inputs"])
    comparison = {
        "schema_version": "pion_hgcer_ab_comparison/v1",
        "method": "non_authoritative_ab_comparison",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source,
        "phase_a_contract_fingerprint": phase,
        "coordinate_fingerprint": coordinate,
        "method_a_comparison_fingerprint": method_a["fingerprint"],
        "method_b_comparison_fingerprint": method_b["fingerprint"],
        "source_method_a_comparison_payload_fingerprint": _fingerprint(method_a),
        "source_method_b_comparison_payload_fingerprint": _fingerprint(method_b),
        "canonical_t_edges": t_edges, "delta_edges": delta_edges,
        "host_state": "proton_cleaned", "source_target_state": "post_proton_noRF",
        "cells": cells,
        "fingerprint_inputs": {"representation": "ab"},
        "non_authoritative": True, "comparison_performed": True,
        "classification_performed": True,
        "classification_scope": "availability_only_non_prescriptive",
        "decision_performed": False, "statistical_compatibility_claimed": False,
        "production_objects_mutated": False, "refinement_applied": False,
    }
    comparison["fingerprint"] = _fingerprint(comparison["fingerprint_inputs"])
    return {
        "schema_version": "pion_hgcer_phase_d_checkpoint/v1",
        "status": "available", "available": True,
        "source_checkpoint_payload_fingerprint": source,
        "non_authoritative": True, "comparison_performed": True,
        "classification_performed": True,
        "classification_scope": "availability_only_non_prescriptive",
        "decision_performed": False,
        "statistical_compatibility_claimed": False,
        "production_objects_mutated": False, "refinement_applied": False,
        "method_a_comparison": method_a,
        "method_b_comparison": method_b,
        "ab_comparison": comparison,
    }


class PionHGCerABCombinationPrototypeTests(unittest.TestCase):
    def test_equal_log_formula_states_and_marginal_input_are_copied(self):
        result = prototype.build_pion_hgcer_ab_combination_prototype(_checkpoint())

        self.assertTrue(result["available"])
        self.assertEqual(result["status"], "available")
        self.assertEqual(result["schema_version"], prototype.PION_HGCER_AB_COMBINATION_PROTOTYPE_SCHEMA_VERSION)
        self.assertEqual(result["method"], "equal_weight_log_geometric_mean")
        self.assertEqual((result["log_weight_method_a"], result["log_weight_method_b"]), (0.5, 0.5))
        self.assertFalse(result["prototype_statistical_optimality_claimed"])
        self.assertEqual(result["prototype_uncertainty_model"], "not_defined")

        comparable = result["cells"][0]
        self.assertEqual(comparable["method_a"]["status"], "marginal")
        self.assertEqual(comparable["comparison"]["diagnostic_interval_relation"], "disjoint")
        self.assertEqual(comparable["prototype_status"], "available")
        self.assertIsNone(comparable["prototype_reason"])
        self.assertAlmostEqual(comparable["prototype_log_scale"], 0.5 * (math.log(4.0) + math.log(9.0)))
        self.assertAlmostEqual(comparable["prototype_relative_scale"], 6.0)

        expected = (
            ("both_present_not_comparable", "frozen_ab_not_comparable"),
            ("a_only", "method_b_unavailable"),
            ("b_only", "method_a_unavailable"),
            ("neither_available", "both_methods_unavailable"),
        )
        for cell, (availability, reason) in zip(result["cells"][1:], expected):
            self.assertEqual(cell["comparison"]["availability"], availability)
            self.assertEqual(cell["prototype_status"], "unavailable")
            self.assertEqual(cell["prototype_reason"], reason)
            self.assertIsNone(cell["prototype_log_scale"])
            self.assertIsNone(cell["prototype_relative_scale"])
        self.assertNotIn("prototype_uncertainty", comparable)
        self.assertNotIn("combined_error", comparable)

        symmetric = _checkpoint()
        symmetric_cell = symmetric["ab_comparison"]["cells"][0]
        symmetric_cell["method_a"] = _method_a(True, 9.0, 8.0, 10.0, "marginal")
        symmetric_cell["method_b"] = _method_b(True, 4.0, 1.0, "available_multi_region")
        symmetric_cell["comparison"].update(
            ratio_B_over_A=4.0 / 9.0,
            log_ratio_B_over_A=math.log(4.0 / 9.0),
        )
        symmetric_result = prototype.build_pion_hgcer_ab_combination_prototype(symmetric)
        self.assertTrue(symmetric_result["available"])
        self.assertAlmostEqual(
            symmetric_result["cells"][0]["prototype_relative_scale"], 6.0
        )

    def test_contract_provenance_geometry_and_grid_failures_are_structured(self):
        cases = (
            (lambda value: value.update(schema_version="wrong"), "phase_d_checkpoint_schema_invalid"),
            (lambda value: value.update(available=False), "phase_d_checkpoint_unavailable"),
            (lambda value: value.update(non_authoritative=False), "phase_d_checkpoint_authority_invalid"),
            (lambda value: value.update(production_objects_mutated=True), "phase_d_checkpoint_authority_invalid"),
            (lambda value: value["ab_comparison"].update(coordinate_fingerprint="wrong"), "ab_comparison_provenance_mismatch"),
            (lambda value: value["ab_comparison"].update(method_a_comparison_fingerprint="wrong"), "ab_comparison_representation_fingerprint_mismatch"),
            (lambda value: value["ab_comparison"].update(source_method_b_comparison_payload_fingerprint="wrong"), "ab_comparison_representation_payload_fingerprint_mismatch"),
            (lambda value: value["ab_comparison"].update(delta_edges=[-5.0, 5.0]), "ab_comparison_geometry_mismatch"),
            (lambda value: value["ab_comparison"].update(cells=value["ab_comparison"]["cells"][:-1]), "ab_comparison_cell_grid_invalid"),
            (lambda value: value["ab_comparison"]["cells"].append(deepcopy(value["ab_comparison"]["cells"][0])), "ab_comparison_cell_grid_invalid"),
        )
        for mutation, expected_reason in cases:
            with self.subTest(expected_reason=expected_reason):
                checkpoint = _checkpoint()
                mutation(checkpoint)
                result = prototype.build_pion_hgcer_ab_combination_prototype(checkpoint)
                self.assertFalse(result["available"])
                self.assertEqual(result["status"], "unavailable")
                self.assertEqual(result["reason"], expected_reason)

    def test_input_is_immutable_result_is_detached_and_fingerprint_is_deterministic(self):
        checkpoint = _checkpoint()
        before = deepcopy(checkpoint)
        first = prototype.build_pion_hgcer_ab_combination_prototype(checkpoint)
        second = prototype.build_pion_hgcer_ab_combination_prototype(deepcopy(checkpoint))

        self.assertEqual(checkpoint, before)
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        first["cells"][0]["method_a"]["candidate"] = 77.0
        self.assertEqual(checkpoint["ab_comparison"]["cells"][0]["method_a"]["comparison_candidate"], 4.0)
        first["fingerprint_inputs"]["cells"][0]["method_a"]["candidate"] = 88.0
        self.assertEqual(first["cells"][0]["method_a"]["candidate"], 77.0)

        changed = _checkpoint()
        changed["ab_comparison"]["cells"][0]["method_a"]["comparison_candidate"] = 5.0
        changed["ab_comparison"]["cells"][0]["method_a"]["comparison_candidate_high"] = 6.0
        changed_result = prototype.build_pion_hgcer_ab_combination_prototype(changed)
        self.assertTrue(changed_result["available"])
        self.assertNotEqual(second["fingerprint"], changed_result["fingerprint"])

    def test_public_signature_and_scientific_boundaries(self):
        self.assertEqual(
            tuple(inspect.signature(prototype.build_pion_hgcer_ab_combination_prototype).parameters),
            ("phase_d_checkpoint",),
        )
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_ab_combination_prototype.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("math.log", source)
        self.assertIn("math.exp", source)
        for forbidden in (
            "find_canonical_bin(", "build_pion_hgcer_method_a(",
            "build_pion_hgcer_method_b(", "build_pion_hgcer_ab_comparison(",
            "numpy.corrcoef", "pearson", "spearman", "scipy.optimize", "curve_fit",
            "minimize", "interpolate", "interp1d", "polyfit", "spline",
            "preferred_method", "selected_method", "use_A", "use_B",
            "acceptance_correction",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
