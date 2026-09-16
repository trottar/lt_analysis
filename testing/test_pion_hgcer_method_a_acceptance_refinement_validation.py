"""Focused coverage for detached Phase F.6.2 refinement validation."""

from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_acceptance_refinement_validation as f62
import pion_hgcer_method_a_reweighting_validation as f61
import test_pion_hgcer_method_a_reweighting_validation as f61_fixtures


KINEMATIC = f61_fixtures.KINEMATIC


def _chain() -> tuple[object, ...]:
    artifacts, hashes, f3, f4_artifact, f5_artifact, authority, f3_authority, f4_authority = f61_fixtures._chain()
    f6_artifact = f61.build_pion_hgcer_method_a_reweighting_validation_artifact(artifacts, f3, f4_artifact, f5_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, input_paths={"f1": {}, "f3": "three", "f4": "four", "f5": "five"}, accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
    f6_authority = {KINEMATIC: {"source_file_sha256": "e" * 64, "validation_fingerprint": f6_artifact["validation"]["fingerprint"], "artifact_fingerprint": f6_artifact["artifact_fingerprint"]}}
    return artifacts, hashes, f3, f4_artifact, f5_artifact, authority, f3_authority, f4_authority, f6_artifact, f6_authority


class AcceptanceRefinementValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (cls.artifacts, cls.hashes, cls.f3, cls.f4, cls.f5, cls.authority, cls.f3_authority, cls.f4_authority, cls.f6, cls.f6_authority) = _chain()

    def _build(self, **changes):
        values = {
            "f1_artifacts": self.artifacts, "f3_artifact": self.f3, "f4_artifact": self.f4, "f5_artifact": self.f5, "f6_1_artifact": self.f6,
            "f1_input_file_hashes": self.hashes, "f3_input_file_sha256": f61_fixtures.F3_SHA, "f4_input_file_sha256": f61_fixtures.F4_SHA, "f5_input_file_sha256": f61_fixtures.F5_SHA, "f6_1_input_file_sha256": "e" * 64,
            "accepted_runtime_authority_by_kinematic": self.authority, "accepted_f4_runtime_authority_by_kinematic": self.f4_authority, "accepted_f3_runtime_authority_by_kinematic": self.f3_authority, "accepted_f6_1_artifact_authority_by_kinematic": self.f6_authority, "bootstrap_test_config": {"replicas": 4},
        }
        values.update(changes)
        return f62.build_pion_hgcer_method_a_acceptance_refinement_validation(**values)

    def test_valid_chain_builds_all_children_and_only_aggregate_content(self):
        result = self._build()
        self.assertTrue(result["available"]); self.assertEqual(len(result["parents"]), 15)
        children = [child for parent in result["parents"] for child in parent["children"]]
        self.assertEqual(len(children), 135)
        self.assertTrue(result["f6_1_reproduction"]["exact_payload_match"])
        self.assertTrue(result["f4_reproduction"]["exact_payload_match"])
        self.assertTrue(result["f5_reproduction"]["exact_payload_match"])
        self.assertEqual(result["f6_1_authority"]["accepted"]["source_file_sha256"], "e" * 64)
        self.assertEqual(result["parents"][0]["children"][0]["joint_distributions"]["SHMS_delta__SHMS_xptar"]["x_edges"], self.artifacts[0]["contract"]["delta_edges"])
        for child in children:
            self.assertEqual(set(child["availability"]), {"has_low_response", "has_prompt_control", "has_full_application", "completely_empty"})
            self.assertEqual(len(child["one_dimensional"]), 6)
            self.assertEqual(set(child["joint_distributions"]), {"analysis_MM__SHMS_xptar", "analysis_MM__SHMS_yptar", "SHMS_delta__SHMS_xptar", "SHMS_delta__SHMS_yptar"})
            self.assertEqual(set(child["population_counts"]["full_by_source"]), {"prompt", "rand", "dummy", "dummy_rand"})
            for shape in child["one_dimensional"].values():
                self.assertEqual(len(shape["edges"]), 41)
                for population in ("L", "B", "A"):
                    if shape[population]["available"]:
                        self.assertAlmostEqual(sum(shape[population]["unit_area"]), 1.0)
        serialized = json.dumps(result, sort_keys=True)
        for forbidden in ("entry_index", "correction_factors", "raw_shape_factors", "in_support_mask", "application_records", "method_a_training_records"):
            self.assertNotIn(forbidden, serialized)
        self.assertFalse(result["production_application_performed"]); self.assertFalse(result["method_b_numerical_dependency"])

    def test_authority_f6_reproduction_and_transient_factor_alignment_fail_closed(self):
        for field in ("source_file_sha256", "validation_fingerprint", "artifact_fingerprint"):
            wrong = deepcopy(self.f6_authority); wrong[KINEMATIC][field] = "f" * 64
            with self.subTest(field=field), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_authority_{}_mismatch".format(field)):
                self._build(accepted_f6_1_artifact_authority_by_kinematic=wrong)
        changed_upstream = deepcopy(self.authority); changed_upstream[KINEMATIC]["f3_source_file_sha256"] = "f" * 64
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
            self._build(accepted_runtime_authority_by_kinematic=changed_upstream)
        corrupt = deepcopy(self.f6); corrupt["validation"]["parents"] = []
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_validation_fingerprint_invalid|f6_1_exact"):
            self._build(f6_1_artifact=corrupt)
        original = f62._f61._validate_and_pair_rows
        with mock.patch.object(f62._f61, "_validate_and_pair_rows", side_effect=ValueError("factor alignment")), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
            self._build()
        self.assertIsNotNone(original)

    def test_metrics_states_bins_and_signed_window_are_explicit(self):
        shape = f62._one_dimensional_payload([], [], [], [], [0.0, 1.0], "empty")
        self.assertFalse(shape["metrics"]["kappa"]["available"]); self.assertEqual(shape["metrics"]["kappa"]["reason"], "low_population_empty")
        correction_zero = f62._one_dimensional_payload([0.1, 0.1], [0.9, 0.9], [1.0, 1.0], [1.0, 1.0], [0.0, 0.5, 1.0], "zero")
        self.assertEqual(correction_zero["metrics"]["kappa"]["reason"], "correction_norm_zero")
        self.assertEqual(f62._coarsen_edges(list(range(41)), "edges"), [float(item) for item in range(0, 41, 4)])
        window = {"mm_min": 1.10, "mm_max": 1.16}
        rows = [{"row": {"source_label": "prompt", "analysis_MM": 1.10, "signed_baseline_event_contribution": 2.0}, "factor": 1.5}, {"row": {"source_label": "rand", "analysis_MM": 1.16, "signed_baseline_event_contribution": -3.0}, "factor": 2.0}]
        signed, _ = f62._signed_payload(rows, [1.0, 1.2], window, "window")
        self.assertEqual(signed["kaon_window"]["P_B_K"]["value"], 2.0)
        self.assertEqual(signed["kaon_window"]["P_A_K"]["value"], 3.0)
        self.assertEqual(signed["kaon_window"]["DeltaP_K"]["value"], 1.0)
        self.assertEqual(signed["variance_proxy"]["V_B"]["value"], 13.0)
        self.assertEqual(signed["variance_proxy"]["V_A"]["value"], 45.0)
        self.assertAlmostEqual(signed["variance_proxy"]["R_V"]["value"], 45.0 / 13.0)
        self.assertTrue(any(value < 0.0 for value in signed["analysis_MM"]["baseline_signed_contents"]))
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "signed_source"):
            f62._signed_payload([{"row": {"source_label": "unknown", "analysis_MM": 1.12, "signed_baseline_event_contribution": 1.0}, "factor": 1.0}], [1.0, 1.2], window, "unknown")
        self.assertAlmostEqual(f62._effective_sample_size([1.0, 2.0], "neff")["value"], 1.8)
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "baseline_weight_invalid"):
            f62._one_dimensional_payload([0.1], [0.2], [-1.0], [1.0], [0.0, 1.0], "nonpositive")

    def test_hand_metrics_and_bootstrap_contract_states(self):
        low = {"available": True, "unit_area": [1.0, 0.0], "reason": None}
        baseline = {"available": True, "unit_area": [0.0, 1.0], "reason": None}
        method_a = {"available": True, "unit_area": [0.5, 0.5], "reason": None}
        metrics = f62._comparison_metrics(low, baseline, method_a)
        self.assertAlmostEqual(metrics["H_B"]["value"], 1.0)
        self.assertAlmostEqual(metrics["TV_LB"]["value"], 1.0)
        self.assertAlmostEqual(metrics["R"]["value"], 0.5)
        self.assertAlmostEqual(metrics["kappa"]["value"], 1.0)
        self.assertAlmostEqual(metrics["rho"]["value"], 0.5)
        self.assertLess(metrics["DeltaH"]["value"], 0.0)
        zero_residual = f62._comparison_metrics(low, low, method_a)
        self.assertEqual(zero_residual["kappa"]["reason"], "residual_norm_zero")
        matrix_low = {"available": True, "unit_area": [[1.0, 0.0], [0.0, 0.0]], "reason": None}
        matrix_baseline = {"available": True, "unit_area": [[0.0, 1.0], [0.0, 0.0]], "reason": None}
        matrix_a = {"available": True, "unit_area": [[0.5, 0.5], [0.0, 0.0]], "reason": None}
        self.assertAlmostEqual(f62._joint_metrics(matrix_low, matrix_baseline, matrix_a)["kappa"]["value"], 1.0)
        policy = dict(f62._BOOTSTRAP_POLICY)
        insufficient = f62._bootstrap_summary([1.0] * 8, 10, "kappa", policy)
        self.assertFalse(insufficient["interval_available"]); self.assertEqual(insufficient["reason"], "insufficient_valid_bootstrap_replicas")
        self.assertEqual(set(insufficient), {"requested_replica_count", "valid_replica_count", "invalid_replica_count", "interval_available", "reason", "ci_low", "ci_high"})

    def test_bootstrap_seed_percentile_and_pairing_are_deterministic(self):
        self.assertEqual(f62._linear_percentile([0.0, 10.0], 25.0), 2.5)
        self.assertEqual(f62._child_seed("Left-lowe", 0, 1, 20260916), f62._child_seed("Left-lowe", 0, 1, 20260916))
        first = self._build(); second = self._build()
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        child = first["parents"][0]["children"][0]
        bootstrap = child["bootstrap"]
        self.assertEqual(bootstrap["policy"]["replicas"], 4)
        self.assertEqual(bootstrap["kaon_window"]["DeltaP_K"]["requested_replica_count"], 4)
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "bootstrap_test_config_invalid"):
            self._build(bootstrap_test_config={"seed": 3})

    def test_aggregate_guard_and_artifact_fingerprint_are_deterministic(self):
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "forbidden_event_persistence:entry_index"):
            f62._assert_aggregate_only_persistence({"entry_index": 1})
        first = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(self.artifacts, self.f3, self.f4, self.f5, self.f6, f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        second = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(deepcopy(self.artifacts), deepcopy(self.f3), deepcopy(self.f4), deepcopy(self.f5), deepcopy(self.f6), f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        self.assertEqual(first["artifact_fingerprint"], second["artifact_fingerprint"])
        for name, expected in {"non_authoritative": True, "validation_only": True, "manual_review_required": True, "accepted_f6_1_consumed": True, "accepted_f6_1_modified": False, "event_correction_evaluated_for_detached_validation": True, "event_correction_persisted": False, "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "method_b_numerical_dependency": False, "automatic_case_classification": False, "case_thresholds_defined": False, "final_yield_uncertainty_claimed": False}.items():
            self.assertEqual(first[name], expected)

    def test_module_has_no_production_or_method_b_import(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_method_a_acceptance_refinement_validation.py").read_text(encoding="utf-8")
        for forbidden in ("import ROOT", "rand_sub", "calculate_yield", "pion_component_subtraction", "pion_hgcer_refinement_method_b"):
            self.assertNotIn(forbidden, source)
        self.assertIn('"method_b_numerical_dependency": False', source)


if __name__ == "__main__":
    unittest.main()
