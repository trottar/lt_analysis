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
        for child in children:
            self.assertEqual(len(child["one_dimensional"]), 6)
            self.assertEqual(set(child["joint_distributions"]), {"analysis_MM__SHMS_xptar", "analysis_MM__SHMS_yptar", "SHMS_delta__SHMS_xptar", "SHMS_delta__SHMS_yptar"})
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
        wrong = deepcopy(self.f6_authority); wrong[KINEMATIC]["source_file_sha256"] = "f" * 64
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_authority_source_file_sha256_mismatch"):
            self._build(accepted_f6_1_artifact_authority_by_kinematic=wrong)
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
        self.assertTrue(any(value < 0.0 for value in signed["analysis_MM"]["baseline_signed_contents"]))

    def test_bootstrap_seed_percentile_and_pairing_are_deterministic(self):
        self.assertEqual(f62._linear_percentile([0.0, 10.0], 25.0), 2.5)
        self.assertEqual(f62._child_seed("Left-lowe", 0, 1, 20260916), f62._child_seed("Left-lowe", 0, 1, 20260916))
        first = self._build(); second = self._build()
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        child = first["parents"][0]["children"][0]
        bootstrap = child["bootstrap"]
        self.assertEqual(bootstrap["policy"]["replicas"], 4)
        self.assertEqual(bootstrap["kaon_window"]["DeltaP_K"]["requested_replicas"], 4)
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "bootstrap_test_config_invalid"):
            self._build(bootstrap_test_config={"seed": 3})

    def test_aggregate_guard_and_artifact_fingerprint_are_deterministic(self):
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "forbidden_event_persistence:entry_index"):
            f62._assert_aggregate_only_persistence({"entry_index": 1})
        first = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(self.artifacts, self.f3, self.f4, self.f5, self.f6, f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        second = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(deepcopy(self.artifacts), deepcopy(self.f3), deepcopy(self.f4), deepcopy(self.f5), deepcopy(self.f6), f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        self.assertEqual(first["artifact_fingerprint"], second["artifact_fingerprint"])

    def test_module_has_no_production_or_method_b_import(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_method_a_acceptance_refinement_validation.py").read_text(encoding="utf-8")
        for forbidden in ("import ROOT", "rand_sub", "calculate_yield", "pion_component_subtraction", "pion_hgcer_refinement_method_b"):
            self.assertNotIn(forbidden, source)
        self.assertIn('"method_b_numerical_dependency": False', source)


if __name__ == "__main__":
    unittest.main()
