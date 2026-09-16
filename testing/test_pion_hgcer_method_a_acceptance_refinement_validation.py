"""Focused coverage for detached Phase F.6.2 refinement validation."""

from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import tempfile
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

    def _prompt_child_inputs(self):
        pairs_by_parent, raw, _ = f62._reproduce_and_pair(self.artifacts, self.f3, self.hashes, f61_fixtures.F3_SHA, self.f4["correction"], self.f3_authority)
        for (setting_id, t_index), pairs in sorted(pairs_by_parent.items()):
            state = deepcopy(raw[setting_id]); state["setting_id"] = setting_id; _t_geometry, phi_edges, _delta_edges = f62._geometry(state["contract"], setting_id)
            for phi_index in range(9):
                low, control, _full = f62._child_rows(state, pairs, t_index, phi_index, phi_edges, "prompt_child")
                if low and control:
                    return state, deepcopy(pairs), t_index, phi_index, phi_edges, low, control
        self.fail("fixture has no populated prompt child")

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

    def test_all_frozen_upstream_authorities_fail_closed(self):
        changed_hashes = dict(self.hashes); first_hash = next(iter(changed_hashes)); changed_hashes[first_hash] = "f" * 64
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
            self._build(f1_input_file_hashes=changed_hashes)
        for field in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint", "f4_source_file_sha256", "f4_correction_fingerprint", "f4_artifact_fingerprint", "f5_source_file_sha256", "f5_propagation_fingerprint", "f5_artifact_fingerprint"):
            wrong = deepcopy(self.authority); wrong[KINEMATIC][field] = "f" * 64
            with self.subTest(authority=field), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
                self._build(accepted_runtime_authority_by_kinematic=wrong)
        for field in ("source_file_sha256", "map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"):
            wrong = deepcopy(self.f3_authority); wrong[KINEMATIC][field] = "f" * 64
            with self.subTest(f3_authority=field), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
                self._build(accepted_f3_runtime_authority_by_kinematic=wrong)
        for field in ("source_file_sha256", "correction_fingerprint", "artifact_fingerprint"):
            wrong = deepcopy(self.f4_authority); wrong[KINEMATIC][field] = "f" * 64
            with self.subTest(f4_authority=field), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
                self._build(accepted_f4_runtime_authority_by_kinematic=wrong)

    def test_f4_transient_review_factor_and_support_contract_fail_closed(self):
        original = f62._f61._f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data

        def broken_review(kind):
            def wrapped(*args, **kwargs):
                correction, review = original(*args, **kwargs)
                review = deepcopy(review)
                if kind == "factor_length":
                    review[0]["correction_factors"] = review[0]["correction_factors"][:-1]
                elif kind == "support_length":
                    review[0]["in_support_mask"] = review[0]["in_support_mask"][:-1]
                elif kind == "nonfinite_factor":
                    review[0]["correction_factors"][0] = float("nan")
                else:
                    review[0]["correction_factors"][0] = 0.0
                return correction, review
            return wrapped

        for kind in ("factor_length", "support_length", "nonfinite_factor", "nonpositive_factor"):
            with self.subTest(kind=kind), mock.patch.object(f62._f61._f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=broken_review(kind)):
                with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
                    self._build()

    def test_geometry_identity_and_phi_semantics_are_directly_closed(self):
        state, pairs, t_index, phi_index, phi_edges, low, control = self._prompt_child_inputs()
        for row in low:
            self.assertEqual(f62._phi_index(math.degrees(float(row["phi"])), phi_edges), phi_index)
        baseline, adjusted = f62._weights(control, "population")
        self.assertEqual(len(baseline), len(adjusted)); self.assertEqual(len(control), len(f62._values(control, "analysis_MM", training=False, label="population")))
        self.assertTrue(all(value > 0.0 for value in baseline)); self.assertTrue(all(value > 0.0 for value in adjusted))
        self.assertEqual(state["contract"]["delta_edges"], self.artifacts[[item["setting"]["phi_setting"] + "-" + item["setting"]["epsilon_filename_token"] for item in self.artifacts].index(state["setting_id"])]["contract"]["delta_edges"])
        prompt_pair = next(item for item in pairs if item["row"]["source_label"] == "prompt" and item["row"]["t_index"] == t_index and item["row"]["phi_index"] == phi_index)
        identity = ("prompt", prompt_pair["row"]["entry_index"])
        missing_training = deepcopy(state); del missing_training["training"][identity]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "prompt_identity_missing_training_control"):
            f62._child_rows(missing_training, pairs, t_index, phi_index, phi_edges, "missing_training")
        phi_mismatch = deepcopy(state); phi_mismatch["training"][identity]["phi"] = float(phi_mismatch["training"][identity]["phi"]) + 0.01
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "phi_semantic_parity_mismatch"):
            f62._child_rows(phi_mismatch, pairs, t_index, phi_index, phi_edges, "phi_mismatch")
        invalid_geometry_pairs = deepcopy(pairs)
        for item in invalid_geometry_pairs:
            if item["row"]["source_label"] == "prompt" and item["row"]["entry_index"] == identity[1]:
                item["row"]["phi_degrees"] = float(phi_edges[-1]) + 1.0
                break
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "application_phi_geometry_invalid"):
            f62._child_rows(state, invalid_geometry_pairs, t_index, phi_index, phi_edges, "invalid_geometry")

    def test_application_identity_and_prompt_training_mismatches_fail_closed(self):
        duplicate = deepcopy(self.artifacts)
        records = duplicate[0]["contract"]["application_records"]
        by_source: dict[str, list[dict[str, object]]] = {}
        for row in records:
            by_source.setdefault(str(row["source_label"]), []).append(row)
        same_source = next(group for group in by_source.values() if len(group) >= 2)
        same_source[1]["entry_index"] = same_source[0]["entry_index"]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "raw_f1_records_invalid"):
            f62._raw_state(duplicate)
        missing = deepcopy(self.artifacts)
        application = missing[0]["contract"]["application_records"]
        prompt_index = next(index for index, row in enumerate(application) if row["source_label"] == "prompt")
        del application[prompt_index]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f4_transient_review_reproduction_failed"):
            f62._reproduce_and_pair(missing, self.f3, self.hashes, f61_fixtures.F3_SHA, self.f4["correction"], self.f3_authority)

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

    def test_zero_control_weights_are_valid_but_negative_or_nonpositive_factors_are_not(self):
        one_dimensional = f62._one_dimensional_payload([0.1, 0.8], [0.2, 0.8], [0.0, 1.0], [0.0, 1.5], [0.0, 0.5, 1.0], "zero_weight_1d")
        for population in ("B", "A"):
            self.assertTrue(one_dimensional[population]["available"])
            self.assertAlmostEqual(sum(one_dimensional[population]["unit_area"]), 1.0)
        two_dimensional = f62._two_dimensional_payload([0.1, 0.8], [0.1, 0.8], [0.2, 0.8], [0.2, 0.8], [0.0, 1.0], [0.0, 1.5], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0], "zero_weight_2d")
        for population in ("B", "A"):
            self.assertTrue(two_dimensional[population]["available"])
            self.assertAlmostEqual(sum(sum(row) for row in two_dimensional[population]["unit_area"]), 1.0)
        all_zero = f62._one_dimensional_payload([0.1], [0.2, 0.8], [0.0, 0.0], [0.0, 0.0], [0.0, 0.5, 1.0], "all_zero_weight")
        self.assertFalse(all_zero["B"]["available"]); self.assertFalse(all_zero["A"]["available"])
        self.assertEqual(all_zero["B"]["reason"], "normalization_invalid")
        self.assertEqual(all_zero["A"]["reason"], "normalization_invalid")
        control = [{"row": {"baseline_pion_weight_w0": 0.0}, "factor": 1.5}]
        self.assertEqual(f62._weights(control, "zero_weight"), ([0.0], [0.0]))
        for bad_control in (
            [{"row": {"baseline_pion_weight_w0": -1.0}, "factor": 1.5}],
            [{"row": {"baseline_pion_weight_w0": 1.0}, "factor": 0.0}],
            [{"row": {"baseline_pion_weight_w0": 1.0}, "factor": -1.0}],
        ):
            with self.subTest(control=bad_control), self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "weight_invalid"):
                f62._weights(bad_control, "invalid_weight")
        self.assertEqual(f62._effective_sample_size([0.0, 1.0, 1.0], "mixed_zero") ["value"], 2.0)
        self.assertEqual(f62._effective_sample_size([0.0, 0.0], "all_zero_neff"), {"available": False, "value": None, "reason": "normalization_invalid"})
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "negative_neff_invalid"):
            f62._effective_sample_size([-1.0, 1.0], "negative_neff")
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "baseline_weight_invalid"):
            f62._two_dimensional_payload([0.1], [0.1], [0.2], [0.2], [-1.0], [1.0], [0.0, 1.0], [0.0, 1.0], "negative_2d")

    def test_zero_weight_prompt_control_is_retained_through_the_child_payload(self):
        state, pairs, t_index, phi_index, phi_edges, _low, original_control = self._prompt_child_inputs()
        zero_pairs = deepcopy(pairs)
        target = next(item for item in zero_pairs if item["row"]["source_label"] == "prompt" and item["row"]["t_index"] == t_index and item["row"]["phi_index"] == phi_index)
        target["row"]["baseline_pion_weight_w0"] = 0.0
        target["row"]["signed_baseline_event_contribution"] = 0.0
        low, controls, _full = f62._child_rows(state, zero_pairs, t_index, phi_index, phi_edges, "zero_weight_child")
        self.assertTrue(low); self.assertEqual(len(controls), len(original_control))
        baseline, adjusted = f62._weights(controls, "zero_weight_child")
        zero_index = next(index for index, item in enumerate(controls) if item["row"]["entry_index"] == target["row"]["entry_index"])
        self.assertEqual(baseline[zero_index], 0.0); self.assertEqual(adjusted[zero_index], 0.0)
        t_geometry, _phi_geometry, delta_edges = f62._geometry(state["contract"], state["setting_id"])
        parent_edges = f62._parent_display_edges(self.f6["validation"])[(state["setting_id"], t_index)]
        child = f62._child_payload(state, zero_pairs, parent_edges, t_index, phi_index, t_geometry[t_index], phi_edges, delta_edges, f62._bootstrap_policy({"replicas": 2}), f62._runtime_kaon_window())
        self.assertEqual(child["population_counts"]["N_control"], len(original_control))
        self.assertTrue(any(value == 0.0 for value in baseline))
        self.assertTrue(all(child["one_dimensional"][name]["B"]["available"] for name in f62._ONE_DIMENSIONAL_VARIABLES))
        self.assertTrue(all(child["one_dimensional"][name]["A"]["available"] for name in f62._ONE_DIMENSIONAL_VARIABLES))
        negative_pairs = deepcopy(zero_pairs); negative_target = next(item for item in negative_pairs if item["row"]["entry_index"] == target["row"]["entry_index"])
        negative_target["row"]["baseline_pion_weight_w0"] = -1.0
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "baseline_weight_invalid"):
            f62._child_rows(state, negative_pairs, t_index, phi_index, phi_edges, "negative_weight_child")

    def test_frozen_edges_finite_values_and_signed_inputs_fail_closed(self):
        display = deepcopy(self.f6["validation"])
        comparison = display["parents"][0]["prompt_shape_comparisons"]["analysis_MM"]
        comparison["edges"] = comparison["edges"][:-1]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_edges_analysis_MM"):
            f62._parent_display_edges(display)
        display = deepcopy(self.f6["validation"])
        comparison = display["parents"][0]["prompt_shape_comparisons"]["analysis_MM"]
        comparison["edges"][2] = comparison["edges"][1]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_edges_analysis_MM"):
            f62._parent_display_edges(display)
        contract = deepcopy(self.artifacts[0]["contract"]); self.assertEqual(f62._geometry(contract, "frozen")[2], self.artifacts[0]["contract"]["delta_edges"])
        contract["delta_edges"][1] = contract["delta_edges"][0]
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "frozen_delta_edges"):
            f62._geometry(contract, "frozen")
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "nonfinite_analysis_MM"):
            f62._value({"analysis_MM": float("nan")}, "analysis_MM", "nonfinite")
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "signed_source_or_factor_invalid"):
            f62._signed_payload([{"row": {"source_label": "prompt", "analysis_MM": 1.12, "signed_baseline_event_contribution": 1.0}, "factor": -1.0}], [1.0, 1.2], {"mm_min": 1.10, "mm_max": 1.16}, "signed_factor")
        changed = deepcopy(self.artifacts)
        row = next(item for item in changed[0]["contract"]["application_records"] if item["source_label"] == "prompt")
        row["signed_baseline_event_contribution"] = float(row["signed_baseline_event_contribution"]) + 1.0
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "f6_1_reproduction_failed"):
            self._build(f1_artifacts=changed)

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

    def test_bootstrap_uses_separate_low_paired_control_and_source_strata_draws(self):
        def training(value: float) -> dict[str, object]:
            return {name: value for name in f62._ONE_DIMENSIONAL_VARIABLES}

        def pair(source: str, value: float, factor: float) -> dict[str, object]:
            return {"row": {**training(value), "source_label": source, "baseline_pion_weight_w0": value, "signed_baseline_event_contribution": value}, "factor": factor}

        low = [training(1.0), training(2.0)]
        control = [pair("prompt", 1.0, 1.1), pair("prompt", 2.0, 1.2), pair("prompt", 3.0, 1.3)]
        full = [pair("prompt", 1.0, 1.1), pair("rand", 2.0, 1.2), pair("rand", 3.0, 1.3)]
        variable_edges = {name: [float(index) for index in range(41)] for name in f62._ONE_DIMENSIONAL_VARIABLES}
        policy = f62._bootstrap_policy({"replicas": 1})

        class FakeGenerator:
            def __init__(self):
                self.draws: list[tuple[int, int]] = []

            def integers(self, low_bound, high_bound, size):
                self.draws.append((int(high_bound), int(size)))
                return np.zeros(int(size), dtype=int)

        fake = FakeGenerator()
        with mock.patch.object(f62.np.random, "default_rng", return_value=fake):
            payload = f62._bootstrap_child(low, control, full, variable_edges, variable_edges["SHMS_delta"], policy, "Left-lowe", 0, 0, {"mm_min": 1.10, "mm_max": 1.16})
        self.assertEqual(fake.draws, [(2, 2), (3, 3), (1, 1), (2, 2)])
        self.assertEqual(payload["policy"]["replicas"], 1)
        self.assertEqual(payload["policy"]["global_seed"], 20260916)
        reordered = self._build(f1_artifacts=list(reversed(self.artifacts)))
        self.assertEqual(self._build()["fingerprint"], reordered["fingerprint"])

    def test_aggregate_guard_and_artifact_fingerprint_are_deterministic(self):
        with self.assertRaisesRegex(f62.MethodAAcceptanceRefinementValidationError, "forbidden_event_persistence:entry_index"):
            f62._assert_aggregate_only_persistence({"entry_index": 1})
        first = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(self.artifacts, self.f3, self.f4, self.f5, self.f6, f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        second = f62.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(deepcopy(self.artifacts), deepcopy(self.f3), deepcopy(self.f4), deepcopy(self.f5), deepcopy(self.f6), f1_input_file_hashes=self.hashes, f3_input_file_sha256=f61_fixtures.F3_SHA, f4_input_file_sha256=f61_fixtures.F4_SHA, f5_input_file_sha256=f61_fixtures.F5_SHA, f6_1_input_file_sha256="e" * 64, input_paths={"f1": {}, "f3": "3", "f4": "4", "f5": "5", "f6_1": "6"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority, accepted_f6_1_artifact_authority_by_kinematic=self.f6_authority, bootstrap_test_config={"replicas": 3})
        self.assertEqual(first["artifact_fingerprint"], second["artifact_fingerprint"])
        for name, expected in {"non_authoritative": True, "validation_only": True, "manual_review_required": True, "accepted_f6_1_consumed": True, "accepted_f6_1_modified": False, "event_correction_evaluated_for_detached_validation": True, "event_correction_persisted": False, "production_application_performed": False, "production_objects_mutated": False, "yield_constructed": False, "cross_section_constructed": False, "root_object_constructed": False, "child_renormalization_performed": False, "smoothing_or_interpolation_performed": False, "absolute_probability_constructed": False, "method_b_numerical_dependency": False, "automatic_case_classification": False, "case_thresholds_defined": False, "final_yield_uncertainty_claimed": False}.items():
            self.assertEqual(first[name], expected)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "nonfinite.json"
            with self.assertRaises((ValueError, f62.MethodAAcceptanceRefinementValidationError)):
                f62.write_pion_hgcer_method_a_acceptance_refinement_validation_json(path, {"value": float("inf")})

    def test_module_has_no_production_or_method_b_import(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_method_a_acceptance_refinement_validation.py").read_text(encoding="utf-8")
        for forbidden in ("import ROOT", "rand_sub", "calculate_yield", "pion_component_subtraction", "pion_hgcer_refinement_method_b"):
            self.assertNotIn(forbidden, source)
        self.assertIn('"method_b_numerical_dependency": False', source)


if __name__ == "__main__":
    unittest.main()
