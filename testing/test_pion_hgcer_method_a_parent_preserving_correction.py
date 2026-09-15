"""Focused pure-Python tests for detached Phase F.4 correction construction."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_acceptance_map as f3_map
import pion_hgcer_method_a_parent_preserving_correction as f4
import test_pion_hgcer_method_a_acceptance_map as f3_tests
import test_pion_hgcer_method_a_acceptance_representation as f1_fixtures


def _artifacts() -> list[dict[str, object]]:
    artifacts = f1_fixtures.make_five_artifacts()
    for artifact in artifacts:
        for row in artifact["contract"]["application_records"]:
            row["phi_status"] = "inside_phi"
            row["baseline_pion_weight_w0"] = 1.0 + 0.001 * (row["entry_index"] % 100)
            row["signed_source_coefficient"] = 1.0 if row["source_label"] == "prompt" else -0.01
            row["signed_baseline_event_contribution"] = row["baseline_pion_weight_w0"] * row["signed_source_coefficient"]
        f1_fixtures._seal_f1_artifact(artifact)
    return artifacts


def _f3(artifacts: list[dict[str, object]], hashes: dict[str, str]) -> dict[str, object]:
    f2 = f3_tests._accepted_f2(artifacts, hashes)
    return f3_map.build_pion_hgcer_method_a_acceptance_map_artifact(artifacts, f2, f1_input_file_hashes=hashes, f2_input_file_sha256="a" * 64, input_paths={"f1": {}, "f2": "two"})


class ParentPreservingCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = _artifacts(); cls.hashes = f1_fixtures.input_hashes(); cls.f3 = _f3(cls.artifacts, cls.hashes)

    def _build(self, artifacts=None, f3=None, **kwargs):
        return f4.build_pion_hgcer_method_a_parent_preserving_correction(self.artifacts if artifacts is None else artifacts, self.f3 if f3 is None else f3, f1_input_file_hashes=self.hashes, f3_input_file_sha256="b" * 64, **kwargs)

    def test_valid_constructs_fifteen_parent_preserving_aggregate_only_corrections(self):
        result = self._build()
        self.assertTrue(result["available"]); self.assertEqual(len(result["parents"]), 15)
        self.assertTrue(result["parent_normalization_constructed"]); self.assertFalse(result["event_correction_persisted"])
        for parent in result["parents"]:
            self.assertTrue(parent["closure_passed"])
            self.assertAlmostEqual(parent["baseline_parent_sum"], parent["adjusted_parent_sum"], places=10)
            self.assertAlmostEqual(parent["ood_final_correction"], 1.0 / parent["parent_normalization"], places=14)
            self.assertNotIn("application_records", parent); self.assertNotIn("event_corrections", parent)
            self.assertAlmostEqual(sum(row["baseline_signed_sum"] for row in parent["canonical_phi_diagnostics"]), parent["baseline_parent_sum"], places=10)
            self.assertAlmostEqual(sum(row["adjusted_signed_sum"] for row in parent["canonical_phi_diagnostics"]), parent["adjusted_parent_sum"], places=10)

    def test_signed_baseline_and_phi_authority_fail_closed(self):
        changed = deepcopy(self.artifacts)
        changed[0]["contract"]["application_records"][0]["signed_baseline_event_contribution"] += 1.0
        f1_fixtures._seal_f1_artifact(changed[0])
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "signed_baseline"):
            self._build(changed)
        changed = deepcopy(self.artifacts); changed[0]["contract"]["application_records"][0]["phi_status"] = "outside_phi"; f1_fixtures._seal_f1_artifact(changed[0])
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "phi_status"):
            self._build(changed)

    def test_stale_f3_fingerprint_and_wrong_basis_fail(self):
        changed = deepcopy(self.f3); changed["acceptance_map"]["models"][0]["logistic_fit"]["coefficients"][0] += 1.0
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "fingerprint"):
            self._build(f3=changed)
        changed = deepcopy(self.f3); changed["acceptance_map"]["accepted_basis"] = "track3"
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "accepted_basis"):
            self._build(f3=changed)
        changed = deepcopy(self.f3); config = changed["acceptance_map"]["algorithm_config"]
        config["regularization_lambda"] = 2.0e-3
        map_value = changed["acceptance_map"]
        map_value["algorithm_fingerprint"] = f3_map._sha256({"map_schema_version": map_value["schema_version"], "fingerprint_schema_version": map_value["fingerprint_schema_version"], "accepted_basis": map_value["accepted_basis"], "ordered_features": map_value["ordered_features"], "algorithm_config": config})
        map_value["fingerprint_inputs"]["algorithm_config"] = config
        map_value["fingerprint_inputs"]["algorithm_fingerprint"] = map_value["algorithm_fingerprint"]
        map_value["fingerprint"] = f3_map._sha256(map_value["fingerprint_inputs"])
        changed["artifact_fingerprint"] = f3_map._sha256({"schema_version": changed["schema_version"], "map_fingerprint": map_value["fingerprint"], "input_paths": changed["provenance"]["input_paths"]})
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "f3_algorithm_config"):
            self._build(f3=changed)

    def test_reconstructed_f3_support_mismatch_makes_global_result_unavailable(self):
        changed = deepcopy(self.f3); model = changed["acceptance_map"]["models"][0]
        model["application_support"]["support_distance_threshold"] += 1.0e-5
        changed["acceptance_map"]["fingerprint_inputs"]["models"] = changed["acceptance_map"]["models"]
        changed["acceptance_map"]["fingerprint"] = f3_map._sha256(changed["acceptance_map"]["fingerprint_inputs"])
        changed["artifact_fingerprint"] = f3_map._sha256({"schema_version": changed["schema_version"], "map_fingerprint": changed["acceptance_map"]["fingerprint"], "input_paths": changed["provenance"]["input_paths"]})
        result = self._build(f3=changed)
        self.assertFalse(result["available"]); self.assertEqual(result["parents"], [])

    def test_nonpositive_signed_parent_sum_makes_global_result_unavailable(self):
        changed = deepcopy(self.artifacts)
        for artifact in changed:
            for row in artifact["contract"]["application_records"]:
                row["signed_source_coefficient"] = -1.0; row["signed_baseline_event_contribution"] = -row["baseline_pion_weight_w0"]
            f1_fixtures._seal_f1_artifact(artifact)
        changed_f3 = _f3(changed, self.hashes)
        result = self._build(changed, changed_f3)
        self.assertFalse(result["available"])
        self.assertIn("baseline_parent_sum_invalid", {row["invalid_reason"] for row in result["invalid_parents"]})

    def test_fingerprint_order_invariance_and_input_sensitivity(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            artifact["contract"]["application_records"].reverse(); artifact["contract"]["method_a_training_records"].reverse(); f1_fixtures._seal_f1_artifact(artifact)
        self.assertEqual(self._build()["fingerprint"], self._build(reordered, _f3(reordered, self.hashes))["fingerprint"])
        alternative_hashes = {name: "{:064x}".format(index + 31) for index, name in enumerate(sorted(self.hashes))}
        alternative = f4.build_pion_hgcer_method_a_parent_preserving_correction(self.artifacts, _f3(self.artifacts, alternative_hashes), f1_input_file_hashes=alternative_hashes, f3_input_file_sha256="b" * 64)
        self.assertNotEqual(self._build()["fingerprint"], alternative["fingerprint"])


if __name__ == "__main__": unittest.main()
