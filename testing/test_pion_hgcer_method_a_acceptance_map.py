"""Focused pure-Python tests for the detached Phase F.3 relative map."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_acceptance_map as acceptance_map
import pion_hgcer_method_a_acceptance_representation as representation
import test_pion_hgcer_method_a_acceptance_representation as f2_fixtures


def _accepted_f2(artifacts: list[dict[str, object]], hashes: dict[str, str]) -> dict[str, object]:
    """Make the frozen synthetic F.2 fixture represent the accepted hgcer3 decision."""
    artifact = representation.build_pion_hgcer_method_a_acceptance_representation_artifact(
        artifacts, input_file_hashes=hashes, input_paths={"f1": {}},
    )
    result = artifact["representation"]
    assert isinstance(result, dict)
    summaries = result["candidate_summaries"]
    assert isinstance(summaries, list)
    for summary in summaries:
        if isinstance(summary, dict) and summary["candidate_id"] in ("delta_only", "track3"):
            summary["information_gate_passed"] = False
            summary["overall_candidate_passed"] = False
    result["recommendation"] = {
        "recommendation_status": "unique_supported_reduced_basis",
        "recommended_basis": "hgcer3",
        "basis_frozen": False,
        "manual_review_required": True,
    }
    fingerprints = result["fingerprint_inputs"]
    assert isinstance(fingerprints, dict)
    fingerprints["candidate_summaries"] = summaries
    fingerprints["recommendation"] = result["recommendation"]
    result["fingerprint"] = representation._sha256(fingerprints)
    artifact["artifact_fingerprint"] = representation._sha256({
        "schema_version": artifact["schema_version"],
        "representation_fingerprint": result["fingerprint"],
        "input_paths": artifact["provenance"]["input_paths"],
    })
    return artifact


class MethodAAcceptanceMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = f2_fixtures.make_five_artifacts()
        cls.hashes = f2_fixtures.input_hashes()
        cls.f2 = _accepted_f2(cls.artifacts, cls.hashes)

    def _build(self, artifacts: list[dict[str, object]] | None = None, f2: dict[str, object] | None = None, **kwargs: object) -> dict[str, object]:
        return acceptance_map.build_pion_hgcer_method_a_acceptance_map(
            self.artifacts if artifacts is None else artifacts,
            self.f2 if f2 is None else f2,
            f1_input_file_hashes=self.hashes,
            f2_input_file_sha256="a" * 64,
            **kwargs,
        )

    def test_valid_map_has_exactly_fifteen_detached_models(self):
        result = self._build()
        self.assertTrue(result["available"])
        self.assertEqual(result["accepted_basis"], "hgcer3")
        self.assertEqual(result["ordered_features"], list(acceptance_map.ACCEPTED_FEATURES))
        self.assertEqual(len(result["models"]), 15)
        self.assertTrue(result["basis_frozen"])
        for forbidden in ("event_records", "event_probabilities", "weights", "correction", "normalization", "grid"):
            self.assertNotIn(forbidden, result)
        self.assertTrue(all(model["f2_support_continuity"]["passed"] for model in result["models"]))
        for model in result["models"]:
            self.assertNotIn("training", model)
            self.assertNotIn("application", model)
            self.assertNotIn("event_records", model)

    def test_reordered_f1_records_preserve_map_fingerprint(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            contract = artifact["contract"]
            contract["method_a_training_records"].reverse()
            contract["application_records"].reverse()
            f2_fixtures._seal_f1_artifact(artifact)
        reordered_f2 = _accepted_f2(reordered, self.hashes)
        changed = self._build(reordered, reordered_f2)
        self.assertEqual(changed["fingerprint"], self._build()["fingerprint"])
        self.assertEqual(changed["models"], self._build()["models"])

    def test_unsealed_f1_tampering_fails_authority(self):
        changed = deepcopy(self.artifacts)
        changed[0]["contract"]["method_a_training_records"][0]["SHMS_delta"] += 0.1
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "fingerprint"):
            self._build(changed)

    def test_resealed_wrong_f1_metadata_fails_authority(self):
        changed = deepcopy(self.artifacts)
        changed[0]["contract"]["feature_metadata"]["phi_is_training_feature"] = True
        f2_fixtures._seal_f1_artifact(changed[0])
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "feature_metadata"):
            self._build(changed)

    def test_stale_f2_representation_or_wrapper_fingerprint_fails(self):
        stale = deepcopy(self.f2)
        stale["representation"]["groups"][0]["candidate_metrics"]["hgcer3"]["application_support"]["application_ood_count"] += 1
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "fingerprint_content"):
            self._build(f2=stale)
        stale = deepcopy(self.f2)
        stale["artifact_fingerprint"] = "0" * 64
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "artifact_fingerprint"):
            self._build(f2=stale)

    def test_wrong_f2_basis_or_feature_order_fails(self):
        changed = deepcopy(self.f2)
        changed["representation"]["recommendation"]["recommended_basis"] = "track3"
        changed["representation"]["fingerprint_inputs"]["recommendation"] = changed["representation"]["recommendation"]
        changed["representation"]["fingerprint"] = representation._sha256(changed["representation"]["fingerprint_inputs"])
        changed["artifact_fingerprint"] = representation._sha256({"schema_version": changed["schema_version"], "representation_fingerprint": changed["representation"]["fingerprint"], "input_paths": changed["provenance"]["input_paths"]})
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "recommendation"):
            self._build(f2=changed)
        changed = deepcopy(self.f2)
        candidate = next(item for item in changed["representation"]["candidate_definitions"] if item["candidate_id"] == "hgcer3")
        candidate["ordered_features"].reverse()
        changed["representation"]["fingerprint_inputs"]["candidate_definitions"] = changed["representation"]["candidate_definitions"]
        changed["representation"]["algorithm_fingerprint"] = representation._sha256({"representation_schema_version": representation.METHOD_A_ACCEPTANCE_REPRESENTATION_SCHEMA_VERSION, "fingerprint_schema_version": representation.METHOD_A_ACCEPTANCE_REPRESENTATION_FINGERPRINT_SCHEMA_VERSION, "candidate_definitions": changed["representation"]["candidate_definitions"], "algorithm_config": changed["representation"]["algorithm_config"]})
        changed["representation"]["fingerprint_inputs"]["algorithm_fingerprint"] = changed["representation"]["algorithm_fingerprint"]
        changed["representation"]["fingerprint"] = representation._sha256(changed["representation"]["fingerprint_inputs"])
        changed["artifact_fingerprint"] = representation._sha256({"schema_version": changed["schema_version"], "representation_fingerprint": changed["representation"]["fingerprint"], "input_paths": changed["provenance"]["input_paths"]})
        with self.assertRaisesRegex(acceptance_map.MethodAAcceptanceMapError, "hgcer3_definition"):
            self._build(f2=changed)

    def test_f2_support_continuity_failure_makes_global_map_unavailable(self):
        changed = deepcopy(self.f2)
        support = changed["representation"]["groups"][0]["candidate_metrics"]["hgcer3"]["application_support"]
        support["training_nn_p50"] += 1.0e-6
        changed["representation"]["fingerprint_inputs"]["groups"] = changed["representation"]["groups"]
        changed["representation"]["fingerprint"] = representation._sha256(changed["representation"]["fingerprint_inputs"])
        changed["artifact_fingerprint"] = representation._sha256({"schema_version": changed["schema_version"], "representation_fingerprint": changed["representation"]["fingerprint"], "input_paths": changed["provenance"]["input_paths"]})
        result = self._build(f2=changed)
        self.assertFalse(result["available"])
        self.assertEqual(result["models"], [])
        self.assertIn("f2_application_support_continuity_failed", {item["invalid_reason"] for item in result["invalid_parents"]})

    def test_insufficient_support_scaling_and_optimizer_fail_closed(self):
        result = self._build(algorithm_config={"minimum_low_count": 26})
        self.assertFalse(result["available"])
        with patch.object(acceptance_map, "_robust_scale", return_value=None):
            result = self._build()
        self.assertFalse(result["available"])
        with patch.object(acceptance_map, "minimize", return_value=type("Result", (), {"success": False, "x": np.zeros(4), "fun": 0.0})()):
            result = self._build()
        self.assertFalse(result["available"])
        self.assertIn("optimizer_not_converged", {item["invalid_reason"] for item in result["invalid_parents"]})

    def test_relative_response_is_one_at_median_and_ood_grid_is_masked(self):
        result = self._build()
        model = result["models"][0]
        training = [row for artifact in acceptance_map._validate_f1_artifacts(self.artifacts) if artifact["setting_id"] == model["setting_id"] for row in artifact["training"] if row["t_index"] == model["canonical_t_index"]]
        median = np.asarray(model["scaler"]["median"], dtype=float)
        values = np.vstack((median, median + 1.0e6))
        response, supported = acceptance_map.evaluate_relative_response_grid(model, values, np.asarray([[row[name] for name in acceptance_map.ACCEPTED_FEATURES] for row in training], dtype=float))
        self.assertAlmostEqual(float(np.exp(np.dot(np.zeros(3), np.asarray(model["logistic_fit"]["coefficients"])))), 1.0, places=14)
        self.assertEqual(model["relative_response_at_training_median"], 1.0)
        self.assertFalse(bool(supported[1])); self.assertTrue(np.isnan(response[1]))

    def test_fingerprint_is_sensitive_to_config_and_artifact_wrapper_is_deterministic(self):
        baseline = self._build()
        changed = self._build(algorithm_config={"regularization_lambda": 2.0e-3})
        self.assertNotEqual(baseline["fingerprint"], changed["fingerprint"])
        alternate_hashes = {name: "{:064x}".format(index + 11) for index, name in enumerate(sorted(self.hashes))}
        alternate_f2 = _accepted_f2(self.artifacts, alternate_hashes)
        alternate = acceptance_map.build_pion_hgcer_method_a_acceptance_map(self.artifacts, alternate_f2, f1_input_file_hashes=alternate_hashes, f2_input_file_sha256="a" * 64)
        self.assertNotEqual(baseline["fingerprint"], alternate["fingerprint"])
        wrapper = acceptance_map.build_pion_hgcer_method_a_acceptance_map_artifact(self.artifacts, self.f2, f1_input_file_hashes=self.hashes, f2_input_file_sha256="a" * 64, input_paths={"f1": "one", "f2": "two"})
        self.assertTrue(wrapper["basis_frozen"])
        self.assertEqual(len(wrapper["artifact_fingerprint"]), 64)


if __name__ == "__main__":
    unittest.main()
