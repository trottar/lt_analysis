"""Focused authority, formula, and aggregate-persistence tests for F.4.Fix.1."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_acceptance_map as f3_map
import pion_hgcer_method_a_parent_preserving_correction as f4
import test_pion_hgcer_method_a_acceptance_map as f3_tests
import test_pion_hgcer_method_a_acceptance_representation as f1_fixtures


KINEMATIC = "Q4p4W2p74"
F3_FILE_SHA = "b" * 64


def _artifacts() -> list[dict[str, object]]:
    artifacts = f1_fixtures.make_five_artifacts()
    for artifact in artifacts:
        for row in artifact["contract"]["application_records"]:
            row["phi_status"] = "inside_phi"
            row["phi_degrees"] = 0.0
            row["baseline_pion_weight_w0"] = 1.0 + 0.001 * (row["entry_index"] % 100)
            row["signed_source_coefficient"] = 1.0 if row["source_label"] == "prompt" else -0.01
            row["signed_baseline_event_contribution"] = row["baseline_pion_weight_w0"] * row["signed_source_coefficient"]
        f1_fixtures._seal_f1_artifact(artifact)
    return artifacts


def _f3(artifacts: list[dict[str, object]], hashes: dict[str, str]) -> dict[str, object]:
    f2 = f3_tests._accepted_f2(artifacts, hashes)
    return f3_map.build_pion_hgcer_method_a_acceptance_map_artifact(artifacts, f2, f1_input_file_hashes=hashes, f2_input_file_sha256="a" * 64, input_paths={"f1": {}, "f2": "two"})


def _authority(f3: dict[str, object], source_file_sha256: str = F3_FILE_SHA) -> dict[str, object]:
    map_value = f3["acceptance_map"]
    return {KINEMATIC: {"source_file_sha256": source_file_sha256, "map_fingerprint": map_value["fingerprint"], "algorithm_fingerprint": map_value["algorithm_fingerprint"], "artifact_fingerprint": f3["artifact_fingerprint"], "farm_source_head": "1" * 40}}


def _reseal_f3(f3: dict[str, object]) -> None:
    map_value = f3["acceptance_map"]
    map_value["fingerprint_inputs"]["models"] = map_value["models"]
    map_value["fingerprint"] = f3_map._sha256(map_value["fingerprint_inputs"])
    f3["artifact_fingerprint"] = f3_map._sha256({"schema_version": f3["schema_version"], "map_fingerprint": map_value["fingerprint"], "input_paths": f3["provenance"]["input_paths"]})


class ParentPreservingCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = _artifacts()
        cls.hashes = f1_fixtures.input_hashes()
        cls.f3 = _f3(cls.artifacts, cls.hashes)
        cls.authority = _authority(cls.f3)

    def _build(self, artifacts=None, f3=None, *, authority=None, f3_sha=F3_FILE_SHA, review=False):
        method = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data if review else f4.build_pion_hgcer_method_a_parent_preserving_correction
        return method(self.artifacts if artifacts is None else artifacts, self.f3 if f3 is None else f3, f1_input_file_hashes=self.hashes, f3_input_file_sha256=f3_sha, accepted_f3_runtime_authority_by_kinematic=self.authority if authority is None else authority)

    def _changed_inputs(self, mutate) -> tuple[list[dict[str, object]], dict[str, object], dict[str, object]]:
        changed = deepcopy(self.artifacts)
        mutate(changed)
        for artifact in changed:
            f1_fixtures._seal_f1_artifact(artifact)
        changed_f3 = _f3(changed, self.hashes)
        return changed, changed_f3, _authority(changed_f3)

    def test_valid_accepted_authority_constructs_fifteen_aggregate_only_parents(self):
        result, review = self._build(review=True)
        self.assertTrue(result["available"])
        self.assertEqual(len(result["parents"]), 15)
        self.assertEqual(len(review), 15)
        authority = result["f3_runtime_authority"]
        self.assertTrue(authority["accepted_authority_match"])
        for key in ("source_file_sha256", "map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"):
            self.assertEqual(authority["accepted"][key], authority["observed"][key])
        serialized = json.dumps(result, sort_keys=True)
        for forbidden in ("application_records", "event_corrections", "correction_factors", "raw_shape_factors", "in_support_mask"):
            self.assertNotIn(forbidden, serialized)
        for parent in result["parents"]:
            self.assertTrue(parent["closure_passed"])
            self.assertAlmostEqual(parent["baseline_parent_sum"], parent["adjusted_parent_sum"], places=10)
            self.assertAlmostEqual(sum(row["baseline_signed_sum"] for row in parent["canonical_phi_diagnostics"]), parent["baseline_parent_sum"], places=10)
            self.assertAlmostEqual(sum(row["adjusted_signed_sum"] for row in parent["canonical_phi_diagnostics"]), parent["adjusted_parent_sum"], places=10)
            self.assertAlmostEqual(sum(row["signed_delta"] for row in parent["canonical_phi_diagnostics"]), 0.0, places=10)

    def test_accepted_authority_rejects_each_observed_identity_and_unknown_kinematic(self):
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "f3_runtime_authority_source_file_sha256_mismatch"):
            self._build(f3_sha="c" * 64)
        for key in ("map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"):
            authority = deepcopy(self.authority)
            authority[KINEMATIC][key] = "c" * 64
            with self.subTest(key=key), self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "f3_runtime_authority_{}_mismatch".format(key)):
                self._build(authority=authority)
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "kinematic_unsupported"):
            self._build(authority={})

    def test_resealed_alternate_f3_is_internally_valid_but_rejected_by_accepted_authority(self):
        alternate = deepcopy(self.f3)
        alternate["acceptance_map"]["models"][0]["logistic_fit"]["coefficients"][0] += 0.25
        _reseal_f3(alternate)
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "f3_runtime_authority_map_fingerprint_mismatch"):
            self._build(f3=alternate)
        accepted_alternate = self._build(f3=alternate, authority=_authority(alternate))
        self.assertTrue(accepted_alternate["available"])
        self.assertNotEqual(self._build()["fingerprint"], accepted_alternate["fingerprint"])

    def test_internal_f3_tampering_wrong_basis_and_inventory_still_fail(self):
        stale = deepcopy(self.f3)
        stale["acceptance_map"]["models"][0]["logistic_fit"]["coefficients"][0] += 1.0
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "fingerprint"):
            self._build(f3=stale)
        wrong_basis = deepcopy(self.f3)
        wrong_basis["acceptance_map"]["accepted_basis"] = "track3"
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "accepted_basis"):
            self._build(f3=wrong_basis)
        missing_model = deepcopy(self.f3)
        missing_model["acceptance_map"]["models"].pop()
        missing_model["acceptance_map"]["fingerprint_inputs"]["models"] = missing_model["acceptance_map"]["models"]
        _reseal_f3(missing_model)
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "f3_model_count"):
            self._build(f3=missing_model, authority=_authority(missing_model))

    def test_canonical_phi_semantics_fail_after_resealed_f1_edits(self):
        cases = {
            "phi_index": (lambda row: row.__setitem__("phi_index", 1), "phi_geometry"),
            "phi_low": (lambda row: row.__setitem__("phi_low", -179.0), "phi_geometry"),
            "phi_high": (lambda row: row.__setitem__("phi_high", 179.0), "phi_geometry"),
            "phi_degrees": (lambda row: row.__setitem__("phi_degrees", -181.0), "phi_assignment"),
            "outside_phi": (lambda row: row.__setitem__("phi_status", "outside_phi"), "phi_status"),
        }
        for name, (mutate_row, expected) in cases.items():
            def mutate(artifacts, apply=mutate_row):
                apply(artifacts[0]["contract"]["application_records"][0])
            artifacts, alternate_f3, authority = self._changed_inputs(mutate)
            with self.subTest(name=name), self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, expected):
                self._build(artifacts, alternate_f3, authority=authority)
        nonfinite = deepcopy(self.artifacts)
        nonfinite[0]["contract"]["application_records"][0]["phi_degrees"] = float("nan")
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "phi_degrees"):
            self._build(nonfinite)

    def test_ood_formula_and_large_finite_in_support_factor_are_not_clipped(self):
        def make_ood(artifacts):
            row = next(row for row in artifacts[0]["contract"]["application_records"] if row["source_label"] == "rand" and row["t_index"] == 0)
            row["SHMS_delta"] = 1.0e6
            row["P_hgcer_xAtCer"] = 1.0e6
            row["P_hgcer_yAtCer"] = -1.0e6
        artifacts, alternate_f3, authority = self._changed_inputs(make_ood)
        result, review = self._build(artifacts, alternate_f3, authority=authority, review=True)
        parent = next(row for row in result["parents"] if row["setting_id"] == "Left-lowe" and row["canonical_t_index"] == 0)
        values = next(row for row in review if row["setting_id"] == "Left-lowe" and row["canonical_t_index"] == 0)
        ood = ~values["in_support_mask"]
        self.assertGreater(np.count_nonzero(ood), 0)
        np.testing.assert_array_equal(values["raw_shape_factors"][ood], np.ones(np.count_nonzero(ood)))
        np.testing.assert_allclose(values["correction_factors"][ood], np.full(np.count_nonzero(ood), 1.0 / parent["parent_normalization"]), rtol=0.0, atol=1.0e-14)
        large = deepcopy(self.f3)
        large["acceptance_map"]["models"][0]["logistic_fit"]["coefficients"] = [8.0, 0.0, 0.0]
        _reseal_f3(large)
        _, large_review = self._build(f3=large, authority=_authority(large), review=True)
        self.assertGreater(float(np.max(large_review[0]["raw_shape_factors"])), 10.0)
        self.assertTrue(np.all(np.isfinite(large_review[0]["correction_factors"])))

    def test_zero_weight_signed_sources_and_multiple_children_remain_parent_only(self):
        def mutate(artifacts):
            contract = artifacts[0]["contract"]
            contract["phi_edges"] = [-180.0, 0.0, 180.0]
            for index, row in enumerate(contract["application_records"]):
                row["phi_index"] = index % 2
                row["phi_low"], row["phi_high"], row["phi_degrees"] = (-180.0, 0.0, -90.0) if index % 2 == 0 else (0.0, 180.0, 90.0)
            row = contract["application_records"][0]
            row["baseline_pion_weight_w0"] = 0.0
            row["signed_baseline_event_contribution"] = 0.0
        artifacts, alternate_f3, authority = self._changed_inputs(mutate)
        result = self._build(artifacts, alternate_f3, authority=authority)
        parent = next(row for row in result["parents"] if row["setting_id"] == "Left-lowe" and row["canonical_t_index"] == 0)
        self.assertEqual(parent["zero_baseline_weight_count"], 1)
        self.assertEqual(len(parent["canonical_phi_diagnostics"]), 2)
        self.assertGreaterEqual(len(parent["source_diagnostics"]), 2)
        self.assertNotIn("source_normalization", parent)
        self.assertNotIn("child_normalization", parent)

    def test_row_order_invariance_and_source_hash_sensitivity(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            artifact["contract"]["application_records"].reverse()
            artifact["contract"]["method_a_training_records"].reverse()
            f1_fixtures._seal_f1_artifact(artifact)
        reordered_f3 = _f3(reordered, self.hashes)
        self.assertEqual(self._build()["fingerprint"], self._build(reordered, reordered_f3, authority=_authority(reordered_f3))["fingerprint"])
        alternate_hashes = {name: "{:064x}".format(index + 31) for index, name in enumerate(sorted(self.hashes))}
        alternate_f3 = _f3(self.artifacts, alternate_hashes)
        alternate = f4.build_pion_hgcer_method_a_parent_preserving_correction(self.artifacts, alternate_f3, f1_input_file_hashes=alternate_hashes, f3_input_file_sha256=F3_FILE_SHA, accepted_f3_runtime_authority_by_kinematic=_authority(alternate_f3))
        self.assertNotEqual(self._build()["fingerprint"], alternate["fingerprint"])


if __name__ == "__main__":
    unittest.main()
