"""Focused pure-Python coverage for detached Phase F.6.1 validation."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_reweighting_validation as f6
import pion_hgcer_method_a_tphi_propagation as f5
import test_pion_hgcer_method_a_acceptance_representation as f1_fixtures
import test_pion_hgcer_method_a_parent_preserving_correction as f4_fixtures
import test_pion_hgcer_method_a_tphi_propagation as f5_fixtures


KINEMATIC = "Q4p4W2p74"
F3_SHA = "b" * 64
F4_SHA = "c" * 64
F5_SHA = "d" * 64


def _enriched_artifacts() -> list[dict[str, object]]:
    artifacts = f5_fixtures._artifacts()
    for artifact in artifacts:
        contract = artifact["contract"]
        training = contract["method_a_training_records"]
        application = contract["application_records"]
        controls: dict[tuple[str, int], dict[str, object]] = {}
        for row in training:
            index, t_index = int(row["entry_index"]), int(row["t_index"])
            phi_edges = contract["phi_edges"]; phi_index = index % 9
            row.update({"allcuts": True, "analysis_t": 0.45 + 0.1 * t_index, "analysis_MM": -0.12 + 0.001 * (index % 30), "Q2": 4.4, "W": 2.74, "epsilon": 0.35, "phi": (phi_edges[phi_index] + phi_edges[phi_index + 1]) / 2.0})
            if row["response_class"] == "control": controls[(str(row["source_label"]), index)] = row
        for row in application:
            index, t_index = int(row["entry_index"]), int(row["t_index"])
            row.update({"allcuts": True, "nommcuts": True, "analysis_t": 0.45 + 0.1 * t_index, "analysis_MM": -0.12 + 0.001 * (index % 30), "Q2": 4.4, "W": 2.74, "epsilon": 0.35})
            matched = controls.get((str(row["source_label"]), index))
            if matched is not None:
                for name in ("SHMS_delta", "P_hgcer_npeSum", "P_hgcer_xAtCer", "P_hgcer_yAtCer", "SHMS_xptar", "SHMS_yptar", "analysis_t", "analysis_MM", "Q2", "W", "epsilon", "allcuts", "nommcuts"):
                    row[name] = matched[name]
                row["phi_degrees"] = matched["phi"]
            else:
                # Keep a signed nonprompt-only MM region so the focused
                # fixture confirms that negative aggregate bins survive.
                row["analysis_MM"] = 0.5 + 0.001 * (index % 20)
        f1_fixtures._seal_f1_artifact(artifact)
    return artifacts


def _chain(artifacts: list[dict[str, object]] | None = None) -> tuple[list[dict[str, object]], dict[str, str], dict[str, object], dict[str, object], dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    artifacts = _enriched_artifacts() if artifacts is None else artifacts
    hashes = f1_fixtures.input_hashes(); f3 = f4_fixtures._f3(artifacts, hashes)
    f3_authority = f4_fixtures._authority(f3, F3_SHA)
    f4_artifact = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, f3, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, input_paths={"f1": {}, "f3": "three"}, accepted_f3_runtime_authority_by_kinematic=f3_authority)
    correction = f4_artifact["correction"]
    assert isinstance(correction, dict)
    f4_authority = {KINEMATIC: {"source_file_sha256": F4_SHA, "correction_fingerprint": correction["fingerprint"], "artifact_fingerprint": f4_artifact["artifact_fingerprint"], "farm_source_head": "2" * 40}}
    f5_artifact = f5.build_pion_hgcer_method_a_tphi_propagation_artifact(artifacts, f3, f4_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, input_paths={"f1": {}, "f3": "three", "f4": "four"}, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
    propagation = f5_artifact["propagation"]
    assert isinstance(propagation, dict)
    authority = {KINEMATIC: {"f1_source_file_sha256": hashes, "f3_source_file_sha256": F3_SHA, "f3_map_fingerprint": f3["acceptance_map"]["fingerprint"], "f3_algorithm_fingerprint": f3["acceptance_map"]["algorithm_fingerprint"], "f3_artifact_fingerprint": f3["artifact_fingerprint"], "f4_source_file_sha256": F4_SHA, "f4_correction_fingerprint": correction["fingerprint"], "f4_artifact_fingerprint": f4_artifact["artifact_fingerprint"], "f5_source_file_sha256": F5_SHA, "f5_propagation_fingerprint": propagation["fingerprint"], "f5_artifact_fingerprint": f5_artifact["artifact_fingerprint"], "farm_source_head": "3" * 40}}
    return artifacts, hashes, f3, f4_artifact, f5_artifact, authority, f3_authority, f4_authority


class ReweightingValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts, cls.hashes, cls.f3, cls.f4, cls.f5, cls.authority, cls.f3_authority, cls.f4_authority = _chain()

    def _build(self, artifacts=None, f3=None, f4_artifact=None, f5_artifact=None, authority=None, *, f5_sha=F5_SHA):
        return f6.build_pion_hgcer_method_a_reweighting_validation(artifacts or self.artifacts, f3 or self.f3, f4_artifact or self.f4, f5_artifact or self.f5, f1_input_file_hashes=self.hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, f5_input_file_sha256=f5_sha, accepted_runtime_authority_by_kinematic=self.authority if authority is None else authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority)

    def test_valid_five_setting_chain_is_aggregate_only_and_exact(self):
        result = self._build()
        self.assertTrue(result["available"]); self.assertEqual(len(result["parents"]), 15)
        self.assertTrue(result["f4_reproduction"]["exact_payload_match"]); self.assertTrue(result["f5_reproduction"]["exact_payload_match"])
        self.assertEqual(result["f5_continuity"]["cell_count"], 135)
        for parent in result["parents"]:
            self.assertTrue(parent["identity_and_phi_parity"]["phi_semantics_closed_without_conversion"])
            self.assertEqual(parent["population_counts"]["prompt_physical_control"], 100)
            for comparison in parent["prompt_shape_comparisons"].values():
                self.assertAlmostEqual(sum(comparison["low_response_unit_area"]), 1.0)
                self.assertAlmostEqual(sum(comparison["baseline_prompt_control_unit_area"]), 1.0)
                self.assertAlmostEqual(sum(comparison["method_a_prompt_control_unit_area"]), 1.0)
                self.assertTrue(np.isfinite(comparison["delta_hellinger"]))
            for matrix in parent["hgcer_xy"].values():
                if isinstance(matrix, list) and matrix and isinstance(matrix[0], list): self.assertAlmostEqual(sum(sum(row) for row in matrix), 1.0)
        serialized = json.dumps(result, sort_keys=True)
        for forbidden in ("entry_index", "correction_factors", "raw_shape_factors", "in_support_mask", "application_records"):
            self.assertNotIn(forbidden, serialized)
        self.assertTrue(any(value < 0.0 for parent in result["parents"] for value in parent["signed_background"]["analysis_MM"]["baseline_signed_contents"]))

    def test_authority_factor_identity_phi_and_f5_fail_closed(self):
        changed_authority = deepcopy(self.authority); changed_authority[KINEMATIC]["f5_source_file_sha256"] = "e" * 64
        with self.assertRaisesRegex(f6.MethodAReweightingValidationError, "f5_source_file_sha256"):
            self._build(authority=changed_authority)
        changed = deepcopy(self.artifacts); changed[0]["contract"]["method_a_training_records"][100]["phi"] += 1.0; f1_fixtures._seal_f1_artifact(changed[0])
        artifacts, hashes, f3, f4_artifact, f5_artifact, authority, f3_authority, f4_authority = _chain(changed)
        with self.assertRaisesRegex(f6.MethodAReweightingValidationError, "phi_semantic"):
            f6.build_pion_hgcer_method_a_reweighting_validation(artifacts, f3, f4_artifact, f5_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, f5_input_file_sha256=F5_SHA, accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        original = f6._f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        def shortened(*args, **kwargs):
            correction, review = original(*args, **kwargs); altered = deepcopy(review); altered[0]["correction_factors"] = altered[0]["correction_factors"][:-1]; return correction, altered
        with mock.patch.object(f6._f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=shortened), self.assertRaisesRegex(f6.MethodAReweightingValidationError, "factor_alignment"):
            self._build()
        corrupt = deepcopy(self.f5); corrupt["propagation"]["setting_templates"][0]["event_counts"][0][0] += 1
        with self.assertRaisesRegex(f6.MethodAReweightingValidationError, "fingerprint|reproduction"):
            self._build(f5_artifact=corrupt)

    def test_artifact_fingerprint_is_deterministic_and_sensitive_to_authority(self):
        first = f6.build_pion_hgcer_method_a_reweighting_validation_artifact(self.artifacts, self.f3, self.f4, self.f5, f1_input_file_hashes=self.hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, f5_input_file_sha256=F5_SHA, input_paths={"f1": {}, "f3": "three", "f4": "four", "f5": "five"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority)
        second = f6.build_pion_hgcer_method_a_reweighting_validation_artifact(deepcopy(self.artifacts), deepcopy(self.f3), deepcopy(self.f4), deepcopy(self.f5), f1_input_file_hashes=self.hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, f5_input_file_sha256=F5_SHA, input_paths={"f1": {}, "f3": "three", "f4": "four", "f5": "five"}, accepted_runtime_authority_by_kinematic=self.authority, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority)
        self.assertEqual(first["artifact_fingerprint"], second["artifact_fingerprint"])
        self.assertEqual(first["validation"]["fingerprint"], second["validation"]["fingerprint"])

    def test_paired_factor_association_survives_outer_record_reordering(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            artifact["contract"]["method_a_training_records"].reverse()
            artifact["contract"]["application_records"].reverse()
            f1_fixtures._seal_f1_artifact(artifact)
        artifacts, hashes, f3, f4_artifact, f5_artifact, authority, f3_authority, f4_authority = _chain(reordered)
        result = f6.build_pion_hgcer_method_a_reweighting_validation(artifacts, f3, f4_artifact, f5_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, f5_input_file_sha256=F5_SHA, accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        self.assertTrue(result["f5_continuity"]["checked"])
        self.assertEqual(result["f5_continuity"]["cell_count"], 135)

    def test_module_has_no_production_or_method_b_dependency(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_method_a_reweighting_validation.py").read_text(encoding="utf-8")
        for forbidden in ("import ROOT", "rand_sub", "calculate_yield", "pion_component_subtraction", "pion_hgcer_refinement_method_b"):
            self.assertNotIn(forbidden, source)
        self.assertIn('"method_b_numerical_dependency": False', source)


if __name__ == "__main__":
    unittest.main()
