"""Focused tests for detached Phase F.5 signed template propagation."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5
import collect_pion_hgcer_validation_bundle as collector
import test_pion_hgcer_method_a_parent_preserving_correction as f4_fixtures
import test_pion_hgcer_method_a_acceptance_representation as f1_fixtures


KINEMATIC = "Q4p4W2p74"
F3_SHA = "b" * 64
F4_SHA = "c" * 64
T_EDGES = [0.4, 0.5666666666666667, 0.7333333333333334, 0.9]
PHI_EDGES = [-180.0, -140.0, -100.0, -60.0, -20.0, 20.0, 60.0, 100.0, 140.0, 180.0]


def _artifacts(*, nine_phi_children: bool = True) -> list[dict[str, object]]:
    artifacts = f4_fixtures._artifacts()
    for artifact in artifacts:
        contract = artifact["contract"]
        contract["t_edges"] = list(T_EDGES); contract["phi_edges"] = list(PHI_EDGES)
        for row in contract["method_a_training_records"] + contract["application_records"]:
            index = int(row["t_index"]); row["t_low"] = T_EDGES[index]; row["t_high"] = T_EDGES[index + 1]
        for row in contract["application_records"]:
            phi_index = int(row["entry_index"]) % 9 if nine_phi_children else 0
            row["phi_index"] = phi_index; row["phi_low"] = PHI_EDGES[phi_index]; row["phi_high"] = PHI_EDGES[phi_index + 1]
            row["phi_degrees"] = (PHI_EDGES[phi_index] + PHI_EDGES[phi_index + 1]) / 2.0; row["phi_status"] = "inside_phi"
        f1_fixtures._seal_f1_artifact(artifact)
    return artifacts


def _f4_artifact(artifacts: list[dict[str, object]]) -> tuple[dict[str, str], dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    hashes = f1_fixtures.input_hashes(); f3 = f4_fixtures._f3(artifacts, hashes); f3_authority = f4_fixtures._authority(f3, F3_SHA)
    artifact = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, f3, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, input_paths={"f1": {}, "f3": "three"}, generated_at_utc="test", git_head="1" * 40, git_status_short="", accepted_f3_runtime_authority_by_kinematic=f3_authority)
    correction = artifact["correction"]
    assert isinstance(correction, dict)
    f4_authority = {KINEMATIC: {"source_file_sha256": F4_SHA, "correction_fingerprint": correction["fingerprint"], "artifact_fingerprint": artifact["artifact_fingerprint"], "farm_source_head": "2" * 40}}
    return hashes, f3, artifact, f3_authority, f4_authority


class TphiPropagationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = _artifacts(); cls.hashes, cls.f3, cls.f4, cls.f3_authority, cls.f4_authority = _f4_artifact(cls.artifacts)

    def _build(self, artifacts=None, f3=None, f4_artifact=None, *, f4_sha=F4_SHA, f4_authority=None):
        return f5.build_pion_hgcer_method_a_tphi_propagation(self.artifacts if artifacts is None else artifacts, self.f3 if f3 is None else f3, self.f4 if f4_artifact is None else f4_artifact, f1_input_file_hashes=self.hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=f4_sha, accepted_f4_runtime_authority_by_kinematic=self.f4_authority if f4_authority is None else f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority)

    def test_builds_one_hundred_thirty_five_explicit_signed_cells_without_event_persistence(self):
        result = self._build()
        self.assertTrue(result["available"]); self.assertEqual(len(result["setting_templates"]), 5); self.assertEqual(len(result["parent_closure"]), 15)
        serialized = json.dumps(result, sort_keys=True)
        for forbidden in ("application_records", "event_corrections", "correction_factors", "raw_shape_factors", "in_support_mask", "event_index", "grid"):
            self.assertNotIn(forbidden, serialized)
        for setting in result["setting_templates"]:
            self.assertEqual(setting["t_edges"], T_EDGES); self.assertEqual(setting["phi_edges"], PHI_EDGES)
            self.assertEqual(sum(sum(row) for row in setting["event_counts"]), sum(parent["application_event_count"] for parent in self.f4["correction"]["parents"] if parent["setting_id"] == setting["setting_id"]))
            for name in ("event_counts", "baseline_signed_contents", "adjusted_signed_contents", "signed_delta_contents", "baseline_share_of_parent", "adjusted_share_of_parent", "redistribution_fraction_of_parent"):
                self.assertEqual((len(setting[name]), len(setting[name][0])), (3, 9))
            for parent in setting["parents"]:
                self.assertAlmostEqual(parent["closure_residual"], 0.0, places=10)

    def test_empty_children_are_explicit_zeroes_but_do_not_require_f4_empty_diagnostics(self):
        artifacts = _artifacts(nine_phi_children=False); hashes, f3, f4_artifact, f3_authority, f4_authority = _f4_artifact(artifacts)
        result = f5.build_pion_hgcer_method_a_tphi_propagation(artifacts, f3, f4_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        for setting in result["setting_templates"]:
            for t_index in range(3):
                self.assertEqual(setting["event_counts"][t_index][1:], [0] * 8)
                self.assertEqual(setting["baseline_signed_contents"][t_index][1:], [0.0] * 8)
                self.assertEqual(setting["adjusted_signed_contents"][t_index][1:], [0.0] * 8)
                self.assertEqual(setting["signed_delta_contents"][t_index][1:], [0.0] * 8)

    def test_reordering_f1_records_preserves_f5_fingerprint(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            artifact["contract"]["method_a_training_records"].reverse(); artifact["contract"]["application_records"].reverse(); f1_fixtures._seal_f1_artifact(artifact)
        hashes, f3, f4_artifact, f3_authority, f4_authority = _f4_artifact(reordered)
        changed = f5.build_pion_hgcer_method_a_tphi_propagation(reordered, f3, f4_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        self.assertEqual(self._build()["fingerprint"], changed["fingerprint"])

    def test_f4_authority_and_reproduction_tampering_fail_closed(self):
        with self.assertRaisesRegex(f5.MethodATPhiPropagationError, "source_file_sha256_mismatch"):
            self._build(f4_sha="d" * 64)
        tampered = deepcopy(self.f4); tampered["correction"]["parents"][0]["baseline_parent_sum"] += 1.0
        with self.assertRaisesRegex(f5.MethodATPhiPropagationError, "fingerprint|reproduction"):
            self._build(f4_artifact=tampered)

    def test_wrong_geometry_and_occupied_f4_child_diagnostic_fail(self):
        wrong = deepcopy(self.artifacts); wrong[0]["contract"]["phi_edges"][1] = -139.0; f1_fixtures._seal_f1_artifact(wrong[0])
        with self.assertRaisesRegex(f5.MethodATPhiPropagationError, "geometry"):
            self._build(artifacts=wrong)
        corrupted = deepcopy(self.f4); corrupted["correction"]["parents"][0]["canonical_phi_diagnostics"].pop()
        correction = corrupted["correction"]; correction["fingerprint"] = f4._sha256(correction["fingerprint_inputs"])
        corrupted["artifact_fingerprint"] = f4._sha256({"schema_version": corrupted["schema_version"], "correction_fingerprint": correction["fingerprint"], "input_paths": corrupted["provenance"]["input_paths"]})
        authority = deepcopy(self.f4_authority); authority[KINEMATIC]["correction_fingerprint"] = correction["fingerprint"]; authority[KINEMATIC]["artifact_fingerprint"] = corrupted["artifact_fingerprint"]
        with self.assertRaisesRegex(f5.MethodATPhiPropagationError, "reproduction|phi"):
            self._build(f4_artifact=corrupted, f4_authority=authority)

    def test_artifact_wrapper_and_fingerprint_are_deterministic_and_sensitive(self):
        artifact = f5.build_pion_hgcer_method_a_tphi_propagation_artifact(self.artifacts, self.f3, self.f4, f1_input_file_hashes=self.hashes, f3_input_file_sha256=F3_SHA, f4_input_file_sha256=F4_SHA, input_paths={"f1": {}, "f3": "three", "f4": "four"}, accepted_f4_runtime_authority_by_kinematic=self.f4_authority, accepted_f3_runtime_authority_by_kinematic=self.f3_authority)
        self.assertEqual(artifact["artifact_fingerprint"], f5._sha256({"schema_version": artifact["schema_version"], "propagation_fingerprint": artifact["propagation"]["fingerprint"], "input_paths": artifact["provenance"]["input_paths"]}))
        self.assertNotEqual(self._build()["fingerprint"], self._build(f4_authority={KINEMATIC: {**self.f4_authority[KINEMATIC], "farm_source_head": "3" * 40}})["fingerprint"])

    def test_f5_profile_is_generic_and_declares_only_f5_f4_f3_and_f1_inputs(self):
        profile = collector.load_validation_profile(REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_f5.json")
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], "67e0298c51759c7a5ba693464d2c2655bf39250d")
        self.assertEqual({item["key"] for item in profile["artifacts"]["global"]}, {"f5_tphi_propagation_json", "f5_tphi_propagation_pdf", "f4_parent_preserving_correction_json", "f3_acceptance_map_json"})
        self.assertEqual(len(profile["artifacts"]["settings"]), 1)


if __name__ == "__main__":
    unittest.main()
