"""Focused pure-Python tests for the detached F.2 representation audit."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import pion_hgcer_method_a_acceptance_representation as representation


SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)
FEATURES = (
    "SHMS_delta", "SHMS_xptar", "SHMS_yptar", "P_hgcer_xAtCer", "P_hgcer_yAtCer",
)


def _features(label: str, index: int, t_index: int) -> dict[str, float]:
    signal = -1.0 if label == "low" else 1.0
    spread = 0.01 * (index % 25) + 0.001 * t_index
    return {
        "SHMS_delta": signal + spread,
        "SHMS_xptar": 0.1 * signal + spread,
        "SHMS_yptar": -0.1 * signal + 0.5 * spread,
        "P_hgcer_xAtCer": 0.2 * signal + 0.3 * spread,
        "P_hgcer_yAtCer": -0.2 * signal + 0.2 * spread,
    }


def _training_record(t_index: int, entry_index: int, response_class: str) -> dict[str, object]:
    npe = 1.0 if response_class == "low" else 3.0
    return {
        "source_label": "prompt", "entry_index": entry_index,
        "response_class": response_class, "nommcuts": True,
        "P_hgcer_npeSum": npe, "t_index": t_index,
        "t_low": float(t_index), "t_high": float(t_index + 1),
        **_features(response_class, entry_index, t_index),
    }


def _application_record(t_index: int, entry_index: int, source: str, index: int) -> dict[str, object]:
    return {
        "source_label": source, "entry_index": entry_index,
        "P_hgcer_npeSum": 3.0, "t_index": t_index,
        "t_low": float(t_index), "t_high": float(t_index + 1),
        "phi_index": 0, "phi_low": -180.0, "phi_high": 180.0,
        "phi_status": "in_range",
        **_features("control", index, t_index),
    }


def make_f1_artifact(phi: str, epsilon: str) -> dict[str, object]:
    training, application = [], []
    for t_index in range(3):
        base = t_index * 1000
        training.extend(_training_record(t_index, base + index, "low") for index in range(25))
        training.extend(_training_record(t_index, base + 100 + index, "control") for index in range(100))
        application.extend(
            _application_record(t_index, base + 100 + index, "prompt", index)
            for index in range(100)
        )
        application.extend(
            _application_record(t_index, base + 10000 + index, "rand", index)
            for index in range(20)
        )
    metadata = {
        "primary_acceptance_features": list(FEATURES),
        "training_population": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
        "training_low_definition": "0_lt_P_hgcer_npeSum_le_2",
        "training_control_definition": "P_hgcer_npeSum_gt_2",
        "application_population": "authoritative_physical_pion_control_P_hgcer_npeSum_gt_2",
        "parent_coordinate": "canonical_t",
        "downstream_yield_coordinates": ["canonical_t", "canonical_phi"],
        "phi_is_training_feature": False, "probability_map_constructed": False,
        "absolute_leakage_probability_claimed": False, "weight_adjustment_constructed": False,
        "method_b_numerical_dependency": False,
        "future_normalization_policy": "future_parent_t_only_no_tphi_child_renormalization",
    }
    contract = {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract/v2",
        "fingerprint_schema_version": "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2",
        "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete",
        "non_authoritative": True, "method_b_numerical_dependency": False,
        "event_application_performed": False, "production_application_performed": False,
        "production_objects_mutated": False, "future_weight_adjustment_constructed": False,
        "refinement_applied": False, "feature_metadata": metadata,
        "method_a_training_records": training, "application_records": application,
        "method_a_training_summary": {"by_t_delta": []},
        "phase_a_contract_fingerprint": "phase-a-contract-{}-{}".format(phi, epsilon),
        "phase_a_pion_event_population_fingerprint": "phase-a-population-{}-{}".format(phi, epsilon),
        "method_a_fingerprint": "method-a-{}-{}".format(phi, epsilon),
        "method_a_event_population_fingerprint": "method-a-population-{}-{}".format(phi, epsilon),
        "part1_config_fingerprint": "part1-{}-{}".format(phi, epsilon),
        "coordinate_fingerprint": "coordinate-{}-{}".format(phi, epsilon),
        "host_state": "post_proton", "source_target_state": "post_proton_noRF",
        "t_edges": [0.0, 1.0, 2.0, 3.0],
        "delta_edges": [-10.0, 0.0, 10.0],
        "phi_edges": [-180.0, 180.0],
    }
    artifact = {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract_artifact/v2",
        "setting": {
            "phi_setting": phi, "epsilon_filename_token": epsilon,
            "epsilon_setting": "low" if epsilon == "lowe" else "high",
            "particle_type": "kaon", "kinematic_token": "Q4p4W2p74",
            "Q2": 4.4, "W": 2.74,
        },
        "contract": contract, "non_authoritative": True,
        "production_objects_mutated": False, "refinement_applied": False,
        "production_application_performed": False, "event_application_performed": False,
    }
    _seal_f1_artifact(artifact)
    return artifact


def _seal_f1_artifact(artifact: dict[str, object]) -> None:
    """Rebuild the real frozen v2 fingerprints after an intentional test edit."""
    contract = artifact["contract"]
    assert isinstance(contract, dict)
    training = contract["method_a_training_records"]
    application = contract["application_records"]
    metadata = contract["feature_metadata"]
    assert isinstance(application, list)
    child_projection = [{
        "source_label": row["source_label"], "entry_index": row["entry_index"],
        "t_index": row["t_index"], "phi_index": row["phi_index"],
        "phi_low": row["phi_low"], "phi_high": row["phi_high"],
        "phi_status": row["phi_status"],
    } for row in application]
    contract["method_a_training_population_fingerprint"] = representation._sha256(training)
    contract["application_population_fingerprint"] = representation._sha256(application)
    contract["acceptance_feature_metadata_fingerprint"] = representation._sha256(metadata)
    contract["application_child_assignment_projection_fingerprint"] = representation._sha256(child_projection)
    fingerprint_inputs = {
        "schema_version": contract["schema_version"],
        "fingerprint_schema_version": contract["fingerprint_schema_version"],
        "phase_a_contract_fingerprint": contract["phase_a_contract_fingerprint"],
        "phase_a_pion_event_population_fingerprint": contract["phase_a_pion_event_population_fingerprint"],
        "method_a_fingerprint": contract["method_a_fingerprint"],
        "method_a_event_population_fingerprint": contract["method_a_event_population_fingerprint"],
        "part1_config_fingerprint": contract["part1_config_fingerprint"],
        "coordinate_fingerprint": contract["coordinate_fingerprint"],
        "host_state": contract["host_state"], "source_target_state": contract["source_target_state"],
        "t_edges": contract["t_edges"], "delta_edges": contract["delta_edges"], "phi_edges": contract["phi_edges"],
        "method_a_training_population_fingerprint": contract["method_a_training_population_fingerprint"],
        "application_population_fingerprint": contract["application_population_fingerprint"],
        "acceptance_feature_metadata_fingerprint": contract["acceptance_feature_metadata_fingerprint"],
        "application_child_assignment_projection_fingerprint": contract["application_child_assignment_projection_fingerprint"],
        "method_a_closure": contract["method_a_training_summary"]["by_t_delta"],
        "feature_metadata": metadata,
    }
    contract["fingerprint_inputs"] = fingerprint_inputs
    contract["fingerprint"] = representation._sha256(fingerprint_inputs)


def make_five_artifacts() -> list[dict[str, object]]:
    return [make_f1_artifact(phi, epsilon) for phi, epsilon in SETTINGS]


def input_hashes() -> dict[str, str]:
    return {"{}-{}".format(phi, epsilon): ("{:064x}".format(index + 1)) for index, (phi, epsilon) in enumerate(SETTINGS)}


class MethodAAcceptanceRepresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = make_five_artifacts()
        cls.hashes = input_hashes()

    def test_valid_five_setting_construction_and_input_immutability(self):
        original = deepcopy(self.artifacts)
        result = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=self.hashes,
        )
        self.assertEqual(self.artifacts, original)
        self.assertTrue(result["available"])
        self.assertEqual(len(result["groups"]), 15)
        self.assertEqual({summary["candidate_id"] for summary in result["candidate_summaries"]}, {
            "delta_only", "track3", "hgcer3", "full5_reference",
        })
        reference = next(summary for summary in result["candidate_summaries"] if summary["candidate_id"] == "full5_reference")
        self.assertFalse(reference["automatic_recommendation_eligible"])
        self.assertFalse(result["basis_frozen"])
        self.assertTrue(result["manual_review_required"])
        self.assertEqual(len(result["algorithm_fingerprint"]), 64)

    def test_row_reorder_is_scientifically_deterministic(self):
        reordered = deepcopy(self.artifacts)
        for artifact in reordered:
            artifact["contract"]["method_a_training_records"].reverse()
            artifact["contract"]["application_records"].reverse()
            _seal_f1_artifact(artifact)
        baseline = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=self.hashes,
        )
        changed = representation.build_pion_hgcer_method_a_acceptance_representation(
            reordered, input_file_hashes=self.hashes,
        )
        self.assertEqual(changed["fingerprint"], baseline["fingerprint"])
        self.assertEqual(changed["candidate_summaries"], baseline["candidate_summaries"])
        self.assertEqual(changed["recommendation"], baseline["recommendation"])

    def test_strict_schema_setting_and_population_validation(self):
        cases = []
        bad_schema = deepcopy(self.artifacts)
        bad_schema[0]["schema_version"] = "pion_hgcer_method_a_acceptance_event_contract_artifact/v1"
        cases.append((bad_schema, "f1_artifact_schema_version"))
        duplicate = deepcopy(self.artifacts)
        duplicate[1]["setting"]["phi_setting"] = "Left"
        duplicate[1]["setting"]["epsilon_filename_token"] = "lowe"
        cases.append((duplicate, "f1_setting_duplicate"))
        mixed = deepcopy(self.artifacts)
        mixed[1]["setting"]["kinematic_token"] = "Q3p0W2p32"
        cases.append((mixed, "f1_kinematic_mixed"))
        malformed_training = deepcopy(self.artifacts)
        malformed_training[0]["contract"]["method_a_training_records"][0]["nommcuts"] = False
        _seal_f1_artifact(malformed_training[0])
        cases.append((malformed_training, "selection_invalid"))
        malformed_application = deepcopy(self.artifacts)
        malformed_application[0]["contract"]["application_records"][0]["P_hgcer_npeSum"] = 2.0
        _seal_f1_artifact(malformed_application[0])
        cases.append((malformed_application, "npe_not_physical_control"))
        missing_prompt = deepcopy(self.artifacts)
        missing_prompt[0]["contract"]["application_records"][0]["entry_index"] = 999999
        _seal_f1_artifact(missing_prompt[0])
        cases.append((missing_prompt, "prompt_application_identity_missing_training_control"))
        for artifacts, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(representation.MethodAAcceptanceRepresentationError, error):
                    representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)

    def test_stale_serialized_f1_fingerprints_fail_before_metrics(self):
        cases = []
        stale_training = deepcopy(self.artifacts)
        stale_training[0]["contract"]["method_a_training_records"][0]["P_hgcer_npeSum"] = 1.5
        cases.append((stale_training, "method_a_training_population_fingerprint_mismatch"))
        stale_application = deepcopy(self.artifacts)
        stale_application[0]["contract"]["application_records"][0]["P_hgcer_npeSum"] = 4.0
        cases.append((stale_application, "application_population_fingerprint_mismatch"))
        stale_metadata = deepcopy(self.artifacts)
        stale_metadata[0]["contract"]["feature_metadata"]["parent_coordinate"] = "t"
        cases.append((stale_metadata, "feature_metadata_invalid"))
        stale_projection = deepcopy(self.artifacts)
        original_projection_fingerprint = stale_projection[0]["contract"]["application_child_assignment_projection_fingerprint"]
        stale_projection[0]["contract"]["application_records"][0]["phi_status"] = "outside"
        _seal_f1_artifact(stale_projection[0])
        stale_projection[0]["contract"]["application_child_assignment_projection_fingerprint"] = original_projection_fingerprint
        cases.append((stale_projection, "application_child_assignment_projection_fingerprint_mismatch"))
        stale_inputs = deepcopy(self.artifacts)
        stale_inputs[0]["contract"]["fingerprint_inputs"]["host_state"] = "tampered"
        cases.append((stale_inputs, "fingerprint_inputs_mismatch"))
        stale_final = deepcopy(self.artifacts)
        stale_final[0]["contract"]["fingerprint"] = "0" * 64
        cases.append((stale_final, "fingerprint_mismatch"))
        for artifacts, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(representation.MethodAAcceptanceRepresentationError, error):
                    representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)

    def test_exact_f1_authority_metadata_and_provenance_are_required(self):
        cases = []
        reason = deepcopy(self.artifacts)
        reason[0]["contract"]["reason"] = "not-available"
        cases.append((reason, "f1_contract_reason"))
        missing_reason = deepcopy(self.artifacts)
        del missing_reason[0]["contract"]["reason"]
        cases.append((missing_reason, "f1_contract_reason"))
        target = deepcopy(self.artifacts)
        target[0]["contract"]["source_target_state"] = "post_proton"
        _seal_f1_artifact(target[0])
        cases.append((target, "source_target_state"))
        missing_provenance = deepcopy(self.artifacts)
        missing_provenance[0]["contract"]["phase_a_contract_fingerprint"] = ""
        _seal_f1_artifact(missing_provenance[0])
        cases.append((missing_provenance, "phase_a_contract_fingerprint_invalid"))
        wrong_metadata = deepcopy(self.artifacts)
        wrong_metadata[0]["contract"]["feature_metadata"]["downstream_yield_coordinates"] = ["canonical_t"]
        _seal_f1_artifact(wrong_metadata[0])
        cases.append((wrong_metadata, "feature_metadata_invalid"))
        wrong_policy = deepcopy(self.artifacts)
        wrong_policy[0]["contract"]["feature_metadata"]["future_normalization_policy"] = "per_child"
        _seal_f1_artifact(wrong_policy[0])
        cases.append((wrong_policy, "feature_metadata_invalid"))
        for artifacts, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(representation.MethodAAcceptanceRepresentationError, error):
                    representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)

    def test_edge_arrays_and_canonical_t_geometry_are_frozen(self):
        cases = []
        nonfinite = deepcopy(self.artifacts)
        nonfinite[0]["contract"]["t_edges"][1] = float("nan")
        cases.append((nonfinite, "f1_t_edges_1_nonfinite"))
        nonincreasing = deepcopy(self.artifacts)
        nonincreasing[0]["contract"]["delta_edges"] = [-10.0, -10.0, 10.0]
        _seal_f1_artifact(nonincreasing[0])
        cases.append((nonincreasing, "f1_delta_edges_not_strictly_increasing"))
        inconsistent = deepcopy(self.artifacts)
        inconsistent[0]["contract"]["method_a_training_records"][0]["t_high"] = 9.0
        _seal_f1_artifact(inconsistent[0])
        cases.append((inconsistent, "t_geometry_mismatch"))
        mixed = deepcopy(self.artifacts)
        mixed_contract = mixed[1]["contract"]
        mixed_contract["t_edges"] = [0.0, 1.0, 2.0, 4.0]
        for row in mixed_contract["method_a_training_records"] + mixed_contract["application_records"]:
            if row["t_index"] == 2:
                row["t_high"] = 4.0
        _seal_f1_artifact(mixed[1])
        cases.append((mixed, "f1_canonical_t_geometry_mixed"))
        for artifacts, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(representation.MethodAAcceptanceRepresentationError, error):
                    representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)

    def test_population_selection_and_identity_remain_strict_after_resealing(self):
        cases = []
        wrong_selection = deepcopy(self.artifacts)
        wrong_selection[0]["contract"]["method_a_training_records"][0]["nommcuts"] = False
        _seal_f1_artifact(wrong_selection[0])
        cases.append((wrong_selection, "selection_invalid"))
        wrong_class = deepcopy(self.artifacts)
        wrong_class[0]["contract"]["method_a_training_records"][0]["response_class"] = "control"
        _seal_f1_artifact(wrong_class[0])
        cases.append((wrong_class, "response_classification_invalid"))
        nonpositive = deepcopy(self.artifacts)
        nonpositive[0]["contract"]["method_a_training_records"][0]["P_hgcer_npeSum"] = 0.0
        _seal_f1_artifact(nonpositive[0])
        cases.append((nonpositive, "response_classification_invalid"))
        nonphysical_application = deepcopy(self.artifacts)
        nonphysical_application[0]["contract"]["application_records"][0]["P_hgcer_npeSum"] = 2.0
        _seal_f1_artifact(nonphysical_application[0])
        cases.append((nonphysical_application, "npe_not_physical_control"))
        missing_prompt = deepcopy(self.artifacts)
        missing_prompt[0]["contract"]["application_records"][0]["entry_index"] = 999999
        _seal_f1_artifact(missing_prompt[0])
        cases.append((missing_prompt, "prompt_application_identity_missing_training_control"))
        duplicate = deepcopy(self.artifacts)
        duplicate[0]["contract"]["application_records"][1]["entry_index"] = duplicate[0]["contract"]["application_records"][0]["entry_index"]
        _seal_f1_artifact(duplicate[0])
        cases.append((duplicate, "record_identity_duplicate"))
        for artifacts, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(representation.MethodAAcceptanceRepresentationError, error):
                    representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)

    def test_insufficient_support_does_not_pool_canonical_t(self):
        artifacts = deepcopy(self.artifacts)
        rows = artifacts[0]["contract"]["method_a_training_records"]
        artifacts[0]["contract"]["method_a_training_records"] = [
            row for row in rows if not (row["t_index"] == 0 and row["response_class"] == "low" and row["entry_index"] % 1000 >= 5)
        ]
        _seal_f1_artifact(artifacts[0])
        result = representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)
        self.assertEqual(result["recommendation"]["recommendation_status"], "no_supported_reduced_basis")
        self.assertIsNone(result["recommendation"]["recommended_basis"])
        self.assertFalse(result["response_support"]["all_groups_satisfy_minimum"])

    def test_recommendation_states_do_not_promote_full5_or_break_ties(self):
        def summary(candidate_id, dimension, passed, eligible=None):
            return {
                "candidate_id": candidate_id, "dimension": dimension,
                "overall_candidate_passed": passed,
                "automatic_recommendation_eligible": (
                    candidate_id != "full5_reference" if eligible is None else eligible
                ),
            }
        self.assertEqual(
            representation._recommendation([], [summary("full5_reference", 5, False)], True)["recommendation_status"],
            "no_supported_reduced_basis",
        )
        self.assertEqual(
            representation._recommendation([], [summary("delta_only", 1, True)], True)["recommendation_status"],
            "unique_supported_reduced_basis",
        )
        self.assertEqual(
            representation._recommendation([], [summary("track3", 3, True), summary("hgcer3", 3, True)], True)["recommendation_status"],
            "multiple_supported_reduced_bases",
        )
        self.assertEqual(
            representation._recommendation([], [summary("delta_only", 1, True), summary("track3", 3, True)], True)["recommendation_status"],
            "unique_minimum_dimension_supported_basis",
        )
        self.assertEqual(
            representation._recommendation([], [summary("full5_reference", 5, True, False)], True)["recommendation_status"],
            "no_supported_reduced_basis",
        )

    def test_support_ood_and_zero_variance_fail_closed(self):
        far = deepcopy(self.artifacts)
        for row in far[0]["contract"]["application_records"]:
            if row["source_label"] != "prompt":
                row["SHMS_delta"] = 1.0e9
        _seal_f1_artifact(far[0])
        result = representation.build_pion_hgcer_method_a_acceptance_representation(far)
        delta = next(summary for summary in result["candidate_summaries"] if summary["candidate_id"] == "delta_only")
        self.assertFalse(delta["application_support_gate_passed"])
        zero = deepcopy(self.artifacts)
        for artifact in zero:
            for row in artifact["contract"]["method_a_training_records"]:
                row["SHMS_delta"] = 0.0
            _seal_f1_artifact(artifact)
        result = representation.build_pion_hgcer_method_a_acceptance_representation(zero)
        metric = result["groups"][0]["candidate_metrics"]["delta_only"]
        self.assertFalse(metric["valid"])
        self.assertEqual(metric["invalid_reason"], "fold_scaling_invalid")

    def test_nonprompt_sparse_flags_include_zero_and_do_not_fabricate_distances(self):
        zero = deepcopy(self.artifacts)
        zero[0]["contract"]["application_records"] = [
            row for row in zero[0]["contract"]["application_records"]
            if row["t_index"] != 0 or row["source_label"] == "prompt"
        ]
        _seal_f1_artifact(zero[0])
        zero_result = representation.build_pion_hgcer_method_a_acceptance_representation(zero)
        zero_support = zero_result["groups"][0]["candidate_metrics"]["delta_only"]["application_support"]
        self.assertEqual(zero_support["nonprompt_application_count"], 0)
        self.assertTrue(zero_support["statistically_sparse"])
        self.assertIsNone(zero_support["application_nn_p50"])

        sparse = deepcopy(self.artifacts)
        sparse[0]["contract"]["application_records"] = [
            row for row in sparse[0]["contract"]["application_records"]
            if row["t_index"] != 0 or row["source_label"] == "prompt" or row["entry_index"] % 10000 == 0
        ]
        _seal_f1_artifact(sparse[0])
        sparse_result = representation.build_pion_hgcer_method_a_acceptance_representation(sparse)
        sparse_support = sparse_result["groups"][0]["candidate_metrics"]["delta_only"]["application_support"]
        self.assertEqual(sparse_support["nonprompt_application_count"], 1)
        self.assertTrue(sparse_support["statistically_sparse"])

        full = representation.build_pion_hgcer_method_a_acceptance_representation(self.artifacts)
        full_support = full["groups"][0]["candidate_metrics"]["delta_only"]["application_support"]
        self.assertEqual(full_support["nonprompt_application_count"], 20)
        self.assertFalse(full_support["statistically_sparse"])

    def test_optimizer_failure_is_an_explicit_closed_failure(self):
        failed = SimpleNamespace(success=False, x=None)
        with patch.object(representation, "minimize", return_value=failed):
            result = representation.build_pion_hgcer_method_a_acceptance_representation(self.artifacts)
        metric = result["groups"][0]["candidate_metrics"]["delta_only"]
        self.assertFalse(metric["valid"])
        self.assertEqual(metric["invalid_reason"], "optimizer_not_converged")

    def test_fingerprint_sensitivity_and_no_event_level_outputs(self):
        baseline = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=self.hashes,
        )
        altered_hashes = dict(self.hashes)
        altered_hashes["Left-lowe"] = "f" * 64
        changed_hash = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=altered_hashes,
        )
        changed_config = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=self.hashes, algorithm_config={"regularization_lambda": 2.0e-3},
        )
        changed_threshold = representation.build_pion_hgcer_method_a_acceptance_representation(
            self.artifacts, input_file_hashes=self.hashes,
            algorithm_config={"information_max_auc_loss_max": 0.04},
        )
        self.assertNotEqual(baseline["fingerprint"], changed_hash["fingerprint"])
        self.assertNotEqual(baseline["fingerprint"], changed_config["fingerprint"])
        self.assertNotEqual(baseline["fingerprint"], changed_threshold["fingerprint"])

        forbidden = {
            "event_score", "prediction", "predicted_probability", "probability",
            "correction", "C_A", "adjusted_weight", "map", "method_b_numerical_input",
        }
        def walk(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    self.assertNotIn(key, forbidden)
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)
        walk(baseline)
        self.assertFalse(baseline["method_b_numerical_dependency"])
        self.assertFalse(baseline["basis_frozen"])
        self.assertTrue(baseline["manual_review_required"])

    def test_artifact_writer_and_filename(self):
        artifact = representation.build_pion_hgcer_method_a_acceptance_representation_artifact(
            self.artifacts, input_file_hashes=self.hashes, input_paths={key: key + ".json" for key in self.hashes},
        )
        self.assertEqual(artifact["schema_version"], "pion_hgcer_method_a_acceptance_representation_artifact/v1")
        self.assertEqual(
            representation.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74"),
            "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-representation.json",
        )


if __name__ == "__main__":
    unittest.main()
