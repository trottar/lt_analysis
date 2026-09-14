"""Focused pure-Python tests for the detached F.2 representation audit."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest


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
        "phi_is_training_feature": False, "probability_map_constructed": False,
        "absolute_leakage_probability_claimed": False, "weight_adjustment_constructed": False,
        "method_b_numerical_dependency": False,
    }
    contract = {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract/v2",
        "fingerprint_schema_version": "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2",
        "status": "available", "available": True, "diagnostic_stage": "complete",
        "non_authoritative": True, "method_b_numerical_dependency": False,
        "event_application_performed": False, "production_application_performed": False,
        "production_objects_mutated": False, "future_weight_adjustment_constructed": False,
        "refinement_applied": False, "feature_metadata": metadata,
        "method_a_training_records": training, "application_records": application,
        "method_a_training_population_fingerprint": "training-{}-{}".format(phi, epsilon),
        "application_population_fingerprint": "application-{}-{}".format(phi, epsilon),
        "acceptance_feature_metadata_fingerprint": "metadata-{}-{}".format(phi, epsilon),
        "application_child_assignment_projection_fingerprint": "child-{}-{}".format(phi, epsilon),
        "coordinate_fingerprint": "coordinate-{}-{}".format(phi, epsilon),
        "fingerprint": "f1-{}-{}".format(phi, epsilon),
    }
    return {
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
        cases.append((malformed_training, "selection_invalid"))
        malformed_application = deepcopy(self.artifacts)
        malformed_application[0]["contract"]["application_records"][0]["P_hgcer_npeSum"] = 2.0
        cases.append((malformed_application, "npe_not_physical_control"))
        missing_prompt = deepcopy(self.artifacts)
        missing_prompt[0]["contract"]["application_records"][0]["entry_index"] = 999999
        cases.append((missing_prompt, "prompt_application_identity_missing_training_control"))
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
        result = representation.build_pion_hgcer_method_a_acceptance_representation(artifacts)
        self.assertEqual(result["recommendation"]["recommendation_status"], "insufficient_response_support")
        self.assertIsNone(result["recommendation"]["recommended_basis"])

    def test_recommendation_states_do_not_promote_full5_or_break_ties(self):
        def summary(candidate_id, dimension, passed):
            return {"candidate_id": candidate_id, "dimension": dimension, "overall_candidate_passed": passed}
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

    def test_support_ood_and_zero_variance_fail_closed(self):
        far = deepcopy(self.artifacts)
        for row in far[0]["contract"]["application_records"]:
            if row["source_label"] != "prompt":
                row["SHMS_delta"] = 1.0e9
        result = representation.build_pion_hgcer_method_a_acceptance_representation(far)
        delta = next(summary for summary in result["candidate_summaries"] if summary["candidate_id"] == "delta_only")
        self.assertFalse(delta["application_support_gate_passed"])
        zero = deepcopy(self.artifacts)
        for artifact in zero:
            for row in artifact["contract"]["method_a_training_records"]:
                row["SHMS_delta"] = 0.0
        result = representation.build_pion_hgcer_method_a_acceptance_representation(zero)
        metric = result["groups"][0]["candidate_metrics"]["delta_only"]
        self.assertFalse(metric["valid"])
        self.assertEqual(metric["invalid_reason"], "fold_scaling_invalid")

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

        forbidden = {"event_score", "prediction", "predicted_probability", "probability", "correction", "C_A", "adjusted_weight"}
        def walk(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    self.assertNotIn(key, forbidden)
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)
        walk(baseline)

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
