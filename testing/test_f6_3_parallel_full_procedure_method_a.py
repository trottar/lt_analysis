"""Deterministic F.6.3 authority and multiplier regressions (no ROOT/farm)."""

from __future__ import annotations

import copy
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import inspect
import math
import sys
import unittest
import tempfile
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src" / "cuts"), str(ROOT / "src" / "utility")]

import pion_component_subtraction as pion
import pion_hgcer_method_a_parallel_full_procedure as f63
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5
import pion_hgcer_method_a_acceptance_representation as f2
import pion_hgcer_method_a_acceptance_map as f3_builder
from testing.test_e8_2_baseline_stage_audit import _Histogram as _YieldHistogram
from testing.test_e8_2_baseline_stage_audit import _load_calculate_yield_module
from testing.test_e8_2_baseline_stage_audit import _public_yield_fixture
from testing import test_pion_hgcer_method_a_tphi_propagation as f5_fixtures


class _Axis:
    def FindBin(self, value):
        return 1


class _Reference:
    def GetXaxis(self):
        return _Axis()


class _Histogram:
    def __init__(self):
        self.fills = []

    def Fill(self, *values):
        self.fills.append(tuple(float(value) for value in values))


def _row(identity=("prompt", 7), *, coefficient=2.0, weight=3.0):
    return {
        "source_label": identity[0], "entry_index": identity[1],
        "t_index": 0, "phi_index": 0, "analysis_MM": 1.12,
        "analysis_t": 0.5, "signed_source_coefficient": coefficient,
        "baseline_pion_weight_w0": weight,
        "signed_baseline_event_contribution": coefficient * weight,
    }


def _accepted_rows(*rows):
    return {
        (str(row["source_label"]), int(row["entry_index"])): {
            name: row[name]
            for name in (
                "t_index", "phi_index", "analysis_MM", "analysis_t",
                "signed_source_coefficient", "baseline_pion_weight_w0",
                "signed_baseline_event_contribution",
            )
        }
        for row in rows
    }


def _raw_f1(*rows):
    return [{"setting": {"phi_setting": "Left", "epsilon_filename_token": "lowe"},
             "contract": {"application_records": list(rows)}}]


def _cache():
    return {
        "entry_index": np.asarray([7], dtype=np.int32),
        "adj_MM": np.asarray([1.12]), "adj_t": np.asarray([0.5]),
        "Q2": np.asarray([1.0]), "W": np.asarray([2.0]),
        "epsilon": np.asarray([0.6]), "theta_cm_deg": np.asarray([30.0]),
        "ssxptar": np.asarray([0.0]), "ssyptar": np.asarray([0.0]),
        "hsxptar": np.asarray([0.0]), "hsyptar": np.asarray([0.0]),
        "allcuts": np.asarray([True]), "nommcuts": np.asarray([True]),
        "t_index": np.asarray([0], dtype=np.int32), "phi_index": np.asarray([0], dtype=np.int32),
        "coefficient": np.asarray([2.0]),
        "allcut_bin_index": {(0, 0): np.asarray([0], dtype=np.int32)},
        "nommcut_bin_index": {(0, 0): np.asarray([0], dtype=np.int32)},
    }


class CandidateLineageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Synthetic inputs exercise the real unchanged F.4/F.5 validators and
        # calculator. They do not stand in for the farm's candidate artifacts.
        cls.artifacts = f5_fixtures._artifacts()
        for artifact in cls.artifacts:
            # Only detector coordinates distinguish response classes. The
            # unchanged F.2 builder therefore selects hgcer3 without overrides.
            for row in artifact["contract"]["method_a_training_records"] + artifact["contract"]["application_records"]:
                spread = 0.01 * (row["entry_index"] % 25) + 0.001 * row["t_index"]
                row.update(SHMS_delta=spread, SHMS_xptar=spread, SHMS_yptar=0.5 * spread)
            for row in artifact["contract"]["application_records"]:
                row["analysis_MM"] = 1.12
                row["analysis_t"] = (row["t_low"] + row["t_high"]) / 2.0
            f5_fixtures.f1_fixtures._seal_f1_artifact(artifact)
        cls.hashes = f5_fixtures.f1_fixtures.input_hashes()
        cls.f2_artifact = f2.build_pion_hgcer_method_a_acceptance_representation_artifact(
            cls.artifacts, input_file_hashes=cls.hashes, input_paths={})
        cls.f2_sha = f63._sha256_bytes(f63._writer_bytes(cls.f2_artifact))
        cls.f2_authority = {"source_file_sha256": cls.f2_sha,
            "representation_fingerprint": cls.f2_artifact["representation"]["fingerprint"],
            "artifact_fingerprint": cls.f2_artifact["artifact_fingerprint"]}
        cls.f3 = f3_builder.build_pion_hgcer_method_a_acceptance_map_artifact(
            cls.artifacts, cls.f2_artifact, f1_input_file_hashes=cls.hashes,
            f2_input_file_sha256=cls.f2_sha, input_paths={})
        cls.f3_sha = hashlib.sha256(json.dumps(cls.f3, sort_keys=True).encode()).hexdigest()
        cls.f3_authority = f5_fixtures.f4_fixtures._authority(cls.f3, cls.f3_sha)
        cls.f3_authority["Q4p4W2p74"]["farm_source_head"] = "0" * 40
        cls.f4_artifact = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(
            cls.artifacts, cls.f3, f1_input_file_hashes=cls.hashes,
            f3_input_file_sha256=cls.f3_sha,
            accepted_f3_runtime_authority_by_kinematic=cls.f3_authority,
            input_paths={"f1": {}, "f3": "candidate-three"},
        )
        cls.f4_sha = hashlib.sha256(json.dumps(cls.f4_artifact, sort_keys=True).encode()).hexdigest()
        correction = cls.f4_artifact["correction"]
        cls.f4_authority = {"Q4p4W2p74": {
            "source_file_sha256": cls.f4_sha,
            "correction_fingerprint": correction["fingerprint"],
            "artifact_fingerprint": cls.f4_artifact["artifact_fingerprint"],
            "farm_source_head": f63.F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
            **{name: correction[name] for name in (
                "f3_source_file_sha256", "f3_map_fingerprint",
                "f3_algorithm_fingerprint", "f3_artifact_fingerprint",
            )},
            "f1_source_file_sha256": cls.hashes,
        }}

    def _reconstruct(self, *, f4_sha=None, f3=None, f3_sha=None, hashes=None, artifacts=None, f2_value=None, f2_sha=None):
        with mock.patch.object(f63, "F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC", self.f3_authority), \
             mock.patch.object(f63, "F6_3_CANDIDATE_F2_AUTHORITY", self.f2_authority), \
             mock.patch.object(f63, "F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC", self.f4_authority):
            return f63.reconstruct_transient_factor_map(
                self.artifacts if artifacts is None else artifacts, self.f3 if f3 is None else f3, self.f4_artifact,
                f2_artifact=self.f2_artifact if f2_value is None else f2_value,
                f2_input_file_sha256=self.f2_sha if f2_sha is None else f2_sha,
                f1_input_file_hashes=self.hashes if hashes is None else hashes,
                f3_input_file_sha256=self.f3_sha if f3_sha is None else f3_sha,
                f4_input_file_sha256=self.f4_sha if f4_sha is None else f4_sha,
                setting_id="Left-lowe",
            )

    def test_exact_candidate_paths_and_canonical_f1_inventory(self):
        paths = f63.accepted_f6_3_artifact_paths(ROOT / "OUTPUT", "Q4p4W2p74")
        self.assertEqual(Path(paths["f2"]).name,
                         "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json")
        self.assertEqual(Path(paths["f3"]).name,
                         "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json")
        self.assertEqual(Path(paths["f4"]).name,
                         "Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json")
        _, _, f1_filename, _, _ = f63._runtime_dependencies()
        self.assertEqual(paths["f1"], {
            "{}-{}".format(phi, epsilon): str(ROOT / "OUTPUT" / f1_filename(phi, "kaon", "Q4p4W2p74", epsilon))
            for phi, epsilon in (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))
        })
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "kinematic_unsupported"):
            f63.accepted_f6_3_artifact_paths(ROOT, "Q3p0W2p32")

    def test_exact_farm_candidate_identity_pins_and_reconstruction_sentinel(self):
        f3_record = f63.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC
        self.assertEqual(f3_record, {"Q4p4W2p74": {
            "source_file_sha256": "c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228",
            "map_fingerprint": "6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548",
            "algorithm_fingerprint": "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912",
            "artifact_fingerprint": "e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302",
            "farm_source_head": "0" * 40,
        }})
        hashes = {
            "Left-lowe": "eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07",
            "Left-highe": "544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e",
            "Center-lowe": "2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16",
            "Center-highe": "c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941",
            "Right-highe": "e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652",
        }
        self.assertEqual(f63.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256, hashes)
        self.assertEqual(f63.F6_3_CANDIDATE_F2_AUTHORITY, {
            "source_file_sha256": "2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e",
            "representation_fingerprint": "e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216",
            "artifact_fingerprint": "87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6"})
        self.assertEqual(f63.F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
                         "463d2657f696ecee33113edc3393ac51083a8944")
        self.assertEqual(f63.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC, {"Q4p4W2p74": {
            "source_file_sha256": "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7",
            "correction_fingerprint": "71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98",
            "artifact_fingerprint": "0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4",
            "farm_source_head": "463d2657f696ecee33113edc3393ac51083a8944",
            "f3_source_file_sha256": f3_record["Q4p4W2p74"]["source_file_sha256"],
            "f3_map_fingerprint": f3_record["Q4p4W2p74"]["map_fingerprint"],
            "f3_algorithm_fingerprint": f3_record["Q4p4W2p74"]["algorithm_fingerprint"],
            "f3_artifact_fingerprint": f3_record["Q4p4W2p74"]["artifact_fingerprint"],
            "f1_source_file_sha256": hashes,
        }})

    def test_real_candidate_reproduction_uses_explicit_records_and_preserves_historical_authorities(self):
        historical_f3 = copy.deepcopy(f4.ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC)
        historical_f4 = copy.deepcopy(f5.ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC)
        inputs = copy.deepcopy((self.artifacts, self.f3, self.f4_artifact))
        with mock.patch.object(f5, "_validate_f4_artifact", wraps=f5._validate_f4_artifact) as validate, \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data",
                               wraps=f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data) as build:
            factors, provenance, rows = self._reconstruct()
        self.assertIs(validate.call_args.args[2], self.f4_authority)
        self.assertIs(build.call_args.kwargs["accepted_f3_runtime_authority_by_kinematic"], self.f3_authority)
        self.assertEqual(self.f3_authority["Q4p4W2p74"]["farm_source_head"], "0" * 40)
        self.assertEqual(provenance["candidate_f3_reconstruction_role"], "candidate_construction_sentinel_only")
        self.assertEqual(provenance["branch_role"], "parallel_nonproduction_method_a_full_analysis")
        self.assertTrue(provenance["current_baseline_candidate_lineage"])
        self.assertFalse(provenance["production_promotion_performed"])
        self.assertFalse(provenance["event_correction_persisted"])
        self.assertEqual(provenance["candidate_f3_source_file_sha256"], self.f3_sha)
        self.assertEqual(provenance["candidate_f4_source_file_sha256"], self.f4_sha)
        self.assertEqual(provenance["scientific_equivalence"]["mode"], "exact-lineage")
        self.assertTrue(provenance["scientific_equivalence"]["all_stages_passed"])
        self.assertEqual(provenance["reviewed_candidate"]["f2"], self.f2_authority)
        self.assertEqual(provenance["current_runtime_lineage"]["f1_source_file_sha256"], self.hashes)
        self.assertEqual(set(factors), set(rows))
        self.assertTrue(all(math.isfinite(value) and value > 0.0 for value in factors.values()))
        live_rows = [{"source_label": source, "entry_index": entry, **row}
                     for (source, entry), row in rows.items()]
        self.assertTrue(f63.validate_live_cache_parity(factors, live_rows, rows)["live_cache_parity_passed"])
        raw_rows = {(row["source_label"], row["entry_index"]): row
                    for row in self.artifacts[0]["contract"]["application_records"]}
        for identity, row in rows.items():
            for field in ("analysis_MM", "analysis_t"):
                self.assertEqual(row[field], raw_rows[identity][field])
        live_rows[0]["analysis_MM"] += 0.01
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f1_analysis_MM_mismatch"):
            f63.validate_live_cache_parity(factors, live_rows, rows)
        self.assertEqual((self.artifacts, self.f3, self.f4_artifact), inputs)
        self.assertEqual(f4.ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC, historical_f3)
        self.assertEqual(f5.ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC, historical_f4)
        # Parent-t closure is preserved by the actual unchanged calculator.
        for t_index in range(3):
            identities = [key for key, row in rows.items() if row["t_index"] == t_index]
            baseline = sum(rows[key]["signed_baseline_event_contribution"] for key in identities)
            adjusted = sum(rows[key]["signed_baseline_event_contribution"] * factors[key] for key in identities)
            self.assertAlmostEqual(adjusted, baseline, places=10)
        self.assertNotIn("correction_factors", json.dumps(provenance))

    def test_wrong_candidate_raw_hashes_fail_closed(self):
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f4_runtime_authority_source_file_sha256_mismatch"):
            self._reconstruct(f4_sha="d" * 64)
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f3_runtime_authority_source_file_sha256_mismatch"):
            self._reconstruct(f3_sha="d" * 64)
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f2_runtime_authority_source_file_sha256_mismatch"):
            self._reconstruct(f2_sha="d" * 64)

    def test_wrong_candidate_f3_fingerprints_fail_closed(self):
        for field in ("map_fingerprint", "algorithm_fingerprint", "artifact_fingerprint"):
            authority = copy.deepcopy(self.f3_authority)
            authority["Q4p4W2p74"][field] = "d" * 64
            with self.subTest(field=field), mock.patch.object(self, "f3_authority", authority), \
                 self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f3_runtime_authority_{}_mismatch".format(field)):
                self._reconstruct()

    def test_wrong_inherited_candidate_f3_and_f1_pins_fail_f5_validation(self):
        for field in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint", "f1_source_file_sha256"):
            authority = copy.deepcopy(self.f4_authority)
            authority["Q4p4W2p74"][field] = {"Left-lowe": "d" * 64} if field == "f1_source_file_sha256" else "d" * 64
            with self.subTest(field=field), mock.patch.object(self, "f4_authority", authority), \
                 self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f4_runtime_authority_{}_mismatch".format(field)):
                self._reconstruct()

    def test_nonzero_reconstruction_head_cannot_reproduce_candidate(self):
        authority = copy.deepcopy(self.f3_authority)
        authority["Q4p4W2p74"]["farm_source_head"] = f63.F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD
        with mock.patch.object(self, "f3_authority", authority), \
             self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "shared_reproduction_mismatch"):
            self._reconstruct()

    def _current_population(self):
        current = copy.deepcopy(self.artifacts)
        for artifact in current:
            artifact["contract"]["phase_a_contract_fingerprint"] += "-new-lineage"
            f5_fixtures.f1_fixtures._seal_f1_artifact(artifact)
        hashes = {f"{a['setting']['phi_setting']}-{a['setting']['epsilon_filename_token']}":
                  f63._sha256_bytes(f63._writer_bytes(a)) for a in current}
        return current, hashes

    def test_real_equivalent_new_raw_and_stable_lineage_uses_current_builders_and_factors(self):
        current, hashes = self._current_population()
        old_stable = f3_builder._validate_f1_artifacts(self.artifacts)
        with mock.patch.object(f2, "build_pion_hgcer_method_a_acceptance_representation_artifact",
                               wraps=f2.build_pion_hgcer_method_a_acceptance_representation_artifact) as build2, \
             mock.patch.object(f3_builder, "build_pion_hgcer_method_a_acceptance_map_artifact",
                               wraps=f3_builder.build_pion_hgcer_method_a_acceptance_map_artifact) as build3, \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data",
                               wraps=f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data) as build4, \
             mock.patch.object(f4, "_validate_f3_artifact", wraps=f4._validate_f3_artifact) as strict:
            factors, provenance, rows = self._reconstruct(artifacts=current, hashes=hashes)
        self.assertEqual((build2.call_count, build3.call_count, build4.call_count, strict.call_count), (1, 1, 1, 1))
        self.assertIs(build2.call_args.args[0], current)
        self.assertIs(build3.call_args.args[0], current)
        self.assertIs(build4.call_args.args[0], current)
        self.assertIsNot(build4.call_args.args[1], self.f3)
        authority = build4.call_args.kwargs["accepted_f3_runtime_authority_by_kinematic"]["Q4p4W2p74"]
        self.assertEqual(authority["source_file_sha256"], build4.call_args.kwargs["f3_input_file_sha256"])
        self.assertNotEqual(authority["map_fingerprint"], self.f3["acceptance_map"]["fingerprint"])
        self.assertEqual(provenance["scientific_equivalence"]["mode"], "equivalent-new-lineage")
        self.assertTrue(provenance["scientific_equivalence"]["all_stages_passed"])
        self.assertEqual(provenance["current_runtime_lineage"]["f1_source_file_sha256"], hashes)
        self.assertEqual(provenance["reviewed_candidate"]["f1_source_file_sha256"], self.hashes)
        stable = provenance["current_runtime_lineage"]["f1_stable_content_fingerprints"]
        self.assertTrue(all(stable[p["setting_id"]] != p["stable_content_fingerprint"] for p in old_stable))
        old_factors, _, _ = self._reconstruct()
        self.assertEqual(factors, old_factors)
        self.assertEqual(set(rows), set(factors))
        self.assertNotIn("correction_factors", json.dumps(provenance))
        # The original ordinary validator still rejects stale F.3 with new F.1.
        with self.assertRaisesRegex(f4.MethodAParentPreservingCorrectionError, "input_content_mismatch"):
            f4.build_pion_hgcer_method_a_parent_preserving_correction(current, self.f3,
                f1_input_file_hashes=hashes, f3_input_file_sha256=self.f3_sha,
                accepted_f3_runtime_authority_by_kinematic=self.f3_authority)

    def test_every_required_scientific_difference_stops_before_later_stages(self):
        current, hashes = self._current_population()
        current2 = f2.build_pion_hgcer_method_a_acceptance_representation_artifact(current, input_file_hashes=hashes, input_paths={})
        current3 = f3_builder.build_pion_hgcer_method_a_acceptance_map_artifact(current, current2,
            f1_input_file_hashes=hashes, f2_input_file_sha256=f63._sha256_bytes(f63._writer_bytes(current2)), input_paths={})
        sha3 = f63._sha256_bytes(f63._writer_bytes(current3))
        authority3 = f5_fixtures.f4_fixtures._authority(current3, sha3)
        authority3["Q4p4W2p74"]["farm_source_head"] = "0" * 40
        current4, review = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data(
            current, current3, f1_input_file_hashes=hashes, f3_input_file_sha256=sha3,
            accepted_f3_runtime_authority_by_kinematic=authority3, input_paths={})
        fields = {
            "f2": ("candidate_definitions", "algorithm_config", "response_support", "groups", "candidate_summaries", "recommendation", "new_scientific_field"),
            "f3": ("models", "accepted_basis", "ordered_features", "algorithm_config", "algorithm_fingerprint"),
            "f4": ("parents", "method_b_numerical_dependency", "production_application_performed", "child_renormalization_performed"),
        }
        bodies = {"f2": "representation", "f3": "acceptance_map", "f4": "correction"}
        base = {"f2": current2, "f3": current3, "f4": current4}
        cases = [(stage, field, None) for stage, names in fields.items() for field in names]
        # Real model/scaler/support shapes and nested parent diagnostics stay compared.
        cases += [("f3", "models", field) for field in ("logistic_fit", "scaler", "application_support", "ood_test_field")]
        cases += [("f4", "parents", field) for field in ("parent_normalization", "correction_factor_summary", "source_diagnostics", "canonical_phi_diagnostics")]
        cases += [("f4", "parents", mode) for mode in ("missing", "extra")]
        for stage, field, nested in cases:
            changed = copy.deepcopy(base)
            body = changed[stage][bodies[stage]]
            if nested in ("missing", "extra"):
                if nested == "missing": body[field].pop()
                else: body[field].append(copy.deepcopy(body[field][0]))
            elif stage == "f3" and nested == "logistic_fit":
                body[field][0][nested]["coefficients"][0] += 0.01
            elif stage == "f3" and nested == "scaler":
                body[field][0][nested]["median"][0] += 0.01
            elif stage == "f3" and nested == "application_support":
                body[field][0][nested]["application_ood_fraction"] += 0.01
            elif nested:
                body[field][0][nested] = {"synthetic_scientific_change": True}
            else:
                body[field] = {"synthetic_scientific_change": True}
            with self.subTest(stage=stage, field=field, nested=nested), \
                 mock.patch.object(f2, "build_pion_hgcer_method_a_acceptance_representation_artifact", return_value=changed["f2"]) as b2, \
                 mock.patch.object(f3_builder, "build_pion_hgcer_method_a_acceptance_map_artifact", return_value=changed["f3"]) as b3, \
                 mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data", return_value=(changed["f4"], review)) as b4, \
                 self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, f"equivalence_mismatch:{stage}:\\$\\.{field}"):
                self._reconstruct(artifacts=current, hashes=hashes)
            self.assertEqual(b2.call_count, 1)
            self.assertEqual(b3.call_count, 0 if stage == "f2" else 1)
            self.assertEqual(b4.call_count, 1 if stage == "f4" else 0)

    def test_candidate_wrappers_flags_fingerprints_and_setting_inventory_fail_closed(self):
        for stage, field in (("f2", "schema_version"), ("f2", "artifact_fingerprint"),
                ("f2", "non_authoritative"), ("f2", "method_b_numerical_dependency"),
                ("f3", "production_application_performed"), ("f3", "map_applied")):
            candidate = copy.deepcopy(self.f2_artifact if stage == "f2" else self.f3)
            candidate[field] = "invalid"
            with self.subTest(stage=stage, field=field), mock.patch.object(f2,
                    "build_pion_hgcer_method_a_acceptance_representation_artifact", side_effect=AssertionError("must not build")), \
                    self.assertRaises(f63.MethodAParallelFullProcedureError):
                self._reconstruct(**({"f2_value": candidate} if stage == "f2" else {"f3": candidate}))
        with self.assertRaises(f63.MethodAParallelFullProcedureError):
            self._reconstruct(artifacts=self.artifacts[:-1])
        wrong = copy.deepcopy(self.artifacts); wrong[0]["setting"]["phi_setting"] = "Right"
        with self.assertRaises(f63.MethodAParallelFullProcedureError):
            self._reconstruct(artifacts=wrong)

    def test_excluded_identity_differences_require_unchanged_science_and_direct_strict_checks(self):
        current, hashes = self._current_population()
        parsed = f3_builder._validate_f1_artifacts(current)
        for excluded in f63.SCIENTIFIC_PROVENANCE_EXCLUSIONS.values():
            for name in excluded:
                self.assertIsNone(f63.first_mismatch(
                    f63.scientific_projection({"science": 1, name: "old"}, excluded),
                    f63.scientific_projection({"science": 1, name: "new"}, excluded)))
        # Ordinary validator rejection covers raw, stable and contract lineage
        # independently, not just a changed aggregate fingerprint.
        cases = []
        raw = dict(self.hashes); raw["Left-lowe"] = "e" * 64
        cases.append((f3_builder._validate_f1_artifacts(self.artifacts), raw))
        stable = copy.deepcopy(f3_builder._validate_f1_artifacts(self.artifacts))
        stable[0]["stable_content_fingerprint"] = "e" * 64
        cases.append((stable, self.hashes))
        contract = copy.deepcopy(f3_builder._validate_f1_artifacts(self.artifacts))
        contract[0]["fingerprints"]["fingerprint"] = "e" * 64
        cases.append((contract, self.hashes))
        for population, raw_hashes in cases:
            with self.assertRaises(f4.MethodAParentPreservingCorrectionError):
                f4._validate_f3_artifact(self.f3, population, raw_hashes, f63._TOLERANCE)
        self.assertNotEqual(parsed[0]["stable_content_fingerprint"], stable[0]["stable_content_fingerprint"])

    def test_explicit_f2_loader_requires_exact_path_without_breaking_read_only_api(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = f63.accepted_f6_3_artifact_paths(directory, "Q4p4W2p74")
            for sid, path in paths["f1"].items():
                artifact = next(a for a in self.artifacts if f"{a['setting']['phi_setting']}-{a['setting']['epsilon_filename_token']}" == sid)
                Path(path).write_bytes(f63._writer_bytes(artifact))
            for stage, artifact in (("f2", self.f2_artifact), ("f3", self.f3), ("f4", self.f4_artifact)):
                Path(paths[stage]).write_bytes(f63._writer_bytes(artifact))
            loaded = f63.load_accepted_f6_3_authority(paths, include_f2=True)
            self.assertEqual(loaded[1], self.f2_artifact)
            self.assertEqual(loaded[-1]["f2"], self.f2_sha)
            self.assertEqual(len(f63.load_accepted_f6_3_authority(paths)), 4)
            Path(paths["f2"]).unlink()
            with self.assertRaises(f63.MethodAParallelFullProcedureError):
                f63.load_accepted_f6_3_authority(paths, include_f2=True)
            self.assertEqual(len(f63.load_accepted_f6_3_authority(paths)), 4)
            for malformed in (b'{"x":1,"x":2}', b'{"x":NaN}'):
                Path(paths["f2"]).write_bytes(malformed)
                with self.assertRaises(f63.MethodAParallelFullProcedureError):
                    f63.load_accepted_f6_3_authority(paths, include_f2=True)


class ParallelAuthorityTests(unittest.TestCase):
    def test_live_cache_parity_passes_without_retaining_factor_values(self):
        result = f63.validate_live_cache_parity(
            {("prompt", 7): 1.25}, [_row()], _accepted_rows(_row()),
        )
        self.assertTrue(result["live_cache_parity_passed"])
        self.assertNotIn("1.25", str(result))

    def test_live_cache_identity_and_baseline_mismatches_fail_closed(self):
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "inventory"):
            f63.validate_live_cache_parity(
                {("prompt", 8): 1.0}, [_row()], _accepted_rows(_row()),
            )
        changed = _row(); changed["baseline_pion_weight_w0"] = 4.0
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "baseline_identity"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [changed], _accepted_rows(_row()),
            )

    def test_live_cache_duplicate_nonfinite_and_nonpositive_fail_closed(self):
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "duplicate"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [_row(), _row()], _accepted_rows(_row()),
            )
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "nonfinite_or_nonpositive"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 0.0}, [_row()], _accepted_rows(_row()),
            )

    def test_reconstruction_rejects_length_duplicate_and_nonpositive_review_factors(self):
        fixture = CandidateLineageTests('test_wrong_candidate_raw_hashes_fail_closed')
        if not hasattr(fixture, 'f2_artifact'):
            CandidateLineageTests.setUpClass()
        real = f63._lineage_reconstruction
        for defect in ('length', 'duplicate', 'nonpositive'):
            def corrupted(*args, **kwargs):
                correction, review, lineage, decision = real(*args, **kwargs)
                review = copy.deepcopy(review)
                selected = [r for r in review if r['setting_id'] == 'Left-lowe']
                if defect == 'length':
                    selected[0]['correction_factors'] = np.asarray([1.0])
                elif defect == 'duplicate':
                    review.append(copy.deepcopy(selected[0]))
                else:
                    selected[0]['correction_factors'][0] = 0.0
                return correction, review, lineage, decision
            reason = {'length': 'row_alignment', 'duplicate': 'identity_duplicate',
                      'nonpositive': 'nonfinite_or_nonpositive'}[defect]
            with self.subTest(defect=defect), mock.patch.object(f63, '_lineage_reconstruction', side_effect=corrupted), self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, reason):
                fixture._reconstruct()

    def test_live_cache_requires_every_f1_baseline_field_to_match(self):
        reference = _row()
        for field, changed_value in (
            ("signed_source_coefficient", 3.0),
            ("baseline_pion_weight_w0", 4.0),
            ("signed_baseline_event_contribution", 7.0),
            ("t_index", 1), ("phi_index", 1),
            ("analysis_MM", 1.13), ("analysis_t", 0.51),
        ):
            with self.subTest(field=field):
                changed = _row(); changed[field] = changed_value
                if field == "signed_source_coefficient":
                    changed["signed_baseline_event_contribution"] = changed_value * changed["baseline_pion_weight_w0"]
                with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "(f1_.*_mismatch|baseline_identity_mismatch)"):
                    f63.validate_live_cache_parity(
                        {("prompt", 7): 1.0}, [changed], _accepted_rows(reference),
                    )

    def test_live_cache_uses_accepted_f4_scaled_tolerance(self):
        accepted = _row()
        accepted["analysis_MM"] = 1.0e6
        within = _row()
        within["analysis_MM"] = 1.0e6 + 5.0e-7
        f63.validate_live_cache_parity(
            {("prompt", 7): 1.0}, [within], _accepted_rows(accepted),
        )
        outside = _row()
        outside["analysis_MM"] = 1.0e6 + 2.0e-6
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "analysis_MM"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [outside], _accepted_rows(accepted),
            )


class TemplateMultiplierTests(unittest.TestCase):
    def _fill(self, multipliers=None):
        templates = {"mm": _Histogram(), "mm_nosub": _Histogram()}
        pion.fill_simc_shape_pion_subtraction_templates(
            templates,
            [{"label": "prompt", "cache_section": _cache(), "coefficient": 2.0, "base_coefficient": 2.0}],
            _Reference(), np.asarray([0.0, 3.0]), {"t_index": 0, "phi_index": 0}, "kaon", "unpolarized",
            method_a_event_multipliers=multipliers,
        )
        return templates

    def test_omitted_and_all_one_multiplier_preserve_baseline(self):
        baseline = self._fill()
        ones = self._fill({("prompt", 7): 1.0})
        self.assertEqual(baseline["mm"].fills, ones["mm"].fills)
        self.assertEqual(baseline["mm_nosub"].fills, ones["mm_nosub"].fills)

    def test_nontrivial_multiplier_scales_allcuts_and_nommcuts_only(self):
        templates = self._fill({("prompt", 7): 1.5})
        self.assertEqual(templates["mm"].fills[0][-1], 9.0)
        self.assertEqual(templates["mm_nosub"].fills[0][-1], 9.0)

    def test_multiplier_lookup_miss_and_invalid_value_are_errors_not_c_one(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            self._fill({})
        with self.assertRaisesRegex(ValueError, "invalid"):
            self._fill({("prompt", 7): math.nan})

    def test_template_extension_neither_normalizes_children_nor_reads_method_b(self):
        text = inspect.getsource(pion.fill_simc_shape_pion_subtraction_templates)
        self.assertNotIn("Scale(", text)
        self.assertNotIn("method_b", text.lower())


class FullParallelBranchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calculate_yield = _load_calculate_yield_module()

    @staticmethod
    def _clone(histogram, *_args, **kwargs):
        clone = histogram.Clone(kwargs.get("name", "parallel_clone"))
        if kwargs.get("reset"):
            clone.Reset()
        return clone

    def _branch_fixture(self, factor=1.0, *, source_coefficient=2.0, base_coefficient=2.0, accepted_coefficient=2.0):
        calc = self.calculate_yield
        pion_input = _YieldHistogram((20.0,), edges=(1.0, 1.3), errors=(0.0,))
        baseline_template = _YieldHistogram((6.0,), edges=(1.0, 1.3), errors=(0.0,))
        baseline_final = _YieldHistogram((14.0,), edges=(1.0, 1.3), errors=(0.0,))
        dummy = _YieldHistogram((5.0,), edges=(1.0, 1.3), errors=(0.0,))
        payload = {
            "accepted": True,
            "H_pion_control_model": object(),
            "weights": object(),
            "H_pion_subtraction_template_MM": baseline_template,
            "H_pion_subtraction_template_MM_nosub": baseline_template.Clone("baseline_nosub"),
            "H_MM_before_pion_subtraction": pion_input,
            "H_MM_after_pion_subtraction": baseline_final,
        }
        entry = {
            "child_valid": True,
            "particle_subtraction_component_payload": payload,
            "H_MM_DUMMY_NORM": dummy,
            "bg_fit1_frac_err": 0.0,
            "bg_fit2_frac_err": 0.0,
        }
        cache = _cache()
        spec = {
            "label": "prompt", "cache_section": cache,
            "coefficient": source_coefficient,
            "base_coefficient": base_coefficient,
        }
        accepted = _accepted_rows(_row(coefficient=accepted_coefficient))
        patches = (
            mock.patch.object(calc, "OUTPATH", "unused"),
            mock.patch.object(calc, "get_active_bg_profile_name", return_value="no_empirical_residual"),
            mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="simc_shape_components"),
            mock.patch.object(calc, "resolve_pion_subtraction_scope", return_value="t_bin"),
            mock.patch.object(calc, "get_particle_subtraction_setting_key", return_value="Q4p4W2p74"),
            mock.patch.object(calc, "accepted_f6_3_artifact_paths", return_value={}),
            mock.patch.object(calc, "load_accepted_f6_3_authority", return_value=([{}], {}, {}, {}, {"Left-lowe": "a", "f2": "a2", "f3": "b", "f4": "c"})),
            mock.patch.object(calc, "reconstruct_transient_factor_map", return_value=({("prompt", 7): factor}, {"accepted": True}, accepted)),
            mock.patch.object(calc, "iter_component_control_source_specs", return_value=[spec]),
            mock.patch.object(calc, "simc_shape_pion_weight_from_value", return_value=3.0),
            mock.patch.object(
                calc,
                "_component_cache_event_coefficient",
                side_effect=pion._component_cache_event_coefficient,
            ),
            mock.patch.object(calc, "clone_root_histogram", side_effect=self._clone),
        )
        def fill(templates, _specs, _reference, _weights, _context, _particle, _pol, *, method_a_event_multipliers=None):
            multiplier = method_a_event_multipliers[("prompt", 7)]
            templates["mm"].Fill(1.12, 6.0 * multiplier)
            templates["mm_nosub"].Fill(1.12, 6.0 * multiplier)
        patches += (mock.patch.object(calc, "fill_simc_shape_pion_subtraction_templates", side_effect=fill),)
        return {
            "calc": calc, "entry": entry, "cache": {"prompt": cache}, "patches": patches,
            "hist": {"phi_setting": "Left"},
            "processed": {"t_bin1phi_bin1": entry},
            "inp": {
                "ParticleType": "kaon", "EPSSET": "lowe", "POL": "unpolarized",
                "mm_min": 1.08, "mm_max": 1.18,
                "data_charge_err_left": 0.0, "dummy_charge_err_left": 0.0,
            },
            "measurements": {(0, 0): {"yield": 14.0, "statistical_error": 0.0, "total_error": 0.0}},
            "baseline_objects": (pion_input, baseline_template, baseline_final),
        }

    def _build(self, factor=1.0, *, fixture_kwargs=None, extra_patches=()):
        fixture = self._branch_fixture(factor, **(fixture_kwargs or {}))
        with ExitStack() as stack:
            for patcher in (*fixture["patches"], *extra_patches):
                stack.enter_context(patcher)
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        return source, fixture

    def test_all_one_parallel_branch_reproduces_baseline_without_mutating_public_objects(self):
        source, fixture = self._build(1.0)
        self.assertTrue(source["available"])
        child = source["children"][0]
        self.assertEqual(child["B_pi_A"].contents, child["B_pi_0"].contents)
        self.assertEqual(child["MM_A"].contents, child["MM_0"].contents)
        self.assertEqual(child["YA"], child["Y0"])
        self.assertEqual(child["YA_statistical_error"], child["Y0_statistical_error"])
        self.assertEqual(child["YA_total_error"], child["Y0_total_error"])
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])
        self.assertEqual(source["lambda_integration_window"], [1.08, 1.18])
        self.assertNotIn("transient_factor_map", str(source))

    def test_nontrivial_parallel_branch_changes_only_the_private_pion_subtraction(self):
        source, fixture = self._build(1.5)
        child = source["children"][0]
        self.assertEqual(child["pion_input"].contents, [20.0])
        self.assertEqual(child["B_pi_0"].contents, [6.0])
        self.assertEqual(child["B_pi_A"].contents, [9.0])
        self.assertEqual(child["MM_0"].contents, [14.0])
        self.assertEqual(child["MM_A"].contents, [11.0])
        self.assertEqual(child["YA"], 11.0)
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])

    def test_effective_filler_coefficient_not_raw_cache_field_controls_live_parity(self):
        source, fixture = self._build(
            fixture_kwargs={
                "source_coefficient": 4.0,
                "base_coefficient": 2.0,
                "accepted_coefficient": 2.0,
            },
        )
        self.assertFalse(source["available"])
        self.assertIn("signed_source_coefficient", source["reason"])
        self.assertEqual(
            pion._component_cache_event_coefficient(
                {"cache_section": _cache(), "coefficient": 4.0, "base_coefficient": 2.0}, 0,
            ),
            4.0,
        )
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])

    def test_branch_runtime_failures_are_unavailable_without_partial_children(self):
        failure_patches = (
            lambda fixture: mock.patch.object(
                fixture["calc"], "clone_root_histogram", side_effect=RuntimeError("clone failure"),
            ),
            lambda fixture: mock.patch.object(
                fixture["calc"], "fill_simc_shape_pion_subtraction_templates", side_effect=RuntimeError("fill failure"),
            ),
            lambda _fixture: mock.patch.object(
                _YieldHistogram, "Add", side_effect=RuntimeError("add failure"),
            ),
        )
        for make_patch in failure_patches:
            with self.subTest(failure=make_patch.__code__.co_firstlineno):
                fixture = self._branch_fixture()
                original = [list(item.contents) for item in fixture["baseline_objects"]]
                with ExitStack() as stack:
                    for patcher in (*fixture["patches"], make_patch(fixture)):
                        stack.enter_context(patcher)
                    source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                        fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                        fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
                    )
                self.assertFalse(source["available"])
                self.assertIn("f6_3_branch_exception:RuntimeError", source["reason"])
                self.assertNotIn("children", source)
                self.assertEqual([item.contents for item in fixture["baseline_objects"]], original)

    def test_unsupported_profile_or_mode_leaves_baseline_untouched_and_is_unavailable(self):
        fixture = self._branch_fixture()
        original = [list(item.contents) for item in fixture["baseline_objects"]]
        with mock.patch.object(fixture["calc"], "get_active_bg_profile_name", return_value="legacy"):
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        self.assertFalse(source["available"])
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], original)
        with mock.patch.object(fixture["calc"], "get_active_bg_profile_name", return_value="no_empirical_residual"), \
             mock.patch.object(fixture["calc"], "resolve_particle_subtraction_mode", return_value="legacy"):
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        self.assertFalse(source["available"])

    def test_parallel_builder_uses_existing_cache_only_and_has_no_method_b_dependency(self):
        text = inspect.getsource(self.calculate_yield._build_f6_3_parallel_method_a_source).lower()
        self.assertNotIn("bin_data(", text)
        self.assertNotIn("pion_hgcer_method_b", text)
        self.assertNotIn("bg_fit(", text)

    @staticmethod
    def _snapshot_public_groups(groups):
        """Return a JSON-serializable complete snapshot of public yield groups."""
        child_keys = sorted((int(t_index), int(phi_index)) for t_index, phi_index in groups)
        return {
            "container_type": type(groups).__name__,
            "top_level_child_keys": [list(key) for key in child_keys],
            "children": {
                "{},{}".format(t_index, phi_index): {
                    "field_keys": sorted(str(field) for field in groups[(t_index, phi_index)]),
                    "values": {
                        str(field): float(value)
                        for field, value in groups[(t_index, phi_index)].items()
                    },
                }
                for t_index, phi_index in child_keys
            },
        }

    def _public_baseline_run(self, *, branch_runtime_failure=False):
        calc = self.calculate_yield
        fixture = _public_yield_fixture(calc)
        fixture["binned_dict"]["yield"] = fixture["binned_dict"].pop("kaon")
        component_payload = {
            "accepted": True,
            "fixture_schema": "f6_3_component_payload_sentinel/v1",
            "component_roles": ["pion", {"source": "prompt", "coefficient": -0.25}],
            "provenance": {"revision": 7, "labels": ["baseline", "preserve"]},
        }
        fixture["processed_entry"]["particle_subtraction_component_payload"] = component_payload
        component_payload_snapshot = copy.deepcopy(component_payload)
        component_payload_json_snapshot = json.dumps(
            component_payload, sort_keys=True, separators=(",", ":"),
        )
        integration = mock.Mock(side_effect=((10.0, 0.3), (5.0, 0.4)))
        common = (
            mock.patch.object(calc, "_resolve_bg_opt_prepass_data_cache", return_value=None),
            mock.patch.object(calc, "_get_cached_kaon_signal_shape_payload", return_value=None),
            mock.patch.object(calc, "_get_cached_kaon_sigma0_shape_payload", return_value=None),
            mock.patch.object(calc, "bin_data", return_value=(fixture["binned_dict"], None, None)),
            mock.patch.object(calc, "get_active_bg_profile_name", return_value="no_empirical_residual"),
            mock.patch.object(calc, "integral_with_stat_error", integration),
        )
        if branch_runtime_failure:
            branch = (
                mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="simc_shape_components"),
                mock.patch.object(calc, "resolve_pion_subtraction_scope", return_value="t_bin"),
                mock.patch.object(calc, "accepted_f6_3_artifact_paths", side_effect=RuntimeError("branch authority failure")),
            )
        else:
            branch = (mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="legacy"),)
        with ExitStack() as stack:
            for patcher in (*common, *branch):
                stack.enter_context(patcher)
            groups = calc.calculate_yield_data(
                "yield", fixture["hist"], [0.0, 1.0], [-180.0, 180.0], fixture["inp_dict"],
            )
        child = fixture["hist"]["_e8_2_baseline_stage_source"]["children"][0]
        return {
            "groups": self._snapshot_public_groups(groups),
            "stage_window_yields": dict(fixture["processed_entry"]["stage_window_yields"]),
            "scale_factor": fixture["binned_dict"]["yield"]["scale_factor"],
            "final_hist": list(fixture["final_hist"].contents),
            "final_errors": list(fixture["final_hist"].errors),
            "e8_2_yield": (
                child["final_yield"], child["statistical_error"], child["total_error"],
            ),
            "component_payload": fixture["processed_entry"]["particle_subtraction_component_payload"],
            "component_payload_original": component_payload,
            "component_payload_snapshot": component_payload_snapshot,
            "component_payload_json_snapshot": component_payload_json_snapshot,
            "f6_3": fixture["hist"]["_f6_3_parallel_method_a_source"],
        }

    def test_public_calculate_yield_baseline_outputs_survive_unavailable_and_runtime_failed_f6_3(self):
        unavailable = self._public_baseline_run()
        failed = self._public_baseline_run(branch_runtime_failure=True)
        self.assertEqual(unavailable["groups"], {
            "container_type": "defaultdict",
            "top_level_child_keys": [[0, 0]],
            "children": {
                "0,0": {
                    "field_keys": ["yield", "yield_err"],
                    "values": {"yield": 10.0, "yield_err": 0.3905124837953327},
                },
            },
        })
        self.assertEqual(unavailable["groups"], failed["groups"])
        self.assertEqual(unavailable["stage_window_yields"], {
            "raw_prompt": 21.0,
            "after_random_subtraction": 18.0,
            "after_dummy_subtraction": 16.0,
            "after_pion_subtraction": 10.0,
        })
        self.assertEqual(unavailable["stage_window_yields"], failed["stage_window_yields"])
        self.assertEqual(unavailable["scale_factor"], [[0.0]])
        self.assertEqual(unavailable["scale_factor"], failed["scale_factor"])
        self.assertEqual(unavailable["final_hist"], [10.0])
        self.assertEqual(unavailable["final_hist"], failed["final_hist"])
        self.assertEqual(unavailable["final_errors"], [0.3])
        self.assertEqual(unavailable["final_errors"], failed["final_errors"])
        self.assertEqual(unavailable["e8_2_yield"], (10.0, 0.3, 0.3905124837953327))
        self.assertEqual(unavailable["e8_2_yield"], failed["e8_2_yield"])
        for result in (unavailable, failed):
            self.assertIs(result["component_payload"], result["component_payload_original"])
            self.assertEqual(result["component_payload"], result["component_payload_snapshot"])
            self.assertEqual(
                json.dumps(result["component_payload"], sort_keys=True, separators=(",", ":")),
                result["component_payload_json_snapshot"],
            )
        self.assertFalse(unavailable["f6_3"]["available"])
        self.assertFalse(failed["f6_3"]["available"])
        self.assertIn("f6_3_branch_exception:RuntimeError", failed["f6_3"]["reason"])


class BranchBoundaryTests(unittest.TestCase):
    def test_unavailable_source_preserves_baseline_and_all_forbidden_boundaries(self):
        source = f63.unavailable_parallel_source("authority_missing", setting_id="Left-lowe")
        self.assertFalse(source["available"])
        self.assertTrue(source["baseline_public_output_unchanged"])
        self.assertFalse(source["production_promotion_performed"])
        self.assertFalse(source["method_b_numerical_dependency"])
        self.assertFalse(source["empirical_residual_used"])
        self.assertFalse(source["event_correction_persisted"])

    def test_helper_is_pure_python_and_does_not_traverse_trees(self):
        text = (ROOT / "src" / "cuts" / "pion_hgcer_method_a_parallel_full_procedure.py").read_text(encoding="utf-8")
        self.assertNotIn("import ROOT", text)
        self.assertNotIn("GetEntries", text)


if __name__ == "__main__":
    unittest.main()
