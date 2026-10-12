"""Deterministic synthetic authority/alignment/arithmetic/owner tests, no farm.

Fixture pin patches affect ONLY the new detached diagnostic, never live F.6.3
or production authority. The unpatched entry point rejects these fixtures.
"""

from contextlib import ExitStack, redirect_stdout, redirect_stderr
import copy
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/"src"/"cuts"), str(ROOT/"testing")]
import pion_hgcer_method_a_reweighted_contribution_debug as debug
import analyze_pion_hgcer_method_a_reweighted_contribution_debug as cli
import test_pion_hgcer_method_a_current_lineage_diagnostic as fixtures
import pion_hgcer_method_a_parallel_full_procedure as lineage


class NeffMomentParityTests(unittest.TestCase):
    def test_one_event_multiplication_and_power_are_discriminating(self):
        supplied = 55.83414573673653
        s = debug.summary([supplied], [1.0], [supplied])
        debug.validate_summary(s)
        self.assertEqual(s["neff"]["w0"]["value"], 1.0)
        # The supplied sample is not discriminating on Windows CPython 3.12.10.
        # Verified alternate: multiply=2917.54439736457, power=2917.5443973645697.
        x = supplied if supplied * supplied != supplied**2 else 54.014298082679645
        self.assertNotEqual(x * x, x**2)
        s = debug.summary([x], [1.0], [x])
        debug.validate_summary(s)
        for name in ("w0", "w_rw", "abs_b", "abs_b_rw"):
            self.assertEqual(s["neff"][name]["value"], 1.0)
            self.assertEqual(s["neff"][name]["numerator"], x * x)

    def test_positive_reweighted_support(self):
        s = debug.summary([1.0, 2.0, 3.0], [2.0, 0.5, 3.0], [2.0, 4.0, 1.0])
        debug.validate_summary(s)
        self.assertEqual(s["neff"]["w_rw"], debug.effective([2.0, 1.0, 9.0]))
        self.assertEqual(s["neff"]["abs_b_rw"], debug.effective([4.0, 2.0, 3.0]))

    def test_zero_support_and_signed_cancellation(self):
        for weights in ([], [0.0, 0.0]):
            with self.subTest(weights=weights):
                s = debug.summary(weights, [1.0] * len(weights), weights)
                debug.validate_summary(s)
                for measure in s["neff"].values():
                    self.assertFalse(measure["available"])
                    self.assertIsNone(measure["value"])
                    self.assertEqual(measure["reason"], "zero_support")
        s = debug.summary([1.0, 1.0], [2.0, 2.0], [1.0, -1.0])
        debug.validate_summary(s)
        self.assertEqual((s["B"], s["B_rw"], s["negative"]), (0.0, 0.0, -1.0))
        self.assertFalse(s["signed_relative_delta"]["available"])
        for measure in s["neff"].values():
            self.assertEqual(measure["value"], 2.0)

    def test_one_ulp_metric_tampering_is_rejected_exactly(self):
        s = debug.summary([1.0, 2.0], [2.0, 3.0], [1.0, -2.0])
        debug.validate_summary(s)
        for name in s["neff"]:
            for field in ("value", "numerator", "denominator"):
                changed = copy.deepcopy(s)
                changed["neff"][name][field] = math.nextafter(s["neff"][name][field], math.inf)
                with self.subTest(name=name, field=field), self.assertRaisesRegex(
                        debug.ContributionDebugError, "^neff_moment_mismatch$"):
                    debug.validate_summary(changed)

    def test_negative_and_nonfinite_moments_remain_rejected(self):
        s = debug.summary([1.0], [1.0], [1.0])
        for field in ("sum_b2", "sum_b_rw2", "sum_w0", "sum_w_rw", "sum_w02", "sum_w_rw2"):
            changed = copy.deepcopy(s)
            changed[field] = -1.0
            with self.subTest(field=field), self.assertRaisesRegex(
                    debug.ContributionDebugError, "^nonnegative_moment_invalid$"):
                debug.validate_summary(changed)
            for value in (math.nan, math.inf, -math.inf):
                changed[field] = value
                # Persisted diagnostics guard every moment before validation.
                with self.subTest(field=field, value=value), self.assertRaisesRegex(
                        debug.ContributionDebugError, "^persisted_nonfinite$"):
                    debug.seal(changed)
        for value in (-1.0, math.nan, math.inf):
            with self.subTest(value=value), self.assertRaisesRegex(
                    debug.ContributionDebugError, "^neff_measure_invalid$"):
                debug.effective([value])

    def test_neff_exceeding_count_is_rejected_even_with_matching_metric(self):
        s = debug.summary([1.0], [1.0], [1.0])
        s["sum_w02"] = 0.25
        s["neff"]["w0"] = {"count": 1, **debug.metric(1.0, 0.25)}
        with self.assertRaisesRegex(debug.ContributionDebugError, "^neff_exceeds_count$"):
            debug.validate_summary(s)


class ContributionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixtures.DiagnosticTests.setUpClass()
        cls.f = fixtures.DiagnosticTests
        cls.pins = {
            "F1_HASHES": cls.f.hashes, "F3_HASH": cls.f.map_sha,
            "F3_MAP": cls.f.map["acceptance_map"]["fingerprint"],
            "F3_ARTIFACT": cls.f.map["artifact_fingerprint"],
            "F3_ALGORITHM": cls.f.map["acceptance_map"]["algorithm_fingerprint"],
            "F4_HASH": cls.f.correction_sha,
            "F4_CORRECTION": cls.f.correction["correction"]["fingerprint"],
            "F4_ARTIFACT": cls.f.correction["artifact_fingerprint"],
            "RECORDED_PARENTS": tuple((p["application_event_count"], p["baseline_parent_sum"], p["parent_normalization"]) for p in cls.f.correction["correction"]["parents"] if p["setting_id"] == "Left-lowe"),
        }
        cls.live_snapshot = copy.deepcopy((lineage.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256,
            lineage.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC,
            lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC))
        with cls.authority():
            cls.result = cls.build()

    @classmethod
    def authority(cls):
        stack = ExitStack()
        for key, value in cls.pins.items():
            stack.enter_context(mock.patch.object(debug, key, copy.deepcopy(value)))
        return stack

    @classmethod
    def build(cls, **overrides):
        values = {"f1_artifacts": cls.f.inputs, "f3_artifact": cls.f.map,
                  "f4_artifact": cls.f.correction,
                  "hashes": {**cls.f.hashes, "f3": cls.f.map_sha, "f4": cls.f.correction_sha},
                  "provenance": {"git_head": "0"*40, "git_status_short": "", "input_basenames": {}}}
        values.update(overrides)
        return debug.build_diagnostic(**values)

    def test_exact_unpatched_historical_pins(self):
        self.assertEqual(debug.F1_HASHES["Left-lowe"], "10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95")
        self.assertEqual(debug.F4_HASH, "1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902")
        self.assertEqual(debug.F4_CORRECTION, "bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368")
        self.assertEqual(debug.REFERENCE_HASH, "ecd547de15aa19a50b73bfd8d9dd1c159027c6438ca38b79bc45ff69a8d0e8c7")
        self.assertNotEqual(debug.F4_HASH, lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC[debug.KINEMATIC]["source_file_sha256"])
        with self.assertRaisesRegex(debug.ContributionDebugError, "raw_sha256"):
            self.build()

    def test_current_candidate_hashes_rejected_before_builder(self):
        hashes = {**lineage.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256, "f3": lineage.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC[debug.KINEMATIC]["source_file_sha256"], "f4": lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC[debug.KINEMATIC]["source_file_sha256"]}
        with mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data") as builder, self.assertRaises(debug.ContributionDebugError):
            self.build(hashes=hashes)
        builder.assert_not_called()

    def test_all_seven_hashes_individually_fail(self):
        with self.authority():
            for role in self.result["input_file_hashes"]:
                hashes = {**self.result["input_file_hashes"], role: "a"*64}
                with self.subTest(role=role), self.assertRaisesRegex(debug.ContributionDebugError, role):
                    self.build(hashes=hashes)

    def test_determinism_input_objects_and_live_pins_preserved(self):
        original = copy.deepcopy((self.f.inputs, self.f.map, self.f.correction))
        with self.authority():
            self.assertEqual(self.build(), self.result)
            self.assertTrue(debug.validate_diagnostic(self.result))
        self.assertEqual((self.f.inputs, self.f.map, self.f.correction), original)
        self.assertEqual((lineage.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256,
            lineage.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC,
            lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC), self.live_snapshot)

    def test_intrinsic_f3_and_f4_fingerprints(self):
        for label, artifact, path in (("f3_artifact", self.f.map, ("artifact_fingerprint",)), ("f3_map", self.f.map, ("acceptance_map", "fingerprint")), ("f3_algorithm", self.f.map, ("acceptance_map", "algorithm_fingerprint")), ("f4_artifact", self.f.correction, ("artifact_fingerprint",)), ("f4_correction", self.f.correction, ("correction", "fingerprint"))):
            altered = copy.deepcopy(artifact); node = altered
            for key in path[:-1]: node = node[key]
            node[path[-1]] = "a"*64
            with self.subTest(label=label), self.authority(), self.assertRaises(debug.ContributionDebugError):
                self.build(**{"f3_artifact" if label.startswith("f3") else "f4_artifact": altered})

    def test_embedded_original_f3_authority_must_be_proven(self):
        altered = copy.deepcopy(self.f.correction)
        altered["correction"]["f3_runtime_authority"]["accepted"]["source_file_sha256"] = "a"*64
        with self.authority(), mock.patch.object(debug.f5, "_validate_f4_artifact", return_value=(altered, altered["correction"], {})), self.assertRaisesRegex(debug.ContributionDebugError, "f3_runtime_authority"):
            self.build(f4_artifact=altered)

    def test_exact_reproduction_required(self):
        shared = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        def changed(*args, **kwargs):
            core, reviews = shared(*args, **kwargs)
            core["parents"][0]["baseline_parent_sum"] = math.nextafter(core["parents"][0]["baseline_parent_sum"], math.inf)
            return core, reviews
        with self.authority(), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=changed), mock.patch.object(debug, "measure_parent") as measure, self.assertRaisesRegex(debug.ContributionDebugError, "exact_reproduction"):
            self.build()
        measure.assert_not_called()

    def test_duplicate_missing_f1_identity_and_contract_fingerprint(self):
        for action in ("duplicate", "missing", "fingerprint"):
            altered = copy.deepcopy(self.f.inputs)
            if action == "duplicate": altered[0]["contract"]["application_records"].append(copy.deepcopy(altered[0]["contract"]["application_records"][0]))
            elif action == "missing": altered[0]["contract"]["application_records"].pop()
            else: altered[0]["contract"]["fingerprint"] = "a"*64
            with self.subTest(action=action), self.authority(), self.assertRaises(debug.ContributionDebugError):
                self.build(f1_artifacts=altered)

    def test_missing_duplicate_parent_or_factor_shape(self):
        shared = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        for action in ("missing", "duplicate", "short", "nonfinite"):
            def changed(*args, **kwargs):
                core, reviews = shared(*args, **kwargs)
                if action == "missing": reviews.pop()
                elif action == "duplicate": reviews.append(reviews[0])
                elif action == "short": reviews[0]["correction_factors"] = reviews[0]["correction_factors"][:-1]
                else: reviews[0]["correction_factors"][0] = float("nan")
                return core, reviews
            with self.subTest(action=action), self.authority(), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=changed), self.assertRaises(debug.ContributionDebugError):
                self.build()

    def measurement_inputs(self):
        parsed = debug.f3._validate_f1_artifacts(self.f.inputs)
        rows = debug.f4._raw_application_rows(self.f.inputs, parsed, 1e-12)[("Left-lowe", 0)]
        raw = {(r["source_label"], r["entry_index"]): copy.deepcopy(r) for r in self.f.inputs[0]["contract"]["application_records"]}
        _, reviews = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            self.f.inputs, self.f.map, f1_input_file_hashes=self.f.hashes,
            f3_input_file_sha256=self.f.map_sha,
            accepted_f3_runtime_authority_by_kinematic=self.f.f3_authority)
        return rows, raw, reviews[0], self.f.correction["correction"]["parents"][0]

    def test_original_mm_join_source_sign_and_factor_order(self):
        rows, raw, review, parent = self.measurement_inputs()
        for action in ("missing_id", "duplicate_id", "assignment", "sign", "nonfinite_mm", "reordered_factors"):
            changed_rows, changed_raw, changed_review = copy.deepcopy((rows, raw, review))
            key = (rows[0]["source_label"], rows[0]["entry_index"])
            if action == "missing_id": changed_raw.pop(key)
            elif action == "duplicate_id": changed_rows[1] = copy.deepcopy(changed_rows[0])
            elif action == "assignment": changed_raw[key]["t_index"] = 2
            elif action == "sign": changed_rows[0]["signed_source_coefficient"] *= -1
            elif action == "nonfinite_mm": changed_raw[key]["analysis_MM"] = float("nan")
            else: changed_review["correction_factors"] = np.roll(changed_review["correction_factors"], 1)
            with self.subTest(action=action), self.assertRaises(debug.ContributionDebugError):
                debug.measure_parent(changed_rows, changed_raw, changed_review, parent, list(debug.f5._PHI_EDGES))

    def test_all_mm_flows_and_boundaries_accumulate_in_joint_cells(self):
        rows, raw, review, parent = self.measurement_inputs()
        coordinates = [.79, .8, debug.MM_EDGES[1], 1.4, 1.41]
        expected = [0]*42
        for i, row in enumerate(rows):
            value = coordinates[i % 5]
            raw[(row["source_label"], row["entry_index"])]["analysis_MM"] = value
            expected[debug.bin_index(value, debug.MM_EDGES)] += 1
        result = debug.measure_parent(rows, raw, review, parent, list(debug.f5._PHI_EDGES))
        self.assertEqual([s["count"] for s in result["mm"]], expected)
        for m in range(42):
            self.assertEqual(sum(row[m]["count"] for row in result["joint"]), expected[m])
        self.assertGreater(expected[0], 0); self.assertGreater(expected[-1], 0)

    def test_shared_producer_input_mutation_is_a_hard_failure(self):
        inputs = copy.deepcopy(self.f.inputs)
        shared = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        def mutating(*args, **kwargs):
            result = shared(*args, **kwargs)
            args[0][0]["unapproved_metadata"] = True
            return result
        with self.authority(), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=mutating), self.assertRaisesRegex(debug.ContributionDebugError, "objects_changed"):
            self.build(f1_artifacts=inputs)

    def test_persisted_scope_authority_and_check_flags_are_strict(self):
        for action in ("scope", "authority", "f3", "check", "reference"):
            changed = copy.deepcopy(self.result); changed.pop("fingerprint")
            if action == "scope": changed["scope"]["setting_id"] = "Center-lowe"
            elif action == "authority": changed["input_authority"]["accepted"]["source_file_sha256"] = "a"*64
            elif action == "f3": changed["f3_map_fingerprint"] = "a"*64
            elif action == "check": changed["checks"]["parent_closure"] = False
            else: changed["historical_evidence"]["sha256"] = "a"*64
            with self.subTest(action=action), self.authority(), self.assertRaises(debug.ContributionDebugError):
                debug.validate_diagnostic(debug.seal(changed))

    def test_optional_hash_pinned_reference_checks_exact_1d_aggregates(self):
        parsed = debug.f3._validate_f1_artifacts(self.f.inputs)
        rows = debug.f4._raw_application_rows(self.f.inputs, parsed, 1e-12)
        raw = {(r["source_label"], r["entry_index"]): r for r in self.f.inputs[0]["contract"]["application_records"]}
        _, reviews = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            self.f.inputs, self.f.map, f1_input_file_hashes=self.f.hashes,
            f3_input_file_sha256=self.f.map_sha, accepted_f3_runtime_authority_by_kinematic=self.f.f3_authority)
        parents = []
        for t in range(3):
            parents.append({"canonical_t_index": t, **debug.legacy._sensitivity_parent(rows[("Left-lowe", t)], raw, reviews[t], self.f.correction["correction"]["parents"][t], list(debug.f5._PHI_EDGES))})
        diagnostic = {"input_file_hashes": {"f1": self.f.hashes, "f3": self.f.map_sha, "f4": self.f.correction_sha}, "measurement_b": {"parents": parents}}
        diagnostic["fingerprint"] = debug.fingerprint(diagnostic)
        reference = debug.legacy.build_artifact(diagnostic, provenance={})
        digest = fixtures.digest(reference)
        with self.authority(), mock.patch.multiple(debug, REFERENCE_HASH=digest, REFERENCE_FP=diagnostic["fingerprint"], REFERENCE_WRAPPER=reference["artifact_fingerprint"]):
            result = self.build(reference=reference, reference_sha=digest)
            self.assertTrue(result["historical_evidence"]["reference_json_checked"])
            changed = copy.deepcopy(reference)
            changed["diagnostic"]["measurement_b"]["parents"][0]["phi"]["bins"][0]["count"] += 1
            with self.assertRaisesRegex(debug.ContributionDebugError, "wrapper_mismatch"):
                self.build(reference=changed, reference_sha=digest)

    def test_signed_absolute_cancellation_and_squared_support(self):
        s = debug.summary([10, 9, 0], [2/11, 1/11, 1000], [10, -9, 0])
        self.assertEqual(s["B"], 1); self.assertEqual(s["A"], 19)
        self.assertEqual(s["positive"], 10); self.assertEqual(s["negative"], -9)
        self.assertAlmostEqual(s["B_rw"], 1); self.assertEqual(s["sum_b2"], 181)
        self.assertEqual(s["zero_weight_count"], 1)
        self.assertAlmostEqual(s["neff"]["abs_b"]["value"], 361/181)
        s = debug.summary([1, 1], [1, 1], [1, -1])
        self.assertFalse(s["signed_relative_delta"]["available"])
        self.assertIsNone(s["signed_relative_delta"]["value"])

    def test_neff_uniform_concentrated_zero_and_empty(self):
        self.assertEqual(debug.effective([1]*10)["value"], 10)
        self.assertEqual(debug.effective([0, 0, 7])["value"], 1)
        for weights in ([], [0, 0]):
            result = debug.effective(weights)
            self.assertFalse(result["available"]); self.assertIsNone(result["value"])
        with self.assertRaises(debug.ContributionDebugError): debug.effective([-1])

    def test_zero_huge_and_invalid_weights_and_products(self):
        s = debug.summary([0, 1e100], [1e100, 1e20], [0, -1e100])
        self.assertTrue(math.isfinite(s["sum_b_rw2"]))
        for w, c, b in (([1], [0], [1]), ([-1], [1], [-1]), ([1], [float("inf")], [1]), ([float("nan")], [1], [1]), ([1e200], [1], [1e200]), ([1], [], [1])):
            with self.subTest(w=w, c=c), self.assertRaises((debug.ContributionDebugError, OverflowError)):
                debug.summary(w, c, b)

    def test_tail_actual_contribution_rank_and_ties(self):
        rows = [{"source_label": "prompt", "entry_index": i} for i in range(100)]
        w = [0., 100.] + [1.]*98; c = [1000., 1.] + [1.]*98; b = list(w)
        t = debug.tails(rows, w, c, b)
        self.assertEqual(t["C"][0]["selected_count"], 1)
        self.assertEqual(t["C"][0]["reweighted_abs_share"]["value"], 0)
        self.assertEqual(t["abs_b_rw"][0]["support"]["A_rw"], 100)
        self.assertEqual(t["C"][1]["selected_count"], 5)
        self.assertEqual(t, debug.tails(rows, w, c, b))
        self.assertNotIn("entry_index", json.dumps(t))
        self.assertEqual(debug.tails([], [], [], [])["C"][0]["selected_count"], 0)

    def test_mm_boundaries_and_flows_last_edge_is_inclusive(self):
        edges = debug.MM_EDGES
        values = [edges[0]-.1, edges[0], edges[1], edges[-1], edges[-1]+.1]
        self.assertEqual([debug.bin_index(v, edges) for v in values], [0, 1, 2, 40, 41])
        self.assertEqual(debug.bin_index(180, debug.f5._PHI_EDGES)-1, 8)
        self.assertEqual(debug.bin_index(-140, debug.f5._PHI_EDGES)-1, 1)
        with self.assertRaises(debug.ContributionDebugError): debug.bin_index(float("nan"), edges)

    def test_joint_3_by_9_by_42_and_all_moments_reconstruct(self):
        with self.authority():
            debug.validate_diagnostic(self.result)
        self.assertEqual(sum(len(row) for p in self.result["parents"] for row in p["joint"]), 1134)
        for p in self.result["parents"]:
            self.assertEqual(sum(c["count"] for row in p["joint"] for c in row), p["totals"]["count"])
            for field in debug.SUM_FIELDS:
                self.assertAlmostEqual(math.fsum(c[field] for row in p["joint"] for c in row), p["totals"][field])
            self.assertEqual(set(p["sources"]), set(debug.SOURCES))

    def test_signed_and_marginal_tamper_rejected_even_after_reseal(self):
        for action in ("count", "moment", "neff", "phi", "source", "tail"):
            changed = copy.deepcopy(self.result)
            p = changed["parents"][0]
            if action == "count": p["joint"][0][0]["count"] += 1
            elif action == "moment": p["joint"][0][0]["sum_b2"] += 1
            elif action == "neff": p["totals"]["neff"]["abs_b"]["value"] += 1
            elif action == "phi": p["phi"].pop()
            elif action == "source": p["sources"].pop("rand")
            else: p["tails"]["C"][0]["reweighted_abs_share"]["value"] = 0.999
            changed.pop("fingerprint")
            with self.subTest(action=action), self.authority(), self.assertRaises(debug.ContributionDebugError):
                debug.validate_diagnostic(debug.seal(changed))

    def test_tampered_fingerprint_and_forbidden_nested_fields(self):
        changed = copy.deepcopy(self.result); changed["fingerprint"] = "a"*64
        with self.authority(), self.assertRaisesRegex(debug.ContributionDebugError, "fingerprint"):
            debug.validate_diagnostic(changed)
        for key in ("entry_index", "event_rows", "event_ids", "correction_factors", "analysis_MM", "weight_tuples", "top_events", "input_paths"):
            with self.subTest(key=key), self.assertRaises(debug.ContributionDebugError):
                debug.guard({"nested": [{key: [1, 2]}]})
        with self.assertRaises(debug.ContributionDebugError): debug.guard({"nested": np.ones(3)})

    def test_production_status_mutation_rejected(self):
        for key in ("production_objects_mutated", "method_b_numerical_dependency", "child_renormalization_performed", "correction_modified"):
            changed = copy.deepcopy(self.result); changed.pop("fingerprint"); changed[key] = True
            with self.authority(), self.assertRaises(debug.ContributionDebugError):
                debug.validate_diagnostic(debug.seal(changed))

    def write_locator(self, root):
        paths = {"f1": {}, "f3": "map.json", "f4": "correction.json"}
        for setting, artifact in zip(self.f.hashes, self.f.inputs):
            name = setting+".json"; (root/name).write_bytes(fixtures.encoded(artifact)); paths["f1"][setting] = name
        (root/"map.json").write_bytes(fixtures.encoded(self.f.map))
        (root/"correction.json").write_bytes(fixtures.encoded(self.f.correction))
        locator = root/"locator.json"; locator.write_bytes(cli.json_bytes(paths))
        return locator, paths

    def args(self, root, locator):
        return ["--inputs-manifest", str(locator), "--outdir", str(root/"artifacts"), "--bundle-dir", str(root/"bundles"), "--kinematic", debug.KINEMATIC, "--phi", "Left", "--epsilon", "lowe"]

    def prepare_root(self, root):
        (root/"artifacts").mkdir(); (root/"bundles").mkdir()
        return self.write_locator(root)

    def test_missing_duplicate_unknown_manifest_roles_and_paths(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, paths = self.write_locator(root)
            for action in ("missing_role", "extra", "missing_file", "duplicate_file"):
                changed = copy.deepcopy(paths)
                if action == "missing_role": changed["f1"].pop("Left-highe")
                elif action == "extra": changed["authority"] = "override"
                elif action == "missing_file": changed["f3"] = "absent.json"
                else: changed["f3"] = changed["f4"]
                locator.write_bytes(cli.json_bytes(changed))
                with self.subTest(action=action), self.assertRaises(debug.ContributionDebugError): cli.load_inputs(locator)
            locator.write_text('{"f1":{},"f1":{},"f3":"map.json","f4":"correction.json"}')
            with self.assertRaisesRegex(debug.ContributionDebugError, "duplicate_json_key"): cli.load_inputs(locator)

    def test_wrong_raw_bytes_rejected_before_scientific_builder(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.write_locator(root)
            (root/"map.json").write_text("{}")
            with mock.patch.object(debug, "build_diagnostic") as builder, self.assertRaisesRegex(debug.ContributionDebugError, "sha256_mismatch:f3"):
                cli.load_diagnostic(locator)
            builder.assert_not_called()

    def test_reference_optional_cannot_weaken_pins_and_wrong_reference_fails(self):
        with self.assertRaises(debug.ContributionDebugError): debug.validate_reference({}, "a"*64)
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.write_locator(root)
            (root/"reference.json").write_text("{}")
            with self.assertRaisesRegex(debug.ContributionDebugError, "reference_raw_sha256"):
                cli.load_diagnostic(locator, root/"reference.json")

    def test_preflight_writes_nothing_and_has_no_renderer(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.prepare_root(root)
            with mock.patch.object(cli, "write_review_pdf") as renderer, redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(self.args(root, locator)+["--preflight-only"]), 0)
            renderer.assert_not_called()
            self.assertEqual(list((root/"artifacts").iterdir()), []); self.assertEqual(list((root/"bundles").iterdir()), [])
            self.assertFalse(json.loads(stdout.getvalue())["outputs_written"])

    def test_owner_real_persisted_only_pdf_package_and_success_stdout(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.prepare_root(root)
            original_inputs = {p.name: p.read_bytes() for p in root.glob("*.json")}
            with redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(self.args(root, locator)), 0)
            path = Path(stdout.getvalue().strip())
            self.assertEqual(path.parent, root/"bundles")
            manifest = cli.verify_zip(path)
            self.assertEqual(manifest["pdf"]["pages"], 8)
            self.assertEqual(list((root/"artifacts").iterdir()), [])
            self.assertEqual({p.name: p.read_bytes() for p in root.glob("*.json")}, original_inputs)
            with zipfile.ZipFile(path) as archive:
                raw = archive.read("diagnostic.json"); pdf = archive.read("diagnostic.pdf")
                self.assertNotIn(str(root).encode(), raw)
                persisted = root/"persisted.json"; persisted.write_bytes(raw)
                repeat = root/"repeat.pdf"
                with mock.patch.object(debug, "build_diagnostic", side_effect=AssertionError("renderer must not compute")):
                    cli.write_review_pdf(repeat, persisted)
                self.assertEqual(pdf, repeat.read_bytes())
                members = {name: archive.read(name) for name in cli.MEMBERS}
            deterministic = root/"repeat.zip"; cli.create_zip(deterministic, members)
            self.assertEqual(path.read_bytes(), deterministic.read_bytes())
            self.assertEqual(cli.check_pdf(pdf), 8)

    def test_failed_pdf_gate_cleans_only_owned_outputs(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.prepare_root(root)
            unrelated = root/"artifacts"/"preserve.txt"; unrelated.write_text("keep")
            def bad_renderer(path, persisted):
                path.write_bytes(b"invalid"); return {"pages": 8}
            with mock.patch.object(cli, "write_review_pdf", side_effect=bad_renderer), redirect_stderr(io.StringIO()) as stderr:
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            self.assertIn("pdf_container_invalid", stderr.getvalue())
            self.assertEqual(list((root/"artifacts").iterdir()), [unrelated]); self.assertEqual(unrelated.read_text(), "keep")
            self.assertEqual(list((root/"bundles").iterdir()), [])

    def test_existing_output_collision_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.prepare_root(root)
            collision = root/"bundles"/"kaonlt-historical-contribution-fixed.zip"; collision.write_text("keep")
            with mock.patch.object(cli.uuid, "uuid4", return_value=mock.Mock(hex="fixed")), mock.patch.object(cli, "load_diagnostic") as load, redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            load.assert_not_called(); self.assertEqual(collision.read_text(), "keep")

    def test_cli_missing_inputs_and_checkout_output_no_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/"artifacts").mkdir(); (root/"bundles").mkdir()
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(self.args(root, root/"absent.json")), 1)
            self.assertEqual(stdout.getvalue(), "")
            with self.assertRaises(debug.ContributionDebugError): cli.roots(ROOT, root/"bundles")

    def test_pdf_missing_page_and_invalid_container(self):
        for raw in (b"", b"%PDF-1.4\n/Type /Page\n/Type /Pages /Count 1\n%%EOF"):
            with self.assertRaises(debug.ContributionDebugError): cli.check_pdf(raw)

    def test_zip_missing_member_and_modified_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root/"bad.zip"
            with zipfile.ZipFile(path, "w") as archive: archive.writestr("diagnostic.json", b"{}")
            with self.assertRaisesRegex(debug.ContributionDebugError, "inventory"): cli.verify_zip(path)
            manifest = debug.seal({"schema_version": "pion_hgcer_method_a_reweighted_contribution_debug_bundle/v1", "status": "non_authoritative", "members": {"diagnostic.json": {"bytes": 2, "sha256": "a"*64}, "diagnostic.pdf": {"bytes": 2, "sha256": "a"*64}}})
            members = {"diagnostic.json": b"{}", "diagnostic.pdf": b"{}", "manifest.json": cli.json_bytes(manifest)}
            path = root/"tampered.zip"; cli.create_zip(path, members)
            with self.assertRaisesRegex(debug.ContributionDebugError, "digest"): cli.verify_zip(path)
            with self.assertRaisesRegex(debug.ContributionDebugError, "bytes"): cli.verify_zip(path, {**members, "diagnostic.pdf": b"changed"})

    def test_final_zip_failure_removes_owner_created_bundle(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory); locator, _ = self.prepare_root(root)
            # Exercise publication cleanup without a repeated expensive render.
            def renderer(path, persisted):
                path.write_bytes(b"fake"); return {"pages": 8, "page_titles": list(cli.PAGE_TITLES), "renderer_input": "persisted aggregates only"}
            calls = []
            def verify(path, expected=None):
                calls.append(path)
                if path.parent == root/"bundles": raise debug.ContributionDebugError("final_zip_tampered")
            with mock.patch.object(cli, "write_review_pdf", side_effect=renderer), mock.patch.object(cli, "check_pdf", return_value=8), mock.patch.object(cli, "verify_zip", side_effect=verify), redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            self.assertEqual(len(calls), 2)
            self.assertEqual(list((root/"artifacts").iterdir()), []); self.assertEqual(list((root/"bundles").iterdir()), [])


if __name__ == "__main__":
    unittest.main()
