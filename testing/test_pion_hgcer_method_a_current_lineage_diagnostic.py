"""Synthetic current-lineage authority, arithmetic and persistence tests; no farm."""

from __future__ import annotations

from contextlib import ExitStack
import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src" / "cuts"), str(ROOT / "testing")]
import pion_hgcer_method_a_current_lineage_diagnostic as diagnostic
import pion_hgcer_method_a_parallel_full_procedure as lineage
import pion_hgcer_method_a_parent_preserving_correction as f4
import test_pion_hgcer_method_a_tphi_propagation as fixtures
import analyze_pion_hgcer_method_a_current_lineage_diagnostic as analyzer


def encoded(value):
    return json.dumps(value, sort_keys=True).encode("utf-8")


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


class DiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = fixtures._artifacts()
        for artifact in cls.inputs:
            for row in artifact["contract"]["application_records"]:
                row["analysis_MM"] = 1.12
                row["analysis_t"] = (row["t_low"] + row["t_high"]) / 2
            fixtures.f1_fixtures._seal_f1_artifact(artifact)
        cls.hashes = {f"{a['setting']['phi_setting']}-{a['setting']['epsilon_filename_token']}": digest(a) for a in cls.inputs}
        cls.map = fixtures.f4_fixtures._f3(cls.inputs, cls.hashes)
        cls.map_sha = digest(cls.map)
        cls.f3_authority = fixtures.f4_fixtures._authority(cls.map, cls.map_sha)
        cls.f3_authority["Q4p4W2p74"]["farm_source_head"] = "0" * 40
        cls.correction = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(
            cls.inputs, cls.map, f1_input_file_hashes=cls.hashes, f3_input_file_sha256=cls.map_sha,
            accepted_f3_runtime_authority_by_kinematic=cls.f3_authority, input_paths={"f1": {}, "f3": "synthetic"})
        assert cls.correction["correction"]["available"]
        cls.correction_sha = digest(cls.correction)
        core = cls.correction["correction"]
        cls.f4_authority = {"Q4p4W2p74": {
            "source_file_sha256": cls.correction_sha, "correction_fingerprint": core["fingerprint"],
            "artifact_fingerprint": cls.correction["artifact_fingerprint"],
            "farm_source_head": lineage.F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
            **{key: core[key] for key in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint")},
            "f1_source_file_sha256": cls.hashes}}

    def authority(self):
        # Synthetic pinned fixture only; production APIs expose no authority override.
        stack = ExitStack()
        stack.enter_context(mock.patch.object(lineage, "F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC", self.f3_authority))
        stack.enter_context(mock.patch.object(lineage, "F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC", self.f4_authority))
        return stack

    def build(self, **overrides):
        kwargs = {"f1_input_file_hashes": self.hashes, "f3_input_file_sha256": self.map_sha,
                  "f4_input_file_sha256": self.correction_sha}
        kwargs.update(overrides)
        with self.authority():
            return diagnostic.build_diagnostic(self.inputs, self.map, self.correction, **kwargs)

    def test_same_inputs_are_deterministic_and_not_mutated(self):
        snapshot = copy.deepcopy((self.inputs, self.map, self.correction))
        a, b = self.build(), self.build()
        self.assertEqual(a, b)
        self.assertEqual((self.inputs, self.map, self.correction), snapshot)
        self.assertTrue(a["f4_exact_reproduction"]["payload_identical"])
        self.assertEqual(a["scope"]["parents"], [0, 1, 2])
        core = dict(a); sha = core.pop("fingerprint")
        self.assertEqual(sha, diagnostic.fingerprint(core))
        a["display_policy"]["phi_edges"][0] = -999
        a["display_policy"]["fixed_edges"]["SHMS_delta"][0] = -999
        self.assertEqual((self.inputs, self.map, self.correction), snapshot)
        self.assertEqual(diagnostic.DISPLAY_EDGES["SHMS_delta"][0], -15)

    def test_scope_and_all_wrong_input_hashes_fail_closed(self):
        for kwargs in ({"kinematic": "Q3p0W2p32"}, {"setting_id": "Left-highe"},
                       {"f3_input_file_sha256": "a" * 64}, {"f4_input_file_sha256": "a" * 64},
                       {"f1_input_file_hashes": {**self.hashes, "Left-lowe": "a" * 64}}):
            with self.subTest(kwargs=kwargs), self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                self.build(**kwargs)

    def test_historical_authority_is_not_accepted(self):
        with self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
            diagnostic.build_diagnostic(self.inputs, self.map, self.correction,
                f1_input_file_hashes=self.hashes, f3_input_file_sha256=self.map_sha,
                f4_input_file_sha256=self.correction_sha)

    def test_tampered_map_and_f4_fingerprints_rejected(self):
        for field in ("fingerprint", "algorithm_fingerprint"):
            changed = copy.deepcopy(self.map)
            changed["acceptance_map"][field] = "a" * 64
            with self.authority(), self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                diagnostic.build_diagnostic(self.inputs, changed, self.correction,
                    f1_input_file_hashes=self.hashes, f3_input_file_sha256=self.map_sha, f4_input_file_sha256=self.correction_sha)
        for field in ("artifact_fingerprint",):
            changed = copy.deepcopy(self.correction); changed[field] = "a" * 64
            with self.authority(), self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                diagnostic.build_diagnostic(self.inputs, self.map, changed,
                    f1_input_file_hashes=self.hashes, f3_input_file_sha256=self.map_sha, f4_input_file_sha256=self.correction_sha)

    def test_exact_f4_reproduction_is_mandatory_before_measurement(self):
        shared = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        def mismatch(*args, **kwargs):
            core, rows = shared(*args, **kwargs)
            core["parents"][0]["parent_normalization"] += .01
            return core, rows
        with mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=mismatch), \
             mock.patch.object(diagnostic, "_proxy_parent") as measure, \
             self.assertRaisesRegex(diagnostic.CurrentLineageDiagnosticError, "exact_reproduction"):
            self.build()
        measure.assert_not_called()

    def test_pid_semantics_and_truth_limitations(self):
        result = self.build()
        self.assertEqual(result["pid_semantics"]["kaon_pid_category"], "P_hgcer_npeSum == 0")
        self.assertEqual(result["pid_semantics"]["weak_positive"], "0 < P_hgcer_npeSum <= 2")
        flags = result["measurement_a"]["limitations"]
        for name in ("direct_hgc_free_true_pion_tag_present_in_consumed_artifacts", "direct_pion_to_kaon_misid_calibration_performed", "absolute_misid_probability_constructed", "proxy_validity_established"):
            self.assertIs(flags[name], False)
        self.assertEqual(flags["proxy_validity_evidence_label"], "NOT VERIFIED")
        self.assertNotIn("NPE=0 pion sample", json.dumps(result))

    def test_weak_positive_counts_bands_and_population_match(self):
        result = self.build()
        training = self.inputs[0]["contract"]["method_a_training_records"]
        application = self.inputs[0]["contract"]["application_records"]
        for p in result["measurement_a"]["parents"]:
            index = p["canonical_t_index"]
            selected = [r for r in training if r["t_index"] == index]
            self.assertEqual(p["classes"]["weak_positive"]["count"], sum(0 < r["P_hgcer_npeSum"] <= 2 for r in selected))
            self.assertEqual(sum(b["count"] for b in p["positive_npe_bands"]), len(selected))
            physical = sum(r["source_label"] == "prompt" and r["t_index"] == index for r in application)
            self.assertEqual(p["training_physical_population_shift"]["matched_count"], physical)
            self.assertEqual(p["training_physical_population_shift"]["physical_control_only_count"], 0)
            self.assertEqual(set(p["training_physical_population_shift"]["physical_minus_training_control_quantile_differences"]),
                             set(diagnostic.FEATURES) | {"raw_relative_response"})
            for c in p["classes"].values():
                self.assertEqual(c["in_support_count"] + c["ood_count"], c["count"])

    def test_npe_band_boundaries_and_zero_rejection(self):
        self.assertEqual([diagnostic._band(n) for n in [.01, 1, 1.01, 2, 2.01, 4, 4.01]], [0, 0, 1, 1, 2, 2, 3])
        with self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
            diagnostic._band(0)

    def test_no_cancellation_and_hand_calculable_cancellation(self):
        self.assertEqual(diagnostic.signed_measure([2, 3], [1, 2])["cancellation_ratio"], 1)
        totals = diagnostic.signed_measure([10, -9], [2, 1])
        for name, value in {"B_signed": 1, "A_abs": 19, "U_signed": 11, "V_abs": 29,
                            "N_signed": 11, "N_abs": 29/19, "cancellation_ratio": 1/19,
                            "normalization_ratio": 209/29, "relative_normalization_difference": 180/29}.items():
            self.assertAlmostEqual(totals[name], value)
        rows = [{"source_label": "prompt"}, {"source_label": "rand"}]
        split = diagnostic.source_decomposition(rows, np.array([10., -9.]), np.array([2., 1.]), totals)
        self.assertEqual(split[0]["signed_fraction"], 10)
        self.assertEqual(split[1]["signed_fraction"], -9)
        self.assertAlmostEqual(sum(s["absolute_support_fraction"] for s in split), 1)
        # Strong redistribution, unchanged signed closure; no N_abs factor returned.
        self.assertAlmostEqual(sum(np.array([10., -9.]) * np.array([2., 1.]) / totals["N_signed"]), 1)
        self.assertNotIn("correction_factors", totals)

    def test_invalid_denominators_and_nonfinite_inputs_rejected(self):
        for b, r in (([0], [1]), ([1, -1], [1, 1]), ([1, -2], [1, 1]),
                     ([2, -1], [1, 3]), ([1], [float("nan")]), ([1], [0])):
            with self.subTest(b=b, r=r), self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                diagnostic.signed_measure(b, r)

    def test_four_authoritative_sources_and_zero_support_are_explicit(self):
        rows = [{"source_label": name} for name in ("prompt", "rand", "dummy", "dummy_rand")]
        b, r = np.asarray([10., -2., -1., 0.]), np.asarray([2., 1., 3., 4.])
        totals = diagnostic.signed_measure(b, r)
        split = diagnostic.source_decomposition(rows, b, r, totals)
        self.assertEqual({row["source_label"] for row in split}, {"prompt", "rand", "dummy", "dummy_rand"})
        self.assertEqual(sum(row["B_signed"] for row in split), 7)
        self.assertEqual(sum(row["A_abs"] for row in split), 13)
        self.assertEqual(sum(row["U_signed"] for row in split), 15)
        zero = next(row for row in split if row["source_label"] == "dummy_rand")
        self.assertIsNone(zero["cancellation_ratio"])
        self.assertFalse(zero["cancellation_available"])

    def test_application_aggregates_preserve_exact_parent_without_child_normalization(self):
        for p in self.build()["measurement_b"]["parents"]:
            totals = p["normalization"]
            phi = p["phi"]["bins"]
            self.assertAlmostEqual(sum(r["baseline_signed_sum"] for r in phi), totals["B_signed"])
            self.assertAlmostEqual(sum(r["adjusted_signed_sum"] for r in phi), p["accepted_adjusted_sum"])
            self.assertAlmostEqual(p["accepted_adjusted_sum"], totals["B_signed"])
            for table in p["coordinate_bins"].values():
                self.assertEqual(sum(r["count"] for r in table["bins"]), totals["count"])
                self.assertAlmostEqual(sum(r["absolute_baseline_support"] for r in table["bins"]), totals["A_abs"])

    def test_display_underflow_overflow_and_edges_preserve_all_support(self):
        table = diagnostic.binned([-5, 0, 1, 2, 8], [0, 1, 2], np.ones(5), np.ones(5))
        self.assertEqual([r["count"] for r in table["bins"]], [1, 1, 2, 1])
        self.assertEqual(sum(r["baseline_signed_sum"] for r in table["bins"]), 5)

    def test_recursive_persistence_guard_and_no_numerical_method_b_dependency(self):
        for key in ("entry_index", "correction_factors", "raw_shape_factors", "production_weight", "identities"):
            with self.subTest(key=key), self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                diagnostic.guard_aggregate_persistence({"nested": [{key: [1, 2]}]})
        result = self.build()
        diagnostic.guard_aggregate_persistence(result)
        self.assertFalse(result["method_b_numerical_dependency"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["alternative_correction_constructed"])
        self.assertNotRegex(Path(diagnostic.__file__).read_text(), r"(?:import|from)\s+(?:ROOT|pion_hgcer_refinement_method_b)")

    def test_deterministic_filenames_overwrite_refusal_and_wrong_scope(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            js, pdf = analyzer.output_paths(root, "Q4p4W2p74")
            self.assertEqual(js.name, diagnostic.STEM + ".json")
            self.assertEqual(pdf.name, diagnostic.STEM + ".pdf")
            for args in ((root, "Q3p0W2p32"), (root, "Q4p4W2p74", "Right", "lowe")):
                with self.assertRaises(diagnostic.CurrentLineageDiagnosticError):
                    analyzer.output_paths(*args)
            js.write_text("preserve")
            with mock.patch.object(analyzer, "load_diagnostic") as loader:
                self.assertEqual(analyzer.main(["--outdir", directory, "--kinematic", "Q4p4W2p74"]), 1)
            loader.assert_not_called()
            self.assertEqual(js.read_text(), "preserve")

    def test_analyzer_real_file_load_pure_python_pdf_and_json(self):
        with tempfile.TemporaryDirectory() as directory, self.authority():
            root = Path(directory)
            paths = lineage.accepted_f6_3_artifact_paths(root, "Q4p4W2p74")
            for artifact in self.inputs:
                setting = artifact["setting"]
                Path(paths["f1"][f"{setting['phi_setting']}-{setting['epsilon_filename_token']}"]).write_bytes(encoded(artifact))
            Path(paths["f3"]).write_bytes(encoded(self.map))
            Path(paths["f4"]).write_bytes(encoded(self.correction))
            before = set(root.iterdir())
            with mock.patch.object(analyzer, "_git_value", return_value="synthetic provenance"):
                self.assertEqual(analyzer.main(["--outdir", directory, "--kinematic", "Q4p4W2p74"]), 0)
            js, pdf = analyzer.output_paths(root, "Q4p4W2p74")
            self.assertEqual(set(root.iterdir()) - before, {js, pdf})
            result = json.loads(js.read_text())
            diagnostic.guard_aggregate_persistence(result)
            self.assertEqual(len(re.findall(rb"/Type /Page\b", pdf.read_bytes())), 8)
            other = root / "repeat.pdf"
            analyzer.write_review_pdf(other, result)
            self.assertEqual(pdf.read_bytes(), other.read_bytes())
            # Wrong candidate bytes must fail before creating either output.
            js.unlink(); pdf.unlink(); other.unlink()
            Path(paths["f3"]).write_text("{}")
            self.assertEqual(analyzer.main(["--outdir", directory, "--kinematic", "Q4p4W2p74"]), 1)
            self.assertFalse(js.exists()); self.assertFalse(pdf.exists())


if __name__ == "__main__":
    unittest.main()
