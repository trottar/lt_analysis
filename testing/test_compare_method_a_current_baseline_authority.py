"""Focused orchestration and comparison tests; no ROOT or farm artifacts."""

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from testing import compare_method_a_current_baseline_authority as audit


def payload(name, scientific):
    return {name: {**scientific, "fingerprint": "old", "fingerprint_inputs": {}, "input_fingerprints": []}}


def parents(normalization=1.0):
    result = []
    for alias in audit.ALIASES:
        for index in range(3):
            row = {metric: 1.0 for metric in audit.PARENT_METRICS}
            row.update({"setting_id": alias, "canonical_t_index": index, "parent_normalization": normalization,
                        "raw_shape_factor_summary": {"p50": 1.0}, "correction_factor_summary": {"p50": 1.0},
                        "source_diagnostics": [{"source_label": "prompt", "event_count": 2, "baseline_signed_sum": 2.0,
                                                "adjusted_signed_sum": 2.0, "signed_delta": 0.0}],
                        "canonical_phi_diagnostics": [{"phi_index": 0, "phi_low": 0.0, "phi_high": 1.0,
                                                       "event_count": 2, "baseline_signed_sum": 2.0,
                                                       "adjusted_signed_sum": 2.0, "signed_delta": 0.0}]})
            result.append(row)
    return result


class ComparisonTests(unittest.TestCase):
    def test_projection_is_owned_by_runtime_and_preserves_types_and_unknown_fields(self):
        import pion_hgcer_method_a_parallel_full_procedure as runtime
        self.assertIs(audit.scientific_projection, runtime.scientific_projection)
        self.assertIs(audit.first_mismatch, runtime.first_mismatch)
        for name in ("F2_PROVENANCE", "F3_PROVENANCE", "F4_PROVENANCE"):
            self.assertIs(getattr(audit, name), getattr(runtime, name))
        self.assertEqual(runtime.first_mismatch({"new_field": 1}, {"new_field": True})["path"], "$.new_field")

    def test_inventory(self):
        entries = [f"{name}=x{i}.json" for i, name in enumerate(audit.ALIASES)]
        self.assertEqual(len(audit._f1_inputs(entries)), 5)
        for bad in (entries[:-1], entries + [entries[0]], entries[:-1] + ["Wrong-low=x"]):
            with self.assertRaises(ValueError):
                audit._f1_inputs(bad)

    def test_f2_projection_and_numerics(self):
        body = {"candidate_definitions": [1], "algorithm_config": {}, "algorithm_fingerprint": "a",
                "response_support": {}, "groups": [0] * 15, "candidate_summaries": [], "recommendation": "hgcer3",
                "fingerprint": "a", "fingerprint_inputs": {"input_content": "old"}, "input_fingerprints": [{"hash": "old"}]}
        changed = copy.deepcopy(body)
        changed.update(fingerprint="b", fingerprint_inputs={"input_content": "new"}, input_fingerprints=[{"hash": "new"}])
        self.assertIsNone(audit.first_mismatch(audit.scientific_projection(body, audit.F2_PROVENANCE), audit.scientific_projection(changed, audit.F2_PROVENANCE)))
        changed["groups"][0] = 2
        self.assertEqual(audit.first_mismatch(audit.scientific_projection(body, audit.F2_PROVENANCE), audit.scientific_projection(changed, audit.F2_PROVENANCE))["path"], "$.groups[0]")
        changed["new_scientific_field"] = 9
        self.assertIn("new_scientific_field", audit.scientific_projection(changed, audit.F2_PROVENANCE))

    def test_f3_projection_coefficients_scaler_support(self):
        body = {"accepted_basis": "hgcer3", "ordered_features": ["x"], "algorithm_config": {}, "algorithm_fingerprint": "a",
                "models": [{"coefficients": [1], "scaler": {"median": 1}, "support": {"fraction": 1}}] * 15,
                "f2_representation_fingerprint": "old", "f2_source_file_sha256": "old", "input_fingerprints": [], "fingerprint_inputs": {}, "fingerprint": "old"}
        changed = copy.deepcopy(body)
        changed["f2_representation_fingerprint"] = "new"
        self.assertIsNone(audit.first_mismatch(audit.scientific_projection(body, audit.F3_PROVENANCE), audit.scientific_projection(changed, audit.F3_PROVENANCE)))
        for field, value in (("coefficients", [2]), ("scaler", {"median": 2}), ("support", {"fraction": 0.5})):
            other = copy.deepcopy(changed)
            other["models"][0][field] = value
            self.assertIn(field, audit.first_mismatch(audit.scientific_projection(body, audit.F3_PROVENANCE), audit.scientific_projection(other, audit.F3_PROVENANCE))["path"])

    def test_f4_and_stage(self):
        a = {"accepted_basis": "hgcer3", "algorithm_config": {}, "parents": parents()}
        c = copy.deepcopy(a)
        c["parents"][0]["parent_normalization"] = 1.25
        details = audit.compare_f4(a, c)
        self.assertEqual(details["global_maxima"]["parent_normalization"]["absolute_difference"], 0.25)
        self.assertEqual(details["setting_maxima"][audit.ALIASES[0]]["parent_normalization"]["relative_difference"], 0.25)
        for values, expected in (((False, False, False), "F2"), ((True, False, False), "F3"), ((True, True, False), "F4"), ((True, True, True), "none")):
            self.assertEqual(audit.first_changed_stage(*values), expected)
        with self.assertRaises(ValueError):
            audit.compare_f4(a, {**c, "parents": c["parents"][:-1]})

    def test_nested_f4_diagnostic_deltas(self):
        accepted = {"accepted_basis": "hgcer3", "algorithm_config": {}, "parents": parents()}
        candidate = copy.deepcopy(accepted)
        candidate["parents"][0]["source_diagnostics"][0]["baseline_signed_sum"] = 2.5
        candidate["parents"][0]["canonical_phi_diagnostics"][0]["adjusted_signed_sum"] = 1.5
        metrics = next(row["metrics"] for row in audit.compare_f4(accepted, candidate)["parents"]
                       if row["setting_id"] == audit.ALIASES[0] and row["canonical_t_index"] == 0)
        source = metrics["source_diagnostics"][0]
        phi = metrics["canonical_phi_diagnostics"][0]
        self.assertEqual(source["baseline_signed_sum"], {"accepted": 2.0, "candidate": 2.5,
                                                         "absolute_difference": 0.5, "relative_difference": 0.25})
        self.assertEqual(phi["adjusted_signed_sum"], {"accepted": 2.0, "candidate": 1.5,
                                                       "absolute_difference": 0.5, "relative_difference": 0.25})
        self.assertEqual(source["source_label"], {"accepted": "prompt", "candidate": "prompt", "match": True})
        self.assertEqual(phi["phi_index"], {"accepted": 0, "candidate": 0,
                                             "absolute_difference": 0, "relative_difference": None})
        leaves = list(audit._numeric_deltas(metrics["source_diagnostics"]))
        self.assertIn(source["baseline_signed_sum"], leaves)
        self.assertIn(phi["adjusted_signed_sum"], list(audit._numeric_deltas(metrics["canonical_phi_diagnostics"])))

    def test_unequal_nested_list_lengths(self):
        accepted = [{"source_label": "prompt", "baseline_signed_sum": 2.0}]
        candidate = accepted + [{"source_label": "random", "baseline_signed_sum": -1.0}]
        added = audit._delta(accepted, candidate)
        self.assertEqual(len(added), 2)
        self.assertEqual(added[1], {"accepted": None, "candidate": candidate[1], "missing_side": "accepted"})
        removed = audit._delta(candidate, accepted)
        self.assertEqual(len(removed), 2)
        self.assertEqual(removed[1], {"accepted": candidate[1], "candidate": None, "missing_side": "candidate"})

    def test_strict_json_and_writer_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.json"
            path.write_bytes(b'{"x":1, "x":2}')
            with self.assertRaises(ValueError):
                audit._read(path)
            path.write_bytes(b'{"x":NaN}')
            with self.assertRaises(ValueError):
                audit._read(path)
            artifact = {"z": [1], "a": "text"}
            out = Path(tmp) / "writer.json"
            audit.f2.write_pion_hgcer_method_a_acceptance_representation_json(out, artifact)
            self.assertEqual(out.read_bytes(), audit._writer_bytes(artifact))
            audit.f3.write_pion_hgcer_method_a_acceptance_map_json(out, artifact)
            self.assertEqual(out.read_bytes(), audit._writer_bytes(artifact))

    def test_orchestration_sha_override_write_safety(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f1_files = []
            original = {}
            for i, alias in enumerate(audit.ALIASES):
                phi, epsilon = alias.split("-")
                path = root / f"f1_{i}.json"
                path.write_bytes((json.dumps({"setting": {"phi_setting": phi, "epsilon_filename_token": epsilon, "kinematic_token": "Q4p4W2p74"}, "contract": {"method_a_training_records": [1], "application_records": [2]}}) + "\n").encode())
                f1_files.append(f"{alias}={path}")
                original[path] = path.read_bytes()
            f2_body = {"candidate_definitions": [], "algorithm_config": {}, "algorithm_fingerprint": "a", "response_support": {}, "groups": [0] * 15, "candidate_summaries": [], "recommendation": "hgcer3"}
            f3_body = {"accepted_basis": "hgcer3", "ordered_features": ["x"], "algorithm_config": {}, "algorithm_fingerprint": "b", "models": [0] * 15, "fingerprint": "map"}
            f4_body = {"accepted_basis": "hgcer3", "algorithm_config": {}, "parents": parents()}
            expected_inputs = [{"setting_id": alias, "setting": {"phi_setting": alias.split("-")[0], "epsilon_filename_token": alias.split("-")[1]}, "source_file_sha256": hashlib.sha256(original[root / f"f1_{i}.json"]).hexdigest(), "stable_f1_content_fingerprint": "stable", "fingerprint": "contract", "method_a_training_population_fingerprint": "train", "application_population_fingerprint": "application"} for i, alias in enumerate(audit.ALIASES)]
            candidate_f2 = payload("representation", {**f2_body, "input_fingerprints": expected_inputs})
            candidate_f3 = {**payload("acceptance_map", f3_body), "artifact_fingerprint": "artifact"}
            candidate_f4 = payload("correction", f4_body)
            accepted_paths = {}
            for stage, artifact in (("f2", candidate_f2), ("f3", candidate_f3), ("f4", candidate_f4)):
                path = root / f"accepted_{stage}.json"
                path.write_bytes(audit._writer_bytes(artifact))
                original[path] = path.read_bytes()
                accepted_paths[stage] = path
            output = root / "comparison.json"
            args = SimpleNamespace(f1=f1_files, accepted_f2=accepted_paths["f2"], accepted_f3=accepted_paths["f3"], accepted_f4=accepted_paths["f4"], output=output, overwrite=False)
            with mock.patch.object(audit.f2, "build_pion_hgcer_method_a_acceptance_representation_artifact", return_value=candidate_f2) as build2, mock.patch.object(audit.f3, "build_pion_hgcer_method_a_acceptance_map_artifact", return_value=candidate_f3) as build3, mock.patch.object(audit.f4, "build_pion_hgcer_method_a_parent_preserving_correction_artifact", return_value=candidate_f4) as build4:
                result = audit.run(args)
            self.assertEqual(build2.call_count, 1)
            self.assertEqual(build3.call_count, 1)
            self.assertEqual(build4.call_count, 1)
            f2_sha = hashlib.sha256(audit._writer_bytes(candidate_f2)).hexdigest()
            f3_sha = hashlib.sha256(audit._writer_bytes(candidate_f3)).hexdigest()
            self.assertEqual(build3.call_args.kwargs["f2_input_file_sha256"], f2_sha)
            self.assertEqual(build4.call_args.kwargs["f3_input_file_sha256"], f3_sha)
            override = build4.call_args.kwargs["accepted_f3_runtime_authority_by_kinematic"]["Q4p4W2p74"]
            self.assertEqual(override["source_file_sha256"], f3_sha)
            self.assertEqual(override["map_fingerprint"], candidate_f3["acceptance_map"]["fingerprint"])
            self.assertEqual(override["algorithm_fingerprint"], "b")
            self.assertEqual(override["artifact_fingerprint"], "artifact")
            self.assertEqual(override["farm_source_head"], "0" * 40)
            self.assertEqual(result["summary"]["first_changed_stage"], "none")
            self.assertEqual(result["f1_inputs"][0]["source_file_sha256"], hashlib.sha256(original[root / "f1_0.json"]).hexdigest())
            with self.assertRaises(FileExistsError):
                audit.run(args)
            args.overwrite = True
            with mock.patch.object(audit.f2, "build_pion_hgcer_method_a_acceptance_representation_artifact", return_value=candidate_f2), mock.patch.object(audit.f3, "build_pion_hgcer_method_a_acceptance_map_artifact", return_value=candidate_f3), mock.patch.object(audit.f4, "build_pion_hgcer_method_a_parent_preserving_correction_artifact", return_value=candidate_f4):
                audit.run(args)
            self.assertTrue(all(path.read_bytes() == content for path, content in original.items()))
            accepted_paths["f2"].write_bytes(b'{"bad":')
            with self.assertRaises(json.JSONDecodeError):
                audit.run(args)
            cli = [item for spec in f1_files for item in ("--f1", spec)] + ["--accepted-f2", str(accepted_paths["f2"]), "--accepted-f3", str(accepted_paths["f3"]), "--accepted-f4", str(accepted_paths["f4"]), "--output", str(output), "--overwrite"]
            with mock.patch("sys.stderr"), self.assertRaises(SystemExit) as failed:
                audit.main(cli)
            self.assertEqual(failed.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
