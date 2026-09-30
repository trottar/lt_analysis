"""Pure-Python materializer gates and publication tests; synthetic inputs only."""

import ast
import copy
import hashlib
import json
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from testing import materialize_method_a_current_baseline_authority as materialize
from testing.test_compare_method_a_current_baseline_authority import parents, payload


class MaterializerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "output"
        self.output.mkdir()
        self.original_inputs = {}
        self.f1_paths = []
        self.f1_rows = []
        self.f2_fp_rows = []
        for index, alias in enumerate(materialize.comparison.ALIASES):
            phi, epsilon = alias.split("-")
            setting = {"phi_setting": phi, "epsilon_filename_token": epsilon, "kinematic_token": "Q4p4W2p74"}
            path = self.root / f"f1_{index}.json"
            path.write_bytes(materialize.comparison._writer_bytes({"setting": setting, "contract": {"method_a_training_records": [1], "application_records": [2]}}))
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            self.original_inputs[path] = path.read_bytes()
            self.f1_paths.append(f"{alias}={path}")
            self.f1_rows.append({"alias": alias, "setting_id": alias, "setting": setting, "source_file_sha256": sha,
                                 "stable_f1_content_fingerprint": "stable", "f1_contract_fingerprint": "contract",
                                 "training_record_count": 1, "application_record_count": 1,
                                 "training_population_fingerprint": "train", "application_population_fingerprint": "application"})
            self.f2_fp_rows.append({"setting_id": alias, "setting": setting, "source_file_sha256": sha,
                                    "stable_f1_content_fingerprint": "stable", "fingerprint": "contract",
                                    "method_a_training_population_fingerprint": "train", "application_population_fingerprint": "application"})
        f2_body = {"candidate_definitions": [], "algorithm_config": {}, "algorithm_fingerprint": "f2_algorithm",
                   "response_support": {}, "groups": [0] * 15, "candidate_summaries": [], "recommendation": "hgcer3"}
        self.candidate_f2 = {**payload("representation", f2_body), "artifact_fingerprint": "f2_artifact", "non_authoritative": True}
        self.candidate_f2["representation"]["input_fingerprints"] = self.f2_fp_rows
        self.candidate_f2["representation"]["fingerprint"] = "candidate_f2_provenance"
        f3_body = {"accepted_basis": "hgcer3", "ordered_features": ["x"], "algorithm_config": {},
                   "algorithm_fingerprint": "f3_algorithm", "models": [0] * 15}
        self.candidate_f3 = {**payload("acceptance_map", f3_body), "artifact_fingerprint": "f3_artifact", "non_authoritative": True}
        self.candidate_f3["acceptance_map"]["fingerprint"] = "candidate_f3_provenance"
        self.f2_sha = hashlib.sha256(materialize.comparison._writer_bytes(self.candidate_f2)).hexdigest()
        self.f3_sha = hashlib.sha256(materialize.comparison._writer_bytes(self.candidate_f3)).hexdigest()
        self.override = {"Q4p4W2p74": {"source_file_sha256": self.f3_sha, "map_fingerprint": "candidate_f3_provenance",
                                         "algorithm_fingerprint": "f3_algorithm", "artifact_fingerprint": "f3_artifact",
                                         "farm_source_head": "0" * 40}}
        f4_body = {"accepted_basis": "hgcer3", "algorithm_config": {}, "parents": parents(1.25),
                   "f3_source_file_sha256": self.f3_sha, "f3_map_fingerprint": "old",
                   "f3_algorithm_fingerprint": "f3_algorithm", "f3_artifact_fingerprint": "f3_artifact"}
        self.candidate_f4 = {**payload("correction", f4_body), "artifact_fingerprint": "f4_artifact", "non_authoritative": True}
        self.accepted = {"f2": copy.deepcopy(self.candidate_f2), "f3": copy.deepcopy(self.candidate_f3),
                         "f4": copy.deepcopy(self.candidate_f4)}
        self.accepted["f2"]["representation"]["fingerprint"] = "accepted_f2_provenance"
        self.accepted["f3"]["acceptance_map"]["fingerprint"] = "accepted_f3_provenance"
        self.accepted["f4"]["correction"]["parents"] = parents(1.0)
        self.accepted_paths = {}
        self.accepted_sha = {}
        for stage, artifact in self.accepted.items():
            path = self.root / f"accepted_{stage}.json"
            path.write_bytes(materialize.comparison._writer_bytes(artifact))
            self.original_inputs[path] = path.read_bytes()
            self.accepted_paths[stage] = path
            self.accepted_sha[stage] = hashlib.sha256(path.read_bytes()).hexdigest()
        f2_result = materialize._scientific_result(self.accepted["f2"], self.candidate_f2, "representation", materialize.comparison.F2_PROVENANCE)
        f3_result = materialize._scientific_result(self.accepted["f3"], self.candidate_f3, "acceptance_map", materialize.comparison.F3_PROVENANCE)
        f4_result = materialize._scientific_result(self.accepted["f4"], self.candidate_f4, "correction", materialize.comparison.F4_PROVENANCE)
        f4_result.update(materialize.comparison.compare_f4(self.accepted["f4"]["correction"], self.candidate_f4["correction"]))
        self.comparison = {"schema_version": "method_a_current_baseline_authority_comparison/v1", "non_authoritative": True,
                           "f1_inputs": self.f1_rows, "accepted_file_sha256": copy.deepcopy(self.accepted_sha),
                           "candidate_serialized_sha256": {"f2": self.f2_sha, "f3": self.f3_sha},
                           "diagnostic_f3_authority_override": self.override,
                           "f2": f2_result, "f3": f3_result, "f4": f4_result,
                           "summary": {"f2_scientific_payload_match": True, "f3_scientific_payload_match": True,
                                       "f4_scientific_payload_match": False, "first_changed_stage": "F4"}}
        self.comparison_path = self.root / "reviewed_comparison.json"
        self._write_comparison()
        self.args = SimpleNamespace(f1=self.f1_paths, accepted_f2=self.accepted_paths["f2"],
                                    accepted_f3=self.accepted_paths["f3"], accepted_f4=self.accepted_paths["f4"],
                                    comparison=self.comparison_path, expected_comparison_sha256=self.comparison_sha,
                                    source_head="A" * 40, output_dir=self.output, overwrite=False)

    def _write_comparison(self):
        self.comparison_path.write_bytes(materialize.comparison._writer_bytes(self.comparison))
        self.comparison_sha = hashlib.sha256(self.comparison_path.read_bytes()).hexdigest()
        self.original_inputs[self.comparison_path] = self.comparison_path.read_bytes()
        if hasattr(self, "args"):
            self.args.expected_comparison_sha256 = self.comparison_sha

    @contextmanager
    def _builders(self, f2=None, f3=None, f4=None):
        order = []
        def builder(stage, value):
            def call(*args, **kwargs):
                order.append((stage, args, kwargs))
                return copy.deepcopy(value)
            return call
        with mock.patch.object(materialize.comparison.f2, "build_pion_hgcer_method_a_acceptance_representation_artifact", side_effect=builder("f2", self.candidate_f2 if f2 is None else f2)) as b2, \
             mock.patch.object(materialize.comparison.f3, "build_pion_hgcer_method_a_acceptance_map_artifact", side_effect=builder("f3", self.candidate_f3 if f3 is None else f3)) as b3, \
             mock.patch.object(materialize.comparison.f4, "build_pion_hgcer_method_a_parent_preserving_correction_artifact", side_effect=builder("f4", self.candidate_f4 if f4 is None else f4)) as b4:
            yield order, (b2, b3, b4)

    def _no_outputs(self):
        self.assertEqual(list(self.output.iterdir()), [])

    def test_inventory_and_input_preflight(self):
        with self._builders():
            for invalid in (self.f1_paths[:-1], self.f1_paths + [self.f1_paths[0]], self.f1_paths[:-1] + ["Wrong-low=x"], self.f1_paths[:-1] + ["Left-lowe"]):
                with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                    materialize.run(SimpleNamespace(**{**vars(self.args), "f1": invalid}))
            for invalid in ("nothex", "G" * 64):
                with self.subTest(sha=invalid), self.assertRaises(ValueError):
                    materialize.run(SimpleNamespace(**{**vars(self.args), "expected_comparison_sha256": invalid}))
            with self.assertRaises(ValueError):
                materialize.run(SimpleNamespace(**{**vars(self.args), "source_head": "bad"}))
            with self.assertRaises(ValueError):
                materialize.run(SimpleNamespace(**{**vars(self.args), "output_dir": self.root / "missing"}))
            with self.assertRaises(ValueError):
                materialize.run(SimpleNamespace(**{**vars(self.args), "accepted_f2": self.args.accepted_f3}))
        self._no_outputs()

    def test_raw_comparison_accepted_and_f1_hash_gates(self):
        with self._builders():
            with self.assertRaisesRegex(ValueError, "comparison raw SHA mismatch"):
                materialize.run(SimpleNamespace(**{**vars(self.args), "expected_comparison_sha256": "0" * 64}))
            for stage in ("f2", "f3", "f4"):
                path = self.accepted_paths[stage]
                path.write_bytes(path.read_bytes() + b" ")
                with self.subTest(stage=stage), self.assertRaisesRegex(ValueError, f"accepted {stage} raw SHA mismatch"):
                    materialize.run(self.args)
                path.write_bytes(self.original_inputs[path])
            f1_path = self.root / "f1_0.json"
            f1_path.write_bytes(f1_path.read_bytes() + b" ")
            with self.assertRaisesRegex(ValueError, "current F.1 SHA mismatch"):
                materialize.run(self.args)
        self._no_outputs()

    def test_comparison_schema_and_f1_provenance_fail_closed(self):
        self.comparison["non_authoritative"] = False
        self._write_comparison()
        with self._builders(), self.assertRaisesRegex(ValueError, "comparison schema/non-authoritative status invalid"):
            materialize.run(self.args)
        self.comparison["non_authoritative"] = True
        self.comparison["f1_inputs"][0]["stable_f1_content_fingerprint"] = "wrong"
        self._write_comparison()
        with self._builders(), self.assertRaisesRegex(ValueError, "current F.1 fingerprint/provenance mismatch"):
            materialize.run(self.args)
        self._no_outputs()

    def test_rebuild_hashes_and_scientific_gates(self):
        with self._builders() as (order, mocks):
            manifest = materialize.run(self.args)
        self.assertEqual([item[0] for item in order], ["f2", "f3", "f4"])
        self.assertTrue(all(builder.call_count == 1 for builder in mocks))
        self.assertEqual(order[1][2]["f2_input_file_sha256"], self.f2_sha)
        self.assertEqual(order[2][2]["f3_input_file_sha256"], self.f3_sha)
        self.assertEqual(order[2][2]["accepted_f3_runtime_authority_by_kinematic"], self.override)
        self.assertEqual(manifest["scientific_gate"]["first_changed_stage"], "F4")
        self.assertEqual(manifest["candidate_outputs"]["f2"]["raw_sha256"], self.f2_sha)
        self.assertEqual(manifest["candidate_outputs"]["f3"]["raw_sha256"], self.f3_sha)
        self.assertTrue(manifest["scientific_gate"]["f2_scientific_payload_match"])
        self.assertTrue(manifest["scientific_gate"]["f3_scientific_payload_match"])
        for stage, body_name, key in (("f2", "representation", "representation_fingerprint"),
                                      ("f3", "acceptance_map", "map_fingerprint"),
                                      ("f4", "correction", "correction_fingerprint")):
            self.assertEqual(manifest["accepted_inputs"][stage][key], self.accepted[stage][body_name]["fingerprint"])
            self.assertEqual(manifest["candidate_outputs"][stage][key], getattr(self, f"candidate_{stage}")[body_name]["fingerprint"])
        for stage in ("f2", "f3"):
            key = {"f2": "representation_fingerprint", "f3": "map_fingerprint"}[stage]
            self.assertNotEqual(manifest["accepted_inputs"][stage][key], manifest["candidate_outputs"][stage][key])
        for stage in ("f2", "f3"):
            for group, artifact in (("accepted_inputs", self.accepted[stage]), ("candidate_outputs", getattr(self, f"candidate_{stage}"))):
                body = artifact[{"f2": "representation", "f3": "acceptance_map"}[stage]]
                self.assertEqual(manifest[group][stage]["algorithm_fingerprint"], body["algorithm_fingerprint"])
        for group in ("accepted_inputs", "candidate_outputs"):
            self.assertEqual(manifest[group]["f3"]["accepted_basis"], "hgcer3")
        self.assertNotIn('"scientific_payload_fingerprint"', json.dumps(manifest))

    def test_candidate_f2_f3_sha_drift_fails_before_publication(self):
        for stage in ("f2", "f3"):
            with self.subTest(stage=stage):
                self.comparison["candidate_serialized_sha256"][stage] = "0" * 64
                self._write_comparison()
                with self._builders(), self.assertRaisesRegex(ValueError, f"candidate F\\.{stage[-1]} writer SHA mismatch"):
                    materialize.run(self.args)
                self._no_outputs()
                self.comparison["candidate_serialized_sha256"][stage] = {"f2": self.f2_sha, "f3": self.f3_sha}[stage]
                self._write_comparison()

    def test_independent_f2_f3_scientific_equality(self):
        for stage, field in (("f2", "recommendation"), ("f3", "accepted_basis")):
            with self.subTest(stage=stage):
                path = self.accepted_paths[stage]
                bad = copy.deepcopy(self.accepted[stage])
                bad[{"f2": "representation", "f3": "acceptance_map"}[stage]][field] = "wrong"
                path.write_bytes(materialize.comparison._writer_bytes(bad))
                self.comparison["accepted_file_sha256"][stage] = hashlib.sha256(path.read_bytes()).hexdigest()
                self._write_comparison()
                with self._builders(), self.assertRaisesRegex(ValueError, f"independent F\\.{stage[-1]} scientific equality failed"):
                    materialize.run(self.args)
                self._no_outputs()
                path.write_bytes(self.original_inputs[path])
                self.comparison["accepted_file_sha256"][stage] = self.accepted_sha[stage]
                self._write_comparison()

    def test_f4_reproduction_and_summary_fail_closed(self):
        self.comparison["f4"]["parents"][0]["metrics"]["baseline_parent_sum"]["candidate"] = 999
        self._write_comparison()
        with self._builders(), self.assertRaisesRegex(ValueError, "F.4 comparison reproduction mismatch: parents"):
            materialize.run(self.args)
        self._no_outputs()
        self.comparison["f4"] = materialize._scientific_result(self.accepted["f4"], self.candidate_f4, "correction", materialize.comparison.F4_PROVENANCE)
        self.comparison["f4"].update(materialize.comparison.compare_f4(self.accepted["f4"]["correction"], self.candidate_f4["correction"]))
        self.comparison["summary"]["first_changed_stage"] = "F3"
        self._write_comparison()
        with self._builders(), self.assertRaisesRegex(ValueError, "comparison summary reproduction mismatch"):
            materialize.run(self.args)
        self._no_outputs()

    def test_public_writers_copy_manifest_last_and_inputs_unchanged(self):
        write_order, publish_order = [], []
        original_writers = [(module, name, getattr(module, name)) for module, name in ((materialize.comparison.f2, materialize.WRITER_NAMES["f2"]), (materialize.comparison.f3, materialize.WRITER_NAMES["f3"]), (materialize.comparison.f4, materialize.WRITER_NAMES["f4"]))]
        original_replace = materialize.os.replace
        def replace(source, target):
            publish_order.append(Path(target).name)
            return original_replace(source, target)
        with self._builders(), mock.patch.object(materialize.os, "replace", side_effect=replace):
            with mock.patch.object(original_writers[0][0], original_writers[0][1], side_effect=lambda p, a: (write_order.append("f2"), original_writers[0][2](p, a))[1]), \
                 mock.patch.object(original_writers[1][0], original_writers[1][1], side_effect=lambda p, a: (write_order.append("f3"), original_writers[1][2](p, a))[1]), \
                 mock.patch.object(original_writers[2][0], original_writers[2][1], side_effect=lambda p, a: (write_order.append("f4"), original_writers[2][2](p, a))[1]):
                manifest = materialize.run(self.args)
        self.assertEqual(write_order, ["f2", "f3", "f4"])
        self.assertEqual(publish_order[-1], materialize._names("Q4p4W2p74")["manifest"])
        self.assertEqual(len(publish_order), 5)
        names = materialize._names("Q4p4W2p74")
        self.assertEqual({path.name for path in self.output.iterdir()}, set(names.values()))
        self.assertFalse({path.name for path in self.accepted_paths.values()} & set(names.values()))
        self.assertEqual((self.output / names["comparison"]).read_bytes(), self.original_inputs[self.comparison_path])
        self.assertEqual(manifest["schema_version"], materialize.SCHEMA)
        self.assertEqual(manifest["source_head"], "a" * 40)
        self.assertTrue(manifest["complete"] and manifest["non_authoritative"])
        self.assertEqual(manifest["errors"], [])
        for flag in ("accepted_authority_mutated", "production_objects_mutated", "production_application_performed", "method_a_promoted"):
            self.assertIs(manifest[flag], False)
        self.assertEqual(manifest["diagnostic_f3_authority_override"], {"accepted_authority": False, "purpose": "diagnostic_candidate_construction_only", "record": self.override})
        self.assertEqual(manifest["comparison_input"]["raw_sha256"], self.comparison_sha)
        self.assertEqual(manifest["accepted_inputs"]["f4"]["raw_sha256"], self.accepted_sha["f4"])
        self.assertEqual(manifest["candidate_outputs"]["f4"]["f3_source_file_sha256"], self.f3_sha)
        self.assertEqual(manifest["scientific_gate"]["f4_global_maxima"], self.comparison["f4"]["global_maxima"])
        self.assertEqual(json.loads((self.output / names["manifest"]).read_bytes()), manifest)
        self.assertTrue(all(path.read_bytes() == content for path, content in self.original_inputs.items()))
        with self._builders(), self.assertRaisesRegex(ValueError, "output already exists"):
            materialize.run(self.args)
        self.args.overwrite = True
        overwrite_order = []
        original_replace = materialize.os.replace
        def overwrite_replace(source, target):
            overwrite_order.append((Path(source), Path(target)))
            return original_replace(source, target)
        with self._builders(), mock.patch.object(materialize.os, "replace", side_effect=overwrite_replace):
            materialize.run(self.args)
        self.assertEqual(overwrite_order[0][0], self.output / names["manifest"])
        self.assertEqual(overwrite_order[-1][1], self.output / names["manifest"])
        self.assertFalse(any(source == self.output / names["manifest"] for source, _target in overwrite_order[1:]))
        self.assertTrue(all(path.read_bytes() == content for path, content in self.original_inputs.items()))

    def test_failed_overwrite_after_candidate_replacement_leaves_no_manifest(self):
        with self._builders():
            materialize.run(self.args)
        names = materialize._names("Q4p4W2p74")
        marker = self.output / names["manifest"]
        self.args.overwrite = True
        original_replace = materialize.os.replace
        replaced = []
        def fail_after_first_candidate(source, target):
            source, target = Path(source), Path(target)
            if target == self.output / names["f2"]:
                self.assertFalse(marker.exists())
                self.assertEqual(replaced, [names["comparison"]])
                raise OSError("candidate publication failure")
            result = original_replace(source, target)
            if target.parent == self.output and target.name in names.values():
                replaced.append(target.name)
            return result
        with self._builders(), mock.patch.object(materialize.os, "replace", side_effect=fail_after_first_candidate):
            with self.assertRaisesRegex(OSError, "candidate publication failure"):
                materialize.run(self.args)
        self.assertFalse(marker.exists())
        self.assertEqual(replaced, [names["comparison"]])
        self.assertEqual((self.output / names["comparison"]).read_bytes(), self.original_inputs[self.comparison_path])
        self.assertTrue(all(path.read_bytes() == content for path, content in self.original_inputs.items()))

    def test_failed_overwrite_before_candidate_replacement_restores_manifest(self):
        with self._builders():
            materialize.run(self.args)
        names = materialize._names("Q4p4W2p74")
        marker = self.output / names["manifest"]
        prior = marker.read_bytes()
        self.args.overwrite = True
        original_replace = materialize.os.replace
        def fail_first_candidate(source, target):
            if Path(target) == self.output / names["comparison"]:
                self.assertFalse(marker.exists())
                raise OSError("publication did not begin")
            return original_replace(source, target)
        with self._builders(), mock.patch.object(materialize.os, "replace", side_effect=fail_first_candidate):
            with self.assertRaisesRegex(OSError, "publication did not begin"):
                materialize.run(self.args)
        self.assertEqual(marker.read_bytes(), prior)
        self.assertTrue(all(path.read_bytes() == content for path, content in self.original_inputs.items()))

    def test_writer_failure_leaves_no_final_outputs(self):
        with self._builders(), mock.patch.object(materialize.comparison.f4, materialize.WRITER_NAMES["f4"], side_effect=OSError("writer failure")):
            with self.assertRaisesRegex(OSError, "writer failure"):
                materialize.run(self.args)
        self._no_outputs()

    def test_public_writer_byte_drift_leaves_no_final_outputs(self):
        def corrupt(path, artifact):
            Path(path).write_bytes(b'{"corrupt": true}\n')
        with self._builders(), mock.patch.object(materialize.comparison.f3, materialize.WRITER_NAMES["f3"], side_effect=corrupt):
            with self.assertRaisesRegex(ValueError, "candidate f3 public-writer bytes differ"):
                materialize.run(self.args)
        self._no_outputs()

    def test_no_root_import(self):
        source = ast.parse(Path(materialize.__file__).read_text(encoding="utf-8"))
        imported = [alias.name for node in ast.walk(source) if isinstance(node, ast.Import) for alias in node.names]
        self.assertNotIn("ROOT", imported)


if __name__ == "__main__":
    unittest.main()
