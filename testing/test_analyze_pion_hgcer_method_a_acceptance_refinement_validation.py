"""CLI coverage for detached Phase F.6.2 refinement validation."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import inspect
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import analyze_pion_hgcer_method_a_acceptance_refinement_validation as analyzer
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_reweighting_validation as f61
import pion_hgcer_method_a_tphi_propagation as f5
import test_pion_hgcer_method_a_reweighting_validation as fixtures


KINEMATIC = fixtures.KINEMATIC


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalyzeAcceptanceRefinementValidationTests(unittest.TestCase):
    def _write_inputs(self, directory: Path) -> tuple[Path, Path, dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
        artifacts = fixtures._enriched_artifacts(); hashes: dict[str, str] = {}; paths: dict[str, str] = {}
        for artifact in artifacts:
            setting = artifact["setting"]; phi, epsilon = setting["phi_setting"], setting["epsilon_filename_token"]
            path = directory / analyzer._f1_filename(phi, KINEMATIC, epsilon)
            path.write_text(json.dumps(artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            key = "{}-{}".format(phi, epsilon); hashes[key] = _sha256(path); paths[key] = str(path.resolve())
        f3_payload = fixtures.f4_fixtures._f3(artifacts, hashes); f3_path = directory / analyzer._f3_filename(KINEMATIC)
        f3_path.write_text(json.dumps(f3_payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"); f3_sha = _sha256(f3_path)
        f3_authority = fixtures.f4_fixtures._authority(f3_payload, f3_sha)
        f4_payload = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, f3_payload, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve())}, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f4_path = directory / analyzer._f4_filename(KINEMATIC)
        f4_path.write_text(json.dumps(f4_payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"); f4_sha = _sha256(f4_path)
        correction = f4_payload["correction"]; assert isinstance(correction, dict)
        f4_authority = {KINEMATIC: {"source_file_sha256": f4_sha, "correction_fingerprint": correction["fingerprint"], "artifact_fingerprint": f4_payload["artifact_fingerprint"], "farm_source_head": "2" * 40}}
        f5_payload = f5.build_pion_hgcer_method_a_tphi_propagation_artifact(artifacts, f3_payload, f4_payload, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve()), "f4": str(f4_path.resolve())}, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f5_path = directory / analyzer._f5_filename(KINEMATIC)
        f5_path.write_text(json.dumps(f5_payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"); f5_sha = _sha256(f5_path)
        propagation = f5_payload["propagation"]; assert isinstance(propagation, dict)
        authority = {KINEMATIC: {"f1_source_file_sha256": hashes, "f3_source_file_sha256": f3_sha, "f3_map_fingerprint": f3_payload["acceptance_map"]["fingerprint"], "f3_algorithm_fingerprint": f3_payload["acceptance_map"]["algorithm_fingerprint"], "f3_artifact_fingerprint": f3_payload["artifact_fingerprint"], "f4_source_file_sha256": f4_sha, "f4_correction_fingerprint": correction["fingerprint"], "f4_artifact_fingerprint": f4_payload["artifact_fingerprint"], "f5_source_file_sha256": f5_sha, "f5_propagation_fingerprint": propagation["fingerprint"], "f5_artifact_fingerprint": f5_payload["artifact_fingerprint"], "farm_source_head": "3" * 40}}
        f6_payload = f61.build_pion_hgcer_method_a_reweighting_validation_artifact(artifacts, f3_payload, f4_payload, f5_payload, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, f5_input_file_sha256=f5_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve()), "f4": str(f4_path.resolve()), "f5": str(f5_path.resolve())}, accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f6_path = directory / analyzer._f6_1_filename(KINEMATIC)
        f6_path.write_text(json.dumps(f6_payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"); f6_sha = _sha256(f6_path)
        f6_authority = {KINEMATIC: {"source_file_sha256": f6_sha, "validation_fingerprint": f6_payload["validation"]["fingerprint"], "artifact_fingerprint": f6_payload["artifact_fingerprint"]}}
        return directory / analyzer.validation.pion_hgcer_method_a_acceptance_refinement_validation_filename(KINEMATIC), directory / analyzer._pdf_filename(KINEMATIC), authority, f4_authority, f3_authority, f6_authority

    def _run(self, directory: Path, json_path: Path, pdf_path: Path, authority: dict[str, object], f4_authority: dict[str, object], f3_authority: dict[str, object], f6_authority: dict[str, object], captured_artifacts: list[tuple[dict[str, object], dict[str, object]]] | None = None) -> int:
        real_builder = analyzer.validation.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact

        def build_with_test_bootstrap(*args, **kwargs):
            self.assertNotIn("bootstrap_test_config", kwargs)
            kwargs["bootstrap_test_config"] = {"replicas": 2}
            artifact = real_builder(*args, **kwargs)
            if captured_artifacts is not None:
                captured_artifacts.append((artifact, deepcopy(artifact)))
            return artifact

        with mock.patch.object(analyzer.validation, "build_pion_hgcer_method_a_acceptance_refinement_validation_artifact", side_effect=build_with_test_bootstrap):
            return analyzer.main(["--outdir", str(directory), "--kinematic", KINEMATIC, "--output-json", str(json_path), "--output-pdf", str(pdf_path)], accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority, accepted_f6_1_artifact_authority_by_kinematic=f6_authority)

    def test_production_analyzer_has_no_bootstrap_test_override(self):
        self.assertNotIn("bootstrap_test_config", inspect.signature(analyzer.main).parameters)
        self.assertNotIn("bootstrap_test_config", inspect.getsource(analyzer.main))
        self.assertNotIn("bootstrap_test_config", inspect.getsource(analyzer.build_argument_parser))

    def test_cli_writes_deterministic_aggregate_artifact_and_review_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            captured: list[tuple[dict[str, object], dict[str, object]]] = []
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority, captured), 0)
            payload = json.loads(output_json.read_text(encoding="utf-8")); self.assertTrue(payload["validation"]["available"])
            self.assertEqual(len(captured), 1); self.assertEqual(captured[0][0], captured[0][1])
            details = sum(1 for parent in payload["validation"]["parents"] for child in parent["children"] if any(child["population_counts"][name] > 0 for name in ("N_low", "N_control", "N_full_application")))
            self.assertEqual(len(re.findall(rb"/Type /Page(?!s)", output_pdf.read_bytes())), 1 + 15 + details + 1)
            serialized = json.dumps(payload, sort_keys=True); self.assertNotIn("correction_factors", serialized); self.assertNotIn("entry_index", serialized)
            self.assertEqual(len(analyzer._OVERVIEW_HEADINGS), 11); self.assertEqual(len(analyzer._OVERVIEW_COLUMN_WIDTHS), 11)
            self.assertEqual(analyzer._OVERVIEW_HEADINGS, ("phi\nrange", "N low/\ncontrol/full", "N_eff\nB/A", "OOD frac\nctrl/full", "DeltaH MM\nvalue / 95% CI", "R MM", "kappa MM\nvalue / 95% CI", "kappa MM×xptar\nvalue / 95% CI", "kappa MM×yptar\nvalue / 95% CI", "DeltaP K\nvalue / 95% CI", "R_V"))
            self.assertTrue(all(sum(1 for child in parent["children"]) == 9 for parent in payload["validation"]["parents"]))
            self.assertTrue(all(len(row) == 11 for parent in payload["validation"]["parents"] for row in analyzer._overview_rows(parent)))
            interval_cell = analyzer._overview_metric_with_interval({"available": True, "value": -0.5258}, {"interval_available": True, "ci_low": -0.6660, "ci_high": -0.42})
            self.assertEqual(interval_cell, "-0.5258\n[-0.666, -0.42]"); self.assertNotIn("; 95% CI", interval_cell)
            unavailable_cell = analyzer._overview_metric_with_interval({"available": True, "value": 0.25}, {"interval_available": False, "reason": "insufficient_valid_bootstrap_replicas"})
            self.assertEqual(unavailable_cell, "0.25\nunavailable:insufficient_valid_bootstrap_replicas")

    def test_cli_rejects_paths_overwrite_and_missing_required_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            self.assertEqual(self._run(directory, directory / "wrong.json", output_pdf, authority, f4_authority, f3_authority, f6_authority), 2)
            self.assertEqual(self._run(directory, output_json, directory / "wrong.pdf", authority, f4_authority, f3_authority, f6_authority), 2)
            output_json.write_text("collision", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 2)
            output_json.unlink(); output_pdf.write_text("collision", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 2)
        for filename in (
            analyzer._f1_filename("Left", KINEMATIC, "lowe"),
            analyzer._f3_filename(KINEMATIC), analyzer._f4_filename(KINEMATIC),
            analyzer._f5_filename(KINEMATIC), analyzer._f6_1_filename(KINEMATIC),
        ):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
                (directory / filename).unlink()
                self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)

    def test_filenames_invalid_inputs_and_f6_authority_are_fail_closed(self):
        self.assertEqual(analyzer._f1_filename("Left", KINEMATIC, "lowe"), "Left_kaon_pion-background_hgcer_method-a-acceptance-contract_Q4p4W2p74_lowe.json")
        self.assertEqual(analyzer.validation.pion_hgcer_method_a_acceptance_refinement_validation_filename(KINEMATIC), "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            wrong = {KINEMATIC: dict(f6_authority[KINEMATIC])}; wrong[KINEMATIC]["source_file_sha256"] = "f" * 64
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, wrong), 1)
            self.assertFalse(any(path.name.startswith(".f6_2_") for path in directory.iterdir()))
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            missing = directory / "missing"
            self.assertEqual(analyzer.main(["--outdir", str(missing), "--kinematic", KINEMATIC, "--output-json", str(missing / analyzer.validation.pion_hgcer_method_a_acceptance_refinement_validation_filename(KINEMATIC)), "--output-pdf", str(missing / analyzer._pdf_filename(KINEMATIC))]), 1)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            (directory / analyzer._f3_filename(KINEMATIC)).write_text("{", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            f1_path = directory / analyzer._f1_filename("Left", KINEMATIC, "lowe")
            payload = json.loads(f1_path.read_text(encoding="utf-8")); payload["setting"]["phi_setting"] = "Wrong"
            f1_path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            (directory / analyzer._f1_filename("Left", KINEMATIC, "lowe")).write_text("{", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)

    def test_cli_rollback_and_unsupported_kinematic(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            original_replace, calls = analyzer.os.replace, 0
            def fail_second(source, destination):
                nonlocal calls
                calls += 1
                if calls == 2: raise OSError("f6_2_test_second_promotion_failure")
                return original_replace(source, destination)
            with mock.patch.object(analyzer, "write_review_pdf", side_effect=lambda path, artifact: path.write_bytes(b"pdf")), mock.patch.object(analyzer.os, "replace", side_effect=fail_second):
                self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())
            self.assertEqual(analyzer.main(["--outdir", str(directory), "--kinematic", "Q9p9W9p99", "--output-json", str(directory / "x.json"), "--output-pdf", str(directory / "x.pdf")]), 1)

    def test_pdf_generation_and_first_promotion_failures_roll_back(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            with mock.patch.object(analyzer, "write_review_pdf", side_effect=OSError("f6_2_pdf_generation_failure")):
                self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists()); self.assertFalse(any(path.name.startswith(".f6_2_") for path in directory.iterdir()))
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority = self._write_inputs(directory)
            with mock.patch.object(analyzer, "write_review_pdf", side_effect=lambda path, artifact: path.write_bytes(b"pdf")), mock.patch.object(analyzer.os, "replace", side_effect=OSError("f6_2_first_promotion_failure")):
                self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority, f6_authority), 1)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())


if __name__ == "__main__":
    unittest.main()
