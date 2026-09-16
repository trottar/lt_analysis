"""CLI coverage for detached Phase F.6.1 Method-A reweighting validation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import analyze_pion_hgcer_method_a_reweighting_validation as analyzer
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5
import test_pion_hgcer_method_a_reweighting_validation as fixtures


KINEMATIC = fixtures.KINEMATIC


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalyzeReweightingValidationTests(unittest.TestCase):
    def _write_inputs(self, directory: Path) -> tuple[Path, Path, dict[str, object], dict[str, object], dict[str, object]]:
        artifacts = fixtures._enriched_artifacts(); hashes: dict[str, str] = {}; paths: dict[str, str] = {}
        for artifact in artifacts:
            setting = artifact["setting"]; phi, epsilon = setting["phi_setting"], setting["epsilon_filename_token"]
            path = directory / analyzer._f1_filename(phi, KINEMATIC, epsilon)
            path.write_text(json.dumps(artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            key = "{}-{}".format(phi, epsilon); hashes[key] = _sha256(path); paths[key] = str(path.resolve())
        f3 = fixtures.f4_fixtures._f3(artifacts, hashes); f3_path = directory / analyzer._f3_filename(KINEMATIC)
        f3_path.write_text(json.dumps(f3, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        f3_sha = _sha256(f3_path); f3_authority = fixtures.f4_fixtures._authority(f3, f3_sha)
        f4_artifact = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, f3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve())}, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f4_path = directory / analyzer._f4_filename(KINEMATIC)
        f4_path.write_text(json.dumps(f4_artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        correction = f4_artifact["correction"]; assert isinstance(correction, dict)
        f4_sha = _sha256(f4_path); f4_authority = {KINEMATIC: {"source_file_sha256": f4_sha, "correction_fingerprint": correction["fingerprint"], "artifact_fingerprint": f4_artifact["artifact_fingerprint"], "farm_source_head": "2" * 40}}
        f5_artifact = f5.build_pion_hgcer_method_a_tphi_propagation_artifact(artifacts, f3, f4_artifact, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve()), "f4": str(f4_path.resolve())}, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f5_path = directory / analyzer._f5_filename(KINEMATIC)
        f5_path.write_text(json.dumps(f5_artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        propagation = f5_artifact["propagation"]; assert isinstance(propagation, dict)
        f5_sha = _sha256(f5_path)
        authority = {KINEMATIC: {"f1_source_file_sha256": hashes, "f3_source_file_sha256": f3_sha, "f3_map_fingerprint": f3["acceptance_map"]["fingerprint"], "f3_algorithm_fingerprint": f3["acceptance_map"]["algorithm_fingerprint"], "f3_artifact_fingerprint": f3["artifact_fingerprint"], "f4_source_file_sha256": f4_sha, "f4_correction_fingerprint": correction["fingerprint"], "f4_artifact_fingerprint": f4_artifact["artifact_fingerprint"], "f5_source_file_sha256": f5_sha, "f5_propagation_fingerprint": propagation["fingerprint"], "f5_artifact_fingerprint": f5_artifact["artifact_fingerprint"], "farm_source_head": "3" * 40}}
        return directory / analyzer.validation.pion_hgcer_method_a_reweighting_validation_filename(KINEMATIC), directory / analyzer._pdf_filename(KINEMATIC), authority, f4_authority, f3_authority

    def _run(self, directory: Path, json_path: Path, pdf_path: Path, authority: dict[str, object], f4_authority: dict[str, object], f3_authority: dict[str, object]) -> int:
        return analyzer.main(["--outdir", str(directory), "--kinematic", KINEMATIC, "--output-json", str(json_path), "--output-pdf", str(pdf_path)], accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)

    def test_cli_writes_deterministic_aggregate_json_and_twenty_seven_page_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority = self._write_inputs(directory)
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority), 0)
            self.assertTrue(output_json.is_file()); self.assertTrue(output_pdf.is_file())
            payload = json.loads(output_json.read_text(encoding="utf-8")); self.assertTrue(payload["validation"]["available"])
            self.assertEqual(len(re.findall(rb"/Type /Page(?!s)", output_pdf.read_bytes())), 27)
            serialized = json.dumps(payload, sort_keys=True)
            self.assertNotIn("correction_factors", serialized); self.assertNotIn("entry_index", serialized)

    def test_cli_rejects_wrong_paths_collisions_and_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority = self._write_inputs(directory)
            self.assertEqual(self._run(directory, directory / "wrong.json", output_pdf, authority, f4_authority, f3_authority), 2)
            output_json.write_text("collision", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority), 2)

    def test_cli_rejects_unsupported_kinematic(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority = self._write_inputs(directory)
            unsupported = "Q9p9W9p99"
            result = analyzer.main(["--outdir", str(directory), "--kinematic", unsupported, "--output-json", str(directory / analyzer.validation.pion_hgcer_method_a_reweighting_validation_filename(unsupported)), "--output-pdf", str(directory / analyzer._pdf_filename(unsupported))], accepted_runtime_authority_by_kinematic=authority, accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)
            self.assertEqual(result, 1)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())

    def test_late_json_promotion_failure_rolls_back_already_promoted_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); output_json, output_pdf, authority, f4_authority, f3_authority = self._write_inputs(directory)
            original_replace = analyzer.os.replace; calls = 0
            def fail_second_replace(source, destination):
                nonlocal calls
                calls += 1
                if calls == 2: raise OSError("f6_1_test_second_promotion_failed")
                return original_replace(source, destination)
            with mock.patch.object(analyzer, "write_review_pdf", side_effect=lambda path, artifact: path.write_bytes(b"pdf")), mock.patch.object(analyzer.os, "replace", side_effect=fail_second_replace):
                self.assertEqual(self._run(directory, output_json, output_pdf, authority, f4_authority, f3_authority), 1)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())

    def test_parser_requires_all_output_arguments(self):
        with self.assertRaises(SystemExit) as captured:
            analyzer.main(["--outdir", "anywhere", "--kinematic", KINEMATIC, "--output-json", "one.json"])
        self.assertEqual(captured.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
