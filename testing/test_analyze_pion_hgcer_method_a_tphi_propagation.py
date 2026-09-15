"""CLI coverage for detached Phase F.5 signed template propagation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import analyze_pion_hgcer_method_a_tphi_propagation as analyzer
import pion_hgcer_method_a_parent_preserving_correction as f4
import test_pion_hgcer_method_a_parent_preserving_correction as f4_fixtures
import test_pion_hgcer_method_a_tphi_propagation as fixtures


KINEMATIC = fixtures.KINEMATIC


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalyzeTphiPropagationTests(unittest.TestCase):
    def _write_inputs(self, directory: Path) -> tuple[Path, Path, dict[str, object], dict[str, object]]:
        artifacts = fixtures._artifacts(); hashes: dict[str, str] = {}; paths: dict[str, str] = {}
        for artifact in artifacts:
            setting = artifact["setting"]; phi, epsilon = setting["phi_setting"], setting["epsilon_filename_token"]
            path = directory / analyzer._f1_filename(phi, KINEMATIC, epsilon)
            path.write_text(json.dumps(artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            key = "{}-{}".format(phi, epsilon); hashes[key] = _sha256(path); paths[key] = str(path.resolve())
        f3 = f4_fixtures._f3(artifacts, hashes); f3_path = directory / analyzer._f3_filename(KINEMATIC)
        f3_path.write_text(json.dumps(f3, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        f3_sha = _sha256(f3_path); f3_authority = f4_fixtures._authority(f3, f3_sha)
        f4_artifact = f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(artifacts, f3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={"f1": paths, "f3": str(f3_path.resolve())}, accepted_f3_runtime_authority_by_kinematic=f3_authority)
        f4_path = directory / analyzer._f4_filename(KINEMATIC)
        f4_path.write_text(json.dumps(f4_artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        correction = f4_artifact["correction"]
        assert isinstance(correction, dict)
        f4_authority = {KINEMATIC: {"source_file_sha256": _sha256(f4_path), "correction_fingerprint": correction["fingerprint"], "artifact_fingerprint": f4_artifact["artifact_fingerprint"], "farm_source_head": "2" * 40}}
        return directory / analyzer.propagation.pion_hgcer_method_a_tphi_propagation_filename(KINEMATIC), directory / analyzer._pdf_filename(KINEMATIC), f4_authority, f3_authority

    def _run(self, directory: Path, json_path: Path, pdf_path: Path, f4_authority: dict[str, object], f3_authority: dict[str, object]) -> int:
        return analyzer.main(["--outdir", str(directory), "--kinematic", KINEMATIC, "--output-json", str(json_path), "--output-pdf", str(pdf_path)], accepted_f4_runtime_authority_by_kinematic=f4_authority, accepted_f3_runtime_authority_by_kinematic=f3_authority)

    def test_cli_writes_exact_named_json_and_twelve_page_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); _, output_pdf, f4_authority, f3_authority = self._write_inputs(directory)
            output_json = directory / "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.json".format(KINEMATIC)
            self.assertEqual(self._run(directory, output_json, output_pdf, f4_authority, f3_authority), 0)
            self.assertTrue(output_json.is_file()); self.assertTrue(output_pdf.is_file())
            payload = json.loads(output_json.read_text(encoding="utf-8")); self.assertTrue(payload["propagation"]["available"])
            self.assertEqual(len(re.findall(rb"/Type /Page(?!s)", output_pdf.read_bytes())), 12)

    def test_cli_requires_exact_paths_and_rejects_collision(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); _, output_pdf, f4_authority, f3_authority = self._write_inputs(directory)
            output_json = directory / "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.json".format(KINEMATIC)
            self.assertEqual(self._run(directory, directory / "wrong.json", output_pdf, f4_authority, f3_authority), 2)
            output_json.write_text("collision", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf, f4_authority, f3_authority), 2)

    def test_cli_missing_required_output_option_exits_with_parser_error(self):
        with self.assertRaises(SystemExit) as captured:
            analyzer.main(["--outdir", "anywhere", "--kinematic", KINEMATIC, "--output-json", "one.json"])
        self.assertEqual(captured.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
