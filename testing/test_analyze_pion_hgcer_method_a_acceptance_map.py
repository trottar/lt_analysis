"""CLI coverage for the detached Phase F.3 acceptance-map analyzer."""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "cuts"), str(REPO_ROOT / "testing")]

import pion_hgcer_method_a_acceptance_representation as representation
import test_pion_hgcer_method_a_acceptance_representation as f2_fixtures


KINEMATIC = "Q4p4W2p74"
SCRIPT = REPO_ROOT / "testing" / "analyze_pion_hgcer_method_a_acceptance_map.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _accepted_f2(artifacts: list[dict[str, object]], hashes: dict[str, str]) -> dict[str, object]:
    artifact = representation.build_pion_hgcer_method_a_acceptance_representation_artifact(artifacts, input_file_hashes=hashes, input_paths={"f1": {}})
    result = artifact["representation"]
    for summary in result["candidate_summaries"]:
        if summary["candidate_id"] in ("delta_only", "track3"):
            summary["information_gate_passed"] = False; summary["overall_candidate_passed"] = False
    result["recommendation"] = {"recommendation_status": "unique_supported_reduced_basis", "recommended_basis": "hgcer3", "basis_frozen": False, "manual_review_required": True}
    result["fingerprint_inputs"]["candidate_summaries"] = result["candidate_summaries"]
    result["fingerprint_inputs"]["recommendation"] = result["recommendation"]
    result["fingerprint"] = representation._sha256(result["fingerprint_inputs"])
    artifact["artifact_fingerprint"] = representation._sha256({"schema_version": artifact["schema_version"], "representation_fingerprint": result["fingerprint"], "input_paths": artifact["provenance"]["input_paths"]})
    return artifact


def _f1_filename(phi: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, KINEMATIC, epsilon)


class AnalyzeMethodAAcceptanceMapTests(unittest.TestCase):
    def _write_inputs(self, directory: Path) -> tuple[Path, Path]:
        artifacts = f2_fixtures.make_five_artifacts()
        hashes: dict[str, str] = {}
        for artifact in artifacts:
            setting = artifact["setting"]
            phi, epsilon = setting["phi_setting"], setting["epsilon_filename_token"]
            path = directory / _f1_filename(phi, epsilon)
            path.write_text(json.dumps(artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            hashes["{}-{}".format(phi, epsilon)] = _sha256(path)
        f2 = _accepted_f2(artifacts, hashes)
        f2_path = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-representation.json".format(KINEMATIC)
        f2_path.write_text(json.dumps(f2, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        return f2_path, directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(KINEMATIC)

    def _run(self, directory: Path, output_json: Path | None = None, output_pdf: Path | None = None) -> subprocess.CompletedProcess[str]:
        command = [sys.executable, str(SCRIPT), "--outdir", str(directory), "--kinematic", KINEMATIC]
        if output_json is not None:
            command.extend(("--output-json", str(output_json)))
        if output_pdf is not None:
            command.extend(("--output-pdf", str(output_pdf)))
        return subprocess.run(command, cwd=REPO_ROOT, text=True, capture_output=True, check=False)

    def test_cli_creates_deterministic_json_and_seven_page_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, output_json = self._write_inputs(directory)
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.pdf".format(KINEMATIC)
            completed = self._run(directory, output_json, output_pdf)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(output_json.is_file()); self.assertTrue(output_pdf.is_file())
            payload = json.loads(output_json.read_text(encoding="utf-8"))
            self.assertTrue(payload["acceptance_map"]["available"])
            # Matplotlib emits one /Type /Page object per rendered review page.
            self.assertEqual(len(re.findall(rb"/Type /Page(?!s)", output_pdf.read_bytes())), 7)

    def test_cli_defaults_outputs_under_outdir(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, output_json = self._write_inputs(directory)
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.pdf".format(KINEMATIC)
            completed = self._run(directory)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(output_json.is_file()); self.assertTrue(output_pdf.is_file())

    def test_cli_rejects_non_deterministic_or_colliding_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, output_json = self._write_inputs(directory)
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.pdf".format(KINEMATIC)
            wrong = directory / "wrong.json"
            self.assertEqual(self._run(directory, wrong, output_pdf).returncode, 2)
            output_json.write_text("collision", encoding="utf-8")
            self.assertEqual(self._run(directory, output_json, output_pdf).returncode, 2)

    def test_cli_rejects_missing_f2_input_without_writing_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, output_json = self._write_inputs(directory)
            (directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-representation.json".format(KINEMATIC)).unlink()
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.pdf".format(KINEMATIC)
            completed = self._run(directory, output_json, output_pdf)
            self.assertEqual(completed.returncode, 1)
            self.assertIn("f3_f2_input_missing", completed.stderr)
            self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())


if __name__ == "__main__":
    unittest.main()
