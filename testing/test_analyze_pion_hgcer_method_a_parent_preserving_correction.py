"""CLI coverage for detached Phase F.4 parent-preserving correction."""

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
import test_pion_hgcer_method_a_parent_preserving_correction as fixtures


KINEMATIC = "Q4p4W2p74"
SCRIPT = REPO_ROOT / "testing" / "analyze_pion_hgcer_method_a_parent_preserving_correction.py"


def _sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()


class AnalyzeParentPreservingCorrectionTests(unittest.TestCase):
    def _write_inputs(self, directory: Path) -> tuple[Path, Path]:
        artifacts = fixtures._artifacts(); hashes: dict[str, str] = {}
        for artifact in artifacts:
            setting = artifact["setting"]; phi, epsilon = setting["phi_setting"], setting["epsilon_filename_token"]
            path = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, KINEMATIC, epsilon)
            path.write_text(json.dumps(artifact, sort_keys=True, indent=2) + "\n", encoding="utf-8"); hashes["{}-{}".format(phi, epsilon)] = _sha256(path)
        f3 = fixtures._f3(artifacts, hashes); f3_path = directory / "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(KINEMATIC)
        f3_path.write_text(json.dumps(f3, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        return f3_path, directory / "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json".format(KINEMATIC)

    def _run(self, directory: Path, output_json: Path | None = None, output_pdf: Path | None = None) -> subprocess.CompletedProcess[str]:
        command = [sys.executable, str(SCRIPT), "--outdir", str(directory), "--kinematic", KINEMATIC]
        if output_json is not None: command += ["--output-json", str(output_json)]
        if output_pdf is not None: command += ["--output-pdf", str(output_pdf)]
        return subprocess.run(command, cwd=REPO_ROOT, text=True, capture_output=True, check=False)

    def test_cli_writes_deterministic_json_and_seven_page_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); _, output_json = self._write_inputs(directory)
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.pdf".format(KINEMATIC)
            completed = self._run(directory)
            self.assertEqual(completed.returncode, 0, completed.stderr); self.assertTrue(output_json.is_file()); self.assertTrue(output_pdf.is_file())
            self.assertTrue(json.loads(output_json.read_text(encoding="utf-8"))["correction"]["available"])
            self.assertEqual(len(re.findall(rb"/Type /Page(?!s)", output_pdf.read_bytes())), 7)

    def test_cli_rejects_wrong_destination_and_collision(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); _, output_json = self._write_inputs(directory)
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.pdf".format(KINEMATIC)
            self.assertEqual(self._run(directory, directory / "wrong.json", output_pdf).returncode, 2)
            output_json.write_text("collision", encoding="utf-8"); self.assertEqual(self._run(directory, output_json, output_pdf).returncode, 2)

    def test_cli_rejects_missing_f3_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); f3_path, output_json = self._write_inputs(directory); f3_path.unlink()
            output_pdf = directory / "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.pdf".format(KINEMATIC)
            completed = self._run(directory, output_json, output_pdf)
            self.assertEqual(completed.returncode, 1); self.assertIn("f4_f3_input_missing", completed.stderr); self.assertFalse(output_json.exists()); self.assertFalse(output_pdf.exists())


if __name__ == "__main__": unittest.main()
