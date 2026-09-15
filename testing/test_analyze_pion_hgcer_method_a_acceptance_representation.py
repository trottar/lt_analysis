"""CLI tests for the detached F.2 representation analyzer."""

from __future__ import annotations

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_PATH = REPO_ROOT / "testing" / "analyze_pion_hgcer_method_a_acceptance_representation.py"
sys.path.insert(0, str(REPO_ROOT / "testing"))
from test_pion_hgcer_method_a_acceptance_representation import SETTINGS, make_five_artifacts


SPEC = importlib.util.spec_from_file_location("_f2_representation_cli", CLI_PATH)
cli = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cli
SPEC.loader.exec_module(cli)


def _write_inputs(directory: Path, artifacts: list[dict[str, object]] | None = None) -> list[dict[str, object]]:
    payloads = make_five_artifacts() if artifacts is None else artifacts
    for artifact in payloads:
        setting = artifact["setting"]
        path = directory / cli._f1_filename(
            setting["phi_setting"], setting["kinematic_token"], setting["epsilon_filename_token"],
        )
        path.write_text(json.dumps(artifact, sort_keys=True), encoding="utf-8")
    return payloads


class AnalyzeMethodAAcceptanceRepresentationTests(unittest.TestCase):
    def _arguments(self, source: Path, target: Path, pdf: Path) -> list[str]:
        return [
            "--outdir", str(source), "--kinematic", "Q4p4W2p74",
            "--output-json", str(target), "--output-pdf", str(pdf),
        ]

    def test_valid_five_file_directory_writes_json_and_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_inputs(source)
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            pdf = root / "explicit-review.pdf"
            self.assertEqual(cli.main(self._arguments(source, target, pdf)), 0)
            payload = json.loads(target.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema_version"], "pion_hgcer_method_a_acceptance_representation_artifact/v1")
            self.assertEqual(len(payload["representation"]["input_fingerprints"]), 5)
            self.assertTrue(all(row["source_file_sha256"] for row in payload["representation"]["input_fingerprints"]))
            self.assertGreater(pdf.stat().st_size, 0)
            self.assertNotIn("import ROOT", CLI_PATH.read_text(encoding="utf-8"))

    def test_missing_wrong_schema_and_conflicting_setting_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_inputs(source)
            missing = source / cli._f1_filename("Right", "Q4p4W2p74", "highe")
            missing.unlink()
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            self.assertEqual(cli.main(self._arguments(source, target, root / "missing.pdf")), 1)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            artifacts = make_five_artifacts()
            artifacts[0]["schema_version"] = "pion_hgcer_method_a_acceptance_event_contract_artifact/v1"
            _write_inputs(source, artifacts)
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            self.assertEqual(cli.main(self._arguments(source, target, root / "schema.pdf")), 1)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            artifacts = make_five_artifacts()
            _write_inputs(source, artifacts)
            artifacts[-1]["setting"]["phi_setting"] = "Left"
            artifacts[-1]["setting"]["epsilon_filename_token"] = "lowe"
            right_path = source / cli._f1_filename("Right", "Q4p4W2p74", "highe")
            right_path.write_text(json.dumps(artifacts[-1], sort_keys=True), encoding="utf-8")
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            self.assertEqual(cli.main(self._arguments(source, target, root / "duplicate.pdf")), 1)

    def test_deterministic_filenames_bind_payload_identity_before_hashing(self):
        def assert_identity_failure(mutator, error):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                source = root / "source"
                source.mkdir()
                artifacts = make_five_artifacts()
                _write_inputs(source, artifacts)
                mutator(artifacts)
                for artifact, (phi, epsilon) in zip(artifacts, cli.CANONICAL_SETTINGS):
                    path = source / cli._f1_filename(phi, "Q4p4W2p74", epsilon)
                    path.write_text(json.dumps(artifact, sort_keys=True), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, error):
                    cli.load_f1_inputs(source, "Q4p4W2p74")

        assert_identity_failure(
            lambda artifacts: artifacts.__setitem__(0, deepcopy(artifacts[2])),
            "f2_f1_input_identity_phi_setting:Left-lowe",
        )
        assert_identity_failure(
            lambda artifacts: artifacts[0]["setting"].__setitem__("phi_setting", "Center"),
            "f2_f1_input_identity_phi_setting:Left-lowe",
        )
        assert_identity_failure(
            lambda artifacts: artifacts[0]["setting"].__setitem__("epsilon_filename_token", "highe"),
            "f2_f1_input_identity_epsilon_filename_token:Left-lowe",
        )
        assert_identity_failure(
            lambda artifacts: [artifact["setting"].__setitem__("kinematic_token", "Q3p0W2p32") for artifact in artifacts],
            "f2_f1_input_identity_kinematic_token:Left-lowe",
        )
        assert_identity_failure(
            lambda artifacts: artifacts[0]["setting"].__setitem__("particle_type", "pion"),
            "f2_f1_input_identity_particle_type:Left-lowe",
        )

    def test_collision_and_existing_outputs_fail_before_any_write(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_inputs(source)
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            self.assertEqual(cli.main(self._arguments(source, target, target)), 2)
            self.assertFalse(target.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_inputs(source)
            target = root / cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74")
            target.write_text("pre-existing", encoding="utf-8")
            pdf = root / "review.pdf"
            self.assertEqual(cli.main(self._arguments(source, target, pdf)), 2)
            self.assertEqual(target.read_text(encoding="utf-8"), "pre-existing")
            self.assertFalse(pdf.exists())

    def test_json_basename_is_deterministic_and_pdf_path_is_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_inputs(source)
            self.assertEqual(
                cli.main(self._arguments(source, root / "wrong.json", root / "review.pdf")),
                2,
            )
            self.assertEqual(
                cli.pion_hgcer_method_a_acceptance_representation_filename("Q4p4W2p74"),
                "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-representation.json",
            )


if __name__ == "__main__":
    unittest.main()
