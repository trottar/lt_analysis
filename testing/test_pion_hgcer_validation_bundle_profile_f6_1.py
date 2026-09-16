"""Focused contract coverage for the detached F.6.1 farm-review bundle."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
COLLECTOR_PATH = REPO_ROOT / "testing" / "collect_pion_hgcer_validation_bundle.py"
PROFILE_PATH = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_f6_1.json"
KINEMATIC = "Q4p4W2p74"
REVIEWED_SOURCE = "bfc4fe421f9fc9139a1992a0ec92e31aa101b86c"
SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)

SPEC = importlib.util.spec_from_file_location("_f6_1_bundle_collector", COLLECTOR_PATH)
collector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = collector
assert SPEC.loader is not None
SPEC.loader.exec_module(collector)


def _clean_command_runner(command, _cwd):
    command = [str(item) for item in command]
    stdout = ""
    if command == ["git", "rev-parse", "HEAD"]:
        stdout = "test-head\n"
    elif command == ["git", "rev-parse", "HEAD^"]:
        stdout = "test-parent\n"
    elif command[-1:] == ["--version"]:
        stdout = "Python fake\n"
    return {"command": command, "returncode": 0, "stdout": stdout, "stderr": ""}


class F61ValidationBundleProfileTests(unittest.TestCase):
    def _declaration(self, profile, scope, key):
        return next(item for item in profile["artifacts"][scope] if item["key"] == key)

    def _artifact_path(self, root, declaration, *, phi="global", epsilon="global"):
        basename = collector._format_profile_basename(
            declaration["basename_template"], phi, KINEMATIC, epsilon,
        )
        return Path(root) / basename

    def _write_inputs(self, root):
        root = Path(root)
        profile = collector.load_validation_profile(PROFILE_PATH)
        paths = {}
        for declaration in profile["artifacts"]["global"]:
            path = self._artifact_path(root, declaration)
            if declaration["kind"] == "json":
                path.write_text(json.dumps({"artifact": declaration["key"]}), encoding="utf-8")
            else:
                path.write_bytes(b"f6_1-review-pdf")
            paths[declaration["key"]] = path
        acceptance = self._declaration(profile, "settings", "method_a_acceptance_contract")
        for phi, epsilon in SETTINGS:
            path = self._artifact_path(root, acceptance, phi=phi, epsilon=epsilon)
            path.write_text(json.dumps({"setting": {"phi_setting": phi, "epsilon_setting": epsilon}}), encoding="utf-8")
            paths["{}-{}".format(phi, epsilon)] = path
        return profile, paths

    def _collect(self, temporary, *, command_runner=_clean_command_runner, mutate=None):
        source = Path(temporary) / "source"
        source.mkdir()
        profile, paths = self._write_inputs(source)
        if mutate is not None:
            mutate(paths)
        output = Path(temporary) / "f6_1-validation.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic=KINEMATIC, output=output,
            profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
            command_runner=command_runner,
        )
        return result, output, profile, paths

    def test_profile_declares_only_the_f6_1_review_package(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_f6_1_method_a_reweighting_validation_farm_review/v1")
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(tuple((item["phi"], item["epsilon"]) for item in profile["settings"]), SETTINGS)
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], REVIEWED_SOURCE)
        self.assertEqual(profile["source_identity"]["allowed_committed_files"], [
            "testing/pion_hgcer_validation_bundle_profile_f6_1.json",
            "testing/test_pion_hgcer_validation_bundle_profile_f6_1.py",
        ])
        self.assertEqual(profile["source_identity"]["allowed_non_analysis_path_prefixes"], ["docs/memory/"])
        self.assertEqual([item["basename_template"] for item in profile["artifacts"]["global"]], [
            "{kinematic}_kaon_pion-background_hgcer_method-a-reweighting-validation.json",
            "{kinematic}_kaon_pion-background_hgcer_method-a-reweighting-validation.pdf",
            "{kinematic}_kaon_pion-background_hgcer_method-a-tphi-propagation.json",
            "{kinematic}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json",
            "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-map.json",
        ])
        self.assertEqual(profile["artifacts"]["settings"], [{
            "key": "method_a_acceptance_contract",
            "basename_template": "{phi}_kaon_pion-background_hgcer_method-a-acceptance-contract_{kinematic}_{epsilon}.json",
            "kind": "json",
            "required": True,
        }])

    def test_complete_synthetic_package_is_atomic_and_reports_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, profile, paths = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            manifest = result["manifest"]
            self.assertTrue(output.exists())
            self.assertTrue(manifest["complete"])
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(manifest["unexpected_committed_files_after_required_analysis_commit"], [])
            self.assertEqual(set(manifest["global_artifacts"]), {
                "f6_1_reweighting_validation_json", "f6_1_reweighting_validation_pdf",
                "f5_tphi_propagation_json", "f4_parent_preserving_correction_json",
                "f3_acceptance_map_json",
            })
            self.assertEqual(len(manifest["settings"]), 5)
            for record in manifest["global_artifacts"].values():
                self.assertEqual(record["status"], "exists")
                if record["artifact"].endswith("_json"):
                    self.assertEqual(record["json_status"], "valid")
            for setting in manifest["settings"]:
                record = setting["artifacts"]["method_a_acceptance_contract"]
                self.assertEqual(record["status"], "exists")
                self.assertEqual(record["json_status"], "valid")
            with zipfile.ZipFile(output) as archive:
                names = archive.namelist()
                self.assertIn("manifest.json", names)
                self.assertIn("source_state.txt", names)
                self.assertIn("source_checks.txt", names)
                for path in paths.values():
                    self.assertEqual(sum(name.endswith(path.name) for name in names), 1)
                self.assertEqual(archive.read(manifest["global_artifacts"]["f6_1_reweighting_validation_json"]["archive_path"]), paths["f6_1_reweighting_validation_json"].read_bytes())

    def test_declared_artifact_failures_are_complete_false_and_best_effort(self):
        cases = (
            ("missing_f6_json", lambda paths: paths["f6_1_reweighting_validation_json"].unlink(), "f6_1_reweighting_validation_json", "missing"),
            ("missing_f6_pdf", lambda paths: paths["f6_1_reweighting_validation_pdf"].unlink(), "f6_1_reweighting_validation_pdf", "missing"),
            ("missing_f5_json", lambda paths: paths["f5_tphi_propagation_json"].unlink(), "f5_tphi_propagation_json", "missing"),
            ("missing_f1_json", lambda paths: paths["Left-lowe"].unlink(), "method_a_acceptance_contract", "missing"),
            ("invalid_f6_json", lambda paths: paths["f6_1_reweighting_validation_json"].write_bytes(b"not-json"), "f6_1_reweighting_validation_json", "invalid"),
        )
        for label, mutation, key, status in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, mutate=mutation)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                if key == "method_a_acceptance_contract":
                    record = result["manifest"]["settings"][0]["artifacts"][key]
                else:
                    record = result["manifest"]["global_artifacts"][key]
                self.assertEqual(record["status"], status)

    def test_source_identity_and_output_collision_fail_closed(self):
        def missing_ancestor(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "merge-base", "--is-ancestor"]:
                result["returncode"] = 1
            return result

        def unexpected_science(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result["stdout"] = "src/cuts/f6_1_unreviewed_science.py\n"
            return result

        for label, runner, error_code in (
            ("required_analysis_commit_not_ancestor", missing_ancestor, "required_analysis_commit_not_present"),
            ("unexpected_committed_science", unexpected_science, "unexpected_committed_files_after_required_analysis_commit"),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, command_runner=runner)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                self.assertIn(error_code, {entry["code"] for entry in result["manifest"]["errors"]})
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            self._write_inputs(source)
            output = Path(temporary) / "f6_1-validation.zip"
            output.write_bytes(b"existing")
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(
                    outdir=source, kinematic=KINEMATIC, output=output,
                    profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
                    command_runner=_clean_command_runner,
                )


if __name__ == "__main__":
    unittest.main()
