"""Focused contract coverage for the detached F.6.2 farm-review bundle."""

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
PROFILE_PATH = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_f6_2.json"
KINEMATIC = "Q4p4W2p74"
REVIEWED_SOURCE = "c88ed65cb18ba6a37358897292b77696016312d1"
SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)

SPEC = importlib.util.spec_from_file_location("_f6_2_bundle_collector", COLLECTOR_PATH)
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


class F62ValidationBundleProfileTests(unittest.TestCase):
    def _declaration(self, profile, scope, key):
        return next(item for item in profile["artifacts"][scope] if item["key"] == key)

    def _artifact_path(self, root, declaration, *, phi="global", epsilon="global"):
        return Path(root) / collector._format_profile_basename(
            declaration["basename_template"], phi, KINEMATIC, epsilon,
        )

    def _write_inputs(self, root):
        root = Path(root)
        profile = collector.load_validation_profile(PROFILE_PATH)
        paths = {}
        for declaration in profile["artifacts"]["global"]:
            path = self._artifact_path(root, declaration)
            if declaration["kind"] == "json":
                path.write_text(json.dumps({"artifact": declaration["key"]}), encoding="utf-8")
            else:
                path.write_bytes(b"f6_2-review-pdf")
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
        output = Path(temporary) / "f6_2-validation.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic=KINEMATIC, output=output,
            profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
            command_runner=command_runner,
        )
        return result, output, profile, paths

    def test_profile_declares_exact_f6_2_review_package(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_f6_2_method_a_acceptance_refinement_validation_farm_review/v1")
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(tuple((item["phi"], item["epsilon"]) for item in profile["settings"]), SETTINGS)
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], REVIEWED_SOURCE)
        self.assertEqual(profile["source_identity"]["allowed_committed_files"], [
            "testing/collect_pion_hgcer_validation_bundle.py",
            "testing/test_collect_pion_hgcer_validation_bundle.py",
            "testing/pion_hgcer_validation_bundle_profile_f6_2.json",
            "testing/test_pion_hgcer_validation_bundle_profile_f6_2.py",
        ])
        self.assertEqual(profile["source_identity"]["allowed_non_analysis_path_prefixes"], ["docs/memory/"])
        self.assertEqual(profile["artifacts"]["global"], [
            {"key": "f6_2_acceptance_refinement_validation_json", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json", "kind": "json", "required": True},
            {"key": "f6_2_acceptance_refinement_validation_pdf", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.pdf", "kind": "file", "required": True},
            {"key": "f6_1_reweighting_validation_json", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-reweighting-validation.json", "kind": "json", "required": True},
            {"key": "f5_tphi_propagation_json", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-tphi-propagation.json", "kind": "json", "required": True},
            {"key": "f4_parent_preserving_correction_json", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json", "kind": "json", "required": True},
            {"key": "f3_acceptance_map_json", "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-map.json", "kind": "json", "required": True},
        ])
        self.assertEqual(profile["artifacts"]["settings"], [{
            "key": "method_a_acceptance_contract",
            "basename_template": "{phi}_kaon_pion-background_hgcer_method-a-acceptance-contract_{kinematic}_{epsilon}.json",
            "kind": "json", "required": True,
        }])

    def test_complete_synthetic_package_is_atomic_and_reports_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, profile, paths = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            manifest = result["manifest"]
            self.assertTrue(output.exists()); self.assertTrue(manifest["complete"])
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(manifest["unexpected_committed_files_after_required_analysis_commit"], [])
            self.assertEqual(set(manifest["global_artifacts"]), {item["key"] for item in profile["artifacts"]["global"]})
            self.assertEqual(len(manifest["settings"]), 5)
            for declaration in profile["artifacts"]["global"]:
                record = manifest["global_artifacts"][declaration["key"]]
                self.assertEqual(record["status"], "exists")
                if declaration["kind"] == "json":
                    self.assertEqual(record["json_status"], "valid")
            for setting in manifest["settings"]:
                record = setting["artifacts"]["method_a_acceptance_contract"]
                self.assertEqual(record["status"], "exists"); self.assertEqual(record["json_status"], "valid")
            with zipfile.ZipFile(output) as archive:
                names = archive.namelist()
                self.assertIn("manifest.json", names); self.assertIn("source_state.txt", names); self.assertIn("source_checks.txt", names)
                for path in paths.values():
                    self.assertEqual(sum(name.endswith(path.name) for name in names), 1)
                for key, path in paths.items():
                    if key in manifest["global_artifacts"]:
                        archive_path = manifest["global_artifacts"][key]["archive_path"]
                    else:
                        setting = next(item for item in manifest["settings"] if "{}-{}".format(item["phi"], item["epsilon"]) == key)
                        archive_path = setting["artifacts"]["method_a_acceptance_contract"]["archive_path"]
                    self.assertEqual(archive.read(archive_path), path.read_bytes())

    def test_declared_artifact_failures_are_best_effort_and_incomplete(self):
        cases = (
            ("missing_f6_2_json", lambda paths: paths["f6_2_acceptance_refinement_validation_json"].unlink(), "global", "f6_2_acceptance_refinement_validation_json", "missing"),
            ("missing_f6_2_pdf", lambda paths: paths["f6_2_acceptance_refinement_validation_pdf"].unlink(), "global", "f6_2_acceptance_refinement_validation_pdf", "missing"),
            ("invalid_f6_2_json", lambda paths: paths["f6_2_acceptance_refinement_validation_json"].write_bytes(b"not-json"), "global", "f6_2_acceptance_refinement_validation_json", "invalid"),
            ("missing_f6_1_json", lambda paths: paths["f6_1_reweighting_validation_json"].unlink(), "global", "f6_1_reweighting_validation_json", "missing"),
            ("missing_f5_json", lambda paths: paths["f5_tphi_propagation_json"].unlink(), "global", "f5_tphi_propagation_json", "missing"),
            ("missing_f4_json", lambda paths: paths["f4_parent_preserving_correction_json"].unlink(), "global", "f4_parent_preserving_correction_json", "missing"),
            ("missing_f3_json", lambda paths: paths["f3_acceptance_map_json"].unlink(), "global", "f3_acceptance_map_json", "missing"),
            ("missing_f1_json", lambda paths: paths["Left-lowe"].unlink(), "setting", "method_a_acceptance_contract", "missing"),
            ("invalid_f1_json", lambda paths: paths["Left-lowe"].write_bytes(b"not-json"), "setting", "method_a_acceptance_contract", "invalid"),
        )
        for label, mutation, scope, key, status in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, mutate=mutation)
                self.assertEqual(result["returncode"], 1); self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                record = result["manifest"]["global_artifacts"][key] if scope == "global" else result["manifest"]["settings"][0]["artifacts"][key]
                self.assertEqual(record["status"], status)

    def test_source_identity_and_output_collision_fail_closed(self):
        def missing_ancestor(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "merge-base", "--is-ancestor"]:
                result["returncode"] = 1
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _profile, _paths = self._collect(temporary, command_runner=missing_ancestor)
            self.assertEqual(result["returncode"], 1); self.assertTrue(output.exists())
            self.assertFalse(result["manifest"]["complete"])
            self.assertIn("required_analysis_commit_not_present", {entry["code"] for entry in result["manifest"]["errors"]})
        for changed_path in (
            "src/cuts/unreviewed_f6_2_change.py",
            "testing/analyze_pion_hgcer_method_a_acceptance_refinement_validation.py",
            "testing/unrelated_validation_test.py",
        ):
            def unexpected_change(command, cwd, path=changed_path):
                result = _clean_command_runner(command, cwd)
                if list(command)[:3] == ["git", "diff", "--name-only"]:
                    result["stdout"] = path + "\n"
                return result

            with self.subTest(path=changed_path), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, command_runner=unexpected_change)
                self.assertEqual(result["returncode"], 1); self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                self.assertEqual(result["manifest"]["unexpected_committed_files_after_required_analysis_commit"], [changed_path])
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"; source.mkdir(); self._write_inputs(source)
            output = Path(temporary) / "f6_2-validation.zip"; output.write_bytes(b"existing")
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(outdir=source, kinematic=KINEMATIC, output=output, profile_path=PROFILE_PATH, repo_root=REPO_ROOT, command_runner=_clean_command_runner)


if __name__ == "__main__":
    unittest.main()
