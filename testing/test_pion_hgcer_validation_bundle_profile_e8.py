"""Focused contract coverage for the E.8 global-only farm-review bundle."""

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
RENDERER_PATH = REPO_ROOT / "testing" / "render_pion_hgcer_method_a_acceptance_refinement_figure_library.py"
PROFILE_PATH = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_e8.json"
KINEMATIC = "Q4p4W2p74"
RENDERER_SOURCE = "921511e1e8ea5350218d7b493a5ea69fa422a3f3"
SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)
FROZEN_JSON_BASENAME = (
    "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json"
)
ALLOWED_COMMITTED_FILES = [
    "testing/collect_pion_hgcer_validation_bundle.py",
    "testing/test_collect_pion_hgcer_validation_bundle.py",
    "testing/pion_hgcer_validation_bundle_profile_e8.json",
    "testing/test_pion_hgcer_validation_bundle_profile_e8.py",
]
GLOBAL_ARTIFACTS = [
    {
        "key": "f6_2_acceptance_refinement_validation_json",
        "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json",
        "kind": "json",
        "required": True,
    },
    {
        "key": "e8_figure_library_pdf",
        "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-figure-library.pdf",
        "kind": "file",
        "required": True,
    },
    {
        "key": "e8_figure_library_manifest_json",
        "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-figure-library-manifest.json",
        "kind": "json",
        "required": True,
    },
]


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


collector = _load_module("_e8_bundle_collector", COLLECTOR_PATH)
renderer = _load_module("_e8_figure_library_renderer", RENDERER_PATH)


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


class E8ValidationBundleProfileTests(unittest.TestCase):
    def _declaration(self, profile, key):
        return next(
            declaration for declaration in profile["artifacts"]["global"]
            if declaration["key"] == key
        )

    def _artifact_path(self, root, declaration):
        return Path(root) / collector._format_profile_basename(
            declaration["basename_template"], "global", KINEMATIC, "global",
        )

    def _write_inputs(self, root):
        root = Path(root)
        profile = collector.load_validation_profile(PROFILE_PATH)
        paths = {}
        for declaration in profile["artifacts"]["global"]:
            path = self._artifact_path(root, declaration)
            if declaration["key"] == "e8_figure_library_pdf":
                path.write_bytes(b"e8-figure-library-pdf-fixture")
            elif declaration["key"] == "e8_figure_library_manifest_json":
                path.write_text(json.dumps({
                    "schema_version": renderer.MANIFEST_SCHEMA,
                    "fixture_provenance": "synthetic-e8-profile-test",
                }, sort_keys=True), encoding="utf-8")
            else:
                path.write_text(json.dumps({
                    "schema_version": renderer.ARTIFACT_SCHEMA,
                    "fixture_provenance": "synthetic-f6-2-input",
                }, sort_keys=True), encoding="utf-8")
            paths[declaration["key"]] = path
        return profile, paths

    def _collect(self, temporary, *, command_runner=_clean_command_runner, mutate=None):
        source = Path(temporary) / "source"
        source.mkdir()
        profile, paths = self._write_inputs(source)
        if mutate is not None:
            mutate(paths)
        output = Path(temporary) / "e8-validation.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic=KINEMATIC, output=output,
            profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
            command_runner=command_runner,
        )
        return result, output, profile, paths

    def test_profile_declares_exact_e8_global_only_package(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_e8_f6_2_figure_library_farm_review/v1")
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(tuple((item["phi"], item["epsilon"]) for item in profile["settings"]), SETTINGS)
        self.assertEqual(profile["artifacts"]["global"], GLOBAL_ARTIFACTS)
        self.assertEqual(profile["artifacts"]["settings"], [])
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], RENDERER_SOURCE)
        self.assertEqual(profile["source_identity"]["allowed_committed_files"], ALLOWED_COMMITTED_FILES)
        self.assertEqual(profile["source_identity"]["allowed_non_analysis_path_prefixes"], ["docs/memory/"])

    def test_profile_basenames_match_e8_renderer_constants(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(renderer.KINEMATIC, KINEMATIC)
        self.assertEqual(
            self._artifact_path(
                Path("."), self._declaration(profile, "f6_2_acceptance_refinement_validation_json"),
            ).name,
            FROZEN_JSON_BASENAME,
        )
        self.assertEqual(
            self._artifact_path(Path("."), self._declaration(profile, "e8_figure_library_pdf")).name,
            renderer.PDF_BASENAME,
        )
        self.assertEqual(
            self._artifact_path(Path("."), self._declaration(profile, "e8_figure_library_manifest_json")).name,
            renderer.MANIFEST_BASENAME,
        )

    def test_complete_synthetic_package_is_global_only_and_reports_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, profile, paths = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            manifest = result["manifest"]
            self.assertTrue(output.exists())
            self.assertTrue(manifest["complete"])
            self.assertEqual(manifest["required_analysis_commit"], RENDERER_SOURCE)
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(manifest["unexpected_committed_files_after_required_analysis_commit"], [])
            self.assertEqual(set(manifest["global_artifacts"]), {
                declaration["key"] for declaration in profile["artifacts"]["global"]
            })
            self.assertEqual(len(manifest["settings"]), len(SETTINGS))
            self.assertTrue(all(setting["artifacts"] == {} for setting in manifest["settings"]))
            for declaration in profile["artifacts"]["global"]:
                key = declaration["key"]
                record = manifest["global_artifacts"][key]
                self.assertEqual(record["status"], "exists")
                self.assertEqual(record["sha256"], collector.sha256_file(paths[key]))
                if declaration["kind"] == "json":
                    self.assertEqual(record["json_status"], "valid")
            with zipfile.ZipFile(output) as archive:
                names = archive.namelist()
                self.assertIn("manifest.json", names)
                self.assertIn("source_state.txt", names)
                self.assertIn("source_checks.txt", names)
                setting_payloads = [
                    name for name in names
                    if name.split("/", 1)[0] in {
                        "Left_lowe", "Left_highe", "Center_lowe", "Center_highe", "Right_highe",
                    } and not name.endswith("/")
                ]
                self.assertEqual(setting_payloads, [])
                for key, path in paths.items():
                    archive_path = manifest["global_artifacts"][key]["archive_path"]
                    self.assertEqual(sum(name == archive_path for name in names), 1)
                    self.assertEqual(archive.read(archive_path), path.read_bytes())

    def test_required_artifact_failures_remain_best_effort_and_incomplete(self):
        cases = (
            ("missing_frozen_json", lambda paths: paths["f6_2_acceptance_refinement_validation_json"].unlink(), "f6_2_acceptance_refinement_validation_json", "missing"),
            ("invalid_frozen_json", lambda paths: paths["f6_2_acceptance_refinement_validation_json"].write_bytes(b"not-json"), "f6_2_acceptance_refinement_validation_json", "invalid"),
            ("missing_e8_pdf", lambda paths: paths["e8_figure_library_pdf"].unlink(), "e8_figure_library_pdf", "missing"),
            ("missing_e8_manifest", lambda paths: paths["e8_figure_library_manifest_json"].unlink(), "e8_figure_library_manifest_json", "missing"),
            ("invalid_e8_manifest", lambda paths: paths["e8_figure_library_manifest_json"].write_bytes(b"not-json"), "e8_figure_library_manifest_json", "invalid"),
        )
        for label, mutation, key, status in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, mutate=mutation)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                record = result["manifest"]["global_artifacts"][key]
                self.assertEqual(record["status"], status)
                if status == "invalid":
                    self.assertEqual(record["json_status"], "invalid")
                self.assertIn(
                    "missing_source_artifact" if status == "missing" else "source_artifact_json_invalid",
                    {entry["code"] for entry in result["manifest"]["errors"]},
                )

    def test_source_identity_is_fail_closed_and_allows_only_declared_paths(self):
        def missing_ancestor(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "merge-base", "--is-ancestor"]:
                result["returncode"] = 1
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _profile, _paths = self._collect(
                temporary, command_runner=missing_ancestor,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertFalse(result["manifest"]["complete"])
            self.assertIn(
                "required_analysis_commit_not_present",
                {entry["code"] for entry in result["manifest"]["errors"]},
            )

        def allowed_changes(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result["stdout"] = "\n".join(ALLOWED_COMMITTED_FILES) + "\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _profile, _paths = self._collect(
                temporary, command_runner=allowed_changes,
            )
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.exists())
            self.assertTrue(result["manifest"]["complete"])
            self.assertEqual(result["manifest"]["unexpected_committed_files_after_required_analysis_commit"], [])

        for changed_path in ("src/cuts/rand_sub.py", "testing/unrelated_e8_change.py"):
            def unexpected_change(command, cwd, path=changed_path):
                result = _clean_command_runner(command, cwd)
                if list(command)[:3] == ["git", "diff", "--name-only"]:
                    result["stdout"] = path + "\n"
                return result

            with self.subTest(path=changed_path), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(
                    temporary, command_runner=unexpected_change,
                )
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                self.assertEqual(
                    result["manifest"]["unexpected_committed_files_after_required_analysis_commit"],
                    [changed_path],
                )
                self.assertIn(
                    "unexpected_committed_files_after_required_analysis_commit",
                    {entry["code"] for entry in result["manifest"]["errors"]},
                )

    def test_output_collision_fails_without_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            self._write_inputs(source)
            output = Path(temporary) / "e8-validation.zip"
            output.write_bytes(b"existing-output")
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(
                    outdir=source, kinematic=KINEMATIC, output=output,
                    profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
                    command_runner=_clean_command_runner,
                )
            self.assertEqual(output.read_bytes(), b"existing-output")


if __name__ == "__main__":
    unittest.main()
