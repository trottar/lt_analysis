"""Focused collector coverage for the E.8.1 five-setting procedure-PDF gate."""

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
PROFILE_PATH = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_e8_1.json"
KINEMATIC = "Q4p4W2p74"
REVIEWED_SOURCE = "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"
SETTINGS = (
    ("Left", "lowe"),
    ("Left", "highe"),
    ("Center", "lowe"),
    ("Center", "highe"),
    ("Right", "highe"),
)
FROZEN_JSON_BASENAME = (
    "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json"
)
ALLOWED_COMMITTED_FILES = [
    "testing/pion_hgcer_validation_bundle_profile_e8_1.json",
    "testing/test_pion_hgcer_validation_bundle_profile_e8_1.py",
]
GLOBAL_ARTIFACTS = [
    {
        "key": "f6_2_acceptance_refinement_validation_json",
        "basename_template": "{kinematic}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json",
        "kind": "json",
        "required": True,
    },
]
SETTING_ARTIFACTS = [
    {
        "key": "full_background_subtraction_pdf",
        "basename_template": "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction.pdf",
        "kind": "file",
        "required": True,
    },
    {
        "key": "full_background_subtraction_manifest_json",
        "basename_template": "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction-manifest.json",
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


collector = _load_module("_e8_1_bundle_collector", COLLECTOR_PATH)


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


class E81ValidationBundleProfileTests(unittest.TestCase):
    @staticmethod
    def _declaration(profile, scope, key):
        return next(
            declaration for declaration in profile["artifacts"][scope]
            if declaration["key"] == key
        )

    @staticmethod
    def _artifact_path(root, declaration, *, phi="global", epsilon="global"):
        return Path(root) / collector._format_profile_basename(
            declaration["basename_template"], phi, KINEMATIC, epsilon,
        )

    def _write_inputs(self, root):
        root = Path(root)
        profile = collector.load_validation_profile(PROFILE_PATH)
        paths = {"global": {}, "settings": {}}
        for declaration in profile["artifacts"]["global"]:
            path = self._artifact_path(root, declaration)
            path.write_bytes(json.dumps({
                "schema_version": "synthetic_f6_2_validation_artifact/v1",
                "fixture_provenance": "e8_1_profile_test",
            }, sort_keys=True, separators=(",", ":")).encode("utf-8"))
            paths["global"][declaration["key"]] = path
        for phi, epsilon in SETTINGS:
            setting_paths = {}
            for declaration in profile["artifacts"]["settings"]:
                path = self._artifact_path(
                    root, declaration, phi=phi, epsilon=epsilon,
                )
                if declaration["kind"] == "json":
                    path.write_bytes(json.dumps({
                        "schema_version": "full_background_subtraction_page_manifest/v1",
                        "setting": {"phi": phi, "epsilon": epsilon},
                    }, sort_keys=True, separators=(",", ":")).encode("utf-8"))
                else:
                    path.write_bytes(
                        "synthetic-procedure-pdf:{}:{}".format(phi, epsilon).encode("ascii")
                    )
                setting_paths[declaration["key"]] = path
            paths["settings"][(phi, epsilon)] = setting_paths
        return profile, paths

    def _collect(self, temporary, *, command_runner=_clean_command_runner, mutate=None):
        source = Path(temporary) / "source"
        source.mkdir()
        profile, paths = self._write_inputs(source)
        if mutate is not None:
            mutate(paths)
        output = Path(temporary) / "e8-1-validation.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic=KINEMATIC, output=output,
            profile_path=PROFILE_PATH, repo_root=REPO_ROOT,
            command_runner=command_runner,
        )
        return result, output, profile, paths

    def test_profile_declares_exact_five_setting_e8_1_package(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(
            profile["validation_profile"],
            "phase_e8_1_full_background_procedure_pdf_farm_review/v1",
        )
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(
            tuple((item["phi"], item["epsilon"]) for item in profile["settings"]),
            SETTINGS,
        )
        self.assertEqual(profile["artifacts"]["global"], GLOBAL_ARTIFACTS)
        self.assertEqual(profile["artifacts"]["settings"], SETTING_ARTIFACTS)
        self.assertEqual(
            profile["source_identity"]["required_analysis_commit"], REVIEWED_SOURCE,
        )
        self.assertEqual(
            profile["source_identity"]["allowed_committed_files"],
            ALLOWED_COMMITTED_FILES,
        )
        self.assertEqual(
            profile["source_identity"]["allowed_non_analysis_path_prefixes"],
            ["docs/memory/"],
        )
        self.assertEqual(
            self._artifact_path(
                Path("."),
                self._declaration(
                    profile, "global", "f6_2_acceptance_refinement_validation_json",
                ),
            ).name,
            FROZEN_JSON_BASENAME,
        )
        self.assertEqual(
            self._artifact_path(
                Path("."),
                self._declaration(profile, "settings", "full_background_subtraction_pdf"),
                phi="Left", epsilon="highe",
            ).name,
            "Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction.pdf",
        )
        self.assertEqual(
            self._artifact_path(
                Path("."),
                self._declaration(
                    profile, "settings", "full_background_subtraction_manifest_json",
                ),
                phi="Left", epsilon="highe",
            ).name,
            "Left_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction-manifest.json",
        )

    def test_complete_synthetic_package_archives_all_five_settings(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, profile, paths = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.exists())
            manifest = result["manifest"]
            self.assertTrue(manifest["complete"])
            self.assertEqual(manifest["required_analysis_commit"], REVIEWED_SOURCE)
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(
                manifest["unexpected_committed_files_after_required_analysis_commit"], [],
            )
            self.assertEqual(
                manifest["requested_settings"],
                [{"phi": phi, "epsilon": epsilon} for phi, epsilon in SETTINGS],
            )
            frozen = manifest["global_artifacts"]["f6_2_acceptance_refinement_validation_json"]
            self.assertEqual(frozen["status"], "exists")
            self.assertEqual(frozen["json_status"], "valid")
            self.assertEqual(
                frozen["sha256"],
                collector.sha256_file(paths["global"]["f6_2_acceptance_refinement_validation_json"]),
            )
            self.assertEqual(len(manifest["settings"]), len(SETTINGS))
            with zipfile.ZipFile(output) as archive:
                names = archive.namelist()
                self.assertIn("global/" + FROZEN_JSON_BASENAME, names)
                self.assertIn("manifest.json", names)
                self.assertIn("source_state.txt", names)
                self.assertIn("source_checks.txt", names)
                self.assertEqual(
                    archive.read(frozen["archive_path"]),
                    paths["global"]["f6_2_acceptance_refinement_validation_json"].read_bytes(),
                )
                for setting_manifest, (phi, epsilon) in zip(manifest["settings"], SETTINGS):
                    self.assertEqual((setting_manifest["phi"], setting_manifest["epsilon"]), (phi, epsilon))
                    self.assertEqual(
                        set(setting_manifest["artifacts"]),
                        {declaration["key"] for declaration in profile["artifacts"]["settings"]},
                    )
                    for declaration in profile["artifacts"]["settings"]:
                        key = declaration["key"]
                        record = setting_manifest["artifacts"][key]
                        source = paths["settings"][(phi, epsilon)][key]
                        self.assertEqual(record["status"], "exists")
                        self.assertEqual(record["sha256"], collector.sha256_file(source))
                        if declaration["kind"] == "json":
                            self.assertEqual(record["json_status"], "valid")
                        self.assertIn(record["archive_path"], names)
                        self.assertEqual(archive.read(record["archive_path"]), source.read_bytes())

    def test_required_artifact_failures_are_best_effort_and_incomplete(self):
        cases = (
            (
                "missing_frozen_json",
                lambda paths: paths["global"]["f6_2_acceptance_refinement_validation_json"].unlink(),
                "global", None, "f6_2_acceptance_refinement_validation_json", "missing",
            ),
            (
                "invalid_frozen_json",
                lambda paths: paths["global"]["f6_2_acceptance_refinement_validation_json"].write_bytes(b"not-json"),
                "global", None, "f6_2_acceptance_refinement_validation_json", "invalid",
            ),
            (
                "missing_setting_pdf",
                lambda paths: paths["settings"][("Left", "highe")]["full_background_subtraction_pdf"].unlink(),
                "setting", ("Left", "highe"), "full_background_subtraction_pdf", "missing",
            ),
            (
                "missing_setting_manifest",
                lambda paths: paths["settings"][("Left", "highe")]["full_background_subtraction_manifest_json"].unlink(),
                "setting", ("Left", "highe"), "full_background_subtraction_manifest_json", "missing",
            ),
            (
                "invalid_setting_manifest",
                lambda paths: paths["settings"][("Left", "highe")]["full_background_subtraction_manifest_json"].write_bytes(b"not-json"),
                "setting", ("Left", "highe"), "full_background_subtraction_manifest_json", "invalid",
            ),
        )
        for label, mutation, scope, setting, key, status in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(temporary, mutate=mutation)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                if scope == "global":
                    record = result["manifest"]["global_artifacts"][key]
                else:
                    setting_manifest = next(
                        value for value in result["manifest"]["settings"]
                        if (value["phi"], value["epsilon"]) == setting
                    )
                    record = setting_manifest["artifacts"][key]
                self.assertEqual(record["status"], status)
                if status == "invalid":
                    self.assertEqual(record["json_status"], "invalid")
                self.assertIn(
                    "missing_source_artifact" if status == "missing" else "source_artifact_json_invalid",
                    {entry["code"] for entry in result["manifest"]["errors"]},
                )

    def test_source_identity_is_fail_closed_with_exact_exceptions(self):
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

        for changed_path in ALLOWED_COMMITTED_FILES + ["docs/memory/CURRENT.md"]:
            def allowed_change(command, cwd, path=changed_path):
                result = _clean_command_runner(command, cwd)
                if list(command)[:3] == ["git", "diff", "--name-only"]:
                    result["stdout"] = path + "\n"
                return result

            with self.subTest(allowed_path=changed_path), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(
                    temporary, command_runner=allowed_change,
                )
                self.assertEqual(result["returncode"], 0)
                self.assertTrue(output.exists())
                self.assertTrue(result["manifest"]["complete"])
                self.assertEqual(
                    result["manifest"]["unexpected_committed_files_after_required_analysis_commit"], [],
                )

        for changed_path in (
            "src/cuts/full_background_subtraction_plots.py",
            "src/cuts/rand_sub.py",
            "src/utility/unrelated_e8_change.py",
            "testing/unrelated_e8_change.py",
        ):
            def rejected_change(command, cwd, path=changed_path):
                result = _clean_command_runner(command, cwd)
                if list(command)[:3] == ["git", "diff", "--name-only"]:
                    result["stdout"] = path + "\n"
                return result

            with self.subTest(rejected_path=changed_path), tempfile.TemporaryDirectory() as temporary:
                result, output, _profile, _paths = self._collect(
                    temporary, command_runner=rejected_change,
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
            output = Path(temporary) / "e8-1-validation.zip"
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
