"""Focused synthetic checks for the detached F.4.Refresh.2 bundle profile."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from testing import collect_pion_hgcer_validation_bundle as collector


REPO_ROOT = Path(__file__).resolve().parents[1]
PROFILE = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile_f4_refresh2.json"
KINEMATIC = "Q4p4W2p74"
SOURCE = "141a3d04f9e5d07be21dba14e0e63212c3990bf1"
SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
            ("Center", "highe"), ("Right", "highe"))
DECLARATIONS = (
    ("f4_refresh1_comparison_input", "{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-comparison-input.json"),
    ("f4_refresh2_candidate_f2", "{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json"),
    ("f4_refresh2_candidate_f3", "{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json"),
    ("f4_refresh2_candidate_f4", "{kinematic}_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json"),
    ("f4_refresh2_materialization_manifest", "{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-materialization-manifest.json"),
)
ALLOWED = ["testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json",
           "testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py",
           "testing/run_f4_refresh2_materialize_verify_package.py",
           "testing/test_run_f4_refresh2_materialize_verify_package.py",
           "testing/test_memory_health.py",
           "tools/check_memory_health.py"]


def _runner(command, _cwd):
    command = [str(item) for item in command]
    stdout = ""
    if command == ["git", "rev-parse", "HEAD"]:
        stdout = "profile-test-head\n"
    elif command == ["git", "rev-parse", "HEAD^"]:
        stdout = "profile-test-parent\n"
    elif command[-1:] == ["--version"]:
        stdout = "Python synthetic\n"
    return {"command": command, "returncode": 0, "stdout": stdout, "stderr": ""}


class F4Refresh2BundleProfileTests(unittest.TestCase):
    def _inputs(self, root):
        source = Path(root) / "source"
        source.mkdir()
        paths = {}
        for key, template in DECLARATIONS:
            path = source / template.format(kinematic=KINEMATIC)
            path.write_text(json.dumps({"synthetic_artifact": key}), encoding="utf-8")
            paths[key] = path
        return source, paths

    def _collect(self, root, *, mutate=None, runner=_runner):
        source, paths = self._inputs(root)
        if mutate is not None:
            mutate(paths)
        output = Path(root) / "f4-refresh2-synthetic.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic=KINEMATIC, output=output,
            profile_path=PROFILE, repo_root=REPO_ROOT, command_runner=runner,
        )
        return result, output, paths

    def test_declarations_and_source_identity_are_exact(self):
        profile = collector.load_validation_profile(PROFILE)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_f4_refresh2_current_baseline_candidate_materialization_farm_review/v1")
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual(tuple((row["phi"], row["epsilon"]) for row in profile["settings"]), SETTINGS)
        self.assertEqual(profile["artifacts"]["settings"], [])
        self.assertEqual(profile["artifacts"]["global"], [
            {"key": key, "basename_template": template, "kind": "json", "required": True}
            for key, template in DECLARATIONS
        ])
        self.assertEqual(profile["source_identity"], {
            "required_analysis_commit": SOURCE,
            "allowed_committed_files": ALLOWED,
            "allowed_non_analysis_path_prefixes": ["docs/memory/"],
        })

    def test_complete_package_has_only_five_global_payloads(self):
        with tempfile.TemporaryDirectory() as temporary:
            def add_undeclared(paths):
                next(iter(paths.values())).parent.joinpath("undeclared.json").write_text("{}", encoding="utf-8")
            result, output, paths = self._collect(temporary, mutate=add_undeclared)
            manifest = result["manifest"]
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(manifest["complete"])
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(manifest["unexpected_committed_files_after_required_analysis_commit"], [])
            self.assertEqual(set(manifest["global_artifacts"]), set(paths))
            self.assertEqual([(row["phi"], row["epsilon"]) for row in manifest["settings"]], list(SETTINGS))
            self.assertTrue(all(row["artifacts"] == {} for row in manifest["settings"]))
            with zipfile.ZipFile(output) as archive:
                names = archive.namelist()
                expected_files = {"manifest.json", "source_state.txt", "source_checks.txt"}
                expected_files.update("global/" + path.name for path in paths.values())
                self.assertEqual({name for name in names if not name.endswith("/")}, expected_files)
                for key, path in paths.items():
                    record = manifest["global_artifacts"][key]
                    self.assertEqual(record["status"], "exists")
                    self.assertEqual(record["json_status"], "valid")
                    self.assertEqual(record["archive_path"], "global/" + path.name)
                    self.assertEqual(names.count(record["archive_path"]), 1)
                    self.assertEqual(archive.read(record["archive_path"]), path.read_bytes())

    def test_each_required_global_artifact_fails_when_missing_or_invalid(self):
        for key, _template in DECLARATIONS:
            for label, mutation, expected_status in (
                ("missing", lambda paths, item=key: paths[item].unlink(), "missing"),
                ("malformed_json", lambda paths, item=key: paths[item].write_bytes(b"not-json"), "invalid"),
            ):
                with self.subTest(key=key, case=label), tempfile.TemporaryDirectory() as temporary:
                    result, output, _paths = self._collect(temporary, mutate=mutation)
                    self.assertEqual(result["returncode"], 1)
                    self.assertTrue(output.exists())
                    self.assertFalse(result["manifest"]["complete"])
                    record = result["manifest"]["global_artifacts"][key]
                    self.assertEqual(record["status"], expected_status)
                    if expected_status == "invalid":
                        self.assertEqual(record["json_status"], "invalid")

    def test_source_provenance_and_output_collision_fail_closed(self):
        def runner_with(*, ancestor=True, changed=""):
            def run(command, cwd):
                result = _runner(command, cwd)
                if list(command)[:3] == ["git", "merge-base", "--is-ancestor"] and not ancestor:
                    result["returncode"] = 1
                if list(command)[:3] == ["git", "diff", "--name-only"]:
                    result["stdout"] = changed
                return result
            return run

        for label, runner, error in (
            ("missing_ancestor", runner_with(ancestor=False), "required_analysis_commit_not_present"),
            ("unreviewed_testing_helper", runner_with(changed="testing/unreviewed_helper.py\n"), "unexpected_committed_files_after_required_analysis_commit"),
            ("unreviewed_tools_helper", runner_with(changed="tools/unreviewed_helper.py\n"), "unexpected_committed_files_after_required_analysis_commit"),
            ("unreviewed_materializer", runner_with(changed="testing/materialize_method_a_current_baseline_authority.py\n"), "unexpected_committed_files_after_required_analysis_commit"),
            ("unreviewed_science", runner_with(changed="src/cuts/unreviewed.py\n"), "unexpected_committed_files_after_required_analysis_commit"),
        ):
            with self.subTest(case=label), tempfile.TemporaryDirectory() as temporary:
                result, output, _paths = self._collect(temporary, runner=runner)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertFalse(result["manifest"]["complete"])
                self.assertIn(error, {item["code"] for item in result["manifest"]["errors"]})
        with tempfile.TemporaryDirectory() as temporary:
            changed = "\n".join((*ALLOWED, "docs/memory/CURRENT.md")) + "\n"
            result, _output, _paths = self._collect(temporary, runner=runner_with(changed=changed))
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(result["manifest"]["complete"])
            self.assertEqual(result["manifest"]["unexpected_committed_files_after_required_analysis_commit"], [])
        with tempfile.TemporaryDirectory() as temporary:
            source, _paths = self._inputs(temporary)
            output = Path(temporary) / "f4-refresh2-synthetic.zip"
            output.write_bytes(b"already exists")
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(
                    outdir=source, kinematic=KINEMATIC, output=output,
                    profile_path=PROFILE, repo_root=REPO_ROOT, command_runner=_runner,
                )
            self.assertEqual(output.read_bytes(), b"already exists")


if __name__ == "__main__":
    unittest.main()
