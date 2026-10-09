"""Distinct generic profile and unchanged-collector delivery fixtures."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from testing import run_e8_2_canonical_five_baseline_yield_gate as owner
from testing.test_run_e8_2_canonical_five_baseline_yield_gate import fixture, write_json, SHA, ROOT


def source_runner(command, cwd):
    command = [str(x) for x in command]
    return {"command": command, "returncode": 0, "stderr": "",
            "stdout": SHA + "\n" if command == ["git", "rev-parse", "HEAD"] else ""}


class ProfileTests(unittest.TestCase):
    def test_exact_inventory_and_no_shared_artifact_aliases(self):
        output = ROOT / "synthetic-attempt.zip"
        profile = owner.resolved_profile(ROOT, SHA, output)
        self.assertEqual(len(profile["settings"]), 5)
        self.assertNotIn({"phi": "Right", "epsilon": "lowe"}, profile["settings"])
        declared = owner.names(profile)
        self.assertEqual(len(declared), 29)
        self.assertEqual(len({Path(name).name for name in declared}), 29)
        self.assertEqual(len(profile["artifacts"]["global"]), 9)
        self.assertTrue(all(e["required"] for e in declared.values()))
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], SHA)

    def test_profile_tampering_rejected(self):
        output = ROOT / "synthetic-attempt.zip"
        profile = json.loads((ROOT / owner.PROFILE).read_text())
        from unittest.mock import patch
        for key in ("settings", "artifacts", "source_identity", "validation_profile"):
            wrong = deepcopy(profile)
            if key == "settings":
                wrong[key].append({"phi": "Right", "epsilon": "lowe"})
            elif key == "artifacts":
                wrong[key]["settings"].pop()
            elif key == "source_identity":
                wrong[key]["allowed_committed_files"] = ["src/main.py"]
            else:
                wrong[key] = "unrelated/v1"
            with self.subTest(key=key), patch.object(owner.collector, "load_validation_profile", return_value=wrong):
                with self.assertRaises(Exception):
                    owner.resolved_profile(ROOT, SHA, output)

    def test_unchanged_collector_delivers_both_shared_epsilons_and_ten_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            worktree, outdir = root / "analysis", root / "out"
            worktree.mkdir()
            outdir.mkdir()
            fixture(worktree, outdir)
            output = root / "collected.zip"
            profile = owner.resolved_profile(ROOT, SHA, output)
            records, cells, metadata = owner.audit_outputs(worktree, outdir, profile, 0, {})
            write_json(outdir / owner.artifact_name(profile, "run_summary"), {"cells": cells, "epsilon_metadata": metadata, "boundaries": owner.BOUNDARIES})
            records[owner.artifact_name(profile, "run_summary")] = {"sha256": owner.sha256(outdir / owner.artifact_name(profile, "run_summary")), "bytes": (outdir / owner.artifact_name(profile, "run_summary")).stat().st_size}
            profile_path = root / "profile.json"
            write_json(profile_path, profile)
            output = root / "collected.zip"
            result = owner.collector.collect_validation_bundle(outdir=outdir, kinematic=owner.KINEMATIC,
                output=output, profile_path=profile_path, repo_root=ROOT, command_runner=source_runner)
            self.assertEqual(result["returncode"], 0, result["manifest"]["errors"])
            owner.verify_zip(output, SHA, records, profile)
            original_zip = output.read_bytes()
            failed_output = root / "missing-attempt.zip"
            failed_profile = owner.resolved_profile(ROOT, SHA, failed_output)
            failed_profile_path = root / "missing-profile.json"
            write_json(failed_profile_path, failed_profile)
            failed = owner.collector.collect_validation_bundle(outdir=outdir, kinematic=owner.KINEMATIC,
                output=failed_output, profile_path=failed_profile_path, repo_root=ROOT, command_runner=source_runner)
            self.assertEqual(failed["returncode"], 1)
            self.assertFalse(failed["manifest"]["complete"])
            with self.assertRaisesRegex(Exception, "bundle_incomplete"):
                owner.verify_zip(failed_output, SHA, records, failed_profile)
            self.assertEqual(output.read_bytes(), original_zip)


if __name__ == "__main__":
    unittest.main()
