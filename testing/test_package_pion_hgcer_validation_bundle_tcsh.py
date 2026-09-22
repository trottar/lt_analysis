"""Static and syntax checks for the farm-only validation-bundle wrapper."""

from __future__ import annotations

import shutil
import subprocess
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("package_pion_hgcer_validation_bundle.tcsh")


class PackageValidationBundleTcshTests(unittest.TestCase):
    def setUp(self) -> None:
        self.source = SCRIPT_PATH.read_text(encoding="utf-8")

    def test_contract_is_input_driven_and_uses_the_generic_collector(self) -> None:
        self.assertTrue(self.source.startswith("#!/bin/tcsh -f\n"))
        for option in (
            "--python", "--bundle-commit", "--profile", "--artifact-dir",
            "--kinematic", "--output", "--phi", "--epsilon", "--immutable",
        ):
            self.assertIn("case {}:".format(option), self.source)
        self.assertIn(
            "testing/collect_pion_hgcer_validation_bundle.py", self.source,
        )
        self.assertIn("git worktree add --detach", self.source)
        self.assertIn("git worktree remove --force \"$worktree\"", self.source)
        self.assertIn("sha256sum \"$immutable_path\"", self.source)
        self.assertIn("unzip -t \"$zip\"", self.source)
        self.assertIn(
            "artifact directory must be the canonical KaonLT artifact root",
            self.source,
        )
        self.assertIn(
            "ZIP output must be in the canonical Globus transfer directory",
            self.source,
        )

    def test_contract_cannot_run_or_mutate_analysis(self) -> None:
        forbidden = (
            "run_Prod_Analysis.sh",
            "main.py",
            "main_iter.py",
            "main_auto.py",
            "git clean",
            "git reset",
            "git stash",
            "eval ",
            "sudo ",
            "qsub ",
            "swif ",
        )
        for token in forbidden:
            self.assertNotIn(token, self.source)

    def test_tcsh_parse_when_available(self) -> None:
        tcsh = shutil.which("tcsh")
        if tcsh is None:
            self.skipTest("tcsh is unavailable on this local host")
        result = subprocess.run(
            [tcsh, "-n", str(SCRIPT_PATH)],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
