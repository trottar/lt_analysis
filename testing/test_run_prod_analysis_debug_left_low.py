"""Static contract coverage for the narrow Left/lowe launcher debug mode."""

from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER_PATH = REPO_ROOT / "run_Prod_Analysis.sh"
MAIN_PATH = REPO_ROOT / "src" / "main.py"


class RunProdAnalysisDebugLeftLowTests(unittest.TestCase):
    """Verify source-level orchestration without starting the production launcher."""

    @classmethod
    def setUpClass(cls):
        cls.launcher = LAUNCHER_PATH.read_text(encoding="utf-8")
        cls.main = MAIN_PATH.read_text(encoding="utf-8")

    def _full_low_high_loop(self):
        preflight_start = self.launcher.index('preflight_dir=$(mktemp -d')
        loop_start = self.launcher.index('for j in "low" "high"; do', preflight_start)
        loop_end = self.launcher.index('elif [[ $i_flag = "true"', loop_start)
        return self.launcher[loop_start:loop_end]

    @staticmethod
    def _bash():
        if os.name == "nt":
            git_bash = Path(r"C:\Program Files\Git\bin\bash.exe")
            return str(git_bash) if git_bash.is_file() else None
        return shutil.which("bash")

    def test_shell_syntax_is_valid(self):
        bash = self._bash()
        if bash is None:
            self.skipTest("a local Bash executable is required for shell syntax coverage")
        result = subprocess.run(
            [bash, "-n", str(LAUNCHER_PATH)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_debug_option_is_documented_and_uses_its_own_state(self):
        self.assertRegex(self.launcher, r"while getopts '[^']*d[^']*' flag; do")
        self.assertIn("d) d_flag='true' ;;", self.launcher)
        self.assertIn("fast kaon Left / lowe debug mode", self.launcher)
        self.assertIn("$d_flag = \"true\" && $ParticleType != \"kaon\"", self.launcher)
        self.assertNotRegex(self.launcher, r"d\)\s+DEBUG=")

    def test_debug_flag_uses_existing_flagged_q2_w_convention(self):
        self.assertRegex(
            self.launcher,
            r'if \[\[ (?=[^\n]*\$d_flag = "true")[^\n]* \]\]; then\s+'
            r'EPSILON=\$j\s+Q2=\$2\s+W=\$3',
        )

    def test_preflight_stays_paired_and_unrestricted_before_debug_full_low(self):
        preflight_start = self.launcher.index('preflight_dir=$(mktemp -d')
        full_loop_start = self.launcher.index('for j in "low" "high"; do', preflight_start)
        preflight = self.launcher[preflight_start:full_loop_start]
        low_capture = preflight.index('LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE="${low_capture}"')
        high_capture = preflight.index('LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE="${high_capture}"')
        shared_preflight = preflight.index("shared_canonical_binning_preflight.py")
        self.assertLess(low_capture, high_capture)
        self.assertLess(high_capture, shared_preflight)
        self.assertNotIn("LT_ANALYSIS_DEBUG_LEFT_LOW", preflight)
        self.assertLess(shared_preflight, self.launcher.index("LT_ANALYSIS_DEBUG_LEFT_LOW=1", full_loop_start))

    def test_debug_control_applies_only_to_full_low_and_stops_before_full_high(self):
        full_loop = self._full_low_high_loop()
        debug_low = 'LT_ANALYSIS_DEBUG_LEFT_LOW=1 python3 main.py "${main_args_low[@]}" || exit 1'
        full_high = 'python3 main.py "${main_args_high[@]}" || exit 1'
        self.assertIn(debug_low, full_loop)
        self.assertNotIn("LT_ANALYSIS_DEBUG_LEFT_LOW=1 python3 main.py \"${main_args_high[@]}\"", full_loop)
        self.assertIn(full_high, full_loop)
        self.assertRegex(
            full_loop,
            re.escape(debug_low)
            + r"\s+echo\s+echo \"Left / lowe debug analysis completed\.\""
            + r"\s+echo \"Full high-epsilon processing is intentionally skipped in -d debug mode\.\""
            + r"\s+exit 0",
        )
        self.assertLess(full_loop.index(debug_low), full_loop.index(full_high))

    def test_debug_low_failure_is_not_converted_to_success(self):
        full_loop = self._full_low_high_loop()
        self.assertIn(
            'LT_ANALYSIS_DEBUG_LEFT_LOW=1 python3 main.py "${main_args_low[@]}" || exit 1',
            full_loop,
        )
        self.assertNotIn("|| true", full_loop)

    def test_ordinary_full_low_and_high_paths_remain_present(self):
        full_loop = self._full_low_high_loop()
        self.assertRegex(
            full_loop,
            r"fi\s+python3 main\.py \"\$\{main_args_low\[@\]\}\" \|\| exit 1",
        )
        self.assertIn('python3 main.py "${main_args_high[@]}" || exit 1', full_loop)

    def test_main_debug_selector_is_left_only_and_leaves_binning_config_untouched(self):
        ast.parse(self.main, filename=str(MAIN_PATH))
        self.assertIn('phisetlist = ["Center", "Left", "Right"]', self.main)
        selector_start = self.main.index('debug_left_low = os.environ.get("LT_ANALYSIS_DEBUG_LEFT_LOW"')
        selector_end = self.main.index("ROOT.gROOT.SetBatch", selector_start)
        selector = self.main[selector_start:selector_end]
        self.assertIn('if debug_left_low:', selector)
        self.assertIn('ParticleType != "kaon"', selector)
        self.assertIn('str(EPSSET).strip().lower() != "low"', selector)
        self.assertIn('LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE', selector)
        self.assertIn('phisetlist = ["Left"]', selector)
        self.assertNotIn("inpDict", selector)
        self.assertNotIn("canonical_t_binning", selector)


if __name__ == "__main__":
    unittest.main()
