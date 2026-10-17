"""Local contract coverage for all-low launcher debug and legacy Left-only control."""

from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER_PATH = REPO_ROOT / "run_Prod_Analysis.sh"
MAIN_PATH = REPO_ROOT / "src" / "main.py"
DIAMOND_PATH = REPO_ROOT / "src" / "cuts" / "diamond.py"


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
        self.assertIn("fast kaon debug mode for all available low-epsilon settings", self.launcher)
        self.assertIn("$d_flag = \"true\" && $ParticleType != \"kaon\"", self.launcher)
        self.assertNotRegex(self.launcher, r"d\)\s+DEBUG=")

    def _cleanup_block(self):
        start = self.launcher.index('if [[ $i_flag != "true" && $a_flag != "true" ]]; then')
        end = self.launcher.index('\n# Efficiency csv file', start)
        return self.launcher[start:end]

    def _exercise_cleanup_block(self, debug, validation_status=0):
        """Execute only the real cleanup block with harmless command stubs."""
        bash = self._bash()
        if bash is None:
            self.skipTest("a local Bash executable is required for cleanup-path coverage")
        with tempfile.TemporaryDirectory() as directory:
            symlinks = Path(directory) / "set_SymLinks.sh"
            symlinks.write_text('#!/bin/bash\nprintf "symlinks:%s\\n" "$1"\n', encoding="utf-8")
            symlinks.chmod(0o755)
            harness = (
                'i_flag=false; a_flag=false; ParticleType=kaon\n'
                'd_flag="$1"\n'
                'validate_external_sigma0_paths_before_cleanup() { '
                'printf "validate\\n"; return ' + str(validation_status) + '; }\n'
                'git() { printf "git:%s\\n" "$*"; }\n'
            )
            return subprocess.run(
                [bash, "-c", harness + self._cleanup_block(), "cleanup-test", debug],
                cwd=directory, text=True, capture_output=True, check=False,
            )

    def test_debug_bypasses_only_cleanup_and_ordinary_cleanup_is_preserved(self):
        self.assertEqual(len(re.findall(r'^\s*git clean -fdx\s*$', self.launcher, re.MULTILINE)), 1)
        self.assertRegex(
            self._cleanup_block(),
            r'if \[\[ \$d_flag != "true" \]\]; then\s+git clean -fdx\s+fi\s+'
            r'\./set_SymLinks\.sh \$ParticleType\s+fi',
        )
        for debug in ("true", "false", ""):
            with self.subTest(debug=debug):
                result = self._exercise_cleanup_block(debug)
                self.assertEqual(result.returncode, 0, result.stderr)
                expected = ["validate", "symlinks:kaon"] if debug == "true" else ["validate", "git:clean -fdx", "symlinks:kaon"]
                self.assertEqual(result.stdout.splitlines(), expected)

    def test_external_sigma0_failure_blocks_cleanup_and_symlinks_in_both_modes(self):
        for debug in ("true", "false"):
            with self.subTest(debug=debug):
                result = self._exercise_cleanup_block(debug, validation_status=1)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(result.stdout.splitlines(), ["validate"])

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
        self.assertLess(shared_preflight, self.launcher.index("LT_ANALYSIS_DEBUG_LEFT_LOW=", full_loop_start))

    def test_debug_control_applies_only_to_full_low_and_stops_before_full_high(self):
        full_loop = self._full_low_high_loop()
        debug_low = 'LT_ANALYSIS_DEBUG_LEFT_LOW= python3 main.py "${main_args_low[@]}" || exit 1'
        full_high = 'python3 main.py "${main_args_high[@]}" || exit 1'
        self.assertIn(debug_low, full_loop)
        self.assertNotIn("LT_ANALYSIS_DEBUG_LEFT_LOW=1 python3 main.py \"${main_args_high[@]}\"", full_loop)
        self.assertIn(full_high, full_loop)
        self.assertRegex(
            full_loop,
            re.escape(debug_low)
            + r"\s+echo\s+echo \"All available low-epsilon settings completed in debug mode\.\""
            + r"\s+echo \"Full high-epsilon processing is intentionally skipped in -d debug mode\.\""
            + r"\s+exit 0",
        )
        self.assertLess(full_loop.index(debug_low), full_loop.index(full_high))

    def test_debug_low_failure_is_not_converted_to_success(self):
        full_loop = self._full_low_high_loop()
        self.assertIn(
            'LT_ANALYSIS_DEBUG_LEFT_LOW= python3 main.py "${main_args_low[@]}" || exit 1',
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

    def _exercise_full_loop(self, debug, inherited_selector, low_status=0, high_status=0):
        """Run only the final orchestration loop; every external call is a stub."""
        bash = self._bash()
        if bash is None:
            self.skipTest("a local Bash executable is required for loop coverage")
        harness = (
            'd_flag="$1"; export LT_ANALYSIS_DEBUG_LEFT_LOW="$2"\n'
            'main_args_low=("low input"); main_args_high=("high input")\n'
            'export_background_sample_paths() { :; }\n'
            'python3() {\n'
            '  printf "CALL:%s:%s:%s\\n" "$EPSILON" "$LT_ANALYSIS_DEBUG_LEFT_LOW" "$2"\n'
            '  if [[ "$EPSILON" = low ]]; then return ' + str(low_status) + '; fi\n'
            '  return ' + str(high_status) + '\n'
            '}\n'
        )
        return subprocess.run(
            [bash, "-c", harness + self._full_low_high_loop(), "loop-test", debug, inherited_selector],
            text=True, capture_output=True, check=False,
        )

    def test_debug_neutralizes_inherited_selector_and_runs_full_low_once(self):
        for inherited in ("", "1", "true", "invalid"):
            with self.subTest(inherited=inherited):
                result = self._exercise_full_loop("true", inherited)
                self.assertEqual(result.returncode, 0, result.stderr)
                calls = [line for line in result.stdout.splitlines() if line.startswith("CALL:")]
                self.assertEqual(calls, ["CALL:low::low input"])
                self.assertIn("All available low-epsilon settings completed", result.stdout)

    def test_full_loop_failures_stop_before_success_or_next_epsilon(self):
        for debug in ("true", "false"):
            with self.subTest(debug=debug):
                result = self._exercise_full_loop(debug, "", low_status=7)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(result.stdout.splitlines(), ["CALL:low::low input"])
        result = self._exercise_full_loop("false", "", high_status=9)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertNotIn("High Epsilon Completed!", result.stdout)

    def test_ordinary_loop_keeps_both_calls_and_inherited_environment(self):
        for inherited in ("", "1"):
            with self.subTest(inherited=inherited):
                result = self._exercise_full_loop("false", inherited)
                self.assertEqual(result.returncode, 0, result.stderr)
                calls = [line for line in result.stdout.splitlines() if line.startswith("CALL:")]
                self.assertEqual(calls, [
                    "CALL:low:{}:low input".format(inherited),
                    "CALL:high:{}:high input".format(inherited),
                ])

    def test_main_ordinary_population_uses_existing_availability_gate(self):
        """Evaluate the actual selector/filter AST without importing ROOT or main."""
        nodes = ast.parse(self.main).body
        ordinary = next(node for node in nodes if isinstance(node, ast.Assign)
                        and any(isinstance(t, ast.Name) and t.id == "phisetlist" for t in node.targets))
        restriction = next(node for node in nodes if isinstance(node, ast.If)
                           and isinstance(node.test, ast.Name) and node.test.id == "debug_left_low"
                           and any(isinstance(n, ast.Assign) and any(
                               isinstance(t, ast.Name) and t.id == "phisetlist" for t in n.targets
                           ) for n in node.body))
        availability = next(node for node in nodes if isinstance(node, ast.Assign)
                            and any(isinstance(t, ast.Name) and t.id == "available_phisetlist"
                                    for t in node.targets))
        code = compile(ast.Module(body=[ordinary, restriction, availability], type_ignores=[]),
                       str(MAIN_PATH), "exec")
        for populated in (["Center", "Left"], ["Center", "Left", "Right"]):
            with self.subTest(populated=populated):
                scope = {"debug_left_low": "", "mm_shift_summary": dict.fromkeys(populated),
                         "kaon_data_coordinate_summary": dict.fromkeys(populated)}
                exec(code, scope)
                self.assertEqual(scope["available_phisetlist"], populated)

    def test_main_debug_selector_preserves_diamond_before_left_only_downstream_work(self):
        ast.parse(self.main, filename=str(MAIN_PATH))
        ordinary_settings = self.main.index('phisetlist = ["Center", "Left", "Right"]')
        validation_start = self.main.index(
            'debug_left_low = os.environ.get("LT_ANALYSIS_DEBUG_LEFT_LOW"'
        )
        diamond_import = self.main.index("from diamond import DiamondPlot")
        diamond_loop = self.main.index("for phiset in phisetlist:", diamond_import)
        diamond_call = self.main.index("inpDict.update(DiamondPlot", diamond_loop)
        diamond_registration = self.main.index(
            'output_file_lst.append(OUTPATH+"/{}_{}_diamond_{}', diamond_call
        )
        debug_restriction = self.main.rindex("if debug_left_low:")
        shift_import = self.main.index("from shift_prep import shift_prep", debug_restriction)
        shift_loop = self.main.index("for phiset in phisetlist:", shift_import)
        early_validation = self.main[validation_start:diamond_import]
        downstream_selector = self.main[debug_restriction:shift_import]

        self.assertLess(ordinary_settings, validation_start)
        self.assertIn('if debug_left_low:', early_validation)
        self.assertIn('ParticleType != "kaon"', early_validation)
        self.assertIn('str(EPSSET).strip().lower() != "low"', early_validation)
        self.assertIn('LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE', early_validation)
        self.assertNotIn('phisetlist = ["Left"]', early_validation)
        self.assertLess(diamond_loop, diamond_call)
        self.assertLess(diamond_call, diamond_registration)
        self.assertLess(diamond_registration, debug_restriction)
        self.assertLess(debug_restriction, shift_loop)
        self.assertIn('phisetlist = ["Left"]', downstream_selector)
        self.assertNotIn("inpDict", downstream_selector)
        self.assertNotIn("canonical_t_binning", downstream_selector)

    def test_diamond_source_is_unchanged(self):
        result = subprocess.run(
            [
                "git",
                "-c",
                "core.autocrlf=false",
                "-c",
                "core.safecrlf=false",
                "diff",
                "--quiet",
                "HEAD",
                "--",
                str(DIAMOND_PATH.relative_to(REPO_ROOT)),
            ],
            cwd=REPO_ROOT,
            check=False,
        )
        self.assertEqual(result.returncode, 0, "src/cuts/diamond.py must remain unchanged")


if __name__ == "__main__":
    unittest.main()
