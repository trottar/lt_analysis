"""Regression coverage for the external-only K-Sigma0 SIMC source."""

from __future__ import annotations

import copy
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from array import array
from pathlib import Path
from unittest import mock


ROOT = None
try:
    import ROOT  # type: ignore
except ImportError:
    ROOT = None


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

from background_sample_manifest import (
    BACKGROUND_SAMPLE_SUFFIXES,
    build_background_sample_manifest,
    sigma0_environment_variable,
)


class Sigma0ManifestTests(unittest.TestCase):
    def test_unconfigured_sigma0_keeps_external_manifest_schema(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "low")

        entry = manifest["by_phi"]["Left"]["sigma0"]
        self.assertEqual(entry["source_strategy"], "external_required")
        self.assertFalse(entry["configured"])
        self.assertIsNone(entry["root"])
        self.assertEqual(
            entry["source_identity"],
            {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "phi_setting": "Left"},
        )
        self.assertEqual(
            set(BACKGROUND_SAMPLE_SUFFIXES)
            | {"source_strategy", "configured", "environment_variable", "source_identity"},
            set(entry),
        )
        self.assertEqual(entry["environment_variable"], "LT_BG_SIGMA0_LEFT_LOW_ROOT")
        self.assertTrue(all(entry[key] is None for key in BACKGROUND_SAMPLE_SUFFIXES))

    def test_low_and_high_phi_roots_are_independent(self):
        with mock.patch.dict(
            os.environ,
            {
                "LT_BG_SIGMA0_LEFT_LOW_ROOT": "D:/sigma0/Q4p4W2p74left_low.root",
                "LT_BG_SIGMA0_LEFT_HIGH_ROOT": "D:/sigma0/Q4p4W2p74left_high.root",
                "LT_BG_SIGMA0_RIGHT_LOW_ROOT": "D:/sigma0/Q4p4W2p74right_low.root",
            },
            clear=True,
        ):
            low_manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "low")
            high_manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "high")

        low_left = low_manifest["by_phi"]["Left"]["sigma0"]
        high_left = high_manifest["by_phi"]["Left"]["sigma0"]
        low_right = low_manifest["by_phi"]["Right"]["sigma0"]
        self.assertEqual(low_left["root"], "D:/sigma0/Q4p4W2p74left_low.root")
        self.assertEqual(low_left["environment_variable"], "LT_BG_SIGMA0_LEFT_LOW_ROOT")
        self.assertEqual(high_left["root"], "D:/sigma0/Q4p4W2p74left_high.root")
        self.assertEqual(high_left["environment_variable"], "LT_BG_SIGMA0_LEFT_HIGH_ROOT")
        self.assertEqual(low_right["root"], "D:/sigma0/Q4p4W2p74right_low.root")
        self.assertEqual(low_right["source_identity"]["phi_setting"], "Right")
        self.assertEqual(high_left["source_identity"]["EPSSET"], "high")

    def test_low_only_configuration_leaves_high_unconfigured(self):
        with mock.patch.dict(
            os.environ,
            {"LT_BG_SIGMA0_LEFT_LOW_ROOT": "D:/sigma0/Q4p4W2p74left_low.root"},
            clear=True,
        ):
            high_manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "high")

        high_left = high_manifest["by_phi"]["Left"]["sigma0"]
        self.assertFalse(high_left["configured"])
        self.assertIsNone(high_left["root"])
        self.assertEqual(high_left["environment_variable"], "LT_BG_SIGMA0_LEFT_HIGH_ROOT")

    def test_generated_epsilon_override_does_not_retarget_sigma0(self):
        with mock.patch.dict(
            os.environ,
            {
                "LT_BG_SAMPLE_EPSILON": "high",
                "LT_BG_SIGMA0_CENTER_LOW_ROOT": "D:/sigma0/Q4p4W2p74center_low.root",
                "LT_BG_SIGMA0_CENTER_HIGH_ROOT": "D:/sigma0/Q4p4W2p74center_high.root",
            },
            clear=True,
        ):
            manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "low")

        center = manifest["by_phi"]["Center"]["sigma0"]
        self.assertEqual(manifest["epsilon"], "high")
        self.assertEqual(
            manifest["by_phi"]["Center"]["neutron"]["source_identity"]["EPSSET"],
            "high",
        )
        self.assertEqual(center["environment_variable"], "LT_BG_SIGMA0_CENTER_LOW_ROOT")
        self.assertEqual(center["root"], "D:/sigma0/Q4p4W2p74center_low.root")
        self.assertEqual(center["source_identity"]["EPSSET"], "low")

    def test_analysis_identity_isolated_from_generated_overrides(self):
        with mock.patch.dict(
            os.environ,
            {
                "LT_BG_SAMPLE_Q2": "3p0",
                "LT_BG_SAMPLE_W": "3p14",
                "LT_BG_SAMPLE_EPSILON": "low",
                "LT_BG_SIGMA0_LEFT_LOW_ROOT": "D:/sigma0/Q3p0W3p14left_low.root",
                "LT_BG_SIGMA0_LEFT_HIGH_ROOT": "D:/sigma0/Q4p4W2p74left_high.root",
            },
            clear=True,
        ):
            manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "high")

        left_sigma0 = manifest["by_phi"]["Left"]["sigma0"]
        self.assertEqual(manifest["q2"], "3p0")
        self.assertEqual(manifest["w"], "3p14")
        self.assertEqual(manifest["epsilon"], "low")
        for background in ("neutron", "delta", "sidis"):
            generated = manifest["by_phi"]["Left"][background]
            self.assertEqual(
                generated["source_identity"],
                {"Q2": "3p0", "W": "3p14", "EPSSET": "low", "phi_setting": "Left"},
            )
            self.assertIn("Prod_Coin_Q3p0W3p14left_lowe.root", generated["root"])
        self.assertEqual(left_sigma0["environment_variable"], "LT_BG_SIGMA0_LEFT_HIGH_ROOT")
        self.assertEqual(left_sigma0["root"], "D:/sigma0/Q4p4W2p74left_high.root")
        self.assertEqual(
            left_sigma0["source_identity"],
            {"Q2": "4p4", "W": "2p74", "EPSSET": "high", "phi_setting": "Left"},
        )

    def test_generic_phi_only_variable_is_ignored(self):
        with mock.patch.dict(
            os.environ,
            {"LT_BG_SIGMA0_LEFT_ROOT": "D:/sigma0/legacy_left.root"},
            clear=True,
        ):
            manifest = build_background_sample_manifest("C:/analysis", "4p4", "2p74", "low")

        left = manifest["by_phi"]["Left"]["sigma0"]
        self.assertFalse(left["configured"])
        self.assertIsNone(left["root"])
        self.assertEqual(left["environment_variable"], "LT_BG_SIGMA0_LEFT_LOW_ROOT")

    def test_only_low_and_high_epsilon_tokens_are_accepted(self):
        self.assertEqual(sigma0_environment_variable("Center", "LOW"), "LT_BG_SIGMA0_CENTER_LOW_ROOT")
        with self.assertRaises(ValueError):
            sigma0_environment_variable("Left", "mid")

    @staticmethod
    def _launcher_bash():
        if os.name == "nt":
            git_bash = Path(r"C:\Program Files\Git\bin\bash.exe")
            return str(git_bash) if git_bash.is_file() else None
        return shutil.which("bash")

    def _run_launcher_cleanup_guard(self, worktree_root, configured):
        bash = self._launcher_bash()
        if bash is None:
            self.skipTest("a local Bash executable is required for launcher guard coverage")
        environment = os.environ.copy()
        for epsilon in ("LOW", "HIGH"):
            for phi in ("RIGHT", "LEFT", "CENTER"):
                environment.pop("LT_BG_SIGMA0_{}_{}_ROOT".format(phi, epsilon), None)
        environment.update(configured)
        shell_driver = r'''
to_shell_path() {
    if command -v cygpath >/dev/null 2>&1; then
        cygpath -u "$1"
    else
        printf '%s\n' "$1"
    fi
}
source <(
    sed -n \
        -e '/^resolve_real_path()/,/^}/p' \
        -e '/^sigma0_external_environment_variable()/,/^}/p' \
        -e '/^path_is_within_worktree()/,/^}/p' \
        -e '/^validate_external_sigma0_paths_before_cleanup()/,/^}/p' \
        "$1"
)
LTANAPATH="$(to_shell_path "$2")"
for sigma0_variable in LT_BG_SIGMA0_RIGHT_LOW_ROOT LT_BG_SIGMA0_LEFT_LOW_ROOT LT_BG_SIGMA0_CENTER_LOW_ROOT LT_BG_SIGMA0_RIGHT_HIGH_ROOT LT_BG_SIGMA0_LEFT_HIGH_ROOT LT_BG_SIGMA0_CENTER_HIGH_ROOT; do
    if [[ -n "${!sigma0_variable:-}" ]]; then
        printf -v "${sigma0_variable}" '%s' "$(to_shell_path "${!sigma0_variable}")"
    fi
done
validate_external_sigma0_paths_before_cleanup
'''
        return subprocess.run(
            [bash, "-c", shell_driver, "bash", str(REPO_ROOT / "run_Prod_Analysis.sh"), str(worktree_root)],
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_launcher_cleanup_guard_rejects_in_tree_and_allows_external_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            worktree_root = Path(directory) / "repo"
            worktree_root.mkdir()
            blocked = self._run_launcher_cleanup_guard(
                worktree_root,
                {"LT_BG_SIGMA0_LEFT_LOW_ROOT": str(worktree_root / "sigma0" / "left.root")},
            )
            outside_link = Path(directory) / "external_link_to_repo"
            try:
                outside_link.symlink_to(worktree_root, target_is_directory=True)
            except OSError:
                symlink_blocked = None
            else:
                symlink_blocked = self._run_launcher_cleanup_guard(
                    worktree_root,
                    {"LT_BG_SIGMA0_RIGHT_LOW_ROOT": str(outside_link / "sigma0" / "right.root")},
                )
            volatile = self._run_launcher_cleanup_guard(
                worktree_root,
                {"LT_BG_SIGMA0_RIGHT_HIGH_ROOT": "/volatile/sigma0/right.root"},
            )
            work = self._run_launcher_cleanup_guard(
                worktree_root,
                {"LT_BG_SIGMA0_CENTER_HIGH_ROOT": "/work/sigma0/center.root"},
            )
            legacy = self._run_launcher_cleanup_guard(
                worktree_root,
                {"LT_BG_SIGMA0_LEFT_ROOT": str(worktree_root / "legacy.root")},
            )

        self.assertNotEqual(blocked.returncode, 0)
        self.assertIn("LT_BG_SIGMA0_LEFT_LOW_ROOT", blocked.stderr)
        if symlink_blocked is not None:
            self.assertNotEqual(symlink_blocked.returncode, 0)
            self.assertIn("LT_BG_SIGMA0_RIGHT_LOW_ROOT", symlink_blocked.stderr)
        self.assertEqual(volatile.returncode, 0, volatile.stderr)
        self.assertEqual(work.returncode, 0, work.stderr)
        self.assertEqual(legacy.returncode, 0, legacy.stderr)

    def test_launcher_preserves_external_sigma0_contract_and_cleanup_order(self):
        launcher = (REPO_ROOT / "run_Prod_Analysis.sh").read_text(encoding="utf-8")
        function_body = launcher.split("export_background_sample_paths()", 1)[1].split("SKIM_OUTPUT_DIR=", 1)[0]
        self.assertIn("for background in neutron delta sidis; do", function_body)
        self.assertNotRegex(function_body, r"for background in .*sigma0")
        self.assertIn('report_external_sigma0_sources "${eps_token}"', function_body)
        self.assertNotRegex(function_body, r"export .*LT_BG_SIGMA0")
        self.assertEqual(launcher.count('export_background_sample_paths "${Q2}" "${W}" "${EPSILON}"'), 2)

        cleanup_block = launcher.split("# Clean all untracked files and recreate symlinks", 1)[1]
        self.assertLess(
            cleanup_block.index("validate_external_sigma0_paths_before_cleanup"),
            cleanup_block.index("git clean -fdx"),
        )
        guard_body = launcher.split("validate_external_sigma0_paths_before_cleanup()", 1)[1].split(
            "export_background_sample_paths()", 1
        )[0]
        self.assertIn("for epsilon in low high; do", guard_body)
        self.assertIn("for phi_setting in Right Left Center; do", guard_body)
        self.assertIn('resolve_real_path "${LTANAPATH}"', guard_body)
        self.assertIn('resolve_real_path "${configured_root}"', guard_body)


@unittest.skipUnless(ROOT is not None, "PyROOT is required for K-Sigma0 ROOT source tests")
class Sigma0ResolverAndLoaderTests(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        import pion_component_shapes as shapes

        self.shapes = shapes
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.hole_cut = type("HoleCut", (), {"IsInside": staticmethod(lambda *_args: False)})()

    def tearDown(self):
        self.tempdir.cleanup()

    def _inp(self, root_filename=None):
        return {
            "Q2": "4p4",
            "W": "2p74",
            "EPSSET": "low",
            "POL": "1",
            "ParticleType": "kaon",
            "tmin": 0.0,
            "tmax": 1.0,
            "mm_min": 0.70,
            "mm_max": 1.30,
            "bg_opt_mm_plot_min": 0.70,
            "bg_opt_mm_plot_max": 1.30,
            "bg_opt_mm_plot_nbins": 100,
            "background_samples": {
                "by_phi": {
                    "Left": {
                        "sigma0": {
                            "source_strategy": "external_required",
                            "configured": bool(root_filename),
                            "environment_variable": "LT_BG_SIGMA0_LEFT_LOW_ROOT",
                            "root": root_filename,
                            "source_identity": {
                                "Q2": "4p4",
                                "W": "2p74",
                                "EPSSET": "low",
                                "phi_setting": "Left",
                            },
                        }
                    }
                }
            },
        }

    def _write_root(self, name, include_tree=True, missing_branch=None, entries=120):
        filename = self.tmp / name
        filename.parent.mkdir(parents=True, exist_ok=True)
        root_file = ROOT.TFile(str(filename), "RECREATE")
        if include_tree:
            tree = ROOT.TTree("h10", "h10")
            values = {
                "missmass": 1.195,
                "t": -0.50,
                "ssdelta": 0.0,
                "ssxptar": 0.0,
                "ssyptar": 0.0,
                "hsdelta": 0.0,
                "hsxptar": 0.0,
                "hsyptar": 0.0,
                "Q2": 4.4,
                "W": 2.74,
                "phgcer_x_det": 100.0,
                "phgcer_y_det": 100.0,
                "phipq": 0.0,
            }
            holders = {}
            for branch_name, value in values.items():
                if branch_name == missing_branch:
                    continue
                holders[branch_name] = array("d", [value])
                tree.Branch(branch_name, holders[branch_name], "{}/D".format(branch_name))
            for index in range(entries):
                if "missmass" in holders:
                    x_value = 1.195 + 0.030 * math.sin(float(index) * 0.23)
                    holders["missmass"][0] = x_value
                tree.Fill()
            tree.Write()
        root_file.Close()
        return str(filename)

    def _load(self, root_filename, **kwargs):
        with mock.patch.object(self.shapes, "set_val"), mock.patch.object(
            self.shapes, "apply_simc_cuts", return_value=True
        ), mock.patch.object(self.shapes, "apply_simc_sub_cuts", return_value=True):
            return self.shapes.load_kaon_simc_sigma0_shape(
                None,
                self._inp(root_filename),
                "Left",
                hgcer_cutg=self.hole_cut,
                context="sigma0_source_test",
                **kwargs,
            )

    def test_explicit_valid_source_loads_and_reports_identity(self):
        root_filename = self._write_root("valid_sigma0.root")
        payload = self._load(root_filename)
        diagnostics = payload["diagnostics"]

        self.assertEqual(diagnostics["requested_root"], root_filename)
        self.assertEqual(diagnostics["resolved_root"], root_filename)
        self.assertEqual(diagnostics["resolution_source"], "explicit_configured_root")
        self.assertTrue(diagnostics["path_exists"])
        self.assertTrue(diagnostics["root_open_success"])
        self.assertTrue(diagnostics["tree_exists"])
        self.assertGreater(diagnostics["tree_entries"], 0)
        self.assertGreater(diagnostics["n_events_seen"], 0)
        self.assertGreater(diagnostics["n_events_passed"], 0)
        self.assertGreater(diagnostics["weighted_integral_before_norm"], 0.0)
        self.assertTrue(diagnostics["normalized"])
        self.assertFalse(diagnostics["fallback_used"])
        self.assertGreater(payload["setting_shape_full"].Integral(), 0.0)
        self.assertEqual(diagnostics["source_identity"]["phi_setting"], "Left")
        self.assertEqual(diagnostics["source_identity"]["EPSSET"], "low")
        self.assertEqual(diagnostics["requested_environment_variable"], "LT_BG_SIGMA0_LEFT_LOW_ROOT")

    def test_unconfigured_and_missing_sources_are_distinct(self):
        unconfigured = self._load(None)["diagnostics"]
        self.assertEqual(unconfigured["fallback_reason"], "no_source_configured")
        self.assertFalse(unconfigured["path_exists"])

        missing = self._load(str(self.tmp / "missing_sigma0.root"))["diagnostics"]
        self.assertEqual(missing["fallback_reason"], "configured_path_does_not_exist")
        self.assertFalse(missing["path_exists"])

    def test_missing_tree_zero_entry_tree_and_missing_branch_are_distinct(self):
        missing_tree = self._load(self._write_root("no_tree.root", include_tree=False))["diagnostics"]
        self.assertTrue(missing_tree["root_open_success"])
        self.assertFalse(missing_tree["tree_exists"])
        self.assertEqual(missing_tree["fallback_reason"], "missing_simc_tree")

        zero_tree = self._load(self._write_root("zero_tree.root", entries=0))["diagnostics"]
        self.assertTrue(zero_tree["tree_exists"])
        self.assertEqual(zero_tree["tree_entries"], 0)
        self.assertEqual(zero_tree["fallback_reason"], "zero_entry_tree")

        missing_branch = self._load(self._write_root("missing_branch.root", missing_branch="hsdelta"))["diagnostics"]
        self.assertTrue(missing_branch["tree_exists"])
        self.assertIn("hsdelta", missing_branch["missing_required_branches"])
        self.assertEqual(missing_branch["fallback_reason"], "incompatible_tree_missing_branches")

    def test_resolver_never_discovers_a_valid_unconfigured_generated_path(self):
        generated_like_root = self._write_root(
            "background_samples/OUTPUTS/sigma0/Prod_Coin_Q4p4W2p74left_lowe.root"
        )
        resolved, provenance = self.shapes.resolve_kaon_simc_sigma0_root_filename(
            None,
            self._inp(None),
            "Left",
        )
        self.assertTrue(os.path.isfile(generated_like_root))
        self.assertIsNone(resolved)
        self.assertEqual(provenance["candidate_roots"], [])
        self.assertEqual(provenance["resolution_source"], "no_source_configured")

    @staticmethod
    def _shape(name, mean, width):
        hist = ROOT.TH1D(name, name, 100, 0.70, 1.30)
        for index in range(1, hist.GetNbinsX() + 1):
            x_value = hist.GetBinCenter(index)
            hist.SetBinContent(index, math.exp(-0.5 * ((x_value - mean) / width) ** 2))
            hist.SetBinError(index, 0.02)
        hist.SetDirectory(0)
        return hist

    def test_explicit_source_reaches_component_builder_and_protected_fit(self):
        import background_config as bgcfg
        import pion_component_fits as fits

        root_filename = self._write_root("plumbing_sigma0.root")
        with mock.patch.dict(
            os.environ,
            {"LT_BG_SIGMA0_LEFT_LOW_ROOT": root_filename},
            clear=True,
        ):
            input_dict = self._inp(None)
            input_dict["background_samples"] = build_background_sample_manifest(
                "C:/analysis",
                "4p4",
                "2p74",
                "low",
            )
            with mock.patch.object(self.shapes, "set_val"), mock.patch.object(
                self.shapes, "apply_simc_cuts", return_value=True
            ), mock.patch.object(self.shapes, "apply_simc_sub_cuts", return_value=True):
                sigma0_payload = self.shapes.load_kaon_simc_sigma0_shape(
                    None,
                    input_dict,
                    "Left",
                    hgcer_cutg=self.hole_cut,
                    context="sigma0_full_plumbing",
                )
        sigma0 = sigma0_payload["setting_shape_full"]
        pi_n = self._shape("plumbing_pi_n", 0.90, 0.018)
        pi_sidis = self._shape("plumbing_pi_sidis", 1.07, 0.030)
        pi_delta = self._shape("plumbing_pi_delta", 1.225, 0.050)
        k_lambda = self._shape("plumbing_k_lambda", 1.115, 0.012)
        target = self._shape("plumbing_target", 0.70, 1.0)
        target.Reset()
        for hist, amplitude in (
            (pi_n, 1.10),
            (pi_sidis, 0.60),
            (k_lambda, 2.00),
            (sigma0, 0.80),
            (pi_delta, 0.25),
        ):
            target.Add(hist, amplitude)
        for index in range(1, target.GetNbinsX() + 1):
            target.SetBinError(index, 0.02)

        original = copy.deepcopy(bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"])
        try:
            config = bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"]
            config["joint_refinement_enabled"] = False
            config["residual_component_shifts_enabled"] = False
            config["pi_delta_signal_protected_fit"]["fit_window"] = None
            result = fits.build_particle_subtraction_component_result(
                target,
                target,
                {"pi_n": pi_n, "pi_delta": pi_delta, "pi_sidis": pi_sidis},
                {
                    "particle_subtraction_mode": "simc_shape_components",
                    "bg_opt_mm_plot_min": 0.70,
                    "bg_opt_mm_plot_max": 1.30,
                },
                "setting-wide",
                kaon_signal_shape=k_lambda,
                kaon_sigma0_shape=sigma0,
                phi_setting="Left",
                context="sigma0_full_plumbing",
            )
        finally:
            bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"] = original

        protected = result["diagnostics"]["kaon"]["pi_delta_signal_protected_fit"]
        self.assertEqual(protected["status"], "success")
        self.assertTrue(protected["template_availability"]["k_lambda_signal"])
        self.assertTrue(protected["template_availability"]["k_sigma0_signal"])
        self.assertTrue(protected["template_availability"]["pi_delta"])
        self.assertTrue(math.isfinite(result["A_delta"]))
        self.assertGreater(result["A_delta"], 0.0)
        self.assertEqual(
            sigma0_payload["diagnostics"]["requested_environment_variable"],
            "LT_BG_SIGMA0_LEFT_LOW_ROOT",
        )


if __name__ == "__main__":
    unittest.main()
