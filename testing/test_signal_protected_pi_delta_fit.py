"""Focused regressions for the signal-protected kaon-side pi-delta fit."""

from __future__ import annotations

import copy
import math
import os
import sys
import unittest


ROOT = None
try:
    import ROOT  # type: ignore
except ImportError:
    ROOT = None


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for relative_path in ("src/cuts", "src/utility"):
    path = os.path.join(REPO_ROOT, relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import background_config as bgcfg


class SignalProtectedPiDeltaConfigTests(unittest.TestCase):
    def test_default_and_partial_override_deep_merge(self):
        default = bgcfg.get_particle_subtraction_component_fit_window_config("kaon_nosub")
        protected = default["pi_delta_signal_protected_fit"]
        self.assertTrue(protected["enabled"])
        self.assertTrue(protected["require_k_lambda_template"])
        self.assertTrue(protected["require_k_sigma0_template"])
        self.assertEqual(protected["failure_policy"], "zero_pi_delta")

        merged = bgcfg._deep_merge_particle_subtraction_config(
            default,
            {"pi_delta_signal_protected_fit": {"template_corr_warn": 0.91}},
        )
        merged_protected = merged["pi_delta_signal_protected_fit"]
        self.assertAlmostEqual(merged_protected["template_corr_warn"], 0.91)
        self.assertTrue(merged_protected["require_k_lambda_template"])
        self.assertEqual(merged_protected["failure_policy"], "zero_pi_delta")


@unittest.skipUnless(ROOT is not None, "PyROOT is required for protected-fit histogram tests")
class SignalProtectedPiDeltaFitTests(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        import pion_component_fits as fits

        self.fits = fits
        self.config = copy.deepcopy(
            bgcfg.get_particle_subtraction_component_fit_window_config("kaon_nosub")
        )
        self.config["pi_delta_signal_protected_fit"].update(
            {"fit_window": None, "template_corr_warn": 0.95}
        )
        self.inp = {"bg_opt_mm_plot_min": 0.70, "bg_opt_mm_plot_max": 1.30}
        self.pi_n = self._shape("pi_n", 0.90, 0.018)
        self.pi_sidis = self._shape("pi_sidis", 1.07, 0.030)
        self.k_lambda = self._shape("k_lambda", 1.115, 0.012)
        self.k_sigma = self._shape("k_sigma", 1.195, 0.018)
        self.pi_delta = self._shape("pi_delta", 1.225, 0.050)

    def _shape(self, name, mean, width):
        hist = ROOT.TH1D(name, name, 120, 0.70, 1.30)
        for index in range(1, hist.GetNbinsX() + 1):
            x_value = hist.GetBinCenter(index)
            value = math.exp(-0.5 * ((x_value - mean) / width) ** 2)
            hist.SetBinContent(index, value)
            hist.SetBinError(index, 0.0)
        hist.SetDirectory(0)
        return hist

    def _target(self, k_sigma=None):
        target = self.pi_n.Clone("kaon_target")
        target.SetDirectory(0)
        target.Reset()
        for hist, amplitude in (
            (self.pi_n, 1.10),
            (self.pi_sidis, 0.60),
            (self.k_lambda, 2.00),
            (self.k_sigma if k_sigma is None else k_sigma, 0.80),
            (self.pi_delta, 0.25),
        ):
            target.Add(hist, amplitude)
        for index in range(1, target.GetNbinsX() + 1):
            target.SetBinError(index, 0.02)
        return target

    def _legacy_payload(self):
        return {
            "A_n": 1.10,
            "A_sidis": 0.60,
            "A_delta": 4.00,
            "k_lambda_reference_scale": 7.50,
            "H_k_lambda_simc_reference": self.k_lambda.Clone("lambda_reference"),
            "diagnostics": {"refined_amplitudes": {"pi_n": 1.10, "pi_sidis": 0.60, "pi_delta": 4.00}},
        }

    def _apply(self, target, sigma=None):
        return self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(),
            target,
            self.pi_n,
            self.pi_delta,
            self.pi_sidis,
            self.k_lambda,
            self.k_sigma if sigma is None else sigma,
            self.config,
            0.70,
            1.30,
            [],
            self.inp,
            "Center",
            0.0,
            context="unit",
        )

    def test_recovers_three_template_solution_and_subtracts_only_pi_delta(self):
        result = self._apply(self._target())
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertAlmostEqual(result["A_n"], 1.10, places=10)
        self.assertAlmostEqual(result["A_sidis"], 0.60, places=10)
        self.assertAlmostEqual(result["A_delta"], 0.25, places=7)
        self.assertAlmostEqual(result["S_lambda"], 2.00, places=7)
        self.assertAlmostEqual(result["S_sigma0"], 0.80, places=7)
        self.assertAlmostEqual(result["k_lambda_reference_scale"], 7.50, places=12)
        self.assertTrue(diagnostic["lambda_reference_integrity"]["shape_identical"])
        self.assertTrue(diagnostic["fit_residual_diagnostic_only"])
        self.assertEqual(diagnostic["physics_output"], "H_pi_delta_protected_after_subtraction")

        for index in range(1, self.pi_delta.GetNbinsX() + 1):
            expected = (
                2.00 * self.k_lambda.GetBinContent(index)
                + 0.80 * self.k_sigma.GetBinContent(index)
            )
            self.assertAlmostEqual(
                result["H_pi_delta_protected_after_subtraction"].GetBinContent(index),
                expected,
                places=7,
            )
            self.assertAlmostEqual(
                result["H_pi_delta_protected_fit_residual"].GetBinContent(index),
                0.0,
                places=7,
            )
        self.assertAlmostEqual(
            result["diagnostics"]["refined_amplitudes"]["pi_delta"], result["A_delta"], places=12
        )

    def test_missing_sigma_zeroes_only_final_pi_delta(self):
        # The data can contain Sigma strength while the required SIMC template
        # is genuinely absent; this must not revive the legacy pi-delta fit.
        result = self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(), self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            self.k_lambda, None, self.config, 0.70, 1.30, [], self.inp, "Center", 0.0, context="missing",
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "missing_required_template")
        self.assertEqual(result["A_delta"], 0.0)
        self.assertAlmostEqual(result["A_n"], 1.10, places=12)
        self.assertAlmostEqual(result["A_sidis"], 0.60, places=12)
        self.assertEqual(diagnostic["legacy_staged_A_delta"], 4.00)

    def test_rank_deficiency_is_not_reported_as_high_correlation_warning_only(self):
        result = self._apply(self._target(k_sigma=self.pi_delta), sigma=self.pi_delta)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "rank_deficient")
        self.assertEqual(result["A_delta"], 0.0)

    def test_strong_sigma_delta_overlap_remains_a_solvable_protected_fit(self):
        overlapping_sigma = self._shape("k_sigma_overlap", 1.210, 0.045)
        result = self._apply(self._target(k_sigma=overlapping_sigma), sigma=overlapping_sigma)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertAlmostEqual(result["A_delta"], 0.25, places=6)
        self.assertAlmostEqual(result["S_sigma0"], 0.80, places=6)


if __name__ == "__main__":
    unittest.main()
