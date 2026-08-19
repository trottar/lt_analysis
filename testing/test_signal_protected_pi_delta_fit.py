"""Focused regressions for the signal-protected kaon-side pi-delta fit."""

from __future__ import annotations

import copy
import math
import os
import sys
import unittest
from unittest import mock


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
        self.assertTrue(protected["allow_lambda_only_fallback"])
        self.assertEqual(protected["failure_policy"], "zero_pi_delta")
        self.assertEqual(protected["fit_window"], (1.10, 1.30))
        self.assertIsNone(protected["lambda_gauge_window"])
        self.assertEqual(protected["lambda_gauge_constraint_mode"], "gaussian")
        self.assertAlmostEqual(protected["lambda_gauge_min_relative_uncertainty"], 0.05)
        self.assertAlmostEqual(protected["lambda_gauge_poor_relative_uncertainty"], 0.35)
        self.assertAlmostEqual(protected["lambda_gauge_min_retention_fraction"], 0.90)
        self.assertAlmostEqual(protected["maximum_chi2_ndf"], 5.0)
        self.assertAlmostEqual(protected["minimum_p_value"], 1.0e-6)
        self.assertAlmostEqual(protected["lambda_gauge_maximum_chi2_ndf"], 10.0)
        self.assertIsNone(protected["lambda_gauge_minimum_p_value"])
        self.assertEqual(protected["lambda_gauge_min_fit_bins"], 2)
        self.assertFalse(bgcfg.resolve_particle_subtraction_root_ownership_debug())
        self.assertTrue(bgcfg.resolve_particle_subtraction_root_ownership_debug(
            {"particle_subtraction_root_ownership_debug": True}
        ))

        merged = bgcfg._deep_merge_particle_subtraction_config(
            default,
            {"pi_delta_signal_protected_fit": {"template_corr_warn": 0.91}},
        )
        merged_protected = merged["pi_delta_signal_protected_fit"]
        self.assertAlmostEqual(merged_protected["template_corr_warn"], 0.91)
        self.assertTrue(merged_protected["require_k_lambda_template"])
        self.assertTrue(merged_protected["allow_lambda_only_fallback"])
        self.assertEqual(merged_protected["failure_policy"], "zero_pi_delta")


@unittest.skipUnless(ROOT is not None, "PyROOT is required for protected-fit histogram tests")
class SignalProtectedPiDeltaFitTests(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        import pion_component_fits as fits
        import pion_component_subtraction as subtraction

        self.fits = fits
        self.subtraction = subtraction
        self.config = copy.deepcopy(
            bgcfg.get_particle_subtraction_component_fit_window_config("kaon_nosub")
        )
        self.config["pi_delta_signal_protected_fit"].update(
            {"fit_window": (1.10, 1.30), "template_corr_warn": 0.95}
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

    def test_single_pair_protected_fit_window_is_normalized(self):
        resolved = self.fits._resolve_pi_delta_signal_protected_fit_config(
            {"pi_delta_signal_protected_fit": {"fit_window": (1.10, 1.30)}}
        )
        self.assertEqual(resolved["fit_window_collection"], [(1.10, 1.30)])

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

    def _apply(self, target, sigma=None, source_diagnostics=None):
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
            kaon_sigma0_source_diagnostics=source_diagnostics,
        )

    def test_recovers_three_template_solution_and_subtracts_only_pi_delta(self):
        result = self._apply(self._target())
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertEqual(diagnostic["fit_variant"], "lambda_sigma0_protected")
        self.assertFalse(diagnostic["fallback_attempted"])
        self.assertFalse(diagnostic["fallback_used"])
        self.assertTrue(diagnostic["sigma0_fitted"])
        self.assertAlmostEqual(result["A_n"], 1.10, places=10)
        self.assertAlmostEqual(result["A_sidis"], 0.60, places=10)
        self.assertAlmostEqual(result["A_delta"], 0.25, places=7)
        self.assertAlmostEqual(result["S_lambda"], 2.00, places=7)
        self.assertAlmostEqual(result["S_sigma0"], 0.80, places=7)
        self.assertAlmostEqual(result["k_lambda_reference_scale"], 7.50, places=12)
        self.assertTrue(diagnostic["lambda_reference_integrity"]["shape_identical"])
        self.assertTrue(diagnostic["fit_residual_diagnostic_only"])
        self.assertEqual(diagnostic["physics_output"], "H_pi_delta_protected_after_subtraction")
        self.assertTrue(diagnostic["closure"]["full_five_component_model"]["passed"])
        self.assertTrue(diagnostic["closure"]["protected_three_component_model"]["passed"])
        self.assertTrue(diagnostic["closure"]["delta_only_physics_output"]["passed"])
        self.assertTrue(diagnostic["early_amplitudes_frozen_integrity"]["unchanged"])
        self.assertTrue(diagnostic["lambda_gauge_solver_success"])
        self.assertTrue(diagnostic["lambda_gauge_quality_passed"])
        self.assertEqual(diagnostic["lambda_gauge_status"], "success")
        self.assertEqual(diagnostic["constraint_metrics"]["mode"], "gaussian")
        self.assertEqual(diagnostic["constraint_metrics"]["source"], "protected_lambda_gauge")
        self.assertTrue(diagnostic["k_lambda_source_available"])
        self.assertTrue(diagnostic["k_lambda_scope_template_available"])
        self.assertTrue(diagnostic["protected_fit_attempted"])
        self.assertTrue(diagnostic["protected_fit_succeeded"])
        self.assertEqual(diagnostic["fit_metrics"]["n_free_spectrum_parameters"], 3)
        self.assertTrue(diagnostic["lambda_preservation"]["gate_passed"])
        self.assertIsNotNone(result["H_pi_delta_lambda_gauge"])
        for component_name in ("k_lambda", "k_sigma0", "pi_delta"):
            self.assertTrue(diagnostic["no_double_scaling"][component_name]["passed"])

        for index in range(1, self.pi_delta.GetNbinsX() + 1):
            self.assertAlmostEqual(
                result["H_pi_delta_lambda_gauge"].GetBinContent(index),
                diagnostic["lambda_gauge"]["amplitude"]
                * self.k_lambda.GetBinContent(index),
                places=10,
            )
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
            full_expected = (
                1.10 * self.pi_n.GetBinContent(index)
                + 0.60 * self.pi_sidis.GetBinContent(index)
                + 2.00 * self.k_lambda.GetBinContent(index)
                + 0.80 * self.k_sigma.GetBinContent(index)
                + 0.25 * self.pi_delta.GetBinContent(index)
            )
            self.assertAlmostEqual(result["fit_hist"].GetBinContent(index), full_expected, places=7)
            self.assertAlmostEqual(
                result["full_fit_residual_hist"].GetBinContent(index),
                self._target().GetBinContent(index) - full_expected,
                places=7,
            )
        self.assertAlmostEqual(
            result["diagnostics"]["refined_amplitudes"]["pi_delta"], result["A_delta"], places=12
        )

    def test_inner_fitter_owns_local_sigma_reference_and_full_fit_contract(self):
        """Regression for the former outer-scope K-Sigma0 reference leak."""
        original = copy.deepcopy(bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"])
        try:
            local = bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"]
            local["residual_component_shifts"] = {"enabled": False}
            local["joint_refinement_enabled"] = False
            local["pi_delta_signal_protected_fit"] = copy.deepcopy(self.config["pi_delta_signal_protected_fit"])
            result = self.fits.fit_kaon_nosub_with_simc_pion_shapes(
                self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
                self.k_lambda, self.k_sigma, self.inp, phi_setting="Center", context="inner_scope",
            )
        finally:
            bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"] = original

        protected = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertIsNotNone(result["H_k_sigma0_simc_reference"])
        self.assertEqual(protected["status"], "success")
        self.assertTrue(protected["closure"]["full_five_component_model"]["passed"])
        self.assertIsNotNone(result["full_fit_residual_hist"])

    def test_shared_builder_publishes_fitter_sigma_reference_with_provenance(self):
        """The public K-Sigma0 artifact must originate in the inner fitter."""
        original = copy.deepcopy(bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"])
        try:
            local = bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"]
            local["residual_component_shifts"] = {"enabled": False}
            local["joint_refinement_enabled"] = False
            local["pi_delta_signal_protected_fit"] = copy.deepcopy(self.config["pi_delta_signal_protected_fit"])
            inp = dict(self.inp, particle_subtraction_mode="simc_shape_components")
            result = self.fits.build_particle_subtraction_component_result(
                self._target(),
                self._target(),
                {"pi_n": self.pi_n, "pi_delta": self.pi_delta, "pi_sidis": self.pi_sidis},
                inp,
                "setting-wide",
                kaon_signal_shape=self.k_lambda,
                kaon_sigma0_shape=self.k_sigma,
                phi_setting="Center",
                context="builder_scope",
            )
        finally:
            bgcfg.PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG["kaon_nosub"] = original

        provenance = result["diagnostics"]["kaon"]["k_sigma0_template_provenance"]
        self.assertIsNotNone(result["H_simc_shape_k_sigma0"])
        self.assertTrue(provenance["simc_input_available"])
        self.assertTrue(provenance["aligned_template_available"])
        self.assertTrue(provenance["fitter_reference_available"])
        self.assertEqual(provenance["published_source"], "inner_fitter_reference")
        self.assertTrue(provenance["published_template_available"])

        weights = self.subtraction.build_simc_shape_pion_control_weights(
            result, model_variant="staged"
        )
        protected_amplitudes = (
            result["diagnostics"]["kaon"]["pi_delta_signal_protected_fit"]
            ["protected_applied_amplitudes"]
        )
        for index in range(1, weights["H_kaon_pion_model"].GetNbinsX() + 1):
            expected = sum(
                float(protected_amplitudes[name])
                * float(result["H_simc_shape_{}".format(name)].GetBinContent(index))
                for name in ("pi_n", "pi_delta", "pi_sidis")
            )
            self.assertAlmostEqual(
                weights["H_kaon_pion_model"].GetBinContent(index), expected, places=9
            )
        self.assertTrue(
            weights["diagnostics"]["signal_templates_excluded_from_subtraction_weight"]
        )

    def test_protected_mode_off_leaves_the_legacy_payload_untouched(self):
        disabled_config = copy.deepcopy(self.config)
        disabled_config["pi_delta_signal_protected_fit"]["enabled"] = False
        payload = self._legacy_payload()
        result = self.fits._apply_signal_protected_pi_delta_fit(
            payload, self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            self.k_lambda, self.k_sigma, disabled_config, 0.70, 1.30, [],
            self.inp, "Center", 0.0, context="disabled",
        )
        self.assertIs(result, payload)
        self.assertEqual(result["A_delta"], 4.00)
        self.assertEqual(
            result["diagnostics"]["pi_delta_signal_protected_fit"]["fit_variant"],
            "disabled",
        )

    def test_missing_sigma_uses_lambda_only_protected_fallback(self):
        # The data can contain Sigma strength while its required external
        # template is absent.  Lambda and pi-delta still form the explicit
        # protected fallback; no staged pi-delta value is revived.
        result = self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(), self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            self.k_lambda, None, self.config, 0.70, 1.30, [], self.inp, "Center", 0.0,
            context="missing",
            kaon_sigma0_source_diagnostics={
                "fallback_reason": "no_source_configured",
                "source_identity": {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "phi_setting": "Center"},
            },
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertEqual(diagnostic["fit_variant"], "lambda_only_protected_fallback")
        self.assertTrue(diagnostic["fallback_attempted"])
        self.assertTrue(diagnostic["fallback_used"])
        self.assertEqual(diagnostic["fallback_reason"], "no_source_configured")
        self.assertEqual(diagnostic["sigma0_source_availability"]["status"], "unavailable")
        self.assertEqual(diagnostic["sigma0_scope_template_availability"]["reason"], "k_sigma0_scope_template_missing")
        self.assertIsNone(result["S_sigma0"])
        self.assertIsNone(result["H_pi_delta_protected_k_sigma0"])
        self.assertGreater(result["A_delta"], 0.0)
        self.assertAlmostEqual(result["A_n"], 1.10, places=12)
        self.assertAlmostEqual(result["A_sidis"], 0.60, places=12)
        self.assertEqual(diagnostic["legacy_staged_A_delta"], 4.00)

    def test_scope_local_sigma_failure_preserves_valid_source_provenance(self):
        empty_sigma = self.k_sigma.Clone("empty_sigma")
        empty_sigma.Reset()
        result = self._apply(
            self._target(k_sigma=empty_sigma),
            sigma=empty_sigma,
            source_diagnostics={
                "normalized": True,
                "configured": True,
                "source_identity": {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "phi_setting": "Center"},
            },
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["fit_variant"], "lambda_only_protected_fallback")
        self.assertTrue(diagnostic["fallback_used"])
        self.assertEqual(diagnostic["sigma0_source_availability"]["status"], "available")
        self.assertEqual(
            diagnostic["sigma0_scope_template_availability"]["reason"],
            "k_sigma0_scope_template_non_positive",
        )
        self.assertEqual(diagnostic["fallback_reason"], "k_sigma0_scope_template_non_positive")

    def test_scope_local_sigma_binning_mismatch_uses_fallback(self):
        mismatched_sigma = ROOT.TH1D("mismatched_sigma", "mismatched_sigma", 60, 0.70, 1.30)
        for index in range(1, mismatched_sigma.GetNbinsX() + 1):
            mismatched_sigma.SetBinContent(index, 1.0)
        mismatched_sigma.SetDirectory(0)
        result = self._apply(self._target(k_sigma=self.k_sigma), sigma=mismatched_sigma)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["fit_variant"], "lambda_only_protected_fallback")
        self.assertEqual(
            diagnostic["sigma0_scope_template_availability"]["reason"],
            "k_sigma0_scope_template_binning_mismatch",
        )

    def test_disabled_lambda_only_policy_retains_zero_pi_delta_safety(self):
        config = copy.deepcopy(self.config)
        config["pi_delta_signal_protected_fit"]["allow_lambda_only_fallback"] = False
        result = self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(), self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            self.k_lambda, None, config, 0.70, 1.30, [], self.inp, "Center", 0.0,
            context="fallback_disabled",
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["fit_variant"], "zero_pi_delta_failure")
        self.assertFalse(diagnostic["fallback_attempted"])
        self.assertFalse(diagnostic["fallback_used"])
        self.assertEqual(result["A_delta"], 0.0)

    def test_lambda_or_pi_delta_remains_mandatory(self):
        missing_payload = self._legacy_payload()
        missing_payload.pop("H_k_lambda_simc_reference")
        missing_lambda = self.fits._apply_signal_protected_pi_delta_fit(
            missing_payload, self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            None, self.k_sigma, self.config, 0.70, 1.30, [], self.inp, "Center", 0.0,
            context="missing_lambda",
        )
        diagnostic = missing_lambda["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["fit_variant"], "zero_pi_delta_failure")
        self.assertFalse(diagnostic["k_lambda_source_available"])
        self.assertFalse(diagnostic["fallback_attempted"])
        self.assertEqual(missing_lambda["A_delta"], 0.0)

    def test_rank_deficiency_is_not_reported_as_high_correlation_warning_only(self):
        result = self._apply(self._target(k_sigma=self.pi_delta), sigma=self.pi_delta)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "rank_deficient")
        self.assertEqual(diagnostic["selected_fit_variant"], "lambda_sigma0_protected")
        self.assertEqual(diagnostic["fit_variant"], "zero_pi_delta_failure")
        self.assertFalse(diagnostic["fallback_attempted"])
        self.assertEqual(result["A_delta"], 0.0)

    def test_lambda_only_rank_failure_remains_zeroed(self):
        payload = self._legacy_payload()
        payload["H_k_lambda_simc_reference"] = self.pi_delta.Clone("rank_lambda_reference")
        result = self.fits._apply_signal_protected_pi_delta_fit(
            payload, self._target(), self.pi_n, self.pi_delta, self.pi_sidis,
            self.pi_delta, None, self.config, 0.70, 1.30, [], self.inp, "Center", 0.0,
            context="fallback_rank",
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["selected_fit_variant"], "lambda_only_protected_fallback")
        self.assertEqual(diagnostic["fit_variant"], "zero_pi_delta_failure")
        self.assertTrue(diagnostic["fallback_attempted"])
        self.assertFalse(diagnostic["fallback_used"])
        self.assertEqual(result["A_delta"], 0.0)

    def test_strong_sigma_delta_overlap_remains_a_solvable_protected_fit(self):
        overlapping_sigma = self._shape("k_sigma_overlap", 1.210, 0.045)
        result = self._apply(self._target(k_sigma=overlapping_sigma), sigma=overlapping_sigma)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertAlmostEqual(result["A_delta"], 0.25, places=6)
        self.assertAlmostEqual(result["S_sigma0"], 0.80, places=6)

    def test_fixed_gauge_excludes_lambda_from_spectrum_ndf_and_rank(self):
        config = copy.deepcopy(self.config)
        config["pi_delta_signal_protected_fit"]["lambda_gauge_constraint_mode"] = "fixed"
        result = self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(), self._target(), self.pi_n, self.pi_delta,
            self.pi_sidis, self.k_lambda, self.k_sigma, config, 0.70, 1.30,
            [], self.inp, "Center", 0.0, context="fixed_gauge",
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertEqual(diagnostic["constraint_metrics"]["mode"], "fixed")
        self.assertEqual(diagnostic["fit_metrics"]["n_free_spectrum_parameters"], 2)
        self.assertEqual(diagnostic["matrix_diagnostics"]["weighted_design_effective_rank"], 2)
        self.assertAlmostEqual(
            diagnostic["proposed_amplitudes"]["k_lambda_signal"],
            diagnostic["lambda_gauge"]["amplitude"],
            places=12,
        )

    def test_poor_gauge_quality_uses_an_inflated_lambda_constraint(self):
        config = copy.deepcopy(self.config)
        config["pi_delta_signal_protected_fit"]["lambda_gauge_maximum_chi2_ndf"] = 1.0e-12
        result = self.fits._apply_signal_protected_pi_delta_fit(
            self._legacy_payload(), self._target(), self.pi_n, self.pi_delta,
            self.pi_sidis, self.k_lambda, self.k_sigma, config, 0.70, 1.30,
            [], self.inp, "Center", 0.0, context="poor_gauge",
        )
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertFalse(diagnostic["lambda_gauge_quality_passed"])
        self.assertEqual(diagnostic["constraint_metrics"]["mode"], "gaussian_inflated")
        self.assertEqual(
            diagnostic["lambda_constraint"]["source"],
            "protected_lambda_gauge_poor_quality",
        )
        self.assertEqual(diagnostic["fit_metrics"]["n_free_spectrum_parameters"], 3)
        self.assertGreater(result["A_delta"], 0.0)
        self.assertGreater(result["S_lambda"], 0.0)

    def test_failed_gauge_uses_canonical_pre_delta_lambda_constraint(self):
        failed_gauge = {
            "status": "insufficient_lambda_gauge_bins",
            "solver_success": False,
            "quality_passed": False,
            "window": (1.105, 1.125),
            "amplitude": None,
            "amplitude_sigma": None,
            "effective_sigma": None,
            "failure_reason": "insufficient_lambda_gauge_bins",
        }
        with mock.patch.object(self.fits, "_fit_lambda_gauge", return_value=failed_gauge):
            result = self._apply(self._target())
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "success")
        self.assertEqual(diagnostic["constraint_metrics"]["mode"], "gaussian_inflated")
        self.assertEqual(
            diagnostic["lambda_constraint"]["source"],
            "canonical_k_lambda_simc_pre_delta",
        )
        self.assertFalse(diagnostic["lambda_gauge_solver_success"])
        self.assertGreater(diagnostic["lambda_constraint"]["amplitude"], 0.0)
        self.assertGreater(result["S_lambda"], 0.0)

    def test_pre_delta_lambda_deficit_is_explicit_and_never_attempts_delta(self):
        target = self._target()
        target.Add(self.k_lambda, -1.5)
        forced_gauge = {
            "status": "success",
            "solver_success": True,
            "quality_passed": True,
            "window": (1.105, 1.125),
            "amplitude": 2.0,
            "amplitude_sigma": 0.10,
            "effective_sigma": 0.10,
            "chi2": 0.0,
            "ndf": 3,
            "chi2_ndf": 0.0,
            "p_value": 1.0,
            "fit_bins": 4,
            "template_integral_window": 1.0,
            "data_integral_window": 0.5,
            "gauge_predicted_yield_window": 2.0,
            "failure_reason": None,
        }
        with mock.patch.object(self.fits, "_fit_lambda_gauge", return_value=forced_gauge):
            result = self._apply(target)
        diagnostic = result["diagnostics"]["pi_delta_signal_protected_fit"]
        self.assertEqual(diagnostic["status"], "lambda_pre_delta_deficit")
        self.assertEqual(diagnostic["lambda_preservation"]["status"], "lambda_pre_delta_deficit")
        self.assertIsNone(diagnostic["proposed_A_delta"])
        self.assertEqual(result["A_delta"], 0.0)


if __name__ == "__main__":
    unittest.main()
