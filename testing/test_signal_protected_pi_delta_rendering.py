"""Non-ROOT PDF-routing contracts for signal-protected kaon pi-delta fits."""

from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)


def _import_fits_without_root():
    try:
        import ROOT  # noqa: F401
        import pion_component_fits as fits
        return fits
    except ImportError:
        root_stub = types.ModuleType("ROOT")
        class _RootImportStub:
            def __init__(self, *_args, **_kwargs):
                pass

        for name in (
            "TFile", "TNtuple", "TText", "TGraph", "TGraphErrors",
            "TCanvas", "TF1", "TFitResultPtr",
        ):
            setattr(root_stub, name, _RootImportStub)
        for name in (
            "kBlack", "kBlue", "kGray", "kGreen", "kMagenta", "kOrange",
            "kRed", "kViolet", "kAzure", "kCyan",
        ):
            setattr(root_stub, name, 1)
        utility_stub = types.ModuleType("utility")
        utility_stub.normalize_hist_to_unit_area = lambda histogram, *_args, **_kwargs: histogram
        with mock.patch.dict(sys.modules, {"ROOT": root_stub, "utility": utility_stub}):
            sys.modules.pop("pion_component_fits", None)
            import pion_component_fits as fits
        return fits


fits = _import_fits_without_root()


class SignalProtectedPiDeltaRenderingTests(unittest.TestCase):
    def setUp(self):
        self.objects = {key: object() for key in (
            "pion_input", "kaon_input", "pi_n", "pi_sidis", "pi_delta",
            "protected_input", "protected_gauge", "protected_lambda", "protected_sigma",
            "protected_delta", "protected_total", "protected_after",
        )}

    def _protected_diagnostic(self, status="success"):
        return {
            "enabled": True,
            "status": status,
            "fit_variant": "lambda_sigma0_protected",
            "selected_fit_variant": "lambda_sigma0_protected",
            "failure_reason": "missing K-Sigma0" if status != "success" else None,
            "failure_policy": "zero_pi_delta",
            "applied_A_delta": 0.25 if status == "success" else 0.0,
            "physics_acceptance_passed": status == "success",
            "solver_success": status == "success",
            "fit_quality_passed": status == "success",
            "lambda_gauge_solver_success": status == "success",
            "lambda_gauge_quality_passed": status == "success",
            "lambda_gauge_status": "success" if status == "success" else "poor_lambda_gauge_quality",
            "legacy_staged_A_delta": 3.0,
            "signal_amplitudes": {"k_lambda_signal": 2.0, "k_sigma0_signal": 0.8},
            "proposed_amplitudes": {"k_lambda_signal": 2.0, "k_sigma0_signal": 0.8, "pi_delta": 0.25},
            "fit_metrics": {"chi2_ndf": 1.2, "fit_p_value": 0.4, "n_fit_bins": 25, "n_free_spectrum_parameters": 3},
            "constraint_metrics": {"mode": "gaussian", "prior_chi2": 0.1, "total_chi2": 2.5, "total_ndf": 23},
            "matrix_diagnostics": {
                "weighted_design_effective_rank": 3,
                "weighted_design_condition_number": 12.0,
            },
            "template_availability": {
                "k_lambda_signal": status == "success",
                "k_sigma0_signal": status == "success",
                "pi_delta": True,
            },
            "lambda_reference_integrity": {"shape_identical": True},
            "lambda_gauge": {
                "status": "success" if status == "success" else "poor_lambda_gauge_quality",
                "solver_success": status == "success",
                "quality_passed": status == "success",
                "window": [1.105, 1.125],
                "window_source": "proton_cleaning.validation_windows.lambda_peak",
                "amplitude": 2.0,
                "amplitude_sigma": 0.10,
                "effective_sigma": 0.10,
                "chi2_ndf": 1.0,
                "p_value": 0.5,
                "fit_bins": 4,
                "data_integral_window": 2.0,
                "gauge_predicted_yield_window": 2.0,
            },
            "lambda_preservation": {
                "status": "bounded",
                "gate_passed": status == "success",
                "gate_reason": None if status == "success" else "fit_not_accepted",
                "lambda_gauge_predicted_yield": 2.0,
                "lambda_pre_delta_yield": 2.1,
                "lambda_pi_delta_removed_yield": 0.05,
                "lambda_after_delta_yield": 2.05,
                "lambda_after_over_gauge": 1.025,
                "lambda_removed_fraction_of_gauge": 0.025,
                "minimum_required_retention": 1.8,
                "a_delta_max": 0.5,
            },
            "early_amplitudes_frozen_integrity": {"unchanged": True},
            "signal_preservation": {
                "k_lambda_signal": {"pi_delta_removed_fraction": 0.01},
                "k_sigma0_signal": {"pi_delta_removed_fraction": 0.02},
            },
            "closure": {
                "protected_three_component_model": {"passed": True},
                "delta_only_physics_output": {"passed": True},
            },
        }

    def _fit_payload(self, enabled=True, diagnostic=True, status="success"):
        payload = {
            "analysis_scope": "setting-wide",
            "resolved_subtraction_config": {
                "kaon_nosub": {
                    "pi_delta_signal_protected_fit": {
                        "enabled": enabled,
                        "failure_policy": "zero_pi_delta",
                    },
                },
            },
            "diagnostics": {"pion": {}, "kaon": {}},
            "H_pion_control_input": self.objects["pion_input"],
            "H_kaon_nosub_input": self.objects["kaon_input"],
            "H_pion_fit_pi_n_scaled": self.objects["pi_n"],
            "H_pion_fit_pi_sidis_scaled": self.objects["pi_sidis"],
            "H_pion_fit_pi_delta_scaled": self.objects["pi_delta"],
            "H_pion_fit_total": object(),
            "H_pion_fit_total_refined": object(),
            "H_kaon_fit_pi_n_scaled": self.objects["pi_n"],
            "H_kaon_fit_pi_sidis_scaled": self.objects["pi_sidis"],
            "H_kaon_fit_pi_delta_scaled": self.objects["pi_delta"],
            "H_kaon_pion_bg_fit_total": object(),
            "H_kaon_fit_total": object(),
            "H_kaon_fit_total_refined": object(),
            "H_kaon_pion_bg_fit_total_refined": object(),
            "H_kaon_pion_bg_fit_total_refined_pre_postrefine": object(),
            "H_pion_fit_step_overlays": [
                {"component_name": "pi_n"},
                {"component_name": "pi_sidis"},
                {"component_name": "pi_delta"},
            ],
            "H_kaon_fit_step_overlays": [
                {"component_name": "pi_n"},
                {"component_name": "pi_sidis"},
                {"component_name": "pi_delta"},
            ],
        }
        if diagnostic:
            payload["diagnostics"]["kaon"]["pi_delta_signal_protected_fit"] = (
                self._protected_diagnostic(status)
            )
        if enabled and diagnostic and status == "success":
            payload.update({
                "H_pi_delta_protected_fit_input": self.objects["protected_input"],
                "H_pi_delta_lambda_gauge": self.objects["protected_gauge"],
                "H_pi_delta_protected_k_lambda": self.objects["protected_lambda"],
                "H_pi_delta_protected_k_sigma0": self.objects["protected_sigma"],
                "H_pi_delta_protected_pi_delta": self.objects["protected_delta"],
                "H_pi_delta_protected_fit_total": self.objects["protected_total"],
                "H_pi_delta_protected_after_subtraction": self.objects["protected_after"],
            })
        return payload

    def _lambda_only_fallback_payload(self):
        payload = self._fit_payload()
        protected = payload["diagnostics"]["kaon"]["pi_delta_signal_protected_fit"]
        protected.update(
            {
                "fit_variant": "lambda_only_protected_fallback",
                "selected_fit_variant": "lambda_only_protected_fallback",
                "fallback_attempted": True,
                "fallback_used": True,
                "fallback_reason": "no_source_configured",
                "sigma0_availability_reason": "no_source_configured",
                "sigma0_fitted": False,
                "signal_amplitudes": {"k_lambda_signal": 2.0, "k_sigma0_signal": None},
                "template_availability": {
                    "k_lambda_signal": True,
                    "k_sigma0_signal": False,
                    "pi_delta": True,
                },
                "closure": {
                    "protected_two_component_model": {"passed": True},
                    "delta_only_physics_output": {"passed": True},
                },
            }
        )
        payload.pop("H_pi_delta_protected_k_sigma0")
        return payload

    def _fit_page_recorders(self):
        overlay_pages = []
        text_pages = []
        step_calls = []
        amplitude_calls = []
        joint_pages = []
        comparison_pages = []

        def overlay(_pdf, base, _label, title, overlays, lines, **_kwargs):
            if base is not None and overlays:
                overlay_pages.append((title, base, list(overlays), list(lines)))

        def protected_overlay(_pdf, base, _label, title, overlays, lines, **_kwargs):
            if base is not None and overlays:
                overlay_pages.append((title, base, list(overlays), list(lines)))

        def text(_pdf, title, header, body):
            text_pages.append((title, list(header), list(body)))

        def joint(_pdf, data, *_args, **kwargs):
            if data is not None:
                joint_pages.append(kwargs.get("title") or (_args[3] if len(_args) > 3 else ""))

        def comparison(_pdf, data, *_args, **kwargs):
            if data is not None:
                comparison_pages.append(kwargs.get("title") or (_args[4] if len(_args) > 4 else ""))

        patches = {
            "_print_component_overlay_page": overlay,
            "_print_protected_overlay_page": protected_overlay,
            "_print_component_text_page": text,
            "_print_joint_refinement_overlay_page": joint,
            "_print_kaon_pion_bg_comparison_page": comparison,
            "_print_residual_shift_template_pages": lambda *_args, **_kwargs: None,
            "_print_residual_shift_scan_pages": lambda *_args, **_kwargs: None,
            "_print_component_step_pages": lambda _pdf, _target, steps, _prefix, label: step_calls.append((label, list(steps or []))),
            "_print_component_amplitude_pages": lambda _pdf, _target, steps, _prefix, label: amplitude_calls.append((label, list(steps or []))),
            "_resolve_kaon_lambda_reference_for_plot": lambda *_args, **_kwargs: (None, None, "test"),
        }
        return patches, overlay_pages, text_pages, step_calls, amplitude_calls, joint_pages, comparison_pages

    def test_success_uses_gauge_fit_preservation_pages_and_filters_kaon_steps(self):
        payload = self._fit_payload()
        patches, overlay_pages, text_pages, step_calls, amplitude_calls, joint_pages, comparison_pages = self._fit_page_recorders()
        with mock.patch.multiple(fits, **patches):
            fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)

        titles = [page[0] for page in overlay_pages]
        protected_titles = [
            title for title in titles
            if "K-Lambda pre-pi-delta gauge" in title
            or "K-Lambda gauged" in title
            or "Lambda preservation" in title
        ]
        self.assertEqual(protected_titles, [
            "K-Lambda pre-pi-delta gauge",
            "Signal-protected final pi-delta fit - K-Lambda gauged + K-Sigma0",
            "Protected pi-delta subtraction - Lambda preservation",
        ])
        self.assertFalse(text_pages)
        self.assertFalse(any("kaon no-sub staged" in title for title in titles))
        self.assertFalse(any("kaon no-sub staged vs refined" in title for title in joint_pages))
        self.assertFalse(comparison_pages)
        first_protected = overlay_pages[titles.index(protected_titles[0])]
        self.assertIs(first_protected[1], self.objects["protected_input"])
        self.assertEqual(
            [overlay[0] for overlay in first_protected[2]],
            [self.objects["protected_gauge"]],
        )
        second_protected = overlay_pages[titles.index(protected_titles[1])]
        self.assertEqual(
            [overlay[0] for overlay in second_protected[2]],
            [self.objects["protected_lambda"], self.objects["protected_sigma"], self.objects["protected_delta"], self.objects["protected_total"]],
        )
        third_protected = overlay_pages[titles.index(protected_titles[2])]
        self.assertEqual(
            [overlay[0] for overlay in third_protected[2]],
            [self.objects["protected_gauge"], self.objects["protected_delta"], self.objects["protected_after"]],
        )
        pion_steps = next(steps for label, steps in step_calls if label == "pion-control")
        kaon_steps = next(steps for label, steps in step_calls if label == "kaon no-sub")
        self.assertEqual([entry["component_name"] for entry in pion_steps], ["pi_n", "pi_sidis", "pi_delta"])
        self.assertEqual([entry["component_name"] for entry in kaon_steps], ["pi_n", "pi_sidis"])
        self.assertIsNot(kaon_steps[0], payload["H_kaon_fit_step_overlays"][0])
        self.assertEqual(
            [entry["component_name"] for entry in payload["H_kaon_fit_step_overlays"]],
            ["pi_n", "pi_sidis", "pi_delta"],
        )
        self.assertEqual(
            [entry["component_name"] for label, entries in amplitude_calls if label == "kaon no-sub" for entry in entries],
            ["pi_n", "pi_sidis"],
        )

    def test_resolved_configuration_overrides_diagnostic_enablement(self):
        payload = self._fit_payload(enabled=False)
        state = fits._resolve_protected_pi_delta_render_state(payload)
        self.assertEqual(state["state"], "disabled")
        self.assertFalse(state["suppress_deprecated_kaon_pages"])

        payload = self._fit_payload(enabled=True, diagnostic=True)
        payload["resolved_subtraction_config"]["kaon_nosub"].pop(
            "pi_delta_signal_protected_fit"
        )
        state = fits._resolve_protected_pi_delta_render_state(payload)
        self.assertEqual(state["state"], "success")
        self.assertTrue(state["suppress_deprecated_kaon_pages"])

    def test_lambda_only_fallback_renders_without_a_sigma0_curve(self):
        payload = self._lambda_only_fallback_payload()
        state = fits._resolve_protected_pi_delta_render_state(payload)
        self.assertEqual(state["state"], "success")
        self.assertEqual(state["fit_variant"], "lambda_only_protected_fallback")

        patches, overlay_pages, text_pages, *_unused = self._fit_page_recorders()
        with mock.patch.multiple(fits, **patches):
            fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)

        protected_pages = [page for page in overlay_pages if "K-Lambda" in page[0] or "Lambda preservation" in page[0]]
        self.assertEqual(len(protected_pages), 3)
        self.assertFalse(text_pages)
        self.assertEqual(
            [overlay[0] for overlay in protected_pages[1][2]],
            [
                self.objects["protected_lambda"],
                self.objects["protected_delta"],
                self.objects["protected_total"],
            ],
        )
        self.assertTrue(
            any("A_delta proposed/applied" in line for line in protected_pages[1][3])
        )

    def test_enabled_failure_or_missing_payload_gets_status_not_legacy_kaon_pages(self):
        for payload, expected in (
            (self._fit_payload(diagnostic=True, status="missing_required_template"), "PROTECTED FIT UNAVAILABLE"),
            (self._fit_payload(diagnostic=False), "PROTECTED RENDER-CONTRACT ERROR"),
        ):
            with self.subTest(expected=expected):
                patches, overlay_pages, text_pages, *_unused = self._fit_page_recorders()
                with mock.patch.multiple(fits, **patches):
                    fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)
                self.assertEqual(len(text_pages), 1)
                self.assertIn(expected, text_pages[0][1][1])
                self.assertFalse(any("kaon no-sub staged" in page[0] for page in overlay_pages))
                self.assertFalse(any("signal-protected final #pi#Delta fit" == page[0] for page in overlay_pages))

    def test_protected_final_comparison_reuses_gauge_without_postsubtraction_scaling(self):
        payload = self._fit_payload()
        with mock.patch.object(fits, "_clone_hist", side_effect=lambda hist, _name, **_kwargs: hist), \
             mock.patch.object(fits, "_build_scaled_reference_hist_with_fallback", side_effect=AssertionError):
            histogram, scale, source = fits._resolve_kaon_lambda_reference_for_plot(
                payload, object(), (1.105, 1.125), "setting-wide", "gauge_check"
            )
        self.assertIs(histogram, self.objects["protected_gauge"])
        self.assertEqual(scale, 2.0)
        self.assertEqual(source, "authoritative pre-delta K-Lambda gauge")

        payload.pop("H_pi_delta_lambda_gauge")
        with mock.patch.object(fits, "_build_scaled_reference_hist_with_fallback", side_effect=AssertionError):
            histogram, scale, source = fits._resolve_kaon_lambda_reference_for_plot(
                payload, object(), (1.105, 1.125), "setting-wide", "gauge_missing"
            )
        self.assertIsNone(histogram)
        self.assertEqual(scale, 2.0)
        self.assertIn("unavailable", source)

    def test_unavailable_status_includes_sigma0_source_provenance_when_present(self):
        payload = self._fit_payload(diagnostic=True, status="missing_required_template")
        payload["kaon_sigma0_source_diagnostics"] = {
            "resolution_source": "configured_path_does_not_exist",
            "requested_environment_variable": "LT_BG_SIGMA0_LEFT_LOW_ROOT",
            "requested_root": "/external/sigma0/left.root",
            "resolved_root": "/external/sigma0/left.root",
            "path_exists": False,
            "tree_exists": None,
            "tree_entries": None,
            "missing_required_branches": [],
            "fallback_reason": "configured_path_does_not_exist",
            "source_identity": {
                "Q2": "4p4",
                "W": "2p74",
                "EPSSET": "low",
                "phi_setting": "Left",
            },
        }
        patches, _overlay_pages, text_pages, *_unused = self._fit_page_recorders()
        with mock.patch.multiple(fits, **patches):
            fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)

        body = text_pages[0][2]
        self.assertTrue(any("K-Sigma0 source: configured_path_does_not_exist Q4p4 W2p74 low Left" in line for line in body))
        self.assertTrue(any("K-Sigma0 environment: LT_BG_SIGMA0_LEFT_LOW_ROOT" in line for line in body))
        self.assertTrue(any("K-Sigma0 requested: /external/sigma0/left.root" in line for line in body))

    def test_failed_lambda_only_attempt_reports_attempt_without_use(self):
        payload = self._fit_payload(diagnostic=True, status="rank_deficient")
        protected = payload["diagnostics"]["kaon"]["pi_delta_signal_protected_fit"]
        protected.update(
            {
                "fit_variant": "zero_pi_delta_failure",
                "selected_fit_variant": "lambda_only_protected_fallback",
                "fallback_attempted": True,
                "fallback_used": False,
                "fallback_reason": "k_sigma0_scope_template_non_positive",
                "sigma0_source_availability": {"status": "available", "reason": None},
                "sigma0_scope_template_availability": {
                    "status": "unavailable",
                    "reason": "k_sigma0_scope_template_non_positive",
                },
            }
        )
        patches, _overlay_pages, text_pages, *_unused = self._fit_page_recorders()
        with mock.patch.multiple(fits, **patches):
            fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)

        body = text_pages[0][2]
        self.assertTrue(any("attempted=True used=False" in line for line in body))
        self.assertTrue(any("scope-template availability: unavailable" in line for line in body))

    def test_disabled_mode_retains_legacy_kaon_sequence(self):
        payload = self._fit_payload(enabled=False, diagnostic=False)
        patches, overlay_pages, text_pages, step_calls, _amplitude_calls, joint_pages, comparison_pages = self._fit_page_recorders()
        with mock.patch.multiple(fits, **patches):
            fits.print_particle_subtraction_component_fit_pages("ignored.pdf", payload)

        self.assertFalse(text_pages)
        self.assertTrue(any("kaon no-sub staged SIMC" in page[0] for page in overlay_pages))
        self.assertTrue(any("kaon no-sub staged vs refined" in title for title in joint_pages))
        self.assertTrue(comparison_pages)
        kaon_steps = next(steps for label, steps in step_calls if label == "kaon no-sub")
        self.assertEqual([entry["component_name"] for entry in kaon_steps], ["pi_n", "pi_sidis", "pi_delta"])

    def test_application_pages_label_protected_model_and_hide_staged_comparisons(self):
        fit_result = self._fit_payload()
        payload = {
            "accepted": True,
            "fit_mode": "staged_plus_joint",
            "H_pion_weight_vs_MM": object(),
            "H_pion_weight_vs_MM_stage": object(),
            "H_kaon_pion_model": object(),
            "H_kaon_pion_model_stage": object(),
            "H_weighted_pion_control_model": object(),
            "H_pion_subtraction_template_MM_nosub": object(),
            "H_MM_nosub_before_pion_subtraction": object(),
            "H_MM_nosub_after_pion_subtraction": object(),
            "diagnostics": {"pion": {}, "kaon": {}},
        }
        self.assertFalse(any(
            key.startswith("H_pi_delta_protected_")
            or key in {"H_pi_delta_lambda_gauge", "H_k_lambda_simc_reference", "H_simc_shape_k_lambda"}
            for key in payload
        ))
        titles = []
        labels = []
        with mock.patch.object(fits, "_resolve_protected_pi_delta_render_state", wraps=fits._resolve_protected_pi_delta_render_state) as protected_state, \
             mock.patch.object(fits, "_resolve_kaon_lambda_reference_for_plot", return_value=(None, None, "test")) as lambda_reference, \
             mock.patch.object(fits, "_print_single_hist_page", side_effect=lambda _pdf, _hist, _label, title, *_args, **_kwargs: titles.append(title)), \
             mock.patch.object(fits, "_print_component_overlay_page", side_effect=lambda _pdf, _hist, label, title, *_args, **_kwargs: (labels.append(label), titles.append(title))):
            fits.print_particle_subtraction_component_application_pages(
                "ignored.pdf", payload, component_fit_result=fit_result
            )

        self.assertTrue(any("protected applied" in title for title in titles))
        self.assertIn("protected applied pion-background model", labels)
        self.assertFalse(any("staged vs refined" in title for title in titles))
        self.assertIs(protected_state.call_args.args[0], fit_result)
        self.assertIs(lambda_reference.call_args.args[0], fit_result)


if __name__ == "__main__":
    unittest.main()
