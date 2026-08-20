"""Non-ROOT contracts for diagnostic-only pion pages and Lambda comparison."""

from __future__ import annotations

import ast
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


def _load_fits_without_root():
    try:
        import ROOT  # noqa: F401
        import pion_component_fits as fits
        return fits
    except ImportError:
        root_stub = types.ModuleType("ROOT")

        class _RootImportStub:
            def __init__(self, *_args, **_kwargs):
                pass

        for name in ("TFile", "TNtuple", "TText", "TGraph", "TGraphErrors", "TCanvas", "TF1", "TFitResultPtr"):
            setattr(root_stub, name, _RootImportStub)
        for name in ("kBlack", "kBlue", "kGray", "kGreen", "kMagenta", "kOrange", "kRed", "kViolet", "kAzure", "kCyan"):
            setattr(root_stub, name, 1)
        utility_stub = types.ModuleType("utility")
        utility_stub.normalize_hist_to_unit_area = lambda histogram, *_args, **_kwargs: histogram
        with mock.patch.dict(sys.modules, {"ROOT": root_stub, "utility": utility_stub}):
            sys.modules.pop("pion_component_fits", None)
            import pion_component_fits as fits
        return fits


fits = _load_fits_without_root()
RAND_SUB_PATH = REPO_ROOT / "src" / "cuts" / "rand_sub.py"
AVE_PER_BIN_PATH = REPO_ROOT / "src" / "binning" / "ave_per_bin.py"


class PionDiagnosticPlotContractTests(unittest.TestCase):
    def test_setting_wide_reconstruction_is_explicitly_cloned_from_active_targets(self):
        source = RAND_SUB_PATH.read_text(encoding="utf-8")
        tree = ast.parse(source)
        function_sources = {
            node.name: ast.get_source_segment(source, node)
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef)
        }
        self.assertIn("active_component_targets", function_sources["rand_sub"])
        self.assertIn(
            "_clone_component_targets_for_setting_wide_diagnostic(\n                            active_component_targets",
            function_sources["rand_sub"],
        )
        diagnostic_builder = function_sources["_clone_component_targets_for_setting_wide_diagnostic"]
        self.assertIn('scope="setting_wide_diagnostic"', diagnostic_builder)
        self.assertIn('role="application_target_{}".format(key)', diagnostic_builder)
        application_builder = function_sources["_apply_component_pion_subtraction_setting"]
        self.assertIn('"input_selection": str(input_selection)', application_builder)
        self.assertIn('"source_target_state": str(source_target_state)', application_builder)
        self.assertIn('"production_applied": not bool(diagnostic_only)', application_builder)

    def test_setting_wide_emission_is_gated_without_suppressing_per_t_pages(self):
        rand_source = RAND_SUB_PATH.read_text(encoding="utf-8")
        ave_source = AVE_PER_BIN_PATH.read_text(encoding="utf-8")
        self.assertIn("setting_wide_pages_enabled", rand_source)
        self.assertIn('emit_setting_wide_pion_diagnostic", True', rand_source)
        self.assertIn("render_setting_wide_template_diagnostics", ave_source)
        self.assertIn('page_id_prefix="pion.t_bin.{}".format(j)', ave_source)
        self.assertIn("print_particle_subtraction_kaon_lambda_comparison_page(", ave_source)

    def test_parent_builder_owns_the_manifest_not_process_hist_data(self):
        source = AVE_PER_BIN_PATH.read_text(encoding="utf-8")
        tree = ast.parse(source)
        function_sources = {
            node.name: ast.get_source_segment(source, node)
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef)
        }
        process_source = function_sources["process_hist_data"]
        self.assertIn("component_page_manifest=None", process_source)
        self.assertNotIn("hist.setdefault", process_source)
        self.assertIn("component_page_manifest=hist.setdefault", function_sources["build_t_bin_pion_parents"])

    def test_page_manifest_rejects_duplicates_and_preserves_scope_authority(self):
        manifest = []
        fits.record_particle_subtraction_page(
            manifest,
            "pion.setting_wide.lambda_comparison",
            scope="setting-wide",
            authoritative=False,
        )
        fits.record_particle_subtraction_page(
            manifest,
            "pion.t_bin.0.lambda_comparison",
            scope="t_bin1",
            authoritative=True,
        )
        self.assertEqual([entry["page_id"] for entry in manifest], [
            "pion.setting_wide.lambda_comparison",
            "pion.t_bin.0.lambda_comparison",
        ])
        self.assertFalse(manifest[0]["authoritative"])
        self.assertTrue(manifest[1]["authoritative"])
        with self.assertRaisesRegex(RuntimeError, "duplicate_particle_subtraction_page_id"):
            fits.record_particle_subtraction_page(
                manifest,
                "pion.t_bin.0.lambda_comparison",
                scope="t_bin1",
                authoritative=True,
            )

    def test_lambda_page_uses_scope_after_spectrum_and_records_once(self):
        after_hist = object()
        before_hist = object()
        canonical_hist = object()
        fit_result = {
            "analysis_scope": "t_bin1",
            "H_kaon_nosub_input": object(),
            "diagnostics": {},
        }
        payload = {
            "analysis_scope": "t_bin1",
            "accepted": True,
            "H_MM_nosub_before_pion_subtraction": before_hist,
            "H_MM_nosub_after_pion_subtraction": after_hist,
        }
        manifest = []
        with mock.patch.object(
            fits,
            "_resolve_kaon_lambda_reference_for_plot",
            return_value=(canonical_hist, 1.0, "canonical", "historical"),
        ) as resolver, mock.patch.object(
            fits,
            "_print_component_overlay_page",
            return_value=True,
        ):
            emitted = fits.print_particle_subtraction_kaon_lambda_comparison_page(
                "ignored.pdf",
                fit_result,
                payload,
                page_manifest=manifest,
                page_id_prefix="pion.t_bin.0",
                authoritative=True,
            )
        self.assertTrue(emitted)
        self.assertIs(resolver.call_args.args[1], after_hist)
        self.assertEqual(
            [entry["page_id"] for entry in manifest],
            ["pion.t_bin.0.lambda_comparison"],
        )

    def test_rejected_application_is_status_only_and_can_still_render_lambda(self):
        rejected = {"accepted": False, "analysis_scope": "setting-wide"}
        fit_result = {"analysis_scope": "setting-wide"}
        with mock.patch.object(
            fits,
            "_print_component_application_status_page",
            return_value=True,
        ) as status_page, mock.patch.object(
            fits,
            "print_particle_subtraction_kaon_lambda_comparison_page",
            return_value=True,
        ) as lambda_page:
            fits.print_particle_subtraction_component_application_pages(
                "ignored.pdf",
                rejected,
                component_fit_result=fit_result,
                page_manifest=[],
                page_id_prefix="pion.setting_wide",
            )
        status_page.assert_called_once()
        lambda_page.assert_called_once()


if __name__ == "__main__":
    unittest.main()
