"""Focused regression checks for proton-cleaning runtime-state handling."""

from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

try:
    import ROOT  # noqa: F401
    import proton_contamination_weights as proton_cleaning
except ImportError:
    # The focused helpers below do not need ROOT numerics.  A narrow import
    # stub keeps their regression coverage available in non-PyROOT CI.
    class _RootImportStub:
        def __init__(self, *_args, **_kwargs):
            pass

    root_stub = types.ModuleType("ROOT")
    for name in (
        "TCanvas", "TLatex", "TLegend", "TLine", "TPaveText", "TH1D", "TH2D",
        "TF1", "TGraphErrors",
    ):
        setattr(root_stub, name, _RootImportStub)
    root_stub.gPad = types.SimpleNamespace()
    root_stub.gStyle = types.SimpleNamespace()
    for name in ("kBlack", "kBlue", "kGray", "kGreen", "kMagenta", "kOrange", "kRed", "kViolet"):
        setattr(root_stub, name, 1)
    with mock.patch.dict(sys.modules, {"ROOT": root_stub}):
        sys.modules.pop("proton_contamination_weights", None)
        import proton_contamination_weights as proton_cleaning


class _FakeAxis:
    def FindFixBin(self, _value):
        return 1

    def GetBinCenter(self, _index):
        return 0.0


class _FakeHistogram:
    def GetMaximum(self):
        return 10.0

    def GetEntries(self):
        return 100.0

    def GetEffectiveEntries(self):
        return 100.0

    def GetNbinsX(self):
        return 10

    def GetBinContent(self, _index):
        return 1.0

    def Integral(self, *_args):
        return 10.0

    def GetXaxis(self):
        return _FakeAxis()

    def Fit(self, _function, _options):
        return 0


class _FakeTF1:
    def __init__(self, *_args):
        self.parameters = {}

    def SetParName(self, *_args):
        pass

    def SetParameter(self, index, value):
        self.parameters[int(index)] = float(value)

    def SetParLimits(self, *_args):
        pass

    def FixParameter(self, index, value):
        self.parameters[int(index)] = float(value)

    def GetParameter(self, index):
        return self.parameters.get(int(index), 1.0)

    def GetParError(self, _index):
        return 0.25

    def GetChisquare(self):
        return 1.0

    def GetNDF(self):
        return 4


class ProtonCleaningRuntimeStatusTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "ctime_hist_range": (-2.0, 2.0),
            "slice_fit": {"minimum_entries": 30},
            "low_aero_offset": {"minimum_prompt_events": 20},
            "tof_offset_validation": dict(proton_cleaning.PROTON_CLEANING_TOF_OFFSET_VALIDATION),
        }
        self.global_shape = {
            "valid": True,
            "fit_min": -1.0,
            "fit_max": 1.0,
            "kaon_mean": 0.0,
            "proton_mean": 0.5,
            "kaon_sigma": 0.15,
            "proton_sigma": 0.15,
        }
        self.timing_constraint = {
            "valid": True,
            "timing_center_model_valid": True,
            "timing_center_source": "low_aero_0_5_fit",
            "offset_refinement_valid": True,
            "offset_refinement_applied": True,
            "predicted_kaon_mean": 0.0,
            "predicted_proton_mean": 0.5,
            "kaon_sigma": 0.15,
            "proton_sigma": 0.15,
        }

    def _slice_prefit(self, histogram, global_shape=None, timing_constraint=None, support_entries=100):
        return proton_cleaning._fit_delta_timing_slice(
            histogram,
            self.global_shape if global_shape is None else global_shape,
            self.config,
            "f_runtime_status",
            support_entries=support_entries,
            timing_constraint=self.timing_constraint if timing_constraint is None else timing_constraint,
        )

    def test_delta_offset_reads_only_current_tof_summary(self):
        histogram = _FakeHistogram()
        tof_summary = {"valid": True, "mean_delta_t_pk_ns": 0.5, "prompt_event_count": 100}
        with (
            mock.patch.object(proton_cleaning.ROOT, "TF1", _FakeTF1),
            mock.patch.object(proton_cleaning, "_find_peak_seed", return_value=(5.0, 0.0)),
            mock.patch.object(proton_cleaning, "_compute_poisson_goodness_of_fit", return_value={"deviance": 1.0, "ndf": 1, "deviance_ndf": 1.0, "deviance_per_entry": 0.1}),
            mock.patch.object(proton_cleaning, "_is_near_bound", return_value=False),
        ):
            result = proton_cleaning._fit_delta_common_timing_offset(
                histogram,
                self.global_shape,
                tof_summary,
                self.config,
                "f_delta_offset_runtime",
                False,
                "ct",
                (-2.0, 2.0),
                2.0,
                support_entries=100,
            )
        self.assertTrue(result["fit_attempted"])

    def test_prefit_status_categories_are_distinct(self):
        self.assertEqual(self._slice_prefit(None)["fit_status"], "missing_histogram")
        self.assertEqual(
            self._slice_prefit(object(), global_shape={"valid": False})["fit_status"],
            "invalid_global_shape",
        )
        self.assertEqual(
            self._slice_prefit(object(), timing_constraint={"valid": False})["fit_status"],
            "invalid_timing_center_model",
        )
        low_support = self._slice_prefit(object(), support_entries=1)
        self.assertEqual(low_support["fit_status"], "insufficient_support")
        self.assertIn("insufficient_entries", low_support["rejection_reasons"])

    def test_selected_tof_summary_tracks_actual_offset_source(self):
        primary = {"offset_fit_aero_mode": "low_aero_0_5", "valid": True}
        fallback = {"offset_fit_aero_mode": "low_aero_0_6_fallback", "valid": True}
        selected = proton_cleaning._build_selected_tof_summary(
            {"selected_timing_center_source": "low_aero_0_6_fit"},
            primary,
            fallback,
        )
        self.assertEqual(selected["offset_fit_aero_mode"], "low_aero_0_6_fallback")
        self.assertEqual(selected["tof_summary_role"], "offset_fit_input")

    def test_stable_fallback_keeps_supported_cell_eligible_for_fit(self):
        stable_offset = {
            "selected_timing_center_source": "stable_global_center_fallback",
            "timing_center_model_valid": True,
            "offset_refinement_valid": False,
            "offset_refinement_applied": False,
            "delta_offset": 0.0,
        }
        constraint = proton_cleaning._build_timing_constraint_for_cell(
            self.global_shape,
            stable_offset,
            {},
            False,
            "ct",
            (-2.0, 2.0),
            2.0,
        )
        self.assertTrue(constraint["valid"])
        self.assertEqual(constraint["timing_center_source"], "stable_global_center_fallback")
        with (
            mock.patch.object(proton_cleaning.ROOT, "TF1", _FakeTF1),
            mock.patch.object(proton_cleaning, "_find_peak_seed", return_value=(5.0, 0.0)),
            mock.patch.object(proton_cleaning, "_extract_root_fit_matrices", return_value=({}, {}, {})),
            mock.patch.object(proton_cleaning, "_compute_poisson_goodness_of_fit", return_value={"deviance": 1.0, "ndf": 1, "deviance_ndf": 1.0, "deviance_per_entry": 0.1}),
        ):
            result = self._slice_prefit(_FakeHistogram(), timing_constraint=constraint, support_entries=100)
        self.assertTrue(result["fit_attempted"])

    def test_stable_fallback_can_require_a_local_timing_offset(self):
        stable_offset = {
            "selected_timing_center_source": "stable_global_center_fallback",
            "timing_center_model_valid": True,
            "offset_refinement_valid": False,
            "offset_refinement_applied": False,
            "delta_offset": 0.0,
        }
        constraint = proton_cleaning._build_timing_constraint_for_cell(
            self.global_shape,
            stable_offset,
            {},
            False,
            "ct",
            (-2.0, 2.0),
            2.0,
            allow_stable_global_center_fallback=False,
        )
        self.assertFalse(constraint["valid"])
        self.assertEqual(constraint["reason"], "local_low_aero_timing_offset_required")

    def test_skip_counter_mapping_is_exact(self):
        self.assertEqual(
            proton_cleaning._cell_fit_skip_counter_key("invalid_global_shape"),
            "cell_fit_skipped_invalid_global_shape_count",
        )
        self.assertEqual(
            proton_cleaning._cell_fit_skip_counter_key("insufficient_support"),
            "cell_fit_skipped_insufficient_support_count",
        )

    def test_timing_t_support_hierarchy_rejects_isolated_cell(self):
        thresholds = {
            "minimum_supported_t_slices": 2,
            "minimum_marginal_t_slices": 1,
            "minimum_supported_coverage": 0.35,
            "minimum_marginal_coverage": 0.15,
            "minimum_modeled_yield": 5.0,
            "minimum_setting_valid_t_cells": 2,
        }
        self.assertEqual(
            proton_cleaning._classify_timing_t_support(1, 0.2, 10.0, thresholds),
            proton_cleaning.SUPPORT_MARGINAL,
        )
        gate = proton_cleaning._classify_timing_t_setting_support(
            [{"support_label": proton_cleaning.SUPPORT_MARGINAL, "valid_t_cells": 1,
              "data_total": 10.0, "fitted_data_total": 2.0, "model_total": 10.0}],
            thresholds,
        )
        self.assertFalse(gate["accepted"])
        self.assertEqual(gate["support_label"], proton_cleaning.SUPPORT_UNSUPPORTED)

    def test_delta_center_interpolation_is_a_pure_timing_center_payload(self):
        left = {
            "valid": True, "kaon_mean": -0.1, "proton_mean": 0.6,
            "kaon_sigma": 0.1, "proton_sigma": 0.2, "proton_yield": 9.0,
        }
        right = {
            "valid": True, "kaon_mean": 0.1, "proton_mean": 1.0,
            "kaon_sigma": 0.3, "proton_sigma": 0.4, "proton_yield": 99.0,
        }
        interpolated = proton_cleaning._interpolate_delta_timing_center_shape(
            left, right, -1.0, 0.0, 1.0,
        )
        self.assertAlmostEqual(interpolated["kaon_mean"], 0.0)
        self.assertAlmostEqual(interpolated["proton_mean"], 0.8)
        self.assertAlmostEqual(interpolated["kaon_sigma"], 0.2)
        self.assertAlmostEqual(interpolated["proton_sigma"], 0.3)
        self.assertEqual(interpolated["timing_center_source"], "neighbor_delta_interpolation")
        self.assertNotIn("proton_yield", interpolated)
        self.assertNotIn("proton_amplitude", interpolated)

    def test_nearest_timing_center_is_pure_and_has_provenance(self):
        nearest = proton_cleaning._nearest_delta_timing_center_shape(
            {
                "valid": True, "fit_min": -2.0, "fit_max": 2.0,
                "kaon_mean": -0.1, "proton_mean": 0.8,
                "kaon_sigma": 0.2, "proton_sigma": 0.3,
                "proton_yield": 25.0,
            },
            2, 0.25, 0.75,
        )
        self.assertEqual(nearest["timing_center_source"], "neighbor_delta_nearest_fallback")
        self.assertEqual(nearest["nearest_neighbor_index"], 2)
        self.assertNotIn("proton_yield", nearest)

    def test_applied_cell_map_zeroes_unsupported_and_weak_cells(self):
        cells = [[
            {
                "valid": True, "cell_fit_valid": True,
                "proton_component_detected": True,
                "raw_proton_yield": 9.0, "proton_yield": 7.0,
            },
            {
                "valid": True, "cell_fit_valid": True,
                "proton_component_detected": False,
                "raw_proton_yield": 11.0, "proton_yield": 0.0,
            },
        ]]
        proton_cleaning._apply_timing_t_cell_map(
            cells, [proton_cleaning.SUPPORT_UNSUPPORTED],
            {"support_label": proton_cleaning.SUPPORT_UNSUPPORTED, "accepted": False},
        )
        self.assertEqual(cells[0][0]["raw_proton_yield"], 9.0)
        self.assertEqual(cells[0][0]["applied_proton_yield"], 0.0)
        self.assertEqual(cells[0][0]["applied_zero_reason"], "unsupported_delta")
        self.assertEqual(cells[0][1]["applied_zero_reason"], "weak_proton_component")

    def test_candidate_score_prefers_production_support_over_global_fit(self):
        rf_score = proton_cleaning._timing_t_candidate_selection_score(
            {
                "accepted": False, "support_label": proton_cleaning.SUPPORT_UNSUPPORTED,
                "supported_delta_bins": 0, "marginal_delta_bins": 0,
                "valid_t_cells": 1, "coverage": 0.2,
            },
            {"valid": True, "proton_component_significance": 20.0, "poisson_deviance_per_entry": 0.01},
            2,
        )
        ct_score = proton_cleaning._timing_t_candidate_selection_score(
            {
                "accepted": True, "support_label": proton_cleaning.SUPPORT_MARGINAL,
                "supported_delta_bins": 0, "marginal_delta_bins": 1,
                "valid_t_cells": 2, "coverage": 0.2,
            },
            {"valid": True, "proton_component_significance": 2.0, "poisson_deviance_per_entry": 0.10},
            1,
        )
        self.assertGreater(ct_score, rf_score)
        ranked = proton_cleaning._rank_timing_t_candidate_evaluations([
            {"candidate": {"timing_branch": "RF"}, "candidate_selection_score": rf_score},
            {"candidate": {"timing_branch": "CTime_ROC1"}, "candidate_selection_score": ct_score},
        ])
        self.assertEqual(ranked[0]["candidate"]["timing_branch"], "CTime_ROC1")

    def test_ct_candidate_uses_exact_epsilon_ranges(self):
        low = proton_cleaning.resolve_timing_t_candidate_configuration(
            {"epsset": "low"}, self.config, "CTime_ROC1", None, None, 0.7, 1.5,
        )
        high = proton_cleaning.resolve_timing_t_candidate_configuration(
            {"epsset": "high"}, self.config, "CTime_ROC1", None, None, 0.7, 1.5,
        )
        self.assertEqual(low["display_range"], (-2.0, 2.0))
        self.assertEqual(low["histogram_bins"], 131)
        self.assertEqual(high["display_range"], (-4.0, 4.0))
        self.assertEqual(high["histogram_bins"], 262)
        self.assertFalse(low["proton_peak_is_lower"])

    def test_aerogel_validation_is_observational(self):
        result = {"t_edges": [0.0, 1.0]}
        rows = [{
            "t_index": 0, "aero_index": 0, "physical_weight": 2.0,
            "proton_weight": 0.25, "adj_mm": 0.85,
        }]
        config = {
            "aerogel_validation": {
                "enabled": True, "slice_edges": (0.0, 5.0, 10.0),
                "affects_event_weights": False, "affects_fit_acceptance": False,
            },
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, config)
        self.assertFalse(validation["affects_event_weights"])
        self.assertFalse(validation["affects_fit_acceptance"])
        self.assertEqual(validation["event_count_by_t_aero"][0][0], 1)
        self.assertAlmostEqual(validation["average_proton_probability_by_t_aero"][0][0], 0.25)

    def test_aerogel_validation_records_signed_support_and_tbin_summary(self):
        result = {"t_edges": [0.0, 1.0, 2.0]}
        rows = [
            {
                "t_index": 0, "aero_index": 0, "source_label": "prompt",
                "is_prompt_source": True, "physical_weight": 2.0,
                "proton_weight": 0.50, "cleaned_factor": 0.50, "adj_mm": 0.85,
            },
            {
                "t_index": 0, "aero_index": 1, "source_label": "rand",
                "is_prompt_source": False, "physical_weight": -0.5,
                "proton_weight": 0.20, "cleaned_factor": 0.80, "adj_mm": 1.115,
            },
        ]
        config = {
            "aerogel_validation": {
                "enabled": True,
                "summary_slice_edges": (0.0, 5.0, 10.0),
                "minimum_events_per_t_bin": 3,
                "affects_event_weights": False,
                "affects_fit_acceptance": False,
            },
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, config)
        self.assertEqual(validation["raw_prompt_event_count_by_t_aero"][0][0], 1)
        self.assertAlmostEqual(validation["signed_event_weight_sum_by_t_aero"][0][1], -0.5)
        self.assertAlmostEqual(validation["proton_probability_sum_by_t_aero"][0][0], 1.0)
        self.assertAlmostEqual(validation["low_mm_removed_yield_by_t_aero"][0][0], 1.0)
        self.assertAlmostEqual(validation["lambda_removed_yield_by_t_aero"][0][1], -0.1)
        self.assertEqual(validation["signed_weight_diagnostics"]["positive_event_count"], 1)
        self.assertEqual(validation["signed_weight_diagnostics"]["negative_event_count"], 1)
        self.assertIn("insufficient_aerogel_support", validation["per_t_bin_summary"][0]["warnings"])

    def test_aerogel_summary_edges_preserve_legacy_slice_override(self):
        self.assertEqual(
            proton_cleaning._resolve_aerogel_summary_edges({"slice_edges": (0.0, 5.0, 10.0)}),
            [0.0, 5.0, 10.0],
        )
        self.assertEqual(
            proton_cleaning._resolve_aerogel_summary_edges({
                "summary_slice_edges": (0.0, 3.0, 9.0),
                "slice_edges": (0.0, 5.0, 10.0),
            }),
            [0.0, 3.0, 9.0],
        )


if __name__ == "__main__":
    unittest.main()
