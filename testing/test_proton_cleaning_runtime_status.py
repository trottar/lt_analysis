"""Focused regression checks for proton-cleaning runtime-state handling."""

from __future__ import annotations

import sys
import types
import unittest
from copy import deepcopy
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
        "TCanvas", "TLatex", "TLegend", "TLine", "TPad", "TPaveText", "TH1D", "TH2D",
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


class _RecordingHistogram:
    """Minimal ROOT-histogram stand-in that records physical Fill coordinates."""

    def __init__(self):
        self.fills = []

    def Fill(self, *values):
        self.fills.append(tuple(float(value) for value in values))


class _ApplicationHistogram(_RecordingHistogram):
    """Enough TH1/TH2 behavior to exercise the application loop without ROOT."""

    instances = []

    def __init__(self, name, _title, bins_x, *_axis_args):
        super().__init__()
        self.name = str(name)
        self.bins_x = int(bins_x)
        self.bins_y = int(_axis_args[2]) if len(_axis_args) >= 3 else 0
        self.contents = {}
        self.__class__.instances.append(self)

    def SetDirectory(self, _directory):
        pass

    def Sumw2(self):
        pass

    def Clone(self, name):
        cloned = object.__new__(self.__class__)
        cloned.fills = []
        cloned.name = str(name)
        cloned.bins_x = self.bins_x
        cloned.bins_y = self.bins_y
        cloned.contents = {}
        self.__class__.instances.append(cloned)
        return cloned

    def Reset(self):
        self.fills = []
        self.contents = {}

    def GetNbinsX(self):
        return self.bins_x

    def GetNbinsY(self):
        return self.bins_y

    def GetBinContent(self, *indices):
        return self.contents.get(tuple(int(index) for index in indices), 0.0)

    def SetBinContent(self, *args):
        *indices, value = args
        self.contents[tuple(int(index) for index in indices)] = float(value)


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

    def test_timing_t_pdf_helpers_import_direct_th1d(self):
        # Cross-stage closure pages use the direct ROOT alias rather than
        # ROOT.TH1D, so retain this import-level regression guard.
        self.assertTrue(hasattr(proton_cleaning, "TH1D"))

    def test_post_proton_rf_requires_explicit_enablement(self):
        config = {
            "apply_post_proton_rf": False,
            "apply_only_low_epsilon_rf": True,
        }
        policy, enabled, applied = proton_cleaning._resolve_post_proton_rf_application(
            config, {"EPSSET": "low"}, None
        )
        self.assertEqual(policy, "epsset_default_after_cleaning")
        self.assertFalse(enabled)
        self.assertFalse(applied)

        config["apply_post_proton_rf"] = True
        _policy, enabled, applied = proton_cleaning._resolve_post_proton_rf_application(
            config, {"EPSSET": "low"}, None
        )
        self.assertTrue(enabled)
        self.assertTrue(applied)

    def test_proton_weight_diagnostics_fill_physical_delta_and_t_coordinates(self):
        delta_sum = _RecordingHistogram()
        delta_norm = _RecordingHistogram()
        delta_t_sum = _RecordingHistogram()
        delta_t_norm = _RecordingHistogram()
        proton_cleaning._fill_proton_weight_coordinate_histograms(
            h_weight_sum_delta=delta_sum,
            h_weight_norm_delta=delta_norm,
            h_weight_sum_delta_secondary=delta_t_sum,
            h_weight_norm_delta_secondary=delta_t_norm,
            physical_delta_value=3.75,
            physical_secondary_value=0.642,
            delta_index=1,
            secondary_index=2,
            delta_edges=(-10.0, 0.0, 10.0),
            secondary_edges=(0.0, 0.4, 0.6, 0.8),
            coefficient=-2.0,
            proton_weight=0.3,
        )
        self.assertEqual(delta_sum.fills, [(3.75, 0.6)])
        self.assertEqual(delta_norm.fills, [(3.75, 2.0)])
        self.assertEqual(delta_t_sum.fills, [(3.75, 0.642, 0.6)])
        self.assertEqual(delta_t_norm.fills, [(3.75, 0.642, 2.0)])
        self.assertNotIn((1.0, 0.6), delta_sum.fills)
        self.assertNotIn((1.0, 2.0, 0.6), delta_t_sum.fills)

    def test_proton_weight_diagnostic_indices_only_control_lookup_eligibility(self):
        delta_sum = _RecordingHistogram()
        delta_norm = _RecordingHistogram()
        delta_t_sum = _RecordingHistogram()
        delta_t_norm = _RecordingHistogram()
        proton_cleaning._fill_proton_weight_coordinate_histograms(
            h_weight_sum_delta=delta_sum,
            h_weight_norm_delta=delta_norm,
            h_weight_sum_delta_secondary=delta_t_sum,
            h_weight_norm_delta_secondary=delta_t_norm,
            physical_delta_value=-4.5,
            physical_secondary_value=0.8,
            delta_index=-1,
            secondary_index=0,
            delta_edges=(-10.0, 0.0, 10.0),
            secondary_edges=(0.0, 1.0),
            coefficient=1.0,
            proton_weight=0.4,
        )
        self.assertEqual(delta_sum.fills, [])
        self.assertEqual(delta_t_sum.fills, [])

    def test_frozen_coordinate_audit_uses_lookup_edges_and_skips_invalid_indices(self):
        integrity = proton_cleaning._empty_frozen_coordinate_integrity(1.0e-9)
        proton_cleaning._record_frozen_coordinate_membership(
            integrity,
            coordinate_name="delta",
            physical_value=10.0 + 5.0e-10,
            stored_index=1,
            edges=(-10.0, 0.0, 10.0),
        )
        proton_cleaning._record_frozen_coordinate_membership(
            integrity,
            coordinate_name="secondary",
            physical_value=0.5,
            stored_index=-1,
            edges=(0.0, 0.4, 0.8),
        )
        proton_cleaning._record_frozen_coordinate_membership(
            integrity,
            coordinate_name="delta",
            physical_value=-0.1,
            stored_index=1,
            edges=(-10.0, 0.0, 10.0),
        )
        self.assertEqual(integrity["delta_checked_count"], 2)
        self.assertEqual(integrity["secondary_checked_count"], 0)
        self.assertEqual(integrity["skipped_invalid_index_count"], 1)
        self.assertEqual(integrity["skipped_invalid_secondary_index_count"], 1)
        self.assertEqual(integrity["delta_mismatch_count"], 1)
        self.assertGreater(integrity["maximum_delta_edge_violation"], 0.0)

    def test_application_coordinate_fix_leaves_committed_lookup_fields_unchanged(self):
        class _Event:
            pass

        event = _Event()
        event.P_hgcer_npeSum = 0.0
        event.P_hgcer_xAtCer = 0.0
        event.P_hgcer_yAtCer = 0.0
        frozen_payload = {
            "delta_index": 0,
            "t_index": 1,
            "support_label": proton_cleaning.SUPPORT_SUPPORTED,
            "proton_weight": 0.25,
            "cleaned_factor": 0.75,
            "final_cleaned_factor": 0.60,
            "rf_accept": True,
            "lambda_gate_status": "pass",
            "lambda_gate_passed": True,
        }
        lookup_before = deepcopy(frozen_payload)
        cleaning_result = {
            "accepted": True,
            "method": proton_cleaning.PROTON_CONTAMINATION_CLEANING_METHOD_TIMING_T_EVENT_WEIGHT,
            "t_edges": [0.0, 0.5, 1.0],
            "delta_edges": [-10.0, 0.0, 10.0],
            "settings": {"strict_mode": True, "t_binning": {"edge_tolerance": 1.0e-9}},
            "selected_timing_branch": "CTime_ROC1",
            "diagnostics": {},
            "_prepared_event_weight_lookup": {"prompt:0": frozen_payload},
        }
        source_bundle = {
            "sources": {"prompt": {"tree": [event], "coefficient": 2.0}},
            "prepared_sources": {
                "prompt": {
                    "entries": {
                        0: {
                            "allcuts": False,
                            "nommcuts": True,
                            "noholecuts": False,
                            "adj_hsdelta": 0.0,
                            "adj_mm": 1.115,
                            "adj_t": 0.75,
                            "delta_value": -4.25,
                            "aero_value": 8.0,
                            "timing_values": {"CTime_ROC1": 0.1},
                        }
                    }
                }
            },
        }
        _ApplicationHistogram.instances = []
        with (
            mock.patch.object(proton_cleaning.ROOT, "TH1D", _ApplicationHistogram),
            mock.patch.object(proton_cleaning.ROOT, "TH2D", _ApplicationHistogram),
        ):
            application = proton_cleaning.apply_kaon_proton_cleaning_to_targets(
                cleaning_result,
                source_bundle,
                {},
                lambda *_args: (True, True, 0.0),
                lambda _evt: 1.115,
                lambda _evt: 0.75,
                None,
                0.9,
                1.2,
            )

        self.assertTrue(application["accepted"])
        self.assertEqual(frozen_payload, lookup_before)
        self.assertEqual(frozen_payload["proton_weight"], 0.25)
        self.assertEqual(frozen_payload["cleaned_factor"], 0.75)
        self.assertEqual(frozen_payload["final_cleaned_factor"], 0.60)
        self.assertTrue(frozen_payload["rf_accept"])
        self.assertEqual(frozen_payload["lambda_gate_status"], "pass")
        self.assertTrue(frozen_payload["lambda_gate_passed"])
        self.assertAlmostEqual(2.0 * frozen_payload["proton_weight"], 0.5)
        self.assertAlmostEqual(2.0 * frozen_payload["final_cleaned_factor"], 1.2)
        histogram_by_name = {hist.name: hist for hist in _ApplicationHistogram.instances}
        self.assertEqual(
            histogram_by_name["H_proton_weight_sum_delta"].fills,
            [(-4.25, 0.5)],
        )
        self.assertEqual(
            histogram_by_name["H_proton_weight_sum_delta_t"].fills,
            [(-4.25, 0.75, 0.5)],
        )
        integrity = application["diagnostics"]["frozen_coordinate_integrity"]
        self.assertTrue(integrity["passed"])
        self.assertEqual(integrity["checked_events"], 1)

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
            "proton_weight": 0.25, "aero_value": 2.0, "adj_mm": 0.85,
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
                "proton_weight": 0.50, "cleaned_factor": 0.50, "aero_value": 2.0, "adj_mm": 0.85,
            },
            {
                "t_index": 0, "aero_index": 1, "source_label": "rand",
                "is_prompt_source": False, "physical_weight": -0.5,
                "proton_weight": 0.20, "cleaned_factor": 0.80, "aero_value": 7.0, "adj_mm": 1.115,
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

    def test_aerogel_thresholds_use_event_npe_not_coarse_bins(self):
        result = {"t_edges": [0.0, 1.0]}
        rows = [
            {"t_index": 0, "aero_value": 4.0, "is_prompt_source": True,
             "physical_weight": 2.0, "proton_weight": 0.5, "cleaned_factor": 0.5, "adj_mm": 0.85},
            {"t_index": 0, "aero_value": 5.0, "is_prompt_source": True,
             "physical_weight": 3.0, "proton_weight": 0.4, "cleaned_factor": 0.6, "adj_mm": 0.85},
            {"t_index": 0, "aero_value": 10.0, "is_prompt_source": True,
             "physical_weight": 4.0, "proton_weight": 0.25, "cleaned_factor": 0.75, "adj_mm": 1.115},
        ]
        base_config = {
            "aerogel_validation": {
                "enabled": True, "summary_slice_edges": (0.0, 3.0, 6.0, 10.0, 25.0),
                "low_reference_max_npe": 5.0, "high_reference_min_npe": 10.0,
            },
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, base_config)
        self.assertEqual(validation["low_aero_raw_prompt_count"], 1)
        self.assertEqual(validation["high_aero_raw_prompt_count"], 1)
        self.assertAlmostEqual(validation["low_aero_signed_yield"], 2.0)
        self.assertAlmostEqual(validation["high_aero_signed_yield"], 4.0)
        changed = dict(base_config)
        changed["aerogel_validation"] = dict(base_config["aerogel_validation"], summary_slice_edges=(0.0, 2.0, 8.0, 25.0))
        changed_validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, changed)
        self.assertEqual(changed_validation["low_aero_raw_prompt_count"], 1)
        self.assertEqual(changed_validation["high_aero_raw_prompt_count"], 1)

    def test_aerogel_ratio_preserves_undefined_and_measured_zero(self):
        result = {"t_edges": [0.0, 1.0]}
        rows = [{
            "t_index": 0, "aero_value": 1.0, "is_prompt_source": True,
            "physical_weight": 2.0, "proton_weight": 0.0, "cleaned_factor": 1.0, "adj_mm": 1.0,
        }]
        config = {
            "aerogel_validation": {"enabled": True, "summary_slice_edges": (0.0, 5.0, 10.0)},
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, config)
        self.assertEqual(validation["average_proton_probability_by_t_aero"][0][0], 0.0)
        self.assertTrue(validation["cells"][0][0]["average_proton_probability_valid"])
        self.assertIsNone(validation["average_proton_probability_by_t_aero"][0][1])
        self.assertFalse(validation["cells"][0][1]["average_proton_probability_valid"])

    def test_aerogel_matrix_payload_uses_named_prompt_and_physical_metrics(self):
        result = {"t_edges": [0.0, 1.0], "diagnostics": {"selected_timing_candidate": {"selected": True}}}
        rows = [
            {"t_index": 0, "aero_value": 2.0, "is_prompt_source": True,
             "physical_weight": 2.0, "proton_weight": 0.25, "cleaned_factor": 0.75, "adj_mm": 0.85},
            {"t_index": 0, "aero_value": 7.0, "is_prompt_source": False,
             "physical_weight": -0.5, "proton_weight": 0.20, "cleaned_factor": 0.80, "adj_mm": 1.115},
        ]
        config = {
            "aerogel_validation": {"enabled": True, "summary_slice_edges": (0.0, 5.0, 10.0)},
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, rows, config)
        matrix = validation["matrix_payload"]
        self.assertEqual(matrix["source"], "frozen_timing_t_lookup_rows")
        self.assertEqual(matrix["metrics"]["selected_prompt_count"][0], [1, 0])
        self.assertEqual(matrix["metrics"]["signed_physical_yield"][0], [2.0, -0.5])
        self.assertEqual(matrix["metrics"]["estimated_proton_yield"][0], [0.5, -0.1])
        self.assertEqual(
            validation["raw_prompt_event_count_by_t_aero"],
            matrix["metrics"]["selected_prompt_count"],
        )
        self.assertAlmostEqual(
            matrix["metrics"]["cleaned_yield"][0][0],
            matrix["metrics"]["signed_physical_yield"][0][0]
            - matrix["metrics"]["estimated_proton_yield"][0][0],
        )

    def test_integrity_rejects_nonempty_source_with_no_coarse_matrix_content(self):
        result = {"t_edges": [0.0, 1.0], "diagnostics": {"selected_timing_candidate": {"selected": True}}}
        config = {
            "aerogel_validation": {"enabled": True, "summary_slice_edges": (0.0, 5.0, 10.0)},
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        validation = proton_cleaning._build_t_aerogel_validation(result, {}, [{
            "t_index": 0, "aero_value": 50.0, "is_prompt_source": True,
            "physical_weight": 1.0, "proton_weight": 0.2, "adj_mm": 0.85,
        }], config)
        self.assertIn(
            "nonempty_frozen_source_has_matrix_content",
            validation["diagnostic_integrity"]["failures"],
        )

    def test_integrity_strict_mode_raises_for_invalid_diagnostic_contract(self):
        result = {"t_edges": [0.0, 1.0], "diagnostics": {"selected_timing_candidate": {"selected": True}}}
        config = {
            "aerogel_validation": {
                "enabled": True, "summary_slice_edges": (0.0, 5.0, 10.0),
                "diagnostic_strict": True,
            },
            "validation_windows": {"low_mm": (0.8, 0.9), "lambda_peak": (1.105, 1.125)},
        }
        with self.assertRaisesRegex(RuntimeError, "timing-t diagnostic integrity failed"):
            proton_cleaning._build_t_aerogel_validation(result, {}, [{
                "t_index": 0, "aero_value": 30.0, "physical_weight": 1.0,
                "proton_weight": 0.2, "adj_mm": 0.85,
            }], config)

    def test_timing_t_summary_uses_selected_candidate_state_not_legacy_arrays(self):
        result = {
            "delta_edges": [-1.0, 0.0, 1.0],
            "diagnostics": {
                "selected_timing_candidate": {"selected": True, "timing_branch": "P_RF_Dist", "candidate_selection_rank": [1]},
                "candidate_selection_tuple": ["setting_accepted"],
                "setting_support": {"accepted": True, "support_label": "supported"},
                "delta_support": [{
                    "delta_index": 0, "support_label": "supported", "data_total": 10.0,
                    "fitted_data_total": 9.0, "model_total": 9.0, "proton_total": 3.0,
                    "kaon_total": 5.0, "other_total": 1.0, "valid_t_cells": 2, "coverage": 0.9,
                }],
                "applied_timing_t_cell_map": [{"delta_index": 0, "applied_proton_yield": 2.0}],
                "event_weight_closure_by_delta": [{
                    "delta_index": 0, "summed_event_proton_probability": 1.0, "event_count": 4,
                }],
                "proton_yield_by_delta": [999.0],
            },
        }
        summary = proton_cleaning._build_timing_t_summary(
            result, {"source": "frozen_timing_t_lookup_rows"},
            {"raw_missing_mass_yield": 7.0, "estimated_proton_missing_mass_yield": 2.0, "cleaned_missing_mass_yield": 5.0},
            {
                "source": "post_lambda_gate_frozen_lookup_rows",
                "global": {"committed_applied_proton_yield": 0.0},
                "per_delta": [{
                    "delta_index": 0,
                    "proposed_frozen_lookup_proton_yield": 1.0,
                    "proposed_frozen_lookup_proton_fraction": 0.1,
                    "proposed_absolute_support_mean_probability": 0.25,
                    "committed_applied_proton_yield": 0.0,
                    "committed_applied_proton_fraction": 0.0,
                    "committed_absolute_support_mean_probability": 0.0,
                }],
            },
        )
        self.assertEqual(summary["selected_candidate"]["timing_branch"], "P_RF_Dist")
        delta = summary["per_delta"][0]
        self.assertEqual(delta["raw_fitted_proton_yield"], 3.0)
        self.assertEqual(delta["proposed_fit_proton_yield"], 2.0)
        self.assertEqual(delta["committed_applied_proton_yield"], 0.0)
        self.assertAlmostEqual(
            delta["mean_proposed_frozen_wp_signed_coefficient_contribution_per_event"], 0.25
        )
        self.assertEqual(
            delta["mean_proposed_frozen_wp_label"],
            proton_cleaning.TIMING_T_PROPOSED_FROZEN_WP_SIGNED_CONTRIBUTION_LABEL,
        )

    def test_frozen_lookup_stage_summary_keeps_proposed_and_committed_yields_distinct(self):
        summary = proton_cleaning._build_timing_t_frozen_lookup_stage_summary([
            {
                "delta_index": 0, "physical_weight": 2.0,
                "proposed_proton_probability": 0.25,
                "applied_proton_probability": 0.0,
            },
            {
                "delta_index": 0, "physical_weight": -0.5,
                "proposed_proton_probability": 0.20,
                "applied_proton_probability": 0.0,
            },
        ])
        global_summary = summary["global"]
        self.assertAlmostEqual(global_summary["raw_signed_yield"], 1.5)
        self.assertAlmostEqual(global_summary["proposed_frozen_lookup_proton_yield"], 0.4)
        self.assertEqual(global_summary["committed_applied_proton_yield"], 0.0)
        self.assertAlmostEqual(global_summary["proposed_absolute_support_mean_probability"], 0.24)
        self.assertEqual(global_summary["committed_absolute_support_mean_probability"], 0.0)

    def test_timing_t_mm_diagnostics_use_only_frozen_lookup_rows(self):
        result = {"t_edges": [0.0, 1.0, 2.0]}
        rows = [
            {
                "t_index": 0, "is_prompt_source": True,
                "physical_weight": 2.0, "proton_weight": 0.25,
                "cleaned_factor": 0.75, "final_cleaned_factor": 0.75,
                "adj_mm": 0.85,
            },
            {
                "t_index": 0, "is_prompt_source": False,
                "physical_weight": -0.5, "proton_weight": 0.20,
                "cleaned_factor": 0.80, "final_cleaned_factor": 0.80,
                "adj_mm": 1.115,
            },
        ]
        config = {
            "mm_diagnostics": {
                "enabled": True,
                "display_range": (0.70, 1.50),
                "display_bins": 80,
                "affects_event_weights": False,
                "affects_fit_acceptance": False,
            },
            "validation_windows": {
                "low_mm": (0.80, 0.90),
                "lambda_peak": (1.105, 1.125),
            },
        }
        payload = proton_cleaning._build_timing_t_mm_diagnostics(result, rows, config)
        self.assertEqual(payload["source"], "frozen_timing_t_lookup_rows")
        self.assertFalse(payload["affects_event_weights"])
        self.assertFalse(payload["affects_fit_acceptance"])
        self.assertEqual(payload["aggregate"]["event_count"], 2)
        self.assertEqual(payload["aggregate"]["raw_prompt_event_count"], 1)
        self.assertAlmostEqual(payload["aggregate"]["raw_missing_mass_yield"], 1.5)
        self.assertAlmostEqual(payload["aggregate"]["estimated_proton_missing_mass_yield"], 0.4)
        self.assertAlmostEqual(payload["aggregate"]["cleaned_missing_mass_yield"], 1.1)
        self.assertAlmostEqual(payload["aggregate"]["pre_rf_cleaning_closure_difference"], 0.0)
        low = payload["per_t_bin_summary"][0]["windows"]["low_mm"]
        self.assertAlmostEqual(low["removed_fraction"], 0.25)
        self.assertTrue(low["removed_fraction_valid"])
        lambda_window = payload["per_t_bin_summary"][0]["windows"]["lambda_peak"]
        self.assertAlmostEqual(lambda_window["removed_fraction"], 0.20)

    def test_timing_t_mm_diagnostics_preserve_undefined_window_fraction(self):
        payload = proton_cleaning._build_timing_t_mm_diagnostics(
            {"t_edges": [0.0, 1.0]},
            [{
                "t_index": 0, "is_prompt_source": True,
                "physical_weight": 1.0, "proton_weight": 0.0,
                "cleaned_factor": 1.0, "adj_mm": 1.0,
            }],
            {
                "mm_diagnostics": {"enabled": True},
                "validation_windows": {"low_mm": (0.80, 0.90)},
            },
        )
        low = payload["per_t_bin_summary"][0]["windows"]["low_mm"]
        self.assertIsNone(low["removed_fraction"])
        self.assertFalse(low["removed_fraction_valid"])

    def test_timing_t_summary_prefers_frozen_mm_diagnostic_totals(self):
        result = {
            "delta_edges": [0.0, 1.0],
            "diagnostics": {
                "timing_t_mm_diagnostics": {
                    "aggregate": {
                        "raw_missing_mass_yield": 11.0,
                        "estimated_proton_missing_mass_yield": 4.0,
                        "cleaned_missing_mass_yield": 7.0,
                    },
                },
            },
        }
        summary = proton_cleaning._build_timing_t_summary(
            result, {},
            {"raw_missing_mass_yield": 99.0, "estimated_proton_missing_mass_yield": 99.0, "cleaned_missing_mass_yield": 99.0},
        )
        self.assertEqual(summary["missing_mass_totals"]["raw"], 11.0)
        self.assertEqual(summary["missing_mass_totals"]["estimated_proton"], 4.0)
        self.assertEqual(summary["missing_mass_totals"]["cleaned"], 7.0)

    def test_timing_t_layout_cross_stage_and_per_t_eligibility_helpers(self):
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(1), (1, 1))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(2), (2, 1))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(4), (2, 2))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(6), (3, 2))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(8), (4, 2))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(10), (5, 2))
        self.assertEqual(proton_cleaning._timing_t_layout_for_panel_count(12), (4, 3))
        self.assertTrue(proton_cleaning._cross_stage_visual_state([
            {"maximum_absolute_difference": 0.0},
        ])["render_pass_status"])
        self.assertFalse(proton_cleaning._cross_stage_visual_state([
            {"maximum_absolute_difference": 0.2},
        ], {"cross_stage_visual_threshold": 0.1})["render_pass_status"])
        self.assertTrue(proton_cleaning._timing_t_per_t_pid_eligible({
            "raw_prompt_event_count": 1, "absolute_event_weight_support": 0.2,
        }, {"per_t_absolute_support_tolerance": 0.1})["eligible"])
        self.assertFalse(proton_cleaning._timing_t_per_t_pid_eligible({
            "raw_prompt_event_count": 0, "absolute_event_weight_support": 10.0,
        })["eligible"])
        self.assertEqual(
            proton_cleaning._timing_t_signed_display_limits(0.0, 0.25),
            (0.0, 0.28),
        )
        signed_limits = proton_cleaning._timing_t_signed_display_limits(-0.10, 0.25)
        self.assertLess(signed_limits[0], -0.10)
        self.assertGreater(signed_limits[1], 0.25)
        self.assertEqual(
            proton_cleaning._timing_t_proposed_model_annotation({"status": "pass"}),
            "PROPOSED MODEL - COMMITTED",
        )
        self.assertEqual(
            proton_cleaning._timing_t_proposed_model_annotation({"status": "fail"}),
            "PROPOSED TIMING MODEL - NOT APPLIED TO PRODUCTION",
        )
        for annotation in (
            proton_cleaning._timing_t_proposed_model_annotation({"status": "pass"}),
            proton_cleaning._timing_t_proposed_model_annotation({"status": "fail"}),
        ):
            self.assertNotIn("\N{EM DASH}", annotation)
            self.assertNotIn("\N{EN DASH}", annotation)
            annotation.encode("ascii")
        pair = proton_cleaning._cross_stage_pair_state(
            [{"prepass_t": 0.4, "prepared_proton_cleaning_adj_t": 0.4}],
            "prepass_t", "prepared_proton_cleaning_adj_t", 1.0e-10,
        )
        self.assertEqual(pair["state"], "PASS")

    def test_timing_t_status_page_fallback_and_strict_rendering_integrity(self):
        class FakeCanvas:
            def __init__(self):
                self.printed = []

            def Modified(self):
                return None

            def Update(self):
                return None

            def Print(self, path):
                self.printed.append(path)

        canvas = FakeCanvas()

        def draw_status(_message, *, color, canvas=None):
            panel = object()
            proton_cleaning._retain_timing_t_root_objects(canvas, panel, status=True)
            return panel

        with mock.patch.object(proton_cleaning, "_draw_timing_t_status_panel", side_effect=draw_status):
            proton_cleaning._print_timing_t_page(
                canvas, "status.pdf", required_status=True,
                fallback_message="missing retained panel",
            )
        self.assertEqual(canvas.printed, ["status.pdf"])
        self.assertNotIn(id(canvas), proton_cleaning._TIMING_T_ROOT_RETAINED_OBJECTS)

        strict_canvas = FakeCanvas()
        with self.assertRaisesRegex(RuntimeError, "missing retained panel"):
            proton_cleaning._print_timing_t_page(
                strict_canvas, "strict.pdf",
                config={"aerogel_validation": {"diagnostic_strict": True}},
                required_status=True, fallback_message="missing retained panel",
            )
        self.assertEqual(strict_canvas.printed, [])

    def test_hgcer_display_audit_preserves_signed_content(self):
        class FakeHistogram:
            def GetNbinsX(self):
                return 2

            def GetNbinsY(self):
                return 2

            def GetBinContent(self, x_bin, y_bin):
                return {(1, 1): 2.0, (1, 2): -3.0}.get((x_bin, y_bin), 0.0)

        result = {
            "method": "timing_t_event_weight",
            "application": {"generic_hgcer_fill_counters": {"final_cleaned": {}}},
            "diagnostics": {},
        }
        audit = proton_cleaning.audit_timing_t_hgcer_display_targets(
            result, {"hgcer_x_mm": FakeHistogram()},
        )
        display = audit["final_display_histograms"]["hgcer_x_mm"]
        self.assertEqual(display["positive_bin_count"], 1)
        self.assertEqual(display["negative_bin_count"], 1)
        self.assertAlmostEqual(display["signed_integral"], -1.0)
        self.assertAlmostEqual(display["absolute_integral"], 5.0)

    def test_hgcer_fill_counter_records_selection_range_and_signed_support(self):
        class Axis:
            def GetXmin(self):
                return -1.0

            def GetXmax(self):
                return 1.0

        class FakeHistogram:
            def GetXaxis(self):
                return Axis()

            def GetYaxis(self):
                return Axis()

        counters = {}
        proton_cleaning._record_hgcer_fill_counter(
            counters, "hgcer_x_mm", FakeHistogram(), 0.25, 0.50, -2.0, True,
        )
        counter = counters["hgcer_x_mm"]
        self.assertEqual(counter["seen"], 1)
        self.assertEqual(counter["selected"], 1)
        self.assertEqual(counter["finite"], 1)
        self.assertEqual(counter["in_range"], 1)
        self.assertEqual(counter["filled"], 1)
        self.assertEqual(counter["nonzero"], 1)
        self.assertEqual(counter["signed_weight_sum"], -2.0)
        self.assertEqual(counter["absolute_weight_sum"], 2.0)


if __name__ == "__main__":
    unittest.main()
