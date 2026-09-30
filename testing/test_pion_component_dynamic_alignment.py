"""Focused regressions for the parent-to-fine pion alignment contract."""

from __future__ import annotations

import copy
import inspect
import math
import os
import sys
import tempfile
import types
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

from background_config import get_pion_component_dynamic_alignment_config

if ROOT is not None:
    import pion_component_fits as fits
else:
    class _RootImportStub(types.ModuleType):
        """Allow deterministic provenance tests to import the ROOT-bound module."""

        def __getattr__(self, name):
            return type(name, (), {})


    def _load_root_free_alignment_module():
        """Load the production predicate without leaking a fake ROOT module."""
        originals = {name: sys.modules.get(name) for name in ("ROOT", "utility", "pion_component_fits")}
        root_stub = _RootImportStub("ROOT")
        for color_name, color_value in {
            "kRed": 632,
            "kAzure": 860,
            "kMagenta": 616,
            "kBlue": 600,
            "kCyan": 432,
        }.items():
            setattr(root_stub, color_name, color_value)
        utility_stub = types.ModuleType("utility")
        utility_stub.normalize_hist_to_unit_area = lambda histogram: histogram
        sys.modules["ROOT"] = root_stub
        sys.modules["utility"] = utility_stub
        sys.modules.pop("pion_component_fits", None)
        try:
            import pion_component_fits as root_free_fits
        finally:
            for name, original in originals.items():
                if original is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = original
        return root_free_fits


    fits = _load_root_free_alignment_module()


class DynamicPionAlignmentConfigTests(unittest.TestCase):
    def test_default_schema_and_inherited_sidis_configuration(self):
        config = get_pion_component_dynamic_alignment_config(
            inp_dict={}, phi_setting="center", setting_key="Q2W-test"
        )
        self.assertEqual(config["alignment_schema_version"], 2)
        self.assertTrue(config["components"]["pi_n"]["enabled"])
        self.assertTrue(config["components"]["pi_delta"]["enabled"])
        self.assertFalse(config["components"]["pi_sidis"]["enabled"])
        self.assertTrue(config["components"]["pi_sidis"]["inherit_parent_by_default"])
        self.assertEqual(config["acceptance"]["data_support_metric"], "positive_integral")
        self.assertTrue(config["acceptance"]["reject_offset_boundary"])
        self.assertFalse(config["acceptance"]["reject_expansion_boundary"])

    def test_fine_workflows_share_the_cache_before_resolve_helper(self):
        for relative_path in ("src/binning/calculate_yield.py", "src/binning/ave_per_bin.py"):
            with open(os.path.join(REPO_ROOT, relative_path), "r") as handle:
                source = handle.read()
            self.assertIn("load_or_resolve_pion_component_alignment", source)


class _AxisStub:
    def __init__(self, edges):
        self._edges = list(edges)

    def GetXmin(self):
        return self._edges[0]

    def GetXmax(self):
        return self._edges[-1]

    def GetBinLowEdge(self, index):
        return self._edges[index - 1]


class _HistogramStub:
    def __init__(self, name, edges, contents):
        self._name = name
        self._axis = _AxisStub(edges)
        self._contents = [0.0, *contents, 0.0]
        self._errors = [0.0, *(math.sqrt(abs(value)) for value in contents), 0.0]

    def GetName(self):
        return self._name

    def GetXaxis(self):
        return self._axis

    def GetNbinsX(self):
        return len(self._contents) - 2

    def GetBinContent(self, index):
        return self._contents[index]

    def GetBinError(self, index):
        return self._errors[index]


class DynamicPionAlignmentPersistenceAndBoundaryTests(unittest.TestCase):
    def _identity(self):
        return {
            "kinematic_setting": "Q2W-test",
            "epsilon": "lowe",
            "phi_setting": "center",
            "analysis_scope": "setting-wide",
            "t_bin": None,
            "phi_bin": None,
        }

    def _pion_control(self, name, edges=(0.70, 0.90, 1.10, 1.30), contents=(12.0, 18.0, 7.0)):
        return _HistogramStub(name, edges, contents)

    def _cache_payload(self, pion_control):
        payload = fits.build_expected_pion_alignment_metadata(
            "Q2W-test",
            "setting-wide",
            self._identity(),
            pion_control,
            {name: None for name in fits.COMPONENT_NAMES},
            inp_dict={},
            phi_setting="center",
            common_setting_shift_gev=0.005,
        )
        payload.update({"accepted": True, "selected_score": 1.0})
        return payload

    def test_cache_reuses_checksum_and_axis_when_histogram_name_changes(self):
        first_control = self._pion_control("pion_control_20260924")
        second_control = self._pion_control("pion_control_20260929")
        payload = self._cache_payload(first_control)
        self.assertEqual(
            payload["alignment_semantics_version"],
            fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION,
        )
        with tempfile.TemporaryDirectory() as directory:
            fits.persist_pion_component_alignment(
                directory, "Q2W-test", "center", "lowe", payload
            )
            original_resolver = fits.resolve_pion_component_alignment
            fits.resolve_pion_component_alignment = lambda *args, **kwargs: self.fail(
                "compatible cache hit must not resolve again"
            )
            try:
                cached, status, reasons, _ = fits.load_or_resolve_pion_component_alignment(
                    directory,
                    "Q2W-test",
                    "center",
                    "lowe",
                    "setting-wide",
                    self._identity(),
                    second_control,
                    {name: None for name in fits.COMPONENT_NAMES},
                    inp_dict={},
                    common_setting_shift_gev=0.005,
                )
            finally:
                fits.resolve_pion_component_alignment = original_resolver
        self.assertEqual(status, "reused")
        self.assertEqual(reasons, [])
        self.assertEqual(cached["pion_control_histogram_identifier"]["hist_name"], "pion_control_20260924")
        self.assertEqual(cached["alignment_semantics_version"], fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION)

    def test_missing_semantics_cache_resolves_once_and_then_reuses(self):
        control = self._pion_control("pion_control")
        old_payload = self._cache_payload(control)
        del old_payload["alignment_semantics_version"]
        calls = []
        with tempfile.TemporaryDirectory() as directory:
            fits.persist_pion_component_alignment(directory, "Q2W-test", "center", "lowe", old_payload)
            original_resolver = fits.resolve_pion_component_alignment

            def resolve_once(*args, **kwargs):
                calls.append(True)
                return self._cache_payload(control)

            fits.resolve_pion_component_alignment = resolve_once
            try:
                resolved, status, reasons, _ = fits.load_or_resolve_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", "setting-wide",
                    self._identity(), control, {name: None for name in fits.COMPONENT_NAMES},
                    inp_dict={}, common_setting_shift_gev=0.005,
                )
                reused, second_status, second_reasons, _ = fits.load_or_resolve_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", "setting-wide",
                    self._identity(), control, {name: None for name in fits.COMPONENT_NAMES},
                    inp_dict={}, common_setting_shift_gev=0.005,
                )
            finally:
                fits.resolve_pion_component_alignment = original_resolver
        self.assertEqual(calls, [True])
        self.assertEqual(status, "rejected_stale_then_created")
        self.assertIn("alignment_semantics_version mismatch", reasons)
        self.assertEqual(resolved["alignment_semantics_version"], fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION)
        self.assertEqual(second_status, "reused")
        self.assertEqual(second_reasons, [])
        self.assertEqual(reused["alignment_semantics_version"], fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION)

    def test_wrong_and_malformed_semantics_cache_records_fail_closed(self):
        control = self._pion_control("pion_control")
        expected = self._cache_payload(control)
        for value in ("pion_component_dynamic_alignment_semantics/v1", None, 2, {}):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as directory:
                stored = copy.deepcopy(expected)
                stored["alignment_semantics_version"] = value
                fits.persist_pion_component_alignment(directory, "Q2W-test", "center", "lowe", stored)
                cached, reasons, _ = fits.load_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", expected
                )
                self.assertIsNone(cached)
                self.assertIn("alignment_semantics_version mismatch", reasons)

    def test_stale_parent_semantics_fails_closed_before_fine_scan(self):
        config = get_pion_component_dynamic_alignment_config(
            inp_dict={}, phi_setting="center", setting_key="Q2W-test"
        )
        config_hash = fits._alignment_hash(config)
        parent = {
            "accepted": True,
            "alignment_schema_version": fits.ALIGNMENT_SCHEMA_VERSION,
            "resolved_configuration_hash": config_hash,
            "parent_setting_key": "Q2W-test",
            "components": {},
        }
        parent["parent_alignment_hash"] = fits._alignment_parent_hash(
            "Q2W-test", 0.005, parent["components"], config_hash
        )
        for version in (None, "pion_component_dynamic_alignment_semantics/v1"):
            with self.subTest(version=version):
                candidate = copy.deepcopy(parent)
                if version is not None:
                    candidate["alignment_semantics_version"] = version
                result = fits.resolve_pion_component_alignment(
                    "Q2W-test", "yield_t_phi", self._identity(),
                    self._pion_control("fine_control"),
                    {name: None for name in fits.COMPONENT_NAMES},
                    parent_alignment=candidate, inp_dict={}, phi_setting="center",
                    common_setting_shift_gev=0.005,
                )
                self.assertEqual(result["source"], "current_common_shift_fallback")
                self.assertEqual(result["alignment_semantics_version"], fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION)

    def test_pion_control_semantic_mismatches_fail_closed(self):
        payload = self._cache_payload(self._pion_control("pion_control_original"))
        with tempfile.TemporaryDirectory() as directory:
            fits.persist_pion_component_alignment(
                directory, "Q2W-test", "center", "lowe", payload
            )
            for expected in (
                self._cache_payload(self._pion_control("same_axis_new_content", contents=(12.0, 19.0, 7.0))),
                self._cache_payload(self._pion_control("changed_axis", edges=(0.70, 0.85, 1.10, 1.30))),
            ):
                cached, reasons, _ = fits.load_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", expected
                )
                self.assertIsNone(cached)
                self.assertIn("pion control histogram identifier mismatch", reasons)
            malformed = copy.deepcopy(payload)
            del malformed["pion_control_histogram_identifier"]["checksum"]
            cached, reasons, _ = fits.load_pion_component_alignment(
                directory, "Q2W-test", "center", "lowe", malformed
            )
        self.assertIsNone(cached)
        self.assertIn("pion control histogram identifier mismatch", reasons)

    def test_template_integral_minimum_accepts_only_machine_scale_roundoff(self):
        self.assertIn(
            "_template_integral_meets_minimum",
            inspect.getsource(fits.scan_pion_component_alignment),
        )
        self.assertTrue(fits._template_integral_meets_minimum(1.0, 1.0))
        self.assertTrue(fits._template_integral_meets_minimum(1.0 - 1.0e-13, 1.0))
        self.assertFalse(fits._template_integral_meets_minimum(1.0 - 1.0e-6, 1.0))


@unittest.skipUnless(ROOT is not None, "PyROOT is required for histogram alignment tests")
class DynamicPionAlignmentTests(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        self.config_override = {
            "enabled": True,
            "components": {
                "pi_n": {
                    "enabled": True,
                    "global_offset_scan_gev": {"minimum": -0.004, "maximum": 0.004, "step": 0.004},
                    "fine_bin_offset_scan_gev": {
                        "minimum_relative_to_parent": 0.0,
                        "maximum_relative_to_parent": 0.0,
                        "step": 0.001,
                    },
                    "window_expansion_candidates_gev": [0.0],
                },
                "pi_delta": {"enabled": False},
                "pi_sidis": {"enabled": False},
            },
            "acceptance": {
                "minimum_active_fit_bins": 2,
                "minimum_evaluation_bins": 2,
                "minimum_data_integral": 0.0,
                "minimum_template_integral": 0.0,
                "minimum_relative_score_improvement": 0.01,
                "maximum_lost_template_integral_fraction": 1.0,
                "reject_offset_boundary": False,
                "reject_expansion_boundary": False,
                "require_localized_minimum": False,
            },
        }
        self.inp_dict = {"pion_component_dynamic_alignment": self.config_override}
        self.source = self._histogram("immutable_pi_n")
        self.sources = {"pi_n": self.source, "pi_delta": None, "pi_sidis": None}

    def _histogram(self, name):
        histogram = ROOT.TH1D(name, name, 120, 0.70, 1.30)
        for index in range(1, histogram.GetNbinsX() + 1):
            x_value = histogram.GetBinCenter(index)
            value = 200.0 * math.exp(-0.5 * ((x_value - 0.91) / 0.018) ** 2)
            histogram.SetBinContent(index, value)
            histogram.SetBinError(index, max(math.sqrt(value), 1.0))
        histogram.SetDirectory(0)
        return histogram

    def _identity(self, scope="yield_t_phi"):
        return {
            "kinematic_setting": "Q2W-test",
            "epsilon": "lowe",
            "phi_setting": "center",
            "analysis_scope": scope,
            "t_bin": {"index": 0, "edges": [-0.20, -0.10]},
            "phi_bin": {"index": 0, "edges": [-10.0, 10.0]} if scope == "yield_t_phi" else None,
            "active_dimensions": {"particle_type": "kaon", "polarization": "unpolarized", "target": "LH2"},
        }

    def _valid_parent(self, residual=0.004, expansion=0.015, global_score=9.0):
        config = get_pion_component_dynamic_alignment_config(
            inp_dict=self.inp_dict, phi_setting="center", setting_key="Q2W-test"
        )
        config_hash = fits._alignment_hash(config)
        canonical = fits._alignment_canonical_windows(self.inp_dict, "center")
        limits = (0.70, 1.30)
        components = {}
        for name in fits.COMPONENT_NAMES:
            source = self.sources.get(name)
            component_residual = residual if name == "pi_n" else 0.0
            component_expansion = expansion if name == "pi_n" else 0.0
            total = 0.005 + component_residual
            component = {
                "accepted": True,
                "residual_shift_gev": component_residual,
                "total_shift_gev": total,
                "window_expansion_gev": component_expansion,
                "effective_windows": fits.transform_component_windows(
                    canonical.get(name), total, component_expansion, limits
                )["effective_windows"],
                "evaluation_envelope": fits.transform_component_windows(
                    canonical.get(name), total, component_expansion, limits
                )["effective_windows"],
            }
            component.update(fits._immutable_source_metadata(source, name))
            components[name] = component
        parent_hash = fits._alignment_parent_hash("Q2W-test", 0.005, components, config_hash)
        return {
            "alignment_schema_version": fits.ALIGNMENT_SCHEMA_VERSION,
            "alignment_semantics_version": fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION,
            "accepted": True,
            "parent_setting_key": "Q2W-test",
            "resolved_configuration_hash": config_hash,
            "parent_alignment_hash": parent_hash,
            "selected_score": global_score,
            "components": components,
        }

    def test_fine_baseline_replays_parent_expansion_and_rescores(self):
        parent = self._valid_parent()
        result = fits.resolve_pion_component_alignment(
            "Q2W-test", "yield_t_phi", self._identity(), self._histogram("fine_pion_control"),
            self.sources, parent_alignment=parent, inp_dict=self.inp_dict,
            phi_setting="center", common_setting_shift_gev=0.005,
        )
        component = result["components"]["pi_n"]
        parent_component = parent["components"]["pi_n"]
        self.assertEqual(result["source"], "parent_setting_fallback")
        self.assertEqual(component["source"], "parent_setting_fallback")
        self.assertTrue(component["parent_baseline_candidate_included"])
        self.assertAlmostEqual(component["residual_shift_gev"], parent_component["residual_shift_gev"])
        self.assertAlmostEqual(component["total_shift_gev"], parent_component["total_shift_gev"])
        self.assertAlmostEqual(component["window_expansion_gev"], parent_component["window_expansion_gev"])
        self.assertEqual(component["effective_windows"], parent_component["effective_windows"])
        self.assertEqual(result["parent_global_score"], 9.0)
        self.assertIsNotNone(result["parent_score_in_current_bin"])
        self.assertNotEqual(result["parent_score_in_current_bin"], result["parent_global_score"])

    def test_missing_parent_uses_common_shift_without_a_local_scan(self):
        result = fits.resolve_pion_component_alignment(
            "Q2W-test", "yield_t_phi", self._identity(), self._histogram("fine_no_parent"),
            self.sources, parent_alignment=None, inp_dict=self.inp_dict,
            phi_setting="center", common_setting_shift_gev=0.005,
        )
        component = result["components"]["pi_n"]
        self.assertEqual(result["source"], "current_common_shift_fallback")
        self.assertEqual(component["source"], "current_common_shift_fallback")
        self.assertNotIn("candidate_summaries", component)

    def test_stale_parent_semantics_uses_common_shift_fallback(self):
        for version in (None, "pion_component_dynamic_alignment_semantics/v1"):
            with self.subTest(version=version):
                parent = self._valid_parent()
                if version is None:
                    del parent["alignment_semantics_version"]
                else:
                    parent["alignment_semantics_version"] = version
                result = fits.resolve_pion_component_alignment(
                    "Q2W-test", "yield_t_phi", self._identity(),
                    self._histogram("fine_stale_parent"), self.sources,
                    parent_alignment=parent, inp_dict=self.inp_dict,
                    phi_setting="center", common_setting_shift_gev=0.005,
                )
                self.assertEqual(result["source"], "current_common_shift_fallback")
                self.assertEqual(result["alignment_semantics_version"], fits.PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION)

    def test_local_correction_is_accepted_only_when_it_beats_the_parent_score(self):
        self.config_override["components"]["pi_n"]["fine_bin_offset_scan_gev"] = {
            "minimum_relative_to_parent": -0.004,
            "maximum_relative_to_parent": 0.004,
            "step": 0.004,
        }
        parent = self._valid_parent(residual=0.0, expansion=0.0)
        target, _ = fits.build_shifted_template_histogram(
            self.source, 0.009, "positive_moves_peak_higher_mm", "fine_local_target"
        )
        result = fits.resolve_pion_component_alignment(
            "Q2W-test", "yield_t_phi", self._identity(), target,
            self.sources, parent_alignment=parent, inp_dict=self.inp_dict,
            phi_setting="center", common_setting_shift_gev=0.005,
        )
        component = result["components"]["pi_n"]
        self.assertTrue(result["accepted"])
        self.assertEqual(result["source"], "bin_local_scan")
        self.assertAlmostEqual(component["local_bin_correction_gev"], 0.004)
        self.assertGreater(component["relative_improvement_over_parent"], 0.01)

    def test_each_alignment_build_uses_the_raw_immutable_template_once(self):
        calls = []
        original_builder = fits.build_shifted_template_histogram

        def recording_builder(source_hist, *args, **kwargs):
            calls.append(source_hist)
            return original_builder(source_hist, *args, **kwargs)

        fits.build_shifted_template_histogram = recording_builder
        try:
            result = fits.resolve_pion_component_alignment(
                "Q2W-test", "setting-wide", self._identity("setting-wide"), self._histogram("parent_control"),
                self.sources, inp_dict=self.inp_dict, phi_setting="center", common_setting_shift_gev=0.005,
            )
        finally:
            fits.build_shifted_template_histogram = original_builder
        self.assertTrue(calls)
        self.assertTrue(all(call is self.source for call in calls))
        provenance = result["components"]["pi_n"]
        self.assertEqual(provenance["immutable_source_hist_name"], "immutable_pi_n")
        self.assertTrue(provenance["candidate_built_from_immutable_source"])
        self.assertTrue(provenance["applied_in_single_shift_operation"])
        self.assertEqual(provenance["shift_operation_count"], 1)

    def test_boundary_penalties_follow_their_independent_rejection_policies(self):
        config = get_pion_component_dynamic_alignment_config(
            inp_dict=self.inp_dict, phi_setting="center", setting_key="Q2W-test"
        )
        config["components"]["pi_n"]["global_offset_scan_gev"] = {
            "minimum": -0.004, "maximum": 0.004, "step": 0.004,
        }
        config["components"]["pi_n"]["window_expansion_candidates_gev"] = [0.0, 0.010]
        config["acceptance"]["reject_offset_boundary"] = False
        config["acceptance"]["reject_expansion_boundary"] = False
        scan = fits.scan_pion_component_alignment(
            self._histogram("boundary_target"), self._histogram("boundary_residual"),
            self.source, "pi_n", fits._alignment_canonical_windows(self.inp_dict, "center")["pi_n"],
            0.005, config, "setting-wide", self._identity("setting-wide"), scan_level="parent",
        )
        expansion_edge_candidates = [
            candidate for candidate in scan["candidate_summaries"]
            if candidate["expansion_boundary_hit"] and not candidate["is_baseline"]
        ]
        self.assertTrue(expansion_edge_candidates)
        self.assertTrue(all(candidate["expansion_boundary_penalty_applied"] == 0.0 for candidate in expansion_edge_candidates))

        config["acceptance"]["reject_offset_boundary"] = True
        config["ranking"]["offset_boundary_penalty"] = 123.0
        scan = fits.scan_pion_component_alignment(
            self._histogram("offset_boundary_target"), self._histogram("offset_boundary_residual"),
            self.source, "pi_n", fits._alignment_canonical_windows(self.inp_dict, "center")["pi_n"],
            0.005, config, "setting-wide", self._identity("setting-wide"), scan_level="parent",
        )
        offset_edge_candidates = [
            candidate for candidate in scan["candidate_summaries"]
            if candidate["offset_boundary_hit"] and not candidate["is_baseline"]
        ]
        self.assertTrue(offset_edge_candidates)
        self.assertTrue(all(candidate["offset_boundary_penalty_applied"] == 123.0 for candidate in offset_edge_candidates))

    def test_support_metrics_do_not_cancel_positive_and_negative_residuals(self):
        histogram = self._histogram("support_metrics")
        histogram.SetBinContent(10, 8.0)
        histogram.SetBinContent(11, -8.0)
        metrics = fits._alignment_support_metrics(histogram, [10, 11], "positive_integral")
        self.assertEqual(metrics["signed_support_integral"], 0.0)
        self.assertEqual(metrics["absolute_support_integral"], 16.0)
        self.assertEqual(metrics["positive_support_integral"], 8.0)
        self.assertEqual(metrics["support_integral_for_acceptance"], 8.0)

    def test_accepted_pi_n_can_be_applied_with_pi_delta_common_shift_fallback(self):
        self.config_override["components"]["pi_delta"]["enabled"] = True
        sources = dict(self.sources)
        sources["pi_delta"] = self._histogram("immutable_pi_delta")
        canonical = fits._alignment_canonical_windows(self.inp_dict, "center")
        original_scan = fits.scan_pion_component_alignment
        original_score = fits._score_complete_alignment

        def fake_scan(*args, **kwargs):
            component_name = args[3]
            effective_windows = fits.transform_component_windows(
                canonical[component_name], 0.005, 0.0, (0.70, 1.30)
            )["effective_windows"]
            candidate = {
                "residual_shift_gev": 0.0, "total_shift_gev": 0.005,
                "window_expansion_gev": 0.0, "local_bin_correction_gev": 0.0,
                "effective_windows": effective_windows, "evaluation_envelope": effective_windows,
                "score": 1.0, "rejection_reasons": ([] if component_name == "pi_n" else ["insufficient positive support"]),
                "warnings": [], "offset_boundary_hit": False, "expansion_boundary_hit": False,
            }
            return {
                "accepted": component_name == "pi_n", "source": "setting_global_scan",
                "baseline_source": "current_common_shift", "baseline_score": 2.0,
                "baseline_residual_shift_gev": 0.0, "baseline_total_shift_gev": 0.005,
                "baseline_window_expansion_gev": 0.0, "canonical_windows": canonical[component_name],
                "evaluation_envelope": effective_windows, "selected_candidate": candidate,
                "candidate_summaries": [candidate], "candidate_count": 1,
                "near_minimum_candidate_count": 1, "near_minimum_candidate_fraction": 1.0,
                "minimum_localized": True,
            }

        def fake_score(_target, _sources, components, _config, _windows=None):
            pi_n_source = (components.get("pi_n") or {}).get("source")
            pi_delta_source = (components.get("pi_delta") or {}).get("source")
            if pi_n_source == "current_common_shift_fallback":
                score = 10.0
            elif pi_delta_source == "current_common_shift_fallback":
                score = 5.0
            else:
                score = 6.0
            return {"score": score, "chi2": score, "ndf": 1, "evaluation_bin_indices": [1]}

        fits.scan_pion_component_alignment = fake_scan
        fits._score_complete_alignment = fake_score
        try:
            result = fits.resolve_pion_component_alignment(
                "Q2W-test", "setting-wide", self._identity("setting-wide"),
                self._histogram("mixed_map_target"), sources, inp_dict=self.inp_dict,
                phi_setting="center", common_setting_shift_gev=0.005,
            )
        finally:
            fits.scan_pion_component_alignment = original_scan
            fits._score_complete_alignment = original_score
        self.assertTrue(result["accepted"])
        self.assertEqual(result["components"]["pi_n"]["component_applied_source"], "setting_global_scan")
        self.assertEqual(result["components"]["pi_delta"]["component_applied_source"], "current_common_shift_fallback")
        self.assertTrue(result["components"]["pi_delta"]["component_fallback_used"])
        self.assertEqual(result["component_counts"]["globally_accepted"], 1)
        self.assertEqual(result["component_counts"]["common_shift_fallback"], 1)

    def test_compatible_cache_is_loaded_without_a_second_scan(self):
        with tempfile.TemporaryDirectory() as directory:
            first, first_status, _, _ = fits.load_or_resolve_pion_component_alignment(
                directory, "Q2W-test", "center", "lowe", "yield_t_phi", self._identity(),
                self._histogram("cached_control"), self.sources, parent_alignment=None,
                inp_dict=self.inp_dict, common_setting_shift_gev=0.005,
            )
            self.assertEqual(first_status, "created")
            original_resolve = fits.resolve_pion_component_alignment
            fits.resolve_pion_component_alignment = lambda *args, **kwargs: self.fail("cache hit must not scan")
            try:
                second, second_status, reasons, _ = fits.load_or_resolve_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", "yield_t_phi", self._identity(),
                    self._histogram("cached_control"), self.sources, parent_alignment=None,
                    inp_dict=self.inp_dict, common_setting_shift_gev=0.005,
                )
            finally:
                fits.resolve_pion_component_alignment = original_resolve
            self.assertEqual(second_status, "reused")
            self.assertEqual(reasons, [])
            self.assertEqual(second["selected_score"], first["selected_score"])

    def test_kaon_side_data_cannot_enter_the_alignment_api(self):
        parameter_names = set(inspect.signature(fits.resolve_pion_component_alignment).parameters)
        self.assertFalse(any("kaon" in name.lower() or name.lower().startswith("a_") for name in parameter_names))
        first = fits.resolve_pion_component_alignment(
            "Q2W-test", "setting-wide", self._identity("setting-wide"), self._histogram("pion_control_a"),
            self.sources, inp_dict=self.inp_dict, phi_setting="center", common_setting_shift_gev=0.005,
        )
        second = fits.resolve_pion_component_alignment(
            "Q2W-test", "setting-wide", self._identity("setting-wide"), self._histogram("pion_control_b"),
            self.sources, inp_dict=self.inp_dict, phi_setting="center", common_setting_shift_gev=0.005,
        )
        self.assertEqual(first["components"], second["components"])

    def test_persistence_rejects_parent_config_source_axis_and_identity_mismatches(self):
        payload = fits.resolve_pion_component_alignment(
            "Q2W-test", "yield_t_phi", self._identity(), self._histogram("persist_control"),
            self.sources, parent_alignment=None, inp_dict=self.inp_dict,
            phi_setting="center", common_setting_shift_gev=0.005,
        )
        with tempfile.TemporaryDirectory() as directory:
            fits.persist_pion_component_alignment(directory, "Q2W-test", "center", "lowe", payload)
            loaded, reasons, _ = fits.load_pion_component_alignment(directory, "Q2W-test", "center", "lowe", payload)
            self.assertIsNotNone(loaded)
            self.assertEqual(reasons, [])
            for field, value in (
                ("parent_alignment_hash", "changed-parent"),
                ("resolved_configuration_hash", "changed-config"),
                ("histogram_axis_specification", {"nbins": 1}),
                ("complete_physical_bin_identity", self._identity("average_t_phi_integrated")),
            ):
                expected = copy.deepcopy(payload)
                expected[field] = value
                rejected, reasons, _ = fits.load_pion_component_alignment(
                    directory, "Q2W-test", "center", "lowe", expected
                )
                self.assertIsNone(rejected)
                self.assertTrue(reasons)
            expected = copy.deepcopy(payload)
            expected["components"]["pi_n"]["immutable_source_hist_checksum"] = "changed-source"
            rejected, reasons, _ = fits.load_pion_component_alignment(
                directory, "Q2W-test", "center", "lowe", expected
            )
            self.assertIsNone(rejected)
            self.assertTrue(reasons)


if __name__ == "__main__":
    unittest.main()
