"""Regression coverage for authoritative per-t pion-parent integrity."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
PION_SUBTRACTION_PATH = REPO_ROOT / "src" / "cuts" / "pion_component_subtraction.py"
CALCULATE_YIELD_PATH = REPO_ROOT / "src" / "binning" / "calculate_yield.py"
CANONICAL_BINNING_PATH = REPO_ROOT / "src" / "utility" / "canonical_binning.py"


class _FakeHistogram:
    def __init__(self):
        self.entries = []

    def Fill(self, value, weight=1.0):
        self.entries.append((float(value), float(weight)))

    @property
    def integral(self):
        return sum(weight for _, weight in self.entries)


class _FakeRootHistogram:
    """Small TH1-compatible test double for structural proposal validation."""

    class _Axis:
        def __init__(self, low, high):
            self._low = float(low)
            self._high = float(high)

        def GetXmin(self):
            return self._low

        def GetXmax(self):
            return self._high

    def __init__(self, contents, errors=None, low=1.0, high=1.3):
        self._contents = [float(value) for value in contents]
        self._errors = [float(value) for value in (errors or [0.0] * len(contents))]
        self._axis = self._Axis(low, high)

    def InheritsFrom(self, name):
        return name in ("TH1", "TObject")

    def GetNbinsX(self):
        return len(self._contents)

    def GetXaxis(self):
        return self._axis

    def Integral(self):
        return sum(self._contents)

    def GetBinContent(self, index):
        if 1 <= int(index) <= len(self._contents):
            return self._contents[int(index) - 1]
        return 0.0

    def GetBinError(self, index):
        if 1 <= int(index) <= len(self._errors):
            return self._errors[int(index) - 1]
        return 0.0


def _load_pion_component_subtraction():
    fake_background_config = types.ModuleType("background_config")
    fake_background_config.PARTICLE_SUBTRACTION_MODE_COMPONENTS = "simc_shape_components"
    fake_background_config.resolve_particle_subtraction_mode = lambda *_args, **_kwargs: "simc_shape_components"
    fake_background_config.resolve_particle_subtraction_fallback_mode = lambda inp=None, **_kwargs: str((inp or {}).get("particle_subtraction_fallback_mode", "error"))
    fake_background_config.resolve_particle_subtraction_weight_clip_bounds = lambda *_args, **_kwargs: (0.0, None)
    fake_background_config.resolve_particle_subtraction_weight_denominator_floor = lambda inp=None, **_kwargs: float((inp or {}).get("particle_subtraction_weight_denominator_floor", 1.0e-12))
    fake_background_config.get_pion_component_dynamic_alignment_config = lambda *_args, **_kwargs: {}
    for name, value in {
        "resolve_particle_subtraction_component_fit_windows": {},
        "resolve_particle_subtraction_component_fit_excluded_windows": [],
        "resolve_particle_subtraction_component_stage_amplitude_windows": {},
        "resolve_particle_subtraction_component_stage_amplitude_modes": {},
        "resolve_particle_subtraction_component_prior_scales": {},
        "resolve_particle_subtraction_component_fit_mode": "staged_plus_joint",
        "resolve_particle_subtraction_component_postfit_scales": {},
        "resolve_particle_subtraction_component_postrefine_scales": {},
        "resolve_particle_subtraction_component_residual_shift_settings": {},
        "resolve_particle_subtraction_component_cleanup_validation_mm_max": None,
    }.items():
        setattr(fake_background_config, name, lambda *_args, _value=value, **_kwargs: _value)
    fake_background_config.get_proton_contamination_cleaning_config = lambda **_kwargs: {
        "t_binning": {"edge_tolerance": 1.0e-9}
    }
    fake_mm_background = types.ModuleType("mm_background_subtraction")
    fake_mm_background.mm_background_weight_from_value = lambda *_args, **_kwargs: 1.0
    fake_ownership = types.ModuleType("root_histogram_ownership")
    fake_ownership.clone_root_histogram = lambda histogram, **_kwargs: histogram
    def _fingerprint(histogram):
        values = [
            float(histogram.GetNbinsX()),
            float(histogram.GetXaxis().GetXmin()),
            float(histogram.GetXaxis().GetXmax()),
        ]
        for index in range(0, int(histogram.GetNbinsX()) + 2):
            values.extend((float(histogram.GetBinContent(index)), float(histogram.GetBinError(index))))
        return hashlib.sha256(repr(values).encode("utf-8")).hexdigest()
    fake_ownership.fingerprint_histogram_content_error = _fingerprint

    with mock.patch.dict(
        sys.modules,
        {
            "background_config": fake_background_config,
            "mm_background_subtraction": fake_mm_background,
            "root_histogram_ownership": fake_ownership,
        },
    ):
        spec = importlib.util.spec_from_file_location(
            "pion_component_subtraction_parent_integrity", PION_SUBTRACTION_PATH
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def _inp():
    return {
        "Q2": "4p4",
        "W": "2p74",
        "EPSSET": "low",
        "pion_subtraction_scope": "t_bin",
        "require_shared_canonical_preflight": True,
        "canonical_t_binning": {
            "canonical_interval_pair_id": "pair-q44w274",
            "canonical_interval_pair_hash": "pair-hash-q44w274",
        },
    }


def _configured_inp():
    inp_dict = _inp()
    inp_dict.update({
        "analysis_runtime_config_hash": "analysis-config-hash",
        "analysis_runtime_config": {"config_hash": "analysis-config-hash"},
    })
    inp_dict["canonical_t_binning"] = {
        **inp_dict["canonical_t_binning"],
        "t_edges": [0.0, 0.4, 0.8],
        "phi_edges": [-180.0, 0.0, 180.0],
        "requested_num_t_bins": 2,
        "actual_num_t_bins": 2,
        "requested_num_phi_bins": 2,
        "actual_num_phi_bins": 2,
    }
    return inp_dict


def _load_canonical_binning():
    spec = importlib.util.spec_from_file_location("canonical_binning_test", CANONICAL_BINNING_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _parent(module, inp_dict, t_index=0, edges=(0.0, 0.4)):
    identity = module.build_t_bin_pion_parent_identity(
        inp_dict, "Left", t_index, edges
    )
    return {
        **identity,
        "analysis_scope": "t_bin{}".format(t_index + 1),
        "fit_result": {"analysis_scope": "t_bin{}".format(t_index + 1)},
    }


class TBinPionParentIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = _load_pion_component_subtraction()

    def test_parent_identity_is_deterministic_and_t_scoped(self):
        first = self.module.build_t_bin_pion_parent_identity(
            _inp(), "Left", 0, (0.0, 0.4)
        )
        repeat = self.module.build_t_bin_pion_parent_identity(
            _inp(), "Left", 0, (0.0, 0.4)
        )
        second_t = self.module.build_t_bin_pion_parent_identity(
            _inp(), "Left", 1, (0.4, 0.8)
        )
        self.assertEqual(first["pion_parent_id"], repeat["pion_parent_id"])
        self.assertNotEqual(first["pion_parent_id"], second_t["pion_parent_id"])

    def test_parent_identity_separates_parent_physics_from_phi_run_provenance(self):
        configured = _configured_inp()
        baseline = self.module.build_t_bin_pion_parent_identity(
            configured, "Left", 0, (0.0, 0.4)
        )
        changed_config = _configured_inp()
        changed_config["analysis_runtime_config_hash"] = "other-config-hash"
        changed_config["analysis_runtime_config"]["config_hash"] = "other-config-hash"
        changed_pair = _configured_inp()
        changed_pair["canonical_t_binning"]["canonical_interval_pair_hash"] = "other-pair"
        changed_phi = _configured_inp()
        changed_phi["canonical_t_binning"]["phi_edges"] = [-180.0, -30.0, 180.0]

        self.assertEqual(
            baseline["pion_parent_id"],
            self.module.build_t_bin_pion_parent_identity(changed_config, "Left", 0, (0.0, 0.4))["pion_parent_id"],
        )
        self.assertEqual(
            baseline["pion_parent_id"],
            self.module.build_t_bin_pion_parent_identity(changed_pair, "Left", 0, (0.0, 0.4))["pion_parent_id"],
        )
        self.assertEqual(
            baseline["pion_parent_id"],
            self.module.build_t_bin_pion_parent_identity(changed_phi, "Left", 0, (0.0, 0.4))["pion_parent_id"],
        )
        changed_fit = _configured_inp()
        changed_fit["particle_subtraction_weight_denominator_floor"] = 1.0e-6
        self.assertNotEqual(
            baseline["pion_parent_id"],
            self.module.build_t_bin_pion_parent_identity(changed_fit, "Left", 0, (0.0, 0.4))["pion_parent_id"],
        )

    def test_histogram_fingerprint_includes_axis_contents_and_errors(self):
        baseline = _FakeRootHistogram((1.0, 2.0, 3.0), errors=(0.1, 0.2, 0.3))
        same = _FakeRootHistogram((1.0, 2.0, 3.0), errors=(0.1, 0.2, 0.3))
        changed_content = _FakeRootHistogram((1.0, 2.1, 3.0), errors=(0.1, 0.2, 0.3))
        changed_error = _FakeRootHistogram((1.0, 2.0, 3.0), errors=(0.1, 0.25, 0.3))
        changed_axis = _FakeRootHistogram((1.0, 2.0, 3.0), errors=(0.1, 0.2, 0.3), high=1.4)

        baseline_fingerprint = self.module.fingerprint_histogram_content_error(baseline)
        self.assertEqual(
            baseline_fingerprint,
            self.module.fingerprint_histogram_content_error(same),
        )
        for histogram in (changed_content, changed_error, changed_axis):
            self.assertNotEqual(
                baseline_fingerprint,
                self.module.fingerprint_histogram_content_error(histogram),
            )

    def test_parent_validation_rejects_wrong_identity_or_edges(self):
        parent = _parent(self.module, _inp())
        resolved = self.module.validate_authoritative_t_bin_pion_parent(
            parent, _inp(), "Left", 0, (0.0, 0.4)
        )
        self.assertEqual(resolved["analysis_scope"], "t_bin1")

        wrong_setting = dict(parent, phi_setting="Center")
        with self.assertRaisesRegex(RuntimeError, "phi_setting"):
            self.module.validate_authoritative_t_bin_pion_parent(
                wrong_setting, _inp(), "Left", 0, (0.0, 0.4)
            )

        wrong_epsilon = dict(parent, epsilon="high")
        with self.assertRaisesRegex(RuntimeError, "epsilon"):
            self.module.validate_authoritative_t_bin_pion_parent(
                wrong_epsilon, _inp(), "Left", 0, (0.0, 0.4)
            )

        wrong_pair = dict(parent, canonical_interval_pair_hash="other-pair")
        self.module.validate_authoritative_t_bin_pion_parent(
            wrong_pair, _inp(), "Left", 0, (0.0, 0.4)
        )
        self.assertNotIn("provenance_validation", wrong_pair)

        with self.assertRaisesRegex(RuntimeError, "missing_authoritative"):
            self.module.validate_authoritative_t_bin_pion_parent(
                None, _inp(), "Left", 0, (0.0, 0.4)
            )

        with self.assertRaisesRegex(RuntimeError, "canonical t-edge mismatch"):
            self.module.validate_authoritative_t_bin_pion_parent(
                parent, _inp(), "Left", 0, (0.0, 0.5)
            )

    def test_parent_validation_rejects_active_coordinate_mismatch(self):
        inp_dict = _configured_inp()
        inp_dict["kaon_data_coordinate_summary"] = {
            "Left": {"coordinate_fingerprint": "kaon-coordinate-a"}
        }
        parent = _parent(self.module, inp_dict)
        parent["coordinate_fingerprint"] = "kaon-coordinate-a"
        parent["kaon_data_coordinate"] = {"coordinate_fingerprint": "kaon-coordinate-a"}
        self.module.validate_authoritative_t_bin_pion_parent(
            parent, inp_dict, "Left", 0, (0.0, 0.4)
        )

        parent["coordinate_fingerprint"] = "kaon-coordinate-b"
        with self.assertRaisesRegex(RuntimeError, "coordinate_fingerprint"):
            self.module.validate_authoritative_t_bin_pion_parent(
                parent, inp_dict, "Left", 0, (0.0, 0.4)
            )

    def test_frozen_collection_requires_exact_canonical_parent_set(self):
        first = _parent(self.module, _inp(), 0, (0.0, 0.4))
        second = _parent(self.module, _inp(), 1, (0.4, 0.8))
        hist = {
            "phi_setting": "Left",
            "_pion_t_bin_parent_results": [first, second],
        }
        resolved = self.module.validate_frozen_t_bin_pion_parent_collection(
            hist, _inp(), (0.0, 0.4, 0.8)
        )
        self.assertIs(resolved[0], first)
        with self.assertRaisesRegex(RuntimeError, "collection:count"):
            self.module.validate_frozen_t_bin_pion_parent_collection(
                {**hist, "_pion_t_bin_parent_results": [first]},
                _inp(),
                (0.0, 0.4, 0.8),
            )

    def test_proposal_validation_is_structural_not_the_production_quality_gate(self):
        template = _FakeRootHistogram((1.0, 2.0, 3.0))
        fit_result = {
            "particle_subtraction_mode": "simc_shape_components",
            "fit_status_pion": "quality_rejected",
            "fit_status_kaon": "quality_rejected",
            "fallback_used": False,
            "A_n": 0.4,
            "A_delta": 0.2,
            "A_sidis": 0.1,
            "B_n": 0.0,
            "B_delta": 0.0,
            "B_sidis": 0.0,
            "H_simc_shape_pi_n": template,
            "H_simc_shape_pi_delta": template,
            "H_simc_shape_pi_sidis": template,
            "diagnostics": {
                "pion": {"validation": {"accepted": False}},
                "kaon": {"validation": {"accepted": False}},
            },
        }
        production_gate = self.module.evaluate_particle_subtraction_component_fit_result(
            fit_result, {"particle_subtraction_fallback_mode": "zero"}
        )
        proposal_gate = self.module.evaluate_component_pion_application_proposal(
            fit_result
        )
        self.assertFalse(production_gate["accepted"])
        self.assertTrue(proposal_gate["available"])
        self.assertEqual(
            proposal_gate["diagnostics"]["active_component_names"],
            ["pi_n", "pi_delta", "pi_sidis"],
        )

    def test_event_level_children_keep_their_own_control_support(self):
        weights = {1.11: 2.0, 1.13: 3.0, 1.15: 5.0}
        self.module.simc_shape_pion_weight_from_value = (
            lambda adj_mm, _reference, _weights: weights[round(float(adj_mm), 2)]
        )
        cache = {
            "allcuts_bin_index": {
                (0, 0): np.asarray([0, 1], dtype=np.int32),
                (0, 1): np.asarray([2], dtype=np.int32),
            },
            "nommcuts_bin_index": {
                (0, 0): np.asarray([0, 1], dtype=np.int32),
                (0, 1): np.asarray([2], dtype=np.int32),
            },
            "adj_MM": np.asarray([1.11, 1.13, 1.15]),
            "adj_t": np.asarray([0.10, 0.20, 0.30]),
            "Q2": np.asarray([4.4, 4.4, 4.4]),
            "W": np.asarray([2.74, 2.74, 2.74]),
            "epsilon": np.asarray([0.45, 0.45, 0.45]),
            "theta_cm_deg": np.asarray([20.0, 20.0, 20.0]),
            "ssxptar": np.asarray([0.0, 0.0, 0.0]),
            "ssyptar": np.asarray([0.0, 0.0, 0.0]),
            "hsxptar": np.asarray([0.0, 0.0, 0.0]),
            "hsyptar": np.asarray([0.0, 0.0, 0.0]),
        }
        specs = [{"cache_section": cache, "coefficient": 1.0}]

        first = {"mm": _FakeHistogram(), "mm_nosub": _FakeHistogram()}
        second = {"mm": _FakeHistogram(), "mm_nosub": _FakeHistogram()}
        first_stats = self.module.fill_simc_shape_pion_subtraction_templates(
            first, specs, object(), object(), {"t_index": 0, "phi_index": 0}, "kaon", "0"
        )
        second_stats = self.module.fill_simc_shape_pion_subtraction_templates(
            second, specs, object(), object(), {"t_index": 0, "phi_index": 1}, "kaon", "0"
        )

        self.assertEqual(first_stats["n_events_allcuts"], 2)
        self.assertEqual(second_stats["n_events_allcuts"], 1)
        self.assertAlmostEqual(first["mm"].integral, 5.0)
        self.assertAlmostEqual(second["mm"].integral, 5.0)
        self.assertNotEqual(first["mm"].entries, second["mm"].entries)
        self.assertAlmostEqual(12.0 - first["mm"].integral, 7.0)

    def test_authoritative_cache_coefficients_are_used_per_record(self):
        source_spec = {
            "cache_section": {"coefficient": np.asarray([2.0, -3.0])},
            "coefficient": -5.0,
            "base_coefficient": 5.0,
        }
        self.assertAlmostEqual(
            self.module._component_cache_event_coefficient(source_spec, 0), -2.0
        )
        self.assertAlmostEqual(
            self.module._component_cache_event_coefficient(source_spec, 1), 3.0
        )

    def test_frozen_parent_policy_keeps_zero_and_skip_bin_distinct(self):
        template = _FakeRootHistogram((1.0, 2.0, 3.0))
        fit = {
            "particle_subtraction_mode": "simc_shape_components",
            "fit_status_pion": "rejected",
            "fit_status_kaon": "rejected",
            "fallback_used": False,
            "A_n": 0.0, "A_delta": 0.0, "A_sidis": 0.0,
            "B_n": 0.0, "B_delta": 0.0, "B_sidis": 0.0,
            "H_simc_shape_pi_n": template,
            "H_simc_shape_pi_delta": template,
            "H_simc_shape_pi_sidis": template,
            "diagnostics": {"pion": {"validation": {"accepted": False}}, "kaon": {"validation": {"accepted": False}}},
        }
        parent = {"pion_parent_id": "parent", "fit_result": fit}
        zero = self.module.resolve_frozen_parent_application_policy(
            parent, {"particle_subtraction_fallback_mode": "zero"}
        )
        skipped = self.module.resolve_frozen_parent_application_policy(
            parent, {"particle_subtraction_fallback_mode": "skip_bin"}
        )
        self.assertEqual(zero["action"], "zero")
        self.assertTrue(zero["child_valid"])
        self.assertEqual(skipped["action"], "skip_bin")
        self.assertFalse(skipped["child_valid"])


class CanonicalBinningTests(unittest.TestCase):
    def test_final_endpoint_belongs_to_final_canonical_bin(self):
        module = _load_canonical_binning()
        edges = (0.0, 0.4, 0.8)
        self.assertEqual(module.find_canonical_bin(0.0, edges), 0)
        self.assertEqual(module.find_canonical_bin(0.4, edges), 1)
        self.assertEqual(module.find_canonical_bin(0.8, edges), 1)
        self.assertEqual(module.find_canonical_bin(0.800001, edges), -1)


class CalculateYieldPayloadScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = CALCULATE_YIELD_PATH.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    @classmethod
    def _function_source(cls, name):
        for node in ast.walk(cls.tree):
            if isinstance(node, ast.FunctionDef) and node.name == name:
                lines = cls.source.splitlines()
                return "\n".join(lines[node.lineno - 1 : node.end_lineno])
        raise AssertionError("function not found: {}".format(name))

    def test_direct_path_has_one_local_payload_and_no_stale_entry(self):
        source = self._function_source("process_hist_data")
        self.assertIn("active_component_payload = _accepted_component_payload(", source)
        self.assertEqual(source.count("active_component_payload = _accepted_component_payload("), 1)
        self.assertNotIn('entry.get("particle_subtraction_component_payload")', source)
        self.assertIn("free_t_phi_pion_fit_forbidden_in_t_bin_scope", source)

    def test_cached_path_resolves_payload_before_fit1_and_reuses_it_for_fit2(self):
        source = self._function_source("_process_hist_data_from_base_cache")
        assignment = source.index("active_component_payload = _accepted_component_payload(")
        fit1 = source.index("fit_name=\"Fit 1\"")
        fit2 = source.index("fit_name=\"Fit 2\"")
        self.assertEqual(source.count("active_component_payload = _accepted_component_payload("), 1)
        self.assertLess(assignment, fit1)
        self.assertIn("pion_component_payload=active_component_payload", source[fit1:fit2])
        self.assertIn("pion_component_payload=active_component_payload", source[fit2:])

    def test_child_application_is_pion_only_and_subtracts_once(self):
        source = self._function_source("_apply_component_pion_subtraction_for_bin")
        self.assertNotIn("k_lambda", source.lower())
        self.assertNotIn("k_sigma0", source.lower())
        subtraction = (
            'hist_bin_dict["H_MM_DATA_{}_{}".format(j, k)].Add(template_hists["mm"], -1.0)'
        )
        self.assertEqual(source.count(subtraction), 1)

    def test_bg_optimizer_cache_is_explicit_single_scale_prepass(self):
        source = self._function_source("prepare_bg_opt_data_base_cache")
        config_builder = source.index("build_bg_optimization_prepass_config(inpDict)")
        process_call = source.index("processed_dict, _, ave_event_cache, sub_event_cache = process_hist_data(")
        self.assertLess(config_builder, process_call)
        self.assertIn('"_cache_contract": cache_contract', source)


if __name__ == "__main__":
    unittest.main()
