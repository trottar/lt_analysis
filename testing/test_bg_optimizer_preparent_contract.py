"""Regression coverage for the non-authoritative BG optimizer prepass."""

from __future__ import annotations

from copy import deepcopy
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
UTILITY_DIR = REPO_ROOT / "src" / "utility"
BINNING_DIR = REPO_ROOT / "src" / "binning"
CALCULATE_YIELD_PATH = BINNING_DIR / "calculate_yield.py"
BG_OPTIMIZATION_PATH = UTILITY_DIR / "bg_optimization.py"
if str(UTILITY_DIR) not in sys.path:
    sys.path.insert(0, str(UTILITY_DIR))

from background_config import (  # noqa: E402
    BG_OPT_PREPASS_CONTEXT,
    build_bg_optimization_prepass_config,
)


class _FallbackModule(types.ModuleType):
    def __getattr__(self, _name):
        value = lambda *_args, **_kwargs: None
        setattr(self, _name, value)
        return value


class _FakeRootContext:
    USER = "test"
    HOST = "test"
    REPLAYPATH = ""
    UTILPATH = ""
    LTANAPATH = ""
    ANATYPE = ""
    OUTPATH = ""


def _module_with_fallbacks(name):
    return _FallbackModule(name)


def _load_calculate_yield_module():
    fake_root = _module_with_fallbacks("ROOT")
    fake_ltsep = _module_with_fallbacks("ltsep")
    fake_ltsep.Root = lambda *_args, **_kwargs: _FakeRootContext()
    fake_ltsep.Misc = object()
    fake_modules = {
        "uproot": _module_with_fallbacks("uproot"),
        "root_numpy": _module_with_fallbacks("root_numpy"),
        "ROOT": fake_root,
        "scipy": _module_with_fallbacks("scipy"),
        "scipy.integrate": _module_with_fallbacks("scipy.integrate"),
        "matplotlib": _module_with_fallbacks("matplotlib"),
        "matplotlib.pyplot": _module_with_fallbacks("matplotlib.pyplot"),
        "ltsep": fake_ltsep,
        "binning_helpers": _module_with_fallbacks("binning_helpers"),
        "theta_cm": _module_with_fallbacks("theta_cm"),
        "utility": _module_with_fallbacks("utility"),
        "prompt_trees": _module_with_fallbacks("prompt_trees"),
        "pion_component_shapes": _module_with_fallbacks("pion_component_shapes"),
        "pion_component_fits": _module_with_fallbacks("pion_component_fits"),
        "pion_component_subtraction": _module_with_fallbacks("pion_component_subtraction"),
        "root_histogram_ownership": _module_with_fallbacks("root_histogram_ownership"),
        "proton_contamination_weights": _module_with_fallbacks("proton_contamination_weights"),
        "mm_background_subtraction": _module_with_fallbacks("mm_background_subtraction"),
    }
    with mock.patch.dict(sys.modules, fake_modules, clear=False):
        spec = importlib.util.spec_from_file_location(
            "calculate_yield_bg_opt_prepass_test", CALCULATE_YIELD_PATH
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def _load_bg_optimization_module():
    fake_frozen_manifest = _module_with_fallbacks("frozen_manifest")
    fake_frozen_manifest.load_zeroth_iteration_input_bundle = lambda *_args, **_kwargs: {}
    fake_find_bins = _module_with_fallbacks("find_bins")
    fake_find_bins.apply_bin_proposal = lambda *_args, **_kwargs: None
    fake_find_bins.propose_bins = lambda *_args, **_kwargs: None
    fake_find_bins.write_bin_interval_files = lambda *_args, **_kwargs: None
    with mock.patch.dict(
        sys.modules,
        {
            "frozen_manifest": fake_frozen_manifest,
            "find_bins": fake_find_bins,
        },
        clear=False,
    ):
        spec = importlib.util.spec_from_file_location(
            "bg_optimization_prepass_test", BG_OPTIMIZATION_PATH
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def _production_input():
    return {
        "ParticleType": "kaon",
        "Q2": "4p4",
        "W": "2p74",
        "EPSSET": "low",
        "POL": "unpolarized",
        "shift_mode": "raw",
        "particle_subtraction_mode": "simc_shape_components",
        "particle_subtraction_fallback_mode": "error",
        "pion_subtraction_scope": "t_bin",
        "canonical_t_binning": {"pair": {"id": "frozen-pair"}},
    }


class BGOptimizerPrepassContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calculate_yield = _load_calculate_yield_module()
        cls.bg_optimization = _load_bg_optimization_module()
        cls.t_edges = np.asarray([0.4, 0.6], dtype=float)
        cls.phi_edges = np.asarray([-180.0, 180.0], dtype=float)

    def _prepass_cache(self, production_input, phi_setting="Center"):
        prepass_input = build_bg_optimization_prepass_config(production_input)
        return {
            "_test_cache_sentinel": "cached-scalar-prepass",
            "_cache_contract": self.calculate_yield._build_bg_opt_data_cache_contract(
                prepass_input,
                phi_setting,
                self.t_edges,
                self.phi_edges,
            ),
        }

    def test_prepass_config_is_top_level_only_and_keeps_production_input_unchanged(self):
        production_input = _production_input()
        original = deepcopy(production_input)
        prepass_input = build_bg_optimization_prepass_config(production_input)
        prepass_input["bg_stat_scale1_by_setting"] = {"low_center": 0.0}

        self.assertEqual(production_input, original)
        self.assertEqual(prepass_input["particle_subtraction_mode"], "single_scale")
        self.assertEqual(prepass_input["particle_subtraction_fallback_mode"], "single_scale")
        self.assertTrue(prepass_input["bg_opt_prepass_mode"])
        self.assertEqual(prepass_input["analysis_subtraction_context"], BG_OPT_PREPASS_CONTEXT)
        self.assertEqual(prepass_input["pion_subtraction_scope"], "t_bin")

    def test_marked_t_bin_prepass_uses_only_matching_cache(self):
        production_input = _production_input()
        prepass_input = build_bg_optimization_prepass_config(production_input)
        prepass_input["bg_opt_use_data_cache"] = True
        cache = self._prepass_cache(production_input)
        hist = {"_bg_opt_data_base_cache": cache}

        resolved = self.calculate_yield._resolve_bg_opt_prepass_data_cache(
            hist,
            prepass_input,
            self.t_edges,
            self.phi_edges,
            "Center",
        )
        self.assertIs(resolved, cache)
        self.assertEqual(resolved["_test_cache_sentinel"], "cached-scalar-prepass")

        with self.assertRaisesRegex(RuntimeError, "cache_contract_mismatch: phi_edges"):
            self.calculate_yield._resolve_bg_opt_prepass_data_cache(
                hist,
                prepass_input,
                self.t_edges,
                np.asarray([-180.0, 0.0, 180.0]),
                "Center",
            )

        unmarked = dict(production_input, bg_opt_use_data_cache=True)
        with self.assertRaisesRegex(RuntimeError, "requires_prepass_mode"):
            self.calculate_yield._resolve_bg_opt_prepass_data_cache(
                hist,
                unmarked,
                self.t_edges,
                self.phi_edges,
                "Center",
            )

        self.assertIsNone(
            self.calculate_yield._resolve_bg_opt_prepass_data_cache(
                hist,
                production_input,
                self.t_edges,
                self.phi_edges,
                "Center",
            )
        )

    def test_calculate_yield_data_delegates_cache_selection_to_prepass_guard(self):
        source = CALCULATE_YIELD_PATH.read_text(encoding="utf-8")
        start = source.index("def calculate_yield_data(")
        self.assertIn(
            "data_base_cache = _resolve_bg_opt_prepass_data_cache(",
            source[start:],
        )

    def test_optimizer_candidate_is_scalar_and_consumes_prepared_cache(self):
        production_input = _production_input()
        original = deepcopy(production_input)
        cache = self._prepass_cache(production_input)
        captured = {"parent_resolver_called": False, "free_component_fit_called": False}

        fake_calculate_yield = types.ModuleType("calculate_yield")

        def prepare_bg_opt_data_base_cache(_hist, inp_dict, t_edges, phi_edges):
            self.assertEqual(inp_dict, production_input)
            self.assertTrue(np.array_equal(t_edges, self.t_edges))
            self.assertTrue(np.array_equal(phi_edges, self.phi_edges))
            captured["cache_built_from"] = build_bg_optimization_prepass_config(inp_dict)
            return cache

        def validate_cache(data_base_cache, inp_dict, phi_setting, t_edges, phi_edges):
            self.assertIs(data_base_cache, cache)
            self.assertEqual(inp_dict["particle_subtraction_mode"], "single_scale")
            self.assertEqual(inp_dict["particle_subtraction_fallback_mode"], "single_scale")
            self.assertTrue(inp_dict["bg_opt_prepass_mode"])
            self.assertEqual(inp_dict["analysis_subtraction_context"], BG_OPT_PREPASS_CONTEXT)
            self.assertEqual(phi_setting, "Center")
            self.assertTrue(np.array_equal(t_edges, self.t_edges))
            self.assertTrue(np.array_equal(phi_edges, self.phi_edges))
            return True

        def find_yield_data(histograms, candidate_input):
            candidate_hist = histograms[0]
            captured["candidate_input"] = deepcopy(candidate_input)
            self.assertIs(candidate_hist["_bg_opt_data_base_cache"], cache)
            self.assertEqual(candidate_hist["_bg_opt_data_base_cache"]["_test_cache_sentinel"], "cached-scalar-prepass")
            captured["cached_path_called"] = True
            return {}

        fake_calculate_yield.prepare_bg_opt_data_base_cache = prepare_bg_opt_data_base_cache
        fake_calculate_yield._validate_bg_opt_data_cache_contract = validate_cache
        fake_calculate_yield.find_yield_data = find_yield_data
        fake_calculate_yield.find_yield_simc = lambda *_args, **_kwargs: {}
        fake_calculate_yield._resolve_required_t_bin_parent = lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("parent resolver called")
        )
        fake_calculate_yield.build_particle_subtraction_component_result = lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("free component fitter called")
        )

        fake_calculate_ratio = types.ModuleType("calculate_ratio")
        fake_calculate_ratio.find_ratio = lambda *_args, **_kwargs: {
            "binned": {"Center": {"ratio": {}}}
        }

        optimizer = self.bg_optimization
        with mock.patch.dict(
            sys.modules,
            {
                "calculate_yield": fake_calculate_yield,
                "calculate_ratio": fake_calculate_ratio,
            },
            clear=False,
        ), mock.patch.object(optimizer, "get_forced_bg_scale1", return_value=0.0), mock.patch.object(
            optimizer, "get_forced_bg_scale2", return_value=0.0
        ):
            result = optimizer._optimize_phi_scale(
                {"phi_setting": "Center"},
                production_input,
                self.t_edges,
                self.phi_edges,
            )

        self.assertEqual(production_input, original)
        self.assertEqual(captured["cache_built_from"]["particle_subtraction_mode"], "single_scale")
        self.assertTrue(captured["cached_path_called"])
        self.assertEqual(captured["candidate_input"]["particle_subtraction_mode"], "single_scale")
        self.assertTrue(captured["candidate_input"]["bg_opt_use_data_cache"])
        self.assertTrue(captured["candidate_input"]["bg_opt_prepass_mode"])
        self.assertFalse(captured["parent_resolver_called"])
        self.assertFalse(captured["free_component_fit_called"])
        self.assertEqual(result["prepass_contract"]["context"], BG_OPT_PREPASS_CONTEXT)


if __name__ == "__main__":
    unittest.main()
