"""Focused deterministic contracts for the E.8.2 baseline-audit sidecar."""

from __future__ import annotations

import importlib.util
import inspect
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
UTILITY_DIR = REPO_ROOT / "src" / "utility"
CUTS_DIR = REPO_ROOT / "src" / "cuts"
BINNING_DIR = REPO_ROOT / "src" / "binning"
CALCULATE_YIELD_PATH = BINNING_DIR / "calculate_yield.py"
for directory in (UTILITY_DIR, CUTS_DIR, BINNING_DIR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

import full_background_subtraction_plots as plots


class _FallbackModule(types.ModuleType):
    def __getattr__(self, name):
        value = lambda *_args, **_kwargs: None
        setattr(self, name, value)
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


def _clone_root_histogram(histogram, *_args, **kwargs):
    if histogram is None:
        return None
    clone = histogram.Clone(kwargs.get("name", "detached_histogram"))
    if kwargs.get("reset"):
        clone.Reset()
    return clone


def _load_calculate_yield_module():
    """Load the producer with lightweight dependency stubs, never a fake producer."""
    fake_root = _module_with_fallbacks("ROOT")
    fake_ltsep = _module_with_fallbacks("ltsep")
    fake_ltsep.Root = lambda *_args, **_kwargs: _FakeRootContext()
    fake_ltsep.Misc = object()
    fake_utility = _module_with_fallbacks("utility")
    fake_utility.is_hist = lambda value: hasattr(value, "Clone")
    fake_utility.integrate_hist_range = lambda histogram, _low, _high: histogram.Integral()
    fake_ownership = _module_with_fallbacks("root_histogram_ownership")
    fake_ownership.clone_root_histogram = _clone_root_histogram
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
        "utility": fake_utility,
        "prompt_trees": _module_with_fallbacks("prompt_trees"),
        "pion_component_shapes": _module_with_fallbacks("pion_component_shapes"),
        "pion_component_fits": _module_with_fallbacks("pion_component_fits"),
        "pion_component_subtraction": _module_with_fallbacks("pion_component_subtraction"),
        "root_histogram_ownership": fake_ownership,
        "proton_contamination_weights": _module_with_fallbacks("proton_contamination_weights"),
        "mm_background_subtraction": _module_with_fallbacks("mm_background_subtraction"),
    }
    with mock.patch.dict(sys.modules, fake_modules, clear=False):
        spec = importlib.util.spec_from_file_location(
            "calculate_yield_e8_2_baseline_stage_audit_test", CALCULATE_YIELD_PATH
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


class _Axis:
    def __init__(self, edges):
        self._edges = list(edges)

    def GetBinLowEdge(self, index):
        return self._edges[int(index) - 1]


class _Histogram:
    """Minimal signed TH1 stand-in with the producer operations under test."""

    def __init__(self, contents, edges=(1.0, 1.1, 1.2, 1.3), errors=None):
        self.contents = [float(value) for value in contents]
        self.edges = [float(value) for value in edges]
        self.errors = list(errors or [0.0 for _unused in self.contents])
        self.directory = object()
        self.fill_calls = []

    def Clone(self, _name):
        clone = type(self)(self.contents, self.edges, self.errors)
        clone.fill_calls = list(self.fill_calls)
        return clone

    def SetDirectory(self, directory):
        self.directory = directory

    def Reset(self):
        self.contents = [0.0 for _unused in self.contents]
        self.errors = [0.0 for _unused in self.errors]
        self.fill_calls = []

    def Scale(self, factor):
        factor = float(factor)
        self.contents = [value * factor for value in self.contents]
        self.errors = [abs(factor) * value for value in self.errors]

    def Add(self, other, factor=1.0):
        factor = float(factor)
        self.contents = [
            left + factor * right
            for left, right in zip(self.contents, other.contents)
        ]

    def Fill(self, value, weight=1.0):
        value = float(value)
        index = next(
            (
                item
                for item, (low, high) in enumerate(zip(self.edges, self.edges[1:]))
                if low <= value <= high
            ),
            None,
        )
        if index is None:
            return
        weight = float(weight)
        self.contents[index] += weight
        self.errors[index] = math.sqrt(self.errors[index] ** 2 + weight ** 2)
        self.fill_calls.append((value, weight))

    def GetNbinsX(self):
        return len(self.contents)

    def GetBinContent(self, index):
        return self.contents[int(index) - 1]

    def GetBinError(self, index):
        return self.errors[int(index) - 1]

    def GetXaxis(self):
        return _Axis(self.edges)

    def Integral(self):
        return sum(self.contents)


class _SinglePassTree:
    def __init__(self, events):
        self.events = tuple(events)
        self.traversals = 0

    def GetEntries(self):
        return len(self.events)

    def __iter__(self):
        self.traversals += 1
        return iter(self.events)


class _Event:
    P_hgcer_xAtCer = 0.0
    P_hgcer_yAtCer = 0.0
    ph_q = 0.0
    MM = 1.15
    Q2 = 4.4
    W = 2.74
    epsilon = 0.5
    ssxptar = 0.0
    ssyptar = 0.0
    hsxptar = 0.0
    hsyptar = 0.0


def _histogram(value, *, errors=None):
    return _Histogram((value, value, value), errors=errors)


def _single_bin_histogram(value, *, error=0.0):
    return _Histogram((value,), edges=(1.0, 1.3), errors=(error,))


def _stage_histograms():
    return {
        "prompt_pre_proton": _single_bin_histogram(10.0),
        "random_component_pre_proton": _single_bin_histogram(2.0),
        "after_random_pre_proton": _single_bin_histogram(8.0),
        "dummy_component_pre_proton": _single_bin_histogram(1.0),
        "after_dummy_pre_proton": _single_bin_histogram(7.0),
        "proton_component_removed": _single_bin_histogram(3.0),
        "after_proton_pre_prune": _single_bin_histogram(4.0),
        "after_proton_post_prune": _single_bin_histogram(4.0),
        "pion_input": _single_bin_histogram(4.0),
        "pion_component_removed": _single_bin_histogram(1.0),
        "after_pion_final": _single_bin_histogram(3.0),
    }


def _source(*, active_profile="no_empirical_residual", broken_closure=False):
    stages = {
        "prompt_pre_proton": _histogram(10.0),
        "random_component_pre_proton": _histogram(2.0),
        "after_random_pre_proton": _histogram(8.0),
        "dummy_component_pre_proton": _histogram(1.0),
        "after_dummy_pre_proton": _histogram(7.0),
        "proton_component_removed": _histogram(3.0),
        "after_proton_pre_prune": _histogram(4.0),
        "after_proton_post_prune": _histogram(4.0),
        "pion_input": _histogram(4.0),
        "pion_component_removed": _histogram(1.0),
        "after_pion_final": _histogram(3.0),
    }
    if broken_closure:
        stages["after_proton_pre_prune"] = _histogram(5.0)
    return {
        "schema_version": plots.E8_2_SOURCE_SCHEMA_VERSION,
        "available": True,
        "reason": None,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "active_profile": active_profile,
        "setting": "Left",
        "epsilon": "lowe",
        "t_edges": [0.0, 1.0],
        "phi_edges": [-180.0, 180.0],
        "mm_edges": [1.0, 1.1, 1.2, 1.3],
        "lambda_window": [1.08, 1.18],
        "children": [{
            "t_index": 0,
            "t_low": 0.0,
            "t_high": 1.0,
            "phi_index": 0,
            "phi_low": -180.0,
            "phi_high": 180.0,
            "valid": True,
            "reason": None,
            "pion_application_status": "accepted",
            "pion_application_reason": None,
            "stages": stages,
            "stage_window_integrals": [
                {"stage": "prompt_pre_proton", "value": 10.0},
                {"stage": "after_random_pre_proton", "value": 8.0},
                {"stage": "after_dummy_pre_proton", "value": 7.0},
                {"stage": "after_proton_pre_prune", "value": 4.0},
                {"stage": "after_proton_post_prune", "value": 4.0},
                {"stage": "after_pion_final", "value": 3.0},
            ],
            "final_yield": 3.0,
            "statistical_error": 0.5,
            "total_error": 0.75,
        }],
    }


def _render_payloads():
    return {
        key: {"available": False, "reason": "fixture"}
        for key in (
            "d6", "d7", "d8", "d9", "d10", "d11", "e2", "e3", "e4",
            "e6", "e7", "f1", "e72", "e8",
        )
    }


def _manifest_setting():
    return {
        "kinematic_token": "Q4p4W2p74",
        "epsilon_filename_token": "lowe",
        "phi_setting": "Left",
        "particle_type": "kaon",
    }


def _public_yield_fixture(calculate_yield):
    """Build the smallest real calculate_yield_data input around bin_data."""
    stage_window_yields = {
        "raw_prompt": 21.0,
        "after_random_subtraction": 18.0,
        "after_dummy_subtraction": 16.0,
        "after_pion_subtraction": 10.0,
    }
    processed_entry = {
        "child_valid": True,
        "stage_window_yields": stage_window_yields,
        "_e8_2_baseline_stage_capture": {
            "available": True,
            "reason": None,
            "early_stages": _stage_histograms(),
            "pion_stage": calculate_yield._build_e8_2_pion_stage_capture({
                "accepted": True,
                "child_valid": True,
                "H_MM_nosub_before_pion_subtraction": _single_bin_histogram(4.0),
                "H_pion_subtraction_template_MM_nosub": _single_bin_histogram(1.0),
                "H_MM_nosub_after_pion_subtraction": _single_bin_histogram(3.0),
            }),
        },
    }
    final_hist = _single_bin_histogram(10.0, error=0.3)
    dummy_hist = _single_bin_histogram(5.0, error=0.4)
    binned_dict = {
        "kaon": {
            "processed_dict": {"t_bin1phi_bin1": processed_entry},
            "support_hist_dict": {},
            "binned_t_data": [(0.0, 1.0)],
            "binned_hist_data": [((0.0, 1.0), (10.0,))],
            "binned_hist_sub": [],
            "mm_hist_data": [final_hist],
            "mm_hist_dummy_norm": [dummy_hist],
            "mm_hist_sub": [],
            "child_validity": [[True]],
            "scale_factor": [[0.0]],
        }
    }
    hist = {
        "InFile_DATA": None,
        "InFile_DUMMY": None,
        "normfac_data": 1.0,
        "normfac_dummy": 1.0,
        "nWindows": 4,
        "phi_setting": "Left",
    }
    inp_dict = {
        "ParticleType": "kaon",
        "EPSSET": "lowe",
        "data_charge_err_left": 0.02,
        "dummy_charge_err_left": 0.03,
        "mm_min": 1.08,
        "mm_max": 1.18,
    }
    return {
        "binned_dict": binned_dict,
        "dummy_hist": dummy_hist,
        "final_hist": final_hist,
        "hist": hist,
        "inp_dict": inp_dict,
        "processed_entry": processed_entry,
        "stage_window_yields": stage_window_yields,
    }


class E82BaselineStageAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calculate_yield = _load_calculate_yield_module()

    def _early_capture(self):
        return {
            "pre_proton": _single_bin_histogram(0.0),
            "post_proton": _single_bin_histogram(0.0),
            "proton_removed": _single_bin_histogram(0.0),
        }

    def _production_group(self):
        return {
            name: [[_single_bin_histogram(0.0)]]
            for name in ("fit1sub", "pisub", "nosub")
        }

    def _process_one_event(self, proton_factor):
        capture = {key: [[histogram]] for key, histogram in self._early_capture().items()}
        production = self._production_group()
        tree = _SinglePassTree([_Event()])
        payload = {"final_cleaned_factor": float(proton_factor), "rf_accept": True}
        with patch.object(self.calculate_yield, "get_kaon_proton_cleaning_event_payload", return_value=payload), patch.object(
            self.calculate_yield, "find_2d_bin_indices", return_value=(0, 0)
        ), patch.object(
            self.calculate_yield, "calculate_theta_cm_deg", return_value=0.0
        ), patch.object(
            self.calculate_yield, "_append_ave_event"
        ):
            self.calculate_yield._process_yield_data_tree(
                tree, [], production, lambda _evt: 1.15, lambda _evt: 0.5,
                [0.0, 1.0], [-180.0, 180.0], "kaon", "unpolarized", None,
                lambda _evt, _low, _high: (False, True, None), 1.08, 1.18,
                lambda *_args, **_kwargs: None, source_label="prompt",
                proton_cleaning_result={"accepted": True}, e8_2_capture=capture,
            )
        return tree, capture, production

    def _raw_capture(self, values):
        return {
            role: {
                channel: [[_single_bin_histogram(contents)]]
                for channel, contents in channels.items()
            }
            for role, channels in values.items()
        }

    def _accepted_pion_payload(self):
        return {
            "accepted": True,
            "child_valid": True,
            "H_MM_nosub_before_pion_subtraction": _single_bin_histogram(4.0),
            "H_pion_subtraction_template_MM_nosub": _single_bin_histogram(1.0),
            "H_MM_nosub_after_pion_subtraction": _single_bin_histogram(3.0),
        }

    def _source_from_producer(self, *, capture=None, pion_stage=None,
                              child_valid=True, measurement=None):
        capture = capture or {
            "available": True,
            "reason": None,
            "early_stages": _stage_histograms(),
        }
        pion_stage = pion_stage or self.calculate_yield._build_e8_2_pion_stage_capture(
            self._accepted_pion_payload()
        )
        processed = {
            "t_bin1phi_bin1": {
                "child_valid": bool(child_valid),
                "_e8_2_baseline_stage_capture": {
                    **capture,
                    "pion_stage": pion_stage,
                },
            }
        }
        if measurement is None:
            measurement = {"yield": 3.0, "statistical_error": 0.5, "total_error": 0.75}
        with patch.object(self.calculate_yield, "get_active_bg_profile_name", return_value="no_empirical_residual"):
            return self.calculate_yield._build_e8_2_baseline_stage_source(
                {"phi_setting": "Left"}, processed, [0.0, 1.0], [-180.0, 180.0],
                {"EPSSET": "lowe", "mm_min": 1.08, "mm_max": 1.18},
                {} if measurement is False else {(0, 0): measurement},
            )

    def _finalization_hist(self, pdf_path, manifest_path):
        return {
            "_e8_2_baseline_stage_source": _source(),
            "_e8_2_full_background_render_state": (
                plots.capture_full_background_subtraction_e8_2_render_state(
                    pdf_path=str(pdf_path), page_manifest_path=str(manifest_path),
                    page_manifest_setting=_manifest_setting(), payloads=_render_payloads(),
                )
            ),
        }

    def _run_public_yield_fixture(self, *, build_e8_2_source):
        fixture = _public_yield_fixture(self.calculate_yield)
        integration = mock.Mock(side_effect=((10.0, 0.3), (5.0, 0.4)))
        common_patches = (
            patch.object(self.calculate_yield, "_resolve_bg_opt_prepass_data_cache", return_value=None),
            patch.object(self.calculate_yield, "_get_cached_kaon_signal_shape_payload", return_value=None),
            patch.object(self.calculate_yield, "_get_cached_kaon_sigma0_shape_payload", return_value=None),
            patch.object(self.calculate_yield, "bin_data", return_value=(
                fixture["binned_dict"], None, None,
            )),
            patch.object(self.calculate_yield, "get_active_bg_profile_name", return_value="no_empirical_residual"),
            patch.object(self.calculate_yield, "integral_with_stat_error", integration),
            patch.object(self.calculate_yield, "resolve_particle_subtraction_mode", return_value="legacy"),
        )
        with common_patches[0], common_patches[1], common_patches[2], common_patches[3], common_patches[4], common_patches[5], common_patches[6]:
            if build_e8_2_source:
                groups = self.calculate_yield.calculate_yield_data(
                    "kaon", fixture["hist"], [0.0, 1.0], [-180.0, 180.0],
                    fixture["inp_dict"],
                )
            else:
                with patch.object(
                    self.calculate_yield,
                    "_build_e8_2_baseline_stage_source",
                    return_value={"available": False, "reason": "fixture_sidecar_disabled"},
                ):
                    groups = self.calculate_yield.calculate_yield_data(
                        "kaon", fixture["hist"], [0.0, 1.0], [-180.0, 180.0],
                        fixture["inp_dict"],
                    )
        return fixture, groups, integration

    def test_same_traversal_capture_factor_one_preserves_production_fill(self):
        tree, capture, production = self._process_one_event(1.0)

        self.assertEqual(tree.traversals, 1)
        self.assertEqual(capture["pre_proton"][0][0].contents, [1.0])
        self.assertEqual(capture["post_proton"][0][0].contents, [1.0])
        self.assertEqual(capture["proton_removed"][0][0].contents, [0.0])
        for name in ("fit1sub", "pisub", "nosub"):
            self.assertEqual(production[name][0][0].contents, [1.0])

    def test_same_traversal_capture_nontrivial_factor_preserves_production_fill(self):
        tree, capture, production = self._process_one_event(0.4)

        self.assertEqual(tree.traversals, 1)
        self.assertEqual(capture["pre_proton"][0][0].contents, [1.0])
        self.assertAlmostEqual(capture["post_proton"][0][0].contents[0], 0.4)
        self.assertAlmostEqual(capture["proton_removed"][0][0].contents[0], 0.6)
        for name in ("fit1sub", "pisub", "nosub"):
            self.assertAlmostEqual(production[name][0][0].contents[0], 0.4)

    def test_early_stage_algebra_uses_existing_normalizations_and_closes(self):
        captures = self._raw_capture({
            "prompt": {"pre_proton": 5.0, "post_proton": 2.0, "proton_removed": 3.0},
            "random": {"pre_proton": 4.0, "post_proton": 2.0, "proton_removed": 2.0},
            "dummy_prompt": {"pre_proton": 5.0, "post_proton": 3.0, "proton_removed": 2.0},
            "dummy_random": {"pre_proton": 4.0, "post_proton": 4.0, "proton_removed": 0.0},
        })
        stages, reason = self.calculate_yield._build_e8_2_early_stage_capture(
            captures, 0, 0, 4, 2.0, 0.5,
            _single_bin_histogram(2.0), _single_bin_histogram(2.0),
        )

        self.assertIsNone(reason)
        self.assertEqual(stages["prompt_pre_proton"].contents, [10.0])
        self.assertEqual(stages["random_component_pre_proton"].contents, [2.0])
        self.assertEqual(stages["after_random_pre_proton"].contents, [8.0])
        self.assertEqual(stages["dummy_component_pre_proton"].contents, [2.0])
        self.assertEqual(stages["after_dummy_pre_proton"].contents, [6.0])
        self.assertEqual(stages["proton_component_removed"].contents, [4.0])
        self.assertEqual(stages["after_proton_pre_prune"].contents, [2.0])
        self.assertEqual(stages["after_proton_post_prune"].contents, [2.0])

    def test_production_pruning_boundary_is_observed_without_sidecar_pruning(self):
        captures = self._raw_capture({
            "prompt": {"pre_proton": 5.0, "post_proton": 2.0, "proton_removed": 3.0},
            "random": {"pre_proton": 4.0, "post_proton": 2.0, "proton_removed": 2.0},
            "dummy_prompt": {"pre_proton": 5.0, "post_proton": 3.0, "proton_removed": 2.0},
            "dummy_random": {"pre_proton": 4.0, "post_proton": 4.0, "proton_removed": 0.0},
        })
        production = _Histogram((4.0, 2.0), edges=(1.0, 1.1, 1.2))
        prune_calls = []

        def existing_production_prune(histogram, event_threshold):
            prune_calls.append((histogram, event_threshold))
            histogram.contents[1] = 0.0

        with patch.object(
            self.calculate_yield, "prune_hist", side_effect=existing_production_prune,
        ):
            pre_prune = self.calculate_yield._clone_hist_for_plot(production)
            self.calculate_yield.prune_hist(production, 1)
            post_prune = self.calculate_yield._clone_hist_for_plot(production)
            stages, reason = self.calculate_yield._build_e8_2_early_stage_capture(
                captures, 0, 0, 4, 2.0, 0.5, pre_prune, post_prune,
            )

        self.assertIsNone(reason)
        self.assertEqual(prune_calls, [(production, 1)])
        self.assertEqual(stages["after_proton_pre_prune"].contents, [4.0, 2.0])
        self.assertEqual(stages["after_proton_post_prune"].contents, [4.0, 0.0])
        self.assertEqual(production.contents, [4.0, 0.0])

        producer_source = inspect.getsource(self.calculate_yield.process_hist_data)
        pre_prune_capture = producer_source.index("e8_2_after_proton_pre_prune = _clone_hist_for_plot")
        production_prune = producer_source.index(
            'prune_hist(\n                hist_bin_dict["H_MM_nosub_DATA_{}_{}".format(j, k)]'
        )
        post_prune_capture = producer_source.index("e8_2_after_proton_post_prune = _clone_hist_for_plot")
        pion_subtraction = producer_source.index(
            "# Pion subtraction by scaling pion background to peak size",
            production_prune,
        )
        self.assertLess(pre_prune_capture, production_prune)
        self.assertLess(production_prune, post_prune_capture)
        self.assertLess(post_prune_capture, pion_subtraction)
        self.assertEqual(producer_source.count("prune_hist("), 15)

    def test_proton_closure_uses_pre_prune_and_post_prune_may_differ(self):
        source = _source()
        stages = source["children"][0]["stages"]
        stages["after_proton_post_prune"] = _histogram(2.0)
        stages["pion_input"] = _histogram(2.0)
        stages["after_pion_final"] = _histogram(1.0)
        source["children"][0]["stage_window_integrals"][4]["value"] = 2.0
        source["children"][0]["stage_window_integrals"][5]["value"] = 1.0

        payload = plots.build_full_background_subtraction_e8_2_payload(source)

        self.assertTrue(payload["available"])
        child = payload["per_t"][0]["children"][0]
        self.assertEqual(child["stages"]["after_proton_pre_prune"].contents, [4.0, 4.0, 4.0])
        self.assertEqual(child["stages"]["after_proton_post_prune"].contents, [2.0, 2.0, 2.0])

    def test_post_prune_pion_handoff_mismatch_fails_closed_in_source_and_reader(self):
        binning_payload = self._accepted_pion_payload()
        binning_payload["H_MM_nosub_before_pion_subtraction"] = _Histogram(
            (4.0,), edges=(1.0, 1.2),
        )
        binning_payload["H_pion_subtraction_template_MM_nosub"] = _Histogram(
            (1.0,), edges=(1.0, 1.2),
        )
        binning_payload["H_MM_nosub_after_pion_subtraction"] = _Histogram(
            (3.0,), edges=(1.0, 1.2),
        )
        binning_source = self._source_from_producer(
            pion_stage=self.calculate_yield._build_e8_2_pion_stage_capture(binning_payload)
        )
        content_payload = self._accepted_pion_payload()
        content_payload["H_MM_nosub_before_pion_subtraction"] = _single_bin_histogram(3.0)
        content_payload["H_pion_subtraction_template_MM_nosub"] = _single_bin_histogram(1.0)
        content_payload["H_MM_nosub_after_pion_subtraction"] = _single_bin_histogram(2.0)
        content_source = self._source_from_producer(
            pion_stage=self.calculate_yield._build_e8_2_pion_stage_capture(content_payload)
        )

        self.assertFalse(binning_source["available"])
        self.assertEqual(
            binning_source["reason"],
            "e8_2_authority_failure_post_prune_pion_input_binning_inconsistent",
        )
        self.assertFalse(content_source["available"])
        self.assertEqual(
            content_source["reason"],
            "e8_2_authority_failure_post_prune_pion_input_content_mismatch",
        )

        reader_source = _source()
        reader_stages = reader_source["children"][0]["stages"]
        reader_stages["pion_input"] = _histogram(2.0)
        reader_stages["after_pion_final"] = _histogram(1.0)
        reader_source["children"][0]["stage_window_integrals"][5]["value"] = 1.0
        reader = plots.build_full_background_subtraction_e8_2_payload(reader_source)
        self.assertFalse(reader["available"])
        self.assertEqual(reader["reason"], "e8_2_post_prune_pion_input_mismatch")

    def test_pion_authority_and_yield_authority_are_exact_and_separate(self):
        payload = self._accepted_pion_payload()
        original_template_contents = list(
            payload["H_pion_subtraction_template_MM_nosub"].contents
        )
        pion_stage = self.calculate_yield._build_e8_2_pion_stage_capture(payload)
        source = self._source_from_producer(pion_stage=pion_stage)

        self.assertTrue(pion_stage["available"])
        self.assertIsNot(
            pion_stage["stages"]["pion_component_removed"],
            payload["H_pion_subtraction_template_MM_nosub"],
        )
        self.assertEqual(pion_stage["stages"]["pion_component_removed"].contents, [1.0])
        self.assertEqual(
            pion_stage["stages"]["pion_input"].contents[0]
            - pion_stage["stages"]["pion_component_removed"].contents[0],
            pion_stage["stages"]["after_pion_final"].contents[0],
        )
        self.assertEqual(
            payload["H_pion_subtraction_template_MM_nosub"].contents,
            original_template_contents,
        )
        self.assertTrue(source["available"])
        child = source["children"][0]
        self.assertEqual(child["stages"]["after_pion_final"].contents, [3.0])
        self.assertEqual(child["final_yield"], 3.0)
        self.assertEqual(child["statistical_error"], 0.5)
        self.assertEqual(child["total_error"], 0.75)
        self.assertEqual(
            tuple(row["stage"] for row in child["stage_window_integrals"]),
            (
                "prompt_pre_proton", "after_random_pre_proton",
                "after_dummy_pre_proton", "after_proton_pre_prune",
                "after_proton_post_prune", "after_pion_final",
            ),
        )

    def test_current_chain_authority_failures_close_the_entire_source(self):
        missing_early = self._source_from_producer(capture={
            "available": False, "reason": "e8_2_auxiliary_source_histogram_missing",
            "early_stages": {},
        })
        missing_yield = self._source_from_producer(measurement=False)
        nonfinite_yield = self._source_from_producer(measurement={
            "yield": 3.0, "statistical_error": float("nan"), "total_error": 0.75,
        })
        malformed_yield = self._source_from_producer(measurement={
            "yield": 3.0, "statistical_error": 0.5, "total_error": "not-a-number",
        })
        missing_template_payload = self._accepted_pion_payload()
        missing_template_payload["H_pion_subtraction_template_MM_nosub"] = None
        missing_template = self._source_from_producer(
            pion_stage=self.calculate_yield._build_e8_2_pion_stage_capture(
                missing_template_payload
            )
        )
        with patch.object(self.calculate_yield, "_clone_hist_for_plot", return_value=None):
            clone_failure_stage = self.calculate_yield._build_e8_2_pion_stage_capture(
                self._accepted_pion_payload()
            )
        clone_failure = self._source_from_producer(pion_stage=clone_failure_stage)

        self.assertEqual(missing_early["reason"], "e8_2_authority_failure_early_stage_capture_missing")
        self.assertEqual(missing_yield["reason"], "e8_2_authority_failure_final_yield_missing")
        self.assertEqual(nonfinite_yield["reason"], "e8_2_authority_failure_final_yield_nonfinite")
        self.assertEqual(malformed_yield["reason"], "e8_2_authority_failure_final_yield_malformed")
        self.assertEqual(
            missing_template["reason"],
            "e8_2_authority_failure_accepted_pion_application_missing_exact_wide_object",
        )
        self.assertEqual(
            clone_failure["reason"],
            "e8_2_authority_failure_accepted_pion_application_clone_failed",
        )
        for source in (
            missing_early, missing_yield, nonfinite_yield, malformed_yield,
            missing_template, clone_failure,
        ):
            self.assertFalse(source["available"])

    def test_authoritative_skip_state_remains_child_level_unavailable(self):
        source = self._source_from_producer(
            child_valid=False,
            pion_stage={
                "available": False,
                "reason": "frozen_parent_skip_bin",
                "pion_application_status": "skip_bin",
                "pion_application_reason": "frozen_parent_skip_bin",
                "stages": {},
            },
        )

        self.assertFalse(source["available"])
        self.assertEqual(source["reason"], "no_valid_e8_2_canonical_child")
        child = source["children"][0]
        self.assertEqual(child["reason"], "frozen_parent_skip_bin")
        self.assertEqual(child["pion_application_status"], "skip_bin")

    def test_e8_2_helpers_keep_only_accepted_stages_and_never_call_bg_fit(self):
        captures = self._raw_capture({
            "prompt": {"pre_proton": 5.0, "post_proton": 2.0, "proton_removed": 3.0},
            "random": {"pre_proton": 4.0, "post_proton": 2.0, "proton_removed": 2.0},
            "dummy_prompt": {"pre_proton": 5.0, "post_proton": 3.0, "proton_removed": 2.0},
            "dummy_random": {"pre_proton": 4.0, "post_proton": 4.0, "proton_removed": 0.0},
        })
        bg_fit = mock.Mock(side_effect=AssertionError("bg_fit must stay dormant"))
        with patch.object(self.calculate_yield, "bg_fit", bg_fit, create=True):
            early_stages, early_reason = self.calculate_yield._build_e8_2_early_stage_capture(
                captures, 0, 0, 4, 2.0, 0.5,
                _single_bin_histogram(2.0), _single_bin_histogram(2.0),
            )
            pion_stage = self.calculate_yield._build_e8_2_pion_stage_capture(
                self._accepted_pion_payload()
            )
            source = self._source_from_producer(pion_stage=pion_stage)
            payload = plots.build_full_background_subtraction_e8_2_payload(source)

        self.assertIsNone(early_reason)
        self.assertEqual(early_stages["after_proton_pre_prune"].contents, [2.0])
        self.assertTrue(pion_stage["available"])
        self.assertTrue(source["available"])
        self.assertTrue(payload["available"])
        child = payload["per_t"][0]["children"][0]
        self.assertEqual(
            tuple(row["stage"] for row in child["stage_window_integrals"]),
            (
                "prompt_pre_proton", "after_random_pre_proton",
                "after_dummy_pre_proton", "after_proton_pre_prune",
                "after_proton_post_prune", "after_pion_final",
            ),
        )
        labels = " ".join(child["stages"]).lower()
        self.assertNotIn("fit", labels)
        self.assertNotIn("method_a", child["stages"])
        self.assertNotIn("method_b", child["stages"])
        self.assertEqual(child["final_yield"], 3.0)
        self.assertEqual(child["statistical_error"], 0.5)
        self.assertEqual(child["total_error"], 0.75)
        bg_fit.assert_not_called()

    def test_non_current_profile_and_bad_closure_fail_closed(self):
        profile_failure = plots.build_full_background_subtraction_e8_2_payload(
            _source(active_profile="nominal_weighted")
        )
        closure_failure = plots.build_full_background_subtraction_e8_2_payload(
            _source(broken_closure=True)
        )

        self.assertFalse(profile_failure["available"])
        self.assertEqual(profile_failure["reason"], "active_profile_not_no_empirical_residual")
        self.assertFalse(closure_failure["available"])
        self.assertEqual(closure_failure["reason"], "e8_2_stage_closure_failed")

    def test_valid_child_schema_rejects_fit_and_method_stage_injection(self):
        for forbidden_stage in (
            "fit_1_empirical_residual",
            "fit_2_empirical_residual",
            "method_a_weight_factor",
            "method_b_template",
        ):
            with self.subTest(forbidden_stage=forbidden_stage):
                source = _source()
                source["children"][0]["stages"][forbidden_stage] = _histogram(0.0)
                payload = plots.build_full_background_subtraction_e8_2_payload(source)
                self.assertFalse(payload["available"])
                self.assertEqual(payload["reason"], "e8_2_valid_child_stage_inventory_invalid")

    def test_malformed_geometry_binning_schema_and_missing_step3_state_fail_closed(self):
        malformed_geometry = _source()
        malformed_geometry["t_edges"] = [0.0, 0.0]
        malformed_binning = _source()
        malformed_binning["mm_edges"] = [1.0, 1.15, 1.3]
        malformed_schema = _source()
        malformed_schema["schema_version"] = "e8_2_baseline_stage_source/v0"

        self.assertEqual(
            plots.build_full_background_subtraction_e8_2_payload(malformed_geometry)["reason"],
            "e8_2_canonical_geometry_invalid",
        )
        self.assertEqual(
            plots.build_full_background_subtraction_e8_2_payload(malformed_binning)["reason"],
            "e8_2_stage_histogram_binning_inconsistent",
        )
        self.assertEqual(
            plots.build_full_background_subtraction_e8_2_payload(malformed_schema)["reason"],
            "e8_2_source_schema_invalid",
        )
        status = plots.finalize_full_background_subtraction_e8_2({
            "_e8_2_baseline_stage_source": _source(),
        })
        self.assertEqual(status, {
            "status": "unavailable",
            "reason": "e8_2_step3_render_state_missing",
        })

    def test_private_integrals_and_statistical_error_remain_producer_owned(self):
        source = self._source_from_producer(measurement={
            "yield": 3.0, "statistical_error": 0.625, "total_error": 0.75,
        })
        source_text = inspect.getsource(self.calculate_yield.calculate_yield_data)
        reader_integral = mock.Mock(side_effect=AssertionError("reader must not integrate yields"))
        with patch.object(plots, "integral_with_stat_error", reader_integral, create=True):
            payload = plots.build_full_background_subtraction_e8_2_payload(source)

        self.assertTrue(payload["available"])
        child = payload["per_t"][0]["children"][0]
        self.assertEqual(child["statistical_error"], 0.625)
        self.assertEqual(
            tuple(row["stage"] for row in child["stage_window_integrals"]),
            (
                "prompt_pre_proton", "after_random_pre_proton",
                "after_dummy_pre_proton", "after_proton_pre_prune",
                "after_proton_post_prune", "after_pion_final",
            ),
        )
        self.assertIn("yld, yld_stat_err = integral_with_stat_error(final_hist)", source_text)
        self.assertIn('"statistical_error": float(yld_stat_err)', source_text)
        reader_integral.assert_not_called()

    def test_public_yield_outputs_and_stage_windows_survive_sidecar_assembly(self):
        baseline_fixture, baseline_groups, baseline_integration = self._run_public_yield_fixture(
            build_e8_2_source=False,
        )
        observed_fixture, observed_groups, observed_integration = self._run_public_yield_fixture(
            build_e8_2_source=True,
        )
        public_cell = observed_groups[(0, 0)]
        public_before_reader = dict(observed_groups[(0, 0)])
        stage_windows_before = dict(observed_fixture["stage_window_yields"])

        self.assertEqual(observed_groups, baseline_groups)
        self.assertEqual(public_cell["kaon"], 10.0)
        self.assertEqual(public_cell["kaon_err"], baseline_groups[(0, 0)]["kaon_err"])
        self.assertNotEqual(public_cell["kaon_err"], 0.3)
        self.assertEqual(baseline_integration.call_count, 2)
        self.assertEqual(observed_integration.call_count, 2)
        self.assertEqual(
            observed_integration.call_args_list,
            [
                mock.call(observed_fixture["final_hist"]),
                mock.call(observed_fixture["dummy_hist"]),
            ],
        )

        source = observed_fixture["hist"]["_e8_2_baseline_stage_source"]
        self.assertTrue(source["available"])
        child = source["children"][0]
        self.assertEqual(child["final_yield"], public_cell["kaon"])
        self.assertEqual(child["statistical_error"], 0.3)
        self.assertEqual(child["total_error"], public_cell["kaon_err"])
        self.assertEqual(
            observed_fixture["processed_entry"]["stage_window_yields"],
            stage_windows_before,
        )
        self.assertIs(
            observed_fixture["processed_entry"]["stage_window_yields"],
            observed_fixture["stage_window_yields"],
        )

        presentation = plots.build_full_background_subtraction_e8_2_payload(source)
        self.assertTrue(presentation["available"])
        self.assertEqual(observed_groups[(0, 0)], public_before_reader)
        self.assertEqual(
            observed_fixture["processed_entry"]["stage_window_yields"],
            stage_windows_before,
        )

    def test_step3_render_state_detaches_histograms(self):
        original = _histogram(4.0)
        payloads = _render_payloads()
        payloads["d6"] = {"available": True, "histogram": original}
        state = plots.capture_full_background_subtraction_e8_2_render_state(
            pdf_path="preliminary.pdf", page_manifest_path="preliminary-manifest.json",
            page_manifest_setting=_manifest_setting(), payloads=payloads,
        )

        original.contents[0] = 999.0
        retained = state["payloads"]["d6"]["histogram"]
        self.assertIsNot(retained, original)
        self.assertEqual(retained.contents[0], 4.0)

    def test_e82_pages_follow_existing_parent_pages_and_precede_handoff(self):
        order = []

        with patch.object(plots, "_import_root", return_value=object()), patch.object(
            plots, "_render_full_background_subtraction_e8_context_page",
            side_effect=lambda *_args: order.append("context") or True,
        ), patch.object(
            plots, "_render_full_background_subtraction_e8_parent_pages",
            side_effect=lambda *_args: order.append("parent"),
        ), patch.object(
            plots, "_render_full_background_subtraction_e8_2_pages",
            side_effect=lambda *_args: order.append("e8_2"),
        ), patch.object(
            plots, "_render_full_background_subtraction_e8_handoff_page",
            side_effect=lambda *_args: order.append("handoff") or True,
        ):
            result = plots.render_full_background_subtraction_procedure_pages(
                "ignored.pdf", {"available": False, "reason": "fixture"},
                {"available": False, "reason": "fixture"},
                e8_payload={"available": True, "parents": ({"t_index": 0},)},
                e8_2_payload={"available": True},
            )

        self.assertEqual(order, ["context", "parent", "e8_2", "handoff"])
        self.assertEqual(result["manifest"][-1]["page_id"], "full_background.e8.handoff")

    def test_proton_page_shows_the_explicit_pre_and_post_prune_handoff(self):
        payload = plots.build_full_background_subtraction_e8_2_payload(_source())
        observed = []

        def record_page(*_args, **kwargs):
            observed.append(kwargs)
            return True

        with patch.object(plots, "_e8_2_render_subtraction_page", side_effect=record_page), patch.object(
            plots, "_e8_2_render_final_mm_page", return_value=True
        ), patch.object(plots, "_e8_2_render_stage_yield_page", return_value=True):
            plots._render_full_background_subtraction_e8_2_pages(
                object(), "ignored.pdf", payload, [], [],
            )

        proton_page = next(page for page in observed if page["semantic_stage"] == "slow_proton_cleaning")
        self.assertIn("prune_hist lies between columns 3 and 4", proton_page["title"])
        self.assertEqual(
            proton_page["columns"],
            (
                ("after dummy pre-proton", "after_dummy_pre_proton"),
                ("proton component removed", "proton_component_removed"),
                ("production after proton, pre-prune", "after_proton_pre_prune"),
                ("after existing production prune_hist; pion input state", "after_proton_post_prune"),
            ),
        )

    def test_finalization_replaces_only_a_complete_temporary_pair(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pdf_path = root / "preliminary.pdf"
            manifest_path = root / "preliminary-manifest.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary-manifest", encoding="utf-8")
            hist = self._finalization_hist(pdf_path, manifest_path)

            def open_temporary(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def render_temporary(*_args, **kwargs):
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            with patch.object(plots, "open_full_background_subtraction_pdf", side_effect=open_temporary), patch.object(
                plots, "close_full_background_subtraction_pdf", return_value=True
            ), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=render_temporary
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)

            self.assertEqual(status["status"], "available")
            self.assertEqual(status["renderer_failures"], [])
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "replacement")
            self.assertEqual(
                json.loads(manifest_path.read_text(encoding="utf-8"))["pages"][0]["page_id"],
                "fixture",
            )
            self.assertNotIn("_e8_2_full_background_render_state", hist)
            self.assertEqual(hist["full_background_subtraction_page_manifest"][0]["page_id"], "fixture")

    def test_renderer_failures_preserve_preliminary_pair_and_render_state(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pdf_path = root / "preliminary.pdf"
            manifest_path = root / "preliminary-manifest.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary-manifest", encoding="utf-8")
            hist = self._finalization_hist(pdf_path, manifest_path)

            def open_temporary(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            with patch.object(plots, "open_full_background_subtraction_pdf", side_effect=open_temporary), patch.object(
                plots, "close_full_background_subtraction_pdf", return_value=True
            ), patch.object(
                plots, "render_full_background_subtraction_procedure_pages",
                return_value={"failures": ["fixture renderer failure"]},
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)

            self.assertEqual(status["status"], "unavailable")
            self.assertEqual(status["reason"], "e8_2_renderer_failures")
            self.assertEqual(status["renderer_failures"], ["fixture renderer failure"])
            self.assertTrue(status["preliminary_artifact_preserved"])
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "preliminary")
            self.assertEqual(manifest_path.read_text(encoding="utf-8"), "preliminary-manifest")
            self.assertIn("_e8_2_full_background_render_state", hist)

    def test_second_install_failure_restores_the_preliminary_pair(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pdf_path = root / "preliminary.pdf"
            manifest_path = root / "preliminary-manifest.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary-manifest", encoding="utf-8")
            hist = self._finalization_hist(pdf_path, manifest_path)
            real_replace = os.replace

            def open_temporary(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def render_temporary(*_args, **kwargs):
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            def fail_manifest_install(source, destination):
                if os.fspath(destination) == os.fspath(manifest_path):
                    raise OSError("fixture manifest installation failure")
                return real_replace(source, destination)

            with patch.object(plots, "open_full_background_subtraction_pdf", side_effect=open_temporary), patch.object(
                plots, "close_full_background_subtraction_pdf", return_value=True
            ), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=render_temporary
            ), patch.object(plots.os, "replace", side_effect=fail_manifest_install):
                status = plots.finalize_full_background_subtraction_e8_2(hist)

            self.assertEqual(status["status"], "unavailable")
            self.assertTrue(status["preliminary_artifact_preserved"])
            self.assertEqual(status["recovery_errors"], [])
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "preliminary")
            self.assertEqual(manifest_path.read_text(encoding="utf-8"), "preliminary-manifest")
            self.assertIn("_e8_2_full_background_render_state", hist)
            self.assertFalse(list(root.glob("*e8_2.recovery*")))

    def test_recovery_failure_never_claims_preliminary_pair_preserved(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pdf_path = root / "preliminary.pdf"
            manifest_path = root / "preliminary-manifest.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary-manifest", encoding="utf-8")
            hist = self._finalization_hist(pdf_path, manifest_path)
            real_replace = os.replace
            real_copy = plots._e8_2_copy_recovery_file

            def open_temporary(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def render_temporary(*_args, **kwargs):
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            def fail_manifest_install(source, destination):
                if os.fspath(destination) == os.fspath(manifest_path):
                    raise OSError("fixture manifest installation failure")
                return real_replace(source, destination)

            def fail_pdf_recovery(source, destination):
                if os.fspath(destination) == os.fspath(pdf_path):
                    raise OSError("fixture recovery failure")
                return real_copy(source, destination)

            with patch.object(plots, "open_full_background_subtraction_pdf", side_effect=open_temporary), patch.object(
                plots, "close_full_background_subtraction_pdf", return_value=True
            ), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=render_temporary
            ), patch.object(plots.os, "replace", side_effect=fail_manifest_install), patch.object(
                plots, "_e8_2_copy_recovery_file", side_effect=fail_pdf_recovery
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)

            self.assertEqual(status["status"], "unavailable")
            self.assertFalse(status["preliminary_artifact_preserved"])
            self.assertTrue(status["recovery_errors"])
            self.assertIn("pdf_restore", status["recovery_errors"][0])
            self.assertIn("_e8_2_full_background_render_state", hist)

    def test_failed_finalization_preserves_preliminary_outputs_before_installation(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            pdf_path = root / "preliminary.pdf"
            manifest_path = root / "preliminary-manifest.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary-manifest", encoding="utf-8")
            hist = self._finalization_hist(pdf_path, manifest_path)

            def open_temporary(path):
                Path(path).write_text("incomplete", encoding="utf-8")
                return True

            with patch.object(plots, "open_full_background_subtraction_pdf", side_effect=open_temporary), patch.object(
                plots, "render_full_background_subtraction_procedure_pages",
                side_effect=RuntimeError("fixture failure"),
            ), patch.object(plots, "close_full_background_subtraction_pdf", return_value=True):
                status = plots.finalize_full_background_subtraction_e8_2(hist)

            self.assertEqual(status["status"], "unavailable")
            self.assertTrue(status["preliminary_artifact_preserved"])
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "preliminary")
            self.assertEqual(manifest_path.read_text(encoding="utf-8"), "preliminary-manifest")
            self.assertIn("_e8_2_full_background_render_state", hist)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
