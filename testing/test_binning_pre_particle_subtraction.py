"""Regression checks for pre-particle-subtraction binning inputs."""

from __future__ import annotations

import importlib.util
import math
import sys
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
FIND_BINS_PATH = REPO_ROOT / "src" / "binning" / "find_bins.py"


class _FakeHistogram:
    def __init__(self, values):
        self.values = list(values)

    def Clone(self):
        return _FakeHistogram(self.values)

    def Scale(self, scale):
        self.values = [float(value) * float(scale) for value in self.values]


def _load_find_bins():
    fake_ltsep = types.ModuleType("ltsep")
    fake_ltsep.Misc = types.SimpleNamespace()
    fake_ltsep.Root = lambda *_args, **_kwargs: types.SimpleNamespace(LTANAPATH="")
    fake_utility = types.ModuleType("utility")
    fake_utility.flatten_hist = lambda histogram: list(histogram.values)
    fake_background_config = types.ModuleType("background_config")
    for name, value in {
        "MIN_PHI_BINS": 1,
        "MIN_T_BINS": 1,
        "PHI_BIN_MAX_DEG": 180.0,
        "PHI_BIN_MIN_DEG": -180.0,
        "PHI_BIN_MIN_EVENTS": 1,
        "T_BIN_ADJUST_MAX_ITERATIONS": 1,
        "T_BIN_ADJUST_TOLERANCE": 0.01,
        "T_BIN_EDGE_BIAS": 1.0,
        "T_BIN_MIN_EVENTS": 1,
    }.items():
        setattr(fake_background_config, name, value)

    with mock.patch.dict(
        sys.modules,
        {
            "ltsep": fake_ltsep,
            "utility": fake_utility,
            "background_config": fake_background_config,
        },
    ):
        spec = importlib.util.spec_from_file_location("find_bins_test_module", FIND_BINS_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


class PreParticleSubtractionBinningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.find_bins = _load_find_bins()

    def test_pre_particle_subtraction_histograms_are_preferred(self):
        hist = {
            "phi_setting": "Center",
            "normfac_data": 1.0,
            "H_t_DATA": _FakeHistogram((10.0, 20.0, 30.0)),
            "H_ph_q_DATA": _FakeHistogram((1.0, 1.1, 1.2)),
            self.find_bins.PRE_PARTICLE_SUBTRACTION_T_HIST_KEY: _FakeHistogram((0.1, 0.2, 0.3)),
            self.find_bins.PRE_PARTICLE_SUBTRACTION_PHI_HIST_KEY: _FakeHistogram((0.1, 0.2, 0.3)),
        }

        t_values, phi_values = self.find_bins._collect_bin_samples(
            [hist],
            {"tmin": 0.0, "tmax": 0.4, "normfac_data": 1.0},
        )

        self.assertEqual(t_values.tolist(), [0.0, 0.2, 0.4])
        self.assertAlmostEqual(phi_values[1], 0.2 * (180.0 / math.pi))

    def test_final_histograms_remain_the_legacy_fallback(self):
        final_t = _FakeHistogram((0.1, 0.2, 0.3))
        final_phi = _FakeHistogram((0.1, 0.2, 0.3))
        selected_t, selected_phi = self.find_bins._get_binning_histograms(
            {"H_t_DATA": final_t, "H_ph_q_DATA": final_phi}
        )

        self.assertIs(selected_t, final_t)
        self.assertIs(selected_phi, final_phi)


if __name__ == "__main__":
    unittest.main()
