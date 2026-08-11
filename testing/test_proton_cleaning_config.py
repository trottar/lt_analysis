"""Regression checks for scoped proton-cleaning configuration."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))

import background_config


class ProtonCleaningConfigTests(unittest.TestCase):
    def test_q3p0_w2p32_requires_local_low_aero_timing_offset(self):
        config = background_config.get_proton_contamination_cleaning_config(
            inp_dict={"Q2": "3p0", "W": "2p32"},
        )

        self.assertFalse(
            config["low_aero_offset"]["allow_stable_global_center_fallback"]
        )
        self.assertIn(
            "Q3p0W2p32:default",
            config["proton_contamination_override_layers"],
        )

    def test_other_settings_preserve_the_stable_center_fallback(self):
        config = background_config.get_proton_contamination_cleaning_config(
            inp_dict={"Q2": "3p0", "W": "3p14"},
        )

        self.assertTrue(
            config["low_aero_offset"]["allow_stable_global_center_fallback"]
        )
        self.assertNotIn(
            "Q3p0W2p32:default",
            config["proton_contamination_override_layers"],
        )

    def test_partial_timing_t_nested_overrides_preserve_defaults(self):
        override = {
            "t_binning": {"edge_tolerance": 2.5e-10},
            "t_cell_fit": {"minimum_entries": 41},
            "t_support_thresholds": {"minimum_setting_valid_t_cells": 3},
            "aerogel_validation": {"enabled": False},
            "mm_diagnostics": {"display_bins": 64},
            "lambda_preservation_gate": {"maximum_lambda_removed_fraction": 0.08},
        }
        with mock.patch.dict(
            background_config.PROTON_CONTAMINATION_CLEANING_RUNTIME_OVERRIDES,
            override,
            clear=True,
        ):
            config = background_config.get_proton_contamination_cleaning_config(
                inp_dict={"Q2": "3p0", "W": "3p14"},
            )
        self.assertEqual(config["t_binning"]["edge_tolerance"], 2.5e-10)
        self.assertTrue(config["t_binning"]["require_metadata_sidecar"])
        self.assertEqual(config["t_cell_fit"]["minimum_entries"], 41)
        self.assertEqual(config["t_cell_fit"]["maximum_chi2_ndf"], 5.0)
        self.assertEqual(config["t_support_thresholds"]["minimum_setting_valid_t_cells"], 3)
        self.assertEqual(config["t_support_thresholds"]["minimum_modeled_yield"], 5.0)
        self.assertFalse(config["aerogel_validation"]["enabled"])
        self.assertFalse(config["aerogel_validation"]["affects_event_weights"])
        self.assertEqual(config["aerogel_validation"]["display_bins"], 100)
        self.assertEqual(
            tuple(config["aerogel_validation"]["summary_slice_edges"]),
            (0.0, 3.0, 6.0, 10.0, 15.0, 25.0),
        )
        self.assertEqual(config["mm_diagnostics"]["display_bins"], 64)
        self.assertTrue(config["mm_diagnostics"]["write_overview_page"])
        self.assertEqual(config["mm_diagnostics"]["max_t_panels_per_page"], 1)
        self.assertFalse(config["mm_diagnostics"]["affects_event_weights"])
        self.assertEqual(
            config["lambda_preservation_gate"]["maximum_lambda_removed_fraction"],
            0.08,
        )
        self.assertEqual(
            config["lambda_preservation_gate"]["validation_window_key"],
            "lambda_peak",
        )
        self.assertEqual(
            config["lambda_preservation_gate"]["minimum_raw_prompt_events"], 20
        )
        self.assertIsNone(config["lambda_preservation_gate"]["minimum_absolute_support"])


if __name__ == "__main__":
    unittest.main()
