"""Regression checks for scoped proton-cleaning configuration."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
