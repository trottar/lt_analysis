"""Central production-analysis configuration regressions."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
UTILITY_DIR = REPO_ROOT / "src" / "utility"
if str(UTILITY_DIR) not in sys.path:
    sys.path.insert(0, str(UTILITY_DIR))

from background_config import (  # noqa: E402
    ANALYSIS_RUNTIME_CONFIG_SENTINEL,
    resolve_analysis_runtime_config,
)


class AnalysisRuntimeConfigTests(unittest.TestCase):
    def test_production_settings_preserve_current_values(self):
        expected = {
            ("0p4", "2p20"): (5, 16, 0.001, 0.035),
            ("2p1", "2p95"): (4, 8, 0.150, 0.400),
            ("3p0", "2p32"): (4, 9, 0.400, 0.800),
            ("3p0", "3p14"): (4, 10, 0.180, 0.600),
            ("4p4", "2p74"): (4, 9, 0.400, 0.900),
            ("5p5", "3p02"): (2, 7, 0.400, 0.900),
        }
        for key, values in expected.items():
            with self.subTest(key=key):
                config = resolve_analysis_runtime_config(*key)
                self.assertEqual(
                    (config["num_t_bins"], config["num_phi_bins"], config["tmin"], config["tmax"]),
                    values,
                )
                self.assertEqual(config["mm_min"], 1.10)
                self.assertEqual(config["mm_max"], 1.16)
                self.assertEqual(config["pion_subtraction_scope"], "t_bin")
                self.assertFalse(config["pion_parent_diagnostic_strict"])

    def test_test_setting_is_explicit_and_unknown_setting_is_rejected(self):
        test_config = resolve_analysis_runtime_config("0p5", "2p40")
        self.assertEqual(test_config["source"], "explicit_test_setting")
        with self.assertRaisesRegex(ValueError, "No authoritative analysis runtime configuration"):
            resolve_analysis_runtime_config("9p9", "9p9")

    def test_launcher_compatibility_sentinel_is_stable(self):
        self.assertEqual(ANALYSIS_RUNTIME_CONFIG_SENTINEL, "__CONFIG__")


if __name__ == "__main__":
    unittest.main()
