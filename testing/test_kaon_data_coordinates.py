"""Unit coverage for the shared kaon-derived experimental coordinate frame."""

from __future__ import annotations

import importlib.util
import sys
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CUTS_PATH = REPO_ROOT / "src" / "cuts"
BINNING_PATH = REPO_ROOT / "src" / "binning"
COORDINATES_PATH = CUTS_PATH / "data_coordinates.py"
APPLY_CUTS_PATH = CUTS_PATH / "apply_cuts.py"
CANONICAL_BINNING_PATH = REPO_ROOT / "src" / "utility" / "canonical_binning.py"


class _Event:
    MM = 1.08
    MandelT = -0.40
    # Deliberately wrong legacy branches: they must never select the frame.
    MM_shift = 9.99
    t_shift = 8.88


def _load_coordinates():
    spec = importlib.util.spec_from_file_location("data_coordinates", COORDINATES_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["data_coordinates"] = module
    spec.loader.exec_module(module)
    return module


def _load_apply_cuts():
    sys.path.insert(0, str(CUTS_PATH))
    sys.path.insert(0, str(BINNING_PATH))
    fake_theta = types.ModuleType("theta_cm")
    fake_theta.calculate_tmin = lambda *_args, **_kwargs: 0.0
    previous = sys.modules.get("theta_cm")
    sys.modules["theta_cm"] = fake_theta
    try:
        spec = importlib.util.spec_from_file_location("apply_cuts_coordinate_test", APPLY_CUTS_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if previous is None:
            sys.modules.pop("theta_cm", None)
        else:
            sys.modules["theta_cm"] = previous


class KaonDataCoordinateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates = _load_coordinates()
        cls.apply_cuts = _load_apply_cuts()

    def _contract(self):
        return self.coordinates.build_kaon_data_coordinate_contract(
            "Left",
            {"shift": 0.012, "data_peak": 1.104, "simc_peak": 1.116},
            {"shift": -0.031, "theta_cm_deg": 18.0},
        )

    def test_contract_is_deterministic_and_rejects_tampering(self):
        first = self._contract()
        second = self._contract()
        self.assertEqual(first["coordinate_fingerprint"], second["coordinate_fingerprint"])

        tampered = dict(first, mm_shift=0.013)
        with self.assertRaisesRegex(RuntimeError, "fingerprint"):
            self.coordinates.validate_kaon_data_coordinate_contract(tampered)

    def test_raw_and_analysis_coordinates_ignore_legacy_root_branches(self):
        result = self.coordinates.analysis_event_coordinates(_Event(), self._contract())
        self.assertEqual(result["raw_mm"], 1.08)
        self.assertEqual(result["raw_t"], 0.40)
        self.assertAlmostEqual(result["analysis_mm"], 1.092)
        self.assertAlmostEqual(result["analysis_t"], 0.369)

    def test_apply_cuts_uses_one_contract_and_disables_legacy_offset(self):
        contract = self._contract()
        self.apply_cuts.set_shift_context(
            phi_setting="Left",
            shift_mode="raw",
            kaon_data_coordinate=contract,
            kaon_data_coordinate_summary={"Left": contract},
        )
        self.assertEqual(self.apply_cuts.get_shift_mode(), "shifted")
        self.assertEqual(self.apply_cuts.get_effective_mm_offset(0.25), 0.0)
        self.assertAlmostEqual(self.apply_cuts.get_shifted_mm(_Event(), mm_offset=0.25), 1.092)
        self.assertAlmostEqual(self.apply_cuts.get_shifted_t(_Event()), 0.369)

    def test_final_canonical_edge_is_inclusive(self):
        spec = importlib.util.spec_from_file_location("canonical_binning_coordinate_test", CANONICAL_BINNING_PATH)
        canonical = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(canonical)
        self.assertEqual(canonical.find_canonical_bin(0.4, (0.0, 0.2, 0.4)), 1)

    def test_missing_t_shift_is_fatal_for_canonical_coordinate_contract(self):
        with self.assertRaisesRegex(RuntimeError, "missing_t_shift"):
            self.coordinates.build_kaon_data_coordinate_contract(
                "Left", {"shift": 0.01}, None, require_t_shift=True
            )


if __name__ == "__main__":
    unittest.main()
