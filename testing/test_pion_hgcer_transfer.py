"""Pure-Python contracts for the Part-2 analysis-stage zero-PE map."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

from background_config import (
    fingerprint_hgcer_pid_contract,
    get_pion_hgcer_transfer_config,
    hgcer_mask_accepts,
    normalize_hgcer_mask,
)
import pion_hgcer_transfer as transfer

try:
    import ROOT
except ImportError:  # pragma: no cover - host capability
    ROOT = None


def _config(**override):
    config = get_pion_hgcer_transfer_config({
        "pion_hgcer_transfer_config": {
            "toy_count": 40,
            "minimum_prompt_pion_records": 20,
            **override,
        }
    })
    return config


def _record(side, npe, t=0, delta=0, source="prompt"):
    return {
        "side": side, "P_hgcer_npeSum": float(npe), "canonical_t_index": t,
        "delta_index": delta, "source_label": source, "allcuts": True,
        "coefficient": 1.0, "adj_MM": 1.1,
    }


def _part1(records, cells=None, no_rf=True):
    source_name = "T_pion_noRF" if no_rf else "T_pion_RF"
    return {
        "status": "available", "coordinate_fingerprint": "coordinate-test",
        "t_edges": [0.1, 0.2], "delta_edges": [-10.0, 0.0, 10.0],
        "source_provenance": {
            "kaon": {"prompt": {"tree_name": "T_kaon_noRF"}},
            "pion": {"prompt": {"tree_name": source_name}},
        },
        "records": records,
        "cells": cells or [
            {"t_index": 0, "delta_index": 0, "support_class": "supported"},
            {"t_index": 0, "delta_index": 1, "support_class": "supported"},
        ],
    }


class PionHGCerTransferPurePythonTests(unittest.TestCase):
    def test_mask_contract_boundaries_and_fingerprint_are_stable(self):
        contract = fingerprint_hgcer_pid_contract()
        self.assertTrue(hgcer_mask_accepts(0.0, contract["masks"]["kaon_tree"]))
        self.assertFalse(hgcer_mask_accepts(1.0e-9, contract["masks"]["kaon_tree"]))
        self.assertTrue(hgcer_mask_accepts(1.0e-9, contract["masks"]["pion_tree"]))
        self.assertFalse(hgcer_mask_accepts(2.0, contract["masks"]["physical_pion_control"]))
        self.assertTrue(hgcer_mask_accepts(2.000001, contract["masks"]["physical_pion_control"]))
        self.assertEqual(contract, fingerprint_hgcer_pid_contract())
        self.assertEqual(normalize_hgcer_mask({"field": "P_hgcer_npeSum", "operator": "<=", "value": 2}), {"field": "P_hgcer_npeSum", "operator": "<=", "value": 2.0})
        with self.assertRaises(ValueError):
            fingerprint_hgcer_pid_contract({
                "kaon_tree": {"field": "P_hgcer_npeSum", "operator": "==", "value": 0.0},
                "pion_tree": {"field": "P_hgcer_npeSum", "operator": ">", "value": 0.0},
                "physical_pion_control": {"field": "P_hgcer_npeSum", "operator": ">", "value": 1.0},
            })

    def test_setting_wide_fallback_is_rejected_by_configuration(self):
        with self.assertRaisesRegex(ValueError, "setting-wide fallback"):
            get_pion_hgcer_transfer_config({
                "pion_hgcer_transfer_config": {"fallback": {"setting_wide_enabled": True}}
            })

    def test_poisson_recovery_and_actual_control_denominator(self):
        rng = np.random.default_rng(20260827)
        true_mu = 3.2
        values = rng.poisson(true_mu, size=3000)
        positive = values[values > 0]
        config = _config()
        fit = transfer.fit_zero_photoelectron_response(
            positive, response_mask=config["pid_contract"]["masks"]["pion_tree"],
            contract=config["pid_contract"], config=config,
        )
        self.assertIn(fit["fit_status"], {"fit_valid", "fit_valid_but_P0_weakly_constrained"})
        self.assertAlmostEqual(fit["mu"], true_mu, delta=0.20)
        self.assertAlmostEqual(fit["P0"], np.exp(-true_mu), delta=0.02)
        above_zero = transfer.zero_truncated_poisson_transfer(
            true_mu, config["pid_contract"]["masks"]["kaon_tree"], config["pid_contract"]["masks"]["pion_tree"],
        )
        above_two = transfer.zero_truncated_poisson_transfer(
            true_mu, config["pid_contract"]["masks"]["kaon_tree"], config["pid_contract"]["masks"]["physical_pion_control"],
        )
        self.assertGreater(above_two, above_zero)
        self.assertGreater(fit["accepted_toy_fraction"], 0.0)

    def test_positive_only_gamma_is_rejected_for_zero_mask_numerator(self):
        with self.assertRaisesRegex(ValueError, "cannot determine"):
            transfer.validate_zero_photoelectron_response_family("positive_gamma")
        self.assertEqual(
            transfer.validate_zero_photoelectron_response_family("zero_truncated_poisson"),
            "zero_truncated_poisson",
        )

    def test_norf_record_and_cache_contract_fail_closed_without_filtering(self):
        config = _config()
        records = {
            "kaon": [_record("kaon", 0.0)],
            "pion": [_record("pion", 1.0)],
        }
        bad_rf = transfer.audit_pion_hgcer_transfer_inputs(
            _part1(records, no_rf=False),
            {"physical_pion_control_mask_fingerprint": config["pid_contract"]["fingerprint"]},
            config["pid_contract"],
        )
        self.assertEqual(bad_rf["status"], "unavailable")
        self.assertEqual(bad_rf["reason"], "norf_provenance_failed")
        bad_content = transfer.audit_pion_hgcer_transfer_inputs(
            _part1({"kaon": [_record("kaon", 1.0)], "pion": [_record("pion", 1.0)]}),
            {"physical_pion_control_mask_fingerprint": config["pid_contract"]["fingerprint"]},
            config["pid_contract"],
        )
        self.assertEqual(bad_content["reason"], "pid_tree_event_content_audit_failed")
        bad_cache = transfer.audit_pion_hgcer_transfer_inputs(
            _part1(records), {"physical_pion_control_mask_fingerprint": "wrong"}, config["pid_contract"],
        )
        self.assertEqual(bad_cache["reason"], "pion_control_cache_mask_contract_mismatch")

    def test_same_t_only_fallback_levels_and_unsupported_cells(self):
        config = _config()
        cells = {
            (0, 0): {"support_class": "supported"}, (0, 1): {"support_class": "supported"},
            (0, 2): {"support_class": "supported"}, (0, 3): {"support_class": "supported"},
            (0, 4): {"support_class": "unsupported"}, (1, 0): {"support_class": "supported"},
        }
        fits = {
            (0, 1): {"fit_status": "fit_valid", "transfer": 1.0, "transfer_total_uncertainty": 0.1},
            (0, 3): {"fit_status": "fit_valid", "transfer": 3.0, "transfer_total_uncertainty": 0.2},
        }
        samples = {key: [1.0] * 25 for key in cells}
        pooled = {"fit_status": "fit_valid", "transfer": 4.0, "transfer_total_uncertainty": 0.4}
        with patch.object(transfer, "fit_zero_photoelectron_response", return_value=pooled):
            result = transfer._solve_fallbacks(cells, fits, samples, config["pid_contract"], config)
        self.assertEqual(result[(0, 2)]["solution_source"], "same_t_delta_bracket")
        self.assertEqual(result[(0, 0)]["solution_source"], "same_t_delta_edge")
        self.assertEqual(result[(0, 4)]["solution_source"], "unsupported")
        self.assertEqual(result[(1, 0)]["solution_source"], "same_t_pooled")
        self.assertNotEqual(result[(1, 0)]["solution_source"], "same_t_delta_edge")

    def test_map_is_frozen_and_weak_direct_fit_does_not_seed_direct_solution(self):
        config = _config()
        rng = np.random.default_rng(9)
        pions = [_record("pion", value, delta=0) for value in rng.poisson(3.0, size=80) if value > 0]
        part1 = _part1({"kaon": [_record("kaon", 0.0)], "pion": pions}, cells=[
            {"t_index": 0, "delta_index": 0, "support_class": "supported"},
        ])
        cache = {"physical_pion_control_mask_fingerprint": config["pid_contract"]["fingerprint"]}
        payload = transfer.build_pion_hgcer_zerope_transfer_map(part1, cache, config=config)
        self.assertTrue(payload["frozen"])
        self.assertTrue(payload["map_fingerprint"])
        self.assertIn(payload["cells"][0]["fit"]["fit_status"], {"fit_valid", "fit_valid_but_P0_weakly_constrained", "fit_invalid"})
        if payload["cells"][0]["fit"]["fit_status"] != "fit_valid":
            self.assertNotEqual(payload["cells"][0]["solution"]["solution_source"], "direct")

    def test_display_catalog_is_ascii_and_calls_host_norf(self):
        for value in transfer.ROOT_SAFE_TRANSFER_LABELS.values():
            self.assertTrue(value.isascii())
            self.assertNotIn("pre-RF", value)
        unavailable = transfer.apply_frozen_pion_hgcer_transfer_map({}, {}, ())
        self.assertEqual(unavailable["status"], "unavailable")


@unittest.skipUnless(ROOT is not None, "PyROOT unavailable")
class PionHGCerTransferPyROOTTests(unittest.TestCase):
    def test_frozen_application_uses_detached_norf_host_and_serializes_manifest(self):
        host = ROOT.TH1D("part2_host", "noRF host", 20, 0.7, 1.5)
        host.SetDirectory(0)
        host.Fill(1.10, 10.0)
        payload = {
            "status": "available", "frozen": True, "map_fingerprint": "map-test",
            "Part3_eligibility": "review_required_incomplete_uncertainty_or_solution",
            "cells": [{"t_index": 0, "delta_index": 0, "solution": {"transfer": 0.25, "solution_source": "direct"}}],
        }
        cache = {"by_t": ({"records": ({
            "allcuts": True, "delta_index": 0, "adj_MM": 1.10, "coefficient": 2.0,
        },)},)}
        before = float(host.Integral())
        application = transfer.apply_frozen_pion_hgcer_transfer_map(
            payload, cache, ({"cleaned_targets_pre_rf": {"h_mm_nosub": host}},),
        )
        self.assertEqual(application["status"], "available")
        self.assertEqual(application["host_label"], "proton-cleaned noRF host")
        self.assertEqual(float(host.Integral()), before)
        self.assertAlmostEqual(application["t_products"][0]["pion_integral"], 0.5)
        payload["application"] = application
        manifest = transfer.expected_pion_hgcer_transfer_page_manifest(payload)
        self.assertEqual(manifest[-1], "part2_proposed_noRF_host_mm")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "part2.json"
            transfer.write_pion_hgcer_zerope_transfer_json(path, payload)
            self.assertIn("noRF", path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
