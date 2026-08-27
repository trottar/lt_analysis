"""Pure-Python contracts for the Part-2 zero-photoelectron transfer map."""

from __future__ import annotations

import math
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_transfer as transfer
from background_config import get_pion_hgcer_transfer_config


def _mask(operator, value):
    return {"field": "P_hgcer_npeSum", "operator": operator, "value": value}


def _config(**overrides):
    config = {
        "response_family": "auto",
        "fit_range": (0.0, 20.0),
        "integer_tolerance": 1.0e-6,
        "minimum_prompt_pion_records": 20,
        "toy_count": 160,
        "minimum_accepted_toy_fraction": 0.50,
        "profile_delta_nll": 0.5,
        "profile_scan_log_mu_half_width": 3.0,
        "profile_scan_points": 25,
        "correlation_condition_number_max": 1.0e8,
        "pair_correlation_abs_max": 0.995,
        "fallback": {
            "delta_interpolation_relative_uncertainty": 0.20,
            "edge_neighbor_relative_uncertainty": 0.30,
            "minimum_pooled_prompt_pion_records": 20,
            "allow_setting_wide": False,
        },
    }
    config.update(overrides)
    return config


def _manifest(tree_name, expression):
    pid_role = "kaon_pid" if "Kaon" in tree_name else "pion_pid"
    return {
        "trees": {
            tree_name: {
                "tree_name": tree_name,
                "cut_name": "synthetic",
                "cut_file": "synthetic.cuts",
                "expression": expression,
                "pid_role": pid_role,
                "rf_state": "noRF",
                "fingerprint": "synthetic-{}".format(tree_name),
            }
        }
    }


def _source_entries():
    kaon_tree = "Cut_Kaon_Events_prompt_noRF"
    pion_tree = "Cut_Pion_Events_prompt_noRF"
    return {
        "kaon": [{"source_label": "prompt", "tree_name": kaon_tree, "manifest": _manifest(kaon_tree, "P_hgcer_npeSum == 0.0"), "rf_state": "noRF", "pid_role": "kaon_pid", "signed_coefficient": 1.0, "coordinate_fingerprint": "coords", "proton_factor_scope": "kaon_host_only"}],
        "pion": [{"source_label": "prompt", "tree_name": pion_tree, "manifest": _manifest(pion_tree, "P_hgcer_npeSum > 0.0"), "rf_state": "noRF", "pid_role": "pion_pid", "signed_coefficient": 1.0, "coordinate_fingerprint": "coords", "proton_factor_scope": "never_pion"}],
    }


class PionHGCerTransferTests(unittest.TestCase):
    def test_config_is_diagnostic_and_validated(self):
        config = get_pion_hgcer_transfer_config({"pion_hgcer_transfer_config": {"toy_count": 31}})
        self.assertTrue(config["enabled"])
        self.assertEqual(config["toy_count"], 31)
        self.assertEqual(config["minimum_prompt_pion_records"], 20)
        with self.assertRaisesRegex(ValueError, "response_family"):
            get_pion_hgcer_transfer_config({"pion_hgcer_transfer_config": {"response_family": "gamma"}})

    def test_exact_tree_masks_and_control_mask_are_distinct(self):
        resolved = transfer.resolve_hgcer_transfer_masks(_source_entries(), _mask(">", 2.0))
        self.assertEqual(resolved["masks"]["S_K_tree"]["operator"], "==")
        self.assertEqual(resolved["masks"]["S_K_tree"]["value"], 0.0)
        self.assertEqual(resolved["masks"]["S_pi_tree"]["operator"], ">")
        self.assertEqual(resolved["masks"]["S_pi_physical_control"]["value"], 2.0)
        self.assertEqual(resolved["response_sample_provenance"]["source_label"], "prompt")
        self.assertEqual(resolved["response_sample_provenance"]["weighting"], "unweighted")
        self.assertTrue(resolved["pid_selection_fingerprint"])

    def test_missing_or_mismatched_provenance_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "missing persisted"):
            transfer.resolve_hgcer_transfer_masks(
                {"kaon": [{"tree_name": "missing_noRF", "manifest": {}, "rf_state": "noRF", "pid_role": "kaon_pid"}], "pion": []},
                _mask(">", 2.0),
            )
        entries = _source_entries()
        pion_tree = "Cut_Pion_Events_prompt_noRF"
        entries["pion"].append({"source_label": "rand", "tree_name": pion_tree, "manifest": _manifest(pion_tree, "P_hgcer_npeSum > 2.0"), "rf_state": "noRF", "pid_role": "pion_pid"})
        with self.assertRaisesRegex(ValueError, "disagree"):
            transfer.resolve_hgcer_transfer_masks(entries, _mask(">", 2.0))

    def test_non_norf_or_wrong_role_provenance_fails_closed(self):
        entries = _source_entries()
        entries["pion"][0]["rf_state"] = "RF"
        with self.assertRaisesRegex(ValueError, "noRF"):
            transfer.resolve_hgcer_transfer_masks(entries, _mask(">", 2.0))
        entries = _source_entries()
        entries["kaon"][0]["pid_role"] = "pion_pid"
        with self.assertRaisesRegex(ValueError, "PID role"):
            transfer.resolve_hgcer_transfer_masks(entries, _mask(">", 2.0))

    def test_zero_truncated_poisson_recovers_p0_and_transfer(self):
        rng = np.random.default_rng(24)
        truth_mu = 1.8
        values = rng.poisson(truth_mu, 6000)
        values = values[values > 0]
        result = transfer.fit_zero_truncated_pion_response(
            values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            physical_control_mask=_mask(">", 2.0),
            config=_config(),
            fingerprint_context="poisson_recovery",
        )
        self.assertEqual(result["fit_status"], "fit_valid")
        self.assertEqual(result["response_family"], "zero_truncated_poisson")
        self.assertAlmostEqual(result["parameters"]["mu"], truth_mu, delta=0.12)
        self.assertAlmostEqual(result["P0"], math.exp(-truth_mu), delta=0.025)
        truth_transfer = math.exp(-truth_mu) / float(1.0 - poisson_cdf(2, truth_mu))
        self.assertAlmostEqual(result["R_pi_to_K"], truth_transfer, delta=0.06)
        self.assertEqual(result["profile_likelihood"]["status"], "two_sided_68pct")
        self.assertGreaterEqual(result["accepted_toy_fraction"], 0.50)

    def test_control_masks_produce_different_transfers_and_host_counts_do_not_enter(self):
        rng = np.random.default_rng(25)
        values = rng.poisson(1.6, 4000)
        values = values[values > 0]
        common = dict(
            values=values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            config=_config(),
            fingerprint_context="mask_difference",
        )
        over_zero = transfer.fit_zero_truncated_pion_response(physical_control_mask=_mask(">", 0.0), **common)
        over_two = transfer.fit_zero_truncated_pion_response(physical_control_mask=_mask(">", 2.0), **common)
        self.assertNotEqual(over_zero["R_pi_to_K"], over_two["R_pi_to_K"])
        self.assertEqual(over_zero["R_pi_to_K"], over_zero["epsilon_pi_to_K"] / over_zero["epsilon_pi_to_pi_control"])
        self.assertEqual(over_two["R_pi_to_K"], over_two["epsilon_pi_to_K"] / over_two["epsilon_pi_to_pi_control"])

    def test_p0_identifiability_weak_status_does_not_need_relative_threshold(self):
        rng = np.random.default_rng(26)
        values = rng.poisson(1.4, 2000)
        values = values[values > 0]
        result = transfer.fit_zero_truncated_pion_response(
            values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            physical_control_mask=_mask(">", 2.0),
            config=_config(profile_scan_log_mu_half_width=0.001),
            fingerprint_context="weak_profile",
        )
        self.assertEqual(result["fit_status"], "fit_valid_but_P0_weakly_constrained")
        self.assertIsNotNone(result["P0_relative_uncertainty"])
        self.assertNotEqual(result["profile_likelihood"]["status"], "two_sided_68pct")

    def test_continuous_response_has_a_structural_zero_atom(self):
        rng = np.random.default_rng(27)
        truth_mu = 1.35
        photoelectron_counts = rng.poisson(truth_mu, 700)
        values = np.asarray([
            rng.gamma(shape=2.0 * count, scale=0.45)
            for count in photoelectron_counts if count > 0
        ])
        result = transfer.fit_zero_truncated_pion_response(
            values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            physical_control_mask=_mask(">", 2.0),
            config=_config(toy_count=60, profile_scan_points=15),
            fingerprint_context="compound_response",
        )
        self.assertIn(result["response_family"], {"zero_truncated_compound_poisson_gamma", "zero_truncated_compound_poisson_exponential"})
        self.assertAlmostEqual(result["parameters"]["mu"], truth_mu, delta=0.45)
        self.assertAlmostEqual(result["P0"], math.exp(-truth_mu), delta=0.14)
        self.assertNotEqual(result["response_family"], "gamma")

    def test_positive_only_gamma_is_rejected_as_a_leakage_model(self):
        values = np.asarray([0.5, 1.0, 1.5] * 8)
        result = transfer.fit_zero_truncated_pion_response(
            values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            physical_control_mask=_mask(">", 2.0),
            config=_config(response_family="gamma"),
            fingerprint_context="reject_positive_only_gamma",
        )
        self.assertEqual(result["fit_status"], "fit_invalid")
        self.assertIn("response_family", result["reason"])

    def test_negative_binomial_alternate_preserves_the_zero_atom_parameterization(self):
        parameters = transfer._model_kind_parameters(
            np.log([1.7, 2.5]), "zero_truncated_negative_binomial"
        )
        self.assertAlmostEqual(
            transfer._model_p0(parameters, "zero_truncated_negative_binomial"),
            math.exp(-1.7),
        )

    def test_covariance_degeneracy_marks_a_fit_weak(self):
        rng = np.random.default_rng(28)
        counts = rng.poisson(1.3, 300)
        values = np.asarray([
            rng.gamma(shape=2.0 * count, scale=0.5)
            for count in counts if count > 0
        ])
        result = transfer.fit_zero_truncated_pion_response(
            values,
            response_mask=_mask(">", 0.0),
            kaon_mask=_mask("==", 0.0),
            physical_control_mask=_mask(">", 2.0),
            config=_config(toy_count=60, profile_scan_points=15, pair_correlation_abs_max=0.10),
            fingerprint_context="weak_covariance",
        )
        self.assertEqual(result["fit_status"], "fit_valid_but_P0_weakly_constrained")
        self.assertTrue(result["covariance_diagnostics"]["near_degenerate"])

    def test_low_toy_acceptance_is_invalid(self):
        rng = np.random.default_rng(29)
        values = rng.poisson(1.5, 1000)
        values = values[values > 0]
        with mock.patch.object(
            transfer, "_response_toys", return_value=np.asarray([[0.1, 0.1, 0.5, 0.2]])
        ):
            result = transfer.fit_zero_truncated_pion_response(
                values,
                response_mask=_mask(">", 0.0),
                kaon_mask=_mask("==", 0.0),
                physical_control_mask=_mask(">", 2.0),
                config=_config(toy_count=60),
                fingerprint_context="low_toy_acceptance",
            )
        self.assertEqual(result["fit_status"], "fit_invalid")
        self.assertEqual(result["reason"], "accepted_toy_fraction_below_minimum")

    def test_map_freezes_and_interpolates_only_within_t(self):
        cells = []
        records = []
        for delta_index in range(3):
            cells.append({
                "t_index": 0, "t_low": 0.4, "t_high": 0.9,
                "delta_index": delta_index, "delta_low": float(delta_index), "delta_high": float(delta_index + 1),
                "support_class": "supported", "kaon": {}, "pion": {},
            })
            if delta_index != 1:
                records.extend({"source_label": "prompt", "allcuts": True, "P_hgcer_npeSum": 1.0, "canonical_t_index": 0, "delta_index": delta_index} for _ in range(20))
        part1 = {
            "status": "available", "coordinate_fingerprint": "coords", "t_edges": [0.4, 0.9], "delta_edges": [0.0, 1.0, 2.0, 3.0],
            "cells": cells, "records": {"pion": records},
        }
        direct = {
            "fit_status": "fit_valid", "R_pi_to_K": 0.2, "sigma_stat": 0.01, "sigma_model": 0.02,
            "P0": 0.1, "P0_relative_uncertainty": 0.1, "epsilon_pi_to_K": 0.1, "epsilon_pi_to_pi_control": 0.5,
        }
        with mock.patch.object(transfer, "fit_zero_truncated_pion_response", side_effect=lambda values, **kwargs: dict(direct) if len(values) >= 20 else {"fit_status": "fit_invalid", "reason": "short"}):
            result = transfer.build_pion_hgcer_transfer_map(
                part1_payload=part1,
                pion_control_cache={"physical_control_hgcer_mask": _mask(">", 2.0)},
                source_entries=_source_entries(),
                config=_config(),
            )
        middle = [cell for cell in result["cells"] if cell["delta_index"] == 1][0]
        self.assertEqual(middle["solution_source"], "delta_interpolation")
        self.assertAlmostEqual(middle["R_pi_to_K"], 0.2)
        with self.assertRaises(TypeError):
            result["status"] = "mutated"

    def test_t_pooled_fallback_never_crosses_canonical_t_or_setting(self):
        cells = []
        records = []
        for t_index in range(2):
            for delta_index in range(2):
                cells.append({
                    "t_index": t_index, "t_low": 0.4 + 0.5 * t_index, "t_high": 0.9 + 0.5 * t_index,
                    "delta_index": delta_index, "delta_low": float(delta_index), "delta_high": float(delta_index + 1),
                    "support_class": "supported", "kaon": {}, "pion": {},
                })
                count = 10 if t_index == 0 else 5
                records.extend({"source_label": "prompt", "allcuts": True, "P_hgcer_npeSum": 1.0, "canonical_t_index": t_index, "delta_index": delta_index} for _ in range(count))
        part1 = {
            "status": "available", "coordinate_fingerprint": "coords", "t_edges": [0.4, 0.9, 1.4], "delta_edges": [0.0, 1.0, 2.0],
            "cells": cells, "records": {"pion": records},
        }
        direct = {"fit_status": "fit_valid", "R_pi_to_K": 0.3, "sigma_stat": 0.01, "sigma_model": 0.02, "P0": 0.1, "P0_relative_uncertainty": 0.1, "epsilon_pi_to_K": 0.1, "epsilon_pi_to_pi_control": 0.4}
        with mock.patch.object(transfer, "fit_zero_truncated_pion_response", side_effect=lambda values, **kwargs: dict(direct) if len(values) >= 20 else {"fit_status": "fit_invalid", "reason": "short"}):
            result = transfer.build_pion_hgcer_transfer_map(
                part1_payload=part1,
                pion_control_cache={"physical_control_hgcer_mask": _mask(">", 2.0)},
                source_entries=_source_entries(),
                config=_config(),
            )
        t0_sources = {cell["solution_source"] for cell in result["cells"] if cell["t_index"] == 0}
        t1_sources = {cell["solution_source"] for cell in result["cells"] if cell["t_index"] == 1}
        self.assertEqual(t0_sources, {"t_pooled"})
        self.assertEqual(t1_sources, {"unsupported"})

    def test_proposed_application_never_uses_proton_factor(self):
        transfer_map = {
            "status": "available", "transfer_map_fingerprint": "frozen", "delta_edges": [0.0, 1.0],
            "cells": [{"t_index": 0, "delta_index": 0, "solution_source": "direct", "R_pi_to_K": 0.25}],
        }
        application = transfer.apply_pion_hgcer_transfer_diagnostic(
            transfer_map,
            {"records": [{"allcuts": True, "t_index": 0, "ssdelta": 0.5, "adj_MM": 1.1, "coefficient": -4.0, "proton_cleaning_factor": 0.01}]},
        )
        self.assertEqual(application["per_t"][0]["pion_estimate_integral"], -1.0)
        self.assertIsNone(application["records"][0]["proton_cleaning_factor"])
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / "transfer.json"
            transfer.write_pion_hgcer_transfer_json(
                sidecar, transfer_map, application=application,
                page_manifest=[{"page_id": "pion_hgcer.part2.status"}],
            )
            payload = json.loads(sidecar.read_text(encoding="utf-8"))
        self.assertEqual(payload["proposed_application"]["per_t"]["0"]["pion_estimate_integral"], -1.0)
        self.assertEqual(payload["page_manifest"][0]["page_id"], "pion_hgcer.part2.status")

    @unittest.skipIf(transfer.ROOT is None, "PyROOT is unavailable")
    def test_pyroot_detached_source_to_map_to_render_contract(self):
        root = transfer.ROOT
        host = root.TH1D("H_part2_host", "host", 10, 0.8, 1.3)
        pion_control = root.TH1D("H_part2_pion_control", "pion", 10, 0.8, 1.3)
        host.Fill(1.1, 10.0)
        part1 = {
            "status": "available", "coordinate_fingerprint": "coords", "t_edges": [0.4, 0.9], "delta_edges": [0.0, 1.0],
            "cells": [{"t_index": 0, "t_low": 0.4, "t_high": 0.9, "delta_index": 0, "delta_low": 0.0, "delta_high": 1.0, "support_class": "supported", "kaon": {}, "pion": {}}],
            "records": {"pion": [{"source_label": "prompt", "allcuts": True, "P_hgcer_npeSum": 1.0, "canonical_t_index": 0, "delta_index": 0} for _ in range(20)]},
        }
        direct = {"fit_status": "fit_valid", "response_family": "zero_truncated_poisson", "R_pi_to_K": 0.25, "sigma_stat": 0.01, "sigma_model": 0.02, "P0": 0.1, "P0_relative_uncertainty": 0.2, "epsilon_pi_to_K": 0.1, "epsilon_pi_to_pi_control": 0.4}
        with mock.patch.object(transfer, "fit_zero_truncated_pion_response", return_value=direct):
            transfer_map = transfer.build_pion_hgcer_transfer_map(
                part1_payload=part1,
                pion_control_cache={"physical_control_hgcer_mask": _mask(">", 2.0)},
                source_entries=_source_entries(),
                config=_config(),
            )
        application = transfer.apply_pion_hgcer_transfer_diagnostic(
            transfer_map,
            {"records": [{"allcuts": True, "t_index": 0, "ssdelta": 0.5, "adj_MM": 1.1, "coefficient": 4.0}], "by_t": [{"t_index": 0, "H_pion_control_cut": pion_control}]},
            {"canonical_t_products": ({"t_index": 0, "final_targets": {"h_mm": host}},), "final_targets": {"h_mm": host}},
        )
        self.assertAlmostEqual(host.Integral(), 10.0)
        self.assertIn("_root_products", application)
        self.assertTrue(application["mm_products"]["global"]["strict_global_sums"])
        with tempfile.TemporaryDirectory() as directory:
            pdf = str(Path(directory) / "part2.pdf")
            root.TCanvas("C_part2_open", "open", 1, 1).Print(pdf + "(")
            manifest = transfer.render_pion_hgcer_transfer_pages(pdf, transfer_map, application=application, close_pdf=True)
        page_ids = [page["page_id"] for page in manifest]
        self.assertEqual(page_ids, [page["page_id"] for page in transfer.expected_pion_hgcer_transfer_page_manifest(transfer_map)])


def poisson_cdf(value, mu):
    return sum(math.exp(-mu) * mu ** count / math.factorial(count) for count in range(int(value) + 1))


if __name__ == "__main__":
    unittest.main()
