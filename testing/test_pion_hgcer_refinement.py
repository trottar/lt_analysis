"""Pure-Python contracts for detached Part-2B HGCer relative refinement."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

from background_config import (  # noqa: E402
    fingerprint_hgcer_pid_contract,
    get_pion_hgcer_refinement_config,
    get_pion_hgcer_transfer_config,
)
import pion_hgcer_refinement as refinement  # noqa: E402
import pion_hgcer_transfer as transfer  # noqa: E402

try:  # pragma: no cover - host capability
    import ROOT
except ImportError:  # pragma: no cover - host capability
    ROOT = None


def _config(**override):
    return get_pion_hgcer_refinement_config({"pion_hgcer_refinement_config": override})


def _raw_transfer(t_count=1, delta_count=3):
    contract = fingerprint_hgcer_pid_contract()
    cells = []
    for t_index in range(t_count):
        for delta_index in range(delta_count):
            cells.append({
                "t_index": t_index, "delta_index": delta_index,
                "support_class": "supported",
                "fit": {"fit_status": "fit_valid", "model_family": "zero_truncated_poisson", "fit_parameters_log": [1.0], "response_mask": contract["masks"]["pion_tree"]},
                "solution": {"solution_source": "direct", "transfer": 0.2 + 0.1 * delta_index, "transfer_statistical_uncertainty": 0.02, "transfer_model_uncertainty": 0.01},
            })
    return {
        "status": "available", "frozen": True, "map_fingerprint": "raw-map-fingerprint",
        "pid_contract": contract, "coordinate_fingerprint": "coordinate-fingerprint",
        "t_edges": list(range(t_count + 1)), "delta_edges": list(range(delta_count + 1)),
        "config": get_pion_hgcer_transfer_config({"pion_hgcer_transfer_config": {"toy_count": 20}}),
        "cells": cells,
    }


def _part1_for_raw(raw):
    records = []
    for cell in raw["cells"]:
        for number in range(25):
            records.append({
                "source_label": "prompt", "allcuts": True,
                "canonical_t_index": cell["t_index"], "delta_index": cell["delta_index"],
                "P_hgcer_npeSum": float(max(1, (number % 5) + 1)),
            })
    return {"records": {"pion": records}}


def _render_payload(t_count):
    raw = _raw_transfer(t_count)
    cells = []
    for cell in raw["cells"]:
        cells.append({
            "t_index": cell["t_index"], "delta_index": cell["delta_index"],
            "B": 10.0 + cell["delta_index"], "L": cell["solution"]["transfer"],
            "C": 0.9 + 0.1 * cell["delta_index"], "source": "stable_direct",
            "sigma_C_total": 0.02,
        })
    return {
        "status": "available", "frozen": True, "raw_hgcer_score_fingerprint": "raw",
        "refinement_fingerprint": "part2b", "t_edges": raw["t_edges"],
        "delta_edges": raw["delta_edges"], "cells": cells,
        "t_results": [{"t_index": index, "alpha": 1.0, "f_N": 1.0, "f_D": 1.0, "status": "available"} for index in range(t_count)],
        "histograms": {index: {"host": object(), "baseline_pion": object(), "refined_pion": object(), "baseline_clean": object(), "refined_clean": object()} for index in range(t_count)},
        "global_histograms": {"host": object(), "baseline_pion": object(), "refined_pion": object(), "baseline_clean": object(), "refined_clean": object()},
    }


class PionHGCerRefinementPurePythonTests(unittest.TestCase):
    def test_config_is_permanently_non_authoritative_and_has_fixed_defaults(self):
        config = _config()
        self.assertEqual(config["signed_support_ratio_floor"], 1.0e-4)
        self.assertEqual(config["seed_tie_delta_nll"], 0.5)
        self.assertEqual(config["seed_ambiguity_relative_score_spread"], 0.10)
        for label in refinement.ROOT_SAFE_REFINEMENT_LABELS.values():
            self.assertTrue(label.isascii())
            self.assertNotIn("pre-RF", label)
        with self.assertRaisesRegex(ValueError, "cannot be authoritative"):
            _config(authoritative=True)
        with self.assertRaisesRegex(ValueError, "canonical t"):
            _config(normalization_scope="setting_wide")

    def test_signed_normalization_preserves_parent_and_never_clips(self):
        result = refinement.calculate_part2b_relative_corrections([
            {"t_index": 0, "delta_index": 0, "B": 20.0, "L": 0.5, "source": "stable_direct"},
            {"t_index": 0, "delta_index": 1, "B": 80.0, "L": 2.0, "source": "stable_direct"},
            {"t_index": 0, "delta_index": 2, "B": 9.0, "L": None, "source": "identity_no_hgcer_refinement"},
        ], _config())
        self.assertEqual(result["status"], "available")
        self.assertAlmostEqual(result["alpha"], 100.0 / 170.0)
        self.assertAlmostEqual(sum(cell["B"] * cell["C"] for cell in result["cells"]), 109.0)
        self.assertGreater(result["cells"][1]["C"], 1.0)  # No arbitrary clipping.
        self.assertEqual(result["cells"][2]["C"], 1.0)
        self.assertGreater(result["f_N"], 1.0e-4)
        self.assertGreater(result["f_D"], 1.0e-4)

    def test_signed_cancellation_floor_uses_identity_and_reports_continuous_ratios(self):
        result = refinement.calculate_part2b_relative_corrections([
            {"t_index": 0, "delta_index": 0, "B": 1.0, "L": 1.0},
            {"t_index": 0, "delta_index": 1, "B": -1.0, "L": 2.0},
        ], _config())
        self.assertEqual(result["status"], "identity_fallback")
        self.assertEqual(result["reason"], "signed_normalization_unstable")
        self.assertEqual(result["f_N"], 0.0)
        self.assertEqual([cell["C"] for cell in result["cells"]], [1.0, 1.0])

    def test_part2b_rebuilds_scores_without_ambiguous_direct_seed(self):
        raw, original = _raw_transfer(), None
        raw["cells"][1]["solution"]["transfer"] = 10.0
        original = deepcopy(raw)
        part1 = _part1_for_raw(raw)
        for record in part1["records"]["pion"]:
            record["P_hgcer_npeSum"] = float(record["delta_index"] + 1)

        def audit(values, *_args, **_kwargs):
            # The middle cell is a tied-seed ambiguity.  The stable neighbors
            # are sufficient for a newly-computed same-t bracket.
            return {"status": "ambiguous" if values and int(values[0]) == 2 else "stable", "reason": "ambiguous"}

        with patch.object(refinement, "audit_part2b_response_score", side_effect=audit), patch.object(refinement, "fit_zero_photoelectron_response", return_value={"fit_status": "fit_invalid", "reason": "not_needed"}):
            result = refinement.build_part2b_usable_score_map(raw, part1, _config())
        middle = result["scores"][(0, 1)]
        self.assertEqual(middle["source"], "same_t_delta_bracket")
        self.assertAlmostEqual(middle["L"], 0.3)
        self.assertEqual(raw, original)  # Frozen Part-2 input was never changed.

    def test_control_cache_contract_rejects_wrong_mask_coordinate_or_source_algebra(self):
        raw, contract = _raw_transfer(), fingerprint_hgcer_pid_contract()
        record = {
            "allcuts": True, "source_label": "prompt", "coefficient": 1.0,
            "rf_state": "noRF", "coordinate_fingerprint": "coordinate-fingerprint",
            "P_hgcer_npeSum": 3.0, "t_index": 0, "adj_t": 0.5,
            "ssdelta": 0.5, "delta_index": 0,
        }
        cache = {
            "physical_pion_control_mask": contract["masks"]["physical_pion_control"],
            "physical_pion_control_mask_fingerprint": contract["fingerprint"],
            "coordinate_fingerprint": "coordinate-fingerprint",
            "by_t": ({"records": (record,), "source_accounting": {"prompt": {"allcuts_records": 1, "coefficient": 1.0}}},),
        }
        self.assertEqual(refinement._audit_control_records(raw, cache)["status"], "available")
        cache["by_t"][0]["records"][0]["P_hgcer_npeSum"] = 2.0
        self.assertEqual(refinement._audit_control_records(raw, cache)["status"], "unavailable")
        cache["by_t"][0]["records"][0]["P_hgcer_npeSum"] = 3.0
        cache["by_t"][0]["source_accounting"]["prompt"]["allcuts_records"] = 2
        self.assertEqual(refinement._audit_control_records(raw, cache)["status"], "unavailable")

    def test_phi_invariant_fingerprint_input_and_dynamic_manifest(self):
        raw = _raw_transfer()
        parents = {0: {"pion_parent_id": "parent"}}
        cache = {"physical_pion_control_mask_fingerprint": "mask"}
        cells = [{"t_index": 0, "delta_index": 0, "B": 1.0, "L": 0.2, "C": 1.0}]
        first = refinement.pion_hgcer_refinement_fingerprint_payload(raw, parents, cache, {"phi_edges": [0, 1], "value": 1}, cells)
        second = refinement.pion_hgcer_refinement_fingerprint_payload(raw, parents, cache, {"phi_edges": [0, 9, 18], "value": 1}, cells)
        self.assertEqual(refinement._fingerprint(first), refinement._fingerprint(second))
        for count in (1, 2, 3, 5):
            manifest = refinement.expected_pion_hgcer_refinement_page_manifest(_render_payload(count))
            page_ids = [page["page_id"] for page in manifest]
            for t_index in range(count):
                for stem in ("raw_L", "C", "B", "pion_mm", "mm_products", "clean_mm", "resonance_closure", "signal_closure"):
                    exact = "part2b_{}_t{}".format(stem, t_index + 1)
                    matching = [page for page in manifest if page["page_id"] in {exact, exact + "_unavailable"}]
                    self.assertEqual(len(matching), 1, (count, stem, t_index))
                    self.assertEqual(matching[0]["t_index"], t_index)
            self.assertEqual(len(page_ids), len(set(page_ids)))

    def test_graphical_specs_require_semantic_primitives_and_old_part2_prefix_is_untouched(self):
        for key, spec in refinement.PART2B_PAGE_SPECS.items():
            if spec["kind"] == "graphical":
                self.assertTrue(spec["required_roles"], key)
                self.assertNotEqual(spec["renderer"], "status")
        old = transfer.expected_pion_hgcer_transfer_page_manifest({"status": "unavailable", "reason": "test"})
        new = old + refinement.expected_pion_hgcer_refinement_page_manifest(_render_payload(2))
        self.assertEqual(new[:len(old)], old)
        self.assertTrue(all(page["page_id"].startswith("part2b_") for page in new[len(old):]))


@unittest.skipIf(ROOT is None, "PyROOT is unavailable")
class PionHGCerRefinementPyROOTTests(unittest.TestCase):
    @staticmethod
    def _histogram(name, fills=()):
        histogram = ROOT.TH1D(name, name, 20, 1.0, 1.4)
        histogram.SetDirectory(0); histogram.Sumw2()
        for value, weight in fills:
            histogram.Fill(value, weight)
        return histogram

    def _products(self, *, coefficient=1.0, expected_weight=1.0):
        raw = _raw_transfer(delta_count=1)
        raw["delta_edges"] = [0.0, 1.0]
        expected = self._histogram("expected_part2b", [(1.12, expected_weight)])
        host = self._histogram("host_part2b", [(1.12, 50.0)])
        contract = raw["pid_contract"]
        record = {
            "allcuts": True, "source_label": "prompt", "coefficient": coefficient,
            "rf_state": "noRF", "coordinate_fingerprint": "coordinate-fingerprint",
            "P_hgcer_npeSum": 3.0, "t_index": 0, "adj_t": 0.5,
            "ssdelta": 0.5, "delta_index": 0, "adj_MM": 1.12,
        }
        cache = {
            "physical_pion_control_mask": contract["masks"]["physical_pion_control"],
            "physical_pion_control_mask_fingerprint": contract["fingerprint"],
            "coordinate_fingerprint": "coordinate-fingerprint",
            "by_t": ({"records": (record,), "source_accounting": {"prompt": {"allcuts_records": 1, "coefficient": coefficient}}},),
        }
        parents = ({
            "t_bin_index": 0, "pion_parent_id": "parent-1",
            "final_diagnostic_application_result": {
                # This is intentionally the current event-template target,
                # never the component-fit/SIMC model histogram.
                "H_pion_subtraction_template_MM": expected,
                "H_pion_control_model": object(), "weights": object(),
            },
        },)
        proton_products = ({"cleaned_targets_pre_rf": {"h_mm": host}},)
        return raw, cache, parents, proton_products, expected, host

    def test_event_template_closure_and_detached_products_preserve_authoritative_inputs(self):
        raw, cache, parents, products, expected, host = self._products()
        raw_before = deepcopy(raw)
        score_view = {"raw_map_fingerprint": raw["map_fingerprint"], "audits": {}, "scores": {(0, 0): {"L": 2.0, "source": "stable_direct"}}}
        with patch.object(refinement, "build_part2b_usable_score_map", return_value=score_view), patch.object(refinement, "simc_shape_pion_weight_from_value", return_value=1.0):
            result = refinement.build_pion_hgcer_relative_refinement(raw, {}, cache, parents, products, config=_config())
        self.assertEqual(result["status"], "available")
        self.assertTrue(result["t_results"][0]["baseline_event_template_closure"]["passed"])
        self.assertTrue(result["t_results"][0]["parent_integral_closure"]["passed"])
        detached = result["histograms"][0]
        self.assertAlmostEqual(detached["baseline_pion"].Integral(0, 21), expected.Integral(0, 21))
        self.assertAlmostEqual(detached["refined_pion"].Integral(0, 21), expected.Integral(0, 21))
        self.assertAlmostEqual(detached["refined_clean"].GetBinContent(detached["refined_clean"].FindBin(1.12)), 49.0)
        self.assertEqual(raw, raw_before)
        self.assertAlmostEqual(host.Integral(0, 21), 50.0)  # Authoritative host unchanged.

        # Detached products intentionally retain negative clean bins.
        empty_host = self._histogram("empty_part2b_host")
        with patch.object(refinement, "build_part2b_usable_score_map", return_value=score_view), patch.object(refinement, "simc_shape_pion_weight_from_value", return_value=1.0):
            negative = refinement.build_pion_hgcer_relative_refinement(
                raw, {}, cache, parents, ({"cleaned_targets_pre_rf": {"h_mm": empty_host}},), config=_config(),
            )
        self.assertLess(float(negative["histograms"][0]["refined_clean"].GetMinimum()), 0.0)

    def test_perturbed_baseline_coefficient_fails_closed_and_renderer_has_semantic_primitives(self):
        raw, cache, parents, products, _expected, _host = self._products(coefficient=2.0, expected_weight=1.0)
        score_view = {"raw_map_fingerprint": raw["map_fingerprint"], "audits": {}, "scores": {(0, 0): {"L": 2.0, "source": "stable_direct"}}}
        with patch.object(refinement, "build_part2b_usable_score_map", return_value=score_view), patch.object(refinement, "simc_shape_pion_weight_from_value", return_value=1.0):
            result = refinement.build_pion_hgcer_relative_refinement(raw, {}, cache, parents, products, config=_config())
        self.assertEqual(result["t_results"][0]["status"], "unavailable")
        self.assertEqual(result["t_results"][0]["correction_reason"], "baseline_pion_event_template_closure_failed")

        # Use a valid payload for a real PDF/semantic rendering check.
        raw, cache, parents, products, _expected, _host = self._products()
        with patch.object(refinement, "build_part2b_usable_score_map", return_value=score_view), patch.object(refinement, "simc_shape_pion_weight_from_value", return_value=1.0):
            valid = refinement.build_pion_hgcer_relative_refinement(raw, {}, cache, parents, products, config=_config())
        reference = self._histogram("part2b_closure", [(1.15, 2.0)])
        inputs = {"by_t": {0: {"simc": {"pi-n": reference, "pi-SIDIS": reference, "pi-delta": reference}, "signal": {"K-Lambda": {"hist": reference}, "K-Sigma0": {"hist": None, "status": "unavailable"}}}}}
        with tempfile.TemporaryDirectory() as temporary_directory:
            pdf = Path(temporary_directory) / "part2b.pdf"
            opener = ROOT.TCanvas("part2b_open", "part2b_open", 1, 1)
            opener.Print(str(pdf) + "[")
            manifest = refinement.render_pion_hgcer_refinement_pages(str(pdf), valid, page_manifest=[], renderer_inputs=inputs, close_pdf=True)
            self.assertTrue(pdf.exists() and pdf.stat().st_size > 0)
            # The existing Part-2 page list is the immutable PDF prefix when
            # Part 2B is appended to the same multipage file.
            old_payload = {"status": "unavailable", "reason": "synthetic Part-2 unavailable"}
            old_pdf = Path(temporary_directory) / "part2_old.pdf"
            old_emitted = transfer.render_pion_hgcer_zerope_transfer_pages(
                old_pdf, old_payload, close_pdf=True,
            )
            prefix_pdf = Path(temporary_directory) / "part2_plus_part2b.pdf"
            prefix_emitted = transfer.render_pion_hgcer_zerope_transfer_pages(
                prefix_pdf, old_payload, close_pdf=False,
            )
            new_emitted = prefix_emitted + refinement.render_pion_hgcer_refinement_pages(
                prefix_pdf, valid, renderer_inputs=inputs, close_pdf=True,
            )
            self.assertEqual(new_emitted[:len(old_emitted)], old_emitted)
            self.assertTrue(all(page["page_id"].startswith("part2b_") for page in new_emitted[len(old_emitted):]))
        for page in manifest:
            if page["page_kind"] == "graphical":
                roles = {entry["role"] for entry in page["semantic_primitives"]}
                self.assertTrue(set(page["required_roles"]).issubset(roles), page["page_id"])
        self.assertEqual(refinement.serialize_pion_hgcer_refinement(valid).get("histograms"), None)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
