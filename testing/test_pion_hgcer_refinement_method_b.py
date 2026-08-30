"""Pure-Python Phase-C Method-B closure contracts."""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_refinement_method_b as method_b


MM_EDGES = (0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)


def _config(*, regions=None, protected=None):
    return {
        "mm_regions": regions or [
            {
                "region_name": "pi_n",
                "mm_low": 0.9,
                "mm_high": 1.0,
                "region_role": "pion_sensitive",
                "window_source": "pion_control.windows.pi_n",
                "mm_offset_applied": 0.0,
            },
            {
                "region_name": "pi_sidis",
                "mm_low": 1.0,
                "mm_high": 1.1,
                "region_role": "pion_sensitive",
                "window_source": "pion_control.windows.pi_sidis",
                "mm_offset_applied": 0.0,
            },
        ],
        "protected_regions": protected or [
            {
                "region_name": "KLambda",
                "mm_low": 1.3,
                "mm_high": 1.35,
                "region_role": "protected_signal",
                "window_source": "proton_cleaning.validation_windows.lambda_peak",
                "mm_offset_applied": 0.0,
            },
            {
                "region_name": "KSigma0",
                "mm_low": 1.35,
                "mm_high": 1.4,
                "region_role": "protected_signal",
                "window_source": "kaon_nosub.windows.k_sigma0_signal",
                "mm_offset_applied": 0.0,
            },
        ],
    }


def _phase_a(
    *,
    t_edges=(0.0, 1.0),
    delta_edges=(0.0, 10.0, 20.0),
    host_state="proton_cleaned",
):
    closure = {
        "passed": True,
        "global_full": {"authoritative": {"edges": list(MM_EDGES), "nbins": len(MM_EDGES) - 1}},
    }
    return {
        "schema_version": "pion_hgcer_event_contract/v1",
        "status": "available",
        "available": True,
        "contract_fingerprint": "phase-a-contract",
        "coordinate_fingerprint": "phase-a-coordinates",
        "pion_event_population_fingerprint": "pion-events",
        "kaon_host_population_fingerprint": "host-events",
        "baseline_weight_provenance": {"state": "frozen"},
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
        "rf_restoration_applied": False,
        "canonical_t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "pion_closure": copy.deepcopy(closure),
        "host_closure": copy.deepcopy(closure),
        "pion_records": [],
        "kaon_host_records": [],
    }


def _index(value, edges):
    for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        if low <= value < high:
            return index
    raise AssertionError("fixture coordinate outside canonical geometry")


def _record(
    population,
    contribution,
    *,
    t,
    delta,
    mm,
    t_edges,
    delta_edges,
    host_state="proton_cleaned",
    source_coefficient=1.0,
    proton_factor=1.0,
):
    record = {
        "analysis_abs_t": float(t),
        "SHMS_delta": float(delta),
        "analysis_MM": float(mm),
        "canonical_t_index": _index(t, t_edges),
        "delta_index": _index(delta, delta_edges),
        "nommcuts": True,
        "signed_source_coefficient": float(source_coefficient),
    }
    if population == "pion":
        record.update({
            "baseline_pion_weight_w0": 99.0,
            "signed_baseline_event_contribution": float(contribution),
        })
    else:
        record.update({
            "host_state": host_state,
            "proton_cleaning_factor": float(proton_factor),
            "signed_host_event_contribution": float(contribution),
        })
    return record


def _add_region_events(
    phase,
    *,
    t=0.5,
    delta=5.0,
    mm=0.95,
    pion_weights=(1.0,) * 10,
    host_weights=(1.0,) * 10,
    source_coefficient=1.0,
    host_state=None,
):
    t_edges = phase["canonical_t_edges"]
    delta_edges = phase["delta_edges"]
    state = host_state or phase["host_state"]
    for weight in pion_weights:
        phase["pion_records"].append(_record(
            "pion", weight, t=t, delta=delta, mm=mm,
            t_edges=t_edges, delta_edges=delta_edges,
            source_coefficient=source_coefficient,
        ))
    for weight in host_weights:
        phase["kaon_host_records"].append(_record(
            "host", weight, t=t, delta=delta, mm=mm,
            t_edges=t_edges, delta_edges=delta_edges, host_state=state,
            source_coefficient=source_coefficient,
        ))


def _supported_phase(*, host_state="proton_cleaned", t_edges=(0.0, 1.0), delta_edges=(0.0, 10.0, 20.0)):
    phase = _phase_a(t_edges=t_edges, delta_edges=delta_edges, host_state=host_state)
    for delta in ((delta_edges[0] + delta_edges[1]) / 2.0, (delta_edges[1] + delta_edges[2]) / 2.0):
        for mm in (0.95, 1.05):
            _add_region_events(phase, delta=delta, mm=mm)
    return phase


def _cell(result, t_index=0, delta_index=0):
    return next(
        entry for entry in result["cells"]
        if entry["t_index"] == t_index and entry["delta_index"] == delta_index
    )


class PionHGCerMethodBTests(unittest.TestCase):
    def test_phase_a_event_contributions_are_used_exactly_once(self):
        phase = _phase_a()
        for delta in (5.0, 15.0):
            _add_region_events(
                phase, delta=delta, mm=0.95,
                pion_weights=(2.0,) * 10, host_weights=(3.0,) * 10,
                source_coefficient=-7.0,
            )
        result = method_b.build_pion_hgcer_method_b(phase, config=_config(regions=[_config()["mm_regions"][0]]))
        row = _cell(result)["regions"][0]
        self.assertAlmostEqual(row["baseline_pion_yield"], 20.0)
        self.assertAlmostEqual(row["host_yield"], 30.0)
        self.assertAlmostEqual(row["raw_ratio"], 1.5)

    def test_identity_host_and_unity_cleaned_host_are_numerically_equivalent(self):
        cleaned = _supported_phase(host_state="proton_cleaned")
        identity = copy.deepcopy(cleaned)
        identity["host_state"] = "identity_no_proton_cleaning"
        for record in identity["kaon_host_records"]:
            record["host_state"] = "identity_no_proton_cleaning"
        first = method_b.build_pion_hgcer_method_b(cleaned, config=_config())
        second = method_b.build_pion_hgcer_method_b(identity, config=_config())
        self.assertEqual(first["parent_region_references"], second["parent_region_references"])
        for left, right in zip(first["cells"], second["cells"]):
            self.assertEqual(left["regions"], right["regions"])
            self.assertEqual(left["shape_chi2"], right["shape_chi2"])
            self.assertEqual(left["candidate_L_B_status"], right["candidate_L_B_status"])

    def test_signed_support_uses_absolute_neff_and_rejects_bad_denominators(self):
        scenarios = {
            "cancellation": ((2.0,) * 10 + (-1.0,) * 10, "baseline_cancellation_dominated"),
            "negative": ((-1.0,) * 10, "baseline_signed_yield_nonpositive"),
            "zero": ((1.0,) * 10 + (-1.0,) * 10, "baseline_signed_yield_nonpositive"),
        }
        for name, (weights, reason) in scenarios.items():
            with self.subTest(name=name):
                phase = _phase_a()
                for delta in (5.0, 15.0):
                    _add_region_events(phase, delta=delta, mm=0.95, pion_weights=weights, host_weights=(1.0,) * 20)
                result = method_b.build_pion_hgcer_method_b(phase, config=_config(regions=[_config()["mm_regions"][0]]))
                row = _cell(result)["regions"][0]
                self.assertEqual(row["support_status"], "unavailable")
                self.assertEqual(row["support_reason"], reason)
                self.assertIsNone(row["raw_ratio"])
        metric = method_b._finish_metric({"record_count": 2, "signed_yield": 1.0, "absolute_weight_support": 3.0, "sumw2": 5.0})
        self.assertAlmostEqual(metric["effective_entries"], 9.0 / 5.0)

    def test_raw_ratio_residual_and_uncertainty_are_exact(self):
        phase = _phase_a()
        for delta in (5.0, 15.0):
            _add_region_events(phase, delta=delta, mm=0.95, pion_weights=(1.0,) * 10, host_weights=(2.0,) * 10)
        result = method_b.build_pion_hgcer_method_b(phase, config=_config(regions=[_config()["mm_regions"][0]]))
        row = _cell(result)["regions"][0]
        self.assertAlmostEqual(row["raw_ratio"], 2.0)
        self.assertAlmostEqual(row["raw_ratio_sigma"], math.sqrt(0.8))
        self.assertAlmostEqual(row["residual"], 10.0)
        self.assertAlmostEqual(row["residual_sigma"], math.sqrt(50.0))

    def test_protected_signal_overlap_excludes_lambda_and_sigma0_regions(self):
        protected = _config()["protected_regions"]
        regions = [
            {**_config()["mm_regions"][0], "region_name": "overlap_lambda", "mm_low": 1.30, "mm_high": 1.34},
            {**_config()["mm_regions"][1], "region_name": "overlap_sigma0", "mm_low": 1.35, "mm_high": 1.39},
        ]
        phase = _supported_phase()
        result = method_b.build_pion_hgcer_method_b(phase, config=_config(regions=regions, protected=protected))
        for region in result["mm_regions"]:
            self.assertFalse(region["available"])
            self.assertTrue(region["protected_signal_overlap"])
        self.assertEqual(result["parent_region_references"], [
            reference for reference in result["parent_region_references"]
            if reference["parent_reference_status"] == "unavailable"
        ])

    def test_parent_references_are_same_t_only_and_require_two_delta_cells(self):
        phase = _phase_a(t_edges=(0.0, 1.0, 2.0))
        for t, ratio in ((0.5, 2.0), (1.5, 3.0)):
            for delta in (5.0, 15.0):
                for mm in (0.95, 1.05):
                    _add_region_events(phase, t=t, delta=delta, mm=mm, host_weights=(ratio,) * 10)
        result = method_b.build_pion_hgcer_method_b(phase, config=_config())
        references = {(entry["t_index"], entry["region_name"]): entry for entry in result["parent_region_references"]}
        self.assertAlmostEqual(references[(0, "pi_n")]["parent_reference_ratio"], 2.0)
        self.assertAlmostEqual(references[(1, "pi_n")]["parent_reference_ratio"], 3.0)
        self.assertAlmostEqual(_cell(result, 1, 0)["regions"][0]["parent_relative_ratio"], 1.0)
        one_cell = _phase_a()
        _add_region_events(one_cell, delta=5.0, mm=0.95)
        unavailable = method_b.build_pion_hgcer_method_b(one_cell, config=_config(regions=[_config()["mm_regions"][0]]))
        reference = unavailable["parent_region_references"][0]
        self.assertEqual(reference["parent_reference_status"], "unavailable")
        self.assertIsNone(reference["parent_reference_ratio"])

    def test_zero_and_one_supported_region_statuses_are_distinct(self):
        empty = _phase_a()
        for delta in (5.0, 15.0):
            _add_region_events(empty, delta=delta, mm=1.2)
        zero = method_b.build_pion_hgcer_method_b(empty, config=_config(regions=[_config()["mm_regions"][0]]))
        zero_cell = _cell(zero)
        self.assertEqual(zero_cell["region_consistency_reason"], "no_supported_parent_relative_regions")
        self.assertEqual(zero_cell["candidate_L_B_status"], "unavailable")
        self.assertEqual(zero_cell["method_B_status"], "unavailable")
        single = _supported_phase()
        one = method_b.build_pion_hgcer_method_b(single, config=_config(regions=[_config()["mm_regions"][0]]))
        one_cell = _cell(one)
        self.assertEqual(one_cell["region_consistency_status"], "insufficient_regions")
        self.assertEqual(one_cell["candidate_L_B_status"], "single_region_only")
        self.assertEqual(one_cell["method_B_status"], "unavailable")

    def test_regional_consistency_and_log_space_candidate_rules(self):
        def row(value, sigma):
            return {"parent_relative_status": "available", "parent_relative_ratio": value, "parent_relative_sigma": sigma}
        consistent, _, rows = method_b._region_consistency([row(1.0, 0.1), row(1.05, 0.1)], _config()["region_consistency"] if "region_consistency" in _config() else method_b.DEFAULT_METHOD_B_CONFIG["region_consistency"])
        self.assertEqual(consistent, "region_consistent")
        self.assertEqual(method_b._region_consistency([row(1.0, 0.1), row(1.25, 0.1)], method_b.DEFAULT_METHOD_B_CONFIG["region_consistency"])[0], "region_marginal")
        self.assertEqual(method_b._region_consistency([row(1.0, 0.1), row(1.6, 0.1)], method_b.DEFAULT_METHOD_B_CONFIG["region_consistency"])[0], "region_inconsistent")
        candidate, uncertainty, status = method_b._candidate(rows, consistent, "good")
        weights = [1.0 / (0.1 / 1.0) ** 2, 1.0 / (0.1 / 1.05) ** 2]
        expected = math.exp((math.log(1.0) * weights[0] + math.log(1.05) * weights[1]) / sum(weights))
        self.assertEqual(status, "available_multi_region")
        self.assertAlmostEqual(candidate, expected)
        self.assertGreater(uncertainty, 0.0)

    def test_poor_shape_veto_is_explicit_and_shape_closure_is_unscaled(self):
        phase = _supported_phase()
        poor_shape = {
            "shape_status": "poor", "shape_reason": "synthetic_shape_failure",
            "shape_chi2": 99.0, "shape_ndf": 2, "shape_chi2_ndf": 49.5,
            "shape_max_abs_pull": 9.0, "shape_usable_bin_count": 2,
            "mm_edges": list(MM_EDGES), "bins": [], "underflow": {}, "overflow": {},
        }
        with mock.patch.object(method_b, "_shape_payload", return_value=poor_shape):
            result = method_b.build_pion_hgcer_method_b(phase, config=_config())
        cell = _cell(result)
        self.assertIsNone(cell["candidate_L_B"])
        self.assertEqual(cell["candidate_L_B_status"], "shape_poor_veto")
        self.assertEqual(cell["method_B_status"], "shape_inconsistent")
        regions = method_b._annotate_regions(
            _config()["mm_regions"], _config()["protected_regions"]
        )
        source = {
            "host_events": [(0.95, 1.0)] * 10 + [(1.05, 1.0)] * 10,
            "pion_events": [(0.95, 1.0)] * 10 + [(1.05, 1.0)] * 10,
        }
        shape = method_b._shape_payload(source, MM_EDGES, regions, _config()["protected_regions"], method_b.DEFAULT_METHOD_B_CONFIG["shape"], 1.0e-12)
        self.assertEqual(shape["shape_status"], "good")
        self.assertEqual(shape["shape_chi2"], 0.0)
        scaled = method_b._shape_payload({"host_events": [(0.95, 2.0)] * 10 + [(1.05, 2.0)] * 10, "pion_events": [(0.95, 1.0)] * 10 + [(1.05, 1.0)] * 10}, MM_EDGES, regions, _config()["protected_regions"], method_b.DEFAULT_METHOD_B_CONFIG["shape"], 1.0e-12)
        self.assertNotEqual(scaled["bins"][2]["residual"], 0.0)

    def test_shape_mismatch_and_protected_bins_are_excluded(self):
        regions = [{**_config()["mm_regions"][0], "mm_low": 0.9, "mm_high": 1.4}]
        protected = _config()["protected_regions"]
        annotated_regions = [{**regions[0], "available": True}]
        mismatch = method_b._shape_payload(
            {"host_events": [(1.05, 1.0)] * 100, "pion_events": [(0.95, 1.0)] * 100},
            MM_EDGES, annotated_regions, protected, method_b.DEFAULT_METHOD_B_CONFIG["shape"], 1.0e-12,
        )
        self.assertEqual(mismatch["shape_status"], "poor")
        protected_only = method_b._shape_payload(
            {"host_events": [(1.32, 100.0), (1.37, 100.0)], "pion_events": [(1.32, 1.0), (1.37, 1.0)]},
            MM_EDGES, annotated_regions, protected, method_b.DEFAULT_METHOD_B_CONFIG["shape"], 1.0e-12,
        )
        self.assertEqual(protected_only["shape_usable_bin_count"], 0)
        self.assertEqual(protected_only["shape_chi2"], None)

    def test_method_a_and_phi_metadata_do_not_affect_method_b(self):
        phase = _supported_phase()
        phase["method_a"] = {"candidate": 0.1, "support": "unsupported", "hgcer_x": 4.0}
        first = method_b.build_pion_hgcer_method_b(phase, config=_config())
        changed = copy.deepcopy(phase)
        changed["method_a"] = {"candidate": 999.0, "support": "supported", "hgcer_x": -999.0}
        changed["phi_setting"] = "Right"
        changed["phi_edges"] = [0.0, 180.0]
        second = method_b.build_pion_hgcer_method_b(changed, config=_config())
        self.assertEqual(first["cells"], second["cells"])
        self.assertEqual(first["parent_region_references"], second["parent_region_references"])
        self.assertEqual(first["fingerprint"], second["fingerprint"])

    def test_missing_delta_is_not_interpolated_and_inputs_are_immutable(self):
        phase = _phase_a(delta_edges=(0.0, 10.0, 20.0, 30.0))
        for delta in (5.0, 25.0):
            for mm in (0.95, 1.05):
                _add_region_events(phase, delta=delta, mm=mm)
        before = copy.deepcopy(phase)
        result = method_b.build_pion_hgcer_method_b(phase, config=_config())
        self.assertEqual(phase, before)
        missing = _cell(result, 0, 1)
        self.assertEqual(missing["regions"][0]["support_status"], "unavailable")
        self.assertEqual(missing["regions"][0]["parent_relative_status"], "unavailable")

    def test_unavailable_result_keeps_available_phase_a_provenance_and_forbids_application(self):
        phase = _phase_a()
        result = method_b.build_pion_hgcer_method_b(phase, config=_config())
        self.assertFalse(result["available"])
        self.assertEqual(result["phase_a_contract_fingerprint"], "phase-a-contract")
        self.assertEqual(result["coordinate_fingerprint"], "phase-a-coordinates")
        self.assertEqual(result["host_state"], "proton_cleaned")
        self.assertEqual(result["t_edges"], [0.0, 1.0])
        self.assertEqual(result["delta_edges"], [0.0, 10.0, 20.0])
        self.assertEqual(result["mm_regions"], _config()["mm_regions"])
        self.assertEqual(result["protected_regions"], _config()["protected_regions"])
        json.dumps(result, allow_nan=False)
        phase_unavailable = _supported_phase()
        phase_unavailable["available"] = False
        no_phase_provenance = method_b.build_pion_hgcer_method_b(
            phase_unavailable, config=_config()
        )
        self.assertNotIn("phase_a_contract_fingerprint", no_phase_provenance)
        self.assertNotIn("host_state", no_phase_provenance)
        forbidden = {"C_B", "C_final", "refined_pion_weight", "applied_refinement_weight"}
        def walk(value):
            if isinstance(value, dict):
                self.assertTrue(forbidden.isdisjoint(value))
                for child in value.values():
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)
        walk(result)
        source = (REPO_ROOT / "src/cuts/pion_hgcer_refinement_method_b.py").read_text(encoding="utf-8")
        for token in ("TCanvas", ".Print(", "apply_refined"):
            self.assertNotIn(token, source)

    def test_runtime_handoff_keeps_method_b_detached_and_after_phase_a(self):
        source = (REPO_ROOT / "src/cuts/rand_sub.py").read_text(encoding="utf-8")
        event_position = source.index("build_pion_hgcer_event_contract(")
        method_a_position = source.index("build_pion_hgcer_method_a(")
        method_b_position = source.index("build_pion_hgcer_method_b(")
        self.assertLess(event_position, method_a_position)
        self.assertLess(event_position, method_b_position)
        self.assertIn('histDict["_pion_hgcer_method_b"]', source)
        self.assertIn("resolve_pion_hgcer_method_b_config(", source)


if __name__ == "__main__":
    unittest.main()
