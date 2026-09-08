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
                "region_name": "pion_sensitive_low",
                "mm_low": 0.8,
                "mm_high": 1.10,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
            },
            {
                "region_name": "pion_sensitive_high",
                "mm_low": 1.23,
                "mm_high": 1.5,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
            },
        ],
        "protected_regions": protected or [
            {
                "region_name": "KLambdaSigma0",
                "mm_low": 1.10,
                "mm_high": 1.23,
                "region_role": "protected_signal",
                "window_source": "method_b.fixed_lambda_sigma_protected",
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
        for mm in (0.95, 1.35):
            _add_region_events(phase, delta=delta, mm=mm)
    return phase


def _adaptive_supported_phase(*, t_edges=(0.0, 1.0), delta_edges=(0.0, 10.0, 20.0)):
    """Populate every frozen Phase-A atomic MM interval with two delta parents."""
    phase = _phase_a(t_edges=t_edges, delta_edges=delta_edges)
    for t_low, t_high in zip(t_edges, t_edges[1:]):
        t = 0.5 * (t_low + t_high)
        for delta_low, delta_high in zip(delta_edges, delta_edges[1:]):
            delta = 0.5 * (delta_low + delta_high)
            for mm in (0.85, 0.95, 1.05, 1.25, 1.35, 1.45):
                _add_region_events(phase, t=t, delta=delta, mm=mm)
    return phase


def _cell(result, t_index=0, delta_index=0):
    return next(
        entry for entry in result["cells"]
        if entry["t_index"] == t_index and entry["delta_index"] == delta_index
    )


def _adaptive_partition_cells(delta_count, events_by_delta):
    """Build the baseline-only source cells consumed by the adaptive partitioner."""
    return {
        (0, delta_index): {
            "t_index": 0,
            "delta_index": delta_index,
            "pion_events": list(events_by_delta.get(delta_index, ())),
        }
        for delta_index in range(delta_count)
    }


class PionHGCerMethodBTests(unittest.TestCase):
    def test_runtime_resolver_uses_only_the_fixed_phase_a_complement_partition(self):
        phase = _phase_a()
        expected_regions = [
            {
                "region_name": "pion_sensitive_low",
                "mm_low": 0.8,
                "mm_high": 1.10,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
            },
            {
                "region_name": "pion_sensitive_high",
                "mm_low": 1.23,
                "mm_high": 1.5,
                "region_role": "pion_sensitive",
                "window_source": "phase_a.mm_binning_complement_of_lambda_sigma_protected",
                "mm_offset_applied": 0.0,
            },
        ]
        expected_protected = [{
            "region_name": "KLambdaSigma0",
            "mm_low": 1.10,
            "mm_high": 1.23,
            "region_role": "protected_signal",
            "window_source": "method_b.fixed_lambda_sigma_protected",
            "mm_offset_applied": 0.0,
        }]
        resolved = []
        for phi_setting, mm_offset_data in (("Left", 0.0), ("Center", 0.015), ("Right", -0.020)):
            config = method_b.resolve_pion_hgcer_method_b_config(
                phase_a_contract=phase,
                inp_dict={"pion_component_windows": {"must_not": "matter"}},
                phi_setting=phi_setting,
                mm_offset_data=mm_offset_data,
            )
            self.assertEqual(config["mm_regions"], expected_regions)
            self.assertEqual(config["protected_regions"], expected_protected)
            resolved.append(config)
        self.assertEqual(resolved[0], resolved[1])
        self.assertEqual(resolved[1], resolved[2])

        low, high = resolved[0]["mm_regions"]
        protected = resolved[0]["protected_regions"][0]
        self.assertTrue(method_b._in_region(0.8, low))
        self.assertFalse(method_b._in_region(1.10, low))
        self.assertTrue(method_b._in_region(1.10, protected))
        self.assertFalse(method_b._in_region(1.23, protected))
        self.assertTrue(method_b._in_region(1.23, high))
        self.assertFalse(method_b._in_region(1.5, high))

        outputs = []
        for phi_setting, mm_offset_data in (("Left", 0.0), ("Center", 0.015), ("Right", -0.020)):
            config = method_b.resolve_pion_hgcer_method_b_config(
                phase_a_contract=_supported_phase(),
                inp_dict={"pion_component_windows": {"must_not": "matter"}},
                phi_setting=phi_setting,
                mm_offset_data=mm_offset_data,
            )
            result = method_b.build_pion_hgcer_method_b(_supported_phase(), config=config)
            self.assertTrue(result["available"])
            outputs.append(result)
        self.assertEqual(outputs[0]["cells"], outputs[1]["cells"])
        self.assertEqual(outputs[1]["cells"], outputs[2]["cells"])
        self.assertEqual(outputs[0]["parent_region_references"], outputs[2]["parent_region_references"])
        self.assertEqual(outputs[0]["fingerprint"], outputs[2]["fingerprint"])
        first_cell = _cell(outputs[0])
        self.assertEqual(
            [row["region_name"] for row in first_cell["regions"]],
            ["pion_sensitive_low", "pion_sensitive_high"],
        )
        self.assertEqual(
            {
                (reference["t_index"], reference["region_name"])
                for reference in outputs[0]["parent_region_references"]
            },
            {(0, "pion_sensitive_low"), (0, "pion_sensitive_high")},
        )
        usable = [
            row for row in first_cell["regions"]
            if row["parent_relative_status"] == "available"
        ]
        expected_candidate, expected_sigma, expected_status = method_b._candidate(
            usable,
            first_cell["region_consistency_status"],
            first_cell["shape_status"],
        )
        self.assertEqual(first_cell["candidate_L_B_status"], expected_status)
        self.assertAlmostEqual(first_cell["candidate_L_B"], expected_candidate)
        self.assertAlmostEqual(first_cell["candidate_L_B_uncertainty"], expected_sigma)

    def test_runtime_resolver_rejects_phase_a_domains_that_do_not_straddle_window(self):
        phase = _phase_a()
        narrow_edges = [0.8, 0.9, 1.0, 1.10, 1.23]
        for closure_name in ("pion_closure", "host_closure"):
            phase[closure_name]["global_full"]["authoritative"] = {
                "edges": narrow_edges,
                "nbins": len(narrow_edges) - 1,
            }
        with self.assertRaises(method_b.MethodBUnavailable) as raised:
            method_b.resolve_pion_hgcer_method_b_config(phase_a_contract=phase)
        self.assertEqual(
            raised.exception.reason,
            "method_b_protected_region_outside_phase_a_domain",
        )

    def test_shape_bins_must_be_fully_inside_the_fixed_sensitive_complement(self):
        config = method_b.resolve_pion_hgcer_method_b_config(
            phase_a_contract=_phase_a()
        )
        regions = method_b._annotate_regions(
            config["mm_regions"], config["protected_regions"]
        )
        protected = config["protected_regions"]
        self.assertTrue(method_b._bin_allowed(0.8, 1.10, regions, protected))
        self.assertFalse(method_b._bin_allowed(1.05, 1.15, regions, protected))
        self.assertFalse(method_b._bin_allowed(1.20, 1.30, regions, protected))
        self.assertTrue(method_b._bin_allowed(1.23, 1.5, regions, protected))

    def test_runtime_resolver_has_no_shared_subtraction_window_dependency(self):
        source = (REPO_ROOT / "src/cuts/pion_hgcer_refinement_method_b.py").read_text(
            encoding="utf-8"
        )
        resolver = source[
            source.index("def resolve_pion_hgcer_method_b_config("):
            source.index("\ndef _normalize_regions", source.index("def resolve_pion_hgcer_method_b_config("))
        ]
        for forbidden in (
            "get_particle_subtraction_component_fit_window_config",
            "get_proton_contamination_cleaning_config",
            "resolve_particle_subtraction_component_fit_windows",
        ):
            self.assertNotIn(forbidden, resolver)

    def test_method_b_partition_identifiers_do_not_reach_production_pion_weights(self):
        production = (
            REPO_ROOT / "src/cuts/pion_component_subtraction.py"
        ).read_text(encoding="utf-8")
        production_weight_application = production[
            production.index("def fill_simc_shape_pion_subtraction_templates("):
        ]
        for method_b_only_identifier in (
            "METHOD_B_PROTECTED_MM_LOW",
            "METHOD_B_PROTECTED_MM_HIGH",
            "METHOD_B_PROTECTED_REGION_NAME",
            "KLambdaSigma0",
        ):
            with self.subTest(identifier=method_b_only_identifier):
                self.assertNotIn(method_b_only_identifier, production_weight_application)

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

    def test_protected_signal_overlap_excludes_fixed_lambda_sigma_window(self):
        protected = _config()["protected_regions"]
        regions = [
            {**_config()["mm_regions"][0], "region_name": "overlap_low", "mm_low": 1.05, "mm_high": 1.15},
            {**_config()["mm_regions"][1], "region_name": "overlap_high", "mm_low": 1.20, "mm_high": 1.25},
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
        self.assertAlmostEqual(references[(0, "pion_sensitive_low")]["parent_reference_ratio"], 2.0)
        self.assertAlmostEqual(references[(1, "pion_sensitive_low")]["parent_reference_ratio"], 3.0)
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
            {"host_events": [(1.15, 100.0), (1.20, 100.0)], "pion_events": [(1.15, 1.0), (1.20, 1.0)]},
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

    def test_adaptive_slice_payload_is_additive_and_legacy_serialization_is_exactly_stable(self):
        phase = _adaptive_supported_phase()
        phase_before = copy.deepcopy(phase)
        with mock.patch.object(
            method_b,
            "_build_adaptive_slice_diagnostic",
            return_value={"synthetic": "adaptive-only"},
        ):
            baseline = method_b.build_pion_hgcer_method_b(phase, config=_config())
        result = method_b.build_pion_hgcer_method_b(phase, config=_config())

        baseline.pop("adaptive_slice_diagnostic")
        adaptive = result.pop("adaptive_slice_diagnostic")
        self.assertEqual(
            json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False),
            json.dumps(baseline, sort_keys=True, separators=(",", ":"), allow_nan=False),
        )
        self.assertEqual(phase, phase_before)
        self.assertTrue(adaptive["available"])
        self.assertEqual(
            adaptive["schema_version"],
            method_b.METHOD_B_ADAPTIVE_SLICE_SCHEMA_VERSION,
        )
        self.assertTrue(adaptive["non_authoritative"])
        self.assertFalse(adaptive["candidate_replaces_legacy_method_b"])
        self.assertEqual(adaptive["derived_parent_baseline_neff_reference"], 20.0)
        self.assertEqual(
            adaptive["partition_minimum_baseline_usable_delta_cells"], 2
        )
        self.assertEqual(adaptive["partition_cell_minimum_baseline_neff"], 10.0)
        self.assertEqual(len(adaptive["t_partitions"][0]["low_slices"]), 3)
        self.assertEqual(len(adaptive["t_partitions"][0]["high_slices"]), 3)
        repeat = method_b.build_pion_hgcer_method_b(phase, config=_config())
        self.assertEqual(
            adaptive["fingerprint"],
            repeat["adaptive_slice_diagnostic"]["fingerprint"],
        )
        self.assertTrue(all(
            cell["adaptive_candidate_status"] == "available_multi_slice"
            for cell in adaptive["cells"]
        ))

    def test_adaptive_partition_uses_baseline_only_and_is_shared_by_delta_cells(self):
        phase = _adaptive_supported_phase()
        baseline = method_b.build_pion_hgcer_method_b(phase, config=_config())
        changed_host = copy.deepcopy(phase)
        for record in changed_host["kaon_host_records"]:
            factor = 37.0 if record["analysis_MM"] < 1.1 else 0.25
            record["signed_host_event_contribution"] *= factor
        changed = method_b.build_pion_hgcer_method_b(changed_host, config=_config())
        base_adaptive = baseline["adaptive_slice_diagnostic"]
        changed_adaptive = changed["adaptive_slice_diagnostic"]
        self.assertEqual(
            base_adaptive["t_partitions"], changed_adaptive["t_partitions"]
        )
        expected_ids = [
            row["slice_id"] for row in base_adaptive["t_partitions"][0]["slices"]
        ]
        self.assertTrue(all(
            [row["slice_id"] for row in cell["slices"]] == expected_ids
            for cell in base_adaptive["cells"]
        ))
        partition = base_adaptive["t_partitions"][0]
        low_slices = sorted(partition["low_slices"], key=lambda row: row["mm_low"])
        high_slices = sorted(partition["high_slices"], key=lambda row: row["mm_low"])
        self.assertEqual(low_slices[0]["mm_low"], MM_EDGES[0])
        self.assertEqual(low_slices[-1]["mm_high"], method_b.METHOD_B_PROTECTED_MM_LOW)
        self.assertEqual(high_slices[0]["mm_low"], method_b.METHOD_B_PROTECTED_MM_HIGH)
        self.assertEqual(high_slices[-1]["mm_high"], MM_EDGES[-1])

    def test_adaptive_partition_requires_two_independently_baseline_supported_delta_cells(self):
        """Aggregate parent Neff cannot close a slice with no locally usable deltas."""
        resolved = method_b._resolved_config(_config())
        t_edges = (0.0, 1.0)
        events = {
            delta_index: [(1.05, 1.0), (1.05, 1.0)]
            for delta_index in range(10)
        }
        partitions, reference = method_b._adaptive_t_partitions(
            _adaptive_partition_cells(10, events), t_edges, MM_EDGES,
            resolved["support"], resolved["parent_reference"],
        )
        low = partitions[0]["low_slices"]
        self.assertEqual(reference, 20.0)
        self.assertEqual(len(low), 1)
        self.assertEqual(low[0]["parent_baseline_neff"], 20.0)
        self.assertEqual(low[0]["baseline_supported_delta_cell_count"], 0)
        self.assertEqual(low[0]["partition_support_status"], "unavailable")
        self.assertEqual(
            low[0]["partition_support_reason"],
            "insufficient_baseline_supported_delta_cells",
        )
        self.assertTrue(all(
            row["reason"] == "baseline_effective_entries_below_minimum"
            for row in low[0]["baseline_delta_support_reasons"]
        ))

        repaired = copy.deepcopy(events)
        repaired[0].extend([(1.05, 1.0)] * 8)
        repaired[1].extend([(0.95, 1.0)] * 8)
        repaired[2] = [(0.85, 1.0)] * 10
        repaired[3] = [(0.85, 1.0)] * 10
        partitions, _ = method_b._adaptive_t_partitions(
            _adaptive_partition_cells(10, repaired), t_edges, MM_EDGES,
            resolved["support"], resolved["parent_reference"],
        )
        low = sorted(partitions[0]["low_slices"], key=lambda row: row["mm_low"])
        self.assertEqual([(row["mm_low"], row["mm_high"]) for row in low], [
            (0.8, 0.9), (0.9, 1.1),
        ])
        inner = low[1]
        self.assertEqual(inner["baseline_supported_delta_indices"], [0, 1])
        self.assertEqual(inner["baseline_supported_delta_cell_count"], 2)
        self.assertEqual(inner["partition_support_status"], "available")

    def test_adaptive_partition_support_is_baseline_only_and_uses_legacy_reason_vocabulary(self):
        resolved = method_b._resolved_config(_config())
        support = resolved["support"]
        self.assertEqual(
            method_b._adaptive_baseline_support_reason(
                method_b._adaptive_metric([(1.0, 1.0)] * 9, 0.8, 1.1), support
            ),
            "baseline_effective_entries_below_minimum",
        )
        self.assertEqual(
            method_b._adaptive_baseline_support_reason(
                method_b._adaptive_metric([(1.0, -1.0)] * 10, 0.8, 1.1), support
            ),
            "baseline_signed_yield_nonpositive",
        )
        self.assertEqual(
            method_b._adaptive_baseline_support_reason(
                method_b._adaptive_metric(
                    [(1.0, 1.0)] * 10 + [(1.0, -1.0)] * 9, 0.8, 1.1
                ), support
            ),
            "baseline_cancellation_dominated",
        )

        phase = _adaptive_supported_phase()
        baseline = method_b.build_pion_hgcer_method_b(phase, config=_config())
        changed_host = copy.deepcopy(phase)
        for index, record in enumerate(changed_host["kaon_host_records"]):
            record["signed_host_event_contribution"] = (
                -100.0 if index % 2 else 0.001
            )
        changed = method_b.build_pion_hgcer_method_b(changed_host, config=_config())
        base_partitions = baseline["adaptive_slice_diagnostic"]["t_partitions"]
        changed_partitions = changed["adaptive_slice_diagnostic"]["t_partitions"]
        self.assertEqual(base_partitions, changed_partitions)
        self.assertEqual(
            baseline["adaptive_slice_diagnostic"]["fingerprint_inputs"]["t_partitions"],
            changed["adaptive_slice_diagnostic"]["fingerprint_inputs"]["t_partitions"],
        )

    def test_adaptive_partition_refines_as_baseline_delta_support_allows_and_keeps_sides_asymmetric(self):
        high = method_b.build_pion_hgcer_method_b(
            _adaptive_supported_phase(), config=_config()
        )["adaptive_slice_diagnostic"]["t_partitions"][0]
        low_phase = _phase_a()
        for delta in (5.0, 15.0):
            for mm in (0.85, 0.95, 1.05):
                _add_region_events(low_phase, delta=delta, mm=mm)
            for mm in (1.25, 1.35, 1.45):
                _add_region_events(
                    low_phase, delta=delta, mm=mm,
                    pion_weights=(1.0,) * 5, host_weights=(1.0,) * 10,
                )
        low = method_b.build_pion_hgcer_method_b(
            low_phase, config=_config()
        )["adaptive_slice_diagnostic"]["t_partitions"][0]
        self.assertGreaterEqual(
            len(high["low_slices"]), len(low["low_slices"])
        )
        self.assertGreater(len(low["low_slices"]), len(low["high_slices"]))
        self.assertEqual(len(low["high_slices"]), 1)

    def test_adaptive_outer_remainder_recomputes_signed_delta_support_recursively(self):
        resolved = method_b._resolved_config(_config())
        events = {index: [] for index in range(4)}
        events[0].extend([(1.05, 1.0)] * 10)
        events[1].extend([(1.05, 1.0)] * 10)
        events[2].extend([(0.95, 1.0)] * 10)
        events[3].extend([(0.95, 1.0)] * 10)
        events[2].extend([(0.85, -1.0)] * 9)
        events[3].extend([(0.85, -1.0)] * 9)
        partitions, _ = method_b._adaptive_t_partitions(
            _adaptive_partition_cells(4, events), (0.0, 1.0), MM_EDGES,
            resolved["support"], resolved["parent_reference"],
        )
        low = partitions[0]["low_slices"]
        self.assertEqual(len(low), 1)
        self.assertEqual((low[0]["mm_low"], low[0]["mm_high"]), (0.8, 1.1))
        self.assertEqual(low[0]["baseline_supported_delta_indices"], [0, 1])
        self.assertEqual(low[0]["partition_support_status"], "available")

    def test_adaptive_partition_repair_does_not_change_local_support_evaluation(self):
        resolved = method_b._resolved_config(_config())
        definition = {
            "slice_id": "t0_low_slice0", "side": "low", "slice_index": 0,
            "mm_low": 0.8, "mm_high": 1.1,
        }
        baseline_failure = method_b._adaptive_slice_row(
            definition, [(0.95, 1.0)] * 10, [(0.95, 1.0)] * 9,
            resolved["support"],
        )
        host_failure = method_b._adaptive_slice_row(
            definition, [(0.95, 1.0)] * 9, [(0.95, 1.0)] * 10,
            resolved["support"],
        )
        usable = method_b._adaptive_slice_row(
            definition, [(0.95, 1.0)] * 10, [(0.95, 1.0)] * 10,
            resolved["support"],
        )
        self.assertEqual(
            baseline_failure["support_reason"],
            "baseline_effective_entries_below_minimum",
        )
        self.assertEqual(host_failure["support_reason"], "host_effective_entries_below_minimum")
        self.assertEqual(usable["support_status"], "usable")

    def test_adaptive_parent_references_are_numerically_equivalent_to_legacy_equations(self):
        resolved = method_b._resolved_config(_config())
        legacy_cells = []
        adaptive_cells = {}
        for delta_index, host_weight in enumerate((1.0, 1.5, 2.0)):
            metrics = {
                "host": method_b._adaptive_metric(
                    [(0.95, host_weight)] * 12, 0.8, 1.1
                ),
                "baseline": method_b._adaptive_metric(
                    [(0.95, 1.0)] * 12, 0.8, 1.1
                ),
            }
            legacy_row = method_b._region_row(
                {"region_name": "equivalent", "available": True, "reason": None},
                metrics, resolved["support"],
            )
            adaptive_row = {
                **copy.deepcopy(legacy_row),
                "slice_id": "t0_low_slice0",
                "side": "low",
                "slice_index": 0,
                "mm_low": 0.8,
                "mm_high": 1.1,
            }
            legacy_cells.append({
                "t_index": 0, "delta_index": delta_index,
                "regions": [legacy_row],
            })
            adaptive_cells[(0, delta_index)] = {
                "t_index": 0, "delta_index": delta_index,
                "slices": [adaptive_row],
            }
        legacy, _ = method_b._parent_references(
            legacy_cells,
            [{"region_name": "equivalent", "available": True}],
            (0.0, 1.0), resolved["support"], resolved["parent_reference"],
        )
        adaptive, _ = method_b._adaptive_parent_references(
            adaptive_cells, (0.0, 1.0), resolved["support"],
            resolved["parent_reference"],
        )
        for field in (
            "parent_reference_ratio", "parent_reference_uncertainty",
            "usable_delta_cell_count", "contributing_delta_indices",
            "combined_host_neff", "combined_baseline_neff",
            "parent_reference_status", "parent_reference_reason",
        ):
            self.assertEqual(legacy[0][field], adaptive[0][field])

    def test_adaptive_candidates_and_consistency_keep_single_and_inconsistent_evidence(self):
        def row(slice_id, value, sigma, side="low"):
            return {
                "slice_id": slice_id,
                "side": side,
                "parent_relative_status": "available",
                "parent_relative_ratio": value,
                "parent_relative_sigma": sigma,
            }

        empty = method_b._adaptive_candidate_and_consistency([])
        self.assertEqual(empty["adaptive_candidate_status"], "unavailable")
        self.assertIsNone(empty["adaptive_candidate"])
        self.assertEqual(empty["slice_consistency_status"], "not_evaluable")
        single = method_b._adaptive_candidate_and_consistency([row("one", 1.60, 0.25)])
        self.assertEqual(single["adaptive_candidate_status"], "available_single_slice")
        self.assertEqual(single["adaptive_candidate"], 1.60)
        self.assertEqual(single["slice_consistency_status"], "not_evaluable")
        self.assertIsNone(single["slice_consistency_ndf"])

        multi_rows = [row("low", 0.40, 0.04), row("high", 2.50, 0.25, "high")]
        multi = method_b._adaptive_candidate_and_consistency(multi_rows)
        self.assertEqual(multi["adaptive_candidate_status"], "available_multi_slice")
        self.assertTrue(math.isfinite(multi["adaptive_candidate"]))
        self.assertEqual(multi["slice_consistency_status"], "evaluated")
        self.assertEqual(multi["slice_consistency_ndf"], 1)
        self.assertGreater(multi["slice_consistency_chi2"], 1.0)
        self.assertEqual(len(multi["slice_consistency_log_pulls"]), 2)
        small = method_b._adaptive_candidate_and_consistency([row("small", 0.05, 0.01)])
        large = method_b._adaptive_candidate_and_consistency([row("large", 25.0, 2.0)])
        self.assertEqual(small["adaptive_candidate"], 0.05)
        self.assertEqual(large["adaptive_candidate"], 25.0)

    def test_adaptive_partition_merges_low_support_and_retains_unsupported_side(self):
        weak = _phase_a()
        for delta in (5.0, 15.0):
            for mm in (0.85, 0.95, 1.05, 1.25, 1.35, 1.45):
                _add_region_events(
                    weak, delta=delta, mm=mm,
                    pion_weights=(1.0,) * 5, host_weights=(1.0,) * 10,
                )
        weak_result = method_b.build_pion_hgcer_method_b(weak, config=_config())
        weak_partition = weak_result["adaptive_slice_diagnostic"]["t_partitions"][0]
        self.assertEqual(len(weak_partition["low_slices"]), 1)
        self.assertEqual(len(weak_partition["high_slices"]), 1)
        self.assertEqual(weak_partition["low_slices"][0]["mm_low"], MM_EDGES[0])
        self.assertEqual(
            weak_partition["low_slices"][0]["mm_high"],
            method_b.METHOD_B_PROTECTED_MM_LOW,
        )

        unsupported = _phase_a()
        for delta in (5.0, 15.0):
            for mm in (1.25, 1.35, 1.45):
                _add_region_events(unsupported, delta=delta, mm=mm)
        unsupported_result = method_b.build_pion_hgcer_method_b(
            unsupported, config=_config()
        )
        adaptive = unsupported_result["adaptive_slice_diagnostic"]
        partition = adaptive["t_partitions"][0]
        self.assertTrue(adaptive["available"])
        self.assertEqual(len(partition["low_slices"]), 1)
        self.assertEqual(
            partition["low_slices"][0]["partition_support_status"], "unavailable"
        )
        self.assertTrue(all(
            row["mm_high"] <= method_b.METHOD_B_PROTECTED_MM_LOW
            or row["mm_low"] >= method_b.METHOD_B_PROTECTED_MM_HIGH
            for row in partition["slices"]
        ))

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
        self.assertFalse(result["adaptive_slice_diagnostic"]["available"])
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
        resolver_call = source[
            source.index("resolve_pion_hgcer_method_b_config("):
            source.index("build_pion_hgcer_method_b(", source.index("resolve_pion_hgcer_method_b_config("))
        ]
        self.assertIn("phase_a_contract=pion_hgcer_event_contract", resolver_call)
        runtime_failure = source[
            source.index("except Exception as exc:", source.index("resolve_pion_hgcer_method_b_config(")):
            source.index('histDict["_pion_hgcer_method_b"]', source.index("resolve_pion_hgcer_method_b_config("))
        ]
        self.assertIn('"reason": getattr(', runtime_failure)
        self.assertIn('"runtime_method_b_configuration_exception"', runtime_failure)


if __name__ == "__main__":
    unittest.main()
