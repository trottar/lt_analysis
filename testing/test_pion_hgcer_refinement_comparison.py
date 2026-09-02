"""Pure-Python regression coverage for Phase-D.1/D.1.1, D.2, and D.3 contracts."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)


import pion_hgcer_refinement_checkpoint as checkpoint
import pion_hgcer_refinement_comparison as comparison


T_EDGES = [0.0, 1.0, 2.0]
DELTA_EDGES = [0.0, 10.0, 20.0]
METHOD_A_SUPPORT_THRESHOLDS = {
    "supported_positive_count": 25,
    "supported_low_count": 5,
    "supported_control_count": 5,
    "marginal_positive_count": 10,
    "minimum_control_count_for_ratio": 1,
    "minimum_low_count_for_ratio": 0,
}
WILSON_Z_95 = 1.959963984540054


def _cells(method, host_state):
    return [
        {
            "t_index": t_index,
            "delta_index": delta_index,
            "t_low": T_EDGES[t_index],
            "t_high": T_EDGES[t_index + 1],
            "delta_low": DELTA_EDGES[delta_index],
            "delta_high": DELTA_EDGES[delta_index + 1],
            "host_state": host_state,
            "method": method,
            "synthetic_observation": 10 * t_index + delta_index,
        }
        for t_index in range(len(T_EDGES) - 1)
        for delta_index in range(len(DELTA_EDGES) - 1)
    ]


def _phase_a(host_state):
    return {
        "status": "available",
        "available": True,
        "contract_fingerprint": "phase-a-contract-fingerprint",
        "coordinate_fingerprint": "phase-a-coordinate-fingerprint",
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
        "canonical_t_edges": list(T_EDGES),
        "delta_edges": list(DELTA_EDGES),
    }


def _phase_a_summary(host_state):
    return {
        "schema_version": "pion_hgcer_event_contract/v1",
        "fingerprint_schema_version": "pion_hgcer_event_contract_fingerprint/v2",
        "status": "available",
        "available": True,
        "contract_fingerprint": "phase-a-contract-fingerprint",
        "coordinate_fingerprint": "phase-a-coordinate-fingerprint",
        "canonical_t_edges": list(T_EDGES),
        "delta_edges": list(DELTA_EDGES),
        "pion_closure_passed": True,
        "host_closure_passed": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_restoration_applied": False,
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
    }


def _wilson_interval(successes, total):
    fraction = float(successes) / float(total)
    z2 = WILSON_Z_95 * WILSON_Z_95
    denominator = 1.0 + z2 / total
    center = (fraction + z2 / (2.0 * total)) / denominator
    spread = (
        WILSON_Z_95
        * (
            fraction * (1.0 - fraction) / total
            + z2 / (4.0 * total * total)
        )
        ** 0.5
        / denominator
    )
    return (
        0.0 if successes == 0 else max(0.0, center - spread),
        1.0 if successes == total else min(1.0, center + spread),
    )


def _method_a_support_class(positive, low, control):
    if (
        positive >= METHOD_A_SUPPORT_THRESHOLDS["supported_positive_count"]
        and low >= METHOD_A_SUPPORT_THRESHOLDS["supported_low_count"]
        and control >= METHOD_A_SUPPORT_THRESHOLDS["supported_control_count"]
    ):
        return "supported"
    if (
        positive >= METHOD_A_SUPPORT_THRESHOLDS["marginal_positive_count"]
        and control >= METHOD_A_SUPPORT_THRESHOLDS["minimum_control_count_for_ratio"]
        and low >= METHOD_A_SUPPORT_THRESHOLDS["minimum_low_count_for_ratio"]
    ):
        return "marginal"
    return "unsupported"


def _method_a_cell(t_index, delta_index, positive, low, control, host_state):
    support_class = _method_a_support_class(positive, low, control)
    available = support_class != "unsupported"
    f_low = float(low) / float(positive) if positive else None
    f_low_low, f_low_high = _wilson_interval(low, positive) if positive else (None, None)
    ratio = float(low) / float(control) if control else None
    ratio_low = (
        f_low_low / (1.0 - f_low_low)
        if f_low_low is not None and f_low_low < 1.0
        else None
    )
    ratio_high = (
        f_low_high / (1.0 - f_low_high)
        if f_low_high is not None and f_low_high < 1.0
        else None
    )
    return {
        "t_index": t_index,
        "delta_index": delta_index,
        "t_low": T_EDGES[t_index],
        "t_high": T_EDGES[t_index + 1],
        "delta_low": DELTA_EDGES[delta_index],
        "delta_high": DELTA_EDGES[delta_index + 1],
        "host_state": host_state,
        "method": "method_a",
        "synthetic_observation": 10 * t_index + delta_index,
        "prompt_positive_count": positive,
        "prompt_low_count": low,
        "prompt_control_count": control,
        "partition_closure_passed": positive == low + control,
        "f_low": f_low,
        "f_low_low": f_low_low if available else None,
        "f_low_high": f_low_high if available else None,
        "R_low_control": ratio,
        "R_low_control_low": ratio_low if available else None,
        "R_low_control_high": ratio_high if available else None,
        "signed_R_low_control": ratio + 0.125 if ratio is not None else None,
        "signed_R_low_control_sigma": 0.25 if ratio is not None else None,
        "prompt_vs_signed_status": "consistent",
        "nommcuts_vs_allcuts_status": "consistent",
        "support_class": support_class,
        "method_A_status": "available" if available else "unavailable",
        "method_A_reason": None if available else "support_insufficient",
    }


def _method_a_cells(host_state, counts=None):
    counts = counts or (
        (30, 10, 20),
        (20, 0, 20),
        (30, 15, 15),
        (20, 5, 15),
    )
    return [
        _method_a_cell(t_index, delta_index, *counts[2 * t_index + delta_index], host_state)
        for t_index in range(len(T_EDGES) - 1)
        for delta_index in range(len(DELTA_EDGES) - 1)
    ]


def _method_a(status, host_state):
    available = status == "available"
    return {
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_a_unavailable",
        "fingerprint": "method-a-fingerprint" if available else None,
        "host_state": host_state,
        "cells": _method_a_cells(host_state) if available else [],
    }


def _method_a_summary(status, host_state):
    available = status == "available"
    return {
        "schema_version": "pion_hgcer_method_a/v1",
        "method": "observed_positive_hgcer_response",
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_a_unavailable",
        "fingerprint": "method-a-fingerprint" if available else None,
        "phase_a_contract_fingerprint": "phase-a-contract-fingerprint",
        "coordinate_fingerprint": "phase-a-coordinate-fingerprint",
        "uncertainty_method": "wilson_95_percent",
        "support_thresholds": dict(METHOD_A_SUPPORT_THRESHOLDS),
        "t_edges": list(T_EDGES) if available else [],
        "delta_edges": list(DELTA_EDGES) if available else [],
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "host_state": host_state,
    }


METHOD_B_REGIONS = ("pi_n", "pi_delta")


def _method_b_parent_reference(t_index, region_name):
    base = {
        (0, "pi_n"): 1.10,
        (0, "pi_delta"): 0.90,
        (1, "pi_n"): 1.05,
        (1, "pi_delta"): 0.95,
    }[(t_index, region_name)]
    return {
        "t_index": t_index,
        "t_low": T_EDGES[t_index],
        "t_high": T_EDGES[t_index + 1],
        "region_name": region_name,
        "usable_delta_cell_count": 2,
        "contributing_delta_indices": [0, 1],
        "combined_host_abs_support": 40.0,
        "combined_host_sumw2": 16.0,
        "combined_host_neff": 20.0,
        "combined_baseline_abs_support": 40.0,
        "combined_baseline_sumw2": 16.0,
        "combined_baseline_neff": 20.0,
        "parent_reference_ratio": base,
        "parent_reference_uncertainty": 0.10,
        "parent_reference_status": "available",
        "parent_reference_reason": None,
        "weighting": "inverse_variance",
    }


def _method_b_region(t_index, delta_index, region_name):
    parent = _method_b_parent_reference(t_index, region_name)
    raw_ratio = parent["parent_reference_ratio"] * (1.0 + 0.05 * delta_index)
    return {
        "region_name": region_name,
        "available": True,
        "reason": None,
        "region_role": "pion_background",
        "window_source": "synthetic_frozen",
        "mm_low": -0.10 if region_name == "pi_n" else 0.10,
        "mm_high": 0.10 if region_name == "pi_n" else 0.30,
        "protected_signal_overlap": False,
        "support_status": "usable",
        "support_reason": None,
        "host_record_count": 20,
        "host_yield": 10.0,
        "host_sumw2": 4.0,
        "host_neff": 20.0,
        "baseline_record_count": 20,
        "baseline_pion_yield": 10.0,
        "baseline_pion_sumw2": 4.0,
        "baseline_pion_neff": 20.0,
        "baseline_pion_significance": 5.0,
        "residual": 0.0,
        "residual_sigma": 2.0,
        "fractional_residual": 0.0,
        "raw_ratio": raw_ratio,
        "raw_ratio_sigma": 0.20,
        "parent_reference_ratio": parent["parent_reference_ratio"],
        "parent_reference_sigma": parent["parent_reference_uncertainty"],
        "parent_relative_ratio": raw_ratio / parent["parent_reference_ratio"],
        "parent_relative_sigma": 0.20,
        "parent_relative_status": "available",
        "parent_relative_reason": None,
    }


def _method_b_cell(t_index, delta_index, host_state):
    return {
        "t_index": t_index,
        "t_low": T_EDGES[t_index],
        "t_high": T_EDGES[t_index + 1],
        "delta_index": delta_index,
        "delta_low": DELTA_EDGES[delta_index],
        "delta_high": DELTA_EDGES[delta_index + 1],
        "host_state": host_state,
        "host_record_count": 40,
        "baseline_record_count": 40,
        "regions": [
            _method_b_region(t_index, delta_index, region_name)
            for region_name in METHOD_B_REGIONS
        ],
        "region_consistency_status": "region_consistent",
        "region_consistency_reason": None,
        "shape_chi2": 1.0,
        "shape_ndf": 2,
        "shape_chi2_ndf": 0.5,
        "shape_max_abs_pull": 0.5,
        "shape_usable_bin_count": 2,
        "shape_status": "good",
        "shape_reason": None,
        "candidate_L_B": 1.0 + 0.05 * delta_index,
        "candidate_L_B_uncertainty": 0.15,
        "candidate_L_B_status": "available_multi_region",
        "method_B_status": "available",
        "method_B_reason": None,
        "bins": [{"not": "carried_by_d3"}],
        "underflow": {"not": "carried_by_d3"},
        "overflow": {"not": "carried_by_d3"},
    }


def _method_b_cells(host_state):
    return [
        _method_b_cell(t_index, delta_index, host_state)
        for t_index in range(len(T_EDGES) - 1)
        for delta_index in range(len(DELTA_EDGES) - 1)
    ]


def _method_b(status, host_state):
    available = status == "available"
    return {
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_b_unavailable",
        "fingerprint": "method-b-fingerprint" if available else None,
        "host_state": host_state,
        "cells": _method_b_cells(host_state) if available else [],
        "parent_region_references": (
            [
                _method_b_parent_reference(t_index, region_name)
                for t_index in range(len(T_EDGES) - 1)
                for region_name in METHOD_B_REGIONS
            ]
            if available
            else []
        ),
    }


def _method_b_summary(status, host_state):
    available = status == "available"
    return {
        "schema_version": "pion_hgcer_method_b/v1",
        "method": "same_t_pion_background_closure",
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_b_unavailable",
        "fingerprint": "method-b-fingerprint" if available else None,
        "phase_a_contract_fingerprint": "phase-a-contract-fingerprint",
        "coordinate_fingerprint": "phase-a-coordinate-fingerprint",
        "host_state": host_state,
        "t_edges": list(T_EDGES) if available else [],
        "delta_edges": list(DELTA_EDGES) if available else [],
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "interpolation_used": False,
    }


def _checkpoint(
    *,
    host_state="proton_cleaned",
    method_a_status="available",
    method_b_status="available",
    Q2="4p4",
    W="2p74",
):
    return checkpoint.build_pion_hgcer_refinement_checkpoint(
        setting={
            "kinematic_token": "Q4p4W2p74",
            "Q2": Q2,
            "W": W,
            "epsilon_setting": "low",
            "epsilon_filename_token": "lowe",
            "phi_setting": "Left",
            "particle_type": "kaon",
        },
        phase_a=_phase_a(host_state),
        phase_a_summary=_phase_a_summary(host_state),
        method_a=_method_a(method_a_status, host_state),
        method_a_summary=_method_a_summary(method_a_status, host_state),
        method_b=_method_b(method_b_status, host_state),
        method_b_summary=_method_b_summary(method_b_status, host_state),
    )


class PionHGCerRefinementComparisonTests(unittest.TestCase):
    def _assert_unavailable(self, payload, reason, stage=None):
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        self.assertEqual(result["schema_version"], comparison.COMPARISON_INPUT_SCHEMA_VERSION)
        self.assertEqual(result["status"], "unavailable")
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], reason)
        self.assertTrue(result["non_authoritative"])
        self.assertFalse(result["comparison_performed"])
        self.assertFalse(result["classification_performed"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["refinement_applied"])
        if stage is not None:
            self.assertEqual(result["diagnostic_stage"], stage)
        return result

    def test_valid_proton_cleaned_checkpoint_returns_complete_detached_contract(self):
        payload = _checkpoint(host_state="proton_cleaned")
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        self.assertEqual(
            set(result),
            {
                "schema_version",
                "status",
                "available",
                "reason",
                "diagnostic_stage",
                "source_checkpoint_schema_version",
                "source_checkpoint_payload_fingerprint",
                "setting",
                "canonical_t_edges",
                "delta_edges",
                "phase_a",
                "method_a",
                "method_b",
                "host_state_summary",
                "non_authoritative",
                "comparison_performed",
                "classification_performed",
                "production_objects_mutated",
                "refinement_applied",
            },
        )
        self.assertEqual(result["status"], "available")
        self.assertTrue(result["available"])
        self.assertEqual(result["reason"], None)
        self.assertEqual(result["diagnostic_stage"], "complete")
        self.assertEqual(
            result["source_checkpoint_schema_version"],
            "pion_hgcer_refinement_checkpoint/v1",
        )
        self.assertEqual(result["setting"], payload["setting"])
        self.assertEqual(result["setting"]["Q2"], "4p4")
        self.assertEqual(result["setting"]["W"], "2p74")
        self.assertIsInstance(result["setting"]["Q2"], str)
        self.assertIsInstance(result["setting"]["W"], str)
        self.assertEqual(result["canonical_t_edges"], T_EDGES)
        self.assertEqual(result["delta_edges"], DELTA_EDGES)
        self.assertEqual(result["phase_a"], payload["phase_a"])
        self.assertEqual(result["method_a"], payload["method_a"])
        self.assertEqual(result["method_b"], payload["method_b"])
        self.assertEqual(result["host_state_summary"], payload["host_state_summary"])
        self.assertEqual(len(result["source_checkpoint_payload_fingerprint"]), 64)
        self.assertTrue(result["non_authoritative"])
        self.assertFalse(result["comparison_performed"])
        self.assertFalse(result["classification_performed"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["refinement_applied"])

    def test_identity_no_proton_cleaning_host_is_a_valid_input(self):
        payload = _checkpoint(host_state="identity_no_proton_cleaning")
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        self.assertTrue(result["available"])
        self.assertEqual(result["phase_a"]["host_state"], "identity_no_proton_cleaning")
        self.assertEqual(
            result["host_state_summary"]["phase_a_host_state"],
            "identity_no_proton_cleaning",
        )

    def test_finite_numeric_setting_metadata_remains_valid_and_exact(self):
        payload = _checkpoint(Q2=4.4, W=2.74)
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        self.assertTrue(result["available"])
        self.assertEqual(result["setting"]["Q2"], 4.4)
        self.assertEqual(result["setting"]["W"], 2.74)
        self.assertIsInstance(result["setting"]["Q2"], float)
        self.assertIsInstance(result["setting"]["W"], float)
        token_result = comparison.build_pion_hgcer_comparison_input_contract(
            _checkpoint(Q2="4p4", W="2p74")
        )
        self.assertNotEqual(
            result["source_checkpoint_payload_fingerprint"],
            token_result["source_checkpoint_payload_fingerprint"],
        )
        integer_result = comparison.build_pion_hgcer_comparison_input_contract(
            _checkpoint(Q2=4, W=2)
        )
        self.assertTrue(integer_result["available"])
        self.assertIsInstance(integer_result["setting"]["Q2"], int)
        self.assertIsInstance(integer_result["setting"]["W"], int)

    def test_native_method_unavailability_is_not_a_contract_failure(self):
        for method_a_status, method_b_status in (
            ("unavailable", "available"),
            ("available", "unavailable"),
            ("unavailable", "unavailable"),
        ):
            with self.subTest(method_a=method_a_status, method_b=method_b_status):
                result = comparison.build_pion_hgcer_comparison_input_contract(
                    _checkpoint(
                        method_a_status=method_a_status,
                        method_b_status=method_b_status,
                    )
                )
                self.assertTrue(result["available"])
                self.assertEqual(result["method_a"]["status"], method_a_status)
                self.assertEqual(result["method_b"]["status"], method_b_status)
                self.assertFalse(result["comparison_performed"])
                self.assertFalse(result["classification_performed"])

    def test_unavailable_methods_may_retain_valid_diagnostic_cells(self):
        for method_name in ("method_a", "method_b"):
            with self.subTest(method=method_name):
                payload = _checkpoint(
                    method_a_status="unavailable" if method_name == "method_a" else "available",
                    method_b_status="unavailable" if method_name == "method_b" else "available",
                )
                payload[method_name]["cells"] = [
                    {
                        "t_index": 0,
                        "delta_index": 0,
                        "t_low": T_EDGES[0],
                        "t_high": T_EDGES[1],
                        "delta_low": DELTA_EDGES[0],
                        "delta_high": DELTA_EDGES[1],
                        "host_state": "proton_cleaned",
                        "diagnostic_only": True,
                    }
                ]
                result = comparison.build_pion_hgcer_comparison_input_contract(payload)
                self.assertTrue(result["available"])
                self.assertEqual(result[method_name]["status"], "unavailable")
                self.assertEqual(result[method_name]["cells"], payload[method_name]["cells"])
                self.assertFalse(result["comparison_performed"])
                self.assertFalse(result["classification_performed"])

    def test_source_checkpoint_is_immutable_and_returned_snapshot_has_no_aliases(self):
        payload = _checkpoint()
        before = copy.deepcopy(payload)
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        self.assertTrue(result["available"])
        self.assertEqual(payload, before)
        result["phase_a"]["summary"]["host_state"] = "changed_only_in_result"
        result["method_a"]["cells"][0]["synthetic_observation"] = -1
        result["method_b"]["parent_region_references"][0]["region_name"] = "changed"
        self.assertEqual(payload, before)

    def test_snapshot_fingerprint_is_canonical_deterministic_and_covers_all_input(self):
        payload = _checkpoint()
        result = comparison.build_pion_hgcer_comparison_input_contract(payload)
        expected = hashlib.sha256(
            json.dumps(
                payload,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            ).encode("ascii")
        ).hexdigest()
        self.assertEqual(result["source_checkpoint_payload_fingerprint"], expected)
        repeated = comparison.build_pion_hgcer_comparison_input_contract(
            copy.deepcopy(payload)
        )
        self.assertEqual(
            repeated["source_checkpoint_payload_fingerprint"],
            result["source_checkpoint_payload_fingerprint"],
        )
        changed = copy.deepcopy(payload)
        changed["method_a"]["cells"][0]["synthetic_observation"] += 1
        changed_result = comparison.build_pion_hgcer_comparison_input_contract(changed)
        self.assertTrue(changed_result["available"])
        self.assertNotEqual(
            changed_result["source_checkpoint_payload_fingerprint"],
            result["source_checkpoint_payload_fingerprint"],
        )

    def test_checkpoint_schema_authority_and_snapshot_failures_are_structured(self):
        payload = _checkpoint()
        payload["schema_version"] = "other/v1"
        self._assert_unavailable(payload, "checkpoint_schema_mismatch", "checkpoint_validation")
        payload = _checkpoint()
        payload["non_authoritative"] = False
        self._assert_unavailable(
            payload, "checkpoint_authority_contract_invalid", "checkpoint_validation"
        )
        payload = _checkpoint()
        del payload["method_a"]
        self._assert_unavailable(payload, "checkpoint_contract_invalid", "checkpoint_validation")
        self._assert_unavailable(
            {"opaque": object()}, "checkpoint_snapshot_invalid", "snapshot_validation"
        )

    def test_setting_contract_rejects_invalid_semantic_or_scalar_metadata(self):
        for key, value in (
            ("particle_type", "pion"),
            ("epsilon_filename_token", "highe"),
            ("epsilon_setting", ["low"]),
        ):
            with self.subTest(key=key, value=repr(value)):
                payload = _checkpoint()
                payload["setting"][key] = value
                self._assert_unavailable(payload, "setting_contract_invalid", "setting_validation")
        invalid_scalars = (
            None,
            True,
            False,
            ["4p4"],
            {"value": "4p4"},
            (),
            set(),
            "",
            "  ",
            " 4p4",
            "4p4 ",
            "../4p4",
            "4p4/extra",
            "4p4\\extra",
            "4p4:extra",
            "4p4\x00",
            float("nan"),
            float("inf"),
            -float("inf"),
        )
        for key in ("Q2", "W"):
            for value in invalid_scalars:
                with self.subTest(key=key, value=repr(value)):
                    payload = _checkpoint()
                    payload["setting"][key] = value
                    self._assert_unavailable(
                        payload, "setting_contract_invalid", "setting_validation"
                    )

    def test_phase_a_provenance_failures_are_narrow(self):
        cases = (
            ("contract_fingerprint", "other-contract", "phase_a_contract_fingerprint_mismatch"),
            ("coordinate_fingerprint", "other-coordinate", "coordinate_fingerprint_mismatch"),
            ("schema_version", "other/v1", "phase_a_schema_mismatch"),
            (
                "fingerprint_schema_version",
                "other/v1",
                "phase_a_fingerprint_schema_mismatch",
            ),
        )
        for key, value, reason in cases:
            with self.subTest(key=key):
                payload = _checkpoint()
                payload["phase_a"]["summary"][key] = value
                self._assert_unavailable(payload, reason, "phase_a_provenance")

    def test_canonical_and_delta_geometry_mismatches_are_never_repaired(self):
        payload = _checkpoint()
        payload["method_a"]["summary"]["t_edges"][1] = 1.25
        self._assert_unavailable(
            payload, "canonical_t_edges_mismatch", "canonical_geometry"
        )
        payload = _checkpoint()
        payload["method_b"]["summary"]["delta_edges"][1] = 12.5
        self._assert_unavailable(payload, "delta_edges_mismatch", "canonical_geometry")
        payload = _checkpoint()
        payload["canonical_t_edges"] = [0.0, 0.0]
        self._assert_unavailable(payload, "canonical_t_edges_invalid", "canonical_geometry")
        payload = _checkpoint()
        payload["delta_edges"] = [0.0, float("inf")]
        self._assert_unavailable(payload, "checkpoint_snapshot_invalid", "snapshot_validation")

    def test_host_and_source_target_provenance_failures_are_narrow(self):
        payload = _checkpoint()
        payload["host_state_summary"]["phase_a_host_state"] = "identity_no_proton_cleaning"
        self._assert_unavailable(payload, "host_state_mismatch", "host_provenance")
        payload = _checkpoint()
        payload["host_state_summary"]["source_target_state"] = "pre_proton_noRF"
        self._assert_unavailable(
            payload, "source_target_state_mismatch", "host_provenance"
        )
        payload = _checkpoint()
        payload["method_b"]["cells"][0]["host_state"] = "identity_no_proton_cleaning"
        self._assert_unavailable(
            payload, "method_b_cell_host_state_mismatch", "method_b_geometry"
        )

    def test_available_method_fingerprint_and_upstream_provenance_are_required(self):
        payload = _checkpoint()
        payload["method_a"]["fingerprint"] = None
        self._assert_unavailable(payload, "method_a_fingerprint_missing", "method_a_provenance")
        payload = _checkpoint()
        payload["method_b"]["fingerprint"] = None
        self._assert_unavailable(payload, "method_b_fingerprint_missing", "method_b_provenance")
        payload = _checkpoint()
        payload["method_a"]["summary"]["phase_a_contract_fingerprint"] = "other"
        self._assert_unavailable(
            payload, "method_a_phase_a_fingerprint_mismatch", "method_a_provenance"
        )
        payload = _checkpoint()
        payload["method_b"]["summary"]["coordinate_fingerprint"] = "other"
        self._assert_unavailable(
            payload, "method_b_coordinate_fingerprint_mismatch", "method_b_provenance"
        )

    def test_available_method_cell_grids_must_be_complete_unique_and_exact(self):
        cases = (
            ("removed", lambda cells: cells.pop(), "cell_grid_mismatch"),
            (
                "duplicated",
                lambda cells: cells.append(copy.deepcopy(cells[0])),
                "cell_grid_mismatch",
            ),
            (
                "t_index_out_of_range",
                lambda cells: cells[0].update(t_index=len(T_EDGES) - 1),
                "cell_grid_mismatch",
            ),
            (
                "delta_index_out_of_range",
                lambda cells: cells[0].update(delta_index=len(DELTA_EDGES) - 1),
                "cell_grid_mismatch",
            ),
            ("t_low", lambda cells: cells[0].update(t_low=0.25), "cell_geometry_mismatch"),
            ("t_high", lambda cells: cells[0].update(t_high=0.75), "cell_geometry_mismatch"),
            (
                "delta_low",
                lambda cells: cells[0].update(delta_low=0.25),
                "cell_geometry_mismatch",
            ),
            (
                "delta_high",
                lambda cells: cells[0].update(delta_high=0.75),
                "cell_geometry_mismatch",
            ),
        )
        for method_name in ("method_a", "method_b"):
            for case_name, mutate, suffix in cases:
                with self.subTest(method=method_name, case=case_name):
                    payload = _checkpoint()
                    mutate(payload[method_name]["cells"])
                    self._assert_unavailable(
                        payload,
                        "{}_{}".format(method_name, suffix),
                        "{}_geometry".format(method_name),
                    )

    def test_method_b_parent_reference_index_is_bounded_without_recalculation(self):
        payload = _checkpoint()
        payload["method_b"]["parent_region_references"][0]["t_index"] = len(T_EDGES) - 1
        self._assert_unavailable(
            payload,
            "method_b_parent_reference_geometry_mismatch",
            "method_b_geometry",
        )

    def test_comparison_module_has_no_runtime_or_analysis_dependencies(self):
        source = (
            REPO_ROOT / "src/cuts/pion_hgcer_refinement_comparison.py"
        ).read_text(encoding="utf-8")
        for forbidden in (
            "import ROOT",
            "import uproot",
            "root_numpy",
            "proton_contamination_weights",
            "pion_component_subtraction",
            "pion_component_fits",
            "pion_t_bin_parents",
            "import rand_sub",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)


class PionHGCerMethodAComparisonTests(unittest.TestCase):
    def _comparison_input(self, **checkpoint_kwargs):
        result = comparison.build_pion_hgcer_comparison_input_contract(
            _checkpoint(**checkpoint_kwargs)
        )
        self.assertTrue(result["available"])
        return result

    def _assert_unavailable(self, payload, reason, stage=None):
        result = comparison.build_pion_hgcer_method_a_comparison(payload)
        self.assertEqual(
            result["schema_version"], comparison.METHOD_A_COMPARISON_SCHEMA_VERSION
        )
        self.assertEqual(result["method"], "method_a_same_t_comparison_representation")
        self.assertEqual(result["status"], "unavailable")
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], reason)
        self.assertTrue(result["non_authoritative"])
        self.assertFalse(result["method_b_numerical_dependency"])
        self.assertFalse(result["comparison_performed"])
        self.assertFalse(result["classification_performed"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["refinement_applied"])
        if stage is not None:
            self.assertEqual(result["diagnostic_stage"], stage)
        return result

    @staticmethod
    def _cell(result, t_index, delta_index):
        return next(
            cell
            for cell in result["cells"]
            if cell["t_index"] == t_index and cell["delta_index"] == delta_index
        )

    @staticmethod
    def _replace_cell(contract, t_index, delta_index, positive, low, control):
        cells = contract["method_a"]["cells"]
        index = next(
            index
            for index, cell in enumerate(cells)
            if cell["t_index"] == t_index and cell["delta_index"] == delta_index
        )
        cells[index] = _method_a_cell(
            t_index, delta_index, positive, low, control, "proton_cleaned"
        )

    @staticmethod
    def _calculation_view(result):
        return {
            "parents": [
                {
                    key: parent[key]
                    for key in (
                        "t_index",
                        "contributing_delta_indices",
                        "contributing_delta_cell_count",
                        "prompt_positive_count",
                        "prompt_low_count",
                        "prompt_control_count",
                        "f_low",
                        "f_low_low",
                        "f_low_high",
                        "R_low_control",
                        "R_low_control_low",
                        "R_low_control_high",
                        "support_class",
                        "parent_reference_status",
                        "parent_reference_reason",
                    )
                }
                for parent in result["parent_references"]
            ],
            "cells": [
                {
                    key: cell[key]
                    for key in (
                        "t_index",
                        "delta_index",
                        "parent_reference_R_low_control",
                        "parent_reference_R_low_control_low",
                        "parent_reference_R_low_control_high",
                        "parent_reference_status",
                        "method_a_comparison_candidate",
                        "method_a_comparison_candidate_low",
                        "method_a_comparison_candidate_high",
                        "method_a_comparison_candidate_status",
                        "method_a_comparison_candidate_reason",
                    )
                }
                for cell in result["cells"]
            ],
        }

    def test_successful_same_t_representation_has_required_detached_shape(self):
        for host_state in ("proton_cleaned", "identity_no_proton_cleaning"):
            with self.subTest(host_state=host_state):
                input_contract = self._comparison_input(host_state=host_state)
                result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
                self.assertEqual(
                    set(result),
                    {
                        "schema_version",
                        "method",
                        "status",
                        "available",
                        "reason",
                        "diagnostic_stage",
                        "source_checkpoint_payload_fingerprint",
                        "source_method_a_payload_fingerprint",
                        "phase_a_contract_fingerprint",
                        "coordinate_fingerprint",
                        "method_a_fingerprint",
                        "canonical_t_edges",
                        "delta_edges",
                        "support_thresholds",
                        "parent_definition",
                        "uncertainty_definition",
                        "candidate_definition",
                        "parent_references",
                        "cells",
                        "fingerprint_inputs",
                        "fingerprint",
                        "non_authoritative",
                        "method_b_numerical_dependency",
                        "comparison_performed",
                        "classification_performed",
                        "production_objects_mutated",
                        "refinement_applied",
                    },
                )
                self.assertTrue(result["available"])
                self.assertEqual(result["reason"], None)
                self.assertEqual(result["diagnostic_stage"], "complete")
                self.assertEqual(result["support_thresholds"], METHOD_A_SUPPORT_THRESHOLDS)
                self.assertEqual(result["parent_definition"], "same_t_aggregate_prompt_low_control_counts")
                self.assertEqual(
                    result["uncertainty_definition"],
                    "wilson_95_percent_aggregated_prompt_counts",
                )
                self.assertEqual(
                    result["candidate_definition"],
                    "same_t_parent_ratio_with_ratio_envelope_bounds",
                )
                self.assertFalse(result["method_b_numerical_dependency"])
                self.assertEqual(len(result["parent_references"]), 2)
                self.assertEqual(len(result["cells"]), 4)
                self.assertEqual(len(result["source_method_a_payload_fingerprint"]), 64)
                self.assertEqual(len(result["fingerprint"]), 64)

    def test_input_and_nested_result_values_are_isolated(self):
        input_contract = self._comparison_input()
        before = copy.deepcopy(input_contract)
        result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        self.assertEqual(input_contract, before)
        result["parent_references"][0]["contributing_delta_indices"].append(99)
        result["cells"][0]["signed_R_low_control"] = -1.0
        result["support_thresholds"]["supported_positive_count"] = -1
        self.assertEqual(input_contract, before)

    def test_method_b_content_is_not_read_or_fingerprinted(self):
        input_contract = self._comparison_input()
        baseline = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        changed = copy.deepcopy(input_contract)
        changed["method_b"] = {"arbitrary": object(), "cells": ["anything"]}
        result = comparison.build_pion_hgcer_method_a_comparison(changed)
        for key in (
            "source_method_a_payload_fingerprint",
            "parent_references",
            "cells",
            "fingerprint",
        ):
            self.assertEqual(result[key], baseline[key])

    def test_same_t_parent_uses_aggregate_counts_not_ratio_averaging(self):
        input_contract = self._comparison_input()
        self._replace_cell(input_contract, 0, 0, 30, 5, 25)
        self._replace_cell(input_contract, 0, 1, 20, 10, 10)
        result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        parent = result["parent_references"][0]
        self.assertEqual(parent["prompt_low_count"], 15)
        self.assertEqual(parent["prompt_control_count"], 35)
        self.assertEqual(parent["R_low_control"], 15.0 / 35.0)
        self.assertNotEqual(parent["R_low_control"], (0.2 + 1.0) / 2.0)
        self.assertEqual(parent["contributing_delta_indices"], [0, 1])

    def test_parent_wilson_bounds_and_asymmetric_candidates_match_definition(self):
        result = comparison.build_pion_hgcer_method_a_comparison(self._comparison_input())
        parent = result["parent_references"][0]
        expected_low, expected_high = _wilson_interval(10, 50)
        self.assertAlmostEqual(parent["f_low"], 10.0 / 50.0, places=15)
        self.assertAlmostEqual(parent["f_low_low"], expected_low, places=15)
        self.assertAlmostEqual(parent["f_low_high"], expected_high, places=15)
        self.assertAlmostEqual(
            parent["R_low_control_low"], expected_low / (1.0 - expected_low), places=15
        )
        self.assertAlmostEqual(
            parent["R_low_control_high"], expected_high / (1.0 - expected_high), places=15
        )
        cell = self._cell(result, 0, 0)
        self.assertAlmostEqual(
            cell["method_a_comparison_candidate"],
            cell["R_low_control"] / parent["R_low_control"],
            places=15,
        )
        self.assertAlmostEqual(
            cell["method_a_comparison_candidate_low"],
            cell["R_low_control_low"] / parent["R_low_control_high"],
            places=15,
        )
        self.assertAlmostEqual(
            cell["method_a_comparison_candidate_high"],
            cell["R_low_control_high"] / parent["R_low_control_low"],
            places=15,
        )
        self.assertLessEqual(
            cell["method_a_comparison_candidate_low"],
            cell["method_a_comparison_candidate"],
        )
        self.assertLessEqual(
            cell["method_a_comparison_candidate"],
            cell["method_a_comparison_candidate_high"],
        )
        self.assertEqual(
            cell["candidate_interval_method"], "ratio_envelope_from_wilson_bounds"
        )
        self.assertEqual(
            cell["candidate_covariance_treatment"], "shared_parent_covariance_not_modeled"
        )

    def test_supported_marginal_zero_low_and_unsupported_cells_follow_native_status(self):
        input_contract = self._comparison_input()
        result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        self.assertEqual(result["parent_references"][0]["parent_reference_status"], "available")
        supported = self._cell(result, 0, 0)
        zero_low = self._cell(result, 0, 1)
        self.assertEqual(supported["method_a_comparison_candidate_status"], "available")
        self.assertEqual(zero_low["method_a_comparison_candidate"], 0.0)
        self.assertEqual(zero_low["method_a_comparison_candidate_low"], 0.0)
        self.assertGreater(zero_low["method_a_comparison_candidate_high"], 0.0)
        self.assertEqual(zero_low["method_a_comparison_candidate_status"], "marginal")

        unsupported_input = self._comparison_input()
        self._replace_cell(unsupported_input, 1, 1, 9, 1, 8)
        unsupported_result = comparison.build_pion_hgcer_method_a_comparison(
            unsupported_input
        )
        unsupported = self._cell(unsupported_result, 1, 1)
        self.assertEqual(unsupported["support_class"], "unsupported")
        self.assertIsNone(unsupported["method_a_comparison_candidate"])
        self.assertIsNone(unsupported["method_a_comparison_candidate_low"])
        self.assertIsNone(unsupported["method_a_comparison_candidate_high"])
        self.assertEqual(unsupported["method_a_comparison_candidate_status"], "unavailable")
        self.assertEqual(
            unsupported["method_a_comparison_candidate_reason"], "cell_support_insufficient"
        )

    def test_unusable_parents_preserve_diagnostics_without_candidates(self):
        input_contract = self._comparison_input()
        self._replace_cell(input_contract, 0, 1, 9, 1, 8)
        result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        parent = result["parent_references"][0]
        self.assertEqual(parent["contributing_delta_cell_count"], 1)
        self.assertEqual(parent["prompt_positive_count"], 30)
        self.assertEqual(parent["parent_reference_status"], "unavailable")
        self.assertEqual(
            parent["parent_reference_reason"], "insufficient_contributing_delta_cells"
        )
        cell = self._cell(result, 0, 0)
        self.assertIsNone(cell["method_a_comparison_candidate"])
        self.assertEqual(
            cell["method_a_comparison_candidate_reason"],
            "insufficient_contributing_delta_cells",
        )

        all_unusable_input = self._comparison_input()
        for t_index in range(2):
            for delta_index in range(2):
                self._replace_cell(all_unusable_input, t_index, delta_index, 9, 1, 8)
        all_unusable = comparison.build_pion_hgcer_method_a_comparison(
            all_unusable_input
        )
        self.assertTrue(all_unusable["available"])
        self.assertTrue(
            all(
                parent["support_class"] == "unsupported"
                and parent["parent_reference_status"] == "unavailable"
                for parent in all_unusable["parent_references"]
            )
        )
        self.assertTrue(
            all(
                cell["method_a_comparison_candidate"] is None
                and cell["method_a_comparison_candidate_status"] == "unavailable"
                for cell in all_unusable["cells"]
            )
        )

    def test_marginal_parent_remains_usable_and_marginal(self):
        input_contract = self._comparison_input()
        self._replace_cell(input_contract, 0, 0, 10, 1, 9)
        self._replace_cell(input_contract, 0, 1, 10, 1, 9)
        result = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        parent = result["parent_references"][0]
        self.assertEqual(parent["support_class"], "marginal")
        self.assertEqual(parent["parent_reference_status"], "marginal")
        self.assertEqual(
            self._cell(result, 0, 0)["method_a_comparison_candidate_status"], "marginal"
        )

    def test_changing_one_t_bin_does_not_change_another(self):
        input_contract = self._comparison_input()
        baseline = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        changed = copy.deepcopy(input_contract)
        self._replace_cell(changed, 1, 0, 40, 20, 20)
        result = comparison.build_pion_hgcer_method_a_comparison(changed)
        self.assertEqual(result["parent_references"][0], baseline["parent_references"][0])
        self.assertEqual(
            [cell for cell in result["cells"] if cell["t_index"] == 0],
            [cell for cell in baseline["cells"] if cell["t_index"] == 0],
        )
        self.assertNotEqual(result["parent_references"][1], baseline["parent_references"][1])

    def test_cross_checks_and_signed_fields_do_not_change_calculations(self):
        input_contract = self._comparison_input()
        baseline = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        for key, value in (
            ("prompt_vs_signed_status", "inconsistent"),
            ("nommcuts_vs_allcuts_status", "inconsistent"),
            ("signed_R_low_control", 987.0),
            ("signed_R_low_control_sigma", 654.0),
        ):
            with self.subTest(key=key):
                changed = copy.deepcopy(input_contract)
                changed["method_a"]["cells"][0][key] = value
                result = comparison.build_pion_hgcer_method_a_comparison(changed)
                self.assertEqual(
                    self._calculation_view(result), self._calculation_view(baseline)
                )
                self.assertNotEqual(
                    result["source_method_a_payload_fingerprint"],
                    baseline["source_method_a_payload_fingerprint"],
                )
                self.assertNotEqual(result["fingerprint"], baseline["fingerprint"])
                self.assertEqual(result["cells"][0][key], value)

    def test_method_a_fingerprints_are_deterministic_and_method_a_specific(self):
        input_contract = self._comparison_input()
        first = comparison.build_pion_hgcer_method_a_comparison(input_contract)
        second = comparison.build_pion_hgcer_method_a_comparison(copy.deepcopy(input_contract))
        for key in (
            "source_method_a_payload_fingerprint",
            "fingerprint",
            "parent_references",
            "cells",
        ):
            self.assertEqual(first[key], second[key])
        changed = copy.deepcopy(input_contract)
        self._replace_cell(changed, 0, 0, 32, 8, 24)
        changed_result = comparison.build_pion_hgcer_method_a_comparison(changed)
        self.assertNotEqual(
            changed_result["source_method_a_payload_fingerprint"],
            first["source_method_a_payload_fingerprint"],
        )
        self.assertNotEqual(changed_result["fingerprint"], first["fingerprint"])

    def test_structured_top_level_and_provenance_failures(self):
        self._assert_unavailable([], "comparison_input_contract_invalid")
        unavailable_input = comparison.build_pion_hgcer_comparison_input_contract(
            _checkpoint(method_a_status="unavailable")
        )
        self._assert_unavailable(unavailable_input, "method_a_unavailable", "method_a_provenance")
        malformed_d1 = copy.deepcopy(self._comparison_input())
        malformed_d1["status"] = "unavailable"
        malformed_d1["available"] = False
        self._assert_unavailable(malformed_d1, "comparison_input_unavailable")
        cases = (
            (
                "authority",
                lambda value: value.update(non_authoritative=False),
                "comparison_input_authority_contract_invalid",
            ),
            (
                "phase_fingerprint",
                lambda value: value["phase_a"].update(contract_fingerprint=None),
                "phase_a_fingerprint_missing",
            ),
            (
                "coordinate_fingerprint",
                lambda value: value["phase_a"].update(coordinate_fingerprint=None),
                "coordinate_fingerprint_missing",
            ),
            (
                "method_fingerprint",
                lambda value: value["method_a"].update(fingerprint=None),
                "method_a_fingerprint_missing",
            ),
            (
                "summary",
                lambda value: value["method_a"]["summary"].update(schema_version="other/v1"),
                "method_a_summary_invalid",
            ),
            (
                "uncertainty",
                lambda value: value["method_a"]["summary"].update(uncertainty_method="other"),
                "method_a_uncertainty_contract_invalid",
            ),
            (
                "thresholds",
                lambda value: value["method_a"]["summary"]["support_thresholds"].pop("supported_low_count"),
                "method_a_support_thresholds_invalid",
            ),
        )
        for name, mutate, reason in cases:
            with self.subTest(name=name):
                payload = self._comparison_input()
                mutate(payload)
                self._assert_unavailable(payload, reason)

    def test_structured_cell_and_geometry_failures_are_not_repaired(self):
        cases = (
            (
                "missing_required_field",
                lambda value: value["method_a"]["cells"][0].pop("signed_R_low_control"),
                "method_a_cells_invalid",
            ),
            (
                "missing_cell",
                lambda value: value["method_a"]["cells"].pop(),
                "method_a_cell_count_contract_invalid",
            ),
            (
                "duplicate_cell",
                lambda value: value["method_a"]["cells"].append(copy.deepcopy(value["method_a"]["cells"][0])),
                "method_a_cell_count_contract_invalid",
            ),
            (
                "invalid_count",
                lambda value: value["method_a"]["cells"][0].update(prompt_low_count=-1),
                "method_a_cell_count_contract_invalid",
            ),
            (
                "invalid_status",
                lambda value: value["method_a"]["cells"][0].update(method_A_status="unavailable"),
                "method_a_cell_status_contract_invalid",
            ),
            (
                "invalid_ratio",
                lambda value: value["method_a"]["cells"][0].update(R_low_control=99.0),
                "method_a_cell_ratio_contract_invalid",
            ),
            (
                "cell_geometry",
                lambda value: value["method_a"]["cells"][0].update(t_high=0.5),
                "canonical_geometry_mismatch",
            ),
            (
                "summary_geometry",
                lambda value: value["method_a"]["summary"]["t_edges"].__setitem__(1, 1.5),
                "canonical_geometry_mismatch",
            ),
            (
                "invalid_edges",
                lambda value: value.update(canonical_t_edges=[0.0, 0.0]),
                "canonical_geometry_invalid",
            ),
        )
        for name, mutate, reason in cases:
            with self.subTest(name=name):
                payload = self._comparison_input()
                mutate(payload)
                self._assert_unavailable(payload, reason)


class PionHGCerMethodBComparisonTests(unittest.TestCase):
    def _comparison_input(self, **checkpoint_kwargs):
        result = comparison.build_pion_hgcer_comparison_input_contract(
            _checkpoint(**checkpoint_kwargs)
        )
        self.assertTrue(result["available"])
        return result

    def _assert_unavailable(self, payload, reason, stage=None):
        result = comparison.build_pion_hgcer_method_b_comparison(payload)
        self.assertEqual(
            result["schema_version"], comparison.METHOD_B_COMPARISON_SCHEMA_VERSION
        )
        self.assertEqual(result["method"], "method_b_comparison_representation")
        self.assertEqual(result["status"], "unavailable")
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], reason)
        self.assertTrue(result["non_authoritative"])
        self.assertFalse(result["method_a_numerical_dependency"])
        self.assertFalse(result["comparison_performed"])
        self.assertFalse(result["classification_performed"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["refinement_applied"])
        if stage is not None:
            self.assertEqual(result["diagnostic_stage"], stage)
        return result

    @staticmethod
    def _cell(result, t_index, delta_index):
        return next(
            cell
            for cell in result["cells"]
            if cell["t_index"] == t_index and cell["delta_index"] == delta_index
        )

    @staticmethod
    def _make_candidate_unavailable(cell, candidate_status, method_status, reason):
        cell["candidate_L_B"] = None
        cell["candidate_L_B_uncertainty"] = None
        cell["candidate_L_B_status"] = candidate_status
        cell["method_B_status"] = method_status
        cell["method_B_reason"] = reason
        if candidate_status == "region_marginal":
            cell["region_consistency_status"] = "region_marginal"
            cell["region_consistency_reason"] = reason
        elif candidate_status == "region_inconsistent":
            cell["region_consistency_status"] = "region_inconsistent"
            cell["region_consistency_reason"] = reason
        elif candidate_status == "shape_poor_veto":
            cell["region_consistency_status"] = "region_consistent"
            cell["shape_status"] = "poor"
            cell["shape_reason"] = reason
        else:
            cell["region_consistency_status"] = "insufficient_regions"
            cell["region_consistency_reason"] = reason

    def test_successful_method_b_projection_has_compact_detached_shape(self):
        for host_state in ("proton_cleaned", "identity_no_proton_cleaning"):
            with self.subTest(host_state=host_state):
                input_contract = self._comparison_input(host_state=host_state)
                result = comparison.build_pion_hgcer_method_b_comparison(input_contract)
                self.assertEqual(
                    set(result),
                    {
                        "schema_version",
                        "method",
                        "status",
                        "available",
                        "reason",
                        "diagnostic_stage",
                        "source_checkpoint_payload_fingerprint",
                        "source_method_b_payload_fingerprint",
                        "phase_a_contract_fingerprint",
                        "coordinate_fingerprint",
                        "method_b_fingerprint",
                        "canonical_t_edges",
                        "delta_edges",
                        "host_state",
                        "source_target_state",
                        "parent_region_references",
                        "cells",
                        "fingerprint_inputs",
                        "fingerprint",
                        "non_authoritative",
                        "method_a_numerical_dependency",
                        "comparison_performed",
                        "classification_performed",
                        "production_objects_mutated",
                        "refinement_applied",
                    },
                )
                self.assertTrue(result["available"])
                self.assertEqual(len(result["cells"]), 4)
                self.assertEqual(len(result["parent_region_references"]), 4)
                self.assertFalse(result["method_a_numerical_dependency"])
                expected_source_fingerprint = hashlib.sha256(
                    json.dumps(
                        input_contract["method_b"],
                        sort_keys=True,
                        separators=(",", ":"),
                        ensure_ascii=True,
                        allow_nan=False,
                    ).encode("ascii")
                ).hexdigest()
                self.assertEqual(
                    result["source_method_b_payload_fingerprint"],
                    expected_source_fingerprint,
                )
                self.assertEqual(len(result["source_method_b_payload_fingerprint"]), 64)
                self.assertEqual(len(result["fingerprint"]), 64)
                cell = self._cell(result, 0, 0)
                self.assertNotIn("bins", cell)
                self.assertNotIn("underflow", cell)
                self.assertNotIn("overflow", cell)
                self.assertEqual(
                    cell["method_b_comparison_candidate"], cell["candidate_L_B"]
                )
                self.assertEqual(
                    cell["method_b_comparison_candidate_uncertainty"],
                    cell["candidate_L_B_uncertainty"],
                )
                self.assertEqual(
                    cell["method_b_comparison_candidate_status"],
                    "available_multi_region",
                )

    def test_input_and_nested_return_values_are_isolated(self):
        input_contract = self._comparison_input()
        before = copy.deepcopy(input_contract)
        result = comparison.build_pion_hgcer_method_b_comparison(input_contract)
        self.assertEqual(input_contract, before)
        result["cells"][0]["regions"][0]["raw_ratio"] = -999.0
        result["parent_region_references"][0]["contributing_delta_indices"].append(99)
        self.assertEqual(input_contract, before)

    def test_method_a_is_completely_opaque_to_d3(self):
        input_contract = self._comparison_input()
        baseline = comparison.build_pion_hgcer_method_b_comparison(input_contract)
        changed = copy.deepcopy(input_contract)
        changed["method_a"] = {"opaque": object(), "cells": [object()]}
        result = comparison.build_pion_hgcer_method_b_comparison(changed)
        for key in (
            "source_method_b_payload_fingerprint",
            "parent_region_references",
            "cells",
            "fingerprint",
        ):
            self.assertEqual(result[key], baseline[key])
        source = (REPO_ROOT / "src/cuts/pion_hgcer_refinement_comparison.py").read_text(
            encoding="utf-8"
        )
        d3_source = source[
            source.index("def _method_b_comparison_unavailable") : source.index(
                "\n__all__"
            )
        ]
        for forbidden in (
            "METHOD_A_COMPARISON_SCHEMA_VERSION",
            "build_pion_hgcer_method_a_comparison",
            "method_a_comparison_candidate",
        ):
            self.assertNotIn(forbidden, d3_source)

    def test_candidate_and_native_statuses_are_preserved_without_promotion(self):
        cases = (
            ("single_region_only", "unavailable", "only_one_region"),
            ("unavailable", "unavailable", "no_usable_regions"),
            ("region_marginal", "marginal", "marginal_agreement"),
            ("region_inconsistent", "internally_inconsistent", "disjoint_regions"),
            ("shape_poor_veto", "shape_inconsistent", "poor_shape"),
        )
        for candidate_status, method_status, reason in cases:
            with self.subTest(candidate_status=candidate_status):
                input_contract = self._comparison_input()
                cell = input_contract["method_b"]["cells"][0]
                self._make_candidate_unavailable(
                    cell, candidate_status, method_status, reason
                )
                result = comparison.build_pion_hgcer_method_b_comparison(input_contract)
                projected = self._cell(result, 0, 0)
                self.assertIsNone(projected["method_b_comparison_candidate"])
                self.assertIsNone(
                    projected["method_b_comparison_candidate_uncertainty"]
                )
                self.assertEqual(
                    projected["method_b_comparison_candidate_status"], candidate_status
                )
                self.assertEqual(
                    projected["method_b_comparison_candidate_reason"], reason
                )
                self.assertEqual(projected["method_B_status"], method_status)

    def test_valid_top_level_method_b_may_have_all_local_candidates_unavailable(self):
        input_contract = self._comparison_input()
        for cell in input_contract["method_b"]["cells"]:
            self._make_candidate_unavailable(
                cell, "single_region_only", "unavailable", "only_one_region"
            )
        result = comparison.build_pion_hgcer_method_b_comparison(input_contract)
        self.assertTrue(result["available"])
        self.assertTrue(
            all(cell["method_b_comparison_candidate"] is None for cell in result["cells"])
        )

    def test_method_b_unavailable_is_structured_without_changing_d1(self):
        input_contract = self._comparison_input(method_b_status="unavailable")
        self._assert_unavailable(input_contract, "method_b_unavailable", "method_b_provenance")

    def test_parent_reference_row_linkage_is_exact_for_ratio_and_sigma(self):
        for key, value in (
            ("parent_reference_ratio", 999.0),
            ("parent_reference_sigma", 999.0),
        ):
            with self.subTest(key=key):
                input_contract = self._comparison_input()
                input_contract["method_b"]["cells"][0]["regions"][0][key] = value
                self._assert_unavailable(
                    input_contract,
                    "method_b_parent_reference_linkage_mismatch",
                    "method_b_parent_references",
                )

    def test_parent_reference_geometry_and_status_contracts_are_not_repaired(self):
        cases = (
            (
                "missing",
                lambda value: value["method_b"]["parent_region_references"].pop(),
                "method_b_parent_references_invalid",
            ),
            (
                "duplicate",
                lambda value: value["method_b"]["parent_region_references"].append(
                    copy.deepcopy(value["method_b"]["parent_region_references"][0])
                ),
                "method_b_parent_references_invalid",
            ),
            (
                "out_of_range_t",
                lambda value: value["method_b"]["parent_region_references"][0].update(t_index=2),
                "method_b_parent_reference_geometry_invalid",
            ),
            (
                "wrong_t_high",
                lambda value: value["method_b"]["parent_region_references"][0].update(t_high=0.5),
                "method_b_parent_reference_geometry_invalid",
            ),
            (
                "unknown_region",
                lambda value: value["method_b"]["parent_region_references"][0].update(region_name="unknown"),
                "method_b_parent_reference_geometry_invalid",
            ),
            (
                "duplicate_delta",
                lambda value: value["method_b"]["parent_region_references"][0].update(contributing_delta_indices=[0, 0]),
                "method_b_parent_reference_geometry_invalid",
            ),
            (
                "count_mismatch",
                lambda value: value["method_b"]["parent_region_references"][0].update(usable_delta_cell_count=1),
                "method_b_parent_reference_geometry_invalid",
            ),
            (
                "available_null_ratio",
                lambda value: value["method_b"]["parent_region_references"][0].update(parent_reference_ratio=None),
                "method_b_parent_reference_status_invalid",
            ),
            (
                "unavailable_missing_reason",
                lambda value: value["method_b"]["parent_region_references"][0].update(parent_reference_status="unavailable", parent_reference_reason=None),
                "method_b_parent_reference_status_invalid",
            ),
        )
        for name, mutate, reason in cases:
            with self.subTest(name=name):
                input_contract = self._comparison_input()
                mutate(input_contract)
                self._assert_unavailable(input_contract, reason)

    def test_region_shape_candidate_and_geometry_failures_are_structured(self):
        cases = (
            (
                "missing_region_field",
                lambda value: value["method_b"]["cells"][0]["regions"][0].pop("raw_ratio"),
                "method_b_region_contract_invalid",
            ),
            (
                "duplicate_region",
                lambda value: value["method_b"]["cells"][0]["regions"].append(
                    copy.deepcopy(value["method_b"]["cells"][0]["regions"][0])
                ),
                "method_b_region_contract_invalid",
            ),
            (
                "region_set_mismatch",
                lambda value: value["method_b"]["cells"][0]["regions"][1].update(region_name="other"),
                "method_b_region_contract_invalid",
            ),
            (
                "invalid_support",
                lambda value: value["method_b"]["cells"][0]["regions"][0].update(support_status="other"),
                "method_b_region_status_invalid",
            ),
            (
                "available_parent_relative_null",
                lambda value: value["method_b"]["cells"][0]["regions"][0].update(parent_relative_ratio=None),
                "method_b_parent_relative_contract_invalid",
            ),
            (
                "invalid_shape",
                lambda value: value["method_b"]["cells"][0].update(shape_status="other"),
                "method_b_status_contract_invalid",
            ),
            (
                "negative_ndf",
                lambda value: value["method_b"]["cells"][0].update(shape_ndf=-1),
                "method_b_shape_contract_invalid",
            ),
            (
                "candidate_value",
                lambda value: value["method_b"]["cells"][0].update(candidate_L_B=None),
                "method_b_candidate_contract_invalid",
            ),
            (
                "cell_host",
                lambda value: value["method_b"]["cells"][0].update(host_state="identity_no_proton_cleaning"),
                "method_b_cell_host_state_mismatch",
            ),
            (
                "cell_geometry",
                lambda value: value["method_b"]["cells"][0].update(delta_high=15.0),
                "method_b_cell_geometry_invalid",
            ),
        )
        for name, mutate, reason in cases:
            with self.subTest(name=name):
                input_contract = self._comparison_input()
                mutate(input_contract)
                self._assert_unavailable(input_contract, reason)

    def test_provenance_geometry_and_authority_failures_are_structured(self):
        self._assert_unavailable([], "comparison_input_contract_invalid")
        unavailable_d1 = copy.deepcopy(self._comparison_input())
        unavailable_d1.update(status="unavailable", available=False)
        self._assert_unavailable(unavailable_d1, "comparison_input_unavailable")
        cases = (
            (
                "authority",
                lambda value: value.update(non_authoritative=False),
                "comparison_input_authority_contract_invalid",
            ),
            (
                "phase_fingerprint",
                lambda value: value["phase_a"].update(contract_fingerprint=None),
                "phase_a_fingerprint_missing",
            ),
            (
                "coordinate",
                lambda value: value["phase_a"].update(coordinate_fingerprint=None),
                "coordinate_fingerprint_missing",
            ),
            (
                "host",
                lambda value: value["host_state_summary"].update(method_b_host_state="other"),
                "host_state_invalid",
            ),
            (
                "source_target",
                lambda value: value["phase_a"].update(source_target_state="other"),
                "source_target_state_mismatch",
            ),
            (
                "method_fingerprint",
                lambda value: value["method_b"].update(fingerprint=None),
                "method_b_fingerprint_missing",
            ),
            (
                "summary",
                lambda value: value["method_b"]["summary"].update(interpolation_used=True),
                "method_b_summary_invalid",
            ),
            (
                "edges",
                lambda value: value.update(canonical_t_edges=[0.0, 0.0]),
                "canonical_geometry_invalid",
            ),
        )
        for name, mutate, reason in cases:
            with self.subTest(name=name):
                input_contract = self._comparison_input()
                mutate(input_contract)
                self._assert_unavailable(input_contract, reason)

    def test_method_b_fingerprints_are_deterministic_and_same_t_projection_is_local(self):
        input_contract = self._comparison_input()
        baseline = comparison.build_pion_hgcer_method_b_comparison(input_contract)
        repeated = comparison.build_pion_hgcer_method_b_comparison(copy.deepcopy(input_contract))
        for key in (
            "source_method_b_payload_fingerprint",
            "fingerprint",
            "parent_region_references",
            "cells",
        ):
            self.assertEqual(repeated[key], baseline[key])
        changed = copy.deepcopy(input_contract)
        changed["method_b"]["cells"][2]["shape_chi2"] = 2.0
        result = comparison.build_pion_hgcer_method_b_comparison(changed)
        self.assertNotEqual(
            result["source_method_b_payload_fingerprint"],
            baseline["source_method_b_payload_fingerprint"],
        )
        self.assertNotEqual(result["fingerprint"], baseline["fingerprint"])
        self.assertEqual(
            [cell for cell in result["cells"] if cell["t_index"] == 0],
            [cell for cell in baseline["cells"] if cell["t_index"] == 0],
        )
        self.assertEqual(
            [parent for parent in result["parent_region_references"] if parent["t_index"] == 0],
            [parent for parent in baseline["parent_region_references"] if parent["t_index"] == 0],
        )


if __name__ == "__main__":
    unittest.main()
