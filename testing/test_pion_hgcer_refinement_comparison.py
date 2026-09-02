"""Pure-Python regression coverage for the frozen Phase-D.1 input contract."""

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


def _method_a(status, host_state):
    available = status == "available"
    return {
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_a_unavailable",
        "fingerprint": "method-a-fingerprint" if available else None,
        "host_state": host_state,
        "cells": _cells("method_a", host_state) if available else [],
    }


def _method_a_summary(status, host_state):
    available = status == "available"
    return {
        "schema_version": "pion_hgcer_method_a/v1",
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_a_unavailable",
        "fingerprint": "method-a-fingerprint" if available else None,
        "phase_a_contract_fingerprint": "phase-a-contract-fingerprint",
        "coordinate_fingerprint": "phase-a-coordinate-fingerprint",
        "t_edges": list(T_EDGES) if available else [],
        "delta_edges": list(DELTA_EDGES) if available else [],
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "host_state": host_state,
    }


def _method_b(status, host_state):
    available = status == "available"
    return {
        "status": status,
        "available": available,
        "reason": None if available else "synthetic_method_b_unavailable",
        "fingerprint": "method-b-fingerprint" if available else None,
        "host_state": host_state,
        "cells": _cells("method_b", host_state) if available else [],
        "parent_region_references": (
            [{"t_index": 0, "region_name": "pi_n"}] if available else []
        ),
    }


def _method_b_summary(status, host_state):
    available = status == "available"
    return {
        "schema_version": "pion_hgcer_method_b/v1",
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
):
    return checkpoint.build_pion_hgcer_refinement_checkpoint(
        setting={
            "kinematic_token": "Q4p4W2p74",
            "Q2": 4.4,
            "W": 2.74,
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
        cases = (
            ("particle_type", "pion"),
            ("epsilon_filename_token", "highe"),
            ("epsilon_setting", ["low"]),
            ("Q2", [4.4]),
            ("W", "2.74"),
        )
        for key, value in cases:
            with self.subTest(key=key):
                payload = _checkpoint()
                payload["setting"][key] = value
                self._assert_unavailable(payload, "setting_contract_invalid", "setting_validation")

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


if __name__ == "__main__":
    unittest.main()
